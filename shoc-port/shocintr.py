import sys, os, math, ctypes, numpy as npy

import matplotlib
#matplotlib.use('AGG')
import matplotlib.pyplot as pl

def is_coll(e):
    try:
        iter(e)
    except:
        return False
    return True

def s_all_rev(): return slice(None, None, -1)
def cc(x): return 0.5*(x[:-1] + x[1:])

class Struct:
    pass

def dispfig(fn_prefix=None, format='pdf', tight=True):
    if tight: pl.tight_layout()
    if not fn_prefix or len(fn_prefix) == 0:
        return pl.show()
    else:
        pl.savefig(fn_prefix + '.' + format, format=format, bbox_inches='tight')

def my_grid():
    pl.grid(True, lw=0.5, ls='-', color=(0.8, 0.8, 0.8), zorder=-1, which='both')
    return pl.gca().set_axisbelow(True)

def pad_lim(lim, pad=0.05, mult=False):
    if mult:
        v = lim[0] * (1 - pad), lim[1] * (1 + pad)
    else:
        d = lim[1] - lim[0]
        delta = pad * d
        v = lim[0] - delta, lim[1] + delta
    return v

def axis_tight_pad(pad=0.05, mult=False):
    pl.axis('tight')
    xl = pl.xlim()
    yl = pl.ylim()
    pl.xlim(pad_lim(xl, pad, mult))
    return pl.ylim(pad_lim(yl, pad, mult))

class pl_plot:
    def __init__(me, figsize, filename, format=None, tight=True):
        me.filename = filename
        me.format = 'pdf' if not None else format
        me.tight = tight
        pl.close()
        pl.figure(num=1, figsize=figsize)
    def cleanup(me):
        dispfig(me.filename, format=me.format, tight=me.tight)
    def __enter__(me): return me
    def __exit__(me, *args): pass
    def __del__(me): return me.cleanup()

# Direct pointer from the numpy array.
as_ctypes = npy.ctypeslib.as_ctypes

# Make a ctypes array.
def ctype_vec(B):
    if  B.dtype == npy.dtype("float64"): f = ctypes.c_double*B.size
    elif B.dtype == npy.dtype("int"): f = ctypes.int*B.size
    else: raise('ctype_vec: not a valid dtype: {}'.format(B.dtype))
    return f(*list(B.reshape(B.size)))

def vec(a): return a.reshape(a.size)

class PhysConsts:
    gravit = 9.80616
    rair   = 287.042
    rh2o   = 461.505
    cpair  = 1004.64
    latvap = 2501000.0
    p0     = 1e5

# Class to call Fortran SHOC.
class Shoc:
    def __init__(me, shcol, nlev, nqtracers):
        me.lib = me.load_lib()
        if not me.lib: return
        me.shoc_init()
        me.set_defaults(shcol, nlev, nqtracers)
        
    def load_lib(me):
        try:
            lib = npy.ctypeslib.load_library("libshoc", ".")
        except Exception as e:
            print(e)
            lib = None
        return lib        

    def shoc_init(me):
        c = PhysConsts
        d = ctypes.c_double
        me.lib.shoc_c_init(d(c.gravit), d(c.rair), d(c.rh2o), d(c.cpair), d(c.latvap))

    def set_defaults(me, shcol, nlev, nqtracers):
        me.shcol = shcol
        me.nlev = nlev
        me.nqtracers = nqtracers
        me.dtime = 300

        aset = npy.zeros
        zc = lambda: aset((me.shcol))
        zcm = lambda: aset((me.shcol, me.nlev))
        zci = lambda: aset((me.shcol, me.nlev+1))
        zcq = lambda: aset((me.shcol, me.nqtracers))
        zcmq = lambda: aset((me.shcol, me.nlev, me.nqtracers))
        def zero(zeroer, names):
            for e in names:
                me.__dict__[e] = zeroer()

        zero(zc, ['host_dx', 'host_dy', 'wthl_sfc', 'wqw_sfc', 'uw_sfc', 'vw_sfc'])
        zero(zcm, ['zt_grid', 'pres', 'pdel', 'thv', 'cldliq', 'w_field', 'tke', 'thetal',
                   'qw', 'u_wind', 'v_wind', 'wthv_sec', 'tk', 'tkh', 'shoc_cldfrac',
                   'shoc_ql', 'shoc_mix', 'w_sec', 'wqls_sec', 'brunt', 'isotropy'])
        zero(zci, ['zi_grid', 'thl_sec', 'qw_sec', 'qwthl_sec', 'wthl_sec', 'wqw_sec',
                   'wtke_sec', 'uw_sec', 'vw_sec', 'w3'])
        zero(zcq, ['wtracer_sfc'])
        zero(zcmq, ['qtracers'])

    def test_lib(me):
        c = ctypes
        i = c.c_int
        d = c.c_double
        v = lambda e: ctype_vec(npy.transpose(e))
        arr = npy.zeros((me.shcol, me.nlev))
        varr = v(arr)
        me.lib.hi(i(me.shcol), i(me.nlev), i(me.nlev+1), i(me.nqtracers),
                  d(me.dtime), varr)
        arr = npy.transpose(npy.array(varr).reshape(arr.shape[s_all_rev()]))
        print(arr)

    def call(me, ntimes=-1, verbose=False):
        def flipvert(e):
            # Flip vertical after recent SHOC index flip.
            if len(e.shape) >= 2 and e.shape[1] >= me.nlev:
                if len(e.shape) == 2: e = e[:,s_all_rev()]
                else: e = e[:,s_all_rev(),:]
            return e
        def v(e):
            return ctype_vec(npy.transpose(flipvert(e)))
        di = me.__dict__
        c = ctypes
        i = c.c_int
        d = c.c_double
        inout = ['tke', 'thetal', 'qw', 'u_wind', 'v_wind', 'qtracers', 'wthv_sec', 'tkh',
                 'tk', 'shoc_cldfrac', 'shoc_ql', 'shoc_mix', 'isotropy', 'w_sec', 'thl_sec',
                 'qw_sec', 'qwthl_sec', 'wthl_sec', 'wqw_sec', 'wtke_sec', 'uw_sec', 'vw_sec',
                 'w3', 'wqls_sec', 'brunt']
        v_inout = [v(di[e]) for e in inout]
        shoc_main = me.lib.shoc_c_main_ntimes if ntimes > 0 else me.lib.shoc_c_main
        if ntimes > 0: v_inout.append(i(ntimes))
        shoc_main(
            #in
            i(me.shcol), i(me.nlev), i((me.nlev+1),), d(me.dtime),
            v(me.host_dx), v(me.host_dy), v(me.thv), v(me.cldliq), v(me.zt_grid), v(me.zi_grid),
            v(me.pres), v(me.pdel), v(me.wthl_sfc), v(me.wqw_sfc), v(me.uw_sfc), v(me.vw_sfc),
            v(me.wtracer_sfc), i(me.nqtracers), v(me.w_field),
            #inout, out
            *v_inout)
        for i, e in enumerate(v_inout):
            if ntimes > 0 and i == len(v_inout)-1: break
            name = str(inout[i])
            if verbose: print(name)
            f = di[name]
            if verbose:
                print(f.shape)
                print(len(v_inout[i]))
                print(npy.array(v_inout[i]).shape)
            di[name] = flipvert(npy.transpose(npy.array(v_inout[i]).reshape(f.shape[s_all_rev()])))

    def get_state(me):
        state = ['thv', 'cldliq', 'zt_grid', 'zi_grid', 'pres', 'pdel',
                 'wthl_sfc', 'wqw_sfc', 'uw_sfc', 'vw_sfc', 'wtracer_sfc', 'w_field',
                 'tke', 'thetal', 'qw', 'u_wind', 'v_wind', 'qtracers', 'wthv_sec', 'tkh',
                 'tk', 'shoc_cldfrac', 'shoc_ql', 'shoc_mix', 'isotropy', 'w_sec', 'thl_sec',
                 'qw_sec', 'qwthl_sec', 'wthl_sec', 'wqw_sec', 'wtke_sec', 'uw_sec', 'vw_sec',
                 'w3', 'wqls_sec', 'brunt']
        d = Struct()
        for s in state:
            d.__dict__[s] = me.__dict__[s]
        return d

## Some case utilities.

def calc_hydrostatic_p(ps, th0, z, theta):
    assert len(z) == len(theta)
    pc = PhysConsts
    k = pc.rair / pc.cpair
    c = -pc.gravit * pc.p0**k / pc.rair
    p = npy.zeros(len(z))
    for i in range(len(z)):
        z0 = 0.0 if i == 0 else z[i-1]
        z1 = z[i]
        th0 = th0 if i == 0 else theta[i-1]
        th1 = theta[i]
        p0 = ps if i == 0 else p[i-1]
        if abs(th0 - th1) < 1e-14 * th0:
            p[i] = (p0**k + k*c*(z1 - z0)/th0)**(1/k)
        else:
            ra = (z1 - z0)/(th1 - th0)
            p[i] = (p0**k + k*c*ra*math.log(th1/th0))**(1/k)
    return p

def p_to_exner(p):
    return (p / PhysConsts.p0)**(PhysConsts.rair / PhysConsts.cpair)

def theta_to_T(theta, p):
    return p_to_exner(p)*theta

def theta_to_thetav(theta, qv):
    zvir = PhysConsts.rh2o / PhysConsts.rair - 1
    return theta * (1 + zvir*qv)

def linterp(x, y, xi):
    assert xi[0] >= x[0] and xi[-1] <= x[-1]
    def _linterp(i, j):
        a = (xi[j] - x[i]) / (x[i+1] - x[i])
        return (1 - a) * y[i] + a * y[i+1]
    yi = npy.zeros(xi.shape)
    j = 0
    for i in range(len(x)-1):
        while j < len(xi) and xi[j] <= x[i+1]:
            yi[j] = _linterp(i, j)
            j = j+1
    return yi

def set_uni_mesh(s, c, ztop):
    s.zi_grid[c] = npy.linspace(0, ztop, s.nlev+1)
    s.zt_grid[c] = cc(s.zi_grid[c])

# A 'case': particular set of initial conditions, surface fluxes, and
# forcings. Copy this class and fill in the functions with different numbers to
# make other cases.
#   WARNING: This case does *not* necessarily exercise every part of SHOC.
class ExampleCase:
    def get_zref(me): return npy.array([0., 520, 1480, 2000, 3000])

    def get_ps(me): return 1015e2

    def get_ic_qv(me, z):
        zref = me.get_zref()
        qvref = 1e-3*npy.array([17., 16.3, 10.7, 4.2, 3])
        return linterp(zref, qvref, z)

    def get_ic_theta(me, z):
        zref = me.get_zref()
        thrf = npy.array([299.7, 298.7, 302.4, 308.2, 312.85])
        return linterp(zref, thrf, z)

    def get_ls_winds(me, z):
        zref = npy.array([0., 700, 3000])
        uref = npy.array([-7.75, -8.75, -4.61])
        vref = npy.array([0., 0, 0])
        wref = npy.array([0., 0, 0])
        return linterp(zref, uref, z), linterp(zref, vref, z), linterp(zref, wref, z)

    def set_theta_pressure(me, s, c, zi, zt):
        theta_zi = me.get_ic_theta(zi)
        theta_zt = me.get_ic_theta(zt)
        ps = me.get_ps()
        s.pdel[c] = npy.abs(npy.diff(calc_hydrostatic_p(ps, theta_zi[0], zi, theta_zi)))
        s.pres[c] = calc_hydrostatic_p(ps, theta_zi[0], zt, theta_zt)
        s.qw[c] = me.get_ic_qv(zt)
        s.thv[c] = theta_to_thetav(theta_zt, s.qw[c])
        s.thetal[c] = theta_zt

    def set_ls_winds(me, s, c, zt):
        u, v, w = me.get_ls_winds(zt)
        s.u_wind[c] = u
        s.v_wind[c] = v
        s.w_field[c] = w

    def set_surface_fluxes(me, s, c):
        s.wthl_sfc[c] = 0.008
        s.wqw_sfc[c] = 5.2e-05
        ustar2 = 0.28**2
        u = s.u_wind[c][0]
        v = s.v_wind[c][0]
        speed = math.sqrt(u*u + v*v)
        s.uw_sfc[c] = -ustar2 * (u / speed)
        s.vw_sfc[c] = -ustar2 * (v / speed)

    def set_forcings(me, s, c): pass

    def set_tke(me, s, c, zt): pass

    def set_tracers(me, s, c, zt):
        for q in range(s.qtracers.shape[2]):
            s.qtracers[c,:,q] = npy.sin(q/3. + 1e-3*zt*(0.2*(c+1)))
            s.wtracer_sfc[c,q] = 0.1*c

    def set_time_dep(me, s):
        for c in range(s.shcol):
            me.set_surface_fluxes(s, c)
            me.set_forcings(s, c)

# Given a case object and number of levels, initialize it with initial
# conditions.
def get_ics(case, nz):
    tc = case
    s = Shoc(1,nz,1)
    s.dtime = 12
    s.host_dx[:] = 5300
    s.host_dy[:] = s.host_dx[0]
    for c in range(s.shcol):
        set_uni_mesh(s, c, 2400)
        zt = s.zt_grid[c]
        tc.set_theta_pressure(s, c, s.zi_grid[c], zt)
        tc.set_ls_winds(s, c, zt)
        tc.set_surface_fluxes(s, c)
        tc.set_tke(s, c, zt)
        tc.set_tracers(s, c, zt)
    tc.set_time_dep(s)
    return s

# Given a Shoc object, a Shoc::get_state object, or a collection of these, make
# a grid of figures showing the basic quantities.
def plot_basics(ss, filename):
    axs = []
    def _plot_basics(s, first):
        plotno = [1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
        grids  = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  1,  1,  1,  0,  1,  1 , 1,  0,  0,  0 ]
        zs = s.zt_grid[0], s.zi_grid[0]
        fields = ['pres', 'pdel', 'thv', 'thetal', 'cldliq', 'shoc_ql',
                  'qw', 'u_wind', 'v_wind', 'w_field', 'tke', 'tkh',
                  'shoc_mix', 'isotropy', 'wtke_sec', 'uw_sec', 'wthl_sec',
                  'wqw_sec', 'w_sec', 'thl_sec', 'qw_sec',
                  'w3', 'wqls_sec', 'brunt', 'qtracers']
        for i in range(len(plotno)):
            if i == 0 or plotno[i] != plotno[i-1]:
                if first: axs.append(pl.subplot(5, 4, plotno[i]))
                else: pl.sca(axs[plotno[i]-1])
            name = fields[i]
            if name == 'qtracers': y = s.qtracers[0,:,0]
            else: y = s.__dict__[name]
            pl.plot(vec(y), 0.001 * zs[grids[i]], '-',
                    label=name if first else None)
            if first: pl.legend(loc='best', fontsize=10)
            axis_tight_pad()
    
    if not is_coll(ss): ss = [ss]
    with pl_plot((16, 20), filename):
        for i, s in enumerate(ss):
            _plot_basics(s, i == 0)

# Given a case, set it up and run it. This function has hardcoded things in it;
# it's just meant to show how to do these steps.
def example_run_case(tc):
    nz = 160
    s = get_ics(tc, nz)
    plot_basics(s, "fig/ics")
    s.dtime = 10
    # Advance in time.
    for ti in range(1,101):
        tc.set_time_dep(s)
        s.call()
    # Save the current state.
    s1 = s.get_state()
    # Then advance further.
    for ti in range(1,21):
        tc.set_time_dep(s)
        s.call()
    # Plot both snapshots.
    plot_basics([s1, s], "fig/fin")

def devtest():
    c = PhysConsts
    print(c.gravit)
    s = Shoc(3,72,10)
    s.test_lib()
    
if __name__ == '__main__':
    example_run_case(ExampleCase())
