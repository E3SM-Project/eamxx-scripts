;;; Parse and plot tridiag performance data.

(require [amb [*]])
(import [amb [*]]
        pickle)

(sv dec (Box)
    dec.machines (, 'v100 'knl 'skx)
    dec.name {'skx "SKX" 'knl "KNL" 'v100 "V100"}
    dec.clr {'skx "g" 'knl "b" 'v100 "r"}
    dec.line {'cr "-" 'thomas "-" 'cusparse ":" 'gttr ":" 'dttr ":"
              's-cr "--" 's-thomas "--"}
    dec.marker-seq "ovs*px^.")

(defn base-solver [solver]
  (sv ns (name solver)
      prefix (if (= (cut solver 0 2) "s-") "s-" ""))
  (cond [(in "cr" ns) (symbol (+ prefix "cr"))]
        [(in "thomas" ns) (symbol (+ prefix "thomas"))]
        [:else solver]))

(defn get-field [ln field]
  (sv pos (.find ln field))
  (when (< pos 0) (return -1))
  (sv pos (+ pos (len field) 1)
      end (.find (cut ln pos) " "))
  (cut ln pos (+ pos end)))

(defclass Data []
  (defn --init-- [me]
    (sv me.d {}))

  (defn pickle [me filename]
    (with [f (open filename "w")] (pickle.dump me.d f)))

  (defn unpickle [me filename]
    (sv me.d (with [f (open filename "r")] (pickle.load f))))

  (defn read [me filename machine precision]
    (sv lns (.split (readall filename) "\n")
        impl None)
    (for [ln lns]
      (if (< (len ln) 4) (continue))
      (sv key (cut ln 0 3))
      (cond [(= key "avx")
             (sv avx (cut ln 4)
                 impl 'standalone)]
            [(= key "siz")
             (sv impl 'scream
                 omp (int (last (.split ln)))
                 avx (get-field ln "avx"))]
            [(= key "omp")
             (sv omp (int (cut ln 5)))]
            [(= key "run")
             (cond [(in " solver " ln)
                    (cond [(= impl 'standalone)
                           (sv (, - - solver - nprob - nrow - nrhs - nwarp)
                               (sscanf ln "s,s,s,s,i,s,i,s,i,s,i"))]
                          [(= impl 'scream)
                           (sv (, - - solver - pack - nprob - nrow - nA - nrhs - nwarp)
                               (sscanf ln "s,s,s,s,i,s,i,s,i,s,i,s,i,s,i")
                               solver (+ "s-" solver
                                         (if pack "_pack" "")
                                         (if (> nA 1) "_amxm" "_a1xm")))
                           (assert (in nrow (, 72 128)))
                           (unless (= machine 'v100)
                             (assert (<= nwarp 1))
                             (sv nwarp 1))])]
                   [(in " et " ln)
                    (sv (, - - et - et-per-datum)
                        (sscanf ln "s,s,f,s,f"))
                    (assoc-nested-append
                     me.d
                     (, machine precision (symbol solver) nrow nrhs nprob omp nwarp)
                     (, et et-per-datum))])])))
  
  (defn get-data [me machine solver d0]
    (cond [(= machine 'v100)
           (sv d1 (get d0 1))
           (if (in (base-solver solver) (, 'cr 'thomas 's-cr 's-thomas))
             (do (sv epd-min 1e3
                     nw-best -1)
                 (for [nw (.keys d1)]
                   (sv epd (second (first (get d1 nw))))
                   (when (< epd epd-min)
                     (sv epd-min epd
                         nw-best nw)))
                 (get d1 nw-best))
             (get d1 (first (.keys d1))))]
          [(= machine 'knl)
           (do (sv epd-min 1e3
                   nthr-best -1)
               (for [nthr (.keys d0)]
                 (sv epd (second (first (get d0 nthr 1))))
                 (when (< epd epd-min)
                   (sv epd-min epd
                       nthr-best nthr)))
               (get d0 nthr-best 1))]
          [:else
           (sv d1 (get d0 (first (.keys d0)))
               d1 (get d1 (first (.keys d1))))
           d1]))

  (defn get-y-value [me np unzipped-data probpersec]
    (median (if probpersec
              (/ np (npy.array (first unzipped-data)))
              (second unzipped-data))))

  (defn ylabel [me probpersec]
    (pl.ylabel (if probpersec
                 "Problems per second"
                 "Seconds per RHS datum")))

  (defn plot-baseline [me precision &optional [all False] [probpersec False]
                       [simple False]]
    ;;(sv nprobs (sort (.keys (get me.d 'v100 precision 'cr_a1x1p 128 1))))
    (sv nprobs (sort (.keys (get me.d 'skx precision 'thomas 128 1))))
    (for [nlev (, 72 128 256)]
      (with [(pl-plot (, 7 9) (.format "baseline-prec{}-nlev{}{}"
                                       precision nlev
                                       (cond [all "-all"] [simple "-simple"] [:else ""])))]
        (sv ctr -1)
        (for [machine dec.machines
              solver (cond [all
                            (, 'cusparse 'gttr 'dttr 'cr_a1x1 'cr_a1x1p 'cr_a1xm
                               'cr_amxm 'thomas 'thomas_pack_a1xm
                               's-thomas_a1xm 's-thomas_pack_a1xm 's-cr_a1xm)]
                           [simple (, 'cusparse 'dttr 'cr_a1x1p 'thomas
                                      's-thomas_a1xm 's-cr_a1xm)]
                           [:else
                            (, 'cusparse 'gttr 'dttr 'cr_a1x1 'cr_a1x1p 'thomas
                               's-thomas_a1xm 's-cr_a1xm)])]
          (sv d (geton me.d machine precision solver))
          (when (none? d) (continue))
          (sv x [] y [])
          (for [np nprobs]
            (sv d0 (geton d nlev 1 np))
            (when (none? d0) (continue))
            (sv data (me.get-data machine solver d0))
            (.append x np)
            (.append y (me.get-y-value np (unzip data) probpersec)))
          (pl.loglog
           x y (+ (get dec.clr machine)
                  (get dec.marker-seq (% (inc! ctr) (len dec.marker-seq)))
                  (get dec.line (base-solver solver)))
           :label (+ (get dec.name machine)
                     " " (name solver))
           :lw 2))
        (my-grid)
        (axis-tight-pad :mult True :pad 0.2)
        (pl.xlabel "Number of columns")
        (me.ylabel probpersec)
        (pl.title (+ "Baseline against cusparseTgtsv2StridedBatch\nand MKL gttr, dttr\n"
                     (.format "{} levels, 1 RHS/column" nlev))
                  :fontsize 12)
        (pl.legend :loc "best" :fontsize 10))))

  (defn plot-hommexx [me precision &optional [probpersec False] [simple False]]
    (sv nprobs (sort (.keys (get me.d 'skx precision 'thomas_pack_amxm 128 16))))
    (for [nlev (, 72 128 256)
          nrhs (, 10 16)]
      (with [(pl-plot (, 7 9) (.format "hommexx-prec{}-nlev{}-nrhs{}{}"
                                       precision nlev nrhs
                                       (cond [simple "-simple"] [:else ""])))]
        (sv ctr -1)
        (for [machine dec.machines
              solver (cond [simple
                            (, 'cusparse 'cr_amxm
                               'dttr
                               'thomas_pack_amxm
                               's-thomas_amxm 's-thomas_pack_amxm 's-cr_amxm)]
                           [:else
                            (, 'cusparse 'cr_amxm 'thomas
                               'gttr 'dttr
                               'thomas_amxm 'thomas_pack_amxm
                               's-thomas_amxm 's-thomas_pack_amxm 's-cr_amxm)])]
          (when (and (!= machine 'v100) (= solver 'thomas)) (continue))
          (sv d (geton me.d machine precision solver))
          (when (none? d) (continue))
          (sv x [] y []
              mult (in solver (, 'gttr 'dttr)))
          (for [np nprobs]
            (sv nrhs-key (if mult 1 nrhs)
                np-key (if mult (* nrhs np) np)
                d0 (geton d nlev nrhs-key np-key))
            (when (none? d0) (continue))
            (sv data (me.get-data machine solver d0))
            (.append x np)
            (.append y (me.get-y-value np (unzip data) probpersec)))
          (pl.loglog
           x y (+ (get dec.clr machine)
                  (get dec.marker-seq (% (inc! ctr) (len dec.marker-seq)))
                  (get dec.line (base-solver solver)))
           :label (+ (get dec.name machine)
                     " " (.replace (name solver) "-amxm" ""))
           :lw 2))
        (my-grid)
        (axis-tight-pad :mult True :pad 0.2)
        (pl.xlabel "Number of elements")
        (me.ylabel probpersec)
        (pl.title (+ "HOMME IMEX Problem:\n"
                     (.format "{} levels, {} L&RHS/column" nlev nrhs))
                  :fontsize 12)
        (pl.legend :loc "best" :fontsize 10))))

  (defn plot-shoc [me precision &optional [probpersec False] [simple False]]
    (sv nprobs (sort (.keys (get me.d 'skx precision 'thomas 128 43))))
    (for [nlev (, 72 128 256)
          nrhs (, 2 13 43)]
      (with [(pl-plot (, 7 9) (.format "shoc-prec{}-nlev{}-nrhs{}{}"
                                       precision nlev nrhs
                                       (cond [simple "-simple"] [:else ""])))]
        (sv ctr -1)
        (for [machine dec.machines
              solver (cond [simple
                            (, 'cusparse 'cr_a1xm 'thomas
                               'dttr
                               'thomas_pack_a1xm
                               's-thomas_a1xm 's-thomas_pack_a1xm 's-cr_a1xm)]
                           [:else
                            (, 'cusparse 'cr_a1xm 'thomas
                               'gttr 'dttr
                               'thomas_pack_a1xm
                               's-thomas_a1xm 's-thomas_pack_a1xm 's 's-cr_a1xm)])]
          (when (and simple (!= machine 'v100) (= solver 'thomas)) (continue))
          (sv d (geton me.d machine precision solver))
          (when (none? d) (continue))
          (sv x [] y [])
          (for [np nprobs]
            (sv d0 (geton d nlev nrhs np))
            (when (none? d0) (continue))
            (sv data (me.get-data machine solver d0))
            (.append x np)
            (.append y (me.get-y-value np (unzip data) probpersec)))
          (pl.loglog
           x y (+ (get dec.clr machine)
                  (get dec.marker-seq (% (inc! ctr) (len dec.marker-seq)))
                  (get dec.line (base-solver solver)))
           :label (+ (get dec.name machine)
                     " " (.replace (name solver) "-a1xm" ""))
           :lw 2))
        (my-grid)
        (axis-tight-pad :mult True :pad 0.2)
        (pl.xlabel "Number of columns")
        (me.ylabel probpersec)
        (pl.title (+ "SHOC Problem:\n"
                     (.format "{} levels, {} RHS/column" nlev nrhs))
                  :fontsize 12)
        (pl.legend :loc "best" :fontsize 10)))))

;;; Data parsing and plots for the report.

(when-inp ["read"]
  (sv d (Data)
      set1 (, "cdata-skx-dp-1.txt" "cdata-knl-dp-1.txt" "cdata-v100-dp-2.txt")
      set2 (, "cdata-skx-dp-2.txt" "cdata-knl-dp-2.txt" "cdata-v100-dp-3.txt")
      set3 (, "cdata-skx-dp-3.txt"
              "cdata-knl-dp-3.txt" ; 136 threads
              "cdata-knl-dp-5.txt" ; 272 threads (4 was killed in power outage)
              "cdata-v100-dp-3.txt")
      set4 (, "cdata-skx-sp-0.txt"
              "cdata-knl-sp-0.txt" ; 272
              "cdata-knl-sp-1.txt" ; 136
              "cdata-v100-sp-0.txt"))
  (for [filename set3]
    (sv tokens (.split filename "-")
        machine (symbol (second tokens))
        precision (if (= (get tokens 2) "dp") 2 1))
    (d.read filename machine precision))
  (d.pickle "tdp-read.pickle"))

(when-inp ["plot"]
  (sv d (Data))
  (d.unpickle "tdp-read.pickle")
  (sv pps True prec 1)
  (d.plot-baseline prec :probpersec pps)
  (d.plot-baseline prec :probpersec pps :simple True)
  (d.plot-baseline prec :probpersec pps :all True)
  (d.plot-shoc prec :probpersec pps)
  (d.plot-shoc prec :probpersec pps :simple True)
  (d.plot-hommexx prec :probpersec pps)
  (d.plot-hommexx prec :probpersec pps :simple True))

;;; Data parsing and plots to check performance of scream's integrated tridiag
;;; solvers against scream-doc's standalone ones.

(when-inp ["perf-check-read"]
  (sv d (Data)
      set1 (, "v100-1.txt" "skx-3.txt" "knl-0.txt" "knl-1.txt")
      set2 (, "v100-2.txt" "skx-6.txt" "knl-2.txt" "knl-3.txt")
      precision 2)
  (for [filename set2]
    (sv tokens (.split filename "-")
        machine (symbol (first tokens)))
    (d.read filename machine precision))
  (d.pickle "tdp-perf-check.pickle"))

(when-inp ["perf-check-plot"]
  (sv d (Data))
  (d.unpickle "tdp-perf-check.pickle")
  (sv pps True prec 2)
  (d.plot-hommexx prec :probpersec pps)
  (d.plot-baseline prec :probpersec pps :all True)
  (d.plot-shoc prec :probpersec pps))
