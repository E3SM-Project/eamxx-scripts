import sys, argparse, os, textwrap, time, glob, subprocess, datetime

TESTING = False # If true, use test.txt to act as squeue.

class Struct:
    pass

class JobState:
    none = 'NONE'
    pend = 'PEND'
    run = 'RUN'
    completing = 'CG'
    done = 'DONE'
    exited = 'EXIT'
    unknown = 'UNKNOWN'

def get_file_time_delta(fname):
    clock_time = time.clock_gettime(0)
    file_time = None
    try:
        file_time = os.path.getmtime(fname)
    except:
        print('Could not get timestamp of {}'.format(fname))
    if file_time is None:
        # The likeliest issue is the file hasn't been created yet or the network
        # file system temporarily has lost access. Neither of these is
        # associated with a hang condition, so return 0.
        return 0
    else:
        return clock_time - file_time

class E3smLogReader:
    nofile = 'NOFILE'
    ok = 'OK'
    failed = 'FAILED'
    done = 'DONE'
    unknown = 'UNKNOWN'

    def __init__(self, rundir, jobid, atm_log_timeout=7200):
        self.rundir = rundir
        self.jobid = jobid
        self.seek = 0 # At each call to study(), read only new text.
        self.e3sm_log_fname = None
        self.atm_log_fname = None
        self.atm_log_timeout = atm_log_timeout
        print(('Initializing E3smLogReader with run directory {}\n' +
               '  job ID {}\n' +
               '  atm.log timeout {}').format(
                   self.rundir, self.jobid, self.atm_log_timeout))

    def study(self):
        if self.e3sm_log_fname is None:
            # Construct the e3sm.log file name.
            fnames = glob.glob(self.rundir + '/' + 'e3sm.log.' + str(self.jobid) + '*')
            if len(fnames) > 1:
                print('Run dir {} and job id {} yield multiple files: {}'.
                      format(self.rundir, self.jobid, fnames))
                return self.unknown, ''
            elif len(fnames) == 0:
                return self.nofile, ''
            else:
                self.e3sm_log_fname = fnames[0]
                self.atm_log_fname = self.e3sm_log_fname.replace('e3sm.log.', 'atm.log.')
                print('Found {}'.format(self.e3sm_log_fname))
        # Before proceeding, check if the e3sm.log file has been gzipped.
        try:
            os.stat(self.e3sm_log_fname + '.gz')
            # It has. The job is nearly done.
            return self.done, ''
        except:
            # It has not.
            pass
        # Main work:
        # 1. Check for strings indicating a failed job.
        try:
            with open(self.e3sm_log_fname, 'r') as f:
                f.seek(self.seek)
                for ln in f:
                    self.seek += len(ln)
                    bad_text = find_bad_text(ln)
                    if bad_text is not None:
                        return self.failed, bad_text
        except:
            print('Could not read {}'.format(self.e3sm_log_fname))
            return self.unknown, ''
        # 2. Check that atm.log has been written sufficiently recently.
        dt = get_file_time_delta(self.atm_log_fname)
        if dt > self.atm_log_timeout:
            print('{} has not been written for {} seconds but timeout is {}'.
                  format(self.atm_log_fname, dt, self.atm_log_timeout))
            return self.failed, 'atm.log is not being written'
        # Everything is fine.
        return self.ok, ''

def find_bad_text(ln):
    "Check for text indicating a failed run in 'ln'."
    bad_strs = ['MPI_ABORT was invoked', 'stack frames', 'FATAL ERROR: Aborting',
                'Error configuring interconnect',
                'ERROR: ld.so', 'cannot open shared object file',
                'error while loading shared libraries',
                'check_dim ERROR: mismatch of input dimension',
                'ERROR:  One or more process',
                'Program received signal SIGABRT',
                'terminated with signal', '(core dumped)', 'Backtrace for this error:',
                'An error occurred in MPI_Init',
                'MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort',
                'MPIDI_OFI_handle_cq_error',
                'Segmentation fault',
                'Aborting since the error handler was set to PIO_INTERNAL_ERROR']
    for b in bad_strs:
        if b in ln:
            return b
    return None

def sp_run(cmd_list):
    "Wrapper to subprocess.run."
    return subprocess.run(cmd_list, capture_output=True)
    
def call_job_lister(jobid):
    "Call squeue."
    if TESTING:
        sp = Struct()
        with open('test.txt', 'r') as f:
            sp.stdout = f.read()
            sp.stderr = sp.stdout
    else:
        sp = sp_run(['squeue', '--job', str(jobid)])
    if 'Invalid job id specified' in str(sp.stderr): return JobState.none
    so = str(sp.stdout)
    if ' PD ' in so: return JobState.pend
    if ' R ' in so: return JobState.run
    if ' CG ' in so: return JobState.completing
    return JobState.unknown

def call_job_killer(jobid):
    "Call scancel."
    # Run a sequence of kill commands to be maximally convincing.
    for cmd in (['scancel', str(jobid)],):
        sp = sp_run(cmd)
        print('{}: {} {}'.format(' '.join(cmd), str(sp.stdout, 'utf-8'),
                                 str(sp.stderr, 'utf-8')))

def monitor(rundir, jobid, sleep=300, atm_log_timeout=7200):
    "Run the job monitoring loop."
    print(f'Monitoring {rundir:s} {jobid:s} with sleep {sleep:d}')
    kill_job = False
    elr = E3smLogReader(rundir, jobid, atm_log_timeout=atm_log_timeout)
    while True:
        job_state = call_job_lister(jobid)
        print('{} {} {}'.format(str(datetime.datetime.now())[0:-7], jobid, job_state))
        if not kill_job and job_state == JobState.run:
            e3sm_log_state, evidence = elr.study()
            if e3sm_log_state == E3smLogReader.failed:
                print(f'Job {jobid:s} failed; evidence: "{evidence:s}"')
                kill_job = True
        elif job_state in (JobState.none, JobState.unknown):
            print(f'Job {jobid:s} is not in the list of jobs.')
            break
        elif job_state in (JobState.done, JobState.exited):
            print('Job {} is in completed state {}.'.format(jobid, job_state))
            break
        if kill_job:
            print(f'Attempting to kill job {jobid:s}')
            if TESTING: break
            call_job_killer(jobid)
        sys.stdout.flush()
        time.sleep(sleep)

def main():
    help = textwrap.dedent("""
    Monitor an E3SM CIME-configured job. The goal is to kill jobs if E3SM has
    aborted rather than have jobs hang with nothing happening.

    This script works as follows. In a sleep loop, check the e3sm.log file for
    evidence of a run that has aborted. If such evidence is found, issue
    commands to kill the job. Once the job is no longer visible in the queue,
    the script exits.

    Run this script on a login node. The sleep command is extremely efficient,
    so the script has essentially no impact on other users on the login node. To
    use this script effectively, you should prefix it with 'nohup' so it
    continues to run if your connection gets interrupted. Alternatively, run
    your terminal inside a 'screen' session, which also is persistent when a
    connection is interrupted.

    Example: Suppose you want to monitor two jobs. In a shell, issue these
    commands in sequence:
       nohup python jobmonitor_frontier.py ${path_to_run_dir1} ${jobid1} > job.1.txt &
       nohup python jobmonitor_frontier.py ${path_to_run_dir2} ${jobid2} > job.2.txt &
    These jobs run in the background, and you can look at job.1.txt and
    job.2.txt for messages about what the job monitor has done so far.
    ${path_to_run_dir1} is something like
       /gpfs/alpine/cli115/proj-shared/username/e3sm_scratch/ne4pg2_ne4pg2.F2010-SCREAMv1/run
    and ${jobid} is printed at ./case.submit and is also available from
       squeue -u username
    """)

    p = argparse.ArgumentParser(description=help, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('rundir', type=str,
                   help='Path to directory containing (future) e3sm.log file.')
    p.add_argument('jobid', type=str,
                   help='Job ID: slurm job ID field.')
    p.add_argument('--sleep', type=int, default='60',
                   help='Time length in seconds to sleep between checks.')
    p.add_argument('--atmlog-timeout', type=int, default='7200',
                   help=textwrap.dedent("""\
                   Maximum permitted time, in seconds, between atm.log writes;
                   if exceeded, assume the job has hung.
                   """))
    p.add_argument('--test', default=False, action='store_true',
                   help=textwrap.dedent("""\
                   Test mode. Set up an e3sm.log.${jobid}.* file and test.txt.
                   Put slurm-like messages in test.txt for jobmonitor_frontier.py
                   to read.
                   """))
    o = p.parse_args()

    global TESTING
    TESTING = o.test

    monitor(o.rundir, o.jobid, sleep=o.sleep, atm_log_timeout=o.atmlog_timeout)

if __name__ == '__main__':
    main()
