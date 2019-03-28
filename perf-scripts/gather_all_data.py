from utils import run_cmd_no_fail

import os
import concurrent.futures as threading3

# MACHINE -> (env_setup, repo-loc, kokkos-install-loc, compiler, batch submit prefix)
MACHINE_METADATA = {
    "melvin"   : (["module purge && module load sems-env && module load sems-gcc/7.3.0 sems-openmpi/1.10.1 sems-gcc/7.3.0 sems-git/2.10.1 sems-cmake/3.10.3 sems-python/3.5.2"], "$(which mpicxx)", ""),
    "bowman"   : (["module load openmpi/1.10.6/intel/17.2.174 git/2.8.2 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-bowman/bin:$PATH"], "$(which mpicxx)", "srun"),
    "blake"    : (["module load openmpi/2.1.5/intel/19.1.144 git/2.9.4 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-blake/bin:$PATH"],  "$(which mpicxx)", "srun"),
    "waterman" : (["module load devpack/latest/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 git/2.10.1", "module switch cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-waterman/bin:$PATH"],
                  "$(which mpicxx)", "bsub -I -q rhel7W"),
    "white"    : (["module load devpack/20181011/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 git/2.10.1 cmake/3.12.3", "export PATH=/ascldap/users/jgfouca/packages/Python-3.6.8-white/bin:$PATH"], "$(which mpicxx)", "bsub -I -q rhel7G"),
}

###############################################################################
class GatherAllData(object):
###############################################################################

    ###########################################################################
    def __init__(self, perf_analysis_args, commit, machines, test_all, scream, local, kokkos, submit):
    ###########################################################################
        self._perf_analysis_args = perf_analysis_args
        self._commit             = commit
        self._machines           = machines
        self._test_all           = test_all
        self._scream             = scream
        self._local              = local
        self._kokkos             = kokkos
        self._submit             = submit

    ###########################################################################
    def formulate_command(self, machine):
    ###########################################################################
        env_setup, compiler, batch = MACHINE_METADATA[machine]

        if self._local:
            scream_docs_repo = os.path.abspath("./scream-docs/micro-apps")
            repo             = os.path.abspath("./scream/components/scream") if self._scream else scream_docs_repo
        else:
            scream_docs_repo = "~/scream-docs-perf-{}/micro-apps".format(machine)
            repo             = "~/scream-perf-{}/components/scream".format(machine) if self._scream else scream_docs_repo

        test_all_script = "../perf-scripts/test-all-scream" if self._scream else "test-all"
        submit_args = "-s" if self._submit else ""

        kokkos_arg = "-k {}".format(self._kokkos) if self._kokkos else ""
        kokkos_loc = self._kokkos if self._kokkos else os.path.join(os.path.dirname(os.path.dirname(repo)), "externals", "kokkos")

        local_cmd = "{}/{} {} {} {} -m {}".format(scream_docs_repo, test_all_script, compiler, kokkos_arg, submit_args, machine) if self._test_all else "../perf-scripts/perf_analysis {} -p".format(self._perf_analysis_args)
        local_cmd = local_cmd.replace("$compiler", compiler)

        setup = "cd {} && git fetch && git reset --hard origin/master && ".format(scream_docs_repo) if (self._scream and not self._local) else ""
        extra_env = ""
        if machine in ["waterman", "white"]:
            extra_env = "OMPI_CXX={}/bin/nvcc_wrapper ".format(kokkos_loc)
        else:
            extra_env = "OMP_PROC_BIND=FALSE "

        repo_setup = "true" if (self._local or not self._commit) else "git fetch && git checkout {} && git submodule update --init".format(self._commit)

        cmd = "{}cd {} && {} && {} && {}{} {}".format(setup, repo, " && ".join(env_setup), repo_setup, extra_env, batch, local_cmd)

        return cmd

    ###########################################################################
    def run_on_machine(self, machine):
    ###########################################################################
        cmd = self.formulate_command(machine)
        print("Starting {} analysis on {} with cmd: {}".format("test-all" if self._test_all else "performance", machine, cmd))

        if self._local:
            run_cmd_no_fail(cmd, arg_stdout=None, arg_stderr=None, verbose=True, exc_type=RuntimeError)
        else:
            try:
                output = run_cmd_no_fail("ssh -o StrictHostKeyChecking=no {} '{}'".format(machine, cmd), exc_type=RuntimeError, combine_output=True)
            except RuntimeError as e:
                output = str(e)
                raise
            finally:
                with open(os.path.join("test-all-results" if self._test_all else "perf-results", self._commit, machine), "w") as fd:
                    fd.write(output)

        print ("Completed {} analysis on {}".format("test-all" if self._test_all else "performance", machine))

    ###########################################################################
    def gather_all_data(self):
    ###########################################################################
        if not self._local:
            os.makedirs(os.path.join("test-all-results" if self._test_all else "perf-results", self._commit))

        success = True

        with threading3.ThreadPoolExecutor(max_workers=len(self._machines)) as executor:
            future_to_machine = {executor.submit(self.run_on_machine, machine): machine for machine in self._machines}
            for future in threading3.as_completed(future_to_machine):
                machine = future_to_machine[future]
                try:
                    future.result()
                except RuntimeError:
                    print('{} failed'.format(machine))
                    success = False

        return success
