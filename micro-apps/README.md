# micro-apps
==================================
This directory contains code for each micro-app we've generated before it gets migrated to
the (scream fork of the) E3SM repo.

## Directory structure as of 1/11/19:
==================================
* micro-sed/:   This contains the infrastructure for running the micro-sed
	        micro-app in standalone mode
* p3/:          Contains the micro-app code that will be eventually used in E3SM's
	        p3 parameterization
* share/:	Contains code likely to be needed by multiple micro-apps. This
	        also happens to be the most technical bits of code.
* nano/:        This was the precursor to the micro-sed micro-app. Old/unused now.
* perf-results: where perf data is archived whenever standalone micro-apps are run.
