1. ln -s ${SCREAMROOT}/components/cam/src/physics/cam/shoc.F90 .

2. make
2a. Optionally, edit Makefile to build debug.

3. mkdir fig

4. Run with a scientific python 3, such as anaconda:
4a. Short run just to see that things are working. Output: fig/ics.pdf and fig/fin.pdf
      python shocinter.py
4b. Long run, convergence test. Output: fig/conv.pdf
      python shocintr.py conv
