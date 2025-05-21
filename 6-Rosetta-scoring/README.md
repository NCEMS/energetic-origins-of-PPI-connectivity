### Run Rosetta FastRelax and scoring

* Uses the protein structures from AF2 as inputs to a Rosetta FastRelax protocol with the ref2015 scoring function

* Note that this step take significant time; uses Python multiprocessing to distribute calculations over N CPUs as specified in the .config file for the run.

* With N = 48 CPUs, this step required ~4 days. 
