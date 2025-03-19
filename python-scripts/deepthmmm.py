import biolib

deeptmhmm = biolib.load('DTU/DeepTMHMM')

biolib.utils.STREAM_STDOUT = True
deeptmhmm_job = deeptmhmm.cli(args='--fasta query.fasta')
deeptmhmm_job.save_files('test_deepthmmm') # Saves all results to `result` dir
