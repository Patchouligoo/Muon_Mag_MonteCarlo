import pycondor, argparse, sys, os.path
from glob import glob
import numpy as np
import pandas as pd

error = '/scratch/runzeli/condor/error'
output = '/scratch/runzeli/condor/output'
log = '/scratch/runzeli/condor/log'
submit = '/scratch/runzeli/condor/submit'

job = pycondor.Job('Runze_muons','./F_B_test.py',
                        error=error,
                        output=output,
                        log=log,
                        submit=submit,
            getenv=True,
            universe='vanilla',
                        verbose=2,
                        request_memory=4000,
                        extra_lines= ['should_transfer_files = YES', 'when_to_transfer_output = ON_EXIT', 'Requirements =  (Machine != "node128.icecube.wisc.edu")']
                        )

for index in range(500):
  for index_2 in [1, 2, 3, 4]:
        job.add_arg(str(index) + " " + str(index_2))

dagman = pycondor.Dagman('lee_muon_test_100B', submit=submit, verbose=2)
dagman.add_job(job)
dagman.build_submit()
