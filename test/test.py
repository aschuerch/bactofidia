#!/bin/bash python

import shutil
import os
import sys
import subprocess


cwd = os.getcwd()
shutil.copy (cwd+"/test/Test_R1.fastq.gz", cwd)
shutil.copy (cwd+"/test/Test_R2.fastq.gz", cwd)

process = subprocess.Popen(['bash','bactofidia.sh', 'Test'], stdout=subprocess.PIPE)
stdout = process.communicate()[0]
print ('STDOUT:{}'.format(stdout))



#can I catch the exit status of the subprocess?

os.remove ("Test_R1.fastq.gz")
os.remove ("Test_R2.fastq.gz")


