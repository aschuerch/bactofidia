#!/bin/bash python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 07:24:26 2017

@author: aschurch
"""
import argparse
import glob
import os
import sys
import logging
import errno
import subprocess
import re



#logging function

#command line arguments: sample names

# check if the samples are present as paired-end sequencing reads

# check if conda is installed

# check if there is a virtual environment 'snakemake'

# Determine read length/ config files

# concatenate for/rev and put into data folder

# check if it is on our hpc -> email and hpc command, else

##parsing##########################################################################
parser = argparse.ArgumentParser(description='''###########################################################################
############      Basic microbial WGS analysis pipeline    ################
##                                                                       ##
## for all samples in this folder. Please use samples with the same read ##
## length and preferably from the same species.                          ## 
##                                                                       ##
## Compressed sequencing files (fastq.gz)                                ##
## must be present in the same folder from where the script is called.   ##
##                                                                       ##
## Use only the sample names to call the script                          ##
##                                                                       ##
## Example:                                                              ##
##                                                                       ##
## python bactofidia.py ECO-RES-PR1-00001  ECO-RES-PR1-00002             ##
##                                                                       ##
##                                                                       ##
## Before running the pipeline for the first time, a virtual             ##
## environment needs to be created. Packages and versions are specified  ##
## in package-list.txt. See bioconda.github.io for available packages.   ##
##                                                                       ##
## Create the environment with                                           ##
##                                                                       ##
## conda create --file package-list.txt -n bactofidia_standard201709     ##
##                                                                       ##
##                                                                       ##
## Anita Schurch Aug 2017                                                ##
###########################################################################''')
parser.add_argument('-path', default='', help='Full path to the folder containing the sequencing files if not in this folder')
parser.add_argument('file', nargs='+', help='Sample name of the sequencing file (R1 and R2) as fastq.gz files')
argument = parser.parse_args()
path = argument.path
if not os.path.exists(path):
    path = os.getcwd()
#    sys.exit ('Please provide the full path to the folder containing the sequencing files with -path [/home/user/data/]')
seqfiles=[]
names= ''
for f in argument.file:
    for name in f:
        for seqfile in glob.glob(path+"/"+name+"*fastq.gz"):
            seqfiles.append(seqfile)

    names += f
    names += ", "

if (len(seqfiles) % 2)  != 0:
    sys.exit('Please provide the sample name to an R1 and and R2 sequencing file. ')



#####################################################################################
def mkdir_p(path):
    try:
        os.makedirs(path)
        logger.debug(path + " directory created")  
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
####################################################################################
def log():
    global logger 
    mkdir_p("log")

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(message)s')

    fh = logging.FileHandler('log/log_'+ '%(asctime)s'+'.txt')
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.debug('A new analysis was started with sample %s', names)
#####################################################################################
def check_conda():
    process = subprocess.Popen(["conda", "create", "-n", "snakemake", "-y", "snakemake", "python=3.5"], stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    
    if 'already exists' in stdout:
        logger.debug('snakemake environment exists') 
    elif 'conda is not' in stdout:
        logger.debug('Please see README.md to install bioconda')
        sys.exit()
    else:
        logger.debug('snakemake environment succesfully installed.')
####################################################################################
def activate_snakemake():
    process = subprocess.Popen(["source", "activate", "snakemake"], stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    print(stdout)
    
    
####################################################################################
def determine_read_length(fqfile):
    print (fqfile)
    process = subprocess.Popen ( ["zcat",fqfile,"|awk'{{if(NR%4==2) print length($1)}}' | sort | uniq -c | sort -rn | head -n 1 | rev | cut -f 1,1 -d " "| rev".format()], stdout=subprocess.PIPE)
    stdout = process.communicate()[0]
    print(stdout)    
    

  #  F=open(fqfile,'r')
   # re_seq=re.compile('^@')
    #numlist=[]
    #lenlist=[]

#    for occurence in re_seq.finditer(F.read()):
 #       gnom,num,seq=occurence.groups()
  #      numlist.append(num)
   #     lenlist.append(len(seq))
    #    print gnom,num,len(seq)
    
   # maxl=max(lenlist)

#####################################################################################
if __name__ == '__main__':
    log()
    logger.debug('Start of the analysis.')
    logger.debug('The results will be generated in results/')
    logger.debug('The logfiles will be generated in log/')
    check_conda()
   # activate_snakemake()
    determine_read_length(seqfiles[0])