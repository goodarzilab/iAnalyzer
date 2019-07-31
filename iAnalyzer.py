#!/usr/bin/python                                                                                                                           
import sys
import subprocess
import shlex
from multiprocessing import Pool
import glob
import os
import re
import argparse
import pandas as pd
import numpy as np
from process import process

#python iAnalyzer.py ./data/umi/metadata.txt ~sample.type hi_v_lo_logFC.txt                          
if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--printMode", help="Only print commands", dest='runMode', action='store_false')
  parser.add_argument("-r", "--runMode", help="Run commands (default)", dest='runMode', action='store_true')
  parser.add_argument("-a", "--aligner", help="Choose aligner package (default BWA)", type=str)
  parser.add_argument("-w", "--weighting", help="Choose weighting method for combining z-score. Options are 'n' (sample size), 'SE' (standard error), and 'SES' (standardized effect size) (default SES)", type=str)
  parser.add_argument("--hasUMI", help="UMI version (default)", dest='umi', action='store_true')
  parser.add_argument("--noUMI", help="Without UMI", dest='umi', action='store_false')
  parser.add_argument("--reflib", help="The sample that is the reference in univariate analysis (default=default).", type=str)
  parser.add_argument("--ref", help="The sample that is the reference in univariate analysis.", type=str)
  parser.add_argument("metadata", help="Specifications of the samples", type=str)
  parser.add_argument("formula", help="Design formula: e.g. type~cell", type=str)
  parser.add_argument("outfile", help="Output file format.", type=str)
  parser.set_defaults(aligner="bowtie2", weighting="SES", umi=True, runMode=True, reflib="default")
  args = parser.parse_args()

  iAnalyzeRdir='../'
  if ('iAnalyzeRdir' in os.environ):
    iAnalyzeRdir = os.getcwd()

  reflib = "{}/reflib/CRISPRi_v2_human/CRISPRi_v2_human_library".format(iAnalyzeRdir)
  if (args.reflib!="default"):
    reflib = args.reflib

  meta = pd.read_csv(args.metadata, sep="\t", header=0, index_col=0)
  print(meta)
  indir = os.path.dirname(args.metadata)
  os.chdir(indir)
  log= open("pipeline.txt", "wt")
  runner = process(metadata=meta, iAnalyzeRdir=iAnalyzeRdir, reffile=reflib, log=log, hasUMI=args.umi, aligner=args.aligner, runMode=args.runMode)
  runner.rm=False
  runner.collpase()
  runner.umi_extract()
  runner.trim()
  runner.align()
  runner.dedup()
  runner.count()
  analysis_type = "univariate"
  runner.rm=True
  dummy = args.formula.split('~')
  if (len(dummy[0])>0 and len(dummy[1])>0):
    analysis_type = "logit"
  print("Based on the formula, the analysis type is identified as {}".format(analysis_type))

  if (analysis_type == "univariate" and runner.rm==True):
    runner.univariate(os.path.basename(args.metadata),args.formula,args.ref,args.weighting,args.outfile)

  log.close()
