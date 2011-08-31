#!/usr/bin/python
"""Tool to copy latest output to baseline directories to create a new baseline for regression testing of AERO codes"""

__author__ = "Mark A. Potts (mpotts@hpti.com)"
__version__ = "$Revision: 0.1 $"
__date__ = "2011/08/11"

import sys, os, re, glob
#import argparse


def newBaseline(params):
  mycwd = os.getcwd()
  if(mycwd.find("Regression.d") < 0):
    os.chdir("Regression.d")  

  if(params[1] == 'ALL'):
    PROBLEM_NAMES=['statics','nlstatics','eigen','freqsweep','dynamics','nldynamics','impe','tempstatics','tempdynamics','tempnldynamics','tempnlstatics','test1','test11','test31']
  else:
    del params[0]
    PROBLEM_NAMES = params

  os.system("cp ../.hg/*.cache baseline")
  for problem_type in PROBLEM_NAMES:
    baselinepath = "baseline/"+problem_type
    if(os.path.exists(baselinepath) == 0):
      os.mkdir(baselinepath)
    command = "cp " + problem_type +"/* baseline/" + problem_type
    print "%s" % command
    os.system(command)
    
if __name__ == "__main__":

# parser = argparse.ArgumentParser(description='Build Regression Testing Inputs.')
# parser.add_argument('--c')
# basepath = parser.parse_args(sys.argv.split())

  params = sys.argv
  newBaseline(params)
