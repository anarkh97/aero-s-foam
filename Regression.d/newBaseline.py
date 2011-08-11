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
    PROBLEM_NAMES=['statics','nlstatics','eigen','freqsweep','dynamics','nldynamics','impe','tempstatics','tempdynamics','tempnldynamics','tempnlstatics']
  else:
    PROBLEM_NAMES = [params[1]]

  for problem_type in PROBLEM_NAMES:
    command = "cp " + problem_type +"/* baseline/" + problem_type
    os.system(command)
    print "%s" % command
    
if __name__ == "__main__":

# parser = argparse.ArgumentParser(description='Build Regression Testing Inputs.')
# parser.add_argument('--c')
# basepath = parser.parse_args(sys.argv.split())

  params = sys.argv
  newBaseline(params)
