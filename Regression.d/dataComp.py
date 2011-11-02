#!/usr/bin/python
"""Tool designed to compare data files produced by regression testing and report discrepancies"""

__author__ = "Mark A. Potts (mpotts@hpti.com)"
__version__ = "$Revision: 0.1 $"
__date__ = "2011/02/17"

import sys, os, re, md5, subprocess, math, glob, datetime, time

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


def directComp(basefile,file,SUMMARY_FILE,outstring):
  COMP = list(open(file).read().splitlines())
  BASE = list(open(basefile,"r").read().splitlines())
  if(file.find("eigen") != -1):
    BASE.sort()
    COMP.sort()
  MaxDiff = 0.0
  RelDiff = 0.0
  TotDiff = 0.0
  nSample = 0
  MaxDiffLine = -1
  MaxDiffLoc = -1
  MaxDiffVals = (-1,-1)
  if(len(BASE) != len(COMP)):
    print bcolors.WARNING + "WARNING-- %s and %s have different lengths!" % (file,basefile) + bcolors.ENDC
    print len(BASE)
    print len(COMP)
    outstring.append( "\tBaseline file length = %d, New file length = %d\n" %(len(BASE),len(COMP)))
    SUMMARY_FILE.write(outstring[0])
    result = 1
    return result
   
  for i in range(len(BASE)):
    if(i > 0):
      BASE[i].replace(","," ");
      COMP[i].replace(","," ");
      BASE[i].strip()
      COMP[i].strip()
      basewords = BASE[i].split()
      compwords = COMP[i].split()
      if(len(basewords) != len(compwords)):
#       print "WARNING--Different line lengths at line %d of %f!!" %(i,file)
        outstring.append( "\tWARNING--Different line lengths at line %d of %f!!" %(i,file) )
        result = 1
      else:
        for j in range(len(basewords)):
          if(not(basewords[j].isalpha())and not(compwords[j].isalpha())):
            diff = float(basewords[j]) - float(compwords[j])
            TotDiff = math.fabs(diff) + TotDiff
            nSample = nSample + 1
            if((0.5*(float(basewords[j]) + float(compwords[j])))>= 1.0e-5):
              if(diff/(0.5*(float(basewords[j]) + float(compwords[j])))> RelDiff):
                RelDiff = diff/(0.5*(float(basewords[j]) + float(compwords[j])))
            if(math.fabs(diff) >= MaxDiff):
              MaxDiff = math.fabs(diff)
              MaxDiffLine = i
              MaxDiffLoc = j
              MaxDiffVals = (basewords[j],compwords[j])


#  print bcolors.OKBLUE + "\t\t\t\tRel diff is %e on line %d at location %d" %(RelDiff,MaxDiffLine,MaxDiffLoc) + bcolors.ENDC
  outstring.append( "\tRel diff is %e on line %d at location %d\n" %(RelDiff,MaxDiffLine,MaxDiffLoc))
# SUMMARY_FILE.write(outstring[0])

#  print bcolors.OKBLUE + "\t\t\t\tvals were %s " % (MaxDiffVals,) + bcolors.ENDC
  outstring.append( "\tvals were %s \n" % (MaxDiffVals,))
# SUMMARY_FILE.write(outstring[0])
  if((MaxDiff > 1.0e-8)&(RelDiff > 1.0e-3)):
#   TotDiff = TotDiff/nSample
#   print "Average Difference = %e %e %d \n" % (TotDiff,TotDiff/nSample,nSample)
    result = 1
  else:
    result = 0
  return result

def dComp(params):
  import time, glob
  import subprocess

  summary_filename = "reg_test_summary"
  mycwd = os.getcwd()
  if(mycwd.find("Regression.d") == -1):
    os.chdir("Regression.d")
  if(os.path.exists(summary_filename)):
    last_mod = os.path.getmtime(summary_filename)
    time_now = time.mktime(time.localtime())
    print "time now = %d and last_mod = %d, diff = %d \n" % (time_now,last_mod,time_now - last_mod)
    if((time_now - last_mod) > 300):
      SUMMARY_FILE = open(summary_filename,"w")
    else:
      SUMMARY_FILE = open(summary_filename,"a")
  else:
    SUMMARY_FILE = open(summary_filename,"w")
  outstring = "\n\nSummary of regression test run as:\n%s\nfrom:%s" % (params,os.getcwd())
  SUMMARY_FILE.write(outstring)
  now = datetime.datetime.now()
  time_now = now.time()
  date_now = now.date()
  outstring = "\n\nTest started at %s on %s\n"% (time_now,date_now)
  SUMMARY_FILE.write(outstring)
  exactMatches = 0
  Total = 0
  loc = -1
  sloc = -1
  rloc = -1
  i = 0 
  pattern = re.compile("\-r")
  runLocal = re.compile("\-l")
  sendMail = re.compile("\-s")

  for s in params:
    if(re.search(pattern,s)):
       loc = i
    if(re.search(sendMail,s)):
       sloc = i
    if(re.search(runLocal,s)):
       rloc = i
    i=i+1

  if(sloc != -1):
    sendMail = 1
    del params[sloc]
  else:
     sendMail = 0
 
  if(loc != -1):
    run = 1
    del params[loc]
  else:
    run = 0 

  if(rloc != -1):
    lrun = 1
    del params[rloc]
  else:
    lrun = 0 

  files = [] 
  if((params[1] == 'ALL')|(params[1] == 'short')):
     
    if(params[1] == 'ALL'):
      PROBLEM_NAMES=['statics','nlstatics','eigen','dynamics','nldynamics',\
                   'freqsweep','impe','tempstatics','tempnlstatics',\
                   'tempdynamics','tempnldynamics','dsvm1','dsvm11','dsvm31',
                   'dsvm2','dsvm13','dsvm15','dsvm19','dsvm20','dsvm21','dsvm22',\
                   'dsvm23','dsvm24','dsvm25','dsvm27a','dsvm27b','dsvm29','dsvm30',\
                   'dsvm32','dsvm34','dsvm35a','dsvm35b','dsvm37','dsvm38','dsvm39',\
                   'dsvm40']
    else:
      PROBLEM_NAMES=['nlstatics','freqsweep','impe','tempstatics','tempnlstatics',\
                   'tempdynamics','tempnldynamics','dsvm1','dsvm31',
                   'dsvm20','dsvm21','dsvm22',\
                   'dsvm23','dsvm24','dsvm25','dsvm27a','dsvm27b','dsvm30',\
                   'dsvm35b','dsvm40']
    del params[0:1]

    for names in PROBLEM_NAMES:
      indir = names+"/"
      os.chdir(indir)
      if(run == 1):
        os.system("rm *.dat")
        qsubScript = "scp."+names
        retcode = subprocess.Popen(["qsub", qsubScript ], stdout=subprocess.PIPE).communicate()[0]
        numbers = filter(lambda x: x.isdigit(), retcode)
        qsubRetfile = "test.o"+numbers
        while True:
          if(os.path.exists(qsubRetfile)):
            break
          time.sleep(10)
      if(lrun == 1):
        os.system("rm *.dat")
        command = "./run."+names+" >& reg.out"
#       command = "sh ./scp."+names+" >& reg.out"
        os.system(command)
      os.chdir('../') 
      for infile in glob.glob( os.path.join(indir, '*.dat') ):
        files.append(infile)
  else: 
    del params[0]
    indirs = params
    for indir in indirs:
      indirp = indir+"/"
      print "%s" % indir
      os.chdir(indir)
      if(run == 1):
        os.system("rm *.dat")
        qsubScript = "scp."+indir
        retcode = subprocess.Popen(["qsub", qsubScript ], stdout=subprocess.PIPE).communicate()[0]
        numbers = filter(lambda x: x.isdigit(), retcode)
        qsubRetfile = "test.o"+numbers
        while True:
          if(os.path.exists(qsubRetfile)):
            break
          time.sleep(10)
      if(lrun == 1):
        os.system("rm *.dat")
        print "current directory is %s\n"% os.getcwd()
        command = "./run."+indir+" >& reg.out"
#       command = "sh ./scp."+indir+" >& reg.out"
        print "command is %s\n"% command
        os.system(command)
      os.chdir('../') 
      for infile in glob.glob( os.path.join(indirp, '*.dat') ):
        files.append(infile)

  result = 0
  
  for file in files:
     basefile = "baseline/"+file
     p = subprocess.Popen("md5sum " + file, shell = True, stdout=subprocess.PIPE) # don't forget to "import subprocess"
     newmd5 = p.stdout.readline().split() # now the md5 variable contains the MD5 sum
     p.wait() # some clean up
     p = subprocess.Popen("md5sum " + basefile, shell = True, stdout=subprocess.PIPE) # don't forget to "import subprocess"
     oldmd5 = p.stdout.readline().split() # now the md5 variable contains the MD5 sum
     p.wait() # some clean up
     Total += 1
     if(newmd5[0] == oldmd5[0]):
       print bcolors.OKGREEN + " \tExact match "+ bcolors.ENDC, file
       exactMatches += 1
       outstring = "\texact match " + file + "\n"
       SUMMARY_FILE.write(outstring)
     else:
       compstring = []
       result = directComp(basefile,file,SUMMARY_FILE,compstring)
       if(result == 1):
         print bcolors.FAIL + " \tDiscrepancy " + bcolors.ENDC, file
         outstring = "\tDiscrepancy " + file + "\n"
         SUMMARY_FILE.write(outstring)
         SUMMARY_FILE.write(compstring[0])
       else:
         exactMatches += 1
         print bcolors.OKGREEN + " \tMatch" + bcolors.ENDC, file
         print bcolors.OKBLUE    + "\t" +compstring[0] + bcolors.ENDC 
         outstring = "\tclose match " + file + "\n"
         SUMMARY_FILE.write(outstring)
         SUMMARY_FILE.write(compstring[0])
  outstring = "%s--%d matches found out of %d total cases\n" %(indir,exactMatches,Total)
  print outstring
  SUMMARY_FILE.write(outstring)
  now = datetime.datetime.now()
  time = now.time()
  date = now.date()
  outstring = "Test completed at %s on %s\n"% (time,date)
  SUMMARY_FILE.write(outstring)

  if(sendMail == 1):
    mail_file = open(".mail_command","w")
    next_build_num = list(open("/lustre/home/hudson/jobs/FEM build/nextBuildNumber","r").read().splitlines())
#   command = "mail -s \"Regression Test Summary\" mpotts@hpti.com < reg_test_summary"
    build_num = int(next_build_num[0])-1
    log_file = "/lustre/home/hudson/jobs/FEM\ build/builds/" + str(build_num) +"/log"
    command = "tail -90 " + log_file + "| mail -s \"Regression Test Summary\" mpotts@hpti.com "
    mail_file.write(command)
    print(command)
    os.system(command)

  sys.exit(result)

if __name__ == "__main__":
  params = sys.argv
  dComp(params)


