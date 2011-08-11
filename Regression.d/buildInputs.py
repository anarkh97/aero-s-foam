#!/usr/bin/python
"""Tool to build input files used for regression testing of AERO codes"""

__author__ = "Mark A. Potts (mpotts@hpti.com)"
__version__ = "$Revision: 0.1 $"
__date__ = "2011/02/14"

import sys, os, re, glob
#import argparse


def buildInputs(params):
  mycwd = os.getcwd()
  print "at start, mycwd is %s \n"% mycwd
  if(mycwd.find("Regression.d") < 0):
    if(os.path.exists("Regression.d")==0):
      os.mkdir("Regression.d")  
    os.chdir("Regression.d")  

  ploc = -1
  bpath = re.compile("\-c")
  i = 0
  copyInputs = 0

  for s in params:
    if(re.search(bpath,s)):
       ploc = i
       copyInputs = 1
    i=i+1

  if(ploc != -1):
    basepath = params[ploc+1]
    del params[ploc+1]
    del params[ploc]

  if(copyInputs == 1): # just copy the all the subdirectories and input files from the basepath parameter
    import shutil
    print "listing for " + basepath
    listdir = os.listdir(basepath)
    mycwd = os.getcwd()
    for infile in listdir:
      fullpath = os.path.join(basepath,infile)
      if (os.path.isdir(fullpath)==True)&(infile[0] != "."):
        if(os.path.exists(os.path.join(mycwd,infile)) != True):
          print "creating directory : " + infile
          newdir = os.mkdir(infile)
        
        print fullpath, os.path.exists(fullpath) 
        origfiles = fullpath + "/*.inp"
        print origfiles, infile
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
        origfiles = fullpath + "/*.include"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
        origfiles = fullpath + "/run.*"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
        origfiles = fullpath + "/scp.*"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
  else:
    if(params[1] == 'ALL'):
      PROBLEM_NAMES=['statics','nlstatics','eigen','freqsweep','dynamics','nldynamics','impe','tempstatics','tempdynamics','tempnldynamics','tempnlstatics']
    else:
      PROBLEM_NAMES = [params[1]]

    for problem_type in PROBLEM_NAMES:
  
      if(os.path.exists(problem_type)==0):
        os.mkdir(problem_type)  

      os.chdir(problem_type)
      dirname = os.getcwd()
      runfilename = "run."+problem_type 
      qsubfilename = "scp."+problem_type 
      RUNFILE = open(runfilename,"w")
      MPIFILE = open(qsubfilename,"w")
      MPIFILE.write("#!/bin/bash\n#PBS -N test\n#PBS -l nodes=1:ppn=8,walltime=1:00:00\n\n")
      MPIFILE.write("cd %s\n" % dirname)
     
      command = "chmod +x " + runfilename
      os.system(command)
      command = "cp ../*.include ."
      os.system(command)

      OUTPUT = ["gdisplac","stressvm","strainxx","strainxy","strainxz",\
                "strainxx","strainxy","strainxz","strainxx","strainxy", \
                "strainxz","inxforce","inyforce","inzforce","axmoment", \
                "aymoment","azmoment","energies","gvelocit","gacceler", \
                "displmod","rotatmod","gdispmod"]

      STATICS = ["sparse","skyline","mumps","spooles","gmres","direct",\
                 "spooles pivot","mumps pivot",\
                 "frontal","pcg","bcg","cr","FETI",\
                 "FETI 1","FETI 2 OLD","FETI 2 NEW","FETI DP","FETI DPH"]

      INCLUDE = ["\"mesh.include\""]

      DYNAMICS = ["time\t0.0\t3.0e-5\t3.0e-2",\
                  "time\t0.0\t1.0e-4\t3.0e-2",\
                  "time\t0.0\t3.0e-4\t3.0e-2"]

      IMPE =      ["freq 10.0\ndamp 1e-6 1.0",\
                   "freq 10.0\ndamp 1e-7 1.0",\
                   "freq 10.0\ndamp 1e-5 1.0"]

      NONLINEAR = ["maxitr 10\nnltol 1.0e-6\nrebuild 1",\
                   "maxitr 10\nnltol 1.0e-7\nrebuild 1",\
                   "maxitr 10\nnltol 1.0e-5\nrebuild 1"]

      EIGEN = ["arpack\nnsbspv 20\nneigpa 12\ntoleig 1.0e-10\ntoljac 1.0e-05",\
               "arpack\nnsbspv 20\nneigpa 12\ntoleig 1.0e-10\ntoljac 1.0e-04",\
               "arpack\nnsbspv 20\nneigpa 12\ntoleig 1.0e-10\ntoljac 1.0e-6"]

      SHIFT = ["0","10","100","1000"]

      if(problem_type == "tempnlstatics"):
        EXTRAS = ["*","*","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*"]
        INCLUDE = ["\"mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 0.0 0.0 202.4 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempnldynamics"):
        DYNAMICS = ["heat\t0.5\ntime\t1.0\t1.0\t200.0"]
        EXTRAS = ["*","newmark","*","*","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*"]
        INCLUDE = ["\"mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 2719.0  0.0 202.4 0.0 0.0 0.0 871.0   0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","NONLINEAR\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,NONLINEAR,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempdynamics"):
        DYNAMICS = ["mech\t0.25\t0.5",\
                    "time\t1.0\t1.0\t200.0"]
        EXTRAS = ["*","newmark","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*","*"]
        INCLUDE = ["\"mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 2719.0  0.0 202.4 0.0 0.0 0.0 871.0   0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempstatics"):
        NAMELIST = ["STATICS\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        EXTRAS = ["*","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc..."]
        INCLUDE = ["\"mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 0.0    0.0 202.4 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0"]
        OPTIONSLIST=[STATICS,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "dynamics"):
        EXTRAS = ["*","newmark\nmech\t0.25000\t0.5000\n*\ttime step\ttotal time","*","*","*"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,INCLUDE]
  
      if(problem_type == "nldynamics"):
        EXTRAS = ["*","newmark\nmech\t0.25000\t0.5000\n*\ttime step\ttotal time","*","*","*","*","*"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,NONLINEAR,OUTPUT,INCLUDE]

      if(problem_type == "impe"):
        NAMELIST = ["IMPE\n","STATICS\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [IMPE,STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*","*"]
  
      if(problem_type == "freqsweep"):
        NAMELIST = ["STATICS\n","IMPE\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,IMPE,OUTPUT,INCLUDE]
        IMPE =    ["freqsweep 1. 3. 3 10\ndamp 1e-6 1.0",\
                   "freqsweep 1. 3. 3 10\ndamp 1e-7 1.0",\
                   "freqsweep 1. 3. 3 10\ndamp 1e-5 1.0"]
        EXTRAS = ["*","recons pade 2 9 10","*","*","*"]

      if(problem_type == "nlstatics"):
        EXTRAS = ["include \"feti.include\"\n*","*","*","*","*"]
  
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
      if(problem_type == "statics"):
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        EXTRAS = ["*","*","*"]
        OPTIONSLIST=[STATICS,OUTPUT,INCLUDE]

      if(problem_type == "eigen"):
        EXTRAS = ["*","*","*","*","*","*"]
        NAMELIST = ["STATICS\n","EIGEN\n","SHIFT ","OUTPUT\n","INCLUDE "]
        OPTIONSLIST=[STATICS,EIGEN,SHIFT,OUTPUT,INCLUDE]
  
      options = 0
      for i in range(len(OPTIONSLIST)):
        if(len(OPTIONSLIST[i]) >  options):
          options = len(OPTIONSLIST[i])
      for i in range(options):
        idname = problem_type
        for j in range(len(OPTIONSLIST)-1):
          if((NAMELIST[j] != "")&(OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find(" ") != 1)&
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("-") == -1 )&
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("\t") == -1 )):
            idname = idname + "_" + OPTIONSLIST[j][i % len(OPTIONSLIST[j])]

        idname = idname.replace(" ","_",10)
  
        OUTPUT_FILENAME = idname+".dat" + " 1"
   
        filename = idname+".inp"
        FILE = open(filename,"w")
        FILE.write("CONTROL\n");
        FILE.write(idname);
        if(idname.find("temp") != -1):
          FILE.write("\n2\n\"nodeset\"\n\"elemset\"\n*\n");
        else:
          FILE.write("\n1\n\"nodeset\"\n\"elemset\"\n*\n");
        for j in range(len(OPTIONSLIST)):
          FILE.write(NAMELIST[j])
          FILE.write(OPTIONSLIST[j][i  % len(OPTIONSLIST[j])])
          if(NAMELIST[j].find("OUTPUT") != -1 ):
            FILE.write(" %s" % OUTPUT_FILENAME);
          FILE.write("\n")
          if(NAMELIST[j].find("FETI") != -1):
            FILE.write("include \"feti.include\"\n*")
          else:
            FILE.write(EXTRAS[j])
          FILE.write("\n")
        FILE.write("END\n")
        FILE.close()
        command = "../../bin/aeros "
        re.compile("FETI");
        if(idname.find("FETI") != -1 ):
          command = command + "-n 2 --dec --nsub 4"
        print "Creating %s" % filename
        MPIFILE.write("echo mpirun -n 4 %s %s\n" % (command,filename.replace(" ","_")))
        MPIFILE.write("mpirun -n 4 %s %s\n" % (command,filename.replace(" ","_")))
        RUNFILE.write("echo %s %s\n" % (command,filename.replace(" ","_")))
        RUNFILE.write("%s %s\n" % (command,filename.replace(" ","_")))
      os.chdir('../')

if __name__ == "__main__":

# parser = argparse.ArgumentParser(description='Build Regression Testing Inputs.')
# parser.add_argument('--c')
# basepath = parser.parse_args(sys.argv.split())

  params = sys.argv
  buildInputs(params)
