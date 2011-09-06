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
        origfiles = fullpath + "/run.*"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
        origfiles = fullpath + "/scp.*"
        command = "cp "+ origfiles+ " "+infile 
        os.system(command)
  else:
    if(params[1] == 'ALL'):
      PROBLEM_NAMES=['statics','nlstatics','eigen','freqsweep','dynamics','nldynamics',\
                     'impe','tempstatics','tempdynamics','tempnldynamics',\
                     'tempnlstatics','dsvm1','dsvm11','dsvm31','dsvm13',\
                     'dsvm2','dsvm15','dsvm19','dsvm20','dsvm21','dsvm22',\
                     'dsvm23','dsvm24','dsvm25','dsvm27a','dsvm27b','dsvm29','dsvm30',\
                     'dsvm31','dsvm32','dsvm34','dsvm35','dsvm37','dsvm38',\
                     'dsvm39','dsvm40']
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
      MPIFILE.write("#!/bin/bash\n#PBS -N test\n#PBS -l nodes=1:ppn=8,walltime=3:00:00\n\n")
      MPIFILE.write("cd %s\n" % dirname)
     
      command = "chmod +x " + runfilename
      os.system(command)
#     command = "cp ../*.include ."
#     os.system(command)

      OUTPUT = ["gdisplac","stressvm","strainxx","strainxy","strainxz",\
                "strainxx","strainxy","strainxz","strainxx","strainxy", \
                "strainxz","inxforce","inyforce","inzforce","axmoment", \
                "aymoment","azmoment","energies","gvelocit","gacceler", \
                "displmod","rotatmod","gdispmod"]

      OUTPUT_EXTRAS = ["1"]
      OUTPUT2 = ""
      STATICS = ["sparse","skyline","mumps","spooles","gmres","direct",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr","FETI",\
                 "FETI DP","FETI DPH"]

      INCLUDE = ["\"../mesh.include\""]

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

      if(problem_type == "dsvm30"):
        OUTPUT = ["displmod"]
        OUTPUT2 = ["stressxx","reaction"]
        OUTPUT_EXTRAS = [" 1","1 NG 2","1 NG 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres","direct",\
                 "spooles pivot","mumps pivot","pcg","bcg","cr",\
                 "FETI DP","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm29"):
        OUTPUT = ["stressvm"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        STATICS = ["mumps pivot","gmres","spooles pivot"]
        NONLINEAR = ["maxitr 60\nnltol 1.0e-10\ndlambda 0.25 1.5"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm27b"):
        OUTPUT = ["displacz"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","skyline",\
                   "spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "dsvm27a"):
        OUTPUT = ["gtempera"]
        OUTPUT2 = ["heatflxz","heatreac"]
        OUTPUT_EXTRAS = [" 1"," 1"," 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","skyline",\
                   "spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "dsvm25"):
        OUTPUT = ["stressp1"]
        OUTPUT2 = ["stressp1"]
        OUTPUT_EXTRAS = ["1 NG 1","1 NG 2"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","skyline",\
                   "spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["global_cor_rbm_tol 1e-5\nconstraints multipliers\n*\nGRBM\n*","*","*","*"]

      if(problem_type == "dsvm24"):
        OUTPUT = ["displacz"]
        OUTPUT_EXTRAS = [" 1 NG 1","1 NG 1 modphase"]
        IMPE = ["freq 500.0","freq 500.0\ndamp 3.18309886e-5 0"]
        NAMELIST = ["IMPE\n","STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot",\
                   "mumps","spooles","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [IMPE,STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm23"):
        OUTPUT = ["displacy"]
        OUTPUT_EXTRAS = [" 1 11"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "dsvm22"):
        OUTPUT = ["displacz"]
        OUTPUT2 = ["displacz","stresszz","stresszz"]
        OUTPUT_EXTRAS = [" 1 NG 1"," 1 NG 2"," 1 NG 1"," 1 NG 2"]
        NAMELIST = ["STATICS\n","","OUTPUT\n","INCLUDE ","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        STATICS_OPTS = ["constraints penalty 1e12","constraints multipliers"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        CONTACT = ["../dsvm22.contact1","../dsvm22.contact2","../dsvm22.contact3",\
                   "../dsvm22.contact4"]
        OPTIONSLIST = [STATICS,STATICS_OPTS,OUTPUT,INCLUDE,CONTACT]
#        EXTRAS = ["constraints penalty 1e12","*","*","*"]
        EXTRAS = ["*","*","*","*","*"]


      if(problem_type == "dsvm21"):
        OUTPUT = ["displacy"]
        NAMELIST = ["STATICS\n","EIGEN\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr"]
        EIGEN = ["arpack\nnsbspv 3\nneigpa 1","arpack\nnsbspv 3\nneigpa 1"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,EIGEN,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*\nGEPS\nBUCKLE\n*","*"]

      if(problem_type == "dsvm20"):
        OUTPUT = ["gdispmod"]
        NAMELIST = ["STATICS\n","EIGEN\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        EIGEN = ["arpack\nnsbspv 8\nneigpa 4"]
        INCLUDE_FILE = "../" + problem_type + ".include"
        INCLUDE = [INCLUDE_FILE]
        OPTIONSLIST = [STATICS,EIGEN,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm19"):
        OUTPUT = ["displacy"]
        OUTPUT_EXTRAS = ["1 27"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","spooles pivot","mumps pivot","skyline",\
                   "mumps","spooles","direct","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE = ["../dsvm19.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "dsvm15"):
        OUTPUT = ["displacy"]
        OUTPUT2 = ["displacy"]
        OUTPUT_EXTRAS = ["1 38046","1 37735"]
        IMPE =    ["freqsweep 0 500 11 50\nrecons pade 2 4 5"]
        NAMELIST = ["IMPE\n","STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["spooles pivot","sparse","mumps pivot",\
                   "spooles","FETI DP","FETI DPH"]
        INCLUDE = ["../dsvm15.include"]
        OPTIONSLIST = [IMPE,STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm2"):
        OUTPUT = ["stressxx","strainxx"]
        OUTPUT_EXTRAS = ["1 NG 1"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps pivot","mumps","spooles",\
                   "spooles pivot","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE = ["../dsvm2.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]
      if(problem_type == "dsvm13"):
        OUTPUT = ["geigenpa"]
        NAMELIST = ["STATICS\n","EIGEN\n","GEPS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["mumps pivot","sparse","skyline","mumps","spooles",\
                   "direct","spooles pivot","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        GEPS = ["buckle"]
        EIGEN = ["arpack LA 4\nnsbspv 3\nneigpa 1\ntoleig 1.0e-6\nshift 100"]
        INCLUDE = ["../dsvm13.include"]
        OPTIONSLIST = [STATICS,EIGEN,GEPS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*","*"]

      if(problem_type == "dsvm31"):
        OUTPUT = ["stressvm","strainvm","strainxx","sp3direc",\
                "stressp3","ep3direc","stressxx","strainp3"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres","direct",\
                   "spooles pivot","mumps pivot","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE = ["../dsvm31.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "dsvm11"):
        OUTPUT = ["displmod","gdisplac","displacz"]
        OUTPUT_EXTRAS = ["1 5"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","mumps",\
                 "spooles pivot","mumps pivot",\
                 "FETI DP","FETI DPH","spooles"]
        NONLINEAR = ["maxitr 40\nnltol 1.0e-6\nrebuild 1",\
                   "maxitr 10\nnltol 1.0e-5\nrebuild 1"]
        INCLUDE = ["../dsvm11.include"]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*","*"]

      if(problem_type == "dsvm1"):
        OUTPUT = ["reaction"]
        NAMELIST = ["STATICS\n","OUTPUT\n","INCLUDE "]
        STATICS = ["sparse","skyline","mumps","spooles","gmres","direct",\
                   "spooles pivot","mumps pivot","pcg","bcg","cr",\
                   "FETI DP","FETI DPH"]
        INCLUDE = ["../dsvm1.include"]
        OPTIONSLIST = [STATICS,OUTPUT,INCLUDE]
        EXTRAS = ["*","*","*"]

      if(problem_type == "tempnlstatics"):
        EXTRAS = ["*","*","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*"]
        INCLUDE = ["\"../mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 0.0 0.0 202.4 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","NONLINEAR\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,NONLINEAR,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempnldynamics"):
        DYNAMICS = ["heat\t0.5\ntime\t1.0\t1.0\t200.0"]
        EXTRAS = ["*","newmark","*","*","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*"]
        INCLUDE = ["\"../mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 2719.0  0.0 202.4 0.0 0.0 0.0 871.0   0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","NONLINEAR\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,NONLINEAR,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempdynamics"):
        DYNAMICS = ["mech\t0.25\t0.5",\
                    "time\t1.0\t1.0\t200.0"]
        EXTRAS = ["*","newmark","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc...","*","*"]
        INCLUDE = ["\"../mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 2719.0  0.0 202.4 0.0 0.0 0.0 871.0   0.0 0.0 0.0 0.0"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "tempstatics"):
        NAMELIST = ["STATICS\n","OUTPUT\n","MATERIALS\n","INCLUDE "]
        EXTRAS = ["*","*","*","*id A   E   nu  rho    c   k     h   P   Ta  q     w   etc..."]
        INCLUDE = ["\"../mesh_temp.include\""]
        OUTPUT = ["gtempera","gtempvel","heatflxx","heatflxy","heatflxz","grdtempx","grdtempy","grdtempz"]
        MATERIALS = ["1   0.0 0.0 0.0 0.0    0.0 202.4 0.0 0.0 0.0 0.0   0.0 0.0 0.0 0.0"]
        OPTIONSLIST=[STATICS,OUTPUT,MATERIALS,INCLUDE]

      if(problem_type == "dynamics"):
        EXTRAS = ["*","newmark\nmech\t0.25000\t0.5000\n*\ttime step\ttotal time","*","*","*"]
        NAMELIST = ["STATICS\n","DYNAMICS\n","OUTPUT\n","INCLUDE "]
        OPTIONSLIST = [STATICS,DYNAMICS,OUTPUT,INCLUDE]
  
      if(problem_type == "nldynamics"):
        STATICS = ["sparse","spooles","spooles pivot","pcg",\
                   "FETI DP","FETI DPH"]
        OUTPUT = ["gdisplac","stressvm","strainxx","strainxz",\
                  "stressxx","stressxz","inzforce","axmoment",  "energies" ]
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
        EXTRAS = ["include \"../feti.include\"\n*","*","*","*","*"]
  
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
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("\n") == -1 )&
                  (NAMELIST[j].find("INCLUDE") == -1 )&
                  (OPTIONSLIST[j][i % len(OPTIONSLIST[j])].find("\t") == -1 )):
            idname = idname + "_" + OPTIONSLIST[j][i % len(OPTIONSLIST[j])]

        idname = idname.replace(" ","_",10)
  

        if(problem_type == "dsvm1"): 
          OUTPUT_FILENAME = idname+"_1.dat 1" + " NG 1\n" 
          OUTPUT_FILENAME = OUTPUT_FILENAME + "reaction "+ idname+"_2.dat 1" + " NG 2" 
        OUTPUT_FILENAME = idname+".dat" 
 
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
            FILE.write(" %s" % OUTPUT_FILENAME)
            FILE.write(" %s" % OUTPUT_EXTRAS[0])
            if(OUTPUT2 != ""):
              for jj in range(len(OUTPUT2)):
                FILE.write("\n")
                OUTPUT_FILENAME = idname + "_" + "%d.dat" % ( jj + 2)
                FILE.write("%s %s" % (OUTPUT2[jj],OUTPUT_FILENAME))
                FILE.write(" %s" % OUTPUT_EXTRAS[jj+1])
          FILE.write("\n")
          if(NAMELIST[j].find("FETI") != -1):
            FILE.write("include \"../feti.include\"\n*")
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
