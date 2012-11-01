# CMake generated Testfile for 
# Source directory: /home/tac688/Research/FEM2/FEM
# Build directory: /home/tac688/Research/FEM2/FEM
# 
# This file includes the relevent testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
ADD_TEST(copyfiles "/home/tac688/Research/FEM2/FEM/Regression.d/copy_files.pl" "/home/tac688/Research/FEM2/FEM/Regression.d" "/home/tac688/Research/FEM2/FEM/Regression.d")
ADD_TEST(buildInputs "/home/tac688/Research/FEM2/FEM/Regression.d/buildInputs.py" "ALL")
ADD_TEST(testStatics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "-n" "statics")
ADD_TEST(testNLstatics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "nlstatics")
ADD_TEST(testEigen "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "eigen")
ADD_TEST(testDynamics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dynamics")
ADD_TEST(testNLDynamics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "nldynamics")
ADD_TEST(testIMPE "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "impe")
ADD_TEST(testFreqSweep "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "freqsweep")
ADD_TEST(testTempStatics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "tempstatics")
ADD_TEST(testTempNLStatics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "tempnlstatics")
ADD_TEST(testTempDynamics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "tempdynamics")
ADD_TEST(testTempNLDynamics "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "tempnldynamics")
ADD_TEST(testRtest1 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm1")
ADD_TEST(testRtest2 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-r" "dsvm2")
ADD_TEST(testRtest11 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-r" "dsvm11")
ADD_TEST(testRtest13 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm13")
ADD_TEST(testRtest19 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm19")
ADD_TEST(testRtest20 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm20")
ADD_TEST(testRtest21 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm21")
ADD_TEST(testRtest22 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm22")
ADD_TEST(testRtest23 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm23")
ADD_TEST(testRtest24 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-r" "dsvm24")
ADD_TEST(testRtest25 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm25")
ADD_TEST(testRtest27a "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm27a")
ADD_TEST(testRtest27b "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm27b")
ADD_TEST(testRtest29 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-r" "dsvm29")
ADD_TEST(testRtest30 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm30")
ADD_TEST(testRtest31 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm31")
ADD_TEST(testRtest32 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm32")
ADD_TEST(testRtest34 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm34")
ADD_TEST(testRtest35a "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm35a")
ADD_TEST(testRtest35b "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm35b")
ADD_TEST(testRtest37 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-r" "dsvm37")
ADD_TEST(testRtest39 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "dsvm39")
ADD_TEST(testRtest40 "/home/tac688/Research/FEM2/FEM/Regression.d/dataComp.py" "-l" "-s" "dsvm40")
ADD_TEST(sendResults "/home/tac688/Research/FEM2/FEM/Regression.d/sendResults.py")
SUBDIRS(Element.d)
SUBDIRS(Feti.d)
SUBDIRS(Driver.d)
SUBDIRS(Comm.d)
SUBDIRS(Corotational.d)
SUBDIRS(Dec.d)
SUBDIRS(HelmAxi.d)
SUBDIRS(Solvers.d)
SUBDIRS(Utils.d)
SUBDIRS(Parser.d)
SUBDIRS(Timers.d)
SUBDIRS(Threads.d)
SUBDIRS(Mortar.d)
SUBDIRS(Math.d)
SUBDIRS(Linpack.d)
SUBDIRS(Sfem.d)
SUBDIRS(Paral.d)
SUBDIRS(Problems.d)
SUBDIRS(GNU-getopt.d)
SUBDIRS(Hetero.d)
SUBDIRS(Material.d)
SUBDIRS(Rom.d)
SUBDIRS(Regression.d)
SUBDIRS(Acme.d)
SUBDIRS(Pita.d)
