#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
using namespace std;

#include <memory>
#include <Utils.d/dbg_alloca.h>

#include <GNU-getopt.d/getopt.h>

#include <Driver.d/Domain.h>
#include <Driver.d/StaticProbType.h>
#include <Driver.d/DynamProbType.h>
#include <Driver.d/EigenProbType.h>
#include <Driver.d/CondProbType.h>
#include <Driver.d/NLStaticProbType.h>
#include <Driver.d/NLDynamProbType.h>
#include <Driver.d/NLMat.h>
#include <Driver.d/TempProbType.h>
#include <Solvers.d/Solver.h>
#include <Problems.d/StaticDescr.h>
#include <Problems.d/DynamDescr.h>
#include <Problems.d/EigenDescr.h>
#include <Problems.d/CondDescr.h>
#include <Problems.d/NonLinStatic.h>
#include <Problems.d/NonLinDynam.h>
#include <Problems.d/TempDescr.h>
#include <Problems.d/ModalDescr.h>
#include <Problems.d/NLModalDescr.h>
#include <Problems.d/DEMProblem.h>
#include <Paral.d/MDStatic.h>
#include <Paral.d/MDInpcStatic.h>
#include <Paral.d/MDDynam.h>
#include <Paral.d/MDEigen.h>
#include <Paral.d/MDNLStatic.h>
#include <Paral.d/MDNLDynam.h>
#include <HelmAxi.d/FourierDescrip.h>
#include <HelmAxi.d/FourierProbTyp.h>
#include <HelmAxi.d/FourierHelmBCs.h>
#include <HelmAxi.d/MDAxiDesc.h>
#include <Driver.d/GeoSource.h>
#include <Dec.d/dec.h>
#include <Parser.d/DecInit.h>
#include <Sfem.d/Sfem.h>
#ifdef DISTRIBUTED
  #include <Pita.d/PitaNonLinDynam.h>
  #include <Pita.d/NLDistrTimeDecompSolver.h>
  #include <OOPita.d/PitaNonLinDynam.h>
  #include <OOPita.d/NlDriver.h>
  #include <OOPita.d/LinearDriver.h>
#endif
#include <Comm.d/Communicator.h>


// .... for different problems and hardware
void writeOptionsToScreen();

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef STRUCTOPT
#include <Structopt.d/Structopt_base.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#include <Structopt.d/Driver_opt.d/SubDomainFactory_opt.h>
#include <Structopt.d/Element_opt.d/Element_opt.h>

#include <Structopt.d/Paral_opt.d/MDStatic_opt.h>
#include <Structopt.d/Structopt_dec.h>
#endif

DecInit * decInit=0;

#ifdef TFLOP
extern map<int,double> weightList;

extern int optind;
extern "C" int getopt (
             int argc,
             const char *argv[],
             const char *optstring );
#endif

// .... global variables

#ifdef STRUCTOPT
Domain *domain = new Domain_opt();
std::auto_ptr<ElementFactory> elemFact(new ElementFactory_opt());
std::auto_ptr<GenSubDomainFactory<double> >   subDomainFactory(new GenSubDomainFactory_opt<double>());
std::auto_ptr<GenSubDomainFactory<DComplex> > subDomainFactoryC(new GenSubDomainFactory_opt<DComplex>());;
#else
Domain *domain = new Domain();
std::auto_ptr<ElementFactory> elemFact(new ElementFactory());
std::auto_ptr<GenSubDomainFactory<double> >   subDomainFactory(new GenSubDomainFactory<double>());
std::auto_ptr<GenSubDomainFactory<DComplex> > subDomainFactoryC(new GenSubDomainFactory<DComplex>());;
#endif

Sfem *sfem = new Sfem();

Connectivity *procMpcToMpc;

long totMemSky       = 0;
long totMemSparse    = 0;
long totMemSpooles   = 0;
long totMemMumps     = 0;

ThreadManager *threadManager = 0;

extern int yyparse(void);
extern FILE* yyin;

bool estFlag=false;
bool weightOutFlag=false;
bool nosa=false;
bool useFull=false;

int verboseFlag = 0;
int salinasFlag = 0;
int totalNewtonIter = 0;
int iterTotal = 0;

SysCom *syscom = 0;
Communicator *structCom = 0;
Communicator *heatStructCom;
Communicator *fluidCom;

// ... main program

#ifdef CREATE_DSO
extern "C"
int entrypoint(int argc, char** argv)
#else
int main(int argc, char** argv)
#endif
{
#ifdef __GNUC__
//  std::set_new_handler(&print_trace_handler);
#endif
 double initTime = getTime();
 double totalMemoryUsed = 0.0;

 int c;
 extern char *optarg;

 /**** DEFAULT ELEMENT WEIGHT VALUES  ***/
 /**** IF DIFFERENT FROM 1.0 add the weight of your object below
       weightList[objectNumber] = makePair(my_weight, my_true_weight)  ***/
 weightList[2] = 2.0;   // FourNodeQuad
 weightList[3] = 2.0;   // Therm3DQuad
 weightList[4] = 2.0;   // Triangle3
 weightList[8] = 3.0;   // ThreeNodeShell
 weightList[10] = 2.0;  // ThermQuadGal
 weightList[17] = 3.0;  // EightNodeBrick
 weightList[19] = 3.0;  // Membrane
 weightList[20] = 3.0;  // Compo3NodeShell
 weightList[2020] = 4.0;  // Compo4NodeShell
 weightList[23] = 3.0;  // Tetrahedral
 weightList[24] = 3.0;	// Pentahedral
 weightList[25] = 4.0;  // TenNodeTetrahedral
 weightList[32] = 3.0;  // HelmQuad8Gal
 weightList[38] = 3.0;  // HelmTri6Gal
 weightList[40] = 2.0;  // TetraHelmGal
 weightList[41] = 2.0;  // TetraHelmGLS
 weightList[42] = 3.0;  // Tetra10HelmGal
 weightList[44] = 3.0;  // HelmBrickGLS
 weightList[45] = 1.0;  // HelmBrick
 weightList[46] = 2.0;  // Therm3NoShell
 weightList[4646] = 2.0;  // Therm4NoShell
 weightList[48] = 2.0;  // QuadConvec
 weightList[49] = 2.0;  // TriangleConvec
 weightList[50] = 3.0;  // TetraTherm
 weightList[51] = 4.0;  // ThermBrick
 weightList[72] = 4.0;  // Brick20
 weightList[88] = 4.0;  // FourNodeShell
 weightList[90] = 3.0;  // HelmPenta
 weightList[91] = 6.0;  // 3d 32 node serendipity brick
 weightList[92] = 5.0;  // 3d 26 node serendipity wedge
 weightList[93] = 5.0;  // 3d 32 node serendipity brick
 weightList[94] = 4.0;  // 3d 26 node serendipity wedge
 weightList[95] = 3.0;  // 3d HelmIsoParamHexa
 weightList[96] = 5.0;  // 3d HelmIsoParamTetra
 weightList[97] = 4.0;  // 3d 15 node serendipity wedge
 weightList[104]= 4.0;  // 3d 18 node Lagrange wedge
 weightList[201] = 3.0; // 3d 8 node brick
 weightList[202] = 3.0; // green lagrange
 weightList[102] = 3.0; // 3d LEIsoParamHexa
 weightList[301] = 1.0; // 2d 4-node sloshing (fluid) quadrilateral
 weightList[302] = 1.0; // 2d 2-node free-surface (fluid)
 weightList[311] = 3.0; // 3d 4-node sloshing (fluid) tetrahedron
 weightList[312] = 2.0; // 3d 3-node free-surface (fluid) triangle
 weightList[321] = 2.0; // 2d 4-node hydroelastic vibration (fluid) quadrilateral
 weightList[331] = 3.0; // 3d 4-node hydroelastic vibration (fluid) tetrahedral

#if defined(USE_MPI)
 SysCom theCom(argc,argv);
 syscom = &theCom;
 structCom = syscom;
#endif

 // Default number of threads equal to one
 int numThreads =  1;
 int numSubdomains = 1;
 int numProcessors = 1;
 int topFlag    = -1;
 int numClusters = 1;

 bool callDec = false; // invoque Dec ? (-D or --call-dec)
 bool exitAfterDec = false;
 bool callSower = false;
 bool exitAfterSower = false;
 // Process command line arguments
 //
 // -n = number of threads lok
 // -d = decomposition file name lok
 // -v = verbose output lok
 // -c = verbose output plus contact status change info
 // -t = output topdomdec file and exit lok
 // -T = output topdomdec file with gaps renumbered sequentially lok
 // -m = output topdomdec file with number of element sets equal to
 //      the number of materials in the input file.lok
 // -r = output topdomdec file for axisymmetric mesh and exit
 // -p = output primal residual at each FETI iteration to a file

 // Process command line arguments for DEC
 //
 // -n <number of threads> OK
 // -d <decomposition file name> OK
 // -s <number of subdomains> OK
 // -t = output topdomdec file and exit OK
 // -T = output topdomdec file with gaps renumbered sequentially OK
 // -m = output topdomdec file with number of element sets equal to OK
 //      the number of materials in the input file.
 // -e = output memory estimate in a file with a .memory extension OK
 // -w = output weights with a .weight extension

 // -e = output memory estimate in a file with a .memory extension OK
 // -w = output weights with a .weight extension

 if(argc == 1) {
   writeOptionsToScreen();
 }

 // getopt_long
 int option_index = 0; // will hold index for long options
 static struct option long_options[] = {
   {"with-dec", 0, 0, 1000},
   {"dec", 0, 0, 1000},
   {"exit", 0, 0, 1002},
   {"deter", 0, 0, 1005},
   {"use-weight-from", 1, 0, 1004},
   {"threads-number", 1, 0, 'n'},
   {"decomposition-filename", 1, 0, 'd'},
   {"nsub", 1, 0, 1001},
   {"subdomains-number", 1, 0, 1001},
   {"processors-number", 1, 0, 1003},
   {"output-topdomdec", 0, 0, 't'},
   {"output-topdomdec-asymetric-mesh", 1, 0, 'r'},
   {"output-primal", 1, 0, 'p'},
   {"contact-status-verbose", 1, 0, 'c'},
   {"output-topdomdec-gaps-renumbered-sequentially", 0, 0, 'T'},
   {"output-topdomdec-element-equals-materials-number", 0, 0, 'm'},
   {"output-topdomdec-element-equals-materials-numberc-gaps-renumbered-sequentially", 0, 0, 'M'},
   {"output-match", 0, 0, 'P'},
   {"output-memory-estimate", 0, 0, 'e'},
   {"mem", 0, 0, 'e'},
   {"output-weights", 0, 0, 'w'},
   {"load", 0, 0, 'w'},
   {"verbose", 1, 0, 'v'},
   {"with-sower",0,0, 1010},
   {"sower",0,0, 1010},
   {"nclus", 1, 0, 1012},
   {0, 0, 0, 0}
 };
 // end getopt_long

 filePrint(stderr,"\n --------- R U N  PARAMETERS ----------\n");
 FILE * weightFile;
 while ((c = getopt_long(argc, argv, "n:d:p:v:c:DVtTPmMr:Pfs:",long_options, &option_index)) != -1)
      switch (c) {
      case 1000 :  // call dec from FEM
	callDec = true;
	break;
      case 1001 :
	numSubdomains = atoi(optarg);
        if(numSubdomains <= 0) numSubdomains = 1;
	break;
      case 1002 :  // exit after dec
	exitAfterDec = true;
	geoSource->setExitAfterDec(true);

	exitAfterSower = true; //TG just in case
	break;
      case 1003 :  // number of processors for dec
	numProcessors = atoi(optarg);
	break;
      case 1004 :
	weightFile = fopen(optarg, "r");
	double w;
	int k;
	if(weightFile)
	  {
	    int i=1;
	    char c;
	    fpos_t position;
	    while(!feof(weightFile))
	      {
		int res=fscanf(weightFile,"%d",&k);
		if(res == 0 || res == EOF)
		  {
		    cerr << "*** WEIGHT FILE CORRUPTED AT LINE " << i << " bad object ID." <<endl;
		    exit(1);
		  }
		fgetpos (weightFile, &position);
		while((c = fgetc(weightFile))==' ')
		  ;
		if(c=='\n')
		  {
		    cerr << "*** WEIGHT FILE CORRUPTED AT LINE " << i << " : no weight specified !" << endl;
		    exit(1);
		  }
		fsetpos (weightFile, &position);
		res=fscanf(weightFile,"%lf",&w);
		if(res == 0 || res == EOF)
		   {
		     cerr << "*** WEIGHT FILE CORRUPTED AT LINE " << i << " : no weight specified !" <<endl;
		     exit(1);
		   }
		weightList[k]=w;
		i++;
	      }
	  }
	else
	  {
	    filePrint(stderr," *******************************************\n");
	    filePrint(stderr," *** ERROR: Cannot open weight file %s ***\n",
           optarg );
           filePrint(stderr," *******************************************\n");
           exit(-1);
	  }
	break;
      case 1005 :
	nosa=true;
	break;
      case 1010 :
	callSower = true;
	domain->setSowering(true);
	break;
      case 1012 :
        numClusters = atoi(optarg);
        if(numClusters <= 0) numClusters = 1;
        break;
      case 'w':
	weightOutFlag = true;
        break;
      case 'e':
        estFlag = true;
	useFull = true;
        break;
      case 'n':
        numThreads = atoi(optarg);
        if(numThreads <= 0) numThreads = 1;
#ifdef USE_OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(numThreads);
#endif
        break;
      case 'd': {
          geoSource->getCheckFileInfo()->decomposition = optarg;
          FILE *f;
          if((f=fopen(optarg,"r"))==(FILE *) NULL ) {
            filePrint(stderr," *******************************************\n");
            filePrint(stderr," *** ERROR: Cannot open decomposition %s ***\n", optarg);
            filePrint(stderr," *******************************************\n");
            exit(-1);
          }
          geoSource->getCheckFileInfo()->decPtr = f;
        }
        break;
      case 'v':
        verboseFlag = 1;
        domain->setVerbose();
        domain->solInfo().fetiInfo.printNumber = atoi(optarg);
        filePrint(stderr," ... Setting Verbose Output Mode    ... \n");
        break;
      case 'V':
        verboseFlag = 1;
        domain->setVerbose();
        domain->solInfo().fetiInfo.printNumber = 10;
        filePrint(stderr," ... Setting Verbose Output Mode    ... \n");
        break;
      case 't':
        topFlag = 0;
        domain->solInfo().setProbType(SolverInfo::Top);
	break;
      case 'T':
        topFlag = 1;
        domain->solInfo().setProbType(SolverInfo::Top);
        break;
      case 'm':
        topFlag = 2;
        domain->solInfo().setProbType(SolverInfo::Top);
        break;
      case 'M':
        topFlag = 7;
        domain->solInfo().setProbType(SolverInfo::Top);
        break;
      case 'x':
	domain->setOutputMatchInTop(true);
	break;
      case 'r':
        topFlag = 3 + atoi(optarg);
        domain->solInfo().setProbType(SolverInfo::Top);
        break;
      case 'p':
        domain->solInfo().fetiInfo.primalFlag = 1;
        break;
      case 'c':
        domain->solInfo().fetiInfo.contactPrintFlag = atoi(optarg);
        break;
      case '?':
        {
          filePrint(stderr," *******************************************\n");
          filePrint(stderr," *** ERROR: Check command line argument  ***\n");
          filePrint(stderr," *******************************************\n");
          exit(-1);
        }
      }

 if(optind < argc - 1) {
   filePrint(stderr," ******************************************************\n");
   filePrint(stderr," *** ERROR: Command line contained errors. Aborting ***\n");
   filePrint(stderr," ******************************************************\n");
   exit(-1);
 }

 if(optind == argc) {
   filePrint(stderr," **************************************\n");
   filePrint(stderr," *** ERROR: No input file specified ***\n");
   filePrint(stderr," **************************************\n");
   exit(-1);
 }

 // Open input file
 if(optind == argc-1) {
   geoSource->getCheckFileInfo()->checkfile = argv[argc-1];
   //FILE *nin = freopen(argv[argc-1],"r",stdin);
   FILE *nin = yyin = fopen(argv[argc-1],"r");
   if(nin == 0) {
     filePrint(stderr," *******************************************\n");
     filePrint(stderr," *** ERROR: Could not open input file: %s\n",argv[argc-1]);
     filePrint(stderr," *******************************************\n");
     exit(-1);
   }
 }
 MatrixTimers &times = domain->getTimers();
 double t1 = getTime();

 // Read input file and time reading
 if(verboseFlag) filePrint(stderr," ... Reading Input File             ...\n");
 startTimerMemory(times.readTime, times.memoryParse);
 int error = yyparse();
 stopTimerMemory(times.readTime, times.memoryParse);
 if(verboseFlag) filePrint(stderr," ... Parsed Input File In %8.2e sec and %8.2e Mb ...\n",
	                   (getTime() - t1)/1000.0, times.memoryParse/oneMegaByte);

 // Check if Input file had any errors
 if(error) {
   filePrint(stderr," ****************************************************\n");
   filePrint(stderr," *** ERROR: Input file contained errors. Aborting ***\n");
   filePrint(stderr," ****************************************************\n");
   exit(error);
 }
 if(decInit != 0 && decInit->skip==false) { // dec initializers in parser !
   if(numThreads == 1) // command line option prevail !!!
     numThreads = decInit->nthreads;
   if(numSubdomains == 1)
     numSubdomains = decInit->nsubs;
   if(numProcessors == 1)
     numProcessors = decInit->nproc;
   if(exitAfterDec==false)
     exitAfterDec = decInit->exitAfterDec;
   if(decInit->weight)
     weightOutFlag = true;
   if(decInit->memory) {
     estFlag = true;
     useFull = true;
   }
   if(!nosa) nosa = decInit->nosa; // PJSA
   if(topFlag < 0) {
     callDec=true;
     if(geoSource->getCheckFileInfo()->decPtr == 0 && decInit->file !=0)
       geoSource->getCheckFileInfo()->checkfile = decInit->file;
   }
 }

#define MAX_CODES 4
#define FLUID_ID 0
#define STRUC_ID 1
#define HEAT_ID  2

 // We do a split
 Communicator* allCom[MAX_CODES];
 if(domain->solInfo().aeroheatFlag >= 0 || domain->solInfo().thermohFlag >= 0) {
   syscom->split(HEAT_ID, MAX_CODES, allCom);
   structCom = allCom[HEAT_ID];
   heatStructCom = allCom[STRUC_ID]; // comunicator between the thermal and
                                     // mechanical structure codes
 }
 else {
   syscom->split(STRUC_ID, MAX_CODES, allCom);
   structCom = allCom[STRUC_ID];
   heatStructCom = allCom[HEAT_ID]; // comunicator between the thermal and
                                    // mechanical structure codes
 }
 fluidCom = allCom[FLUID_ID];

 threadManager = new ThreadManager(numThreads);
 if(threadManager->numThr() != numThreads) { //HB: for checking purpose
   filePrint(stderr," *** WARNING: number of threads requested: %d\n",numThreads);
   filePrint(stderr,"              number of threads created  : %d\n",threadManager->numThr());
   filePrint(stderr," -> tip: if you are running on a Linux platform\n");
   filePrint(stderr,"         you may need to activate OpenMP and compile with an OpenMP compliant\n");
   filePrint(stderr,"         compiler (for instance, icpc or g++ version 4.2)\n");
 }
 fprintf(stderr,"KW: OK.\n");
 if(geoSource->binaryInput) geoSource->readGlobalBinaryData(); // SOWERX
#ifdef SOWER_SURFS
 else {
#endif
   fprintf(stderr,"KW: OK. I'm here.\n");
   // HB for checking Mortar & generating LMPCs from Mortar tied conditions
   if(topFlag < 0) {
     domain->SetUpSurfaces(&(geoSource->GetNodes()));
     if(!callSower) {
       domain->ProcessSurfaceBCs();
       domain->SetMortarPairing();
       if(!domain->tdenforceFlag()) {
         if(domain->solInfo().isNonLin()) { // for nonlinear statics and dynamics just process the tied surfaces here
           domain->InitializeStaticContactSearch(MortarHandler::TIED);
           domain->PerformStaticContactSearch(MortarHandler::TIED);
           domain->ExpComputeMortarLMPC(MortarHandler::TIED);
           domain->CreateMortarToMPC();
         }
         else {
           domain->ComputeMortarLMPC();
           domain->computeMatchingWetInterfaceLMPC();
           domain->CreateMortarToMPC();
         }
       }
     }
//#ifdef MORTAR_DEBUG
     fprintf(stderr,"KW: going to print surface entities.\n");
     domain->PrintSurfaceEntities();
//     domain->PrintMortarConds();
//#endif
     //domain->printLMPC();
     //if(domain->solInfo().fetiInfo.c_normalize) domain->normalizeLMPC(); // PJSA 5-24-06
   }
#ifdef SOWER_SURFS
 }
#endif

 if(domain->solInfo().type != 2 && !geoSource->getDirectMPC())
   geoSource->addMpcElements(domain->getNumLMPC(), *(domain->getLMPC()));

 if((domain->solInfo().type != 2 || (!domain->solInfo().isMatching && (domain->solInfo().fetiInfo.fsi_corner != 0))) && !domain->solInfo().HEV)
   geoSource->addFsiElements(domain->getNumFSI(), domain->getFSI());

 if(!geoSource->binaryInput) {
   domain->setUpData();
 }

 if(callDec) {
//   if(domain->solInfo().type == 2 || domain->solInfo().type == 3) { // DEC requires FETI or BLOCKDIAG to be activated (see below)
     Dec::dec(numProcessors, numThreads, numSubdomains, topFlag);
     if(exitAfterDec && !callSower) {
       filePrint(stderr," ... Exiting after Dec run          ...\n");
       filePrint(stderr," --------------------------------------\n");
       delete threadManager; // PJSA
       exit(0);
     }
//   }
//   else {
//     //  It doesn't make sense to use DEC outside of FETI. However when preparing files, user might not specify
//     //  the FETI solver in the input file as they are not solving yet. In that case proper FETI initializations do
//     //  not occur and the decomposition is different than what it would have been if the FETI keyword was here.
//     //  Thus as of 09/14/06 FETI has to be there to be able to use the decomposer.
//     filePrint(stderr,"*******************************************\n");
//     filePrint(stderr,"*** ERROR: FETI command missing            \n");
//     filePrint(stderr,"*******************************************\n");
//     delete threadManager; // PJSA
//     exit(0);
//   }
 }
 useFull = true; // or TenNodeTetraHedral will crush ! (bad design not from me !)

// if(!geoSource->binaryInput) {
//   domain->setUpData();
// }

 if(callSower) {
   filePrint(stderr," ... Writing Distributed Binary Input Files ... \n");
   geoSource->writeDistributedInputFiles(numClusters, domain); //add domain as argument for surfaces
   if(exitAfterSower) {
     filePrint(stderr," ... Exiting after Sower run        ...\n");
     filePrint(stderr," --------------------------------------\n");
     exit(0);
   }
 }

 // Make TOPDOM/DEC input file and exit from fem
 if(topFlag >= 0) {
   if(topFlag == 5 || topFlag == 6)
     domain->makeAxiTopFile(topFlag,fourHelmBC->numSlices);
   else if(domain->probType() == SolverInfo::Top) {
     fprintf(stderr," ... Memory to Parse       %14.3f Mb\n",
             times.memoryParse/oneMegaByte);
     fprintf(stderr," ... Memory to Set up Data %14.3f Mb\n",
             times.memorySetUp/oneMegaByte);
     domain->makeTopFile(topFlag);
   }
   exit(-1);
 }

 if(domain->solInfo().noninpc || domain->solInfo().inpc) domain->initSfem();

 // 1. check to see if a decomposition has been provided or requested (three options: -d, --dec, or DECOMP)
 bool domain_decomp = (geoSource->getCheckFileInfo()->decPtr || callDec || decInit || geoSource->binaryInput);
 // 2. check to see how many cpus are available (mpi processes and threads)
#ifdef USE_MPI
 bool parallel_proc = (structCom->numCPUs()*threadManager->numThr() > 1 || geoSource->binaryInput);
#else
 bool parallel_proc = (threadManager->numThr() > 1);
#endif
 // 3. choose lumped mass (also pressure and gravity) and diagonal "solver" for explicit dynamics 
 if(domain->solInfo().newmarkBeta == 0) {
   domain->solInfo().subtype = 10;
   domain->solInfo().getFetiInfo().solvertype = FetiInfo::diagonal;
   if(parallel_proc || domain_decomp) domain->solInfo().type = 3;
   geoSource->setMRatio(0.0);
   geoSource->setConsistentQFlag(false);
   geoSource->setConsistentPFlag(false);
   if(verboseFlag) filePrint(stderr, " ... Explicit Dynamics: lumped mass matrix, gravity and pressure will be used ... \n");
 }

 if(domain->solInfo().aeroFlag >= 0)
   filePrint(stderr," ... AeroElasticity Flag   = %d\n", domain->solInfo().aeroFlag);
 if(domain->solInfo().thermoeFlag >= 0)
   filePrint(stderr," ... ThermoElasticity Flag = %d\n", domain->solInfo().thermoeFlag);
 if(domain->solInfo().aeroheatFlag >= 0)
   filePrint(stderr," ... AeroThermo Flag       = %d\n", domain->solInfo().aeroheatFlag);
 if(domain->solInfo().thermohFlag >= 0)
   filePrint(stderr," ... ThermoElasticity Flag = %d\n", domain->solInfo().thermohFlag);

 // Domain Decomposition tasks
 //   type == 2 (FETI) and type == 3 (BLOCKDIAG) are always Domain Decomposition methods
 //   type == 1 && iterType == 0 (PCG) is a Domain Decomposition method only if a decomposition is provided or requested
 //   type == 0 && subtype == 9 (MUMPS) is a Domain Decomposition method only if a decomposition is provided or requested
 if(domain->solInfo().type == 2 || domain->solInfo().type == 3
    || (domain->solInfo().type == 1 && domain->solInfo().iterType == 0 && domain_decomp)
    || (domain->solInfo().type == 0 && domain->solInfo().subtype == 9 && domain_decomp)) {

   if(parallel_proc) {
#ifdef USE_MPI
     if(structCom->numCPUs() > 1) {
       if(threadManager->numThr() == 1)
         filePrint(stderr, " ... Parallel processing: Launching %d MPI processes ...\n", structCom->numCPUs());
       else
         filePrint(stderr, " ... Parallel processing: Launching %d MPI processes and forking %d threads per MPI process ...\n", structCom->numCPUs(), threadManager->numThr());
     } else
#endif
     filePrint(stderr, " ... Parallel processing: Forking %d Threads ...\n", threadManager->numThr());
   }

   switch(domain->probType()) {
     case SolverInfo::Dynamic: 
     case SolverInfo::TempDynamic:
       {
        if(domain->solInfo().mdPita) { // Not implemented yet
          filePrint(stderr, " ... PITA does not support multidomain - Aborting...\n");
        } else {
          MultiDomainDynam dynamProb(domain);
          DynamicSolver < MDDynamMat, DistrVector, MultiDomDynPostProcessor,
       		MultiDomainDynam, double > dynamSolver(&dynamProb);
          dynamSolver.solve();
          fflush(stderr);
        }
       }
       break;

      case SolverInfo::NonLinStatic: {
        MDNLStatic nlstatic(domain);
        NLStaticSolver < ParallelSolver,DistrVector,MultiDomainPostProcessor,
                      MDNLStatic, DistrGeomState >
                      nlsolver(&nlstatic);
        nlsolver.solve();
	}
       break;
     case SolverInfo::ArcLength: {
	MDNLStatic nlstatic(domain);
	NLStaticSolver < ParallelSolver,DistrVector,MultiDomainPostProcessor,
                      MDNLStatic, DistrGeomState >
                      nlsolver(&nlstatic);
	nlsolver.arclength();
	}
       break;
     case SolverInfo::Modal: { // CBM
        GenMultiDomainEigen<double> eigenProb(domain);
        EigenSolver<MDDynamMat, GenDistrVector<double>, GenDistrVectorSet<double>, GenMultiDomainEigenPostProcessor<double>,
                    GenMultiDomainEigen<double> > * eigenSolver = 0;
        switch(domain->solInfo().eigenSolverType) {
           case SolverInfo::Arpack :
             {
#ifdef USE_ARPACK
	       eigenSolver = new SymArpackSolver<MDDynamMat, GenDistrVector<double>, GenDistrVectorSet<double>,
                                                GenMultiDomainEigenPostProcessor<double>, GenMultiDomainEigen<double> > (&eigenProb);
#else
              filePrint(stderr," *** ERROR: executable not linked with ARPACK. See flag USE_ARPACK in Makefile.\n");
              exit(-1);
#endif
             }
             break;
           default :
             {
              filePrint(stderr,"ERROR: Eigensolver type is not implemented for multiple domain. Please use ARPACK.\n");
              exit(-1);
             }
           }
           eigenSolver->solve();
       }
       break;
     case SolverInfo::AxiHelm: {
         MDAxiDesc fDescr(domain, fourHelmBC, globalMPCs, globalScatter);
         fDescr.solve();
       }
       break;
     case SolverInfo::Static: {
       if(geoSource->isShifted()) filePrint(stderr, " ... Frequency Response Analysis ");
       if(domain->isComplex()) {
         if(domain->solInfo().inpc) cerr << "inpc not implemented for complex domain \n";
         if(geoSource->isShifted()) filePrint(stderr, "in Complex Domain ...\n");
#ifdef STRUCTOPT
	 // quick-fix, for debugging SOpt
	 if (dynamic_cast<Domain_opt*>(domain)->getStructoptFlag())
	   {
	     // structopt requested
#ifdef DISTRIBUTED
	     GenDistrDomain_opt<DComplex> dd(dynamic_cast<Domain_opt*>(domain));
#else
	     GenDecDomain_opt<DComplex> dd(dynamic_cast<Domain_opt*>(domain));
#endif
	     dd.structoptSolve();
	   }
	 else
#endif
	   {
         GenMultiDomainStatic<complex<double> > statProb(domain);
         StaticSolver<complex<double>, GenMDDynamMat<complex<double> >, GenDistrVector<complex<double> >,
                      GenMultiDomainPostProcessor<complex<double> >, GenMultiDomainStatic<complex<double> >,
                      GenDistrVector<complex<double> > >
           statSolver(&statProb);
         statSolver.solve();
	   }
       }
       else {
         if(geoSource->isShifted()) filePrint(stderr, "in Real Domain ...\n");
         if(domain->solInfo().inpc) {
           GenMultiDomainInpcStatic<double> statProb(domain);
           StaticSolver<double, AllOps<double>, DistrBlockVector<double>,
                        GenMultiDomainInpcPostProcessor<double>, GenMultiDomainInpcStatic<double>,
                        DistrBlockVector<complex<double> > >
                statSolver(&statProb);
           statSolver.solve();
         }
         else {
#ifdef STRUCTOPT
	 // quick-fix, for debugging SOpt
	 if (dynamic_cast<Domain_opt*>(domain)->getStructoptFlag())
	   {
	     // structopt requested
#ifdef DISTRIBUTED
	     GenDistrDomain_opt<double> dd(dynamic_cast<Domain_opt*>(domain));
#else
	     GenDecDomain_opt<double> dd(dynamic_cast<Domain_opt*>(domain));
#endif
	     dd.structoptSolve();
	   }
	 else
#endif
	   {
         GenMultiDomainStatic<double> statProb(domain);
         StaticSolver<double, GenMDDynamMat<double>, GenDistrVector<double>,
                      GenMultiDomainPostProcessor<double>, GenMultiDomainStatic<double>,
                      GenDistrVector<complex<double> > >
              statSolver(&statProb);
         statSolver.solve();
	   }
         }
       }
       break;
     }
     case SolverInfo::HelmholtzFreqSweep:
     case SolverInfo::Helmholtz: {
       filePrint(stderr, " ... Acoustic Scattering Helmholtz Analysis ");
       if(domain->isComplex()) {
         filePrint(stderr, "in Complex Domain ...\n");
         GenMultiDomainStatic<complex<double> > FAProb(domain);
         StaticSolver<complex<double>, GenMDDynamMat<complex<double> >, GenDistrVector<complex<double> >,
                      GenMultiDomainPostProcessor<complex<double> >, GenMultiDomainStatic<complex<double> >,
                      GenDistrVector<complex<double> > >
            FASolver(&FAProb);
         FASolver.solve();
       }
       else {
         filePrint(stderr,"in Real Domain ...\n");
         GenMultiDomainStatic<double> FAProb(domain);
         StaticSolver<double, GenMDDynamMat<double>, GenDistrVector<double>,
                      GenMultiDomainPostProcessor<double>, GenMultiDomainStatic<double>,
                      GenDistrVector<complex<double> > >
              FASolver(&FAProb);
         FASolver.solve();
       }
     } break;
     case SolverInfo::NonLinDynam: {
       if(domain->solInfo().newmarkBeta == 0) { // explicit
         MultiDomainDynam dynamProb(domain);
         DynamicSolver < MDDynamMat, DistrVector, MultiDomDynPostProcessor,
               MultiDomainDynam, double > dynamSolver(&dynamProb);
         dynamSolver.solve();
       }
       else {
         MDNLDynamic nldynamic(domain);
         NLDynamSolver <ParallelSolver, DistrVector, MultiDomainPostProcessor,
                        MDNLDynamic, DistrGeomState> nldynamicSolver(&nldynamic);
         nldynamicSolver.solve();
       }
     } break;
     default:
       filePrint(stderr,"*** ERROR: Problem type %d is not supported multi-domain mode\n", domain->probType());
   }

   totalMemoryUsed = double(memoryUsed()+totMemSpooles+totMemMumps)/oneMegaByte;//CBM
   delete threadManager;

 }
 else {

   //--------------------------------
#ifdef STRUCTOPT
   // print what type of problem we are doing
   if (dynamic_cast<Domain_opt*>(domain)->getStructoptFlag()) {
     if (domain->solInfo().aeroFlag > 0) {
       fprintf(stderr," ... Multiphysics structural opt: %d ...\n",
	       domain->solInfo().aeroFlag);
     }
     else {
       fprintf(stderr," ... Structural optimization: %d	...\n",
	       dynamic_cast<Domain_opt*>(domain)->getStructoptFlag());
     }
   }
   else if (dynamic_cast<Domain_opt*>(domain)->getReliabilityFlag() > 0) {
     fprintf(stderr," ... Structural reliability: %d    ...\n",
	     dynamic_cast<Domain_opt*>(domain)->getReliabilityFlag());
   }
   /*
   else if (domain->solInfo().aeroFlag > 0) {
     fprintf(stderr," ... Multiphysics structural: %d   ...\n",
	     domain->solInfo().aeroFlag);
   }
   else if (domain->solInfo().eigenEmFlag > 0) {
     fprintf(stderr," ... MP Structural eigenvalue: %d  ...\n",
	     domain->solInfo().eigenEmFlag);
   }
   else if (domain->solInfo().podBasEmFlag > 0) {
     fprintf(stderr," ... MP Structural POD Decomp: %d  ...\n",
	     domain->solInfo().podBasEmFlag);
   }
   else if (domain->solInfo().elecstrcFlag > 0) {
     fprintf(stderr," ... Multiphysics Electrostatic: %d...\n",
	     domain->solInfo().elecstrcFlag );
   }
   else if (domain->solInfo().optelecstrcFlag > 0) {
     fprintf(stderr," ... MP Electrostatic optimiz: %d  ...\n",
	     domain->solInfo().optelecstrcFlag);
   }
   else if (domain->solInfo().eigenEmFlag > 0) {
     fprintf(stderr," ... MP Structural eigenvalue: %d  ...\n",
	     domain->solInfo().eigenEmFlag);
   }
   else if (domain->solInfo().eigelecstrcFlag > 0) {
     fprintf(stderr," ... MP Electrostatic eigenval: %d ...\n",
	     domain->solInfo().eigelecstrcFlag);
   }
   else if (domain->solInfo().podelecstrcFlag > 0) {
     fprintf(stderr," ... MP Electrostatic POD Dec: %d  ...\n",
	     domain->solInfo().podelecstrcFlag);
   }
   else if (domain->solInfo().meshmotFlag > 0) {
     fprintf(stderr," ... Multiphysics Mesh-motion: %d  ...\n",
	     domain->solInfo().meshmotFlag );
   }
   else if (domain->solInfo().optmeshmotFlag > 0) {
     fprintf(stderr," ... MP Mesh-motion optimiz: %d    ...\n",
	     domain->solInfo().optmeshmotFlag);
   }
   else if (domain->solInfo().eigmeshmotFlag > 0) {
     fprintf(stderr," ... MP Mesh-motion eigenvalue: %d ...\n",
	     domain->solInfo().eigmeshmotFlag);
   }
   else if (domain->solInfo().podmeshmotFlag > 0) {
     fprintf(stderr," ... MP Mesh-motion POD Decomp: %d ...\n",
	     domain->solInfo().podmeshmotFlag);
   }
   else if (domain->solInfo().thermoeFlag > 0) {
     fprintf(stderr," ... MP Thermoelastic: %d ...\n",
	     domain->solInfo().thermoeFlag);
   }
   */

   if ( dynamic_cast<Domain_opt*>(domain)->getStructoptFlag()        > 0 ||
	dynamic_cast<Domain_opt*>(domain)->getReliabilityFlag()      > 0 /*||
	domain->solInfo().aeroFlag        > 0 ||
	domain->solInfo().eigenEmFlag     > 0 ||
	domain->solInfo().podBasEmFlag    > 0 ||
	domain->solInfo().elecstrcFlag    > 0 ||
	domain->solInfo().optelecstrcFlag > 0 ||
	domain->solInfo().eigelecstrcFlag > 0 ||
	domain->solInfo().podelecstrcFlag > 0 ||
	domain->solInfo().meshmotFlag     > 0 ||
	domain->solInfo().optmeshmotFlag  > 0 ||
	domain->solInfo().eigelecstrcFlag > 0 ||
	domain->solInfo().eigmeshmotFlag  > 0 ||
	domain->solInfo().podmeshmotFlag  > 0 ||
	domain->solInfo().thermoeFlag    == 2 */ )
     {
       dynamic_cast<Domain_opt*>(domain)->copyNodes();	//copy nodes; always needed
       if (dynamic_cast<Domain_opt*>(domain)->getStructoptFlag() > 0 )
	 { dynamic_cast<Domain_opt*>(domain)->structoptSolve(); }
       else
	 {
	   if(dynamic_cast<Domain_opt*>(domain)->getReliabilityFlag() > 0)
	     { dynamic_cast<Domain_opt*>(domain)->reliabilitySolve(); }
	 }
     /* Execute other types of analysis such as AEROELASTIC or MESHMOTION */
     }
   else
#endif
   switch(domain->probType()) {
     case SolverInfo::DisEnrM: {
        filePrint(stderr, " ... DEM Problem ...\n");
        DEM dem;
        dem.run(domain,geoSource);
        break;
     }
     case SolverInfo::Static:
       {
         if(geoSource->isShifted()) filePrint(stderr, " ... Frequency Response Helmholtz Analysis ");
         if(domain->isComplex()) {
           if(geoSource->isShifted()) filePrint(stderr, "in Complex Domain ...\n");
           SingleDomainStatic<DComplex, GenVector<DComplex>, GenSolver<DComplex> >
             statProb(domain);
           StaticSolver<DComplex, AllOps<DComplex>, /*GenSolver<DComplex>,*/ GenVector<DComplex>,
                        SingleDomainPostProcessor<DComplex, GenVector<DComplex>, GenSolver<DComplex> >,
                        SingleDomainStatic<DComplex, GenVector<DComplex>, GenSolver<DComplex> >, GenVector<DComplex> >
             statSolver(&statProb);
           statSolver.solve();
         }
         else {
           if(geoSource->isShifted()) filePrint(stderr, "in Real Domain ...\n");
           SingleDomainStatic<double, Vector, Solver> statProb(domain);
           StaticSolver<double, AllOps<double>, /*Solver,*/ Vector,
	     	        SingleDomainPostProcessor<double, Vector, Solver>,
		        SingleDomainStatic<double, Vector, Solver>, GenVector<DComplex> >
             statSolver(&statProb);
           statSolver.solve();
         }
       }
       break;

     case SolverInfo::Dynamic:
     //case SolverInfo::TempDynamic:
       {
        if(domain->solInfo().modal) {
	  fprintf(stderr," ... Modal Method  ...\n");
          ModalDescr<double> * modalProb = new ModalDescr<double>(domain);
          DynamicSolver<ModalOps, Vector, ModalDescr<double>, ModalDescr<double>, double>
              modalSolver(modalProb);
          modalSolver.solve();
        }
        else {
          if (domain->solInfo().activatePita) {
#ifdef DISTRIBUTED
             SingleDomainDynamic dynamProb(domain);
             Pita::LinearDriver::Ptr driver;
             if (!domain->solInfo().pitaTimeReversible) {
                 fprintf(stderr," ... Linear PITA ...\n");
                 driver = linearPitaDriverNew(&dynamProb);
             } else {
                 fprintf(stderr," ... Time-reversible linear PITA ...\n");
                 driver = linearReversiblePitaDriverNew(&dynamProb);
             }
             driver->solve();
#else
             fprintf(stderr," ... PITA requires distributed version ...\n");
#endif
	  } else {
            if (domain->solInfo().ATDARBFlag>=1.5) {
/*              SingleDomainDynamic dynamProb(domain);
              DynamicSolver <GenDynamMat<DComplex>, Vector,
                    SDDynamPostProcessor, SingleDomainDynamic, DComplex>
                dynaSolver(&dynamProb);
              dynaSolver.solve();
*/            }
            else {
              SingleDomainDynamic dynamProb(domain);
              DynamicSolver <GenDynamMat<double>, Vector,
                    SDDynamPostProcessor, SingleDomainDynamic, double>
                dynaSolver(&dynamProb);
              dynaSolver.solve();
            }
	  }
	}
      }
      break;
     case SolverInfo::TempDynamic:
       {
         SingleDomainTemp tempProb(domain);
         TempSolver<DynamMat, Vector, SDTempDynamPostProcessor,
                    SingleDomainTemp> dynaSolver(&tempProb);
         dynaSolver.solve();
       }
       break;
     case SolverInfo::Modal:
       { //CBM
 	 SingleDomainEigen eigenProb(domain);
         EigenSolver<DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> *eigenSolver;
          switch(domain->solInfo().eigenSolverType) {
            case SolverInfo::SubSpace :
              eigenSolver = new SubSpaceSolver <DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> (&eigenProb);
              break;
            case SolverInfo::LobPcg :
              eigenSolver = new LOBPCGSolver<DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> (&eigenProb);
              break;
            case SolverInfo::Arpack :
              {
#ifdef USE_ARPACK
		eigenSolver = new SymArpackSolver <DynamMat, Vector, VectorSet, SDEigenPostProcessor, SingleDomainEigen> (&eigenProb);
#else
                  filePrint(stderr," *** ERROR: executable not linked with ARPACK. See flag USE_ARPACK in Makefile.\n");
                  exit(-1);
#endif
              }
              break;
            default :
                filePrint(stderr," *** ERROR: unknown eigensolver.\n");
                exit(-1);
            }
	  eigenSolver->solve();
       }
       break;
     case SolverInfo::NonLinStatic:
       {
	 NonLinStatic nlstatic(domain);
         NLStaticSolver<Solver, Vector, SingleDomainPostProcessor<double, Vector, Solver>,
		        NonLinStatic, GeomState>
           nlsolver(&nlstatic);
	 nlsolver.solve();
       }
       break;
     case SolverInfo::MatNonLinStatic:
       {
         NLMatProbDesc nlstatic(domain);
         NLStaticSolver<Solver, Vector, NLMatProbDesc, NLMatProbDesc,
                        NLState, TotalUpdater<NLMatProbDesc, Vector, NLState> >
           nlsolver(&nlstatic);
         nlsolver.solve();
       }
       break;
     case SolverInfo::NonLinDynam: {
       if(domain->solInfo().modal) {
         NLModalDescr nlModalDescr(domain);
         NLDynamSolver< NLModalOpSolver, Vector, SDDynamPostProcessor, NLModalDescr,
           ModalGeomState, NLModalUpdater< NLModalDescr, Vector, ModalGeomState > >
           nlmodalsolver(&nlModalDescr);
         nlmodalsolver.solve();
       }
       else {
         if(domain->solInfo().activatePita) {
#ifdef DISTRIBUTED
           if (domain->solInfo().pitaTimeReversible) {
             filePrint(stderr, " ... Time-reversible nonlinear PITA ...\n");
             Pita::PitaNonLinDynamic pitaProblem(domain);
             Pita::NlDriver::Ptr pitaDriver = nlPitaDriverNew(&pitaProblem);
             pitaDriver->solve();
           } else {
             filePrint(stderr, " ... Nonlinear PITA ...\n");
             PitaNonLinDynamic pitaProblem(domain);
             NLDistrTimeDecompSolver pitaSolver(&pitaProblem);
             pitaSolver.solve();
           }
           fprintf(stderr, "End NlPita\n");
#else
           fprintf(stderr," ... PITA requires distributed version ...\n");
#endif
         }
         else {
           if(domain->solInfo().newmarkBeta == 0) { // explicit
             SingleDomainDynamic dynamProb(domain);
             DynamicSolver <GenDynamMat<double>, Vector,
                   SDDynamPostProcessor, SingleDomainDynamic, double>
               dynaSolver(&dynamProb);
             dynaSolver.solve();
           }
           else { // implicit
             NonLinDynamic nldynamic(domain);
             NLDynamSolver <Solver, Vector, SDDynamPostProcessor, NonLinDynamic, GeomState> nldynamicSolver(&nldynamic);
             nldynamicSolver.solve();
           }
         }
       }
     }
     break;
     case SolverInfo::MatNonLinDynam: {
         NLMatProbDesc nldynamic(domain);
         NLDynamSolver < Solver, Vector, NLMatProbDesc,
             NLMatProbDesc, NLState, TotalUpdater<NLMatProbDesc, Vector, NLState> >
                 nldynamicSolver(&nldynamic);
         nldynamicSolver.solve();
       }
       break;
     case SolverInfo::ArcLength:
	// KHP: Keep both methods around until satisfied with one or the other
	// KHP: because of difference in convergence criteria when doing
	// KHP: single domain vs. multiple domain.
        /*{
        NonLinStatic nlstatic(domain);
        NLStaticSolver <Solver,Vector,SingleDomainPostProcessor<double,Vector,Solver>,
                      NonLinStatic, GeomState  > nlsolver(&nlstatic);
        nlsolver.arclength();
        }*/
        domain->arcLength();
        break;
     case SolverInfo::ConditionNumber: {
        SingleDomainCond condProb(domain);
        CondSolver< DynamMat,SDCondPostProcessor, SingleDomainCond >
		condSolver(&condProb);
        condSolver.solve();
        }
        break;
     case SolverInfo::HelmholtzFreqSweep:
     case SolverInfo::Helmholtz:
       {
         filePrint(stderr, " ... Acoustic Helmholtz problem ");
         if(domain->isComplex() > 0) {
           filePrint(stderr, "in complex domain ...\n");
           SingleDomainStatic<DComplex, ComplexVector, ComplexSolver>
             FAProb(domain);
           StaticSolver<DComplex, AllOps<DComplex>, /*ComplexSolver,*/ ComplexVector,
                        SingleDomainPostProcessor<DComplex, ComplexVector, ComplexSolver>,
                        SingleDomainStatic<DComplex, ComplexVector, ComplexSolver>, ComplexVector >
             FASolver(&FAProb);
           FASolver.solve();
         }
         else {
           filePrint(stderr, "in real domain ...\n");
           SingleDomainStatic<double, Vector, Solver> FAProb(domain);
           StaticSolver<double, AllOps<double>, /*Solver,*/ Vector,
                        SingleDomainPostProcessor<double, Vector, Solver>,
                        SingleDomainStatic<double, Vector, Solver>, ComplexVector >
             FASolver(&FAProb);
           FASolver.solve();
         }
       }
       break;
     case SolverInfo::Top:
       {
         double mass = domain->computeStructureMass();
         fprintf(stderr," ... Structure mass = %10.4f    ...\n",mass);
	 domain->makeTopFile(topFlag);
       }
       break;
     case SolverInfo::AxiHelm:
       {
         FourierStatic fDescr(domain, fourHelmBC);
         FourierSolver fSolver(&fDescr);
         fSolver.solve();
       }
       break;
     case SolverInfo::None:
       {
         if(domain->solInfo().massFlag) {
           domain->preProcessing();
           double mass = domain->computeStructureMass();
           fprintf(stderr," ... Structure mass = %10.4f    ...\n",mass);
         }
         fprintf(stderr," ... No Analysis Type selected      ...\n");
       }
       break;
   }
   totalMemoryUsed = double(memoryUsed()+totMemSpooles+totMemMumps)/oneMegaByte;//CBM
 }

#ifdef DISTRIBUTED
// double totMem = (double) totalMemoryUsed;//CBM
// if(structCom) totMem = structCom->globalSum(totMem);
// totalMemoryUsed = (long) totMem;
 if(structCom) totalMemoryUsed = structCom->globalSum(totalMemoryUsed);
#endif

 domain->printStatistics(); // PJSA 4-2-08
 filePrint(stderr," --------------------------------------\n");
 filePrint(stderr," ... Total Time           = %.2e s\n",
         (getTime() - initTime)/1000.0);
 filePrint(stderr," ... Total Memory Used    = %.2e Mb\n", totalMemoryUsed);
 if(domain->solInfo().isNonLin() && domain->solInfo().newmarkBeta != 0.0)
   filePrint(stderr," ... Total Newton Iterations = %4d \n", totalNewtonIter);
 if(iterTotal > 0)
   filePrint(stderr," ... Total Krylov Iterations = %4d \n", iterTotal);
 filePrint(stderr," --------------------------------------\n");

 if(geoSource) { delete geoSource; geoSource = 0; }
 //if(communicator) { delete communicator; communicator = 0; }
 if(domain) { delete domain; domain = 0; }
#ifdef USE_MPI
 // if(fetiCom) { delete fetiCom; fetiCom = 0; }
 //delete syscom;
#endif
}

void
writeOptionsToScreen()
{
 fprintf(stderr,"\n*********** FEM Options **********************************************************\n");
 fprintf(stderr," -n [number]                   = specified number of threads (FETI)\n");
 fprintf(stderr," -d [decomposition_filename]   = specified decomposition file (FETI)\n");
 fprintf(stderr," -v [verbose_frequency]        = verbose is turned on; specified frequency\n"
	        "                                 for output to screen of FETI iteration count (FETI)\n"
                "                                 and subspace iteration convergence\n");
 fprintf(stderr," -c                            = contact status is outputted to screen (FETI)\n");
/*
 fprintf(stderr," -p                            = primal residual is outputted to screen (FETI)\n");
*/
 fprintf(stderr," -t                            = input file is converted to XPost format\n");
 fprintf(stderr," -T                            = input file is converted to XPost format;\n");
 fprintf(stderr,"                                 all numbering gaps are removed\n");
 fprintf(stderr," -m                            = input file is converted to XPost format;\n"
 	        "                                 each material is gathered in a separate element set\n");
 fprintf(stderr," -M                            = input file is converted to XPost format;\n"
 	        "                                 all numbering gaps are removed and each material is\n"
                "                                 gathered in a separate element set\n");
/*
 fprintf(stderr," -r                            = axisymmetric geometry contained in input file is\n"
                "                                 converted to XPost format\n");
*/
 fprintf(stderr," -P                            = Xpost patterns are automatically generated for the\n");
 fprintf(stderr,"                                 the various Xpost element sets; useful only\n");
 fprintf(stderr,"                                 in conjonction with the -m and -M options\n");
 fprintf(stderr,"                                 which can generate multiple Xpost element sets;\n");
 fprintf(stderr,"                                 automatically generates a global element set");
 fprintf(stderr,"\n*********** DOMDEC Options (FETI) *************************************************\n");
 fprintf(stderr," --dec                         = embedded DOMDEC module is applied to input file\n"
	        "                                 to generate a mesh decomposition\n");
 fprintf(stderr," --deter                       = optimization of mesh decomposition is performed using a\n"
	        "                                 deterministic algorithm\n"
	        "                                 to generate a mesh decomposition\n");
 fprintf(stderr," --exit                        = run is normally terminated after mesh decomposition\n"
	        "                                 is generated\n");
 fprintf(stderr," --nsub [number]               = specified number of subdomains\n");
 fprintf(stderr," --mem                         = memory estimate associated with generated mesh\n"
	        "                                 decomposition is outputted in file xxxxx.mem\n");
 fprintf(stderr," --load                        = load and load balance information associated with\n");

 fprintf(stderr," --sower                       = embedded SOWER module is applied to input file to\n"
	        "                                 generate binary distributed data\n");
 fprintf(stderr," --exit                        = run is normally terminated after binary distributed\n"
	        "                                 data is generated\n");

 fprintf(stderr," --nclus [number]              = specified number of clusters\n");

 fprintf(stderr,"************************************************************************************\n");
 exit(-1);


}

