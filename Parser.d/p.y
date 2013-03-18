%{
#include <iostream>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <map>
#include <cstdlib>
#include <Parser.d/AuxDefs.h>
#include <Driver.d/Domain.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#ifdef STRUCTOPT
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#endif

 int numColumns = 3;
 double amplitude = 1.0;
 int PitaTS = 1;         //CD: Pita
%}

%union
{
 int ival;
 double fval;
 char * strval;
 NumedNode nval;
 NumedElem eval;
 NumList nl;
 BCond bcval;
 BCList *bclist;
 LMPCons *lmpcons;
 MFTTData *ymtt;
 MFTTData *ctett;
 ComplexBCList *cxbclist;
 ComplexBCond cxbcval;
 FrameData frame;
 NodalFrameData nframe;
 MFTTData *mftval;
 MFTTData *mptval;
 MFTTData *hftval;
 LayerData ldata;
 LayInfo *linfo;
 CoefData coefdata;
 ComplexFDBC complexFDBC;
 ComplexFNBC complexFNBC;
 AxiMPC axiMPC;
 DoubleList dlist;
 SurfaceEntity* SurfObj;
 MortarHandler* MortarCondObj;
 LMPCTerm* mpcterm;
 GeoSource::Rprop rprop;
 OutputInfo oinfo;
 ConstraintOptions copt;
}

%expect 6

%token ACTUATORS AERO AEROH AEROTYPE AMAT ANALYSIS ARCLENGTH ATTRIBUTES ANGULAROUTTYPE
%token AUGMENT AUGMENTTYPE AVERAGED ATDARB ACOU ATDDNB ATDROB ARPACK ATDDIR ATDNEU
%token AXIHDIR AXIHNEU AXINUMMODES AXINUMSLICES AXIHSOMMER AXIMPC AUXCOARSESOLVER ACMECNTL ADDEDMASS AEROEMBED AUGMENTED
%token BLOCKDIAG BOFFSET BUCKLE BGTL BMPC BINARYINPUT BINARYOUTPUT
%token CHECKTOKEN COARSESOLVER COEF CFRAMES COLLOCATEDTYPE CONVECTION COMPOSITE CONDITION
%token CONTROL CORNER CORNERTYPE CURVE CCTTOL CCTSOLVER CRHS COUPLEDSCALE CONTACTSURFACES CMPC CNORM
%token COMPLEXOUTTYPE CONSTRMAT CASES CONSTRAINEDSURFACES CSFRAMES CSTYPE
%token CONWEP
%token DAMPING DblConstant DEM DIMASS DISP DIRECT DLAMBDA DP DYNAM DETER DECOMPOSE DECOMPFILE DMPC DEBUGCNTL DEBUGICNTL
%token CONSTRAINTS MULTIPLIERS PENALTY
%token EIGEN EFRAMES ELSCATTERER END ELHSOMMERFELD EXPLICIT EPSILON ELEMENTARYFUNCTIONTYPE
%token FABMAT FACOUSTICS FETI FETI2TYPE FETIPREC FFP FFPDIR FITALG FNAME FLUX FORCE FRONTAL FETIH FILTEREIG
%token FREQSWEEP FREQSWEEP1 FREQSWEEP2 FREQSWEEPA FSINTERFACE FSISCALING FSIELEMENT NOLOCALFSISPLITING FSICORNER FFIDEBUG FAILSAFE FRAMETYPE
%token GEPS GLOBALTOL GRAVITY GRBM GTGSOLVER GLOBALCRBMTOL GROUP GROUPTYPE GOLDFARBTOL GOLDFARBCHECK
%token HDIRICHLET HEAT HFETI HNEUMAN HSOMMERFELD HFTT
%token HELMHOLTZ HNBO HELMMF HELMSO HSCBO HWIBO HZEM HZEMFILTER HLMPC 
%token HERMITIAN HESSIAN
%token IACC IDENTITY IDIS IDIS6 IntConstant INTERFACELUMPED ITEMP ITERTYPE IVEL 
%token INCIDENCE IHDIRICHLET IHDSWEEP IHNEUMANN ISOLVERTYPE INPC INFINTY
%token JACOBI KRYLOVTYPE KIRLOC
%token LAYC LAYN LAYD LAYO LAYMAT LFACTOR LMPC LOAD LOBPCG LOCALSOLVER LINESEARCH LUMPED
%token MASS MATERIALS MATLAB MAXITR MAXORTHO MAXVEC MODAL MPCPRECNO MPCPRECNOID MPCTYPE MPCTYPEID MPCSCALING MPCELEMENT MPCBLOCKID 
%token MPCBLK_OVERLAP MFTT MPTT MRHS MPCCHECK MUMPSICNTL MUMPSCNTL MECH MODEFILTER MOMENTTYPE
%token NDTYPE NEIGPA NEWMARK NewLine NL NLMAT NLPREC NOCOARSE NODETOKEN NONINPC
%token NSBSPV NLTOL NUMCGM NOSECONDARY NFRAMES
%token OPTIMIZATION OUTPUT OUTPUT6 OUTPUTFRAME
%token QSTATIC QLOAD
%token PITA PITADISP6 PITAVEL6 NOFORCE MDPITA GLOBALBASES LOCALBASES TIMEREVERSIBLE REMOTECOARSE ORTHOPROJTOL READINITSEED JUMPCVG JUMPOUTPUT
%token PRECNO PRECONDITIONER PRELOAD PRESSURE PRINTMATLAB PROJ PIVOT PRECTYPE PRECTYPEID PICKANYCORNER PADEPIVOT PROPORTIONING PLOAD PADEPOLES POINTSOURCE PLANEWAVE PTOL PLANTOL PMAXIT
%token RADIATION RBMFILTER RBMSET READMODE REBUILD RENUM RENUMBERID REORTHO RESTART RECONS RECONSALG REBUILDCCT RANDOM RPROP RNORM REVERSENORMALS RIGID
%token SCALING SCALINGTYPE SENSORS SOLVERTYPE SHIFT
%token SPOOLESTAU SPOOLESSEED SPOOLESMAXSIZE SPOOLESMAXDOMAINSIZE SPOOLESMAXZEROS SPOOLESMSGLVL SPOOLESSCALE SPOOLESPIVOT SPOOLESRENUM SPARSEMAXSUP SPARSEDEFBLK
%token STATS STRESSID SUBSPACE SURFACE SAVEMEMCOARSE SPACEDIMENSION SCATTERER STAGTOL SCALED SWITCH STABLE SUBTYPE STEP SOWER SHELLTHICKNESS SURF SPRINGMAT
%token TANGENT TEMP TIME TOLEIG TOLFETI TOLJAC TOLPCG TOPFILE TOPOLOGY TRBM THERMOE THERMOH 
%token TETT TOLCGM TURKEL TIEDSURFACES THETA HRC THIRDNODE THERMMAT TDENFORC TESTULRICH THRU TOPFLAG
%token USE USERDEFINEDISP USERDEFINEFORCE UPROJ UNSYMMETRIC USING
%token VERSION WETCORNERS XPOST YMTT 
%token ZERO BINARY GEOMETRY DECOMPOSITION GLOBAL MATCHER CPUMAP
%token NODALCONTACT MODE FRIC GAP
%token OUTERLOOP EDGEWS WAVETYPE ORTHOTOL IMPE FREQ DPH WAVEMETHOD
%token MATSPEC MATUSAGE BILINEARPLASTIC FINITESTRAINPLASTIC LINEARELASTIC STVENANTKIRCHHOFF LINPLSTRESS READ OPTCTV ISOTROPICLINEARELASTIC NEOHOOKEAN ISOTROPICLINEARELASTICJ2PLASTIC ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS HYPERELASTIC MOONEYRIVLIN HENCKY LOGSTRAINPLASTIC SVKPLSTRESS
%token SURFACETOPOLOGY MORTARTIED MORTARSCALING MORTARINTEGRATIONRULE SEARCHTOL STDMORTAR DUALMORTAR WETINTERFACE
%token NSUBS EXITAFTERDEC SKIP OUTPUTMEMORY OUTPUTWEIGHT
%token WEIGHTLIST GMRESRESIDUAL 
%token SLOSH SLGRAV SLZEM SLZEMFILTER 
%token PDIR HEFSB HEFRS HEINTERFACE  // Added for HEV Problem, EC, 20080512
%token SNAPFI PODROB TRNVCT OFFSET ORTHOG SVDTOKEN CONVERSIONTOKEN CONVFI SAMPLING PODSIZEMAX REFSUBSTRACT TOLER

%type <complexFDBC> AxiHD
%type <complexFNBC> AxiHN
%type <axiMPC>   AxiLmpc
%type <bclist>   BCDataList IDisp6 TBCDataList PBCDataList AtdDirScatterer AtdNeuScatterer IDisp6Pita IVel6Pita
%type <bclist>   DirichletBC NeumanBC TempDirichletBC TempNeumanBC TempConvection TempRadiation ModalValList
%type <bclist>   HEVDirichletBC HEVDBCDataList HEVFRSBCList HEVFRSBC HEVFRSBCElem 
%type <bcval>    BC_Data TBC_Data ModalVal PBC_Data HEVDBC_Data
%type <coefdata> CoefList
%type <cxbcval>  ComplexBC_Data ComplexMPCHeader
%type <mpcterm>  MPCLine ComplexMPCLine
%type <cxbclist> ComplexBCDataList ComplexNeumanBC ComplexDirichletBC 
%type <frame>    Frame
%type <nframe>   NodalFrame
%type <fval>     Float DblConstant
%type <ival>     AEROTYPE Attributes AUGMENTTYPE AVERAGED 
%type <ival>     COLLOCATEDTYPE CORNERTYPE COMPLEXOUTTYPE TDENFORC CSTYPE ANGULAROUTTYPE
%type <ival>     ELEMENTARYFUNCTIONTYPE FETIPREC FETI2TYPE FRAMETYPE
%type <ival>     GTGSOLVER Integer IntConstant ITERTYPE
%type <ival>     RBMSET RENUMBERID OPTCTV
%type <rprop>    RPROP
%type <ival>     WAVETYPE WAVEMETHOD
%type <ival>     SCALINGTYPE SOLVERTYPE STRESSID SURFACE MOMENTTYPE
%type <ldata>    LayData LayoData LayMatData
%type <linfo>    LaycInfo LaynInfo LaydInfo LayoInfo
%type <mftval>   MFTTInfo
%type <mptval>   MPTTInfo
%type <hftval>   HFTTInfo
%type <ival>     NDTYPE OUTPUTFRAME
%type <ival>     GROUPTYPE
%type <nl>       NodeNums SommNodeNums 
%type <nval>     Node
%type <lmpcons>  MPCList ComplexMPCList MPCHeader
%type <strval>   FNAME 
%type <ymtt>     YMTTList
%type <ctett>    TETTList
%type <dlist>    FloatList
%type <SurfObj>  FaceSet
%type <MortarCondObj> MortarCondition TiedSurfaces ContactSurfaces
%type <ival>     MPCTYPEID MPCPRECNOID MPCBLOCKID
%type <ival>     ISOLVERTYPE RECONSALG
%type <ival>     PRECTYPEID SWITCH TOPFLAG
%type <oinfo>    OutInfo
%type <copt>     ConstraintOptionsData
%%
FinalizedData:
	All END
	 { 
          return 0;
         }
	;
All:
	Component
	| All Component
	;
Component:
	NodeSet 
        | DirichletBC
        { if(geoSource->setDirichlet($1->n,$1->d) < 0) return -1; delete $1; }
        | NeumanBC
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
        | LMPConstrain 
        | ComplexLMPConstrain 
	| ElemSet
	| FrameDList
        | NodalFrameDList
        | ConstrainedSurfaceFrameDList
	| Attributes
	{}
	| Materials
        | Statics
	| Pressure
	| Lumped
        {}
        | Preload
        {}
	| Renumbering
	| IDisp
	| Mode
        | IDisp6Pita
        | IVel6Pita
	| IDisp6
	{}
	| IVel
	| ITemp
	| SensorLocations
	| ActuatorLocations
	| UsddLocations
	| UsdfLocations
        | DynInfo
	| SloshInfo 
	| HEVibInfo 
        | QstaticInfo
	| Output
        | RenumberOutput
        | NLInfo
	| DiscrMasses
	| Composites
	| Cframes
        | LayMat
	| MFTTInfo
	  { domain->setMFTT($1); }
	| MPTTInfo
	  { domain->setMPTT($1); }
        | HFTTInfo
          { domain->setHFTT($1); }
        | YMTTable
        | TETTable
	| RbmTolerance
        | ToleranceInfo
        | Gravity
	| RbmFilterInfo
	| HzemFilterInfo
	| SlzemFilterInfo
        | ModeFilterInfo
        | Restart
	| LoadCInfo
	| UseCInfo
        | Control
        | Optimization
        | AnalysisInfo
        | OrthoInfo
//        | NewtonInfo
	| AeroInfo
        | AeroHeatInfo
        | ThermohInfo
        | ThermoeInfo
        | HzemInfo
        | SlzemInfo
	| CondInfo
	| MassInfo
	| TempDirichletBC
        { if(geoSource->setDirichlet($1->n,$1->d) < 0) return -1; }
        | TempNeumanBC
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
        | TempConvection
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
	| HEVDirichletBC
	{ if(geoSource->setDirichletFluid($1->n,$1->d) < 0) return -1; }
	| HEVFRSBC
	{ if(geoSource->setDirichletFluid($1->n,$1->d) < 0) return -1; }
	| TempRadiation
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
        | TopInfo
        | OldHelmInfo
        | DEMInfo
	| HelmInfo
        | HelmMFInfo
        | HelmSOInfo
        | Scatterer
        | FarFieldPattern
        | FarFieldPatternDirs
        | KirchhoffLocations
        | HelmHoltzBC
        | EleHelmHoltzBC
        | EleScatterer
        | NeumScatterer
        | WetScatterer
        | HelmScatterer
        | ComplexNeumanBC
        { if(domain->setComplexNeuman($1->n,$1->d) < 0) return -1; }
        | IComplexNeumannBC
        | ComplexDirichletBC
        { if(domain->setComplexDirichlet($1->n,$1->d) < 0) return -1; }
        | IComplexDirichletBC
        | IComplexDirichletBCSweep
        | AxiHDir
        | AxiHNeu
        | AxiNumModes
        | AxiNumSlices
        | AxiHSommer
        | AxiMPC
        | Sower
        | BinarySpec
        | AtdDirScatterer
        { if(geoSource->setDirichlet($1->n,$1->d) < 0) return -1; }
        | AtdNeuScatterer
        { if(geoSource->setNeuman($1->n,$1->d) < 0) return -1; }
	| AtdArbScatterer
	| AtdNeumScatterer
	| AtdRobinScatterer
        | Decompose
	| WeightList
	| NodalContact
	{}
        | ModeInfo
        | Noninpc
        | Inpc
        | Group
        | Random
	//| Shift
        | Impe
	| MatSpec
	| MatUsage
	| FaceSet
	{}
	| MortarCondition
	{}
        | TiedSurfaces
        {}
        | WetInterface
        {}
        | FSInterface
        {}
        | HEInterface
        {}
        | ContactSurfaces
        {}
        | ConstrainedSurfaces
        {}
	| BoffsetList
	| ParallelInTimeInfo 
        | AcmeControls
        | Constraints
	| SvdToken
	| Sampling
        | ConversionToken
        ;
Noninpc:
        NONINPC NewLine Integer Integer NewLine
          { domain->solInfo().noninpc = true;
            sfem->setOrder($3); 
            domain->solInfo().nsample = $4;
          }
         ;
Inpc:
        INPC NewLine Integer NewLine
          { domain->solInfo().inpc = true;
            sfem->setOrder($3);
          }
        ;
Group:
        GROUP NewLine
        | Group GROUPTYPE Integer Integer NewLine
        { if ($2 == OutputInfo::Attribute)  geoSource->setGroupAttribute($3-1, $4-1);
          else if ($2 == OutputInfo::Nodal)  geoSource->setNodeGroup($3-1, $4);
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", $2);  exit(-1); }
        }
        | Group GROUPTYPE Integer Integer Integer NewLine
        { int i;
          if ($2 == OutputInfo::Attribute)  {
            for(i=$3; i<$4+1; ++i)
              geoSource->setGroupAttribute(i-1,$5-1);
          }
          else if ($2 == OutputInfo::Nodal)  {
            for(i=$3; i<$4+1; ++i)
              geoSource->setNodeGroup(i-1, $5);
          }
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", $2);  exit(-1); }
        }
        | Group GROUPTYPE SURF Integer Integer NewLine
        { if ($2 == OutputInfo::Nodal) geoSource->setSurfaceGroup($4-1, $5);
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", $2);  exit(-1); }
        }
        ;
Random:
        RANDOM NewLine
        | Random Integer RPROP Float Float NewLine
        { geoSource->setGroupRandomProperty($2-1,$3,$4,$5); }
        ; 
Impe:
        IMPE NewLine FREQ Float NewLine
          { geoSource->setImpe($4); }
        //| IMPE NewLine SHIFT Float NewLine
        //  { geoSource->setShift($4); }
        | IMPE NewLine FREQSWEEP1 Float Float Integer NewLine
          { geoSource->setImpe($4); domain->addFrequencies1(2.0*PI*$4, 2.0*PI*$5, $6); }
        | IMPE NewLine FREQSWEEP2 Float Float Integer NewLine
          { geoSource->setImpe($4); domain->addFrequencies2(2.0*PI*$4, 2.0*PI*$5, $6); }
        | IMPE NewLine FREQSWEEP Float Float Integer Integer NewLine
          { geoSource->setImpe($4); domain->addFrequencies(2.0*PI*$4, 2.0*PI*$5, $6, $7); }
        | IMPE NewLine FREQSWEEPA Integer Integer Integer Float Float Float Integer Integer Integer NewLine
        {
          geoSource->setImpe($7);
          domain->addFrequencies($7, $8, 2,$5);
          domain->solInfo().isAdaptSweep = true;
          domain->solInfo().adaptSweep.maxP = $4;
          domain->solInfo().adaptSweep.numS = $5;
          domain->solInfo().adaptSweep.dgp_flag = bool($6);
          domain->solInfo().adaptSweep.w1 = 2.0*PI*$7;
          domain->solInfo().adaptSweep.w2 = 2.0*PI*$8;
          domain->solInfo().adaptSweep.atol = $9;
          domain->solInfo().adaptSweep.minRHS = $10;
          domain->solInfo().adaptSweep.maxRHS = $11;
          domain->solInfo().adaptSweep.deltaRHS = $12;
          domain->solInfo().nFreqSweepRHS = $11;
        }
        | IMPE NewLine FreqSweep 
        | Impe ReconsInfo
        | Impe DampInfo
        | Impe PadePivotInfo
        | Impe PadePolesInfo
        ;
PadePivotInfo:
        PADEPIVOT Float NewLine
        { domain->solInfo().pade_pivot = true; domain->solInfo().pade_tol = $2; }
        ;
PadePolesInfo:
        PADEPOLES NewLine
        { domain->solInfo().pade_poles = true; }
        | PADEPOLES Float Float NewLine
        { domain->solInfo().pade_poles = true; 
          domain->solInfo().pade_poles_sigmaL = $2; domain->solInfo().pade_poles_sigmaU = $3; }
FreqSweep:  
        /* FREQSWEEP f_0 */
        FREQSWEEP Float NewLine
        { geoSource->setImpe($2); domain->addCoarseFrequency(2.0*PI*$2); }
        /* f_i num_fine_intervals (i-1 to i) */
        | FreqSweep Float Integer NewLine
        { domain->addFrequencies(2.0*PI*$2, $3); }
        ;
ReconsInfo:
        RECONS RECONSALG Integer NewLine 
        { domain->solInfo().freqSweepMethod = $2; 
          int &l = domain->solInfo().padeL, &m = domain->solInfo().padeM, &n = domain->solInfo().padeN;
          switch($2) {
            case SolverInfo::Taylor:
              domain->solInfo().nFreqSweepRHS = $3+1; // taylor
              break;
            case SolverInfo::Pade1:
              n = 1;
              domain->solInfo().nFreqSweepRHS = l+m+1;
              break;
            case SolverInfo::Pade:
            case SolverInfo::Fourier:
              n = $3;
              domain->solInfo().nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SolverInfo::PadeLanczos:
              n = $3;
              if(m%n != 0) m = m/n*(n+1)-m%n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().nFreqSweepRHS = m/n;
              break;
            case SolverInfo::GalProjection:
              n = $3;
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::KrylovGalProjection:
              n = $3;
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::QRGalProjection:
              n = $3;
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
          }
        }
        | RECONS RECONSALG Integer Integer Integer NewLine  
        { domain->solInfo().freqSweepMethod = $2;
          int &l = domain->solInfo().padeL, &m = domain->solInfo().padeM, &n = domain->solInfo().padeN;
          switch($2) {
            case SolverInfo::Taylor:
              domain->solInfo().nFreqSweepRHS = $3+1; // taylor
              break;
            case SolverInfo::Pade1:
              n = 1;
              l = $4;
              m = $5;
              domain->solInfo().nFreqSweepRHS = l+m+1;
              break;
            case SolverInfo::Pade:
            case SolverInfo::Fourier:
              n = $3;
              l = $4; 
              m = $5;
              domain->solInfo().nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SolverInfo::PadeLanczos:
              n = $3;
              m = $5;
              if(m%n != 0) m = m/n*(n+1)-m%n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().nFreqSweepRHS = m/n;
              break;
            case SolverInfo::GalProjection:
              n = $3;
              l = $4;
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::KrylovGalProjection:
              n = $3;
              l = $4;
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::QRGalProjection:
              n = $3;
              l = $4;
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
          }
        }
        ;
Sower:
	SOWER NewLine
          { geoSource->binaryInput = true; geoSource->binaryOutput = true; }
        | BINARYINPUT SWITCH NewLine 
          { geoSource->binaryInput = bool($2); }
        | BINARYOUTPUT SWITCH NewLine
          { geoSource->binaryOutput = bool($2); }
        ;
BinarySpec:
        BINARY NewLine
        | BinarySpec GEOMETRY FNAME NewLine
          { geoSource->setGeo($3); }
        | BinarySpec DECOMPOSITION FNAME NewLine
          { geoSource->setDecomp($3); }
        | BinarySpec GLOBAL   FNAME NewLine
          { geoSource->setGlob($3); }
        | BinarySpec MATCHER  FNAME NewLine
          { geoSource->setMatch($3); }
        | BinarySpec CPUMAP  FNAME NewLine
          { geoSource->setCpuMap($3); }
        ;
AnalysisInfo:
        ANALYSIS Integer NewLine
        { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis($2); 
#endif
	}
        ;
Decompose :
       DECOMPOSE NewLine
         {if(decInit==0) decInit = new DecInit(); }
       | Decompose DECOMPFILE FNAME NewLine
         {decInit->file = strdup($3);}
       | Decompose NSUBS Integer NewLine
         {decInit->nsubs = $3; }
       | Decompose OUTPUTWEIGHT NewLine
         {decInit->weight = true; }
       | Decompose OUTPUTMEMORY NewLine
         {decInit->memory = true; }
       | Decompose EXITAFTERDEC NewLine
         {decInit->exitAfterDec = true;}
       | Decompose SKIP NewLine
         {decInit->skip = true;}
       | Decompose DETER NewLine
         {decInit->nosa = true; }
       ;
WeightList :
       WEIGHTLIST NewLine
        {}
       | WeightList Integer Float NewLine
         {
	   // map<int,double >::iterator it = weightList.find($2);
	   //if(it == weightList.end())
	     weightList[$2] = $3;
	 }
        ;
MFTTInfo:
	MFTT NewLine
	{ $$ = new MFTTData; }
	| MFTTInfo Float Float NewLine
	{ $$ = $1; $$->add($2,$3); }
	;
MPTTInfo:
        MPTT NewLine
        { $$ = new MFTTData; }
        | MPTTInfo Float Float NewLine
        { $$ = $1; $$->add($2,$3); }
        ;
HFTTInfo:
        HFTT NewLine
        { $$ = new MFTTData; }
        | HFTTInfo Float Float NewLine
        { $$ = $1; $$->add($2,$3); }
        ;
Composites:
	COMPOSITE NewLine
	| Composites CoefInfo
	| Composites LaycInfo
	| Composites LaynInfo
        | Composites LaydInfo
        | Composites LayoInfo
	;
Cframes:
	CFRAMES NewLine
	| Cframes Frame
	{ geoSource->addCFrame($2.num,$2.d); }
	;
CoefInfo:
	COEF Integer NewLine CoefList
	{ geoSource->addCoefInfo($2-1,$4); }
	;
CoefList:
	Integer Integer Float NewLine
	{ $$.zero(); $$.setCoef($1-1,$2-1,$3); }
	| CoefList Integer Integer Float NewLine
	{ $$.setCoef($2-1,$3-1,$4); }
	;
LaycInfo:
	LAYC Integer NewLine
	{ $$ = new LayInfo(0); geoSource->addLay($2-1,$$); }
	| LaycInfo LayData
	{ $1->add($2.lnum,$2.d,$2.matid); }
	;
LaynInfo:
	LAYN Integer NewLine 
	{ $$ = new LayInfo(1); geoSource->addLay($2-1,$$); }
	| LaynInfo LayData
	{ $1->add($2.lnum,$2.d,$2.matid); }
	;
LaydInfo:
        LAYD Integer NewLine
        { $$ = new LayInfo(0); geoSource->addLay($2-1,$$); }
        | LaydInfo LayoData
        { $1->add($2.lnum,$2.d,$2.matid); }
        ;
LayoInfo:
        LAYO Integer NewLine
        { $$ = new LayInfo(1); geoSource->addLay($2-1,$$); }
        | LayoInfo LayoData
        { $1->add($2.lnum,$2.d,$2.matid); }
        ;
LayData:
	Integer Float Float Float Float Float Float Float Float Float NewLine
	{ $$.lnum = $1-1;
          $$.matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
	  $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
	  $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10; }
        | Integer Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.lnum = $1-1;
          $$.matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10;
          $$.d[9] = $11;$$.d[10]= $12; }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.lnum = $1-1;
          $$.matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10;
          $$.d[9] = $11;$$.d[10]= $12; $$.d[11] = $13; } 
	;
LayoData:
        Integer Integer Float Float NewLine // PJSA 3-30-05: elastic constants to be read later from LAYMAT
        { $$.lnum = $1-1;  $$.matid = $2-1; $$.d[7] = $3; $$.d[8] = $4; }
        ;
LayMat:
        LAYMAT NewLine
        | LayMat LayMatData
        { geoSource->addLayMat($2.matid, $2.d); }
        ;
LayMatData:
        Integer Float Float Float Float Float NewLine // E1 E2 nu12 G12 rho
                                                      // note: coefficients of mutual influence are zero
                                                      //       coefficients of thermal expansion are zero
                                                      //       reference temperature is zero
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = 0.0; $$.d[5] = 0.0; $$.d[6] = $6; 
          $$.d[7] = 0; $$.d[8] = 0; $$.d[9] = 0; }
        | Integer Float Float Float Float Float Float Float NewLine // E1 E2 nu12 G12 mu1,12 mu2,12 rho
                                                                    // note: coefficients of thermal expansion are zero
                                                                    //       reference temperature is zero
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7; $$.d[6] = $8;
          $$.d[7] = 0; $$.d[8] = 0; $$.d[9] = 0; }
        | Integer Float Float Float Float Float Float Float Float Float NewLine // E1 E2 nu12 G12 mu1,12 mu2,12 rho, cte1, cte2
                                                                                // note: reference temperature is zero
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7; $$.d[6] = $8;
          $$.d[7] = $9; $$.d[8] = $10; $$.d[9] = 0; }
        | Integer Float Float Float Float Float Float Float Float Float Float NewLine // E1 E2 nu12 G12 mu1,12 mu2,12 rho, cte1, cte2, ta
        { $$.matid = $1-1; $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7; $$.d[6] = $8; 
          $$.d[7] = $9; $$.d[8] = $10; $$.d[9] = $11; }
        ;
DiscrMasses:
	DIMASS NewLine
	| DiscrMasses Integer Integer Float NewLine
	  { domain->addDMass($2-1,$3-1,$4); }
        | DiscrMasses Integer Integer Integer Float NewLine
          { domain->addDMass($2-1,$3-1,$5,$4-1); }
	;
Gravity:
	GRAVITY NewLine
	| Gravity Float Float Float NewLine
	  { domain->setGravity($2,$3,$4); }
	;
Restart:
        RESTART NewLine
        | Restart FNAME FNAME FNAME NewLine
        { geoSource->getCheckFileInfo()->lastRestartFile = $2;
          geoSource->getCheckFileInfo()->outputExt = $3;
          geoSource->getCheckFileInfo()->FlagRST = $4; }
        | Restart FNAME FNAME NewLine
        { geoSource->getCheckFileInfo()->lastRestartFile = $2;
          geoSource->getCheckFileInfo()->outputExt = $3;}
        | Restart FNAME Integer NewLine
        { geoSource->getCheckFileInfo()->currentRestartFile = $2;
          domain->solInfo().nRestart = $3; }
        ;
LoadCInfo:
	LOAD FNAME NewLine
       { geoSource->setControlFile($2);
         geoSource->setControlRoutine((char *) "controlObj");}
	;
UseCInfo:
	USE FNAME NewLine
	{ geoSource->setControlRoutine($2); }
	;
SensorLocations:
        SENSORS NewLine BCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations($3->n,$3->d) < 0) return -1; }
	;
ActuatorLocations:
        ACTUATORS NewLine BCDataList
        { for(int i=0; i<$3->n; ++i) { $3->d[i].type = BCond::Actuators; $3->d[i].caseid = 0; }
          if(geoSource->setActuatorLocations($3->n,$3->d) < 0) return -1; 
          if(geoSource->setNeuman($3->n,$3->d) < 0)            return -1; }
	;
UsdfLocations:
	USERDEFINEFORCE NewLine BCDataList
        { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<$3->n; ++i) { $3->d[i].type = BCond::Usdf; $3->d[i].caseid = 0; }
          if(geoSource->setUsdfLocation($3->n,$3->d) < 0) return -1;
          if(geoSource->setNeuman($3->n,$3->d) < 0)       return -1; } 
	;
UsddLocations:
	USERDEFINEDISP NewLine BCDataList
	{ geoSource->binaryInputControlLeft = true;
          for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Usdd;
          if(geoSource->setUsddLocation($3->n,$3->d) < 0) return -1;
          if(geoSource->setDirichlet($3->n,$3->d) < 0)    return -1; }
	;
Output:
	OUTPUT NewLine
        { numColumns = 3; } // set number of output columns to 3 
	| OUTPUT6 NewLine
	{ numColumns = 6; } // set number of output columns to 6 
        | OUTPUT Integer NewLine
        { numColumns = 3; geoSource->setOutLimit($2); }
        | OUTPUT6 Integer NewLine
        { numColumns = 6; geoSource->setOutLimit($2); }
        | Output OutInfo NewLine
        { $2.finalize(numColumns); geoSource->addOutput($2); }
        ;
OutInfo:
          STRESSID FNAME Integer // unformatted output for all nodes
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.filename = $2; $$.interval = $3; }
        | STRESSID Integer Integer FNAME Integer // formatted output for all nodes
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; }
        | STRESSID FNAME Integer Integer // unformatted output for one node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.filename = $2; $$.interval = $3; $$.nodeNumber = $4-1; }
        | STRESSID FNAME Integer GROUPTYPE Integer // TDL: unformatted output for node group or node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.filename = $2; $$.interval = $3; 
          if ($4 == OutputInfo::NodeGroup) $$.groupNumber = $5; else $$.nodeNumber = $5-1;}
        | STRESSID Integer Integer FNAME Integer Integer // formatted output for one node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; $$.nodeNumber = $6-1; }
        | STRESSID Integer Integer FNAME Integer GROUPTYPE Integer // TDL: formatted output for node group or node
        { $$.initialize(); $$.type = (OutputInfo::Type) $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; if ($6 == OutputInfo::NodeGroup) $$.groupNumber = $7; else $$.nodeNumber = $7-1; }
        | TDENFORC FNAME Integer // unformatted output for all nodes (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.filename = $2; $$.interval = $3; }
        | TDENFORC Integer Integer FNAME Integer // formatted output for all nodes (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; }
        | TDENFORC FNAME Integer Integer // unformatted output for one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.filename = $2; $$.interval = $3; $$.nodeNumber = $4-1; }
        | TDENFORC FNAME Integer GROUPTYPE Integer // TDL: unformatted output for group or one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.filename = $2; $$.interval = $3; if ($4 == OutputInfo::NodeGroup) $$.groupNumber = $5; else $$.nodeNumber = $5-1; }
        | TDENFORC Integer Integer FNAME Integer Integer // formatted output for one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; $$.nodeNumber = $6-1; }
        | TDENFORC Integer Integer FNAME Integer GROUPTYPE Integer // TDL: formatted output for group or one node (for explicit dynamics tied/contact surfaces)
        { $$.initialize(); $$.type = OutputInfo::TDEnforcement; $$.tdenforc_var = $1; $$.width = $2; $$.precision = $3; $$.filename = $4; $$.interval = $5; if ($6 == OutputInfo::NodeGroup) $$.groupNumber = $7; else $$.nodeNumber = $7-1; }

        | OutInfo NODETOKEN Integer
        { $$.nodeNumber = $3-1; }
        | OutInfo SURFACE
        { $$.surface = $2; }
        | OutInfo Float Float
        { $$.ylayer = $2; $$.zlayer = $3; }
        | OutInfo AVERAGED
        { $$.averageFlg = $2; }
        | OutInfo COMPLEXOUTTYPE
        { $$.complexouttype = $2; }
        | OutInfo COMPLEXOUTTYPE Integer
        { $$.complexouttype = $2; $$.ncomplexout = $3; }
        | OutInfo ANGULAROUTTYPE
        { $$.angularouttype = $2; }
        | OutInfo NDTYPE
        { $$.ndtype = $2; }
        | OutInfo NDTYPE Integer
        { $$.ndtype = $2; sfem->setnsamp_out($3); }
        | OutInfo OUTPUTFRAME
        { $$.oframe = (OutputInfo::FrameType) $2; }
        | OutInfo MATLAB 
        { $$.matlab = true; }
        ;
RenumberOutput:
        XPOST TOPFLAG NewLine
        { domain->outFlag = $2; }
        ;
DynInfo:
        DynamInfo
	| EIGEN NewLine
	{ domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; }
	| EIGEN Integer NewLine
	{ domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = $2;}
	| NEIGPA Integer NewLine
	{ domain->solInfo().nEig = $2; }
	| SUBSPACE NewLine
	{ domain->solInfo().eigenSolverType = SolverInfo::SubSpace;}
        | SUBSPACE Integer Float Float NewLine
        { domain->solInfo().setSubSpaceInfo($2,$3,$4); }
	| NSBSPV Integer NewLine
	{ domain->solInfo().subspaceSize = $2;}
	| TOLEIG Float NewLine
	{ domain->solInfo().tolEig = $2; }
	| TOLJAC Float NewLine
	{ domain->solInfo().tolJac = $2; }
        | EXPLICIT NewLine
        { domain->solInfo().explicitK = true; }
        | SHIFT Float NewLine
        { geoSource->setShift($2); }
        | ARPACK NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack; }
        | ARPACK FNAME NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = $2; }
        | ARPACK FNAME Integer NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = $2; 
          domain->solInfo().arpack_mode = $3; }
        | ARPACK Float Float NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue($2, int($3)); }
        | ARPACK Float Float Integer NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues($2, $3, $4);}
        | FILTEREIG SWITCH NewLine
        { domain->solInfo().filtereig = bool($2); }
        | JACOBI Integer NewLine
        { domain->solInfo().eigenSolverSubType = $2; } 
        | LOBPCG NewLine
        { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;}
        | DynInfo MAXITR Integer NewLine
        { domain->solInfo().maxitEig = $3; }
        | TESTULRICH NewLine
          { domain->solInfo().test_ulrich = true; }
        | ADDEDMASS Integer NewLine
        { domain->solInfo().addedMass = $2; }
	;
SloshInfo:
        SLOSH NewLine
        { domain->solInfo().sloshing = 1; }
	| SLGRAV Float NewLine
	{ domain->setGravitySloshing($2); }
	;
MassInfo:
	MASS NewLine
        { domain->solInfo().massFlag = 1; }
	;
CondInfo:
	CONDITION NewLine Float Integer NewLine
        { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol($3, $4); }
	| CONDITION NewLine
        { domain->solInfo().setProbType(SolverInfo::ConditionNumber);}
	;
TopInfo:
	TOPFILE NewLine
	{ domain->solInfo().setProbType(SolverInfo::Top); }
	;
DynamInfo:
        DYNAM NewLine 
        { domain->solInfo().setProbType(SolverInfo::Dynamic); }
	| DynamInfo TimeIntegration
        | DynamInfo TimeInfo
        | DynamInfo DampInfo
        | DynamInfo MODAL NewLine
        { domain->solInfo().modal = true; }
        | DynamInfo STABLE Integer NewLine
        { domain->solInfo().stable = $3; }
        | DynamInfo STABLE SWITCH NewLine
        { domain->solInfo().stable = $3; }
        | DynamInfo STABLE Integer Float Float Integer Integer NewLine
        { domain->solInfo().stable = $3;
          domain->solInfo().stable_cfl = $4;
          domain->solInfo().stable_tol = $5;
          domain->solInfo().stable_maxit = $6;
          domain->solInfo().stable_freq = $7;
        }
        | DynamInfo IACC SWITCH NewLine
        { domain->solInfo().iacc_switch = bool($3); }
        | DynamInfo ZERO SWITCH NewLine
        { domain->solInfo().zeroRot = bool($3); }
        | DynamInfo NOSECONDARY NewLine
        { domain->solInfo().no_secondary = true; }
        | DynamInfo CHECKTOKEN NewLine
        { domain->solInfo().check_energy_balance = true; }
        | DynamInfo CHECKTOKEN Float Float NewLine
        { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = $3; 
          domain->solInfo().epsilon2 = $4; }
        | DynamInfo CONWEP Float Float Float Float Float NewLine
        { domain->solInfo().ConwepOnOff = true;
          // Note: chargeWeight must be entered in the units of mass of the problem, not units of force.
          BlastLoading::InputFileData.ConwepGlobalOnOff = true;
          BlastLoading::InputFileData.ExplosivePosition[0]    = $3;
          BlastLoading::InputFileData.ExplosivePosition[1]    = $4;
          BlastLoading::InputFileData.ExplosivePosition[2]    = $5;
          BlastLoading::InputFileData.ExplosiveDetonationTime = $7;
          BlastLoading::InputFileData.BlastType               = BlastLoading::BlastData::AirBurst; // ($7 == 0 ? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst);
          BlastLoading::InputFileData.ScaleLength             = 1.0;
          BlastLoading::InputFileData.ScaleTime               = 1.0;
          BlastLoading::InputFileData.ScaleMass               = 1.0;
          BlastLoading::InputFileData.ExplosiveWeight         = $6*2.2; // The 2.2 factor is to convert from kilograms to pounds force.
          BlastLoading::InputFileData.ExplosiveWeightCubeRoot = pow(BlastLoading::InputFileData.ExplosiveWeight,1.0/3.0);
        }
        ;
TimeIntegration:
        NEWMARK NewLine
        { domain->solInfo().timeIntegration = SolverInfo::Newmark; }
        | MECH NewmarkSecondOrder NewLine
        | ACOU NewmarkSecondOrder NewLine
        { domain->solInfo().acoustic = true; }
        | HEAT NewmarkFirstOrder NewLine
        ;
NewmarkSecondOrder:
        Float Float // beta  gamma
        { domain->solInfo().setNewmarkSecondOrderInfo($1,$2); }
        | Float Float Float Float // beta  gamma alphaf alpham
        { domain->solInfo().setNewmarkSecondOrderInfo($1,$2,$3,$4); }
        | Float // rhoinfty
        { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,$1); }
        | Float Float Float // beta  gamma coef
        { domain->solInfo().setNewmarkSecondOrderInfo($1,$2);
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = $3; }
        ;
NewmarkFirstOrder:
        Float // epsilon
        { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
          domain->solInfo().order = 1;
         }
          else 
          domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo($1); 
        }
        ;
QstaticInfo:
        QSTATIC NewLine
        { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; }
        | QstaticInfo QMechInfo
        | QstaticInfo QHeatInfo
        | QstaticInfo MODAL NewLine
        { domain->solInfo().modal = true; }
	;
QMechInfo:
        MECH Float Float Integer NewLine
        { domain->solInfo().setQuasistaticInfo($2, 0, $3, $4); }
        | MECH Float Float Integer Float NewLine
        { domain->solInfo().setQuasistaticInfo($2, 0, $3, $4, $5); }
	;
QHeatInfo:
        HEAT Float Float Float Integer NewLine
        { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo($2, $3, $4, $5); }
/* need to check this
        | MECH Float Float Integer Float Float NewLine
        { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo($2, $3, $4, $5, $6); }
*/
	;
AeroInfo:
        AERO AEROTYPE NewLine
        { domain->solInfo().setAero($2); 
          domain->solInfo().isCollocated = 0; }
	| AERO NewLine AEROTYPE NewLine
        { domain->solInfo().setAero($3); 
          domain->solInfo().isCollocated = 0; }
	| AERO NewLine AEROTYPE Float Float NewLine
        { domain->solInfo().setAero($3);
          domain->solInfo().isCollocated = 0;
          if($3 < 6 || $3 == 20) {
              domain->solInfo().alphas[0] = $4+$5;
              domain->solInfo().alphas[1] = -$5;
          }
        }
        | AERO NewLine AEROTYPE Float NewLine
        { domain->solInfo().setAero($3);
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = $4;
        }
	| AeroInfo COLLOCATEDTYPE NewLine
	{ domain->solInfo().isCollocated = $2; }
        | AeroInfo MATCHER FNAME NewLine
        { geoSource->setMatch($3); }
        | AeroInfo AeroEmbeddedSurfaceInfo NewLine
        {}
	;
AeroEmbeddedSurfaceInfo:
        AEROEMBED
        {}
        | AeroEmbeddedSurfaceInfo Integer
        { domain->AddAeroEmbedSurfaceId($2); }
	;
AeroHeatInfo:
        AEROH NewLine AEROTYPE Float Float NewLine
        { domain->solInfo().setAeroHeat($3, $4, $5); }
	;
ThermohInfo:
        THERMOH NewLine 
        { domain->solInfo().setThermoh(1); }
	;
ThermoeInfo:
        THERMOE NewLine 
        { domain->solInfo().setThermoe(1); }
	;
ModeInfo:
        MODE NewLine 
        { domain->solInfo().setModeDecomp(1); }
	;
HzemInfo:
        HZEM NewLine
        { domain->solInfo().hzemFlag=1; }
	;
SlzemInfo:
        SLZEM NewLine
        { domain->solInfo().slzemFlag=1; }
	;
RbmTolerance:
        TRBM NewLine Float NewLine
        { domain->solInfo().setTrbm($3); }
/* 
        | TRBM NewLine Float Float NewLine
        { domain->solInfo().setTrbm($3); }
*/
	;
ToleranceInfo:
        GRBM NewLine Float Float NewLine
        { domain->solInfo().setGrbm($3,$4); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");}
        | GRBM NewLine Float NewLine
        { domain->solInfo().setGrbm($3); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");}
        | GRBM NewLine 
        { domain->solInfo().setGrbm();
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");}
	;
ModeFilterInfo:
        MODEFILTER Integer NewLine
        { domain->solInfo().modeFilterFlag = $2; }
        | MODEFILTER NewLine
        { domain->solInfo().modeFilterFlag = 1; }
        ;
RbmFilterInfo:
	RBMFILTER Integer NewLine
	{ domain->solInfo().useRbmFilter($2); }
	| RBMFILTER NewLine
	{ domain->solInfo().useRbmFilter(1); }
        | RBMFILTER NewLine RbmList NewLine
	;
RbmList:
    Integer
    { if($1 < 1 || $1 > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", $1);
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[$1-1] = 1;
    }
    | RbmList Integer
    { if($2 < 1 || $2 > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", $2);
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[$2-1] = 1;
    }
    ;
HzemFilterInfo:
        HZEMFILTER NewLine
        { domain->solInfo().hzemFilterFlag=1; }
	;
SlzemFilterInfo:
        SLZEMFILTER NewLine
        { domain->solInfo().slzemFilterFlag=1; }
	;
TimeInfo:
	TIME Float Float Float NewLine
	{ domain->solInfo().setTimes($4,$3,$2); }
	;
ParallelInTimeInfo:
        PITA NewLine Integer Integer ParallelInTimeOptions
        {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime($3,$4,1);
        }
        |
        PITA NewLine Integer Integer Integer ParallelInTimeOptions
        {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime($3,$4,$5);
        }
        | MDPITA NewLine Integer Integer Integer Integer ParallelInTimeOptions
        {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime($3,$4,$5); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        }
        ;
ParallelInTimeOptions:
        NewLine ParallelInTimeKeyWord ParallelInTimeOptions
        | NewLine
        ;
ParallelInTimeKeyWord:    
        NOFORCE
        { domain->solInfo().pitaNoForce = true; }
        | GLOBALBASES Integer
        { domain->solInfo().pitaGlobalBasisImprovement = $2; }
        | LOCALBASES
        { domain->solInfo().pitaLocalBasisImprovement = 1; }
        | TIMEREVERSIBLE
        { domain->solInfo().pitaTimeReversible = true; }
        | REMOTECOARSE
        { domain->solInfo().pitaRemoteCoarse = true; }
        | ORTHOPROJTOL Float
        { domain->solInfo().pitaProjTol = $2; }
        | READINITSEED
        { domain->solInfo().pitaReadInitSeed = true; }
        | JUMPCVG
        { domain->solInfo().pitaJumpCvgRatio = 0.0; }
        | JUMPCVG Float
        { domain->solInfo().pitaJumpCvgRatio = $2; }
        | JUMPOUTPUT
        { domain->solInfo().pitaJumpMagnOutput = true; }
        ;
DampInfo:
	DAMPING Float Float NewLine
	{ domain->solInfo().setDamping($2,$3); }
	| DAMPING MODAL NewLine ModalValList
	{ if(geoSource->setModalDamping($4->n, $4->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
	;
ComplexDirichletBC:
	HDIRICHLET NewLine ComplexBCDataList
	{ $$ = $3; }
	;
IComplexDirichletBC:
        IHDIRICHLET NewLine Float Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, $3, $4, $5);
        }
        | IHDIRICHLET NewLine Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, $3, $4, 0.0);
        }
        ;
IComplexDirichletBCSweep:
        IHDSWEEP NewLine Integer Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections($3, $4, $5);
        }
        | IHDSWEEP Integer NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections($2,0.0,0.0,0.0);
        }
        | IComplexDirichletBCSweep DirectionVector
        ;
DirectionVector:
        Integer Float Float Float NewLine
        {
           domain->setWaveDirections($1, $2, $3, $4);
        }
        ;
IComplexNeumannBC:
        IHNEUMANN NewLine Float Float Float NewLine
        {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, $3, $4, $5);
        }
        ;
DirichletBC:
        DISP NewLine
        { $$ = new BCList; }
        | DirichletBC BC_Data
        { $2.type = BCond::Displacements; $$->add($2); }
        | DirichletBC Integer THRU Integer Integer Float NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6, BCond::Displacements); $$->add(bc); } }
        | DirichletBC Integer THRU Integer STEP Integer Integer Float NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8, BCond::Displacements); $$->add(bc); } }
        | DirichletBC SURF BC_Data
        { BCond *surf_bc = new BCond[1];
          surf_bc[0] = $3;
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
/* TODO | DirichletBC SURF Integer Integer Float USING ConstraintOptionsData NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1; surf_bc[0].dofnum = $4-1; surf_bc[0].val = $5;
          surf_bc[0].type = BCond::Lmpc;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
*/
        ;
ConstrainedSurfaces:
        CONSTRAINEDSURFACES NewLine
        /*                    surface_id type   attribute_id csframe_id */
        | ConstrainedSurfaces Integer    CSTYPE Integer      Integer      NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $2-1;
          surf_bc[0].type = (BCond::BCType) $3; //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = $4-1;
          surf_bc[0].val = $5-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
        }
HEVDirichletBC:
        PDIR NewLine HEVDBCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Pdir; $$ = $3; }
	;
HEVDBCDataList:
	HEVDBC_Data
	{ $$ = new BCList; $$->add($1); }
	| HEVDBCDataList HEVDBC_Data
	{ $$ = $1; $$->add($2); }
	;
HEVDBC_Data:
	Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = 10; $$.val = $2; }
	| Integer NewLine
	{ $$.nnum = $1-1; $$.dofnum = 10; $$.val = 0.0; }
	;
HEVFRSBC:  // Added for HEV problem, EC, 20080512
        HEFRS NewLine HEVFRSBCList
        { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; 

          int* allPDirNodes = new int[$3->n];

          for (int ii=0; ii < $3->n; ii++)
            allPDirNodes[ii]=($3->d[ii]).nnum;
          sort(allPDirNodes,allPDirNodes + $3->n);

          int maxFSNodes = 32;
          int* allHEVFSNodes = new int[maxFSNodes];
          allHEVFSNodes[0] = allPDirNodes[0];

          int numHEVFSNodes = 1;
          for (int ii = 1; ii < $3->n; ++ii)  {
            if (numHEVFSNodes == maxFSNodes)  {
              int nMaxFSNodes = maxFSNodes*3/2;
              int* nAllHEVFSNodes = new int[nMaxFSNodes];
              for (int kk = 0; kk < maxFSNodes; kk++)
                nAllHEVFSNodes[kk] = allHEVFSNodes[kk];
              delete [] allHEVFSNodes;
              allHEVFSNodes = nAllHEVFSNodes;
              maxFSNodes = nMaxFSNodes;
            }

            if (allPDirNodes[ii]!=allPDirNodes[ii-1])
              allHEVFSNodes[numHEVFSNodes++] = allPDirNodes[ii];
          }
          delete [] allPDirNodes;
          
          $$ = new BCList;

          for (int ii=0; ii<numHEVFSNodes; ii++)
            $$->add(allHEVFSNodes[ii],10,0.0); 
          delete [] allHEVFSNodes;

          for(int i=0; i<$$->n; ++i) $$->d[i].type = BCond::Hefrs;
        }
        ;
HEVFRSBCList:
        HEVFRSBCElem
        { $$ = new BCList; 
          for (int ii = 0; ii < $1->n; ii++) 
           $$->add(($1->d)[ii]); }
        | HEVFRSBCList HEVFRSBCElem
        { $$ = $1; 
          for (int ii = 0; ii < $2->n; ii++) 
           $$->add(($2->d)[ii]); }
        ;
HEVFRSBCElem: // Added for HEV problem, EC, 20080512
        Integer Integer SommNodeNums NewLine
        { $$ = new BCList;
          //BCond *bclist = new BCond[$3.num];
          for(int i=0; i<$3.num; ++i) 
          { $$->add($3.nd[i],10,0.0); } }
          //geoSource->addDirichlet($3.num, bclist); }
        ;
TempDirichletBC:
        TEMP NewLine
        { $$ = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;}
        | TempDirichletBC Integer Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = $3; bc.type = BCond::Temperatures; $$->add(bc); }
        | TempDirichletBC SURF Integer Float NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1;
          surf_bc[0].val = $4;
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
	;
TempNeumanBC:
        FLUX NewLine
        { $$ = new BCList; }
        | FLUX Integer NewLine
        { $$ = new BCList($2); }
        | TempNeumanBC Integer Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = $3; bc.type = BCond::Flux; bc.caseid = $$->caseid; $$->add(bc); }
        | TempNeumanBC SURF Integer Float NewLine
        { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = $3-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = $4;
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].caseid = $$->caseid;
          geoSource->addSurfaceNeuman(1,surf_bc); }
	;
TempConvection:
        CONVECTION NewLine
        { $$ = new BCList; }
        | CONVECTION Integer NewLine
        { $$ = new BCList($2); }
        | TempConvection Integer Float Float Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = $3*$4*$5; bc.type = BCond::Convection; bc.caseid = $$->caseid; $$->add(bc); }
	;
TempRadiation:
        RADIATION NewLine
        { $$ = new BCList; }
        | RADIATION Integer NewLine
        { $$ = new BCList($2); }
        | TempRadiation Integer Float Float Float NewLine
        { $$ = $1; BCond bc; bc.nnum = $2-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*$3*$4*$5*$5*$5*$5; bc.type = BCond::Radiation; bc.caseid = $$->caseid; $$->add(bc); }
        ;
HelmHoltzBC:
        HSOMMERFELD NewLine SommerfeldBCDataList
	;
SommerfeldBCDataList:
        SommerfeldBC_Data
        | SommerfeldBCDataList SommerfeldBC_Data
	;
SommerfeldBC_Data:
        Integer Integer Float NewLine
          { domain->addSommer(new LineSommerBC($1-1, $2-1)); }
        | Integer Integer Integer Float NewLine
          { domain->addSommer(new TriangleSommerBC($1-1,$2-1,$3-1)); }
        | Integer Integer Integer Integer Float NewLine
          { domain->addSommer(new QuadSommerBC($1-1,$2-1,$3-1, $4-1)); }
	;
EleHelmHoltzBC:
        ELHSOMMERFELD NewLine SommerElement
        | EleHelmHoltzBC SommerElement
        ;
SommerElement:
        Integer Integer SommNodeNums NewLine
        { domain->addSommerElem($1-1, $2, 1.0, $3.num, $3.nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        }
        ;
SommNodeNums:
        Integer
        { $$.num = 1; $$.nd[0] = $1-1; }
        | SommNodeNums Integer
        { if($$.num == 64) return -1;
          $$.nd[$$.num] = $2-1; $$.num++; }
        ;
Scatterer:
        SCATTERER NewLine ScattererEleList
        ;
ScattererEleList:
        ScattererEle
        | ScattererEleList ScattererEle
        ;
ScattererEle:
        Integer Integer Float NewLine
        { domain->addScatter(new LineSommerBC($1-1, $2-1));
          domain->addNeum(new LineSommerBC($1-1, $2-1)); }
        | Integer Integer Integer Float NewLine
        { domain->addScatter(new TriangleSommerBC($1-1,$2-1,$3-1));
          domain->addNeum(new TriangleSommerBC($1-1,$2-1,$3-1)); }
        | Integer Integer Integer Integer Float NewLine
        { domain->addScatter(new QuadSommerBC($1-1,$2-1,$3-1, $4-1));
          domain->addNeum(new QuadSommerBC($1-1,$2-1,$3-1, $4-1)); }
        ;
EleScatterer:
        ELSCATTERER NewLine ScatterElement
        | EleScatterer ScatterElement
        ;
ScatterElement:
        Integer Integer SommNodeNums NewLine
        { domain->addScatterElem($1-1, $2, 1.0, $3.num, $3.nd);
          domain->addNeumElem($1-1, $2, 1.0, $3.num, $3.nd); }
        ;
NeumScatterer:
        HNBO NewLine NeumElement
        | NeumScatterer NeumElement
        ;
NeumElement:
        Integer Integer SommNodeNums NewLine
        { domain->addNeumElem($1-1, $2, 1.0, $3.num, $3.nd); }
        ;
WetScatterer:
        HWIBO NewLine WetInterfaceElement
        | WetScatterer WetInterfaceElement
        ;
WetInterfaceElement:
        Integer Integer SommNodeNums NewLine
        { domain->addWetElem($1-1, $2, 1.0, $3.num, $3.nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; }
        ;
HelmScatterer:
        HSCBO NewLine HelmScattererElement
        | HelmScatterer HelmScattererElement
        ;
HelmScattererElement:
        Integer Integer SommNodeNums NewLine
        { domain->addScatterElem($1-1, $2, 1.0, $3.num, $3.nd);}
        ;
PBC_Data:
	Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = 7; $$.val = $2; }
	;
PBCDataList:
	PBC_Data
	{ $$ = new BCList; $$->add($1); }
	| PBCDataList PBC_Data
	{ $$ = $1; $$->add($2); }
	;
AtdDirScatterer:
        ATDDIR NewLine PBCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Atddir; $$ = $3; }
        ;
AtdNeuScatterer:
        ATDNEU NewLine PBCDataList
        { for(int i=0; i<$3->n; ++i) { $3->d[i].type = BCond::Atdneu; $3->d[i].caseid = 0; } $$ = $3; }
        ;
AtdArbScatterer:
	ATDARB Float NewLine
        { domain->solInfo().ATDARBFlag = $2;}
        | AtdArbScatterer SommerElement
		;
AtdNeumScatterer:
	ATDDNB Float NewLine
        { domain->solInfo().ATDDNBVal = $2;}
        | AtdNeumScatterer NeumElement
        ;
AtdRobinScatterer:
	ATDROB Float Float Float NewLine
        { domain->solInfo().ATDROBVal = $2;
          domain->solInfo().ATDROBalpha = $3;
          domain->solInfo().ATDROBbeta = $4;}
        | AtdRobinScatterer HelmScattererElement
        ;
FarFieldPattern:
        FFP Integer NewLine
        { domain->setFFP($2); }
        | FFP Integer Integer NewLine
        { domain->setFFP($2,$3); }
        ;
FarFieldPatternDirs:
        FFPDIR Integer NewLine FFPDirList
        {
           domain->setFFP($2);
        }
        ;
AxiHDir:
        AXIHDIR NewLine Float Float NewLine Float Float Float NewLine
        { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex($3,$4));
          fourHelmBC->setDir($6, $7, $8);
        }
        | AxiHDir AxiHD
        { fourHelmBC->addDirichlet($2); }
	;
AxiHD:
        Integer NewLine
        { $$ = FDBC($1-1); }
	;
AxiHNeu:
        AXIHNEU NewLine Float Float NewLine Float Float Float NewLine
        { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex($3,$4));
          fourHelmBC->setDir($6, $7, $8);
        }
        | AxiHNeu AxiHN
        { fourHelmBC->addNeuman($2); }
	;
AxiHN:
        Integer Integer NewLine
        { $$ = FNBC($1-1, $2-1); }
        | Integer Integer Integer NewLine
        { $$ = FNBC($1-1, $2-1, $3-1); }
	;
AxiNumModes:
        AXINUMMODES Integer NewLine
        {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setModes($2);
          domain->solInfo().setProbType(SolverInfo::AxiHelm);
        }
	;
AxiNumSlices:
        AXINUMSLICES Integer NewLine
        {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSlices($2);
        }
	;
AxiHSommer:
        AXIHSOMMER NewLine Integer NewLine Float Float NewLine
        { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSomType($3);
          fourHelmBC->setSurf($5, $6);
        }
        | AxiHSommer AxiHSData
	;
AxiHSData:
        Integer Integer NewLine
        { fourHelmBC->addSommer(new LineAxiSommer($1-1, $2-1)); }
        | Integer Integer Integer NewLine
        { fourHelmBC->addSommer(new Line2AxiSommer($1-1, $2-1, $3-1)); }
	;
AxiMPC:
        AXIMPC NewLine
        { if( globalMPCs== NULL) globalMPCs = new MPCData(); }
        | AxiMPC AxiLmpc
        { globalMPCs->addMPC($2); }
	;
AxiLmpc:
        Integer Float Float Integer NewLine
        { $$ = MPC($1-1, $2, $3, $4, DComplex(1.0,0.0), 0.0, 0.0, 0.0); }
        | Integer Float Float Integer Float Float NewLine
        { $$ = MPC($1-1, $2, $3, $4, DComplex($5,$6) , 0.0, 0.0, 0.0); }
        | Integer Float Float Integer Float Float Float Float Float NewLine
        { $$ = MPC($1-1, $2, $3, $4, DComplex($5,$6) , $7, $8, $9); }
	;
Mode:
	READMODE FNAME NewLine
	{ domain->solInfo().readInROBorModes = $2;
	  domain->solInfo().readmodeCalled = true; }
	| READMODE FNAME Integer NewLine
	{ domain->solInfo().readInROBorModes = $2;
          domain->solInfo().readmodeCalled = true; 
 	  domain->solInfo().maxSizePodRom = $3; }	
	;
IDisp:
        IDIS NewLine
        { }
        | IDIS ZERO NewLine
        { domain->solInfo().zeroInitialDisp = 1; }
	| IDIS NewLine BCDataList
	{ for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis($3->n,$3->d) < 0) return -1; }
	| IDisp MODAL NewLine ModalValList
	{ for(int i=0; i<$4->n; ++i) $4->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal($4->n, $4->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
	;
IDisp6:
	IDIS6 Float NewLine
	{ $$ = new BCList; amplitude = $2;  }
	| IDIS6 NewLine
	{ $$ = new BCList; amplitude = 1.0; }
	| IDisp6 Integer Float Float Float Float Float Float NewLine
        { BCond bc; /* add 6 boundary conditions */
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = amplitude*$3; $$->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*$4; $$->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*$5; $$->add(bc);
                          bc.dofnum = 3; bc.val = amplitude*$6; $$->add(bc);
                          bc.dofnum = 4; bc.val = amplitude*$7; $$->add(bc);
                          bc.dofnum = 5; bc.val = amplitude*$8; $$->add(bc);
          for(int i=0; i<$$->n; ++i) $$->d[i].type = BCond::Idisp6;
          geoSource->setIDis6($$->n, $$->d);
        }
        | IDisp6 Integer Float Float Float NewLine
        { BCond bc; /* add 6 boundary conditions */
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = amplitude*$3; $$->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*$4; $$->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*$5; $$->add(bc);
                          bc.dofnum = 3; bc.val = 0.0         ; $$->add(bc);
                          bc.dofnum = 4; bc.val = 0.0         ; $$->add(bc);
                          bc.dofnum = 5; bc.val = 0.0         ; $$->add(bc);
          for(int i=0; i<$$->n; ++i) $$->d[i].type = BCond::Idisp6;
          geoSource->setIDis6($$->n, $$->d);
        }
	| GEPS NewLine
	{ fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); }
	| BUCKLE NewLine
	{ domain->solInfo().buckling = 1; }
	;
IDisp6Pita:
        /* PITA: Get initial seed displacement for each time-slice. */
        PITADISP6 Integer NewLine
        { $$ = new BCList; PitaTS = $2; }
        | IDisp6Pita Integer Float Float Float Float Float Float NewLine
        { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = $3; $$->add(bc);
                          bc.dofnum = 1; bc.val = $4; $$->add(bc);
                          bc.dofnum = 2; bc.val = $5; $$->add(bc);
                          bc.dofnum = 3; bc.val = $6; $$->add(bc);
                          bc.dofnum = 4; bc.val = $7; $$->add(bc);
                          bc.dofnum = 5; bc.val = $8; $$->add(bc);
          geoSource->setPitaIDis6($$->n, $$->d, PitaTS);
        } 
        ;
IVel6Pita:
        /* PITA: Get initial seed velocity for each time-slice. */
        PITAVEL6 Integer NewLine
        { $$ = new BCList; PitaTS = $2; }
        | IVel6Pita Integer Float Float Float Float Float Float NewLine
        { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = $2-1; bc.dofnum = 0; bc.val = $3; $$->add(bc);
                          bc.dofnum = 1; bc.val = $4; $$->add(bc);
                          bc.dofnum = 2; bc.val = $5; $$->add(bc);
                          bc.dofnum = 3; bc.val = $6; $$->add(bc);
                          bc.dofnum = 4; bc.val = $7; $$->add(bc);
                          bc.dofnum = 5; bc.val = $8; $$->add(bc);
          geoSource->setPitaIVel6($$->n, $$->d, PitaTS);
        }
        ;
IVel:
        IVEL NewLine
        { }
        | IVEL NewLine BCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel($3->n,$3->d) < 0) return -1; }
        | IVel MODAL NewLine ModalValList
        { for(int i=0; i<$4->n; ++i) $4->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal($4->n, $4->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
	;
ITemp:
        ITEMP NewLine TBCDataList
        { for(int i=0; i<$3->n; ++i) $3->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis($3->n,$3->d) < 0) return -1; }
	;
NeumanBC:
        FORCE NewLine
        { $$ = new BCList; }
        | FORCE Integer NewLine
        { $$ = new BCList($2); }
        | NeumanBC BC_Data
        { $2.type = BCond::Forces; $2.caseid = $$->caseid; $$->add($2); }
        | NeumanBC Integer THRU Integer Integer Float NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6, BCond::Forces, $$->caseid); $$->add(bc); } }
        | NeumanBC Integer THRU Integer STEP Integer Integer Float NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8, BCond::Forces, $$->caseid); $$->add(bc); } }
        | NeumanBC Integer THRU Integer Integer Float MOMENTTYPE NewLine
        { for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6, BCond::Forces, $$->caseid, (BCond::MomentType) $7); $$->add(bc); } }
        | NeumanBC Integer THRU Integer STEP Integer Integer Float MOMENTTYPE NewLine
        { for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8, BCond::Forces, $$->caseid, (BCond::MomentType) $7); $$->add(bc); } }
        | NeumanBC SURF BC_Data
        { BCond *surf_bc = new BCond[1];
          surf_bc[0] = $3;
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].caseid = $$->caseid;
          geoSource->addSurfaceNeuman(1,surf_bc); }
        ;
BCDataList:
	BC_Data
	{ $$ = new BCList; $$->add($1); }
	| BCDataList BC_Data
	{ $$ = $1; $$->add($2); }
        | Integer THRU Integer Integer Float NewLine
        { $$ = new BCList; for(int i=$1; i<=$3; ++i) { BCond bc; bc.setData(i-1, $4-1, $5); $$->add(bc); } }
        | BCDataList Integer THRU Integer Integer Float NewLine
        { $$ = $1; for(int i=$2; i<=$4; ++i) { BCond bc; bc.setData(i-1, $5-1, $6); $$->add(bc); } }
        | Integer THRU Integer STEP Integer Integer Float NewLine
        { $$ = new BCList; for(int i=$1; i<=$3; i+=$5) { BCond bc; bc.setData(i-1, $6-1, $7); $$->add(bc); } }
        | BCDataList Integer THRU Integer STEP Integer Integer Float NewLine
        { $$ = $1; for(int i=$2; i<=$4; i+=$6) { BCond bc; bc.setData(i-1, $7-1, $8); $$->add(bc); } }
	;
ModalValList:
	ModalVal
	{ $$ = new BCList; $$->add($1); }
	| ModalValList ModalVal
	{ $$ = $1; $$->add($2); }
	;
TBCDataList:
	TBC_Data
	{ $$ = new BCList; $$->add($1); }
	| TBCDataList TBC_Data
	{ $$ = $1; $$->add($2); }
	;
YMTTable:
        YMTT NewLine
        | YMTT NewLine YMTTList
        ;
YMTTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addYMTT($$);}
        | YMTTList Float Float NewLine
        { $$->add($2, $3); }
        | YMTTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addYMTT($$);}
        ;
TETTable:
        TETT NewLine
        | TETT NewLine TETTList
        ;
TETTList:
        CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($2); $$->add($4, $5); domain->addCTETT($$);}
        | TETTList Float Float NewLine
        { $$->add($2, $3); }
        | TETTList CURVE Integer NewLine Float Float NewLine
        { $$ = new MFTTData($3); $$->add($5, $6); domain->addCTETT($$);}
        ;
LMPConstrain:
        LMPC NewLine
        | LMPC NewLine MPCList
	;
MPCList:
        MPCHeader MPCLine
        { $$ = $1;
          $$->addterm($2);
          domain->addLMPC($$); }
        | MPCList MPCLine
        { $$->addterm($2); }
        | MPCList MPCHeader MPCLine
        { $$ = $2;
          $$->addterm($3);
          domain->addLMPC($$); }
	;
MPCHeader:
        Integer NewLine
        { $$ = new LMPCons($1, 0.0); 
          $$->setSource(mpc::Lmpc); }
        | Integer Float NewLine
        { $$ = new LMPCons($1, $2); 
          $$->setSource(mpc::Lmpc); }
        | Integer Float MODE Integer NewLine
        { $$ = new LMPCons($1, $2);
          $$->type = $4; 
          $$->setSource(mpc::Lmpc); }
        | Integer Float ConstraintOptionsData NewLine
        { $$ = new LMPCons($1, $2);
          $$->lagrangeMult = $3.lagrangeMult;
          $$->penalty = $3.penalty; 
          $$->setSource(mpc::Lmpc); }
	;
MPCLine:
        Integer Integer Float NewLine
        { if($3 == 0.0) {
            fprintf(stderr," *** WARNING: zero coefficient in LMPC\n");
            fprintf(stderr," ***          node %d dof %d\n",$1,$2);
          }
          $$ = new LMPCTerm();
          $$->nnum = $1-1;
          $$->dofnum = $2-1;
          $$->coef.r_value = $3;
        }
	;
ComplexLMPConstrain:
        HLMPC NewLine
        | HLMPC NewLine ComplexMPCList
        ;
ComplexMPCList:
        ComplexMPCHeader ComplexMPCLine
        { $$ = new LMPCons($1.nnum,$1.reval,$1.imval,$2); domain->addLMPC($$); }
        | ComplexMPCList ComplexMPCLine
        { $$->addterm($2); }
        | ComplexMPCList ComplexMPCHeader ComplexMPCLine
        { $$ = new LMPCons($2.nnum,$2.reval,$2.imval,$3); domain->addLMPC($$); }
        ;
ComplexMPCHeader:
        Integer CRHS Float Float NewLine // PJSA: added CRHS keyword to eliminate conflict
        { $$.nnum=$1; $$.reval=$3; $$.imval=$4; }
        | Integer Float NewLine
        { $$.nnum=$1; $$.reval=$2; $$.imval=0.0; }
        | Integer NewLine
        { $$.nnum=$1; $$.reval=0.0; $$.imval=0.0; }
        ;
ComplexMPCLine:
        Integer Integer Float Float NewLine
        { if(($3==0.0) && ($4==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",$1,$2);
          return -1;
          }
          else { $$ = new LMPCTerm(true); $$->nnum=($1-1); $$->dofnum=($2-1); $$->coef.c_value=DComplex($3, $4); }
        }
        | Integer Integer Float NewLine
        { if($3==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",$1,$2);
          return -1;
          }
          else { $$ = new LMPCTerm(true); $$->nnum=($1-1); $$->dofnum=($2-1); $$->coef.c_value=DComplex($3,0.0); }
        }
        ;
ComplexNeumanBC:
	HNEUMAN NewLine ComplexBCDataList
        { $$ = $3; }
        | HNEUMAN Integer NewLine ComplexBCDataList
        { for(int i=0; i<$4->n; ++i) $4->d[i].caseid = $2;
          $$ = $4; }
        ;
ComplexBCDataList:
	ComplexBC_Data
	{ $$ = new ComplexBCList; $$->add($1); }
        | ComplexBCDataList ComplexBC_Data
        { $$ = $1; $$->add($2); }
	;
Materials:
	MATERIALS NewLine MatData
	| Materials MatData
	;
MatData:
	Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
	{ StructProp sp; 
	  sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10; 
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float DAMPING Float Float NewLine
        { StructProp sp;
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.betaDamp = $17; sp.alphaDamp = $18;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float RIGID Integer Float NewLine
        { StructProp sp;
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.lagrangeMult = bool($17);
          sp.initialPenalty = sp.penalty = $18;
          sp.type = StructProp::Constraint;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
	| Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
	{ StructProp sp; 
	  sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10; 
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
	  sp.ymin = $16; sp.ymax = $17; sp.zmin = $18; sp.zmax = $19;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float DAMPING Float Float NewLine
        { StructProp sp;
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.ymin = $16; sp.ymax = $17; sp.zmin = $18; sp.zmax = $19;
          sp.betaDamp = $21; sp.alphaDamp = $22;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float RIGID Integer Float NewLine
        { StructProp sp;
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.Iyy = $14; sp.Izz = $15;
          sp.ymin = $16; sp.ymax = $17; sp.zmin = $18; sp.zmax = $19;
          sp.lagrangeMult = bool($21);
          sp.initialPenalty = sp.penalty = $22;
          sp.type = StructProp::Constraint;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.A = $2; sp.E = $3; sp.nu = $4; sp.rho = $5;
          sp.c = $6; sp.k = $7; sp.eh = $8;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp ); 
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13;  
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float DAMPING Float NewLine
        { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = $2;  sp.E = $3;  sp.nu  = $4;  sp.rho = $5;
          sp.c = $6;  sp.k = $7;  sp.eh  = $8;  sp.P   = $9;  sp.Ta  = $10;
          sp.Q = $11; sp.W = $12; sp.Ixx = $13; sp.betaDamp = $15;
          sp.isReal = true;
          geoSource->addMat( $1-1, sp );
        }
        | Integer AMAT Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,0.0);
          sp.fp.PMLtype = int($4);
          sp.fp.gamma = $5;
          sp.fp.Rx = $6;
          sp.fp.Sx = $7;
          sp.fp.Ry = $8;
          sp.fp.Sy = $9;
          sp.fp.Rz = $10;
          sp.fp.Sz = $11;
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,0.0);
          sp.rho = $4;
          sp.fp.PMLtype = int($5);
          sp.fp.gamma = $6;
          sp.fp.Rx = $7;
          sp.fp.Sx = $8;
          sp.fp.Ry = $9;
          sp.fp.Sy = $10;
          sp.fp.Rz = $11;
          sp.fp.Sz = $12;
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,$4);
          sp.rho = $5;
          sp.fp.PMLtype = int($6);
          sp.fp.gamma = $7;
          sp.fp.Rx = $8;
          sp.fp.Sx = $9;
          sp.fp.Ry = $10;
          sp.fp.Sy = $11;
          sp.fp.Rz = $12;
          sp.fp.Sz = $13;
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,0.0);
          sp.rho = $4;
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->solInfo().acoustic = true;
        }
        | Integer AMAT Float Float Float NewLine
        { StructProp sp;
          sp.soundSpeed = complex<double>($3,$4);
          sp.rho = $5;
          sp.type = StructProp::Fluid;
          geoSource->addMat( $1-1, sp );
          domain->solInfo().acoustic = true;
        }
	| Integer FABMAT Integer Float Float Float Float Float Float Float Float Float Integer Integer Integer NewLine
        { StructProp sp;
          sp.F_op = $3;
          sp.E = $4;
          sp.rho = $5;
	  sp.A = $6;
          sp.F_Uc = $7;
          sp.F_Uf = $8;
          sp.lambda = $9;
          sp.F_h = $10;
          sp.F_d = $11;
          sp.F_dlambda = $12;
          sp.F_np = $13;
          sp.F_Nf = $14;
	  sp.Seed = $15;
	  sp.isReal = true;
          sp.type = StructProp::Fabric;
          geoSource->addMat( $1-1, sp );
        }
        | Integer THERMMAT Float Float Float Float Float Float Float Float Float NewLine
        { StructProp sp; 
          sp.A = $3;  sp.rho = $4; sp.Q = $5; sp.c = $6; 
          sp.sigma = $7;  sp.k = $8;  sp.eh  = $9;  sp.P   = $10;  sp.Ta  = $11;
          sp.isReal = true;
          sp.type = StructProp::Thermal;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float NewLine
        { StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float Float Float NewLine
        { StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.amplitude = $5;
          sp.omega = $6;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float Float Float Float NewLine
        { StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.amplitude = $5;
          sp.omega = $6;
          sp.phase = $7;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float Float Float Float Float NewLine
        { StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.amplitude = $5;
          sp.omega = $6;
          sp.phase = $7;
          sp.offset = $8;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float Float Float Float Float Float Integer NewLine
        { StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.amplitude = $5;
          sp.omega = $6;
          sp.phase = $7;
          sp.B = $8;
          sp.C = $9;
          sp.relop = $10;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float ELEMENTARYFUNCTIONTYPE Float Float Float Float NewLine
        { // new style for joints with prescribed motion by 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.funtype = $5;
          sp.amplitude = $6;
          sp.offset = $7;
          sp.c1 = $8;
          sp.c2 = $9;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float ELEMENTARYFUNCTIONTYPE Float Float Float Float Float NewLine
        { // new style for joints with prescribed motion by 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.funtype = $5;
          sp.amplitude = $6;
          sp.offset = $7;
          sp.c1 = $8;
          sp.c2 = $9;
          sp.c3 = $10;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float ELEMENTARYFUNCTIONTYPE Float Float Float Float Float Float NewLine
        { // new style for joints with prescribed motion by 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.funtype = $5;
          sp.amplitude = $6;
          sp.offset = $7;
          sp.c1 = $8;
          sp.c2 = $9;
          sp.c3 = $10;
          sp.c4 = $11;
          sp.type = StructProp::Constraint;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float SPRINGMAT Float NewLine
        { // use for RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.type = StructProp::Constraint;
          sp.k1 = $6;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float SPRINGMAT Float Float NewLine
        { // use for UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.penalty = $4;
          sp.type = StructProp::Constraint;
          sp.k1 = $6;
          sp.k2 = $7;
          geoSource->addMat( $1-1, sp );
        }
        | Integer CONSTRMAT Integer Float SPRINGMAT Float Float Float NewLine
        { // use for SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool($3);
          sp.initialPenalty = sp.penalty = $4;
          sp.type = StructProp::Constraint;
          sp.k1 = $6;
          sp.k2 = $7;
          sp.k3 = $8;
          geoSource->addMat( $1-1, sp );
        }
        | Integer SPRINGMAT Float NewLine
        { // use for TorsionalSpringType1 or TranslationalSpring
          StructProp sp;
          sp.k1 = $3;
          geoSource->addMat( $1-1, sp );
        }
	;
ElemSet:
	TOPOLOGY NewLine Element
	| ElemSet Element
	;
FaceSet:
	SURFACETOPOLOGY Integer NewLine 
	{ if($2 == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetReverseNormals(false);
          domain->AddSurfaceEntity($$);
        }
        | SURFACETOPOLOGY Integer REVERSENORMALS NewLine 
        { if($2 == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetReverseNormals(true);
          domain->AddSurfaceEntity($$);
        }
        | SURFACETOPOLOGY Integer SHELLTHICKNESS Float NewLine 
        { if($2 == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetIsShellFace(true);
          $$->SetShellThickness($4);
          domain->AddSurfaceEntity($$);
        }
        | SURFACETOPOLOGY Integer SHELLTHICKNESS Float REVERSENORMALS NewLine 
        { if($2 == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          $$ = new SurfaceEntity($2);
          $$->SetIsShellFace(true);
          $$->SetShellThickness($4);
          $$->SetReverseNormals(true);
          domain->AddSurfaceEntity($$);
        }
        | FaceSet Integer Integer NodeNums NewLine 
        { if($$->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[$4.num];
            for(int i=0; i<$4.num; ++i) nodes[$4.num-1-i] = $4.nd[i];
            $$->AddFaceElement($2-1, $3, $4.num, nodes);
            delete [] nodes;
          }
/* moved to  Domain::SetUpSurfaces
          else if($$->GetIsShellFace()) { // for acme shell it is necessary to include both sides of the element in the face block
            int etype;
            if($3 == 1) etype = 5; // SHELLQUADFACEL4
            else if($3 == 3) etype = 6; // SHELLTRIFACEL3
            else { cerr << " *** ERROR: Surface element type " << $3 << " not supported with SHELL_THICKNESS option\n"; exit(-1); }
            $$->AddFaceElement(2*($2-1), etype, $4.num, $4.nd);
            int *nodes = new int[$4.num];
            for(int i=0; i<$4.num; ++i) nodes[$4.num-1-i] = $4.nd[i];
            $$->AddFaceElement(2*($2-1)+1, etype, $4.num, nodes);
            delete [] nodes;
          }
*/
          else $$->AddFaceElement($2-1, $3, $4.num, $4.nd);
        }
	; 
MortarCondition:
	MORTARTIED Integer Integer NewLine
	{ $$ = new MortarHandler($2, $3); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer STDMORTAR NewLine
	{ $$ = new MortarHandler($2, $3); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer DUALMORTAR NewLine
	{ $$ = new MortarHandler($2, $3); $$->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond($$); 
        }
	| MORTARTIED Integer Integer SEARCHTOL Float NewLine
	{ $$ = new MortarHandler($2, $3, $5); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer SEARCHTOL Float Float NewLine
	{ $$ = new MortarHandler($2, $3, $5, $6); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer STDMORTAR SEARCHTOL Float NewLine
	{ $$ = new MortarHandler($2, $3, $6); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer STDMORTAR SEARCHTOL Float Float NewLine
	{ $$ = new MortarHandler($2, $3, $6, $7); domain->AddMortarCond($$); }
	| MORTARTIED Integer Integer DUALMORTAR SEARCHTOL Float NewLine
	{ $$ = new MortarHandler($2, $3, $6); $$->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond($$); 
        }
	| MORTARTIED Integer Integer DUALMORTAR SEARCHTOL Float Float NewLine
	{ $$ = new MortarHandler($2, $3, $6, $7); $$->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond($$); 
        }
	;
WetInterface:
        WETINTERFACE Integer Integer NewLine
        { domain->addWetInterface($2, $3); domain->solInfo().isCoupled = true; }
        | WETINTERFACE Integer NewLine
        { domain->addWetInterface($2, $2); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; }
        ;
TiedSurfaces:
        // $2 = pair id, $3 = master surface, $4 = slave surface, $5 = mortar type, $6 = normal search tolerace, $7 = tangential search tolerance
        // $8 = number of iterations (TD enforcement), $9 = convergence tolerance (TD enforcement), $10 = friction coefficient
        TIEDSURFACES NewLine
        { }
        | TiedSurfaces Integer Integer Integer NewLine
        {
          $$ = new MortarHandler($3, $4);
          $$->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond($$); 
        }
        | TiedSurfaces Integer Integer Integer Integer NewLine
        { 
          $$ = new MortarHandler($3, $4); 
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          domain->AddMortarCond($$); 
        }
        | TiedSurfaces Integer Integer Integer Integer Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float Integer Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetTDEnfParams($8, $9);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float Integer Float Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetMortarType($5);
          $$->SetTDEnfParams($8, $9);
          $$->SetFrictionCoef($10);
          $$->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetConstraintOptions($5);
          domain->AddMortarCond($$);
        }
        | TiedSurfaces Integer Integer Integer Integer Float Float ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::TIED);
          $$->SetMortarType($5);
          $$->SetConstraintOptions($8);
          domain->AddMortarCond($$);
        }
        ;
FSInterface:
        FSINTERFACE NewLine
        { }
        | FSInterface Integer Integer Integer NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4); domain->solInfo().isCoupled = true; 
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        | FSInterface Integer Integer Integer Float Float NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4, $5, $6); domain->solInfo().isCoupled = true;
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        ;
HEVibInfo:
        HEFSB NewLine HEVInterfaceElement
        { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
        | HEVibInfo HEVInterfaceElement
        ;
HEVInterfaceElement:
        Integer Integer SommNodeNums NewLine
        { domain->addWetElem($1-1, $2, 1.0, $3.num, $3.nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
        ;
HEInterface:
        HEINTERFACE NewLine
        { }
        | HEInterface Integer Integer Integer NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4); domain->solInfo().HEV = 1;
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        | HEInterface Integer Integer Integer Float Float NewLine
        { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface($3, $4, $5, $6); domain->solInfo().HEV = 1;
          if($3 == $4) domain->solInfo().isMatching = true;
        }
        ;
ContactSurfaces:
        // $2 = pair id, $3 = master surface, $4 = slave surface, $5 = mortar type, $6 = normal search tolerace, $7 = tangential search tolerance
        // $8 = number of iterations (TD enforcement), $9 = convergence tolerance (TD enforcement), $10 = friction coefficient
        CONTACTSURFACES NewLine
        { }
        | ContactSurfaces Integer Integer Integer NewLine
        {
          $$ = new MortarHandler($3, $4); 
          $$->SetInteractionType(MortarHandler::CTC); 
          $$->SetMortarType(MortarHandler::STD); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer NewLine
        {
          $$ = new MortarHandler($3, $4); 
          $$->SetInteractionType(MortarHandler::CTC); 
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float Float NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float Float Integer Float NewLine
        { /* this one is for frictionless */
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($8, $9);
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float Float Integer Float Float NewLine
        { /* this one is for constant friction */
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($8, $9);
          $$->SetFrictionCoef($10);
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float Float Integer Float Float Float Float NewLine
        { /* this one is for velocity dependent friction */
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($8, $9);
          $$->SetFrictionCoef($10, $11, $12);
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float Float Integer Float Float Float Float Float NewLine
        { /* this one is for pressure dependent friction */
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetTDEnfParams($8, $9);
          $$->SetFrictionCoef($10, $11, $12, $13);
          $$->SetMortarType($5); 
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4); 
          $$->SetInteractionType(MortarHandler::CTC); 
          $$->SetMortarType(MortarHandler::STD); 
          $$->SetConstraintOptions($5);
          domain->AddMortarCond($$);
        }
        | ContactSurfaces Integer Integer Integer Integer Float Float ConstraintOptionsData NewLine
        {
          $$ = new MortarHandler($3, $4, $6, $7);
          $$->SetInteractionType(MortarHandler::CTC);
          $$->SetMortarType($5); 
          $$->SetConstraintOptions($8);
          domain->AddMortarCond($$);
        }
        ;
AcmeControls:
        ACMECNTL Integer NewLine
        { domain->solInfo().dist_acme = $2; }
        | FFIDEBUG Integer NewLine
        { domain->solInfo().ffi_debug = bool($2); }
        | MORTARSCALING Float NewLine
        { domain->solInfo().mortar_scaling = $2; }
        | MORTARINTEGRATIONRULE Integer NewLine
        { domain->solInfo().mortar_integration_rule = $2; }
NodeSet:
	NODETOKEN NewLine Node
	{ geoSource->addNode($3.num, $3.xyz, $3.cp, $3.cd); }
	| NodeSet Node
	{ geoSource->addNode($2.num, $2.xyz, $2.cp, $2.cd); }
	;
Node:
	Integer Float Float Float NewLine
	{ $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = $4;  $$.cp = 0;  $$.cd = 0; }
	| Integer Float Float NewLine
	{ $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = 0.0; $$.cp = 0;  $$.cd = 0; }
	| Integer Float NewLine
	{ $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = 0.0; $$.xyz[2] = 0.0; $$.cp = 0;  $$.cd = 0; }
        | Integer Float Float Float Integer Integer NewLine
        { $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = $4;  $$.cp = $5; $$.cd = $6;
          if($5 != 0) domain->solInfo().basicPosCoords = false;
          if($6 != 0) domain->solInfo().basicDofCoords = false; }
        | Integer Float Float Float Integer NewLine
        { $$.num = $1-1; $$.xyz[0] = $2; $$.xyz[1] = $3;  $$.xyz[2] = $4;  $$.cp = $5; $$.cd = $5;
          if($5 != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } }
	;
Element:
	Integer Integer NodeNums NewLine
	{ /* Define each Element */
          geoSource->addElem($1-1, $2, $3.num, $3.nd);}
	;
NodeNums:
	Integer
	{ $$.num = 1; $$.nd[0] = $1-1;}
	| NodeNums Integer
	{ if($$.num == 125) return -1; 
          $$.nd[$$.num] = $2-1; $$.num++;}
	;
BC_Data:
	Integer Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = $3; $$.mtype = BCond::Axial; }
	| Integer Integer NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = 0.0; $$.mtype = BCond::Axial; }
        | Integer Integer Float MOMENTTYPE NewLine
        { $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = $3; $$.mtype = (BCond::MomentType) $4; }
        | Integer Integer MOMENTTYPE NewLine
        { $$.nnum = $1-1; $$.dofnum = $2-1; $$.val = 0.0; $$.mtype = (BCond::MomentType) $3; }
	;
ModalVal:
	Integer Float NewLine
	{ $$.nnum = $1-1;  $$.dofnum = -1;  $$.val = $2; }
	;
TBC_Data:
	Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = 6; $$.val = $2; }
	;
ComplexBC_Data:
	Integer Integer Float Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.reval = $3; $$.imval = $4;  }
	| Integer Integer Float NewLine
	{ $$.nnum = $1-1; $$.dofnum = $2-1; $$.reval = $3; $$.imval = 0.0; }
	;
ConstrainedSurfaceFrameDList:
        CSFRAMES NewLine
        | ConstrainedSurfaceFrameDList Frame
        { geoSource->setCSFrame($2.num,$2.d); }
        ;
FrameDList:
        EFRAMES NewLine
	| FrameDList Frame
	{ geoSource->setFrame($2.num,$2.d); }
	;
Frame:
	Integer Float Float Float Float Float Float Float Float Float NewLine
	{ $$.num = $1-1; 
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10; }
        | Integer THIRDNODE Integer NewLine
        { $$.num = $1-1;
          geoSource->makeEframe($1-1, $3, $$.d); }
	;
NodalFrameDList:
        NFRAMES NewLine
        | NodalFrameDList NodalFrame
        { geoSource->setNodalFrame($2.id,$2.o,$2.d,$2.type); }
        ;
NodalFrame:
        Integer Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = NFrameData::Rectangular;
          $$.o[0] = 0;  $$.o[1] = 0;  $$.o[2] = 0;
          $$.d[0] = $2; $$.d[1] = $3; $$.d[2] = $4;
          $$.d[3] = $5; $$.d[4] = $6; $$.d[5] = $7;
          $$.d[6] = $8; $$.d[7] = $9; $$.d[8] = $10; }
        | Integer Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = NFrameData::Rectangular;
          $$.o[0] = $2;  $$.o[1] = $3;  $$.o[2] = $4;
          $$.d[0] = $5;  $$.d[1] = $6;  $$.d[2] = $7;
          $$.d[3] = $8;  $$.d[4] = $9;  $$.d[5] = $10;
          $$.d[6] = $11; $$.d[7] = $12; $$.d[8] = $13; }
        | Integer FRAMETYPE Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = $2;
          $$.o[0] = 0;  $$.o[1] = 0;   $$.o[2] = 0;
          $$.d[0] = $3; $$.d[1] = $4;  $$.d[2] = $5;
          $$.d[3] = $6; $$.d[4] = $7;  $$.d[5] = $8;
          $$.d[6] = $9; $$.d[7] = $10; $$.d[8] = $11; }
        | Integer FRAMETYPE Float Float Float Float Float Float Float Float Float Float Float Float NewLine
        { $$.id = $1;
          $$.type = $2;
          $$.o[0] = $3;  $$.o[1] = $4;  $$.o[2] = $5;
          $$.d[0] = $6;  $$.d[1] = $7;  $$.d[2] = $8;
          $$.d[3] = $9;  $$.d[4] = $10; $$.d[5] = $11;
          $$.d[6] = $12; $$.d[7] = $13; $$.d[8] = $14; }
        ;
BoffsetList:
	BOFFSET NewLine
	| BoffsetList Integer Integer Float Float Float NewLine
	{ OffsetData od;
	  od.first = $2-1; od.last = $3-1;
	  od.o[0] = $4; od.o[1] = $5; od.o[2] = $6; 
	  geoSource->addOffset(od); }
	;
Attributes:
	ATTRIBUTES NewLine 
        { $$ = 0; }
	| Attributes Integer Integer NewLine
	{ geoSource->setAttrib($2-1,$3-1); }
        | Attributes Integer Integer HRC Float NewLine // added HRC keyword for Hyper Reduction Coefficient
        { geoSource->setAttrib($2-1,$3-1); 
	  geoSource->setElementLumpingWeight($2 - 1, $5);
	  domain->solInfo().elemLumpPodRom = true; }
	| Attributes Integer Integer Integer Integer NewLine
	{ geoSource->setAttrib($2-1,$3-1,$4-1,$5-1); }
        | Attributes Integer Integer Integer Integer HRC Float NewLine // added HRC keyword for Hyper Reduction Coefficient
        { geoSource->setAttrib($2-1,$3-1,$4-1,$5-1);
	  geoSource->setElementLumpingWeight($2 - 1, $7); 
	  domain->solInfo().elemLumpPodRom = true; }
        | Attributes Integer Integer Integer THETA Float NewLine // PJSA: added THETA keyword to eliminate conflict
        { geoSource->setAttrib($2-1,$3-1,$4-1,-1,$6); }
	| Attributes Integer Integer IDENTITY NewLine
        { int i;
          for(i=$2; i<$3+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        }
	| Attributes Integer Integer Integer NewLine
	{ int i;
	  for(i=$2; i<$3+1; ++i)
 	    geoSource->setAttrib(i-1,$4-1);
	}
	| Attributes Integer Integer Integer Integer Integer NewLine
	{ int i;
	  for(i=$2; i<$3+1; ++i)
	    geoSource->setAttrib(i-1, $4-1, $5-1, $6-1);
	}
        | Attributes Integer Integer Integer Integer THETA Float NewLine // PJSA: added THETA keyword to eliminate conflict
        { int i;
          for(i=$2; i<$3+1; ++i)
            geoSource->setAttrib(i-1, $4-1, $5-1, -1, $7);
        }
	;
Pressure:
	PRESSURE NewLine
	| Pressure Integer Float NewLine
	{ geoSource->setElementPressure( $2-1, $3 ); }
        // Element-by-element Conwep OnOff flag:
        //| Pressure Integer Float String NewLine
        //{ geoSource->setElementPressure( $2-1, $3, ($4 == "On" ? true : false) ); }
	| Pressure Integer Integer Float NewLine
	{ int i; 
          for(i=$2; i<($3+1); ++i)
            geoSource->setElementPressure( i-1, $4 );  
        }
        | Pressure SURF PBC_Data
        { BCond *surf_pbc = new BCond[1];
          surf_pbc[0] = $3;
          geoSource->addSurfacePressure(1,surf_pbc); }
	;
Lumped:
	LUMPED NewLine
	{ geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        }
        | LUMPED Integer NewLine
        { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, $2);
          geoSource->setConsistentPFlag(false);
        }
	;
Preload:
        PRELOAD NewLine
        { }
        | Preload Integer Float NewLine
        { geoSource->setElementPreLoad( $2-1, $3 ); }
        | Preload Integer THRU Integer Float NewLine
        { int i;
          for(i=$2; i<($4+1); ++i)
            geoSource->setElementPreLoad( i-1, $5 );
        }
        | Preload Integer Float Float Float NewLine
        { double load[3] = { $3, $4, $5 };
          geoSource->setElementPreLoad( $2-1, load ); }
        | Preload Integer THRU Integer Float Float Float NewLine
        { double load[3] = { $5, $6, $7 };
          int i;
          for(i=$2; i<($4+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        }
	;
Statics:
        Solver
        | IterSolver
        | Statics CASES CasesList NewLine
        ;
CasesList:
        Integer
        { domain->solInfo().loadcases.push_back($1); }
        | CasesList Integer
        { domain->solInfo().loadcases.push_back($2); }
        ;
IterSolver:
        STATS NewLine ITERTYPE NewLine
        { domain->solInfo().type = 1;
          domain->solInfo().iterType = $3;
          domain->solInfo().setProbType(SolverInfo::Static); }
        | IterSolver PRECNO Integer NewLine
        { domain->solInfo().precond = $3; }
        | IterSolver MAXITR Integer NewLine
        { domain->solInfo().maxit = $3; }
        | IterSolver TOLPCG Float NewLine
        { domain->solInfo().tol = $3; }
        | IterSolver MAXORTHO Integer NewLine
        { domain->solInfo().maxvecsize = $3; }
        | IterSolver SUBTYPE Integer NewLine
        { domain->solInfo().iterSubtype = $3; }
        ;
Solver:
	STATS NewLine DIRECT NewLine
        { domain->solInfo().setSolver(0); 
          domain->solInfo().setProbType(SolverInfo::Static); }
	| STATS NewLine DIRECT Integer NewLine
	{ domain->solInfo().setSolver($4); 
          domain->solInfo().setProbType(SolverInfo::Static); }
	| STATS NewLine SOLVERTYPE NewLine
	{ domain->solInfo().setSolver($3); 
          domain->solInfo().setProbType(SolverInfo::Static); }
        | STATS NewLine SOLVERTYPE PIVOT NewLine
        { domain->solInfo().setSolver($3);
          domain->solInfo().setProbType(SolverInfo::Static);
          if($3 < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
        | STATS NewLine SOLVERTYPE UNSYMMETRIC NewLine
        { domain->solInfo().setSolver($3);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().getNLInfo().unsymmetric = true; }
        | STATS NewLine ITERTYPE Integer NewLine
        { domain->solInfo().setSolver($3,$4);    
          domain->solInfo().setProbType(SolverInfo::Static); }
        | STATS NewLine ITERTYPE Integer Float NewLine
        { domain->solInfo().setSolver($3,$4,$5);    
          domain->solInfo().setProbType(SolverInfo::Static); }
	| STATS NewLine ITERTYPE Integer Float Integer NewLine
	{ domain->solInfo().setSolver($3,$4,$5,$6); 
          domain->solInfo().setProbType(SolverInfo::Static); }
	| STATS NewLine ITERTYPE Integer Float Integer Integer NewLine
	{ domain->solInfo().setSolver($3,$4,$5,$6,$7); 
          domain->solInfo().setProbType(SolverInfo::Static); }
        | STATS NewLine ITERTYPE Integer Float Integer Integer Integer NewLine
        { domain->solInfo().setSolver($3,$4,$5,$6,$7,$8);
          domain->solInfo().setProbType(SolverInfo::Static); }
/*
        | STATS NewLine ITERTYPE NewLine PRECNO Integer NewLine TOLPCG Float NewLine MAXITR Integer NewLine
        { domain->solInfo().setSolver($6,$9,$12,$3,3);
          domain->solInfo().setProbType(SolverInfo::Static); }
*/
        | STATS NewLine FETI Integer Float NewLine
        { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = $4;
          domain->solInfo().fetiInfo.tol      = $5;
          domain->solInfo().fetiInfo.maxortho = $4;
          domain->solInfo().type =(2); }
        | STATS NewLine FETI Integer NewLine
        { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ($4-1); } 
	| STATS NewLine FETI DP NewLine
	{ domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; }
        | STATS NewLine FETI DPH NewLine
        { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; }
        | STATS NewLine FETI Integer FETI2TYPE NewLine
        { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ($4-1); 
          domain->solInfo().fetiInfo.feti2version 
                  = (FetiInfo::Feti2Version) $5; } 
        | STATS NewLine FETI Integer Float Integer NewLine
        { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = $4;
          domain->solInfo().fetiInfo.tol      = $5;
          domain->solInfo().fetiInfo.maxortho = $6;
          domain->solInfo().type =(2); }
        | FETI Integer Float Integer NewLine
        { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = $2;
          domain->solInfo().fetiInfo.tol      = $3;
          domain->solInfo().fetiInfo.maxortho = $4;
          domain->solInfo().type =(2); }
	| STATS NewLine FETI NewLine
	{ domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static); }
	| FETI NewLine
        { domain->solInfo().type =(2);}
	| FETI Integer NewLine
	{ domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ($2-1);}
	| FETI DP NewLine
	{ domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; }
        | FETI DPH NewLine
        { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; }
	| FETI Integer FETI2TYPE NewLine
	{ domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ($2-1);
          domain->solInfo().fetiInfo.feti2version = (FetiInfo::Feti2Version) $3; 
        }
	| BLOCKDIAG NewLine SOLVERTYPE NewLine
	{
	  domain->solInfo().type = 3;
          domain->solInfo().subtype = $3;
          domain->solInfo().getFetiInfo().solvertype = (FetiInfo::Solvertype)($3);
	}
        | SPARSEMAXSUP Integer NewLine
        { domain->solInfo().sparse_maxsup  = $2; }
        | SPARSEDEFBLK Integer NewLine
        { domain->solInfo().sparse_defblk  = $2; }
        | SPOOLESTAU Float NewLine
        { domain->solInfo().spooles_tau  = $2; }
        | SPOOLESMAXSIZE Integer NewLine
        { domain->solInfo().spooles_maxsize = $2; }
        | SPOOLESMAXDOMAINSIZE Integer NewLine
        { if($2 < 0) {
            $2 = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          domain->solInfo().spooles_maxdomainsize = $2; }
        | SPOOLESSEED Integer NewLine
        { domain->solInfo().spooles_seed = $2; }
        | SPOOLESMAXZEROS Float NewLine
        { if(($2 < 0.0) || ($2 > 1.0)) {
            $2 = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          domain->solInfo().spooles_maxzeros = $2; }
        | SPOOLESMSGLVL Integer NewLine
        { domain->solInfo().spooles_msglvl = $2; }
        | SPOOLESPIVOT SWITCH NewLine
        { domain->solInfo().pivot = bool($2); }
        | SPOOLESSCALE Integer NewLine
        { domain->solInfo().spooles_scale = $2; }
        | SPOOLESRENUM Integer NewLine
        { domain->solInfo().spooles_renum = $2; }
	| MUMPSICNTL Integer Integer NewLine
	{ domain->solInfo().mumps_icntl[$2] = $3; }
	| MUMPSCNTL Integer Float NewLine
	{ domain->solInfo().mumps_cntl[$2] = $3; }
        | GOLDFARBTOL Float NewLine
        { domain->solInfo().goldfarb_tol = $2; }
        | GOLDFARBCHECK SWITCH NewLine
        { domain->solInfo().goldfarb_check = bool($2); }
	| Solver MAXITR Integer NewLine 
	{ domain->solInfo().fetiInfo.maxit = $3; }
        | DEBUGICNTL Integer Integer NewLine
        { domain->solInfo().debug_icntl[$2] = $3; }
        | DEBUGCNTL Integer Float NewLine
        { domain->solInfo().debug_cntl[$2] = $3; }
/* potential conflict/confusion with LUMPED for mass matrix etc.
        | FETIPREC NewLine
        { domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) $1; }
*/
	| Solver PRECNO FETIPREC NewLine
        { domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) $3; }
        | Solver PRECNO LUMPED NewLine
        { domain->solInfo().fetiInfo.precno = FetiInfo::lumped; }
        | Solver PRECNO Integer NewLine
        { if(($3 < 0) || ($3 > 3)) { 
            $3 = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) $3;
	}
        | PRECTYPE Integer NewLine
        { if(($2 < 0) || ($2 > 1)) {
            $2 = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) $2;
        }
        | PRECTYPE PRECTYPEID NewLine
        { if(($2 < 0) || ($2 > 1)) {
            $2 = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) $2;
        }
        | TOLFETI Float NewLine
  	{ domain->solInfo().fetiInfo.tol = $2; }
        | TOLFETI Float Float NewLine
        { domain->solInfo().fetiInfo.tol = $2; 
          domain->solInfo().fetiInfo.absolute_tol = $3; }
        | STAGTOL Float NewLine
        { domain->solInfo().fetiInfo.stagnation_tol = $2; }
        | STAGTOL Float Float NewLine
        { domain->solInfo().fetiInfo.stagnation_tol = $2;
          domain->solInfo().fetiInfo.absolute_stagnation_tol = $3; }
        | PTOL Float Float NewLine
        { domain->solInfo().fetiInfo.primal_proj_tol = $2;
          domain->solInfo().fetiInfo.dual_proj_tol = $3; }
        | PMAXIT Integer Integer NewLine
        { domain->solInfo().fetiInfo.primal_plan_maxit = $2;
          domain->solInfo().fetiInfo.dual_plan_maxit = $3; }
        | PLANTOL Float Float NewLine
        { domain->solInfo().fetiInfo.primal_plan_tol = $2;
          domain->solInfo().fetiInfo.dual_plan_tol = $3; }
	| Solver MAXORTHO Integer NewLine
	{ domain->solInfo().fetiInfo.maxortho = $3; }
	| NOCOARSE NewLine
	{ domain->solInfo().fetiInfo.noCoarse = 1; }
	| PROJ Integer NewLine
	{ if($2 == 1) 
            domain->solInfo().fetiInfo.nonLocalQ = 0;
          else if($2 == 2) {
            domain->solInfo().fetiInfo.nonLocalQ = 1;
            if(domain->solInfo().fetiInfo.version == FetiInfo::feti2) {
              domain->solInfo().fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used"
                             " with FETI 2\n");
            }
          } else if($2 == 3) {
            domain->solInfo().fetiInfo.nonLocalQ = 1;
            domain->solInfo().fetiInfo.nQ = 3;
          } else if($2 == 4) {
            domain->solInfo().fetiInfo.nonLocalQ = 1;
            domain->solInfo().fetiInfo.nQ = 4;
          } else
            fprintf(stderr," *** WARNING: This projector does not exist,"
                           " using basic projector\n");
        }
	| SCALING Integer NewLine
	{ if(($2 < 0) || ($2 > 2)) $2 = 1; 
          domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) $2; }
	| SCALING SCALINGTYPE NewLine
	{ domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) $2; }
        | MPCSCALING Integer NewLine
        { if(($2 < 0) || ($2 > 2)) $2 = 2;
          domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) $2; }
        | MPCSCALING SCALINGTYPE NewLine
        { domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) $2; }
        | FSISCALING Integer NewLine
        { if(($2 < 0) || ($2 > 2)) $2 = 2;
          domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) $2; }
        | FSISCALING SCALINGTYPE NewLine
        { domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) $2; }
        | MPCELEMENT NewLine
        { domain->solInfo().fetiInfo.mpc_element = true; }
        | FSIELEMENT NewLine
        { domain->solInfo().fetiInfo.fsi_element = true; }
        | FSICORNER Integer NewLine
        { domain->solInfo().fetiInfo.fsi_corner = $2; }
        | NOLOCALFSISPLITING NewLine
        { domain->solInfo().fetiInfo.splitLocalFsi = false; } 
        | COUPLEDSCALE Float NewLine
        { domain->solInfo().coupled_scale = $2; }
        | WETCORNERS NewLine
        { domain->solInfo().fetiInfo.wetcorners = true; }
	| CORNER CORNERTYPE NewLine
        { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) $2; }
        | CORNER CORNERTYPE Integer NewLine
        { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) $2; 
          domain->solInfo().fetiInfo.pick_unsafe_corners = bool($3);
        }
        | CORNER AUGMENTTYPE NewLine
        { if($2 == 0) {
            domain->solInfo().fetiInfo.corners = FetiInfo::noCorners;
            domain->solInfo().fetiInfo.pickAnyCorner = 0; 
            domain->solInfo().fetiInfo.bmpc = true;
            domain->solInfo().fetiInfo.pick_unsafe_corners = false;
            domain->solInfo().fetiInfo.augment = FetiInfo::none;
          }
        }
        | AUGMENT AUGMENTTYPE NewLine
        {
          if(domain->solInfo().fetiInfo.dph_flag && ($2 == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          }
          else domain->solInfo().fetiInfo.augment = (FetiInfo::AugmentType) $2;

          if(domain->solInfo().fetiInfo.augment == FetiInfo::Edges) {
            domain->solInfo().fetiInfo.rbmType = FetiInfo::translation;
            domain->solInfo().fetiInfo.nGs = 3;
          } 
          else if(domain->solInfo().fetiInfo.augment == FetiInfo::WeightedEdges) {
            domain->solInfo().fetiInfo.rbmType = FetiInfo::translation;
            domain->solInfo().fetiInfo.nGs = 3;
          }
          else if(domain->solInfo().fetiInfo.augment == FetiInfo::Gs) {
            domain->solInfo().fetiInfo.rbmType = FetiInfo::all;
            domain->solInfo().fetiInfo.nGs = 6;
          }
        }
        | AUGMENT AUGMENTTYPE RBMSET NewLine
        {
          if(domain->solInfo().fetiInfo.dph_flag && ($2 == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          }
          else domain->solInfo().fetiInfo.augment = (FetiInfo::AugmentType) $2;
          if(domain->solInfo().fetiInfo.dph_flag && ($3 > 2) && ($3 < 6)) {
            std::cerr << "WARNING: Selected rbm type is unsupported for FETI-DPH, set to translation \n";
            domain->solInfo().fetiInfo.rbmType = FetiInfo::translation;
          }
          else domain->solInfo().fetiInfo.rbmType = (FetiInfo::RbmType) $3;

          if(domain->solInfo().fetiInfo.rbmType == FetiInfo::all)
            domain->solInfo().fetiInfo.nGs = 6;
          else if(domain->solInfo().fetiInfo.rbmType != FetiInfo::None)
            domain->solInfo().fetiInfo.nGs = 3;
        }
        | AUGMENT EDGEWS Integer NewLine
        { domain->solInfo().fetiInfo.numdir = $3; 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
        | AUGMENT EDGEWS WAVETYPE Integer NewLine
        { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) $3;
          domain->solInfo().fetiInfo.numdir = $4; 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
        | AUGMENT EDGEWS Integer WAVEMETHOD NewLine
        { domain->solInfo().fetiInfo.numdir = $3;
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) $4;
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
        | AUGMENT EDGEWS WAVETYPE Integer WAVEMETHOD NewLine
        { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) $3;
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) $5;
          domain->solInfo().fetiInfo.numdir = $4;
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
	| ORTHOTOL Float NewLine
	{ domain->solInfo().fetiInfo.orthotol = $2; }
        | ORTHOTOL Float Float NewLine
        { domain->solInfo().fetiInfo.orthotol = $2; 
          domain->solInfo().fetiInfo.orthotol2 = $3; }
	| GLOBALTOL Float NewLine
	{ domain->solInfo().fetiInfo.grbm_tol = $2; }
        | GLOBALCRBMTOL Float NewLine
	{ domain->solInfo().fetiInfo.crbm_tol = $2; }
        | CCTTOL Float NewLine
        { domain->solInfo().fetiInfo.cct_tol = $2; }
        | REBUILDCCT SWITCH NewLine
        { domain->solInfo().fetiInfo.rebuildcct = int($2); }
        | UPROJ Integer NewLine
        { domain->solInfo().fetiInfo.uproj = $2; }
	| PRINTMATLAB NewLine
	{ domain->solInfo().fetiInfo.printMatLab = 1; }
	| LOCALSOLVER SOLVERTYPE NewLine
	{ domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) $2; }
	| COARSESOLVER SOLVERTYPE NewLine
	{ domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) $2; }
        | AUXCOARSESOLVER SOLVERTYPE NewLine
        {  domain->solInfo().fetiInfo.auxCoarseSolver = (FetiInfo::Solvertype) $2; }
        | CCTSOLVER SOLVERTYPE NewLine
        { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) $2; }
        | LOCALSOLVER SOLVERTYPE PIVOT NewLine
        { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) $2;
          if($2 < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; } 
        | LOCALSOLVER SOLVERTYPE SCALED NewLine
        { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) $2;
          domain->solInfo().localScaled = true; }
        | COARSESOLVER SOLVERTYPE PIVOT NewLine
        { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) $2; 
          if($2 < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
        | AUXCOARSESOLVER SOLVERTYPE PIVOT NewLine
        { domain->solInfo().fetiInfo.auxCoarseSolver  = (FetiInfo::Solvertype) $2;
          if($2 < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
        | COARSESOLVER SOLVERTYPE SCALED NewLine
        { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) $2;
          domain->solInfo().coarseScaled = true; }
        | CCTSOLVER SOLVERTYPE PIVOT NewLine
        { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) $2; 
          if($2 < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
        | CCTSOLVER SOLVERTYPE SCALED NewLine
        { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) $2; 
          if($2!=0) fprintf(stderr," *** WARNING: Scaling not supported for this CCt solver \n");
          else domain->solInfo().fetiInfo.cctScaled = true; }
	| VERSION Integer NewLine
	{ 
          if($2 == 1)
            domain->solInfo().fetiInfo.version = FetiInfo::feti1; 
          else if($2 == 2) {
            domain->solInfo().fetiInfo.version = FetiInfo::feti2;
            if(domain->solInfo().fetiInfo.nonLocalQ == 1) {
              domain->solInfo().fetiInfo.nonLocalQ = 0;
              domain->solInfo().fetiInfo.nQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } else if($2 == 3) {
            domain->solInfo().fetiInfo.version = FetiInfo::feti3;
            if(domain->solInfo().fetiInfo.nonLocalQ == 1) {
              domain->solInfo().fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } 
          else {
            domain->solInfo().fetiInfo.version = FetiInfo::feti1;
	    fprintf(stderr," *** WARNING: Version does not exist,"
                           " using FETI 1\n");
          }
	}
	| GTGSOLVER NewLine
	{ domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) $1; }
        | GMRESRESIDUAL NewLine
        { domain->solInfo().fetiInfo.gmresResidual = true; }
        | GMRESRESIDUAL SWITCH NewLine
        { domain->solInfo().fetiInfo.gmresResidual = bool($2); }
        | PICKANYCORNER Integer NewLine
        { domain->solInfo().fetiInfo.pickAnyCorner = $2; }
/* deprecated
        | FRONTAL NewLine
	{ domain->solInfo().fetiInfo.solvertype = FetiInfo::frontal; }
*/
        | NLPREC NewLine
        { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().fetiInfo.nlPrecFlg = 1; 
          domain->solInfo().setKrylov(); 
        }
	| NLPREC Integer NewLine
        { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
	  domain->solInfo().fetiInfo.nlPrecFlg = $2;
	  domain->solInfo().setKrylov();
	}
        | HFETI Integer Float NewLine
        { domain->solInfo().fetiInfo.maxit = $2;
          domain->solInfo().fetiInfo.tol = $3;
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
        | HFETI Integer Float Integer NewLine
        { domain->solInfo().fetiInfo.maxit = $2;
          domain->solInfo().fetiInfo.tol = $3;
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.numcgm = $4;
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
        | HFETI Integer Float Integer Float NewLine
        {
          domain->solInfo().fetiInfo.maxit = $2;
          domain->solInfo().fetiInfo.tol = $3;
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.numcgm = $4;
          domain->solInfo().fetiInfo.tolcgm = $5;
          domain->solInfo().fetiInfo.spaceDimension = 2;
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
        | HFETI Integer Float Integer Float Integer NewLine
        { domain->solInfo().fetiInfo.maxit = $2;
          domain->solInfo().fetiInfo.tol = $3;
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.numcgm = $4;
          domain->solInfo().fetiInfo.tolcgm = $5;
          domain->solInfo().fetiInfo.spaceDimension = $6;
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
        | HFETI NewLine
        { domain->solInfo().type =(2); 
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
        | NUMCGM Integer NewLine
        { domain->solInfo().fetiInfo.numcgm = $2; }
        | NUMCGM Integer Float NewLine
        { domain->solInfo().fetiInfo.numcgm = $2;
          domain->solInfo().fetiInfo.numcgm2 = $3; }
        | TOLCGM Float NewLine
        { domain->solInfo().fetiInfo.tolcgm = $2; }
        | SPACEDIMENSION Integer NewLine
        { domain->solInfo().fetiInfo.spaceDimension = $2; }
        | KRYLOVTYPE Integer NewLine
        { domain->solInfo().fetiInfo.krylovtype = $2; }
        | KRYLOVTYPE ISOLVERTYPE NewLine
        { domain->solInfo().fetiInfo.krylovtype =  $2; }
        | INTERFACELUMPED NewLine
        { domain->solInfo().fetiInfo.lumpedinterface = 1; }
        | SAVEMEMCOARSE NewLine
        { domain->solInfo().fetiInfo.saveMemCoarse = 1; }
        | TURKEL Integer NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureFlag = 0;
        }
        | TURKEL Integer Float NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureConst1 = $3;
          domain->curvatureFlag = 1;
        }
        | TURKEL Integer Float Float NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureConst1 = $3;
          domain->curvatureConst2 = $4;
          domain->curvatureFlag = 2;
        }
        | OUTERLOOP ITERTYPE NewLine
        { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) $2; }
        | OUTERLOOP ITERTYPE HERMITIAN NewLine
        { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) $2;
          domain->solInfo().fetiInfo.complex_hermitian = true; }
        | MPCTYPE Integer NewLine
        { domain->solInfo().fetiInfo.mpcflag = $2; }
        | MPCTYPE MPCTYPEID NewLine
        { domain->solInfo().fetiInfo.mpcflag = $2; }
        | MPCPRECNO Integer NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2; }
        | MPCPRECNO MPCPRECNOID NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2; }
        | MPCPRECNO Integer MPCBLK_OVERLAP Integer NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2;
          domain->solInfo().fetiInfo.mpcBlkOverlap = $4; }
        | MPCPRECNO MPCPRECNOID MPCBLK_OVERLAP Integer NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2;
          domain->solInfo().fetiInfo.mpcBlkOverlap = $4; }
        | MPCPRECNO Integer Integer NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2; 
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) $3; }
        | MPCPRECNO MPCPRECNOID MPCBLOCKID NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2;
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) $3; }
        | MPCPRECNO Integer Integer MPCBLK_OVERLAP Integer NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2;
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) $3; 
          domain->solInfo().fetiInfo.mpcBlkOverlap = $5; }
        | MPCPRECNO MPCPRECNOID MPCBLOCKID MPCBLK_OVERLAP Integer NewLine
        { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) $2;
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) $3; 
          domain->solInfo().fetiInfo.mpcBlkOverlap = $5; }
        | MRHS Integer NewLine
        { if($2 < 1) domain->solInfo().fetiInfo.useMRHS = false; }
        | PROPORTIONING Float NewLine
        { domain->solInfo().fetiInfo.gamma = $2; }
        | LINESEARCH Integer Float NewLine
        { domain->solInfo().fetiInfo.linesearch_maxit = $2;
          domain->solInfo().fetiInfo.linesearch_tau = $3; }
        | BMPC SWITCH NewLine
        { domain->solInfo().fetiInfo.bmpc = bool($2); }
        | DMPC SWITCH NewLine
        { domain->solInfo().fetiInfo.dmpc = bool($2); }
        | CMPC SWITCH NewLine
        { domain->solInfo().fetiInfo.cmpc = bool($2); }
        | CNORM SWITCH NewLine
        { domain->solInfo().fetiInfo.c_normalize = bool($2); }
        | MPCCHECK Integer NewLine
        { domain->solInfo().dbccheck = bool($2); }
	;
OldHelmInfo:
        FACOUSTICS NewLine FAcousticData
	;
FAcousticData:
        Float NewLine
        {
          /*domain->omega = $1;*/ geoSource->setOmega($1);
          StructProp sp; 
          sp.kappaHelm = $1;
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
	;
Constraints:
        CONSTRAINTS ConstraintOptionsData NewLine
        { if(!$2.lagrangeMult && $2.penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = $2.lagrangeMult;
          domain->solInfo().penalty = $2.penalty;
          domain->solInfo().constraint_hess = $2.constraint_hess; 
          domain->solInfo().constraint_hess_eps = $2.constraint_hess_eps; }
        | CONSTRAINTS NewLine ConstraintOptionsData NewLine
        { if(!$3.lagrangeMult && $3.penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = $3.lagrangeMult;
          domain->solInfo().penalty = $3.penalty;
          domain->solInfo().constraint_hess = $3.constraint_hess;
          domain->solInfo().constraint_hess_eps = $3.constraint_hess_eps; }
        ;
ConstraintOptionsData:
        DIRECT
        { // Direct elimination of slave dofs
          $$.lagrangeMult = false;
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
        }
        | DIRECT Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2; }
        | DIRECT Float Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2;
          domain->solInfo().coefFilterTol = $3; }
        | DIRECT Float Float Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2; 
          domain->solInfo().coefFilterTol = $3;
          domain->solInfo().rhsZeroTol = $4; }
        | DIRECT Float Float Float Float
        { $$.lagrangeMult = false; 
          $$.penalty = 0.0;
          $$.constraint_hess = 0;
          $$.constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = $2;
          domain->solInfo().coefFilterTol = $3; 
          domain->solInfo().rhsZeroTol = $4;
          domain->solInfo().inconsistentTol = $5; }
        | MULTIPLIERS
        { // Treatment of constraints through Lagrange multipliers method
          $$.lagrangeMult = true; 
          $$.penalty = 0.0;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | PENALTY Float
        { // Treatment of constraints through penalty method
          $$.lagrangeMult = false;
          $$.penalty = $2;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | MULTIPLIERS PENALTY Float
        { // Treatment of constraints through augmented Lagrangian method
          $$.lagrangeMult = true;
          $$.penalty = $3;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | AUGMENTED Float
        { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          $$.lagrangeMult = true;
          $$.penalty = $2;
          $$.constraint_hess = 1;
          $$.constraint_hess_eps = 0.0; }
        | ConstraintOptionsData HESSIAN Integer
        { $$.constraint_hess = $3;
          $$.constraint_hess_eps = 0; }
        | ConstraintOptionsData HESSIAN Integer Float
        { $$.constraint_hess = $3;
          $$.constraint_hess_eps = $4; }
        ;
HelmInfo:
        HELMHOLTZ NewLine
        { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
        | FETIH NewLine
        { domain->solInfo().type = (2); 
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
        | BGTL Integer NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureFlag = 0;
        }
        | BGTL Integer Float NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureConst1 = $3;
          domain->curvatureFlag = 1;
        }
        | BGTL Integer Float Float NewLine
        {
          domain->sommerfeldType = $2;
          domain->curvatureConst1 = $3;
          domain->curvatureConst2 = $4;
          domain->curvatureFlag = 2;
        }
        | POINTSOURCE Integer NewLine IncidenceList
        {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        }
        | PLANEWAVE Integer NewLine IncidenceList
        {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
        | INCIDENCE Integer NewLine IncidenceList
        {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
        | HelmInfo DampInfo
        | HelmInfo ReconsInfo 
        | HelmInfo PadePivotInfo
        ;
IncidenceList:
        IncidenceVector
        | IncidenceList IncidenceVector
        ;
IncidenceVector:
        Float Float Float NewLine
        { domain->setWaveDirections(0, $1, $2, $3); }
        ;
KirchhoffLocations:
        KIRLOC NewLine Float Float Float NewLine{
          domain->setKirchhoffLocations($3, $4, $5);
        } 
        | KirchhoffLocations Float Float Float NewLine {
          domain->setKirchhoffLocations($2, $3, $4);
        }
        ;
FFPDirList:
        FFPDirVector
        | FFPDirList FFPDirVector
        ;
FFPDirVector:
        Float Float Float NewLine
        { domain->setFFPDirections($1, $2, $3); }
        ;
HelmMFInfo:
        HELMMF NewLine FAcousticDataMF
        ;
FAcousticDataMF:
        Float NewLine
        {
          /*domain->omega = $1;*/ geoSource->setOmega($1);
          StructProp sp;
          sp.kappaHelm = $1;
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        }
        ;
HelmSOInfo:
        HELMSO NewLine FAcousticDataSO
        ;
FAcousticDataSO:
        Float NewLine
        {
          /*domain->omega = $1;*/ geoSource->setOmega($1);
          StructProp sp;
          sp.kappaHelm = $1;
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        }
        ;
DEMInfo:
        DEM NewLine
        {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        }
        ;
NLInfo:
        NL NewLine
        { 
          if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
            domain->solInfo().probType = SolverInfo::NonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::Dynamic)
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
            domain->solInfo().order = 1;
            domain->solInfo().probType = SolverInfo::NonLinDynam;
            domain->solInfo().probType = SolverInfo::TempDynamic;
          }
          domain->solInfo().fetiInfo.type = FetiInfo::nonlinear; // XXXX
        }
        | NLInfo ARCLENGTH NewLine
        { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength;
        }
        | NLInfo NLMAT NewLine
        { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam;
        }
        | NLInfo MAXITR Integer NewLine
        { domain->solInfo().getNLInfo().maxiter = $3; }
        | NLInfo NLTOL Float NewLine
        { domain->solInfo().getNLInfo().tolRes = $3; }
        | NLInfo NLTOL Float Float NewLine
        { domain->solInfo().getNLInfo().tolRes = $3;
          domain->solInfo().getNLInfo().tolInc = $4; }
        | NLInfo NLTOL Float Float Float Float NewLine
        { domain->solInfo().getNLInfo().tolRes = $3;
          domain->solInfo().getNLInfo().tolInc = $4;
          domain->solInfo().getNLInfo().absTolRes = $5;
          domain->solInfo().getNLInfo().absTolInc = $6; }
        | NLInfo DLAMBDA Float Float NewLine
        { domain->solInfo().getNLInfo().dlambda = $3;
          domain->solInfo().getNLInfo().maxLambda = $4; }
        | NLInfo DLAMBDA Float Float Integer Integer NewLine
        { domain->solInfo().getNLInfo().dlambda = $3; 
          domain->solInfo().getNLInfo().maxLambda = $4;
          domain->solInfo().getNLInfo().extMin = $5;
          domain->solInfo().getNLInfo().extMax = $6; }
        | NLInfo FITALG Integer NewLine
        { domain->solInfo().getNLInfo().fitAlgShell = $3;
          domain->solInfo().getNLInfo().fitAlgBeam  = $3; }
        | NLInfo FITALG Integer Integer NewLine
        { domain->solInfo().getNLInfo().fitAlgShell = $3;
          domain->solInfo().getNLInfo().fitAlgBeam  = $4; }
        | NLInfo UNSYMMETRIC NewLine
        { domain->solInfo().getNLInfo().unsymmetric = true; }
        | NLInfo LFACTOR Float NewLine
        { domain->solInfo().getNLInfo().lfactor = $3; }
/* conflict
        | NLInfo LINESEARCH NewLine
        { domain->solInfo().getNLInfo().linesearch = true; }
*/
        | NLInfo FAILSAFE NewLine
        { domain->solInfo().getNLInfo().failsafe = true; }
        | NLInfo FAILSAFE Float NewLine
        { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = $3; }
        | NLInfo PENALTY Integer Float Float NewLine
        { domain->solInfo().num_penalty_its = $3; 
          domain->solInfo().penalty_tol = $4;
          domain->solInfo().penalty_beta = $5; }
        | NLInfo NewtonInfo
        ;
NewtonInfo:
        REBUILD Integer NewLine
        { 
          domain->solInfo().setNewton($2); 
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear; 
        }
	| REBUILD Integer Integer NewLine
	{ 
          domain->solInfo().setNewton($2); 
          int rebuildK    = $2; 
          int rebuildPrec = $3;
          if(rebuildK > 1) rebuildPrec = rebuildK;
          domain->solInfo().fetiInfo.nTang = rebuildK;
          domain->solInfo().fetiInfo.nPrec = rebuildPrec;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        }
	| REBUILD NewLine TANGENT Integer NewLine
	{
          domain->solInfo().setNewton($4);
          domain->solInfo().fetiInfo.nPrec = $4;
          domain->solInfo().fetiInfo.nTang = $4;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        }
	| REBUILD NewLine TANGENT Integer NewLine PRECONDITIONER Integer NewLine
	{
	  domain->solInfo().setNewton($4); 
          int rebuildK    = $4; 
          int rebuildPrec = $7;
          if(rebuildK > 1) rebuildPrec = rebuildK;
          domain->solInfo().fetiInfo.nTang = rebuildK;
          domain->solInfo().fetiInfo.nPrec = rebuildPrec;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
	}
	;
OrthoInfo:
	REORTHO NewLine
        { domain->solInfo().setReOrtho(); }
	;
Control:
	CONTROL NewLine
        | CONTROL NewLine FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($3,$7,$9); domain->solInfo().soltyp = $5; }
        | CONTROL NewLine FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($3,$7,$9,$11); domain->solInfo().soltyp = $5; }
/*
        | Control FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($2,$6,$8); }
        | Control FNAME NewLine Integer NewLine FNAME NewLine FNAME NewLine FNAME NewLine
         { geoSource->setControl($2,$6,$8,$10); }
*/
	;
Optimization:
	OPTIMIZATION FNAME NewLine
        { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = $2;
#endif
        }
        | Optimization NewLine FNAME NewLine
        { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = $3;
#endif
 }
	;
NodalContact: 
        NODALCONTACT NewLine
        | NODALCONTACT MODE Integer NewLine
        { domain->solInfo().contact_mode = $3; }

        // glbNode01 glbNode02 nx ny nz
        | NodalContact Integer Integer Float Float Float NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6); }
        // glbNode01 glbNode02 nx ny nz GAP normalGap
        | NodalContact Integer Integer Float Float Float GAP Float NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, $8);}
        // glbNode01 glbNode02 nx ny nz MODE mode
        | NodalContact Integer Integer Float Float Float MODE Integer NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, 0.0, $8);}
        | NodalContact Integer Integer Float Float Float ConstraintOptionsData NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, 0.0, -1, $7.lagrangeMult, $7.penalty);}
        // glbNode01 glbNode02 nx ny nz GAP normalGap MODE mode
        | NodalContact Integer Integer Float Float Float GAP Float MODE Integer NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, $8, $10);}
        | NodalContact Integer Integer Float Float Float GAP Float MODE Integer ConstraintOptionsData NewLine
        { domain->addNodalCTC($2-1, $3-1, $4, $5, $6, $8, $10, $11.lagrangeMult, $11.penalty);}
	;
MatSpec:
	MATSPEC NewLine
	| MatSpec Integer BILINEARPLASTIC Float Float Float Float Float NewLine
	 { 
           geoSource->addMaterial($2-1, 
             new BilinPlasKinHardMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer BILINEARPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new BilinPlasKinHardMat($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer FINITESTRAINPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new FiniteStrainPlasKinHardMat($4, $5, $6, $7, $8, $9) );
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new LogStrainPlasKinHardMat($4, $5, $6, $7, $8) );
         }
        | MatSpec Integer LOGSTRAINPLASTIC Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new LogStrainPlasKinHardMat($4, $5, $6, $7, $8, $9) );
         }
	| MatSpec Integer LINEARELASTIC Float Float Float NewLine
	 { 
           geoSource->addMaterial($2-1, 
             new ElaLinIsoMat($4, $5, $6));
	 }
        | MatSpec Integer STVENANTKIRCHHOFF Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat($4, $5, $6));
         }
        | MatSpec Integer HENCKY Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new HenckyMat($4, $5, $6));
         }
        | MatSpec Integer LINPLSTRESS Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ElaLinIsoMat2D($4, $5, $6, $7));
         }
        | MatSpec Integer SVKPLSTRESS Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new StVenantKirchhoffMat2D($4, $5, $6, $7));
         }
        | MatSpec Integer ISOTROPICLINEARELASTIC Float Float Float NewLine
          {
            double params[3] = { $4, $5, $6 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          }
        | MatSpec Integer NEOHOOKEAN Float Float Float NewLine
          {
            double params[4] = { $4, $5, $6, -1 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<NeoHookean>(params));
          }
        | MatSpec Integer NEOHOOKEAN Float Float Float Float NewLine
          {
            double params[4] = { $4, $5, $6, $7 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<NeoHookean>(params));
          }
        | MatSpec Integer MOONEYRIVLIN Float Float Float Float NewLine
          {
            double params[5] = { $4, $5, $6, $7, -1 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<MooneyRivlin>(params));
          }
        | MatSpec Integer MOONEYRIVLIN Float Float Float Float Float NewLine
          {
            double params[5] = { $4, $5, $6, $7, $8 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<MooneyRivlin>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTIC Float Float Float Float Float Float NewLine
          {
            double params[6] = { $4, $5, $6, $7, $8, $9 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float NewLine
          {
            double params[8] = { $4, $5, $6, $7, $8, $9, 1.0e-6, -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float Float NewLine
          {
            double params[8] = { $4, $5, $6, $7, $8, $9, $10, -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
        | MatSpec Integer ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS Float Float Float Float Float Float Float Float NewLine
          {
            double params[8] = { $4, $5, $6, $7, $8, $9, $10, $11 };
            geoSource->addMaterial($2-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23));
         }
        | MatSpec Integer OPTCTV Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         }
        | MatSpec Integer OPTCTV Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, $9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         }
        | MatSpec Integer OPTCTV Float Float Float Float Float Float Float NewLine
         {
           geoSource->addMaterial($2-1,
             new ExpMat($3, $4, $5, $6, $7, $8, $9, $10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         }
	| MatSpec READ FNAME FNAME NewLine
	 {
	   geoSource->loadMaterial($3, $4);
	 }
	| MatSpec Integer FNAME FloatList NewLine
	 {
	   geoSource->addMaterial($2-1, $3, $4);
	 }
	;
MatUsage:
	MATUSAGE NewLine
	| MatUsage Integer Integer NewLine
	  { geoSource->setMatUsage($2-1, $3-1); }
	| MatUsage Integer Integer Integer NewLine
	  {
            for(int i = $2-1; i < $3; ++i)
	      geoSource->setMatUsage(i, $4-1);
	  }
	;
FloatList:
	{ $$.nval = 0; }
	| FloatList Float
	{ 
          if($1.nval == 32) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          $$ = $1;
          $$.v[$$.nval++] = $2;
 	}
	;
Renumbering:
	RENUM NewLine RENUMBERID NewLine
	{ domain->solInfo().setRenum($3);
          domain->solInfo().setSparseRenum($3); 
          domain->solInfo().setSpoolesRenum($3); }
        | RENUM NewLine RENUMBERID NewLine RENUMBERID NewLine
        { domain->solInfo().setRenum($3);
          domain->solInfo().setSparseRenum($5); }
        | RENUM NewLine RENUMBERID NewLine RENUMBERID NewLine RENUMBERID NewLine
        { domain->solInfo().setRenum($3);
          domain->solInfo().setSparseRenum($5); 
          domain->solInfo().setSpoolesRenum($7); }
	;

SvdToken:
    SVDTOKEN NewLine
  { domain->solInfo().activatePodRom = true; 
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().svdPodRom = true;}
  | SvdToken SvdOption NewLine
  ;

SvdOption:
    SNAPFI FNAME
  { domain->solInfo().snapfiPodRom = $2; }
  | SNAPFI FNAME Integer 
  { domain->solInfo().snapfiPodRom = $2;
    if ($3 == 1) domain->solInfo().statevectPodRom = true;
    if ($3 == 2) domain->solInfo().residvectPodRom = true;
    if ($3 == 3) domain->solInfo().jacobvectPodRom = true;
    if ($3 == 4) domain->solInfo().forcevectPodRom = true;
    if ($3 == 5) domain->solInfo().accelvectPodRom = true;}
  | PODSIZEMAX Integer
  { domain->solInfo().maxSizePodRom = $2; }
  ;

Sampling:
    SAMPLING NewLine 
  { domain->solInfo().activatePodRom = true; 
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().samplingPodRom = true; }
  | Sampling SamplingOption NewLine
  ;

SamplingOption:
    PODROB FNAME
  { domain->solInfo().readInROBorModes = $2; }
  | TRNVCT FNAME
  { domain->solInfo().statePodRomFile = $2; }
  | TRNVCT FNAME FNAME
  { domain->solInfo().statePodRomFile = $2;
    domain->solInfo().velocPodRomFile = $3; }
  | TOLER Float
  { domain->solInfo().tolPodRom = $2; }
  | SKIP Integer
  { domain->solInfo().skipPodRom = $2; }
  | OFFSET Integer
  { domain->solInfo().skipOffSet = $2; }
  | PODSIZEMAX Integer
  { domain->solInfo().maxSizePodRom = $2; }
  ;

ConversionToken:
    CONVERSIONTOKEN NewLine
  { domain->solInfo().activatePodRom = true;
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().ROMPostProcess = true; }
  | ConversionToken ConversionOption NewLine
  ;

ConversionOption:
    CONVFI FNAME
  { domain->solInfo().RODConversionFiles.push_back($2); 
    domain->solInfo().numRODFile += 1; }
  ;

Integer:
	IntConstant
	{ $$ = $1; }
	;

Float:
	IntConstant
	{ $$ = $1; }
	| DblConstant 
	{ $$ = $1; }
        | INFINTY
        { $$ = std::numeric_limits<double>::infinity(); }
        | EPSILON
        { $$ = std::numeric_limits<double>::epsilon(); }
	;
%%
