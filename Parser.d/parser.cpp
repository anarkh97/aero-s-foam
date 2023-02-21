/* A Bison parser, made by GNU Bison 3.0.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "3.0.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
#line 1 "p.y" /* yacc.c:339  */

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <map>
#include <cstdlib>
#include <Parser.d/AuxDefs.h>
#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.h>
#include <Element.d/NonLinearity.d/CrushableFoam.h>
#include <Element.d/NonLinearity.d/ElaLinIsoMat.h>
#include <Element.d/NonLinearity.d/2DMat.h>
#include <Element.d/NonLinearity.d/ExpMat.h>
#include <Element.d/NonLinearity.d/MaterialWrapper.h>
#include <Element.d/NonLinearity.d/PronyViscoElastic.h>
#include <Element.d/NonLinearity.d/OgdenMat.h>
#include <Element.d/NonLinearity.d/SimoElasticMat.h>
#include <Element.d/NonLinearity.d/SimoPlasticMat.h>
#include <Element.d/NonLinearity.d/NeoHookeanMat.h>
#include <Element.d/NonLinearity.d/MooneyRivlinMat.h>
#include <Element.d/NonLinearity.d/BrittleFractureTB.h>
#include <Element.d/NonLinearity.d/PlaneStressMat.h>
#include <Driver.d/Domain.h>
#include <Sfem.d/Sfem.h>
#include <Utils.d/DistHelper.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#ifdef STRUCTOPT
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#endif
#ifdef USE_EIGEN3
#include <Eigen/Core>
#endif

 int numColumns = 3;
 double amplitude = 1.0;
 int PitaTS = 1;         //CD: Pita
 extern std::string clusterData_;
 extern std::string subdomains_;
 extern std::string decomposition_;
 extern std::string connectivity_;
 extern bool randomShuffle;
 extern bool allowMechanisms;
 extern bool useScotch;

#line 112 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:339  */

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "parser.hpp".  */
#ifndef YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMTESTING_PARSER_D_PARSER_HPP_INCLUDED
# define YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMTESTING_PARSER_D_PARSER_HPP_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token type.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    ACTUATORS = 258,
    ADJOINTBASIS = 259,
    AERO = 260,
    AEROH = 261,
    AEROTYPE = 262,
    AGRESSIVETOLERANCES = 263,
    ALPROC = 264,
    AMAT = 265,
    ANALYSIS = 266,
    ARCLENGTH = 267,
    ATTRIBUTES = 268,
    ANGULAROUTTYPE = 269,
    ARUBBERMAT = 270,
    AUGMENT = 271,
    AUGMENTTYPE = 272,
    AUTOTOL = 273,
    AUXILIARY = 274,
    AVERAGED = 275,
    ATDARB = 276,
    ACOU = 277,
    ATDDNB = 278,
    ATDROB = 279,
    ARPACK = 280,
    ATDDIR = 281,
    ATDNEU = 282,
    ALLOWMECHANISMS = 283,
    AUXCOARSESOLVER = 284,
    ACMECNTL = 285,
    ADDEDMASS = 286,
    AEROEMBED = 287,
    ANDESCLR = 288,
    ANDESCQR = 289,
    ANDESBETAB = 290,
    ANDESALPHA = 291,
    ANDESBETAM = 292,
    AUGMENTED = 293,
    BLOCKDIAG = 294,
    BOFFSET = 295,
    BUCKLE = 296,
    BGTL = 297,
    BMPC = 298,
    BINARYINPUT = 299,
    BINARYOUTPUT = 300,
    BLOCKSIZE = 301,
    CHECKTOKEN = 302,
    COARSESOLVER = 303,
    COEF = 304,
    CFRAMES = 305,
    COLLOCATEDTYPE = 306,
    CONVECTION = 307,
    COMPOSITE = 308,
    CONDITION = 309,
    CONTACT = 310,
    CONTROL = 311,
    CORNER = 312,
    CORNERTYPE = 313,
    CURVE = 314,
    CCTTOL = 315,
    CCTSOLVER = 316,
    CRHS = 317,
    COUPLEDSCALE = 318,
    CONTACTSURFACES = 319,
    CMPC = 320,
    CNORM = 321,
    COMPLEXOUTTYPE = 322,
    CONSTRMAT = 323,
    CASES = 324,
    CONSTRAINEDSURFACES = 325,
    CSFRAMES = 326,
    CSTYPE = 327,
    CONSTANT = 328,
    CONWEP = 329,
    DAMPING = 330,
    DblConstant = 331,
    DEFAULTPENALTY = 332,
    DELETEELEMENTS = 333,
    DEM = 334,
    DIMASS = 335,
    DISP = 336,
    DIRECT = 337,
    DLAMBDA = 338,
    DP = 339,
    DYNAM = 340,
    DETER = 341,
    DECOMPOSE = 342,
    DECOMPFILE = 343,
    DMPC = 344,
    DEBUGCNTL = 345,
    DEBUGICNTL = 346,
    DOCLUSTERING = 347,
    DOROWCLUSTERING = 348,
    ANGLE = 349,
    DUALBASIS = 350,
    DUALRB = 351,
    KMEANS = 352,
    CRANDOM = 353,
    CONSTRAINTS = 354,
    MULTIPLIERS = 355,
    PENALTY = 356,
    ELLUMP = 357,
    EIGEN = 358,
    EFRAMES = 359,
    ELSCATTERER = 360,
    END = 361,
    ELHSOMMERFELD = 362,
    ENGINEERING = 363,
    ETEMP = 364,
    EXPLICIT = 365,
    EXTFOL = 366,
    EPSILON = 367,
    ELEMENTARYFUNCTIONTYPE = 368,
    FABMAT = 369,
    FABRICMAP = 370,
    FABRICMAT = 371,
    FACE = 372,
    FACOUSTICS = 373,
    FETI = 374,
    FETI2TYPE = 375,
    FETIPREC = 376,
    FFP = 377,
    FFPDIR = 378,
    FITALG = 379,
    FNAME = 380,
    FLAGMN = 381,
    FLUX = 382,
    FORCE = 383,
    FRONTAL = 384,
    FETIH = 385,
    FIELDWEIGHTLIST = 386,
    FILTEREIG = 387,
    FLUID = 388,
    FREEPLAY = 389,
    FREQSWEEP = 390,
    FREQSWEEP1 = 391,
    FREQSWEEP2 = 392,
    FREQSWEEPA = 393,
    FREQSWEEPAW = 394,
    FSGL = 395,
    FSINTERFACE = 396,
    FSISCALING = 397,
    FSIELEMENT = 398,
    NOLOCALFSISPLITING = 399,
    FSICORNER = 400,
    FFIDEBUG = 401,
    FAILSAFE = 402,
    FRAMETYPE = 403,
    GEPS = 404,
    GLOBALSEARCHCULL = 405,
    GLOBALTOL = 406,
    GRAVITY = 407,
    GRBM = 408,
    GTGSOLVER = 409,
    GLOBALCRBMTOL = 410,
    GROUP = 411,
    GROUPTYPE = 412,
    GOLDFARBTOL = 413,
    HDIRICHLET = 414,
    HEAT = 415,
    HFETI = 416,
    HNEUMAN = 417,
    HSOMMERFELD = 418,
    HFTT = 419,
    HELMHOLTZ = 420,
    HNBO = 421,
    HELMMF = 422,
    HELMSO = 423,
    HSCBO = 424,
    HWIBO = 425,
    HZEM = 426,
    HZEMFILTER = 427,
    HLMPC = 428,
    HERMITIAN = 429,
    HESSIAN = 430,
    IACC = 431,
    IDENTITY = 432,
    IDIS = 433,
    IDIS6 = 434,
    ILUDROPTOL = 435,
    IntConstant = 436,
    INTERFACELUMPED = 437,
    ITEMP = 438,
    ITERTYPE = 439,
    IVEL = 440,
    IMESH = 441,
    INCIDENCE = 442,
    IHDIRICHLET = 443,
    IHDSWEEP = 444,
    IHNEUMANN = 445,
    ISOLVERTYPE = 446,
    INPC = 447,
    INFINTY = 448,
    JACOBI = 449,
    KEYLETTER = 450,
    KRYLOVTYPE = 451,
    KIRLOC = 452,
    LAYC = 453,
    LAYN = 454,
    LAYD = 455,
    LAYO = 456,
    LAYMAT = 457,
    LFACTOR = 458,
    LISRBM = 459,
    LMPC = 460,
    LOAD = 461,
    LOADCASE = 462,
    LOBPCG = 463,
    LOCALREDUCEDORDERBASES = 464,
    LOCALSOLVER = 465,
    LINESEARCH = 466,
    LUMPED = 467,
    KSPARAM = 468,
    KSMAX = 469,
    MASS = 470,
    MASSAUGMENTATION = 471,
    MATERIALS = 472,
    MATLAB = 473,
    MAXITR = 474,
    MAXELEM = 475,
    MAXORTHO = 476,
    MAXVEC = 477,
    MODAL = 478,
    MPCPRECNO = 479,
    MPCPRECNOID = 480,
    MPCTYPE = 481,
    MPCTYPEID = 482,
    MPCSCALING = 483,
    MPCELEMENT = 484,
    MPCBLOCKID = 485,
    MPCBLK_OVERLAP = 486,
    MFTT = 487,
    MRHS = 488,
    MPCCHECK = 489,
    MUMPSICNTL = 490,
    MUMPSCNTL = 491,
    MUMPSMINEQ = 492,
    MUMPSSTRIDE = 493,
    MECH = 494,
    MODDAMP = 495,
    MODEFILTER = 496,
    MOMENTTYPE = 497,
    MPROJECT = 498,
    MAXIMUM = 499,
    NOWARPEDVOLUME = 500,
    NDTYPE = 501,
    NEIGPA = 502,
    NEWMARK = 503,
    NewLine = 504,
    NEWTON = 505,
    NL = 506,
    NLMAT = 507,
    NLPREC = 508,
    NLMEMPTYPE = 509,
    NOCOARSE = 510,
    NODETOKEN = 511,
    NOMULTIPLEINTERACTIONS = 512,
    NONINPC = 513,
    NSBSPV = 514,
    NLTOL = 515,
    NUMCGM = 516,
    NONORMALSMOOTHING = 517,
    NOSECONDARY = 518,
    NOGHOSTING = 519,
    NFRAMES = 520,
    NORMALSMOOTHINGDISTANCE = 521,
    SENSITIVITY = 522,
    SENSITIVITYMETHOD = 523,
    SKIPPHYSICALFACES = 524,
    OUTPUT = 525,
    OUTPUT6 = 526,
    OUTPUTFRAME = 527,
    OLDDYNAMICSEARCH = 528,
    QSTATIC = 529,
    QLOAD = 530,
    QUASISTATIC = 531,
    PITA = 532,
    PITADISP6 = 533,
    PITAVEL6 = 534,
    NOFORCE = 535,
    MDPITA = 536,
    GLOBALBASES = 537,
    LOCALBASES = 538,
    TIMEREVERSIBLE = 539,
    REMOTECOARSE = 540,
    ORTHOPROJTOL = 541,
    READINITSEED = 542,
    JUMPCVG = 543,
    JUMPOUTPUT = 544,
    PRECNO = 545,
    PRECONDITIONER = 546,
    PRELOAD = 547,
    PRESSURE = 548,
    PRINTMATLAB = 549,
    printmatlab = 550,
    PRINTNUMBER = 551,
    PROJ = 552,
    PIVOT = 553,
    PRECTYPE = 554,
    PRECTYPEID = 555,
    PICKANYCORNER = 556,
    PADEPIVOT = 557,
    PROPORTIONING = 558,
    PLOAD = 559,
    PADEPOLES = 560,
    POINTSOURCE = 561,
    PLANEWAVE = 562,
    PTOL = 563,
    PLANTOL = 564,
    PMAXIT = 565,
    PIECEWISE = 566,
    PARAMETERS = 567,
    RADIATION = 568,
    RBMFILTER = 569,
    RBMSET = 570,
    READMODE = 571,
    READSENSITIVITY = 572,
    REBUILD = 573,
    RESOLUTIONMETHOD = 574,
    REVERSEORDER = 575,
    REDFOL = 576,
    RENUM = 577,
    RENUMBERID = 578,
    REORTHO = 579,
    RESTART = 580,
    RECONS = 581,
    RECONSALG = 582,
    REBUILDCCT = 583,
    RANDOM = 584,
    RPROP = 585,
    RNORM = 586,
    REVERSENORMALS = 587,
    ROBTYPE = 588,
    ROTVECOUTTYPE = 589,
    RESCALING = 590,
    RUBDAFT = 591,
    SCALING = 592,
    SCALINGTYPE = 593,
    SDETAFT = 594,
    SENSORS = 595,
    SHELLFABRICMAP = 596,
    SHELLFABRICMAT = 597,
    SHELLSIMPLELOFTING = 598,
    SOLVERCNTL = 599,
    SOLVERHANDLE = 600,
    SOLVERTYPE = 601,
    SHIFT = 602,
    SHARPNONSHARPANGLE = 603,
    SPOOLESTAU = 604,
    SPOOLESSEED = 605,
    SPOOLESMAXSIZE = 606,
    SPOOLESMAXDOMAINSIZE = 607,
    SPOOLESMAXZEROS = 608,
    SPOOLESMSGLVL = 609,
    SPOOLESSCALE = 610,
    SPOOLESPIVOT = 611,
    SPOOLESRENUM = 612,
    SPARSEMAXSUP = 613,
    SPARSEDEFBLK = 614,
    SPARSERENUM = 615,
    STATS = 616,
    STRESSID = 617,
    SUBCYCLE = 618,
    SUBSPACE = 619,
    SURFACE = 620,
    STR_THERM_OPTION = 621,
    SAVEMEMCOARSE = 622,
    SPACEDIMENSION = 623,
    SCATTERER = 624,
    STAGTOL = 625,
    SCALED = 626,
    SWITCH = 627,
    STABLE = 628,
    SUBTYPE = 629,
    STEP = 630,
    SOWER = 631,
    SHELLTHICKNESS = 632,
    SURF = 633,
    SPRINGMAT = 634,
    SENSITIVITYID = 635,
    TABLE = 636,
    TANGENT = 637,
    TDENFORCE = 638,
    TEMP = 639,
    TIME = 640,
    TOLEIG = 641,
    TOLFETI = 642,
    TOLJAC = 643,
    TOLPCG = 644,
    TOLSEN = 645,
    TOPFILE = 646,
    TOPOLOGY = 647,
    TRBM = 648,
    TRBMlc = 649,
    THERMOE = 650,
    THERMOH = 651,
    RATIOTOLSEN = 652,
    TETT = 653,
    TOLCGM = 654,
    TURKEL = 655,
    TIEDSURFACES = 656,
    THETA = 657,
    PROJSOL = 658,
    CENTER = 659,
    POSELEM = 660,
    HOTSTART = 661,
    HRC = 662,
    THIRDNODE = 663,
    THERMMAT = 664,
    TDENFORC = 665,
    TESTULRICH = 666,
    THRU = 667,
    TRIVIAL = 668,
    NUMROMCPUS = 669,
    THICKNESSGROUPLIST = 670,
    USE = 671,
    USERDEFINEDISP = 672,
    USERDEFINEFORCE = 673,
    UPROJ = 674,
    UNSYMMETRIC = 675,
    USING = 676,
    USESCOTCH = 677,
    VERBOSE = 678,
    VERSION = 679,
    WETCORNERS = 680,
    YMTT = 681,
    SS1DT = 682,
    SS2DT = 683,
    YMST = 684,
    YSST = 685,
    YSSRT = 686,
    ZERO = 687,
    BINARY = 688,
    GEOMETRY = 689,
    DECOMPOSITION = 690,
    GLOBAL = 691,
    MATCHER = 692,
    CPUMAP = 693,
    NODALCONTACT = 694,
    MODE = 695,
    FRIC = 696,
    GAP = 697,
    OUTERLOOP = 698,
    EDGEWS = 699,
    WAVETYPE = 700,
    ORTHOTOL = 701,
    IMPE = 702,
    FREQ = 703,
    DPH = 704,
    WAVEMETHOD = 705,
    MATSPEC = 706,
    MATUSAGE = 707,
    BILINEARPLASTIC = 708,
    FINITESTRAINPLASTIC = 709,
    CRUSHABLEFOAM = 710,
    LINEARELASTIC = 711,
    STVENANTKIRCHHOFF = 712,
    TULERBUTCHER = 713,
    LINPLSTRESS = 714,
    READ = 715,
    OPTCTV = 716,
    ISOTROPICLINEARELASTIC = 717,
    VISCOLINEARELASTIC = 718,
    VISCOSTVENANTKIRCHHOFF = 719,
    NEOHOOKEAN = 720,
    VISCONEOHOOKEAN = 721,
    ISOTROPICLINEARELASTICJ2PLASTIC = 722,
    ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 723,
    HYPERELASTIC = 724,
    MOONEYRIVLIN = 725,
    VISCOMOONEYRIVLIN = 726,
    HENCKY = 727,
    OGDEN = 728,
    SIMOELASTIC = 729,
    SIMOPLASTIC = 730,
    LOGSTRAINPLASTIC = 731,
    SVKPLSTRESS = 732,
    VISCOLINPLSTRESS = 733,
    VISCOSVKPLSTRESS = 734,
    VISCOFABRICMAP = 735,
    SHELLVISCOFABRICMAP = 736,
    VISCOFABRICMAT = 737,
    SHELLVISCOFABRICMAT = 738,
    PARTITIONGAP = 739,
    PLANESTRESSLINEAR = 740,
    PLANESTRESSSTVENANTKIRCHHOFF = 741,
    PLANESTRESSNEOHOOKEAN = 742,
    PLANESTRESSMOONEYRIVLIN = 743,
    PLANESTRESSBILINEARPLASTIC = 744,
    PLANESTRESSFINITESTRAINPLASTIC = 745,
    PLANESTRESSVISCOLINEARELASTIC = 746,
    PLANESTRESSVISCOSTVENANTKIRCHHOFF = 747,
    PLANESTRESSVISCONEOHOOKEAN = 748,
    PLANESTRESSVISCOMOONEYRIVLIN = 749,
    SURFACETOPOLOGY = 750,
    MORTARTIED = 751,
    MORTARSCALING = 752,
    MORTARINTEGRATIONRULE = 753,
    SEARCHTOL = 754,
    STDMORTAR = 755,
    DUALMORTAR = 756,
    WETINTERFACE = 757,
    NSUBS = 758,
    EXITAFTERDEC = 759,
    SKIP = 760,
    ROBCSOLVE = 761,
    RANDOMSAMPLE = 762,
    OUTPUTMEMORY = 763,
    OUTPUTWEIGHT = 764,
    SOLVER = 765,
    SPNNLSSOLVERTYPE = 766,
    MAXSIZE = 767,
    CLUSTERSOLVER = 768,
    CLUSTERSOLVERTYPE = 769,
    WEIGHTLIST = 770,
    GMRESRESIDUAL = 771,
    SLOSH = 772,
    SLGRAV = 773,
    SLZEM = 774,
    SLZEMFILTER = 775,
    PDIR = 776,
    HEFSB = 777,
    HEFRS = 778,
    HEINTERFACE = 779,
    SNAPFI = 780,
    VELSNAPFI = 781,
    ACCSNAPFI = 782,
    DSVSNAPFI = 783,
    MUVSNAPFI = 784,
    PODROB = 785,
    ROMENERGY = 786,
    TRNVCT = 787,
    OFFSET = 788,
    ORTHOG = 789,
    SVDTOKEN = 790,
    CONVERSIONTOKEN = 791,
    CONVFI = 792,
    ROMRES = 793,
    SAMPLING = 794,
    SNAPSHOTPROJECT = 795,
    PODSIZEMAX = 796,
    REFSUBTRACT = 797,
    TOLER = 798,
    NORMALIZETOKEN = 799,
    FNUMBER = 800,
    SNAPWEIGHT = 801,
    ROBFI = 802,
    STAVCT = 803,
    VELVCT = 804,
    ACCVCT = 805,
    CONWEPCFG = 806,
    SCALEPOSCOORDS = 807,
    NODEPOSCOORDS = 808,
    MESHSCALEFACTOR = 809,
    PSEUDOGNAT = 810,
    PSEUDOGNATELEM = 811,
    USENMF = 812,
    USENMFC = 813,
    USEGREEDY = 814,
    USEPQN = 815,
    FILTERROWS = 816,
    VECTORNORM = 817,
    REBUILDFORCE = 818,
    REBUILDCONSTRAINT = 819,
    SAMPNODESLOT = 820,
    REDUCEDSTIFFNESS = 821,
    UDEIMBASIS = 822,
    FORCEROB = 823,
    CONSTRAINTROB = 824,
    DEIMINDICES = 825,
    UDEIMINDICES = 826,
    SVDFORCESNAP = 827,
    SVDCONSTRAINTSNAP = 828,
    USEMASSNORMALIZEDBASIS = 829,
    USECONSTANTMASS = 830,
    ONLINEMASSNORMALIZEBASIS = 831,
    STACKED = 832,
    NUMTHICKNESSGROUP = 833,
    STRESSNODELIST = 834,
    DISPNODELIST = 835,
    RELAXATIONSEN = 836,
    QRFACTORIZATION = 837,
    QMATRIX = 838,
    RMATRIX = 839,
    XMATRIX = 840,
    EIGENVALUE = 841,
    NPMAX = 842,
    BSSPLH = 843,
    PGSPLH = 844,
    LIB = 845
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 48 "p.y" /* yacc.c:355  */

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
 MFTTData *sdetaft;
#ifdef USE_EIGEN3
 GenMFTTData<Eigen::Vector4d> *rubdaft;
#endif
 SS2DTData *ss2dt;
 ComplexBCList *cxbclist;
 ComplexBCond cxbcval;
 FrameData frame;
 NodalFrameData nframe;
 TrivialPair<MFTTData*,int> mftval;
 TrivialPair<MFTTData*,int> hftval;
 LayerData ldata;
 LayInfo *linfo;
 CoefData coefdata;
 DoubleList dlist;
 StringList slist;
 SurfaceEntity* SurfObj;
 MortarHandler* MortarCondObj;
 LMPCTerm* mpcterm;
 GeoSource::Rprop rprop;
 OutputInfo oinfo;
 ConstraintOptions copt;
 BlastLoading::BlastData blastData;
 SolverCntl* scntl;
 FreeplayProps freeplayProps;
 ModalParams::Type mpt;

#line 783 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMTESTING_PARSER_D_PARSER_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 800 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:358  */

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif

#ifndef YY_ATTRIBUTE
# if (defined __GNUC__                                               \
      && (2 < __GNUC__ || (__GNUC__ == 2 && 96 <= __GNUC_MINOR__)))  \
     || defined __SUNPRO_C && 0x5110 <= __SUNPRO_C
#  define YY_ATTRIBUTE(Spec) __attribute__(Spec)
# else
#  define YY_ATTRIBUTE(Spec) /* empty */
# endif
#endif

#ifndef YY_ATTRIBUTE_PURE
# define YY_ATTRIBUTE_PURE   YY_ATTRIBUTE ((__pure__))
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# define YY_ATTRIBUTE_UNUSED YY_ATTRIBUTE ((__unused__))
#endif

#if !defined _Noreturn \
     && (!defined __STDC_VERSION__ || __STDC_VERSION__ < 201112)
# if defined _MSC_VER && 1200 <= _MSC_VER
#  define _Noreturn __declspec (noreturn)
# else
#  define _Noreturn YY_ATTRIBUTE ((__noreturn__))
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN \
    _Pragma ("GCC diagnostic push") \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")\
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif


#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYSIZE_T yynewbytes;                                            \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                          \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  542
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   8575

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  591
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  260
/* YYNRULES -- Number of rules.  */
#define YYNRULES  1430
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  3799

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   845

#define YYTRANSLATE(YYX)                                                \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, without out-of-bounds checking.  */
static const yytype_uint16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   196,   196,   203,   204,   207,   208,   210,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   234,   235,   237,
     238,   239,   240,   241,   242,   244,   245,   246,   247,   248,
     249,   252,   253,   254,   255,   256,   257,   258,   259,   260,
     261,   262,   263,   264,   265,   266,   268,   270,   271,   272,
     273,   274,   275,   276,   277,   278,   279,   280,   281,   282,
     283,   284,   285,   286,   287,   288,   289,   290,   291,   292,
     293,   294,   295,   296,   297,   298,   299,   300,   301,   302,
     304,   306,   308,   310,   312,   314,   315,   316,   317,   318,
     319,   320,   321,   322,   323,   324,   325,   326,   327,   328,
     329,   330,   332,   333,   335,   336,   337,   338,   340,   342,
     343,   344,   345,   346,   347,   348,   350,   351,   352,   353,
     354,   356,   357,   358,   359,   361,   363,   365,   367,   369,
     371,   373,   375,   376,   377,   378,   379,   380,   381,   382,
     383,   384,   385,   388,   395,   401,   402,   407,   419,   425,
     426,   430,   432,   434,   436,   440,   444,   448,   470,   494,
     527,   528,   529,   530,   531,   534,   546,   550,   552,   557,
     562,   566,   608,   660,   661,   663,   671,   673,   681,   683,
     685,   687,   689,   693,   701,   703,   705,   707,   709,   711,
     713,   715,   717,   719,   721,   723,   725,   729,   731,   735,
     737,   739,   741,   743,   746,   748,   750,   754,   756,   758,
     762,   764,   766,   768,   770,   772,   774,   776,   780,   781,
     782,   783,   784,   785,   788,   789,   793,   795,   804,   806,
     817,   819,   823,   825,   829,   831,   835,   837,   841,   843,
     847,   854,   861,   870,   874,   875,   879,   885,   890,   895,
     901,   902,   904,   906,   911,   912,   916,   917,   921,   924,
     929,   934,   938,   943,   949,   956,   959,   961,   963,   965,
     973,   975,   977,   979,   981,   983,   985,   987,   989,   991,
     993,   997,   999,  1001,  1003,  1006,  1008,  1010,  1012,  1014,
    1016,  1018,  1020,  1023,  1025,  1027,  1029,  1031,  1033,  1035,
    1037,  1039,  1041,  1043,  1045,  1047,  1049,  1051,  1053,  1055,
    1057,  1059,  1063,  1064,  1067,  1070,  1072,  1074,  1076,  1078,
    1080,  1082,  1084,  1086,  1088,  1090,  1093,  1097,  1100,  1103,
    1105,  1107,  1110,  1112,  1114,  1118,  1122,  1128,  1130,  1133,
    1135,  1139,  1141,  1145,  1147,  1152,  1155,  1159,  1163,  1165,
    1167,  1170,  1172,  1173,  1174,  1175,  1177,  1179,  1181,  1183,
    1190,  1192,  1194,  1196,  1198,  1202,  1207,  1213,  1215,  1219,
    1223,  1228,  1242,  1244,  1245,  1247,  1250,  1252,  1254,  1256,
    1262,  1273,  1276,  1277,  1278,  1280,  1284,  1286,  1290,  1300,
    1303,  1312,  1329,  1350,  1355,  1357,  1359,  1361,  1365,  1367,
    1371,  1375,  1379,  1383,  1385,  1389,  1393,  1397,  1401,  1403,
    1405,  1407,  1410,  1415,  1420,  1425,  1433,  1435,  1439,  1441,
    1443,  1446,  1449,  1457,  1467,  1471,  1475,  1479,  1485,  1490,
    1499,  1500,  1503,  1505,  1507,  1509,  1511,  1513,  1515,  1517,
    1519,  1521,  1525,  1529,  1536,  1542,  1546,  1552,  1560,  1566,
    1571,  1574,  1580,  1588,  1590,  1592,  1594,  1596,  1610,  1612,
    1622,  1626,  1628,  1632,  1634,  1638,  1679,  1683,  1689,  1695,
    1697,  1700,  1702,  1704,  1713,  1715,  1717,  1720,  1731,  1733,
    1735,  1738,  1749,  1751,  1753,  1758,  1761,  1762,  1765,  1767,
    1769,  1773,  1774,  1777,  1783,  1785,  1790,  1793,  1794,  1797,
    1800,  1803,  1808,  1809,  1812,  1817,  1818,  1821,  1825,  1826,
    1829,  1835,  1836,  1839,  1843,  1847,  1849,  1853,  1857,  1859,
    1863,  1865,  1868,  1870,  1873,  1877,  1880,  1882,  1886,  1892,
    1894,  1896,  1898,  1902,  1904,  1906,  1907,  1908,  1914,  1916,
    1918,  1921,  1925,  1932,  1934,  1936,  1947,  1958,  1961,  1966,
    1968,  1981,  1983,  1995,  1997,  2000,  2004,  2011,  2014,  2020,
    2026,  2028,  2030,  2032,  2034,  2036,  2038,  2040,  2049,  2052,
    2055,  2058,  2063,  2065,  2067,  2069,  2071,  2073,  2077,  2079,
    2083,  2085,  2089,  2090,  2093,  2095,  2097,  2101,  2102,  2105,
    2107,  2109,  2113,  2114,  2117,  2119,  2121,  2125,  2126,  2129,
    2131,  2133,  2135,  2137,  2141,  2142,  2145,  2147,  2149,  2153,
    2154,  2157,  2159,  2161,  2165,  2166,  2169,  2171,  2173,  2177,
    2178,  2181,  2183,  2185,  2189,  2190,  2193,  2201,  2209,  2219,
    2220,  2221,  2227,  2230,  2234,  2238,  2240,  2246,  2249,  2252,
    2256,  2261,  2269,  2288,  2289,  2292,  2294,  2296,  2300,  2302,
    2304,  2308,  2316,  2326,  2328,  2333,  2335,  2339,  2340,  2343,
    2350,  2373,  2381,  2405,  2411,  2418,  2426,  2439,  2455,  2472,
    2489,  2497,  2505,  2523,  2530,  2537,  2548,  2564,  2571,  2582,
    2590,  2602,  2617,  2630,  2647,  2659,  2672,  2686,  2702,  2719,
    2737,  2750,  2764,  2779,  2796,  2814,  2833,  2841,  2850,  2862,
    2875,  2884,  2895,  2908,  2923,  2933,  2946,  2960,  2977,  2984,
    2994,  3002,  3018,  3027,  3028,  3031,  3037,  3043,  3050,  3058,
    3069,  3071,  3073,  3077,  3079,  3081,  3083,  3085,  3089,  3095,
    3097,  3105,  3107,  3113,  3120,  3127,  3134,  3142,  3149,  3157,
    3165,  3173,  3184,  3186,  3192,  3200,  3203,  3206,  3212,  3214,
    3220,  3230,  3232,  3234,  3236,  3238,  3240,  3246,  3252,  3259,
    3266,  3273,  3280,  3288,  3297,  3306,  3317,  3318,  3320,  3322,
    3324,  3326,  3328,  3330,  3332,  3334,  3336,  3338,  3340,  3342,
    3344,  3346,  3348,  3351,  3353,  3355,  3357,  3359,  3361,  3363,
    3365,  3367,  3369,  3372,  3374,  3380,  3384,  3386,  3388,  3390,
    3394,  3399,  3404,  3406,  3411,  3413,  3415,  3417,  3421,  3425,
    3429,  3431,  3435,  3436,  3440,  3441,  3445,  3450,  3455,  3456,
    3460,  3467,  3474,  3481,  3490,  3491,  3498,  3500,  3503,  3506,
    3511,  3513,  3517,  3522,  3528,  3530,  3534,  3540,  3542,  3546,
    3552,  3557,  3562,  3567,  3574,  3577,  3580,  3582,  3584,  3588,
    3589,  3590,  3593,  3597,  3601,  3603,  3607,  3608,  3611,  3615,
    3616,  3619,  3625,  3627,  3629,  3633,  3639,  3644,  3651,  3655,
    3661,  3667,  3672,  3679,  3684,  3693,  3698,  3712,  3716,  3718,
    3720,  3725,  3728,  3736,  3738,  3740,  3743,  3745,  3747,  3749,
    3751,  3753,  3755,  3757,  3761,  3763,  3765,  3767,  3769,  3771,
    3773,  3775,  3777,  3779,  3781,  3783,  3787,  3789,  3793,  3795,
    3799,  3801,  3803,  3804,  3808,  3814,  3816,  3819,  3823,  3828,
    3830,  3834,  3836,  3840,  3844,  3848,  3852,  3857,  3862,  3867,
    3871,  3876,  3882,  3889,  3897,  3906,  3909,  3913,  3919,  3925,
    3930,  3935,  3940,  3946,  3956,  3969,  3983,  3999,  4015,  4018,
    4022,  4026,  4028,  4030,  4032,  4034,  4036,  4038,  4040,  4042,
    4044,  4046,  4053,  4055,  4062,  4064,  4066,  4068,  4070,  4072,
    4074,  4076,  4078,  4080,  4082,  4084,  4087,  4090,  4092,  4094,
    4096,  4103,  4109,  4115,  4117,  4120,  4122,  4125,  4128,  4131,
    4134,  4136,  4155,  4158,  4160,  4163,  4165,  4168,  4170,  4172,
    4174,  4176,  4178,  4180,  4182,  4185,  4193,  4212,  4228,  4232,
    4237,  4242,  4248,  4252,  4254,  4257,  4259,  4261,  4263,  4265,
    4267,  4269,  4272,  4274,  4276,  4278,  4280,  4304,  4306,  4308,
    4310,  4313,  4316,  4318,  4320,  4322,  4324,  4326,  4328,  4330,
    4332,  4335,  4337,  4339,  4341,  4343,  4346,  4349,  4352,  4355,
    4359,  4363,  4365,  4367,  4370,  4372,  4374,  4378,  4381,  4392,
    4398,  4406,  4413,  4420,  4428,  4437,  4447,  4453,  4459,  4465,
    4471,  4474,  4479,  4496,  4501,  4507,  4514,  4519,  4524,  4529,
    4530,  4531,  4534,  4535,  4538,  4542,  4545,  4550,  4551,  4554,
    4558,  4561,  4572,  4575,  4586,  4592,  4594,  4597,  4609,  4611,
    4614,  4619,  4621,  4623,  4625,  4627,  4630,  4635,  4638,  4643,
    4646,  4649,  4651,  4653,  4655,  4658,  4664,  4670,  4677,  4679,
    4682,  4686,  4692,  4694,  4697,  4702,  4708,  4717,  4721,  4722,
    4724,  4734,  4735,  4739,  4742,  4745,  4747,  4750,  4752,  4756,
    4757,  4762,  4767,  4772,  4782,  4793,  4804,  4809,  4814,  4819,
    4829,  4840,  4851,  4855,  4859,  4863,  4874,  4879,  4884,  4889,
    4899,  4910,  4921,  4926,  4931,  4936,  4942,  4948,  4954,  4959,
    4964,  4969,  4975,  4981,  4987,  4992,  4997,  5002,  5008,  5014,
    5020,  5025,  5030,  5035,  5040,  5045,  5050,  5055,  5060,  5065,
    5070,  5075,  5080,  5085,  5090,  5096,  5102,  5107,  5112,  5118,
    5124,  5129,  5135,  5140,  5146,  5151,  5156,  5161,  5171,  5182,
    5193,  5198,  5203,  5208,  5218,  5229,  5240,  5246,  5255,  5264,
    5274,  5284,  5293,  5302,  5311,  5320,  5329,  5338,  5347,  5356,
    5365,  5374,  5383,  5392,  5401,  5410,  5420,  5430,  5439,  5448,
    5458,  5468,  5477,  5486,  5496,  5506,  5511,  5517,  5526,  5536,
    5545,  5555,  5560,  5566,  5575,  5585,  5594,  5604,  5611,  5618,
    5625,  5632,  5639,  5646,  5653,  5660,  5667,  5674,  5681,  5688,
    5695,  5702,  5709,  5716,  5723,  5728,  5733,  5738,  5744,  5750,
    5759,  5768,  5774,  5780,  5789,  5798,  5803,  5808,  5813,  5818,
    5823,  5831,  5840,  5859,  5863,  5869,  5870,  5872,  5879,  5880,
    5891,  5892,  5903,  5907,  5910,  5917,  5921,  5925,  5936,  5938,
    5940,  5942,  5944,  5946,  5948,  5950,  5952,  5955,  5957,  5959,
    5961,  5964,  5966,  5968,  5970,  5983,  5989,  5996,  6000,  6007,
    6011,  6013,  6015,  6017,  6020,  6022,  6025,  6029,  6031,  6033,
    6036,  6040,  6044,  6048,  6052,  6056,  6060,  6061,  6062,  6063,
    6067,  6071,  6075,  6077,  6080,  6083,  6085,  6088,  6092,  6094,
    6096,  6099,  6101,  6103,  6106,  6110,  6112,  6114,  6116,  6118,
    6120,  6122,  6124,  6127,  6131,  6133,  6135,  6137,  6139,  6141,
    6143,  6145,  6147,  6149,  6151,  6154,  6158,  6161,  6165,  6170,
    6172,  6174,  6176,  6178,  6180,  6182,  6184,  6186,  6189,  6192,
    6194,  6201,  6208,  6209,  6214,  6223,  6225,  6229,  6234,  6236,
    6242,  6246,  6250,  6253,  6257,  6262,  6264,  6269,  6271,  6273,
    6275
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ACTUATORS", "ADJOINTBASIS", "AERO",
  "AEROH", "AEROTYPE", "AGRESSIVETOLERANCES", "ALPROC", "AMAT", "ANALYSIS",
  "ARCLENGTH", "ATTRIBUTES", "ANGULAROUTTYPE", "ARUBBERMAT", "AUGMENT",
  "AUGMENTTYPE", "AUTOTOL", "AUXILIARY", "AVERAGED", "ATDARB", "ACOU",
  "ATDDNB", "ATDROB", "ARPACK", "ATDDIR", "ATDNEU", "ALLOWMECHANISMS",
  "AUXCOARSESOLVER", "ACMECNTL", "ADDEDMASS", "AEROEMBED", "ANDESCLR",
  "ANDESCQR", "ANDESBETAB", "ANDESALPHA", "ANDESBETAM", "AUGMENTED",
  "BLOCKDIAG", "BOFFSET", "BUCKLE", "BGTL", "BMPC", "BINARYINPUT",
  "BINARYOUTPUT", "BLOCKSIZE", "CHECKTOKEN", "COARSESOLVER", "COEF",
  "CFRAMES", "COLLOCATEDTYPE", "CONVECTION", "COMPOSITE", "CONDITION",
  "CONTACT", "CONTROL", "CORNER", "CORNERTYPE", "CURVE", "CCTTOL",
  "CCTSOLVER", "CRHS", "COUPLEDSCALE", "CONTACTSURFACES", "CMPC", "CNORM",
  "COMPLEXOUTTYPE", "CONSTRMAT", "CASES", "CONSTRAINEDSURFACES",
  "CSFRAMES", "CSTYPE", "CONSTANT", "CONWEP", "DAMPING", "DblConstant",
  "DEFAULTPENALTY", "DELETEELEMENTS", "DEM", "DIMASS", "DISP", "DIRECT",
  "DLAMBDA", "DP", "DYNAM", "DETER", "DECOMPOSE", "DECOMPFILE", "DMPC",
  "DEBUGCNTL", "DEBUGICNTL", "DOCLUSTERING", "DOROWCLUSTERING", "ANGLE",
  "DUALBASIS", "DUALRB", "KMEANS", "CRANDOM", "CONSTRAINTS", "MULTIPLIERS",
  "PENALTY", "ELLUMP", "EIGEN", "EFRAMES", "ELSCATTERER", "END",
  "ELHSOMMERFELD", "ENGINEERING", "ETEMP", "EXPLICIT", "EXTFOL", "EPSILON",
  "ELEMENTARYFUNCTIONTYPE", "FABMAT", "FABRICMAP", "FABRICMAT", "FACE",
  "FACOUSTICS", "FETI", "FETI2TYPE", "FETIPREC", "FFP", "FFPDIR", "FITALG",
  "FNAME", "FLAGMN", "FLUX", "FORCE", "FRONTAL", "FETIH",
  "FIELDWEIGHTLIST", "FILTEREIG", "FLUID", "FREEPLAY", "FREQSWEEP",
  "FREQSWEEP1", "FREQSWEEP2", "FREQSWEEPA", "FREQSWEEPAW", "FSGL",
  "FSINTERFACE", "FSISCALING", "FSIELEMENT", "NOLOCALFSISPLITING",
  "FSICORNER", "FFIDEBUG", "FAILSAFE", "FRAMETYPE", "GEPS",
  "GLOBALSEARCHCULL", "GLOBALTOL", "GRAVITY", "GRBM", "GTGSOLVER",
  "GLOBALCRBMTOL", "GROUP", "GROUPTYPE", "GOLDFARBTOL", "HDIRICHLET",
  "HEAT", "HFETI", "HNEUMAN", "HSOMMERFELD", "HFTT", "HELMHOLTZ", "HNBO",
  "HELMMF", "HELMSO", "HSCBO", "HWIBO", "HZEM", "HZEMFILTER", "HLMPC",
  "HERMITIAN", "HESSIAN", "IACC", "IDENTITY", "IDIS", "IDIS6",
  "ILUDROPTOL", "IntConstant", "INTERFACELUMPED", "ITEMP", "ITERTYPE",
  "IVEL", "IMESH", "INCIDENCE", "IHDIRICHLET", "IHDSWEEP", "IHNEUMANN",
  "ISOLVERTYPE", "INPC", "INFINTY", "JACOBI", "KEYLETTER", "KRYLOVTYPE",
  "KIRLOC", "LAYC", "LAYN", "LAYD", "LAYO", "LAYMAT", "LFACTOR", "LISRBM",
  "LMPC", "LOAD", "LOADCASE", "LOBPCG", "LOCALREDUCEDORDERBASES",
  "LOCALSOLVER", "LINESEARCH", "LUMPED", "KSPARAM", "KSMAX", "MASS",
  "MASSAUGMENTATION", "MATERIALS", "MATLAB", "MAXITR", "MAXELEM",
  "MAXORTHO", "MAXVEC", "MODAL", "MPCPRECNO", "MPCPRECNOID", "MPCTYPE",
  "MPCTYPEID", "MPCSCALING", "MPCELEMENT", "MPCBLOCKID", "MPCBLK_OVERLAP",
  "MFTT", "MRHS", "MPCCHECK", "MUMPSICNTL", "MUMPSCNTL", "MUMPSMINEQ",
  "MUMPSSTRIDE", "MECH", "MODDAMP", "MODEFILTER", "MOMENTTYPE", "MPROJECT",
  "MAXIMUM", "NOWARPEDVOLUME", "NDTYPE", "NEIGPA", "NEWMARK", "NewLine",
  "NEWTON", "NL", "NLMAT", "NLPREC", "NLMEMPTYPE", "NOCOARSE", "NODETOKEN",
  "NOMULTIPLEINTERACTIONS", "NONINPC", "NSBSPV", "NLTOL", "NUMCGM",
  "NONORMALSMOOTHING", "NOSECONDARY", "NOGHOSTING", "NFRAMES",
  "NORMALSMOOTHINGDISTANCE", "SENSITIVITY", "SENSITIVITYMETHOD",
  "SKIPPHYSICALFACES", "OUTPUT", "OUTPUT6", "OUTPUTFRAME",
  "OLDDYNAMICSEARCH", "QSTATIC", "QLOAD", "QUASISTATIC", "PITA",
  "PITADISP6", "PITAVEL6", "NOFORCE", "MDPITA", "GLOBALBASES",
  "LOCALBASES", "TIMEREVERSIBLE", "REMOTECOARSE", "ORTHOPROJTOL",
  "READINITSEED", "JUMPCVG", "JUMPOUTPUT", "PRECNO", "PRECONDITIONER",
  "PRELOAD", "PRESSURE", "PRINTMATLAB", "printmatlab", "PRINTNUMBER",
  "PROJ", "PIVOT", "PRECTYPE", "PRECTYPEID", "PICKANYCORNER", "PADEPIVOT",
  "PROPORTIONING", "PLOAD", "PADEPOLES", "POINTSOURCE", "PLANEWAVE",
  "PTOL", "PLANTOL", "PMAXIT", "PIECEWISE", "PARAMETERS", "RADIATION",
  "RBMFILTER", "RBMSET", "READMODE", "READSENSITIVITY", "REBUILD",
  "RESOLUTIONMETHOD", "REVERSEORDER", "REDFOL", "RENUM", "RENUMBERID",
  "REORTHO", "RESTART", "RECONS", "RECONSALG", "REBUILDCCT", "RANDOM",
  "RPROP", "RNORM", "REVERSENORMALS", "ROBTYPE", "ROTVECOUTTYPE",
  "RESCALING", "RUBDAFT", "SCALING", "SCALINGTYPE", "SDETAFT", "SENSORS",
  "SHELLFABRICMAP", "SHELLFABRICMAT", "SHELLSIMPLELOFTING", "SOLVERCNTL",
  "SOLVERHANDLE", "SOLVERTYPE", "SHIFT", "SHARPNONSHARPANGLE",
  "SPOOLESTAU", "SPOOLESSEED", "SPOOLESMAXSIZE", "SPOOLESMAXDOMAINSIZE",
  "SPOOLESMAXZEROS", "SPOOLESMSGLVL", "SPOOLESSCALE", "SPOOLESPIVOT",
  "SPOOLESRENUM", "SPARSEMAXSUP", "SPARSEDEFBLK", "SPARSERENUM", "STATS",
  "STRESSID", "SUBCYCLE", "SUBSPACE", "SURFACE", "STR_THERM_OPTION",
  "SAVEMEMCOARSE", "SPACEDIMENSION", "SCATTERER", "STAGTOL", "SCALED",
  "SWITCH", "STABLE", "SUBTYPE", "STEP", "SOWER", "SHELLTHICKNESS", "SURF",
  "SPRINGMAT", "SENSITIVITYID", "TABLE", "TANGENT", "TDENFORCE", "TEMP",
  "TIME", "TOLEIG", "TOLFETI", "TOLJAC", "TOLPCG", "TOLSEN", "TOPFILE",
  "TOPOLOGY", "TRBM", "TRBMlc", "THERMOE", "THERMOH", "RATIOTOLSEN",
  "TETT", "TOLCGM", "TURKEL", "TIEDSURFACES", "THETA", "PROJSOL", "CENTER",
  "POSELEM", "HOTSTART", "HRC", "THIRDNODE", "THERMMAT", "TDENFORC",
  "TESTULRICH", "THRU", "TRIVIAL", "NUMROMCPUS", "THICKNESSGROUPLIST",
  "USE", "USERDEFINEDISP", "USERDEFINEFORCE", "UPROJ", "UNSYMMETRIC",
  "USING", "USESCOTCH", "VERBOSE", "VERSION", "WETCORNERS", "YMTT",
  "SS1DT", "SS2DT", "YMST", "YSST", "YSSRT", "ZERO", "BINARY", "GEOMETRY",
  "DECOMPOSITION", "GLOBAL", "MATCHER", "CPUMAP", "NODALCONTACT", "MODE",
  "FRIC", "GAP", "OUTERLOOP", "EDGEWS", "WAVETYPE", "ORTHOTOL", "IMPE",
  "FREQ", "DPH", "WAVEMETHOD", "MATSPEC", "MATUSAGE", "BILINEARPLASTIC",
  "FINITESTRAINPLASTIC", "CRUSHABLEFOAM", "LINEARELASTIC",
  "STVENANTKIRCHHOFF", "TULERBUTCHER", "LINPLSTRESS", "READ", "OPTCTV",
  "ISOTROPICLINEARELASTIC", "VISCOLINEARELASTIC", "VISCOSTVENANTKIRCHHOFF",
  "NEOHOOKEAN", "VISCONEOHOOKEAN", "ISOTROPICLINEARELASTICJ2PLASTIC",
  "ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS", "HYPERELASTIC",
  "MOONEYRIVLIN", "VISCOMOONEYRIVLIN", "HENCKY", "OGDEN", "SIMOELASTIC",
  "SIMOPLASTIC", "LOGSTRAINPLASTIC", "SVKPLSTRESS", "VISCOLINPLSTRESS",
  "VISCOSVKPLSTRESS", "VISCOFABRICMAP", "SHELLVISCOFABRICMAP",
  "VISCOFABRICMAT", "SHELLVISCOFABRICMAT", "PARTITIONGAP",
  "PLANESTRESSLINEAR", "PLANESTRESSSTVENANTKIRCHHOFF",
  "PLANESTRESSNEOHOOKEAN", "PLANESTRESSMOONEYRIVLIN",
  "PLANESTRESSBILINEARPLASTIC", "PLANESTRESSFINITESTRAINPLASTIC",
  "PLANESTRESSVISCOLINEARELASTIC", "PLANESTRESSVISCOSTVENANTKIRCHHOFF",
  "PLANESTRESSVISCONEOHOOKEAN", "PLANESTRESSVISCOMOONEYRIVLIN",
  "SURFACETOPOLOGY", "MORTARTIED", "MORTARSCALING",
  "MORTARINTEGRATIONRULE", "SEARCHTOL", "STDMORTAR", "DUALMORTAR",
  "WETINTERFACE", "NSUBS", "EXITAFTERDEC", "SKIP", "ROBCSOLVE",
  "RANDOMSAMPLE", "OUTPUTMEMORY", "OUTPUTWEIGHT", "SOLVER",
  "SPNNLSSOLVERTYPE", "MAXSIZE", "CLUSTERSOLVER", "CLUSTERSOLVERTYPE",
  "WEIGHTLIST", "GMRESRESIDUAL", "SLOSH", "SLGRAV", "SLZEM", "SLZEMFILTER",
  "PDIR", "HEFSB", "HEFRS", "HEINTERFACE", "SNAPFI", "VELSNAPFI",
  "ACCSNAPFI", "DSVSNAPFI", "MUVSNAPFI", "PODROB", "ROMENERGY", "TRNVCT",
  "OFFSET", "ORTHOG", "SVDTOKEN", "CONVERSIONTOKEN", "CONVFI", "ROMRES",
  "SAMPLING", "SNAPSHOTPROJECT", "PODSIZEMAX", "REFSUBTRACT", "TOLER",
  "NORMALIZETOKEN", "FNUMBER", "SNAPWEIGHT", "ROBFI", "STAVCT", "VELVCT",
  "ACCVCT", "CONWEPCFG", "SCALEPOSCOORDS", "NODEPOSCOORDS",
  "MESHSCALEFACTOR", "PSEUDOGNAT", "PSEUDOGNATELEM", "USENMF", "USENMFC",
  "USEGREEDY", "USEPQN", "FILTERROWS", "VECTORNORM", "REBUILDFORCE",
  "REBUILDCONSTRAINT", "SAMPNODESLOT", "REDUCEDSTIFFNESS", "UDEIMBASIS",
  "FORCEROB", "CONSTRAINTROB", "DEIMINDICES", "UDEIMINDICES",
  "SVDFORCESNAP", "SVDCONSTRAINTSNAP", "USEMASSNORMALIZEDBASIS",
  "USECONSTANTMASS", "ONLINEMASSNORMALIZEBASIS", "STACKED",
  "NUMTHICKNESSGROUP", "STRESSNODELIST", "DISPNODELIST", "RELAXATIONSEN",
  "QRFACTORIZATION", "QMATRIX", "RMATRIX", "XMATRIX", "EIGENVALUE",
  "NPMAX", "BSSPLH", "PGSPLH", "LIB", "$accept", "FinalizedData", "All",
  "Component", "Noninpc", "Inpc", "Group", "Random", "Impe",
  "ImpeDampInfo", "PadePivotInfo", "PadePolesInfo", "FreqSweep",
  "ReconsInfo", "BinarySpec", "AnalysisInfo", "Decompose", "WeightList",
  "FieldWeightList", "MFTTInfo", "HFTTInfo", "LoadCase", "Composites",
  "Cframes", "CoefInfo", "CoefList", "LaycInfo", "LaynInfo", "LaydInfo",
  "LayoInfo", "LayData", "LayoData", "LayMat", "LayMatData", "DiscrMasses",
  "Gravity", "Restart", "LoadCInfo", "UseCInfo", "SensorLocations",
  "ActuatorLocations", "UsdfLocations", "UsddLocations", "Output",
  "OutInfo", "DynInfo", "PrintMat", "DeleteElements", "DeleteElementsList",
  "SloshInfo", "MassInfo", "CondInfo", "TopInfo", "ModalInfo", "DynamInfo",
  "Conwep", "ConwepData", "TimeIntegration", "NewmarkSecondOrder",
  "NewmarkFirstOrder", "QstaticInfo", "QMechInfo", "QHeatInfo", "AeroInfo",
  "AeroEmbeddedSurfaceInfo", "AeroHeatInfo", "ThermohInfo", "ThermoeInfo",
  "ModeInfo", "HzemInfo", "SlzemInfo", "RbmTolerance", "ToleranceInfo",
  "ModeFilterInfo", "RbmFilterInfo", "RbmList", "HzemFilterInfo",
  "SlzemFilterInfo", "TimeInfo", "ParallelInTimeInfo",
  "ParallelInTimeOptions", "ParallelInTimeKeyWord", "DampInfo",
  "ComplexDirichletBC", "IComplexDirichletBC", "IComplexDirichletBCSweep",
  "DirectionVector", "IComplexNeumannBC", "DirichletBC",
  "ConstrainedSurfaces", "HEVDirichletBC", "HEVDBCDataList", "HEVDBC_Data",
  "HEVFRSBC", "HEVFRSBCList", "HEVFRSBCElem", "TempDirichletBC",
  "TempNeumanBC", "TempConvection", "TempRadiation", "HelmHoltzBC",
  "SommerfeldBCDataList", "SommerfeldBC_Data", "EleHelmHoltzBC",
  "SommerElement", "SommNodeNums", "Scatterer", "ScattererEleList",
  "ScattererEle", "EleScatterer", "ScatterElement", "NeumScatterer",
  "NeumElement", "WetScatterer", "WetInterfaceElement", "HelmScatterer",
  "HelmScattererElement", "PBC_Data", "PBCDataList", "AtdDirScatterer",
  "AtdNeuScatterer", "AtdArbScatterer", "AtdNeumScatterer",
  "AtdRobinScatterer", "FarFieldPattern", "FarFieldPatternDirs",
  "ReadModeInfo", "Mode", "IDisp", "IDisp6", "IDisp6Pita", "IVel6Pita",
  "IVel", "ITemp", "ETemp", "NeumanBC", "ModalNeumanBC", "BCDataList",
  "ModalValList", "TBCDataList", "YMTTable", "YMTTList", "TETTable",
  "TETTList", "SS1DTable", "SS1DTList", "SS2DTable", "SS2DTList",
  "YSSTable", "YSSTList", "YSSRTable", "YSSRTList", "YMSTable", "YMSTList",
  "SDETAFTable", "SDETAFList", "RUBDAFTable", "RUBDAFList", "LMPConstrain",
  "ModalLMPConstrain", "MPCList", "MPCHeader", "MPCLine",
  "ComplexLMPConstrain", "ComplexMPCList", "ComplexMPCHeader",
  "ComplexMPCLine", "ComplexNeumanBC", "ComplexBCDataList", "Materials",
  "MatData", "FreeplayProps", "ElemSet", "FaceSet", "MortarCondition",
  "WetInterface", "TiedSurfaces", "FSInterface", "HEVibInfo",
  "HEVInterfaceElement", "HEInterface", "ContactSurfaces",
  "ContactSurfacesInfo", "Parameters", "NodeSet", "Node", "Element",
  "NodeNums", "BC_Data", "ModalVal", "TBC_Data", "ComplexBC_Data",
  "ConstrainedSurfaceFrameDList", "FrameDList", "Frame", "NodalFrameDList",
  "NodalFrame", "BoffsetList", "Attributes", "Ellump",
  "LocalReducedOrderBases", "LocalBasesAuxi", "LocalBasesCent",
  "ReducedStiffness", "UDeimBasis", "SampNodeSlot", "Pressure", "Lumped",
  "MassAugmentation", "Preload", "Sensitivity", "DispNode", "StressNode",
  "ThicknessGroup", "Statics", "CasesList", "SolverMethod", "Solvercntl",
  "Solver", "OldHelmInfo", "FAcousticData", "Constraints",
  "ConstraintOptionsData", "HelmInfo", "IncidenceList", "IncidenceVector",
  "KirchhoffLocations", "FFPDirList", "FFPDirVector", "HelmMFInfo",
  "FAcousticDataMF", "HelmSOInfo", "FAcousticDataSO", "DEMInfo", "AlProc",
  "NLInfo", "NewtonInfo", "OrthoInfo", "Control", "NodalContact",
  "MatSpec", "MatUsage", "FloatList", "StringList", "Renumbering",
  "SvdToken", "SvdOption", "DeimIndices", "UDeimIndices", "Sampling",
  "SnapshotProject", "SamplingOption", "ConwepConfig", "MeshScaleFactor",
  "ScalePosCoords", "NodePosCoords", "ConversionToken", "ConversionOption",
  "Integer", "Float", YY_NULLPTR
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   614,
     615,   616,   617,   618,   619,   620,   621,   622,   623,   624,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     635,   636,   637,   638,   639,   640,   641,   642,   643,   644,
     645,   646,   647,   648,   649,   650,   651,   652,   653,   654,
     655,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     665,   666,   667,   668,   669,   670,   671,   672,   673,   674,
     675,   676,   677,   678,   679,   680,   681,   682,   683,   684,
     685,   686,   687,   688,   689,   690,   691,   692,   693,   694,
     695,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   706,   707,   708,   709,   710,   711,   712,   713,   714,
     715,   716,   717,   718,   719,   720,   721,   722,   723,   724,
     725,   726,   727,   728,   729,   730,   731,   732,   733,   734,
     735,   736,   737,   738,   739,   740,   741,   742,   743,   744,
     745,   746,   747,   748,   749,   750,   751,   752,   753,   754,
     755,   756,   757,   758,   759,   760,   761,   762,   763,   764,
     765,   766,   767,   768,   769,   770,   771,   772,   773,   774,
     775,   776,   777,   778,   779,   780,   781,   782,   783,   784,
     785,   786,   787,   788,   789,   790,   791,   792,   793,   794,
     795,   796,   797,   798,   799,   800,   801,   802,   803,   804,
     805,   806,   807,   808,   809,   810,   811,   812,   813,   814,
     815,   816,   817,   818,   819,   820,   821,   822,   823,   824,
     825,   826,   827,   828,   829,   830,   831,   832,   833,   834,
     835,   836,   837,   838,   839,   840,   841,   842,   843,   844,
     845
};
# endif

#define YYPACT_NINF -2827

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-2827)))

#define YYTABLE_NINF -1426

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    7993,  -157,    97,  -144,   854,   891,  2286,  2286,  2286,   870,
     -96,   999,   854,   -66,   -63,   854,   117,  1025,   196,   215,
     275,   -14,   355,   369,   393,   416,   446,   466,   472,   476,
     514,   238,  1042,  1063,   523,   559,   576,   605,   641,   661,
     854,   854,  1069,  1095,   667,  -258,   690,   737,   748,   990,
     769,   789,  1151,   795,  1223,   797,   829,   834,   846,   859,
     873,   879,   885,   914,   -87,   991,   916,   919,   854,   921,
    1306,   928,   943,   854,   955,   963,   993,    72,  1319,   968,
     969,  1322,   -68,   106,   972,  1432,  1438,   854,   979,  1000,
     980,   854,   986,   987,   681,  1051,   992,  1002,   854,   854,
    1003,  1009,  1445,  1139,   854,   854,  1029,  1446,  1463,   337,
    1036,  1040,  1047,  1056,  1062,  1071,  1072,   854,  2286,  1074,
    1466,  1081,  1085,  2286,  2286,  1093,  1100,  1102,  1109,  1111,
    1112,  1113,  1122,  1253,  1131,  1141,  1142,  1143,  1144,  1145,
    1204,  1206,  1212,    29,  1473,  1483,  1215,  1216,   854,   854,
     854,  1219,  1227,  2286,  1231,  1232,  1237,  1242,  1248,  1255,
    1261,  1263,  1266,  1267,  2286,  1274,  1276,  1277,  1278,  1282,
     585,  1220,  7423, -2827, -2827, -2827,  1345,   854,     0,    75,
   -2827,    -8,   854,    34,  2286,  2286,   854,   931,   854,   854,
    1080,  2286,  1392, -2827, -2827, -2827, -2827, -2827,   -51,   174,
    1318, -2827,  2286, -2827, -2827, -2827, -2827,   133, -2827,   708,
       4, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
    1336, -2827, -2827, -2827, -2827, -2827,   854, -2827,   232,   854,
   -2827, -2827,   338,   359,   389,   854, -2827,   854, -2827,   854,
     854,   854,   854, -2827, -2827,   854,   854,   854, -2827, -2827,
      85,  1321,   854,   854,   854,  1325,  1330, -2827,   419, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
    2286, -2827, -2827,   854,   854,   854, -2827, -2827,   854,   854,
     854,   854,  -128,   138,   854,   854,   854,   854,   854,   854,
     462,    13,  2286,   898,   854,   424, -2827, -2827,   854,   270,
     558, -2827, -2827, -2827,   166,  2286, -2827, -2827, -2827,    10,
   -2827, -2827,   854,  -147,   854, -2827,   285,  5165,  5165,  1817,
    5165, -2827,  -389,   854,  1296,  1549,  1551, -2827, -2827,  1305,
   -2827,  1320, -2827, -2827, -2827, -2827,  1331,  1342,  2286,  1504,
   -2827,  2286,   854,   854,  1346,  1348, -2827, -2827,  1038, -2827,
   -2827,  1349, -2827,  2286,  1454, -2827,   854, -2827, -2827,  2286,
    2286, -2827, -2827, -2827, -2827, -2827,  2286,  2286,  1499,  2286,
    1268,   376, -2827,  1360, -2827,  1369, -2827,   854,   854,   854,
   -2827,  2286,  1606,  1386, -2827,  1387,  1415,  1390, -2827,  1396,
   -2827, -2827, -2827,  2286,  2286, -2827,   854,   854,  1397,   854,
   -2827,  1400, -2827,   854,  2286,  2286,   854,   854, -2827, -2827,
     854,   854,  1402, -2827,  1406,   854,   854,  1413,  2286,   854,
    1417,  2286,   854,  1420,  2286, -2827,  1285,  2286,  1425, -2827,
      -2, -2827, -2827, -2827,  1426,  1430, -2827,  1435,   854, -2827,
    1437, -2827,  1442,  1443, -2827,   854,  2286,   854,  1444, -2827,
   -2827,  1447,  1617, -2827,  1452,  1455,  1623, -2827,  1467, -2827,
     854,  1469,  1471,   854, -2827, -2827,  1472,  -207,  1477,  1482,
   -2827, -2827,  1489, -2827,  1672,  1676,   854,  1419, -2827, -2827,
   -2827,  1364,  1687,   854,  1505,  1506, -2827, -2827,  2286,   854,
   -2827,  1509,  1514, -2827,   854,  2286, -2827, -2827,  1705, -2827,
   -2827,  1516, -2827,   854,  1708,  1711,  1407,  1747,  1752,  1757,
   -2827, -2827,   854, -2827,  1575, -2827,  1577, -2827, -2827,   194,
     854,  1699, -2827, -2827,  1580, -2827, -2827,   854,   854,   854,
   -2827, -2827, -2827, -2827, -2827,  2286, -2827, -2827, -2827, -2827,
   -2827,  1585, -2827, -2827, -2827,   450,  1511,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  1041,  1512,  2286, -2827, -2827, -2827,
    2286, -2827,  1479,  1486,  1718,  1720,  1729,  1735,  1737,  1610,
    1614,  1740,  1619,  1700,  1620,   854,  1621,  1639,  1642,  1643,
    2286,   854,   854,   854,   854,  2286,  2286,  2385,   854,   854,
     854,   854,   854, -2827,   854,   854,   854,   854, -2827,   -11,
   -2827,  2286,  1649,   854,  2286,   696,   854, -2827,  -120,   756,
     796,  1774,  1776,  1779,  1780,   198,   854,  1703,   854,  2286,
    1044,  2286,  2286,  1541,  1712,  2286,  1673,  1674,  1677,  1562,
    -117,  1567,  2286,  1570,   439, -2827, -2827, -2827,  2286,  1715,
    2286, -2827, -2827, -2827,  1681,   854,  1828,  1719,   854,  1837,
   -2827,  2286,   854, -2827,   205,  1886,   854,   122,   854,  2286,
     854,  2286,  2286, -2827,   854, -2827,   854, -2827,   854, -2827,
     854, -2827,   854, -2827, -2827, -2827,  1835, -2827,   -58,  1854,
    2286,  2286,  2286,  1863,  1713,   854, -2827,   229,  1724, -2827,
      76, -2827,   854,   854,   854,   854, -2827,   854,   854, -2827,
      -1,   854,  1733,  1738,   854,  2286,  2286,  2286,  2286,  2286,
    1739,   854,  1770,  1782,   854,  1783,  1784,  1785,  1790,  2286,
    1793,  1805,   854,  1810,  2286,  1812,  2286,   854, -2827,  2286,
   -2827, -2827, -2827,  1171,   854,   -52,  1814,  1815,  2286,   854,
    1940, -2827, -2827,  1818,   803,   854,  1823,   854,   854,  1078,
     140,  2286,  2286,  1609,  1948,  2286,  2286,   854,   854,  1827,
    2286,  1836,  2286,   854,  1896,  1722,   -60,  1899,  1911,   854,
    1045,   854,   640,   854, -2827,  6402, -2827, -2827, -2827,  2286,
    1841,  2286,   854,   854,  1082,  2286,   113,   854,  1844,  1860,
    2286,   854,   854,  1869,  1919, -2827,   854,  1994,  4736,   854,
     854,   854,   854,  1755,   854,   854,  1758,  1630, -2827, -2827,
   -2827, -2827, -2827,  2286,   854,  2025,   -53,  1781, -2827,  1897,
     854,   854, -2827,   854,  1903,  2286,  2037,  1795,  2286,   854,
    1797,  1800,  1803,  1831,  1834,  1838,  1843,  1847,   854,   854,
    1678,  2286,  2045,  2091,   854,   854,  2286,  1852,  1884,   854,
    1885,  2101,  1887,  1889,  2156,  2157,  1914,  1917,  1923,  1924,
    1925,   854,   854,   854,  1950,  1958,  2053,  2055,  2056,  2286,
    2182,   430,  2057,  2183,   854,  2060,   854, -2827,   278, -2827,
    1084,  2286, -2827, -2827, -2827, -2827,  2286, -2827,  2061,  1943,
   -2827,   854,  2286,   854,   854, -2827, -2827,  1091, -2827,   854,
    2062,  2063,  2064,  2286,  1944, -2827,  2286,  2286, -2827,   433,
     854, -2827, -2827, -2827, -2827, -2827,   854, -2827,  2286, -2827,
    2065, -2827,  2067,  2286, -2827,  1960,  2094, -2827,  1099,  2286,
     854, -2827,   854,   854,   854,   854, -2827,   854, -2827, -2827,
   -2827,  2072, -2827,  2079, -2827, -2827,   854,   854,   772,   854,
   -2827, -2827,   854,   854,  2286,  2286,  2286, -2827,  2286,  2084,
   -2827,  2286,   854,   854,   854,  1121,  2286, -2827,  1964, -2827,
    1966, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,  2286,
     854, -2827, -2827, -2827,  2090, -2827, -2827, -2827,  2092, -2827,
     854, -2827, -2827,   854, -2827, -2827,  2095,  2286,  2286, -2827,
   -2827,  2096, -2827,  2097, -2827,  2102,   854,   178,   854,  1186,
     854,   341, -2827,  2286,   854, -2827,   854, -2827, -2827, -2827,
    2103,   854,  1773, -2827,   854,   854,  1865,   854,  1990,   854,
     346,   854,  1992,   854,  2019,   854,  2281,  2107, -2827, -2827,
   -2827,  2116,  2286,  -183, -2827,  2118, -2827,   854, -2827,  1148,
   -2827,   854, -2827,   854,  2286, -2827,   854,   854,  2286,  2286,
    1154,  2286,  2286,  2286,  2286,  2127, -2827,  2286,   854,  2133,
     854,   482,   525,  2134,  2138,  2152,  2154,  2159, -2827, -2827,
    2160, -2827, -2827,  2167, -2827,  2194, -2827, -2827, -2827, -2827,
    2205,  2210,  2211,  2220,  2221,  2222,  2226,   854,   854,  2228,
     547,  2229,  2232,  2235,  2238, -2827,  2286, -2827, -2827,   854,
   -2827,   854,  2286,  2286,  2218,   898,  2286,   556,  2239, -2827,
     854,   854,   692,   854,   854,   854,   854, -2827, -2827, -2827,
   -2827, -2827, -2827,   854, -2827,   854, -2827,   854, -2827, -2827,
    2117, -2827, -2827, -2827,  2286,  2241, -2827, -2827, -2827,  2247,
    2286, -2827,  2286,  2286,  2250, -2827,  2257, -2827, -2827,  2258,
     854, -2827, -2827,  2261,  2269,  1203,  1985,  2286,  2270,   854,
   -2827, -2827,  2286, -2827,  2276,  2286, -2827,  2283,  2284, -2827,
   -2827, -2827, -2827, -2827,  2286, -2827,   854,   854,  2286,   854,
    2287,  2286,  2290,  2286,  2286,  2286,   854,   854,   854,   854,
     854,   854,  2418,  2419,   854,  2293,  2286,  2286,  2286,   854,
    2302,   854, -2827,   854, -2827,  2286,  2286,    -5,   854,  2286,
    2286,  2286,  2286,   854,   854,   854,   854,   854,   854, -2827,
     854,  -136,   854, -2827, -2827,  2305,  2318,  2324,  2325,  2327,
    2331, -2827,  2335, -2827, -2827,  2337, -2827, -2827, -2827, -2827,
    2338, -2827, -2827,  2342, -2827,  2353, -2827,  2354,  2375,  1205,
    2286,  2286,  2286,  2286,   408, -2827, -2827,  2376,   854,  2377,
   -2827,  2380, -2827,  2074,   155,   854,  2203,   300,   854,  1269,
    2381,  2383,  2387,  2391,  2400,  2404,  2086, -2827,  2093, -2827,
     854,  2420, -2827,  2421,  2104, -2827, -2827,  2423,  2426,  2430,
   -2827,  2441,  2442,   878, -2827,  2286, -2827,  1270,  2444, -2827,
    2286,  2445, -2827,  2447,  2470,  2476,  1289,     6,   341,  2295,
     341,   971,  2286,   341,  2296,  2371,  2477,   517,  2501,  2505,
     854,  2286,  2286,  2286,  2286,  2511,  1298,   341,   854,   854,
     854,   983,  1162,   533,  2512,   854,   854,   854,   854,   854,
    2110,  2513,   854,   612,   592,   854,   854,   712,   854,  2286,
    2286,  2286,   854,  2394,   545,  2286,   854,   854,   854,  2286,
     854,   854,   854,   854,   854,   854,  2532,   854,  2286,   854,
    2286,  2286,  2286,  2286,   854,  2122,   854,  2534,  2601,  2286,
     480,  2286, -2827,  2286,  2286,  2260, -2827,  2538,  2540, -2827,
    2304, -2827,  2544, -2827, -2827,  1290,  2314,  2549, -2827, -2827,
    2550,  2286,  2675,  2286,  2286, -2827,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,  2341, -2827,  2708,  2709,
   -2827, -2827,   854, -2827,  2286,  2679,  2679,  2679,  2679,  2679,
   -2827, -2827, -2827, -2827, -2827, -2827,  2434,  2679, -2827,   898,
    2286,  2286, -2827,  2560,   854, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
     854,  2692, -2827,   854, -2827, -2827, -2827,   898, -2827,  2693,
   -2827, -2827,   854,   854, -2827, -2827, -2827, -2827, -2827, -2827,
     854,   854, -2827, -2827, -2827, -2827, -2827,  2577,  2286,   597,
   -2827,   854, -2827, -2827, -2827,   293,   854, -2827,  1297,  2286,
    2578, -2827, -2827,  2584, -2827,  2589,   854, -2827,  2593,  2596,
     854, -2827, -2827,  2286, -2827,  2286, -2827, -2827,  2286, -2827,
    2599, -2827, -2827,  2286, -2827,  2286,   854,  2608,  2411, -2827,
    2412,  2610, -2827,  2286,   854, -2827,   898, -2827, -2827,   854,
   -2827,   718, -2827,   854,  2286, -2827,  2612,  2286, -2827,  2286,
    1380,  2286,  2286, -2827,  2286,  2613,   854, -2827,   892, -2827,
     854, -2827,     2,  2619,   176,  2621,  2622,  2624, -2827, -2827,
    2413,   854, -2827,  2286,  2286, -2827, -2827,  2527,  2625,   854,
    2286,  2627,   854,  2286,  6402,  2628, -2827,   898, -2827,  2629,
     854,  2286,  2634,   854,  2286,  2636,   854,  2286,    92,   854,
   -2827,  2637,   854,  2286,  2639,   854,  2286,  2640,   854,  2286,
   -2827, -2827,   203, -2827,  2286,  -199,  -198, -2827, -2827, -2827,
    2642, -2827,   854,  2643,   854,  2414,  2286,  2644, -2827,   854,
     854,   854,   854,   854, -2827,  2646,  2415, -2827,  2647,  2651,
   -2827,  2653, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,  2655,  2672,
   -2827,   854,  2677,  2286, -2827, -2827, -2827, -2827,  2286,  2286,
    2682,  2286,  2286,  2685,  2286,  2689,  2719,  2721, -2827, -2827,
     395,  2723, -2827,   694,   652,  2855,   821,  2856, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827,  2286,  2738,   833, -2827, -2827,
   -2827,   854, -2827,  2286, -2827, -2827, -2827,  2286, -2827,  2286,
    2286, -2827, -2827,  2286, -2827,   854, -2827, -2827,  2286,   420,
     854,  2742,   554, -2827,  2746, -2827,  2286,  2286,  2286,  2427,
   -2827,  2500,  2526,  2590,  2597,  2638,   854,   854,   854,   854,
    2286,  2286,  2286,   854,   854,   854,   425,  2286,  1381,  2286,
    2286, -2827,  2286,   854,    74,  2286,  2286,  2286,   218,  2286,
    2286,  2614, -2827,  2616,   490,  1383,  2620,  1384,  2748, -2827,
     854,   854, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827,  1401,  2286,  2286,  2286,
     -44,  2750, -2827,   169,   -17, -2827,  2896, -2827, -2827, -2827,
    2757,   854,   481,  2286,   854,   489, -2827,  2758,  2286, -2827,
    2286, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827,   854,   854, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827,  2774, -2827,  2631,  2635, -2827,  2716, -2827,  2775,
   -2827, -2827, -2827, -2827, -2827,  1408,   734,  2779,  -137, -2827,
    2780, -2827,  2782,  2717,  2789, -2827,  2790,  2799, -2827,  2800,
    2803, -2827, -2827,  2806,  2808,  2829,  2831,  2839, -2827,  2842,
    2863, -2827,  2286,  2865,  2873,   894,   745,  2874,  2875,  2879,
    2883, -2827,  2887,   854,  2286,  2890,  2891, -2827,  2898, -2827,
    2899,  2901,  2902,  2903,  2910, -2827,  2911,  2914,  2915,  2916,
    2922,  2927,  2286,  2286,   854,  2929,  2931,  2932,  2935,  2940,
    2950,  2954,  2956,  2967,  2972,  2984,  2986,  2998,  3007, -2827,
    3016,  1409,  3019,  1412,  3021,  3023,  3026,  3027, -2827,  3028,
    3032, -2827,   522,  1418, -2827,  3034,  3040,  2718,  2286, -2827,
    3043, -2827, -2827, -2827,  1422, -2827, -2827,  1423, -2827,  2798,
   -2827, -2827,  2286,  3045,   854,   854,  1429,   854,   854,  2286,
    2286,  2286,   -50,   -30,  2286,  2286,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,   -27,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,   854,   854,   854,   854,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,  2286, -2827,  3047, -2827,
   -2827, -2827,  2787, -2827, -2827,   854, -2827,  2286,   854, -2827,
   -2827, -2827,  3183,   854,   854, -2827,  3187,   854,   854, -2827,
   -2827, -2827,  2286, -2827, -2827,   854,   492, -2827,  1431,  3065,
   -2827, -2827, -2827, -2827, -2827,  3076,  2286,  2286, -2827, -2827,
   -2827,  2286,   854,   854,   854,  3086, -2827,  3089,  2286,  1460,
     898,  3093, -2827,  2286,  2286, -2827, -2827,  2286, -2827,  3095,
    3101,  3102,  3103, -2827, -2827,  2286, -2827,   854,   532, -2827,
   -2827,  2988, -2827,   854, -2827,  2728, -2827,  3107,   854,  3113,
    2286,  3117,  2286,  2286,  3118,  3119, -2827,   898,  3121,  2286,
    3122,  3123,  2286,  3128,  3129,  2286,  3147,  3154,  3033, -2827,
     222,  1475,  2286,  3157,  3163,  2286,  3164,  3166,  2286,  3168,
    3170, -2827,  3179,  1484, -2827,  2286, -2827,  2286, -2827,  2801,
   -2827,  3185, -2827,  3193,  3194, -2827,   854,  3199,  3203,  3105,
    3127, -2827, -2827,   854, -2827, -2827, -2827, -2827, -2827,   854,
     854,   854,  2286,  2286,  2286, -2827,  2286,  2286, -2827,  3206,
   -2827, -2827, -2827,   854,  3207, -2827,  3210, -2827,   854, -2827,
     854,   854, -2827,   854,  2286, -2827,  3214, -2827, -2827,  3221,
     854,  1530,  3224,   854,  1547,  3226,   854,  2286,  3227, -2827,
     854,  3233, -2827,  2286,  3235,  3236, -2827, -2827, -2827, -2827,
   -2827, -2827,  3239,  1558,  1576,   854,  1587,  2286,  2286,   854,
     854,  2286,  1607, -2827,  2286,  2286,  1622,   903,  3240,  2286,
    2286, -2827,  2286,   854,  2286,  2286,  2286,  2286, -2827,  2286,
    2286, -2827, -2827, -2827, -2827,   610,   832, -2827,  2286, -2827,
   -2827,  2286, -2827,  3247,  2286, -2827,  2802,  2286,  2286,  2286,
    3261, -2827, -2827,  2286,   309, -2827,  2286,   380,  3262, -2827,
    2286, -2827,  3269,   520,  2286, -2827,  3275, -2827,  1624,  3276,
     854,  2809, -2827, -2827,  3284, -2827,  1637, -2827,  2812, -2827,
   -2827,  3288, -2827,  3290, -2827,   854,  -189, -2827, -2827, -2827,
    3291, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827,  3292, -2827, -2827,  -155,   854, -2827,   854,
   -2827,   347, -2827, -2827, -2827, -2827, -2827,  3294,  3296, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827,  3298,  3300,  3306, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827,  3308, -2827, -2827,  3310, -2827, -2827, -2827, -2827,
   -2827, -2827,  3311, -2827, -2827,  3312, -2827, -2827, -2827,   854,
    2815, -2827, -2827,  1638, -2827,  2286, -2827,  3323,  2286, -2827,
     854,   854, -2827, -2827,   854,   854,  2286,  2286,  2286, -2827,
    2286,  2286, -2827,  2286,  2286,  2286,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286, -2827,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,   854,   854,   854,   854,
    2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286,
   -2827, -2827,   854,  2286,  2286, -2827, -2827,   854, -2827, -2827,
   -2827,  3324,   497,   854,  2286, -2827,  3331, -2827,  3458,  2286,
   -2827,  3335,   854,   854,   854, -2827,  1654, -2827,  3337,  2286,
    3342, -2827,  1659,  3343,  3344, -2827, -2827, -2827, -2827,  3345,
     611, -2827,  3346, -2827, -2827,   854, -2827, -2827, -2827,  2286,
   -2827,  2286, -2827,  3107, -2827,  3107,  3274,  2286,  2286,  2286,
    2286,  2286, -2827,  2286,  3350, -2827,  2286,  2286, -2827,  2286,
    2286, -2827,  2286,  2286, -2827,  3352,  1683,  3064, -2827, -2827,
    2286,  2286, -2827,  2286,  2286, -2827,  2286,  2286, -2827, -2827,
   -2827,  3354,  1696,  1697, -2827, -2827, -2827, -2827,  3356, -2827,
   -2827,  1701,  2286,  3358,   854,  2286,   854,  2286,  2286,  3360,
    2286,  2286, -2827,   854, -2827, -2827, -2827,   836, -2827,   842,
   -2827, -2827, -2827,   854, -2827,  1702, -2827,  3368, -2827,  3369,
   -2827,   854,  3370, -2827,  2286, -2827,  3372, -2827, -2827, -2827,
   -2827,  3373, -2827,  3377, -2827,  2286,  2286,  2286,   854,   855,
   -2827,  1716,  2286,  2286, -2827,  3382,  2286, -2827,   930, -2827,
    2286,  1721,   956,  3383,  2286,  2286,  2286,  3384,  2286,  2286,
   -2827, -2827,   760,   869,  3385,  3388, -2827,  2286, -2827,  2828,
    2286,  2286,  3390, -2827,  3391,  3392, -2827,  -201, -2827,  2286,
    2286,  3393, -2827,  3395, -2827, -2827,  3397,   537, -2827, -2827,
    2286, -2827,  2849, -2827,  2881, -2827, -2827,  2893, -2827,  2928,
   -2827, -2827,  -170, -2827,  3400, -2827, -2827,   854, -2827,  3407,
    3408,   854, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827,  3409, -2827,   854, -2827,  2948,  3417, -2827,   -13,
     854,  2286,   854,  2286,  2286,  2286,  2286,  2286,    14,  2286,
      30,  2286,  1726,  3418,  2286,  2286,  -186,  2286,  2286,  2286,
    2286,  2286,  2286,    33,  2286,  3424,  2286,  2286,  2286,  2286,
    2286,   854,   854,  2286,  2286,  2286,  2286,  2286,  2286,  2286,
    2286,  2286,  2286,  2286,  2286,  2286,  2286, -2827,  2286, -2827,
     854,  2286,   854,  3426, -2827,  3427,   854, -2827,   854, -2827,
    2957, -2827,  3428, -2827, -2827,  3432, -2827, -2827, -2827, -2827,
     844, -2827, -2827, -2827, -2827, -2827, -2827,  3434,  2286,  2286,
    2286,  3436,  2286,  3438, -2827,  3441,  2286,  3442,  2286,  3443,
    2286, -2827,  2286,  3446,  1734,  3449,  2286,  3452,  2286,  3457,
    2286, -2827, -2827,  3459, -2827,  3466, -2827, -2827,   854,   854,
   -2827,  2286,  3467,  2286,  2286, -2827,  2286,  1736,  3468,   854,
   -2827,   854, -2827,  3472, -2827,  3476, -2827, -2827,  2286, -2827,
    3478, -2827, -2827, -2827,  2286,  2286,  2286,  2286,  3480, -2827,
   -2827,  2286,  2286,  2286, -2827,  3486,  2286,  2286, -2827,   345,
    2286, -2827,  3487,  2286, -2827,   958, -2827,  2286,  2286,  2286,
   -2827,  2286,  2286, -2827, -2827,   864,   717, -2827, -2827,   854,
   -2827,  3488,  2286,  2286, -2827, -2827, -2827, -2827,  2286,  3492,
     388, -2827, -2827, -2827, -2827,  3494,  3496, -2827,  2974, -2827,
    3001, -2827,  3501, -2827,  3502, -2827,  3505, -2827,  3508, -2827,
   -2827,  3513, -2827,  3515, -2827,  3518, -2827, -2827,   854,  2286,
     866,  2286,  2286,  2286,  2286,  2286,  2286,  2286,  2286, -2827,
    2286,  2286,  2286, -2827,  2286,  2286,  1798, -2827,  1799, -2827,
    2286,  2286, -2827,  2286,  2286,  2286,  2286,  -181,  2286,  2286,
   -2827,  2286,  2286,  1857, -2827,  2286,  2286,  1859,  2286,  2286,
    2286,  2286,  2286,  2286,    47,    49,  -175,  2286,  2286,  2286,
    2286,  2286,  2286,  2286, -2827, -2827, -2827,   854,  3520,  2286,
   -2827,  3648, -2827, -2827,  3530, -2827, -2827, -2827, -2827,  2286,
    2286,  3532, -2827,  3538, -2827, -2827,  3544, -2827,  3545, -2827,
    3547,  1862, -2827, -2827,  2286, -2827,  3549, -2827,  3550, -2827,
    3551, -2827, -2827,   854,   854,  3552, -2827,  2286,  2286,  2286,
   -2827,  2286, -2827, -2827, -2827, -2827, -2827,  3554, -2827,  3555,
    3556,  3557,   884, -2827,  2286,  2286,   183, -2827,  2286,  2286,
    2286, -2827,  2286, -2827,  3558,  2286, -2827,   504,  2286,  2286,
    2286,  2286,  2286, -2827, -2827,   899,  2286,  2286, -2827,  2286,
    2286,   455, -2827,  3560, -2827, -2827, -2827, -2827,  3002, -2827,
    3013, -2827, -2827, -2827, -2827, -2827, -2827, -2827,  3561,   268,
   -2827,  1864,  2286,  1908,  2286,  1941,  1993,  2286,  3562,  2286,
    -151,  3565,  2286,  -134, -2827,  2286, -2827,  1998,  2286,  2286,
    2286,  2286,  2286,  2286, -2827,  2286,  2286,  3566,  2286,  -109,
   -2827,  2286,  2003,  2020, -2827,  2286,  2286,  2286,  2286,  2286,
    2286,  2286, -2827,  2286,  2286, -2827,  2286,  2286, -2827,  2286,
     -89,  2286,  2286,  2286,  2286,  2286,  2286,  2286, -2827,  3567,
    3568, -2827,  2286,  2286, -2827, -2827, -2827, -2827, -2827,  2286,
    2021,  2044, -2827, -2827, -2827, -2827,   854,  2286, -2827,   552,
    2286,  2286,  2052, -2827, -2827, -2827, -2827,  3570, -2827,  2286,
    2286, -2827,  2286,   282,  2286,  3578,  2286,   499, -2827,  2286,
    2286, -2827,  2286,  2286,  2286,  2286,  2082, -2827,   738,  2286,
    2286,  2286,  3580, -2827, -2827, -2827,  3044, -2827,  3056, -2827,
   -2827,   854, -2827,  2286,  2286, -2827,  2286,  2286, -2827,  2083,
   -2827,  2085,  2286, -2827,  2286, -2827,  2286, -2827,  2286, -2827,
    2286,  3582, -2827,  2106,  2286,  2286,  2286,  2286,  2207,  2214,
    2286,  2286, -2827,  2286, -2827,  2286,  2249, -2827,  3586, -2827,
    2301,  3589,  2286,  2286,  2286,  2286,  2286,  2286,  2286,   -25,
    2286,   -12,  2286, -2827,  2286,  2364,  2379,  2286,  2286,  2286,
    2286,  3590, -2827,  3716,  3594,  2286, -2827, -2827,  2286,  2388,
     854,  2286,   854,  3596,  2286,  2286, -2827,  2286, -2827,  2286,
    2286,  3597, -2827,  2286,  -177,  2286, -2827,  2286, -2827,  2286,
     516,  3598,  2286,   854,  2286,  2286,  2286, -2827,  2286, -2827,
     913,  2286,  2286,  2286, -2827, -2827,  3109, -2827,  3600,   765,
    3602,  2389,  3605,  2390, -2827,  2286, -2827,  2286,  2286,  3606,
    2286,  3608,  2286, -2827, -2827,  2396,  2286,  2286,  3610,  2286,
   -2827,  2416, -2827,  2417,  2286,  2286,  3611,  2286, -2827,  2440,
   -2827, -2827,  2286, -2827,  2286,  2286,  2286,  2286,  2286,  2286,
    2286, -2827,  2286,  2286, -2827,  2286,  2286,  2286, -2827,  2492,
   -2827,  2499,  2286,  2286,  2286,  2286, -2827,  3613, -2827,  3614,
    2565, -2827, -2827,  3619,  3620,   854,   854,  2286,  3621,  2566,
    2286,  2286, -2827,  3623, -2827,  2286, -2827,  3624,  3627, -2827,
    2286,   189, -2827,  2286,  3628,  2286,  2286,  2286,  2286, -2827,
    2286,  2286,  2570, -2827,  3631, -2827, -2827,   924, -2827, -2827,
    2286, -2827, -2827,  2286,  2571,  2572,  2574, -2827,  2286, -2827,
    2286, -2827,  3143,  2286,  2286, -2827,  2286, -2827,  3146, -2827,
    3202,  3632,  2286, -2827,  2286, -2827,  2575,  2576,  2286,  2286,
    2286,  2286,  2286,  2286,  3634,  2286,  3640,  2286,  3641,  2286,
   -2827,  2286, -2827,  2286,  2286,  2286,  2286,  2286, -2827, -2827,
   -2827,  2598, -2827, -2827,   854,  2701, -2827, -2827,  3642,  2720,
    3645, -2827,  3646, -2827, -2827,  3649, -2827,  2286,  3650, -2827,
    2286,  2286,  3652,  2286,  2286,  2778, -2827,  2286, -2827, -2827,
    3653,  3654, -2827,  2805, -2827,  2811, -2827,  2882,  3655,  3656,
   -2827,  3249,    57,    65,    19, -2827,  3657, -2827,  3658, -2827,
    2286,  3659, -2827,  2913, -2827,  3067,  2286,  2286,  2286,  2286,
    2286,  2286, -2827,  2286, -2827,  2286, -2827,  3660,  3073,  3110,
    2286,  2286,  2286,  2286, -2827, -2827,  2286, -2827, -2827,  3195,
   -2827, -2827, -2827,  3661, -2827,   854,  3662, -2827,  2286, -2827,
   -2827,  2286,  2286, -2827, -2827, -2827,  3277, -2827,  3285, -2827,
    3208, -2827, -2827, -2827,   854, -2827,  2286,  2286, -2827,  2286,
    2286, -2827,  2286, -2827, -2827,    22, -2827, -2827,  3279, -2827,
    3314,  3281,  3297,  2286,  2286,  2286,  2286,  3664,  3665, -2827,
   -2827,  3299, -2827,  3338,   141,   144,    62,  2286,  3371, -2827,
    3666, -2827,   854, -2827,   782,  2286,  3667, -2827,  3668, -2827,
    3669, -2827,  3670,   854,  2286,    63,  2286,   119,  2286, -2827,
    2286, -2827,  3415, -2827,  3672, -2827,  2286, -2827,  2286,  3430,
    3462,  2286,  2286, -2827, -2827, -2827,  3327, -2827,  3341, -2827,
    2286,  2286, -2827,  2286,  2286, -2827,  2286,   154, -2827,  3673,
   -2827,   854,  2286, -2827,  2286,  3674, -2827, -2827, -2827, -2827,
     854,  2286, -2827,  2286,  2286, -2827,  2286,  2286,  2286, -2827,
    3553, -2827,  3677,  3678, -2827,  2286, -2827,  2286,  2286,  2286,
   -2827,  3682, -2827,  3683,  2286,   188,  2286,   212,  2286, -2827,
    2286, -2827,  3684,  3685,   852, -2827,  2286,  3686,  2286,  3687,
    2286,  3688,  2286, -2827,  3559, -2827, -2827,  3689,  3690,  3577,
    3584, -2827, -2827,  2286, -2827,  2286,  2286, -2827,  2286,  2286,
    2286, -2827, -2827,  2286, -2827,  2286,  2286, -2827,  2286, -2827,
    2286, -2827,  3691, -2827,  3593, -2827, -2827, -2827,  2286, -2827,
    2286,  3692,  2286,  3693,  2286,  3694,  2286,  2286,  2286,  2286,
    3695,  3696, -2827, -2827,  3637,  3697,  3698, -2827,  2286, -2827,
    2286, -2827,  3699,  3700,  2286,   854, -2827, -2827, -2827,  3643,
   -2827, -2827,  3701,  3702, -2827, -2827,   896,  2286, -2827,  3644,
   -2827, -2827,  2286, -2827,  2286,  3703, -2827,  3647,  2286,   854,
   -2827, -2827,  3663,  3704,  2286, -2827,  3671, -2827,  2286, -2827,
    3705,   854, -2827,   854,  2286,  2286,  2286,  3706, -2827
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_uint16 yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     3,   127,   128,   129,   130,   131,   116,
      79,   122,   123,   124,    55,    56,    57,    52,    53,    54,
      51,    69,    74,    75,    76,    38,    39,    41,    40,    49,
      42,    43,    44,    46,    88,    87,    95,   322,    45,    48,
      81,    82,    83,    84,   126,    85,    86,    67,    68,    73,
      70,    71,    72,   143,   113,   114,   115,   112,     6,   141,
      92,    93,    89,    90,    91,    94,   105,   106,   101,   107,
     108,   109,   110,   117,   118,   119,   120,   121,   102,   103,
      31,    30,    34,    32,    33,    35,    36,    37,     7,     8,
      58,    59,    60,    61,    62,    64,    63,    65,    66,    10,
       9,    11,   111,    22,    12,   134,   135,   137,   136,   138,
      47,   139,   140,   144,     5,    15,    13,    14,   142,    16,
      17,    18,    20,    21,    19,    25,    26,    27,    28,    78,
      23,    24,    96,   145,    98,   104,    99,   100,    97,    50,
      80,    77,   125,   132,   133,    29,   146,   147,   148,   149,
     151,   150,   152,     0,     0,     0,     0,  1425,  1426,     0,
     836,     0,  1428,  1430,  1427,  1429,     0,     0,     0,     0,
     334,     0,     0,     0,     0,     0,   834,   558,     0,   234,
     488,     0,   228,   356,  1138,   761,     0,   468,   822,     0,
       0,  1104,   260,   463,   361,   194,     0,  1071,  1076,     0,
       0,     0,   854,     0,   323,     0,   824,     0,     0,     0,
     332,     0,     0,     0,   484,     0,   570,     0,   209,     0,
     752,   557,   264,   420,     0,   155,     0,     0,     0,     0,
     217,     0,  1082,     0,     0,     0,     0,     0,   415,   434,
     653,   548,     0,   554,     0,     0,   563,     0,     0,     0,
       0,     0,     0,     0,     0,   254,   639,     0,     0,   220,
       0,   341,   859,   885,     0,     0,   353,     0,     0,   214,
       0,   427,     0,     0,  1107,     0,     0,     0,     0,   828,
     893,     0,     0,   280,     0,     0,     0,   281,     0,   391,
       0,     0,     0,     0,   888,   872,     0,     0,     0,     0,
     776,   492,     0,   429,     0,     0,     0,     0,  1137,   266,
     159,   634,   629,     0,     0,     0,   920,   327,     0,     0,
     479,     0,     0,   357,     0,     0,   412,   411,   597,   741,
     343,     0,   275,     0,   592,   602,   607,   624,   614,   619,
     183,  1141,     0,   413,     0,   161,     0,  1149,  1305,     0,
       0,     0,   207,   351,     0,   416,   435,     0,     0,     0,
     758,  1315,  1420,  1355,  1360,     0,   869,   864,   866,  1351,
    1353,     0,     1,     2,     4,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   172,   173,   174,
     170,   171,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   229,   230,   231,   232,   233,   235,     0,
     255,     0,     0,     0,     0,     0,     0,   276,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   362,   363,   364,     0,     0,
       0,   392,   393,   408,     0,     0,     0,     0,     0,     0,
     460,     0,     0,   464,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   502,     0,   513,     0,   516,     0,   519,
       0,   522,     0,   531,   533,   535,     0,   546,     0,     0,
       0,     0,     0,     0,     0,     0,   572,     0,     0,   668,
       0,   724,     0,     0,     0,     0,   756,     0,     0,   767,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   805,     0,
     823,   825,   829,     0,     0,     0,     0,     0,     0,     0,
       0,   860,   861,     0,  1427,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   961,   921,  1091,  1090,  1089,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1133,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1310,  1310,
    1310,  1310,  1310,     0,     0,     0,     0,     0,  1310,     0,
       0,     0,  1341,     0,     0,  1350,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1357,
    1358,  1359,     0,     0,     0,     0,   273,   582,     0,   399,
       0,     0,   193,   837,   530,   532,     0,   335,     0,     0,
     525,   527,     0,   528,     0,   344,  1083,     0,   489,     0,
       0,     0,     0,     0,     0,  1079,  1072,     0,  1077,     0,
       0,  1069,   855,   324,   512,   501,   569,   590,     0,  1067,
       0,   536,     0,     0,   485,     0,   571,   339,     0,     0,
     455,   665,     0,   663,     0,   495,   496,     0,   218,   515,
    1100,     0,  1102,     0,   521,   518,   654,     0,     0,   550,
     549,   553,   567,   564,     0,     0,     0,   459,     0,     0,
     340,     0,     0,   640,     0,     0,     0,   270,     0,   221,
       0,   886,   354,   887,   667,   215,   426,   326,   803,     0,
       0,   329,   288,   284,     0,   282,   289,   285,     0,   283,
       0,   559,   561,     0,   873,   345,     0,     0,     0,   493,
     428,     0,   543,     0,   545,     0,     0,   635,     0,   630,
     272,     0,   333,     0,   506,   507,     0,   330,   331,   723,
       0,     0,   598,   271,   274,     0,   593,     0,   603,     0,
     608,     0,   625,     0,   615,     0,   620,     0,   414,   162,
     725,     0,     0,     0,   740,     0,   352,   470,   471,     0,
     755,   475,   476,     0,     0,   325,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   177,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   205,   201,
       0,   204,   202,     0,   206,     0,   199,   200,   198,   197,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   243,     0,   245,   247,     0,
     249,     0,     0,     0,     0,     0,     0,     0,     0,   279,
       0,     0,     0,     0,     0,     0,     0,   319,   320,   318,
     321,   310,   307,   308,   317,   314,   290,     0,   316,   311,
       0,   304,   305,   312,     0,     0,   348,   350,   349,     0,
     388,   377,     0,     0,     0,   390,     0,   365,   358,     0,
       0,   382,   372,     0,     0,     0,     0,     0,     0,     0,
     366,   359,     0,   394,     0,     0,   404,     0,     0,   407,
     409,   432,   431,   433,     0,   467,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   577,     0,   643,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   763,
       0,     0,     0,   792,   791,     0,     0,     0,     0,     0,
       0,   788,     0,   789,   790,     0,   781,   783,   778,   779,
       0,   793,   786,     0,   780,     0,   787,     0,     0,     0,
       0,     0,     0,     0,     0,   858,   857,     0,     0,     0,
     865,     0,   868,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   918,     0,   916,
       0,     0,   958,     0,     0,   931,   933,     0,     0,     0,
     945,     0,     0,     0,   953,     0,   939,     0,     0,   923,
       0,     0,   935,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1109,     0,     0,     0,  1128,     0,     0,  1105,
       0,  1106,     0,  1108,  1110,     0,     0,     0,  1121,  1111,
       0,     0,     0,     0,     0,  1308,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1333,  1342,  1348,
    1347,  1340,  1329,  1331,  1344,  1317,  1318,  1319,  1320,  1321,
    1334,  1322,  1326,  1325,  1324,  1323,  1327,  1332,  1412,     0,
       0,     0,  1316,     0,     0,  1380,  1392,  1391,  1378,  1409,
    1379,  1384,  1386,  1385,  1387,  1388,  1369,  1370,  1389,  1390,
    1362,  1365,  1371,  1372,  1368,  1400,  1401,     0,  1399,  1381,
    1402,  1403,  1393,  1396,  1404,  1405,  1375,  1376,  1377,  1406,
       0,     0,  1352,  1354,  1415,  1418,  1356,     0,     0,     0,
    1361,  1422,  1424,  1421,   583,     0,     0,   400,     0,     0,
       0,   336,   337,     0,   526,     0,   529,  1084,     0,     0,
       0,   762,   380,     0,   347,  1073,  1078,  1070,  1080,   591,
       0,  1068,   537,   538,  1097,     0,     0,     0,     0,   419,
       0,     0,   666,     0,   664,   497,     0,  1101,  1103,     0,
     656,     0,   655,     0,     0,   660,     0,  1088,  1092,     0,
       0,     0,     0,   154,     0,     0,     0,   645,     0,   644,
       0,   647,     0,     0,     0,     0,     0,     0,   286,   287,
       0,     0,   346,  1086,  1087,   430,   544,  1312,     0,     0,
       0,     0,     0,     0,   960,     0,   508,     0,   417,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1308,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1142,   726,     0,   730,     0,     0,     0,   739,   472,   474,
       0,   477,     0,     0,     0,     0,     0,     0,   179,     0,
       0,     0,     0,     0,   176,     0,     0,   163,     0,     0,
     184,     0,   186,   188,   189,   190,   191,   192,   195,   203,
     196,   208,   210,   213,   212,   211,   216,   219,     0,     0,
     225,     0,     0,     0,   242,   244,   246,   248,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   268,   269,
       0,     0,   815,     0,   291,     0,   297,     0,   309,   315,
     303,   313,   306,   342,   384,   386,     0,     0,   385,   370,
     383,   454,   588,     0,   379,   368,   367,     0,   373,     0,
       0,   371,   360,     0,   395,     0,   405,   406,     0,     0,
       0,     0,     0,   480,     0,   486,     0,     0,     0,     0,
     504,     0,     0,     0,     0,     0,     0,     0,   551,     0,
       0,     0,     0,   565,     0,   568,     0,     0,     0,     0,
       0,   684,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   812,     0,     0,     0,     0,     0,     0,   764,
       0,   768,   777,   797,   798,   799,   800,   801,   794,   802,
     784,   785,   782,   795,   796,   808,     0,     0,     0,     0,
       0,     0,   840,     0,     0,   856,     0,   863,   867,   870,
       0,     0,     0,     0,     0,     0,   874,     0,     0,   889,
       0,   899,   900,   894,   895,   896,   898,   901,   919,   902,
     917,   903,     0,   897,   925,   922,   932,   934,   929,   950,
     952,   951,     0,   946,     0,     0,   940,     0,   930,     0,
     959,   936,   937,   938,   926,     0,     0,     0,     0,  1034,
       0,  1033,     0,     0,     0,  1035,     0,     0,   987,     0,
       0,  1009,  1011,     0,     0,     0,     0,     0,  1047,     0,
       0,  1032,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1008,     0,     0,     0,     0,     0,  1040,     0,  1000,
       0,     0,     0,     0,     0,  1030,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1048,
       0,     0,     0,     0,     0,     0,     0,     0,   962,     0,
       0,  1012,     0,     0,  1037,     0,     0,     0,     0,  1119,
       0,  1129,  1122,  1123,     0,  1113,  1114,     0,  1134,     0,
    1132,  1112,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1306,     0,  1343,
    1349,  1330,  1345,  1311,  1328,     0,  1337,  1339,     0,  1413,
    1364,  1363,  1366,  1373,     0,  1411,  1382,  1394,  1397,  1407,
    1408,  1417,     0,  1419,  1423,     0,     0,   403,     0,     0,
     534,   338,   524,  1085,   355,     0,     0,  1074,  1081,   819,
    1098,     0,   578,     0,     0,     0,   418,     0,   424,     0,
       0,     0,   657,     0,     0,   659,  1093,     0,   457,     0,
       0,     0,     0,   642,   646,     0,   648,     0,     0,   641,
     222,     0,   223,     0,   153,   441,   437,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   328,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1308,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   727,     0,     0,   731,     0,   732,     0,   473,     0,
    1414,     0,   156,     0,     0,   175,     0,     0,     0,     0,
       0,   178,   181,     0,   180,   185,   187,   227,   226,   236,
       0,     0,     0,     0,     0,   827,     0,     0,   263,     0,
     261,   265,   267,     0,     0,   817,     0,   814,     0,   293,
       0,     0,   299,     0,   389,   378,     0,   452,   589,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   483,
       0,     0,   487,     0,     0,     0,   503,   505,   514,   517,
     520,   523,     0,     0,     0,   552,     0,     0,     0,   566,
       0,     0,     0,   721,     0,     0,     0,     0,     0,     0,
       0,   685,     0,     0,     0,     0,     0,     0,   718,     0,
       0,   811,   813,   729,   742,     0,     0,   753,     0,   757,
     759,     0,   766,     0,   769,   807,     0,     0,     0,     0,
       0,   838,   850,     0,     0,   851,     0,     0,     0,   871,
       0,   876,     0,     0,     0,   875,     0,   878,     0,     0,
       0,     0,   949,   947,     0,   954,     0,   941,     0,   924,
     927,     0,  1016,     0,  1022,     0,     0,  1064,  1015,  1013,
       0,  1027,  1065,  1066,  1007,  1006,  1010,  1025,  1026,   981,
     984,  1046,  1045,     0,   985,   986,     0,     0,  1054,     0,
    1053,     0,  1052,  1051,  1005,  1004,  1061,     0,     0,   979,
     980,  1041,  1042,   988,   989,   990,  1031,   964,  1001,   992,
     991,  1039,  1062,     0,     0,     0,  1028,  1003,  1002,   969,
     972,   970,   971,   973,   974,   975,   976,   967,   968,   966,
    1044,   995,     0,   983,   993,     0,   982,   965,  1043,  1029,
     963,  1036,     0,  1049,  1023,     0,  1038,  1096,  1117,     0,
       0,  1120,  1124,     0,  1115,     0,  1135,     0,     0,  1303,
       0,     0,  1304,  1309,     0,     0,     0,     0,     0,  1174,
       0,     0,  1180,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1186,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1307,  1346,     0,     0,     0,  1367,  1374,     0,  1383,  1395,
    1398,     0,     0,     0,     0,   401,     0,   410,     0,     0,
    1075,     0,   580,   579,     0,   421,     0,   821,     0,     0,
       0,   498,     0,     0,     0,   456,   458,   462,  1095,     0,
       0,   650,     0,   804,   442,     0,   444,   445,   446,     0,
     448,   449,   451,     0,   438,     0,  1313,     0,     0,     0,
       0,     0,   632,     0,     0,   509,     0,     0,   600,     0,
       0,   595,     0,     0,   605,     0,     0,     0,  1308,   611,
       0,     0,   627,     0,     0,   617,     0,     0,   622,   728,
     733,     0,     0,     0,   478,   158,   157,   160,     0,   164,
     165,     0,     0,     0,     0,     0,   238,     0,     0,     0,
       0,     0,   262,     0,   277,   816,   294,   292,   300,   298,
     387,   453,   818,     0,   374,     0,   436,     0,   396,     0,
     461,     0,     0,   469,     0,   481,     0,   490,   494,   547,
     539,     0,   540,     0,   556,     0,     0,     0,     0,     0,
     680,     0,     0,     0,   687,     0,     0,   706,     0,   692,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     747,   743,     0,     0,     0,     0,   765,   770,   806,     0,
       0,     0,     0,   839,     0,     0,   841,     0,   844,     0,
       0,     0,   862,     0,   880,   881,     0,     0,   879,   890,
       0,   891,     0,   904,     0,   948,   955,     0,   942,     0,
     928,  1017,     0,  1018,     0,  1014,  1063,     0,  1058,     0,
       0,     0,  1057,   977,   978,   997,   999,   998,   996,   994,
    1050,  1024,     0,  1130,     0,  1125,     0,     0,  1136,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1335,     0,  1416,
       0,     0,     0,     0,   402,     0,     0,  1099,   581,   423,
       0,   820,     0,   499,   662,     0,   658,  1094,   652,   649,
       0,   224,   443,   447,   450,   440,   439,     0,     0,     0,
       0,     0,     0,     0,   510,     0,     0,     0,     0,     0,
       0,  1308,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   734,   735,     0,   737,     0,   166,   169,     0,     0,
     182,     0,     0,     0,     0,   253,     0,     0,     0,     0,
     295,     0,   301,     0,   375,     0,   398,   397,     0,   465,
       0,   491,   541,   542,     0,     0,     0,     0,     0,   573,
     681,     0,     0,     0,   689,     0,     0,     0,   710,     0,
       0,   688,     0,     0,   708,     0,   693,     0,     0,     0,
     719,     0,     0,   748,   744,     0,     0,   754,   760,   771,
     810,     0,     0,     0,   835,   842,   843,   847,     0,     0,
       0,   852,   877,   883,   882,     0,     0,   905,     0,   906,
       0,   956,     0,   943,     0,  1019,     0,  1020,     0,  1056,
    1055,     0,  1118,     0,  1126,     0,  1116,  1143,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1173,
       0,     0,     0,  1179,     0,     0,     0,  1295,     0,  1226,
       0,     0,  1255,     0,     0,     0,     0,     0,     0,     0,
    1185,     0,     0,     0,  1284,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1336,  1338,  1410,     0,     0,     0,
     584,     0,   381,   422,     0,   500,   661,   651,  1314,     0,
       0,     0,   631,     0,   511,   599,     0,   594,     0,   604,
       0,     0,  1308,  1308,     0,   626,     0,   616,     0,   621,
       0,   736,   738,     0,     0,     0,   240,     0,     0,     0,
     256,     0,   278,   296,   302,   369,   376,     0,   482,     0,
       0,     0,     0,   575,     0,     0,     0,   707,     0,     0,
       0,   714,     0,   690,     0,     0,   712,     0,     0,     0,
       0,     0,     0,   749,   745,     0,     0,     0,   809,     0,
       0,     0,   853,     0,   845,   884,   892,   907,     0,   908,
       0,   957,   944,  1021,  1060,  1059,  1131,  1127,     0,     0,
    1146,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1190,     0,  1296,     0,     0,     0,
       0,     0,     0,     0,  1261,     0,     0,     0,     0,     0,
    1267,     0,     0,     0,  1192,     0,     0,     0,     0,     0,
       0,     0,  1202,     0,     0,  1206,     0,     0,  1210,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   585,     0,
       0,   425,     0,     0,   637,   633,   601,   596,   606,     0,
       0,     0,  1308,   628,   618,   623,     0,     0,   241,     0,
       0,     0,     0,   466,   555,   560,   562,     0,   574,     0,
       0,   694,     0,     0,     0,     0,     0,     0,   709,     0,
       0,   716,     0,     0,     0,     0,     0,   750,     0,   772,
       0,     0,     0,   848,   846,   909,     0,   910,     0,  1145,
    1144,     0,  1194,     0,     0,  1196,     0,     0,  1150,     0,
    1156,     0,     0,  1177,     0,  1172,     0,  1183,     0,  1178,
       0,     0,  1297,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1189,     0,  1184,     0,     0,  1285,     0,  1166,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1212,     0,     0,     0,     0,     0,     0,
       0,     0,   586,  1139,     0,     0,  1308,   609,     0,     0,
       0,     0,     0,     0,     0,     0,   257,     0,   576,     0,
       0,     0,   695,     0,     0,     0,   711,     0,   697,     0,
       0,     0,     0,     0,     0,     0,     0,   673,     0,   746,
       0,   773,     0,     0,   849,   911,     0,   912,     0,     0,
       0,     0,     0,     0,  1151,     0,  1157,     0,     0,     0,
       0,     0,     0,  1191,  1298,     0,     0,     0,     0,     0,
    1287,     0,  1291,     0,     0,     0,     0,     0,  1268,     0,
    1286,  1167,     0,  1193,     0,     0,     0,     0,     0,     0,
       0,  1203,     0,     0,  1207,     0,     0,     0,  1214,     0,
    1220,     0,     0,     0,     0,     0,   587,     0,   636,     0,
       0,  1308,   612,     0,     0,   237,     0,     0,     0,     0,
       0,     0,   700,     0,   696,     0,   722,     0,     0,   698,
       0,     0,   713,     0,     0,     0,     0,     0,     0,   751,
       0,     0,     0,   913,     0,   914,  1147,     0,  1195,  1198,
       0,  1197,  1200,     0,     0,     0,     0,  1176,     0,  1182,
       0,  1299,     0,     0,     0,  1256,     0,  1288,     0,  1292,
       0,     0,     0,  1188,     0,  1269,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1215,     0,  1221,     0,     0,     0,     0,     0,  1140,   638,
     610,     0,   167,   168,   239,     0,   826,   258,     0,     0,
       0,   701,     0,   715,   703,     0,   699,     0,     0,   691,
       0,     0,     0,     0,   774,     0,   830,     0,   915,  1148,
       0,     0,  1152,     0,  1158,     0,  1162,     0,     0,     0,
    1300,     0,     0,     0,     0,  1289,     0,  1293,     0,  1262,
       0,     0,  1270,     0,  1168,     0,     0,     0,     0,     0,
       0,     0,  1204,     0,  1208,     0,  1211,     0,     0,     0,
       0,     0,     0,     0,   613,   250,     0,   259,   677,     0,
     720,   702,   704,     0,   717,     0,     0,   683,     0,   775,
     832,     0,     0,  1199,  1201,  1153,     0,  1159,     0,  1163,
       0,  1175,  1181,  1301,     0,  1227,     0,     0,  1247,     0,
       0,  1257,     0,  1290,  1294,     0,  1187,  1271,     0,  1169,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1213,
    1216,     0,  1222,     0,     0,     0,     0,     0,     0,   678,
       0,   705,     0,   686,     0,     0,     0,  1154,     0,  1160,
       0,  1164,     0,     0,     0,     0,     0,     0,     0,  1263,
       0,  1272,     0,  1170,     0,  1231,     0,  1233,     0,     0,
       0,     0,     0,  1205,  1209,  1217,     0,  1223,     0,  1243,
       0,     0,  1251,     0,     0,  1259,     0,     0,   251,     0,
     679,     0,     0,   674,     0,     0,   831,  1155,  1161,  1165,
       0,     0,  1228,     0,     0,  1248,     0,     0,     0,  1273,
       0,  1171,     0,     0,  1235,     0,  1237,     0,     0,     0,
    1218,     0,  1224,     0,     0,     0,     0,     0,     0,  1265,
       0,   252,     0,     0,     0,   833,     0,     0,     0,     0,
       0,     0,     0,  1274,     0,  1232,  1234,     0,     0,     0,
       0,  1219,  1225,     0,  1244,     0,     0,  1252,     0,     0,
       0,   682,   675,     0,   669,     0,     0,  1229,     0,  1249,
       0,  1258,     0,  1275,     0,  1236,  1238,  1239,     0,  1241,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1264,  1276,     0,     0,     0,  1245,     0,  1253,
       0,  1260,     0,     0,     0,     0,  1230,  1250,  1277,     0,
    1240,  1242,     0,     0,  1266,   670,     0,     0,  1278,     0,
    1246,  1254,     0,   671,     0,     0,  1279,     0,     0,     0,
    1302,  1280,     0,     0,     0,  1281,     0,   672,     0,  1282,
       0,     0,  1283,     0,     0,     0,     0,     0,   676
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -2827, -2827, -2827,  3784, -2827, -2827, -2827, -2827, -2827, -2827,
    3675, -2827, -2827,  3676, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2175, -2827, -2827, -2827, -2827,
    3363,  3362, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,  3601, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827,  3604, -2827,  3339, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2079, -2827,  3679, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827,  2918, -2827, -2827,  2909, -2827, -2827, -2827, -2827,
   -2827, -2827,  3036, -2827,  -135, -1169, -2827, -2827,  2953, -2827,
    3595, -2827,   307, -2827,  3569, -2827,   -22,  -884,  -325, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827,  3497, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827,   502, -1196,  3571,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827,  3011,  -952, -2827, -2827,  3029,  -931, -2827,  -378, -2827,
    3539, -2826, -2827, -2827, -2827, -2827, -2827, -2827, -2827,  3453,
   -2827, -2827, -2827, -2827, -2827,  -440,  3490,  2761,  -187, -1731,
    -751,  -916, -2827, -2827,   620, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827,   191, -2827,  2971, -2827,
   -2827, -2827,  -335, -2827,   201,  -305, -2827, -2827,  2424, -2827,
   -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827, -2827, -1636,  2166, -2827, -2827, -2827, -2827, -2827, -2827,
   -2827,  1018,  3680, -2827, -2827, -2827, -2827, -2827,  3780,    -6
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   171,   172,   173,   174,   175,   176,   177,   178,   557,
     558,   559,   560,   561,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   593,  2179,   594,   595,   596,   597,
    1105,  1108,   189,   600,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   199,   615,   200,   201,   202,   617,   203,
     204,   205,   206,   634,   207,   208,  1483,   635,  1149,  1154,
     209,   641,   642,   210,   647,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   649,   221,   222,   636,   223,
    2116,  2513,   637,   224,   225,   226,   650,   227,   228,   229,
     230,  1047,  1048,   231,  1051,  1052,   232,   233,   234,   235,
     236,   935,   936,   237,   663,  1769,   238,  1014,  1015,   239,
     665,   240,   667,   241,   669,   242,   671,   890,   891,   243,
     244,   245,   246,   247,   248,   249,   677,   250,   251,   252,
     253,   254,   255,   256,   257,   258,   259,   876,  1741,   916,
     260,  1026,   261,  1022,   262,  1028,   263,  1030,   264,  1034,
     265,  1036,   266,  1032,   267,  1009,   268,  1007,   269,   270,
     963,   964,  1597,   271,   946,   947,  1580,   272,   930,   273,
     689,  2845,   274,   275,   276,   277,   278,   279,   280,   696,
     281,   282,   700,   283,   284,   728,   691,  1801,   877,  1742,
     917,   931,   285,   286,   598,   287,   732,   288,   289,   290,
     291,   741,   742,   292,   293,   294,   295,   296,   297,   298,
     299,  1861,  1288,  1286,   300,  1294,   774,   301,   775,   302,
     919,   303,   371,   304,  1587,  1588,   305,  1563,  1564,   306,
     940,   307,   942,   308,  1400,   309,   795,   310,   311,   312,
     313,   314,  1996,  1465,   315,   316,   824,   317,   318,   319,
     320,   864,   825,   321,   870,   871,   322,   875,  1743,  2846
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     336,   337,   338,   341,  2141,   978,  2566,  1544,  1778,  1544,
    2208,   607,  1599,  1783,  1572,  1785,  1582,  1572,   893,   933,
     569,  3029,   780,  1886,  1299,   366,   332,  3034,  1771,  1772,
    1773,  1774,   739,   366,   327,   909,   643,   366,  2514,   910,
     366,   653,   995,   394,   327,  1202,   332,  2208,  2877,   332,
    2154,  2156,  2208,   327,  2208,   644,   581,   435,  1806,   414,
    2673,   327,   333,  2932,   327,   332,  1653,  2280,  3094,   367,
     427,   686,  3374,  1473,  3118,   547,  2677,   367,   570,  2895,
     571,   367,   333,   446,   367,   333,  1215,   368,   369,   676,
     332,  1216,   323,   781,  2678,   368,   369,   328,  3205,   368,
     369,   333,   368,   369,   324,   326,   332,   328,  1789,   332,
     673,   782,   485,  1809,   389,  3209,   328,   491,   492,   562,
     563,   327,  1399,   332,   328,   332,   333,   328,   327,   327,
     327,   334,   572,   332,   783,   548,   549,   550,   551,   552,
    3224,   332,   333,   335,  1217,   333,   702,   524,   873,   874,
     968,   334,   332,   342,   334,   619,   703,   784,   535,   333,
    3243,   333,   411,   335,   327,  1559,   335,   582,   704,   333,
     334,   705,   706,   707,   708,   709,   327,   333,   585,   586,
     620,   436,   335,   346,   328,   604,   347,  2249,   333,  1300,
    1218,   328,   328,   328,   583,   334,   618,   428,   332,  2419,
    2138,  1559,  3375,  3165,  3166,  2281,  2878,   335,   621,  3169,
    1790,   334,  1131,   785,   334,   710,   332,   332,  1132,  2422,
     332,   786,  2436,   335,  3341,   675,   335,   328,   334,   787,
     334,   332,  2285,  1887,   333,   355,  2907,  3344,   334,   328,
     335,   621,   335,   915,  1791,   332,   334,   969,  1229,   910,
     335,  2106,   333,   333,   332,  1164,   333,   334,   335,   332,
     788,  2674,   789,  2919,   688,  1133,   327,   333,  3561,   335,
     790,  3609,  2933,   584,   332,  1203,   366,  3095,   511,  2923,
    2896,   333,  2940,  3119,   711,   699,   743,   746,   712,  2250,
     333,  1219,  1121,   622,   327,   333,  3112,   996,  3115,   779,
    2155,  2157,   553,   334,  1810,   554,  3555,  3206,  2315,   623,
     333,  3635,  3652,   797,  3558,   335,  1654,  1655,  1656,  1474,
     367,   334,   334,  2251,  3210,   334,   555,   606,   791,   328,
    2537,   800,   886,   335,   335,   889,   334,   335,   368,   369,
    3277,  2139,   897,  3281,  3282,   412,   325,   899,   335,  3225,
     334,  2208,  2257,   903,   618,  1263,   624,   328,   332,   334,
     905,   906,   335,   908,   334,  1231,   349,   645,  3655,  3244,
    2082,   335,   625,   626,  1792,   920,   335,   801,   802,   334,
     761,   627,   970,   713,   944,  2286,   327,   928,   929,  1301,
    3629,   335,   714,  3632,   333,   715,   628,  1111,   941,   943,
     716,   717,   718,  3679,   719,   573,   626,   720,  2420,   629,
     327,   721,   955,   327,   574,   958,  1134,   740,   961,  1119,
    2645,   966,   332,   764,   792,  2110,   356,  2908,  2423,  2909,
     793,  2437,  3161,  3342,  2775,  1793,  2776,  3704,  3466,  1230,
     979,   646,  2107,  1040,  1135,   352,  3345,  1136,   556,   328,
    1888,  3377,  2151,  2252,  1137,  1220,  3383,   722,   333,   327,
     766,  3707,   475,   334,   353,  1185,   794,  2258,   553,   512,
    1138,  2538,  2920,   328,   327,   335,   328,  3562,   437,  3030,
    3610,   723,  1013,   751,   752,  1221,   724,   370,  2924,  1020,
    2283,  2941,   555,  2159,  1169,   575,   576,   577,  1212,  3053,
     578,   579,   767,  2536,  2208,  3113,   630,  3116,  2208,   564,
     565,   566,   567,   568,  2253,  3556,   631,  3190,   632,   327,
    3636,  3653,   328,  3559,   354,   768,  1041,   334,   366,  1054,
    1302,  3272,  1139,  1140,  1189,  2152,   609,   328,   753,   335,
     327,  1059,  1060,  1061,  1062,  1063,  1064,  1065,  1067,  1846,
    1069,   910,  1278,   674,  1070, -1310,  1574,  3468,  2646,  1619,
    2111,   327,  3162,  1141,  1142,   633,  3182,  1841,  3467,  1546,
     327,  1042,   367,   736,  1090,   332,   327,  3656,  2681,  1095,
    1096,  1099,   328,  2235,   610,  1831,   476,   754,  2239,   327,
     368,   369,   332,  1112,  3031,  1113,  2682,   761,  1116,  3630,
     327,   327,  3633,   328,   357,   327,   327,  1679,   910,  1144,
     652,   333,  3680,  1150,  1152,  1153,  1155,  1186,   358,  1150,
     327,   762,   725,   332,   328,   911,  1167,   763,   333,  2648,
     332,   327,  1172,   328,  1175,   726,   727,  3054,  3170,   328,
     764,  1213,   359,   327,  2104,  1184,  3705,   765,  2092,   366,
    1681,  1190,   328,  1192,  1143,  1194,  1195,  1832,  1572,   333,
     755,  3273,  1544,   328,   328,   360,   333,   756,   328,   328,
    3708,   327,  1847,   327,  1206,  1207,  1208,   766,   327, -1310,
     334,  1717,  1557,   328,  1222,   757,   771,   772,  1170,  1534,
    1536,   803,   335,   367,   328,   361,  2392,   334,   327,  1236,
    1237,  1238,  1239,  1240,  3183,  2065,   328,   910,  3191,   335,
     939,   368,   369,  1250,   327,   362,   656,  1934,  1255,   767,
    1257,   363,  2043,  1259,  1639,   364,   327,  1261,   334,  1974,
    2291,  1680,  1267,  1931,   328,   334,   328,   658,  2295,  2264,
     335,   328,   768,  1277,  1279,  1280,  1281,   335,  3278,  1284,
    1285,  2208,  2208,  3171,  1291,   366,  1293,   611,   612,   613,
     614,   328,  1534,   365,  1310,  3379,  1534,   660,   332,  2655,
    2193,  2393,   376,  1391,  1682,  1393,   366,   328,  1397,  1398,
    1584,  2501,  2649,   737,  1405,   910,  2884,  2650,   804,   328,
     805,   806,   769,   327,   332,  2216,  1701,   685,   807,   367,
    2240,  3262,   748,   366,   333,  1718,   451,  1470,   377,  2198,
     808,   809,   810,   811,   812,  1833,   813,   368,   369,   903,
     367,  1117,  1486,  1534,  1932,   378,   814,   815,  1056,   816,
     333,   817,   818,   327,  1584,  1499,   819,  1534,   368,   369,
    1504,  1935,   820,   821,   822,   823,  2063,   367,   332,   758,
     759,   760,  1975,  2292,   379,  1899,   328,  3642,   332,  2630,
    2769,  2296,   327,   903,  1528,   368,   369,  2473,   638,   770,
     366,  1919,  2750,   334,  1538,  1539,   452,   327,  3279, -1425,
    1540,  1123,  1794,  1946,   333,   335,  1545,  2482,  2483,  1312,
     380,  1548,  2656,   327,   333,  3380,   328,  1553,   327,   744,
    1555,  1556,  2794,   771,   772,   730,   731,   366,   332,  2885,
     381,   335,  1560,   949,   367, -1425,   388,  1565,   953,  1702,
   -1425,  1125,  1570,  1571,  3263,   328,   327,  3713,   332,  2220,
     453,   639,   368,   369,  1721,   910,  2196,   327,  1313,   390,
     328,  1722,  1586,  2197,   333,   332,   332,   640,  1589,  1590,
    1591,   367,  1592,   334,   332,  1594,   328,   541,   773,  1602,
    1603,   328,   328,   334,   333,   335,  3044,  1585,   332,   368,
     369,  3772,   332,  1606,   332,   335,  2339,   327,  2201,   332,
     588,   333,   333,  2312, -1425,  1010,   391,  3289,  1892,   328,
     333,  1589,  1589,  2819,  2340,   339, -1425,   392,  1872,  2821,
     328,  1620,   327,  1623,   333,  1024,   332,  1625,   333,  2863,
     333,  1314,  1938,   334,  3396,   333,  1631,   327,   395,   910,
    1634,  1585,  1637,   327,  1640,   335,  1643,  2208,  1646,  1893,
    1649,  3643,   332,   334,   332,   327,  1652,  2616,   396,   910,
     328,   910,   333,  1660,   399,   335,   402, -1425,  1663,  2313,
     334,   334,  1666,  1667,  1669,  1670,  1671,  1672,  1673,   334,
    1315,  1675,   335,   335,  2847,   328,   332,   332,   333,   332,
     333,   335,   327,   744,   910,  2206,   332,   334,   403,   744,
     328,  2631,  2207,   404,   334,   335,   328,  3365,   910,   335,
    2853,   335,  3035,  2977,  1703,   405,   335,  2838,   328,   910,
    1708,  3714,   333,   333,  2839,   333,  1711,  1712,   406,  1715,
    1716,   334,   333,  3043,   332,  3070,  1723,   332,  2864,   340,
     332,   332,   407,   335,  2336,  2337,  3157,  1873,   408,   589,
     590,   591,   592,  3158,   409,   328,   328,   334,  1732,   334,
     330,  1601,   328,  2338,  1735,  3773,  1736,  1737,  3177,   335,
     333,   335,  2617,   333,   332,  2991,   333,   333,   332,  1747,
     332,  1750,  3389,   410,   327,   415,  1753,   332,   416,  1755,
     418,   334,   334,  3479,   334,   332,   455,   421,  1758,  2848,
     327,   334,  1761,   335,   335,  1764,   335,  1766,  1767,  1768,
     333,  3454,   422,   335,   333,  1275,   333,   332,  1613,  1614,
    1780,  1781,  1782,   333,   424,  2854,   327,  3036,  1915,  1787,
    1788,   333,   425,  1797,  1798,  1799,  1800,   431,   432,   334,
     542,   438,   334,   327,   332,   334,   334,   328,   444,   447,
     332,   335,   327,   333,   335,   449,   450,   335,   335,   393,
     413,   459,   426,   328,   327,  1622,   456,   332,   343,   445,
     327,   460,   463,  1826,  1827,  1828,  1829,  1830,   464,   744,
     333,   327,   332,   334,   467,   334,   333,  2108,  1842,   328,
    1845,   335,   334,  1850,   350,   335,   327,   335,   470,   332,
     334,   332,  2096,   333,   335,   477,   328,   896,  2758,   478,
    1066,   372,   335,  1151,  1309,   328,   479,  1874,   333,  1875,
     457,  1877,   334,   602,  1879,   480,   366,   328,  2096,  2096,
    1885,   481,   374,   328,   335,   333,  1894,   333,   384,  1260,
     482,   483,   328,   486,   328,  1904,  1905,  1906,  1907,   334,
     489,  1396,   327,  1537,   490,   334,   865,   868,   872,   328,
    1547,   335,   493,   327,   386,   332,   332,   335,  1569,   494,
     367,   495,   334,  1941,  1942,  1943,  3140,  3141,   496,  1948,
     497,   498,   499,  1952,   335,   332,   332,   334,   368,   369,
    1601,   500,  1961,   332,  1963,  1964,  1965,  1966,   501,   335,
     502,   333,   333,  1973,   334,  1976,   334,  1977,  1978,  1917,
     503,   504,   505,   506,   507,   328,   335,  1659,   335,  1987,
     397,   333,   333,  1668,   327,  1992,   328,  1994,  1995,   333,
    1997,  1998,  1999,  2000,  2001,  2002,  2003,  2004,  2005,  2006,
    2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,  2016,
    2017,  2018,  2019,  2020,  2021,  2022,  2023,  2024,  2025,  2026,
    2027,  2028,  2029,  2030,  2031,  2032,  2033,  2034,  2035,  2036,
     334,   334,  1746,   508,  1825,   509,   332,   332,  2042,   332,
     332,   510,   335,   335,   517,   518,   327,   328,   522,  2265,
     334,   334,   400,  2046,  2047,  2048,   523,   332,   334,   327,
     525,   526,   335,   335,   332,   332,   527,   327,   332,  1909,
     335,   528,   333,   333,   332,   333,   333,   529,   332,   332,
     327,  2055,   545,   327,   530,   332,  3259,   332,   962,  1889,
     531,  1891,   532,   333,  1895,   533,   534,   605,  1849,  1876,
     333,   333,  2062,   536,   333,   537,   538,   539,  1911,   328,
     333,   540,  2068,  2069,   333,   333,   332,   616,  1884,  1986,
     648,   333,   328,   333,   679,   879,  2067,  2076,   683,  2077,
     328,   332,  2078,   684,   882,   419,   880,  1565,   881,  2081,
     332,   334,   334,   328,   334,   334,   328,  2089,   429,   883,
    2091,   433,   333,   335,   335,  1586,   335,   335,  2094,   900,
     884,  1589,   334,  2097,  2099,  2100,  2101,   333,  2102,   334,
     334,   885,  1602,   334,   335,   894,   333,   895,   898,   334,
     907,   335,   335,   334,   334,   335,   332,  1589,  1589,   912,
     334,   335,   334,   327,  2122,   335,   335,  2125,   913,   327,
    3360,  2128,   335,   332,   335,  2131,   327,   327,  2134,  2098,
    2243,  2137,  2267,  2270,   332,   923,   924,  2144,   925,   926,
    2147,   334,   333,  2150,   327,   927,   934,   327,  2153,   938,
    2275,   950,   332,   335,   327,   951,   334,  2310,  2381,   333,
    2164,  2384,   954,   332,   327,   334,   957,  2394,   335,   960,
     333,  2402,  2404,  2503,   967,   971,   328,   335,  2412,   972,
    2475,   439,   328,   332,   973,   327,   975,   441,   333,   328,
     328,   976,   977,   981,   465,   471,   982,  2182,   332,   333,
     332,   985,  2183,  2184,   986,  2186,  2187,   328,  2189,  2487,
     328,   334,   473,   332,   332,   487,   989,   328,   991,   333,
     992,   994,   513,   335,  2539,  3451,   997,   328,   334,  2204,
     332,   998,   515,  2550,   333,   332,   333,  2209,   999,   334,
     335,  2210,  1005,  2211,  2212,  1006,  1008,  2213,   328,   333,
     333,   335,  2215,   887,  1011,  1012,  2221,   334,  1017,   332,
    2223,  2224,  2225,  1018,  1021,  1023,   333,  1025,   334,   335,
    1027,   333,   332,   332,  2236,  2237,  2238,   332,   332,  2584,
     335,  2242,  2244,  2245,  2246,  1029,  2247,   327,   334,  2254,
    2255,  2256,   332,  2259,  2260,   333,  2588,   332,   327,  2268,
     335,  2271,   332,   334,   327,   334,  1031,  2600,   333,   333,
     332,  1033,   332,   333,   333,   335,  1035,   335,   334,   334,
    2276,  2277,  2278,  2279,  1038,  2602,  1039,  2284,   333,  1046,
     335,   335,  1630,   333,  1055,   334,  2604,  2293,   333,  1068,
     334,  1058,  2298,  1073,  2299,  1074,   333,   335,   333,   332,
     328,  1071,   335,   327,  1075,   921,  2610,   327,  1072,  1078,
    1076,   328,  1077,  1079,   334,  1080,   983,   328,  1081,  1084,
    1086,  2614,   987,  2659,   332,   332,   335,   334,   334,  2311,
     327,   327,   334,   334,   327,   333,  2666,  2695,  1087,   335,
     335,  1088,  1089,   327,   335,   335,   327,   334,  1114,  1127,
     327,  1128,   334,  2759,  1129,  1130,  2333,   334,  2764,   335,
     333,   333,   826,  1156,   335,   334,   328,   334,  2348,   335,
     328,  1000,  1160,  1161,  1633,  1002,  1162,   335,   827,   335,
    1176,  2632,  2792,   332,  1163,   332,  2363,  2364,   332,  1166,
     332,   332,  1168,   328,   328,  2802,  2804,   328,  1044,  1082,
    2807,  2824,  1146,  1178,   334,  2382,   328,  2385,  1187,   328,
    1201,  1157,  1211,   328,  1173,  2840,   335,  2395,  1179,   333,
    2851,   333,  2400,  1214,   333,  2927,   333,   333,  2403,   334,
     334,  2405,  1233,  2994,   332,  3010,  2408,  1234,  1241,  1282,
    2413,   335,   335,  2416,  2417,  2418,  2421,  2424,  2425,  2426,
    2427,  2428,  2429,  2430,  2431,  2432,  2433,  2434,  2435,  2438,
    2439,  2440,  2441,  2442,  2443,  2444,  2445,   332,   327,  1243,
     333,  2450,  2451,  2452,  2453,  2454,  2455,  2456,  2457,  2458,
    2459,  1244,  1246,  1247,  1248,   327,   828,   829,   334,  1249,
     334,  2463,  1251,   334,   327,   334,   334,  3084,  3086,  1636,
     335,  1642,   335,   333,  1252,   335,  2471,   335,   335,  1254,
     830,  1256,  2476,  1265,  1266,  1269,   332,  1270,   332,   332,
    2479,  2480,  1272,  1283,   332,  2481,  1290,   327,  1645,   332,
     327,   328,  2486,  2488,  2490,  1292,  1182,  2492,  2493,   334,
    1392,  2494,   327,  1403,  1298,   332,   332,   332,   328,  2499,
     327,   335,   333,  1204,   333,   333,  3100,   328,  3104,  1404,
     333,  3139,  1209,  3192,  2517,   333,  2519,  2520,  1408,  1412,
     332,  2524,   334,  2526,   327,   327,  2529,  1460,   332,  2532,
    1463,   333,   333,   333,   335,  2413,  2540,   831,   832,  2543,
     328,   327,  2546,   328,  1464,  1296,  1478,  2551,  1304,  2552,
    1472,  2553,  1482,  1476,   833,   328,   333,  3195,   332,   332,
    1306,   332,  1484,   328,   333,  2770,   327,  1485,  1409,  1488,
    1500,   334,  1489,   334,   334,  1490,  2567,  2568,  2569,   334,
    2570,  2571,   332,   335,   334,   335,   335,   328,   328,  1498,
    3198,   335,  1542,  1554,   333,   333,   335,   333,  2580,  1522,
     334,   334,   334,  1491,   328,  2585,  1492,  1523,  2589,  1566,
    1493,  2592,   335,   335,   335,  1494,  1501,  2596,   333,  1495,
     834,   835,   836,   837,  1505,   334,  1509,  2601,  2603,   328,
    2605,  2606,  2607,   334,  1748,  2609,  2611,   335,  2612,  2613,
    2615,  2618,  3200,  2620,  2621,   335,  2622,  3212,  2624,  2625,
    2626,  2627,  3227,  2628,  2629,   327,  1506,  1508,   332,  1510,
    2633,  1511,  2634,   334,   334,  2635,   334,   327,  2637,  3229,
    3257,  2640,  2641,  2642,   327,   335,   335,  2644,   335,   332,
    2647,  1512,  1513,   332,  2653,   327,  1514,   334,  2657,  1515,
     332,   327,  2660,  3258,   333,  1516,  1517,  1518,  2865,   335,
    2667,  3266,  1524,   327,  1525,  1526,  1530,  1527,  1531,  1533,
    1541,  1550,  1551,  1552,  1561,   333,  1562,  1568,   328,   333,
    1844,  1577,   838,  1839,   839,   332,   333,   840,  1578,   841,
     328,  3287,  3304,  1593,  3306,  1857,  1604,   328,  1605,  1608,
    1648,  1609,  1859,  1713,  1612,  1615,  1616,   842,   328,   843,
     844,  1617,  1628,  1865,   328,  3314,  1650,   332,   845,  1927,
     846,   333,   332,   334,  2910,  1651,   328,  1657,   819,   866,
     867,  1968,   847,   848,   849,   335,  1674,   332,   850,   851,
     852,   853,  1677,  1683,   334,   854,   855,  1684,   334,   856,
     857,   858,   859,   333,   860,   334,   335,  2696,   333,  2697,
     335,  1685,  2699,  1686,   861,   862,   863,   335,  1687,  1688,
    2704,  2705,  2706,   333,  2707,  2708,  1689,  2709,  2710,  2711,
    2712,  2713,  2714,  2715,  2716,  2717,  2718,  2719,  2720,  2721,
     334,  2722,  2723,  2724,  2725,  2726,  2727,  2728,  2729,  2730,
     332,   327,   335,  1690,  2735,  2736,  2737,  2738,  2739,  2740,
    2741,  2742,  2743,  2744,  1691,   332,  3320,  2746,  2747,  1692,
    1693,   332,   334,  3322,   332,   332,   332,   334,  2753,  1694,
    1695,  1696,   332,  2756,   335,  1697,   333,  1700,  1704,   335,
    2760,  1705,   334,  2762,  1706,   327,  2765,  1707,  1719,  1731,
    1733,   333,   332,   332,   335,   327,  1734,   333,  3328,  1738,
     333,   333,   333,  2773,   328,  2774,  1739,  1740,   333,  1979,
    1744,  2778,  2779,  2780,  2781,  2782,   332,  2783,  1745,  1751,
    2785,  2786,   327,  2787,  2788,  1754,  2789,  2790,   333,   333,
    2413,  3045,  1756,  1757,  2795,  2796,  1763,  2797,  2798,  1765,
    2799,  2800,  1779,  1776,  1777,   334,  2803,  2805,   328,  1097,
    3331,  1784,   333,  1983,  1812,  2808,  2809,   335,   328,  2812,
     334,  2813,  2814,  1988,  2816,  2817,   334,  1813,   332,   334,
     334,   334,   335,  1814,  1815,   332,  1816,   334,   335,  2825,
    1817,   335,   335,   335,  1818,   328,  1819,  1820,  2830,   335,
    2037,  1821,   327,   327,   327,   327,   327,   334,   334,  2834,
    2835,  2836,  1822,  1823,   333,  2841,  2842,  2843,   327,   335,
     335,   333,  2849,  3348,  2850,  2852,  2855,  1098,  2857,  2858,
    2859,   334,  2861,  2862,  1824,  1835,  1837,  2866,  3350,  1838,
    1851,  2869,  1852,   335,  2872,  2873,  1853,  3362,  3399,  3402,
    1854,   332,   332,  2879,  2880,  3411,   332,   332,   332,  1855,
     332,   332,   332,  1856,  2886,   328,   328,   328,   328,   328,
    2084,  2086,  2115,  2162,  2172,  3417,  3419,  1890,  1896,  1863,
    1864,   328,  1867,   334,   332,  1868,  2226,   333,   333,  1869,
     334,   327,   333,   333,   333,   335,   333,   333,   333,  3425,
    1870,  1871,   335,  1878,  1880,  2912,  1881,  2914,  2915,  2916,
    2917,  2918,  2921,  2922,  2925,  2926,  2928,   327,  2930,  2931,
     333,  2934,  2935,  2936,  2937,  2938,  2939,  2942,  2943,  1882,
    2945,  2946,  2947,  2948,  2949,  1883,  1898,  2952,  2953,  2954,
    2955,  2956,  2957,  2958,  2959,  2960,  2961,  2962,  2963,  2964,
    2965,  3440,  2966,  1897,   328,  2968,   334,   334,  3442,  2228,
    1901,   334,   334,   334,  1902,   334,   334,   334,   335,   335,
    1908,  1921,  1929,   335,   335,   335,  1945,   335,   335,   335,
     328,   327,  2979,  2980,  2981,  2229,  2983,   332,   327,   334,
    2986,  1959,  2988,  1971,  2990,  1972,  2992,  1981,  2413,  1982,
    2996,   335,  2998,  1985,  3000,   327,   332,   327,  1990,  1991,
    1993,   327,  2039,  2040,  2043,  3005,  2044,  3007,  3008,  2049,
    3009,  3011,   327,   333,  3450,  3457,   327,  2052,  2056,  3476,
    3482,  3484,  3017,  3486,  3502,  3504,  2061,  2070,  3019,  3020,
    3021,  3022,   333,  2071,   328,  3024,  3025,  3026,  2072,  2230,
    3028,   328,  2073,  3290,  3032,  2074,  2231,  3524,  2079,  3037,
    2119,  3038,  3039,  3040,   332,  3041,  3042,  2083,   328,  2088,
     328,  2095,  2103,  2261,   328,  2263,  3049,  3050,  2109,  2269,
    2112,  2113,  3051,  2114,  2120,   328,  2123,  2126,  2129,   328,
    2303,   332,   334,  2132,  2305,  2135,  2142,   332,  2145,  2148,
     333,  2158,  2160,  2165,   335,  2171,  2174,   327,   327,   327,
    2175,   334,  2176,  3069,  2177,  3071,  3072,  3073,  3074,  3075,
    3076,  3077,  3078,   335,  3079,  3080,  3081,   333,  3082,  3083,
    3085,  2178,  3087,   333,  3088,  3089,  2181,  3090,  3091,  3092,
    3093,  2185,  3096,  3097,  2188,  3098,  3099,  3101,  2190,  3102,
    3103,  3105,  3106,  3107,  3108,  3109,  3110,  3111,  3114,  3117,
    3525,  3120,  3121,  3122,  3123,  3124,  3125,  3126,   332,   334,
     328,   328,   328,  3129,  3397,  2307,  2319,  2398,  2191,  3528,
    2192,   335,  2195,  3132,  3133,  1466,  1467,  1468,  1469,   327,
    2200,  2203,   327,   327,  1477,  2413,   334,  2205,  3142,   332,
     327,  2219,   334,   327,   333,  2222,   327,  2272,   335,  2282,
    2232,  3149,  3150,  3151,   335,  3152,  2289,  2297,  2504,   327,
    2505,  2506,  2507,  2508,  2509,  2510,  2511,  2512,  3159,  3160,
    3163,  2288,  3164,  2302,  2309,   333,  3167,  3540,  2314,  2317,
     327,  2318,  3172,  3173,  3174,  3175,  3176,  1529,  2321,  2322,
    3178,  3179,   328,  3180,  3181,   328,   328,  2406,  2323,  2324,
    2554,  2638,  2325,   328,  3545,  2326,   328,  2327,  2663,   328,
    3547,  2668,   327,   334,  2693,  3193,  3194,  3196,  3197,  3199,
    3201,  3202,   328,  3204,   327,   335,  3208,  2870,  2328,  3211,
    2329,  3213,  3214,  3215,  3216,  3217,  3218,  3219,  2330,  3220,
    3221,  2331,  3223,   328,   334,  3226,  3228,  3230,  2887,  3231,
    3232,  3233,  3234,  3235,  3236,  3237,   335,  3238,  3239,   327,
    3240,  3241,  2332,  3242,  2334,  3245,  3246,  3247,  3248,  3249,
    3250,  3251,  2335,  2342,  2343,   328,  3254,  3255,  2344,   327,
    2889,  3549,  2345,  3256,  2413,  2413,  2346,   328,   327,  2349,
    2350,  3261,  2891,   332,  3264,  3265,  3267,  2351,  2352,   332,
    2353,  2354,  2355,  3269,  3270,   327,  3271,  3274,  3275,  2356,
    2357,  3280,  3567,  2358,  2359,  2360,  3283,  3284,  3285,  3286,
    3288,  2361,   328,  3291,  3292,  3293,  2362,  2893,  2366,   333,
    2367,  2368,   327,   327,  2369,   333,   332,  3300,  3301,  2370,
    3302,  3303,   328,  3305,   327,  3307,  3308,  2904,  3309,  2371,
    3310,   328,  3311,  2372,  3312,  2373,  2973,  3315,  3316,  3317,
    3318,  3319,  3321,  3323,  3324,  3325,  2374,  3326,   328,  3327,
    3329,  2375,   333,  3057,  3332,   327,  3334,  3335,  3336,  3337,
    3338,  3339,  3340,  2376,  3343,  2377,  3346,   327,  3347,  3349,
    3351,  3352,  3353,  3354,  3355,   328,   328,  2378,   334,  3359,
    3059,  3185,  3361,  2413,   334,  3364,  2379,   328,  3367,  3368,
     335,  3369,  3187,  3370,  3371,  2380,   335,  3373,  2383,  3376,
    2386,   332,  2387,  3378,  3381,  2388,  2389,  2390,  3385,  3386,
    3387,  2391,  3388,  2396,   332,  3390,  3391,  3392,   328,  2397,
     327,   334,  2401,  3295,  2409,  3400,  2460,  3403,  2461,  3404,
     328,  3405,  3406,   335,  3408,  3297,  3410,   333,  2465,  3412,
    3413,  3414,  2468,  3416,  2477,  3418,  3569,  3420,  3421,  3422,
     333,  3424,  3580,  3426,   327,  2478,  3427,   327,  3428,  3429,
    3430,  3431,  3432,  3433,  3434,  2484,  3435,  3436,  2485,  3437,
    3438,  3439,  2491,  3441,  2495,  3443,  3444,  3445,  3446,  3447,
    2496,  2497,  2498,   328,  2413,   332,  2115,   332,  3393,  3582,
    2502,  3455,  2516,  3458,  3459,  3460,  2518,  2521,  2522,  3462,
    2525,  2527,  2528,   332,  3465,   332,   334,  2530,  2531,  3470,
    3471,  3472,  3473,   327,  3474,  3475,  3477,   328,   335,   334,
     328,   333,  3490,   333,  3480,  3495,  2533,  3481,  3483,  3485,
    3487,   335,  3488,  2534,  3489,  2535,  2541,  3492,  3493,   333,
    3494,   333,  2542,  2544,   332,  2545,  3500,  2547,  3501,  2548,
    3503,  3505,  3506,  3507,  3508,  3509,  3510,  3511,  2549,  3513,
     327,  3515,  2561,  3517,  2555,  3518,  2793,  3519,  3520,  3521,
    3522,  3523,  2556,  2557,  3589,  2413,   328,   332,  2559,  3526,
     333,  3497,  2560,  3529,  2562,  2572,  2574,  3601,   327,  2575,
     334,  3533,   334,  2581,  3535,  3536,   327,  3538,  3539,  3541,
    2582,  3542,   335,  2586,   335,  2590,  2593,  3546,   334,  3548,
     334,  3550,  2595,   333,  2597,  2598,  3557,  3560,  2599,  2619,
     335,   332,   335,   328,  3565,   327,  2636,  3568,  3553,  3570,
    3571,  3572,  3573,  3574,  3575,  3576,   332,  3577,   327,  3578,
    2643,  2652,  3581,  3583,  3584,  3585,  3586,  3587,  2654,   334,
    3588,   328,   327,  3590,  2658,  2661,  3597,   333,  3611,   328,
    3615,   335,  3594,  2665,  3599,  3595,  3596,  2670,   332,  2671,
    2675,  2676,   333,  2683,  3602,  2684,  3617,  2685,  3625,  2686,
    3604,  3605,   334,  3606,  3607,  2687,  3608,  2688,   328,  2689,
    2690,  2691,  3612,  3613,   335,  3616,  3618,  3619,  3620,  3621,
    3622,   328,  2698,  2749,   333,  3626,  3670,  3628,  3631,  3634,
    2754,  3637,  3639,  2755,  2757,   328,  2761,  3627,  3644,  3645,
    3672,  2763,  2766,  2767,  2768,  2771,   334,  2777,  3651,  2784,
    3654,  2791,  3657,  2801,  3658,  2806,  3660,  2810,   335,  2815,
    3662,   334,  3663,  3665,  3667,  3668,  3669,  2826,  2827,  2829,
    3638,  2831,  2832,   335,  3674,  3675,  2833,  3676,  3677,   332,
    3678,  2844,  2856,  2860,  2867,   332,  3683,  2868,  3684,  2874,
    2875,  2876,  2881,   334,  2882,  3687,  2883,  3688,  3689,  2897,
    3690,  3691,  3692,   332,  3694,   335,  2899,  2900,  2902,  3697,
     332,  3698,  3699,  3700,  3659,   333,  2906,  2929,  3703,   332,
    3706,   333,  3709,  2944,  3710,  2970,  2971,  2975,  3715,  3664,
    3716,  2976,  3718,  2978,  3720,  2982,  3722,  2984,  3724,   333,
    2985,  2987,  2989,  3728,  3730,  2993,   333,  3731,  2995,  3732,
    3733,  2997,  3734,  3735,  3736,   333,  2999,  3737,  3001,  3738,
    3739,  3666,  3740,   332,  3741,  3002,  3006,  3012,  3744,   332,
     332,  3015,  3745,   332,  3746,  3016,  3748,  3018,  3750,  3023,
    3752,  3753,  3754,  3755,   334,  3027,  3033,  3048,  3759,   332,
     334,  3052,  3762,  3055,  3763,  3056,   335,   332,  3766,   333,
    3061,  3062,   335,  3769,  3063,   333,   333,  3064,   334,   333,
    3774,  3775,  3065,  3777,  3066,   334,  3778,  3067,  3779,  3128,
     335,  3782,  3783,  3130,   334,   333,  3786,   335,  3788,  3131,
    3790,  3134,  3791,   333,   329,   331,   335,  3135,  3795,  3796,
    3797,   344,   345,  3136,  3137,   348,  3138,   351,  3143,  3144,
    3145,  3148,  3693,  3153,  3154,  3155,  3156,  3168,  3723,  3184,
    3189,  3203,   373,   375,  3207,  3222,  3252,  3253,   334,  3268,
     382,   383,   385,   387,   334,   334,  3727,  3276,   334,  3294,
     335,  3313,   398,  3729,   401,  3330,   335,   335,  3333,  3356,
     335,  3357,  3743,  3358,   334,  3366,  3372,  3382,   417,  3395,
     420,  3398,   334,   423,  3401,  3407,   335,  3409,   430,  3415,
    3423,   434,  3448,  3449,   335,   440,   442,   443,  3452,  3453,
    3456,   448,  3461,  3463,   454,   458,  3464,  3469,   461,   462,
    3478,  3499,   466,  3512,   468,   469,  3758,   472,   474,  3514,
    3516,  3527,  3768,  3776,  3530,  3531,  3781,   484,  3532,  3534,
     488,  3537,  3543,  3544,  3551,  3552,  3563,  3564,  3566,  3579,
    3591,  3593,  3785,  3623,  3624,  3640,  3646,  3647,  3648,  3649,
    3789,  3661,  3681,  3685,   514,   516,  3695,  3696,   519,   520,
     521,  3701,  3702,  3711,  3712,  3717,  3719,  3721,  3725,  3726,
    3742,  3747,  3749,  3751,  3756,  3757,  3760,  3761,  3764,  3765,
    3770,  3771,  3780,  3787,  3792,  3798,   544,   546,  1107,  1110,
    1661,   904,   580,   902,  1159,  1658,   587,  1626,   599,   601,
     603,  1575,   914,  1004,  1596,  1579,   945,   974,   608,   776,
     777,  1050,  1624,   778,  1019,  1803,   952,  2080,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   869,
       0,     0,     0,     0,     0,     0,   651,     0,   654,   655,
       0,     0,   657,   659,   661,   662,     0,   664,     0,   666,
     668,   670,   672,     0,     0,   664,   668,   672,     0,     0,
     678,     0,   680,   681,   682,     0,     0,     0,   687,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   690,   692,   693,     0,     0,   694,   695,
     697,   698,   701,     0,   729,   599,   599,   733,   734,   735,
     738,     0,     0,   745,   747,   749,     0,     0,   750,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   796,   798,   799,     0,     0,     0,     0,     0,
       0,     0,     0,   878,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   888,
       0,     0,   892,   892,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   901,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   666,   664,   918,
       0,     0,   922,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   932,   932,     0,   937,
       0,     0,     0,   668,     0,     0,   672,   670,     0,     0,
     948,   878,     0,     0,     0,   918,   878,     0,     0,   956,
       0,     0,   959,     0,     0,     0,   965,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   690,     0,
       0,     0,     0,     0,     0,   729,     0,   980,     0,     0,
       0,     0,   984,     0,     0,     0,   988,     0,     0,     0,
     990,     0,     0,   993,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1001,  1003,   678,     0,     0,     0,
       0,     0,     0,   878,     0,     0,     0,     0,     0,  1016,
       0,     0,     0,     0,   692,     0,     0,     0,     0,     0,
       0,     0,     0,   878,     0,     0,     0,     0,     0,     0,
       0,     0,  1037,     0,     0,     0,     0,     0,     0,     0,
    1043,  1045,     0,     0,     0,     0,     0,  1049,   697,  1053,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1057,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1083,     0,  1085,     0,     0,     0,     0,
       0,  1091,  1092,  1093,  1094,     0,     0,     0,  1100,  1101,
    1102,  1103,  1104,     0,  1106,  1106,  1109,  1109,     0,     0,
       0,     0,     0,  1115,     0,  1118,  1120,     0,  1122,  1124,
    1126,     0,     0,     0,     0,     0,  1145,  1147,  1148,     0,
       0,     0,     0,     0,  1158,     0,     0,     0,     0,     0,
    1165,     0,     0,     0,  1171,     0,     0,     0,     0,  1174,
       0,     0,     0,     0,     0,  1177,     0,  1180,  1181,  1183,
       0,     0,  1120,     0,  1122,     0,  1188,     0,  1191,     0,
    1193,     0,     0,     0,  1196,     0,  1197,     0,  1198,     0,
    1199,     0,  1200,     0,     0,     0,     0,     0,     0,  1205,
       0,     0,     0,  1210,     0,  1120,     0,  1122,     0,     0,
       0,     0,  1223,  1224,  1225,  1226,     0,  1227,  1228,     0,
       0,  1232,     0,     0,  1235,     0,     0,     0,     0,     0,
       0,  1242,     0,     0,  1245,     0,     0,     0,     0,     0,
       0,     0,  1253,     0,     0,     0,     0,  1258,     0,     0,
       0,     0,     0,     0,  1262,  1264,     0,     0,     0,  1268,
       0,     0,     0,     0,     0,  1271,     0,  1273,  1274,  1276,
       0,     0,     0,     0,     0,     0,     0,  1287,  1289,     0,
       0,     0,     0,  1295,  1297,     0,  1303,  1305,  1307,  1308,
       0,  1311,     0,  1316,     0,     0,     0,     0,     0,     0,
       0,     0,  1394,  1395,     0,     0,  1401,  1402,     0,     0,
       0,  1406,  1407,     0,  1410,     0,  1411,     0,     0,  1456,
    1457,  1458,  1459,     0,  1461,  1462,     0,     0,     0,     0,
       0,     0,     0,     0,  1471,     0,  1475,     0,     0,     0,
    1479,  1480,     0,  1481,     0,     0,     0,     0,     0,  1487,
       0,     0,     0,     0,     0,     0,     0,     0,  1496,  1497,
       0,     0,     0,     0,  1502,  1503,     0,     0,     0,  1507,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1519,  1520,  1521,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1532,     0,  1535,     0,  1122,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1543,
       0,   892,     0,   892,   892,     0,     0,     0,     0,  1549,
       0,     0,     0,     0,  1147,     0,     0,     0,     0,     0,
    1558,     0,     0,     0,     0,     0,   918,     0,     0,     0,
       0,     0,     0,     0,     0,  1567,     0,     0,     0,     0,
     932,     0,  1573,   932,   932,   937,     0,  1576,     0,     0,
       0,     0,     0,     0,     0,     0,  1581,  1583,     0,  1535,
       0,     0,   918,  1535,     0,     0,     0,     0,     0,     0,
       0,     0,  1595,  1598,  1600,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1607,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1610,     0,     0,  1611,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1618,     0,  1621,     0,
    1535,     0,     0,     0,  1016,     0,  1627,     0,     0,     0,
       0,  1629,     0,     0,  1535,  1632,     0,  1635,     0,  1638,
       0,  1641,     0,  1644,     0,  1647,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1049,     0,     0,
       0,  1053,     0,  1662,     0,     0,  1664,  1665,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1676,     0,
    1678,  1413,  1414,     0,     0,     0,     0,     0,     0,     0,
       0,  1415,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1698,  1699,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1709,
       0,  1710,     0,     0,     0,  1714,     0,     0,     0,     0,
    1122,  1720,     0,  1724,  1725,  1726,  1727,     0,     0,     0,
       0,     0,     0,  1728,     0,  1729,     0,  1730,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1749,     0,     0,  1752,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1759,  1760,     0,  1762,
       0,     0,     0,     0,     0,     0,  1770,  1770,  1770,  1770,
    1770,  1775,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1786,     0,     0,     0,  1795,  1796,     0,
       0,     0,     0,  1802,  1802,  1804,  1805,  1770,  1807,     0,
    1808,     0,  1811,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1834,     0,     0,     0,  1836,     0,
       0,     0,     0,  1840,     0,  1843,     0,     0,  1848,     0,
       0,     0,     0,     0,     0,     0,  1858,     0,  1860,     0,
    1862,     0,     0,     0,  1866,     0,     0,  1416,  1417,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1900,     0,     0,
    1903,     0,     0,     0,     0,     0,  1910,     0,  1912,  1913,
    1914,  1916,  1918,  1920,     0,  1922,  1923,  1924,  1925,  1926,
    1928,     0,  1930,  1933,     0,  1936,  1937,  1939,  1940,     0,
       0,     0,  1944,     0,  1947,     0,  1949,  1950,  1951,     0,
    1953,  1954,  1955,  1956,  1957,  1958,     0,  1960,     0,  1962,
       0,     0,     0,     0,  1967,  1969,  1970,     0,     0,     0,
       0,     0,     0,     0,     0,  1980,     0,     0,     0,     0,
    1984,     0,     0,     0,     0,     0,  1989,     0,     0,  1418,
    1419,  1420,  1421,  1422,     0,  1423,     0,  1424,  1425,  1426,
    1427,  1428,  1429,  1430,  1431,     0,  1432,  1433,  1434,  1435,
    1436,  1437,  1438,  1439,  1440,  1441,  1442,  1443,  1444,  1445,
       0,  1446,  1447,  1448,  1449,  1450,  1451,  1452,  1453,  1454,
    1455,     0,     0,     0,     0,     0,  2038,     0,     0,     0,
       0,     0,  2041,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2045,
     826,     0,     0,     0,  2050,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   827,     0,     0,     0,
    2051,     0,     0,  2053,     0,     0,     0,  2054,     0,     0,
       0,     0,  2057,  2058,     0,     0,     0,     0,     0,     0,
    2059,  2060,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2064,     0,     0,     0,  1122,  2066,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   892,     0,     0,     0,
    2075,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2085,     0,
    2087,     0,     0,     0,   932,     0,  2090,     0,     0,  1583,
       0,  2093,     0,  2093,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1600,     0,  2105,     0,
    2105,     0,     0,     0,   828,   829,     0,     0,     0,     0,
    2117,  2118,     0,     0,     0,     0,     0,     0,     0,  2121,
       0,     0,  2124,     0,     0,     0,     0,  2127,   830,     0,
    2130,     0,     0,  2133,     0,     0,  2136,     0,     0,  2140,
       0,     0,  2143,     0,     0,  2146,     0,     0,  2149,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1770,     0,  2161,  2163,     0,     0,     0,  2166,
    2167,  2168,  2169,  2170,     0,     0,  2173,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2180,     0,     0,     0,   831,   832,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2194,     0,   833,     0,  2199,     0,  2202,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2214,     0,     0,     0,  2217,
    2218,     0,     0,     0,     0,     0,     0,     0,     0,  2227,
       0,  2227,  2227,  2227,  2227,     0,  2233,  2234,     0,     0,
       0,     0,     0,     0,     0,     0,  2241,     0,   834,   835,
     836,   837,     0,  2248,     0,     0,     0,     0,     0,     0,
       0,  2262,     0,  2262,  2266,     0,  2227,     0,     0,     0,
    2273,  2274,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2287,     0,     0,     0,     0,     0,
       0,  2290,     0,     0,  2294,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2300,  2301,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2304,  2306,     0,  2308,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2316,     0,
     838,     0,   839,  2320,     0,   840,     0,   841,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   842,  2341,   843,   844,     0,
       0,     0,     0,  2347,     0,     0,   845,     0,   846,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     847,   848,   849,     0,  2365,     0,   850,   851,   852,   853,
       0,     0,     0,   854,   855,     0,     0,   856,   857,   858,
     859,     0,   860,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   861,   862,   863,     0,     0,  2399,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2407,
       0,     0,     0,     0,  2410,  2411,     0,  2414,  2415,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2446,  2447,  2448,  2449,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2462,     0,     0,  2464,     0,
       0,     0,     0,  2466,  2467,     0,     0,  2469,  2470,     0,
       0,     0,     0,     0,     0,  2472,  2474,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2489,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2500,     0,     0,
       0,     0,     0,   729,     0,     0,     0,     0,  2515,     0,
       0,     0,     0,     0,     0,     0,     0,  2523,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2227,
       0,     0,     0,     0,     0,     0,  2558,     0,     0,     0,
       0,     0,     0,  2563,     0,     0,     0,     0,     0,  2564,
    2565,  2180,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2573,     0,     0,     0,     0,  2576,     0,
    2577,  2578,     0,  2579,     0,     0,     0,     0,     0,     0,
    2583,     0,     0,  2587,     0,     0,  2591,     0,     0,     0,
    2594,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2608,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2623,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2639,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2651,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2662,  2664,     0,     0,     0,     0,     0,     0,  2669,     0,
       0,     0,     0,     0,     0,  2672,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2679,     0,  2680,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2692,
    2694,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2700,  2701,     0,     0,  2702,  2703,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2731,  2732,  2733,  2734,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2745,     0,     0,     0,     0,  2748,     0,     0,
       0,     0,  2751,  2752,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2772,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2811,     0,  2564,     0,     0,     0,
       0,     0,     0,  2818,     0,     0,     0,  2820,     0,  2822,
       0,     0,     0,  2823,     0,     0,     0,     0,     0,     0,
       0,  2828,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2837,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1317,  2871,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1318,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2888,     0,  2890,  1319,     0,  2892,     0,  2894,
    1320,     0,     0,     0,     0,     0,     0,  2898,     0,  1321,
       0,  2901,  1322,  1323,     0,     0,     0,  1324,  1325,     0,
       0,     0,     0,     0,  2903,     0,  2905,     0,     0,     0,
    2911,     0,  2913,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2950,  2951,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1326,     0,     0,     0,     0,     0,     0,
    2967,     0,  2969,     0,     0,     0,  2972,     0,     0,     0,
    2974,     0,     0,     0,  1327,  1328,  1329,  1330,     0,     0,
       0,     0,     0,  1331,     0,     0,     0,  1332,     0,     0,
    1333,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1334,     0,  1335,     0,     0,     0,  3003,  3004,
       0,     0,     0,     0,     0,     0,     0,     0,  1336,  3013,
       0,  3014,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1337,  1338,     0,     0,     0,     0,     0,     0,
       0,  1339,     0,  1340,     0,     0,  1341,     0,  1342,     0,
    1343,  1344,     0,     0,     0,  1345,     0,  1346,  1347,  1348,
    1349,     0,     0,     0,     0,     0,  3046,     0,     0,  3047,
       0,     0,     0,     0,     0,  1350,     0,  1351,     0,     0,
       0,     0,     0,  1352,     0,     0,     0,     0,  3058,     0,
    3060,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3068,     0,
       0,     0,  1353,     0,     0,     0,     0,  1354,  1355,  1356,
       0,  1357,     0,  1358,     0,  1359,     0,     0,     0,     0,
    1360,  1361,  1362,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1363,     0,     0,     0,     0,     0,     0,     0,     0,  1364,
       0,     0,     0,     0,     0,     0,     0,  3127,     0,     0,
       0,  1365,  1366,  1367,  1368,  1369,  1370,  1371,     0,  1372,
    1373,  1374,  1375,     0,     0,     0,     0,     0,     0,  1376,
    1377,     0,  1378,     0,     0,     0,  1379,     0,     0,     0,
       0,     0,     0,  3146,  3147,     0,     0,     0,     0,  1380,
       0,  1381,     0,     0,     0,     0,  1382,     0,     0,     0,
       0,  1383,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1384,     0,     0,     0,  1385,  1386,  1387,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3186,     0,
    3188,     0,     0,     0,     0,  1388,     0,     0,  1389,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1390,     0,
       0,     0,     0,     0,     0,     0,  3260,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  3296,     0,  3298,     0,
       0,  3299,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    3363,     0,  2180,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  3384,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  3394,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2564,  2180,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  3491,     0,     0,     0,     0,     0,  3496,     0,
    3498,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2564,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  3554,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3592,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  3598,     0,  3600,     0,
       0,     0,     0,     0,  3603,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    3614,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  3641,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  3650,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  3671,     0,  3673,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  3682,     0,     0,     0,     0,     1,     0,     2,     3,
    3686,     0,     0,     0,     4,     0,     5,     0,     0,     0,
       0,     0,     0,     0,     6,     0,     7,     8,     9,    10,
      11,     0,     0,     0,    12,     0,     0,     0,     0,     0,
       0,     0,     0,    13,    14,    15,     0,     0,     0,     0,
       0,     0,     0,    16,     0,    17,    18,    19,     0,    20,
       0,     0,     0,     0,     0,     0,     0,    21,     0,     0,
       0,     0,     0,    22,    23,     0,     0,    24,     0,     0,
       0,    25,    26,    27,    28,     0,     0,     0,    29,     0,
      30,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    31,     0,     0,    32,    33,    34,    35,   543,
      36,     0,    37,    38,     0,  3767,     0,     0,     0,     0,
       0,    39,     0,     0,     0,    40,    41,     0,     0,     0,
      42,    43,     0,     0,    44,    45,     0,     0,     0,  3784,
       0,     0,     0,     0,    46,     0,     0,     0,     0,     0,
       0,  3793,    47,  3794,     0,    48,    49,     0,     0,    50,
       0,     0,    51,     0,     0,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,     0,     0,     0,
       0,    64,    65,     0,     0,     0,    66,     0,    67,     0,
      68,    69,    70,    71,     0,    72,     0,    73,     0,     0,
      74,     0,     0,     0,     0,    75,     0,     0,    76,    77,
      78,    79,    80,     0,     0,    81,     0,     0,    82,    83,
      84,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    85,     0,     0,     0,     0,
       0,     0,     0,     0,    86,     0,     0,     0,     0,     0,
      87,     0,     0,     0,    88,     0,     0,     0,     0,    89,
       0,    90,    91,     0,     0,     0,     0,     0,    92,     0,
      93,     0,     0,    94,    95,     0,     0,    96,     0,     0,
      97,    98,    99,     0,   100,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   101,   102,   103,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   104,
     105,     0,     0,     0,     0,   106,   107,   108,     0,   109,
       0,     0,     0,     0,     0,   110,     0,   111,   112,     0,
       0,     0,   113,     0,     0,     0,     0,     0,     0,   114,
       0,     0,   115,   116,     0,     0,     0,   117,     0,     0,
     118,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   119,     0,     0,   120,     0,     0,
       0,     0,   121,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   122,     0,   123,
       0,   124,     0,     0,   125,   126,   127,     0,   128,   129,
       0,   130,     0,     0,   131,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   132,     0,     0,     0,     0,   133,
     134,   135,     0,     0,     0,     0,     0,     0,     0,   136,
     137,   138,   139,   140,   141,     0,   142,     0,     0,     0,
       0,     0,   143,   144,     0,     0,     0,     0,     0,     0,
     145,     0,     0,     0,   146,   147,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   148,   149,
       0,     0,     0,     0,     0,   150,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   151,     0,
     152,   153,   154,   155,   156,   157,   158,   159,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   160,   161,
       0,     0,   162,   163,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   165,   166,
     167,     0,     0,   168,   169,     0,     1,     0,     2,     3,
       0,     0,     0,     0,     4,   170,     5,     0,     0,     0,
       0,     0,     0,     0,     6,     0,     7,     8,     9,    10,
      11,     0,     0,     0,    12,     0,     0,     0,     0,     0,
       0,     0,     0,    13,    14,    15,     0,     0,     0,     0,
       0,     0,     0,    16,     0,    17,    18,    19,     0,    20,
       0,     0,     0,     0,     0,     0,     0,    21,     0,     0,
       0,     0,     0,    22,    23,     0,     0,    24,     0,     0,
       0,    25,    26,    27,    28,     0,     0,     0,    29,     0,
      30,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    31,     0,     0,    32,    33,    34,    35,     0,
      36,     0,    37,    38,     0,     0,     0,     0,     0,     0,
       0,    39,     0,     0,     0,    40,    41,     0,     0,     0,
      42,    43,     0,     0,    44,    45,     0,     0,     0,     0,
       0,     0,     0,     0,    46,     0,     0,     0,     0,     0,
       0,     0,    47,     0,     0,    48,    49,     0,     0,    50,
       0,     0,    51,     0,     0,    52,    53,    54,    55,    56,
      57,    58,    59,    60,    61,    62,    63,     0,     0,     0,
       0,    64,    65,     0,     0,     0,    66,     0,    67,     0,
      68,    69,    70,    71,     0,    72,     0,    73,     0,     0,
      74,     0,     0,     0,     0,    75,     0,     0,    76,    77,
      78,    79,    80,     0,     0,    81,     0,     0,    82,    83,
      84,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    85,     0,     0,     0,     0,
       0,     0,     0,     0,    86,     0,     0,     0,     0,     0,
      87,     0,     0,     0,    88,     0,     0,     0,     0,    89,
       0,    90,    91,     0,     0,     0,     0,     0,    92,     0,
      93,     0,     0,    94,    95,     0,     0,    96,     0,     0,
      97,    98,    99,     0,   100,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   101,   102,   103,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   104,
     105,     0,     0,     0,     0,   106,   107,   108,     0,   109,
       0,     0,     0,     0,     0,   110,     0,   111,   112,     0,
       0,     0,   113,     0,     0,     0,     0,     0,     0,   114,
       0,     0,   115,   116,     0,     0,     0,   117,     0,     0,
     118,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   119,     0,     0,   120,     0,     0,
       0,     0,   121,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   122,     0,   123,
       0,   124,     0,     0,   125,   126,   127,     0,   128,   129,
       0,   130,     0,     0,   131,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   132,     0,     0,     0,     0,   133,
     134,   135,     0,     0,     0,     0,     0,     0,     0,   136,
     137,   138,   139,   140,   141,     0,   142,     0,     0,     0,
       0,     0,   143,   144,     0,     0,     0,     0,     0,     0,
     145,     0,     0,     0,   146,   147,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   148,   149,
       0,     0,     0,     0,     0,   150,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   151,     0,
     152,   153,   154,   155,   156,   157,   158,   159,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   160,   161,
       0,     0,   162,   163,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   164,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   165,   166,
     167,     0,     0,   168,   169,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   170
};

static const yytype_int16 yycheck[] =
{
       6,     7,     8,     9,  1640,   445,  2181,   891,  1204,   893,
    1741,   198,   964,  1209,   930,  1211,   947,   933,   343,   397,
      28,  2847,    12,    17,    84,    38,    76,  2853,  1197,  1198,
    1199,  1200,    19,    38,   181,   370,    32,    38,  2117,   175,
      38,   228,   249,    49,   181,   103,    76,  1778,   249,    76,
     249,   249,  1783,   181,  1785,    51,    22,   125,  1227,    65,
     249,   181,   112,   249,   181,    76,   249,   111,   249,    82,
      76,   258,   249,   126,   249,    75,   231,    82,    86,   249,
      88,    82,   112,    89,    82,   112,    10,   100,   101,     4,
      76,    15,   249,    83,   249,   100,   101,   244,   249,   100,
     101,   112,   100,   101,     7,   249,    76,   244,   113,    76,
     245,   101,   118,   249,   372,   249,   244,   123,   124,    44,
      45,   181,     9,    76,   244,    76,   112,   244,   181,   181,
     181,   181,   140,    76,   124,   135,   136,   137,   138,   139,
     249,    76,   112,   193,    68,   112,     8,   153,   537,   538,
     152,   181,    76,   249,   181,    22,    18,   147,   164,   112,
     249,   112,   249,   193,   181,   916,   193,   133,    30,   112,
     181,    33,    34,    35,    36,    37,   181,   112,   184,   185,
      47,   249,   193,   249,   244,   191,   249,   113,   112,   249,
     114,   244,   244,   244,   160,   181,   202,   125,    76,   249,
     108,   952,   379,  3029,  3030,   249,   407,   193,    75,  3035,
     215,   181,    14,   203,   181,    77,    76,    76,    20,   249,
      76,   211,   249,   193,   249,   247,   193,   244,   181,   219,
     181,    76,   249,   227,   112,   249,   249,   249,   181,   244,
     193,    75,   193,   378,   249,    76,   181,   249,   249,   175,
     193,   249,   112,   112,    76,   372,   112,   181,   193,    76,
     250,   450,   252,   249,   270,    67,   181,   112,   249,   193,
     260,   249,   458,   239,    76,   333,    38,   458,   249,   249,
     450,   112,   249,   458,   146,   413,   292,   293,   150,   215,
     112,   215,   412,   160,   181,   112,   249,   504,   249,   305,
     499,   499,   302,   181,   440,   305,   249,   458,   445,   176,
     112,   249,   249,   460,   249,   193,   499,   500,   501,   372,
      82,   181,   181,   249,   458,   181,   326,   378,   318,   244,
     108,    46,   338,   193,   193,   341,   181,   193,   100,   101,
    3166,   249,   348,  3169,  3170,   432,   249,   353,   193,   458,
     181,  2082,   134,   359,   360,   407,   223,   244,    76,   181,
     366,   367,   193,   369,   181,   700,   249,   363,   249,   458,
    1566,   193,   239,   240,   379,   381,   193,    92,    93,   181,
      39,   248,   384,   245,   406,   402,   181,   393,   394,   449,
     249,   193,   254,   249,   112,   257,   263,   408,   404,   405,
     262,   263,   264,   249,   266,   413,   240,   269,   458,   276,
     181,   273,   418,   181,   422,   421,   218,   404,   424,   606,
     111,   427,    76,    82,   414,   249,   440,   440,   458,   442,
     420,   458,   249,   458,  2513,   440,  2515,   249,   249,   440,
     446,   437,   440,   249,   246,   249,   458,   249,   448,   244,
     444,  3277,   249,   379,   256,   379,  3282,   319,   112,   181,
     119,   249,   125,   181,   249,   652,   456,   249,   302,   440,
     272,   249,   458,   244,   181,   193,   244,   458,   372,   134,
     458,   343,   488,   213,   214,   409,   348,   249,   458,   495,
     321,   458,   326,  1662,    55,   503,   504,   505,   685,   111,
     508,   509,   161,  2139,  2235,   458,   373,   458,  2239,   434,
     435,   436,   437,   438,   440,   458,   383,   249,   385,   181,
     458,   458,   244,   458,   249,   184,   332,   181,    38,   535,
     590,   249,   334,   335,   412,   332,   362,   244,   268,   193,
     181,   547,   548,   549,   550,   551,   552,   553,   554,   249,
     556,   175,   412,   246,   560,   125,   934,  3383,   249,   381,
     384,   181,   379,   365,   366,   432,   111,   412,   379,   894,
     181,   377,    82,   111,   580,    76,   181,   458,   231,   585,
     586,   587,   244,  1779,   410,   177,   249,   317,  1784,   181,
     100,   101,    76,   599,   249,   601,   249,    39,   604,   458,
     181,   181,   458,   244,   249,   181,   181,   125,   175,   615,
     378,   112,   458,   619,   620,   621,   622,   412,   249,   625,
     181,    63,   484,    76,   244,   249,   632,    69,   112,   249,
      76,   181,   638,   244,   640,   497,   498,   249,   134,   244,
      82,   412,   249,   181,  1596,   651,   458,    89,  1579,    38,
     125,   657,   244,   659,   456,   661,   662,   249,  1574,   112,
     390,   379,  1546,   244,   244,   249,   112,   397,   244,   244,
     458,   181,   372,   181,   680,   681,   682,   119,   181,   249,
     181,   125,   249,   244,   690,   415,   345,   346,   249,   876,
     412,   406,   193,    82,   244,   249,   174,   181,   181,   705,
     706,   707,   708,   709,   249,   412,   244,   175,   440,   193,
     403,   100,   101,   719,   181,   249,   378,   125,   724,   161,
     726,   249,   125,   729,   378,   249,   181,   733,   181,   249,
     249,   249,   738,   121,   244,   181,   244,   378,   249,   249,
     193,   244,   184,   749,   750,   751,   752,   193,   249,   755,
     756,  2482,  2483,   249,   760,    38,   762,   583,   584,   585,
     586,   244,   949,   249,   770,   249,   953,   378,    76,   249,
     375,   249,   249,   779,   249,   781,    38,   244,   784,   785,
      62,   249,   402,   321,   790,   175,   249,   407,   503,   244,
     505,   506,   234,   181,    76,   375,   249,   378,   513,    82,
     375,   249,   378,    38,   112,   249,   125,   813,   249,   157,
     525,   526,   527,   528,   529,   407,   531,   100,   101,   825,
      82,   125,   828,  1010,   212,   249,   541,   542,   378,   544,
     112,   546,   547,   181,    62,   841,   551,  1024,   100,   101,
     846,   249,   557,   558,   559,   560,   249,    82,    76,   579,
     580,   581,   372,   372,   249,   338,   244,    75,    76,   249,
     249,   372,   181,   869,   870,   100,   101,   375,   160,   311,
      38,   338,   375,   181,   880,   881,   195,   181,   379,    76,
     886,   125,  1217,   338,   112,   193,   892,  2083,  2084,   249,
     249,   897,   372,   181,   112,   379,   244,   903,   181,   181,
     906,   907,  2538,   345,   346,   285,   286,    38,    76,   372,
     249,   193,   918,   411,    82,   112,   249,   923,   416,   372,
     117,   125,   928,   929,   372,   244,   181,    75,    76,   375,
     249,   223,   100,   101,   242,   175,   242,   181,   298,   249,
     244,   249,   948,   249,   112,    76,    76,   239,   954,   955,
     956,    82,   958,   181,    76,   961,   244,   372,   400,   965,
     966,   244,   244,   181,   112,   193,   249,   249,    76,   100,
     101,    75,    76,   979,    76,   193,   231,   181,   157,    76,
      49,   112,   112,   249,   181,   483,   249,   249,    17,   244,
     112,   997,   998,   157,   249,   125,   193,   249,   120,   157,
     244,  1007,   181,  1009,   112,   503,    76,  1013,   112,   249,
     112,   371,   300,   181,   249,   112,  1022,   181,   249,   175,
    1026,   249,  1028,   181,  1030,   193,  1032,  2758,  1034,    58,
    1036,   249,    76,   181,    76,   181,  1042,   134,   249,   175,
     244,   175,   112,  1049,   249,   193,   249,   244,  1054,   315,
     181,   181,  1058,  1059,  1060,  1061,  1062,  1063,  1064,   181,
     420,  1067,   193,   193,   134,   244,    76,    76,   112,    76,
     112,   193,   181,   181,   175,   242,    76,   181,   249,   181,
     244,   249,   249,   249,   181,   193,   244,  3262,   175,   193,
     134,   193,   134,   249,  1100,   249,   193,   242,   244,   175,
    1106,   249,   112,   112,   249,   112,  1112,  1113,   249,  1115,
    1116,   181,   112,   249,    76,   249,  1122,    76,   249,   249,
      76,    76,   249,   193,   230,   231,   242,   249,   249,   198,
     199,   200,   201,   249,   249,   244,   244,   181,  1144,   181,
     249,   249,   244,   249,  1150,   249,  1152,  1153,   249,   193,
     112,   193,   249,   112,    76,  2791,   112,   112,    76,  1165,
      76,  1167,   249,   249,   181,   249,  1172,    76,   249,  1175,
     249,   181,   181,   249,   181,    76,   125,   249,  1184,   249,
     181,   181,  1188,   193,   193,  1191,   193,  1193,  1194,  1195,
     112,  3366,   249,   193,   112,   117,   112,    76,   997,   998,
    1206,  1207,  1208,   112,   249,   249,   181,   249,   225,  1215,
    1216,   112,   249,  1219,  1220,  1221,  1222,   249,   249,   181,
       0,   249,   181,   181,    76,   181,   181,   244,   249,   249,
      76,   193,   181,   112,   193,   249,   249,   193,   193,   249,
     249,   249,   249,   244,   181,    59,   195,    76,   249,   249,
     181,   249,   249,  1259,  1260,  1261,  1262,  1263,   249,   181,
     112,   181,    76,   181,   125,   181,   112,  1602,  1274,   244,
    1276,   193,   181,  1279,   249,   193,   181,   193,   249,    76,
     181,    76,  1587,   112,   193,   249,   244,   249,  2484,   249,
     249,   249,   193,   249,   249,   244,   249,  1303,   112,  1305,
     249,  1307,   181,   223,  1310,   249,    38,   244,  1613,  1614,
    1316,   249,   249,   244,   193,   112,  1322,   112,   249,   148,
     249,   249,   244,   249,   244,  1331,  1332,  1333,  1334,   181,
     249,   249,   181,   249,   249,   181,   318,   319,   320,   244,
     249,   193,   249,   181,   249,    76,    76,   193,   249,   249,
      82,   249,   181,  1359,  1360,  1361,  2992,  2993,   249,  1365,
     249,   249,   249,  1369,   193,    76,    76,   181,   100,   101,
     249,   249,  1378,    76,  1380,  1381,  1382,  1383,   125,   193,
     249,   112,   112,  1389,   181,  1391,   181,  1393,  1394,   227,
     249,   249,   249,   249,   249,   244,   193,   249,   193,  1405,
     249,   112,   112,   249,   181,  1411,   244,  1413,  1414,   112,
    1416,  1417,  1418,  1419,  1420,  1421,  1422,  1423,  1424,  1425,
    1426,  1427,  1428,  1429,  1430,  1431,  1432,  1433,  1434,  1435,
    1436,  1437,  1438,  1439,  1440,  1441,  1442,  1443,  1444,  1445,
    1446,  1447,  1448,  1449,  1450,  1451,  1452,  1453,  1454,  1455,
     181,   181,   249,   249,   249,   249,    76,    76,  1464,    76,
      76,   249,   193,   193,   249,   249,   181,   244,   249,  1804,
     181,   181,   249,  1479,  1480,  1481,   249,    76,   181,   181,
     249,   249,   193,   193,    76,    76,   249,   181,    76,   191,
     193,   249,   112,   112,    76,   112,   112,   249,    76,    76,
     181,  1507,   157,   181,   249,    76,  3142,    76,   223,  1318,
     249,  1320,   249,   112,  1323,   249,   249,   125,   249,   249,
     112,   112,  1528,   249,   112,   249,   249,   249,  1337,   244,
     112,   249,  1538,  1539,   112,   112,    76,   219,   249,   249,
     204,   112,   244,   112,   223,   249,   249,  1553,   223,  1555,
     244,    76,  1558,   223,   249,   249,     7,  1563,     7,  1565,
      76,   181,   181,   244,   181,   181,   244,  1573,   249,   249,
    1576,   249,   112,   193,   193,  1581,   193,   193,  1584,   125,
     249,  1587,   181,  1589,  1590,  1591,  1592,   112,  1594,   181,
     181,   249,  1598,   181,   193,   249,   112,   249,   249,   181,
     101,   193,   193,   181,   181,   193,    76,  1613,  1614,   249,
     181,   193,   181,   181,  1620,   193,   193,  1623,   249,   181,
    3256,  1627,   193,    76,   193,  1631,   181,   181,  1634,   249,
     249,  1637,   249,   249,    76,   249,   249,  1643,   223,   249,
    1646,   181,   112,  1649,   181,   249,   249,   181,  1654,   249,
     249,   249,    76,   193,   181,   249,   181,   249,   249,   112,
    1666,   249,   249,    76,   181,   181,   249,   249,   193,   249,
     112,   249,   249,  2113,   249,   249,   244,   193,   249,   249,
     249,   249,   244,    76,   249,   181,   249,   249,   112,   244,
     244,   249,   249,   249,   249,   249,   249,  1703,    76,   112,
      76,   249,  1708,  1709,   249,  1711,  1712,   244,  1714,   249,
     244,   181,   249,    76,    76,   249,   249,   244,   249,   112,
     249,   249,   249,   193,   249,  3361,   249,   244,   181,  1735,
      76,   249,   249,   249,   112,    76,   112,  1743,   249,   181,
     193,  1747,   323,  1749,  1750,   381,    59,  1753,   244,   112,
     112,   193,  1758,   249,   249,   249,  1762,   181,   249,    76,
    1766,  1767,  1768,   249,    59,   249,   112,    59,   181,   193,
      59,   112,    76,    76,  1780,  1781,  1782,    76,    76,   249,
     193,  1787,  1788,  1789,  1790,   378,  1792,   181,   181,  1795,
    1796,  1797,    76,  1799,  1800,   112,   249,    76,   181,  1805,
     193,  1807,    76,   181,   181,   181,    59,   249,   112,   112,
      76,    59,    76,   112,   112,   193,    59,   193,   181,   181,
    1826,  1827,  1828,  1829,   249,   249,   249,  1833,   112,   249,
     193,   193,    59,   112,   249,   181,   249,  1843,   112,   327,
     181,   330,  1848,   125,  1850,   125,   112,   193,   112,    76,
     244,   372,   193,   181,   125,   249,   249,   181,   372,   249,
     125,   244,   125,   249,   181,   125,   249,   244,   249,   249,
     249,   249,   249,   249,    76,    76,   193,   181,   181,  1885,
     181,   181,   181,   181,   181,   112,   249,   249,   249,   193,
     193,   249,   249,   181,   193,   193,   181,   181,   249,   125,
     181,   125,   181,   249,   125,   125,  1912,   181,   249,   193,
     112,   112,    95,   372,   193,   181,   244,   181,  1924,   193,
     244,   249,   249,   249,    59,   249,   249,   193,   111,   193,
     249,  2266,   249,    76,   372,    76,  1942,  1943,    76,   372,
      76,    76,   372,   244,   244,   249,   249,   244,   249,   249,
     249,   249,   249,   125,   181,  1961,   244,  1963,    72,   244,
     125,   249,   249,   244,   249,   249,   193,  1973,   249,   112,
     249,   112,  1978,   249,   112,   249,   112,   112,  1984,   181,
     181,  1987,   249,   249,    76,   249,  1992,   249,   249,   380,
    1996,   193,   193,  1999,  2000,  2001,  2002,  2003,  2004,  2005,
    2006,  2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,
    2016,  2017,  2018,  2019,  2020,  2021,  2022,    76,   181,   249,
     112,  2027,  2028,  2029,  2030,  2031,  2032,  2033,  2034,  2035,
    2036,   249,   249,   249,   249,   181,   219,   220,   181,   249,
     181,  2047,   249,   181,   181,   181,   181,   249,   249,    59,
     193,    59,   193,   112,   249,   193,  2062,   193,   193,   249,
     243,   249,  2068,   249,   249,   125,    76,   249,    76,    76,
    2076,  2077,   249,   125,    76,  2081,   249,   181,    59,    76,
     181,   244,  2088,  2089,  2090,   249,   249,  2093,  2094,   181,
     249,  2097,   181,   249,   372,    76,    76,    76,   244,  2105,
     181,   193,   112,   249,   112,   112,   249,   244,   249,   249,
     112,   249,   249,   249,  2120,   112,  2122,  2123,   249,   125,
      76,  2127,   181,  2129,   181,   181,  2132,   372,    76,  2135,
     372,   112,   112,   112,   193,  2141,  2142,   320,   321,  2145,
     244,   181,  2148,   244,   514,   249,   249,  2153,   249,  2155,
     125,  2157,   249,   372,   337,   244,   112,   249,    76,    76,
     249,    76,   125,   244,   112,  2500,   181,   372,   249,   372,
     125,   181,   372,   181,   181,   372,  2182,  2183,  2184,   181,
    2186,  2187,    76,   193,   181,   193,   193,   244,   244,   511,
     249,   193,   249,   249,   112,   112,   193,   112,  2204,   249,
     181,   181,   181,   372,   244,  2211,   372,   249,  2214,   249,
     372,  2217,   193,   193,   193,   372,   125,  2223,   112,   372,
     403,   404,   405,   406,   372,   181,   125,  2233,  2234,   244,
    2236,  2237,  2238,   181,   249,  2241,  2242,   193,  2244,  2245,
    2246,  2247,   249,  2249,  2250,   193,  2252,   249,  2254,  2255,
    2256,  2257,   249,  2259,  2260,   181,   372,   372,    76,   372,
    2266,   372,  2268,   181,   181,  2271,   181,   181,  2274,   249,
     249,  2277,  2278,  2279,   181,   193,   193,  2283,   193,    76,
    2286,   125,   125,    76,  2290,   181,   372,   181,  2294,   372,
      76,   181,  2298,   249,   112,   372,   372,   372,  2633,   193,
    2306,   249,   249,   181,   249,   249,   249,   125,   125,   249,
     249,   249,   249,   249,   249,   112,   249,   223,   244,   112,
     117,   249,   505,   249,   507,    76,   112,   510,   249,   512,
     244,   249,   249,   249,   249,   249,   372,   244,   372,   249,
      59,   249,   249,   125,   249,   249,   249,   530,   244,   532,
     533,   249,   249,   249,   244,   249,   249,    76,   541,   249,
     543,   112,    76,   181,  2699,   249,   244,   249,   551,   552,
     553,   249,   555,   556,   557,   193,   249,    76,   561,   562,
     563,   564,   249,   249,   181,   568,   569,   249,   181,   572,
     573,   574,   575,   112,   577,   181,   193,  2403,   112,  2405,
     193,   249,  2408,   249,   587,   588,   589,   193,   249,   249,
    2416,  2417,  2418,   112,  2420,  2421,   249,  2423,  2424,  2425,
    2426,  2427,  2428,  2429,  2430,  2431,  2432,  2433,  2434,  2435,
     181,  2437,  2438,  2439,  2440,  2441,  2442,  2443,  2444,  2445,
      76,   181,   193,   249,  2450,  2451,  2452,  2453,  2454,  2455,
    2456,  2457,  2458,  2459,   249,    76,   249,  2463,  2464,   249,
     249,    76,   181,   249,    76,    76,    76,   181,  2474,   249,
     249,   249,    76,  2479,   193,   249,   112,   249,   249,   193,
    2486,   249,   181,  2489,   249,   181,  2492,   249,   249,   372,
     249,   112,    76,    76,   193,   181,   249,   112,   249,   249,
     112,   112,   112,  2509,   244,  2511,   249,   249,   112,   249,
     249,  2517,  2518,  2519,  2520,  2521,    76,  2523,   249,   249,
    2526,  2527,   181,  2529,  2530,   249,  2532,  2533,   112,   112,
    2536,  2866,   249,   249,  2540,  2541,   249,  2543,  2544,   249,
    2546,  2547,   249,   125,   125,   181,  2552,  2553,   244,   164,
     249,   249,   112,   249,   249,  2561,  2562,   193,   244,  2565,
     181,  2567,  2568,   249,  2570,  2571,   181,   249,    76,   181,
     181,   181,   193,   249,   249,    76,   249,   181,   193,  2585,
     249,   193,   193,   193,   249,   244,   249,   249,  2594,   193,
     249,   249,   181,   181,   181,   181,   181,   181,   181,  2605,
    2606,  2607,   249,   249,   112,  2611,  2612,  2613,   181,   193,
     193,   112,  2618,   249,  2620,  2621,  2622,   232,  2624,  2625,
    2626,   181,  2628,  2629,   249,   249,   249,  2633,   249,   249,
     249,  2637,   249,   193,  2640,  2641,   249,   249,   249,   249,
     249,    76,    76,  2649,  2650,   249,    76,    76,    76,   249,
      76,    76,    76,   249,  2660,   244,   244,   244,   244,   244,
     249,   249,   249,   249,   249,   249,   249,   372,   372,   249,
     249,   244,   249,   181,    76,   249,   249,   112,   112,   249,
     181,   181,   112,   112,   112,   193,   112,   112,   112,   249,
     249,   249,   193,   249,   249,  2701,   249,  2703,  2704,  2705,
    2706,  2707,  2708,  2709,  2710,  2711,  2712,   181,  2714,  2715,
     112,  2717,  2718,  2719,  2720,  2721,  2722,  2723,  2724,   249,
    2726,  2727,  2728,  2729,  2730,   249,   249,  2733,  2734,  2735,
    2736,  2737,  2738,  2739,  2740,  2741,  2742,  2743,  2744,  2745,
    2746,   249,  2748,   372,   244,  2751,   181,   181,   249,   249,
     249,   181,   181,   181,   249,   181,   181,   181,   193,   193,
     249,   249,   249,   193,   193,   193,   372,   193,   193,   193,
     244,   181,  2778,  2779,  2780,   249,  2782,    76,   181,   181,
    2786,   249,  2788,   249,  2790,   184,  2792,   249,  2794,   249,
    2796,   193,  2798,   249,  2800,   181,    76,   181,   249,   249,
     125,   181,    94,    94,   125,  2811,   372,  2813,  2814,   249,
    2816,  2817,   181,   112,   249,   249,   181,   125,   125,   249,
     249,   249,  2828,   249,   249,   249,   249,   249,  2834,  2835,
    2836,  2837,   112,   249,   244,  2841,  2842,  2843,   249,   249,
    2846,   244,   249,  3178,  2850,   249,   249,   249,   249,  2855,
     323,  2857,  2858,  2859,    76,  2861,  2862,   249,   244,   249,
     244,   249,   249,   249,   244,   249,  2872,  2873,   249,   249,
     249,   249,  2878,   249,   249,   244,   249,   249,   249,   244,
     249,    76,   181,   249,   249,   249,   249,    76,   249,   249,
     112,   249,   249,   249,   193,   249,   249,   181,   181,   181,
     249,   181,   249,  2909,   249,  2911,  2912,  2913,  2914,  2915,
    2916,  2917,  2918,   193,  2920,  2921,  2922,   112,  2924,  2925,
    2926,   249,  2928,   112,  2930,  2931,   249,  2933,  2934,  2935,
    2936,   249,  2938,  2939,   249,  2941,  2942,  2943,   249,  2945,
    2946,  2947,  2948,  2949,  2950,  2951,  2952,  2953,  2954,  2955,
     249,  2957,  2958,  2959,  2960,  2961,  2962,  2963,    76,   181,
     244,   244,   244,  2969,  3299,   249,   249,   249,   249,   249,
     249,   193,   249,  2979,  2980,   809,   810,   811,   812,   181,
     125,   125,   181,   181,   818,  2991,   181,   249,  2994,    76,
     181,   249,   181,   181,   112,   249,   181,   249,   193,   249,
     362,  3007,  3008,  3009,   193,  3011,   249,   249,   280,   181,
     282,   283,   284,   285,   286,   287,   288,   289,  3024,  3025,
    3026,   125,  3028,   249,   249,   112,  3032,   249,   249,   249,
     181,   249,  3038,  3039,  3040,  3041,  3042,   871,   249,   249,
    3046,  3047,   244,  3049,  3050,   244,   244,   249,   249,   249,
     249,   249,   249,   244,   249,   249,   244,   249,   249,   244,
     249,   249,   181,   181,   249,  3071,  3072,  3073,  3074,  3075,
    3076,  3077,   244,  3079,   181,   193,  3082,   249,   249,  3085,
     249,  3087,  3088,  3089,  3090,  3091,  3092,  3093,   249,  3095,
    3096,   249,  3098,   244,   181,  3101,  3102,  3103,   249,  3105,
    3106,  3107,  3108,  3109,  3110,  3111,   193,  3113,  3114,   181,
    3116,  3117,   249,  3119,   249,  3121,  3122,  3123,  3124,  3125,
    3126,  3127,   249,   249,   249,   244,  3132,  3133,   249,   181,
     249,   249,   249,  3139,  3140,  3141,   249,   244,   181,   249,
     249,  3147,   249,    76,  3150,  3151,  3152,   249,   249,    76,
     249,   249,   249,  3159,  3160,   181,  3162,  3163,  3164,   249,
     249,  3167,   249,   249,   249,   249,  3172,  3173,  3174,  3175,
    3176,   249,   244,  3179,  3180,  3181,   249,   249,   249,   112,
     249,   249,   181,   181,   249,   112,    76,  3193,  3194,   249,
    3196,  3197,   244,  3199,   181,  3201,  3202,   249,  3204,   249,
    3206,   244,  3208,   249,  3210,   249,   249,  3213,  3214,  3215,
    3216,  3217,  3218,  3219,  3220,  3221,   249,  3223,   244,  3225,
    3226,   249,   112,   249,  3230,   181,  3232,  3233,  3234,  3235,
    3236,  3237,  3238,   249,  3240,   249,  3242,   181,  3244,  3245,
    3246,  3247,  3248,  3249,  3250,   244,   244,   249,   181,  3255,
     249,   249,  3258,  3259,   181,  3261,   249,   244,  3264,  3265,
     193,  3267,   249,  3269,  3270,   249,   193,  3273,   249,  3275,
     249,    76,   249,  3279,  3280,   249,   249,   249,  3284,  3285,
    3286,   249,  3288,   249,    76,  3291,  3292,  3293,   244,   249,
     181,   181,   249,   249,   249,  3301,   249,  3303,   511,  3305,
     244,  3307,  3308,   193,  3310,   249,  3312,   112,   125,  3315,
    3316,  3317,   125,  3319,   249,  3321,   249,  3323,  3324,  3325,
     112,  3327,   249,  3329,   181,   249,  3332,   181,  3334,  3335,
    3336,  3337,  3338,  3339,  3340,   249,  3342,  3343,   249,  3345,
    3346,  3347,   249,  3349,   249,  3351,  3352,  3353,  3354,  3355,
     249,   249,   249,   244,  3360,    76,   249,    76,   249,   249,
     372,  3367,   249,  3369,  3370,  3371,   249,   249,   249,  3375,
     249,   249,   249,    76,  3380,    76,   181,   249,   249,  3385,
    3386,  3387,  3388,   181,  3390,  3391,  3392,   244,   193,   181,
     244,   112,   249,   112,  3400,   249,   249,  3403,  3404,  3405,
    3406,   193,  3408,   249,  3410,   372,   249,  3413,  3414,   112,
    3416,   112,   249,   249,    76,   249,  3422,   249,  3424,   249,
    3426,  3427,  3428,  3429,  3430,  3431,  3432,  3433,   249,  3435,
     181,  3437,   327,  3439,   249,  3441,   372,  3443,  3444,  3445,
    3446,  3447,   249,   249,   249,  3451,   244,    76,   249,  3455,
     112,   249,   249,  3459,   327,   249,   249,   249,   181,   249,
     181,  3467,   181,   249,  3470,  3471,   181,  3473,  3474,  3475,
     249,  3477,   193,   249,   193,   249,   249,  3483,   181,  3485,
     181,  3487,   249,   112,   249,   249,  3492,  3493,   249,   249,
     193,    76,   193,   244,  3500,   181,   249,  3503,   249,  3505,
    3506,  3507,  3508,  3509,  3510,  3511,    76,  3513,   181,  3515,
     249,   249,  3518,  3519,  3520,  3521,  3522,  3523,   249,   181,
    3526,   244,   181,  3529,   249,   249,   249,   112,   249,   244,
     249,   193,  3538,   249,   249,  3541,  3542,   249,    76,   249,
     249,   249,   112,   249,  3550,   249,   249,   249,   249,   249,
    3556,  3557,   181,  3559,  3560,   249,  3562,   249,   244,   249,
     249,   249,  3568,   249,   193,  3571,  3572,  3573,  3574,  3575,
    3576,   244,   249,   249,   112,  3581,   249,  3583,  3584,  3585,
     249,  3587,  3588,   125,   249,   244,   249,   249,  3594,  3595,
     249,   249,   249,   249,   249,   249,   181,   323,  3604,   249,
    3606,   249,  3608,   249,  3610,   249,  3612,   249,   193,   249,
    3616,   181,  3618,  3619,  3620,  3621,  3622,   249,   249,   249,
     249,   249,   249,   193,  3630,  3631,   249,  3633,  3634,    76,
    3636,   249,   249,   249,   249,    76,  3642,   249,  3644,   249,
     249,   249,   249,   181,   249,  3651,   249,  3653,  3654,   249,
    3656,  3657,  3658,    76,  3660,   193,   249,   249,   249,  3665,
      76,  3667,  3668,  3669,   249,   112,   249,   249,  3674,    76,
    3676,   112,  3678,   249,  3680,   249,   249,   249,  3684,   249,
    3686,   249,  3688,   249,  3690,   249,  3692,   249,  3694,   112,
     249,   249,   249,  3699,  3700,   249,   112,  3703,   249,  3705,
    3706,   249,  3708,  3709,  3710,   112,   249,  3713,   249,  3715,
    3716,   249,  3718,    76,  3720,   249,   249,   249,  3724,    76,
      76,   249,  3728,    76,  3730,   249,  3732,   249,  3734,   249,
    3736,  3737,  3738,  3739,   181,   249,   249,   249,  3744,    76,
     181,   249,  3748,   249,  3750,   249,   193,    76,  3754,   112,
     249,   249,   193,  3759,   249,   112,   112,   249,   181,   112,
    3766,  3767,   249,  3769,   249,   181,  3772,   249,  3774,   249,
     193,  3777,  3778,   125,   181,   112,  3782,   193,  3784,   249,
    3786,   249,  3788,   112,     4,     5,   193,   249,  3794,  3795,
    3796,    11,    12,   249,   249,    15,   249,    17,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,    32,    33,   249,   249,   249,   249,   181,   249,
      40,    41,    42,    43,   181,   181,   249,   249,   181,   249,
     193,   249,    52,   249,    54,   249,   193,   193,   249,   249,
     193,   125,   249,   249,   181,   249,   249,   249,    68,   249,
      70,   249,   181,    73,   249,   249,   193,   249,    78,   249,
     249,    81,   249,   249,   193,    85,    86,    87,   249,   249,
     249,    91,   249,   249,    94,    95,   249,   249,    98,    99,
     249,   249,   102,   249,   104,   105,   249,   107,   108,   249,
     249,   249,   249,   249,   249,   249,   249,   117,   249,   249,
     120,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   144,   145,   249,   249,   148,   149,
     150,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   172,   177,   595,   597,
    1051,   360,   182,   359,   625,  1047,   186,  1014,   188,   189,
     190,   935,   377,   476,   963,   946,   407,   438,   198,   304,
     304,   528,  1011,   304,   494,  1224,   415,  1563,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   319,
      -1,    -1,    -1,    -1,    -1,    -1,   226,    -1,   228,   229,
      -1,    -1,   232,   233,   234,   235,    -1,   237,    -1,   239,
     240,   241,   242,    -1,    -1,   245,   246,   247,    -1,    -1,
     250,    -1,   252,   253,   254,    -1,    -1,    -1,   258,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   273,   274,   275,    -1,    -1,   278,   279,
     280,   281,   282,    -1,   284,   285,   286,   287,   288,   289,
     290,    -1,    -1,   293,   294,   295,    -1,    -1,   298,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   312,   313,   314,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   323,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   339,
      -1,    -1,   342,   343,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   356,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   377,   378,   379,
      -1,    -1,   382,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   396,   397,    -1,   399,
      -1,    -1,    -1,   403,    -1,    -1,   406,   407,    -1,    -1,
     410,   411,    -1,    -1,    -1,   415,   416,    -1,    -1,   419,
      -1,    -1,   422,    -1,    -1,    -1,   426,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   438,    -1,
      -1,    -1,    -1,    -1,    -1,   445,    -1,   447,    -1,    -1,
      -1,    -1,   452,    -1,    -1,    -1,   456,    -1,    -1,    -1,
     460,    -1,    -1,   463,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   474,   475,   476,    -1,    -1,    -1,
      -1,    -1,    -1,   483,    -1,    -1,    -1,    -1,    -1,   489,
      -1,    -1,    -1,    -1,   494,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   503,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   512,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     520,   521,    -1,    -1,    -1,    -1,    -1,   527,   528,   529,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   545,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   573,    -1,   575,    -1,    -1,    -1,    -1,
      -1,   581,   582,   583,   584,    -1,    -1,    -1,   588,   589,
     590,   591,   592,    -1,   594,   595,   596,   597,    -1,    -1,
      -1,    -1,    -1,   603,    -1,   605,   606,    -1,   608,   609,
     610,    -1,    -1,    -1,    -1,    -1,   616,   617,   618,    -1,
      -1,    -1,    -1,    -1,   624,    -1,    -1,    -1,    -1,    -1,
     630,    -1,    -1,    -1,   634,    -1,    -1,    -1,    -1,   639,
      -1,    -1,    -1,    -1,    -1,   645,    -1,   647,   648,   649,
      -1,    -1,   652,    -1,   654,    -1,   656,    -1,   658,    -1,
     660,    -1,    -1,    -1,   664,    -1,   666,    -1,   668,    -1,
     670,    -1,   672,    -1,    -1,    -1,    -1,    -1,    -1,   679,
      -1,    -1,    -1,   683,    -1,   685,    -1,   687,    -1,    -1,
      -1,    -1,   692,   693,   694,   695,    -1,   697,   698,    -1,
      -1,   701,    -1,    -1,   704,    -1,    -1,    -1,    -1,    -1,
      -1,   711,    -1,    -1,   714,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   722,    -1,    -1,    -1,    -1,   727,    -1,    -1,
      -1,    -1,    -1,    -1,   734,   735,    -1,    -1,    -1,   739,
      -1,    -1,    -1,    -1,    -1,   745,    -1,   747,   748,   749,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   757,   758,    -1,
      -1,    -1,    -1,   763,   764,    -1,   766,   767,   768,   769,
      -1,   771,    -1,   773,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   782,   783,    -1,    -1,   786,   787,    -1,    -1,
      -1,   791,   792,    -1,   794,    -1,   796,    -1,    -1,   799,
     800,   801,   802,    -1,   804,   805,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   814,    -1,   816,    -1,    -1,    -1,
     820,   821,    -1,   823,    -1,    -1,    -1,    -1,    -1,   829,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   838,   839,
      -1,    -1,    -1,    -1,   844,   845,    -1,    -1,    -1,   849,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   861,   862,   863,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   874,    -1,   876,    -1,   878,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   889,
      -1,   891,    -1,   893,   894,    -1,    -1,    -1,    -1,   899,
      -1,    -1,    -1,    -1,   904,    -1,    -1,    -1,    -1,    -1,
     910,    -1,    -1,    -1,    -1,    -1,   916,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   925,    -1,    -1,    -1,    -1,
     930,    -1,   932,   933,   934,   935,    -1,   937,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   946,   947,    -1,   949,
      -1,    -1,   952,   953,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   962,   963,   964,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     980,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     990,    -1,    -1,   993,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1006,    -1,  1008,    -1,
    1010,    -1,    -1,    -1,  1014,    -1,  1016,    -1,    -1,    -1,
      -1,  1021,    -1,    -1,  1024,  1025,    -1,  1027,    -1,  1029,
      -1,  1031,    -1,  1033,    -1,  1035,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1047,    -1,    -1,
      -1,  1051,    -1,  1053,    -1,    -1,  1056,  1057,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1068,    -1,
    1070,   115,   116,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   125,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1097,  1098,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1109,
      -1,  1111,    -1,    -1,    -1,  1115,    -1,    -1,    -1,    -1,
    1120,  1121,    -1,  1123,  1124,  1125,  1126,    -1,    -1,    -1,
      -1,    -1,    -1,  1133,    -1,  1135,    -1,  1137,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1166,    -1,    -1,  1169,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1186,  1187,    -1,  1189,
      -1,    -1,    -1,    -1,    -1,    -1,  1196,  1197,  1198,  1199,
    1200,  1201,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1213,    -1,    -1,    -1,  1217,  1218,    -1,
      -1,    -1,    -1,  1223,  1224,  1225,  1226,  1227,  1228,    -1,
    1230,    -1,  1232,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1264,    -1,    -1,    -1,  1268,    -1,
      -1,    -1,    -1,  1273,    -1,  1275,    -1,    -1,  1278,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1286,    -1,  1288,    -1,
    1290,    -1,    -1,    -1,  1294,    -1,    -1,   341,   342,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1327,    -1,    -1,
    1330,    -1,    -1,    -1,    -1,    -1,  1336,    -1,  1338,  1339,
    1340,  1341,  1342,  1343,    -1,  1345,  1346,  1347,  1348,  1349,
    1350,    -1,  1352,  1353,    -1,  1355,  1356,  1357,  1358,    -1,
      -1,    -1,  1362,    -1,  1364,    -1,  1366,  1367,  1368,    -1,
    1370,  1371,  1372,  1373,  1374,  1375,    -1,  1377,    -1,  1379,
      -1,    -1,    -1,    -1,  1384,  1385,  1386,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1395,    -1,    -1,    -1,    -1,
    1400,    -1,    -1,    -1,    -1,    -1,  1406,    -1,    -1,   453,
     454,   455,   456,   457,    -1,   459,    -1,   461,   462,   463,
     464,   465,   466,   467,   468,    -1,   470,   471,   472,   473,
     474,   475,   476,   477,   478,   479,   480,   481,   482,   483,
      -1,   485,   486,   487,   488,   489,   490,   491,   492,   493,
     494,    -1,    -1,    -1,    -1,    -1,  1456,    -1,    -1,    -1,
      -1,    -1,  1462,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1479,
      95,    -1,    -1,    -1,  1484,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   111,    -1,    -1,    -1,
    1500,    -1,    -1,  1503,    -1,    -1,    -1,  1507,    -1,    -1,
      -1,    -1,  1512,  1513,    -1,    -1,    -1,    -1,    -1,    -1,
    1520,  1521,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1531,    -1,    -1,    -1,  1535,  1536,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1546,    -1,    -1,    -1,
    1550,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1568,    -1,
    1570,    -1,    -1,    -1,  1574,    -1,  1576,    -1,    -1,  1579,
      -1,  1581,    -1,  1583,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1596,    -1,  1598,    -1,
    1600,    -1,    -1,    -1,   219,   220,    -1,    -1,    -1,    -1,
    1610,  1611,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1619,
      -1,    -1,  1622,    -1,    -1,    -1,    -1,  1627,   243,    -1,
    1630,    -1,    -1,  1633,    -1,    -1,  1636,    -1,    -1,  1639,
      -1,    -1,  1642,    -1,    -1,  1645,    -1,    -1,  1648,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1662,    -1,  1664,  1665,    -1,    -1,    -1,  1669,
    1670,  1671,  1672,  1673,    -1,    -1,  1676,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1701,    -1,    -1,    -1,   320,   321,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1720,    -1,   337,    -1,  1724,    -1,  1726,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1755,    -1,    -1,    -1,  1759,
    1760,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1769,
      -1,  1771,  1772,  1773,  1774,    -1,  1776,  1777,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1786,    -1,   403,   404,
     405,   406,    -1,  1793,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1801,    -1,  1803,  1804,    -1,  1806,    -1,    -1,    -1,
    1810,  1811,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1834,    -1,    -1,    -1,    -1,    -1,
      -1,  1841,    -1,    -1,  1844,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1861,  1862,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1874,  1875,    -1,  1877,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1888,    -1,
     505,    -1,   507,  1893,    -1,   510,    -1,   512,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   530,  1916,   532,   533,    -1,
      -1,    -1,    -1,  1923,    -1,    -1,   541,    -1,   543,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     555,   556,   557,    -1,  1944,    -1,   561,   562,   563,   564,
      -1,    -1,    -1,   568,   569,    -1,    -1,   572,   573,   574,
     575,    -1,   577,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   587,   588,   589,    -1,    -1,  1977,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1989,
      -1,    -1,    -1,    -1,  1994,  1995,    -1,  1997,  1998,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2023,  2024,  2025,  2026,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2045,    -1,    -1,  2048,    -1,
      -1,    -1,    -1,  2053,  2054,    -1,    -1,  2057,  2058,    -1,
      -1,    -1,    -1,    -1,    -1,  2065,  2066,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2090,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2107,    -1,    -1,
      -1,    -1,    -1,  2113,    -1,    -1,    -1,    -1,  2118,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2127,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2159,
      -1,    -1,    -1,    -1,    -1,    -1,  2166,    -1,    -1,    -1,
      -1,    -1,    -1,  2173,    -1,    -1,    -1,    -1,    -1,  2179,
    2180,  2181,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2193,    -1,    -1,    -1,    -1,  2198,    -1,
    2200,  2201,    -1,  2203,    -1,    -1,    -1,    -1,    -1,    -1,
    2210,    -1,    -1,  2213,    -1,    -1,  2216,    -1,    -1,    -1,
    2220,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2240,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2253,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2276,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2287,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2300,  2301,    -1,    -1,    -1,    -1,    -1,    -1,  2308,    -1,
      -1,    -1,    -1,    -1,    -1,  2315,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2337,    -1,  2339,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2399,
    2400,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2410,  2411,    -1,    -1,  2414,  2415,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2446,  2447,  2448,  2449,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2462,    -1,    -1,    -1,    -1,  2467,    -1,    -1,
      -1,    -1,  2472,  2473,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2505,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2564,    -1,  2566,    -1,    -1,    -1,
      -1,    -1,    -1,  2573,    -1,    -1,    -1,  2577,    -1,  2579,
      -1,    -1,    -1,  2583,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2591,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2608,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    16,  2639,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    29,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2662,    -1,  2664,    43,    -1,  2667,    -1,  2669,
      48,    -1,    -1,    -1,    -1,    -1,    -1,  2677,    -1,    57,
      -1,  2681,    60,    61,    -1,    -1,    -1,    65,    66,    -1,
      -1,    -1,    -1,    -1,  2694,    -1,  2696,    -1,    -1,    -1,
    2700,    -1,  2702,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2731,  2732,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   121,    -1,    -1,    -1,    -1,    -1,    -1,
    2750,    -1,  2752,    -1,    -1,    -1,  2756,    -1,    -1,    -1,
    2760,    -1,    -1,    -1,   142,   143,   144,   145,    -1,    -1,
      -1,    -1,    -1,   151,    -1,    -1,    -1,   155,    -1,    -1,
     158,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   180,    -1,   182,    -1,    -1,    -1,  2808,  2809,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,  2819,
      -1,  2821,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   210,   211,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   219,    -1,   221,    -1,    -1,   224,    -1,   226,    -1,
     228,   229,    -1,    -1,    -1,   233,    -1,   235,   236,   237,
     238,    -1,    -1,    -1,    -1,    -1,  2866,    -1,    -1,  2869,
      -1,    -1,    -1,    -1,    -1,   253,    -1,   255,    -1,    -1,
      -1,    -1,    -1,   261,    -1,    -1,    -1,    -1,  2888,    -1,
    2890,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2908,    -1,
      -1,    -1,   290,    -1,    -1,    -1,    -1,   295,   296,   297,
      -1,   299,    -1,   301,    -1,   303,    -1,    -1,    -1,    -1,
     308,   309,   310,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     328,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   337,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2967,    -1,    -1,
      -1,   349,   350,   351,   352,   353,   354,   355,    -1,   357,
     358,   359,   360,    -1,    -1,    -1,    -1,    -1,    -1,   367,
     368,    -1,   370,    -1,    -1,    -1,   374,    -1,    -1,    -1,
      -1,    -1,    -1,  3003,  3004,    -1,    -1,    -1,    -1,   387,
      -1,   389,    -1,    -1,    -1,    -1,   394,    -1,    -1,    -1,
      -1,   399,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   419,    -1,    -1,    -1,   423,   424,   425,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3058,    -1,
    3060,    -1,    -1,    -1,    -1,   443,    -1,    -1,   446,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3146,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3186,    -1,  3188,    -1,
      -1,  3191,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    3260,    -1,  3262,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  3283,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3296,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3365,  3366,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3412,    -1,    -1,    -1,    -1,    -1,  3418,    -1,
    3420,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  3454,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  3491,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3535,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3546,    -1,  3548,    -1,
      -1,    -1,    -1,    -1,  3554,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    3570,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3592,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  3603,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3626,    -1,  3628,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  3641,    -1,    -1,    -1,    -1,     3,    -1,     5,     6,
    3650,    -1,    -1,    -1,    11,    -1,    13,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    21,    -1,    23,    24,    25,    26,
      27,    -1,    -1,    -1,    31,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    40,    41,    42,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    50,    -1,    52,    53,    54,    -1,    56,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    64,    -1,    -1,
      -1,    -1,    -1,    70,    71,    -1,    -1,    74,    -1,    -1,
      -1,    78,    79,    80,    81,    -1,    -1,    -1,    85,    -1,
      87,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    99,    -1,    -1,   102,   103,   104,   105,   106,
     107,    -1,   109,   110,    -1,  3755,    -1,    -1,    -1,    -1,
      -1,   118,    -1,    -1,    -1,   122,   123,    -1,    -1,    -1,
     127,   128,    -1,    -1,   131,   132,    -1,    -1,    -1,  3779,
      -1,    -1,    -1,    -1,   141,    -1,    -1,    -1,    -1,    -1,
      -1,  3791,   149,  3793,    -1,   152,   153,    -1,    -1,   156,
      -1,    -1,   159,    -1,    -1,   162,   163,   164,   165,   166,
     167,   168,   169,   170,   171,   172,   173,    -1,    -1,    -1,
      -1,   178,   179,    -1,    -1,    -1,   183,    -1,   185,    -1,
     187,   188,   189,   190,    -1,   192,    -1,   194,    -1,    -1,
     197,    -1,    -1,    -1,    -1,   202,    -1,    -1,   205,   206,
     207,   208,   209,    -1,    -1,   212,    -1,    -1,   215,   216,
     217,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   232,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   241,    -1,    -1,    -1,    -1,    -1,
     247,    -1,    -1,    -1,   251,    -1,    -1,    -1,    -1,   256,
      -1,   258,   259,    -1,    -1,    -1,    -1,    -1,   265,    -1,
     267,    -1,    -1,   270,   271,    -1,    -1,   274,    -1,    -1,
     277,   278,   279,    -1,   281,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   292,   293,   294,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   306,
     307,    -1,    -1,    -1,    -1,   312,   313,   314,    -1,   316,
      -1,    -1,    -1,    -1,    -1,   322,    -1,   324,   325,    -1,
      -1,    -1,   329,    -1,    -1,    -1,    -1,    -1,    -1,   336,
      -1,    -1,   339,   340,    -1,    -1,    -1,   344,    -1,    -1,
     347,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   361,    -1,    -1,   364,    -1,    -1,
      -1,    -1,   369,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   384,    -1,   386,
      -1,   388,    -1,    -1,   391,   392,   393,    -1,   395,   396,
      -1,   398,    -1,    -1,   401,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   411,    -1,    -1,    -1,    -1,   416,
     417,   418,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   426,
     427,   428,   429,   430,   431,    -1,   433,    -1,    -1,    -1,
      -1,    -1,   439,   440,    -1,    -1,    -1,    -1,    -1,    -1,
     447,    -1,    -1,    -1,   451,   452,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   495,   496,
      -1,    -1,    -1,    -1,    -1,   502,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   515,    -1,
     517,   518,   519,   520,   521,   522,   523,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   535,   536,
      -1,    -1,   539,   540,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   554,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   565,   566,
     567,    -1,    -1,   570,   571,    -1,     3,    -1,     5,     6,
      -1,    -1,    -1,    -1,    11,   582,    13,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    21,    -1,    23,    24,    25,    26,
      27,    -1,    -1,    -1,    31,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    40,    41,    42,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    50,    -1,    52,    53,    54,    -1,    56,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    64,    -1,    -1,
      -1,    -1,    -1,    70,    71,    -1,    -1,    74,    -1,    -1,
      -1,    78,    79,    80,    81,    -1,    -1,    -1,    85,    -1,
      87,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    99,    -1,    -1,   102,   103,   104,   105,    -1,
     107,    -1,   109,   110,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   118,    -1,    -1,    -1,   122,   123,    -1,    -1,    -1,
     127,   128,    -1,    -1,   131,   132,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   141,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   149,    -1,    -1,   152,   153,    -1,    -1,   156,
      -1,    -1,   159,    -1,    -1,   162,   163,   164,   165,   166,
     167,   168,   169,   170,   171,   172,   173,    -1,    -1,    -1,
      -1,   178,   179,    -1,    -1,    -1,   183,    -1,   185,    -1,
     187,   188,   189,   190,    -1,   192,    -1,   194,    -1,    -1,
     197,    -1,    -1,    -1,    -1,   202,    -1,    -1,   205,   206,
     207,   208,   209,    -1,    -1,   212,    -1,    -1,   215,   216,
     217,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   232,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   241,    -1,    -1,    -1,    -1,    -1,
     247,    -1,    -1,    -1,   251,    -1,    -1,    -1,    -1,   256,
      -1,   258,   259,    -1,    -1,    -1,    -1,    -1,   265,    -1,
     267,    -1,    -1,   270,   271,    -1,    -1,   274,    -1,    -1,
     277,   278,   279,    -1,   281,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   292,   293,   294,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   306,
     307,    -1,    -1,    -1,    -1,   312,   313,   314,    -1,   316,
      -1,    -1,    -1,    -1,    -1,   322,    -1,   324,   325,    -1,
      -1,    -1,   329,    -1,    -1,    -1,    -1,    -1,    -1,   336,
      -1,    -1,   339,   340,    -1,    -1,    -1,   344,    -1,    -1,
     347,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   361,    -1,    -1,   364,    -1,    -1,
      -1,    -1,   369,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   384,    -1,   386,
      -1,   388,    -1,    -1,   391,   392,   393,    -1,   395,   396,
      -1,   398,    -1,    -1,   401,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   411,    -1,    -1,    -1,    -1,   416,
     417,   418,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   426,
     427,   428,   429,   430,   431,    -1,   433,    -1,    -1,    -1,
      -1,    -1,   439,   440,    -1,    -1,    -1,    -1,    -1,    -1,
     447,    -1,    -1,    -1,   451,   452,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   495,   496,
      -1,    -1,    -1,    -1,    -1,   502,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   515,    -1,
     517,   518,   519,   520,   521,   522,   523,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   535,   536,
      -1,    -1,   539,   540,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   554,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   565,   566,
     567,    -1,    -1,   570,   571,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   582
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     3,     5,     6,    11,    13,    21,    23,    24,    25,
      26,    27,    31,    40,    41,    42,    50,    52,    53,    54,
      56,    64,    70,    71,    74,    78,    79,    80,    81,    85,
      87,    99,   102,   103,   104,   105,   107,   109,   110,   118,
     122,   123,   127,   128,   131,   132,   141,   149,   152,   153,
     156,   159,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   172,   173,   178,   179,   183,   185,   187,   188,
     189,   190,   192,   194,   197,   202,   205,   206,   207,   208,
     209,   212,   215,   216,   217,   232,   241,   247,   251,   256,
     258,   259,   265,   267,   270,   271,   274,   277,   278,   279,
     281,   292,   293,   294,   306,   307,   312,   313,   314,   316,
     322,   324,   325,   329,   336,   339,   340,   344,   347,   361,
     364,   369,   384,   386,   388,   391,   392,   393,   395,   396,
     398,   401,   411,   416,   417,   418,   426,   427,   428,   429,
     430,   431,   433,   439,   440,   447,   451,   452,   495,   496,
     502,   515,   517,   518,   519,   520,   521,   522,   523,   524,
     535,   536,   539,   540,   554,   565,   566,   567,   570,   571,
     582,   592,   593,   594,   595,   596,   597,   598,   599,   605,
     606,   607,   608,   609,   610,   611,   612,   613,   614,   623,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     636,   637,   638,   640,   641,   642,   643,   645,   646,   651,
     654,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     665,   667,   668,   670,   674,   675,   676,   678,   679,   680,
     681,   684,   687,   688,   689,   690,   691,   694,   697,   700,
     702,   704,   706,   710,   711,   712,   713,   714,   715,   716,
     718,   719,   720,   721,   722,   723,   724,   725,   726,   727,
     731,   733,   735,   737,   739,   741,   743,   745,   747,   749,
     750,   754,   758,   760,   763,   764,   765,   766,   767,   768,
     769,   771,   772,   774,   775,   783,   784,   786,   788,   789,
     790,   791,   794,   795,   796,   797,   798,   799,   800,   801,
     805,   808,   810,   812,   814,   817,   820,   822,   824,   826,
     828,   829,   830,   831,   832,   835,   836,   838,   839,   840,
     841,   844,   847,   249,     7,   249,   249,   181,   244,   849,
     249,   849,    76,   112,   181,   193,   850,   850,   850,   125,
     249,   850,   249,   249,   849,   849,   249,   249,   849,   249,
     249,   849,   249,   249,   249,   249,   440,   249,   249,   249,
     249,   249,   249,   249,   249,   249,    38,    82,   100,   101,
     249,   813,   249,   849,   249,   849,   249,   249,   249,   249,
     249,   249,   849,   849,   249,   849,   249,   849,   249,   372,
     249,   249,   249,   249,   850,   249,   249,   249,   849,   249,
     249,   849,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   432,   249,   850,   249,   249,   849,   249,   249,
     849,   249,   249,   849,   249,   249,   249,   850,   125,   249,
     849,   249,   249,   249,   849,   125,   249,   372,   249,   249,
     849,   249,   849,   849,   249,   249,   850,   249,   849,   249,
     249,   125,   195,   249,   849,   125,   195,   249,   849,   249,
     249,   849,   849,   249,   249,   249,   849,   125,   849,   849,
     249,   249,   849,   249,   849,   125,   249,   249,   249,   249,
     249,   249,   249,   249,   849,   850,   249,   249,   849,   249,
     249,   850,   850,   249,   249,   249,   249,   249,   249,   249,
     249,   125,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   440,   249,   849,   249,   849,   249,   249,   849,
     849,   849,   249,   249,   850,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   850,   249,   249,   249,   249,
     249,   372,     0,   106,   594,   157,   849,    75,   135,   136,
     137,   138,   139,   302,   305,   326,   448,   600,   601,   602,
     603,   604,    44,    45,   434,   435,   436,   437,   438,    28,
      86,    88,   140,   413,   422,   503,   504,   505,   508,   509,
     849,    22,   133,   160,   239,   850,   850,   849,    49,   198,
     199,   200,   201,   615,   617,   618,   619,   620,   785,   849,
     624,   849,   223,   849,   850,   125,   378,   779,   849,   362,
     410,   583,   584,   585,   586,   635,   219,   639,   850,    22,
      47,    75,   160,   176,   223,   239,   240,   248,   263,   276,
     373,   383,   385,   432,   644,   648,   669,   673,   160,   223,
     239,   652,   653,    32,    51,   363,   437,   655,   204,   666,
     677,   849,   378,   779,   849,   849,   378,   849,   378,   849,
     378,   849,   849,   695,   849,   701,   849,   703,   849,   705,
     849,   707,   849,   695,   703,   707,     4,   717,   849,   223,
     849,   849,   849,   223,   223,   378,   779,   849,   850,   761,
     849,   777,   849,   849,   849,   849,   770,   849,   849,   413,
     773,   849,     8,    18,    30,    33,    34,    35,    36,    37,
      77,   146,   150,   245,   254,   257,   262,   263,   264,   266,
     269,   273,   319,   343,   348,   484,   497,   498,   776,   849,
     785,   785,   787,   849,   849,   849,   111,   321,   849,    19,
     404,   792,   793,   850,   181,   849,   850,   849,   378,   849,
     849,   213,   214,   268,   317,   390,   397,   415,   579,   580,
     581,    39,    63,    69,    82,    89,   119,   161,   184,   234,
     311,   345,   346,   400,   807,   809,   601,   604,   673,   850,
      12,    83,   101,   124,   147,   203,   211,   219,   250,   252,
     260,   318,   414,   420,   456,   827,   849,   460,   849,   849,
      46,    92,    93,   406,   503,   505,   506,   513,   525,   526,
     527,   528,   529,   531,   541,   542,   544,   546,   547,   551,
     557,   558,   559,   560,   837,   843,    95,   111,   219,   220,
     243,   320,   321,   337,   403,   404,   405,   406,   505,   507,
     510,   512,   530,   532,   533,   541,   543,   555,   556,   557,
     561,   562,   563,   564,   568,   569,   572,   573,   574,   575,
     577,   587,   588,   589,   842,   842,   552,   553,   842,   843,
     845,   846,   842,   537,   538,   848,   728,   779,   849,   249,
       7,     7,   249,   249,   249,   249,   850,   249,   849,   850,
     708,   709,   849,   709,   249,   249,   249,   850,   249,   850,
     125,   849,   647,   850,   639,   850,   850,   101,   850,   813,
     175,   249,   249,   249,   701,   695,   730,   781,   849,   811,
     850,   249,   849,   249,   249,   223,   249,   249,   850,   850,
     759,   782,   849,   759,   249,   692,   693,   849,   249,   703,
     821,   850,   823,   850,   707,   705,   755,   756,   849,   728,
     249,   249,   730,   728,   249,   850,   849,   249,   850,   849,
     249,   850,   223,   751,   752,   849,   850,   249,   152,   249,
     384,   249,   249,   249,   761,   249,   249,   249,   776,   850,
     849,   249,   249,   249,   849,   249,   249,   249,   849,   249,
     849,   249,   249,   849,   249,   249,   504,   249,   249,   249,
     249,   849,   249,   849,   717,   323,   381,   748,    59,   746,
     728,   249,   249,   850,   698,   699,   849,   249,   249,   777,
     850,    59,   734,   249,   728,    59,   732,    59,   736,   378,
     738,    59,   744,    59,   740,    59,   742,   849,   249,   249,
     249,   332,   377,   849,   249,   849,   249,   682,   683,   849,
     770,   685,   686,   849,   850,   249,   378,   849,   330,   850,
     850,   850,   850,   850,   850,   850,   249,   850,   327,   850,
     850,   372,   372,   125,   125,   125,   125,   125,   249,   249,
     125,   249,   249,   849,   249,   849,   249,   249,   249,   249,
     850,   849,   849,   849,   849,   850,   850,   164,   232,   850,
     849,   849,   849,   849,   849,   621,   849,   621,   622,   849,
     622,   408,   850,   850,   249,   849,   850,   125,   849,   779,
     849,   412,   849,   125,   849,   125,   849,   125,   125,   125,
     125,    14,    20,    67,   218,   246,   249,   256,   272,   334,
     335,   365,   366,   456,   850,   849,   249,   849,   849,   649,
     850,   249,   850,   850,   650,   850,   372,   249,   849,   649,
     249,   249,   249,   372,   372,   849,   372,   850,   372,    55,
     249,   849,   850,   249,   849,   850,   249,   849,   125,   249,
     849,   849,   249,   849,   850,   779,   412,    72,   849,   412,
     850,   849,   850,   849,   850,   850,   849,   849,   849,   849,
     849,   125,   103,   333,   249,   849,   850,   850,   850,   249,
     849,   249,   779,   412,   249,    10,    15,    68,   114,   215,
     379,   409,   850,   849,   849,   849,   849,   849,   849,   249,
     440,   813,   849,   249,   249,   849,   850,   850,   850,   850,
     850,   249,   849,   249,   249,   849,   249,   249,   249,   249,
     850,   249,   249,   849,   249,   850,   249,   850,   849,   850,
     148,   850,   849,   407,   849,   249,   249,   850,   849,   125,
     249,   849,   249,   849,   849,   117,   849,   850,   412,   850,
     850,   850,   380,   125,   850,   850,   804,   849,   803,   849,
     249,   850,   249,   850,   806,   849,   249,   849,   372,    84,
     249,   449,   590,   849,   249,   849,   249,   849,   849,   249,
     850,   849,   249,   298,   371,   420,   849,    16,    29,    43,
      48,    57,    60,    61,    65,    66,   121,   142,   143,   144,
     145,   151,   155,   158,   180,   182,   196,   210,   211,   219,
     221,   224,   226,   228,   229,   233,   235,   236,   237,   238,
     253,   255,   261,   290,   295,   296,   297,   299,   301,   303,
     308,   309,   310,   328,   337,   349,   350,   351,   352,   353,
     354,   355,   357,   358,   359,   360,   367,   368,   370,   374,
     387,   389,   394,   399,   419,   423,   424,   425,   443,   446,
     516,   850,   249,   850,   849,   849,   249,   850,   850,     9,
     825,   849,   849,   249,   249,   850,   849,   849,   249,   249,
     849,   849,   125,   115,   116,   125,   341,   342,   453,   454,
     455,   456,   457,   459,   461,   462,   463,   464,   465,   466,
     467,   468,   470,   471,   472,   473,   474,   475,   476,   477,
     478,   479,   480,   481,   482,   483,   485,   486,   487,   488,
     489,   490,   491,   492,   493,   494,   849,   849,   849,   849,
     372,   849,   849,   372,   514,   834,   834,   834,   834,   834,
     850,   849,   125,   126,   372,   849,   372,   834,   249,   849,
     849,   849,   249,   647,   125,   372,   850,   849,   372,   372,
     372,   372,   372,   372,   372,   372,   849,   849,   511,   850,
     125,   125,   849,   849,   850,   372,   372,   849,   372,   125,
     372,   372,   125,   125,   372,   372,   372,   372,   372,   849,
     849,   849,   249,   249,   249,   249,   249,   125,   850,   834,
     249,   125,   849,   249,   779,   849,   412,   249,   850,   850,
     850,   249,   249,   849,   708,   850,   709,   249,   850,   849,
     249,   249,   249,   850,   249,   850,   850,   249,   849,   781,
     850,   249,   249,   818,   819,   850,   249,   849,   223,   249,
     850,   850,   782,   849,   759,   693,   849,   249,   249,   756,
     757,   849,   757,   849,    62,   249,   850,   815,   816,   850,
     850,   850,   850,   249,   850,   849,   752,   753,   849,   753,
     849,   249,   850,   850,   372,   372,   850,   849,   249,   249,
     849,   849,   249,   815,   815,   249,   249,   249,   849,   381,
     850,   849,    59,   850,   809,   850,   699,   849,   249,   849,
      59,   850,   849,    59,   850,   849,    59,   850,   849,   378,
     850,   849,    59,   850,   849,    59,   850,   849,    59,   850,
     249,   249,   850,   249,   499,   500,   501,   249,   683,   249,
     850,   686,   849,   850,   849,   849,   850,   850,   249,   850,
     850,   850,   850,   850,   249,   850,   849,   249,   849,   125,
     249,   125,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   849,   849,
     249,   249,   372,   850,   249,   249,   249,   249,   850,   849,
     849,   850,   850,   125,   849,   850,   850,   125,   249,   249,
     849,   242,   249,   850,   849,   849,   849,   849,   849,   849,
     849,   372,   850,   249,   249,   850,   850,   850,   249,   249,
     249,   729,   780,   849,   249,   249,   249,   850,   249,   849,
     850,   249,   849,   850,   249,   850,   249,   249,   850,   849,
     849,   850,   849,   249,   850,   249,   850,   850,   850,   696,
     849,   696,   696,   696,   696,   849,   125,   125,   729,   249,
     850,   850,   850,   729,   249,   729,   849,   850,   850,   113,
     215,   249,   379,   440,   813,   849,   849,   850,   850,   850,
     850,   778,   849,   778,   849,   849,   696,   849,   849,   249,
     440,   849,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   850,   850,   850,   850,
     850,   177,   249,   407,   849,   249,   849,   249,   249,   249,
     849,   412,   850,   849,   117,   850,   249,   372,   849,   249,
     850,   249,   249,   249,   249,   249,   249,   249,   849,   249,
     849,   802,   849,   249,   249,   249,   849,   249,   249,   249,
     249,   249,   120,   249,   850,   850,   249,   850,   249,   850,
     249,   249,   249,   249,   249,   850,    17,   227,   444,   807,
     372,   807,    17,    58,   850,   807,   372,   372,   249,   338,
     849,   249,   249,   849,   850,   850,   850,   850,   249,   191,
     849,   807,   849,   849,   849,   225,   849,   227,   849,   338,
     849,   249,   849,   849,   849,   849,   849,   249,   849,   249,
     849,   121,   212,   849,   125,   249,   849,   849,   300,   849,
     849,   850,   850,   850,   849,   372,   338,   849,   850,   849,
     849,   849,   850,   849,   849,   849,   849,   849,   849,   249,
     849,   850,   849,   850,   850,   850,   850,   849,   249,   849,
     849,   249,   184,   850,   249,   372,   850,   850,   850,   249,
     849,   249,   249,   249,   849,   249,   249,   850,   249,   849,
     249,   249,   850,   125,   850,   850,   833,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   249,   849,    94,
      94,   849,   850,   125,   372,   849,   850,   850,   850,   249,
     849,   849,   125,   849,   849,   850,   125,   849,   849,   849,
     849,   249,   850,   249,   849,   412,   849,   249,   850,   850,
     249,   249,   249,   249,   249,   849,   850,   850,   850,   249,
     819,   850,   729,   249,   249,   849,   249,   849,   249,   850,
     849,   850,   757,   849,   850,   249,   816,   850,   249,   850,
     850,   850,   850,   249,   753,   849,   249,   440,   813,   249,
     249,   384,   249,   249,   249,   249,   671,   849,   849,   323,
     249,   849,   850,   249,   849,   850,   249,   849,   850,   249,
     849,   850,   249,   849,   850,   249,   849,   850,   108,   249,
     849,   833,   249,   849,   850,   249,   849,   850,   249,   849,
     850,   249,   332,   850,   249,   499,   249,   499,   249,   696,
     249,   849,   249,   849,   850,   249,   849,   849,   849,   849,
     849,   249,   249,   849,   249,   249,   249,   249,   249,   616,
     849,   249,   850,   850,   850,   249,   850,   850,   249,   850,
     249,   249,   249,   375,   849,   249,   242,   249,   157,   849,
     125,   157,   849,   125,   850,   249,   242,   249,   780,   850,
     850,   850,   850,   850,   849,   850,   375,   849,   849,   249,
     375,   850,   249,   850,   850,   850,   249,   849,   249,   249,
     249,   249,   362,   849,   849,   729,   850,   850,   850,   729,
     375,   849,   850,   249,   850,   850,   850,   850,   849,   113,
     215,   249,   379,   440,   850,   850,   850,   134,   249,   850,
     850,   249,   849,   249,   249,   813,   849,   249,   850,   249,
     249,   850,   249,   849,   849,   249,   850,   850,   850,   850,
     111,   249,   249,   321,   850,   249,   402,   849,   125,   249,
     849,   249,   372,   850,   849,   249,   372,   249,   850,   850,
     849,   849,   249,   249,   849,   249,   849,   249,   849,   249,
     249,   850,   249,   315,   249,   445,   849,   249,   249,   249,
     849,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   850,   249,   249,   230,   231,   249,   231,
     249,   849,   249,   249,   249,   249,   249,   849,   850,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   850,   850,   849,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   850,   249,   249,   850,   249,   249,   249,   249,
     249,   249,   174,   249,   249,   850,   249,   249,   249,   849,
     850,   249,   249,   850,   249,   850,   249,   849,   850,   249,
     849,   849,   249,   850,   849,   849,   850,   850,   850,   249,
     458,   850,   249,   458,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   249,   458,   850,   850,
     850,   850,   850,   850,   850,   850,   849,   849,   849,   849,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     249,   511,   849,   850,   849,   125,   849,   849,   125,   849,
     849,   850,   849,   375,   849,   249,   850,   249,   249,   850,
     850,   850,   729,   729,   249,   249,   850,   249,   850,   849,
     850,   249,   850,   850,   850,   249,   249,   249,   249,   850,
     849,   249,   372,   776,   280,   282,   283,   284,   285,   286,
     287,   288,   289,   672,   671,   849,   249,   850,   249,   850,
     850,   249,   249,   849,   850,   249,   850,   249,   249,   850,
     249,   249,   850,   249,   249,   372,   833,   108,   249,   249,
     850,   249,   249,   850,   249,   249,   850,   249,   249,   249,
     249,   850,   850,   850,   249,   249,   249,   249,   849,   249,
     249,   327,   327,   849,   849,   849,   616,   850,   850,   850,
     850,   850,   249,   849,   249,   249,   849,   849,   849,   849,
     850,   249,   249,   849,   249,   850,   249,   849,   249,   850,
     249,   849,   850,   249,   849,   249,   850,   249,   249,   249,
     249,   850,   249,   850,   249,   850,   850,   850,   849,   850,
     249,   850,   850,   850,   249,   850,   134,   249,   850,   249,
     850,   850,   850,   849,   850,   850,   850,   850,   850,   850,
     249,   249,   813,   850,   850,   850,   249,   850,   249,   849,
     850,   850,   850,   249,   850,   111,   249,   850,   249,   402,
     407,   849,   249,   850,   249,   249,   372,   850,   249,   249,
     850,   249,   849,   249,   849,   249,   249,   850,   249,   849,
     249,   249,   849,   249,   450,   249,   249,   231,   249,   849,
     849,   231,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   849,   249,   849,   249,   850,   850,   249,   850,
     849,   849,   849,   849,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   849,   849,   849,   849,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   849,   850,   850,   849,   249,
     375,   849,   849,   850,   249,   125,   850,   249,   729,   249,
     850,   249,   850,   249,   249,   850,   249,   249,   249,   249,
     813,   249,   849,   850,   850,   671,   671,   323,   850,   850,
     850,   850,   850,   850,   249,   850,   850,   850,   850,   850,
     850,   249,   249,   372,   833,   850,   850,   850,   850,   850,
     850,   249,   249,   850,   249,   850,   249,   249,   850,   850,
     249,   849,   850,   850,   850,   249,   850,   850,   849,   157,
     849,   157,   849,   849,   249,   850,   249,   249,   849,   249,
     850,   249,   249,   249,   850,   850,   850,   849,   242,   249,
     249,   850,   850,   850,   249,   762,   850,   134,   249,   850,
     850,   249,   850,   134,   249,   850,   249,   850,   850,   850,
     249,   850,   850,   249,   249,   813,   850,   249,   249,   850,
     249,   849,   850,   850,   249,   249,   249,   249,   407,   850,
     850,   249,   249,   249,   249,   372,   850,   249,   849,   249,
     849,   249,   849,   249,   849,   249,   450,   249,   849,   249,
     249,   849,   249,   849,   249,   849,   249,   249,   440,   442,
     813,   849,   850,   849,   850,   850,   850,   850,   850,   249,
     458,   850,   850,   249,   458,   850,   850,   249,   850,   249,
     850,   850,   249,   458,   850,   850,   850,   850,   850,   850,
     249,   458,   850,   850,   249,   850,   850,   850,   850,   850,
     849,   849,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   849,   850,   849,
     249,   249,   849,   249,   849,   249,   249,   249,   249,   850,
     850,   850,   249,   850,   249,   249,   850,   249,   850,   249,
     850,   833,   850,   249,   249,   249,   850,   249,   850,   249,
     850,   249,   249,   849,   849,   850,   249,   850,   850,   850,
     249,   850,   249,   849,   849,   249,   249,   850,   249,   850,
     850,   850,   850,   249,   850,   850,   850,   249,   850,   762,
     134,   249,   850,   249,   762,   134,   249,   850,   850,   850,
     850,   850,   850,   249,   249,   813,   849,   849,   249,   850,
     850,   850,   249,   111,   249,   249,   249,   249,   849,   249,
     849,   249,   249,   249,   249,   249,   249,   249,   849,   850,
     249,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   850,   249,   850,   249,   850,   850,   850,
     850,   850,   850,   850,   249,   458,   850,   850,   850,   850,
     249,   850,   850,   850,   249,   850,   850,   850,   850,   850,
     850,   850,   249,   458,   850,   249,   458,   850,   249,   458,
     850,   850,   850,   850,   850,   850,   850,   849,   249,   850,
     125,   249,   850,   850,   249,   249,   249,   249,   249,   249,
     833,   833,   850,   249,   249,   249,   849,   849,   249,   850,
     850,   850,   850,   249,   249,   249,   249,   242,   249,   850,
     850,   249,   379,   850,   850,   762,   762,   850,   249,   762,
     134,   249,   850,   850,   850,   850,   850,   249,   850,   850,
     850,   850,   111,   249,   249,   249,   849,   249,   849,   249,
     249,   440,   249,   850,   850,   249,   850,   850,   249,   850,
     249,   850,   850,   249,   850,   249,   458,   249,   850,   249,
     458,   850,   249,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   249,   850,   249,   458,   850,   249,   850,   249,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   850,   249,   458,   850,   850,   850,   850,   850,
     850,   850,   249,   249,   850,   850,   850,   249,   249,   833,
     849,   850,   249,   372,   850,   850,   249,   850,   249,   850,
     850,   850,   249,   379,   850,   850,   249,   762,   249,   379,
     850,   762,   762,   850,   850,   850,   850,   249,   850,   249,
     813,   850,   850,   850,   249,   249,   849,   249,   849,   849,
     850,   850,   850,   850,   249,   850,   249,   850,   850,   850,
     850,   850,   850,   249,   249,   850,   850,   850,   850,   850,
     249,   850,   249,   850,   850,   850,   850,   850,   249,   850,
     249,   249,   850,   249,   850,   850,   850,   850,   850,   850,
     850,   249,   458,   850,   249,   458,   850,   850,   249,   850,
     249,   850,   850,   850,   850,   850,   249,   125,   249,   850,
     833,   850,   249,   849,   850,   616,   249,   850,   850,   850,
     850,   850,   249,   850,   249,   379,   850,   762,   850,   249,
     379,   850,   249,   762,   849,   850,   850,   850,   850,   249,
     850,   850,   850,   249,   849,   249,   249,   813,   249,   249,
     850,   249,   249,   850,   850,   850,   850,   249,   850,   249,
     850,   249,   850,   850,   850,   249,   850,   249,   850,   249,
     850,   850,   850,   249,   850,   249,   850,   850,   850,   850,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     249,   850,   249,   850,   850,   850,   850,   850,   249,   249,
     249,   833,   249,   249,   616,   850,   249,   249,   850,   850,
     850,   249,   850,   249,   249,   850,   249,   379,   762,   249,
     850,   850,   850,   850,   850,   850,   249,   850,   249,   249,
     850,   850,   249,   850,   249,   850,   249,   850,   850,   850,
     249,   849,   850,   850,   850,   249,   849,   249,   849,   249,
     850,   850,   249,   850,   249,   850,   850,   850,   850,   850,
     850,   850,   249,   850,   249,   850,   249,   850,   850,   850,
     850,   850,   850,   850,   249,   249,   850,   249,   249,   850,
     249,   249,   249,   850,   249,   850,   850,   249,   850,   850,
     249,   850,   850,   249,   249,   249,   850,   249,   850,   249,
     850,   249,   249,   249,   849,   249,   458,   850,   249,   458,
     850,   249,   458,   249,   249,   850,   249,   249,   850,   249,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   249,
     249,   850,   249,   850,   850,   850,   850,   850,   850,   249,
     850,   249,   849,   249,   850,   850,   850,   249,   849,   249,
     849,   249,   850,   849,   850,   850,   850,   850,   850,   249,
     458,   249,   850,   249,   849,   249,   850,   249,   850,   850,
     850,   850,   850,   249,   249,   249,   850,   249,   850,   249,
     458,   850,   249,   458,   850,   249,   458,   850,   249,   850,
     249,   849,    75,   249,   850,   850,   249,   249,   249,   249,
     849,   850,   249,   458,   850,   249,   458,   850,   850,   249,
     850,   249,   850,   850,   249,   850,   249,   850,   850,   850,
     249,   849,   249,   849,   850,   850,   850,   850,   850,   249,
     458,   249,   849,   850,   850,   249,   849,   850,   850,   850,
     850,   850,   850,   249,   850,   249,   249,   850,   850,   850,
     850,   249,   249,   850,   249,   458,   850,   249,   458,   850,
     850,   249,   249,    75,   249,   850,   850,   249,   850,   249,
     850,   249,   850,   249,   850,   249,   249,   249,   850,   249,
     850,   850,   850,   850,   850,   850,   850,   850,   850,   850,
     850,   850,   249,   249,   850,   850,   850,   249,   850,   249,
     850,   249,   850,   850,   850,   850,   249,   249,   249,   850,
     249,   249,   850,   850,   249,   249,   850,   849,   249,   850,
     249,   249,    75,   249,   850,   850,   249,   850,   850,   850,
     249,   249,   850,   850,   849,   249,   850,   249,   850,   249,
     850,   850,   249,   849,   849,   850,   850,   850,   249
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   591,   592,   593,   593,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   595,   596,   597,   597,   597,   597,   598,
     598,   599,   599,   599,   599,   599,   599,   599,   599,   599,
     599,   599,   599,   599,   599,   600,   601,   602,   602,   603,
     603,   604,   604,   605,   605,   605,   605,   605,   605,   605,
     605,   605,   605,   606,   607,   607,   607,   607,   607,   607,
     607,   607,   607,   607,   607,   607,   607,   608,   608,   609,
     609,   609,   609,   609,   610,   610,   610,   611,   611,   611,
     612,   612,   612,   612,   612,   612,   612,   612,   613,   613,
     613,   613,   613,   613,   614,   614,   615,   615,   615,   615,
     616,   616,   617,   617,   618,   618,   619,   619,   620,   620,
     621,   621,   621,   622,   623,   623,   624,   624,   624,   624,
     625,   625,   625,   625,   626,   626,   627,   627,   627,   627,
     628,   629,   630,   631,   632,   633,   633,   633,   633,   633,
     634,   634,   634,   634,   634,   634,   634,   634,   634,   634,
     634,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   636,   636,   636,   636,   636,   636,   636,   636,
     636,   636,   636,   636,   636,   636,   636,   636,   636,   636,
     636,   636,   636,   636,   636,   637,   637,   638,   638,   639,
     639,   640,   640,   641,   641,   642,   642,   643,   644,   644,
     644,   645,   645,   645,   645,   645,   645,   645,   645,   645,
     645,   645,   645,   645,   645,   645,   645,   645,   645,   645,
     646,   647,   648,   648,   648,   648,   649,   649,   649,   649,
     650,   651,   651,   651,   651,   651,   652,   652,   653,   654,
     654,   654,   654,   654,   654,   654,   654,   654,   655,   655,
     656,   657,   658,   659,   659,   660,   661,   662,   663,   663,
     663,   663,   663,   663,   663,   663,   664,   664,   665,   665,
     665,   665,   666,   666,   667,   668,   669,   670,   670,   670,
     671,   671,   672,   672,   672,   672,   672,   672,   672,   672,
     672,   672,   673,   673,   673,   674,   675,   675,   676,   676,
     676,   677,   678,   679,   679,   679,   679,   679,   680,   680,
     681,   682,   682,   683,   683,   684,   685,   685,   686,   687,
     687,   687,   687,   687,   688,   688,   688,   688,   689,   689,
     689,   689,   690,   690,   690,   691,   692,   692,   693,   693,
     693,   694,   694,   695,   696,   696,   697,   698,   698,   699,
     699,   699,   700,   700,   701,   702,   702,   703,   704,   704,
     705,   706,   706,   707,   708,   709,   709,   710,   711,   711,
     712,   712,   713,   713,   714,   714,   715,   715,   716,   717,
     717,   717,   717,   718,   718,   718,   718,   718,   719,   719,
     719,   719,   719,   720,   720,   720,   720,   720,   720,   721,
     721,   722,   722,   723,   723,   723,   723,   724,   724,   725,
     726,   726,   726,   726,   726,   726,   726,   726,   727,   727,
     727,   727,   728,   728,   728,   728,   728,   728,   729,   729,
     730,   730,   731,   731,   732,   732,   732,   733,   733,   734,
     734,   734,   735,   735,   736,   736,   736,   737,   737,   738,
     738,   738,   738,   738,   739,   739,   740,   740,   740,   741,
     741,   742,   742,   742,   743,   743,   744,   744,   744,   745,
     745,   746,   746,   746,   747,   747,   748,   748,   748,   749,
     749,   749,   750,   750,   751,   751,   751,   752,   752,   752,
     752,   752,   753,   754,   754,   755,   755,   755,   756,   756,
     756,   757,   757,   758,   758,   759,   759,   760,   760,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   761,   761,   761,   761,   761,   761,   761,   761,
     761,   761,   762,   763,   763,   764,   764,   764,   764,   764,
     765,   765,   765,   765,   765,   765,   765,   765,   765,   766,
     766,   767,   767,   767,   767,   767,   767,   767,   767,   767,
     767,   767,   768,   768,   768,   769,   769,   770,   771,   771,
     771,   772,   772,   772,   772,   772,   772,   773,   773,   773,
     773,   773,   773,   773,   773,   773,   774,   774,   774,   774,
     774,   774,   774,   774,   774,   774,   774,   774,   774,   774,
     774,   774,   774,   774,   774,   774,   774,   774,   774,   774,
     774,   774,   774,   775,   775,   775,   776,   776,   776,   776,
     776,   777,   778,   778,   779,   779,   779,   779,   780,   781,
     782,   782,   783,   783,   784,   784,   785,   785,   786,   786,
     787,   787,   787,   787,   788,   788,   789,   789,   789,   789,
     789,   789,   789,   789,   789,   789,   789,   789,   789,   789,
     789,   789,   789,   789,   790,   790,   790,   790,   790,   791,
     791,   791,   792,   793,   794,   794,   795,   795,   795,   796,
     796,   796,   797,   797,   797,   797,   797,   797,   797,   797,
     797,   797,   797,   797,   797,   798,   798,   799,   800,   800,
     800,   800,   800,   801,   801,   801,   801,   801,   801,   801,
     801,   801,   801,   801,   802,   802,   802,   802,   802,   802,
     802,   802,   802,   802,   802,   802,   803,   803,   804,   804,
     805,   805,   805,   805,   805,   805,   805,   805,   805,   805,
     805,   806,   806,   807,   807,   807,   807,   807,   807,   807,
     807,   807,   807,   807,   807,   807,   807,   807,   807,   807,
     807,   807,   807,   807,   807,   807,   807,   807,   807,   807,
     808,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   809,   809,   809,
     809,   809,   809,   809,   809,   809,   809,   810,   811,   812,
     812,   813,   813,   813,   813,   813,   813,   813,   813,   813,
     813,   813,   814,   814,   814,   814,   814,   814,   814,   814,
     814,   814,   815,   815,   816,   817,   817,   818,   818,   819,
     820,   821,   822,   823,   824,   825,   825,   826,   826,   826,
     826,   826,   826,   826,   826,   826,   826,   826,   826,   826,
     826,   826,   826,   826,   826,   826,   826,   826,   826,   826,
     826,   826,   826,   826,   827,   827,   827,   828,   829,   829,
     829,   830,   830,   830,   830,   830,   830,   830,   830,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   831,   831,   831,   831,   831,
     831,   831,   831,   831,   831,   832,   832,   832,   833,   833,
     834,   834,   835,   835,   835,   836,   836,   837,   837,   837,
     837,   837,   837,   837,   837,   837,   837,   837,   837,   837,
     837,   837,   837,   837,   837,   837,   837,   837,   837,   837,
     837,   837,   837,   837,   837,   837,   837,   837,   837,   837,
     837,   838,   838,   839,   839,   840,   840,   840,   840,   840,
     841,   841,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   843,   843,   844,   845,   845,   845,   846,   846,
     847,   847,   848,   848,   848,   849,   849,   850,   850,   850,
     850
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     5,     4,     2,     5,     6,     6,     2,
       6,     2,     3,     4,     6,     6,     7,    12,    12,     7,
       2,     2,     2,     2,     2,     4,     3,     2,     4,     3,
       4,     4,     6,     2,     4,     5,     4,     5,     4,     4,
       4,     4,     4,     3,     2,     4,     4,     3,     3,     3,
       3,     3,     3,     4,     3,     3,     3,     2,     4,     2,
       4,     4,     4,     4,     2,     3,     4,     2,     3,     4,
       2,     3,     5,     5,     7,     4,     5,     5,     2,     2,
       2,     2,     2,     2,     2,     2,     4,    10,     5,    11,
       4,     5,     3,     2,     3,     2,     3,     2,     3,     2,
      11,    13,    14,     5,     2,     2,     7,     9,    11,    12,
       2,     5,     6,     5,     2,     5,     2,     5,     4,     4,
       3,     3,     3,     3,     3,     2,     2,     6,     8,     3,
       2,     2,     3,     3,     3,     3,     4,     4,     3,     3,
       3,     3,     5,     4,     5,     6,     7,     3,     5,     4,
       5,     6,     7,     3,     2,     2,     3,     2,     2,     3,
       2,     2,     2,     3,     2,     3,     2,     2,     2,     2,
       2,     2,     1,     2,     3,     3,     3,     2,     5,     3,
       3,     3,     2,     3,     2,     3,     4,     4,     5,     3,
       3,     2,     4,     2,     3,     3,     4,     4,     3,     2,
       2,     2,     3,     2,     3,     5,     2,     2,     2,     2,
       3,     2,     2,     2,     2,     3,     3,     4,     4,     8,
       4,     4,     3,     4,     6,     7,     8,     3,     5,     4,
       4,     6,     2,     3,     3,     3,     2,     4,     1,     3,
       1,     2,     2,     2,     3,     4,     5,     6,     6,     3,
       4,     6,     7,     5,     3,     4,     4,     3,     1,     2,
       6,     2,     2,     2,     3,     2,     2,     4,     5,     4,
       2,     6,     8,     7,     5,     9,     3,     2,     3,     2,
       4,     3,     2,     2,     2,     2,     5,     5,     6,     7,
       3,     1,     1,     2,     1,     1,     1,     2,     1,     1,
       2,     1,     4,     5,     3,     3,     6,     5,     6,     3,
       2,     5,     6,     2,     2,     7,     9,     3,     2,     6,
       3,     1,     2,     3,     2,     3,     1,     2,     4,     2,
       4,     6,     8,     5,     2,     3,     4,     5,     2,     3,
       6,     7,     2,     3,     6,     3,     1,     2,     4,     5,
       6,     3,     2,     4,     1,     2,     3,     1,     2,     4,
       5,     6,     3,     2,     4,     3,     2,     4,     3,     2,
       4,     3,     2,     4,     3,     1,     2,     3,     3,     4,
       3,     2,     3,     2,     5,     2,     3,     4,     4,     5,
       5,     6,     6,     3,     4,     3,     2,     6,     2,     3,
       3,     4,     5,     3,     2,     9,     6,     2,     2,     3,
       9,     3,     9,     2,     3,     4,     5,     3,     4,     3,
       2,     3,     2,     7,     9,     8,    10,     3,     5,     6,
       6,     7,     1,     2,     6,     7,     8,     9,     1,     2,
       1,     2,     2,     3,     6,     4,     7,     2,     3,     6,
       4,     7,     2,     3,     6,     4,     7,     2,     3,     8,
      10,     4,     9,    11,     2,     3,     6,     4,     7,     2,
       3,     6,     4,     7,     2,     3,     6,     4,     7,     2,
       3,     6,     4,     7,     2,     3,     9,     7,    10,     2,
       3,     5,     5,     3,     2,     2,     3,     2,     3,     5,
       4,     6,     4,     2,     3,     2,     2,     3,     5,     3,
       2,     5,     4,     3,     4,     1,     2,     3,     2,    16,
      19,    20,    23,     9,    14,    16,    30,    12,    13,    14,
       5,     6,    16,    12,     3,     4,    13,     5,     6,     6,
       7,    11,     5,     6,     8,     9,    10,     9,    10,    11,
      10,    11,    12,    11,    12,    13,     5,     7,     6,     8,
       6,     9,     7,    10,     7,    11,     8,    12,     4,     6,
      12,     4,     5,     3,     2,     3,     4,     5,     6,     5,
       4,     5,     5,     6,     7,     7,     8,     7,     8,     4,
       3,     2,     5,     6,     7,     8,    10,     6,     7,     8,
       9,    11,     2,     5,     7,     3,     2,     4,     2,     5,
       7,     2,     4,     3,     4,     6,     5,     1,     3,     4,
       5,     6,     8,     9,    11,    12,     2,     4,     3,     3,
       3,     3,     4,     3,     4,     4,     3,     3,     3,     3,
       3,     3,     3,     3,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     3,     6,     2,     5,     4,     3,     7,
       6,     4,     1,     2,     4,     3,     5,     4,     3,     3,
       5,     4,     2,     2,     2,     2,    11,     4,     2,     2,
      11,    14,    12,    15,     2,     7,     2,     3,     5,     6,
       4,     6,     7,     7,     6,     8,     9,     7,     9,    10,
       5,     5,     7,     8,     2,     3,     4,     3,     3,     2,
       2,     2,     5,     3,     2,     3,     2,     4,     3,     2,
       4,     5,     2,     3,     4,     5,     5,     7,     5,     6,
       6,     6,     7,     7,     8,     2,     3,     3,     2,     4,
       6,     6,     8,     2,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     3,     4,     4,     5,     5,     6,
       6,     7,     7,     8,     8,     9,     1,     2,     1,     2,
       2,     2,     4,     3,     5,     4,     4,     5,     6,     4,
       4,     1,     2,     2,     3,     2,     3,     3,     3,     2,
       3,     4,     5,     6,     7,     2,     3,     4,     5,     4,
       3,     3,     3,     2,     4,     5,     6,     7,     2,     3,
       4,     1,     3,     4,     4,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     5,     5,     4,
       4,     4,     4,     4,     4,     4,     4,     3,     4,     4,
       4,     4,     4,     4,     5,     4,     5,     5,     5,     5,
       3,     4,     4,     4,     4,     4,     4,     4,     3,     3,
       4,     3,     3,     4,     5,     4,     4,     5,     5,     6,
       6,     7,     4,     4,     5,     4,     4,     4,     4,     4,
       3,     4,     3,     3,     3,     3,     4,     3,     4,     4,
       3,     4,     4,     4,     4,     4,     4,     3,     3,     4,
       5,     4,     4,     4,     4,     6,     6,     5,     5,     7,
       7,     4,     4,     5,     4,     4,     4,     3,     2,     3,
       4,     1,     2,     3,     4,     5,     1,     2,     3,     2,
       3,     4,     2,     3,     4,     5,     4,     4,     4,     2,
       2,     2,     1,     2,     4,     6,     5,     1,     2,     4,
       3,     2,     3,     2,     2,     1,     1,     2,     3,     3,
       3,     3,     4,     4,     4,     5,     7,     5,     7,     4,
       5,     3,     4,     4,     5,     6,     7,     8,     3,     4,
       6,     8,     4,     2,     3,     4,     5,     2,     2,    10,
      12,     2,     4,     7,     9,     9,     8,    11,    12,     2,
       9,    10,    12,    13,    14,    15,     9,    10,    12,    13,
      14,    15,    12,    13,    14,    15,     9,    10,    12,    13,
      14,    15,     9,     7,     5,    13,    11,     9,     9,     7,
       5,    13,    11,     9,     9,     7,     5,    13,    11,     9,
       8,    10,     8,    10,     9,    11,     9,    11,    11,    13,
      11,    13,     8,    10,    12,    14,     8,    10,    12,    14,
       8,    12,     9,    13,    10,    11,    13,    14,    15,    16,
      10,    11,    13,    14,    15,    16,     7,    13,    15,    17,
      19,    14,    16,    14,    16,    15,    17,    15,    17,    17,
      19,    17,    19,    14,    16,    18,    20,    13,    15,    17,
      19,    14,    16,    18,    20,     7,    11,    13,    17,    14,
      18,     8,    12,    14,    18,    15,    19,     8,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,     7,     9,    10,    10,    11,    12,
      13,    10,    11,    12,    13,     7,     8,     9,    10,    11,
      12,    13,    22,     5,     5,     2,     4,     5,     0,     2,
       0,     2,     4,     6,     8,     2,     3,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     3,     2,
       3,     2,     2,     2,     2,     5,     6,     3,     6,     3,
       2,     1,     2,     3,     2,     3,     4,     2,     2,     3,
       1,     2,     3,     2,     3,     2,     3,     2,     2,     2,
       2,     3,     2,     3,     3,     2,     3,     4,     2,     2,
       2,     2,     2,     3,     4,     2,     2,     2,     2,     2,
       2,     2,     3,     4,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     3,     4,     2,     3,     4,     2,
       2,     2,     2,     2,     2,     2,     2,     3,     3,     2,
       6,     3,     2,     3,     5,     2,     5,     3,     2,     3,
       2,     3,     2,     3,     2,     1,     1,     1,     1,     1,
       1
};


#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)
#define YYEMPTY         (-2)
#define YYEOF           0

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
do                                                              \
  if (yychar == YYEMPTY)                                        \
    {                                                           \
      yychar = (Token);                                         \
      yylval = (Value);                                         \
      YYPOPSTACK (yylen);                                       \
      yystate = *yyssp;                                         \
      goto yybackup;                                            \
    }                                                           \
  else                                                          \
    {                                                           \
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;                                                  \
    }                                                           \
while (0)

/* Error token number */
#define YYTERROR        1
#define YYERRCODE       256



/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


# define YY_SYMBOL_PRINT(Title, Type, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Type, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*----------------------------------------.
| Print this symbol's value on YYOUTPUT.  |
`----------------------------------------*/

static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  YYUSE (yytype);
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyoutput, "%s %s (",
             yytype < YYNTOKENS ? "token" : "nterm", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yytype_int16 *yyssp, YYSTYPE *yyvsp, int yyrule)
{
  unsigned long int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       yystos[yyssp[yyi + 1 - yynrhs]],
                       &(yyvsp[(yyi + 1) - (yynrhs)])
                                              );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
yystrlen (const char *yystr)
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
yystpcpy (char *yydest, const char *yysrc)
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
        switch (*++yyp)
          {
          case '\'':
          case ',':
            goto do_not_strip_quotes;

          case '\\':
            if (*++yyp != '\\')
              goto do_not_strip_quotes;
            /* Fall through.  */
          default:
            if (yyres)
              yyres[yyn] = *yyp;
            yyn++;
            break;

          case '"':
            if (yyres)
              yyres[yyn] = '\0';
            return yyn;
          }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULLPTR, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULLPTR;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULLPTR, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yytype);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}




/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       'yyss': related to states.
       'yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        yytype_int16 *yyss1 = yyss;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),
                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yytype_int16 *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 197 "p.y" /* yacc.c:1646  */
    { 
          if(domain->solInfo().piecewise || domain->solInfo().freeplay) domain->solInfo().activatePiecewise();
          return 0;
         }
#line 5490 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 209 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; delete (yyvsp[0].bclist); }
#line 5496 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 211 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5502 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 16:
#line 221 "p.y" /* yacc.c:1646  */
    { int j = geoSource->getLocalIndex();
          if(geoSource->elementLumpingWeightLocalSize(j)>0) geoSource->setLocalIndex(j+1); }
#line 5509 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 233 "p.y" /* yacc.c:1646  */
    {}
#line 5515 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 236 "p.y" /* yacc.c:1646  */
    {}
#line 5521 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 34:
#line 243 "p.y" /* yacc.c:1646  */
    {}
#line 5527 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 250 "p.y" /* yacc.c:1646  */
    { if(geoSource->setUsddLocation((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)    return -1; delete (yyvsp[0].bclist); }
#line 5534 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 267 "p.y" /* yacc.c:1646  */
    { domain->setMFTT((yyvsp[0].mftval).first, (yyvsp[0].mftval).second); }
#line 5540 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 269 "p.y" /* yacc.c:1646  */
    { domain->setHFTT((yyvsp[0].hftval).first, (yyvsp[0].hftval).second); }
#line 5546 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 89:
#line 303 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5552 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 90:
#line 305 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5558 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 91:
#line 307 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5564 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 92:
#line 309 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichletFluid((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5570 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 93:
#line 311 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichletFluid((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5576 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 94:
#line 313 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5582 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 111:
#line 331 "p.y" /* yacc.c:1646  */
    { if(domain->setComplexNeuman((yyvsp[0].cxbclist)->n,(yyvsp[0].cxbclist)->d) < 0) return -1; }
#line 5588 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 113:
#line 334 "p.y" /* yacc.c:1646  */
    { if(domain->setComplexDirichlet((yyvsp[0].cxbclist)->n,(yyvsp[0].cxbclist)->d) < 0) return -1; }
#line 5594 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 117:
#line 339 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5600 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 118:
#line 341 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5606 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 125:
#line 349 "p.y" /* yacc.c:1646  */
    {}
#line 5612 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 134:
#line 360 "p.y" /* yacc.c:1646  */
    {}
#line 5618 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 135:
#line 362 "p.y" /* yacc.c:1646  */
    {}
#line 5624 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 136:
#line 364 "p.y" /* yacc.c:1646  */
    {}
#line 5630 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 137:
#line 366 "p.y" /* yacc.c:1646  */
    {}
#line 5636 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 138:
#line 368 "p.y" /* yacc.c:1646  */
    {}
#line 5642 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 139:
#line 370 "p.y" /* yacc.c:1646  */
    {}
#line 5648 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 140:
#line 372 "p.y" /* yacc.c:1646  */
    {}
#line 5654 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 141:
#line 374 "p.y" /* yacc.c:1646  */
    {}
#line 5660 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 153:
#line 389 "p.y" /* yacc.c:1646  */
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[-2].ival)); 
            domain->solInfo().nsample = (yyvsp[-1].ival);
          }
#line 5669 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 154:
#line 396 "p.y" /* yacc.c:1646  */
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[-1].ival));
          }
#line 5677 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 156:
#line 403 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-3].ival) == OutputInfo::Attribute)  geoSource->setAttributeGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1);
          else if ((yyvsp[-3].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[-3].ival));  exit(-1); }
        }
#line 5686 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 157:
#line 408 "p.y" /* yacc.c:1646  */
    { int i;
          if ((yyvsp[-4].ival) == OutputInfo::Attribute)  {
            for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
              geoSource->setAttributeGroup(i-1,(yyvsp[-1].ival)-1);
          }
          else if ((yyvsp[-4].ival) == OutputInfo::Nodal)  {
            for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
              geoSource->setNodeGroup(i-1, (yyvsp[-1].ival));
          }
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[-4].ival));  exit(-1); }
        }
#line 5702 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 158:
#line 420 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-4].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[-4].ival));  exit(-1); }
        }
#line 5710 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 160:
#line 427 "p.y" /* yacc.c:1646  */
    { geoSource->setGroupRandomProperty((yyvsp[-4].ival)-1,(yyvsp[-3].rprop),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 5716 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 161:
#line 431 "p.y" /* yacc.c:1646  */
    { domain->solInfo().curSweepParam = 0; }
#line 5722 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 162:
#line 433 "p.y" /* yacc.c:1646  */
    { domain->solInfo().curSweepParam = (yyvsp[-1].ival); }
#line 5728 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 163:
#line 435 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-1].fval)); }
#line 5734 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 164:
#line 437 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-3].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies1(2.0*PI*(yyvsp[-3].fval), 2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5742 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 165:
#line 441 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-3].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies2(2.0*PI*(yyvsp[-3].fval), 2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5750 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 166:
#line 445 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-4].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies(2.0*PI*(yyvsp[-4].fval), 2.0*PI*(yyvsp[-3].fval), (yyvsp[-2].ival), (yyvsp[-1].ival)); }
#line 5758 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 167:
#line 449 "p.y" /* yacc.c:1646  */
    {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-9].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*(yyvsp[-9].fval), 2.0*PI*(yyvsp[-8].fval), 2, (yyvsp[-7].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = (yyvsp[-4].ival);
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[-7].ival);
          if ((yyvsp[-6].ival) == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
          else if ((yyvsp[-6].ival) == SweepParams::WCAWEGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
          else
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[-9].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[-8].fval);
          domain->solInfo().getSweepParams()->adaptSweep.atol = (yyvsp[-5].fval);
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = (yyvsp[-2].ival);
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = (yyvsp[-1].ival);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-2].ival);
        }
#line 5784 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 168:
#line 471 "p.y" /* yacc.c:1646  */
    {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-9].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*(yyvsp[-9].fval), 2.0*PI*(yyvsp[-8].fval), 2, (yyvsp[-7].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = (yyvsp[-4].ival);
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[-7].ival);
          if ((yyvsp[-6].ival) == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
          else if ((yyvsp[-6].ival) == SweepParams::WCAWEGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
          else
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[-9].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[-8].fval);
          domain->solInfo().getSweepParams()->adaptSweep.atol = (yyvsp[-5].fval);
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = -1;
          domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-3].ival);
          domain->solInfo().getSweepParams()->adaptSweep.ctolf = (yyvsp[-2].fval);
          domain->solInfo().getSweepParams()->adaptSweep.tol1f = (yyvsp[-1].fval);
        }
#line 5812 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 169:
#line 495 "p.y" /* yacc.c:1646  */
    {
          if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-4].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addFrequencies(2.0*PI*(yyvsp[-4].fval), 2.0*PI*(yyvsp[-3].fval), 2, (yyvsp[-2].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = 6;
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[-2].ival);
          if ((yyvsp[-1].ival) == SweepParams::KrylovGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 0; 
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else if ((yyvsp[-1].ival) == SweepParams::WCAWEGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 2;
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = 1;
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 16;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[-4].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[-3].fval);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
        }
#line 5849 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 175:
#line 535 "p.y" /* yacc.c:1646  */
    {
          if((yyvsp[-3].ival) == 1) {
            domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
            domain->solInfo().getSweepParams()->alphaD = (yyvsp[-1].fval);
            domain->solInfo().getSweepParams()->betaD = (yyvsp[-2].fval);
            domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
          }
          else return -1; // only RAYDAMP is allowed here
        }
#line 5863 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 176:
#line 547 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_pivot = true; domain->solInfo().getSweepParams()->pade_tol = (yyvsp[-1].fval); }
#line 5869 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 177:
#line 551 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_poles = true; }
#line 5875 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 178:
#line 553 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_poles = true; 
          domain->solInfo().getSweepParams()->pade_poles_sigmaL = (yyvsp[-2].fval); domain->solInfo().getSweepParams()->pade_poles_sigmaU = (yyvsp[-1].fval); }
#line 5882 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 179:
#line 558 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-1].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addCoarseFrequency(2.0*PI*(yyvsp[-1].fval)); }
#line 5890 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 180:
#line 563 "p.y" /* yacc.c:1646  */
    { domain->addFrequencies(2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5896 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 181:
#line 567 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->freqSweepMethod = (yyvsp[-2].ival); 
          int &l = domain->solInfo().getSweepParams()->padeL,
              &m = domain->solInfo().getSweepParams()->padeM,
              &n = domain->solInfo().getSweepParams()->padeN;
          switch((yyvsp[-2].ival)) {
            case SweepParams::Taylor:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-1].ival)+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = (yyvsp[-1].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = (yyvsp[-1].ival);
              if(m%n != 0) m = (m/n+1)*n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = (yyvsp[-1].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = (yyvsp[-1].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::WCAWEGalProjection:
              n = (yyvsp[-1].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        }
#line 5942 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 182:
#line 609 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->freqSweepMethod = (yyvsp[-4].ival);
          int &l = domain->solInfo().getSweepParams()->padeL,
              &m = domain->solInfo().getSweepParams()->padeM,
              &n = domain->solInfo().getSweepParams()->padeN;
          switch((yyvsp[-4].ival)) {
            case SweepParams::Taylor:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[-3].ival)+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              l = (yyvsp[-2].ival);
              m = (yyvsp[-1].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival); 
              m = (yyvsp[-1].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = (yyvsp[-3].ival);
              m = (yyvsp[-1].ival);
              if(m%n != 0) m = (m/n+1)*n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::WCAWEGalProjection:
              n = (yyvsp[-3].ival);
              l = (yyvsp[-2].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        }
#line 5996 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 184:
#line 662 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInput = bool((yyvsp[-1].ival)); }
#line 6002 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 185:
#line 664 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInput = bool((yyvsp[-2].ival));
            std::string prefix = (yyvsp[-1].strval);
            clusterData_ = prefix + ".msh";
            decomposition_ = prefix + ".dec";
            connectivity_ = prefix + ".con";
            subdomains_ = prefix + ".sub";
          }
#line 6014 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 186:
#line 672 "p.y" /* yacc.c:1646  */
    { geoSource->binaryOutput = bool((yyvsp[-1].ival)); }
#line 6020 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 187:
#line 674 "p.y" /* yacc.c:1646  */
    { geoSource->binaryOutput = bool((yyvsp[-2].ival));
            int len = strlen((yyvsp[-1].strval));
            char *file = new char[len+5];
            strcpy(file, (yyvsp[-1].strval));
            strcat(file,".con");
            geoSource->setGlob(file);
          }
#line 6032 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 188:
#line 682 "p.y" /* yacc.c:1646  */
    { geoSource->setGeo((yyvsp[-1].strval)); }
#line 6038 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 189:
#line 684 "p.y" /* yacc.c:1646  */
    { geoSource->setDecomp((yyvsp[-1].strval)); }
#line 6044 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 190:
#line 686 "p.y" /* yacc.c:1646  */
    { geoSource->setGlob((yyvsp[-1].strval)); }
#line 6050 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 191:
#line 688 "p.y" /* yacc.c:1646  */
    { geoSource->setMatch((yyvsp[-1].strval)); }
#line 6056 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 192:
#line 690 "p.y" /* yacc.c:1646  */
    { geoSource->setCpuMap((yyvsp[-1].strval)); }
#line 6062 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 193:
#line 694 "p.y" /* yacc.c:1646  */
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[-1].ival)); 
#endif
	}
#line 6072 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 194:
#line 702 "p.y" /* yacc.c:1646  */
    {if(decInit==0) decInit = new DecInit(); }
#line 6078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 195:
#line 704 "p.y" /* yacc.c:1646  */
    {decInit->file = strdup((yyvsp[-1].strval));}
#line 6084 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 196:
#line 706 "p.y" /* yacc.c:1646  */
    {decInit->nsubs = (yyvsp[-1].ival); }
#line 6090 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 197:
#line 708 "p.y" /* yacc.c:1646  */
    {decInit->weight = true; }
#line 6096 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 198:
#line 710 "p.y" /* yacc.c:1646  */
    {decInit->memory = true; }
#line 6102 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 199:
#line 712 "p.y" /* yacc.c:1646  */
    {decInit->exitAfterDec = true;}
#line 6108 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 200:
#line 714 "p.y" /* yacc.c:1646  */
    {decInit->skip = true;}
#line 6114 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 201:
#line 716 "p.y" /* yacc.c:1646  */
    {decInit->nosa = true; }
#line 6120 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 202:
#line 718 "p.y" /* yacc.c:1646  */
    {decInit->trivial = true; }
#line 6126 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 203:
#line 720 "p.y" /* yacc.c:1646  */
    {decInit->trivial = true; randomShuffle = bool((yyvsp[-1].ival)); }
#line 6132 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 204:
#line 722 "p.y" /* yacc.c:1646  */
    {decInit->fsgl = true; }
#line 6138 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 205:
#line 724 "p.y" /* yacc.c:1646  */
    {allowMechanisms = true; }
#line 6144 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 206:
#line 726 "p.y" /* yacc.c:1646  */
    {useScotch = true; }
#line 6150 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 207:
#line 730 "p.y" /* yacc.c:1646  */
    {}
#line 6156 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 208:
#line 732 "p.y" /* yacc.c:1646  */
    { weightList[(yyvsp[-2].ival)] = (yyvsp[-1].fval); }
#line 6162 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 209:
#line 736 "p.y" /* yacc.c:1646  */
    {}
#line 6168 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 210:
#line 738 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Acoustic] = (yyvsp[-1].ival); }
#line 6174 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 211:
#line 740 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Structural] = (yyvsp[-1].ival); }
#line 6180 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 212:
#line 742 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Thermal] = (yyvsp[-1].ival); }
#line 6186 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 213:
#line 744 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Fluid] = (yyvsp[-1].ival); }
#line 6192 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 214:
#line 747 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = 0; }
#line 6198 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 215:
#line 749 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = (yyvsp[-1].ival); }
#line 6204 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 216:
#line 751 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first->add((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6210 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 217:
#line 755 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = 0; }
#line 6216 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 218:
#line 757 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = (yyvsp[-1].ival); }
#line 6222 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 219:
#line 759 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first->add((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6228 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 220:
#line 763 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 6234 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 221:
#line 765 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 6240 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 222:
#line 767 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival); domain->setLoadFactorGrav((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6246 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 223:
#line 769 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival); domain->setLoadFactorTemp((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6252 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 224:
#line 771 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-5].ival); domain->setLoadFactorGrav((yyvsp[-5].ival), (yyvsp[-3].ival)); domain->setLoadFactorTemp((yyvsp[-5].ival), (yyvsp[-1].ival)); }
#line 6258 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 225:
#line 773 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactor((yyval.ival), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 6264 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 226:
#line 775 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactorMFTT((yyval.ival), (yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6270 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 227:
#line 777 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactorHFTT((yyval.ival), (yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6276 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 235:
#line 790 "p.y" /* yacc.c:1646  */
    { geoSource->addCFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d); }
#line 6282 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 236:
#line 794 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).coefFlag = false; geoSource->addCoefInfo((yyvsp[-2].ival)-1,(yyvsp[0].coefdata)); }
#line 6288 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 237:
#line 796 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).c[6][0] = (yyvsp[-7].fval);
          (yyvsp[0].coefdata).c[6][1] = (yyvsp[-6].fval);
          (yyvsp[0].coefdata).c[6][2] = (yyvsp[-5].fval);
          (yyvsp[0].coefdata).c[6][3] = (yyvsp[-4].fval);
          (yyvsp[0].coefdata).c[6][4] = (yyvsp[-3].fval);
          (yyvsp[0].coefdata).c[6][5] = (yyvsp[-2].fval);
          (yyvsp[0].coefdata).coefFlag = false;
          geoSource->addCoefInfo((yyvsp[-8].ival)-1,(yyvsp[0].coefdata)); }
#line 6301 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 238:
#line 805 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).coefFlag = (yyvsp[-2].ival); geoSource->addCoefInfo((yyvsp[-3].ival)-1,(yyvsp[0].coefdata)); }
#line 6307 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 239:
#line 807 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).c[6][0] = (yyvsp[-8].fval);
          (yyvsp[0].coefdata).c[6][1] = (yyvsp[-7].fval);
          (yyvsp[0].coefdata).c[6][2] = (yyvsp[-6].fval);
          (yyvsp[0].coefdata).c[6][3] = (yyvsp[-5].fval);
          (yyvsp[0].coefdata).c[6][4] = (yyvsp[-4].fval);
          (yyvsp[0].coefdata).c[6][5] = (yyvsp[-3].fval);
          (yyvsp[0].coefdata).coefFlag = (yyvsp[-2].ival);
          geoSource->addCoefInfo((yyvsp[-9].ival)-1,(yyvsp[0].coefdata)); }
#line 6320 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 240:
#line 818 "p.y" /* yacc.c:1646  */
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6326 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 241:
#line 820 "p.y" /* yacc.c:1646  */
    { (yyval.coefdata).setCoef((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6332 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 242:
#line 824 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6338 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 243:
#line 826 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6344 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 244:
#line 830 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6350 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 245:
#line 832 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6356 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 246:
#line 836 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6362 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 247:
#line 838 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6368 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 248:
#line 842 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6374 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 249:
#line 844 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6380 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 250:
#line 848 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-10].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-9].fval); (yyval.ldata).d[1] = (yyvsp[-8].fval); (yyval.ldata).d[2] = (yyvsp[-7].fval);
	  (yyval.ldata).d[3] = (yyvsp[-6].fval); (yyval.ldata).d[4] = (yyvsp[-5].fval); (yyval.ldata).d[5] = (yyvsp[-4].fval);
	  (yyval.ldata).d[6] = (yyvsp[-3].fval); (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval);
          (yyval.ldata).d[9] = 0;  (yyval.ldata).d[10] = 0; (yyval.ldata).d[11] = 0; }
#line 6391 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 251:
#line 855 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-12].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-11].fval); (yyval.ldata).d[1] = (yyvsp[-10].fval); (yyval.ldata).d[2] = (yyvsp[-9].fval);
          (yyval.ldata).d[3] = (yyvsp[-8].fval); (yyval.ldata).d[4] = (yyvsp[-7].fval); (yyval.ldata).d[5] = (yyvsp[-6].fval);
          (yyval.ldata).d[6] = (yyvsp[-5].fval); (yyval.ldata).d[7] = (yyvsp[-4].fval); (yyval.ldata).d[8] = (yyvsp[-3].fval);
          (yyval.ldata).d[9] = (yyvsp[-2].fval);(yyval.ldata).d[10]= (yyvsp[-1].fval);(yyval.ldata).d[11] = 0; }
#line 6402 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 252:
#line 862 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-13].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-12].fval); (yyval.ldata).d[1] = (yyvsp[-11].fval); (yyval.ldata).d[2] = (yyvsp[-10].fval);
          (yyval.ldata).d[3] = (yyvsp[-9].fval); (yyval.ldata).d[4] = (yyvsp[-8].fval); (yyval.ldata).d[5] = (yyvsp[-7].fval);
          (yyval.ldata).d[6] = (yyvsp[-6].fval); (yyval.ldata).d[7] = (yyvsp[-5].fval); (yyval.ldata).d[8] = (yyvsp[-4].fval);
          (yyval.ldata).d[9] = (yyvsp[-3].fval);(yyval.ldata).d[10]= (yyvsp[-2].fval); (yyval.ldata).d[11] = (yyvsp[-1].fval); }
#line 6413 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 253:
#line 871 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-4].ival)-1;  (yyval.ldata).matid = (yyvsp[-3].ival)-1; (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval); }
#line 6419 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 255:
#line 876 "p.y" /* yacc.c:1646  */
    { geoSource->addLayMat((yyvsp[0].ldata).matid, (yyvsp[0].ldata).d); }
#line 6425 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 256:
#line 882 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-6].ival)-1; (yyval.ldata).d[0] = (yyvsp[-5].fval); (yyval.ldata).d[1] = (yyvsp[-4].fval); (yyval.ldata).d[2] = (yyvsp[-3].fval);
          (yyval.ldata).d[3] = (yyvsp[-2].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[-1].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
#line 6433 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 257:
#line 887 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-8].ival)-1; (yyval.ldata).d[0] = (yyvsp[-7].fval); (yyval.ldata).d[1] = (yyvsp[-6].fval); (yyval.ldata).d[2] = (yyvsp[-5].fval);
          (yyval.ldata).d[3] = (yyvsp[-4].fval); (yyval.ldata).d[4] = (yyvsp[-3].fval); (yyval.ldata).d[5] = (yyvsp[-2].fval); (yyval.ldata).d[6] = (yyvsp[-1].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
#line 6441 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 258:
#line 892 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-10].ival)-1; (yyval.ldata).d[0] = (yyvsp[-9].fval); (yyval.ldata).d[1] = (yyvsp[-8].fval); (yyval.ldata).d[2] = (yyvsp[-7].fval);
          (yyval.ldata).d[3] = (yyvsp[-6].fval); (yyval.ldata).d[4] = (yyvsp[-5].fval); (yyval.ldata).d[5] = (yyvsp[-4].fval); (yyval.ldata).d[6] = (yyvsp[-3].fval);
          (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval); (yyval.ldata).d[9] = 0; }
#line 6449 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 259:
#line 896 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-11].ival)-1; (yyval.ldata).d[0] = (yyvsp[-10].fval); (yyval.ldata).d[1] = (yyvsp[-9].fval); (yyval.ldata).d[2] = (yyvsp[-8].fval);
          (yyval.ldata).d[3] = (yyvsp[-7].fval); (yyval.ldata).d[4] = (yyvsp[-6].fval); (yyval.ldata).d[5] = (yyvsp[-5].fval); (yyval.ldata).d[6] = (yyvsp[-4].fval); 
          (yyval.ldata).d[7] = (yyvsp[-3].fval); (yyval.ldata).d[8] = (yyvsp[-2].fval); (yyval.ldata).d[9] = (yyvsp[-1].fval); }
#line 6457 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 261:
#line 903 "p.y" /* yacc.c:1646  */
    { domain->addDMass((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6463 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 262:
#line 905 "p.y" /* yacc.c:1646  */
    { domain->addDMass((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-1].fval),(yyvsp[-2].ival)-1); }
#line 6469 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 263:
#line 907 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modalDIMASS = true;
          domain->solInfo().reducedMassFile = (yyvsp[-1].strval); }
#line 6476 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 265:
#line 913 "p.y" /* yacc.c:1646  */
    { domain->setGravity((yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6482 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 267:
#line 918 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[-3].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[-2].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[-1].strval); }
#line 6490 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 268:
#line 922 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[-2].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval);}
#line 6497 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 269:
#line 925 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[-2].strval);
          domain->solInfo().nRestart = (yyvsp[-1].ival); }
#line 6504 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 270:
#line 930 "p.y" /* yacc.c:1646  */
    { geoSource->setControlFile((yyvsp[-1].strval));
         geoSource->setControlRoutine((char *) "controlObj");}
#line 6511 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 271:
#line 935 "p.y" /* yacc.c:1646  */
    { geoSource->setControlRoutine((yyvsp[-1].strval)); }
#line 6517 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 272:
#line 939 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 6524 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 273:
#line 944 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Actuators; }
          if(geoSource->setActuatorLocations((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)            return -1; }
#line 6532 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 274:
#line 950 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Usdf; }
          if(geoSource->setUsdfLocation((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)       return -1; }
#line 6541 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 275:
#line 957 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInputControlLeft = true;
          (yyval.bclist) = new BCList; }
#line 6548 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 276:
#line 960 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Usdd; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 6554 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 277:
#line 962 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-4].ival); i<=(yyvsp[-2].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-1].ival)-1, 0., BCond::Usdd); (yyval.bclist)->add(bc); } }
#line 6560 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 278:
#line 964 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); i+=(yyvsp[-2].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-1].ival)-1, 0., BCond::Usdd); (yyval.bclist)->add(bc); } }
#line 6566 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 279:
#line 966 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Usdd;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
#line 6576 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 280:
#line 974 "p.y" /* yacc.c:1646  */
    { numColumns = 3; }
#line 6582 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 281:
#line 976 "p.y" /* yacc.c:1646  */
    { numColumns = 6; }
#line 6588 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 282:
#line 978 "p.y" /* yacc.c:1646  */
    { numColumns = 3; geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6594 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 283:
#line 980 "p.y" /* yacc.c:1646  */
    { numColumns = 6; geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6600 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 284:
#line 982 "p.y" /* yacc.c:1646  */
    { numColumns = 3; domain->outFlag = (yyvsp[-1].ival); }
#line 6606 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 285:
#line 984 "p.y" /* yacc.c:1646  */
    { numColumns = 6; domain->outFlag = (yyvsp[-1].ival); }
#line 6612 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 286:
#line 986 "p.y" /* yacc.c:1646  */
    { numColumns = 3; domain->outFlag = (yyvsp[-2].ival); geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6618 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 287:
#line 988 "p.y" /* yacc.c:1646  */
    { numColumns = 6; domain->outFlag = (yyvsp[-2].ival); geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6624 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 288:
#line 990 "p.y" /* yacc.c:1646  */
    { numColumns = 3; geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval); }
#line 6630 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 289:
#line 992 "p.y" /* yacc.c:1646  */
    { numColumns = 6; geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval); }
#line 6636 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 290:
#line 994 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[-1].oinfo)); }
#line 6642 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 291:
#line 998 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6648 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 292:
#line 1000 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-4].ival); (yyval.oinfo).width = (yyvsp[-3].ival); (yyval.oinfo).precision = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6654 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 293:
#line 1002 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6660 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 294:
#line 1004 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); 
          if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1;}
#line 6667 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 295:
#line 1007 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-5].ival); (yyval.oinfo).width = (yyvsp[-4].ival); (yyval.oinfo).precision = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6673 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 296:
#line 1009 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-6].ival); (yyval.oinfo).width = (yyvsp[-5].ival); (yyval.oinfo).precision = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6679 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 297:
#line 1011 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6685 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 298:
#line 1013 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-4].ival); (yyval.oinfo).width = (yyvsp[-3].ival); (yyval.oinfo).precision = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6691 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 299:
#line 1015 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6697 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 300:
#line 1017 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6703 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 301:
#line 1019 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-5].ival); (yyval.oinfo).width = (yyvsp[-4].ival); (yyval.oinfo).precision = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6709 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 302:
#line 1021 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-6].ival); (yyval.oinfo).width = (yyvsp[-5].ival); (yyval.oinfo).precision = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6715 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 303:
#line 1024 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6721 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 304:
#line 1026 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).surface = (yyvsp[0].ival); }
#line 6727 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 305:
#line 1028 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).str_therm_option = (yyvsp[0].ival); }
#line 6733 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 306:
#line 1030 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ylayer = (yyvsp[-1].fval); (yyval.oinfo).zlayer = (yyvsp[0].fval); }
#line 6739 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 307:
#line 1032 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).averageFlg = (yyvsp[0].ival); }
#line 6745 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 308:
#line 1034 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).complexouttype = (yyvsp[0].ival); }
#line 6751 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 309:
#line 1036 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).complexouttype = (yyvsp[-1].ival); (yyval.oinfo).ncomplexout = (yyvsp[0].ival); }
#line 6757 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 310:
#line 1038 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).angularouttype = (yyvsp[0].ival); }
#line 6763 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 311:
#line 1040 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rotvecouttype = (yyvsp[0].ival); }
#line 6769 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 312:
#line 1042 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rotvecouttype = OutputInfo::Linear; }
#line 6775 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 313:
#line 1044 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rescaling = bool((yyvsp[0].ival)); }
#line 6781 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 314:
#line 1046 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ndtype = (yyvsp[0].ival); }
#line 6787 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 315:
#line 1048 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ndtype = (yyvsp[-1].ival); sfem->setnsamp_out((yyvsp[0].ival)); }
#line 6793 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 316:
#line 1050 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).oframe = (OutputInfo::FrameType) (yyvsp[0].ival); }
#line 6799 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 317:
#line 1052 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).matlab = true; }
#line 6805 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 318:
#line 1054 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xmatrixname = (yyvsp[0].strval); }
#line 6811 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 319:
#line 1056 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qmatrixname = (yyvsp[0].strval); }
#line 6817 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 320:
#line 1058 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rmatrixname = (yyvsp[0].strval); }
#line 6823 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 321:
#line 1060 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenvaluename = (yyvsp[0].strval); }
#line 6829 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 323:
#line 1065 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; }
#line 6836 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 324:
#line 1068 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[-1].ival);}
#line 6843 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 325:
#line 1071 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qrfactorization = (yyvsp[-1].ival);}
#line 6849 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 326:
#line 1073 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nEig = (yyvsp[-1].ival); }
#line 6855 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 327:
#line 1075 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;}
#line 6861 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 328:
#line 1077 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setSubSpaceInfo((yyvsp[-3].ival),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6867 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 329:
#line 1079 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subspaceSize = (yyvsp[-1].ival);}
#line 6873 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 330:
#line 1081 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolEig = (yyvsp[-1].fval); }
#line 6879 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 331:
#line 1083 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolJac = (yyvsp[-1].fval); }
#line 6885 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 332:
#line 1085 "p.y" /* yacc.c:1646  */
    { domain->solInfo().explicitK = true; }
#line 6891 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 333:
#line 1087 "p.y" /* yacc.c:1646  */
    { geoSource->setShift((yyvsp[-1].fval)); }
#line 6897 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 334:
#line 1089 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; }
#line 6903 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 335:
#line 1091 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[-1].strval); }
#line 6910 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 336:
#line 1094 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[-2].strval); 
          domain->solInfo().arpack_mode = (yyvsp[-1].ival); }
#line 6918 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 337:
#line 1098 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[-2].fval), int((yyvsp[-1].fval))); }
#line 6925 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 338:
#line 1101 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival));}
#line 6932 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 339:
#line 1104 "p.y" /* yacc.c:1646  */
    { domain->solInfo().filtereig = bool((yyvsp[-1].ival)); }
#line 6938 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 340:
#line 1106 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverSubType = (yyvsp[-1].ival); }
#line 6944 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 341:
#line 1108 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;}
#line 6951 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 342:
#line 1111 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxitEig = (yyvsp[-1].ival); }
#line 6957 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 343:
#line 1113 "p.y" /* yacc.c:1646  */
    { domain->solInfo().test_ulrich = true; }
#line 6963 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 344:
#line 1115 "p.y" /* yacc.c:1646  */
    { domain->solInfo().addedMass = (yyvsp[-1].ival); }
#line 6969 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 345:
#line 1119 "p.y" /* yacc.c:1646  */
    { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = (yyvsp[-1].strval);
          domain->solInfo().printMatLabExit = true; }
#line 6977 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 346:
#line 1123 "p.y" /* yacc.c:1646  */
    { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = (yyvsp[-2].strval); 
          domain->solInfo().printMatLabExit = true; }
#line 6985 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 347:
#line 1129 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elementDeletion = true; }
#line 6991 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 349:
#line 1134 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[-1].fval); domain->solInfo().deleteElements.insert(std::pair<int,double>((yyvsp[0].ival)-1,(yyvsp[-1].fval))); }
#line 6997 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 350:
#line 1136 "p.y" /* yacc.c:1646  */
    { domain->solInfo().deleteElements.insert(std::pair<int,double>((yyvsp[0].ival)-1,(yyvsp[-1].fval))); }
#line 7003 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 351:
#line 1140 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sloshing = 1; }
#line 7009 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 352:
#line 1142 "p.y" /* yacc.c:1646  */
    { domain->setGravitySloshing((yyvsp[-1].fval)); }
#line 7015 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 353:
#line 1146 "p.y" /* yacc.c:1646  */
    { domain->solInfo().massFlag = 1; }
#line 7021 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 354:
#line 1148 "p.y" /* yacc.c:1646  */
    { domain->solInfo().massFlag = 1;
          domain->solInfo().massFile = std::string((yyvsp[-1].strval)); }
#line 7028 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 355:
#line 1153 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7035 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 356:
#line 1156 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);}
#line 7041 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 357:
#line 1160 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Top); }
#line 7047 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 358:
#line 1164 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal_id.push_back((yyvsp[0].ival)); }
#line 7053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 359:
#line 1166 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal_id.push_back((yyvsp[0].ival)); }
#line 7059 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 360:
#line 1168 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contact_modal_id.push_back((yyvsp[0].ival)); }
#line 7065 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 361:
#line 1171 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Dynamic); }
#line 7071 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 365:
#line 1176 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
#line 7077 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 366:
#line 1178 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; }
#line 7083 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 367:
#line 1180 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-1].ival); }
#line 7089 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 368:
#line 1182 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-1].ival); }
#line 7095 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 369:
#line 1184 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-5].ival);
          domain->solInfo().stable_cfl = (yyvsp[-4].fval);
          domain->solInfo().stable_tol = (yyvsp[-3].fval);
          domain->solInfo().stable_maxit = (yyvsp[-2].ival);
          domain->solInfo().stable_freq = (yyvsp[-1].ival);
        }
#line 7106 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 370:
#line 1191 "p.y" /* yacc.c:1646  */
    { domain->solInfo().iacc_switch = bool((yyvsp[-1].ival)); }
#line 7112 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 371:
#line 1193 "p.y" /* yacc.c:1646  */
    { domain->solInfo().zeroRot = bool((yyvsp[-1].ival)); }
#line 7118 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 372:
#line 1195 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_secondary = true; }
#line 7124 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 373:
#line 1197 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-1].ival); }
#line 7130 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 374:
#line 1199 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-3].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-2].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-1].fval); }
#line 7138 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 375:
#line 1203 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-4].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-3].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-2].fval);
          domain->solInfo().tdenforceInitia = (yyvsp[-1].fval); }
#line 7147 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 376:
#line 1208 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-5].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-4].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-3].fval);
          domain->solInfo().tdenforceInitia = (yyvsp[-2].fval);
          domain->solInfo().tdenforceFinal = (yyvsp[-1].fval); }
#line 7157 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 377:
#line 1214 "p.y" /* yacc.c:1646  */
    { domain->solInfo().check_energy_balance = true; }
#line 7163 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 378:
#line 1216 "p.y" /* yacc.c:1646  */
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[-2].fval); 
          domain->solInfo().epsilon2 = (yyvsp[-1].fval); }
#line 7171 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 379:
#line 1220 "p.y" /* yacc.c:1646  */
    { domain->solInfo().quasistatic = bool((yyvsp[-1].ival)); }
#line 7177 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 380:
#line 1224 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ConwepOnOff = true;
          BlastLoading::InputFileData = (yyvsp[-1].blastData); }
#line 7184 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 381:
#line 1229 "p.y" /* yacc.c:1646  */
    { // Note: chargeWeight must be entered in the units of mass of the problem, not units of force.
          (yyval.blastData).ExplosivePosition[0] = (yyvsp[-5].fval);
          (yyval.blastData).ExplosivePosition[1] = (yyvsp[-4].fval);
          (yyval.blastData).ExplosivePosition[2] = (yyvsp[-3].fval);
          (yyval.blastData).ExplosiveDetonationTime = (yyvsp[-1].fval);
          ((yyvsp[0].ival) == 0.0? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst);
          (yyval.blastData).ScaleLength = 1.0;
          (yyval.blastData).ScaleTime = 1.0;
          (yyval.blastData).ScaleMass = 1.0;
          (yyval.blastData).ExplosiveWeight = (yyvsp[-2].fval)*2.2; // The 2.2 factor is to convert from kilograms to pounds force.
          (yyval.blastData).ExplosiveWeightCubeRoot = pow((yyval.blastData).ExplosiveWeight,1.0/3.0);
        }
#line 7201 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 382:
#line 1243 "p.y" /* yacc.c:1646  */
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; }
#line 7207 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 384:
#line 1246 "p.y" /* yacc.c:1646  */
    { domain->solInfo().acoustic = true; }
#line 7213 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 386:
#line 1251 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-1].fval),(yyvsp[0].fval)); }
#line 7219 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 387:
#line 1253 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval),(yyvsp[0].fval)); }
#line 7225 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 388:
#line 1255 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[0].fval)); }
#line 7231 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 389:
#line 1257 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-2].fval),(yyvsp[-1].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[0].fval); }
#line 7239 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 390:
#line 1263 "p.y" /* yacc.c:1646  */
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
            domain->solInfo().order = 1;
          }
          else 
            domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo((yyvsp[0].fval)); 
        }
#line 7252 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 391:
#line 1274 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; }
#line 7259 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 394:
#line 1279 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
#line 7265 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 395:
#line 1281 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back((yyvsp[-1].ival)); }
#line 7271 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 396:
#line 1285 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setQuasistaticInfo((yyvsp[-3].fval), 0, (yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7277 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 397:
#line 1287 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setQuasistaticInfo((yyvsp[-4].fval), 0, (yyvsp[-3].fval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 7283 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 398:
#line 1291 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7290 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 399:
#line 1301 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-1].ival)); 
          domain->solInfo().isCollocated = 0; }
#line 7297 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 400:
#line 1304 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-1].ival)); 
          domain->solInfo().isCollocated = 0;
          if((yyvsp[-1].ival) == 20 || (yyvsp[-1].ival) == 22) { // set default alphas for C0
            domain->solInfo().alphas[0] = 0.5+0.375;
            domain->solInfo().alphas[1] = -0.375;
            domain->solInfo().alphasv    = 0.0;
          }
        }
#line 7310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 401:
#line 1313 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-3].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[-3].ival) == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = (yyvsp[-2].fval);
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = (yyvsp[-2].fval)+(yyvsp[-1].fval);
            domain->solInfo().alphas[1] = -(yyvsp[-1].fval);
            domain->solInfo().alphasv    = 0.0;
          }
        }
#line 7331 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 402:
#line 1330 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-4].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[-4].ival) == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = (yyvsp[-3].fval);
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = (yyvsp[-3].fval)+(yyvsp[-2].fval);
            domain->solInfo().alphas[1] = -(yyvsp[-2].fval);

            // This option added by Alex Main.  The last number on this line is the value used for velocity prediciton,
            // which is useful for embedded simulations.
            domain->solInfo().alphasv    = (yyvsp[-1].fval);
          }
        }
#line 7355 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 403:
#line 1351 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-2].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[-1].fval);
        }
#line 7364 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 404:
#line 1356 "p.y" /* yacc.c:1646  */
    { domain->solInfo().isCollocated = (yyvsp[-1].ival); }
#line 7370 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 405:
#line 1358 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subcycle = (yyvsp[-1].ival); }
#line 7376 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 406:
#line 1360 "p.y" /* yacc.c:1646  */
    { geoSource->setMatch((yyvsp[-1].strval)); }
#line 7382 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 407:
#line 1362 "p.y" /* yacc.c:1646  */
    {}
#line 7388 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 408:
#line 1366 "p.y" /* yacc.c:1646  */
    {}
#line 7394 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 409:
#line 1368 "p.y" /* yacc.c:1646  */
    { domain->AddAeroEmbedSurfaceId((yyvsp[0].ival)); }
#line 7400 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 410:
#line 1372 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAeroHeat((yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 7406 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 411:
#line 1376 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setThermoh(1); }
#line 7412 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 412:
#line 1380 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setThermoe(1); }
#line 7418 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 413:
#line 1384 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setModeDecomp(1); }
#line 7424 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 414:
#line 1386 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setModeDecomp(1, (yyvsp[-1].ival)); }
#line 7430 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 415:
#line 1390 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hzemFlag=1; }
#line 7436 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 416:
#line 1394 "p.y" /* yacc.c:1646  */
    { domain->solInfo().slzemFlag=1; }
#line 7442 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 417:
#line 1398 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setTrbm((yyvsp[-1].fval)); }
#line 7448 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 418:
#line 1402 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 7454 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 419:
#line 1404 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-1].fval)); }
#line 7460 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 420:
#line 1406 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm(); }
#line 7466 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 421:
#line 1408 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-3].fval),(yyvsp[-2].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[-1].ival)); }
#line 7473 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 422:
#line 1411 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-2].fval),(yyvsp[-1].fval));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-6].fval); grbm_ref[1] = (yyvsp[-5].fval); grbm_ref[2] = (yyvsp[-4].fval);
        }
#line 7482 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 423:
#line 1416 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-1].fval));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-5].fval); grbm_ref[1] = (yyvsp[-4].fval); grbm_ref[2] = (yyvsp[-3].fval);
        }
#line 7491 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 424:
#line 1421 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm();
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-3].fval); grbm_ref[1] = (yyvsp[-2].fval); grbm_ref[2] = (yyvsp[-1].fval);
        }
#line 7500 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 425:
#line 1426 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-3].fval),(yyvsp[-2].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[-1].ival));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-7].fval); grbm_ref[1] = (yyvsp[-6].fval); grbm_ref[2] = (yyvsp[-5].fval);
        }
#line 7510 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 426:
#line 1434 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modeFilterFlag = (yyvsp[-1].ival); }
#line 7516 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 427:
#line 1436 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modeFilterFlag = 1; }
#line 7522 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 428:
#line 1440 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter((yyvsp[-1].ival)); }
#line 7528 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 429:
#line 1442 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter(1); }
#line 7534 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 430:
#line 1444 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter((yyvsp[-2].ival));
          domain->solInfo().filterQ = (yyvsp[-1].ival); }
#line 7541 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 432:
#line 1450 "p.y" /* yacc.c:1646  */
    { if((yyvsp[0].ival) < 1) {
        fprintf(stderr, " *** ERROR: mode %d specified under RBMFILTER is invalid.\n", (yyvsp[0].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters.insert((yyvsp[0].ival)-1);
    }
#line 7553 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 433:
#line 1458 "p.y" /* yacc.c:1646  */
    { if((yyvsp[0].ival) < 1) {
        fprintf(stderr, " *** ERROR: mode %d specified under RBMFILTER is invalid.\n", (yyvsp[0].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters.insert((yyvsp[0].ival)-1);
    }
#line 7565 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 434:
#line 1468 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hzemFilterFlag=1; }
#line 7571 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 435:
#line 1472 "p.y" /* yacc.c:1646  */
    { domain->solInfo().slzemFilterFlag=1; }
#line 7577 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 436:
#line 1476 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setTimes((yyvsp[-1].fval),(yyvsp[-2].fval),(yyvsp[-3].fval)); }
#line 7583 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 437:
#line 1480 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[-2].ival),(yyvsp[-1].ival),1);
        }
#line 7592 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 438:
#line 1486 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[-3].ival),(yyvsp[-2].ival),(yyvsp[-1].ival));
        }
#line 7601 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 439:
#line 1491 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[-4].ival),(yyvsp[-3].ival),(yyvsp[-2].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        }
#line 7612 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 442:
#line 1504 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaNoForce = true; }
#line 7618 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 443:
#line 1506 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[0].ival); }
#line 7624 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 444:
#line 1508 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaLocalBasisImprovement = 1; }
#line 7630 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 445:
#line 1510 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaTimeReversible = true; }
#line 7636 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 446:
#line 1512 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaRemoteCoarse = true; }
#line 7642 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 447:
#line 1514 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaProjTol = (yyvsp[0].fval); }
#line 7648 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 448:
#line 1516 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaReadInitSeed = true; }
#line 7654 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 449:
#line 1518 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpCvgRatio = 0.0; }
#line 7660 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 450:
#line 1520 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[0].fval); }
#line 7666 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 451:
#line 1522 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpMagnOutput = true; }
#line 7672 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 452:
#line 1526 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-3].ival) == 1) domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
          else return -1; // only RAYDAMP is allowed here
        }
#line 7680 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 453:
#line 1530 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-4].ival) == 1) {
            domain->solInfo().setDamping((yyvsp[-3].fval),(yyvsp[-2].fval));
            domain->solInfo().mtypeDamp = (int)(yyvsp[-1].ival);
          }
          else return -1;
        }
#line 7691 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 454:
#line 1537 "p.y" /* yacc.c:1646  */
    { if(geoSource->setModalDamping((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true;
        }
#line 7699 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 455:
#line 1543 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 7705 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 456:
#line 1547 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7715 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 457:
#line 1553 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-2].fval), (yyvsp[-1].fval), 0.0);
        }
#line 7725 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 458:
#line 1561 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7735 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 459:
#line 1567 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[-1].ival),0.0,0.0,0.0);
        }
#line 7744 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 461:
#line 1575 "p.y" /* yacc.c:1646  */
    {
           domain->setWaveDirections((yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7752 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 462:
#line 1581 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7762 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 463:
#line 1589 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7768 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 464:
#line 1591 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7774 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 465:
#line 1593 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
#line 7780 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 466:
#line 1595 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
#line 7786 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 467:
#line 1597 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
#line 7796 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 469:
#line 1613 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-4].ival)-1;
          surf_bc[0].type = (BCond::BCType) (yyvsp[-3].ival); //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = (yyvsp[-2].ival)-1;
          surf_bc[0].val = (yyvsp[-1].ival)-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
          if(geoSource->getNumSurfaceConstraint() > 1) delete [] surf_bc;
        }
#line 7809 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 470:
#line 1623 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[0].bclist); }
#line 7815 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 471:
#line 1627 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7821 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 472:
#line 1629 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7827 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 473:
#line 1633 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 7833 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 474:
#line 1635 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-1].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; }
#line 7839 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 475:
#line 1639 "p.y" /* yacc.c:1646  */
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; 

          int* allPDirNodes = new int[(yyvsp[0].bclist)->n];

          for (int ii=0; ii < (yyvsp[0].bclist)->n; ii++)
            allPDirNodes[ii]=((yyvsp[0].bclist)->d[ii]).nnum;
          std::sort(allPDirNodes,allPDirNodes + (yyvsp[0].bclist)->n);

          int maxFSNodes = 32;
          int* allHEVFSNodes = new int[maxFSNodes];
          allHEVFSNodes[0] = allPDirNodes[0];

          int numHEVFSNodes = 1;
          for (int ii = 1; ii < (yyvsp[0].bclist)->n; ++ii)  {
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
          
          (yyval.bclist) = new BCList;

          for (int ii=0; ii<numHEVFSNodes; ii++)
            (yyval.bclist)->add(allHEVFSNodes[ii],10,0.0); 
          delete [] allHEVFSNodes;

          for(int i=0; i<(yyval.bclist)->n; ++i) (yyval.bclist)->d[i].type = BCond::Hefrs;
        }
#line 7882 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 476:
#line 1680 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[0].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[0].bclist)->d)[ii]); }
#line 7890 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 477:
#line 1684 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); 
          for (int ii = 0; ii < (yyvsp[0].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[0].bclist)->d)[ii]); }
#line 7898 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 478:
#line 1690 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList;
          for(int i=0; i<(yyvsp[-1].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[-1].nl).nd[i],10,0.0); } }
#line 7906 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 479:
#line 1696 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;}
#line 7912 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 480:
#line 1698 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-3].bclist); BCond bc; bc.nnum = (yyvsp[-2].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-1].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); }
#line 7919 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 481:
#line 1701 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-4].ival); i<=(yyvsp[-2].ival); ++i) { BCond bc; bc.setData(i-1, 6, (yyvsp[-1].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } }
#line 7925 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 482:
#line 1703 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); i+=(yyvsp[-2].ival)) { BCond bc; bc.setData(i-1, 6, (yyvsp[-1].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } }
#line 7931 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 483:
#line 1705 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-2].ival)-1;
          surf_bc[0].val = (yyvsp[-1].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
#line 7942 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 484:
#line 1714 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7948 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 485:
#line 1716 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7954 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 486:
#line 1718 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-3].bclist); BCond bc; bc.nnum = (yyvsp[-2].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-1].fval); bc.type = BCond::Flux; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); }
#line 7961 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 487:
#line 1721 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-2].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[-1].fval);
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 7974 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 488:
#line 1732 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7980 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 489:
#line 1734 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7986 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 490:
#line 1736 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-5].bclist); BCond bc; bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval); bc.type = BCond::Convection; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); }
#line 7993 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 491:
#line 1739 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-4].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval);
          surf_bc[0].type = BCond::Convection;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 8006 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 492:
#line 1750 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 8012 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 493:
#line 1752 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 8018 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 494:
#line 1754 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-5].bclist); BCond bc; bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval); bc.type = BCond::Radiation; (yyval.bclist)->add(bc); }
#line 8025 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 498:
#line 1766 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8031 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 499:
#line 1768 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)); }
#line 8037 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 500:
#line 1770 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8043 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 503:
#line 1778 "p.y" /* yacc.c:1646  */
    { domain->addSommerElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        }
#line 8051 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 504:
#line 1784 "p.y" /* yacc.c:1646  */
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[0].ival)-1; }
#line 8057 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 505:
#line 1786 "p.y" /* yacc.c:1646  */
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[0].ival)-1; (yyval.nl).num++; }
#line 8064 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 509:
#line 1798 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8071 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 510:
#line 1801 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)); }
#line 8078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 511:
#line 1804 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8085 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 514:
#line 1813 "p.y" /* yacc.c:1646  */
    { domain->addScatterElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
          domain->addNeumElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); }
#line 8092 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 517:
#line 1822 "p.y" /* yacc.c:1646  */
    { domain->addNeumElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); }
#line 8098 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 520:
#line 1830 "p.y" /* yacc.c:1646  */
    { domain->addWetElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; }
#line 8106 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 523:
#line 1840 "p.y" /* yacc.c:1646  */
    { domain->addScatterElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);}
#line 8112 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 524:
#line 1844 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 8118 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 525:
#line 1848 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8124 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 526:
#line 1850 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8130 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 527:
#line 1854 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[0].bclist); }
#line 8136 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 528:
#line 1858 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[0].bclist); }
#line 8142 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 529:
#line 1860 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[0].bclist); }
#line 8148 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 530:
#line 1864 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDARBFlag = (yyvsp[-1].fval);}
#line 8154 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 532:
#line 1869 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDDNBVal = (yyvsp[-1].fval);}
#line 8160 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 534:
#line 1874 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDROBVal = (yyvsp[-3].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[-2].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[-1].fval);}
#line 8168 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 536:
#line 1881 "p.y" /* yacc.c:1646  */
    { domain->setFFP((yyvsp[-1].ival)); }
#line 8174 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 537:
#line 1883 "p.y" /* yacc.c:1646  */
    { domain->setFFP((yyvsp[-2].ival),(yyvsp[-1].ival)); }
#line 8180 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 538:
#line 1887 "p.y" /* yacc.c:1646  */
    {
           domain->setFFP((yyvsp[-2].ival));
        }
#line 8188 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 539:
#line 1893 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-4].ival)] = ModalParams(ModalParams::Eigen, (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8194 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 540:
#line 1895 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-4].ival)] = ModalParams((yyvsp[-3].mpt), (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8200 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 541:
#line 1897 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-5].ival)] = ModalParams(ModalParams::Eigen, (yyvsp[-3].strval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 8206 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 542:
#line 1899 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-5].ival)] = ModalParams((yyvsp[-4].mpt), (yyvsp[-3].strval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 8212 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 543:
#line 1903 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, (yyvsp[-1].strval)); }
#line 8218 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 544:
#line 1905 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8224 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 547:
#line 1909 "p.y" /* yacc.c:1646  */
    { domain->solInfo().adjointMap[(OutputInfo::Type)(yyvsp[-1].ival)] = domain->solInfo().readInAdjointROB.size();
          domain->solInfo().readInAdjointROB.push_back((yyvsp[-3].strval));
          domain->solInfo().maxSizeAdjointBasis.push_back((yyvsp[-2].ival)); }
#line 8232 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 548:
#line 1915 "p.y" /* yacc.c:1646  */
    { }
#line 8238 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 549:
#line 1917 "p.y" /* yacc.c:1646  */
    { domain->solInfo().zeroInitialDisp = 1; }
#line 8244 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 550:
#line 1919 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)  return -1; }
#line 8251 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 551:
#line 1922 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
#line 8259 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 552:
#line 1926 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().idis_modal_id = (yyvsp[-2].ival); }
#line 8268 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 553:
#line 1933 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[-1].fval);  }
#line 8274 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 554:
#line 1935 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; amplitude = 1.0; }
#line 8280 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 555:
#line 1937 "p.y" /* yacc.c:1646  */
    { BCond bc; /* add 6 boundary conditions */
          bc.type = BCond::Idisp6;
          bc.nnum = (yyvsp[-7].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[-6].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[-5].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[-4].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = amplitude*(yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = amplitude*(yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = amplitude*(yyvsp[-1].fval); (yyval.bclist)->add(bc);
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        }
#line 8295 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 556:
#line 1948 "p.y" /* yacc.c:1646  */
    { BCond bc; /* add 6 boundary conditions */
          bc.type = BCond::Idisp6;
          bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[-1].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = 0.0         ; (yyval.bclist)->add(bc);
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        }
#line 8310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 557:
#line 1959 "p.y" /* yacc.c:1646  */
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); }
#line 8317 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 558:
#line 1962 "p.y" /* yacc.c:1646  */
    { domain->solInfo().buckling = 1; }
#line 8323 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 559:
#line 1967 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[-1].ival); }
#line 8329 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 560:
#line 1969 "p.y" /* yacc.c:1646  */
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[-7].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[-6].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[-5].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[-4].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[-1].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIDis6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        }
#line 8343 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 561:
#line 1982 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[-1].ival); }
#line 8349 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 562:
#line 1984 "p.y" /* yacc.c:1646  */
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[-7].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[-6].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[-5].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[-4].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[-3].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[-2].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[-1].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIVel6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        }
#line 8363 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 563:
#line 1996 "p.y" /* yacc.c:1646  */
    { }
#line 8369 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 564:
#line 1998 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8376 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 565:
#line 2001 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
#line 8384 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 566:
#line 2005 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().ivel_modal_id = (yyvsp[-2].ival); }
#line 8393 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 567:
#line 2012 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8400 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 568:
#line 2015 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true; }
#line 8408 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 569:
#line 2021 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGEPS();
          for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Etemperatures;
          if(geoSource->setIDis6((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8416 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 570:
#line 2027 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 8422 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 571:
#line 2029 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 8428 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 572:
#line 2031 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Forces; (yyvsp[0].bcval).loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8434 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 573:
#line 2033 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } }
#line 8440 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 574:
#line 2035 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } }
#line 8446 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 575:
#line 2037 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[-1].ival)); (yyval.bclist)->add(bc); } }
#line 8452 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 576:
#line 2039 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-8].ival); i<=(yyvsp[-6].ival); i+=(yyvsp[-4].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[-3].ival)); (yyval.bclist)->add(bc); } }
#line 8458 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 577:
#line 2041 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 8469 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 578:
#line 2050 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8476 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 579:
#line 2053 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Forces; (yyvsp[0].bclist)->d[i].loadsetid = (yyvsp[-4].ival); }
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8483 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 580:
#line 2056 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8490 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 581:
#line 2059 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Forces; (yyvsp[0].bclist)->d[i].loadsetid = (yyvsp[-5].ival); }
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8497 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 582:
#line 2064 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8503 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 583:
#line 2066 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8509 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 584:
#line 2068 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); }}
#line 8515 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 585:
#line 2070 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-6].bclist); for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); }}
#line 8521 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 586:
#line 2072 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); } }
#line 8527 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 587:
#line 2074 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-8].bclist); for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); } }
#line 8533 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 588:
#line 2078 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8539 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 589:
#line 2080 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8545 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 590:
#line 2084 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8551 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 591:
#line 2086 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8557 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 594:
#line 2094 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMTT((yyval.ymtt));}
#line 8563 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 595:
#line 2096 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8569 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 596:
#line 2098 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMTT((yyval.ymtt));}
#line 8575 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 599:
#line 2106 "p.y" /* yacc.c:1646  */
    { (yyval.ctett) = new MFTTData((yyvsp[-4].ival)); (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addCTETT((yyval.ctett));}
#line 8581 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 600:
#line 2108 "p.y" /* yacc.c:1646  */
    { (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8587 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 601:
#line 2110 "p.y" /* yacc.c:1646  */
    { (yyval.ctett) = new MFTTData((yyvsp[-4].ival)); (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addCTETT((yyval.ctett));}
#line 8593 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 604:
#line 2118 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSS1DT((yyval.ymtt));}
#line 8599 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 605:
#line 2120 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8605 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 606:
#line 2122 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSS1DT((yyval.ymtt));}
#line 8611 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 609:
#line 2130 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-6].ival), (yyvsp[-4].dlist), true); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8617 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 610:
#line 2132 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-8].ival), (yyvsp[-4].dlist), bool((yyvsp[-6].ival))); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8623 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 611:
#line 2134 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); }
#line 8629 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 612:
#line 2136 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-6].ival), (yyvsp[-4].dlist), true); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8635 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 613:
#line 2138 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-8].ival), (yyvsp[-4].dlist), bool((yyvsp[-6].ival))); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8641 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 616:
#line 2146 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSST((yyval.ymtt));}
#line 8647 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 617:
#line 2148 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8653 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 618:
#line 2150 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSST((yyval.ymtt));}
#line 8659 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 621:
#line 2158 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSSRT((yyval.ymtt));}
#line 8665 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 622:
#line 2160 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8671 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 623:
#line 2162 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSSRT((yyval.ymtt));}
#line 8677 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 626:
#line 2170 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMST((yyval.ymtt));}
#line 8683 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 627:
#line 2172 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8689 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 628:
#line 2174 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMST((yyval.ymtt));}
#line 8695 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 631:
#line 2182 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft) = new MFTTData((yyvsp[-4].ival)); (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSDETAFT((yyval.sdetaft));}
#line 8701 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 632:
#line 2184 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8707 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 633:
#line 2186 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft) = new MFTTData((yyvsp[-4].ival)); (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSDETAFT((yyval.sdetaft));}
#line 8713 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 636:
#line 2194 "p.y" /* yacc.c:1646  */
    { 
#ifdef USE_EIGEN3
          (yyval.rubdaft) = new GenMFTTData<Eigen::Vector4d>((yyvsp[-7].ival)); (yyval.rubdaft)->add((yyvsp[-5].fval), Eigen::Vector4d((yyvsp[-4].fval),(yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval))); domain->addRUBDAFT((yyval.rubdaft));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
#line 8725 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 637:
#line 2202 "p.y" /* yacc.c:1646  */
    {
#ifdef USE_EIGEN3
          (yyval.rubdaft)->add((yyvsp[-5].fval), Eigen::Vector4d((yyvsp[-4].fval),(yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval)));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
#line 8737 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 638:
#line 2210 "p.y" /* yacc.c:1646  */
    {
#ifdef USE_EIGEN3
          (yyval.rubdaft) = new GenMFTTData<Eigen::Vector4d>((yyvsp[-7].ival)); (yyval.rubdaft)->add((yyvsp[-5].fval), Eigen::Vector4d((yyvsp[-4].fval),(yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval))); domain->addRUBDAFT((yyval.rubdaft));
#else
          std::cerr << " *** ERROR: RUBDAFT command requires AERO-S configured with Eigen library. Exiting...\n"; exit(-1);
#endif
        }
#line 8749 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 641:
#line 2222 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xLMPCFactor = (yyvsp[-3].fval);
          domain->solInfo().yLMPCFactor = (yyvsp[-2].fval);
          domain->solInfo().zLMPCFactor = (yyvsp[-1].fval); }
#line 8757 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 642:
#line 2228 "p.y" /* yacc.c:1646  */
    { domain->solInfo().localDualBasisSize.push_back((yyvsp[-1].ival));
          domain->solInfo().modalLMPC = true; }
#line 8764 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 643:
#line 2231 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackROMLMPCVec((yyvsp[-1].fval)); }
#line 8770 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 644:
#line 2235 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = (yyvsp[-1].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[0].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
#line 8778 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 645:
#line 2239 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons)->addterm((yyvsp[0].mpcterm)); }
#line 8784 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 646:
#line 2241 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = (yyvsp[-1].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[0].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
#line 8792 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 647:
#line 2247 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8799 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 648:
#line 2250 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-2].ival), (yyvsp[-1].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8806 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 649:
#line 2253 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-4].ival), (yyvsp[-3].fval));
          (yyval.lmpcons)->type = (yyvsp[-1].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8814 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 650:
#line 2257 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-3].ival), (yyvsp[-2].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[-1].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8823 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 651:
#line 2262 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-5].ival), (yyvsp[-4].fval));
          (yyval.lmpcons)->type = (yyvsp[-2].ival);
          (yyval.lmpcons)->lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[-1].copt).penalty;
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8833 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 652:
#line 2270 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].fval) == 0.0) {
            fprintf(stderr," *** WARNING: zero coefficient in LMPC\n");
            fprintf(stderr," ***          node %d dof %d\n",(yyvsp[-3].ival),(yyvsp[-2].ival));
          }
          (yyval.mpcterm) = new LMPCTerm();
          (yyval.mpcterm)->nnum = (yyvsp[-3].ival)-1;
          (yyval.mpcterm)->dofnum = (yyvsp[-2].ival)-1;
          (yyval.mpcterm)->coef.r_value = (yyvsp[-1].fval);
          if((yyvsp[-2].ival) == 1){
            (yyval.mpcterm)->coef.r_value /= domain->solInfo().xLMPCFactor;
          } else if ((yyvsp[-2].ival) == 2) {
            (yyval.mpcterm)->coef.r_value /= domain->solInfo().yLMPCFactor;
          } else if ((yyvsp[-2].ival) == 3) {
            (yyval.mpcterm)->coef.r_value /= domain->solInfo().zLMPCFactor;
          }
        }
#line 8854 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 655:
#line 2293 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].cxbcval).nnum,(yyvsp[-1].cxbcval).reval,(yyvsp[-1].cxbcval).imval,(yyvsp[0].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
#line 8860 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 656:
#line 2295 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons)->addterm((yyvsp[0].mpcterm)); }
#line 8866 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 657:
#line 2297 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].cxbcval).nnum,(yyvsp[-1].cxbcval).reval,(yyvsp[-1].cxbcval).imval,(yyvsp[0].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
#line 8872 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 658:
#line 2301 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-4].ival); (yyval.cxbcval).reval=(yyvsp[-2].fval); (yyval.cxbcval).imval=(yyvsp[-1].fval); }
#line 8878 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 659:
#line 2303 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-2].ival); (yyval.cxbcval).reval=(yyvsp[-1].fval); (yyval.cxbcval).imval=0.0; }
#line 8884 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 660:
#line 2305 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-1].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; }
#line 8890 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 661:
#line 2309 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-2].fval)==0.0) && ((yyvsp[-1].fval)==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[-4].ival),(yyvsp[-3].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[-4].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[-3].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[-2].fval), (yyvsp[-1].fval)); }
        }
#line 8902 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 662:
#line 2317 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].fval)==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[-3].ival),(yyvsp[-2].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[-3].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[-2].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[-1].fval),0.0); }
        }
#line 8914 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 663:
#line 2327 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 8920 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 664:
#line 2329 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].cxbclist)->n; ++i) (yyvsp[0].cxbclist)->d[i].loadsetid = (yyvsp[-2].ival);
          (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 8927 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 665:
#line 2334 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[0].cxbcval)); }
#line 8933 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 666:
#line 2336 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[-1].cxbclist); (yyval.cxbclist)->add((yyvsp[0].cxbcval)); }
#line 8939 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 669:
#line 2344 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
	  sp.A = (yyvsp[-14].fval);  sp.E = (yyvsp[-13].fval);  sp.nu  = (yyvsp[-12].fval);  sp.rho = (yyvsp[-11].fval);
          sp.c = (yyvsp[-10].fval);  sp.k = (yyvsp[-9].fval);  sp.eh  = (yyvsp[-8].fval);  sp.P   = (yyvsp[-7].fval);  sp.Ta  = (yyvsp[-6].fval); 
          sp.Q = (yyvsp[-5].fval); sp.W = (yyvsp[-4].fval); sp.Ixx = (yyvsp[-3].fval); sp.Iyy = (yyvsp[-2].fval); sp.Izz = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 8950 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 670:
#line 2351 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-17].fval);  sp.E = (yyvsp[-16].fval);  sp.nu  = (yyvsp[-15].fval);  sp.rho = (yyvsp[-14].fval);
          sp.c = (yyvsp[-13].fval);  sp.k = (yyvsp[-12].fval);  sp.eh  = (yyvsp[-11].fval);  sp.P   = (yyvsp[-10].fval);  sp.Ta  = (yyvsp[-9].fval);
          sp.Q = (yyvsp[-8].fval); sp.W = (yyvsp[-7].fval); sp.Ixx = (yyvsp[-6].fval); sp.Iyy = (yyvsp[-5].fval); sp.Izz = (yyvsp[-4].fval);
          switch((yyvsp[-3].ival)) {
            case 1 : // RAYDAMP
              sp.betaDamp = (yyvsp[-2].fval); sp.alphaDamp = (yyvsp[-1].fval);
              break;
            case 2 : // STRDAMP
              sp.etaDamp = (yyvsp[-2].fval); sp.betaDamp = (yyvsp[-1].fval);
              break;
            case 3 : // RUBDAMP
              sp.eta_E = (yyvsp[-2].fval);
              if(sp.eta_E >= 0) sp.eta_mu = (yyvsp[-1].fval);
              sp.E0 = (yyvsp[-16].fval);
              sp.mu0 = (yyvsp[-16].fval)/(2*(1+(yyvsp[-15].fval)));
              break;
            default :
              return -1;
          }
          geoSource->addMat( (yyvsp[-18].ival)-1, sp );
        }
#line 8977 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 671:
#line 2374 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
	  sp.A = (yyvsp[-18].fval);  sp.E = (yyvsp[-17].fval);  sp.nu  = (yyvsp[-16].fval);  sp.rho = (yyvsp[-15].fval);
          sp.c = (yyvsp[-14].fval);  sp.k = (yyvsp[-13].fval);  sp.eh  = (yyvsp[-12].fval);  sp.P   = (yyvsp[-11].fval);  sp.Ta  = (yyvsp[-10].fval); 
          sp.Q = (yyvsp[-9].fval); sp.W = (yyvsp[-8].fval); sp.Ixx = (yyvsp[-7].fval); sp.Iyy = (yyvsp[-6].fval); sp.Izz = (yyvsp[-5].fval);
	  sp.ymin = (yyvsp[-4].fval); sp.ymax = (yyvsp[-3].fval); sp.zmin = (yyvsp[-2].fval); sp.zmax = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-19].ival)-1, sp );
        }
#line 8989 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 672:
#line 2382 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-21].fval);  sp.E = (yyvsp[-20].fval);  sp.nu  = (yyvsp[-19].fval);  sp.rho = (yyvsp[-18].fval);
          sp.c = (yyvsp[-17].fval);  sp.k = (yyvsp[-16].fval);  sp.eh  = (yyvsp[-15].fval);  sp.P   = (yyvsp[-14].fval);  sp.Ta  = (yyvsp[-13].fval);
          sp.Q = (yyvsp[-12].fval); sp.W = (yyvsp[-11].fval); sp.Ixx = (yyvsp[-10].fval); sp.Iyy = (yyvsp[-9].fval); sp.Izz = (yyvsp[-8].fval);
          sp.ymin = (yyvsp[-7].fval); sp.ymax = (yyvsp[-6].fval); sp.zmin = (yyvsp[-5].fval); sp.zmax = (yyvsp[-4].fval);
          switch((yyvsp[-3].ival)) {
            case 1 : // RAYDAMP
              sp.betaDamp = (yyvsp[-2].fval); sp.alphaDamp = (yyvsp[-1].fval);
              break;
            case 2 : // STRDAMP
              sp.etaDamp = (yyvsp[-2].fval); sp.betaDamp = (yyvsp[-1].fval);
              break;
            case 3 : // RUBDAMP
              sp.eta_E = (yyvsp[-2].fval); 
              if(sp.eta_E >= 0) sp.eta_mu = (yyvsp[-1].fval);
              sp.E0 = (yyvsp[-20].fval);
              sp.mu0 = (yyvsp[-20].fval)/(2*(1+(yyvsp[-19].fval)));
              break;
            default :
              return -1;
          }
          geoSource->addMat( (yyvsp[-22].ival)-1, sp );
        }
#line 9017 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 673:
#line 2406 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-7].fval); sp.E = (yyvsp[-6].fval); sp.nu = (yyvsp[-5].fval); sp.rho = (yyvsp[-4].fval);
          sp.c = (yyvsp[-3].fval); sp.k = (yyvsp[-2].fval); sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-8].ival)-1, sp ); 
        }
#line 9027 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 674:
#line 2412 "p.y" /* yacc.c:1646  */
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[-12].fval);  sp.E = (yyvsp[-11].fval);  sp.nu  = (yyvsp[-10].fval);  sp.rho = (yyvsp[-9].fval);
          sp.c = (yyvsp[-8].fval);  sp.k = (yyvsp[-7].fval);  sp.eh  = (yyvsp[-6].fval);  sp.P   = (yyvsp[-5].fval);  sp.Ta  = (yyvsp[-4].fval);
          sp.Q = (yyvsp[-3].fval); sp.W = (yyvsp[-2].fval); sp.Ixx = (yyvsp[-1].fval);  
          geoSource->addMat( (yyvsp[-13].ival)-1, sp );
        }
#line 9038 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 675:
#line 2419 "p.y" /* yacc.c:1646  */
    { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[-14].fval);  sp.E = (yyvsp[-13].fval);  sp.nu  = (yyvsp[-12].fval);  sp.rho = (yyvsp[-11].fval);
          sp.c = (yyvsp[-10].fval);  sp.k = (yyvsp[-9].fval);  sp.eh  = (yyvsp[-8].fval);  sp.P   = (yyvsp[-7].fval);  sp.Ta  = (yyvsp[-6].fval);
          sp.Q = (yyvsp[-5].fval); sp.W = (yyvsp[-4].fval); sp.Ixx = (yyvsp[-3].fval);
          if((yyvsp[-2].ival) == 0) sp.betaDamp = (yyvsp[-1].fval); else return -1;
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 9050 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 676:
#line 2428 "p.y" /* yacc.c:1646  */
    { StructProp sp; // this is used for the reduced mesh file output in Rom.d/MeshOutput.C
                         // all properties relevant to structural nonlinear dynamics should be included
          sp.A = (yyvsp[-28].fval);  sp.E = (yyvsp[-27].fval);  sp.nu  = (yyvsp[-26].fval);  sp.rho = (yyvsp[-25].fval);
          sp.c = (yyvsp[-24].fval);  sp.k = (yyvsp[-23].fval);  sp.eh  = (yyvsp[-22].fval);  sp.P   = (yyvsp[-21].fval);  sp.Ta  = (yyvsp[-20].fval);
          sp.Q = (yyvsp[-19].fval); sp.W = (yyvsp[-18].fval); sp.Ixx = (yyvsp[-17].fval); sp.Iyy = (yyvsp[-16].fval); sp.Izz = (yyvsp[-15].fval);
          sp.ymin = (yyvsp[-14].fval); sp.ymax = (yyvsp[-13].fval); sp.zmin = (yyvsp[-12].fval); sp.zmax = (yyvsp[-11].fval);
          sp.betaDamp = (yyvsp[-10].fval); sp.alphaDamp = (yyvsp[-9].fval); 
          sp.lagrangeMult = bool((yyvsp[-8].ival)); sp.penalty = (yyvsp[-7].fval); sp.initialPenalty = (yyvsp[-6].fval);
          sp.funtype = (yyvsp[-5].ival); sp.type = StructProp::PropType((yyvsp[-4].ival)); sp.k1 = (yyvsp[-3].fval); sp.k2 = (yyvsp[-2].fval); sp.k3 = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-29].ival)-1, sp );
        }
#line 9066 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 677:
#line 2440 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-9].fval),0.0);
          sp.fp.PMLtype = int((yyvsp[-8].fval));
          sp.fp.gamma = (yyvsp[-7].fval);
          sp.fp.Rx = (yyvsp[-6].fval);
          sp.fp.Sx = (yyvsp[-5].fval);
          sp.fp.Ry = (yyvsp[-4].fval);
          sp.fp.Sy = (yyvsp[-3].fval);
          sp.fp.Rz = (yyvsp[-2].fval);
          sp.fp.Sz = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
#line 9086 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 678:
#line 2456 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-10].fval),0.0);
          sp.rho = (yyvsp[-9].fval);
          sp.fp.PMLtype = int((yyvsp[-8].fval));
          sp.fp.gamma = (yyvsp[-7].fval);
          sp.fp.Rx = (yyvsp[-6].fval);
          sp.fp.Sx = (yyvsp[-5].fval);
          sp.fp.Ry = (yyvsp[-4].fval);
          sp.fp.Sy = (yyvsp[-3].fval);
          sp.fp.Rz = (yyvsp[-2].fval);
          sp.fp.Sz = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-12].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
#line 9107 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 679:
#line 2473 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-11].fval),(yyvsp[-10].fval));
          sp.rho = (yyvsp[-9].fval);
          sp.fp.PMLtype = int((yyvsp[-8].fval));
          sp.fp.gamma = (yyvsp[-7].fval);
          sp.fp.Rx = (yyvsp[-6].fval);
          sp.fp.Sx = (yyvsp[-5].fval);
          sp.fp.Ry = (yyvsp[-4].fval);
          sp.fp.Sy = (yyvsp[-3].fval);
          sp.fp.Rz = (yyvsp[-2].fval);
          sp.fp.Sz = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-13].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        }
#line 9128 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 680:
#line 2490 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-2].fval),0.0);
          sp.rho = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
          domain->solInfo().acoustic = true;
        }
#line 9140 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 681:
#line 2498 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[-3].fval),(yyvsp[-2].fval));
          sp.rho = (yyvsp[-1].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
          domain->solInfo().acoustic = true;
        }
#line 9152 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 682:
#line 2506 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.F_op = (yyvsp[-13].ival);
          sp.E = (yyvsp[-12].fval);
          sp.rho = (yyvsp[-11].fval);
	  sp.A = (yyvsp[-10].fval);
          sp.F_Uc = (yyvsp[-9].fval);
          sp.F_Uf = (yyvsp[-8].fval);
          sp.lambda = (yyvsp[-7].fval);
          sp.F_h = (yyvsp[-6].fval);
          sp.F_d = (yyvsp[-5].fval);
          sp.F_dlambda = (yyvsp[-4].fval);
          sp.F_np = (yyvsp[-3].ival);
          sp.F_Nf = (yyvsp[-2].ival);
	  sp.Seed = (yyvsp[-1].ival);
          sp.type = StructProp::Fabric;
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 9174 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 683:
#line 2524 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
          sp.A = (yyvsp[-9].fval); sp.rho = (yyvsp[-8].fval); sp.Q = (yyvsp[-7].fval); sp.c = (yyvsp[-6].fval); 
          sp.sigma = (yyvsp[-5].fval); sp.k = (yyvsp[-4].fval); sp.eh = (yyvsp[-3].fval); sp.P = (yyvsp[-2].fval); sp.Ta = (yyvsp[-1].fval);
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9185 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 684:
#line 2531 "p.y" /* yacc.c:1646  */
    { // rigid element or joint with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-2].ival)-1, sp );
        }
#line 9196 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 685:
#line 2538 "p.y" /* yacc.c:1646  */
    { // rigid element or joint
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-1].copt).penalty;
          sp.constraint_hess = (yyvsp[-1].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9211 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 686:
#line 2549 "p.y" /* yacc.c:1646  */
    { //           m     Ixx   Iyy   Izz   Ixy   Iyz   Ixz   cx    cy    cz
          // discrete mass with offset
          StructProp sp;
          sp.rho = (yyvsp[-10].fval);
          sp.Ixx = (yyvsp[-9].fval);
          sp.Iyy = (yyvsp[-8].fval);
          sp.Izz = (yyvsp[-7].fval);
          sp.Ixy = (yyvsp[-6].fval);
          sp.Iyz = (yyvsp[-5].fval);
          sp.Ixz = (yyvsp[-4].fval);
          sp.cx  = (yyvsp[-3].fval);
          sp.cy  = (yyvsp[-2].fval);
          sp.cz  = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-12].ival)-1, sp );
        }
#line 9231 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 687:
#line 2565 "p.y" /* yacc.c:1646  */
    { // rigid 8-node brick element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9242 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 688:
#line 2572 "p.y" /* yacc.c:1646  */
    { // rigid 8-node brick element with mass
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-3].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-3].copt).penalty;
          sp.constraint_hess = (yyvsp[-3].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-3].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9257 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 689:
#line 2583 "p.y" /* yacc.c:1646  */
    { // rigid beam or shell element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[-2].fval);
          sp.A = sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9269 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 690:
#line 2591 "p.y" /* yacc.c:1646  */
    { // rigid beam or shell element with mass
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-4].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-4].copt).penalty;
          sp.constraint_hess = (yyvsp[-4].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-4].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = (yyvsp[-2].fval);
          sp.A = sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 691:
#line 2603 "p.y" /* yacc.c:1646  */
    { // constraint function element XXX deprecated format
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[-8].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[-7].fval);
          sp.amplitude = (yyvsp[-6].fval);
          sp.omega = (yyvsp[-5].fval);
          sp.phase = (yyvsp[-4].fval);
          sp.B = (yyvsp[-3].fval);
          sp.C = (yyvsp[-2].fval);
          sp.relop = (yyvsp[-1].ival);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9304 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 692:
#line 2618 "p.y" /* yacc.c:1646  */
    { // constraint function element with default constraint options
          StructProp sp;
          sp.amplitude = 0;
          sp.omega = 0;
          sp.phase = 0;
          sp.B = 0;
          sp.C = 0;
          sp.relop = (yyvsp[-1].ival);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9321 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 693:
#line 2631 "p.y" /* yacc.c:1646  */
    { // constraint function element
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-3].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-3].copt).penalty;
          sp.constraint_hess = (yyvsp[-3].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-3].copt).constraint_hess_eps;
          sp.amplitude = 0;
          sp.omega = 0;
          sp.phase = 0;
          sp.B = 0;
          sp.C = 0;
          sp.relop = (yyvsp[-1].ival);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9342 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 694:
#line 2648 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-5].ival);
          sp.amplitude = (yyvsp[-4].fval);
          sp.offset = (yyvsp[-3].fval);
          sp.c1 = (yyvsp[-2].fval);
          sp.c2 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-7].ival)-1, sp );
        }
#line 9358 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 695:
#line 2660 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-6].ival);
          sp.amplitude = (yyvsp[-5].fval);
          sp.offset = (yyvsp[-4].fval);
          sp.c1 = (yyvsp[-3].fval);
          sp.c2 = (yyvsp[-2].fval);
          sp.c3 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-8].ival)-1, sp );
        }
#line 9375 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 696:
#line 2673 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.c3 = (yyvsp[-2].fval);
          sp.c4 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9393 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 697:
#line 2687 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-6].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-6].copt).penalty;
          sp.constraint_hess = (yyvsp[-6].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-6].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-5].ival);
          sp.amplitude = (yyvsp[-4].fval);
          sp.offset = (yyvsp[-3].fval);
          sp.c1 = (yyvsp[-2].fval);
          sp.c2 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-8].ival)-1, sp );
        }
#line 9413 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 698:
#line 2703 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-7].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-7].copt).penalty;
          sp.constraint_hess = (yyvsp[-7].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-7].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-6].ival);
          sp.amplitude = (yyvsp[-5].fval);
          sp.offset = (yyvsp[-4].fval);
          sp.c1 = (yyvsp[-3].fval);
          sp.c2 = (yyvsp[-2].fval);
          sp.c3 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9434 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 699:
#line 2720 "p.y" /* yacc.c:1646  */
    { // joint-with-driver, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-8].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-8].copt).penalty;
          sp.constraint_hess = (yyvsp[-8].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-8].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.c3 = (yyvsp[-2].fval);
          sp.c4 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9456 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 700:
#line 2738 "p.y" /* yacc.c:1646  */
    { // actuated joint, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9473 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 701:
#line 2751 "p.y" /* yacc.c:1646  */
    { // actuated joint, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-8].ival);
          sp.amplitude = (yyvsp[-7].fval);
          sp.offset = (yyvsp[-6].fval);
          sp.c1 = (yyvsp[-5].fval);
          sp.c2 = (yyvsp[-4].fval);
          sp.c3 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9491 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 702:
#line 2765 "p.y" /* yacc.c:1646  */
    { // actuated joints, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[-9].ival);
          sp.amplitude = (yyvsp[-8].fval);
          sp.offset = (yyvsp[-7].fval);
          sp.c1 = (yyvsp[-6].fval);
          sp.c2 = (yyvsp[-5].fval);
          sp.c3 = (yyvsp[-4].fval);
          sp.c4 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9510 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 703:
#line 2780 "p.y" /* yacc.c:1646  */
    { // actuated joint, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-8].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-8].copt).penalty;
          sp.constraint_hess = (yyvsp[-8].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-8].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-7].ival);
          sp.amplitude = (yyvsp[-6].fval);
          sp.offset = (yyvsp[-5].fval);
          sp.c1 = (yyvsp[-4].fval);
          sp.c2 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9531 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 704:
#line 2797 "p.y" /* yacc.c:1646  */
    { // actuated joint, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-9].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-9].copt).penalty;
          sp.constraint_hess = (yyvsp[-9].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-9].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-8].ival);
          sp.amplitude = (yyvsp[-7].fval);
          sp.offset = (yyvsp[-6].fval);
          sp.c1 = (yyvsp[-5].fval);
          sp.c2 = (yyvsp[-4].fval);
          sp.c3 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9553 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 705:
#line 2815 "p.y" /* yacc.c:1646  */
    { // actuated joints, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-10].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-10].copt).penalty;
          sp.constraint_hess = (yyvsp[-10].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-10].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[-9].ival);
          sp.amplitude = (yyvsp[-8].fval);
          sp.offset = (yyvsp[-7].fval);
          sp.c1 = (yyvsp[-6].fval);
          sp.c2 = (yyvsp[-5].fval);
          sp.c3 = (yyvsp[-4].fval);
          sp.c4 = (yyvsp[-3].fval);
          sp.k1 = (yyvsp[-1].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-12].ival)-1, sp );
        }
#line 9576 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 706:
#line 2834 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9588 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 707:
#line 2842 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-3].fval);
          sp.freeplay[0] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9601 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 708:
#line 2851 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-3].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-3].copt).penalty;
          sp.constraint_hess = (yyvsp[-3].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-3].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9617 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 709:
#line 2863 "p.y" /* yacc.c:1646  */
    { // RevoluteJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-5].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-5].copt).penalty;
          sp.constraint_hess = (yyvsp[-5].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-5].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-3].fval);
          sp.freeplay[0] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-7].ival)-1, sp );
        }
#line 9634 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 710:
#line 2876 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-2].fval);
          sp.k2 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9647 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 711:
#line 2885 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-5].fval);
          sp.k2 = (yyvsp[-4].fval);
          sp.freeplay[0] = (yyvsp[-2].freeplayProps);
          sp.freeplay[1] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-8].ival)-1, sp );
        }
#line 9662 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 712:
#line 2896 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-4].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-4].copt).penalty;
          sp.constraint_hess = (yyvsp[-4].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-4].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-2].fval);
          sp.k2 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9679 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 713:
#line 2909 "p.y" /* yacc.c:1646  */
    { // UniversalJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-7].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-7].copt).penalty;
          sp.constraint_hess = (yyvsp[-7].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-7].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-5].fval);
          sp.k2 = (yyvsp[-4].fval);
          sp.freeplay[0] = (yyvsp[-2].freeplayProps);
          sp.freeplay[1] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-9].ival)-1, sp );
        }
#line 9698 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 714:
#line 2924 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-3].fval);
          sp.k2 = (yyvsp[-2].fval);
          sp.k3 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-6].ival)-1, sp );
        }
#line 9712 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 715:
#line 2934 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringComboWithFreeplay with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[-7].fval);
          sp.k2 = (yyvsp[-6].fval);
          sp.k3 = (yyvsp[-5].fval);
          sp.freeplay[0] = (yyvsp[-3].freeplayProps);
          sp.freeplay[1] = (yyvsp[-2].freeplayProps);
          sp.freeplay[2] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-10].ival)-1, sp );
        }
#line 9729 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 716:
#line 2947 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-5].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-5].copt).penalty;
          sp.constraint_hess = (yyvsp[-5].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-5].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-3].fval);
          sp.k2 = (yyvsp[-2].fval);
          sp.k3 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-7].ival)-1, sp );
        }
#line 9747 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 717:
#line 2961 "p.y" /* yacc.c:1646  */
    { // SphericalJointSpringComboWithFreeplay
          StructProp sp;
          sp.lagrangeMult = (yyvsp[-9].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[-9].copt).penalty;
          sp.constraint_hess = (yyvsp[-9].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[-9].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[-7].fval);
          sp.k2 = (yyvsp[-6].fval);
          sp.k3 = (yyvsp[-5].fval);
          sp.freeplay[0] = (yyvsp[-3].freeplayProps);
          sp.freeplay[1] = (yyvsp[-2].freeplayProps);
          sp.freeplay[2] = (yyvsp[-1].freeplayProps);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9768 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 718:
#line 2978 "p.y" /* yacc.c:1646  */
    { // TorsionalSpring or TranslationalSpring (types 200,201,202)
          StructProp sp;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9779 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 719:
#line 2985 "p.y" /* yacc.c:1646  */
    { // 1-sided TorsionalSpring or TranslationalSpring with freeplay (types 203,204,205)
          StructProp sp;
          sp.k1 = (yyvsp[-3].fval);
          sp.rho = 0;
          sp.freeplay[0].ul = (yyvsp[-1].fval);
          sp.freeplay[0].dz = 0.0;
          sp.freeplay[0].uz = 1.0;
          geoSource->addMat( (yyvsp[-5].ival)-1, sp );
        }
#line 9793 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 720:
#line 2995 "p.y" /* yacc.c:1646  */
    { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = (yyvsp[-9].fval); sp.dE = (yyvsp[-8].fval)/(2*M_PI); sp.eta_E = (yyvsp[-7].fval); sp.deta_E = (yyvsp[-6].fval)/(2*M_PI);
          sp.mu0 = (yyvsp[-5].fval); sp.dmu = (yyvsp[-4].fval)/(2*M_PI); sp.eta_mu = (yyvsp[-3].fval); sp.deta_mu = (yyvsp[-2].fval)/(2*M_PI);
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9805 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 721:
#line 3003 "p.y" /* yacc.c:1646  */
    { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9815 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 722:
#line 3019 "p.y" /* yacc.c:1646  */
    { // 5-parameter freeplay model
          (yyval.freeplayProps).ll = (yyvsp[-4].fval);
          (yyval.freeplayProps).ul = (yyvsp[-3].fval);
          (yyval.freeplayProps).lz = (yyvsp[-2].fval);
          (yyval.freeplayProps).dz = (yyvsp[-1].fval);
          (yyval.freeplayProps).uz = (yyvsp[0].fval);
        }
#line 9827 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 725:
#line 3032 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-1].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9837 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 726:
#line 3038 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-2].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-2].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9847 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 727:
#line 3044 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-3].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-3].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[-1].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9858 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 728:
#line 3051 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-4].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-4].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[-2].fval));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9870 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 729:
#line 3059 "p.y" /* yacc.c:1646  */
    { if((yyval.SurfObj)->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[(yyvsp[-1].nl).num];
            for(int i=0; i<(yyvsp[-1].nl).num; ++i) nodes[(yyvsp[-1].nl).num-1-i] = (yyvsp[-1].nl).nd[i];
            (yyval.SurfObj)->AddFaceElement((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, nodes);
            delete [] nodes;
          }
          else (yyval.SurfObj)->AddFaceElement((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
        }
#line 9883 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 730:
#line 3070 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9889 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 731:
#line 3072 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9895 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 732:
#line 3074 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9903 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 733:
#line 3078 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9909 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 734:
#line 3080 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9915 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 735:
#line 3082 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9921 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 736:
#line 3084 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9927 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 737:
#line 3086 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-1].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9935 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 738:
#line 3090 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9943 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 739:
#line 3096 "p.y" /* yacc.c:1646  */
    { domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().isCoupled = true; }
#line 9949 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 740:
#line 3098 "p.y" /* yacc.c:1646  */
    { domain->addWetInterface((yyvsp[-1].ival), (yyvsp[-1].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; }
#line 9957 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 741:
#line 3106 "p.y" /* yacc.c:1646  */
    { }
#line 9963 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 742:
#line 3108 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9973 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 743:
#line 3114 "p.y" /* yacc.c:1646  */
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-1].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9984 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 744:
#line 3121 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 9995 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 745:
#line 3128 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-3].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10006 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 746:
#line 3135 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-7].ival), (yyvsp[-6].ival), (yyvsp[-4].fval), (yyvsp[-3].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-5].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-2].ival), (yyvsp[-1].fval));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10018 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 747:
#line 3143 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10029 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 748:
#line 3150 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10041 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 749:
#line 3158 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-3].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 750:
#line 3166 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-3].fval), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-4].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10065 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 751:
#line 3174 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-5].fval), (yyvsp[-4].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-6].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-3].ival), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 752:
#line 3185 "p.y" /* yacc.c:1646  */
    { }
#line 10084 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 753:
#line 3187 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[-2].ival) == (yyvsp[-1].ival)) domain->solInfo().isMatching = true;
        }
#line 10094 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 754:
#line 3193 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[-4].ival) == (yyvsp[-3].ival)) domain->solInfo().isMatching = true;
        }
#line 10104 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 755:
#line 3201 "p.y" /* yacc.c:1646  */
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
#line 10111 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 757:
#line 3207 "p.y" /* yacc.c:1646  */
    { domain->addWetElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
#line 10119 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 758:
#line 3213 "p.y" /* yacc.c:1646  */
    { }
#line 10125 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 759:
#line 3215 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[-2].ival) == (yyvsp[-1].ival)) domain->solInfo().isMatching = true;
        }
#line 10135 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 760:
#line 3221 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[-4].ival) == (yyvsp[-3].ival)) domain->solInfo().isMatching = true;
        }
#line 10145 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 761:
#line 3231 "p.y" /* yacc.c:1646  */
    { }
#line 10151 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 762:
#line 3233 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contactsurface_mode = (yyvsp[-1].ival); }
#line 10157 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 763:
#line 3235 "p.y" /* yacc.c:1646  */
    { domain->AddMortarCond((yyvsp[-1].MortarCondObj)); }
#line 10163 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 764:
#line 3237 "p.y" /* yacc.c:1646  */
    { (yyvsp[-2].MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt)); domain->AddMortarCond((yyvsp[-2].MortarCondObj)); }
#line 10169 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 765:
#line 3239 "p.y" /* yacc.c:1646  */
    { (yyvsp[-4].MortarCondObj)->SetConstraintOptions((yyvsp[-3].copt)); (yyvsp[-4].MortarCondObj)->SetCtcMode((yyvsp[-1].ival)); domain->AddMortarCond((yyvsp[-4].MortarCondObj)); }
#line 10175 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 766:
#line 3241 "p.y" /* yacc.c:1646  */
    { (yyvsp[-3].MortarCondObj)->SetCtcMode((yyvsp[-1].ival)); domain->AddMortarCond((yyvsp[-3].MortarCondObj)); }
#line 10181 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 767:
#line 3247 "p.y" /* yacc.c:1646  */
    { domain->solInfo().trivial_detection = true;
          (yyval.MortarCondObj) = new MortarHandler(0, 0);
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD);
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode); }
#line 10191 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 768:
#line 3253 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-1].ival), (yyvsp[0].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD);
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10202 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 769:
#line 3260 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[0].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10213 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 770:
#line 3267 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10224 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 771:
#line 3274 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10235 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 772:
#line 3281 "p.y" /* yacc.c:1646  */
    { /* frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-3].fval), (yyvsp[-2].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-1].ival), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-4].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10247 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 773:
#line 3289 "p.y" /* yacc.c:1646  */
    { /* constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-7].ival), (yyvsp[-6].ival), (yyvsp[-4].fval), (yyvsp[-3].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-2].ival), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-5].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 774:
#line 3298 "p.y" /* yacc.c:1646  */
    { /* velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-6].fval), (yyvsp[-5].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-4].ival), (yyvsp[-3].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-7].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10273 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 775:
#line 3307 "p.y" /* yacc.c:1646  */
    { /* pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-7].fval), (yyvsp[-6].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[-5].ival), (yyvsp[-4].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-8].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10286 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 777:
#line 3319 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dist_acme = (yyvsp[-1].ival); }
#line 10292 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 778:
#line 3321 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_secondary = true; }
#line 10298 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 779:
#line 3323 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_ghosting = true; }
#line 10304 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 780:
#line 3325 "p.y" /* yacc.c:1646  */
    { domain->solInfo().shell_simple_lofting = true; }
#line 10310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 781:
#line 3327 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_multiple_interactions = true; }
#line 10316 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 782:
#line 3329 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sharp_non_sharp_angle = (yyvsp[-1].fval); }
#line 10322 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 783:
#line 3331 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normal_smoothing = false; }
#line 10328 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 784:
#line 3333 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normal_smoothing_distance = (yyvsp[-1].fval); }
#line 10334 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 785:
#line 3335 "p.y" /* yacc.c:1646  */
    { domain->solInfo().resolution_method = (yyvsp[-1].ival); }
#line 10340 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 786:
#line 3337 "p.y" /* yacc.c:1646  */
    { domain->solInfo().old_dynamic_search = true; }
#line 10346 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 787:
#line 3339 "p.y" /* yacc.c:1646  */
    { domain->solInfo().partition_gap = true; }
#line 10352 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 788:
#line 3341 "p.y" /* yacc.c:1646  */
    { domain->solInfo().default_penalty = true; }
#line 10358 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 789:
#line 3343 "p.y" /* yacc.c:1646  */
    { domain->solInfo().global_search_cull = true; }
#line 10364 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 790:
#line 3345 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_warped_volume = true; }
#line 10370 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 791:
#line 3347 "p.y" /* yacc.c:1646  */
    { domain->solInfo().auto_tol = true; }
#line 10376 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 792:
#line 3349 "p.y" /* yacc.c:1646  */
    { domain->solInfo().auto_tol = true;
          domain->solInfo().agressive_tolerances = true; }
#line 10383 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 793:
#line 3352 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skip_physical_faces = true; }
#line 10389 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 794:
#line 3354 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ffi_debug = bool((yyvsp[-1].ival)); }
#line 10395 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 795:
#line 3356 "p.y" /* yacc.c:1646  */
    { domain->solInfo().mortar_scaling = (yyvsp[-1].fval); }
#line 10401 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 796:
#line 3358 "p.y" /* yacc.c:1646  */
    { domain->solInfo().mortar_integration_rule = (yyvsp[-1].ival); }
#line 10407 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 797:
#line 3360 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_clr = (yyvsp[-1].fval); }
#line 10413 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 798:
#line 3362 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_cqr = (yyvsp[-1].fval); }
#line 10419 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 799:
#line 3364 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_betab = (yyvsp[-1].fval); }
#line 10425 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 800:
#line 3366 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_alpha = (yyvsp[-1].fval); }
#line 10431 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 801:
#line 3368 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_betam = (yyvsp[-1].fval); }
#line 10437 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 802:
#line 3370 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nlmembrane_pressure_type = (yyvsp[-1].ival); }
#line 10443 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 803:
#line 3373 "p.y" /* yacc.c:1646  */
    { geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10449 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 804:
#line 3375 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scalePosCoords = true;
          domain->solInfo().xScaleFactor = (yyvsp[-4].fval);
          domain->solInfo().yScaleFactor = (yyvsp[-3].fval);
          domain->solInfo().zScaleFactor = (yyvsp[-2].fval);
          geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10459 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 805:
#line 3381 "p.y" /* yacc.c:1646  */
    { geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10465 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 806:
#line 3385 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-4].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-3].fval); (yyval.nval).xyz[1] = (yyvsp[-2].fval);  (yyval.nval).xyz[2] = (yyvsp[-1].fval);  (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10471 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 807:
#line 3387 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-3].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-2].fval); (yyval.nval).xyz[1] = (yyvsp[-1].fval);  (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10477 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 808:
#line 3389 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-2].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-1].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10483 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 809:
#line 3391 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-6].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-5].fval); (yyval.nval).xyz[1] = (yyvsp[-4].fval);  (yyval.nval).xyz[2] = (yyvsp[-3].fval);  (yyval.nval).cp = (yyvsp[-2].ival); (yyval.nval).cd = (yyvsp[-1].ival);
          if((yyvsp[-2].ival) != 0) domain->solInfo().basicPosCoords = false;
          if((yyvsp[-1].ival) != 0) domain->solInfo().basicDofCoords = false; }
#line 10491 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 810:
#line 3395 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-5].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-4].fval); (yyval.nval).xyz[1] = (yyvsp[-3].fval);  (yyval.nval).xyz[2] = (yyvsp[-2].fval);  (yyval.nval).cp = (yyvsp[-1].ival); (yyval.nval).cd = (yyvsp[-1].ival);
          if((yyvsp[-1].ival) != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } }
#line 10498 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 811:
#line 3400 "p.y" /* yacc.c:1646  */
    { /* Define each Element */
          geoSource->addElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);}
#line 10505 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 812:
#line 3405 "p.y" /* yacc.c:1646  */
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[0].ival)-1;}
#line 10511 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 813:
#line 3407 "p.y" /* yacc.c:1646  */
    { if((yyval.nl).num == 500) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[0].ival)-1; (yyval.nl).num++;}
#line 10518 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 814:
#line 3412 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-3].ival)-1; (yyval.bcval).dofnum = (yyvsp[-2].ival)-1; (yyval.bcval).val = (yyvsp[-1].fval); (yyval.bcval).mtype = BCond::Axial; }
#line 10524 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 815:
#line 3414 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = (yyvsp[-1].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = BCond::Axial; }
#line 10530 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 816:
#line 3416 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-4].ival)-1; (yyval.bcval).dofnum = (yyvsp[-3].ival)-1; (yyval.bcval).val = (yyvsp[-2].fval); (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[-1].ival); }
#line 10536 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 817:
#line 3418 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-3].ival)-1; (yyval.bcval).dofnum = (yyvsp[-2].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[-1].ival); }
#line 10542 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 818:
#line 3422 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[-1].fval); }
#line 10548 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 819:
#line 3426 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 10554 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 820:
#line 3430 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum = (yyvsp[-4].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[-3].ival)-1; (yyval.cxbcval).reval = (yyvsp[-2].fval); (yyval.cxbcval).imval = (yyvsp[-1].fval);  }
#line 10560 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 821:
#line 3432 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum = (yyvsp[-3].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[-2].ival)-1; (yyval.cxbcval).reval = (yyvsp[-1].fval); (yyval.cxbcval).imval = 0.0; }
#line 10566 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 823:
#line 3437 "p.y" /* yacc.c:1646  */
    { geoSource->setCSFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d); }
#line 10572 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 825:
#line 3442 "p.y" /* yacc.c:1646  */
    { geoSource->setFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d);  }
#line 10578 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 826:
#line 3446 "p.y" /* yacc.c:1646  */
    { (yyval.frame).num = (yyvsp[-10].ival)-1; 
          (yyval.frame).d[0] = (yyvsp[-9].fval); (yyval.frame).d[1] = (yyvsp[-8].fval); (yyval.frame).d[2] = (yyvsp[-7].fval);
          (yyval.frame).d[3] = (yyvsp[-6].fval); (yyval.frame).d[4] = (yyvsp[-5].fval); (yyval.frame).d[5] = (yyvsp[-4].fval);
          (yyval.frame).d[6] = (yyvsp[-3].fval); (yyval.frame).d[7] = (yyvsp[-2].fval); (yyval.frame).d[8] = (yyvsp[-1].fval); }
#line 10587 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 827:
#line 3451 "p.y" /* yacc.c:1646  */
    { (yyval.frame).num = (yyvsp[-3].ival)-1;
          geoSource->makeEframe((yyvsp[-3].ival)-1, (yyvsp[-1].ival), (yyval.frame).d); }
#line 10594 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 829:
#line 3457 "p.y" /* yacc.c:1646  */
    { geoSource->setNodalFrame((yyvsp[0].nframe).id,(yyvsp[0].nframe).o,(yyvsp[0].nframe).d,(yyvsp[0].nframe).type); }
#line 10600 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 830:
#line 3461 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-10].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;  (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[-9].fval); (yyval.nframe).d[1] = (yyvsp[-8].fval); (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval); (yyval.nframe).d[4] = (yyvsp[-5].fval); (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10611 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 831:
#line 3468 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-13].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = (yyvsp[-12].fval);  (yyval.nframe).o[1] = (yyvsp[-11].fval);  (yyval.nframe).o[2] = (yyvsp[-10].fval);
          (yyval.nframe).d[0] = (yyvsp[-9].fval);  (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval);  (yyval.nframe).d[4] = (yyvsp[-5].fval);  (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10622 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 832:
#line 3475 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-11].ival);
          (yyval.nframe).type = (yyvsp[-10].ival);
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;   (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[-9].fval); (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval); (yyval.nframe).d[4] = (yyvsp[-5].fval);  (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10633 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 833:
#line 3482 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-14].ival);
          (yyval.nframe).type = (yyvsp[-13].ival);
          (yyval.nframe).o[0] = (yyvsp[-12].fval);  (yyval.nframe).o[1] = (yyvsp[-11].fval);  (yyval.nframe).o[2] = (yyvsp[-10].fval);
          (yyval.nframe).d[0] = (yyvsp[-9].fval);  (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval);  (yyval.nframe).d[4] = (yyvsp[-5].fval); (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10644 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 835:
#line 3492 "p.y" /* yacc.c:1646  */
    { OffsetData od;
	  od.first = (yyvsp[-5].ival)-1; od.last = (yyvsp[-4].ival)-1;
	  od.o[0] = (yyvsp[-3].fval); od.o[1] = (yyvsp[-2].fval); od.o[2] = (yyvsp[-1].fval); 
	  geoSource->addOffset(od); }
#line 10653 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 836:
#line 3499 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 10659 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 837:
#line 3501 "p.y" /* yacc.c:1646  */
    { geoSource->setLocalIndex((yyvsp[-1].ival)-1); }
#line 10665 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 838:
#line 3504 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-3].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10672 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 839:
#line 3507 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-4].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10680 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 840:
#line 3512 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-2].ival)-1,(yyvsp[-1].ival)-1); }
#line 10686 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 841:
#line 3514 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1); 
          geoSource->setElementLumpingWeight((yyvsp[-4].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10694 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 842:
#line 3518 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-5].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; 
          domain->solInfo().reduceFollower = true; }
#line 10703 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 843:
#line 3523 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-5].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10712 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 844:
#line 3529 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].ival)-1); }
#line 10718 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 845:
#line 3531 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-6].ival)-1,(yyvsp[-1].fval)); 
          domain->solInfo().elemLumpPodRom = true; }
#line 10726 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 846:
#line 3535 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-7].ival)-1,(yyvsp[-2].fval));   
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10735 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 847:
#line 3541 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,-2,(yyvsp[-1].fval)); }
#line 10741 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 848:
#line 3543 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,-2,(yyvsp[-3].fval));
          geoSource->setElementLumpingWeight((yyvsp[-7].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10749 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 849:
#line 3547 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-8].ival)-1,(yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,-2,(yyvsp[-4].fval)); 
          geoSource->setElementLumpingWeight((yyvsp[-8].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10758 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 850:
#line 3553 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        }
#line 10767 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 851:
#line 3558 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
            geoSource->setAttrib(i-1,(yyvsp[-1].ival)-1);
        }
#line 10776 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 852:
#line 3563 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-5].ival); i<(yyvsp[-4].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1);
        }
#line 10785 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 853:
#line 3568 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-6].ival); i<(yyvsp[-5].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, -2, (yyvsp[-1].fval));
        }
#line 10794 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 854:
#line 3575 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex(0); }
#line 10801 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 855:
#line 3578 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex((yyvsp[-1].ival)-1); }
#line 10808 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 856:
#line 3581 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-2].ival) - 1, (yyvsp[-1].fval)); }
#line 10814 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 857:
#line 3583 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = true;}
#line 10820 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 858:
#line 3585 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = true;}
#line 10826 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 862:
#line 3594 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInLocalBasesAuxi[std::make_pair((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)] = std::string((yyvsp[-1].strval)); }
#line 10832 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 863:
#line 3598 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInLocalBasesCent.push_back(std::string((yyvsp[-1].strval))); }
#line 10838 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 864:
#line 3602 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ReducedStiffness = true;}
#line 10844 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 865:
#line 3604 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackStiffVec((yyvsp[-1].fval));}
#line 10850 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 867:
#line 3609 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodSize = (yyvsp[-2].ival);
          domain->solInfo().maxDeimBasisSize = (yyvsp[-1].ival);}
#line 10857 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 868:
#line 3612 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackUDEIMVec((yyvsp[-1].fval));}
#line 10863 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 870:
#line 3617 "p.y" /* yacc.c:1646  */
    { domain->solInfo().DEIMPodRom = true;
          geoSource->setSampleNodesAndSlots((yyvsp[-2].ival)-1,(yyvsp[-1].ival));}
#line 10870 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 871:
#line 3620 "p.y" /* yacc.c:1646  */
    { geoSource->setSampleElemsAndDOFs((yyvsp[-2].ival)-1,(yyvsp[-1].ival));
          domain->solInfo().UDEIMPodRom = true;}
#line 10877 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 872:
#line 3626 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 10883 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 873:
#line 3628 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 10889 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 874:
#line 3630 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-2].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          geoSource->setElementPressure(pbc); }
#line 10897 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 875:
#line 3634 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-3].ival); i < ((yyvsp[-2].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            geoSource->setElementPressure(pbc);
          } }
#line 10907 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 876:
#line 3640 "p.y" /* yacc.c:1646  */
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[-2].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
#line 10916 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 877:
#line 3645 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-4].ival); i <= (yyvsp[-2].ival); ++i) {
            PressureBCond *pbc = new PressureBCond[1];
            pbc[0].setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            geoSource->addSurfacePressure(1, pbc);
            if(geoSource->getNumSurfacePressure() > 1) delete [] pbc;
          } }
#line 10927 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 878:
#line 3652 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-3].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          geoSource->setElementPressure(pbc); }
#line 10935 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 879:
#line 3656 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-4].ival); i < ((yyvsp[-3].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
            geoSource->setElementPressure(pbc);
          } }
#line 10945 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 880:
#line 3662 "p.y" /* yacc.c:1646  */
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[-3].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
#line 10954 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 881:
#line 3668 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-4].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          pbc.face = (yyvsp[-2].ival)-1;
          geoSource->setElementPressure(pbc); }
#line 10963 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 882:
#line 3673 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-5].ival); i < ((yyvsp[-4].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            pbc.face = (yyvsp[-2].ival)-1;
            geoSource->setElementPressure(pbc);
          } }
#line 10974 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 883:
#line 3680 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-5].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          pbc.face = (yyvsp[-3].ival)-1;
          geoSource->setElementPressure(pbc); }
#line 10983 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 884:
#line 3685 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-6].ival); i < ((yyvsp[-5].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
            pbc.face = (yyvsp[-3].ival)-1;
            geoSource->setElementPressure(pbc);
          } }
#line 10994 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 885:
#line 3694 "p.y" /* yacc.c:1646  */
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        }
#line 11003 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 886:
#line 3699 "p.y" /* yacc.c:1646  */
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[-1].ival));
          geoSource->setConsistentPFlag(false);
        }
#line 11012 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 887:
#line 3713 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassAugmentation = (yyvsp[-1].ival); }
#line 11018 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 888:
#line 3717 "p.y" /* yacc.c:1646  */
    { }
#line 11024 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 889:
#line 3719 "p.y" /* yacc.c:1646  */
    { geoSource->setElementPreLoad( (yyvsp[-2].ival)-1, (yyvsp[-1].fval) ); }
#line 11030 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 890:
#line 3721 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-4].ival); i<((yyvsp[-2].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[-1].fval) );
        }
#line 11039 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 891:
#line 3726 "p.y" /* yacc.c:1646  */
    { double load[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
          geoSource->setElementPreLoad( (yyvsp[-4].ival)-1, load ); }
#line 11046 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 892:
#line 3729 "p.y" /* yacc.c:1646  */
    { double load[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
          int i;
          for(i=(yyvsp[-6].ival); i<((yyvsp[-4].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        }
#line 11056 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 893:
#line 3737 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivity = true; }
#line 11062 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 894:
#line 3739 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivityMethod = (SolverInfo::SensitivityMethod) (yyvsp[-1].ival); }
#line 11068 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 895:
#line 3741 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readShapeSen = true;
          domain->solInfo().readInShapeSen = (yyvsp[-1].strval); }
#line 11075 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 896:
#line 3744 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivityTol = (yyvsp[-1].fval); }
#line 11081 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 897:
#line 3746 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qsMaxvelSen = (yyvsp[-1].fval); }
#line 11087 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 898:
#line 3748 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ratioSensitivityTol = (yyvsp[-1].fval); }
#line 11093 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 899:
#line 3750 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ksParameter = (yyvsp[-1].fval); }
#line 11099 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 900:
#line 3752 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ksMax = (yyvsp[-1].fval); }
#line 11105 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 901:
#line 3754 "p.y" /* yacc.c:1646  */
    { }
#line 11111 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 902:
#line 3756 "p.y" /* yacc.c:1646  */
    { }
#line 11117 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 903:
#line 3758 "p.y" /* yacc.c:1646  */
    { }
#line 11123 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 904:
#line 3762 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11129 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 905:
#line 3764 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11135 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 906:
#line 3766 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11141 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 907:
#line 3768 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11147 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 908:
#line 3770 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11153 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 909:
#line 3772 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11159 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 910:
#line 3774 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11165 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 911:
#line 3776 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11171 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 912:
#line 3778 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11177 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 913:
#line 3780 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11183 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 914:
#line 3782 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11189 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 915:
#line 3784 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11195 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 916:
#line 3788 "p.y" /* yacc.c:1646  */
    { domain->setStressNodes((yyvsp[0].ival)); }
#line 11201 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 917:
#line 3790 "p.y" /* yacc.c:1646  */
    { domain->setStressNodes((yyvsp[0].ival)); }
#line 11207 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 918:
#line 3794 "p.y" /* yacc.c:1646  */
    { domain->setThicknessGroup((yyvsp[0].ival)); }
#line 11213 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 919:
#line 3796 "p.y" /* yacc.c:1646  */
    { domain->setThicknessGroup((yyvsp[0].ival)); }
#line 11219 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 920:
#line 3800 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Static); }
#line 11225 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 921:
#line 3802 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solvercntl = (yyvsp[0].scntl); }
#line 11231 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 923:
#line 3805 "p.y" /* yacc.c:1646  */
    { // activate piecewise constant configuration dependent external forces for a linear dynamic analysis
          domain->solInfo().piecewise = true;
        }
#line 11239 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 924:
#line 3809 "p.y" /* yacc.c:1646  */
    { // activate piecewise constant configuration dependent external forces for a linear static analysis
          domain->solInfo().piecewise = true;
          domain->solInfo().piecewise_dlambda = (yyvsp[-2].fval);
          domain->solInfo().piecewise_maxLambda = (yyvsp[-1].fval);
        }
#line 11249 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 925:
#line 3815 "p.y" /* yacc.c:1646  */
    { domain->solInfo().coupled_scale = (yyvsp[-1].fval); }
#line 11255 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 926:
#line 3817 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-1].ival);
          domain->curvatureFlag = 0; }
#line 11262 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 927:
#line 3820 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-2].ival);
          domain->curvatureConst1 = (yyvsp[-1].fval);
          domain->curvatureFlag = 1; }
#line 11270 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 928:
#line 3824 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-3].ival);
          domain->curvatureConst1 = (yyvsp[-2].fval);
          domain->curvatureConst2 = (yyvsp[-1].fval);
          domain->curvatureFlag = 2; }
#line 11279 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 929:
#line 3829 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dmpc = bool((yyvsp[-1].ival)); }
#line 11285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 930:
#line 3831 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dbccheck = bool((yyvsp[-1].ival)); }
#line 11291 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 931:
#line 3835 "p.y" /* yacc.c:1646  */
    { domain->solInfo().loadcases.push_back((yyvsp[0].ival)); }
#line 11297 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 932:
#line 3837 "p.y" /* yacc.c:1646  */
    { domain->solInfo().loadcases.push_back((yyvsp[0].ival)); }
#line 11303 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 933:
#line 3841 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = 0; }
#line 11311 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 934:
#line 3845 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = (yyvsp[-1].ival); }
#line 11319 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 935:
#line 3849 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = (yyvsp[-1].ival); }
#line 11327 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 936:
#line 3853 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->pivot = true; }
#line 11336 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 937:
#line 3858 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->scaled = true; }
#line 11345 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 938:
#line 3863 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->unsymmetric = true; }
#line 11354 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 939:
#line 3868 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-1].ival); }
#line 11362 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 940:
#line 3872 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-2].ival);
      (yyval.scntl)->precond = (yyvsp[-1].ival); }
#line 11371 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 941:
#line 3877 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-3].ival);
      (yyval.scntl)->precond = (yyvsp[-2].ival);
      (yyval.scntl)->tol=(yyvsp[-1].fval); }
#line 11381 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 942:
#line 3883 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-4].ival);
      (yyval.scntl)->precond = (yyvsp[-3].ival);
      (yyval.scntl)->tol = (yyvsp[-2].fval);
      (yyval.scntl)->maxit = (yyvsp[-1].ival); }
#line 11392 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 943:
#line 3890 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-5].ival);
      (yyval.scntl)->precond = (yyvsp[-4].ival);
      (yyval.scntl)->tol = (yyvsp[-3].fval);
      (yyval.scntl)->maxit = (yyvsp[-2].ival);
      (yyval.scntl)->iterSubtype = (yyvsp[-1].ival); }
#line 11404 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 944:
#line 3898 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-6].ival);
      (yyval.scntl)->precond = (yyvsp[-5].ival);
      (yyval.scntl)->tol = (yyvsp[-4].fval);
      (yyval.scntl)->maxit = (yyvsp[-3].ival);
      (yyval.scntl)->iterSubtype = (yyvsp[-2].ival);
      (yyval.scntl)->maxvecsize = (yyvsp[-1].ival); }
#line 11417 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 945:
#line 3907 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti; }
#line 11424 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 946:
#line 3910 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.version = (FetiInfo::Version) ((yyvsp[-1].ival)-1); }
#line 11432 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 947:
#line 3914 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-2].ival); }
#line 11442 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 948:
#line 3920 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-1].ival); }
#line 11452 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 949:
#line 3926 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.version = (FetiInfo::Version) ((yyvsp[-2].ival)-1);
      (yyval.scntl)->fetiInfo.feti2version = (FetiInfo::Feti2Version) (yyvsp[-1].ival); }
#line 11461 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 950:
#line 3931 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp; }
#line 11470 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 951:
#line 3936 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::FetiLib;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp; }
#line 11479 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 952:
#line 3941 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true; }
#line 11489 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 953:
#line 3947 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners3;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11503 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 954:
#line 3957 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11520 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 955:
#line 3970 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-1].ival);
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11538 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 956:
#line 3984 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-4].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-3].fval);
      (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.spaceDimension = 2;
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11558 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 957:
#line 4000 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-5].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-4].fval);
      (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.spaceDimension = (yyvsp[-1].ival);
      (yyval.scntl)->fetiInfo.krylovtype = 1;
      (yyval.scntl)->fetiInfo.scaling = FetiInfo::tscaling;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true;
      (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
      (yyval.scntl)->fetiInfo.rbmType = FetiInfo::None;
      (yyval.scntl)->fetiInfo.nGs = 0; }
#line 11578 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 958:
#line 4016 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::BlockDiag; }
#line 11585 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 959:
#line 4019 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = &domain->solInfo().solvercntls[(yyvsp[-1].ival)]; }
#line 11591 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 960:
#line 4023 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = &domain->solInfo().solvercntls[(yyvsp[-2].ival)];
          *(yyval.scntl) = *(yyvsp[0].scntl); }
#line 11598 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 961:
#line 4027 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = (yyvsp[0].scntl); }
#line 11604 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 962:
#line 4029 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->verbose = 1; }
#line 11610 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 963:
#line 4031 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->verbose = (yyvsp[-1].ival); }
#line 11616 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 964:
#line 4033 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.printNumber = (yyvsp[-1].ival); }
#line 11622 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 965:
#line 4035 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->trbm = (yyvsp[-1].fval); }
#line 11628 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 966:
#line 4037 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_renum = (yyvsp[-1].ival); }
#line 11634 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 967:
#line 4039 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_maxsup = (yyvsp[-1].ival); }
#line 11640 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 968:
#line 4041 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_defblk = (yyvsp[-1].ival); }
#line 11646 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 969:
#line 4043 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_tau = (yyvsp[-1].fval); }
#line 11652 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 970:
#line 4045 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_maxsize = (yyvsp[-1].ival); }
#line 11658 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 971:
#line 4047 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) < 0) {
            (yyvsp[-1].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          (yyval.scntl)->spooles_maxdomainsize = (yyvsp[-1].ival); }
#line 11669 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 972:
#line 4054 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_seed = (yyvsp[-1].ival); }
#line 11675 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 973:
#line 4056 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].fval) < 0.0) || ((yyvsp[-1].fval) > 1.0)) {
            (yyvsp[-1].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          (yyval.scntl)->spooles_maxzeros = (yyvsp[-1].fval); }
#line 11686 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 974:
#line 4063 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_msglvl = (yyvsp[-1].ival); }
#line 11692 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 975:
#line 4065 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_scale = (yyvsp[-1].ival); }
#line 11698 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 976:
#line 4067 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_renum = (yyvsp[-1].ival); }
#line 11704 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 977:
#line 4069 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_icntl[(yyvsp[-2].ival)] = (yyvsp[-1].ival); }
#line 11710 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 978:
#line 4071 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_cntl[(yyvsp[-2].ival)] = (yyvsp[-1].fval); }
#line 11716 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 979:
#line 4073 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_mineq = (yyvsp[-1].ival); }
#line 11722 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 980:
#line 4075 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_stride = (yyvsp[-1].ival); }
#line 11728 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 981:
#line 4077 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->goldfarb_tol = (yyvsp[-1].fval); }
#line 11734 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 982:
#line 4079 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->tol = (yyvsp[-1].fval); }
#line 11740 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 983:
#line 4081 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->iterSubtype = (yyvsp[-1].ival); }
#line 11746 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 984:
#line 4083 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->ilu_droptol = (yyvsp[-1].fval); }
#line 11752 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 985:
#line 4085 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->maxit = (yyvsp[-1].ival); 
          (yyval.scntl)->fetiInfo.maxit = (yyvsp[-1].ival); }
#line 11759 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 986:
#line 4088 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->maxvecsize = (yyvsp[-1].ival);
          (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-1].ival); }
#line 11766 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 987:
#line 4091 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11772 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 988:
#line 4093 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11778 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 989:
#line 4095 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = FetiInfo::lumped; }
#line 11784 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 990:
#line 4097 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->precond = (yyvsp[-1].ival);
          if(isFeti((yyval.scntl)->type) && (((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 3))) {
            (yyvsp[-1].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11795 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 991:
#line 4104 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 1)) {
            (yyvsp[-1].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          (yyval.scntl)->fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[-1].ival); }
#line 11805 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 992:
#line 4110 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 1)) {
            (yyvsp[-1].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          (yyval.scntl)->fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[-1].ival); }
#line 11815 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 993:
#line 4116 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval); }
#line 11821 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 994:
#line 4118 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval); 
          (yyval.scntl)->fetiInfo.absolute_tol = (yyvsp[-1].fval); }
#line 11828 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 995:
#line 4121 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.stagnation_tol = (yyvsp[-1].fval); }
#line 11834 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 996:
#line 4123 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.stagnation_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.absolute_stagnation_tol = (yyvsp[-1].fval); }
#line 11841 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 997:
#line 4126 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_proj_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.dual_proj_tol = (yyvsp[-1].fval); }
#line 11848 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 998:
#line 4129 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_plan_maxit = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.dual_plan_maxit = (yyvsp[-1].ival); }
#line 11855 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 999:
#line 4132 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_plan_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.dual_plan_tol = (yyvsp[-1].fval); }
#line 11862 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1000:
#line 4135 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.noCoarse = 1; }
#line 11868 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1001:
#line 4137 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 1) 
            (yyval.scntl)->fetiInfo.nonLocalQ = 0;
          else if((yyvsp[-1].ival) == 2) {
            (yyval.scntl)->fetiInfo.nonLocalQ = 1;
            if((yyval.scntl)->fetiInfo.version == FetiInfo::feti2) {
              (yyval.scntl)->fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used"
                             " with FETI 2\n");
            }
          } else if((yyvsp[-1].ival) == 3) {
            (yyval.scntl)->fetiInfo.nonLocalQ = 1;
            (yyval.scntl)->fetiInfo.nQ = 3;
          } else if((yyvsp[-1].ival) == 4) {
            (yyval.scntl)->fetiInfo.nonLocalQ = 1;
            (yyval.scntl)->fetiInfo.nQ = 4;
          } else
            fprintf(stderr," *** WARNING: This projector does not exist,"
                           " using basic projector\n"); }
#line 11891 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1002:
#line 4156 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 1; 
          (yyval.scntl)->fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11898 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1003:
#line 4159 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11904 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1004:
#line 4161 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 2;
          (yyval.scntl)->fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11911 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1005:
#line 4164 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11917 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1006:
#line 4166 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 2;
          (yyval.scntl)->fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11924 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1007:
#line 4169 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11930 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1008:
#line 4171 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_element = true; }
#line 11936 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1009:
#line 4173 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_element = true; }
#line 11942 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1010:
#line 4175 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_corner = (yyvsp[-1].ival); }
#line 11948 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1011:
#line 4177 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.splitLocalFsi = false; }
#line 11954 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1012:
#line 4179 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.wetcorners = true; }
#line 11960 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1013:
#line 4181 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[-1].ival); }
#line 11966 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1014:
#line 4183 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[-2].ival); 
          (yyval.scntl)->fetiInfo.pick_unsafe_corners = bool((yyvsp[-1].ival)); }
#line 11973 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1015:
#line 4186 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 0) {
            (yyval.scntl)->fetiInfo.corners = FetiInfo::noCorners;
            (yyval.scntl)->fetiInfo.pickAnyCorner = 0; 
            (yyval.scntl)->fetiInfo.bmpc = true;
            (yyval.scntl)->fetiInfo.pick_unsafe_corners = false;
            (yyval.scntl)->fetiInfo.augment = FetiInfo::none;
          } }
#line 11985 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1016:
#line 4194 "p.y" /* yacc.c:1646  */
    { if((yyval.scntl)->fetiInfo.dph_flag && ((yyvsp[-1].ival) == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
          }
          else (yyval.scntl)->fetiInfo.augment = (FetiInfo::AugmentType) (yyvsp[-1].ival);

          if((yyval.scntl)->fetiInfo.augment == FetiInfo::Edges) {
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::translation;
            (yyval.scntl)->fetiInfo.nGs = 3;
          } 
          else if((yyval.scntl)->fetiInfo.augment == FetiInfo::WeightedEdges) {
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::translation;
            (yyval.scntl)->fetiInfo.nGs = 3;
          }
          else if((yyval.scntl)->fetiInfo.augment == FetiInfo::Gs) {
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::all;
            (yyval.scntl)->fetiInfo.nGs = 6;
          } }
#line 12008 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1017:
#line 4213 "p.y" /* yacc.c:1646  */
    { if((yyval.scntl)->fetiInfo.dph_flag && ((yyvsp[-2].ival) == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges;
          }
          else (yyval.scntl)->fetiInfo.augment = (FetiInfo::AugmentType) (yyvsp[-2].ival);
          if((yyval.scntl)->fetiInfo.dph_flag && ((yyvsp[-1].ival) > 2) && ((yyvsp[-1].ival) < 6)) {
            std::cerr << "WARNING: Selected rbm type is unsupported for FETI-DPH, set to translation \n";
            (yyval.scntl)->fetiInfo.rbmType = FetiInfo::translation;
          }
          else (yyval.scntl)->fetiInfo.rbmType = (FetiInfo::RbmType) (yyvsp[-1].ival);

          if((yyval.scntl)->fetiInfo.rbmType == FetiInfo::all)
            (yyval.scntl)->fetiInfo.nGs = 6;
          else if((yyval.scntl)->fetiInfo.rbmType != FetiInfo::None)
            (yyval.scntl)->fetiInfo.nGs = 3; }
#line 12028 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1018:
#line 4229 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numdir = (yyvsp[-1].ival); 
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12036 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1019:
#line 4233 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.numdir = (yyvsp[-1].ival); 
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12045 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1020:
#line 4238 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numdir = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[-1].ival);
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12054 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1021:
#line 4243 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[-1].ival);
          (yyval.scntl)->fetiInfo.numdir = (yyvsp[-2].ival);
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12064 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1022:
#line 4249 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-1].ival) == 2)
           (yyval.scntl)->fetiInfo.augmentimpl = FetiInfo::Primal;
        }
#line 12072 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1023:
#line 4253 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.orthotol = (yyvsp[-1].fval); }
#line 12078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1024:
#line 4255 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.orthotol = (yyvsp[-2].fval); 
          (yyval.scntl)->fetiInfo.orthotol2 = (yyvsp[-1].fval); }
#line 12085 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1025:
#line 4258 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.grbm_tol = (yyvsp[-1].fval); }
#line 12091 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1026:
#line 4260 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.crbm_tol = (yyvsp[-1].fval); }
#line 12097 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1027:
#line 4262 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cct_tol = (yyvsp[-1].fval); }
#line 12103 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1028:
#line 4264 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.rebuildcct = int((yyvsp[-1].ival)); }
#line 12109 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1029:
#line 4266 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.uproj = (yyvsp[-1].ival); }
#line 12115 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1030:
#line 4268 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.printMatLab = 1; }
#line 12121 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1031:
#line 4270 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->printMatLab = 1;
          (yyval.scntl)->printMatLabFile = (yyvsp[-1].strval); }
#line 12128 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1032:
#line 4273 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.local_cntl = (yyval.scntl)->fetiInfo.kii_cntl = (yyvsp[0].scntl); }
#line 12134 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1033:
#line 4275 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.coarse_cntl = (yyvsp[0].scntl); }
#line 12140 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1034:
#line 4277 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.auxcoarse_cntl = (yyvsp[0].scntl); }
#line 12146 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1035:
#line 4279 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cct_cntl = (yyvsp[0].scntl); }
#line 12152 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1036:
#line 4281 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 1)
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti1; 
          else if((yyvsp[-1].ival) == 2) {
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti2;
            if((yyval.scntl)->fetiInfo.nonLocalQ == 1) {
              (yyval.scntl)->fetiInfo.nonLocalQ = 0;
              (yyval.scntl)->fetiInfo.nQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } else if((yyvsp[-1].ival) == 3) {
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti3;
            if((yyval.scntl)->fetiInfo.nonLocalQ == 1) {
              (yyval.scntl)->fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } 
          else {
            (yyval.scntl)->fetiInfo.version = FetiInfo::feti1;
	    fprintf(stderr," *** WARNING: Version does not exist,"
                           " using FETI 1\n");
          } }
#line 12180 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1037:
#line 4305 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gmresResidual = true; }
#line 12186 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1038:
#line 4307 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gmresResidual = bool((yyvsp[-1].ival)); }
#line 12192 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1039:
#line 4309 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.pickAnyCorner = (yyvsp[-1].ival); }
#line 12198 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1040:
#line 4311 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.type = FetiInfo::nonlinear;
          (yyval.scntl)->fetiInfo.nlPrecFlg = 1; }
#line 12205 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1041:
#line 4314 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.type = FetiInfo::nonlinear;
	  (yyval.scntl)->fetiInfo.nlPrecFlg = (yyvsp[-1].ival); }
#line 12212 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1042:
#line 4317 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-1].ival); }
#line 12218 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1043:
#line 4319 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-1].fval); }
#line 12224 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1044:
#line 4321 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.spaceDimension = (yyvsp[-1].ival); }
#line 12230 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1045:
#line 4323 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.krylovtype = (yyvsp[-1].ival); }
#line 12236 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1046:
#line 4325 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.krylovtype = (yyvsp[-1].ival); }
#line 12242 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1047:
#line 4327 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.lumpedinterface = 1; }
#line 12248 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1048:
#line 4329 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.saveMemCoarse = 1; }
#line 12254 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1049:
#line 4331 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[-1].ival); }
#line 12260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1050:
#line 4333 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.complex_hermitian = true; }
#line 12267 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1051:
#line 4336 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpcflag = (yyvsp[-1].ival); }
#line 12273 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1052:
#line 4338 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpcflag = (yyvsp[-1].ival); }
#line 12279 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1053:
#line 4340 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-1].ival); }
#line 12285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1054:
#line 4342 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-1].ival); }
#line 12291 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1055:
#line 4344 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12298 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1056:
#line 4347 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12305 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1057:
#line 4350 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-2].ival); 
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-1].ival); }
#line 12312 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1058:
#line 4353 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-1].ival); }
#line 12319 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1059:
#line 4356 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-4].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12327 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1060:
#line 4360 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-4].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-3].ival); 
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12335 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1061:
#line 4364 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) < 1) (yyval.scntl)->fetiInfo.useMRHS = false; }
#line 12341 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1062:
#line 4366 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gamma = (yyvsp[-1].fval); }
#line 12347 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1063:
#line 4368 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.linesearch_maxit = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.linesearch_tau = (yyvsp[-1].fval); }
#line 12354 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1064:
#line 4371 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.bmpc = bool((yyvsp[-1].ival)); }
#line 12360 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1065:
#line 4373 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cmpc = bool((yyvsp[-1].ival)); }
#line 12366 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1066:
#line 4375 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.c_normalize = bool((yyvsp[-1].ival)); }
#line 12372 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1068:
#line 4382 "p.y" /* yacc.c:1646  */
    {
          geoSource->setOmega((yyvsp[-1].fval));
          StructProp sp; 
          sp.kappaHelm = (yyvsp[-1].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
#line 12385 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1069:
#line 4393 "p.y" /* yacc.c:1646  */
    { if(!(yyvsp[-1].copt).lagrangeMult && (yyvsp[-1].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[-1].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[-1].copt).constraint_hess; 
          domain->solInfo().constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps; }
#line 12395 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1070:
#line 4399 "p.y" /* yacc.c:1646  */
    { if(!(yyvsp[-1].copt).lagrangeMult && (yyvsp[-1].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[-1].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[-1].copt).constraint_hess;
          domain->solInfo().constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps; }
#line 12405 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1071:
#line 4407 "p.y" /* yacc.c:1646  */
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        }
#line 12416 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1072:
#line 4414 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[0].fval); }
#line 12427 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1073:
#line 4421 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[-1].fval);
          domain->solInfo().coefFilterTol = (yyvsp[0].fval); }
#line 12439 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1074:
#line 4429 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[-2].fval); 
          domain->solInfo().coefFilterTol = (yyvsp[-1].fval);
          domain->solInfo().rhsZeroTol = (yyvsp[0].fval); }
#line 12452 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1075:
#line 4438 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[-3].fval);
          domain->solInfo().coefFilterTol = (yyvsp[-2].fval); 
          domain->solInfo().rhsZeroTol = (yyvsp[-1].fval);
          domain->solInfo().inconsistentTol = (yyvsp[0].fval); }
#line 12466 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1076:
#line 4448 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12476 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1077:
#line 4454 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12486 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1078:
#line 4460 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12496 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1079:
#line 4466 "p.y" /* yacc.c:1646  */
    { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12506 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1080:
#line 4472 "p.y" /* yacc.c:1646  */
    { (yyval.copt).constraint_hess = (yyvsp[0].ival);
          (yyval.copt).constraint_hess_eps = 0; }
#line 12513 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1081:
#line 4475 "p.y" /* yacc.c:1646  */
    { (yyval.copt).constraint_hess = (yyvsp[-1].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[0].fval); }
#line 12520 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1082:
#line 4480 "p.y" /* yacc.c:1646  */
    { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
#line 12529 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1083:
#line 4497 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-1].ival);
          domain->curvatureFlag = 0;
        }
#line 12538 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1084:
#line 4502 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-2].ival);
          domain->curvatureConst1 = (yyvsp[-1].fval);
          domain->curvatureFlag = 1;
        }
#line 12548 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1085:
#line 4508 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-3].ival);
          domain->curvatureConst1 = (yyvsp[-2].fval);
          domain->curvatureConst2 = (yyvsp[-1].fval);
          domain->curvatureFlag = 2;
        }
#line 12559 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1086:
#line 4515 "p.y" /* yacc.c:1646  */
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        }
#line 12568 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1087:
#line 4520 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
#line 12577 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1088:
#line 4525 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
#line 12586 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1094:
#line 4539 "p.y" /* yacc.c:1646  */
    { domain->setWaveDirections(0, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12592 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1095:
#line 4542 "p.y" /* yacc.c:1646  */
    {
          domain->setKirchhoffLocations((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 12600 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1096:
#line 4545 "p.y" /* yacc.c:1646  */
    {
          domain->setKirchhoffLocations((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 12608 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1099:
#line 4555 "p.y" /* yacc.c:1646  */
    { domain->setFFPDirections((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12614 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1101:
#line 4562 "p.y" /* yacc.c:1646  */
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[-1].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[-1].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        }
#line 12627 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1103:
#line 4576 "p.y" /* yacc.c:1646  */
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[-1].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[-1].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        }
#line 12640 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1104:
#line 4587 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        }
#line 12648 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1105:
#line 4593 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 12654 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1106:
#line 4595 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 12660 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1107:
#line 4598 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
            domain->solInfo().probType = SolverInfo::NonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::Dynamic)
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
            domain->solInfo().order = 1;
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          }
          domain->solInfo().solvercntl->fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().getNLInfo().setDefaults(); /* just in case PIECEWISE is used under statics */
          domain->solInfo().nlFlag = 1; /* can be used for decomposition when a different treatment is required, e.g. contact */ }
#line 12676 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1108:
#line 4610 "p.y" /* yacc.c:1646  */
    {}
#line 12682 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1109:
#line 4612 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength; }
#line 12689 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1110:
#line 4615 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam; }
#line 12698 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1111:
#line 4620 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linearelastic = 1; }
#line 12704 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1112:
#line 4622 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linearelastic = (yyvsp[-1].ival); }
#line 12710 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1113:
#line 4624 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[-1].ival); }
#line 12716 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1114:
#line 4626 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-1].fval); }
#line 12722 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1115:
#line 4628 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[-1].fval); }
#line 12729 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1116:
#line 4631 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-4].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[-1].fval); }
#line 12738 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1117:
#line 4636 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[-1].fval); }
#line 12745 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1118:
#line 4639 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[-4].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[-2].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[-1].ival); }
#line 12754 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1119:
#line 4644 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[-1].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[-1].ival); }
#line 12761 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1120:
#line 4647 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[-2].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[-1].ival); }
#line 12768 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1121:
#line 4650 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().unsymmetric = true; }
#line 12774 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1122:
#line 4652 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[-1].fval); }
#line 12780 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1123:
#line 4654 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-1].ival); }
#line 12786 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1124:
#line 4656 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-2].ival); 
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-1].ival); }
#line 12793 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1125:
#line 4659 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-3].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-2].ival);
          // note: currently we use either c1 or c2, but never both
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-1].fval);
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-1].fval); }
#line 12803 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1126:
#line 4665 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-4].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-3].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-2].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[-1].fval); }
#line 12813 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1127:
#line 4671 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-5].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-4].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-3].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().linesearch.verbose = bool((yyvsp[-1].ival)); }
#line 12824 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1128:
#line 4678 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().failsafe = true; }
#line 12830 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1129:
#line 4680 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[-1].fval); }
#line 12837 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1130:
#line 4683 "p.y" /* yacc.c:1646  */
    { domain->solInfo().num_penalty_its = (yyvsp[-3].ival); 
          domain->solInfo().penalty_tol = (yyvsp[-2].fval);
          domain->solInfo().penalty_beta = (yyvsp[-1].fval); }
#line 12845 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1131:
#line 4687 "p.y" /* yacc.c:1646  */
    { domain->solInfo().num_penalty_its = (yyvsp[-5].ival);
          domain->solInfo().penalty_tol = (yyvsp[-4].fval);
          domain->solInfo().penalty_beta = (yyvsp[-3].fval);
          domain->solInfo().reinit_lm = bool((yyvsp[-2].ival));
          domain->solInfo().lm_update_flag = (yyvsp[-1].ival); }
#line 12855 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1132:
#line 4693 "p.y" /* yacc.c:1646  */
    { domain->solInfo().numberOfRomCPUs = (yyvsp[-1].ival); }
#line 12861 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1134:
#line 4698 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-1].ival)); 
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear; 
        }
#line 12870 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1135:
#line 4703 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-2].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[-1].ival);
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
#line 12880 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1136:
#line 4709 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-3].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[-2].ival);
          domain->solInfo().piecewise_contact = bool((yyvsp[-1].ival));
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
#line 12891 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1137:
#line 4718 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setReOrtho(); }
#line 12897 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1139:
#line 4723 "p.y" /* yacc.c:1646  */
    { geoSource->setControl((yyvsp[-7].strval),(yyvsp[-3].strval),(yyvsp[-1].strval)); domain->solInfo().soltyp = (yyvsp[-5].ival); }
#line 12903 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1140:
#line 4725 "p.y" /* yacc.c:1646  */
    { geoSource->setControl((yyvsp[-9].strval),(yyvsp[-5].strval),(yyvsp[-3].strval),(yyvsp[-1].strval)); domain->solInfo().soltyp = (yyvsp[-7].ival); }
#line 12909 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1142:
#line 4736 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contact_mode = (yyvsp[-1].ival); }
#line 12915 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1143:
#line 4740 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12921 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1144:
#line 4743 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval));}
#line 12927 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1145:
#line 4746 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 0.0, (yyvsp[-1].ival));}
#line 12933 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1146:
#line 4748 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), 0.0, -1, (yyvsp[-1].copt).lagrangeMult, (yyvsp[-1].copt).penalty);}
#line 12939 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1147:
#line 4751 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-9].ival)-1, (yyvsp[-8].ival)-1, (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-1].ival));}
#line 12945 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1148:
#line 4753 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-10].ival)-1, (yyvsp[-9].ival)-1, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-4].fval), (yyvsp[-2].ival), (yyvsp[-1].copt).lagrangeMult, (yyvsp[-1].copt).penalty);}
#line 12951 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1150:
#line 4758 "p.y" /* yacc.c:1646  */
    { 
           geoSource->addMaterial((yyvsp[-7].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12960 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1151:
#line 4763 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new BilinPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12969 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1152:
#line 4768 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new BilinPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12978 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1153:
#line 4773 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new BilinPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new BilinPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 12992 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1154:
#line 4783 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new BilinPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new BilinPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval)) );
           }
         }
#line 13007 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1155:
#line 4794 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new BilinPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new BilinPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
           }
         }
#line 13022 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1156:
#line 4805 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13031 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1157:
#line 4810 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13040 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1158:
#line 4815 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13049 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1159:
#line 4820 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13063 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1160:
#line 4830 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval)) );
           }
         }
#line 13078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1161:
#line 4841 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), 414) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new FiniteStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
           }
         }
#line 13093 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1162:
#line 4852 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-10].ival)-1, new CrushableFoam((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
          }
#line 13101 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1163:
#line 4856 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-11].ival)-1, new CrushableFoam((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
          }
#line 13109 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1164:
#line 4860 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-12].ival)-1, new CrushableFoam((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
          }
#line 13117 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1165:
#line 4864 "p.y" /* yacc.c:1646  */
    {
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              geoSource->addMaterial((yyvsp[-13].ival)-1, new CrushableFoam((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
              domain->solInfo().elementDeletion = true;
            }
            else {
              geoSource->addMaterial((yyvsp[-13].ival)-1, new CrushableFoam((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 
                                     (yyvsp[-2].fval), std::numeric_limits<double>::infinity()) );
            }
          }
#line 13132 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1166:
#line 4875 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13141 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1167:
#line 4880 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13150 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1168:
#line 4885 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13159 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1169:
#line 4890 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13173 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1170:
#line 4900 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval)) );
           }
         }
#line 13188 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1171:
#line 4911 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival)) );
           }
         }
#line 13203 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1172:
#line 4922 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new ElaLinIsoMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13212 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1173:
#line 4927 "p.y" /* yacc.c:1646  */
    { 
           geoSource->addMaterial((yyvsp[-5].ival)-1, 
             new ElaLinIsoMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
	 }
#line 13221 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1174:
#line 4932 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new ElaLinIsoMat((yyvsp[-1].fval)));
         }
#line 13230 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1175:
#line 4937 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13240 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1176:
#line 4943 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13250 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1177:
#line 4949 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1178:
#line 4955 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13269 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1179:
#line 4960 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13278 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1180:
#line 4965 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-1].fval)));
         }
#line 13287 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1181:
#line 4970 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13297 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1182:
#line 4976 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13307 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1183:
#line 4982 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13317 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1184:
#line 4988 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new HenckyMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13326 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1185:
#line 4993 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new HenckyMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13335 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1186:
#line 4998 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new HenckyMat((yyvsp[-1].fval)));
         }
#line 13344 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1187:
#line 5003 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13354 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1188:
#line 5009 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13364 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1189:
#line 5015 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13374 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1190:
#line 5021 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new ElaLinIsoMat2D((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0));
         }
#line 13383 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1191:
#line 5026 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new ElaLinIsoMat2D((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13392 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1192:
#line 5031 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0));
         }
#line 13401 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1193:
#line 5036 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13410 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1194:
#line 5041 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FabricMap((yyvsp[-5].fval), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].fval), 0, 0, FabricMap::GREEN_LAGRANGE));
         }
#line 13419 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1195:
#line 5046 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMap((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMap::GREEN_LAGRANGE));
         }
#line 13428 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1196:
#line 5051 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FabricMap((yyvsp[-5].fval), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].fval), 0, 0, FabricMap::INFINTESIMAL));
         }
#line 13437 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1197:
#line 5056 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMap((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMap::INFINTESIMAL));
         }
#line 13446 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1198:
#line 5061 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMat((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0, FabricMat::GREEN_LAGRANGE));
         }
#line 13455 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1199:
#line 5066 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new FabricMat((yyvsp[-9].fval), (yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMat::GREEN_LAGRANGE));
         }
#line 13464 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1200:
#line 5071 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMat((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0, FabricMat::INFINTESIMAL));
         }
#line 13473 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1201:
#line 5076 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new FabricMat((yyvsp[-9].fval), (yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMat::INFINTESIMAL));
         }
#line 13482 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1202:
#line 5081 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<ElaLinIsoMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13491 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1203:
#line 5086 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<ElaLinIsoMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13500 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1204:
#line 5091 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13510 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1205:
#line 5097 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13520 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1206:
#line 5103 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<StVenantKirchhoffMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13529 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1207:
#line 5108 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<StVenantKirchhoffMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13538 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1208:
#line 5113 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13548 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1209:
#line 5119 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13558 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1210:
#line 5125 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<NeoHookeanMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13567 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1211:
#line 5130 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<NeoHookeanMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13577 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1212:
#line 5136 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new PlaneStressMat<MooneyRivlinMat>((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13586 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1213:
#line 5141 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<BrittleFractureTB<MooneyRivlinMat> >((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13596 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1214:
#line 5147 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13605 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1215:
#line 5152 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13614 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1216:
#line 5157 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13623 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1217:
#line 5162 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval), (yyvsp[-2].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13637 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1218:
#line 5172 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
           }
         }
#line 13652 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1219:
#line 5183 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
           }
         }
#line 13667 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1220:
#line 5194 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13676 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1221:
#line 5199 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13685 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1222:
#line 5204 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13694 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1223:
#line 5209 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval), (yyvsp[-2].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13708 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1224:
#line 5219 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-13].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-1].fval), (yyvsp[-3].fval)) );
           }
         }
#line 13723 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1225:
#line 5230 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-3].fval) > 0 && (yyvsp[-3].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-14].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval),
                                    std::numeric_limits<double>::infinity(), (yyvsp[-2].fval), (yyvsp[-1].ival), (yyvsp[-4].fval)) );
           }
         }
#line 13738 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1226:
#line 5241 "p.y" /* yacc.c:1646  */
    {
            double params[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-5].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          }
#line 13748 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1227:
#line 5247 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<ElaLinIsoMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13761 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1228:
#line 5256 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PronyViscoElastic<ElaLinIsoMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13774 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1229:
#line 5265 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13788 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1230:
#line 5275 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13802 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1231:
#line 5285 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PronyViscoElastic<ElaLinIsoMat2D>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13815 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1232:
#line 5294 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-14].ival)-1,
             new PronyViscoElastic<ElaLinIsoMat2D>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13828 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1233:
#line 5303 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13841 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1234:
#line 5312 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-14].ival)-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13854 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1235:
#line 5321 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-13].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-11].fval), (yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-7].fval), 0, 0, FabricMap::GREEN_LAGRANGE, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13867 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1236:
#line 5330 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMap::GREEN_LAGRANGE, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13880 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1237:
#line 5339 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-13].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-11].fval), (yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-7].fval), 0, 0, FabricMap::INFINTESIMAL, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13893 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1238:
#line 5348 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMap::INFINTESIMAL, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13906 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1239:
#line 5357 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, FabricMat::GREEN_LAGRANGE, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13919 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1240:
#line 5366 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-17].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-15].fval), (yyvsp[-14].ival), (yyvsp[-13].ival), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMat::GREEN_LAGRANGE, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13932 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1241:
#line 5375 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, FabricMat::INFINTESIMAL, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13945 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1242:
#line 5384 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-17].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-15].fval), (yyvsp[-14].ival), (yyvsp[-13].ival), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMat::INFINTESIMAL, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13958 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1243:
#line 5393 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13971 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1244:
#line 5402 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-14].ival)-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13984 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1245:
#line 5411 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval); 
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13998 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1246:
#line 5421 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-18].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >((yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14012 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1247:
#line 5431 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<StVenantKirchhoffMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14025 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1248:
#line 5440 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PronyViscoElastic<StVenantKirchhoffMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14038 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1249:
#line 5449 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14052 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1250:
#line 5459 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14066 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1251:
#line 5469 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14079 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1252:
#line 5478 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-14].ival)-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14092 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1253:
#line 5487 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14106 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1254:
#line 5497 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-18].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >((yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14120 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1255:
#line 5507 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-5].ival)-1,
              new NeoHookeanMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14129 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1256:
#line 5512 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new BrittleFractureTB<NeoHookeanMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14139 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1257:
#line 5518 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<NeoHookeanMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14152 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1258:
#line 5527 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14166 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1259:
#line 5537 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<NeoHookeanMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14179 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1260:
#line 5546 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14193 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1261:
#line 5556 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-6].ival)-1,
              new MooneyRivlinMat((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14202 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1262:
#line 5561 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new BrittleFractureTB<MooneyRivlinMat>((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14212 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1263:
#line 5567 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PronyViscoElastic<MooneyRivlinMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14225 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1264:
#line 5576 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14239 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1265:
#line 5586 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PlaneStressMat<PronyViscoElastic<MooneyRivlinMat> >((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14252 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1266:
#line 5595 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> > >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14266 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1267:
#line 5605 "p.y" /* yacc.c:1646  */
    {
            double mu[1] = { (yyvsp[-3].fval) };
            double alpha[1] = { (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-6].ival)-1, new OgdenMat((yyvsp[-4].fval), mu, alpha, K));
          }
#line 14277 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1268:
#line 5612 "p.y" /* yacc.c:1646  */
    {
            double mu[2] = { (yyvsp[-5].fval), (yyvsp[-4].fval) };
            double alpha[2] = { (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-8].ival)-1, new OgdenMat((yyvsp[-6].fval), mu, alpha, K));
          }
#line 14288 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1269:
#line 5619 "p.y" /* yacc.c:1646  */
    {
            double mu[2] = { (yyvsp[-6].fval), (yyvsp[-5].fval) }; 
            double alpha[2] = { (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-9].ival)-1, new OgdenMat((yyvsp[-7].fval), mu, alpha, K));
          }
#line 14299 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1270:
#line 5626 "p.y" /* yacc.c:1646  */
    {
            double mu[3] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval) };
            double alpha[3] = { (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-10].ival)-1, new OgdenMat((yyvsp[-8].fval), mu, alpha, K));
          }
#line 14310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1271:
#line 5633 "p.y" /* yacc.c:1646  */
    {
            double mu[3] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval) };
            double alpha[3] = { (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-11].ival)-1, new OgdenMat((yyvsp[-9].fval), mu, alpha, K));
          }
#line 14321 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1272:
#line 5640 "p.y" /* yacc.c:1646  */
    {
            double mu[4] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval) };
            double alpha[4] = { (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-12].ival)-1, new OgdenMat((yyvsp[-10].fval), mu, alpha, K));
          }
#line 14332 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1273:
#line 5647 "p.y" /* yacc.c:1646  */
    {
            double mu[4] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval) };
            double alpha[4] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-13].ival)-1, new OgdenMat((yyvsp[-11].fval), mu, alpha, K));
          }
#line 14343 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1274:
#line 5654 "p.y" /* yacc.c:1646  */
    {
            double mu[5] = { (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval) };
            double alpha[5] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-14].ival)-1, new OgdenMat((yyvsp[-12].fval), mu, alpha, K));
          }
#line 14354 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1275:
#line 5661 "p.y" /* yacc.c:1646  */
    {
            double mu[5] = { (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval) };
            double alpha[5] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-15].ival)-1, new OgdenMat((yyvsp[-13].fval), mu, alpha, K));
          }
#line 14365 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1276:
#line 5668 "p.y" /* yacc.c:1646  */
    {
            double mu[6] = { (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval) };
            double alpha[6] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-16].ival)-1, new OgdenMat((yyvsp[-14].fval), mu, alpha, K));
          }
#line 14376 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1277:
#line 5675 "p.y" /* yacc.c:1646  */
    {
            double mu[6] = { (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval) };
            double alpha[6] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-17].ival)-1, new OgdenMat((yyvsp[-15].fval), mu, alpha, K));
          }
#line 14387 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1278:
#line 5682 "p.y" /* yacc.c:1646  */
    {
            double mu[7] = { (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval) };
            double alpha[7] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-18].ival)-1, new OgdenMat((yyvsp[-16].fval), mu, alpha, K));
          }
#line 14398 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1279:
#line 5689 "p.y" /* yacc.c:1646  */
    {
            double mu[7] = { (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval) };
            double alpha[7] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-19].ival)-1, new OgdenMat((yyvsp[-17].fval), mu, alpha, K));
          }
#line 14409 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1280:
#line 5696 "p.y" /* yacc.c:1646  */
    {
            double mu[8] = { (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval) };
            double alpha[8] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-20].ival)-1, new OgdenMat((yyvsp[-18].fval), mu, alpha, K));
          }
#line 14420 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1281:
#line 5703 "p.y" /* yacc.c:1646  */
    {
            double mu[8] = { (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval) };
            double alpha[8] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-21].ival)-1, new OgdenMat((yyvsp[-19].fval), mu, alpha, K));
          }
#line 14431 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1282:
#line 5710 "p.y" /* yacc.c:1646  */
    {
            double mu[9] = { (yyvsp[-19].fval), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval) };
            double alpha[9] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-22].ival)-1, new OgdenMat((yyvsp[-20].fval), mu, alpha, K));
          }
#line 14442 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1283:
#line 5717 "p.y" /* yacc.c:1646  */
    {
            double mu[9] = { (yyvsp[-20].fval), (yyvsp[-19].fval), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval) };
            double alpha[9] = { (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-23].ival)-1, new OgdenMat((yyvsp[-21].fval), mu, alpha, K));
          }
#line 14453 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1284:
#line 5724 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new SimoElasticMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14462 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1285:
#line 5729 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new SimoPlasticMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 14471 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1286:
#line 5734 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new SimoPlasticMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 14480 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1287:
#line 5739 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-8].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
#line 14490 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1288:
#line 5745 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
#line 14500 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1289:
#line 5751 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0. };
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14513 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1290:
#line 5760 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), double((yyvsp[-1].ival)) };
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14526 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1291:
#line 5769 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-8].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
#line 14536 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1292:
#line 5775 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
#line 14546 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1293:
#line 5781 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0. };
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14559 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1294:
#line 5790 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), double((yyvsp[-1].ival)) };
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14572 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1295:
#line 5799 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new ExpMat((yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14581 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1296:
#line 5804 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new ExpMat((yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14590 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1297:
#line 5809 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new ExpMat((yyvsp[-6].ival), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14599 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1298:
#line 5814 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new ExpMat((yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14608 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1299:
#line 5819 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new ExpMat((yyvsp[-8].ival), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14617 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1300:
#line 5824 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new ExpMat((yyvsp[-9].ival), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14629 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1301:
#line 5832 "p.y" /* yacc.c:1646  */
    {
           ExpMat *mat = new ExpMat((yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval));
           mat->yssrtid = (yyvsp[-1].ival);
           geoSource->addMaterial((yyvsp[-11].ival)-1, mat);
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14642 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1302:
#line 5842 "p.y" /* yacc.c:1646  */
    {
           ExpMat *mat = new ExpMat((yyvsp[-19].ival), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval));
           mat->yssrtid = (yyvsp[-10].ival);
           mat->optcor0 = (yyvsp[-9].ival);
           mat->optcor1 = (yyvsp[-8].ival);
           mat->optprj = (yyvsp[-7].ival);
           mat->opthgc = (yyvsp[-6].ival);
           mat->prmhgc[0] = (yyvsp[-5].fval);
           mat->prmhgc[1] = (yyvsp[-4].fval);
           mat->prmhgc[2] = (yyvsp[-3].fval);
           mat->ngqpt2 = (yyvsp[-2].ival);
           mat->ematpro[18] = (yyvsp[-1].fval);
           geoSource->addMaterial((yyvsp[-20].ival)-1, mat);
           if((yyvsp[-11].fval) > 0 && (yyvsp[-11].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14664 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1303:
#line 5860 "p.y" /* yacc.c:1646  */
    {
	   geoSource->loadMaterial((yyvsp[-2].strval), (yyvsp[-1].strval));
	 }
#line 14672 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1304:
#line 5864 "p.y" /* yacc.c:1646  */
    {
	   geoSource->addMaterial((yyvsp[-3].ival)-1, (yyvsp[-2].strval), (yyvsp[-1].dlist));
	 }
#line 14680 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1306:
#line 5871 "p.y" /* yacc.c:1646  */
    { geoSource->setMatUsage((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 14686 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1307:
#line 5873 "p.y" /* yacc.c:1646  */
    {
            for(int i = (yyvsp[-3].ival)-1; i < (yyvsp[-2].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[-1].ival)-1);
	  }
#line 14695 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1308:
#line 5879 "p.y" /* yacc.c:1646  */
    { (yyval.dlist).nval = 0; }
#line 14701 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1309:
#line 5881 "p.y" /* yacc.c:1646  */
    { 
          if((yyvsp[-1].dlist).nval == 64) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[-1].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[0].fval);
 	}
#line 14714 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1310:
#line 5891 "p.y" /* yacc.c:1646  */
    { (yyval.slist).nval = 0; }
#line 14720 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1311:
#line 5893 "p.y" /* yacc.c:1646  */
    { 
          if((yyvsp[-1].slist).nval == 32) {
             fprintf(stderr, "Too many files!\n");
	     exit(-1);
          }
          (yyval.slist) = (yyvsp[-1].slist);
          (yyval.slist).v[(yyval.slist).nval++] = (yyvsp[0].strval);
 	}
#line 14733 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1312:
#line 5904 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-1].ival));
          domain->solInfo().setSparseRenum((yyvsp[-1].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[-1].ival)); }
#line 14741 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1313:
#line 5908 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-3].ival));
          domain->solInfo().setSparseRenum((yyvsp[-1].ival)); }
#line 14748 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1314:
#line 5911 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-5].ival));
          domain->solInfo().setSparseRenum((yyvsp[-3].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[-1].ival)); }
#line 14756 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1315:
#line 5918 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().svdPodRom = true;}
#line 14764 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1317:
#line 5926 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().snapfiPodRom.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14770 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1318:
#line 5937 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().velocPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14776 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1319:
#line 5939 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().accelPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14782 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1320:
#line 5941 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().dsvPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14788 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1321:
#line 5943 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().muvPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14794 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1322:
#line 5945 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[0].ival); }
#line 14800 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1323:
#line 5947 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14806 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1324:
#line 5949 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14812 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1325:
#line 5951 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14818 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1326:
#line 5953 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subtractRefPodRom = true;
    domain->solInfo().readInLocalBasesCent.push_back(std::string((yyvsp[0].strval))); }
#line 14825 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1327:
#line 5956 "p.y" /* yacc.c:1646  */
    { domain->solInfo().flagss = (yyvsp[0].ival); }
#line 14831 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1328:
#line 5958 "p.y" /* yacc.c:1646  */
    { domain->solInfo().flagss = (yyvsp[-1].ival); domain->solInfo().flagrs = (yyvsp[0].ival); }
#line 14837 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1329:
#line 5960 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[0].ival); }
#line 14843 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1330:
#line 5962 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[-1].ival);
    domain->solInfo().skipOffSet = (yyvsp[0].ival); }
#line 14850 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1331:
#line 5965 "p.y" /* yacc.c:1646  */
    { domain->solInfo().robcSolve = bool((yyvsp[0].ival)); }
#line 14856 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1332:
#line 5967 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().robfi.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14862 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1333:
#line 5969 "p.y" /* yacc.c:1646  */
    { domain->solInfo().svdBlockSize = (yyvsp[0].ival); }
#line 14868 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1334:
#line 5971 "p.y" /* yacc.c:1646  */
    { domain->solInfo().romEnergy = (yyvsp[0].fval); }
#line 14874 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1335:
#line 5984 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf            = 3;
    domain->solInfo().nmfMaxIter         = (yyvsp[-3].ival);
    domain->solInfo().nmfTol             = (yyvsp[-2].fval);
    domain->solInfo().nmfPqnNumInnerIter = (yyvsp[-1].ival);
    domain->solInfo().nmfPqnAlpha        = (yyvsp[0].fval); }
#line 14884 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1336:
#line 5990 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = (yyvsp[-4].ival);
    domain->solInfo().nmfDelROBDim = (yyvsp[-3].ival);
    domain->solInfo().nmfRandInit  = (yyvsp[-2].ival);
    domain->solInfo().nmfMaxIter   = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 14895 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1337:
#line 5997 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 1;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol     = (yyvsp[0].fval); }
#line 14903 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1338:
#line 6001 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 4;   
    domain->solInfo().nmfMaxIter = (yyvsp[-4].ival); 
    domain->solInfo().nmfTol     = (yyvsp[-3].fval); 
    domain->solInfo().nmfcAlpha  = (yyvsp[-2].fval);
    domain->solInfo().nmfcBeta   = (yyvsp[-1].fval);
    domain->solInfo().nmfcGamma  = (yyvsp[0].fval);}
#line 14914 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1339:
#line 6008 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 4;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol     = (yyvsp[0].fval); }
#line 14922 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1340:
#line 6012 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nmfNumSub = (yyvsp[0].ival); }
#line 14928 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1341:
#line 6014 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 2; }
#line 14934 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1342:
#line 6016 "p.y" /* yacc.c:1646  */
    { domain->solInfo().clustering = (yyvsp[0].ival); }
#line 14940 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1343:
#line 6018 "p.y" /* yacc.c:1646  */
    { domain->solInfo().clustering = (yyvsp[-1].ival); 
    domain->solInfo().clusterSubspaceAngle = true; }
#line 14947 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1344:
#line 6021 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[0].ival); }
#line 14953 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1345:
#line 6023 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[-1].ival);
    domain->solInfo().tolPodRom = (yyvsp[0].fval);}
#line 14960 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1346:
#line 6026 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[-2].ival);
    domain->solInfo().tolPodRom = (yyvsp[-1].fval);
    domain->solInfo().solverTypeSpnnls = (yyvsp[0].ival); }
#line 14968 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1347:
#line 6030 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hotstartSample = bool((yyvsp[0].ival)); }
#line 14974 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1348:
#line 6032 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rowClustering = (yyvsp[0].ival); }
#line 14980 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1349:
#line 6034 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rowClustering = (yyvsp[-1].ival);
    domain->solInfo().clusterSubspaceAngle = true; }
#line 14987 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1351:
#line 6041 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().DEIMBasisPod = true; }
#line 14995 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1353:
#line 6049 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().UDEIMBasisPod = true; }
#line 15003 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1355:
#line 6057 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().samplingPodRom = true; }
#line 15011 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1360:
#line 6068 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().snapProjPodRom = true; }
#line 15019 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1362:
#line 6076 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInROBorModes.push_back((yyvsp[0].strval)); }
#line 15025 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1363:
#line 6078 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInROBorModes.push_back((yyvsp[-1].strval));
    domain->solInfo().localBasisSize.push_back((yyvsp[0].ival)); }
#line 15032 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1364:
#line 6081 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInDualROB.push_back((yyvsp[-1].strval));
    domain->solInfo().localDualBasisSize.push_back((yyvsp[0].ival)); }
#line 15039 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1365:
#line 6084 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[0].strval)); }
#line 15045 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1366:
#line 6086 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[-1].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[0].strval)); }
#line 15052 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1367:
#line 6089 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[-2].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[-1].strval));
    domain->solInfo().accelPodRomFile.push_back((yyvsp[0].strval)); }
#line 15060 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1368:
#line 6093 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolPodRom = (yyvsp[0].fval); }
#line 15066 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1369:
#line 6095 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[0].ival); }
#line 15072 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1370:
#line 6097 "p.y" /* yacc.c:1646  */
    { domain->solInfo().randomSampleSize = (yyvsp[0].ival); 
    domain->solInfo().randomVecSampling = true; }
#line 15079 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1371:
#line 6100 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipOffSet = (yyvsp[0].ival); }
#line 15085 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1372:
#line 6102 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[0].ival); }
#line 15091 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1373:
#line 6104 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[-1].ival); 
    domain->solInfo().forcePodSize = (yyvsp[0].ival);}
#line 15098 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1374:
#line 6107 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[-2].ival); 
    domain->solInfo().forcePodSize = (yyvsp[-1].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); }
#line 15106 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1375:
#line 6111 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[0].ival)); }
#line 15112 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1376:
#line 6113 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useConstantMassForces = bool((yyvsp[0].ival)); }
#line 15118 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1377:
#line 6115 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stackedElementSampling = bool((yyvsp[0].ival)); }
#line 15124 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1378:
#line 6117 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassOrthogonalProjection = bool((yyvsp[0].ival)); }
#line 15130 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1379:
#line 6119 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = bool((yyvsp[0].ival)); }
#line 15136 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1380:
#line 6121 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = bool((yyvsp[0].ival)); }
#line 15142 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1381:
#line 6123 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15148 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1382:
#line 6125 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[-1].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15155 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1383:
#line 6128 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[-2].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[-1].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15163 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1384:
#line 6132 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useScalingSpnnls = bool((yyvsp[0].ival)); }
#line 15169 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1385:
#line 6134 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useCenterSpnnls = bool((yyvsp[0].ival)); }
#line 15175 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1386:
#line 6136 "p.y" /* yacc.c:1646  */
    { domain->solInfo().projectSolution = bool((yyvsp[0].ival)); }
#line 15181 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1387:
#line 6138 "p.y" /* yacc.c:1646  */
    { domain->solInfo().positiveElements = bool((yyvsp[0].ival)); }
#line 15187 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1388:
#line 6140 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hotstartSample = bool((yyvsp[0].ival)); }
#line 15193 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1389:
#line 6142 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeSpnnls = (yyvsp[0].ival); }
#line 15199 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1390:
#line 6144 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizeSpnnls = (yyvsp[0].fval); }
#line 15205 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1391:
#line 6146 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxElemSpnnls = (yyvsp[0].ival); }
#line 15211 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1392:
#line 6148 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxIterSpnnls = (yyvsp[0].fval); }
#line 15217 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1393:
#line 6150 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[0].strval); }
#line 15223 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1394:
#line 6152 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[-1].strval);
    domain->solInfo().forcePodSize = (yyvsp[0].ival); }
#line 15230 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1395:
#line 6155 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[-2].strval); 
    domain->solInfo().forcePodSize = (yyvsp[-1].ival); 
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); }
#line 15238 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1396:
#line 6159 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[0].strval); 
    domain->solInfo().ConstraintBasisPod = true;}
#line 15245 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1397:
#line 6162 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[-1].strval);
    domain->solInfo().constraintPodSize = (yyvsp[0].ival); 
    domain->solInfo().ConstraintBasisPod = true; }
#line 15253 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1398:
#line 6166 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[-2].strval);
    domain->solInfo().constraintPodSize = (yyvsp[-1].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); 
    domain->solInfo().ConstraintBasisPod = true; }
#line 15262 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1399:
#line 6171 "p.y" /* yacc.c:1646  */
    { domain->solInfo().filterSnapshotRows = bool((yyvsp[0].ival)); }
#line 15268 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1400:
#line 6173 "p.y" /* yacc.c:1646  */
    { domain->solInfo().selectFullNode = bool((yyvsp[0].ival)); }
#line 15274 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1401:
#line 6175 "p.y" /* yacc.c:1646  */
    { domain->solInfo().selectFullElem = bool((yyvsp[0].ival)); }
#line 15280 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1402:
#line 6177 "p.y" /* yacc.c:1646  */
    { domain->solInfo().computeForceSnap = bool((yyvsp[0].ival)); }
#line 15286 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1403:
#line 6179 "p.y" /* yacc.c:1646  */
    { domain->solInfo().computeConstraintSnap = bool((yyvsp[0].ival)); }
#line 15292 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1404:
#line 6181 "p.y" /* yacc.c:1646  */
    { domain->solInfo().orthogForceSnap = bool((yyvsp[0].ival)); }
#line 15298 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1405:
#line 6183 "p.y" /* yacc.c:1646  */
    { domain->solInfo().orthogConstraintSnap = bool((yyvsp[0].ival)); }
#line 15304 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1406:
#line 6185 "p.y" /* yacc.c:1646  */
    { domain->solInfo().npMax = (yyvsp[0].ival); }
#line 15310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1407:
#line 6187 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scpkMB= (yyvsp[-1].ival);
    domain->solInfo().scpkNB= (yyvsp[0].ival); }
#line 15317 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1408:
#line 6190 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scpkMP= (yyvsp[-1].ival);
    domain->solInfo().scpkNP= (yyvsp[0].ival); }
#line 15324 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1409:
#line 6193 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useReverseOrder = bool((yyvsp[0].ival)); }
#line 15330 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1410:
#line 6195 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = (yyvsp[-4].ival);
    domain->solInfo().nmfDelROBDim = (yyvsp[-3].ival);
    domain->solInfo().nmfRandInit = (yyvsp[-2].ival);
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 15341 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1411:
#line 6202 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 15349 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1413:
#line 6210 "p.y" /* yacc.c:1646  */
    { domain->solInfo().conwepConfigurations.push_back((yyvsp[-1].blastData)); }
#line 15355 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1414:
#line 6215 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scalePosCoords = true;
     domain->solInfo().xScaleFactor = (yyvsp[-3].fval);
     domain->solInfo().yScaleFactor = (yyvsp[-2].fval);
     domain->solInfo().zScaleFactor = (yyvsp[-1].fval);}
#line 15364 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1415:
#line 6224 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePOSCFG = true; }
#line 15370 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1416:
#line 6226 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xScaleFactors.push_back((yyvsp[-3].fval));
     domain->solInfo().yScaleFactors.push_back((yyvsp[-2].fval));
     domain->solInfo().zScaleFactors.push_back((yyvsp[-1].fval)); }
#line 15378 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1417:
#line 6230 "p.y" /* yacc.c:1646  */
    { domain->solInfo().MassOrthogonalBasisFiles.push_back((yyvsp[-1].strval)); }
#line 15384 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1418:
#line 6235 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePOSCFG = true; }
#line 15390 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1419:
#line 6237 "p.y" /* yacc.c:1646  */
    { domain->solInfo().NodeTrainingFiles.push_back(std::string((yyvsp[-1].slist).v[0]));
     for(int i=1; i<(yyvsp[-1].slist).nval; ++i) domain->solInfo().MassOrthogonalBasisFiles.push_back(std::string((yyvsp[-1].slist).v[i])); }
#line 15397 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1420:
#line 6243 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().ROMPostProcess = true; }
#line 15405 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1422:
#line 6251 "p.y" /* yacc.c:1646  */
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[0].strval)); 
    domain->solInfo().numRODFile += 1; }
#line 15412 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1423:
#line 6254 "p.y" /* yacc.c:1646  */
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[-1].strval));
    domain->solInfo().numRODFile += 1; 
    domain->solInfo().skipPodRom = (yyvsp[0].ival);}
#line 15420 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1424:
#line 6258 "p.y" /* yacc.c:1646  */
    { domain->solInfo().romresidType = (yyvsp[0].ival); }
#line 15426 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1425:
#line 6263 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 15432 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1426:
#line 6265 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = std::numeric_limits<int>::max(); }
#line 15438 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1427:
#line 6270 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[0].ival); }
#line 15444 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1428:
#line 6272 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[0].fval); }
#line 15450 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1429:
#line 6274 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = std::numeric_limits<double>::infinity(); }
#line 15456 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1430:
#line 6276 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); }
#line 15462 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;


#line 15466 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  return yyresult;
}
#line 6278 "p.y" /* yacc.c:1906  */

