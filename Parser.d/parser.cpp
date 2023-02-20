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

#line 111 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:339  */

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
    LINEARELASTIC = 710,
    STVENANTKIRCHHOFF = 711,
    TULERBUTCHER = 712,
    LINPLSTRESS = 713,
    READ = 714,
    OPTCTV = 715,
    ISOTROPICLINEARELASTIC = 716,
    VISCOLINEARELASTIC = 717,
    VISCOSTVENANTKIRCHHOFF = 718,
    NEOHOOKEAN = 719,
    VISCONEOHOOKEAN = 720,
    ISOTROPICLINEARELASTICJ2PLASTIC = 721,
    ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 722,
    HYPERELASTIC = 723,
    MOONEYRIVLIN = 724,
    VISCOMOONEYRIVLIN = 725,
    HENCKY = 726,
    OGDEN = 727,
    SIMOELASTIC = 728,
    SIMOPLASTIC = 729,
    LOGSTRAINPLASTIC = 730,
    SVKPLSTRESS = 731,
    VISCOLINPLSTRESS = 732,
    VISCOSVKPLSTRESS = 733,
    VISCOFABRICMAP = 734,
    SHELLVISCOFABRICMAP = 735,
    VISCOFABRICMAT = 736,
    SHELLVISCOFABRICMAT = 737,
    PARTITIONGAP = 738,
    PLANESTRESSLINEAR = 739,
    PLANESTRESSSTVENANTKIRCHHOFF = 740,
    PLANESTRESSNEOHOOKEAN = 741,
    PLANESTRESSMOONEYRIVLIN = 742,
    PLANESTRESSBILINEARPLASTIC = 743,
    PLANESTRESSFINITESTRAINPLASTIC = 744,
    PLANESTRESSVISCOLINEARELASTIC = 745,
    PLANESTRESSVISCOSTVENANTKIRCHHOFF = 746,
    PLANESTRESSVISCONEOHOOKEAN = 747,
    PLANESTRESSVISCOMOONEYRIVLIN = 748,
    SURFACETOPOLOGY = 749,
    MORTARTIED = 750,
    MORTARSCALING = 751,
    MORTARINTEGRATIONRULE = 752,
    SEARCHTOL = 753,
    STDMORTAR = 754,
    DUALMORTAR = 755,
    WETINTERFACE = 756,
    NSUBS = 757,
    EXITAFTERDEC = 758,
    SKIP = 759,
    ROBCSOLVE = 760,
    RANDOMSAMPLE = 761,
    OUTPUTMEMORY = 762,
    OUTPUTWEIGHT = 763,
    SOLVER = 764,
    SPNNLSSOLVERTYPE = 765,
    MAXSIZE = 766,
    CLUSTERSOLVER = 767,
    CLUSTERSOLVERTYPE = 768,
    WEIGHTLIST = 769,
    GMRESRESIDUAL = 770,
    SLOSH = 771,
    SLGRAV = 772,
    SLZEM = 773,
    SLZEMFILTER = 774,
    PDIR = 775,
    HEFSB = 776,
    HEFRS = 777,
    HEINTERFACE = 778,
    SNAPFI = 779,
    VELSNAPFI = 780,
    ACCSNAPFI = 781,
    DSVSNAPFI = 782,
    MUVSNAPFI = 783,
    PODROB = 784,
    ROMENERGY = 785,
    TRNVCT = 786,
    OFFSET = 787,
    ORTHOG = 788,
    SVDTOKEN = 789,
    CONVERSIONTOKEN = 790,
    CONVFI = 791,
    ROMRES = 792,
    SAMPLING = 793,
    SNAPSHOTPROJECT = 794,
    PODSIZEMAX = 795,
    REFSUBTRACT = 796,
    TOLER = 797,
    NORMALIZETOKEN = 798,
    FNUMBER = 799,
    SNAPWEIGHT = 800,
    ROBFI = 801,
    STAVCT = 802,
    VELVCT = 803,
    ACCVCT = 804,
    CONWEPCFG = 805,
    SCALEPOSCOORDS = 806,
    NODEPOSCOORDS = 807,
    MESHSCALEFACTOR = 808,
    PSEUDOGNAT = 809,
    PSEUDOGNATELEM = 810,
    USENMF = 811,
    USENMFC = 812,
    USEGREEDY = 813,
    USEPQN = 814,
    FILTERROWS = 815,
    VECTORNORM = 816,
    REBUILDFORCE = 817,
    REBUILDCONSTRAINT = 818,
    SAMPNODESLOT = 819,
    REDUCEDSTIFFNESS = 820,
    UDEIMBASIS = 821,
    FORCEROB = 822,
    CONSTRAINTROB = 823,
    DEIMINDICES = 824,
    UDEIMINDICES = 825,
    SVDFORCESNAP = 826,
    SVDCONSTRAINTSNAP = 827,
    USEMASSNORMALIZEDBASIS = 828,
    USECONSTANTMASS = 829,
    ONLINEMASSNORMALIZEBASIS = 830,
    STACKED = 831,
    NUMTHICKNESSGROUP = 832,
    STRESSNODELIST = 833,
    DISPNODELIST = 834,
    RELAXATIONSEN = 835,
    QRFACTORIZATION = 836,
    QMATRIX = 837,
    RMATRIX = 838,
    XMATRIX = 839,
    EIGENVALUE = 840,
    NPMAX = 841,
    BSSPLH = 842,
    PGSPLH = 843,
    LIB = 844
  };
#endif

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED

union YYSTYPE
{
#line 47 "p.y" /* yacc.c:355  */

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

#line 781 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:355  */
};

typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_HOME_ANARKHEDE_TINKERCLIFFS_FEMTESTING_PARSER_D_PARSER_HPP_INCLUDED  */

/* Copy the second part of user declarations.  */

#line 798 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:358  */

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
#define YYLAST   8491

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  590
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  260
/* YYNRULES -- Number of rules.  */
#define YYNRULES  1426
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  3782

/* YYTRANSLATE[YYX] -- Symbol number corresponding to YYX as returned
   by yylex, with out-of-bounds checking.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   844

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
     585,   586,   587,   588,   589
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   195,   195,   202,   203,   206,   207,   209,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   233,   234,   236,
     237,   238,   239,   240,   241,   243,   244,   245,   246,   247,
     248,   251,   252,   253,   254,   255,   256,   257,   258,   259,
     260,   261,   262,   263,   264,   265,   267,   269,   270,   271,
     272,   273,   274,   275,   276,   277,   278,   279,   280,   281,
     282,   283,   284,   285,   286,   287,   288,   289,   290,   291,
     292,   293,   294,   295,   296,   297,   298,   299,   300,   301,
     303,   305,   307,   309,   311,   313,   314,   315,   316,   317,
     318,   319,   320,   321,   322,   323,   324,   325,   326,   327,
     328,   329,   331,   332,   334,   335,   336,   337,   339,   341,
     342,   343,   344,   345,   346,   347,   349,   350,   351,   352,
     353,   355,   356,   357,   358,   360,   362,   364,   366,   368,
     370,   372,   374,   375,   376,   377,   378,   379,   380,   381,
     382,   383,   384,   387,   394,   400,   401,   406,   418,   424,
     425,   429,   431,   433,   435,   439,   443,   447,   469,   493,
     526,   527,   528,   529,   530,   533,   545,   549,   551,   556,
     561,   565,   607,   659,   660,   662,   670,   672,   680,   682,
     684,   686,   688,   692,   700,   702,   704,   706,   708,   710,
     712,   714,   716,   718,   720,   722,   724,   728,   730,   734,
     736,   738,   740,   742,   745,   747,   749,   753,   755,   757,
     761,   763,   765,   767,   769,   771,   773,   775,   779,   780,
     781,   782,   783,   784,   787,   788,   792,   794,   803,   805,
     816,   818,   822,   824,   828,   830,   834,   836,   840,   842,
     846,   853,   860,   869,   873,   874,   878,   884,   889,   894,
     900,   901,   903,   905,   910,   911,   915,   916,   920,   923,
     928,   933,   937,   942,   948,   955,   958,   960,   962,   964,
     972,   974,   976,   978,   980,   982,   984,   986,   988,   990,
     992,   996,   998,  1000,  1002,  1005,  1007,  1009,  1011,  1013,
    1015,  1017,  1019,  1022,  1024,  1026,  1028,  1030,  1032,  1034,
    1036,  1038,  1040,  1042,  1044,  1046,  1048,  1050,  1052,  1054,
    1056,  1058,  1062,  1063,  1066,  1069,  1071,  1073,  1075,  1077,
    1079,  1081,  1083,  1085,  1087,  1089,  1092,  1096,  1099,  1102,
    1104,  1106,  1109,  1111,  1113,  1117,  1121,  1127,  1129,  1132,
    1134,  1138,  1140,  1144,  1146,  1151,  1154,  1158,  1162,  1164,
    1166,  1169,  1171,  1172,  1173,  1174,  1176,  1178,  1180,  1182,
    1189,  1191,  1193,  1195,  1197,  1201,  1206,  1212,  1214,  1218,
    1222,  1227,  1242,  1244,  1245,  1247,  1250,  1252,  1254,  1256,
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
    4829,  4840,  4851,  4856,  4861,  4866,  4876,  4887,  4898,  4903,
    4908,  4913,  4919,  4925,  4931,  4936,  4941,  4946,  4952,  4958,
    4964,  4969,  4974,  4979,  4985,  4991,  4997,  5002,  5007,  5012,
    5017,  5022,  5027,  5032,  5037,  5042,  5047,  5052,  5057,  5062,
    5067,  5073,  5079,  5084,  5089,  5095,  5101,  5106,  5112,  5117,
    5123,  5128,  5133,  5138,  5148,  5159,  5170,  5175,  5180,  5185,
    5195,  5206,  5217,  5223,  5232,  5241,  5251,  5261,  5270,  5279,
    5288,  5297,  5306,  5315,  5324,  5333,  5342,  5351,  5360,  5369,
    5378,  5387,  5397,  5407,  5416,  5425,  5435,  5445,  5454,  5463,
    5473,  5483,  5488,  5494,  5503,  5513,  5522,  5532,  5537,  5543,
    5552,  5562,  5571,  5581,  5588,  5595,  5602,  5609,  5616,  5623,
    5630,  5637,  5644,  5651,  5658,  5665,  5672,  5679,  5686,  5693,
    5700,  5705,  5710,  5715,  5721,  5727,  5736,  5745,  5751,  5757,
    5766,  5775,  5780,  5785,  5790,  5795,  5800,  5808,  5817,  5836,
    5840,  5846,  5847,  5849,  5856,  5857,  5868,  5869,  5880,  5884,
    5887,  5894,  5898,  5902,  5913,  5915,  5917,  5919,  5921,  5923,
    5925,  5927,  5929,  5932,  5934,  5936,  5938,  5941,  5943,  5945,
    5947,  5960,  5966,  5973,  5977,  5984,  5988,  5990,  5992,  5994,
    5997,  5999,  6002,  6006,  6008,  6010,  6013,  6017,  6021,  6025,
    6029,  6033,  6037,  6038,  6039,  6040,  6044,  6048,  6052,  6054,
    6057,  6060,  6062,  6065,  6069,  6071,  6073,  6076,  6078,  6080,
    6083,  6087,  6089,  6091,  6093,  6095,  6097,  6099,  6101,  6104,
    6108,  6110,  6112,  6114,  6116,  6118,  6120,  6122,  6124,  6126,
    6128,  6131,  6135,  6138,  6142,  6147,  6149,  6151,  6153,  6155,
    6157,  6159,  6161,  6163,  6166,  6169,  6171,  6178,  6185,  6186,
    6191,  6200,  6202,  6206,  6211,  6213,  6219,  6223,  6227,  6230,
    6234,  6239,  6241,  6246,  6248,  6250,  6252
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
  "FINITESTRAINPLASTIC", "LINEARELASTIC", "STVENANTKIRCHHOFF",
  "TULERBUTCHER", "LINPLSTRESS", "READ", "OPTCTV",
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
     835,   836,   837,   838,   839,   840,   841,   842,   843,   844
};
# endif

#define YYPACT_NINF -2823

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-2823)))

#define YYTABLE_NINF -1422

#define yytable_value_is_error(Yytable_value) \
  0

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    7910,  -163,    59,  -149,   857,   740,  2382,  2382,  2382,   888,
    -138,   884,   857,  -116,  -114,   857,  -106,   948,   -76,   -57,
     -27,   -21,    69,   100,   205,   265,   300,   368,   384,   445,
     477,   490,   992,   999,   518,   546,   555,   577,   582,   593,
     857,   857,  1054,  1064,   607,  -212,   669,   690,   725,   975,
     747,   766,  1085,   770,  1088,   782,   784,   839,   842,   853,
     872,   891,   894,   896,  -200,   984,   899,   937,   857,   941,
    1094,   942,   950,   857,   955,   968,  1000,   118,  1137,   972,
     977,  1160,   -26,  -183,   983,  1279,  1298,   857,   991,  1003,
     993,   857,   998,  1009,   846,   979,  1014,  1015,   857,   857,
    1022,  1025,  1308,   394,   857,   857,  1027,  1321,  1397,   333,
    1029,  1034,  1042,  1045,  1062,  1066,  1068,   857,  2382,  1070,
    1411,  1071,  1075,  2382,  2382,  1082,  1084,  1087,  1090,  1098,
    1101,  1102,  1103,   500,  1109,  1113,  1117,  1118,  1120,  1122,
    1124,  1129,  1130,   -16,  1419,  1434,  1135,  1140,   857,   857,
     857,  1142,  1143,  2382,  1144,  1145,  1146,  1154,  1203,  1204,
    1205,  1213,  1215,  1219,  2382,  1222,  1231,  1234,  1236,  1239,
     613,  1490,  7341, -2823, -2823, -2823,  1340,   857,   361,    40,
   -2823,     1,   857,   106,  2382,  2382,   857,   736,   857,   857,
     585,  2382,  1378, -2823, -2823, -2823, -2823, -2823,   339,    20,
    1280, -2823,  2382, -2823, -2823, -2823, -2823,    28, -2823,   400,
     -12, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
    1302, -2823, -2823, -2823, -2823, -2823,   857, -2823,   374,   857,
   -2823, -2823,   387,   392,   489,   857, -2823,   857, -2823,   857,
     857,   857,   857, -2823, -2823,   857,   857,   857, -2823, -2823,
      54,  1287,   857,   857,   857,  1290,  1293, -2823,   530, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
    2382, -2823, -2823,   857,   857,   857, -2823, -2823,   857,   857,
     857,   857,  -130,   175,   857,   857,   857,   857,   857,   857,
     376,    13,  2382,  1772,   857,   533, -2823, -2823,   857,   270,
     699, -2823, -2823, -2823,    41,  2382, -2823, -2823, -2823,    51,
   -2823, -2823,   857,  -140,   857, -2823,   287,  5137,  5137,  1815,
    5137, -2823,   446,   857,  1258,  1501,  1504, -2823, -2823,  1269,
   -2823,  1273, -2823, -2823, -2823, -2823,  1275,  1276,  2382,  1444,
   -2823,  2382,   857,   857,  1277,  1278, -2823, -2823,  1013, -2823,
   -2823,  1292, -2823,  2382,  1424, -2823,   857, -2823, -2823,  2382,
    2382, -2823, -2823, -2823, -2823, -2823,  2382,  2382,  1435,  2382,
    1206,   664, -2823,  1295, -2823,  1304, -2823,   857,   857,   857,
   -2823,  2382,  1462,  1305, -2823,  1313,  1356,  1314, -2823,  1337,
   -2823, -2823, -2823,  2382,  2382, -2823,   857,   857,  1339,   857,
   -2823,  1341, -2823,   857,  2382,  2382,   857,   857, -2823, -2823,
     857,   857,  1345, -2823,  1347,   857,   857,  1348,  2382,   857,
    1353,  2382,   857,  1354,  2382, -2823,  1079,  2382,  1359, -2823,
      22, -2823, -2823, -2823,  1368,  1385, -2823,  1386,   857, -2823,
    1388, -2823,  1389,  1396, -2823,   857,  2382,   857,  1401, -2823,
   -2823,  1405,  1470, -2823,  1416,  1422,  1472, -2823,  1423, -2823,
     857,  1425,  1426,   857, -2823, -2823,  1431,  -213,  1432,  1437,
   -2823, -2823,  1440, -2823,  1503,  1622,   857,  1282, -2823, -2823,
   -2823,  1309,  1632,   857,  1443,  1446, -2823, -2823,  2382,   857,
   -2823,  1454,  1460, -2823,   857,  2382, -2823, -2823,  1651, -2823,
   -2823,  1469, -2823,   857,  1664,  1665,  1366,  1666,  1667,  1672,
   -2823, -2823,   857, -2823,  1488, -2823,  1496, -2823, -2823,   351,
     857,  1675, -2823, -2823,  1499, -2823, -2823,   857,   857,   857,
   -2823, -2823, -2823, -2823, -2823,  2382, -2823, -2823, -2823, -2823,
   -2823,  1505, -2823, -2823, -2823,   538,  1427,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  1037,  1436,  2382, -2823, -2823, -2823,
    2382, -2823,  1381,  1390,  1639,  1652,  1661,  1677,  1680,  1527,
    1560,  1687,  1565,  1688,  1566,   857,  1579,  1580,  1582,  1583,
    2382,   857,   857,   857,   857,  2382,  2382,  1536,   857,   857,
     857,   857,   857, -2823,   857,   857,   857,   857, -2823,    49,
   -2823,  2382,  1591,   857,  2382,   241,   857, -2823,  -120,   671,
     751,  1717,  1722,  1729,  1730,   459,   857,  1702,   857,  2382,
    1041,  2382,  2382,  1485,  1806,  2382,  1610,  1612,  1616,  1502,
     309,  1507,  2382,  1513,   624, -2823, -2823, -2823,  2382,  1881,
    2382, -2823, -2823, -2823,  1619,   857,  1762,  1941,   857,  1958,
   -2823,  2382,   857, -2823,   -89,  1818,   857,   -39,   857,  2382,
     857,  2382,  2382, -2823,   857, -2823,   857, -2823,   857, -2823,
     857, -2823,   857, -2823, -2823, -2823,  1766, -2823,   -23,  2069,
    2382,  2382,  2382,  2091,  1646,   857, -2823,   -60,  1648, -2823,
      25, -2823,   857,   857,   857,   857, -2823,   857,   857, -2823,
      -5,   857,  1650,  1660,   857,  2382,  2382,  2382,  2382,  2382,
    1662,   857,  1663,  1671,   857,  1674,  1684,  1690,  1694,  2382,
    1698,  1699,   857,  1703,  2382,  1706,  2382,   857, -2823,  2382,
   -2823, -2823, -2823,  2180,   857,   -87,  1708,  1709,  2382,   857,
    1776, -2823, -2823,  1710,  1786,   857,  1719,   857,   857,  1988,
     -31,  2382,  2382,  1533,  1791,  2382,  2382,   857,   857,  1723,
    2382,  1724,  2382,   857,  2099,  1549,   -58,  2114,  2136,   857,
    1046,   857,   704,   857, -2823,  6378, -2823, -2823, -2823,  2382,
    1725,  2382,   857,   857,  1086,  2382,   104,   857,  1733,  1735,
    2382,   857,   857,  1739,  2155, -2823,   857,  1904,  4719,   857,
     857,   857,   857,  1668,   857,   857,  1670,  1478, -2823, -2823,
   -2823, -2823, -2823,  2382,   857,  1907,   -28,  1681, -2823,  1784,
     857,   857, -2823,   857,  1787,  2382,  1919,  1712,  2382,   857,
    1716,  1718,  1720,  1726,  1732,  1738,  1744,  1745,   857,   857,
    1542,  2382,  1931,  1936,   857,   857,  2382,  1751,  1756,   857,
    1770,  1945,  1779,  1781,  1946,  1950,  1792,  1793,  1794,  1795,
    1801,   857,   857,   857,  1829,  1834,  1838,  1842,  1845,  2382,
    1748,   363,  1847,  1981,   857,  1858,   857, -2823,    19, -2823,
    1091,  2382, -2823, -2823, -2823, -2823,  2382, -2823,  1860,  2246,
   -2823,   857,  2382,   857,   857, -2823, -2823,  1096, -2823,   857,
    1862,  1889,  1892,  2382,  2257, -2823,  2382,  2382, -2823,   679,
     857, -2823, -2823, -2823, -2823, -2823,   857, -2823,  2382, -2823,
    1895, -2823,  1914,  2382, -2823,  2297,  1964, -2823,  1099,  2382,
     857, -2823,   857,   857,   857,   857, -2823,   857, -2823, -2823,
   -2823,  1948, -2823,  1952, -2823, -2823,   857,   857,   881,   857,
   -2823, -2823,   857,   857,  2382,  2382,  2382, -2823,  2382,  1955,
   -2823,  2382,   857,   857,   857,  1100,  2382, -2823,  1822, -2823,
    1827, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,  2382,
     857, -2823, -2823, -2823,  1956, -2823, -2823, -2823,  1959, -2823,
     857, -2823, -2823,   857, -2823, -2823,  1973,  2382,  2382, -2823,
   -2823,  1975, -2823,  2005, -2823,  2008,   857,    91,   857,  1866,
     857,    63, -2823,  2382,   857, -2823,   857, -2823, -2823, -2823,
    2010,   857,  1987, -2823,   857,   857,  1989,   857,  1990,   857,
     114,   857,  2017,   857,  2688,   857,  2778,  2016, -2823, -2823,
   -2823,  2028,  2382,  -171, -2823,  2032, -2823,   857, -2823,  1119,
   -2823,   857, -2823,   857,  2382, -2823,   857,   857,  2382,  2382,
    1149,  2382,  2382,  2382,  2382,  2034, -2823,  2382,   857,  2035,
     857,   405,   481,  2036,  2038,  2039,  2044,  2048, -2823, -2823,
    2055, -2823, -2823,  2056, -2823,  2057, -2823, -2823, -2823, -2823,
    2058,  2061,  2062,  2063,  2065,  2067,  2071,   857,   857,  2078,
     415,  2082,  2084,  2088,  2089, -2823,  2382, -2823, -2823,   857,
   -2823,   857,  2382,  2382,  2214,  1772,  2382,   504,  2092, -2823,
     857,   857,   878,   857,   857,   857,   857, -2823, -2823, -2823,
   -2823, -2823, -2823,   857, -2823,   857, -2823,   857, -2823, -2823,
    1841, -2823, -2823, -2823,  2382,  2093, -2823, -2823, -2823,  2096,
    2382, -2823,  2382,  2382,  2100, -2823,  2101, -2823, -2823,  2102,
     857, -2823, -2823,  2103,  2104,  1209,  2306,  2382,  2105,   857,
   -2823, -2823,  2382, -2823,  2110,  2382, -2823,  2119,  2123, -2823,
   -2823, -2823, -2823, -2823,  2382, -2823,   857,   857,  2382,   857,
    2130,  2382,  2132,  2382,  2382,  2382,   857,   857,   857,   857,
     857,   857,  2268,  2273,   857,  2156,  2382,  2382,  2382,   857,
    2158,   857, -2823,   857, -2823,  2382,  2382,     6,   857,  2382,
    2382,  2382,  2382,   857,   857,   857,   857,   857,   857, -2823,
     857,  -118,   857, -2823, -2823,  2164,  2188,  2190,  2202,  2204,
    2207, -2823,  2208, -2823, -2823,  2210, -2823, -2823, -2823, -2823,
    2215, -2823, -2823,  2217, -2823,  2218, -2823,  2219,  2222,  1268,
    2382,  2382,  2382,  2382,   215, -2823, -2823,  2223,   857,  2225,
   -2823,  2226, -2823,  2334,   -22,   857,  2364,  -152,   857,  1270,
    2233,  2235,  2236,  2237,  2239,  2243,  2338, -2823,  2361, -2823,
     857,  2244, -2823,  2248,  2410, -2823, -2823,  2249,  2254,  2255,
   -2823,  2258,  2267,   885, -2823,  2382, -2823,  1284,  2279, -2823,
    2382,  2287, -2823,  2290,  2291,  2302,  1285,    11,    63,  2150,
      63,    57,  2382,    63,  2182,  2195,  2311,   289,  2321,  2324,
     857,  2382,  2382,  2382,  2382,  2330,  1038,    63,   857,   857,
     857,   744,  1078,   413,  2332,   857,   857,   857,   857,   857,
    2411,  2339,   857,   806,   537,   857,   857,   269,   857,  2382,
    2382,  2382,   857,  2227,   625,  2382,   857,   857,   857,  2382,
     857,   857,   857,   857,   857,   857,  2340,   857,  2382,   857,
    2382,  2382,  2382,  2382,   857,  2412,   857,  2341,  2417,  2382,
     -83,  2382, -2823,  2382,  2382,  2414, -2823,  2365,  2372, -2823,
    2441, -2823,  2374, -2823, -2823,  1288,  2506,  2376, -2823, -2823,
    2377,  2382,  2502,  2382,  2382, -2823,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  2382,  2565, -2823,  2535,  2540, -2823,
   -2823,   857, -2823,  2382,  2514,  2514,  2514,  2514,  2514, -2823,
   -2823, -2823, -2823, -2823, -2823,  2277,  2514, -2823,  1772,  2382,
    2382, -2823,  2396,   857, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,   857,
    2527, -2823,   857, -2823, -2823, -2823,  1772, -2823,  2539, -2823,
   -2823,   857,   857, -2823, -2823, -2823, -2823, -2823, -2823,   857,
     857, -2823, -2823, -2823, -2823, -2823,  2418,  2382,   646, -2823,
     857, -2823, -2823, -2823,    33,   857, -2823,  1289,  2382,  2419,
   -2823, -2823,  2421, -2823,  2422,   857, -2823,  2423,  2425,   857,
   -2823, -2823,  2382, -2823,  2382, -2823, -2823,  2382, -2823,  2426,
   -2823, -2823,  2382, -2823,  2382,   857,  2429,  2571, -2823,  2585,
    2434, -2823,  2382,   857, -2823,  1772, -2823, -2823,   857, -2823,
     828, -2823,   857,  2382, -2823,  2442,  2382, -2823,  2382,  1294,
    2382,  2382, -2823,  2382,  2466,   857, -2823,   892, -2823,   857,
   -2823,    26,  2473,   194,  2491,  2494,  2508, -2823, -2823,  2606,
     857, -2823,  2382,  2382, -2823, -2823,  2343,  2509,   857,  2382,
    2510,   857,  2382,  6378,  2518, -2823,  1772, -2823,  2522,   857,
    2382,  2524,   857,  2382,  2528,   857,  2382,   -55,   857, -2823,
    2530,   857,  2382,  2532,   857,  2382,  2534,   857,  2382, -2823,
   -2823,   232, -2823,  2382,  -219,  -201, -2823, -2823, -2823,  2536,
   -2823,   857,  2542,   857,  2621,  2382,  2544, -2823,   857,   857,
     857,   857,   857, -2823,  2545,  2627, -2823,  2546,  2547, -2823,
    2548, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823,  2549,  2550, -2823,
     857,  2556,  2382, -2823, -2823, -2823, -2823,  2382,  2382,  2564,
    2382,  2382,  2573,  2382,  2589,  2590,  2592, -2823, -2823,   321,
    2593, -2823,   574,   722,  2664,   805,  2718, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823,  2382,  2597,   643, -2823, -2823, -2823,
     857, -2823,  2382, -2823, -2823, -2823,  2382, -2823,  2382,  2382,
   -2823, -2823,  2382, -2823,   857, -2823, -2823,  2382,   380,   857,
    2604,   202, -2823,  2608, -2823,  2382,  2382,  2382,  2630, -2823,
    2631,  2642,  2712,  2737,  2504,   857,   857,   857,   857,  2382,
    2382,  2382,   857,   857,   857,   554,  2382,  1374,  2382,  2382,
   -2823,  2382,   857,   161,  2382,  2382,  2382,   534,  2382,  2382,
    2745, -2823,  2773,   609,  1379,  2797,  1380,  2618, -2823,   857,
     857, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823,  1383,  2382,  2382,  2382,   145,
    2623, -2823,   404,   219, -2823,  2753, -2823, -2823, -2823,  2628,
     857,   310,  2382,   857,   317, -2823,  2633,  2382, -2823,  2382,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
     857,   857, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823,  2634, -2823,  2799,  2803, -2823,  2806, -2823,  2635, -2823,
   -2823, -2823, -2823, -2823,  1400,   313,  2638,  -139, -2823,  2639,
   -2823,  2640,  2807,  2643, -2823,  2648,  2649, -2823,  2651,  2659,
   -2823, -2823,  2667,  2704,  2708,  2711,  2713, -2823,  2717,  2720,
   -2823,  2382,  2721,  2723,   901,  1065,  2724,  2725,  2734,  2736,
   -2823,  2742,   857,  2382,  2743,  2744, -2823,  2749, -2823,  2754,
    2755,  2757,  2758,  2759, -2823,  2760,  2761,  2762,  2766,  2770,
    2772,  2382,  2382,   857,  2775,  2776,  2791,  2793,  2795,  2796,
    2800,  2808,  2816,  2818,  2821,  2822,  2832,  2835, -2823,  2838,
    1408,  2850,  1417,  2866,  2869,  2873,  2875, -2823,  2876,  2880,
   -2823,  -103,  1418, -2823,  2882,  2883,  2824,  2382, -2823,  2884,
   -2823, -2823, -2823,  1428, -2823, -2823,  1433, -2823,  2842, -2823,
   -2823,  2382,  2887,   857,   857,  1438,   857,   857,  2382,  2382,
     -30,    44,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,    53,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,   857,   857,   857,   857,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  2382, -2823,  2891, -2823, -2823, -2823,
    2632, -2823, -2823,   857, -2823,  2382,   857, -2823, -2823, -2823,
    3016,   857,   857, -2823,  3019,   857,   857, -2823, -2823, -2823,
    2382, -2823, -2823,   857,   556, -2823,  1459,  2896, -2823, -2823,
   -2823, -2823, -2823,  2899,  2382,  2382, -2823, -2823, -2823,  2382,
     857,   857,   857,  2903, -2823,  2904,  2382,  1463,  1772,  2906,
   -2823,  2382,  2382, -2823, -2823,  2382, -2823,  2907,  2908,  2909,
    2915, -2823, -2823,  2382, -2823,   857,   721, -2823, -2823,  2798,
   -2823,   857, -2823,  3499, -2823,  2916,   857,  2920,  2382,  2922,
    2382,  2382,  2923,  2927, -2823,  1772,  2929,  2382,  2933,  2936,
    2382,  2938,  2946,  2382,  2948,  2959,  2843, -2823,   -32,  1483,
    2382,  2961,  2965,  2382,  2976,  2978,  2382,  2980,  2990, -2823,
    2991,  1492, -2823,  2382, -2823,  2382, -2823,  2858, -2823,  2997,
   -2823,  3000,  3005, -2823,   857,  3009,  3013,  2937,  2940, -2823,
   -2823,   857, -2823, -2823, -2823, -2823, -2823,   857,   857,   857,
    2382,  2382,  2382, -2823,  2382,  2382, -2823,  3024, -2823, -2823,
   -2823,   857,  3026, -2823,  3027, -2823,   857, -2823,   857,   857,
   -2823,   857,  2382, -2823,  3035, -2823, -2823,  3037,   857,  1546,
    3039,   857,  1557,  3043,   857,  2382,  3054, -2823,   857,  3058,
   -2823,  2382,  3062,  3064, -2823, -2823, -2823, -2823, -2823, -2823,
    3069,  1585,  1586,   857,  1603,  2382,  2382,   857,   857,  2382,
    1618, -2823,  2382,  2382,  1623,   905,  3080,  2382,  2382, -2823,
    2382,   857,  2382,  2382,  2382,  2382, -2823,  2382,  2382, -2823,
   -2823, -2823, -2823,   731,   765, -2823,  2382, -2823, -2823,  2382,
   -2823,  3086,  2382, -2823,  2872,  2382,  2382,  2382,  3093, -2823,
   -2823,  2382,   204, -2823,  2382,   288,  3098, -2823,  2382, -2823,
    3104,   397,  2382, -2823,  3108, -2823,  1637,  3109,   857,  2949,
   -2823, -2823,  3113, -2823,  1644, -2823,  2993, -2823, -2823,  3114,
   -2823,  3118, -2823,   857,  -111, -2823, -2823, -2823,  3133, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823,  3138, -2823, -2823,   716,   857, -2823,   857, -2823,   745,
   -2823, -2823, -2823, -2823, -2823,  3142,  3144, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823,  3147,  3150,  3152, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
    3153, -2823, -2823,  3169, -2823, -2823, -2823, -2823, -2823, -2823,
    3171, -2823, -2823,  3173, -2823, -2823, -2823,   857,  2994, -2823,
   -2823,  1657, -2823,  2382, -2823,  3175,  2382, -2823,   857,   857,
   -2823, -2823,   857,   857,  2382,  2382, -2823,  2382,  2382, -2823,
    2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382, -2823,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382,   857,   857,   857,   857,  2382,  2382,  2382,
    2382,  2382,  2382,  2382,  2382,  2382,  2382, -2823, -2823,   857,
    2382,  2382, -2823, -2823,   857, -2823, -2823, -2823,  3177,   566,
     857,  2382, -2823,  3183, -2823,  3143,  2382, -2823,  3184,   857,
     857,   857, -2823,  1658, -2823,  3185,  2382,  3188, -2823,  1682,
    3198,  3203, -2823, -2823, -2823, -2823,  3207,   548, -2823,  3219,
   -2823, -2823,   857, -2823, -2823, -2823,  2382, -2823,  2382, -2823,
    2916, -2823,  2916,  3149,  2382,  2382,  2382,  2382,  2382, -2823,
    2382,  3221, -2823,  2382,  2382, -2823,  2382,  2382, -2823,  2382,
    2382, -2823,  3224,  1689,  3105, -2823, -2823,  2382,  2382, -2823,
    2382,  2382, -2823,  2382,  2382, -2823, -2823, -2823,  3230,  1695,
    1696, -2823, -2823, -2823, -2823,  3231, -2823, -2823,  1711,  2382,
    3234,   857,  2382,   857,  2382,  2382,  3235,  2382,  2382, -2823,
     857, -2823, -2823, -2823,   927, -2823,   961, -2823, -2823, -2823,
     857, -2823,  1715, -2823,  3237, -2823,  3245, -2823,   857,  3247,
   -2823,  2382, -2823,  3249, -2823, -2823, -2823, -2823,  3258, -2823,
    3263, -2823,  2382,  2382,  2382,   857,   877, -2823,  1721,  2382,
    2382, -2823,  3266,  2382, -2823,   930, -2823,  2382,  1734,   958,
    3268,  2382,  2382,  2382,  3269,  2382,  2382, -2823, -2823,   749,
     844,  3270,  3277, -2823,  2382, -2823,  3008,  2382,  2382,  3279,
   -2823,  3280,  3282, -2823,  -180, -2823,  2382,  2382,  3285, -2823,
    3294, -2823, -2823,  3295,   485, -2823, -2823,  2382, -2823,  3036,
   -2823,  3100, -2823, -2823,  3102, -2823,  3115, -2823, -2823,   -80,
   -2823,  3310, -2823, -2823,   857, -2823,  3311,  3317,   857, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,  3321,
   -2823,   857, -2823,  3117,  3322, -2823,     9,   857,  2382,   857,
    2382,  2382,  2382,  2382,    60,  2382,    64,  2382,  1737,  3323,
    2382,  2382,  -184,  2382,  2382,  2382,  2382,  2382,  2382,    68,
    2382,  3328,  2382,  2382,  2382,  2382,  2382,   857,   857,  2382,
    2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,  2382,  2382, -2823,  2382, -2823,   857,  2382,   857,  3331,
   -2823,  3338, -2823, -2823,   857, -2823,  3134, -2823,  3340, -2823,
   -2823,  3342, -2823, -2823, -2823, -2823,   763, -2823, -2823, -2823,
   -2823, -2823, -2823,  3343,  2382,  2382,  2382,  3346,  2382,  3353,
   -2823,  3354,  2382,  3361,  2382,  3368,  2382, -2823,  2382,  3370,
    1741,  3374,  2382,  3376,  2382,  3378,  2382, -2823, -2823,  3383,
   -2823,  3390, -2823, -2823,   857,   857, -2823,  2382,  3391,  2382,
    2382, -2823,  2382,  1768,  3398,   857, -2823,   857, -2823,  3399,
   -2823,  3400, -2823, -2823,  2382, -2823,  3409, -2823, -2823, -2823,
    2382,  2382,  2382,  2382,  3410, -2823, -2823,  2382,  2382,  2382,
   -2823,  3413,  2382,  2382, -2823,   587,  2382, -2823,  3429,  2382,
   -2823,   971, -2823,  2382,  2382,  2382, -2823,  2382,  2382, -2823,
   -2823,   797,   640, -2823, -2823,   857, -2823,  3432,  2382,  2382,
   -2823, -2823, -2823, -2823,  2382,  3439,   244, -2823, -2823, -2823,
   -2823,  3440,  3442, -2823,  3160, -2823,  3200, -2823,  3450, -2823,
    3451, -2823,  3459, -2823,  3463, -2823, -2823,  3469, -2823,  3474,
   -2823,  3477, -2823, -2823,   857,  2382,   833,  2382,  2382,  2382,
    2382,  2382,  2382,  2382, -2823,  2382,  2382,  2382, -2823,  2382,
    2382,  1769, -2823,  1782, -2823,  2382,  2382, -2823,  2382,  2382,
    2382,  2382,  -177,  2382,  2382, -2823,  2382,  2382,  1788, -2823,
    2382,  2382,  1864,  2382,  2382,  2382,  2382,  2382,  2382,   101,
     119,  -127,  2382,  2382,  2382,  2382,  2382,  2382,  2382, -2823,
   -2823, -2823,   857,  3479,  2382, -2823,  3605, -2823,  3483, -2823,
   -2823, -2823, -2823,  2382,  2382,  3484, -2823,  3485, -2823, -2823,
    3486, -2823,  3496, -2823,  3498,  1939, -2823, -2823,  2382, -2823,
    3501, -2823,  3507, -2823,  3509, -2823, -2823,   857,   857,  3515,
   -2823,  2382,  2382,  2382, -2823,  2382, -2823, -2823, -2823, -2823,
   -2823,  3520, -2823,  3521,  3527,  3528,   921, -2823,  2382,  2382,
     480, -2823,  2382,  2382,  2382, -2823,  2382, -2823,  3540,  2382,
   -2823,   602,  2382,  2382,  2382,  2382,  2382, -2823, -2823,   848,
    2382,  2382, -2823,  2382,  2382,   273, -2823,  3541, -2823, -2823,
   -2823, -2823,  3222, -2823,  3265, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823,  3542,   148, -2823,  1962,  2382,  1965,  2382,  1967,
    1991,  3543,  2382,  -119,  3544,  2382,  -101, -2823,  2382, -2823,
    1996,  2382,  2382,  2382,  2382,  2382,  2382, -2823,  2382,  2382,
    3545,  2382,  -100, -2823,  2382,  2019,  2042, -2823,  2382,  2382,
    2382,  2382,  2382,  2382,  2382, -2823,  2382,  2382, -2823,  2382,
    2382, -2823,  2382,   -85,  2382,  2382,  2382,  2382,  2382,  2382,
    2382, -2823,  3549,  3550, -2823,  2382,  2382, -2823, -2823, -2823,
   -2823, -2823,  2382,  2050,  2074, -2823, -2823, -2823, -2823,   857,
    2382, -2823,   527,  2382,  2382,  2080, -2823, -2823, -2823, -2823,
    3551, -2823,  2382,  2382, -2823,  2382,   499,  2382,  3552,  2382,
     523, -2823,  2382,  2382, -2823,  2382,  2382,  2382,  2382,  2081,
   -2823,   606,  2382,  2382,  2382,  3553, -2823, -2823, -2823,  3281,
   -2823,  3283, -2823, -2823,   857, -2823,  2382,  2382, -2823,  2382,
    2382, -2823,  2083, -2823,  2115, -2823,  2382, -2823,  2382, -2823,
    2382, -2823,  2382,  3561, -2823,  2141,  2382,  2382,  2382,  2382,
    2203,  2213,  2382,  2382, -2823,  2382, -2823,  2382,  2280, -2823,
    3562, -2823,  2298,  3563,  2382,  2382,  2382,  2382,  2382,  2382,
    2382,   -79,  2382,   -68,  2382, -2823,  2382,  2384,  2387,  2382,
    2382,  2382,  2382,  3564, -2823,  3689,  3569,  2382, -2823, -2823,
    2382,  2393,   857,  2382,   857,  3574,  2382,  2382, -2823,  2382,
   -2823,  2382,  2382,  3575, -2823,  2382,  -182,  2382, -2823,  2382,
   -2823,  2382,   531,  3576,  2382,   857,  2382,  2382,  2382, -2823,
    2382, -2823,   854,  2382,  2382,  2382, -2823, -2823,  3297, -2823,
    3577,   660,  3578,  2413,  3579,  2437, -2823,  2382, -2823,  2382,
    3580,  2382,  3581,  2382, -2823, -2823,  2454,  2382,  2382,  3583,
    2382, -2823,  2457, -2823,  2472,  2382,  2382,  3586,  2382, -2823,
    2488, -2823, -2823,  2382, -2823,  2382,  2382,  2382,  2382,  2382,
    2382,  2382, -2823,  2382,  2382, -2823,  2382,  2382,  2382, -2823,
    2495, -2823,  2496,  2382,  2382,  2382,  2382, -2823,  3588, -2823,
    3589,  2561, -2823, -2823,  3590,  3591,   857,   857,  2382,  3593,
    2567,  2382,  2382, -2823,  3594, -2823,  2382, -2823,  3596,  3597,
   -2823,  2382,   263, -2823,  2382,  3598,  2382,  2382,  2382,  2382,
   -2823,  2382,  2382,  2568, -2823,  3602, -2823, -2823,   865, -2823,
   -2823,  2382, -2823, -2823,  2382,  2570,  2572, -2823,  2382, -2823,
    2382, -2823,  3332,  2382,  2382, -2823,  2382, -2823,  3357, -2823,
    3377,  3603,  2382, -2823,  2382, -2823,  2581,  2663,  2382,  2382,
    2382,  2382,  2382,  2382,  3604,  2382,  3606,  2382,  3607,  2382,
   -2823,  2382, -2823,  2382,  2382,  2382,  2382,  2382, -2823, -2823,
   -2823,  2680, -2823, -2823,   857,  2752, -2823, -2823,  3610,  2782,
    3611, -2823,  3614, -2823, -2823,  3615, -2823,  2382,  3617, -2823,
    2382,  2382,  3620,  2382,  2382,  2783, -2823,  2382, -2823, -2823,
    3623,  3624, -2823,  2784, -2823,  2809,  3625,  3629, -2823,  3392,
     153,   184,   -50, -2823,  3630, -2823,  3632, -2823,  2382,  3633,
   -2823,  2924, -2823,  2942,  2382,  2382,  2382,  2382,  2382,  2382,
   -2823,  2382, -2823,  2382, -2823,  3635,  3067,  3101,  2382,  2382,
    2382,  2382, -2823, -2823,  2382, -2823, -2823,  3187, -2823, -2823,
   -2823,  3636, -2823,   857,  3637, -2823,  2382, -2823, -2823,  2382,
    2382, -2823, -2823, -2823,  3424, -2823,  3426, -2823, -2823, -2823,
     857, -2823,  2382,  2382, -2823,  2382,  2382, -2823,  2382, -2823,
   -2823,    -1, -2823, -2823,  3193, -2823,  3475,  3214,  3257,  2382,
    2382,  2382,  2382,  3638,  3639, -2823, -2823,  3267, -2823,  3272,
     193,   212,    10,  2382,  3329, -2823,  3640, -2823,   857, -2823,
     726,  2382,  3642, -2823,  3643, -2823,  3644,   857,  2382,    67,
    2382,    77,  2382, -2823,  2382, -2823,  3364, -2823,  3645, -2823,
    2382, -2823,  2382,  3369,  3423,  2382,  2382, -2823, -2823, -2823,
    3493, -2823,  3513, -2823,  2382,  2382, -2823,  2382,  2382, -2823,
    2382,    82, -2823,  3646, -2823,   857,  2382, -2823,  2382,  3647,
   -2823, -2823, -2823,   857,  2382, -2823,  2382,  2382, -2823,  2382,
    2382,  2382, -2823,  3457, -2823,  3648,  3649, -2823,  2382, -2823,
    2382,  2382,  2382, -2823,  3650, -2823,  3651,  2382,   110,  2382,
     113,  2382, -2823,  2382, -2823,  3652,  3653,   786, -2823,  2382,
    3654,  2382,  3655,  2382,  3656,  2382, -2823,  3471, -2823, -2823,
    3657,  3660,  3473,  3503, -2823, -2823,  2382, -2823,  2382,  2382,
   -2823,  2382,  2382,  2382, -2823, -2823,  2382, -2823,  2382,  2382,
   -2823,  2382, -2823,  2382, -2823,  3661, -2823,  3517, -2823, -2823,
   -2823,  2382, -2823,  2382,  3665,  2382,  3666,  2382,  3667,  2382,
    2382,  2382,  2382,  3668,  3669, -2823, -2823,  3548,  3670,  3671,
   -2823,  2382, -2823,  2382, -2823,  3672,  3673,  2382,   857, -2823,
   -2823, -2823,  3567, -2823, -2823,  3674,  3675, -2823, -2823,   818,
    2382, -2823,  3626, -2823, -2823,  2382, -2823,  2382,  3676, -2823,
    3627,  2382,   857, -2823, -2823,  3628,  3677,  2382, -2823,  3641,
   -2823,  2382, -2823,  3678,   857, -2823,   857,  2382,  2382,  2382,
    3679, -2823
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
     151,   150,   152,     0,     0,     0,     0,  1421,  1422,     0,
     836,     0,  1424,  1426,  1423,  1425,     0,     0,     0,     0,
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
     183,  1141,     0,   413,     0,   161,     0,  1149,  1301,     0,
       0,     0,   207,   351,     0,   416,   435,     0,     0,     0,
     758,  1311,  1416,  1351,  1356,     0,   869,   864,   866,  1347,
    1349,     0,     1,     2,     4,     0,     0,     0,     0,     0,
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
       0,   860,   861,     0,  1423,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   961,   921,  1091,  1090,  1089,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1133,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1306,  1306,
    1306,  1306,  1306,     0,     0,     0,     0,     0,  1306,     0,
       0,     0,  1337,     0,     0,  1346,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1353,
    1354,  1355,     0,     0,     0,     0,   273,   582,     0,   399,
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
       0,     0,     0,     0,     0,  1304,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1329,  1338,  1344,  1343,
    1336,  1325,  1327,  1340,  1313,  1314,  1315,  1316,  1317,  1330,
    1318,  1322,  1321,  1320,  1319,  1323,  1328,  1408,     0,     0,
       0,  1312,     0,     0,  1376,  1388,  1387,  1374,  1405,  1375,
    1380,  1382,  1381,  1383,  1384,  1365,  1366,  1385,  1386,  1358,
    1361,  1367,  1368,  1364,  1396,  1397,     0,  1395,  1377,  1398,
    1399,  1389,  1392,  1400,  1401,  1371,  1372,  1373,  1402,     0,
       0,  1348,  1350,  1411,  1414,  1352,     0,     0,     0,  1357,
    1418,  1420,  1417,   583,     0,     0,   400,     0,     0,     0,
     336,   337,     0,   526,     0,   529,  1084,     0,     0,     0,
     762,   380,     0,   347,  1073,  1078,  1070,  1080,   591,     0,
    1068,   537,   538,  1097,     0,     0,     0,     0,   419,     0,
       0,   666,     0,   664,   497,     0,  1101,  1103,     0,   656,
       0,   655,     0,     0,   660,     0,  1088,  1092,     0,     0,
       0,     0,   154,     0,     0,     0,   645,     0,   644,     0,
     647,     0,     0,     0,     0,     0,     0,   286,   287,     0,
       0,   346,  1086,  1087,   430,   544,  1308,     0,     0,     0,
       0,     0,     0,   960,     0,   508,     0,   417,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1304,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1142,
     726,     0,   730,     0,     0,     0,   739,   472,   474,     0,
     477,     0,     0,     0,     0,     0,     0,   179,     0,     0,
       0,     0,     0,   176,     0,     0,   163,     0,     0,   184,
       0,   186,   188,   189,   190,   191,   192,   195,   203,   196,
     208,   210,   213,   212,   211,   216,   219,     0,     0,   225,
       0,     0,     0,   242,   244,   246,   248,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   268,   269,     0,
       0,   815,     0,   291,     0,   297,     0,   309,   315,   303,
     313,   306,   342,   384,   386,     0,     0,   385,   370,   383,
     454,   588,     0,   379,   368,   367,     0,   373,     0,     0,
     371,   360,     0,   395,     0,   405,   406,     0,     0,     0,
       0,     0,   480,     0,   486,     0,     0,     0,     0,   504,
       0,     0,     0,     0,     0,     0,     0,   551,     0,     0,
       0,     0,   565,     0,   568,     0,     0,     0,     0,     0,
     684,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   812,     0,     0,     0,     0,     0,     0,   764,     0,
     768,   777,   797,   798,   799,   800,   801,   794,   802,   784,
     785,   782,   795,   796,   808,     0,     0,     0,     0,     0,
       0,   840,     0,     0,   856,     0,   863,   867,   870,     0,
       0,     0,     0,     0,     0,   874,     0,     0,   889,     0,
     899,   900,   894,   895,   896,   898,   901,   919,   902,   917,
     903,     0,   897,   925,   922,   932,   934,   929,   950,   952,
     951,     0,   946,     0,     0,   940,     0,   930,     0,   959,
     936,   937,   938,   926,     0,     0,     0,     0,  1034,     0,
    1033,     0,     0,     0,  1035,     0,     0,   987,     0,     0,
    1009,  1011,     0,     0,     0,     0,     0,  1047,     0,     0,
    1032,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1008,     0,     0,     0,     0,     0,  1040,     0,  1000,     0,
       0,     0,     0,     0,  1030,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1048,     0,
       0,     0,     0,     0,     0,     0,     0,   962,     0,     0,
    1012,     0,     0,  1037,     0,     0,     0,     0,  1119,     0,
    1129,  1122,  1123,     0,  1113,  1114,     0,  1134,     0,  1132,
    1112,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1302,     0,  1339,  1345,  1326,
    1341,  1307,  1324,     0,  1333,  1335,     0,  1409,  1360,  1359,
    1362,  1369,     0,  1407,  1378,  1390,  1393,  1403,  1404,  1413,
       0,  1415,  1419,     0,     0,   403,     0,     0,   534,   338,
     524,  1085,   355,     0,     0,  1074,  1081,   819,  1098,     0,
     578,     0,     0,     0,   418,     0,   424,     0,     0,     0,
     657,     0,     0,   659,  1093,     0,   457,     0,     0,     0,
       0,   642,   646,     0,   648,     0,     0,   641,   222,     0,
     223,     0,   153,   441,   437,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   328,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1304,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   727,
       0,     0,   731,     0,   732,     0,   473,     0,  1410,     0,
     156,     0,     0,   175,     0,     0,     0,     0,     0,   178,
     181,     0,   180,   185,   187,   227,   226,   236,     0,     0,
       0,     0,     0,   827,     0,     0,   263,     0,   261,   265,
     267,     0,     0,   817,     0,   814,     0,   293,     0,     0,
     299,     0,   389,   378,     0,   452,   589,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   483,     0,     0,
     487,     0,     0,     0,   503,   505,   514,   517,   520,   523,
       0,     0,     0,   552,     0,     0,     0,   566,     0,     0,
       0,   721,     0,     0,     0,     0,     0,     0,     0,   685,
       0,     0,     0,     0,     0,     0,   718,     0,     0,   811,
     813,   729,   742,     0,     0,   753,     0,   757,   759,     0,
     766,     0,   769,   807,     0,     0,     0,     0,     0,   838,
     850,     0,     0,   851,     0,     0,     0,   871,     0,   876,
       0,     0,     0,   875,     0,   878,     0,     0,     0,     0,
     949,   947,     0,   954,     0,   941,     0,   924,   927,     0,
    1016,     0,  1022,     0,     0,  1064,  1015,  1013,     0,  1027,
    1065,  1066,  1007,  1006,  1010,  1025,  1026,   981,   984,  1046,
    1045,     0,   985,   986,     0,     0,  1054,     0,  1053,     0,
    1052,  1051,  1005,  1004,  1061,     0,     0,   979,   980,  1041,
    1042,   988,   989,   990,  1031,   964,  1001,   992,   991,  1039,
    1062,     0,     0,     0,  1028,  1003,  1002,   969,   972,   970,
     971,   973,   974,   975,   976,   967,   968,   966,  1044,   995,
       0,   983,   993,     0,   982,   965,  1043,  1029,   963,  1036,
       0,  1049,  1023,     0,  1038,  1096,  1117,     0,     0,  1120,
    1124,     0,  1115,     0,  1135,     0,     0,  1299,     0,     0,
    1300,  1305,     0,     0,     0,     0,  1170,     0,     0,  1176,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1182,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1303,  1342,     0,
       0,     0,  1363,  1370,     0,  1379,  1391,  1394,     0,     0,
       0,     0,   401,     0,   410,     0,     0,  1075,     0,   580,
     579,     0,   421,     0,   821,     0,     0,     0,   498,     0,
       0,     0,   456,   458,   462,  1095,     0,     0,   650,     0,
     804,   442,     0,   444,   445,   446,     0,   448,   449,   451,
       0,   438,     0,  1309,     0,     0,     0,     0,     0,   632,
       0,     0,   509,     0,     0,   600,     0,     0,   595,     0,
       0,   605,     0,     0,     0,  1304,   611,     0,     0,   627,
       0,     0,   617,     0,     0,   622,   728,   733,     0,     0,
       0,   478,   158,   157,   160,     0,   164,   165,     0,     0,
       0,     0,     0,   238,     0,     0,     0,     0,     0,   262,
       0,   277,   816,   294,   292,   300,   298,   387,   453,   818,
       0,   374,     0,   436,     0,   396,     0,   461,     0,     0,
     469,     0,   481,     0,   490,   494,   547,   539,     0,   540,
       0,   556,     0,     0,     0,     0,     0,   680,     0,     0,
       0,   687,     0,     0,   706,     0,   692,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   747,   743,     0,
       0,     0,     0,   765,   770,   806,     0,     0,     0,     0,
     839,     0,     0,   841,     0,   844,     0,     0,     0,   862,
       0,   880,   881,     0,     0,   879,   890,     0,   891,     0,
     904,     0,   948,   955,     0,   942,     0,   928,  1017,     0,
    1018,     0,  1014,  1063,     0,  1058,     0,     0,     0,  1057,
     977,   978,   997,   999,   998,   996,   994,  1050,  1024,     0,
    1130,     0,  1125,     0,     0,  1136,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1331,     0,  1412,     0,     0,     0,     0,
     402,     0,   381,  1099,   581,   423,     0,   820,     0,   499,
     662,     0,   658,  1094,   652,   649,     0,   224,   443,   447,
     450,   440,   439,     0,     0,     0,     0,     0,     0,     0,
     510,     0,     0,     0,     0,     0,     0,  1304,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   734,   735,     0,
     737,     0,   166,   169,     0,     0,   182,     0,     0,     0,
       0,   253,     0,     0,     0,     0,   295,     0,   301,     0,
     375,     0,   398,   397,     0,   465,     0,   491,   541,   542,
       0,     0,     0,     0,     0,   573,   681,     0,     0,     0,
     689,     0,     0,     0,   710,     0,     0,   688,     0,     0,
     708,     0,   693,     0,     0,     0,   719,     0,     0,   748,
     744,     0,     0,   754,   760,   771,   810,     0,     0,     0,
     835,   842,   843,   847,     0,     0,     0,   852,   877,   883,
     882,     0,     0,   905,     0,   906,     0,   956,     0,   943,
       0,  1019,     0,  1020,     0,  1056,  1055,     0,  1118,     0,
    1126,     0,  1116,  1143,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1169,     0,     0,     0,  1175,     0,
       0,     0,  1291,     0,  1222,     0,     0,  1251,     0,     0,
       0,     0,     0,     0,     0,  1181,     0,     0,     0,  1280,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1332,
    1334,  1406,     0,     0,     0,   584,     0,   422,     0,   500,
     661,   651,  1310,     0,     0,     0,   631,     0,   511,   599,
       0,   594,     0,   604,     0,     0,  1304,  1304,     0,   626,
       0,   616,     0,   621,     0,   736,   738,     0,     0,     0,
     240,     0,     0,     0,   256,     0,   278,   296,   302,   369,
     376,     0,   482,     0,     0,     0,     0,   575,     0,     0,
       0,   707,     0,     0,     0,   714,     0,   690,     0,     0,
     712,     0,     0,     0,     0,     0,     0,   749,   745,     0,
       0,     0,   809,     0,     0,     0,   853,     0,   845,   884,
     892,   907,     0,   908,     0,   957,   944,  1021,  1060,  1059,
    1131,  1127,     0,     0,  1146,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1186,     0,  1292,
       0,     0,     0,     0,     0,     0,     0,  1257,     0,     0,
       0,     0,     0,  1263,     0,     0,     0,  1188,     0,     0,
       0,     0,     0,     0,     0,  1198,     0,     0,  1202,     0,
       0,  1206,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   585,     0,     0,   425,     0,     0,   637,   633,   601,
     596,   606,     0,     0,     0,  1304,   628,   618,   623,     0,
       0,   241,     0,     0,     0,     0,   466,   555,   560,   562,
       0,   574,     0,     0,   694,     0,     0,     0,     0,     0,
       0,   709,     0,     0,   716,     0,     0,     0,     0,     0,
     750,     0,   772,     0,     0,     0,   848,   846,   909,     0,
     910,     0,  1145,  1144,     0,  1190,     0,     0,  1192,     0,
       0,  1150,     0,  1156,     0,  1173,     0,  1168,     0,  1179,
       0,  1174,     0,     0,  1293,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1185,     0,  1180,     0,     0,  1281,
       0,  1162,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1208,     0,     0,     0,     0,
       0,     0,     0,     0,   586,  1139,     0,     0,  1304,   609,
       0,     0,     0,     0,     0,     0,     0,     0,   257,     0,
     576,     0,     0,     0,   695,     0,     0,     0,   711,     0,
     697,     0,     0,     0,     0,     0,     0,     0,     0,   673,
       0,   746,     0,   773,     0,     0,   849,   911,     0,   912,
       0,     0,     0,     0,     0,     0,  1151,     0,  1157,     0,
       0,     0,     0,     0,  1187,  1294,     0,     0,     0,     0,
       0,  1283,     0,  1287,     0,     0,     0,     0,     0,  1264,
       0,  1282,  1163,     0,  1189,     0,     0,     0,     0,     0,
       0,     0,  1199,     0,     0,  1203,     0,     0,     0,  1210,
       0,  1216,     0,     0,     0,     0,     0,   587,     0,   636,
       0,     0,  1304,   612,     0,     0,   237,     0,     0,     0,
       0,     0,     0,   700,     0,   696,     0,   722,     0,     0,
     698,     0,     0,   713,     0,     0,     0,     0,     0,     0,
     751,     0,     0,     0,   913,     0,   914,  1147,     0,  1191,
    1194,     0,  1193,  1196,     0,     0,     0,  1172,     0,  1178,
       0,  1295,     0,     0,     0,  1252,     0,  1284,     0,  1288,
       0,     0,     0,  1184,     0,  1265,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1211,     0,  1217,     0,     0,     0,     0,     0,  1140,   638,
     610,     0,   167,   168,   239,     0,   826,   258,     0,     0,
       0,   701,     0,   715,   703,     0,   699,     0,     0,   691,
       0,     0,     0,     0,   774,     0,   830,     0,   915,  1148,
       0,     0,  1152,     0,  1158,     0,     0,     0,  1296,     0,
       0,     0,     0,  1285,     0,  1289,     0,  1258,     0,     0,
    1266,     0,  1164,     0,     0,     0,     0,     0,     0,     0,
    1200,     0,  1204,     0,  1207,     0,     0,     0,     0,     0,
       0,     0,   613,   250,     0,   259,   677,     0,   720,   702,
     704,     0,   717,     0,     0,   683,     0,   775,   832,     0,
       0,  1195,  1197,  1153,     0,  1159,     0,  1171,  1177,  1297,
       0,  1223,     0,     0,  1243,     0,     0,  1253,     0,  1286,
    1290,     0,  1183,  1267,     0,  1165,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1209,  1212,     0,  1218,     0,
       0,     0,     0,     0,     0,   678,     0,   705,     0,   686,
       0,     0,     0,  1154,     0,  1160,     0,     0,     0,     0,
       0,     0,     0,  1259,     0,  1268,     0,  1166,     0,  1227,
       0,  1229,     0,     0,     0,     0,     0,  1201,  1205,  1213,
       0,  1219,     0,  1239,     0,     0,  1247,     0,     0,  1255,
       0,     0,   251,     0,   679,     0,     0,   674,     0,     0,
     831,  1155,  1161,     0,     0,  1224,     0,     0,  1244,     0,
       0,     0,  1269,     0,  1167,     0,     0,  1231,     0,  1233,
       0,     0,     0,  1214,     0,  1220,     0,     0,     0,     0,
       0,     0,  1261,     0,   252,     0,     0,     0,   833,     0,
       0,     0,     0,     0,     0,     0,  1270,     0,  1228,  1230,
       0,     0,     0,     0,  1215,  1221,     0,  1240,     0,     0,
    1248,     0,     0,     0,   682,   675,     0,   669,     0,     0,
    1225,     0,  1245,     0,  1254,     0,  1271,     0,  1232,  1234,
    1235,     0,  1237,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1260,  1272,     0,     0,     0,
    1241,     0,  1249,     0,  1256,     0,     0,     0,     0,  1226,
    1246,  1273,     0,  1236,  1238,     0,     0,  1262,   670,     0,
       0,  1274,     0,  1242,  1250,     0,   671,     0,     0,  1275,
       0,     0,     0,  1298,  1276,     0,     0,     0,  1277,     0,
     672,     0,  1278,     0,     0,  1279,     0,     0,     0,     0,
       0,   676
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -2823, -2823, -2823,  3757, -2823, -2823, -2823, -2823, -2823, -2823,
    3631, -2823, -2823,  3634, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2174, -2823, -2823, -2823, -2823,
    3335,  3334, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,  3572, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823,  3582, -2823,  3308, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2084, -2823,  3658, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823,  2889, -2823, -2823,  2886, -2823, -2823, -2823, -2823,
   -2823, -2823,  2999, -2823,   158, -1175, -2823, -2823,  2925, -2823,
    3565, -2823,  -167, -2823,  3536, -2823,  -191,  -880,  -331, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823,  3468, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823,   652, -1194,  3531,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823,  2984,  -958, -2823, -2823,  3002,  -940, -2823,  -383, -2823,
    3512, -2822, -2823, -2823, -2823, -2823, -2823, -2823, -2823,  3427,
   -2823, -2823, -2823, -2823, -2823,  -441,  3460,  2732,  -190, -1722,
    -784,  -914, -2823, -2823,   952, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823,   178, -2823,  2947, -2823,
   -2823, -2823,  -336, -2823,   253,  -635, -2823, -2823,  2395, -2823,
   -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823, -2823, -1630,  1491, -2823, -2823, -2823, -2823, -2823, -2823,
   -2823,   552,  3659, -2823, -2823, -2823, -2823, -2823,  3763,    -6
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   171,   172,   173,   174,   175,   176,   177,   178,   557,
     558,   559,   560,   561,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   593,  2177,   594,   595,   596,   597,
    1105,  1108,   189,   600,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   199,   615,   200,   201,   202,   617,   203,
     204,   205,   206,   634,   207,   208,  1482,   635,  1149,  1154,
     209,   641,   642,   210,   647,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   649,   221,   222,   636,   223,
    2114,  2510,   637,   224,   225,   226,   650,   227,   228,   229,
     230,  1047,  1048,   231,  1051,  1052,   232,   233,   234,   235,
     236,   935,   936,   237,   663,  1768,   238,  1014,  1015,   239,
     665,   240,   667,   241,   669,   242,   671,   890,   891,   243,
     244,   245,   246,   247,   248,   249,   677,   250,   251,   252,
     253,   254,   255,   256,   257,   258,   259,   876,  1740,   916,
     260,  1026,   261,  1022,   262,  1028,   263,  1030,   264,  1034,
     265,  1036,   266,  1032,   267,  1009,   268,  1007,   269,   270,
     963,   964,  1596,   271,   946,   947,  1579,   272,   930,   273,
     689,  2841,   274,   275,   276,   277,   278,   279,   280,   696,
     281,   282,   700,   283,   284,   728,   691,  1800,   877,  1741,
     917,   931,   285,   286,   598,   287,   732,   288,   289,   290,
     291,   741,   742,   292,   293,   294,   295,   296,   297,   298,
     299,  1860,  1288,  1286,   300,  1294,   774,   301,   775,   302,
     919,   303,   371,   304,  1586,  1587,   305,  1562,  1563,   306,
     940,   307,   942,   308,  1400,   309,   795,   310,   311,   312,
     313,   314,  1995,  1464,   315,   316,   824,   317,   318,   319,
     320,   864,   825,   321,   870,   871,   322,   875,  1742,  2842
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     336,   337,   338,   341,   978,  2563,  1598,  1581,   607,  2139,
    1777,  1543,   893,  1543,   933,  1782,  1571,  1784,  2206,  1571,
     643,  3023,  1770,  1771,  1772,  1773,  1299,  3028,  1885,   569,
    2152,  2511,   739,   366,   909,  1215,   995,   332,   653,   644,
    1216,   327,   327,   394,   366,   332,   332,   366,  2154,   411,
     619,   327,  1805,  2136,   332,  2206,   675,   910,   676,   414,
    2206,   327,  2206,   780,   366,  2927,   324,  3365,   686,  2873,
     427,  2390,  3087,   333,  1891,   620,  2534,   367,  1652,   674,
    1202,   333,   333,   446,   562,   563,   323,   570,   367,   571,
     333,   367,   327,  1217,   327,   368,   369,  1845,  1472,   435,
     326,   332,   761,   621,   328,   328,   368,   369,   367,   368,
     369,   342,   485,  1399,   328,  1892,   621,   491,   492,  1788,
     332,   327,  3111,   327,   328,   332,   368,   369,   581,   332,
    3197,  1808,  1558,   346,   781,   347,   332,   333,  2670,  1218,
     332,   572,   334,   349,   332,   764,  2391,   524,  3201,  3216,
     334,   334,   782,   327,   335,   328,   333,   328,   535,   334,
     389,   333,   335,   335,  3235,   333,  1973,   332,  1558,  2891,
    3332,   335,   333,   352,   968,   783,   333,   332,   585,   586,
     333,  3335,   766,   702,   328,   604,   328,   327,   622,   437,
     332,  1300,   353,   703,  2137,   332,   618,  3366,   784,  3547,
     327,  3158,  3159,   333,   623,   704,   334,  3162,   705,   706,
     707,   708,   709,   333,   327,   944,   328,  2535,   335,  2416,
    1846,  1789,   354,   436,   767,   334,   333,  2874,   355,   332,
     334,   333,   412,   511,   334,   327,   939,   335,  1886,   582,
    1219,   334,   335,   428,  1229,   334,   335,   768,  3593,   334,
     328,   624,   710,   335,   785,  1790,  2278,   335,  2903,  3619,
     332,   335,   786,   328,   688,   333,   583,   625,   626,   332,
     787,   969,   334,  2928,  2247,  2104,   627,   328,   332,  2153,
    3088,   626,   334,   699,   335,   327,   743,   746,   332,  1974,
     996,   628,  1121,  2419,   335,   334,   333,  2155,   328,   779,
     334,   788,  2433,   789,   629,   333,  2313,   335,   325,  2914,
    1203,   790,   335,  2918,   333,  2642,  3635,  2935,   357,   797,
    1263,   711,  1809,  1186,   333,   712,  3638,  1653,  1654,  1655,
    3112,  3662,   886,   800,   334,   889,   910,  3269,  3198,  2671,
    3273,  3274,   897,   553,  1473,   584,   335,   899,   328,   358,
    3105,   645,  1213,   903,   618,  3047,  3202,  3217,  2206,  3687,
     905,   906,  3690,   908,  1231,   334,  1117,   555,  3108,   791,
    2892,  2080,  3236,  1189,   334,   920,  2248,   335,  3333,   801,
     802,  1278,   609,   334,  3175,  1791,   335,   928,   929,  3336,
    1840,  1301,  1830,   334,  2279,   335,   327,  3183,   941,   943,
     327,   630,  3541,   673,  1220,   335,   970,  3548,   771,   772,
    2249,   631,   955,   632,   573,   958,  1119,   740,   961,   356,
     713,   966,   327,   574,   512,   646,  2771,  2417,  2772,   714,
     610,  1535,   715,  3544,  1221,  1230,   547,   716,   717,   718,
     979,   719,  3613,  2108,   720,  2063,  1792,  3368,   721,  2904,
     327,  2905,  3374,  2643,   359,  1887,  3594,  1111,   475,   328,
     633,  3616,  1185,   328,  1831,   792,  2105,  3620,  2283,   327,
     327,   793,  1618,  1131,   564,   565,   566,   567,   568,  1132,
     332,  2149,  1013,   751,   752,   328,  2157,   736, -1306,  1020,
     327,   332,  1638,  3048,   722,  1212,   548,   549,   550,   551,
     552,  2420,   327,   575,   576,   577,   794,  2533,   578,   579,
    2434,  2206,  3456,   328,   360,  2206,   333,  2915,   723,   467,
     327,  2919,  3176,   724,  3636,  2936,  1133,   333,   366,  1054,
    1678,  1302,   328,   328,  3639,   332,   915,  2645,   753,  3663,
    2250,  1059,  1060,  1061,  1062,  1063,  1064,  1065,  1067,   361,
    1069,  1573,  3458,   328,  1070,   327,   332,   327,  3106,  2289,
     638,   327,  2310,  1545,  2150,   328,  2293,  3688,   327,  1937,
    3691,   333,   367,   327,  1090,   332,  3109,  2218,  2109,  1095,
    1096,  1099,   476,   328,  2233,   334,   366,   754,  3184,  2237,
     368,   369,   333,  1112,   327,  1113,   334,   335,  1116,   332,
    1040,  2251,   611,   612,   613,   614,  1680,   332,   335,  1144,
    3542,   333, -1306,  1150,  1152,  1153,  1155,   362,   328,  1150,
     328,  2284,  1832,   639,   328,   501,  1167,  1898,  2311,  1716,
     367,   328,  1172,   363,  1175,   333,   328,  2102,  2090,   640,
     334,  3545,  3457,   333,   366,  1184,  2652,   366,   368,   369,
    3614,  1190,   335,  1192,  1679,  1194,  1195,   328,   725,  1571,
     755,   334,  1933,   553,  1700,  1543,   554,   756,  2255,  3617,
     327,   726,   727,   335,  1206,  1207,  1208,  1134,   366,  1169,
     334,  1164,  2290,  1041,  1222,   757,  1533,   555,   367,  2294,
    2646,   367,   335,   803,   364,  2647,  2191,   737,   366,  1236,
    1237,  1238,  1239,  1240,   334,  1135,   368,   369,  1136,   368,
     369,   327,   334,  1250,   327,  1137,   335,   606,  1255,   327,
    1257,  3024,   367,  1259,   335,  2281,   365,  1261,  1042,  3154,
    1681,  1138,  1267,   328,  2880,   327,  3163,   327,   761,   370,
     368,   369,   367,  1277,  1279,  1280,  1281,   327,  3264,  1284,
    1285,  1918,   652,  1717,  1291,  2214,  1293,  2206,  2206,  1533,
     368,   369,   762,  1533,  1310,   656,   327,   376,   763,  2653,
     658,  2041,  3270,  1391,   328,  1393,  3254,   328,  1397,  1398,
    3370,   764,   328,  2256,  1405,   588,  1934,  1701,   765,   804,
     327,   805,   806,  1139,  1140,   377,  1123,  2765,   328,   807,
     328,  3626,   332,   366,   378,   327,   327,  1469,   602,   556,
     328,   808,   809,   810,   811,   812,  2194,   813,   766,   903,
    1533,   327,  1485,  2195,  1141,  1142,   379,   814,   815,   328,
     816,   380,   817,   818,  1533,  1498,  3025,   819,   333,   910,
    1503,   332,   381,   820,   821,   822,   823,   367,   758,   759,
     760,  3164,   327,   328,   910,  3281,   388,  2881,  2262,  3155,
     767,  3696,   332,   903,  1527,   368,   369,   660,   328,   328,
     865,   868,   872,  1170,  1537,  1538,  1125,   333,  3265,  2196,
    1539,  1793,   366,   768,   328,  2204,  1544,  2479,  2480,  3038,
    1583,  1547,  2205,  3755,   332,  2061,   910,  1552,   333,  3255,
    1554,  1555,  3271,   327,   332,  2790,   910,   334,   685,  3387,
    3371,   748,  1559,   911,  1143,   328,  1056,  1564,   390,   335,
     332,   327,  1569,  1570,   910,   327,   367,  1930,  1556,  2238,
     333,  2470,   327,   769,   589,   590,   591,   592,   910,   391,
     333,  2746,  1585,  1583,   368,   369,   334,  2674,  1588,  1589,
    1590,  2094,  1591,  1312,   332,  1593,   333,   332,   335,  1601,
    1602,   332,  2199,  1945,   332,  2675,   328,   334,   332,  1914,
    2498,   451,   910,  1605,   392,  3627,  2678,  2094,  2094,   335,
    2627,   332,   873,   874,   328,   541,   327,   327,   328,   330,
     333,  1588,  1588,   333,  2679,   328,   395,   333,  2859,   334,
     333,  1619,  1313,  1622,   333,  1871,   332,  1624,   910,   744,
     770,   335,  2971,   339,  2628,   396,  1630,   333,  1931,   399,
    1633,   335,  1636,   910,  1639,   334,  1642,   327,  1645,   910,
    1648,   402,  2206,   403,   332,  3697,  1651,   335,   327,  2613,
     910,   452,   333,  1659,   771,   772,  3037,   332,  1662,   328,
     328,   332,  1665,  1666,  1668,  1669,  1670,  1671,  1672,   334,
     332,  1674,   334,   949,  2843,   327,   334,  3756,   953,   334,
     333,   335,   328,   744,   335,  1314,   332,  1584,   335,   332,
    3356,   335,  3064,   333,  2815,   335,   334,   333,   404,   332,
     328,   405,  2849,  2860,  1702,   453,   333,  3170,   335,   773,
    1707,   328,   406,  3380,   455,  3029,  1710,  1711,   327,  1714,
    1715,   334,   333,   332,  3469,   333,  1722,   332,  2817,  2834,
    1720,   407,   332,   335,  1315,   333,  2835,  1721,   328,   327,
    1584,  2334,  2335,   343,  1872,  1010,   328,   340,  1731,   334,
     408,  1600,   327,   409,  1734,   410,  1735,  1736,   415,   333,
    2336,   335,   334,   333,  2614,  1024,   334,  2985,   333,  1746,
     327,  1749,   332,  3150,   335,   334,  1752,   332,   335,  1754,
    3151,   328,   332,   327,   456,   332,   332,   335,  1757,  2844,
     327,   334,  1760,  3444,   334,  1763,   416,  1765,  1766,  1767,
     418,   421,   328,   335,   334,   332,   335,   350,   333,   422,
    1779,  1780,  1781,   333,   424,   328,   335,  2850,   333,  1786,
    1787,   333,   333,  1796,  1797,  1798,  1799,   425,   334,   327,
    3030,   431,   334,   328,   393,   332,   432,   334,   457,  1908,
     335,   333,   438,   413,   335,   327,   328,   730,   731,   335,
     444,   372,   447,   328,   366,   327,   327,   449,   374,   426,
    1612,  1613,   445,  1825,  1826,  1827,  1828,  1829,   450,   327,
     327,   333,   896,   459,   460,  2106,   327,   334,  1841,   327,
    1844,   463,   334,  1849,   464,   327,   470,   334,   477,   335,
     334,   334,   328,   478,   335,   332,  1066,  2754,   367,   335,
    1151,   479,   335,   335,   480,  1309,  2337,  1873,   328,  1874,
     334,  1876,   962,   384,  1878,  1916,   368,   369,   328,   328,
    1884,   481,   335,   386,  2338,   482,  1893,   483,   327,   486,
     489,   333,   328,   328,   490,  1903,  1904,  1905,  1906,   328,
     334,   493,   328,   494,   397,  1396,   495,   400,   328,   496,
    1536,   327,   335,   419,   332,  1546,   332,   497,  1568,  1600,
     498,   499,   500,  1940,  1941,  1942,  3133,  3134,   502,  1947,
     332,   332,   503,  1951,   332,   332,   504,   505,  1658,   506,
     332,   507,  1960,   508,  1962,  1963,  1964,  1965,   509,   510,
     333,   328,   333,  1972,   517,  1975,   429,  1976,  1977,   518,
     334,   522,   523,   525,   526,   527,   333,   333,  1667,  1986,
     333,   333,   335,   528,   328,  1991,   333,  1993,  1994,   433,
    1996,  1997,  1998,  1999,  2000,  2001,  2002,  2003,  2004,  2005,
    2006,  2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,
    2016,  2017,  2018,  2019,  2020,  2021,  2022,  2023,  2024,  2025,
    2026,  2027,  2028,  2029,  2030,  2031,  2032,  2033,  2034,   334,
     332,   334,   529,   530,   531,   332,   332,  2040,  1745,   332,
     327,   335,   532,   335,   533,   334,   334,  2263,   534,   334,
     334,   536,  2044,  2045,  2046,   334,   332,   335,   335,   327,
     537,   335,   335,   538,   332,   539,   333,   335,   540,   327,
     542,   333,   333,   332,   332,   333,  1888,   545,  1890,   616,
    2053,  1894,   327,   605,   332,  3251,   648,   879,   880,   332,
     679,   881,   333,   683,   332,  1910,   684,  1824,   882,  1848,
     333,  2060,   883,   328,   884,   885,   894,   895,   439,   333,
     333,  2066,  2067,  1875,  1883,   332,   907,  1985,  2065,   332,
     333,   898,   328,  2096,   912,   333,  2074,   441,  2075,   900,
     333,  2076,   328,   913,   923,   334,  1564,   465,  2079,   332,
     334,   334,   924,   926,   334,   328,  2087,   335,   332,  2089,
     471,   333,   335,   335,  1585,   333,   335,  2092,   327,   925,
    1588,   334,  2095,  2097,  2098,  2099,   927,  2100,   934,   334,
     938,  1601,   327,   335,   950,   333,   951,   954,   334,   334,
     327,   335,   957,   960,   333,  1005,  1588,  1588,   967,   334,
     335,   335,   332,  2120,   334,   327,  2123,   971,  3351,   334,
    2126,   335,   332,  2241,  2129,   327,   335,  2132,  2265,  2268,
    2135,   335,  2273,   332,   972,   973,  2142,   975,   976,  2145,
     334,   328,  2148,   327,   334,   977,   473,  2151,   333,  2308,
     981,   327,   335,   327,   982,   328,   335,  2379,   333,  2162,
     487,   332,   332,   328,   334,   985,  2382,  2392,   513,   333,
    2500,   986,   989,   334,   991,   992,   335,  2400,   328,   332,
     994,   997,  2402,   515,   327,   335,   998,  2410,   328,   999,
    1006,  1008,  1011,   887,   332,  1012,  2180,   333,   333,   332,
    1097,  2181,  2182,  1017,  2184,  2185,   328,  2187,  2472,  1018,
    1021,   921,  2484,   332,   328,   333,   328,   334,  1023,   983,
     332,   987,  3441,  1025,  1027,  1031,  1033,   334,  2202,   335,
     333,  1035,  2536,   332,   332,   333,  2207,  1038,   334,   335,
    2208,  2547,  2209,  2210,  1029,  1039,  2211,   328,  1046,   333,
     335,  2213,  1000,  1071,  1055,  2219,   333,  1058,   332,  2221,
    2222,  2223,  1072,  1068,  1073,   332,   334,   334,  1098,   333,
     333,   332,   332,  2234,  2235,  2236,  1078,  1074,   335,   335,
    2240,  2242,  2243,  2244,   334,  2245,  1075,   332,  2252,  2253,
    2254,   332,  2257,  2258,   333,  2581,   335,   332,  2266,   334,
    2269,   333,  1076,   327,   334,  1077,  2585,   333,   333,  1079,
     332,   335,  1080,   332,  1081,  1084,   335,   332,   334,  2274,
    2275,  2276,  2277,   333,   332,   334,  2282,   333,  1086,  1087,
     335,  1088,  1089,   333,  2597,  2599,  2291,   335,   334,   334,
    1114,  2296,  1127,  2297,   332,   332,   333,  1128,   332,   333,
     335,   335,  2601,   333,  1129,  1130,   327,  1156,   332,  1160,
     333,  1161, -1421,   334,   332,  1162,   328,  2607,  1176,   327,
     334,  1002,  2611,  1526,  1163,   335,   334,   334,  2309,  1166,
     333,   333,   335,   327,   333,  1168,  2656,  1178,   335,   335,
    1187,  1201,   334,  2663,   333,  1211,   334,  1214, -1421,  1233,
     333,  1269,   334, -1421,   335,  2331,  2692,  2755,   335,  1234,
     826,  1241,  1243,  1282,   335,   334,  1283,  2346,   334,   328,
    1244,  1298,   334,  1246,  1044,  1621,   827,   335,  2629,   334,
     335,  2760,   328,  1247,   335,  2361,  2362,  1082,  2788,  1248,
     332,   335,   332,  1249,  2798,  2800,   328,  1251,  1252,   334,
     334,  1146,  1254,   744,  2380,  1256,  2383,  1265,  1266,  1270,
    2803,   335,   335,   334,  2820,   335,  2393, -1421,  1272,   334,
    2836,  2398,  1290,  1292,  1392,   335,   333,  2401,   333, -1421,
    2403,   335,  1403,  2847,  1404,  2406,  2922,   327,  1408,  2411,
    2988,  1463,  2414,  2415,  2418,  2421,  2422,  2423,  2424,  2425,
    2426,  2427,  2428,  2429,  2430,  2431,  2432,  2435,  2436,  2437,
    2438,  2439,  2440,  2441,  2442,   332,   328,  3004,  3077,  2447,
    2448,  2449,  2450,  2451,  2452,  2453,  2454,  2455,  2456,  1412,
   -1421,  3079,  1471,  1477,   828,   829,  1481,  3093,   332,  2460,
    1459,   332,  1462,   332,  1483,   334,  1629,   334,  1632,  1635,
     328,   333,  1497,  1475,  2468,  1157,  1499,   335,   830,   335,
    2473,  1500,   327,   332,   332,   332,   332,   332,  2476,  2477,
    1508,  1511,   332,  2478,   333,  1512,  1641,   333,  1521,   333,
    2483,  2485,  2487,  1522,  1484,  2489,  2490,  1523,  1487,  2491,
    1488,  1524,  1489,   332,  1525,   332,  1529,  2496,  1490,   333,
     333,   333,   333,   333,  1491,  1275,  1530,  1532,   333,  1540,
    1492,  1549,  2514,  3097,  2516,  2517,  1493,  1494,   332,  2521,
     334,  2523,   327,  1504,  2526,   328,   332,  2529,  1505,   333,
    1173,   333,   335,  2411,  2537,   831,   832,  2540,  1550,   327,
    2543,  1551,  1507,   334,  1560,  2548,   334,  2549,   334,  2550,
     332,  1509,   833,  1510,   333,   335,   332,   332,   335,   332,
     335,  2766,   333,  1561,  1513,  1514,  1515,  1516,   334,   744,
     334,   334,   334,  1517,  2564,  2565,  2566,   334,  2567,  2568,
     335,   335,   335,   335,   335,   328,   333,  1567,  3132,   335,
    1179,   332,   333,   333,  1603,   333,  2577,  1576,   334,  1604,
     334,  1577,   328,  2582,  1592,  1607,  2586,  1182,  1608,  2589,
     335,  3185,   335,  1730,  3188,  2593,  3191,   332,   834,   835,
     836,   837,  1611,   334,  1614,  2598,  2600,   333,  2602,  2603,
    2604,   334,   328,  2606,  2608,   335,  2609,  2610,  2612,  2615,
    3193,  2617,  2618,   335,  2619,  3204,  2621,  2622,  2623,  2624,
     327,  2625,  2626,   333,  1615,   334,   332,  1616,  2630,  1627,
    2631,   334,   334,  2632,   334,  1649,  2634,   335,  3219,  2637,
    2638,  2639,   327,   335,   335,  2641,   335,  1650,  2644,   332,
     327,  1656,  2650,  1673,  1676,  1682,  2654,  1683,  1684,   332,
    2657,  3221,   333,  1685,  2861,   327,   334,  1686,  2664,  3249,
    1465,  1466,  1467,  1468,  1687,  1688,  1689,  1690,   335,  1476,
    1691,  1692,  1693,   328,  1694,   333,  1695,   327,  1204,   838,
    1696,   839,   334,  3250,   840,   333,   841,  1699,  1260,  3258,
    3279,  1703,  3296,  1704,   335,   328,   327,  1705,  1706,  1712,
    1209,  1718,  1732,   328,   842,  1733,   843,   844,  1296,  1737,
    1738,  1739,  1743,  1744,  1750,   845,   332,   846,   328,  1753,
    2906,   334,  1528,  1304,  3298,   819,   866,   867,  1755,   847,
     848,   849,  1756,   335,   332,   850,   851,   852,   853,  1762,
     328,  1764,   854,   855,   334,  1306,   856,   857,   858,   859,
    3305,   860,   333,  1775,   334,  2693,   335,  2694,  1776,   328,
    2696,   861,   862,   863,  1409,  1778,   335,  1783,  2701,  2702,
     333,  2703,  2704,  1811,  2705,  2706,  2707,  2708,  2709,  2710,
    2711,  2712,  2713,  2714,  2715,  2716,  2717,   327,  2718,  2719,
    2720,  2721,  2722,  2723,  2724,  2725,  2726,  1812,   327,  1813,
     332,  2731,  2732,  2733,  2734,  2735,  2736,  2737,  2738,  2739,
    2740,  1814,  3311,  1815,  2742,  2743,  1816,  1817,   332,  1818,
     332,   334,  3313,   332,  1819,  2749,  1820,  1821,  1822,   332,
    2752,  1823,  1834,   335,  1836,  1837,   333,  2756,   327,   334,
    2758,  1843,  1850,  2761,  1851,  1852,  1853,   327,  1854,   332,
     328,   335,  1855,  1862,   333,  1541,   333,  1863,  1866,   333,
    2769,   328,  2770,  1867,  1868,   333,  1553,  1869,  2774,  2775,
    2776,  2777,  2778,   332,  2779,   327,  1870,  2781,  2782,   327,
    2783,  2784,  1889,  2785,  2786,   333,  3039,  2411,  1877,  3319,
     332,  2791,  2792,   332,  2793,  2794,  1879,  2795,  2796,  1880,
    1881,   328,   327,  2799,  2801,   334,  1565,  3322,   332,   333,
     328,  1882,  2804,  2805,  1895,  1747,  2808,   335,  2809,  2810,
    1897,  2812,  2813,   334,   332,   334,   333,  1896,   334,   333,
    1900,   332,   332,  1901,   334,   335,  2821,   335,   328,  1907,
     335,  1920,   328,  1838,   333,  2826,   335,  1856,  1928,  1958,
    1970,   327,   327,   327,   334,   327,  2830,  2831,  2832,  1944,
     333,  1971,  2837,  2838,  2839,   328,   335,   333,   333,  2845,
    1858,  2846,  2848,  2851,  1980,  2853,  2854,  2855,   334,  2857,
    2858,  1981,   327,  1984,  2862,  1989,  1990,  1992,  2865,  2037,
     335,  2868,  2869,  3339,  2038,   334,  3341,   332,   334,  2041,
    2875,  2876,  3353,   332,   332,  2047,   332,   335,   332,  2042,
     335,  2882,  2050,   334,   328,   328,   328,   332,   328,  1864,
    1926,  1967,  3390,  1978,  2054,   335,  2117,  2059,  2068,   334,
    2069,  2070,  2071,   333,  2072,  2077,   334,   334,  2081,   333,
     333,   335,   333,  2086,   333,   328,  3393,   327,   335,   335,
    1982,  2093,  2908,   333,  2910,  2911,  2912,  2913,  2916,  2917,
    2920,  2921,  2923,  3401,  2925,  2926,  3407,  2929,  2930,  2931,
    2932,  2933,  2934,  2937,  2938,  2101,  2940,  2941,  2942,  2943,
    2944,  3409,  2107,  2947,  2948,  2949,  2950,  2951,  2952,  2953,
    2954,  2955,  2956,  2957,  2958,  2959,  2960,  3415,  2961,   332,
    2110,  2963,   334,  2111,  3430,  3432,   327,  1644,   334,   334,
     328,   334,   327,   334,   335,  1987,   332,  2112,  2118,  2121,
     335,   335,   334,   335,   332,   335,   327,  2124,  2973,  2974,
    2975,  2127,  2977,  2130,   335,   333,  2980,  2133,  2982,  2140,
    2984,  2143,  2986,  2146,  2411,  2156,  2990,   327,  2992,  2198,
    2994,  2158,   333,  2163,  2169,  2172,  2173,  2174,  2175,  2176,
     333,  2999,   327,  3001,  3002,  2179,  3003,  3005,   327,   328,
    3440,   327,   327,  2183,  2035,   328,  3447,  3466,  3011,  3472,
    2082,  3474,  2186,   327,  3013,  3014,  3015,  3016,   332,   328,
    3490,  3018,  3019,  3020,  2084,  3282,  3022,  1647,  2188,  2189,
    3026,  2190,  2193,  2201,   334,  3031,  2203,  3032,  3033,  3034,
     328,  3035,  3036,  2217,   332,  2113,   335,  2220,   332,   332,
     332,   334,  3043,  3044,   333,   328,  2230,  2270,  3045,   334,
    2160,   328,  2280,   335,   328,   328,  2170,  2287,  2286,  2224,
    2226,   335,  2295,  2300,  2307,   332,   328,  2312,  2315,  2316,
     333,  2227,  2319,   327,   333,   333,   333,  2320,  2321,  3063,
    2322,  3065,  3066,  3067,  3068,  3069,  3070,  3071,  2323,  3072,
    3073,  3074,  3492,  3075,  3076,  3078,  2324,  3080,   327,  3081,
    3082,   333,  3083,  3084,  3085,  3086,   327,  3089,  3090,  3512,
    3091,  3092,  3094,   334,  3095,  3096,  3098,  3099,  3100,  3101,
    3102,  3103,  3104,  3107,  3110,   335,  3113,  3114,  3115,  3116,
    3117,  3118,  3119,  2325,   327,  3388,   328,  2326,  3122,   334,
    2327,  2228,  2328,   334,   334,   334,  2329,  3125,  3126,  2330,
    2332,   335,  2333,  2340,  2341,   335,   335,   335,   327,  2411,
     327,   328,  3135,  2342,   327,  2343,  2229,   327,   327,   328,
     334,  2344,  2347,  2348,  2259,  3142,  3143,  3144,  2349,  3145,
     332,  3513,   335,  2350,  2351,   327,  2352,  2353,  2354,  2355,
    2356,  2357,  3152,  3153,  3156,  2358,  3157,   328,   332,  2359,
    3160,  2360,  2261,   327,  2364,  2365,  3165,  3166,  3167,  3168,
    3169,  3516,  3528,  3533,  3171,  3172,   333,  3173,  3174,   327,
    2366,   328,  2367,   328,  2368,  2369,  2267,   328,  2301,  2370,
     328,   328,  2303,   327,   333,  2305,  2317,  2371,  3535,  3186,
    3187,  3189,  3190,  3192,  3194,  2372,  3196,  2373,   328,  3200,
    2374,  2375,  3203,  2396,  3205,  3206,  3207,  3208,  3209,  3210,
    3211,  2376,  3212,  3213,  2377,  3215,   328,  2378,  3218,  3220,
    3222,  2404,  3223,  3224,  3225,  3226,  3227,  3228,  3229,  2381,
    3230,  3231,   328,  3232,  3233,   334,  3234,  2551,  3237,  3238,
    3239,  3240,  3241,  3242,  3243,  2384,   328,   335,  2385,  3246,
    3247,  2635,  2386,   334,  2387,  2388,  3248,  2411,  2411,  2389,
     327,  2394,  2395,  2399,  3253,   335,  2407,  3256,  3257,  3259,
    2457,  2462,  2458,   332,  2465,  2474,  3261,  3262,  2475,  3263,
    3266,  3267,  2481,  2482,  3272,  2488,  2492,  2493,  2494,  3275,
    3276,  3277,  3278,  3280,  2495,  2113,  3283,  3284,  3285,  2513,
    2499,  2515,  2518,  3553,   327,   327,  2519,   332,  2522,   333,
    3292,  3293,  2524,  3294,  3295,  2525,  3297,  2527,  3299,   327,
    3300,  3555,  3301,   328,  3302,  2528,  3303,  2530,  2660,  3306,
    3307,  3308,  3309,  3310,  3312,  3314,  3315,  3316,  2531,  3317,
    2538,  3318,  3320,   333,  2539,  2532,  3323,   327,  3325,  3326,
    3327,  3328,  3329,  3330,  3331,  2541,  3334,  2542,  3337,  2544,
    3338,  3340,  3342,  3343,  3344,  3345,  3346,   328,   328,  2545,
    2546,  3350,  2665,  2690,  3352,  2411,  2552,  3355,   334,  2553,
    3358,  3359,   328,  3360,  2554,  3361,  3362,  2866,  2556,  3364,
     335,  3367,  2557,   332,  2558,  3369,  3372,  2559,  2751,   332,
    3376,  3377,  3378,  2569,  3379,  2571,  2572,  3381,  3382,  3383,
     328,   327,   334,   327,  2578,  2883,  2579,  3391,  2583,  3394,
     332,  3395,  2587,  3396,   335,  3398,   327,  3400,   327,   333,
    3402,  3403,  3404,  2590,  3406,   333,  3408,  2592,  3410,  3411,
    3412,  2594,  3414,  2595,  3416,   327,  3566,  3417,  2596,  3418,
    3419,  3420,  3421,  3422,  3423,  3424,   333,  3425,  3426,  2616,
    3427,  3428,  3429,   332,  3431,  2633,  3433,  3434,  3435,  3436,
    3437,   327,  2640,   332,   328,  2411,   328,  2649,   332,  2885,
    3568,  2887,  3445,  2651,  3448,  3449,  3450,  2655,  2658,   328,
    3452,   328,  2662,  2667,  2889,  3455,  2900,  2668,   334,   333,
    3460,  3461,  3462,  3463,   334,  3464,  3465,  3467,   328,   333,
     335,   327,  2672,  2967,   333,  3470,   335,  2673,  3471,  3473,
    3475,  2680,  3476,  2681,  3477,   334,  2682,  3480,  3481,  2683,
    3482,  2684,  2685,   327,   328,   332,  3488,   335,  3489,  3051,
    3491,  3493,  3494,  3495,  3496,  3497,  3498,  3499,  2686,  3501,
    2687,  3503,  2688,  3505,  2695,  3506,  2745,  3507,  3508,  3509,
    3510,  3511,  2750,  2753,  2757,  2411,  3575,  2759,   334,  3514,
     332,   333,  3595,  3517,   328,   332,   327,  2762,   334,  3053,
     335,  3521,  2763,   334,  3523,  3524,  2764,  3526,  3527,  3529,
     335,  3530,   327,  3599,   327,   335,   328,  3534,  2767,  3536,
    2780,  3178,  2773,  2787,  3543,  3546,   333,  2789,   327,  2797,
    2802,   333,  3551,  2806,  2811,  3554,  2822,  3556,  3557,  3558,
    3559,  3560,  3561,  3562,  2823,  3563,  2825,  3564,  2827,   332,
    3567,  3569,  3570,  3571,  3572,  3573,  3601,  2828,  3574,   328,
     334,  3576,  2829,   327,  3180,  2840,  3609,  2852,  2856,  2863,
    3580,  3611,   335,  3581,  3582,   328,  2864,   328,  2870,  2871,
    3287,  2872,  3289,   332,  2877,   333,  3588,  3589,   327,  3590,
    3591,   328,  3592,  2878,  2879,   334,  3384,   332,  3596,   332,
     334,  3600,  3602,  3603,  3604,  3605,  3606,   335,   327,  2893,
    2895,  3610,   335,  3612,  3615,  3618,  2896,  3621,  3623,   333,
    2898,  2902,  2924,   327,  3628,  3629,   328,  2939,  3622,   332,
    2965,  3478,  3634,   333,  3637,   333,  3640,  2966,  3641,  2969,
    3643,  2970,  2972,   332,  3645,  2976,  3646,  3648,  3650,  3651,
    3652,   328,  2978,  2979,   334,   327,  3483,   327,  3657,  3658,
    2981,  3659,  3660,  3642,  3661,   333,   335,  2983,  3647,  2987,
    3666,   328,  3667,  2989,   332,  2991,  3485,  2993,  3670,   333,
    3671,  3672,  2995,  3673,  3674,  3675,   328,  3677,   334,  2996,
    3000,  3539,  3680,   332,  3681,  3682,  3683,  3006,  3009,  3010,
     335,  3686,   334,  3689,   334,  3692,   327,  3693,  3012,  3017,
     333,  3698,  3021,  3699,   335,  3701,   335,  3703,   328,  3705,
     328,  3707,  3649,  3583,   327,  3585,  3711,  3713,  3027,   333,
    3714,  3042,  3715,  3716,   334,  3717,  3718,  3719,  3046,  3049,
    3720,  3050,  3721,  3722,   327,  3723,   335,  3724,   334,  3055,
    3056,  3727,   332,   332,   332,  3728,  3676,  3729,  3057,  3731,
     335,  3733,  3058,  3735,  3736,  3737,  3738,   332,  3059,   328,
    3706,  3742,  3710,  3060,  3597,  3745,  3061,  3746,  3121,   334,
    3123,  3749,  3124,  3127,  3128,  3129,  3752,   328,   333,   333,
     333,   335,  3653,  3757,  3758,  3130,  3760,  3131,   334,  3761,
    3136,  3762,  3712,   333,  3765,  3766,  3137,   328,  3138,  3769,
     335,  3771,  3655,  3773,  3141,  3774,  3726,   329,   331,  3146,
    3147,  3778,  3779,  3780,   344,   345,  3148,  3149,   348,  2501,
     351,  2502,  2503,  2504,  2505,  2506,  2507,  2508,  2509,  3161,
    3177,  3182,  3195,  3199,  3214,   373,   375,  3741,  3244,  3245,
    3260,  3268,  3286,   382,   383,   385,   387,   334,   334,   334,
    3304,  3321,  3324,  3347,  3348,   398,  3751,   401,  3349,   335,
     335,   335,   334,  3357,  3363,  3373,  3386,  3389,  3392,  3397,
    3399,   417,  3405,   420,   335,  3413,   423,  3438,  3439,  3442,
    3443,   430,  3446,  3451,   434,  3453,  3454,  3459,   440,   442,
     443,  3468,  3487,  3500,   448,  3502,  3504,   454,   458,  3515,
    3518,   461,   462,  3519,  3520,   466,  3522,   468,   469,  3525,
     472,   474,  3531,  3532,  3537,  3759,  3764,  3768,  3538,  3549,
     484,  3550,  3552,   488,  3565,  3577,  3579,  3607,  3608,  3624,
    3772,  3630,  3631,  3632,  3644,  3664,  3668,  3678,  3679,  3684,
    3685,  3694,  3695,  3700,  3702,  3704,  3708,   514,   516,  3709,
    3725,   519,   520,   521,  3730,  3732,  3734,  3739,  3740,  3743,
    3744,  3747,  3748,  3753,  3754,  3763,  3770,  3775,  3781,   544,
    1107,  1110,   904,  1159,  1574,   776,  1657,  1660,   777,  1625,
     546,   902,   914,   945,  1004,   580,   952,  1595,  1578,   587,
     974,   599,   601,   603,  1019,  1050,  1802,  2078,  1623,     0,
       0,   608,   778,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   869,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   651,
       0,   654,   655,     0,     0,   657,   659,   661,   662,     0,
     664,     0,   666,   668,   670,   672,     0,     0,   664,   668,
     672,     0,     0,   678,     0,   680,   681,   682,     0,     0,
       0,   687,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   690,   692,   693,     0,
       0,   694,   695,   697,   698,   701,     0,   729,   599,   599,
     733,   734,   735,   738,     0,     0,   745,   747,   749,     0,
       0,   750,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   796,   798,   799,     0,     0,
       0,     0,     0,     0,     0,     0,   878,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   888,     0,     0,   892,   892,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   901,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     666,   664,   918,     0,     0,   922,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   932,
     932,     0,   937,     0,     0,     0,   668,     0,     0,   672,
     670,     0,     0,   948,   878,     0,     0,     0,   918,   878,
       0,     0,   956,     0,     0,   959,     0,     0,     0,   965,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   690,     0,     0,     0,     0,     0,     0,   729,     0,
     980,     0,     0,     0,     0,   984,     0,     0,     0,   988,
       0,     0,     0,   990,     0,     0,   993,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1001,  1003,   678,
       0,     0,     0,     0,     0,     0,   878,     0,     0,     0,
       0,     0,  1016,     0,     0,     0,     0,   692,     0,     0,
       0,     0,     0,     0,     0,     0,   878,     0,     0,     0,
       0,     0,     0,     0,     0,  1037,     0,     0,     0,     0,
       0,     0,     0,  1043,  1045,     0,     0,     0,     0,     0,
    1049,   697,  1053,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1057,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1083,     0,  1085,     0,
       0,     0,     0,     0,  1091,  1092,  1093,  1094,     0,     0,
       0,  1100,  1101,  1102,  1103,  1104,     0,  1106,  1106,  1109,
    1109,     0,     0,     0,     0,     0,  1115,     0,  1118,  1120,
       0,  1122,  1124,  1126,     0,     0,     0,     0,     0,  1145,
    1147,  1148,     0,     0,     0,     0,     0,  1158,     0,     0,
       0,     0,     0,  1165,     0,     0,     0,  1171,     0,     0,
       0,     0,  1174,     0,     0,     0,     0,     0,  1177,     0,
    1180,  1181,  1183,     0,     0,  1120,     0,  1122,     0,  1188,
       0,  1191,     0,  1193,     0,     0,     0,  1196,     0,  1197,
       0,  1198,     0,  1199,     0,  1200,     0,     0,     0,     0,
       0,     0,  1205,     0,     0,     0,  1210,     0,  1120,     0,
    1122,     0,     0,     0,     0,  1223,  1224,  1225,  1226,     0,
    1227,  1228,     0,     0,  1232,     0,     0,  1235,     0,     0,
       0,     0,     0,     0,  1242,     0,     0,  1245,     0,     0,
       0,     0,     0,     0,     0,  1253,     0,     0,     0,     0,
    1258,     0,     0,     0,     0,     0,     0,  1262,  1264,     0,
       0,     0,  1268,     0,     0,     0,     0,     0,  1271,     0,
    1273,  1274,  1276,     0,     0,     0,     0,     0,     0,     0,
    1287,  1289,     0,     0,     0,     0,  1295,  1297,     0,  1303,
    1305,  1307,  1308,     0,  1311,     0,  1316,     0,     0,     0,
       0,     0,     0,     0,     0,  1394,  1395,     0,     0,  1401,
    1402,     0,     0,     0,  1406,  1407,     0,  1410,     0,  1411,
       0,     0,  1455,  1456,  1457,  1458,     0,  1460,  1461,     0,
       0,     0,     0,     0,     0,     0,     0,  1470,     0,  1474,
       0,     0,     0,  1478,  1479,     0,  1480,     0,     0,     0,
       0,     0,  1486,     0,     0,     0,     0,     0,     0,     0,
       0,  1495,  1496,     0,     0,     0,     0,  1501,  1502,     0,
       0,     0,  1506,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1518,  1519,  1520,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1531,     0,  1534,
       0,  1122,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1542,     0,   892,     0,   892,   892,     0,     0,
       0,     0,  1548,     0,     0,     0,     0,  1147,     0,     0,
       0,     0,     0,  1557,     0,     0,     0,     0,     0,   918,
       0,     0,     0,     0,     0,     0,     0,     0,  1566,     0,
       0,     0,     0,   932,     0,  1572,   932,   932,   937,     0,
    1575,     0,     0,     0,     0,     0,     0,     0,     0,  1580,
    1582,     0,  1534,     0,     0,   918,  1534,     0,     0,     0,
       0,     0,     0,     0,     0,  1594,  1597,  1599,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1606,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1609,     0,     0,  1610,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1617,
       0,  1620,     0,  1534,     0,     0,     0,  1016,     0,  1626,
       0,     0,     0,     0,  1628,     0,     0,  1534,  1631,     0,
    1634,     0,  1637,     0,  1640,     0,  1643,     0,  1646,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1049,     0,     0,     0,  1053,     0,  1661,     0,     0,  1663,
    1664,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1675,     0,  1677,  1413,  1414,     0,     0,     0,     0,
       0,     0,     0,     0,  1415,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1697,  1698,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1708,     0,  1709,     0,     0,     0,  1713,     0,
       0,     0,     0,  1122,  1719,     0,  1723,  1724,  1725,  1726,
       0,     0,     0,     0,     0,     0,  1727,     0,  1728,     0,
    1729,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1748,
       0,     0,  1751,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1758,
    1759,     0,  1761,     0,     0,     0,     0,     0,     0,  1769,
    1769,  1769,  1769,  1769,  1774,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1785,     0,     0,     0,
    1794,  1795,     0,     0,     0,     0,  1801,  1801,  1803,  1804,
    1769,  1806,     0,  1807,     0,  1810,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1833,     0,     0,
       0,  1835,     0,     0,     0,     0,  1839,     0,  1842,     0,
       0,  1847,     0,     0,     0,     0,     0,     0,     0,  1857,
       0,  1859,     0,  1861,     0,     0,     0,  1865,     0,     0,
    1416,  1417,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1899,     0,     0,  1902,     0,     0,     0,     0,     0,  1909,
       0,  1911,  1912,  1913,  1915,  1917,  1919,     0,  1921,  1922,
    1923,  1924,  1925,  1927,     0,  1929,  1932,     0,  1935,  1936,
    1938,  1939,     0,     0,     0,  1943,     0,  1946,     0,  1948,
    1949,  1950,     0,  1952,  1953,  1954,  1955,  1956,  1957,     0,
    1959,     0,  1961,     0,     0,     0,     0,  1966,  1968,  1969,
       0,     0,     0,     0,     0,     0,     0,     0,  1979,     0,
       0,     0,     0,  1983,     0,     0,     0,     0,     0,  1988,
       0,     0,  1418,  1419,  1420,  1421,     0,  1422,     0,  1423,
    1424,  1425,  1426,  1427,  1428,  1429,  1430,     0,  1431,  1432,
    1433,  1434,  1435,  1436,  1437,  1438,  1439,  1440,  1441,  1442,
    1443,  1444,     0,  1445,  1446,  1447,  1448,  1449,  1450,  1451,
    1452,  1453,  1454,     0,     0,     0,     0,     0,  2036,     0,
       0,     0,     0,     0,  2039,     0,     0,     0,     0,     0,
       0,     0,   826,     0,     0,     0,     0,     0,     0,     0,
       0,  2043,     0,     0,     0,     0,  2048,     0,   827,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2049,     0,     0,  2051,     0,     0,     0,  2052,
       0,     0,     0,     0,  2055,  2056,     0,     0,     0,     0,
       0,     0,  2057,  2058,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2062,     0,     0,     0,  1122,  2064,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   892,     0,
       0,     0,  2073,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2083,     0,  2085,     0,     0,     0,   932,     0,  2088,     0,
       0,  1582,     0,  2091,     0,  2091,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   828,   829,  1599,     0,
    2103,     0,  2103,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2115,  2116,     0,     0,     0,     0,     0,     0,
     830,  2119,     0,     0,  2122,     0,     0,     0,     0,  2125,
       0,     0,  2128,     0,     0,  2131,     0,     0,  2134,     0,
       0,  2138,     0,     0,  2141,     0,     0,  2144,     0,     0,
    2147,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1769,     0,  2159,  2161,     0,     0,
       0,  2164,  2165,  2166,  2167,  2168,     0,     0,  2171,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   831,   832,     0,
       0,     0,     0,  2178,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   833,     0,     0,     0,     0,     0,
       0,     0,  2192,     0,     0,     0,  2197,     0,  2200,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2212,     0,     0,
       0,  2215,  2216,     0,     0,     0,     0,     0,     0,     0,
       0,  2225,     0,  2225,  2225,  2225,  2225,     0,  2231,  2232,
     834,   835,   836,   837,     0,     0,     0,     0,  2239,     0,
       0,     0,     0,     0,     0,  2246,     0,     0,     0,     0,
       0,     0,     0,  2260,     0,  2260,  2264,     0,  2225,     0,
       0,     0,  2271,  2272,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2285,     0,     0,     0,
       0,     0,     0,  2288,     0,     0,  2292,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2298,  2299,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2302,  2304,     0,  2306,
       0,   838,     0,   839,     0,     0,   840,     0,   841,     0,
    2314,     0,     0,     0,     0,  2318,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   842,     0,   843,   844,
       0,     0,     0,     0,     0,     0,     0,   845,  2339,   846,
       0,     0,     0,     0,     0,  2345,     0,     0,     0,     0,
       0,   847,   848,   849,     0,     0,     0,   850,   851,   852,
     853,     0,     0,     0,   854,   855,  2363,     0,   856,   857,
     858,   859,     0,   860,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   861,   862,   863,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2397,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2405,     0,     0,     0,     0,  2408,  2409,     0,  2412,
    2413,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2443,  2444,  2445,  2446,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2459,     0,     0,  2461,
       0,     0,     0,     0,  2463,  2464,     0,     0,  2466,  2467,
       0,     0,     0,     0,     0,     0,  2469,  2471,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2486,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2497,     0,
       0,     0,     0,     0,   729,     0,     0,     0,     0,  2512,
       0,     0,     0,     0,     0,     0,     0,     0,  2520,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2225,     0,     0,     0,     0,     0,     0,  2555,     0,     0,
       0,     0,     0,     0,  2560,     0,     0,     0,     0,     0,
    2561,  2562,  2178,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2570,     0,     0,     0,     0,  2573,
       0,  2574,  2575,     0,  2576,     0,     0,     0,     0,     0,
       0,  2580,     0,     0,  2584,     0,     0,  2588,     0,     0,
       0,  2591,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2605,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2620,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2636,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2648,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2659,  2661,     0,     0,     0,     0,     0,     0,  2666,
       0,     0,     0,     0,     0,     0,  2669,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2676,     0,
    2677,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2689,  2691,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2697,  2698,     0,     0,  2699,  2700,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2727,  2728,  2729,  2730,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2741,     0,     0,     0,     0,  2744,     0,     0,
       0,     0,  2747,  2748,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2768,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2807,     0,  2561,     0,     0,     0,
       0,     0,     0,  2814,     0,     0,     0,  2816,     0,  2818,
       0,     0,     0,  2819,     0,     0,     0,     0,     0,     0,
       0,  2824,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2833,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1317,     0,     0,     0,     0,  2867,
       0,     0,     0,     0,     0,     0,     0,  1318,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1319,  2884,     0,  2886,     0,  1320,  2888,     0,  2890,
       0,     0,     0,     0,     0,  1321,     0,  2894,  1322,  1323,
       0,  2897,     0,  1324,  1325,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2899,     0,  2901,     0,     0,     0,
    2907,     0,  2909,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2945,  2946,     0,     0,     0,     0,     0,     0,     0,  1326,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2962,
       0,  2964,     0,     0,     0,     0,     0,     0,     0,  2968,
    1327,  1328,  1329,  1330,     0,     0,     0,     0,     0,  1331,
       0,     0,     0,  1332,     0,     0,  1333,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1334,     0,
    1335,     0,     0,     0,     0,     0,     0,  2997,  2998,     0,
       0,     0,     0,     0,  1336,     0,     0,     0,  3007,     0,
    3008,     0,     0,     0,     0,     0,     0,     0,  1337,  1338,
       0,     0,     0,     0,     0,     0,     0,  1339,     0,  1340,
       0,     0,  1341,     0,  1342,     0,  1343,  1344,     0,     0,
       0,  1345,     0,  1346,  1347,  1348,  1349,     0,     0,     0,
       0,     0,     0,     0,     0,  3040,     0,     0,  3041,     0,
       0,  1350,     0,  1351,     0,     0,     0,     0,     0,  1352,
       0,     0,     0,     0,     0,     0,     0,  3052,     0,  3054,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  3062,  1353,     0,
       0,     0,     0,  1354,  1355,  1356,     0,  1357,     0,  1358,
       0,  1359,     0,     0,     0,     0,  1360,  1361,  1362,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1363,     0,     0,     0,
       0,     0,     0,     0,     0,  1364,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3120,     0,  1365,  1366,  1367,
    1368,  1369,  1370,  1371,     0,  1372,  1373,  1374,  1375,     0,
       0,     0,     0,     0,     0,  1376,  1377,     0,  1378,     0,
       0,     0,  1379,     0,     0,     0,     0,     0,     0,     0,
    3139,  3140,     0,     0,     0,  1380,     0,  1381,     0,     0,
       0,     0,  1382,     0,     0,     0,     0,  1383,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1384,     0,     0,
       0,  1385,  1386,  1387,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3179,     0,  3181,     0,     0,
       0,  1388,     0,     0,  1389,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1390,     0,     0,     0,     0,     0,     0,
       0,     0,  3252,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  3288,     0,  3290,     0,     0,  3291,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3354,     0,  2178,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  3375,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  3385,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2561,
    2178,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  3479,     0,     0,     0,     0,
       0,  3484,     0,  3486,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2561,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  3540,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  3578,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  3584,     0,  3586,
       0,     0,     0,  3587,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  3598,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  3625,     0,     0,     1,     0,     2,     3,     0,     0,
    3633,     0,     4,     0,     5,     0,     0,     0,     0,     0,
       0,     0,     6,     0,     7,     8,     9,    10,    11,     0,
       0,     0,    12,  3654,     0,  3656,     0,     0,     0,     0,
       0,    13,    14,    15,     0,     0,     0,     0,  3665,     0,
       0,    16,     0,    17,    18,    19,  3669,    20,     0,     0,
       0,     0,     0,     0,     0,    21,     0,     0,     0,     0,
       0,    22,    23,     0,     0,    24,     0,     0,     0,    25,
      26,    27,    28,     0,     0,     0,    29,     0,    30,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      31,     0,     0,    32,    33,    34,    35,   543,    36,     0,
      37,    38,     0,     0,     0,     0,     0,     0,     0,    39,
       0,     0,     0,    40,    41,     0,     0,     0,    42,    43,
       0,     0,    44,    45,     0,     0,     0,     0,     0,     0,
       0,     0,    46,     0,     0,     0,     0,     0,     0,     0,
      47,     0,     0,    48,    49,     0,     0,    50,     0,     0,
      51,  3750,     0,    52,    53,    54,    55,    56,    57,    58,
      59,    60,    61,    62,    63,     0,     0,     0,     0,    64,
      65,     0,     0,     0,    66,  3767,    67,     0,    68,    69,
      70,    71,     0,    72,     0,    73,     0,  3776,    74,  3777,
       0,     0,     0,    75,     0,     0,    76,    77,    78,    79,
      80,     0,     0,    81,     0,     0,    82,    83,    84,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    85,     0,     0,     0,     0,     0,     0,
       0,     0,    86,     0,     0,     0,     0,     0,    87,     0,
       0,     0,    88,     0,     0,     0,     0,    89,     0,    90,
      91,     0,     0,     0,     0,     0,    92,     0,    93,     0,
       0,    94,    95,     0,     0,    96,     0,     0,    97,    98,
      99,     0,   100,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   101,   102,   103,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   104,   105,     0,
       0,     0,     0,   106,   107,   108,     0,   109,     0,     0,
       0,     0,     0,   110,     0,   111,   112,     0,     0,     0,
     113,     0,     0,     0,     0,     0,     0,   114,     0,     0,
     115,   116,     0,     0,     0,   117,     0,     0,   118,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   119,     0,     0,   120,     0,     0,     0,     0,
     121,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   122,     0,   123,     0,   124,
       0,     0,   125,   126,   127,     0,   128,   129,     0,   130,
       0,     0,   131,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   132,     0,     0,     0,     0,   133,   134,   135,
       0,     0,     0,     0,     0,     0,     0,   136,   137,   138,
     139,   140,   141,     0,   142,     0,     0,     0,     0,     0,
     143,   144,     0,     0,     0,     0,     0,     0,   145,     0,
       0,     0,   146,   147,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   148,   149,     0,     0,     0,
       0,     0,   150,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   151,     0,   152,   153,   154,
     155,   156,   157,   158,   159,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   160,   161,     0,     0,   162,
     163,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   164,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   165,   166,   167,     0,     0,
     168,   169,     0,     1,     0,     2,     3,     0,     0,     0,
       0,     4,   170,     5,     0,     0,     0,     0,     0,     0,
       0,     6,     0,     7,     8,     9,    10,    11,     0,     0,
       0,    12,     0,     0,     0,     0,     0,     0,     0,     0,
      13,    14,    15,     0,     0,     0,     0,     0,     0,     0,
      16,     0,    17,    18,    19,     0,    20,     0,     0,     0,
       0,     0,     0,     0,    21,     0,     0,     0,     0,     0,
      22,    23,     0,     0,    24,     0,     0,     0,    25,    26,
      27,    28,     0,     0,     0,    29,     0,    30,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    31,
       0,     0,    32,    33,    34,    35,     0,    36,     0,    37,
      38,     0,     0,     0,     0,     0,     0,     0,    39,     0,
       0,     0,    40,    41,     0,     0,     0,    42,    43,     0,
       0,    44,    45,     0,     0,     0,     0,     0,     0,     0,
       0,    46,     0,     0,     0,     0,     0,     0,     0,    47,
       0,     0,    48,    49,     0,     0,    50,     0,     0,    51,
       0,     0,    52,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    63,     0,     0,     0,     0,    64,    65,
       0,     0,     0,    66,     0,    67,     0,    68,    69,    70,
      71,     0,    72,     0,    73,     0,     0,    74,     0,     0,
       0,     0,    75,     0,     0,    76,    77,    78,    79,    80,
       0,     0,    81,     0,     0,    82,    83,    84,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    85,     0,     0,     0,     0,     0,     0,     0,
       0,    86,     0,     0,     0,     0,     0,    87,     0,     0,
       0,    88,     0,     0,     0,     0,    89,     0,    90,    91,
       0,     0,     0,     0,     0,    92,     0,    93,     0,     0,
      94,    95,     0,     0,    96,     0,     0,    97,    98,    99,
       0,   100,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   101,   102,   103,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   104,   105,     0,     0,
       0,     0,   106,   107,   108,     0,   109,     0,     0,     0,
       0,     0,   110,     0,   111,   112,     0,     0,     0,   113,
       0,     0,     0,     0,     0,     0,   114,     0,     0,   115,
     116,     0,     0,     0,   117,     0,     0,   118,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   119,     0,     0,   120,     0,     0,     0,     0,   121,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   122,     0,   123,     0,   124,     0,
       0,   125,   126,   127,     0,   128,   129,     0,   130,     0,
       0,   131,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   132,     0,     0,     0,     0,   133,   134,   135,     0,
       0,     0,     0,     0,     0,     0,   136,   137,   138,   139,
     140,   141,     0,   142,     0,     0,     0,     0,     0,   143,
     144,     0,     0,     0,     0,     0,     0,   145,     0,     0,
       0,   146,   147,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   148,   149,     0,     0,     0,     0,
       0,   150,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   151,     0,   152,   153,   154,   155,
     156,   157,   158,   159,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   160,   161,     0,     0,   162,   163,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   164,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   165,   166,   167,     0,     0,   168,
     169,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   170
};

static const yytype_int16 yycheck[] =
{
       6,     7,     8,     9,   445,  2179,   964,   947,   198,  1639,
    1204,   891,   343,   893,   397,  1209,   930,  1211,  1740,   933,
      32,  2843,  1197,  1198,  1199,  1200,    84,  2849,    17,    28,
     249,  2115,    19,    38,   370,    10,   249,    76,   228,    51,
      15,   181,   181,    49,    38,    76,    76,    38,   249,   249,
      22,   181,  1227,   108,    76,  1777,   247,   175,     4,    65,
    1782,   181,  1784,    12,    38,   249,     7,   249,   258,   249,
      76,   174,   249,   112,    17,    47,   108,    82,   249,   246,
     103,   112,   112,    89,    44,    45,   249,    86,    82,    88,
     112,    82,   181,    68,   181,   100,   101,   249,   126,   125,
     249,    76,    39,    75,   244,   244,   100,   101,    82,   100,
     101,   249,   118,     9,   244,    58,    75,   123,   124,   113,
      76,   181,   249,   181,   244,    76,   100,   101,    22,    76,
     249,   249,   916,   249,    83,   249,    76,   112,   249,   114,
      76,   140,   181,   249,    76,    82,   249,   153,   249,   249,
     181,   181,   101,   181,   193,   244,   112,   244,   164,   181,
     372,   112,   193,   193,   249,   112,   249,    76,   952,   249,
     249,   193,   112,   249,   152,   124,   112,    76,   184,   185,
     112,   249,   119,     8,   244,   191,   244,   181,   160,   372,
      76,   249,   249,    18,   249,    76,   202,   379,   147,   249,
     181,  3023,  3024,   112,   176,    30,   181,  3029,    33,    34,
      35,    36,    37,   112,   181,   406,   244,   249,   193,   249,
     372,   215,   249,   249,   161,   181,   112,   407,   249,    76,
     181,   112,   432,   249,   181,   181,   403,   193,   227,   133,
     215,   181,   193,   125,   249,   181,   193,   184,   249,   181,
     244,   223,    77,   193,   203,   249,   111,   193,   249,   249,
      76,   193,   211,   244,   270,   112,   160,   239,   240,    76,
     219,   249,   181,   457,   113,   249,   248,   244,    76,   498,
     457,   240,   181,   413,   193,   181,   292,   293,    76,   372,
     503,   263,   412,   249,   193,   181,   112,   498,   244,   305,
     181,   250,   249,   252,   276,   112,   445,   193,   249,   249,
     333,   260,   193,   249,   112,   111,   249,   249,   249,   459,
     407,   146,   440,   412,   112,   150,   249,   498,   499,   500,
     457,   249,   338,    46,   181,   341,   175,  3159,   457,   450,
    3162,  3163,   348,   302,   372,   239,   193,   353,   244,   249,
     249,   363,   412,   359,   360,   111,   457,   457,  2080,   249,
     366,   367,   249,   369,   700,   181,   125,   326,   249,   318,
     450,  1565,   457,   412,   181,   381,   215,   193,   457,    92,
      93,   412,   362,   181,   111,   379,   193,   393,   394,   457,
     412,   449,   177,   181,   249,   193,   181,   249,   404,   405,
     181,   373,   249,   245,   379,   193,   384,   457,   345,   346,
     249,   383,   418,   385,   413,   421,   606,   404,   424,   440,
     245,   427,   181,   422,   440,   437,  2510,   457,  2512,   254,
     410,   412,   257,   249,   409,   440,    75,   262,   263,   264,
     446,   266,   249,   249,   269,   412,   440,  3269,   273,   440,
     181,   442,  3274,   249,   249,   444,   457,   408,   125,   244,
     432,   249,   652,   244,   249,   414,   440,   457,   249,   181,
     181,   420,   381,    14,   434,   435,   436,   437,   438,    20,
      76,   249,   488,   213,   214,   244,  1661,   111,   125,   495,
     181,    76,   378,   249,   319,   685,   135,   136,   137,   138,
     139,   457,   181,   502,   503,   504,   455,  2137,   507,   508,
     457,  2233,   249,   244,   249,  2237,   112,   457,   343,   125,
     181,   457,   249,   348,   457,   457,    67,   112,    38,   535,
     125,   589,   244,   244,   457,    76,   378,   249,   268,   457,
     379,   547,   548,   549,   550,   551,   552,   553,   554,   249,
     556,   934,  3374,   244,   560,   181,    76,   181,   457,   249,
     160,   181,   249,   894,   332,   244,   249,   457,   181,   300,
     457,   112,    82,   181,   580,    76,   457,   375,   384,   585,
     586,   587,   249,   244,  1778,   181,    38,   317,   440,  1783,
     100,   101,   112,   599,   181,   601,   181,   193,   604,    76,
     249,   440,   582,   583,   584,   585,   125,    76,   193,   615,
     457,   112,   249,   619,   620,   621,   622,   249,   244,   625,
     244,   402,   407,   223,   244,   125,   632,   338,   315,   125,
      82,   244,   638,   249,   640,   112,   244,  1595,  1578,   239,
     181,   457,   379,   112,    38,   651,   249,    38,   100,   101,
     457,   657,   193,   659,   249,   661,   662,   244,   483,  1573,
     390,   181,   125,   302,   249,  1545,   305,   397,   134,   457,
     181,   496,   497,   193,   680,   681,   682,   218,    38,    55,
     181,   372,   372,   332,   690,   415,   876,   326,    82,   372,
     402,    82,   193,   406,   249,   407,   375,   321,    38,   705,
     706,   707,   708,   709,   181,   246,   100,   101,   249,   100,
     101,   181,   181,   719,   181,   256,   193,   378,   724,   181,
     726,   134,    82,   729,   193,   321,   249,   733,   377,   249,
     249,   272,   738,   244,   249,   181,   134,   181,    39,   249,
     100,   101,    82,   749,   750,   751,   752,   181,   249,   755,
     756,   338,   378,   249,   760,   375,   762,  2479,  2480,   949,
     100,   101,    63,   953,   770,   378,   181,   249,    69,   372,
     378,   125,   249,   779,   244,   781,   249,   244,   784,   785,
     249,    82,   244,   249,   790,    49,   249,   372,    89,   502,
     181,   504,   505,   334,   335,   249,   125,   249,   244,   512,
     244,    75,    76,    38,   249,   181,   181,   813,   223,   448,
     244,   524,   525,   526,   527,   528,   242,   530,   119,   825,
    1010,   181,   828,   249,   365,   366,   249,   540,   541,   244,
     543,   249,   545,   546,  1024,   841,   249,   550,   112,   175,
     846,    76,   249,   556,   557,   558,   559,    82,   578,   579,
     580,   249,   181,   244,   175,   249,   249,   372,   249,   379,
     161,    75,    76,   869,   870,   100,   101,   378,   244,   244,
     318,   319,   320,   249,   880,   881,   125,   112,   379,   157,
     886,  1217,    38,   184,   244,   242,   892,  2081,  2082,   249,
      62,   897,   249,    75,    76,   249,   175,   903,   112,   372,
     906,   907,   379,   181,    76,  2535,   175,   181,   378,   249,
     379,   378,   918,   249,   455,   244,   378,   923,   249,   193,
      76,   181,   928,   929,   175,   181,    82,   121,   249,   375,
     112,   375,   181,   234,   198,   199,   200,   201,   175,   249,
     112,   375,   948,    62,   100,   101,   181,   231,   954,   955,
     956,  1586,   958,   249,    76,   961,   112,    76,   193,   965,
     966,    76,   157,   338,    76,   249,   244,   181,    76,   225,
     249,   125,   175,   979,   249,   249,   231,  1612,  1613,   193,
     249,    76,   536,   537,   244,   372,   181,   181,   244,   249,
     112,   997,   998,   112,   249,   244,   249,   112,   249,   181,
     112,  1007,   298,  1009,   112,   120,    76,  1013,   175,   181,
     311,   193,   249,   125,   249,   249,  1022,   112,   212,   249,
    1026,   193,  1028,   175,  1030,   181,  1032,   181,  1034,   175,
    1036,   249,  2754,   249,    76,   249,  1042,   193,   181,   134,
     175,   195,   112,  1049,   345,   346,   249,    76,  1054,   244,
     244,    76,  1058,  1059,  1060,  1061,  1062,  1063,  1064,   181,
      76,  1067,   181,   411,   134,   181,   181,   249,   416,   181,
     112,   193,   244,   181,   193,   371,    76,   249,   193,    76,
    3254,   193,   249,   112,   157,   193,   181,   112,   249,    76,
     244,   249,   134,   249,  1100,   249,   112,   249,   193,   400,
    1106,   244,   249,   249,   125,   134,  1112,  1113,   181,  1115,
    1116,   181,   112,    76,   249,   112,  1122,    76,   157,   242,
     242,   249,    76,   193,   420,   112,   249,   249,   244,   181,
     249,   230,   231,   249,   249,   483,   244,   249,  1144,   181,
     249,   249,   181,   249,  1150,   249,  1152,  1153,   249,   112,
     249,   193,   181,   112,   249,   503,   181,  2787,   112,  1165,
     181,  1167,    76,   242,   193,   181,  1172,    76,   193,  1175,
     249,   244,    76,   181,   195,    76,    76,   193,  1184,   249,
     181,   181,  1188,  3357,   181,  1191,   249,  1193,  1194,  1195,
     249,   249,   244,   193,   181,    76,   193,   249,   112,   249,
    1206,  1207,  1208,   112,   249,   244,   193,   249,   112,  1215,
    1216,   112,   112,  1219,  1220,  1221,  1222,   249,   181,   181,
     249,   249,   181,   244,   249,    76,   249,   181,   249,   191,
     193,   112,   249,   249,   193,   181,   244,   285,   286,   193,
     249,   249,   249,   244,    38,   181,   181,   249,   249,   249,
     997,   998,   249,  1259,  1260,  1261,  1262,  1263,   249,   181,
     181,   112,   249,   249,   249,  1601,   181,   181,  1274,   181,
    1276,   249,   181,  1279,   249,   181,   249,   181,   249,   193,
     181,   181,   244,   249,   193,    76,   249,  2481,    82,   193,
     249,   249,   193,   193,   249,   249,   231,  1303,   244,  1305,
     181,  1307,   223,   249,  1310,   227,   100,   101,   244,   244,
    1316,   249,   193,   249,   249,   249,  1322,   249,   181,   249,
     249,   112,   244,   244,   249,  1331,  1332,  1333,  1334,   244,
     181,   249,   244,   249,   249,   249,   249,   249,   244,   249,
     249,   181,   193,   249,    76,   249,    76,   249,   249,   249,
     249,   249,   249,  1359,  1360,  1361,  2986,  2987,   249,  1365,
      76,    76,   249,  1369,    76,    76,   249,   249,   249,   249,
      76,   249,  1378,   249,  1380,  1381,  1382,  1383,   249,   249,
     112,   244,   112,  1389,   249,  1391,   249,  1393,  1394,   249,
     181,   249,   249,   249,   249,   249,   112,   112,   249,  1405,
     112,   112,   193,   249,   244,  1411,   112,  1413,  1414,   249,
    1416,  1417,  1418,  1419,  1420,  1421,  1422,  1423,  1424,  1425,
    1426,  1427,  1428,  1429,  1430,  1431,  1432,  1433,  1434,  1435,
    1436,  1437,  1438,  1439,  1440,  1441,  1442,  1443,  1444,  1445,
    1446,  1447,  1448,  1449,  1450,  1451,  1452,  1453,  1454,   181,
      76,   181,   249,   249,   249,    76,    76,  1463,   249,    76,
     181,   193,   249,   193,   249,   181,   181,  1803,   249,   181,
     181,   249,  1478,  1479,  1480,   181,    76,   193,   193,   181,
     249,   193,   193,   249,    76,   249,   112,   193,   249,   181,
       0,   112,   112,    76,    76,   112,  1318,   157,  1320,   219,
    1506,  1323,   181,   125,    76,  3135,   204,   249,     7,    76,
     223,     7,   112,   223,    76,  1337,   223,   249,   249,   249,
     112,  1527,   249,   244,   249,   249,   249,   249,   249,   112,
     112,  1537,  1538,   249,   249,    76,   101,   249,   249,    76,
     112,   249,   244,   249,   249,   112,  1552,   249,  1554,   125,
     112,  1557,   244,   249,   249,   181,  1562,   249,  1564,    76,
     181,   181,   249,   249,   181,   244,  1572,   193,    76,  1575,
     249,   112,   193,   193,  1580,   112,   193,  1583,   181,   223,
    1586,   181,  1588,  1589,  1590,  1591,   249,  1593,   249,   181,
     249,  1597,   181,   193,   249,   112,   249,   249,   181,   181,
     181,   193,   249,   249,   112,   323,  1612,  1613,   249,   181,
     193,   193,    76,  1619,   181,   181,  1622,   249,  3248,   181,
    1626,   193,    76,   249,  1630,   181,   193,  1633,   249,   249,
    1636,   193,   249,    76,   249,   249,  1642,   249,   249,  1645,
     181,   244,  1648,   181,   181,   249,   249,  1653,   112,   249,
     249,   181,   193,   181,   249,   244,   193,   249,   112,  1665,
     249,    76,    76,   244,   181,   249,   249,   249,   249,   112,
    2111,   249,   249,   181,   249,   249,   193,   249,   244,    76,
     249,   249,   249,   249,   181,   193,   249,   249,   244,   249,
     381,    59,   249,   249,    76,   249,  1702,   112,   112,    76,
     164,  1707,  1708,   249,  1710,  1711,   244,  1713,   249,   249,
      59,   249,   249,    76,   244,   112,   244,   181,   249,   249,
      76,   249,  3352,    59,    59,    59,    59,   181,  1734,   193,
     112,    59,   249,    76,    76,   112,  1742,   249,   181,   193,
    1746,   249,  1748,  1749,   378,   249,  1752,   244,   249,   112,
     193,  1757,   249,   372,   249,  1761,   112,   330,    76,  1765,
    1766,  1767,   372,   327,   125,    76,   181,   181,   232,   112,
     112,    76,    76,  1779,  1780,  1781,   249,   125,   193,   193,
    1786,  1787,  1788,  1789,   181,  1791,   125,    76,  1794,  1795,
    1796,    76,  1798,  1799,   112,   249,   193,    76,  1804,   181,
    1806,   112,   125,   181,   181,   125,   249,   112,   112,   249,
      76,   193,   125,    76,   249,   249,   193,    76,   181,  1825,
    1826,  1827,  1828,   112,    76,   181,  1832,   112,   249,   249,
     193,   249,   249,   112,   249,   249,  1842,   193,   181,   181,
     249,  1847,   125,  1849,    76,    76,   112,   125,    76,   112,
     193,   193,   249,   112,   125,   125,   181,   372,    76,   249,
     112,   249,    76,   181,    76,   249,   244,   249,   249,   181,
     181,   249,   249,   125,   372,   193,   181,   181,  1884,   372,
     112,   112,   193,   181,   112,   372,   249,   125,   193,   193,
      72,   125,   181,   249,   112,   249,   181,   249,   112,   249,
     112,   125,   181,   117,   193,  1911,   249,   249,   193,   249,
      95,   249,   249,   380,   193,   181,   125,  1923,   181,   244,
     249,   372,   181,   249,   249,    59,   111,   193,  2264,   181,
     193,   249,   244,   249,   193,  1941,  1942,   249,   249,   249,
      76,   193,    76,   249,   249,   249,   244,   249,   249,   181,
     181,   249,   249,   181,  1960,   249,  1962,   249,   249,   249,
     249,   193,   193,   181,   249,   193,  1972,   181,   249,   181,
     249,  1977,   249,   249,   249,   193,   112,  1983,   112,   193,
    1986,   193,   249,   249,   249,  1991,   249,   181,   249,  1995,
     249,   513,  1998,  1999,  2000,  2001,  2002,  2003,  2004,  2005,
    2006,  2007,  2008,  2009,  2010,  2011,  2012,  2013,  2014,  2015,
    2016,  2017,  2018,  2019,  2020,    76,   244,   249,   249,  2025,
    2026,  2027,  2028,  2029,  2030,  2031,  2032,  2033,  2034,   125,
     244,   249,   125,   249,   219,   220,   249,   249,    76,  2045,
     372,    76,   372,    76,   125,   181,    59,   181,    59,    59,
     244,   112,   510,   372,  2060,   249,   125,   193,   243,   193,
    2066,   125,   181,    76,    76,    76,    76,    76,  2074,  2075,
     125,   125,    76,  2079,   112,   125,    59,   112,   249,   112,
    2086,  2087,  2088,   249,   372,  2091,  2092,   249,   372,  2095,
     372,   249,   372,    76,   249,    76,   249,  2103,   372,   112,
     112,   112,   112,   112,   372,   117,   125,   249,   112,   249,
     372,   249,  2118,   249,  2120,  2121,   372,   372,    76,  2125,
     181,  2127,   181,   372,  2130,   244,    76,  2133,   372,   112,
     249,   112,   193,  2139,  2140,   320,   321,  2143,   249,   181,
    2146,   249,   372,   181,   249,  2151,   181,  2153,   181,  2155,
      76,   372,   337,   372,   112,   193,    76,    76,   193,    76,
     193,  2497,   112,   249,   372,   372,   372,   372,   181,   181,
     181,   181,   181,   372,  2180,  2181,  2182,   181,  2184,  2185,
     193,   193,   193,   193,   193,   244,   112,   223,   249,   193,
     249,    76,   112,   112,   372,   112,  2202,   249,   181,   372,
     181,   249,   244,  2209,   249,   249,  2212,   249,   249,  2215,
     193,   249,   193,   372,   249,  2221,   249,    76,   403,   404,
     405,   406,   249,   181,   249,  2231,  2232,   112,  2234,  2235,
    2236,   181,   244,  2239,  2240,   193,  2242,  2243,  2244,  2245,
     249,  2247,  2248,   193,  2250,   249,  2252,  2253,  2254,  2255,
     181,  2257,  2258,   112,   249,   181,    76,   249,  2264,   249,
    2266,   181,   181,  2269,   181,   249,  2272,   193,   249,  2275,
    2276,  2277,   181,   193,   193,  2281,   193,   249,  2284,    76,
     181,   249,  2288,   249,   249,   249,  2292,   249,   249,    76,
    2296,   249,   112,   249,  2630,   181,   181,   249,  2304,   249,
     809,   810,   811,   812,   249,   249,   249,   249,   193,   818,
     249,   249,   249,   244,   249,   112,   249,   181,   249,   504,
     249,   506,   181,   249,   509,   112,   511,   249,   148,   249,
     249,   249,   249,   249,   193,   244,   181,   249,   249,   125,
     249,   249,   249,   244,   529,   249,   531,   532,   249,   249,
     249,   249,   249,   249,   249,   540,    76,   542,   244,   249,
    2696,   181,   871,   249,   249,   550,   551,   552,   249,   554,
     555,   556,   249,   193,    76,   560,   561,   562,   563,   249,
     244,   249,   567,   568,   181,   249,   571,   572,   573,   574,
     249,   576,   112,   125,   181,  2401,   193,  2403,   125,   244,
    2406,   586,   587,   588,   249,   249,   193,   249,  2414,  2415,
     112,  2417,  2418,   249,  2420,  2421,  2422,  2423,  2424,  2425,
    2426,  2427,  2428,  2429,  2430,  2431,  2432,   181,  2434,  2435,
    2436,  2437,  2438,  2439,  2440,  2441,  2442,   249,   181,   249,
      76,  2447,  2448,  2449,  2450,  2451,  2452,  2453,  2454,  2455,
    2456,   249,   249,   249,  2460,  2461,   249,   249,    76,   249,
      76,   181,   249,    76,   249,  2471,   249,   249,   249,    76,
    2476,   249,   249,   193,   249,   249,   112,  2483,   181,   181,
    2486,   117,   249,  2489,   249,   249,   249,   181,   249,    76,
     244,   193,   249,   249,   112,   249,   112,   249,   249,   112,
    2506,   244,  2508,   249,   249,   112,   249,   249,  2514,  2515,
    2516,  2517,  2518,    76,  2520,   181,   249,  2523,  2524,   181,
    2526,  2527,   372,  2529,  2530,   112,  2862,  2533,   249,   249,
      76,  2537,  2538,    76,  2540,  2541,   249,  2543,  2544,   249,
     249,   244,   181,  2549,  2550,   181,   249,   249,    76,   112,
     244,   249,  2558,  2559,   372,   249,  2562,   193,  2564,  2565,
     249,  2567,  2568,   181,    76,   181,   112,   372,   181,   112,
     249,    76,    76,   249,   181,   193,  2582,   193,   244,   249,
     193,   249,   244,   249,   112,  2591,   193,   249,   249,   249,
     249,   181,   181,   181,   181,   181,  2602,  2603,  2604,   372,
     112,   184,  2608,  2609,  2610,   244,   193,   112,   112,  2615,
     249,  2617,  2618,  2619,   249,  2621,  2622,  2623,   181,  2625,
    2626,   249,   181,   249,  2630,   249,   249,   125,  2634,    94,
     193,  2637,  2638,   249,    94,   181,   249,    76,   181,   125,
    2646,  2647,   249,    76,    76,   249,    76,   193,    76,   372,
     193,  2657,   125,   181,   244,   244,   244,    76,   244,   249,
     249,   249,   249,   249,   125,   193,   323,   249,   249,   181,
     249,   249,   249,   112,   249,   249,   181,   181,   249,   112,
     112,   193,   112,   249,   112,   244,   249,   181,   193,   193,
     249,   249,  2698,   112,  2700,  2701,  2702,  2703,  2704,  2705,
    2706,  2707,  2708,   249,  2710,  2711,   249,  2713,  2714,  2715,
    2716,  2717,  2718,  2719,  2720,   249,  2722,  2723,  2724,  2725,
    2726,   249,   249,  2729,  2730,  2731,  2732,  2733,  2734,  2735,
    2736,  2737,  2738,  2739,  2740,  2741,  2742,   249,  2744,    76,
     249,  2747,   181,   249,   249,   249,   181,    59,   181,   181,
     244,   181,   181,   181,   193,   249,    76,   249,   249,   249,
     193,   193,   181,   193,    76,   193,   181,   249,  2774,  2775,
    2776,   249,  2778,   249,   193,   112,  2782,   249,  2784,   249,
    2786,   249,  2788,   249,  2790,   249,  2792,   181,  2794,   125,
    2796,   249,   112,   249,   249,   249,   249,   249,   249,   249,
     112,  2807,   181,  2809,  2810,   249,  2812,  2813,   181,   244,
     249,   181,   181,   249,   249,   244,   249,   249,  2824,   249,
     249,   249,   249,   181,  2830,  2831,  2832,  2833,    76,   244,
     249,  2837,  2838,  2839,   249,  3171,  2842,    59,   249,   249,
    2846,   249,   249,   125,   181,  2851,   249,  2853,  2854,  2855,
     244,  2857,  2858,   249,    76,   249,   193,   249,    76,    76,
      76,   181,  2868,  2869,   112,   244,   362,   249,  2874,   181,
     249,   244,   249,   193,   244,   244,   249,   249,   125,   249,
     249,   193,   249,   249,   249,    76,   244,   249,   249,   249,
     112,   249,   249,   181,   112,   112,   112,   249,   249,  2905,
     249,  2907,  2908,  2909,  2910,  2911,  2912,  2913,   249,  2915,
    2916,  2917,   249,  2919,  2920,  2921,   249,  2923,   181,  2925,
    2926,   112,  2928,  2929,  2930,  2931,   181,  2933,  2934,   249,
    2936,  2937,  2938,   181,  2940,  2941,  2942,  2943,  2944,  2945,
    2946,  2947,  2948,  2949,  2950,   193,  2952,  2953,  2954,  2955,
    2956,  2957,  2958,   249,   181,  3291,   244,   249,  2964,   181,
     249,   249,   249,   181,   181,   181,   249,  2973,  2974,   249,
     249,   193,   249,   249,   249,   193,   193,   193,   181,  2985,
     181,   244,  2988,   249,   181,   249,   249,   181,   181,   244,
     181,   249,   249,   249,   249,  3001,  3002,  3003,   249,  3005,
      76,   249,   193,   249,   249,   181,   249,   249,   249,   249,
     249,   249,  3018,  3019,  3020,   249,  3022,   244,    76,   249,
    3026,   249,   249,   181,   249,   249,  3032,  3033,  3034,  3035,
    3036,   249,   249,   249,  3040,  3041,   112,  3043,  3044,   181,
     249,   244,   249,   244,   249,   249,   249,   244,   249,   249,
     244,   244,   249,   181,   112,   249,   249,   249,   249,  3065,
    3066,  3067,  3068,  3069,  3070,   249,  3072,   249,   244,  3075,
     249,   249,  3078,   249,  3080,  3081,  3082,  3083,  3084,  3085,
    3086,   249,  3088,  3089,   249,  3091,   244,   249,  3094,  3095,
    3096,   249,  3098,  3099,  3100,  3101,  3102,  3103,  3104,   249,
    3106,  3107,   244,  3109,  3110,   181,  3112,   249,  3114,  3115,
    3116,  3117,  3118,  3119,  3120,   249,   244,   193,   249,  3125,
    3126,   249,   249,   181,   249,   249,  3132,  3133,  3134,   249,
     181,   249,   249,   249,  3140,   193,   249,  3143,  3144,  3145,
     249,   125,   510,    76,   125,   249,  3152,  3153,   249,  3155,
    3156,  3157,   249,   249,  3160,   249,   249,   249,   249,  3165,
    3166,  3167,  3168,  3169,   249,   249,  3172,  3173,  3174,   249,
     372,   249,   249,   249,   181,   181,   249,    76,   249,   112,
    3186,  3187,   249,  3189,  3190,   249,  3192,   249,  3194,   181,
    3196,   249,  3198,   244,  3200,   249,  3202,   249,   249,  3205,
    3206,  3207,  3208,  3209,  3210,  3211,  3212,  3213,   249,  3215,
     249,  3217,  3218,   112,   249,   372,  3222,   181,  3224,  3225,
    3226,  3227,  3228,  3229,  3230,   249,  3232,   249,  3234,   249,
    3236,  3237,  3238,  3239,  3240,  3241,  3242,   244,   244,   249,
     249,  3247,   249,   249,  3250,  3251,   249,  3253,   181,   249,
    3256,  3257,   244,  3259,   249,  3261,  3262,   249,   249,  3265,
     193,  3267,   249,    76,   327,  3271,  3272,   327,   125,    76,
    3276,  3277,  3278,   249,  3280,   249,   249,  3283,  3284,  3285,
     244,   181,   181,   181,   249,   249,   249,  3293,   249,  3295,
      76,  3297,   249,  3299,   193,  3301,   181,  3303,   181,   112,
    3306,  3307,  3308,   249,  3310,   112,  3312,   249,  3314,  3315,
    3316,   249,  3318,   249,  3320,   181,   249,  3323,   249,  3325,
    3326,  3327,  3328,  3329,  3330,  3331,   112,  3333,  3334,   249,
    3336,  3337,  3338,    76,  3340,   249,  3342,  3343,  3344,  3345,
    3346,   181,   249,    76,   244,  3351,   244,   249,    76,   249,
     249,   249,  3358,   249,  3360,  3361,  3362,   249,   249,   244,
    3366,   244,   249,   249,   249,  3371,   249,   249,   181,   112,
    3376,  3377,  3378,  3379,   181,  3381,  3382,  3383,   244,   112,
     193,   181,   249,   249,   112,  3391,   193,   249,  3394,  3395,
    3396,   249,  3398,   249,  3400,   181,   249,  3403,  3404,   249,
    3406,   249,   249,   181,   244,    76,  3412,   193,  3414,   249,
    3416,  3417,  3418,  3419,  3420,  3421,  3422,  3423,   249,  3425,
     249,  3427,   249,  3429,   249,  3431,   249,  3433,  3434,  3435,
    3436,  3437,   249,   249,   249,  3441,   249,   249,   181,  3445,
      76,   112,   249,  3449,   244,    76,   181,   249,   181,   249,
     193,  3457,   249,   181,  3460,  3461,   249,  3463,  3464,  3465,
     193,  3467,   181,   249,   181,   193,   244,  3473,   249,  3475,
     249,   249,   323,   249,  3480,  3481,   112,   372,   181,   249,
     249,   112,  3488,   249,   249,  3491,   249,  3493,  3494,  3495,
    3496,  3497,  3498,  3499,   249,  3501,   249,  3503,   249,    76,
    3506,  3507,  3508,  3509,  3510,  3511,   249,   249,  3514,   244,
     181,  3517,   249,   181,   249,   249,   249,   249,   249,   249,
    3526,   249,   193,  3529,  3530,   244,   249,   244,   249,   249,
     249,   249,   249,    76,   249,   112,  3542,  3543,   181,  3545,
    3546,   244,  3548,   249,   249,   181,   249,    76,  3554,    76,
     181,  3557,  3558,  3559,  3560,  3561,  3562,   193,   181,   249,
     249,  3567,   193,  3569,  3570,  3571,   249,  3573,  3574,   112,
     249,   249,   249,   181,  3580,  3581,   244,   249,   249,    76,
     249,   249,  3588,   112,  3590,   112,  3592,   249,  3594,   249,
    3596,   249,   249,    76,  3600,   249,  3602,  3603,  3604,  3605,
    3606,   244,   249,   249,   181,   181,   249,   181,  3614,  3615,
     249,  3617,  3618,   249,  3620,   112,   193,   249,   249,   249,
    3626,   244,  3628,   249,    76,   249,   249,   249,  3634,   112,
    3636,  3637,   249,  3639,  3640,  3641,   244,  3643,   181,   249,
     249,   249,  3648,    76,  3650,  3651,  3652,   249,   249,   249,
     193,  3657,   181,  3659,   181,  3661,   181,  3663,   249,   249,
     112,  3667,   249,  3669,   193,  3671,   193,  3673,   244,  3675,
     244,  3677,   249,   249,   181,   249,  3682,  3683,   249,   112,
    3686,   249,  3688,  3689,   181,  3691,  3692,  3693,   249,   249,
    3696,   249,  3698,  3699,   181,  3701,   193,  3703,   181,   249,
     249,  3707,    76,    76,    76,  3711,   249,  3713,   249,  3715,
     193,  3717,   249,  3719,  3720,  3721,  3722,    76,   249,   244,
     249,  3727,   249,   249,   249,  3731,   249,  3733,   249,   181,
     125,  3737,   249,   249,   249,   249,  3742,   244,   112,   112,
     112,   193,   249,  3749,  3750,   249,  3752,   249,   181,  3755,
     249,  3757,   249,   112,  3760,  3761,   249,   244,   249,  3765,
     193,  3767,   249,  3769,   249,  3771,   249,     4,     5,   249,
     249,  3777,  3778,  3779,    11,    12,   249,   249,    15,   280,
      17,   282,   283,   284,   285,   286,   287,   288,   289,   249,
     249,   249,   249,   249,   249,    32,    33,   249,   249,   249,
     249,   249,   249,    40,    41,    42,    43,   181,   181,   181,
     249,   249,   249,   249,   125,    52,   249,    54,   249,   193,
     193,   193,   181,   249,   249,   249,   249,   249,   249,   249,
     249,    68,   249,    70,   193,   249,    73,   249,   249,   249,
     249,    78,   249,   249,    81,   249,   249,   249,    85,    86,
      87,   249,   249,   249,    91,   249,   249,    94,    95,   249,
     249,    98,    99,   249,   249,   102,   249,   104,   105,   249,
     107,   108,   249,   249,   249,   249,   249,   249,   249,   249,
     117,   249,   249,   120,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   144,   145,   249,
     249,   148,   149,   150,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   172,
     595,   597,   360,   625,   935,   304,  1047,  1051,   304,  1014,
     177,   359,   377,   407,   476,   182,   415,   963,   946,   186,
     438,   188,   189,   190,   494,   528,  1224,  1562,  1011,    -1,
      -1,   198,   304,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   319,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   226,
      -1,   228,   229,    -1,    -1,   232,   233,   234,   235,    -1,
     237,    -1,   239,   240,   241,   242,    -1,    -1,   245,   246,
     247,    -1,    -1,   250,    -1,   252,   253,   254,    -1,    -1,
      -1,   258,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   273,   274,   275,    -1,
      -1,   278,   279,   280,   281,   282,    -1,   284,   285,   286,
     287,   288,   289,   290,    -1,    -1,   293,   294,   295,    -1,
      -1,   298,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   312,   313,   314,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   323,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   339,    -1,    -1,   342,   343,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   356,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     377,   378,   379,    -1,    -1,   382,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   396,
     397,    -1,   399,    -1,    -1,    -1,   403,    -1,    -1,   406,
     407,    -1,    -1,   410,   411,    -1,    -1,    -1,   415,   416,
      -1,    -1,   419,    -1,    -1,   422,    -1,    -1,    -1,   426,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   438,    -1,    -1,    -1,    -1,    -1,    -1,   445,    -1,
     447,    -1,    -1,    -1,    -1,   452,    -1,    -1,    -1,   456,
      -1,    -1,    -1,   460,    -1,    -1,   463,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   474,   475,   476,
      -1,    -1,    -1,    -1,    -1,    -1,   483,    -1,    -1,    -1,
      -1,    -1,   489,    -1,    -1,    -1,    -1,   494,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   503,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   520,   521,    -1,    -1,    -1,    -1,    -1,
     527,   528,   529,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   545,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   573,    -1,   575,    -1,
      -1,    -1,    -1,    -1,   581,   582,   583,   584,    -1,    -1,
      -1,   588,   589,   590,   591,   592,    -1,   594,   595,   596,
     597,    -1,    -1,    -1,    -1,    -1,   603,    -1,   605,   606,
      -1,   608,   609,   610,    -1,    -1,    -1,    -1,    -1,   616,
     617,   618,    -1,    -1,    -1,    -1,    -1,   624,    -1,    -1,
      -1,    -1,    -1,   630,    -1,    -1,    -1,   634,    -1,    -1,
      -1,    -1,   639,    -1,    -1,    -1,    -1,    -1,   645,    -1,
     647,   648,   649,    -1,    -1,   652,    -1,   654,    -1,   656,
      -1,   658,    -1,   660,    -1,    -1,    -1,   664,    -1,   666,
      -1,   668,    -1,   670,    -1,   672,    -1,    -1,    -1,    -1,
      -1,    -1,   679,    -1,    -1,    -1,   683,    -1,   685,    -1,
     687,    -1,    -1,    -1,    -1,   692,   693,   694,   695,    -1,
     697,   698,    -1,    -1,   701,    -1,    -1,   704,    -1,    -1,
      -1,    -1,    -1,    -1,   711,    -1,    -1,   714,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   722,    -1,    -1,    -1,    -1,
     727,    -1,    -1,    -1,    -1,    -1,    -1,   734,   735,    -1,
      -1,    -1,   739,    -1,    -1,    -1,    -1,    -1,   745,    -1,
     747,   748,   749,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     757,   758,    -1,    -1,    -1,    -1,   763,   764,    -1,   766,
     767,   768,   769,    -1,   771,    -1,   773,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   782,   783,    -1,    -1,   786,
     787,    -1,    -1,    -1,   791,   792,    -1,   794,    -1,   796,
      -1,    -1,   799,   800,   801,   802,    -1,   804,   805,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   814,    -1,   816,
      -1,    -1,    -1,   820,   821,    -1,   823,    -1,    -1,    -1,
      -1,    -1,   829,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   838,   839,    -1,    -1,    -1,    -1,   844,   845,    -1,
      -1,    -1,   849,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   861,   862,   863,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   874,    -1,   876,
      -1,   878,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   889,    -1,   891,    -1,   893,   894,    -1,    -1,
      -1,    -1,   899,    -1,    -1,    -1,    -1,   904,    -1,    -1,
      -1,    -1,    -1,   910,    -1,    -1,    -1,    -1,    -1,   916,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   925,    -1,
      -1,    -1,    -1,   930,    -1,   932,   933,   934,   935,    -1,
     937,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   946,
     947,    -1,   949,    -1,    -1,   952,   953,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   962,   963,   964,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   980,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   990,    -1,    -1,   993,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1006,
      -1,  1008,    -1,  1010,    -1,    -1,    -1,  1014,    -1,  1016,
      -1,    -1,    -1,    -1,  1021,    -1,    -1,  1024,  1025,    -1,
    1027,    -1,  1029,    -1,  1031,    -1,  1033,    -1,  1035,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1047,    -1,    -1,    -1,  1051,    -1,  1053,    -1,    -1,  1056,
    1057,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1068,    -1,  1070,   115,   116,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   125,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1097,  1098,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1109,    -1,  1111,    -1,    -1,    -1,  1115,    -1,
      -1,    -1,    -1,  1120,  1121,    -1,  1123,  1124,  1125,  1126,
      -1,    -1,    -1,    -1,    -1,    -1,  1133,    -1,  1135,    -1,
    1137,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1166,
      -1,    -1,  1169,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1186,
    1187,    -1,  1189,    -1,    -1,    -1,    -1,    -1,    -1,  1196,
    1197,  1198,  1199,  1200,  1201,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1213,    -1,    -1,    -1,
    1217,  1218,    -1,    -1,    -1,    -1,  1223,  1224,  1225,  1226,
    1227,  1228,    -1,  1230,    -1,  1232,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1264,    -1,    -1,
      -1,  1268,    -1,    -1,    -1,    -1,  1273,    -1,  1275,    -1,
      -1,  1278,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1286,
      -1,  1288,    -1,  1290,    -1,    -1,    -1,  1294,    -1,    -1,
     341,   342,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1327,    -1,    -1,  1330,    -1,    -1,    -1,    -1,    -1,  1336,
      -1,  1338,  1339,  1340,  1341,  1342,  1343,    -1,  1345,  1346,
    1347,  1348,  1349,  1350,    -1,  1352,  1353,    -1,  1355,  1356,
    1357,  1358,    -1,    -1,    -1,  1362,    -1,  1364,    -1,  1366,
    1367,  1368,    -1,  1370,  1371,  1372,  1373,  1374,  1375,    -1,
    1377,    -1,  1379,    -1,    -1,    -1,    -1,  1384,  1385,  1386,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1395,    -1,
      -1,    -1,    -1,  1400,    -1,    -1,    -1,    -1,    -1,  1406,
      -1,    -1,   453,   454,   455,   456,    -1,   458,    -1,   460,
     461,   462,   463,   464,   465,   466,   467,    -1,   469,   470,
     471,   472,   473,   474,   475,   476,   477,   478,   479,   480,
     481,   482,    -1,   484,   485,   486,   487,   488,   489,   490,
     491,   492,   493,    -1,    -1,    -1,    -1,    -1,  1455,    -1,
      -1,    -1,    -1,    -1,  1461,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    95,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1478,    -1,    -1,    -1,    -1,  1483,    -1,   111,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1499,    -1,    -1,  1502,    -1,    -1,    -1,  1506,
      -1,    -1,    -1,    -1,  1511,  1512,    -1,    -1,    -1,    -1,
      -1,    -1,  1519,  1520,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1530,    -1,    -1,    -1,  1534,  1535,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1545,    -1,
      -1,    -1,  1549,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1567,    -1,  1569,    -1,    -1,    -1,  1573,    -1,  1575,    -1,
      -1,  1578,    -1,  1580,    -1,  1582,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   219,   220,  1595,    -1,
    1597,    -1,  1599,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1609,  1610,    -1,    -1,    -1,    -1,    -1,    -1,
     243,  1618,    -1,    -1,  1621,    -1,    -1,    -1,    -1,  1626,
      -1,    -1,  1629,    -1,    -1,  1632,    -1,    -1,  1635,    -1,
      -1,  1638,    -1,    -1,  1641,    -1,    -1,  1644,    -1,    -1,
    1647,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1661,    -1,  1663,  1664,    -1,    -1,
      -1,  1668,  1669,  1670,  1671,  1672,    -1,    -1,  1675,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   320,   321,    -1,
      -1,    -1,    -1,  1700,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   337,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1719,    -1,    -1,    -1,  1723,    -1,  1725,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1754,    -1,    -1,
      -1,  1758,  1759,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1768,    -1,  1770,  1771,  1772,  1773,    -1,  1775,  1776,
     403,   404,   405,   406,    -1,    -1,    -1,    -1,  1785,    -1,
      -1,    -1,    -1,    -1,    -1,  1792,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1800,    -1,  1802,  1803,    -1,  1805,    -1,
      -1,    -1,  1809,  1810,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1833,    -1,    -1,    -1,
      -1,    -1,    -1,  1840,    -1,    -1,  1843,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1860,  1861,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1873,  1874,    -1,  1876,
      -1,   504,    -1,   506,    -1,    -1,   509,    -1,   511,    -1,
    1887,    -1,    -1,    -1,    -1,  1892,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   529,    -1,   531,   532,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   540,  1915,   542,
      -1,    -1,    -1,    -1,    -1,  1922,    -1,    -1,    -1,    -1,
      -1,   554,   555,   556,    -1,    -1,    -1,   560,   561,   562,
     563,    -1,    -1,    -1,   567,   568,  1943,    -1,   571,   572,
     573,   574,    -1,   576,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   586,   587,   588,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1976,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1988,    -1,    -1,    -1,    -1,  1993,  1994,    -1,  1996,
    1997,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2021,  2022,  2023,  2024,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2043,    -1,    -1,  2046,
      -1,    -1,    -1,    -1,  2051,  2052,    -1,    -1,  2055,  2056,
      -1,    -1,    -1,    -1,    -1,    -1,  2063,  2064,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2088,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2105,    -1,
      -1,    -1,    -1,    -1,  2111,    -1,    -1,    -1,    -1,  2116,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2125,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2157,    -1,    -1,    -1,    -1,    -1,    -1,  2164,    -1,    -1,
      -1,    -1,    -1,    -1,  2171,    -1,    -1,    -1,    -1,    -1,
    2177,  2178,  2179,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2191,    -1,    -1,    -1,    -1,  2196,
      -1,  2198,  2199,    -1,  2201,    -1,    -1,    -1,    -1,    -1,
      -1,  2208,    -1,    -1,  2211,    -1,    -1,  2214,    -1,    -1,
      -1,  2218,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2238,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2251,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2274,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2285,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2298,  2299,    -1,    -1,    -1,    -1,    -1,    -1,  2306,
      -1,    -1,    -1,    -1,    -1,    -1,  2313,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2335,    -1,
    2337,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2397,  2398,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2408,  2409,    -1,    -1,  2412,  2413,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2443,  2444,  2445,  2446,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2459,    -1,    -1,    -1,    -1,  2464,    -1,    -1,
      -1,    -1,  2469,  2470,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2502,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2561,    -1,  2563,    -1,    -1,    -1,
      -1,    -1,    -1,  2570,    -1,    -1,    -1,  2574,    -1,  2576,
      -1,    -1,    -1,  2580,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2588,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2605,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    16,    -1,    -1,    -1,    -1,  2636,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    29,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    43,  2659,    -1,  2661,    -1,    48,  2664,    -1,  2666,
      -1,    -1,    -1,    -1,    -1,    57,    -1,  2674,    60,    61,
      -1,  2678,    -1,    65,    66,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2691,    -1,  2693,    -1,    -1,    -1,
    2697,    -1,  2699,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2727,  2728,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   121,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2746,
      -1,  2748,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2756,
     142,   143,   144,   145,    -1,    -1,    -1,    -1,    -1,   151,
      -1,    -1,    -1,   155,    -1,    -1,   158,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   180,    -1,
     182,    -1,    -1,    -1,    -1,    -1,    -1,  2804,  2805,    -1,
      -1,    -1,    -1,    -1,   196,    -1,    -1,    -1,  2815,    -1,
    2817,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   210,   211,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   219,    -1,   221,
      -1,    -1,   224,    -1,   226,    -1,   228,   229,    -1,    -1,
      -1,   233,    -1,   235,   236,   237,   238,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2862,    -1,    -1,  2865,    -1,
      -1,   253,    -1,   255,    -1,    -1,    -1,    -1,    -1,   261,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2884,    -1,  2886,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2904,   290,    -1,
      -1,    -1,    -1,   295,   296,   297,    -1,   299,    -1,   301,
      -1,   303,    -1,    -1,    -1,    -1,   308,   309,   310,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   328,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   337,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2962,    -1,   349,   350,   351,
     352,   353,   354,   355,    -1,   357,   358,   359,   360,    -1,
      -1,    -1,    -1,    -1,    -1,   367,   368,    -1,   370,    -1,
      -1,    -1,   374,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2997,  2998,    -1,    -1,    -1,   387,    -1,   389,    -1,    -1,
      -1,    -1,   394,    -1,    -1,    -1,    -1,   399,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   419,    -1,    -1,
      -1,   423,   424,   425,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3052,    -1,  3054,    -1,    -1,
      -1,   443,    -1,    -1,   446,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   515,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3139,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3179,    -1,  3181,    -1,    -1,  3184,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3252,    -1,  3254,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3275,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  3288,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3356,
    3357,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  3402,    -1,    -1,    -1,    -1,
      -1,  3408,    -1,  3410,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  3444,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  3479,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  3523,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  3534,    -1,  3536,
      -1,    -1,    -1,  3540,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  3556,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  3578,    -1,    -1,     3,    -1,     5,     6,    -1,    -1,
    3587,    -1,    11,    -1,    13,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    21,    -1,    23,    24,    25,    26,    27,    -1,
      -1,    -1,    31,  3610,    -1,  3612,    -1,    -1,    -1,    -1,
      -1,    40,    41,    42,    -1,    -1,    -1,    -1,  3625,    -1,
      -1,    50,    -1,    52,    53,    54,  3633,    56,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    64,    -1,    -1,    -1,    -1,
      -1,    70,    71,    -1,    -1,    74,    -1,    -1,    -1,    78,
      79,    80,    81,    -1,    -1,    -1,    85,    -1,    87,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      99,    -1,    -1,   102,   103,   104,   105,   106,   107,    -1,
     109,   110,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   118,
      -1,    -1,    -1,   122,   123,    -1,    -1,    -1,   127,   128,
      -1,    -1,   131,   132,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   141,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     149,    -1,    -1,   152,   153,    -1,    -1,   156,    -1,    -1,
     159,  3738,    -1,   162,   163,   164,   165,   166,   167,   168,
     169,   170,   171,   172,   173,    -1,    -1,    -1,    -1,   178,
     179,    -1,    -1,    -1,   183,  3762,   185,    -1,   187,   188,
     189,   190,    -1,   192,    -1,   194,    -1,  3774,   197,  3776,
      -1,    -1,    -1,   202,    -1,    -1,   205,   206,   207,   208,
     209,    -1,    -1,   212,    -1,    -1,   215,   216,   217,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   232,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   241,    -1,    -1,    -1,    -1,    -1,   247,    -1,
      -1,    -1,   251,    -1,    -1,    -1,    -1,   256,    -1,   258,
     259,    -1,    -1,    -1,    -1,    -1,   265,    -1,   267,    -1,
      -1,   270,   271,    -1,    -1,   274,    -1,    -1,   277,   278,
     279,    -1,   281,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   292,   293,   294,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   306,   307,    -1,
      -1,    -1,    -1,   312,   313,   314,    -1,   316,    -1,    -1,
      -1,    -1,    -1,   322,    -1,   324,   325,    -1,    -1,    -1,
     329,    -1,    -1,    -1,    -1,    -1,    -1,   336,    -1,    -1,
     339,   340,    -1,    -1,    -1,   344,    -1,    -1,   347,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   361,    -1,    -1,   364,    -1,    -1,    -1,    -1,
     369,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   384,    -1,   386,    -1,   388,
      -1,    -1,   391,   392,   393,    -1,   395,   396,    -1,   398,
      -1,    -1,   401,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   411,    -1,    -1,    -1,    -1,   416,   417,   418,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   426,   427,   428,
     429,   430,   431,    -1,   433,    -1,    -1,    -1,    -1,    -1,
     439,   440,    -1,    -1,    -1,    -1,    -1,    -1,   447,    -1,
      -1,    -1,   451,   452,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   494,   495,    -1,    -1,    -1,
      -1,    -1,   501,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   514,    -1,   516,   517,   518,
     519,   520,   521,   522,   523,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   534,   535,    -1,    -1,   538,
     539,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   553,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   564,   565,   566,    -1,    -1,
     569,   570,    -1,     3,    -1,     5,     6,    -1,    -1,    -1,
      -1,    11,   581,    13,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    21,    -1,    23,    24,    25,    26,    27,    -1,    -1,
      -1,    31,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      40,    41,    42,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      50,    -1,    52,    53,    54,    -1,    56,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    64,    -1,    -1,    -1,    -1,    -1,
      70,    71,    -1,    -1,    74,    -1,    -1,    -1,    78,    79,
      80,    81,    -1,    -1,    -1,    85,    -1,    87,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    99,
      -1,    -1,   102,   103,   104,   105,    -1,   107,    -1,   109,
     110,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   118,    -1,
      -1,    -1,   122,   123,    -1,    -1,    -1,   127,   128,    -1,
      -1,   131,   132,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   141,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   149,
      -1,    -1,   152,   153,    -1,    -1,   156,    -1,    -1,   159,
      -1,    -1,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   172,   173,    -1,    -1,    -1,    -1,   178,   179,
      -1,    -1,    -1,   183,    -1,   185,    -1,   187,   188,   189,
     190,    -1,   192,    -1,   194,    -1,    -1,   197,    -1,    -1,
      -1,    -1,   202,    -1,    -1,   205,   206,   207,   208,   209,
      -1,    -1,   212,    -1,    -1,   215,   216,   217,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   232,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   241,    -1,    -1,    -1,    -1,    -1,   247,    -1,    -1,
      -1,   251,    -1,    -1,    -1,    -1,   256,    -1,   258,   259,
      -1,    -1,    -1,    -1,    -1,   265,    -1,   267,    -1,    -1,
     270,   271,    -1,    -1,   274,    -1,    -1,   277,   278,   279,
      -1,   281,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   292,   293,   294,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   306,   307,    -1,    -1,
      -1,    -1,   312,   313,   314,    -1,   316,    -1,    -1,    -1,
      -1,    -1,   322,    -1,   324,   325,    -1,    -1,    -1,   329,
      -1,    -1,    -1,    -1,    -1,    -1,   336,    -1,    -1,   339,
     340,    -1,    -1,    -1,   344,    -1,    -1,   347,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   361,    -1,    -1,   364,    -1,    -1,    -1,    -1,   369,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   384,    -1,   386,    -1,   388,    -1,
      -1,   391,   392,   393,    -1,   395,   396,    -1,   398,    -1,
      -1,   401,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   411,    -1,    -1,    -1,    -1,   416,   417,   418,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   426,   427,   428,   429,
     430,   431,    -1,   433,    -1,    -1,    -1,    -1,    -1,   439,
     440,    -1,    -1,    -1,    -1,    -1,    -1,   447,    -1,    -1,
      -1,   451,   452,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   494,   495,    -1,    -1,    -1,    -1,
      -1,   501,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   514,    -1,   516,   517,   518,   519,
     520,   521,   522,   523,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   534,   535,    -1,    -1,   538,   539,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   553,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   564,   565,   566,    -1,    -1,   569,
     570,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   581
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
     430,   431,   433,   439,   440,   447,   451,   452,   494,   495,
     501,   514,   516,   517,   518,   519,   520,   521,   522,   523,
     534,   535,   538,   539,   553,   564,   565,   566,   569,   570,
     581,   591,   592,   593,   594,   595,   596,   597,   598,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   622,
     624,   625,   626,   627,   628,   629,   630,   631,   632,   633,
     635,   636,   637,   639,   640,   641,   642,   644,   645,   650,
     653,   655,   656,   657,   658,   659,   660,   661,   662,   663,
     664,   666,   667,   669,   673,   674,   675,   677,   678,   679,
     680,   683,   686,   687,   688,   689,   690,   693,   696,   699,
     701,   703,   705,   709,   710,   711,   712,   713,   714,   715,
     717,   718,   719,   720,   721,   722,   723,   724,   725,   726,
     730,   732,   734,   736,   738,   740,   742,   744,   746,   748,
     749,   753,   757,   759,   762,   763,   764,   765,   766,   767,
     768,   770,   771,   773,   774,   782,   783,   785,   787,   788,
     789,   790,   793,   794,   795,   796,   797,   798,   799,   800,
     804,   807,   809,   811,   813,   816,   819,   821,   823,   825,
     827,   828,   829,   830,   831,   834,   835,   837,   838,   839,
     840,   843,   846,   249,     7,   249,   249,   181,   244,   848,
     249,   848,    76,   112,   181,   193,   849,   849,   849,   125,
     249,   849,   249,   249,   848,   848,   249,   249,   848,   249,
     249,   848,   249,   249,   249,   249,   440,   249,   249,   249,
     249,   249,   249,   249,   249,   249,    38,    82,   100,   101,
     249,   812,   249,   848,   249,   848,   249,   249,   249,   249,
     249,   249,   848,   848,   249,   848,   249,   848,   249,   372,
     249,   249,   249,   249,   849,   249,   249,   249,   848,   249,
     249,   848,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   432,   249,   849,   249,   249,   848,   249,   249,
     848,   249,   249,   848,   249,   249,   249,   849,   125,   249,
     848,   249,   249,   249,   848,   125,   249,   372,   249,   249,
     848,   249,   848,   848,   249,   249,   849,   249,   848,   249,
     249,   125,   195,   249,   848,   125,   195,   249,   848,   249,
     249,   848,   848,   249,   249,   249,   848,   125,   848,   848,
     249,   249,   848,   249,   848,   125,   249,   249,   249,   249,
     249,   249,   249,   249,   848,   849,   249,   249,   848,   249,
     249,   849,   849,   249,   249,   249,   249,   249,   249,   249,
     249,   125,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   440,   249,   848,   249,   848,   249,   249,   848,
     848,   848,   249,   249,   849,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   849,   249,   249,   249,   249,
     249,   372,     0,   106,   593,   157,   848,    75,   135,   136,
     137,   138,   139,   302,   305,   326,   448,   599,   600,   601,
     602,   603,    44,    45,   434,   435,   436,   437,   438,    28,
      86,    88,   140,   413,   422,   502,   503,   504,   507,   508,
     848,    22,   133,   160,   239,   849,   849,   848,    49,   198,
     199,   200,   201,   614,   616,   617,   618,   619,   784,   848,
     623,   848,   223,   848,   849,   125,   378,   778,   848,   362,
     410,   582,   583,   584,   585,   634,   219,   638,   849,    22,
      47,    75,   160,   176,   223,   239,   240,   248,   263,   276,
     373,   383,   385,   432,   643,   647,   668,   672,   160,   223,
     239,   651,   652,    32,    51,   363,   437,   654,   204,   665,
     676,   848,   378,   778,   848,   848,   378,   848,   378,   848,
     378,   848,   848,   694,   848,   700,   848,   702,   848,   704,
     848,   706,   848,   694,   702,   706,     4,   716,   848,   223,
     848,   848,   848,   223,   223,   378,   778,   848,   849,   760,
     848,   776,   848,   848,   848,   848,   769,   848,   848,   413,
     772,   848,     8,    18,    30,    33,    34,    35,    36,    37,
      77,   146,   150,   245,   254,   257,   262,   263,   264,   266,
     269,   273,   319,   343,   348,   483,   496,   497,   775,   848,
     784,   784,   786,   848,   848,   848,   111,   321,   848,    19,
     404,   791,   792,   849,   181,   848,   849,   848,   378,   848,
     848,   213,   214,   268,   317,   390,   397,   415,   578,   579,
     580,    39,    63,    69,    82,    89,   119,   161,   184,   234,
     311,   345,   346,   400,   806,   808,   600,   603,   672,   849,
      12,    83,   101,   124,   147,   203,   211,   219,   250,   252,
     260,   318,   414,   420,   455,   826,   848,   459,   848,   848,
      46,    92,    93,   406,   502,   504,   505,   512,   524,   525,
     526,   527,   528,   530,   540,   541,   543,   545,   546,   550,
     556,   557,   558,   559,   836,   842,    95,   111,   219,   220,
     243,   320,   321,   337,   403,   404,   405,   406,   504,   506,
     509,   511,   529,   531,   532,   540,   542,   554,   555,   556,
     560,   561,   562,   563,   567,   568,   571,   572,   573,   574,
     576,   586,   587,   588,   841,   841,   551,   552,   841,   842,
     844,   845,   841,   536,   537,   847,   727,   778,   848,   249,
       7,     7,   249,   249,   249,   249,   849,   249,   848,   849,
     707,   708,   848,   708,   249,   249,   249,   849,   249,   849,
     125,   848,   646,   849,   638,   849,   849,   101,   849,   812,
     175,   249,   249,   249,   700,   694,   729,   780,   848,   810,
     849,   249,   848,   249,   249,   223,   249,   249,   849,   849,
     758,   781,   848,   758,   249,   691,   692,   848,   249,   702,
     820,   849,   822,   849,   706,   704,   754,   755,   848,   727,
     249,   249,   729,   727,   249,   849,   848,   249,   849,   848,
     249,   849,   223,   750,   751,   848,   849,   249,   152,   249,
     384,   249,   249,   249,   760,   249,   249,   249,   775,   849,
     848,   249,   249,   249,   848,   249,   249,   249,   848,   249,
     848,   249,   249,   848,   249,   249,   503,   249,   249,   249,
     249,   848,   249,   848,   716,   323,   381,   747,    59,   745,
     727,   249,   249,   849,   697,   698,   848,   249,   249,   776,
     849,    59,   733,   249,   727,    59,   731,    59,   735,   378,
     737,    59,   743,    59,   739,    59,   741,   848,   249,   249,
     249,   332,   377,   848,   249,   848,   249,   681,   682,   848,
     769,   684,   685,   848,   849,   249,   378,   848,   330,   849,
     849,   849,   849,   849,   849,   849,   249,   849,   327,   849,
     849,   372,   372,   125,   125,   125,   125,   125,   249,   249,
     125,   249,   249,   848,   249,   848,   249,   249,   249,   249,
     849,   848,   848,   848,   848,   849,   849,   164,   232,   849,
     848,   848,   848,   848,   848,   620,   848,   620,   621,   848,
     621,   408,   849,   849,   249,   848,   849,   125,   848,   778,
     848,   412,   848,   125,   848,   125,   848,   125,   125,   125,
     125,    14,    20,    67,   218,   246,   249,   256,   272,   334,
     335,   365,   366,   455,   849,   848,   249,   848,   848,   648,
     849,   249,   849,   849,   649,   849,   372,   249,   848,   648,
     249,   249,   249,   372,   372,   848,   372,   849,   372,    55,
     249,   848,   849,   249,   848,   849,   249,   848,   125,   249,
     848,   848,   249,   848,   849,   778,   412,    72,   848,   412,
     849,   848,   849,   848,   849,   849,   848,   848,   848,   848,
     848,   125,   103,   333,   249,   848,   849,   849,   849,   249,
     848,   249,   778,   412,   249,    10,    15,    68,   114,   215,
     379,   409,   849,   848,   848,   848,   848,   848,   848,   249,
     440,   812,   848,   249,   249,   848,   849,   849,   849,   849,
     849,   249,   848,   249,   249,   848,   249,   249,   249,   249,
     849,   249,   249,   848,   249,   849,   249,   849,   848,   849,
     148,   849,   848,   407,   848,   249,   249,   849,   848,   125,
     249,   848,   249,   848,   848,   117,   848,   849,   412,   849,
     849,   849,   380,   125,   849,   849,   803,   848,   802,   848,
     249,   849,   249,   849,   805,   848,   249,   848,   372,    84,
     249,   449,   589,   848,   249,   848,   249,   848,   848,   249,
     849,   848,   249,   298,   371,   420,   848,    16,    29,    43,
      48,    57,    60,    61,    65,    66,   121,   142,   143,   144,
     145,   151,   155,   158,   180,   182,   196,   210,   211,   219,
     221,   224,   226,   228,   229,   233,   235,   236,   237,   238,
     253,   255,   261,   290,   295,   296,   297,   299,   301,   303,
     308,   309,   310,   328,   337,   349,   350,   351,   352,   353,
     354,   355,   357,   358,   359,   360,   367,   368,   370,   374,
     387,   389,   394,   399,   419,   423,   424,   425,   443,   446,
     515,   849,   249,   849,   848,   848,   249,   849,   849,     9,
     824,   848,   848,   249,   249,   849,   848,   848,   249,   249,
     848,   848,   125,   115,   116,   125,   341,   342,   453,   454,
     455,   456,   458,   460,   461,   462,   463,   464,   465,   466,
     467,   469,   470,   471,   472,   473,   474,   475,   476,   477,
     478,   479,   480,   481,   482,   484,   485,   486,   487,   488,
     489,   490,   491,   492,   493,   848,   848,   848,   848,   372,
     848,   848,   372,   513,   833,   833,   833,   833,   833,   849,
     848,   125,   126,   372,   848,   372,   833,   249,   848,   848,
     848,   249,   646,   125,   372,   849,   848,   372,   372,   372,
     372,   372,   372,   372,   372,   848,   848,   510,   849,   125,
     125,   848,   848,   849,   372,   372,   848,   372,   125,   372,
     372,   125,   125,   372,   372,   372,   372,   372,   848,   848,
     848,   249,   249,   249,   249,   249,   125,   849,   833,   249,
     125,   848,   249,   778,   848,   412,   249,   849,   849,   849,
     249,   249,   848,   707,   849,   708,   249,   849,   848,   249,
     249,   249,   849,   249,   849,   849,   249,   848,   780,   849,
     249,   249,   817,   818,   849,   249,   848,   223,   249,   849,
     849,   781,   848,   758,   692,   848,   249,   249,   755,   756,
     848,   756,   848,    62,   249,   849,   814,   815,   849,   849,
     849,   849,   249,   849,   848,   751,   752,   848,   752,   848,
     249,   849,   849,   372,   372,   849,   848,   249,   249,   848,
     848,   249,   814,   814,   249,   249,   249,   848,   381,   849,
     848,    59,   849,   808,   849,   698,   848,   249,   848,    59,
     849,   848,    59,   849,   848,    59,   849,   848,   378,   849,
     848,    59,   849,   848,    59,   849,   848,    59,   849,   249,
     249,   849,   249,   498,   499,   500,   249,   682,   249,   849,
     685,   848,   849,   848,   848,   849,   849,   249,   849,   849,
     849,   849,   849,   249,   849,   848,   249,   848,   125,   249,
     125,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   848,   848,   249,
     249,   372,   849,   249,   249,   249,   249,   849,   848,   848,
     849,   849,   125,   848,   849,   849,   125,   249,   249,   848,
     242,   249,   849,   848,   848,   848,   848,   848,   848,   848,
     372,   849,   249,   249,   849,   849,   849,   249,   249,   249,
     728,   779,   848,   249,   249,   249,   849,   249,   848,   849,
     249,   848,   849,   249,   849,   249,   249,   849,   848,   848,
     849,   848,   249,   849,   249,   849,   849,   849,   695,   848,
     695,   695,   695,   695,   848,   125,   125,   728,   249,   849,
     849,   849,   728,   249,   728,   848,   849,   849,   113,   215,
     249,   379,   440,   812,   848,   848,   849,   849,   849,   849,
     777,   848,   777,   848,   848,   695,   848,   848,   249,   440,
     848,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   849,   849,   849,   849,   849,
     177,   249,   407,   848,   249,   848,   249,   249,   249,   848,
     412,   849,   848,   117,   849,   249,   372,   848,   249,   849,
     249,   249,   249,   249,   249,   249,   249,   848,   249,   848,
     801,   848,   249,   249,   249,   848,   249,   249,   249,   249,
     249,   120,   249,   849,   849,   249,   849,   249,   849,   249,
     249,   249,   249,   249,   849,    17,   227,   444,   806,   372,
     806,    17,    58,   849,   806,   372,   372,   249,   338,   848,
     249,   249,   848,   849,   849,   849,   849,   249,   191,   848,
     806,   848,   848,   848,   225,   848,   227,   848,   338,   848,
     249,   848,   848,   848,   848,   848,   249,   848,   249,   848,
     121,   212,   848,   125,   249,   848,   848,   300,   848,   848,
     849,   849,   849,   848,   372,   338,   848,   849,   848,   848,
     848,   849,   848,   848,   848,   848,   848,   848,   249,   848,
     849,   848,   849,   849,   849,   849,   848,   249,   848,   848,
     249,   184,   849,   249,   372,   849,   849,   849,   249,   848,
     249,   249,   249,   848,   249,   249,   849,   249,   848,   249,
     249,   849,   125,   849,   849,   832,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   249,   848,    94,    94,   848,
     849,   125,   372,   848,   849,   849,   849,   249,   848,   848,
     125,   848,   848,   849,   125,   848,   848,   848,   848,   249,
     849,   249,   848,   412,   848,   249,   849,   849,   249,   249,
     249,   249,   249,   848,   849,   849,   849,   249,   818,   849,
     728,   249,   249,   848,   249,   848,   249,   849,   848,   849,
     756,   848,   849,   249,   815,   849,   249,   849,   849,   849,
     849,   249,   752,   848,   249,   440,   812,   249,   249,   384,
     249,   249,   249,   249,   670,   848,   848,   323,   249,   848,
     849,   249,   848,   849,   249,   848,   849,   249,   848,   849,
     249,   848,   849,   249,   848,   849,   108,   249,   848,   832,
     249,   848,   849,   249,   848,   849,   249,   848,   849,   249,
     332,   849,   249,   498,   249,   498,   249,   695,   249,   848,
     249,   848,   849,   249,   848,   848,   848,   848,   848,   249,
     249,   848,   249,   249,   249,   249,   249,   615,   848,   249,
     849,   849,   849,   249,   849,   849,   249,   849,   249,   249,
     249,   375,   848,   249,   242,   249,   157,   848,   125,   157,
     848,   125,   849,   249,   242,   249,   779,   849,   849,   849,
     849,   849,   848,   849,   375,   848,   848,   249,   375,   849,
     249,   849,   849,   849,   249,   848,   249,   249,   249,   249,
     362,   848,   848,   728,   849,   849,   849,   728,   375,   848,
     849,   249,   849,   849,   849,   849,   848,   113,   215,   249,
     379,   440,   849,   849,   849,   134,   249,   849,   849,   249,
     848,   249,   249,   812,   848,   249,   849,   249,   249,   849,
     249,   848,   848,   249,   849,   849,   849,   849,   111,   249,
     249,   321,   849,   249,   402,   848,   125,   249,   848,   249,
     372,   849,   848,   249,   372,   249,   849,   849,   848,   848,
     249,   249,   848,   249,   848,   249,   848,   249,   249,   849,
     249,   315,   249,   445,   848,   249,   249,   249,   848,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   849,   249,   249,   230,   231,   249,   231,   249,   848,
     249,   249,   249,   249,   249,   848,   849,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     249,   849,   849,   848,   249,   249,   249,   249,   249,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   249,
     849,   249,   249,   849,   249,   249,   249,   249,   249,   249,
     174,   249,   249,   849,   249,   249,   249,   848,   849,   249,
     249,   849,   249,   849,   249,   848,   849,   249,   848,   848,
     249,   849,   848,   848,   849,   849,   249,   457,   849,   249,
     457,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   249,   457,   849,   849,   849,   849,   849,
     849,   849,   849,   848,   848,   848,   848,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   249,   510,   848,
     849,   848,   125,   848,   848,   125,   848,   848,   849,   848,
     375,   848,   249,   849,   249,   249,   849,   849,   849,   728,
     728,   249,   249,   849,   249,   849,   848,   849,   249,   849,
     849,   849,   249,   249,   249,   249,   849,   848,   249,   372,
     775,   280,   282,   283,   284,   285,   286,   287,   288,   289,
     671,   670,   848,   249,   849,   249,   849,   849,   249,   249,
     848,   849,   249,   849,   249,   249,   849,   249,   249,   849,
     249,   249,   372,   832,   108,   249,   249,   849,   249,   249,
     849,   249,   249,   849,   249,   249,   249,   249,   849,   849,
     849,   249,   249,   249,   249,   848,   249,   249,   327,   327,
     848,   848,   848,   615,   849,   849,   849,   849,   849,   249,
     848,   249,   249,   848,   848,   848,   848,   849,   249,   249,
     848,   249,   849,   249,   848,   249,   849,   249,   848,   849,
     249,   848,   249,   849,   249,   249,   249,   249,   849,   249,
     849,   249,   849,   849,   849,   848,   849,   249,   849,   849,
     849,   249,   849,   134,   249,   849,   249,   849,   849,   849,
     848,   849,   849,   849,   849,   849,   849,   249,   249,   812,
     849,   849,   849,   249,   849,   249,   848,   849,   849,   849,
     249,   849,   111,   249,   849,   249,   402,   407,   848,   249,
     849,   249,   249,   372,   849,   249,   249,   849,   249,   848,
     249,   848,   249,   249,   849,   249,   848,   249,   249,   848,
     249,   450,   249,   249,   231,   249,   848,   848,   231,   249,
     249,   249,   249,   249,   249,   249,   249,   249,   249,   848,
     249,   848,   249,   849,   849,   249,   849,   848,   848,   848,
     848,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   848,   848,   848,
     848,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   848,   849,   849,   848,   249,   375,   848,   848,   849,
     249,   125,   849,   249,   728,   249,   849,   249,   849,   249,
     249,   849,   249,   249,   249,   249,   812,   249,   848,   849,
     849,   670,   670,   323,   849,   849,   849,   849,   849,   849,
     249,   849,   849,   849,   849,   849,   849,   249,   249,   372,
     832,   849,   849,   849,   849,   849,   849,   249,   249,   849,
     249,   849,   249,   249,   849,   849,   249,   848,   849,   849,
     849,   249,   849,   849,   848,   157,   848,   157,   848,   848,
     249,   849,   249,   249,   848,   249,   849,   249,   249,   249,
     849,   849,   849,   848,   242,   249,   249,   849,   849,   849,
     249,   761,   849,   134,   249,   849,   849,   249,   849,   134,
     249,   849,   249,   849,   849,   849,   249,   849,   849,   249,
     249,   812,   849,   249,   249,   849,   249,   848,   849,   849,
     249,   249,   249,   249,   407,   849,   849,   249,   249,   249,
     249,   372,   849,   249,   848,   249,   848,   249,   848,   249,
     848,   249,   450,   249,   848,   249,   249,   848,   249,   848,
     249,   848,   249,   249,   440,   442,   812,   848,   849,   848,
     849,   849,   849,   849,   249,   457,   849,   849,   249,   457,
     849,   849,   249,   849,   249,   849,   849,   249,   457,   849,
     849,   849,   849,   849,   849,   249,   457,   849,   849,   249,
     849,   849,   849,   849,   849,   848,   848,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   848,   849,   848,   249,   249,   249,   848,   249,
     249,   249,   249,   849,   849,   849,   249,   849,   249,   249,
     849,   249,   849,   249,   849,   832,   849,   249,   249,   249,
     849,   249,   849,   249,   849,   249,   249,   848,   848,   849,
     249,   849,   849,   849,   249,   849,   249,   848,   848,   249,
     249,   849,   249,   849,   849,   849,   849,   249,   849,   849,
     849,   249,   849,   761,   134,   249,   849,   249,   761,   134,
     249,   849,   849,   849,   849,   849,   849,   249,   249,   812,
     848,   848,   249,   849,   849,   849,   249,   111,   249,   249,
     249,   249,   848,   249,   848,   249,   249,   249,   249,   249,
     249,   249,   848,   849,   249,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   249,   849,   249,
     849,   849,   849,   849,   849,   849,   849,   249,   457,   849,
     849,   849,   849,   249,   849,   849,   849,   249,   849,   849,
     849,   849,   849,   849,   849,   249,   457,   849,   249,   457,
     849,   249,   457,   849,   849,   849,   849,   849,   849,   849,
     848,   249,   849,   125,   249,   849,   849,   249,   249,   249,
     249,   249,   249,   832,   832,   849,   249,   249,   249,   848,
     848,   249,   849,   849,   849,   849,   249,   249,   249,   249,
     242,   249,   849,   849,   249,   379,   849,   849,   761,   761,
     849,   249,   761,   134,   249,   849,   849,   849,   849,   849,
     249,   849,   849,   849,   849,   111,   249,   249,   249,   848,
     249,   848,   249,   249,   440,   249,   849,   849,   249,   849,
     849,   249,   849,   249,   849,   249,   849,   249,   457,   249,
     849,   249,   457,   849,   249,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   249,   849,   249,   457,   849,   249,
     849,   249,   849,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   249,   457,   849,   849,   849,
     849,   849,   849,   849,   249,   249,   849,   849,   849,   249,
     249,   832,   848,   849,   249,   372,   849,   849,   249,   849,
     249,   849,   849,   849,   249,   379,   849,   849,   249,   761,
     249,   379,   849,   761,   761,   849,   849,   849,   849,   249,
     849,   249,   812,   849,   849,   849,   249,   249,   848,   249,
     848,   848,   849,   849,   849,   849,   249,   849,   249,   849,
     849,   849,   849,   849,   249,   249,   849,   849,   849,   849,
     849,   249,   849,   249,   849,   849,   849,   849,   849,   249,
     849,   249,   249,   849,   249,   849,   849,   849,   849,   849,
     849,   849,   249,   457,   849,   249,   457,   849,   849,   249,
     849,   249,   849,   849,   849,   849,   849,   249,   125,   249,
     849,   832,   849,   249,   848,   849,   615,   249,   849,   849,
     849,   849,   849,   249,   849,   249,   379,   849,   761,   849,
     249,   379,   849,   249,   761,   848,   849,   849,   849,   849,
     249,   849,   849,   849,   249,   848,   249,   249,   812,   249,
     249,   849,   249,   249,   849,   849,   849,   249,   849,   249,
     849,   249,   849,   849,   849,   249,   849,   249,   849,   249,
     849,   849,   849,   249,   849,   249,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   849,   849,   849,   849,   849,
     249,   849,   249,   849,   849,   849,   849,   849,   249,   249,
     249,   832,   249,   249,   615,   849,   249,   249,   849,   849,
     849,   249,   849,   249,   249,   849,   249,   379,   761,   249,
     849,   849,   849,   849,   849,   849,   249,   849,   249,   249,
     849,   849,   249,   849,   249,   849,   849,   849,   249,   848,
     849,   849,   849,   249,   848,   249,   848,   249,   849,   849,
     249,   849,   249,   849,   849,   849,   849,   849,   849,   849,
     249,   849,   249,   849,   249,   849,   849,   849,   849,   849,
     849,   849,   249,   249,   849,   249,   249,   849,   249,   249,
     249,   849,   249,   849,   849,   249,   849,   849,   249,   849,
     849,   249,   249,   249,   849,   249,   849,   249,   249,   249,
     848,   249,   457,   849,   249,   457,   849,   249,   457,   249,
     249,   849,   249,   249,   849,   249,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   249,   249,   849,   249,   849,
     849,   849,   849,   849,   849,   249,   849,   249,   848,   249,
     849,   849,   849,   249,   848,   249,   848,   848,   849,   849,
     849,   849,   849,   249,   457,   249,   849,   249,   848,   249,
     849,   249,   849,   849,   849,   849,   849,   249,   249,   249,
     849,   249,   849,   249,   457,   849,   249,   457,   849,   249,
     457,   849,   249,   849,   249,   848,    75,   249,   849,   849,
     249,   249,   249,   848,   849,   249,   457,   849,   249,   457,
     849,   849,   249,   849,   249,   849,   849,   249,   849,   249,
     849,   849,   849,   249,   848,   249,   848,   849,   849,   849,
     849,   849,   249,   457,   249,   848,   849,   849,   249,   848,
     849,   849,   849,   849,   849,   849,   249,   849,   249,   249,
     849,   849,   849,   849,   249,   249,   849,   249,   457,   849,
     249,   457,   849,   849,   249,   249,    75,   249,   849,   849,
     249,   849,   249,   849,   249,   849,   249,   849,   249,   249,
     249,   849,   249,   849,   849,   849,   849,   849,   849,   849,
     849,   849,   849,   849,   849,   249,   249,   849,   849,   849,
     249,   849,   249,   849,   249,   849,   849,   849,   849,   249,
     249,   249,   849,   249,   249,   849,   849,   249,   249,   849,
     848,   249,   849,   249,   249,    75,   249,   849,   849,   249,
     849,   849,   849,   249,   249,   849,   849,   848,   249,   849,
     249,   849,   249,   849,   849,   249,   848,   848,   849,   849,
     849,   249
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   590,   591,   592,   592,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   593,   594,   595,   596,   596,   596,   596,   597,
     597,   598,   598,   598,   598,   598,   598,   598,   598,   598,
     598,   598,   598,   598,   598,   599,   600,   601,   601,   602,
     602,   603,   603,   604,   604,   604,   604,   604,   604,   604,
     604,   604,   604,   605,   606,   606,   606,   606,   606,   606,
     606,   606,   606,   606,   606,   606,   606,   607,   607,   608,
     608,   608,   608,   608,   609,   609,   609,   610,   610,   610,
     611,   611,   611,   611,   611,   611,   611,   611,   612,   612,
     612,   612,   612,   612,   613,   613,   614,   614,   614,   614,
     615,   615,   616,   616,   617,   617,   618,   618,   619,   619,
     620,   620,   620,   621,   622,   622,   623,   623,   623,   623,
     624,   624,   624,   624,   625,   625,   626,   626,   626,   626,
     627,   628,   629,   630,   631,   632,   632,   632,   632,   632,
     633,   633,   633,   633,   633,   633,   633,   633,   633,   633,
     633,   634,   634,   634,   634,   634,   634,   634,   634,   634,
     634,   634,   634,   634,   634,   634,   634,   634,   634,   634,
     634,   634,   634,   634,   634,   634,   634,   634,   634,   634,
     634,   634,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   636,   636,   637,   637,   638,
     638,   639,   639,   640,   640,   641,   641,   642,   643,   643,
     643,   644,   644,   644,   644,   644,   644,   644,   644,   644,
     644,   644,   644,   644,   644,   644,   644,   644,   644,   644,
     645,   646,   647,   647,   647,   647,   648,   648,   648,   648,
     649,   650,   650,   650,   650,   650,   651,   651,   652,   653,
     653,   653,   653,   653,   653,   653,   653,   653,   654,   654,
     655,   656,   657,   658,   658,   659,   660,   661,   662,   662,
     662,   662,   662,   662,   662,   662,   663,   663,   664,   664,
     664,   664,   665,   665,   666,   667,   668,   669,   669,   669,
     670,   670,   671,   671,   671,   671,   671,   671,   671,   671,
     671,   671,   672,   672,   672,   673,   674,   674,   675,   675,
     675,   676,   677,   678,   678,   678,   678,   678,   679,   679,
     680,   681,   681,   682,   682,   683,   684,   684,   685,   686,
     686,   686,   686,   686,   687,   687,   687,   687,   688,   688,
     688,   688,   689,   689,   689,   690,   691,   691,   692,   692,
     692,   693,   693,   694,   695,   695,   696,   697,   697,   698,
     698,   698,   699,   699,   700,   701,   701,   702,   703,   703,
     704,   705,   705,   706,   707,   708,   708,   709,   710,   710,
     711,   711,   712,   712,   713,   713,   714,   714,   715,   716,
     716,   716,   716,   717,   717,   717,   717,   717,   718,   718,
     718,   718,   718,   719,   719,   719,   719,   719,   719,   720,
     720,   721,   721,   722,   722,   722,   722,   723,   723,   724,
     725,   725,   725,   725,   725,   725,   725,   725,   726,   726,
     726,   726,   727,   727,   727,   727,   727,   727,   728,   728,
     729,   729,   730,   730,   731,   731,   731,   732,   732,   733,
     733,   733,   734,   734,   735,   735,   735,   736,   736,   737,
     737,   737,   737,   737,   738,   738,   739,   739,   739,   740,
     740,   741,   741,   741,   742,   742,   743,   743,   743,   744,
     744,   745,   745,   745,   746,   746,   747,   747,   747,   748,
     748,   748,   749,   749,   750,   750,   750,   751,   751,   751,
     751,   751,   752,   753,   753,   754,   754,   754,   755,   755,
     755,   756,   756,   757,   757,   758,   758,   759,   759,   760,
     760,   760,   760,   760,   760,   760,   760,   760,   760,   760,
     760,   760,   760,   760,   760,   760,   760,   760,   760,   760,
     760,   760,   760,   760,   760,   760,   760,   760,   760,   760,
     760,   760,   760,   760,   760,   760,   760,   760,   760,   760,
     760,   760,   760,   760,   760,   760,   760,   760,   760,   760,
     760,   760,   761,   762,   762,   763,   763,   763,   763,   763,
     764,   764,   764,   764,   764,   764,   764,   764,   764,   765,
     765,   766,   766,   766,   766,   766,   766,   766,   766,   766,
     766,   766,   767,   767,   767,   768,   768,   769,   770,   770,
     770,   771,   771,   771,   771,   771,   771,   772,   772,   772,
     772,   772,   772,   772,   772,   772,   773,   773,   773,   773,
     773,   773,   773,   773,   773,   773,   773,   773,   773,   773,
     773,   773,   773,   773,   773,   773,   773,   773,   773,   773,
     773,   773,   773,   774,   774,   774,   775,   775,   775,   775,
     775,   776,   777,   777,   778,   778,   778,   778,   779,   780,
     781,   781,   782,   782,   783,   783,   784,   784,   785,   785,
     786,   786,   786,   786,   787,   787,   788,   788,   788,   788,
     788,   788,   788,   788,   788,   788,   788,   788,   788,   788,
     788,   788,   788,   788,   789,   789,   789,   789,   789,   790,
     790,   790,   791,   792,   793,   793,   794,   794,   794,   795,
     795,   795,   796,   796,   796,   796,   796,   796,   796,   796,
     796,   796,   796,   796,   796,   797,   797,   798,   799,   799,
     799,   799,   799,   800,   800,   800,   800,   800,   800,   800,
     800,   800,   800,   800,   801,   801,   801,   801,   801,   801,
     801,   801,   801,   801,   801,   801,   802,   802,   803,   803,
     804,   804,   804,   804,   804,   804,   804,   804,   804,   804,
     804,   805,   805,   806,   806,   806,   806,   806,   806,   806,
     806,   806,   806,   806,   806,   806,   806,   806,   806,   806,
     806,   806,   806,   806,   806,   806,   806,   806,   806,   806,
     807,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   808,   808,   808,
     808,   808,   808,   808,   808,   808,   808,   809,   810,   811,
     811,   812,   812,   812,   812,   812,   812,   812,   812,   812,
     812,   812,   813,   813,   813,   813,   813,   813,   813,   813,
     813,   813,   814,   814,   815,   816,   816,   817,   817,   818,
     819,   820,   821,   822,   823,   824,   824,   825,   825,   825,
     825,   825,   825,   825,   825,   825,   825,   825,   825,   825,
     825,   825,   825,   825,   825,   825,   825,   825,   825,   825,
     825,   825,   825,   825,   826,   826,   826,   827,   828,   828,
     828,   829,   829,   829,   829,   829,   829,   829,   829,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   830,   830,   830,   830,   830,   830,   830,   830,   830,
     830,   831,   831,   831,   832,   832,   833,   833,   834,   834,
     834,   835,   835,   836,   836,   836,   836,   836,   836,   836,
     836,   836,   836,   836,   836,   836,   836,   836,   836,   836,
     836,   836,   836,   836,   836,   836,   836,   836,   836,   836,
     836,   836,   836,   836,   836,   836,   836,   837,   837,   838,
     838,   839,   839,   839,   839,   839,   840,   840,   841,   841,
     841,   841,   841,   841,   841,   841,   841,   841,   841,   841,
     841,   841,   841,   841,   841,   841,   841,   841,   841,   841,
     841,   841,   841,   841,   841,   841,   841,   841,   841,   841,
     841,   841,   841,   841,   841,   841,   841,   841,   841,   841,
     841,   841,   841,   841,   841,   841,   841,   841,   842,   842,
     843,   844,   844,   844,   845,   845,   846,   846,   847,   847,
     847,   848,   848,   849,   849,   849,   849
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
       4,     5,     2,     3,     3,     3,     2,     4,     1,     3,
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
      14,    15,     9,    10,    12,    13,    14,    15,     9,     7,
       5,    13,    11,     9,     9,     7,     5,    13,    11,     9,
       9,     7,     5,    13,    11,     9,     8,    10,     8,    10,
       9,    11,     9,    11,    11,    13,    11,    13,     8,    10,
      12,    14,     8,    10,    12,    14,     8,    12,     9,    13,
      10,    11,    13,    14,    15,    16,    10,    11,    13,    14,
      15,    16,     7,    13,    15,    17,    19,    14,    16,    14,
      16,    15,    17,    15,    17,    17,    19,    17,    19,    14,
      16,    18,    20,    13,    15,    17,    19,    14,    16,    18,
      20,     7,    11,    13,    17,    14,    18,     8,    12,    14,
      18,    15,    19,     8,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
       7,     9,    10,    10,    11,    12,    13,    10,    11,    12,
      13,     7,     8,     9,    10,    11,    12,    13,    22,     5,
       5,     2,     4,     5,     0,     2,     0,     2,     4,     6,
       8,     2,     3,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     3,     2,     3,     2,     2,     2,
       2,     5,     6,     3,     6,     3,     2,     1,     2,     3,
       2,     3,     4,     2,     2,     3,     1,     2,     3,     2,
       3,     2,     3,     2,     2,     2,     2,     3,     2,     3,
       3,     2,     3,     4,     2,     2,     2,     2,     2,     3,
       4,     2,     2,     2,     2,     2,     2,     2,     3,     4,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       3,     4,     2,     3,     4,     2,     2,     2,     2,     2,
       2,     2,     2,     3,     3,     2,     6,     3,     2,     3,
       5,     2,     5,     3,     2,     3,     2,     3,     2,     3,
       2,     1,     1,     1,     1,     1,     1
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
#line 196 "p.y" /* yacc.c:1646  */
    { 
          if(domain->solInfo().piecewise || domain->solInfo().freeplay) domain->solInfo().activatePiecewise();
          return 0;
         }
#line 5465 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 6:
#line 208 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; delete (yyvsp[0].bclist); }
#line 5471 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 7:
#line 210 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5477 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 16:
#line 220 "p.y" /* yacc.c:1646  */
    { int j = geoSource->getLocalIndex();
          if(geoSource->elementLumpingWeightLocalSize(j)>0) geoSource->setLocalIndex(j+1); }
#line 5484 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 26:
#line 232 "p.y" /* yacc.c:1646  */
    {}
#line 5490 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 28:
#line 235 "p.y" /* yacc.c:1646  */
    {}
#line 5496 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 34:
#line 242 "p.y" /* yacc.c:1646  */
    {}
#line 5502 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 40:
#line 249 "p.y" /* yacc.c:1646  */
    { if(geoSource->setUsddLocation((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)    return -1; delete (yyvsp[0].bclist); }
#line 5509 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 55:
#line 266 "p.y" /* yacc.c:1646  */
    { domain->setMFTT((yyvsp[0].mftval).first, (yyvsp[0].mftval).second); }
#line 5515 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 56:
#line 268 "p.y" /* yacc.c:1646  */
    { domain->setHFTT((yyvsp[0].hftval).first, (yyvsp[0].hftval).second); }
#line 5521 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 89:
#line 302 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5527 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 90:
#line 304 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5533 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 91:
#line 306 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5539 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 92:
#line 308 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichletFluid((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5545 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 93:
#line 310 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichletFluid((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5551 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 94:
#line 312 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5557 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 111:
#line 330 "p.y" /* yacc.c:1646  */
    { if(domain->setComplexNeuman((yyvsp[0].cxbclist)->n,(yyvsp[0].cxbclist)->d) < 0) return -1; }
#line 5563 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 113:
#line 333 "p.y" /* yacc.c:1646  */
    { if(domain->setComplexDirichlet((yyvsp[0].cxbclist)->n,(yyvsp[0].cxbclist)->d) < 0) return -1; }
#line 5569 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 117:
#line 338 "p.y" /* yacc.c:1646  */
    { if(geoSource->setDirichlet((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5575 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 118:
#line 340 "p.y" /* yacc.c:1646  */
    { if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 5581 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 125:
#line 348 "p.y" /* yacc.c:1646  */
    {}
#line 5587 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 134:
#line 359 "p.y" /* yacc.c:1646  */
    {}
#line 5593 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 135:
#line 361 "p.y" /* yacc.c:1646  */
    {}
#line 5599 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 136:
#line 363 "p.y" /* yacc.c:1646  */
    {}
#line 5605 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 137:
#line 365 "p.y" /* yacc.c:1646  */
    {}
#line 5611 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 138:
#line 367 "p.y" /* yacc.c:1646  */
    {}
#line 5617 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 139:
#line 369 "p.y" /* yacc.c:1646  */
    {}
#line 5623 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 140:
#line 371 "p.y" /* yacc.c:1646  */
    {}
#line 5629 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 141:
#line 373 "p.y" /* yacc.c:1646  */
    {}
#line 5635 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 153:
#line 388 "p.y" /* yacc.c:1646  */
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[-2].ival)); 
            domain->solInfo().nsample = (yyvsp[-1].ival);
          }
#line 5644 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 154:
#line 395 "p.y" /* yacc.c:1646  */
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[-1].ival));
          }
#line 5652 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 156:
#line 402 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-3].ival) == OutputInfo::Attribute)  geoSource->setAttributeGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1);
          else if ((yyvsp[-3].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[-3].ival));  exit(-1); }
        }
#line 5661 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 157:
#line 407 "p.y" /* yacc.c:1646  */
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
#line 5677 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 158:
#line 419 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-4].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[-2].ival)-1, (yyvsp[-1].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[-4].ival));  exit(-1); }
        }
#line 5685 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 160:
#line 426 "p.y" /* yacc.c:1646  */
    { geoSource->setGroupRandomProperty((yyvsp[-4].ival)-1,(yyvsp[-3].rprop),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 5691 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 161:
#line 430 "p.y" /* yacc.c:1646  */
    { domain->solInfo().curSweepParam = 0; }
#line 5697 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 162:
#line 432 "p.y" /* yacc.c:1646  */
    { domain->solInfo().curSweepParam = (yyvsp[-1].ival); }
#line 5703 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 163:
#line 434 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-1].fval)); }
#line 5709 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 164:
#line 436 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-3].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies1(2.0*PI*(yyvsp[-3].fval), 2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5717 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 165:
#line 440 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-3].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies2(2.0*PI*(yyvsp[-3].fval), 2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5725 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 166:
#line 444 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-4].fval));
            domain->setFrequencySet(domain->solInfo().curSweepParam);
            domain->addFrequencies(2.0*PI*(yyvsp[-4].fval), 2.0*PI*(yyvsp[-3].fval), (yyvsp[-2].ival), (yyvsp[-1].ival)); }
#line 5733 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 167:
#line 448 "p.y" /* yacc.c:1646  */
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
#line 5759 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 168:
#line 470 "p.y" /* yacc.c:1646  */
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
#line 5787 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 169:
#line 494 "p.y" /* yacc.c:1646  */
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
#line 5824 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 175:
#line 534 "p.y" /* yacc.c:1646  */
    {
          if((yyvsp[-3].ival) == 1) {
            domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
            domain->solInfo().getSweepParams()->alphaD = (yyvsp[-1].fval);
            domain->solInfo().getSweepParams()->betaD = (yyvsp[-2].fval);
            domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
          }
          else return -1; // only RAYDAMP is allowed here
        }
#line 5838 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 176:
#line 546 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_pivot = true; domain->solInfo().getSweepParams()->pade_tol = (yyvsp[-1].fval); }
#line 5844 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 177:
#line 550 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_poles = true; }
#line 5850 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 178:
#line 552 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getSweepParams()->pade_poles = true; 
          domain->solInfo().getSweepParams()->pade_poles_sigmaL = (yyvsp[-2].fval); domain->solInfo().getSweepParams()->pade_poles_sigmaU = (yyvsp[-1].fval); }
#line 5857 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 179:
#line 557 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().curSweepParam == 0) geoSource->setImpe((yyvsp[-1].fval));
          domain->setFrequencySet(domain->solInfo().curSweepParam);
          domain->addCoarseFrequency(2.0*PI*(yyvsp[-1].fval)); }
#line 5865 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 180:
#line 562 "p.y" /* yacc.c:1646  */
    { domain->addFrequencies(2.0*PI*(yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 5871 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 181:
#line 566 "p.y" /* yacc.c:1646  */
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
#line 5917 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 182:
#line 608 "p.y" /* yacc.c:1646  */
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
#line 5971 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 184:
#line 661 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInput = bool((yyvsp[-1].ival)); }
#line 5977 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 185:
#line 663 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInput = bool((yyvsp[-2].ival));
            std::string prefix = (yyvsp[-1].strval);
            clusterData_ = prefix + ".msh";
            decomposition_ = prefix + ".dec";
            connectivity_ = prefix + ".con";
            subdomains_ = prefix + ".sub";
          }
#line 5989 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 186:
#line 671 "p.y" /* yacc.c:1646  */
    { geoSource->binaryOutput = bool((yyvsp[-1].ival)); }
#line 5995 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 187:
#line 673 "p.y" /* yacc.c:1646  */
    { geoSource->binaryOutput = bool((yyvsp[-2].ival));
            int len = strlen((yyvsp[-1].strval));
            char *file = new char[len+5];
            strcpy(file, (yyvsp[-1].strval));
            strcat(file,".con");
            geoSource->setGlob(file);
          }
#line 6007 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 188:
#line 681 "p.y" /* yacc.c:1646  */
    { geoSource->setGeo((yyvsp[-1].strval)); }
#line 6013 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 189:
#line 683 "p.y" /* yacc.c:1646  */
    { geoSource->setDecomp((yyvsp[-1].strval)); }
#line 6019 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 190:
#line 685 "p.y" /* yacc.c:1646  */
    { geoSource->setGlob((yyvsp[-1].strval)); }
#line 6025 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 191:
#line 687 "p.y" /* yacc.c:1646  */
    { geoSource->setMatch((yyvsp[-1].strval)); }
#line 6031 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 192:
#line 689 "p.y" /* yacc.c:1646  */
    { geoSource->setCpuMap((yyvsp[-1].strval)); }
#line 6037 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 193:
#line 693 "p.y" /* yacc.c:1646  */
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[-1].ival)); 
#endif
	}
#line 6047 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 194:
#line 701 "p.y" /* yacc.c:1646  */
    {if(decInit==0) decInit = new DecInit(); }
#line 6053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 195:
#line 703 "p.y" /* yacc.c:1646  */
    {decInit->file = strdup((yyvsp[-1].strval));}
#line 6059 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 196:
#line 705 "p.y" /* yacc.c:1646  */
    {decInit->nsubs = (yyvsp[-1].ival); }
#line 6065 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 197:
#line 707 "p.y" /* yacc.c:1646  */
    {decInit->weight = true; }
#line 6071 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 198:
#line 709 "p.y" /* yacc.c:1646  */
    {decInit->memory = true; }
#line 6077 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 199:
#line 711 "p.y" /* yacc.c:1646  */
    {decInit->exitAfterDec = true;}
#line 6083 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 200:
#line 713 "p.y" /* yacc.c:1646  */
    {decInit->skip = true;}
#line 6089 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 201:
#line 715 "p.y" /* yacc.c:1646  */
    {decInit->nosa = true; }
#line 6095 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 202:
#line 717 "p.y" /* yacc.c:1646  */
    {decInit->trivial = true; }
#line 6101 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 203:
#line 719 "p.y" /* yacc.c:1646  */
    {decInit->trivial = true; randomShuffle = bool((yyvsp[-1].ival)); }
#line 6107 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 204:
#line 721 "p.y" /* yacc.c:1646  */
    {decInit->fsgl = true; }
#line 6113 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 205:
#line 723 "p.y" /* yacc.c:1646  */
    {allowMechanisms = true; }
#line 6119 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 206:
#line 725 "p.y" /* yacc.c:1646  */
    {useScotch = true; }
#line 6125 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 207:
#line 729 "p.y" /* yacc.c:1646  */
    {}
#line 6131 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 208:
#line 731 "p.y" /* yacc.c:1646  */
    { weightList[(yyvsp[-2].ival)] = (yyvsp[-1].fval); }
#line 6137 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 209:
#line 735 "p.y" /* yacc.c:1646  */
    {}
#line 6143 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 210:
#line 737 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Acoustic] = (yyvsp[-1].ival); }
#line 6149 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 211:
#line 739 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Structural] = (yyvsp[-1].ival); }
#line 6155 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 212:
#line 741 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Thermal] = (yyvsp[-1].ival); }
#line 6161 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 213:
#line 743 "p.y" /* yacc.c:1646  */
    { fieldWeightList[(int)Element::Fluid] = (yyvsp[-1].ival); }
#line 6167 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 214:
#line 746 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = 0; }
#line 6173 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 215:
#line 748 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = (yyvsp[-1].ival); }
#line 6179 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 216:
#line 750 "p.y" /* yacc.c:1646  */
    { (yyval.mftval).first->add((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6185 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 217:
#line 754 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = 0; }
#line 6191 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 218:
#line 756 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = (yyvsp[-1].ival); }
#line 6197 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 219:
#line 758 "p.y" /* yacc.c:1646  */
    { (yyval.hftval).first->add((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6203 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 220:
#line 762 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 6209 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 221:
#line 764 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 6215 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 222:
#line 766 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival); domain->setLoadFactorGrav((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6221 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 223:
#line 768 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-3].ival); domain->setLoadFactorTemp((yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6227 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 224:
#line 770 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-5].ival); domain->setLoadFactorGrav((yyvsp[-5].ival), (yyvsp[-3].ival)); domain->setLoadFactorTemp((yyvsp[-5].ival), (yyvsp[-1].ival)); }
#line 6233 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 225:
#line 772 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactor((yyval.ival), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 6239 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 226:
#line 774 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactorMFTT((yyval.ival), (yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6245 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 227:
#line 776 "p.y" /* yacc.c:1646  */
    { domain->setLoadFactorHFTT((yyval.ival), (yyvsp[-3].ival), (yyvsp[-1].ival)); }
#line 6251 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 235:
#line 789 "p.y" /* yacc.c:1646  */
    { geoSource->addCFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d); }
#line 6257 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 236:
#line 793 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).coefFlag = false; geoSource->addCoefInfo((yyvsp[-2].ival)-1,(yyvsp[0].coefdata)); }
#line 6263 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 237:
#line 795 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).c[6][0] = (yyvsp[-7].fval);
          (yyvsp[0].coefdata).c[6][1] = (yyvsp[-6].fval);
          (yyvsp[0].coefdata).c[6][2] = (yyvsp[-5].fval);
          (yyvsp[0].coefdata).c[6][3] = (yyvsp[-4].fval);
          (yyvsp[0].coefdata).c[6][4] = (yyvsp[-3].fval);
          (yyvsp[0].coefdata).c[6][5] = (yyvsp[-2].fval);
          (yyvsp[0].coefdata).coefFlag = false;
          geoSource->addCoefInfo((yyvsp[-8].ival)-1,(yyvsp[0].coefdata)); }
#line 6276 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 238:
#line 804 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).coefFlag = (yyvsp[-2].ival); geoSource->addCoefInfo((yyvsp[-3].ival)-1,(yyvsp[0].coefdata)); }
#line 6282 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 239:
#line 806 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].coefdata).c[6][0] = (yyvsp[-8].fval);
          (yyvsp[0].coefdata).c[6][1] = (yyvsp[-7].fval);
          (yyvsp[0].coefdata).c[6][2] = (yyvsp[-6].fval);
          (yyvsp[0].coefdata).c[6][3] = (yyvsp[-5].fval);
          (yyvsp[0].coefdata).c[6][4] = (yyvsp[-4].fval);
          (yyvsp[0].coefdata).c[6][5] = (yyvsp[-3].fval);
          (yyvsp[0].coefdata).coefFlag = (yyvsp[-2].ival);
          geoSource->addCoefInfo((yyvsp[-9].ival)-1,(yyvsp[0].coefdata)); }
#line 6295 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 240:
#line 817 "p.y" /* yacc.c:1646  */
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6301 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 241:
#line 819 "p.y" /* yacc.c:1646  */
    { (yyval.coefdata).setCoef((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6307 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 242:
#line 823 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6313 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 243:
#line 825 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6319 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 244:
#line 829 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6325 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 245:
#line 831 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6331 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 246:
#line 835 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6337 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 247:
#line 837 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6343 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 248:
#line 841 "p.y" /* yacc.c:1646  */
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[-1].ival)-1,(yyval.linfo)); }
#line 6349 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 249:
#line 843 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].linfo)->add((yyvsp[0].ldata).lnum,(yyvsp[0].ldata).d,(yyvsp[0].ldata).matid); }
#line 6355 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 250:
#line 847 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-10].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-9].fval); (yyval.ldata).d[1] = (yyvsp[-8].fval); (yyval.ldata).d[2] = (yyvsp[-7].fval);
	  (yyval.ldata).d[3] = (yyvsp[-6].fval); (yyval.ldata).d[4] = (yyvsp[-5].fval); (yyval.ldata).d[5] = (yyvsp[-4].fval);
	  (yyval.ldata).d[6] = (yyvsp[-3].fval); (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval);
          (yyval.ldata).d[9] = 0;  (yyval.ldata).d[10] = 0; (yyval.ldata).d[11] = 0; }
#line 6366 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 251:
#line 854 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-12].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-11].fval); (yyval.ldata).d[1] = (yyvsp[-10].fval); (yyval.ldata).d[2] = (yyvsp[-9].fval);
          (yyval.ldata).d[3] = (yyvsp[-8].fval); (yyval.ldata).d[4] = (yyvsp[-7].fval); (yyval.ldata).d[5] = (yyvsp[-6].fval);
          (yyval.ldata).d[6] = (yyvsp[-5].fval); (yyval.ldata).d[7] = (yyvsp[-4].fval); (yyval.ldata).d[8] = (yyvsp[-3].fval);
          (yyval.ldata).d[9] = (yyvsp[-2].fval);(yyval.ldata).d[10]= (yyvsp[-1].fval);(yyval.ldata).d[11] = 0; }
#line 6377 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 252:
#line 861 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-13].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[-12].fval); (yyval.ldata).d[1] = (yyvsp[-11].fval); (yyval.ldata).d[2] = (yyvsp[-10].fval);
          (yyval.ldata).d[3] = (yyvsp[-9].fval); (yyval.ldata).d[4] = (yyvsp[-8].fval); (yyval.ldata).d[5] = (yyvsp[-7].fval);
          (yyval.ldata).d[6] = (yyvsp[-6].fval); (yyval.ldata).d[7] = (yyvsp[-5].fval); (yyval.ldata).d[8] = (yyvsp[-4].fval);
          (yyval.ldata).d[9] = (yyvsp[-3].fval);(yyval.ldata).d[10]= (yyvsp[-2].fval); (yyval.ldata).d[11] = (yyvsp[-1].fval); }
#line 6388 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 253:
#line 870 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).lnum = (yyvsp[-4].ival)-1;  (yyval.ldata).matid = (yyvsp[-3].ival)-1; (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval); }
#line 6394 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 255:
#line 875 "p.y" /* yacc.c:1646  */
    { geoSource->addLayMat((yyvsp[0].ldata).matid, (yyvsp[0].ldata).d); }
#line 6400 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 256:
#line 881 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-6].ival)-1; (yyval.ldata).d[0] = (yyvsp[-5].fval); (yyval.ldata).d[1] = (yyvsp[-4].fval); (yyval.ldata).d[2] = (yyvsp[-3].fval);
          (yyval.ldata).d[3] = (yyvsp[-2].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[-1].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
#line 6408 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 257:
#line 886 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-8].ival)-1; (yyval.ldata).d[0] = (yyvsp[-7].fval); (yyval.ldata).d[1] = (yyvsp[-6].fval); (yyval.ldata).d[2] = (yyvsp[-5].fval);
          (yyval.ldata).d[3] = (yyvsp[-4].fval); (yyval.ldata).d[4] = (yyvsp[-3].fval); (yyval.ldata).d[5] = (yyvsp[-2].fval); (yyval.ldata).d[6] = (yyvsp[-1].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
#line 6416 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 258:
#line 891 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-10].ival)-1; (yyval.ldata).d[0] = (yyvsp[-9].fval); (yyval.ldata).d[1] = (yyvsp[-8].fval); (yyval.ldata).d[2] = (yyvsp[-7].fval);
          (yyval.ldata).d[3] = (yyvsp[-6].fval); (yyval.ldata).d[4] = (yyvsp[-5].fval); (yyval.ldata).d[5] = (yyvsp[-4].fval); (yyval.ldata).d[6] = (yyvsp[-3].fval);
          (yyval.ldata).d[7] = (yyvsp[-2].fval); (yyval.ldata).d[8] = (yyvsp[-1].fval); (yyval.ldata).d[9] = 0; }
#line 6424 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 259:
#line 895 "p.y" /* yacc.c:1646  */
    { (yyval.ldata).matid = (yyvsp[-11].ival)-1; (yyval.ldata).d[0] = (yyvsp[-10].fval); (yyval.ldata).d[1] = (yyvsp[-9].fval); (yyval.ldata).d[2] = (yyvsp[-8].fval);
          (yyval.ldata).d[3] = (yyvsp[-7].fval); (yyval.ldata).d[4] = (yyvsp[-6].fval); (yyval.ldata).d[5] = (yyvsp[-5].fval); (yyval.ldata).d[6] = (yyvsp[-4].fval); 
          (yyval.ldata).d[7] = (yyvsp[-3].fval); (yyval.ldata).d[8] = (yyvsp[-2].fval); (yyval.ldata).d[9] = (yyvsp[-1].fval); }
#line 6432 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 261:
#line 902 "p.y" /* yacc.c:1646  */
    { domain->addDMass((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].fval)); }
#line 6438 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 262:
#line 904 "p.y" /* yacc.c:1646  */
    { domain->addDMass((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-1].fval),(yyvsp[-2].ival)-1); }
#line 6444 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 263:
#line 906 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modalDIMASS = true;
          domain->solInfo().reducedMassFile = (yyvsp[-1].strval); }
#line 6451 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 265:
#line 912 "p.y" /* yacc.c:1646  */
    { domain->setGravity((yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6457 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 267:
#line 917 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[-3].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[-2].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[-1].strval); }
#line 6465 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 268:
#line 921 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[-2].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval);}
#line 6472 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 269:
#line 924 "p.y" /* yacc.c:1646  */
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[-2].strval);
          domain->solInfo().nRestart = (yyvsp[-1].ival); }
#line 6479 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 270:
#line 929 "p.y" /* yacc.c:1646  */
    { geoSource->setControlFile((yyvsp[-1].strval));
         geoSource->setControlRoutine((char *) "controlObj");}
#line 6486 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 271:
#line 934 "p.y" /* yacc.c:1646  */
    { geoSource->setControlRoutine((yyvsp[-1].strval)); }
#line 6492 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 272:
#line 938 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 6499 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 273:
#line 943 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Actuators; }
          if(geoSource->setActuatorLocations((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)            return -1; }
#line 6507 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 274:
#line 949 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Usdf; }
          if(geoSource->setUsdfLocation((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)       return -1; }
#line 6516 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 275:
#line 956 "p.y" /* yacc.c:1646  */
    { geoSource->binaryInputControlLeft = true;
          (yyval.bclist) = new BCList; }
#line 6523 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 276:
#line 959 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Usdd; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 6529 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 277:
#line 961 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-4].ival); i<=(yyvsp[-2].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-1].ival)-1, 0., BCond::Usdd); (yyval.bclist)->add(bc); } }
#line 6535 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 278:
#line 963 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); i+=(yyvsp[-2].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-1].ival)-1, 0., BCond::Usdd); (yyval.bclist)->add(bc); } }
#line 6541 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 279:
#line 965 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Usdd;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
#line 6551 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 280:
#line 973 "p.y" /* yacc.c:1646  */
    { numColumns = 3; }
#line 6557 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 281:
#line 975 "p.y" /* yacc.c:1646  */
    { numColumns = 6; }
#line 6563 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 282:
#line 977 "p.y" /* yacc.c:1646  */
    { numColumns = 3; geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6569 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 283:
#line 979 "p.y" /* yacc.c:1646  */
    { numColumns = 6; geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6575 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 284:
#line 981 "p.y" /* yacc.c:1646  */
    { numColumns = 3; domain->outFlag = (yyvsp[-1].ival); }
#line 6581 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 285:
#line 983 "p.y" /* yacc.c:1646  */
    { numColumns = 6; domain->outFlag = (yyvsp[-1].ival); }
#line 6587 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 286:
#line 985 "p.y" /* yacc.c:1646  */
    { numColumns = 3; domain->outFlag = (yyvsp[-2].ival); geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6593 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 287:
#line 987 "p.y" /* yacc.c:1646  */
    { numColumns = 6; domain->outFlag = (yyvsp[-2].ival); geoSource->setOutLimit((yyvsp[-1].ival)); }
#line 6599 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 288:
#line 989 "p.y" /* yacc.c:1646  */
    { numColumns = 3; geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval); }
#line 6605 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 289:
#line 991 "p.y" /* yacc.c:1646  */
    { numColumns = 6; geoSource->getCheckFileInfo()->outputExt = (yyvsp[-1].strval); }
#line 6611 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 290:
#line 993 "p.y" /* yacc.c:1646  */
    { (yyvsp[-1].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[-1].oinfo)); }
#line 6617 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 291:
#line 997 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6623 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 292:
#line 999 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-4].ival); (yyval.oinfo).width = (yyvsp[-3].ival); (yyval.oinfo).precision = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6629 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 293:
#line 1001 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6635 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 294:
#line 1003 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); 
          if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1;}
#line 6642 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 295:
#line 1006 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-5].ival); (yyval.oinfo).width = (yyvsp[-4].ival); (yyval.oinfo).precision = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6648 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 296:
#line 1008 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[-6].ival); (yyval.oinfo).width = (yyvsp[-5].ival); (yyval.oinfo).precision = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6654 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 297:
#line 1010 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6660 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 298:
#line 1012 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-4].ival); (yyval.oinfo).width = (yyvsp[-3].ival); (yyval.oinfo).precision = (yyvsp[-2].ival); (yyval.oinfo).filename = (yyvsp[-1].strval); (yyval.oinfo).interval = (yyvsp[0].ival); }
#line 6666 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 299:
#line 1014 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6672 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 300:
#line 1016 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6678 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 301:
#line 1018 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-5].ival); (yyval.oinfo).width = (yyvsp[-4].ival); (yyval.oinfo).precision = (yyvsp[-3].ival); (yyval.oinfo).filename = (yyvsp[-2].strval); (yyval.oinfo).interval = (yyvsp[-1].ival); (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6684 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 302:
#line 1020 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[-6].ival); (yyval.oinfo).width = (yyvsp[-5].ival); (yyval.oinfo).precision = (yyvsp[-4].ival); (yyval.oinfo).filename = (yyvsp[-3].strval); (yyval.oinfo).interval = (yyvsp[-2].ival); if ((yyvsp[-1].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[0].ival); else (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6690 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 303:
#line 1023 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).nodeNumber = (yyvsp[0].ival)-1; }
#line 6696 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 304:
#line 1025 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).surface = (yyvsp[0].ival); }
#line 6702 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 305:
#line 1027 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).str_therm_option = (yyvsp[0].ival); }
#line 6708 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 306:
#line 1029 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ylayer = (yyvsp[-1].fval); (yyval.oinfo).zlayer = (yyvsp[0].fval); }
#line 6714 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 307:
#line 1031 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).averageFlg = (yyvsp[0].ival); }
#line 6720 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 308:
#line 1033 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).complexouttype = (yyvsp[0].ival); }
#line 6726 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 309:
#line 1035 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).complexouttype = (yyvsp[-1].ival); (yyval.oinfo).ncomplexout = (yyvsp[0].ival); }
#line 6732 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 310:
#line 1037 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).angularouttype = (yyvsp[0].ival); }
#line 6738 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 311:
#line 1039 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rotvecouttype = (yyvsp[0].ival); }
#line 6744 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 312:
#line 1041 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rotvecouttype = OutputInfo::Linear; }
#line 6750 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 313:
#line 1043 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).rescaling = bool((yyvsp[0].ival)); }
#line 6756 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 314:
#line 1045 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ndtype = (yyvsp[0].ival); }
#line 6762 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 315:
#line 1047 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).ndtype = (yyvsp[-1].ival); sfem->setnsamp_out((yyvsp[0].ival)); }
#line 6768 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 316:
#line 1049 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).oframe = (OutputInfo::FrameType) (yyvsp[0].ival); }
#line 6774 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 317:
#line 1051 "p.y" /* yacc.c:1646  */
    { (yyval.oinfo).matlab = true; }
#line 6780 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 318:
#line 1053 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xmatrixname = (yyvsp[0].strval); }
#line 6786 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 319:
#line 1055 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qmatrixname = (yyvsp[0].strval); }
#line 6792 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 320:
#line 1057 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rmatrixname = (yyvsp[0].strval); }
#line 6798 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 321:
#line 1059 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenvaluename = (yyvsp[0].strval); }
#line 6804 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 323:
#line 1064 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; }
#line 6811 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 324:
#line 1067 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[-1].ival);}
#line 6818 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 325:
#line 1070 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qrfactorization = (yyvsp[-1].ival);}
#line 6824 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 326:
#line 1072 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nEig = (yyvsp[-1].ival); }
#line 6830 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 327:
#line 1074 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;}
#line 6836 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 328:
#line 1076 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setSubSpaceInfo((yyvsp[-3].ival),(yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 6842 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 329:
#line 1078 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subspaceSize = (yyvsp[-1].ival);}
#line 6848 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 330:
#line 1080 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolEig = (yyvsp[-1].fval); }
#line 6854 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 331:
#line 1082 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolJac = (yyvsp[-1].fval); }
#line 6860 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 332:
#line 1084 "p.y" /* yacc.c:1646  */
    { domain->solInfo().explicitK = true; }
#line 6866 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 333:
#line 1086 "p.y" /* yacc.c:1646  */
    { geoSource->setShift((yyvsp[-1].fval)); }
#line 6872 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 334:
#line 1088 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; }
#line 6878 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 335:
#line 1090 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[-1].strval); }
#line 6885 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 336:
#line 1093 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[-2].strval); 
          domain->solInfo().arpack_mode = (yyvsp[-1].ival); }
#line 6893 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 337:
#line 1097 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[-2].fval), int((yyvsp[-1].fval))); }
#line 6900 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 338:
#line 1100 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival));}
#line 6907 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 339:
#line 1103 "p.y" /* yacc.c:1646  */
    { domain->solInfo().filtereig = bool((yyvsp[-1].ival)); }
#line 6913 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 340:
#line 1105 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverSubType = (yyvsp[-1].ival); }
#line 6919 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 341:
#line 1107 "p.y" /* yacc.c:1646  */
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;}
#line 6926 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 342:
#line 1110 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxitEig = (yyvsp[-1].ival); }
#line 6932 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 343:
#line 1112 "p.y" /* yacc.c:1646  */
    { domain->solInfo().test_ulrich = true; }
#line 6938 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 344:
#line 1114 "p.y" /* yacc.c:1646  */
    { domain->solInfo().addedMass = (yyvsp[-1].ival); }
#line 6944 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 345:
#line 1118 "p.y" /* yacc.c:1646  */
    { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = (yyvsp[-1].strval);
          domain->solInfo().printMatLabExit = true; }
#line 6952 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 346:
#line 1122 "p.y" /* yacc.c:1646  */
    { domain->solInfo().printMatLab = true;
          domain->solInfo().printMatLabFile = (yyvsp[-2].strval); 
          domain->solInfo().printMatLabExit = true; }
#line 6960 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 347:
#line 1128 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elementDeletion = true; }
#line 6966 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 349:
#line 1133 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[-1].fval); domain->solInfo().deleteElements.insert(std::pair<int,double>((yyvsp[0].ival)-1,(yyvsp[-1].fval))); }
#line 6972 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 350:
#line 1135 "p.y" /* yacc.c:1646  */
    { domain->solInfo().deleteElements.insert(std::pair<int,double>((yyvsp[0].ival)-1,(yyvsp[-1].fval))); }
#line 6978 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 351:
#line 1139 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sloshing = 1; }
#line 6984 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 352:
#line 1141 "p.y" /* yacc.c:1646  */
    { domain->setGravitySloshing((yyvsp[-1].fval)); }
#line 6990 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 353:
#line 1145 "p.y" /* yacc.c:1646  */
    { domain->solInfo().massFlag = 1; }
#line 6996 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 354:
#line 1147 "p.y" /* yacc.c:1646  */
    { domain->solInfo().massFlag = 1;
          domain->solInfo().massFile = std::string((yyvsp[-1].strval)); }
#line 7003 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 355:
#line 1152 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7010 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 356:
#line 1155 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);}
#line 7016 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 357:
#line 1159 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Top); }
#line 7022 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 358:
#line 1163 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal_id.push_back((yyvsp[0].ival)); }
#line 7028 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 359:
#line 1165 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal_id.push_back((yyvsp[0].ival)); }
#line 7034 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 360:
#line 1167 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contact_modal_id.push_back((yyvsp[0].ival)); }
#line 7040 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 361:
#line 1170 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Dynamic); }
#line 7046 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 365:
#line 1175 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
#line 7052 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 366:
#line 1177 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; }
#line 7058 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 367:
#line 1179 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-1].ival); }
#line 7064 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 368:
#line 1181 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-1].ival); }
#line 7070 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 369:
#line 1183 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stable = (yyvsp[-5].ival);
          domain->solInfo().stable_cfl = (yyvsp[-4].fval);
          domain->solInfo().stable_tol = (yyvsp[-3].fval);
          domain->solInfo().stable_maxit = (yyvsp[-2].ival);
          domain->solInfo().stable_freq = (yyvsp[-1].ival);
        }
#line 7081 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 370:
#line 1190 "p.y" /* yacc.c:1646  */
    { domain->solInfo().iacc_switch = bool((yyvsp[-1].ival)); }
#line 7087 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 371:
#line 1192 "p.y" /* yacc.c:1646  */
    { domain->solInfo().zeroRot = bool((yyvsp[-1].ival)); }
#line 7093 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 372:
#line 1194 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_secondary = true; }
#line 7099 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 373:
#line 1196 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-1].ival); }
#line 7105 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 374:
#line 1198 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-3].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-2].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-1].fval); }
#line 7113 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 375:
#line 1202 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-4].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-3].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-2].fval);
          domain->solInfo().tdenforceInitia = (yyvsp[-1].fval); }
#line 7122 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 376:
#line 1207 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tdenforceFlag = (yyvsp[-5].ival);
          domain->solInfo().tdenforceMaxItr = (yyvsp[-4].ival);
          domain->solInfo().tdenforceTolAbs = (yyvsp[-3].fval);
          domain->solInfo().tdenforceInitia = (yyvsp[-2].fval);
          domain->solInfo().tdenforceFinal = (yyvsp[-1].fval); }
#line 7132 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 377:
#line 1213 "p.y" /* yacc.c:1646  */
    { domain->solInfo().check_energy_balance = true; }
#line 7138 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 378:
#line 1215 "p.y" /* yacc.c:1646  */
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[-2].fval); 
          domain->solInfo().epsilon2 = (yyvsp[-1].fval); }
#line 7146 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 379:
#line 1219 "p.y" /* yacc.c:1646  */
    { domain->solInfo().quasistatic = bool((yyvsp[-1].ival)); }
#line 7152 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 380:
#line 1223 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ConwepOnOff = true;
          BlastLoading::InputFileData = (yyvsp[-1].blastData); }
#line 7159 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 381:
#line 1228 "p.y" /* yacc.c:1646  */
    { // Note: chargeWeight must be entered in the units of mass of the problem, not units of force.
          (yyval.blastData).ExplosivePosition[0] = (yyvsp[-4].fval);
          (yyval.blastData).ExplosivePosition[1] = (yyvsp[-3].fval);
          (yyval.blastData).ExplosivePosition[2] = (yyvsp[-2].fval);
          (yyval.blastData).ExplosiveDetonationTime = (yyvsp[0].fval);
          (yyval.blastData).BlastType = BlastLoading::BlastData::AirBurst; // ($5 == 0 ? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst);
          (yyval.blastData).ScaleLength = 1.0;
          (yyval.blastData).ScaleTime = 1.0;
          (yyval.blastData).ScaleMass = 1.0;
          (yyval.blastData).ExplosiveWeight = (yyvsp[-1].fval)*2.2; // The 2.2 factor is to convert from kilograms to pounds force.
          (yyval.blastData).ExplosiveWeightCubeRoot = pow((yyval.blastData).ExplosiveWeight,1.0/3.0);
        }
#line 7176 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 382:
#line 1243 "p.y" /* yacc.c:1646  */
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; }
#line 7182 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 384:
#line 1246 "p.y" /* yacc.c:1646  */
    { domain->solInfo().acoustic = true; }
#line 7188 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 386:
#line 1251 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-1].fval),(yyvsp[0].fval)); }
#line 7194 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 387:
#line 1253 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-3].fval),(yyvsp[-2].fval),(yyvsp[-1].fval),(yyvsp[0].fval)); }
#line 7200 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 388:
#line 1255 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[0].fval)); }
#line 7206 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 389:
#line 1257 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[-2].fval),(yyvsp[-1].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[0].fval); }
#line 7214 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7227 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 391:
#line 1274 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; }
#line 7234 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 394:
#line 1279 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back(0); }
#line 7240 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 395:
#line 1281 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modal = true; domain->solInfo().modal_id.push_back((yyvsp[-1].ival)); }
#line 7246 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 396:
#line 1285 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setQuasistaticInfo((yyvsp[-3].fval), 0, (yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7252 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 397:
#line 1287 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setQuasistaticInfo((yyvsp[-4].fval), 0, (yyvsp[-3].fval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 7258 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 398:
#line 1291 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].ival)); }
#line 7265 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 399:
#line 1301 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-1].ival)); 
          domain->solInfo().isCollocated = 0; }
#line 7272 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7306 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7330 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 403:
#line 1351 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAero((yyvsp[-2].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[-1].fval);
        }
#line 7339 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 404:
#line 1356 "p.y" /* yacc.c:1646  */
    { domain->solInfo().isCollocated = (yyvsp[-1].ival); }
#line 7345 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 405:
#line 1358 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subcycle = (yyvsp[-1].ival); }
#line 7351 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 406:
#line 1360 "p.y" /* yacc.c:1646  */
    { geoSource->setMatch((yyvsp[-1].strval)); }
#line 7357 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 407:
#line 1362 "p.y" /* yacc.c:1646  */
    {}
#line 7363 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 408:
#line 1366 "p.y" /* yacc.c:1646  */
    {}
#line 7369 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 409:
#line 1368 "p.y" /* yacc.c:1646  */
    { domain->AddAeroEmbedSurfaceId((yyvsp[0].ival)); }
#line 7375 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 410:
#line 1372 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setAeroHeat((yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 7381 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 411:
#line 1376 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setThermoh(1); }
#line 7387 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 412:
#line 1380 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setThermoe(1); }
#line 7393 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 413:
#line 1384 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setModeDecomp(1); }
#line 7399 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 414:
#line 1386 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setModeDecomp(1, (yyvsp[-1].ival)); }
#line 7405 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 415:
#line 1390 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hzemFlag=1; }
#line 7411 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 416:
#line 1394 "p.y" /* yacc.c:1646  */
    { domain->solInfo().slzemFlag=1; }
#line 7417 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 417:
#line 1398 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setTrbm((yyvsp[-1].fval)); }
#line 7423 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 418:
#line 1402 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-2].fval),(yyvsp[-1].fval)); }
#line 7429 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 419:
#line 1404 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-1].fval)); }
#line 7435 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 420:
#line 1406 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm(); }
#line 7441 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 421:
#line 1408 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-3].fval),(yyvsp[-2].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[-1].ival)); }
#line 7448 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 422:
#line 1411 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-2].fval),(yyvsp[-1].fval));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-6].fval); grbm_ref[1] = (yyvsp[-5].fval); grbm_ref[2] = (yyvsp[-4].fval);
        }
#line 7457 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 423:
#line 1416 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-1].fval));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-5].fval); grbm_ref[1] = (yyvsp[-4].fval); grbm_ref[2] = (yyvsp[-3].fval);
        }
#line 7466 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 424:
#line 1421 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm();
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-3].fval); grbm_ref[1] = (yyvsp[-2].fval); grbm_ref[2] = (yyvsp[-1].fval);
        }
#line 7475 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 425:
#line 1426 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGrbm((yyvsp[-3].fval),(yyvsp[-2].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[-1].ival));
          std::vector<double> &grbm_ref = domain->solInfo().grbm_ref;
          grbm_ref.resize(3); grbm_ref[0] = (yyvsp[-7].fval); grbm_ref[1] = (yyvsp[-6].fval); grbm_ref[2] = (yyvsp[-5].fval);
        }
#line 7485 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 426:
#line 1434 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modeFilterFlag = (yyvsp[-1].ival); }
#line 7491 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 427:
#line 1436 "p.y" /* yacc.c:1646  */
    { domain->solInfo().modeFilterFlag = 1; }
#line 7497 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 428:
#line 1440 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter((yyvsp[-1].ival)); }
#line 7503 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 429:
#line 1442 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter(1); }
#line 7509 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 430:
#line 1444 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useRbmFilter((yyvsp[-2].ival));
          domain->solInfo().filterQ = (yyvsp[-1].ival); }
#line 7516 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7528 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7540 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 434:
#line 1468 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hzemFilterFlag=1; }
#line 7546 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 435:
#line 1472 "p.y" /* yacc.c:1646  */
    { domain->solInfo().slzemFilterFlag=1; }
#line 7552 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 436:
#line 1476 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setTimes((yyvsp[-1].fval),(yyvsp[-2].fval),(yyvsp[-3].fval)); }
#line 7558 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 437:
#line 1480 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[-2].ival),(yyvsp[-1].ival),1);
        }
#line 7567 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 438:
#line 1486 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[-3].ival),(yyvsp[-2].ival),(yyvsp[-1].ival));
        }
#line 7576 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 439:
#line 1491 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[-4].ival),(yyvsp[-3].ival),(yyvsp[-2].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        }
#line 7587 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 442:
#line 1504 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaNoForce = true; }
#line 7593 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 443:
#line 1506 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[0].ival); }
#line 7599 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 444:
#line 1508 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaLocalBasisImprovement = 1; }
#line 7605 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 445:
#line 1510 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaTimeReversible = true; }
#line 7611 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 446:
#line 1512 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaRemoteCoarse = true; }
#line 7617 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 447:
#line 1514 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaProjTol = (yyvsp[0].fval); }
#line 7623 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 448:
#line 1516 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaReadInitSeed = true; }
#line 7629 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 449:
#line 1518 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpCvgRatio = 0.0; }
#line 7635 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 450:
#line 1520 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[0].fval); }
#line 7641 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 451:
#line 1522 "p.y" /* yacc.c:1646  */
    { domain->solInfo().pitaJumpMagnOutput = true; }
#line 7647 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 452:
#line 1526 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-3].ival) == 1) domain->solInfo().setDamping((yyvsp[-2].fval),(yyvsp[-1].fval));
          else return -1; // only RAYDAMP is allowed here
        }
#line 7655 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 453:
#line 1530 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-4].ival) == 1) {
            domain->solInfo().setDamping((yyvsp[-3].fval),(yyvsp[-2].fval));
            domain->solInfo().mtypeDamp = (int)(yyvsp[-1].ival);
          }
          else return -1;
        }
#line 7666 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 454:
#line 1537 "p.y" /* yacc.c:1646  */
    { if(geoSource->setModalDamping((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true;
        }
#line 7674 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 455:
#line 1543 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 7680 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 456:
#line 1547 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7690 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 457:
#line 1553 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-2].fval), (yyvsp[-1].fval), 0.0);
        }
#line 7700 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 458:
#line 1561 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7710 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 459:
#line 1567 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[-1].ival),0.0,0.0,0.0);
        }
#line 7719 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 461:
#line 1575 "p.y" /* yacc.c:1646  */
    {
           domain->setWaveDirections((yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7727 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 462:
#line 1581 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 7737 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 463:
#line 1589 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7743 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 464:
#line 1591 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7749 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 465:
#line 1593 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
#line 7755 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 466:
#line 1595 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
#line 7761 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 467:
#line 1597 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; }
#line 7771 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7784 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 470:
#line 1623 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[0].bclist); }
#line 7790 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 471:
#line 1627 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7796 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 472:
#line 1629 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 7802 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 473:
#line 1633 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 7808 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 474:
#line 1635 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-1].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; }
#line 7814 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7857 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 476:
#line 1680 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[0].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[0].bclist)->d)[ii]); }
#line 7865 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 477:
#line 1684 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); 
          for (int ii = 0; ii < (yyvsp[0].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[0].bclist)->d)[ii]); }
#line 7873 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 478:
#line 1690 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList;
          for(int i=0; i<(yyvsp[-1].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[-1].nl).nd[i],10,0.0); } }
#line 7881 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 479:
#line 1696 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;}
#line 7887 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 480:
#line 1698 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-3].bclist); BCond bc; bc.nnum = (yyvsp[-2].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-1].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); }
#line 7894 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 481:
#line 1701 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-4].ival); i<=(yyvsp[-2].ival); ++i) { BCond bc; bc.setData(i-1, 6, (yyvsp[-1].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } }
#line 7900 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 482:
#line 1703 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); i+=(yyvsp[-2].ival)) { BCond bc; bc.setData(i-1, 6, (yyvsp[-1].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } }
#line 7906 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 483:
#line 1705 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[-2].ival)-1;
          surf_bc[0].val = (yyvsp[-1].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
#line 7917 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 484:
#line 1714 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7923 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 485:
#line 1716 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7929 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 486:
#line 1718 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-3].bclist); BCond bc; bc.nnum = (yyvsp[-2].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-1].fval); bc.type = BCond::Flux; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); }
#line 7936 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7949 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 488:
#line 1732 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7955 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 489:
#line 1734 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7961 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 490:
#line 1736 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-5].bclist); BCond bc; bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval); bc.type = BCond::Convection; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); }
#line 7968 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 7981 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 492:
#line 1750 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 7987 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 493:
#line 1752 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 7993 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 494:
#line 1754 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-5].bclist); BCond bc; bc.nnum = (yyvsp[-4].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[-3].fval)*(yyvsp[-2].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval)*(yyvsp[-1].fval); bc.type = BCond::Radiation; (yyval.bclist)->add(bc); }
#line 8000 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 498:
#line 1766 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8006 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 499:
#line 1768 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)); }
#line 8012 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 500:
#line 1770 "p.y" /* yacc.c:1646  */
    { domain->addSommer(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8018 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 503:
#line 1778 "p.y" /* yacc.c:1646  */
    { domain->addSommerElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        }
#line 8026 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 504:
#line 1784 "p.y" /* yacc.c:1646  */
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[0].ival)-1; }
#line 8032 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 505:
#line 1786 "p.y" /* yacc.c:1646  */
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[0].ival)-1; (yyval.nl).num++; }
#line 8039 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 509:
#line 1798 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8046 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 510:
#line 1801 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)); }
#line 8053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 511:
#line 1804 "p.y" /* yacc.c:1646  */
    { domain->addScatter(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1)); }
#line 8060 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 514:
#line 1813 "p.y" /* yacc.c:1646  */
    { domain->addScatterElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
          domain->addNeumElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); }
#line 8067 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 517:
#line 1822 "p.y" /* yacc.c:1646  */
    { domain->addNeumElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); }
#line 8073 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 520:
#line 1830 "p.y" /* yacc.c:1646  */
    { domain->addWetElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; }
#line 8081 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 523:
#line 1840 "p.y" /* yacc.c:1646  */
    { domain->addScatterElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);}
#line 8087 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 524:
#line 1844 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 8093 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 525:
#line 1848 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8099 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 526:
#line 1850 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8105 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 527:
#line 1854 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[0].bclist); }
#line 8111 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 528:
#line 1858 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[0].bclist); }
#line 8117 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 529:
#line 1860 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[0].bclist); }
#line 8123 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 530:
#line 1864 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDARBFlag = (yyvsp[-1].fval);}
#line 8129 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 532:
#line 1869 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDDNBVal = (yyvsp[-1].fval);}
#line 8135 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 534:
#line 1874 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ATDROBVal = (yyvsp[-3].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[-2].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[-1].fval);}
#line 8143 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 536:
#line 1881 "p.y" /* yacc.c:1646  */
    { domain->setFFP((yyvsp[-1].ival)); }
#line 8149 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 537:
#line 1883 "p.y" /* yacc.c:1646  */
    { domain->setFFP((yyvsp[-2].ival),(yyvsp[-1].ival)); }
#line 8155 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 538:
#line 1887 "p.y" /* yacc.c:1646  */
    {
           domain->setFFP((yyvsp[-2].ival));
        }
#line 8163 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 539:
#line 1893 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-4].ival)] = ModalParams(ModalParams::Eigen, (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8169 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 540:
#line 1895 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-4].ival)] = ModalParams((yyvsp[-3].mpt), (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8175 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 541:
#line 1897 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-5].ival)] = ModalParams(ModalParams::Eigen, (yyvsp[-3].strval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 8181 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 542:
#line 1899 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[(yyvsp[-5].ival)] = ModalParams((yyvsp[-4].mpt), (yyvsp[-3].strval), (yyvsp[-2].ival), (yyvsp[-1].fval)); }
#line 8187 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 543:
#line 1903 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, (yyvsp[-1].strval)); }
#line 8193 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 544:
#line 1905 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInModes[0] = ModalParams(ModalParams::Undefined, (yyvsp[-2].strval), (yyvsp[-1].ival)); }
#line 8199 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 547:
#line 1909 "p.y" /* yacc.c:1646  */
    { domain->solInfo().adjointMap[(OutputInfo::Type)(yyvsp[-1].ival)] = domain->solInfo().readInAdjointROB.size();
          domain->solInfo().readInAdjointROB.push_back((yyvsp[-3].strval));
          domain->solInfo().maxSizeAdjointBasis.push_back((yyvsp[-2].ival)); }
#line 8207 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 548:
#line 1915 "p.y" /* yacc.c:1646  */
    { }
#line 8213 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 549:
#line 1917 "p.y" /* yacc.c:1646  */
    { domain->solInfo().zeroInitialDisp = 1; }
#line 8219 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 550:
#line 1919 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0)  return -1; }
#line 8226 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 551:
#line 1922 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
#line 8234 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 552:
#line 1926 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().idis_modal_id = (yyvsp[-2].ival); }
#line 8243 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 553:
#line 1933 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[-1].fval);  }
#line 8249 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 554:
#line 1935 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; amplitude = 1.0; }
#line 8255 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8270 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 557:
#line 1959 "p.y" /* yacc.c:1646  */
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); }
#line 8292 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 558:
#line 1962 "p.y" /* yacc.c:1646  */
    { domain->solInfo().buckling = 1; }
#line 8298 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 559:
#line 1967 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[-1].ival); }
#line 8304 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8318 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 561:
#line 1982 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[-1].ival); }
#line 8324 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8338 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 563:
#line 1996 "p.y" /* yacc.c:1646  */
    { }
#line 8344 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 564:
#line 1998 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8351 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 565:
#line 2001 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
#line 8359 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 566:
#line 2005 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true;
          domain->solInfo().ivel_modal_id = (yyvsp[-2].ival); }
#line 8368 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 567:
#line 2012 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8375 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 568:
#line 2015 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDisModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1;
          domain->solInfo().modalCalled = true; }
#line 8383 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 569:
#line 2021 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setGEPS();
          for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Etemperatures;
          if(geoSource->setIDis6((yyvsp[0].bclist)->n,(yyvsp[0].bclist)->d) < 0) return -1; }
#line 8391 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 570:
#line 2027 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; }
#line 8397 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 571:
#line 2029 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList((yyvsp[-1].ival)); }
#line 8403 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 572:
#line 2031 "p.y" /* yacc.c:1646  */
    { (yyvsp[0].bcval).type = BCond::Forces; (yyvsp[0].bcval).loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8409 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 573:
#line 2033 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } }
#line 8415 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 574:
#line 2035 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } }
#line 8421 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 575:
#line 2037 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-6].ival); i<=(yyvsp[-4].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[-1].ival)); (yyval.bclist)->add(bc); } }
#line 8427 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 576:
#line 2039 "p.y" /* yacc.c:1646  */
    { for(int i=(yyvsp[-8].ival); i<=(yyvsp[-6].ival); i+=(yyvsp[-4].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[-3].ival)); (yyval.bclist)->add(bc); } }
#line 8433 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 577:
#line 2041 "p.y" /* yacc.c:1646  */
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[0].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; }
#line 8444 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 578:
#line 2050 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8451 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 579:
#line 2053 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Forces; (yyvsp[0].bclist)->d[i].loadsetid = (yyvsp[-4].ival); }
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8458 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 580:
#line 2056 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) (yyvsp[0].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8465 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 581:
#line 2059 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].bclist)->n; ++i) { (yyvsp[0].bclist)->d[i].type = BCond::Forces; (yyvsp[0].bclist)->d[i].loadsetid = (yyvsp[-5].ival); }
          if(geoSource->setNeumanModal((yyvsp[0].bclist)->n, (yyvsp[0].bclist)->d) < 0) return -1; }
#line 8472 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 582:
#line 2064 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8478 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 583:
#line 2066 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8484 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 584:
#line 2068 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); }}
#line 8490 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 585:
#line 2070 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-6].bclist); for(int i=(yyvsp[-5].ival); i<=(yyvsp[-3].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); }}
#line 8496 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 586:
#line 2072 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); } }
#line 8502 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 587:
#line 2074 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-8].bclist); for(int i=(yyvsp[-7].ival); i<=(yyvsp[-5].ival); i+=(yyvsp[-3].ival)) { BCond bc; bc.setData(i-1, (yyvsp[-2].ival)-1, (yyvsp[-1].fval)); (yyval.bclist)->add(bc); } }
#line 8508 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 588:
#line 2078 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8514 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 589:
#line 2080 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8520 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 590:
#line 2084 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8526 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 591:
#line 2086 "p.y" /* yacc.c:1646  */
    { (yyval.bclist) = (yyvsp[-1].bclist); (yyval.bclist)->add((yyvsp[0].bcval)); }
#line 8532 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 594:
#line 2094 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMTT((yyval.ymtt));}
#line 8538 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 595:
#line 2096 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8544 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 596:
#line 2098 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMTT((yyval.ymtt));}
#line 8550 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 599:
#line 2106 "p.y" /* yacc.c:1646  */
    { (yyval.ctett) = new MFTTData((yyvsp[-4].ival)); (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addCTETT((yyval.ctett));}
#line 8556 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 600:
#line 2108 "p.y" /* yacc.c:1646  */
    { (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8562 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 601:
#line 2110 "p.y" /* yacc.c:1646  */
    { (yyval.ctett) = new MFTTData((yyvsp[-4].ival)); (yyval.ctett)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addCTETT((yyval.ctett));}
#line 8568 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 604:
#line 2118 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSS1DT((yyval.ymtt));}
#line 8574 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 605:
#line 2120 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8580 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 606:
#line 2122 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSS1DT((yyval.ymtt));}
#line 8586 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 609:
#line 2130 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-6].ival), (yyvsp[-4].dlist), true); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8592 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 610:
#line 2132 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-8].ival), (yyvsp[-4].dlist), bool((yyvsp[-6].ival))); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8598 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 611:
#line 2134 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); }
#line 8604 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 612:
#line 2136 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-6].ival), (yyvsp[-4].dlist), true); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8610 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 613:
#line 2138 "p.y" /* yacc.c:1646  */
    { (yyval.ss2dt) = new SS2DTData((yyvsp[-8].ival), (yyvsp[-4].dlist), bool((yyvsp[-6].ival))); (yyval.ss2dt)->add((yyvsp[-2].fval), (yyvsp[-1].dlist)); domain->addSS2DT((yyval.ss2dt)); }
#line 8616 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 616:
#line 2146 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSST((yyval.ymtt));}
#line 8622 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 617:
#line 2148 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8628 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 618:
#line 2150 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSST((yyval.ymtt));}
#line 8634 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 621:
#line 2158 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSSRT((yyval.ymtt));}
#line 8640 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 622:
#line 2160 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8646 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 623:
#line 2162 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYSSRT((yyval.ymtt));}
#line 8652 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 626:
#line 2170 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMST((yyval.ymtt));}
#line 8658 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 627:
#line 2172 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8664 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 628:
#line 2174 "p.y" /* yacc.c:1646  */
    { (yyval.ymtt) = new MFTTData((yyvsp[-4].ival)); (yyval.ymtt)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addYMST((yyval.ymtt));}
#line 8670 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 631:
#line 2182 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft) = new MFTTData((yyvsp[-4].ival)); (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSDETAFT((yyval.sdetaft));}
#line 8676 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 632:
#line 2184 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 8682 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 633:
#line 2186 "p.y" /* yacc.c:1646  */
    { (yyval.sdetaft) = new MFTTData((yyvsp[-4].ival)); (yyval.sdetaft)->add((yyvsp[-2].fval), (yyvsp[-1].fval)); domain->addSDETAFT((yyval.sdetaft));}
#line 8688 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8700 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8712 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8724 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 641:
#line 2222 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xLMPCFactor = (yyvsp[-3].fval);
          domain->solInfo().yLMPCFactor = (yyvsp[-2].fval);
          domain->solInfo().zLMPCFactor = (yyvsp[-1].fval); }
#line 8732 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 642:
#line 2228 "p.y" /* yacc.c:1646  */
    { domain->solInfo().localDualBasisSize.push_back((yyvsp[-1].ival));
          domain->solInfo().modalLMPC = true; }
#line 8739 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 643:
#line 2231 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackROMLMPCVec((yyvsp[-1].fval)); }
#line 8745 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 644:
#line 2235 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = (yyvsp[-1].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[0].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
#line 8753 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 645:
#line 2239 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons)->addterm((yyvsp[0].mpcterm)); }
#line 8759 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 646:
#line 2241 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = (yyvsp[-1].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[0].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
#line 8767 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 647:
#line 2247 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8774 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 648:
#line 2250 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-2].ival), (yyvsp[-1].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8781 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 649:
#line 2253 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-4].ival), (yyvsp[-3].fval));
          (yyval.lmpcons)->type = (yyvsp[-1].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8789 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 650:
#line 2257 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-3].ival), (yyvsp[-2].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[-1].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8798 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 651:
#line 2262 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-5].ival), (yyvsp[-4].fval));
          (yyval.lmpcons)->type = (yyvsp[-2].ival);
          (yyval.lmpcons)->lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[-1].copt).penalty;
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
#line 8808 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8829 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 655:
#line 2293 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].cxbcval).nnum,(yyvsp[-1].cxbcval).reval,(yyvsp[-1].cxbcval).imval,(yyvsp[0].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
#line 8835 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 656:
#line 2295 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons)->addterm((yyvsp[0].mpcterm)); }
#line 8841 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 657:
#line 2297 "p.y" /* yacc.c:1646  */
    { (yyval.lmpcons) = new LMPCons((yyvsp[-1].cxbcval).nnum,(yyvsp[-1].cxbcval).reval,(yyvsp[-1].cxbcval).imval,(yyvsp[0].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
#line 8847 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 658:
#line 2301 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-4].ival); (yyval.cxbcval).reval=(yyvsp[-2].fval); (yyval.cxbcval).imval=(yyvsp[-1].fval); }
#line 8853 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 659:
#line 2303 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-2].ival); (yyval.cxbcval).reval=(yyvsp[-1].fval); (yyval.cxbcval).imval=0.0; }
#line 8859 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 660:
#line 2305 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum=(yyvsp[-1].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; }
#line 8865 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8877 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8889 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 663:
#line 2327 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 8895 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 664:
#line 2329 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].cxbclist)->n; ++i) (yyvsp[0].cxbclist)->d[i].loadsetid = (yyvsp[-2].ival);
          (yyval.cxbclist) = (yyvsp[0].cxbclist); }
#line 8902 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 665:
#line 2334 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[0].cxbcval)); }
#line 8908 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 666:
#line 2336 "p.y" /* yacc.c:1646  */
    { (yyval.cxbclist) = (yyvsp[-1].cxbclist); (yyval.cxbclist)->add((yyvsp[0].cxbcval)); }
#line 8914 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 669:
#line 2344 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
	  sp.A = (yyvsp[-14].fval);  sp.E = (yyvsp[-13].fval);  sp.nu  = (yyvsp[-12].fval);  sp.rho = (yyvsp[-11].fval);
          sp.c = (yyvsp[-10].fval);  sp.k = (yyvsp[-9].fval);  sp.eh  = (yyvsp[-8].fval);  sp.P   = (yyvsp[-7].fval);  sp.Ta  = (yyvsp[-6].fval); 
          sp.Q = (yyvsp[-5].fval); sp.W = (yyvsp[-4].fval); sp.Ixx = (yyvsp[-3].fval); sp.Iyy = (yyvsp[-2].fval); sp.Izz = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-15].ival)-1, sp );
        }
#line 8925 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8952 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8964 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 8992 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 673:
#line 2406 "p.y" /* yacc.c:1646  */
    { StructProp sp;
          sp.A = (yyvsp[-7].fval); sp.E = (yyvsp[-6].fval); sp.nu = (yyvsp[-5].fval); sp.rho = (yyvsp[-4].fval);
          sp.c = (yyvsp[-3].fval); sp.k = (yyvsp[-2].fval); sp.eh = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-8].ival)-1, sp ); 
        }
#line 9002 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 674:
#line 2412 "p.y" /* yacc.c:1646  */
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[-12].fval);  sp.E = (yyvsp[-11].fval);  sp.nu  = (yyvsp[-10].fval);  sp.rho = (yyvsp[-9].fval);
          sp.c = (yyvsp[-8].fval);  sp.k = (yyvsp[-7].fval);  sp.eh  = (yyvsp[-6].fval);  sp.P   = (yyvsp[-5].fval);  sp.Ta  = (yyvsp[-4].fval);
          sp.Q = (yyvsp[-3].fval); sp.W = (yyvsp[-2].fval); sp.Ixx = (yyvsp[-1].fval);  
          geoSource->addMat( (yyvsp[-13].ival)-1, sp );
        }
#line 9013 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9025 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9041 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9061 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9082 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9103 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9115 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9127 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9149 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 683:
#line 2524 "p.y" /* yacc.c:1646  */
    { StructProp sp; 
          sp.A = (yyvsp[-9].fval); sp.rho = (yyvsp[-8].fval); sp.Q = (yyvsp[-7].fval); sp.c = (yyvsp[-6].fval); 
          sp.sigma = (yyvsp[-5].fval); sp.k = (yyvsp[-4].fval); sp.eh = (yyvsp[-3].fval); sp.P = (yyvsp[-2].fval); sp.Ta = (yyvsp[-1].fval);
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[-11].ival)-1, sp );
        }
#line 9160 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 684:
#line 2531 "p.y" /* yacc.c:1646  */
    { // rigid element or joint with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-2].ival)-1, sp );
        }
#line 9171 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9186 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9206 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 687:
#line 2565 "p.y" /* yacc.c:1646  */
    { // rigid 8-node brick element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-4].ival)-1, sp );
        }
#line 9217 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9232 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9244 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9279 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9296 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9317 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9333 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9350 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9368 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9388 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9409 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9431 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9448 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9466 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9485 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9506 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9528 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9551 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9563 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9576 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9592 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9609 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9622 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9637 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9654 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9673 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9687 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9704 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9722 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9743 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 718:
#line 2978 "p.y" /* yacc.c:1646  */
    { // TorsionalSpring or TranslationalSpring (types 200,201,202)
          StructProp sp;
          sp.k1 = (yyvsp[-1].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9754 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9768 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9780 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 721:
#line 3003 "p.y" /* yacc.c:1646  */
    { // Acoustic rubber
          StructProp sp;
          sp.E0 = sp.E = (yyvsp[-1].fval);
          geoSource->addMat( (yyvsp[-3].ival)-1, sp );
        }
#line 9790 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9802 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 725:
#line 3032 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-1].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9812 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 726:
#line 3038 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-2].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-2].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9822 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 727:
#line 3044 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-3].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[-3].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[-1].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
#line 9833 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9845 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9858 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 730:
#line 3070 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9864 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 731:
#line 3072 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9870 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 732:
#line 3074 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9878 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 733:
#line 3078 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9884 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 734:
#line 3080 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9890 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 735:
#line 3082 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9896 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 736:
#line 3084 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
#line 9902 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 737:
#line 3086 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-1].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9910 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 738:
#line 3090 "p.y" /* yacc.c:1646  */
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9918 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 739:
#line 3096 "p.y" /* yacc.c:1646  */
    { domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().isCoupled = true; }
#line 9924 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 740:
#line 3098 "p.y" /* yacc.c:1646  */
    { domain->addWetInterface((yyvsp[-1].ival), (yyvsp[-1].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; }
#line 9932 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 741:
#line 3106 "p.y" /* yacc.c:1646  */
    { }
#line 9938 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 742:
#line 3108 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9948 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 743:
#line 3114 "p.y" /* yacc.c:1646  */
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-1].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
#line 9959 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 744:
#line 3121 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 9970 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 745:
#line 3128 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-2].fval), (yyvsp[-1].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-3].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 9981 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 9993 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 747:
#line 3143 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
#line 10004 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10016 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10028 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10040 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 752:
#line 3185 "p.y" /* yacc.c:1646  */
    { }
#line 10059 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 753:
#line 3187 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[-2].ival) == (yyvsp[-1].ival)) domain->solInfo().isMatching = true;
        }
#line 10069 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 754:
#line 3193 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[-4].ival) == (yyvsp[-3].ival)) domain->solInfo().isMatching = true;
        }
#line 10079 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 755:
#line 3201 "p.y" /* yacc.c:1646  */
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
#line 10086 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 757:
#line 3207 "p.y" /* yacc.c:1646  */
    { domain->addWetElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), 1.0, (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
#line 10094 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 758:
#line 3213 "p.y" /* yacc.c:1646  */
    { }
#line 10100 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 759:
#line 3215 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-2].ival), (yyvsp[-1].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[-2].ival) == (yyvsp[-1].ival)) domain->solInfo().isMatching = true;
        }
#line 10110 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 760:
#line 3221 "p.y" /* yacc.c:1646  */
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].fval), (yyvsp[-1].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[-4].ival) == (yyvsp[-3].ival)) domain->solInfo().isMatching = true;
        }
#line 10120 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 761:
#line 3231 "p.y" /* yacc.c:1646  */
    { }
#line 10126 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 762:
#line 3233 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contactsurface_mode = (yyvsp[-1].ival); }
#line 10132 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 763:
#line 3235 "p.y" /* yacc.c:1646  */
    { domain->AddMortarCond((yyvsp[-1].MortarCondObj)); }
#line 10138 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 764:
#line 3237 "p.y" /* yacc.c:1646  */
    { (yyvsp[-2].MortarCondObj)->SetConstraintOptions((yyvsp[-1].copt)); domain->AddMortarCond((yyvsp[-2].MortarCondObj)); }
#line 10144 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 765:
#line 3239 "p.y" /* yacc.c:1646  */
    { (yyvsp[-4].MortarCondObj)->SetConstraintOptions((yyvsp[-3].copt)); (yyvsp[-4].MortarCondObj)->SetCtcMode((yyvsp[-1].ival)); domain->AddMortarCond((yyvsp[-4].MortarCondObj)); }
#line 10150 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 766:
#line 3241 "p.y" /* yacc.c:1646  */
    { (yyvsp[-3].MortarCondObj)->SetCtcMode((yyvsp[-1].ival)); domain->AddMortarCond((yyvsp[-3].MortarCondObj)); }
#line 10156 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 767:
#line 3247 "p.y" /* yacc.c:1646  */
    { domain->solInfo().trivial_detection = true;
          (yyval.MortarCondObj) = new MortarHandler(0, 0);
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD);
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode); }
#line 10166 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 768:
#line 3253 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-1].ival), (yyvsp[0].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD);
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10177 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 769:
#line 3260 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-2].ival), (yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[0].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10188 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 770:
#line 3267 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-1].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10199 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 771:
#line 3274 "p.y" /* yacc.c:1646  */
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-1].fval), (yyvsp[0].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[-2].ival));
          (yyval.MortarCondObj)->SetCtcMode(domain->solInfo().contactsurface_mode);
        }
#line 10210 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10222 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10235 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10248 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 10261 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 777:
#line 3319 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dist_acme = (yyvsp[-1].ival); }
#line 10267 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 778:
#line 3321 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_secondary = true; }
#line 10273 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 779:
#line 3323 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_ghosting = true; }
#line 10279 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 780:
#line 3325 "p.y" /* yacc.c:1646  */
    { domain->solInfo().shell_simple_lofting = true; }
#line 10285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 781:
#line 3327 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_multiple_interactions = true; }
#line 10291 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 782:
#line 3329 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sharp_non_sharp_angle = (yyvsp[-1].fval); }
#line 10297 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 783:
#line 3331 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normal_smoothing = false; }
#line 10303 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 784:
#line 3333 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normal_smoothing_distance = (yyvsp[-1].fval); }
#line 10309 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 785:
#line 3335 "p.y" /* yacc.c:1646  */
    { domain->solInfo().resolution_method = (yyvsp[-1].ival); }
#line 10315 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 786:
#line 3337 "p.y" /* yacc.c:1646  */
    { domain->solInfo().old_dynamic_search = true; }
#line 10321 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 787:
#line 3339 "p.y" /* yacc.c:1646  */
    { domain->solInfo().partition_gap = true; }
#line 10327 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 788:
#line 3341 "p.y" /* yacc.c:1646  */
    { domain->solInfo().default_penalty = true; }
#line 10333 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 789:
#line 3343 "p.y" /* yacc.c:1646  */
    { domain->solInfo().global_search_cull = true; }
#line 10339 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 790:
#line 3345 "p.y" /* yacc.c:1646  */
    { domain->solInfo().no_warped_volume = true; }
#line 10345 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 791:
#line 3347 "p.y" /* yacc.c:1646  */
    { domain->solInfo().auto_tol = true; }
#line 10351 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 792:
#line 3349 "p.y" /* yacc.c:1646  */
    { domain->solInfo().auto_tol = true;
          domain->solInfo().agressive_tolerances = true; }
#line 10358 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 793:
#line 3352 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skip_physical_faces = true; }
#line 10364 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 794:
#line 3354 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ffi_debug = bool((yyvsp[-1].ival)); }
#line 10370 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 795:
#line 3356 "p.y" /* yacc.c:1646  */
    { domain->solInfo().mortar_scaling = (yyvsp[-1].fval); }
#line 10376 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 796:
#line 3358 "p.y" /* yacc.c:1646  */
    { domain->solInfo().mortar_integration_rule = (yyvsp[-1].ival); }
#line 10382 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 797:
#line 3360 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_clr = (yyvsp[-1].fval); }
#line 10388 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 798:
#line 3362 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_cqr = (yyvsp[-1].fval); }
#line 10394 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 799:
#line 3364 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_betab = (yyvsp[-1].fval); }
#line 10400 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 800:
#line 3366 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_alpha = (yyvsp[-1].fval); }
#line 10406 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 801:
#line 3368 "p.y" /* yacc.c:1646  */
    { domain->solInfo().andes_betam = (yyvsp[-1].fval); }
#line 10412 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 802:
#line 3370 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nlmembrane_pressure_type = (yyvsp[-1].ival); }
#line 10418 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 803:
#line 3373 "p.y" /* yacc.c:1646  */
    { geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10424 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 804:
#line 3375 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scalePosCoords = true;
          domain->solInfo().xScaleFactor = (yyvsp[-4].fval);
          domain->solInfo().yScaleFactor = (yyvsp[-3].fval);
          domain->solInfo().zScaleFactor = (yyvsp[-2].fval);
          geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10434 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 805:
#line 3381 "p.y" /* yacc.c:1646  */
    { geoSource->addNode((yyvsp[0].nval).num, (yyvsp[0].nval).xyz, (yyvsp[0].nval).cp, (yyvsp[0].nval).cd); }
#line 10440 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 806:
#line 3385 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-4].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-3].fval); (yyval.nval).xyz[1] = (yyvsp[-2].fval);  (yyval.nval).xyz[2] = (yyvsp[-1].fval);  (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10446 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 807:
#line 3387 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-3].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-2].fval); (yyval.nval).xyz[1] = (yyvsp[-1].fval);  (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10452 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 808:
#line 3389 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-2].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-1].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; }
#line 10458 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 809:
#line 3391 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-6].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-5].fval); (yyval.nval).xyz[1] = (yyvsp[-4].fval);  (yyval.nval).xyz[2] = (yyvsp[-3].fval);  (yyval.nval).cp = (yyvsp[-2].ival); (yyval.nval).cd = (yyvsp[-1].ival);
          if((yyvsp[-2].ival) != 0) domain->solInfo().basicPosCoords = false;
          if((yyvsp[-1].ival) != 0) domain->solInfo().basicDofCoords = false; }
#line 10466 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 810:
#line 3395 "p.y" /* yacc.c:1646  */
    { (yyval.nval).num = (yyvsp[-5].ival)-1; (yyval.nval).xyz[0] = (yyvsp[-4].fval); (yyval.nval).xyz[1] = (yyvsp[-3].fval);  (yyval.nval).xyz[2] = (yyvsp[-2].fval);  (yyval.nval).cp = (yyvsp[-1].ival); (yyval.nval).cd = (yyvsp[-1].ival);
          if((yyvsp[-1].ival) != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } }
#line 10473 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 811:
#line 3400 "p.y" /* yacc.c:1646  */
    { /* Define each Element */
          geoSource->addElem((yyvsp[-3].ival)-1, (yyvsp[-2].ival), (yyvsp[-1].nl).num, (yyvsp[-1].nl).nd);}
#line 10480 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 812:
#line 3405 "p.y" /* yacc.c:1646  */
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[0].ival)-1;}
#line 10486 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 813:
#line 3407 "p.y" /* yacc.c:1646  */
    { if((yyval.nl).num == 500) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[0].ival)-1; (yyval.nl).num++;}
#line 10493 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 814:
#line 3412 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-3].ival)-1; (yyval.bcval).dofnum = (yyvsp[-2].ival)-1; (yyval.bcval).val = (yyvsp[-1].fval); (yyval.bcval).mtype = BCond::Axial; }
#line 10499 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 815:
#line 3414 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = (yyvsp[-1].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = BCond::Axial; }
#line 10505 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 816:
#line 3416 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-4].ival)-1; (yyval.bcval).dofnum = (yyvsp[-3].ival)-1; (yyval.bcval).val = (yyvsp[-2].fval); (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[-1].ival); }
#line 10511 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 817:
#line 3418 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-3].ival)-1; (yyval.bcval).dofnum = (yyvsp[-2].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[-1].ival); }
#line 10517 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 818:
#line 3422 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[-1].fval); }
#line 10523 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 819:
#line 3426 "p.y" /* yacc.c:1646  */
    { (yyval.bcval).nnum = (yyvsp[-2].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[-1].fval); }
#line 10529 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 820:
#line 3430 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum = (yyvsp[-4].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[-3].ival)-1; (yyval.cxbcval).reval = (yyvsp[-2].fval); (yyval.cxbcval).imval = (yyvsp[-1].fval);  }
#line 10535 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 821:
#line 3432 "p.y" /* yacc.c:1646  */
    { (yyval.cxbcval).nnum = (yyvsp[-3].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[-2].ival)-1; (yyval.cxbcval).reval = (yyvsp[-1].fval); (yyval.cxbcval).imval = 0.0; }
#line 10541 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 823:
#line 3437 "p.y" /* yacc.c:1646  */
    { geoSource->setCSFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d); }
#line 10547 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 825:
#line 3442 "p.y" /* yacc.c:1646  */
    { geoSource->setFrame((yyvsp[0].frame).num,(yyvsp[0].frame).d);  }
#line 10553 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 826:
#line 3446 "p.y" /* yacc.c:1646  */
    { (yyval.frame).num = (yyvsp[-10].ival)-1; 
          (yyval.frame).d[0] = (yyvsp[-9].fval); (yyval.frame).d[1] = (yyvsp[-8].fval); (yyval.frame).d[2] = (yyvsp[-7].fval);
          (yyval.frame).d[3] = (yyvsp[-6].fval); (yyval.frame).d[4] = (yyvsp[-5].fval); (yyval.frame).d[5] = (yyvsp[-4].fval);
          (yyval.frame).d[6] = (yyvsp[-3].fval); (yyval.frame).d[7] = (yyvsp[-2].fval); (yyval.frame).d[8] = (yyvsp[-1].fval); }
#line 10562 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 827:
#line 3451 "p.y" /* yacc.c:1646  */
    { (yyval.frame).num = (yyvsp[-3].ival)-1;
          geoSource->makeEframe((yyvsp[-3].ival)-1, (yyvsp[-1].ival), (yyval.frame).d); }
#line 10569 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 829:
#line 3457 "p.y" /* yacc.c:1646  */
    { geoSource->setNodalFrame((yyvsp[0].nframe).id,(yyvsp[0].nframe).o,(yyvsp[0].nframe).d,(yyvsp[0].nframe).type); }
#line 10575 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 830:
#line 3461 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-10].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;  (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[-9].fval); (yyval.nframe).d[1] = (yyvsp[-8].fval); (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval); (yyval.nframe).d[4] = (yyvsp[-5].fval); (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10586 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 831:
#line 3468 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-13].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = (yyvsp[-12].fval);  (yyval.nframe).o[1] = (yyvsp[-11].fval);  (yyval.nframe).o[2] = (yyvsp[-10].fval);
          (yyval.nframe).d[0] = (yyvsp[-9].fval);  (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval);  (yyval.nframe).d[4] = (yyvsp[-5].fval);  (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10597 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 832:
#line 3475 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-11].ival);
          (yyval.nframe).type = (yyvsp[-10].ival);
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;   (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[-9].fval); (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval); (yyval.nframe).d[4] = (yyvsp[-5].fval);  (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10608 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 833:
#line 3482 "p.y" /* yacc.c:1646  */
    { (yyval.nframe).id = (yyvsp[-14].ival);
          (yyval.nframe).type = (yyvsp[-13].ival);
          (yyval.nframe).o[0] = (yyvsp[-12].fval);  (yyval.nframe).o[1] = (yyvsp[-11].fval);  (yyval.nframe).o[2] = (yyvsp[-10].fval);
          (yyval.nframe).d[0] = (yyvsp[-9].fval);  (yyval.nframe).d[1] = (yyvsp[-8].fval);  (yyval.nframe).d[2] = (yyvsp[-7].fval);
          (yyval.nframe).d[3] = (yyvsp[-6].fval);  (yyval.nframe).d[4] = (yyvsp[-5].fval); (yyval.nframe).d[5] = (yyvsp[-4].fval);
          (yyval.nframe).d[6] = (yyvsp[-3].fval); (yyval.nframe).d[7] = (yyvsp[-2].fval); (yyval.nframe).d[8] = (yyvsp[-1].fval); }
#line 10619 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 835:
#line 3492 "p.y" /* yacc.c:1646  */
    { OffsetData od;
	  od.first = (yyvsp[-5].ival)-1; od.last = (yyvsp[-4].ival)-1;
	  od.o[0] = (yyvsp[-3].fval); od.o[1] = (yyvsp[-2].fval); od.o[2] = (yyvsp[-1].fval); 
	  geoSource->addOffset(od); }
#line 10628 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 836:
#line 3499 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 10634 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 837:
#line 3501 "p.y" /* yacc.c:1646  */
    { geoSource->setLocalIndex((yyvsp[-1].ival)-1); }
#line 10640 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 838:
#line 3504 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-3].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10647 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 839:
#line 3507 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-4].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10655 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 840:
#line 3512 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-2].ival)-1,(yyvsp[-1].ival)-1); }
#line 10661 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 841:
#line 3514 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1); 
          geoSource->setElementLumpingWeight((yyvsp[-4].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10669 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 842:
#line 3518 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-5].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; 
          domain->solInfo().reduceFollower = true; }
#line 10678 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 843:
#line 3523 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-5].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10687 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 844:
#line 3529 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1,(yyvsp[-1].ival)-1); }
#line 10693 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 845:
#line 3531 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-6].ival)-1,(yyvsp[-1].fval)); 
          domain->solInfo().elemLumpPodRom = true; }
#line 10701 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 846:
#line 3535 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[-7].ival)-1,(yyvsp[-2].fval));   
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10710 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 847:
#line 3541 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-5].ival)-1,(yyvsp[-4].ival)-1,(yyvsp[-3].ival)-1,-2,(yyvsp[-1].fval)); }
#line 10716 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 848:
#line 3543 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,(yyvsp[-5].ival)-1,-2,(yyvsp[-3].fval));
          geoSource->setElementLumpingWeight((yyvsp[-7].ival)-1,(yyvsp[-1].fval));
          domain->solInfo().elemLumpPodRom = true; }
#line 10724 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 849:
#line 3547 "p.y" /* yacc.c:1646  */
    { geoSource->setAttrib((yyvsp[-8].ival)-1,(yyvsp[-7].ival)-1,(yyvsp[-6].ival)-1,-2,(yyvsp[-4].fval)); 
          geoSource->setElementLumpingWeight((yyvsp[-8].ival)-1,(yyvsp[-2].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; }
#line 10733 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 850:
#line 3553 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        }
#line 10742 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 851:
#line 3558 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-3].ival); i<(yyvsp[-2].ival)+1; ++i)
            geoSource->setAttrib(i-1,(yyvsp[-1].ival)-1);
        }
#line 10751 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 852:
#line 3563 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-5].ival); i<(yyvsp[-4].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1);
        }
#line 10760 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 853:
#line 3568 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-6].ival); i<(yyvsp[-5].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, -2, (yyvsp[-1].fval));
        }
#line 10769 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 854:
#line 3575 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex(0); }
#line 10776 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 855:
#line 3578 "p.y" /* yacc.c:1646  */
    { domain->solInfo().elemLumpPodRom = true;
          geoSource->setLocalIndex((yyvsp[-1].ival)-1); }
#line 10783 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 856:
#line 3581 "p.y" /* yacc.c:1646  */
    { geoSource->setElementLumpingWeight((yyvsp[-2].ival) - 1, (yyvsp[-1].fval)); }
#line 10789 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 857:
#line 3583 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = true;}
#line 10795 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 858:
#line 3585 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = true;}
#line 10801 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 862:
#line 3594 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInLocalBasesAuxi[std::make_pair((yyvsp[-3].ival)-1,(yyvsp[-2].ival)-1)] = std::string((yyvsp[-1].strval)); }
#line 10807 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 863:
#line 3598 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInLocalBasesCent.push_back(std::string((yyvsp[-1].strval))); }
#line 10813 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 864:
#line 3602 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ReducedStiffness = true;}
#line 10819 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 865:
#line 3604 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackStiffVec((yyvsp[-1].fval));}
#line 10825 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 867:
#line 3609 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodSize = (yyvsp[-2].ival);
          domain->solInfo().maxDeimBasisSize = (yyvsp[-1].ival);}
#line 10832 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 868:
#line 3612 "p.y" /* yacc.c:1646  */
    { geoSource->pushBackUDEIMVec((yyvsp[-1].fval));}
#line 10838 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 870:
#line 3617 "p.y" /* yacc.c:1646  */
    { domain->solInfo().DEIMPodRom = true;
          geoSource->setSampleNodesAndSlots((yyvsp[-2].ival)-1,(yyvsp[-1].ival));}
#line 10845 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 871:
#line 3620 "p.y" /* yacc.c:1646  */
    { geoSource->setSampleElemsAndDOFs((yyvsp[-2].ival)-1,(yyvsp[-1].ival));
          domain->solInfo().UDEIMPodRom = true;}
#line 10852 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 872:
#line 3626 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = 0; }
#line 10858 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 873:
#line 3628 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[-1].ival); }
#line 10864 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 874:
#line 3630 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-2].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          geoSource->setElementPressure(pbc); }
#line 10872 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 875:
#line 3634 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-3].ival); i < ((yyvsp[-2].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            geoSource->setElementPressure(pbc);
          } }
#line 10882 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 876:
#line 3640 "p.y" /* yacc.c:1646  */
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[-2].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
#line 10891 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 877:
#line 3645 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-4].ival); i <= (yyvsp[-2].ival); ++i) {
            PressureBCond *pbc = new PressureBCond[1];
            pbc[0].setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            geoSource->addSurfacePressure(1, pbc);
            if(geoSource->getNumSurfacePressure() > 1) delete [] pbc;
          } }
#line 10902 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 878:
#line 3652 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-3].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          geoSource->setElementPressure(pbc); }
#line 10910 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 879:
#line 3656 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-4].ival); i < ((yyvsp[-3].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
            geoSource->setElementPressure(pbc);
          } }
#line 10920 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 880:
#line 3662 "p.y" /* yacc.c:1646  */
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[-3].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; }
#line 10929 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 881:
#line 3668 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-4].ival)-1, (yyvsp[-1].fval), (yyval.ival), true);
          pbc.face = (yyvsp[-2].ival)-1;
          geoSource->setElementPressure(pbc); }
#line 10938 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 882:
#line 3673 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-5].ival); i < ((yyvsp[-4].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-1].fval), (yyval.ival), true);
            pbc.face = (yyvsp[-2].ival)-1;
            geoSource->setElementPressure(pbc);
          } }
#line 10949 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 883:
#line 3680 "p.y" /* yacc.c:1646  */
    { PressureBCond pbc;
          pbc.setData((yyvsp[-5].ival)-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
          pbc.face = (yyvsp[-3].ival)-1;
          geoSource->setElementPressure(pbc); }
#line 10958 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 884:
#line 3685 "p.y" /* yacc.c:1646  */
    { for(int i = (yyvsp[-6].ival); i < ((yyvsp[-5].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[-2].fval), (yyval.ival), (yyvsp[-1].ival));
            pbc.face = (yyvsp[-3].ival)-1;
            geoSource->setElementPressure(pbc);
          } }
#line 10969 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 885:
#line 3694 "p.y" /* yacc.c:1646  */
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        }
#line 10978 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 886:
#line 3699 "p.y" /* yacc.c:1646  */
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[-1].ival));
          geoSource->setConsistentPFlag(false);
        }
#line 10987 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 887:
#line 3713 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassAugmentation = (yyvsp[-1].ival); }
#line 10993 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 888:
#line 3717 "p.y" /* yacc.c:1646  */
    { }
#line 10999 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 889:
#line 3719 "p.y" /* yacc.c:1646  */
    { geoSource->setElementPreLoad( (yyvsp[-2].ival)-1, (yyvsp[-1].fval) ); }
#line 11005 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 890:
#line 3721 "p.y" /* yacc.c:1646  */
    { int i;
          for(i=(yyvsp[-4].ival); i<((yyvsp[-2].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[-1].fval) );
        }
#line 11014 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 891:
#line 3726 "p.y" /* yacc.c:1646  */
    { double load[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
          geoSource->setElementPreLoad( (yyvsp[-4].ival)-1, load ); }
#line 11021 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 892:
#line 3729 "p.y" /* yacc.c:1646  */
    { double load[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
          int i;
          for(i=(yyvsp[-6].ival); i<((yyvsp[-4].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        }
#line 11031 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 893:
#line 3737 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivity = true; }
#line 11037 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 894:
#line 3739 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivityMethod = (SolverInfo::SensitivityMethod) (yyvsp[-1].ival); }
#line 11043 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 895:
#line 3741 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readShapeSen = true;
          domain->solInfo().readInShapeSen = (yyvsp[-1].strval); }
#line 11050 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 896:
#line 3744 "p.y" /* yacc.c:1646  */
    { domain->solInfo().sensitivityTol = (yyvsp[-1].fval); }
#line 11056 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 897:
#line 3746 "p.y" /* yacc.c:1646  */
    { domain->solInfo().qsMaxvelSen = (yyvsp[-1].fval); }
#line 11062 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 898:
#line 3748 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ratioSensitivityTol = (yyvsp[-1].fval); }
#line 11068 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 899:
#line 3750 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ksParameter = (yyvsp[-1].fval); }
#line 11074 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 900:
#line 3752 "p.y" /* yacc.c:1646  */
    { domain->solInfo().ksMax = (yyvsp[-1].fval); }
#line 11080 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 901:
#line 3754 "p.y" /* yacc.c:1646  */
    { }
#line 11086 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 902:
#line 3756 "p.y" /* yacc.c:1646  */
    { }
#line 11092 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 903:
#line 3758 "p.y" /* yacc.c:1646  */
    { }
#line 11098 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 904:
#line 3762 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11104 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 905:
#line 3764 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11110 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 906:
#line 3766 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11116 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 907:
#line 3768 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11122 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 908:
#line 3770 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11128 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 909:
#line 3772 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11134 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 910:
#line 3774 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11140 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 911:
#line 3776 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11146 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 912:
#line 3778 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11152 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 913:
#line 3780 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11158 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 914:
#line 3782 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11164 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 915:
#line 3784 "p.y" /* yacc.c:1646  */
    { domain->setDispNode((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].ival)-1, (yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 11170 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 916:
#line 3788 "p.y" /* yacc.c:1646  */
    { domain->setStressNodes((yyvsp[0].ival)); }
#line 11176 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 917:
#line 3790 "p.y" /* yacc.c:1646  */
    { domain->setStressNodes((yyvsp[0].ival)); }
#line 11182 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 918:
#line 3794 "p.y" /* yacc.c:1646  */
    { domain->setThicknessGroup((yyvsp[0].ival)); }
#line 11188 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 919:
#line 3796 "p.y" /* yacc.c:1646  */
    { domain->setThicknessGroup((yyvsp[0].ival)); }
#line 11194 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 920:
#line 3800 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setProbType(SolverInfo::Static); }
#line 11200 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 921:
#line 3802 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solvercntl = (yyvsp[0].scntl); }
#line 11206 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 923:
#line 3805 "p.y" /* yacc.c:1646  */
    { // activate piecewise constant configuration dependent external forces for a linear dynamic analysis
          domain->solInfo().piecewise = true;
        }
#line 11214 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 924:
#line 3809 "p.y" /* yacc.c:1646  */
    { // activate piecewise constant configuration dependent external forces for a linear static analysis
          domain->solInfo().piecewise = true;
          domain->solInfo().piecewise_dlambda = (yyvsp[-2].fval);
          domain->solInfo().piecewise_maxLambda = (yyvsp[-1].fval);
        }
#line 11224 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 925:
#line 3815 "p.y" /* yacc.c:1646  */
    { domain->solInfo().coupled_scale = (yyvsp[-1].fval); }
#line 11230 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 926:
#line 3817 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-1].ival);
          domain->curvatureFlag = 0; }
#line 11237 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 927:
#line 3820 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-2].ival);
          domain->curvatureConst1 = (yyvsp[-1].fval);
          domain->curvatureFlag = 1; }
#line 11245 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 928:
#line 3824 "p.y" /* yacc.c:1646  */
    { domain->sommerfeldType = (yyvsp[-3].ival);
          domain->curvatureConst1 = (yyvsp[-2].fval);
          domain->curvatureConst2 = (yyvsp[-1].fval);
          domain->curvatureFlag = 2; }
#line 11254 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 929:
#line 3829 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dmpc = bool((yyvsp[-1].ival)); }
#line 11260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 930:
#line 3831 "p.y" /* yacc.c:1646  */
    { domain->solInfo().dbccheck = bool((yyvsp[-1].ival)); }
#line 11266 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 931:
#line 3835 "p.y" /* yacc.c:1646  */
    { domain->solInfo().loadcases.push_back((yyvsp[0].ival)); }
#line 11272 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 932:
#line 3837 "p.y" /* yacc.c:1646  */
    { domain->solInfo().loadcases.push_back((yyvsp[0].ival)); }
#line 11278 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 933:
#line 3841 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = 0; }
#line 11286 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 934:
#line 3845 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = (yyvsp[-1].ival); }
#line 11294 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 935:
#line 3849 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
          (yyval.scntl)->type = SolverSelection::Direct;
          (yyval.scntl)->subtype = (yyvsp[-1].ival); }
#line 11302 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 936:
#line 3853 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->pivot = true; }
#line 11311 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 937:
#line 3858 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->scaled = true; }
#line 11320 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 938:
#line 3863 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Direct;
      (yyval.scntl)->subtype = (yyvsp[-2].ival);
      (yyval.scntl)->unsymmetric = true; }
#line 11329 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 939:
#line 3868 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-1].ival); }
#line 11337 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 940:
#line 3872 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-2].ival);
      (yyval.scntl)->precond = (yyvsp[-1].ival); }
#line 11346 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 941:
#line 3877 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-3].ival);
      (yyval.scntl)->precond = (yyvsp[-2].ival);
      (yyval.scntl)->tol=(yyvsp[-1].fval); }
#line 11356 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 942:
#line 3883 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Iterative;
      (yyval.scntl)->iterType = (yyvsp[-4].ival);
      (yyval.scntl)->precond = (yyvsp[-3].ival);
      (yyval.scntl)->tol = (yyvsp[-2].fval);
      (yyval.scntl)->maxit = (yyvsp[-1].ival); }
#line 11367 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11379 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11392 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 945:
#line 3907 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti; }
#line 11399 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 946:
#line 3910 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.version = (FetiInfo::Version) ((yyvsp[-1].ival)-1); }
#line 11407 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 947:
#line 3914 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-2].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval);
      (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-2].ival); }
#line 11417 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 948:
#line 3920 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.maxit = (yyvsp[-3].ival);
      (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval);
      (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-1].ival); }
#line 11427 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 949:
#line 3926 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.version = (FetiInfo::Version) ((yyvsp[-2].ival)-1);
      (yyval.scntl)->fetiInfo.feti2version = (FetiInfo::Feti2Version) (yyvsp[-1].ival); }
#line 11436 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 950:
#line 3931 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp; }
#line 11445 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 951:
#line 3936 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::FetiLib;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp; }
#line 11454 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 952:
#line 3941 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::Feti;
      (yyval.scntl)->fetiInfo.corners = FetiInfo::allCorners6;
      (yyval.scntl)->fetiInfo.version = FetiInfo::fetidp;
      (yyval.scntl)->fetiInfo.dph_flag = true; }
#line 11464 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11478 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11495 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11513 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11533 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11553 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 958:
#line 4016 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = new SolverCntl(default_cntl);
      (yyval.scntl)->type = SolverSelection::BlockDiag; }
#line 11560 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 959:
#line 4019 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = &domain->solInfo().solvercntls[(yyvsp[-1].ival)]; }
#line 11566 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 960:
#line 4023 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = &domain->solInfo().solvercntls[(yyvsp[-2].ival)];
          *(yyval.scntl) = *(yyvsp[0].scntl); }
#line 11573 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 961:
#line 4027 "p.y" /* yacc.c:1646  */
    { (yyval.scntl) = (yyvsp[0].scntl); }
#line 11579 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 962:
#line 4029 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->verbose = 1; }
#line 11585 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 963:
#line 4031 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->verbose = (yyvsp[-1].ival); }
#line 11591 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 964:
#line 4033 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.printNumber = (yyvsp[-1].ival); }
#line 11597 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 965:
#line 4035 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->trbm = (yyvsp[-1].fval); }
#line 11603 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 966:
#line 4037 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_renum = (yyvsp[-1].ival); }
#line 11609 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 967:
#line 4039 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_maxsup = (yyvsp[-1].ival); }
#line 11615 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 968:
#line 4041 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->sparse_defblk = (yyvsp[-1].ival); }
#line 11621 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 969:
#line 4043 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_tau = (yyvsp[-1].fval); }
#line 11627 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 970:
#line 4045 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_maxsize = (yyvsp[-1].ival); }
#line 11633 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 971:
#line 4047 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) < 0) {
            (yyvsp[-1].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          (yyval.scntl)->spooles_maxdomainsize = (yyvsp[-1].ival); }
#line 11644 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 972:
#line 4054 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_seed = (yyvsp[-1].ival); }
#line 11650 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 973:
#line 4056 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].fval) < 0.0) || ((yyvsp[-1].fval) > 1.0)) {
            (yyvsp[-1].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          (yyval.scntl)->spooles_maxzeros = (yyvsp[-1].fval); }
#line 11661 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 974:
#line 4063 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_msglvl = (yyvsp[-1].ival); }
#line 11667 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 975:
#line 4065 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_scale = (yyvsp[-1].ival); }
#line 11673 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 976:
#line 4067 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->spooles_renum = (yyvsp[-1].ival); }
#line 11679 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 977:
#line 4069 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_icntl[(yyvsp[-2].ival)] = (yyvsp[-1].ival); }
#line 11685 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 978:
#line 4071 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_cntl[(yyvsp[-2].ival)] = (yyvsp[-1].fval); }
#line 11691 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 979:
#line 4073 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_mineq = (yyvsp[-1].ival); }
#line 11697 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 980:
#line 4075 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->mumps_stride = (yyvsp[-1].ival); }
#line 11703 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 981:
#line 4077 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->goldfarb_tol = (yyvsp[-1].fval); }
#line 11709 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 982:
#line 4079 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->tol = (yyvsp[-1].fval); }
#line 11715 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 983:
#line 4081 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->iterSubtype = (yyvsp[-1].ival); }
#line 11721 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 984:
#line 4083 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->ilu_droptol = (yyvsp[-1].fval); }
#line 11727 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 985:
#line 4085 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->maxit = (yyvsp[-1].ival); 
          (yyval.scntl)->fetiInfo.maxit = (yyvsp[-1].ival); }
#line 11734 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 986:
#line 4088 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->maxvecsize = (yyvsp[-1].ival);
          (yyval.scntl)->fetiInfo.maxortho = (yyvsp[-1].ival); }
#line 11741 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 987:
#line 4091 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11747 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 988:
#line 4093 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11753 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 989:
#line 4095 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.precno = FetiInfo::lumped; }
#line 11759 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 990:
#line 4097 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->precond = (yyvsp[-1].ival);
          if(isFeti((yyval.scntl)->type) && (((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 3))) {
            (yyvsp[-1].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          (yyval.scntl)->fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[-1].ival); }
#line 11770 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 991:
#line 4104 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 1)) {
            (yyvsp[-1].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          (yyval.scntl)->fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[-1].ival); }
#line 11780 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 992:
#line 4110 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 1)) {
            (yyvsp[-1].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          (yyval.scntl)->fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[-1].ival); }
#line 11790 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 993:
#line 4116 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tol = (yyvsp[-1].fval); }
#line 11796 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 994:
#line 4118 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tol = (yyvsp[-2].fval); 
          (yyval.scntl)->fetiInfo.absolute_tol = (yyvsp[-1].fval); }
#line 11803 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 995:
#line 4121 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.stagnation_tol = (yyvsp[-1].fval); }
#line 11809 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 996:
#line 4123 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.stagnation_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.absolute_stagnation_tol = (yyvsp[-1].fval); }
#line 11816 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 997:
#line 4126 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_proj_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.dual_proj_tol = (yyvsp[-1].fval); }
#line 11823 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 998:
#line 4129 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_plan_maxit = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.dual_plan_maxit = (yyvsp[-1].ival); }
#line 11830 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 999:
#line 4132 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.primal_plan_tol = (yyvsp[-2].fval);
          (yyval.scntl)->fetiInfo.dual_plan_tol = (yyvsp[-1].fval); }
#line 11837 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1000:
#line 4135 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.noCoarse = 1; }
#line 11843 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11866 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1002:
#line 4156 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 1; 
          (yyval.scntl)->fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11873 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1003:
#line 4159 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11879 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1004:
#line 4161 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 2;
          (yyval.scntl)->fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11886 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1005:
#line 4164 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11892 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1006:
#line 4166 "p.y" /* yacc.c:1646  */
    { if(((yyvsp[-1].ival) < 0) || ((yyvsp[-1].ival) > 2)) (yyvsp[-1].ival) = 2;
          (yyval.scntl)->fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11899 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1007:
#line 4169 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[-1].ival); }
#line 11905 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1008:
#line 4171 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_element = true; }
#line 11911 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1009:
#line 4173 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_element = true; }
#line 11917 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1010:
#line 4175 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.fsi_corner = (yyvsp[-1].ival); }
#line 11923 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1011:
#line 4177 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.splitLocalFsi = false; }
#line 11929 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1012:
#line 4179 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.wetcorners = true; }
#line 11935 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1013:
#line 4181 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[-1].ival); }
#line 11941 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1014:
#line 4183 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[-2].ival); 
          (yyval.scntl)->fetiInfo.pick_unsafe_corners = bool((yyvsp[-1].ival)); }
#line 11948 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11960 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 11983 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12003 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1018:
#line 4229 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numdir = (yyvsp[-1].ival); 
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12011 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1019:
#line 4233 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.numdir = (yyvsp[-1].ival); 
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12020 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1020:
#line 4238 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numdir = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[-1].ival);
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12029 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1021:
#line 4243 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[-1].ival);
          (yyval.scntl)->fetiInfo.numdir = (yyvsp[-2].ival);
          if((yyval.scntl)->fetiInfo.augment == FetiInfo::none)
            (yyval.scntl)->fetiInfo.augment = FetiInfo::Edges; }
#line 12039 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1022:
#line 4249 "p.y" /* yacc.c:1646  */
    { if ((yyvsp[-1].ival) == 2)
           (yyval.scntl)->fetiInfo.augmentimpl = FetiInfo::Primal;
        }
#line 12047 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1023:
#line 4253 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.orthotol = (yyvsp[-1].fval); }
#line 12053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1024:
#line 4255 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.orthotol = (yyvsp[-2].fval); 
          (yyval.scntl)->fetiInfo.orthotol2 = (yyvsp[-1].fval); }
#line 12060 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1025:
#line 4258 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.grbm_tol = (yyvsp[-1].fval); }
#line 12066 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1026:
#line 4260 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.crbm_tol = (yyvsp[-1].fval); }
#line 12072 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1027:
#line 4262 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cct_tol = (yyvsp[-1].fval); }
#line 12078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1028:
#line 4264 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.rebuildcct = int((yyvsp[-1].ival)); }
#line 12084 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1029:
#line 4266 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.uproj = (yyvsp[-1].ival); }
#line 12090 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1030:
#line 4268 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.printMatLab = 1; }
#line 12096 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1031:
#line 4270 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->printMatLab = 1;
          (yyval.scntl)->printMatLabFile = (yyvsp[-1].strval); }
#line 12103 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1032:
#line 4273 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.local_cntl = (yyval.scntl)->fetiInfo.kii_cntl = (yyvsp[0].scntl); }
#line 12109 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1033:
#line 4275 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.coarse_cntl = (yyvsp[0].scntl); }
#line 12115 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1034:
#line 4277 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.auxcoarse_cntl = (yyvsp[0].scntl); }
#line 12121 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1035:
#line 4279 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cct_cntl = (yyvsp[0].scntl); }
#line 12127 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12155 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1037:
#line 4305 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gmresResidual = true; }
#line 12161 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1038:
#line 4307 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gmresResidual = bool((yyvsp[-1].ival)); }
#line 12167 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1039:
#line 4309 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.pickAnyCorner = (yyvsp[-1].ival); }
#line 12173 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1040:
#line 4311 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.type = FetiInfo::nonlinear;
          (yyval.scntl)->fetiInfo.nlPrecFlg = 1; }
#line 12180 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1041:
#line 4314 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.type = FetiInfo::nonlinear;
	  (yyval.scntl)->fetiInfo.nlPrecFlg = (yyvsp[-1].ival); }
#line 12187 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1042:
#line 4317 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.numcgm = (yyvsp[-1].ival); }
#line 12193 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1043:
#line 4319 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.tolcgm = (yyvsp[-1].fval); }
#line 12199 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1044:
#line 4321 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.spaceDimension = (yyvsp[-1].ival); }
#line 12205 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1045:
#line 4323 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.krylovtype = (yyvsp[-1].ival); }
#line 12211 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1046:
#line 4325 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.krylovtype = (yyvsp[-1].ival); }
#line 12217 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1047:
#line 4327 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.lumpedinterface = 1; }
#line 12223 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1048:
#line 4329 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.saveMemCoarse = 1; }
#line 12229 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1049:
#line 4331 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[-1].ival); }
#line 12235 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1050:
#line 4333 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.complex_hermitian = true; }
#line 12242 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1051:
#line 4336 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpcflag = (yyvsp[-1].ival); }
#line 12248 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1052:
#line 4338 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpcflag = (yyvsp[-1].ival); }
#line 12254 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1053:
#line 4340 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-1].ival); }
#line 12260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1054:
#line 4342 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-1].ival); }
#line 12266 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1055:
#line 4344 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12273 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1056:
#line 4347 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12280 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1057:
#line 4350 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-2].ival); 
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-1].ival); }
#line 12287 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1058:
#line 4353 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-1].ival); }
#line 12294 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1059:
#line 4356 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-4].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-3].ival);
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12302 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1060:
#line 4360 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[-4].ival);
          (yyval.scntl)->fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[-3].ival); 
          (yyval.scntl)->fetiInfo.mpcBlkOverlap = (yyvsp[-1].ival); }
#line 12310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1061:
#line 4364 "p.y" /* yacc.c:1646  */
    { if((yyvsp[-1].ival) < 1) (yyval.scntl)->fetiInfo.useMRHS = false; }
#line 12316 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1062:
#line 4366 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.gamma = (yyvsp[-1].fval); }
#line 12322 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1063:
#line 4368 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.linesearch_maxit = (yyvsp[-2].ival);
          (yyval.scntl)->fetiInfo.linesearch_tau = (yyvsp[-1].fval); }
#line 12329 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1064:
#line 4371 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.bmpc = bool((yyvsp[-1].ival)); }
#line 12335 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1065:
#line 4373 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.cmpc = bool((yyvsp[-1].ival)); }
#line 12341 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1066:
#line 4375 "p.y" /* yacc.c:1646  */
    { (yyval.scntl)->fetiInfo.c_normalize = bool((yyvsp[-1].ival)); }
#line 12347 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12360 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1069:
#line 4393 "p.y" /* yacc.c:1646  */
    { if(!(yyvsp[-1].copt).lagrangeMult && (yyvsp[-1].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[-1].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[-1].copt).constraint_hess; 
          domain->solInfo().constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps; }
#line 12370 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1070:
#line 4399 "p.y" /* yacc.c:1646  */
    { if(!(yyvsp[-1].copt).lagrangeMult && (yyvsp[-1].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[-1].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[-1].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[-1].copt).constraint_hess;
          domain->solInfo().constraint_hess_eps = (yyvsp[-1].copt).constraint_hess_eps; }
#line 12380 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1071:
#line 4407 "p.y" /* yacc.c:1646  */
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        }
#line 12391 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1072:
#line 4414 "p.y" /* yacc.c:1646  */
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[0].fval); }
#line 12402 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12414 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12427 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12441 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1076:
#line 4448 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12451 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1077:
#line 4454 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12461 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1078:
#line 4460 "p.y" /* yacc.c:1646  */
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12471 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1079:
#line 4466 "p.y" /* yacc.c:1646  */
    { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[0].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
#line 12481 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1080:
#line 4472 "p.y" /* yacc.c:1646  */
    { (yyval.copt).constraint_hess = (yyvsp[0].ival);
          (yyval.copt).constraint_hess_eps = 0; }
#line 12488 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1081:
#line 4475 "p.y" /* yacc.c:1646  */
    { (yyval.copt).constraint_hess = (yyvsp[-1].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[0].fval); }
#line 12495 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1082:
#line 4480 "p.y" /* yacc.c:1646  */
    { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
#line 12504 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1083:
#line 4497 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-1].ival);
          domain->curvatureFlag = 0;
        }
#line 12513 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1084:
#line 4502 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-2].ival);
          domain->curvatureConst1 = (yyvsp[-1].fval);
          domain->curvatureFlag = 1;
        }
#line 12523 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1085:
#line 4508 "p.y" /* yacc.c:1646  */
    {
          domain->sommerfeldType = (yyvsp[-3].ival);
          domain->curvatureConst1 = (yyvsp[-2].fval);
          domain->curvatureConst2 = (yyvsp[-1].fval);
          domain->curvatureFlag = 2;
        }
#line 12534 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1086:
#line 4515 "p.y" /* yacc.c:1646  */
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        }
#line 12543 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1087:
#line 4520 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
#line 12552 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1088:
#line 4525 "p.y" /* yacc.c:1646  */
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
#line 12561 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1094:
#line 4539 "p.y" /* yacc.c:1646  */
    { domain->setWaveDirections(0, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12567 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1095:
#line 4542 "p.y" /* yacc.c:1646  */
    {
          domain->setKirchhoffLocations((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 12575 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1096:
#line 4545 "p.y" /* yacc.c:1646  */
    {
          domain->setKirchhoffLocations((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval));
        }
#line 12583 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1099:
#line 4555 "p.y" /* yacc.c:1646  */
    { domain->setFFPDirections((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12589 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12602 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12615 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1104:
#line 4587 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        }
#line 12623 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1105:
#line 4593 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 12629 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1106:
#line 4595 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 12635 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12651 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1108:
#line 4610 "p.y" /* yacc.c:1646  */
    {}
#line 12657 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1109:
#line 4612 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength; }
#line 12664 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1110:
#line 4615 "p.y" /* yacc.c:1646  */
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam; }
#line 12673 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1111:
#line 4620 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linearelastic = 1; }
#line 12679 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1112:
#line 4622 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linearelastic = (yyvsp[-1].ival); }
#line 12685 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1113:
#line 4624 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[-1].ival); }
#line 12691 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1114:
#line 4626 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-1].fval); }
#line 12697 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1115:
#line 4628 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[-1].fval); }
#line 12704 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1116:
#line 4631 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[-4].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[-1].fval); }
#line 12713 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1117:
#line 4636 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[-1].fval); }
#line 12720 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1118:
#line 4639 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[-4].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[-2].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[-1].ival); }
#line 12729 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1119:
#line 4644 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[-1].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[-1].ival); }
#line 12736 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1120:
#line 4647 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[-2].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[-1].ival); }
#line 12743 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1121:
#line 4650 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().unsymmetric = true; }
#line 12749 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1122:
#line 4652 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[-1].fval); }
#line 12755 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1123:
#line 4654 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-1].ival); }
#line 12761 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1124:
#line 4656 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-2].ival); 
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-1].ival); }
#line 12768 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1125:
#line 4659 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-3].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-2].ival);
          // note: currently we use either c1 or c2, but never both
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-1].fval);
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-1].fval); }
#line 12778 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1126:
#line 4665 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-4].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-3].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-2].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[-1].fval); }
#line 12788 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1127:
#line 4671 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[-5].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[-4].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[-3].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[-3].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[-2].fval);
          domain->solInfo().getNLInfo().linesearch.verbose = bool((yyvsp[-1].ival)); }
#line 12799 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1128:
#line 4678 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().failsafe = true; }
#line 12805 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1129:
#line 4680 "p.y" /* yacc.c:1646  */
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[-1].fval); }
#line 12812 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1130:
#line 4683 "p.y" /* yacc.c:1646  */
    { domain->solInfo().num_penalty_its = (yyvsp[-3].ival); 
          domain->solInfo().penalty_tol = (yyvsp[-2].fval);
          domain->solInfo().penalty_beta = (yyvsp[-1].fval); }
#line 12820 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1131:
#line 4687 "p.y" /* yacc.c:1646  */
    { domain->solInfo().num_penalty_its = (yyvsp[-5].ival);
          domain->solInfo().penalty_tol = (yyvsp[-4].fval);
          domain->solInfo().penalty_beta = (yyvsp[-3].fval);
          domain->solInfo().reinit_lm = bool((yyvsp[-2].ival));
          domain->solInfo().lm_update_flag = (yyvsp[-1].ival); }
#line 12830 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1132:
#line 4693 "p.y" /* yacc.c:1646  */
    { domain->solInfo().numberOfRomCPUs = (yyvsp[-1].ival); }
#line 12836 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1134:
#line 4698 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-1].ival)); 
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear; 
        }
#line 12845 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1135:
#line 4703 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-2].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[-1].ival);
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
#line 12855 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1136:
#line 4709 "p.y" /* yacc.c:1646  */
    {
          domain->solInfo().setNewton((yyvsp[-3].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[-2].ival);
          domain->solInfo().piecewise_contact = bool((yyvsp[-1].ival));
          domain->solInfo().solvercntl->fetiInfo.type  = FetiInfo::nonlinear;
        }
#line 12866 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1137:
#line 4718 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setReOrtho(); }
#line 12872 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1139:
#line 4723 "p.y" /* yacc.c:1646  */
    { geoSource->setControl((yyvsp[-7].strval),(yyvsp[-3].strval),(yyvsp[-1].strval)); domain->solInfo().soltyp = (yyvsp[-5].ival); }
#line 12878 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1140:
#line 4725 "p.y" /* yacc.c:1646  */
    { geoSource->setControl((yyvsp[-9].strval),(yyvsp[-5].strval),(yyvsp[-3].strval),(yyvsp[-1].strval)); domain->solInfo().soltyp = (yyvsp[-7].ival); }
#line 12884 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1142:
#line 4736 "p.y" /* yacc.c:1646  */
    { domain->solInfo().contact_mode = (yyvsp[-1].ival); }
#line 12890 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1143:
#line 4740 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-5].ival)-1, (yyvsp[-4].ival)-1, (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)); }
#line 12896 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1144:
#line 4743 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval));}
#line 12902 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1145:
#line 4746 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-7].ival)-1, (yyvsp[-6].ival)-1, (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), 0.0, (yyvsp[-1].ival));}
#line 12908 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1146:
#line 4748 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-6].ival)-1, (yyvsp[-5].ival)-1, (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), 0.0, -1, (yyvsp[-1].copt).lagrangeMult, (yyvsp[-1].copt).penalty);}
#line 12914 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1147:
#line 4751 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-9].ival)-1, (yyvsp[-8].ival)-1, (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-1].ival));}
#line 12920 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1148:
#line 4753 "p.y" /* yacc.c:1646  */
    { domain->addNodalCTC((yyvsp[-10].ival)-1, (yyvsp[-9].ival)-1, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-4].fval), (yyvsp[-2].ival), (yyvsp[-1].copt).lagrangeMult, (yyvsp[-1].copt).penalty);}
#line 12926 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1150:
#line 4758 "p.y" /* yacc.c:1646  */
    { 
           geoSource->addMaterial((yyvsp[-7].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12935 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1151:
#line 4763 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new BilinPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12944 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1152:
#line 4768 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new BilinPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 12953 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12967 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12982 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 12997 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1156:
#line 4805 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13006 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1157:
#line 4810 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13015 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1158:
#line 4815 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13024 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 13038 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 13053 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 13068 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1162:
#line 4852 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13077 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1163:
#line 4857 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13086 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1164:
#line 4862 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13095 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1165:
#line 4867 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-11].ival)-1, new LogStrainPlasKinHardMat((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13109 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1166:
#line 4877 "p.y" /* yacc.c:1646  */
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
#line 13124 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1167:
#line 4888 "p.y" /* yacc.c:1646  */
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
#line 13139 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1168:
#line 4899 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new ElaLinIsoMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13148 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1169:
#line 4904 "p.y" /* yacc.c:1646  */
    { 
           geoSource->addMaterial((yyvsp[-5].ival)-1, 
             new ElaLinIsoMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
	 }
#line 13157 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1170:
#line 4909 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new ElaLinIsoMat((yyvsp[-1].fval)));
         }
#line 13166 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1171:
#line 4914 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13176 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1172:
#line 4920 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13186 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1173:
#line 4926 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<ElaLinIsoMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13196 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1174:
#line 4932 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13205 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1175:
#line 4937 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13214 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1176:
#line 4942 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new StVenantKirchhoffMat((yyvsp[-1].fval)));
         }
#line 13223 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1177:
#line 4947 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13233 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1178:
#line 4953 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13243 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1179:
#line 4959 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<StVenantKirchhoffMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13253 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1180:
#line 4965 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new HenckyMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13262 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1181:
#line 4970 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new HenckyMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13271 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1182:
#line 4975 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-3].ival)-1,
             new HenckyMat((yyvsp[-1].fval)));
         }
#line 13280 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1183:
#line 4980 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13290 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1184:
#line 4986 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13300 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1185:
#line 4992 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new BrittleFractureTB<HenckyMat>((yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13310 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1186:
#line 4998 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new ElaLinIsoMat2D((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0));
         }
#line 13319 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1187:
#line 5003 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new ElaLinIsoMat2D((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13328 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1188:
#line 5008 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0));
         }
#line 13337 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1189:
#line 5013 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13346 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1190:
#line 5018 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FabricMap((yyvsp[-5].fval), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].fval), 0, 0, FabricMap::GREEN_LAGRANGE));
         }
#line 13355 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1191:
#line 5023 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMap((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMap::GREEN_LAGRANGE));
         }
#line 13364 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1192:
#line 5028 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new FabricMap((yyvsp[-5].fval), (yyvsp[-4].ival), (yyvsp[-3].ival), (yyvsp[-2].ival), (yyvsp[-1].fval), 0, 0, FabricMap::INFINTESIMAL));
         }
#line 13373 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1193:
#line 5033 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMap((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMap::INFINTESIMAL));
         }
#line 13382 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1194:
#line 5038 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMat((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0, FabricMat::GREEN_LAGRANGE));
         }
#line 13391 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1195:
#line 5043 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new FabricMat((yyvsp[-9].fval), (yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMat::GREEN_LAGRANGE));
         }
#line 13400 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1196:
#line 5048 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new FabricMat((yyvsp[-7].fval), (yyvsp[-6].ival), (yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0, 0, FabricMat::INFINTESIMAL));
         }
#line 13409 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1197:
#line 5053 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new FabricMat((yyvsp[-9].fval), (yyvsp[-8].ival), (yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), FabricMat::INFINTESIMAL));
         }
#line 13418 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1198:
#line 5058 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<ElaLinIsoMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13427 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1199:
#line 5063 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<ElaLinIsoMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13436 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1200:
#line 5068 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13446 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1201:
#line 5074 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PlaneStressMat<BrittleFractureTB<ElaLinIsoMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13456 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1202:
#line 5080 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<StVenantKirchhoffMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13465 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1203:
#line 5085 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<StVenantKirchhoffMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13474 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1204:
#line 5090 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13484 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1205:
#line 5096 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PlaneStressMat<BrittleFractureTB<StVenantKirchhoffMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13494 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1206:
#line 5102 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new PlaneStressMat<NeoHookeanMat>((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13503 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1207:
#line 5107 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new PlaneStressMat<BrittleFractureTB<NeoHookeanMat> >((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13513 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1208:
#line 5113 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new PlaneStressMat<MooneyRivlinMat>((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13522 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1209:
#line 5118 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<BrittleFractureTB<MooneyRivlinMat> >((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
           domain->solInfo().elementDeletion = true;
         }
#line 13532 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1210:
#line 5124 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13541 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1211:
#line 5129 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13550 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1212:
#line 5134 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13559 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1213:
#line 5139 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval), (yyvsp[-2].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<BilinPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13573 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1214:
#line 5149 "p.y" /* yacc.c:1646  */
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
#line 13588 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1215:
#line 5160 "p.y" /* yacc.c:1646  */
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
#line 13603 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1216:
#line 5171 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13612 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1217:
#line 5176 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13621 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1218:
#line 5181 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-11].ival)-1,
             new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 13630 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1219:
#line 5186 "p.y" /* yacc.c:1646  */
    {
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-1].fval), (yyvsp[-2].fval)) );
             domain->solInfo().elementDeletion = true;
           }
           else {
             geoSource->addMaterial((yyvsp[-12].ival)-1, new PlaneStressMat<FiniteStrainPlasKinHardMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval)) );
           }
         }
#line 13644 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1220:
#line 5196 "p.y" /* yacc.c:1646  */
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
#line 13659 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1221:
#line 5207 "p.y" /* yacc.c:1646  */
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
#line 13674 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1222:
#line 5218 "p.y" /* yacc.c:1646  */
    {
            double params[3] = { (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-5].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          }
#line 13684 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1223:
#line 5224 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<ElaLinIsoMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13697 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1224:
#line 5233 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PronyViscoElastic<ElaLinIsoMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13710 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1225:
#line 5242 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13724 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1226:
#line 5252 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13738 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1227:
#line 5262 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PronyViscoElastic<ElaLinIsoMat2D>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13751 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1228:
#line 5271 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-14].ival)-1,
             new PronyViscoElastic<ElaLinIsoMat2D>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13764 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1229:
#line 5280 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-12].ival)-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13777 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1230:
#line 5289 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-14].ival)-1,
             new PronyViscoElastic<StVenantKirchhoffMat2D>((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13790 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1231:
#line 5298 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-13].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-11].fval), (yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-7].fval), 0, 0, FabricMap::GREEN_LAGRANGE, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13803 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1232:
#line 5307 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMap::GREEN_LAGRANGE, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13816 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1233:
#line 5316 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-13].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-11].fval), (yyvsp[-10].ival), (yyvsp[-9].ival), (yyvsp[-8].ival), (yyvsp[-7].fval), 0, 0, FabricMap::INFINTESIMAL, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13829 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1234:
#line 5325 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-2].fval);
           double gtwo   = (yyvsp[-4].fval);
           double gone   = (yyvsp[-6].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMap>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMap::INFINTESIMAL, ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 13842 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1235:
#line 5334 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, FabricMat::GREEN_LAGRANGE, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13855 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1236:
#line 5343 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-17].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-15].fval), (yyvsp[-14].ival), (yyvsp[-13].ival), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMat::GREEN_LAGRANGE, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13868 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1237:
#line 5352 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-15].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-13].fval), (yyvsp[-12].ival), (yyvsp[-11].ival), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), 0, 0, FabricMat::INFINTESIMAL, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13881 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1238:
#line 5361 "p.y" /* yacc.c:1646  */
    {
           double gthree = (yyvsp[-4].fval);
           double gtwo   = (yyvsp[-6].fval);
           double gone   = (yyvsp[-8].fval);
           double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
           geoSource->addMaterial((yyvsp[-17].ival)-1,
             new PronyViscoElastic<FabricMat>((yyvsp[-15].fval), (yyvsp[-14].ival), (yyvsp[-13].ival), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), FabricMat::INFINTESIMAL, ginf, (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval)));
         }
#line 13894 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1239:
#line 5370 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13907 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1240:
#line 5379 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-14].ival)-1,
              new PlaneStressMat<PronyViscoElastic<ElaLinIsoMat> >((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13920 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1241:
#line 5388 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval); 
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13934 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1242:
#line 5398 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-18].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<ElaLinIsoMat> > >((yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13948 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1243:
#line 5408 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<StVenantKirchhoffMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13961 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1244:
#line 5417 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PronyViscoElastic<StVenantKirchhoffMat>((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 13974 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1245:
#line 5426 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 13988 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1246:
#line 5436 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14002 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1247:
#line 5446 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14015 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1248:
#line 5455 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-14].ival)-1,
              new PlaneStressMat<PronyViscoElastic<StVenantKirchhoffMat> >((yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14028 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1249:
#line 5464 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14042 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1250:
#line 5474 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-18].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<StVenantKirchhoffMat> > >((yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14056 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1251:
#line 5484 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-5].ival)-1,
              new NeoHookeanMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14065 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1252:
#line 5489 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new BrittleFractureTB<NeoHookeanMat>((yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14075 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1253:
#line 5495 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new PronyViscoElastic<NeoHookeanMat>((yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14088 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1254:
#line 5504 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-15].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> >((yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14102 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1255:
#line 5514 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PlaneStressMat<PronyViscoElastic<NeoHookeanMat> >((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14115 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1256:
#line 5523 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<NeoHookeanMat> > >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14129 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1257:
#line 5533 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-6].ival)-1,
              new MooneyRivlinMat((yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14138 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1258:
#line 5538 "p.y" /* yacc.c:1646  */
    {
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new BrittleFractureTB<MooneyRivlinMat>((yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14148 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1259:
#line 5544 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-12].ival)-1,
              new PronyViscoElastic<MooneyRivlinMat>((yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14161 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1260:
#line 5553 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-16].ival)-1,
              new BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> >((yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14175 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1261:
#line 5563 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-2].fval);
            double gtwo   = (yyvsp[-4].fval);
            double gone   = (yyvsp[-6].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-13].ival)-1,
              new PlaneStressMat<PronyViscoElastic<MooneyRivlinMat> >((yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), ginf, (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
          }
#line 14188 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1262:
#line 5572 "p.y" /* yacc.c:1646  */
    {
            double gthree = (yyvsp[-6].fval);
            double gtwo   = (yyvsp[-8].fval);
            double gone   = (yyvsp[-10].fval);
            double ginf   = 1.0 - (gone + gtwo + gthree); // use convention that prony series sums to 1
            geoSource->addMaterial((yyvsp[-17].ival)-1,
              new PlaneStressMat<BrittleFractureTB<PronyViscoElastic<MooneyRivlinMat> > >((yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), ginf, (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), (yyvsp[-5].fval)));
            domain->solInfo().elementDeletion = true;
          }
#line 14202 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1263:
#line 5582 "p.y" /* yacc.c:1646  */
    {
            double mu[1] = { (yyvsp[-3].fval) };
            double alpha[1] = { (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-6].ival)-1, new OgdenMat((yyvsp[-4].fval), mu, alpha, K));
          }
#line 14213 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1264:
#line 5589 "p.y" /* yacc.c:1646  */
    {
            double mu[2] = { (yyvsp[-5].fval), (yyvsp[-4].fval) };
            double alpha[2] = { (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-8].ival)-1, new OgdenMat((yyvsp[-6].fval), mu, alpha, K));
          }
#line 14224 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1265:
#line 5596 "p.y" /* yacc.c:1646  */
    {
            double mu[2] = { (yyvsp[-6].fval), (yyvsp[-5].fval) }; 
            double alpha[2] = { (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-9].ival)-1, new OgdenMat((yyvsp[-7].fval), mu, alpha, K));
          }
#line 14235 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1266:
#line 5603 "p.y" /* yacc.c:1646  */
    {
            double mu[3] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval) };
            double alpha[3] = { (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-10].ival)-1, new OgdenMat((yyvsp[-8].fval), mu, alpha, K));
          }
#line 14246 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1267:
#line 5610 "p.y" /* yacc.c:1646  */
    {
            double mu[3] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval) };
            double alpha[3] = { (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-11].ival)-1, new OgdenMat((yyvsp[-9].fval), mu, alpha, K));
          }
#line 14257 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1268:
#line 5617 "p.y" /* yacc.c:1646  */
    {
            double mu[4] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval) };
            double alpha[4] = { (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-12].ival)-1, new OgdenMat((yyvsp[-10].fval), mu, alpha, K));
          }
#line 14268 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1269:
#line 5624 "p.y" /* yacc.c:1646  */
    {
            double mu[4] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval) };
            double alpha[4] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-13].ival)-1, new OgdenMat((yyvsp[-11].fval), mu, alpha, K));
          }
#line 14279 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1270:
#line 5631 "p.y" /* yacc.c:1646  */
    {
            double mu[5] = { (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval) };
            double alpha[5] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-14].ival)-1, new OgdenMat((yyvsp[-12].fval), mu, alpha, K));
          }
#line 14290 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1271:
#line 5638 "p.y" /* yacc.c:1646  */
    {
            double mu[5] = { (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval) };
            double alpha[5] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-15].ival)-1, new OgdenMat((yyvsp[-13].fval), mu, alpha, K));
          }
#line 14301 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1272:
#line 5645 "p.y" /* yacc.c:1646  */
    {
            double mu[6] = { (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval) };
            double alpha[6] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-16].ival)-1, new OgdenMat((yyvsp[-14].fval), mu, alpha, K));
          }
#line 14312 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1273:
#line 5652 "p.y" /* yacc.c:1646  */
    {
            double mu[6] = { (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval) };
            double alpha[6] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-17].ival)-1, new OgdenMat((yyvsp[-15].fval), mu, alpha, K));
          }
#line 14323 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1274:
#line 5659 "p.y" /* yacc.c:1646  */
    {
            double mu[7] = { (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval) };
            double alpha[7] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-18].ival)-1, new OgdenMat((yyvsp[-16].fval), mu, alpha, K));
          }
#line 14334 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1275:
#line 5666 "p.y" /* yacc.c:1646  */
    {
            double mu[7] = { (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval) };
            double alpha[7] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-19].ival)-1, new OgdenMat((yyvsp[-17].fval), mu, alpha, K));
          }
#line 14345 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1276:
#line 5673 "p.y" /* yacc.c:1646  */
    {
            double mu[8] = { (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval), (yyvsp[-10].fval) };
            double alpha[8] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-20].ival)-1, new OgdenMat((yyvsp[-18].fval), mu, alpha, K));
          }
#line 14356 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1277:
#line 5680 "p.y" /* yacc.c:1646  */
    {
            double mu[8] = { (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval) };
            double alpha[8] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-21].ival)-1, new OgdenMat((yyvsp[-19].fval), mu, alpha, K));
          }
#line 14367 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1278:
#line 5687 "p.y" /* yacc.c:1646  */
    {
            double mu[9] = { (yyvsp[-19].fval), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval), (yyvsp[-11].fval) };
            double alpha[9] = { (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval) };
            double K[1] = { (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-22].ival)-1, new OgdenMat((yyvsp[-20].fval), mu, alpha, K));
          }
#line 14378 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1279:
#line 5694 "p.y" /* yacc.c:1646  */
    {
            double mu[9] = { (yyvsp[-20].fval), (yyvsp[-19].fval), (yyvsp[-18].fval), (yyvsp[-17].fval), (yyvsp[-16].fval), (yyvsp[-15].fval), (yyvsp[-14].fval), (yyvsp[-13].fval), (yyvsp[-12].fval) };
            double alpha[9] = { (yyvsp[-11].fval), (yyvsp[-10].fval), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval) };
            double K[2] = { (yyvsp[-2].fval), (yyvsp[-1].fval) };
            geoSource->addMaterial((yyvsp[-23].ival)-1, new OgdenMat((yyvsp[-21].fval), mu, alpha, K));
          }
#line 14389 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1280:
#line 5701 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new SimoElasticMat((yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14398 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1281:
#line 5706 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new SimoPlasticMat((yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 14407 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1282:
#line 5711 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new SimoPlasticMat((yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)) );
         }
#line 14416 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1283:
#line 5716 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-8].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
#line 14426 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1284:
#line 5722 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
#line 14436 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1285:
#line 5728 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0. };
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14449 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1286:
#line 5737 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), double((yyvsp[-1].ival)) };
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
            if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14462 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1287:
#line 5746 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 1.0e-6, std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-8].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
#line 14472 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1288:
#line 5752 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), std::numeric_limits<double>::infinity(), 0. };
            geoSource->addMaterial((yyvsp[-9].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
#line 14482 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1289:
#line 5758 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval), 0. };
            geoSource->addMaterial((yyvsp[-10].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14495 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1290:
#line 5767 "p.y" /* yacc.c:1646  */
    {
            double params[9] = { (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), double((yyvsp[-1].ival)) };
            geoSource->addMaterial((yyvsp[-11].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
            if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
              domain->solInfo().elementDeletion = true;
            }
          }
#line 14508 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1291:
#line 5776 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-5].ival)-1,
             new ExpMat((yyvsp[-4].ival), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14517 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1292:
#line 5781 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-6].ival)-1,
             new ExpMat((yyvsp[-5].ival), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14526 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1293:
#line 5786 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-7].ival)-1,
             new ExpMat((yyvsp[-6].ival), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14535 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1294:
#line 5791 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-8].ival)-1,
             new ExpMat((yyvsp[-7].ival), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14544 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1295:
#line 5796 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-9].ival)-1,
             new ExpMat((yyvsp[-8].ival), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
         }
#line 14553 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1296:
#line 5801 "p.y" /* yacc.c:1646  */
    {
           geoSource->addMaterial((yyvsp[-10].ival)-1,
             new ExpMat((yyvsp[-9].ival), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval), (yyvsp[-1].fval)));
           if((yyvsp[-1].fval) > 0 && (yyvsp[-1].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14565 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1297:
#line 5809 "p.y" /* yacc.c:1646  */
    {
           ExpMat *mat = new ExpMat((yyvsp[-10].ival), (yyvsp[-9].fval), (yyvsp[-8].fval), (yyvsp[-7].fval), (yyvsp[-6].fval), (yyvsp[-5].fval), (yyvsp[-4].fval), (yyvsp[-3].fval), (yyvsp[-2].fval));
           mat->yssrtid = (yyvsp[-1].ival);
           geoSource->addMaterial((yyvsp[-11].ival)-1, mat);
           if((yyvsp[-2].fval) > 0 && (yyvsp[-2].fval) < std::numeric_limits<double>::infinity()) {
             domain->solInfo().elementDeletion = true;
           }
         }
#line 14578 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1298:
#line 5819 "p.y" /* yacc.c:1646  */
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
#line 14600 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1299:
#line 5837 "p.y" /* yacc.c:1646  */
    {
	   geoSource->loadMaterial((yyvsp[-2].strval), (yyvsp[-1].strval));
	 }
#line 14608 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1300:
#line 5841 "p.y" /* yacc.c:1646  */
    {
	   geoSource->addMaterial((yyvsp[-3].ival)-1, (yyvsp[-2].strval), (yyvsp[-1].dlist));
	 }
#line 14616 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1302:
#line 5848 "p.y" /* yacc.c:1646  */
    { geoSource->setMatUsage((yyvsp[-2].ival)-1, (yyvsp[-1].ival)-1); }
#line 14622 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1303:
#line 5850 "p.y" /* yacc.c:1646  */
    {
            for(int i = (yyvsp[-3].ival)-1; i < (yyvsp[-2].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[-1].ival)-1);
	  }
#line 14631 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1304:
#line 5856 "p.y" /* yacc.c:1646  */
    { (yyval.dlist).nval = 0; }
#line 14637 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1305:
#line 5858 "p.y" /* yacc.c:1646  */
    { 
          if((yyvsp[-1].dlist).nval == 64) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[-1].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[0].fval);
 	}
#line 14650 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1306:
#line 5868 "p.y" /* yacc.c:1646  */
    { (yyval.slist).nval = 0; }
#line 14656 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1307:
#line 5870 "p.y" /* yacc.c:1646  */
    { 
          if((yyvsp[-1].slist).nval == 32) {
             fprintf(stderr, "Too many files!\n");
	     exit(-1);
          }
          (yyval.slist) = (yyvsp[-1].slist);
          (yyval.slist).v[(yyval.slist).nval++] = (yyvsp[0].strval);
 	}
#line 14669 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1308:
#line 5881 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-1].ival));
          domain->solInfo().setSparseRenum((yyvsp[-1].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[-1].ival)); }
#line 14677 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1309:
#line 5885 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-3].ival));
          domain->solInfo().setSparseRenum((yyvsp[-1].ival)); }
#line 14684 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1310:
#line 5888 "p.y" /* yacc.c:1646  */
    { domain->solInfo().setRenum((yyvsp[-5].ival));
          domain->solInfo().setSparseRenum((yyvsp[-3].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[-1].ival)); }
#line 14692 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1311:
#line 5895 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().svdPodRom = true;}
#line 14700 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1313:
#line 5903 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().snapfiPodRom.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14706 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1314:
#line 5914 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().velocPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14712 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1315:
#line 5916 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().accelPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14718 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1316:
#line 5918 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().dsvPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14724 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1317:
#line 5920 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().muvPodRomFile.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14730 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1318:
#line 5922 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[0].ival); }
#line 14736 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1319:
#line 5924 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14742 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1320:
#line 5926 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14748 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1321:
#line 5928 "p.y" /* yacc.c:1646  */
    { domain->solInfo().normalize = (yyvsp[0].ival); }
#line 14754 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1322:
#line 5930 "p.y" /* yacc.c:1646  */
    { domain->solInfo().subtractRefPodRom = true;
    domain->solInfo().readInLocalBasesCent.push_back(std::string((yyvsp[0].strval))); }
#line 14761 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1323:
#line 5933 "p.y" /* yacc.c:1646  */
    { domain->solInfo().flagss = (yyvsp[0].ival); }
#line 14767 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1324:
#line 5935 "p.y" /* yacc.c:1646  */
    { domain->solInfo().flagss = (yyvsp[-1].ival); domain->solInfo().flagrs = (yyvsp[0].ival); }
#line 14773 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1325:
#line 5937 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[0].ival); }
#line 14779 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1326:
#line 5939 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[-1].ival);
    domain->solInfo().skipOffSet = (yyvsp[0].ival); }
#line 14786 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1327:
#line 5942 "p.y" /* yacc.c:1646  */
    { domain->solInfo().robcSolve = bool((yyvsp[0].ival)); }
#line 14792 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1328:
#line 5944 "p.y" /* yacc.c:1646  */
    { for(int i=0; i<(yyvsp[0].slist).nval; ++i) domain->solInfo().robfi.push_back(std::string((yyvsp[0].slist).v[i])); }
#line 14798 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1329:
#line 5946 "p.y" /* yacc.c:1646  */
    { domain->solInfo().svdBlockSize = (yyvsp[0].ival); }
#line 14804 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1330:
#line 5948 "p.y" /* yacc.c:1646  */
    { domain->solInfo().romEnergy = (yyvsp[0].fval); }
#line 14810 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1331:
#line 5961 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf            = 3;
    domain->solInfo().nmfMaxIter         = (yyvsp[-3].ival);
    domain->solInfo().nmfTol             = (yyvsp[-2].fval);
    domain->solInfo().nmfPqnNumInnerIter = (yyvsp[-1].ival);
    domain->solInfo().nmfPqnAlpha        = (yyvsp[0].fval); }
#line 14820 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1332:
#line 5967 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = (yyvsp[-4].ival);
    domain->solInfo().nmfDelROBDim = (yyvsp[-3].ival);
    domain->solInfo().nmfRandInit  = (yyvsp[-2].ival);
    domain->solInfo().nmfMaxIter   = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 14831 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1333:
#line 5974 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 1;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol     = (yyvsp[0].fval); }
#line 14839 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1334:
#line 5978 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 4;   
    domain->solInfo().nmfMaxIter = (yyvsp[-4].ival); 
    domain->solInfo().nmfTol     = (yyvsp[-3].fval); 
    domain->solInfo().nmfcAlpha  = (yyvsp[-2].fval);
    domain->solInfo().nmfcBeta   = (yyvsp[-1].fval);
    domain->solInfo().nmfcGamma  = (yyvsp[0].fval);}
#line 14850 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1335:
#line 5985 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf    = 4;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol     = (yyvsp[0].fval); }
#line 14858 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1336:
#line 5989 "p.y" /* yacc.c:1646  */
    { domain->solInfo().nmfNumSub = (yyvsp[0].ival); }
#line 14864 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1337:
#line 5991 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 2; }
#line 14870 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1338:
#line 5993 "p.y" /* yacc.c:1646  */
    { domain->solInfo().clustering = (yyvsp[0].ival); }
#line 14876 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1339:
#line 5995 "p.y" /* yacc.c:1646  */
    { domain->solInfo().clustering = (yyvsp[-1].ival); 
    domain->solInfo().clusterSubspaceAngle = true; }
#line 14883 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1340:
#line 5998 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[0].ival); }
#line 14889 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1341:
#line 6000 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[-1].ival);
    domain->solInfo().tolPodRom = (yyvsp[0].fval);}
#line 14896 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1342:
#line 6003 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeCluster = (yyvsp[-2].ival);
    domain->solInfo().tolPodRom = (yyvsp[-1].fval);
    domain->solInfo().solverTypeSpnnls = (yyvsp[0].ival); }
#line 14904 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1343:
#line 6007 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hotstartSample = bool((yyvsp[0].ival)); }
#line 14910 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1344:
#line 6009 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rowClustering = (yyvsp[0].ival); }
#line 14916 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1345:
#line 6011 "p.y" /* yacc.c:1646  */
    { domain->solInfo().rowClustering = (yyvsp[-1].ival);
    domain->solInfo().clusterSubspaceAngle = true; }
#line 14923 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1347:
#line 6018 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().DEIMBasisPod = true; }
#line 14931 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1349:
#line 6026 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().UDEIMBasisPod = true; }
#line 14939 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1351:
#line 6034 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().samplingPodRom = true; }
#line 14947 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1356:
#line 6045 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().snapProjPodRom = true; }
#line 14955 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1358:
#line 6053 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInROBorModes.push_back((yyvsp[0].strval)); }
#line 14961 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1359:
#line 6055 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInROBorModes.push_back((yyvsp[-1].strval));
    domain->solInfo().localBasisSize.push_back((yyvsp[0].ival)); }
#line 14968 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1360:
#line 6058 "p.y" /* yacc.c:1646  */
    { domain->solInfo().readInDualROB.push_back((yyvsp[-1].strval));
    domain->solInfo().localDualBasisSize.push_back((yyvsp[0].ival)); }
#line 14975 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1361:
#line 6061 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[0].strval)); }
#line 14981 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1362:
#line 6063 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[-1].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[0].strval)); }
#line 14988 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1363:
#line 6066 "p.y" /* yacc.c:1646  */
    { domain->solInfo().statePodRomFile.push_back((yyvsp[-2].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[-1].strval));
    domain->solInfo().accelPodRomFile.push_back((yyvsp[0].strval)); }
#line 14996 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1364:
#line 6070 "p.y" /* yacc.c:1646  */
    { domain->solInfo().tolPodRom = (yyvsp[0].fval); }
#line 15002 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1365:
#line 6072 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipPodRom = (yyvsp[0].ival); }
#line 15008 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1366:
#line 6074 "p.y" /* yacc.c:1646  */
    { domain->solInfo().randomSampleSize = (yyvsp[0].ival); 
    domain->solInfo().randomVecSampling = true; }
#line 15015 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1367:
#line 6077 "p.y" /* yacc.c:1646  */
    { domain->solInfo().skipOffSet = (yyvsp[0].ival); }
#line 15021 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1368:
#line 6079 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[0].ival); }
#line 15027 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1369:
#line 6081 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[-1].ival); 
    domain->solInfo().forcePodSize = (yyvsp[0].ival);}
#line 15034 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1370:
#line 6084 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizePodRom = (yyvsp[-2].ival); 
    domain->solInfo().forcePodSize = (yyvsp[-1].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); }
#line 15042 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1371:
#line 6088 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[0].ival)); }
#line 15048 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1372:
#line 6090 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useConstantMassForces = bool((yyvsp[0].ival)); }
#line 15054 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1373:
#line 6092 "p.y" /* yacc.c:1646  */
    { domain->solInfo().stackedElementSampling = bool((yyvsp[0].ival)); }
#line 15060 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1374:
#line 6094 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useMassOrthogonalProjection = bool((yyvsp[0].ival)); }
#line 15066 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1375:
#line 6096 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = bool((yyvsp[0].ival)); }
#line 15072 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1376:
#line 6098 "p.y" /* yacc.c:1646  */
    { domain->solInfo().reduceFollower = bool((yyvsp[0].ival)); }
#line 15078 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1377:
#line 6100 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15084 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1378:
#line 6102 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[-1].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15091 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1379:
#line 6105 "p.y" /* yacc.c:1646  */
    { domain->solInfo().PODerrornorm.push_back((yyvsp[-2].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[-1].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[0].strval)); }
#line 15099 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1380:
#line 6109 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useScalingSpnnls = bool((yyvsp[0].ival)); }
#line 15105 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1381:
#line 6111 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useCenterSpnnls = bool((yyvsp[0].ival)); }
#line 15111 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1382:
#line 6113 "p.y" /* yacc.c:1646  */
    { domain->solInfo().projectSolution = bool((yyvsp[0].ival)); }
#line 15117 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1383:
#line 6115 "p.y" /* yacc.c:1646  */
    { domain->solInfo().positiveElements = bool((yyvsp[0].ival)); }
#line 15123 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1384:
#line 6117 "p.y" /* yacc.c:1646  */
    { domain->solInfo().hotstartSample = bool((yyvsp[0].ival)); }
#line 15129 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1385:
#line 6119 "p.y" /* yacc.c:1646  */
    { domain->solInfo().solverTypeSpnnls = (yyvsp[0].ival); }
#line 15135 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1386:
#line 6121 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxSizeSpnnls = (yyvsp[0].fval); }
#line 15141 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1387:
#line 6123 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxElemSpnnls = (yyvsp[0].ival); }
#line 15147 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1388:
#line 6125 "p.y" /* yacc.c:1646  */
    { domain->solInfo().maxIterSpnnls = (yyvsp[0].fval); }
#line 15153 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1389:
#line 6127 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[0].strval); }
#line 15159 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1390:
#line 6129 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[-1].strval);
    domain->solInfo().forcePodSize = (yyvsp[0].ival); }
#line 15166 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1391:
#line 6132 "p.y" /* yacc.c:1646  */
    { domain->solInfo().forcePodRomFile = (yyvsp[-2].strval); 
    domain->solInfo().forcePodSize = (yyvsp[-1].ival); 
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); }
#line 15174 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1392:
#line 6136 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[0].strval); 
    domain->solInfo().ConstraintBasisPod = true;}
#line 15181 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1393:
#line 6139 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[-1].strval);
    domain->solInfo().constraintPodSize = (yyvsp[0].ival); 
    domain->solInfo().ConstraintBasisPod = true; }
#line 15189 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1394:
#line 6143 "p.y" /* yacc.c:1646  */
    { domain->solInfo().constraintPodRomFile = (yyvsp[-2].strval);
    domain->solInfo().constraintPodSize = (yyvsp[-1].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[0].ival); 
    domain->solInfo().ConstraintBasisPod = true; }
#line 15198 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1395:
#line 6148 "p.y" /* yacc.c:1646  */
    { domain->solInfo().filterSnapshotRows = bool((yyvsp[0].ival)); }
#line 15204 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1396:
#line 6150 "p.y" /* yacc.c:1646  */
    { domain->solInfo().selectFullNode = bool((yyvsp[0].ival)); }
#line 15210 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1397:
#line 6152 "p.y" /* yacc.c:1646  */
    { domain->solInfo().selectFullElem = bool((yyvsp[0].ival)); }
#line 15216 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1398:
#line 6154 "p.y" /* yacc.c:1646  */
    { domain->solInfo().computeForceSnap = bool((yyvsp[0].ival)); }
#line 15222 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1399:
#line 6156 "p.y" /* yacc.c:1646  */
    { domain->solInfo().computeConstraintSnap = bool((yyvsp[0].ival)); }
#line 15228 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1400:
#line 6158 "p.y" /* yacc.c:1646  */
    { domain->solInfo().orthogForceSnap = bool((yyvsp[0].ival)); }
#line 15234 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1401:
#line 6160 "p.y" /* yacc.c:1646  */
    { domain->solInfo().orthogConstraintSnap = bool((yyvsp[0].ival)); }
#line 15240 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1402:
#line 6162 "p.y" /* yacc.c:1646  */
    { domain->solInfo().npMax = (yyvsp[0].ival); }
#line 15246 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1403:
#line 6164 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scpkMB= (yyvsp[-1].ival);
    domain->solInfo().scpkNB= (yyvsp[0].ival); }
#line 15253 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1404:
#line 6167 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scpkMP= (yyvsp[-1].ival);
    domain->solInfo().scpkNP= (yyvsp[0].ival); }
#line 15260 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1405:
#line 6170 "p.y" /* yacc.c:1646  */
    { domain->solInfo().useReverseOrder = bool((yyvsp[0].ival)); }
#line 15266 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1406:
#line 6172 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfNumROBDim = (yyvsp[-4].ival);
    domain->solInfo().nmfDelROBDim = (yyvsp[-3].ival);
    domain->solInfo().nmfRandInit = (yyvsp[-2].ival);
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 15277 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1407:
#line 6179 "p.y" /* yacc.c:1646  */
    { domain->solInfo().use_nmf = 1;
    domain->solInfo().nmfMaxIter = (yyvsp[-1].ival);
    domain->solInfo().nmfTol = (yyvsp[0].fval); }
#line 15285 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1409:
#line 6187 "p.y" /* yacc.c:1646  */
    { domain->solInfo().conwepConfigurations.push_back((yyvsp[-1].blastData)); }
#line 15291 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1410:
#line 6192 "p.y" /* yacc.c:1646  */
    { domain->solInfo().scalePosCoords = true;
     domain->solInfo().xScaleFactor = (yyvsp[-3].fval);
     domain->solInfo().yScaleFactor = (yyvsp[-2].fval);
     domain->solInfo().zScaleFactor = (yyvsp[-1].fval);}
#line 15300 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1411:
#line 6201 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePOSCFG = true; }
#line 15306 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1412:
#line 6203 "p.y" /* yacc.c:1646  */
    { domain->solInfo().xScaleFactors.push_back((yyvsp[-3].fval));
     domain->solInfo().yScaleFactors.push_back((yyvsp[-2].fval));
     domain->solInfo().zScaleFactors.push_back((yyvsp[-1].fval)); }
#line 15314 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1413:
#line 6207 "p.y" /* yacc.c:1646  */
    { domain->solInfo().MassOrthogonalBasisFiles.push_back((yyvsp[-1].strval)); }
#line 15320 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1414:
#line 6212 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePOSCFG = true; }
#line 15326 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1415:
#line 6214 "p.y" /* yacc.c:1646  */
    { domain->solInfo().NodeTrainingFiles.push_back(std::string((yyvsp[-1].slist).v[0]));
     for(int i=1; i<(yyvsp[-1].slist).nval; ++i) domain->solInfo().MassOrthogonalBasisFiles.push_back(std::string((yyvsp[-1].slist).v[i])); }
#line 15333 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1416:
#line 6220 "p.y" /* yacc.c:1646  */
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().ROMPostProcess = true; }
#line 15341 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1418:
#line 6228 "p.y" /* yacc.c:1646  */
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[0].strval)); 
    domain->solInfo().numRODFile += 1; }
#line 15348 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1419:
#line 6231 "p.y" /* yacc.c:1646  */
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[-1].strval));
    domain->solInfo().numRODFile += 1; 
    domain->solInfo().skipPodRom = (yyvsp[0].ival);}
#line 15356 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1420:
#line 6235 "p.y" /* yacc.c:1646  */
    { domain->solInfo().romresidType = (yyvsp[0].ival); }
#line 15362 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1421:
#line 6240 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = (yyvsp[0].ival); }
#line 15368 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1422:
#line 6242 "p.y" /* yacc.c:1646  */
    { (yyval.ival) = std::numeric_limits<int>::max(); }
#line 15374 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1423:
#line 6247 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[0].ival); }
#line 15380 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1424:
#line 6249 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = (yyvsp[0].fval); }
#line 15386 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1425:
#line 6251 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = std::numeric_limits<double>::infinity(); }
#line 15392 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;

  case 1426:
#line 6253 "p.y" /* yacc.c:1646  */
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); }
#line 15398 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
    break;


#line 15402 "/home/anarkhede/tinkercliffs/FEMTesting/Parser.d/parser.cpp" /* yacc.c:1646  */
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
#line 6255 "p.y" /* yacc.c:1906  */

