
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 1 "p.y"

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
 extern std::string clusterData_;
 extern std::string subdomains_;
 extern std::string decomposition_;
 extern std::string connectivity_;


/* Line 189 of yacc.c  */
#line 100 "/home/tac688/codes/AEROS/Parser.d/parser.cpp"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ACTUATORS = 258,
     AERO = 259,
     AEROH = 260,
     AEROTYPE = 261,
     ALPROC = 262,
     AMAT = 263,
     ANALYSIS = 264,
     ARCLENGTH = 265,
     ATTRIBUTES = 266,
     ANGULAROUTTYPE = 267,
     AUGMENT = 268,
     AUGMENTTYPE = 269,
     AVERAGED = 270,
     ATDARB = 271,
     ACOU = 272,
     ATDDNB = 273,
     ATDROB = 274,
     ARPACK = 275,
     ATDDIR = 276,
     ATDNEU = 277,
     AXIHDIR = 278,
     AXIHNEU = 279,
     AXINUMMODES = 280,
     AXINUMSLICES = 281,
     AXIHSOMMER = 282,
     AXIMPC = 283,
     AUXCOARSESOLVER = 284,
     ACMECNTL = 285,
     ADDEDMASS = 286,
     AEROEMBED = 287,
     AUGMENTED = 288,
     BLOCKDIAG = 289,
     BOFFSET = 290,
     BUCKLE = 291,
     BGTL = 292,
     BMPC = 293,
     BINARYINPUT = 294,
     BINARYOUTPUT = 295,
     BLOCKSIZE = 296,
     CHECKTOKEN = 297,
     COARSESOLVER = 298,
     COEF = 299,
     CFRAMES = 300,
     COLLOCATEDTYPE = 301,
     CONVECTION = 302,
     COMPOSITE = 303,
     CONDITION = 304,
     CONTROL = 305,
     CORNER = 306,
     CORNERTYPE = 307,
     CURVE = 308,
     CCTTOL = 309,
     CCTSOLVER = 310,
     CRHS = 311,
     COUPLEDSCALE = 312,
     CONTACTSURFACES = 313,
     CMPC = 314,
     CNORM = 315,
     COMPLEXOUTTYPE = 316,
     CONSTRMAT = 317,
     CASES = 318,
     CONSTRAINEDSURFACES = 319,
     CSFRAMES = 320,
     CSTYPE = 321,
     CONSTANT = 322,
     CONWEP = 323,
     DAMPING = 324,
     DblConstant = 325,
     DEM = 326,
     DIMASS = 327,
     DISP = 328,
     DIRECT = 329,
     DLAMBDA = 330,
     DP = 331,
     DYNAM = 332,
     DETER = 333,
     DECOMPOSE = 334,
     DECOMPFILE = 335,
     DMPC = 336,
     DEBUGCNTL = 337,
     DEBUGICNTL = 338,
     CONSTRAINTS = 339,
     MULTIPLIERS = 340,
     PENALTY = 341,
     ELLUMP = 342,
     EIGEN = 343,
     EFRAMES = 344,
     ELSCATTERER = 345,
     END = 346,
     ELHSOMMERFELD = 347,
     ETEMP = 348,
     EXPLICIT = 349,
     EXTFOL = 350,
     EPSILON = 351,
     ELEMENTARYFUNCTIONTYPE = 352,
     FABMAT = 353,
     FACE = 354,
     FACOUSTICS = 355,
     FETI = 356,
     FETI2TYPE = 357,
     FETIPREC = 358,
     FFP = 359,
     FFPDIR = 360,
     FITALG = 361,
     FNAME = 362,
     FLUX = 363,
     FORCE = 364,
     FRONTAL = 365,
     FETIH = 366,
     FIELDWEIGHTLIST = 367,
     FILTEREIG = 368,
     FLUID = 369,
     FREQSWEEP = 370,
     FREQSWEEP1 = 371,
     FREQSWEEP2 = 372,
     FREQSWEEPA = 373,
     FSGL = 374,
     FSINTERFACE = 375,
     FSISCALING = 376,
     FSIELEMENT = 377,
     NOLOCALFSISPLITING = 378,
     FSICORNER = 379,
     FFIDEBUG = 380,
     FAILSAFE = 381,
     FRAMETYPE = 382,
     GEPS = 383,
     GLOBALTOL = 384,
     GRAVITY = 385,
     GRBM = 386,
     GTGSOLVER = 387,
     GLOBALCRBMTOL = 388,
     GROUP = 389,
     GROUPTYPE = 390,
     GOLDFARBTOL = 391,
     GOLDFARBCHECK = 392,
     HDIRICHLET = 393,
     HEAT = 394,
     HFETI = 395,
     HNEUMAN = 396,
     HSOMMERFELD = 397,
     HFTT = 398,
     HELMHOLTZ = 399,
     HNBO = 400,
     HELMMF = 401,
     HELMSO = 402,
     HSCBO = 403,
     HWIBO = 404,
     HZEM = 405,
     HZEMFILTER = 406,
     HLMPC = 407,
     HERMITIAN = 408,
     HESSIAN = 409,
     IACC = 410,
     IDENTITY = 411,
     IDIS = 412,
     IDIS6 = 413,
     IntConstant = 414,
     INTERFACELUMPED = 415,
     ITEMP = 416,
     ITERTYPE = 417,
     IVEL = 418,
     INCIDENCE = 419,
     IHDIRICHLET = 420,
     IHDSWEEP = 421,
     IHNEUMANN = 422,
     ISOLVERTYPE = 423,
     INPC = 424,
     INFINTY = 425,
     JACOBI = 426,
     KEYLETTER = 427,
     KRYLOVTYPE = 428,
     KIRLOC = 429,
     LAYC = 430,
     LAYN = 431,
     LAYD = 432,
     LAYO = 433,
     LAYMAT = 434,
     LFACTOR = 435,
     LMPC = 436,
     LOAD = 437,
     LOADCASE = 438,
     LOBPCG = 439,
     LOCALSOLVER = 440,
     LINESEARCH = 441,
     LUMPED = 442,
     MASS = 443,
     MATERIALS = 444,
     MATLAB = 445,
     MAXITR = 446,
     MAXORTHO = 447,
     MAXVEC = 448,
     MODAL = 449,
     MPCPRECNO = 450,
     MPCPRECNOID = 451,
     MPCTYPE = 452,
     MPCTYPEID = 453,
     MPCSCALING = 454,
     MPCELEMENT = 455,
     MPCBLOCKID = 456,
     MPCBLK_OVERLAP = 457,
     MFTT = 458,
     MRHS = 459,
     MPCCHECK = 460,
     MUMPSICNTL = 461,
     MUMPSCNTL = 462,
     MECH = 463,
     MODDAMP = 464,
     MODEFILTER = 465,
     MOMENTTYPE = 466,
     MPROJECT = 467,
     MAXIMUM = 468,
     NDTYPE = 469,
     NEIGPA = 470,
     NEWMARK = 471,
     NewLine = 472,
     NEWTON = 473,
     NL = 474,
     NLMAT = 475,
     NLPREC = 476,
     NOCOARSE = 477,
     NODETOKEN = 478,
     NONINPC = 479,
     NSBSPV = 480,
     NLTOL = 481,
     NUMCGM = 482,
     NOSECONDARY = 483,
     NFRAMES = 484,
     OPTIMIZATION = 485,
     OUTPUT = 486,
     OUTPUT6 = 487,
     OUTPUTFRAME = 488,
     QSTATIC = 489,
     QLOAD = 490,
     PITA = 491,
     PITADISP6 = 492,
     PITAVEL6 = 493,
     NOFORCE = 494,
     MDPITA = 495,
     GLOBALBASES = 496,
     LOCALBASES = 497,
     TIMEREVERSIBLE = 498,
     REMOTECOARSE = 499,
     ORTHOPROJTOL = 500,
     READINITSEED = 501,
     JUMPCVG = 502,
     JUMPOUTPUT = 503,
     PRECNO = 504,
     PRECONDITIONER = 505,
     PRELOAD = 506,
     PRESSURE = 507,
     PRINTMATLAB = 508,
     PROJ = 509,
     PIVOT = 510,
     PRECTYPE = 511,
     PRECTYPEID = 512,
     PICKANYCORNER = 513,
     PADEPIVOT = 514,
     PROPORTIONING = 515,
     PLOAD = 516,
     PADEPOLES = 517,
     POINTSOURCE = 518,
     PLANEWAVE = 519,
     PTOL = 520,
     PLANTOL = 521,
     PMAXIT = 522,
     PIECEWISE = 523,
     RADIATION = 524,
     RAYDAMP = 525,
     RBMFILTER = 526,
     RBMSET = 527,
     READMODE = 528,
     REBUILD = 529,
     REDFOL = 530,
     RENUM = 531,
     RENUMBERID = 532,
     REORTHO = 533,
     RESTART = 534,
     RECONS = 535,
     RECONSALG = 536,
     REBUILDCCT = 537,
     RANDOM = 538,
     RPROP = 539,
     RNORM = 540,
     REVERSENORMALS = 541,
     ROTVECOUTTYPE = 542,
     RESCALING = 543,
     SCALING = 544,
     SCALINGTYPE = 545,
     STRDAMP = 546,
     SDETAFT = 547,
     SENSORS = 548,
     SOLVERTYPE = 549,
     SHIFT = 550,
     SPOOLESTAU = 551,
     SPOOLESSEED = 552,
     SPOOLESMAXSIZE = 553,
     SPOOLESMAXDOMAINSIZE = 554,
     SPOOLESMAXZEROS = 555,
     SPOOLESMSGLVL = 556,
     SPOOLESSCALE = 557,
     SPOOLESPIVOT = 558,
     SPOOLESRENUM = 559,
     SPARSEMAXSUP = 560,
     SPARSEDEFBLK = 561,
     STATS = 562,
     STRESSID = 563,
     SUBSPACE = 564,
     SURFACE = 565,
     SAVEMEMCOARSE = 566,
     SPACEDIMENSION = 567,
     SCATTERER = 568,
     STAGTOL = 569,
     SCALED = 570,
     SWITCH = 571,
     STABLE = 572,
     SUBTYPE = 573,
     STEP = 574,
     SOWER = 575,
     SHELLTHICKNESS = 576,
     SURF = 577,
     SPRINGMAT = 578,
     TANGENT = 579,
     TDENFORCE = 580,
     TEMP = 581,
     TIME = 582,
     TOLEIG = 583,
     TOLFETI = 584,
     TOLJAC = 585,
     TOLPCG = 586,
     TOPFILE = 587,
     TOPOLOGY = 588,
     TRBM = 589,
     THERMOE = 590,
     THERMOH = 591,
     TETT = 592,
     TOLCGM = 593,
     TURKEL = 594,
     TIEDSURFACES = 595,
     THETA = 596,
     HRC = 597,
     THIRDNODE = 598,
     THERMMAT = 599,
     TDENFORC = 600,
     TESTULRICH = 601,
     THRU = 602,
     TRIVIAL = 603,
     USE = 604,
     USERDEFINEDISP = 605,
     USERDEFINEFORCE = 606,
     UPROJ = 607,
     UNSYMMETRIC = 608,
     USING = 609,
     VERSION = 610,
     WETCORNERS = 611,
     YMTT = 612,
     ZERO = 613,
     BINARY = 614,
     GEOMETRY = 615,
     DECOMPOSITION = 616,
     GLOBAL = 617,
     MATCHER = 618,
     CPUMAP = 619,
     NODALCONTACT = 620,
     MODE = 621,
     FRIC = 622,
     GAP = 623,
     OUTERLOOP = 624,
     EDGEWS = 625,
     WAVETYPE = 626,
     ORTHOTOL = 627,
     IMPE = 628,
     FREQ = 629,
     DPH = 630,
     WAVEMETHOD = 631,
     MATSPEC = 632,
     MATUSAGE = 633,
     BILINEARPLASTIC = 634,
     FINITESTRAINPLASTIC = 635,
     LINEARELASTIC = 636,
     STVENANTKIRCHHOFF = 637,
     LINPLSTRESS = 638,
     READ = 639,
     OPTCTV = 640,
     ISOTROPICLINEARELASTIC = 641,
     NEOHOOKEAN = 642,
     ISOTROPICLINEARELASTICJ2PLASTIC = 643,
     ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 644,
     HYPERELASTIC = 645,
     MOONEYRIVLIN = 646,
     HENCKY = 647,
     LOGSTRAINPLASTIC = 648,
     SVKPLSTRESS = 649,
     SURFACETOPOLOGY = 650,
     MORTARTIED = 651,
     MORTARSCALING = 652,
     MORTARINTEGRATIONRULE = 653,
     SEARCHTOL = 654,
     STDMORTAR = 655,
     DUALMORTAR = 656,
     WETINTERFACE = 657,
     NSUBS = 658,
     EXITAFTERDEC = 659,
     SKIP = 660,
     OUTPUTMEMORY = 661,
     OUTPUTWEIGHT = 662,
     SOLVER = 663,
     SPNNLSSOLVERTYPE = 664,
     MAXSIZE = 665,
     WEIGHTLIST = 666,
     GMRESRESIDUAL = 667,
     SLOSH = 668,
     SLGRAV = 669,
     SLZEM = 670,
     SLZEMFILTER = 671,
     PDIR = 672,
     HEFSB = 673,
     HEFRS = 674,
     HEINTERFACE = 675,
     SNAPFI = 676,
     PODROB = 677,
     TRNVCT = 678,
     OFFSET = 679,
     ORTHOG = 680,
     SVDTOKEN = 681,
     CONVERSIONTOKEN = 682,
     CONVFI = 683,
     SAMPLING = 684,
     SNAPSHOTPROJECT = 685,
     PODSIZEMAX = 686,
     REFSUBSTRACT = 687,
     TOLER = 688,
     NORMALIZETOKEN = 689,
     FNUMBER = 690,
     SNAPWEIGHT = 691,
     ROBFI = 692,
     STAVCT = 693,
     VELVCT = 694,
     ACCVCT = 695,
     CONWEPCFG = 696,
     PSEUDOGNAT = 697,
     PSEUDOGNATELEM = 698,
     VECTORNORM = 699,
     REBUILDFORCE = 700,
     SAMPNODESLOT = 701,
     REDUCEDSTIFFNESS = 702,
     UDEIMBASIS = 703,
     FORCEROB = 704,
     DEIMINDICES = 705,
     UDEIMINDICES = 706,
     SVDFORCESNAP = 707,
     USEMASSNORMALIZEDBASIS = 708,
     OPTSENSITIVITY = 709,
     SENSITIVITYID = 710,
     SENSITIVITYTYPE = 711,
     SENSITIVITYMETHOD = 712,
     QRFACTORIZATION = 713,
     QMATRIX = 714,
     RMATRIX = 715,
     XMATRIX = 716
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 28 "p.y"

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
 ComplexBCList *cxbclist;
 ComplexBCond cxbcval;
 FrameData frame;
 NodalFrameData nframe;
 TrivialPair<MFTTData*,int> mftval;
 TrivialPair<MFTTData*,int> hftval;
 LayerData ldata;
 LayInfo *linfo;
 CoefData coefdata;
 ComplexFDBC complexFDBC;
 ComplexFNBC complexFNBC;
 AxiMPC axiMPC;
 DoubleList dlist;
 StringList slist;
 SurfaceEntity* SurfObj;
 MortarHandler* MortarCondObj;
 LMPCTerm* mpcterm;
 GeoSource::Rprop rprop;
 OutputInfo oinfo;
 ConstraintOptions copt;
 BlastLoading::BlastData blastData;
 SensitivityInfo sinfo;



/* Line 214 of yacc.c  */
#line 636 "/home/tac688/codes/AEROS/Parser.d/parser.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 648 "/home/tac688/codes/AEROS/Parser.d/parser.cpp"

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
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
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
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
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
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
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
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  523
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   6306

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  462
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  240
/* YYNRULES -- Number of rules.  */
#define YYNRULES  1126
/* YYNRULES -- Number of states.  */
#define YYNSTATES  2769

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   716

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
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
     455,   456,   457,   458,   459,   460,   461
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     6,     8,    11,    13,    15,    17,    19,
      21,    23,    25,    27,    29,    31,    33,    35,    37,    39,
      41,    43,    45,    47,    49,    51,    53,    55,    57,    59,
      61,    63,    65,    67,    69,    71,    73,    75,    77,    79,
      81,    83,    85,    87,    89,    91,    93,    95,    97,    99,
     101,   103,   105,   107,   109,   111,   113,   115,   117,   119,
     121,   123,   125,   127,   129,   131,   133,   135,   137,   139,
     141,   143,   145,   147,   149,   151,   153,   155,   157,   159,
     161,   163,   165,   167,   169,   171,   173,   175,   177,   179,
     181,   183,   185,   187,   189,   191,   193,   195,   197,   199,
     201,   203,   205,   207,   209,   211,   213,   215,   217,   219,
     221,   223,   225,   227,   229,   231,   233,   235,   237,   239,
     241,   243,   245,   247,   249,   251,   253,   255,   257,   259,
     261,   263,   265,   267,   269,   271,   273,   275,   277,   279,
     281,   283,   285,   287,   289,   291,   293,   295,   301,   306,
     309,   315,   322,   329,   332,   339,   345,   353,   361,   370,
     384,   393,   401,   410,   419,   429,   444,   454,   458,   461,
     464,   467,   470,   474,   477,   482,   486,   491,   496,   503,
     506,   511,   517,   522,   528,   533,   538,   543,   548,   553,
     557,   560,   565,   570,   574,   578,   582,   586,   590,   594,
     598,   601,   606,   609,   614,   619,   624,   629,   632,   636,
     641,   644,   648,   653,   656,   660,   665,   671,   677,   680,
     683,   686,   689,   692,   695,   698,   701,   706,   717,   722,
     728,   732,   735,   739,   742,   746,   749,   753,   756,   768,
     782,   797,   803,   806,   809,   817,   827,   839,   852,   855,
     861,   868,   871,   877,   880,   886,   891,   896,   900,   904,
     908,   912,   916,   920,   923,   927,   931,   934,   937,   940,
     944,   948,   952,   956,   961,   966,   970,   974,   980,   985,
     991,   998,  1006,  1010,  1016,  1020,  1026,  1031,  1037,  1044,
    1052,  1056,  1059,  1063,  1066,  1069,  1073,  1076,  1079,  1082,
    1086,  1089,  1093,  1096,  1099,  1102,  1105,  1108,  1110,  1113,
    1117,  1121,  1125,  1128,  1134,  1138,  1142,  1146,  1149,  1153,
    1156,  1160,  1165,  1170,  1176,  1180,  1184,  1187,  1192,  1195,
    1199,  1202,  1206,  1209,  1215,  1218,  1221,  1224,  1227,  1230,
    1233,  1237,  1242,  1247,  1256,  1261,  1266,  1270,  1275,  1279,
    1285,  1290,  1296,  1299,  1303,  1307,  1311,  1314,  1319,  1321,
    1325,  1327,  1330,  1333,  1336,  1340,  1346,  1353,  1360,  1364,
    1369,  1376,  1382,  1386,  1391,  1395,  1397,  1400,  1407,  1410,
    1413,  1416,  1419,  1422,  1427,  1433,  1438,  1441,  1448,  1452,
    1455,  1459,  1462,  1467,  1469,  1472,  1475,  1478,  1484,  1490,
    1497,  1505,  1509,  1511,  1513,  1516,  1518,  1520,  1522,  1525,
    1527,  1529,  1532,  1534,  1539,  1545,  1549,  1553,  1560,  1566,
    1573,  1577,  1580,  1586,  1593,  1596,  1599,  1607,  1617,  1621,
    1624,  1631,  1635,  1637,  1640,  1644,  1647,  1651,  1653,  1656,
    1661,  1664,  1669,  1676,  1685,  1691,  1694,  1698,  1703,  1709,
    1712,  1716,  1723,  1726,  1730,  1737,  1741,  1743,  1746,  1751,
    1757,  1764,  1768,  1771,  1776,  1778,  1781,  1785,  1787,  1790,
    1795,  1801,  1808,  1812,  1815,  1820,  1824,  1827,  1832,  1836,
    1839,  1844,  1848,  1851,  1856,  1860,  1862,  1865,  1869,  1873,
    1878,  1882,  1885,  1889,  1892,  1898,  1901,  1905,  1910,  1915,
    1925,  1928,  1931,  1941,  1944,  1948,  1953,  1957,  1961,  1969,
    1972,  1976,  1981,  1984,  1987,  1993,  2001,  2012,  2016,  2021,
    2026,  2032,  2037,  2040,  2044,  2048,  2053,  2057,  2060,  2070,
    2077,  2080,  2083,  2087,  2097,  2101,  2111,  2114,  2118,  2123,
    2127,  2131,  2134,  2138,  2141,  2149,  2159,  2168,  2179,  2183,
    2189,  2196,  2198,  2201,  2208,  2216,  2225,  2235,  2237,  2240,
    2242,  2245,  2248,  2252,  2259,  2264,  2272,  2275,  2279,  2286,
    2291,  2299,  2302,  2306,  2313,  2318,  2326,  2329,  2333,  2336,
    2339,  2343,  2346,  2350,  2356,  2361,  2368,  2373,  2376,  2380,
    2383,  2386,  2390,  2396,  2400,  2403,  2409,  2414,  2418,  2423,
    2425,  2428,  2432,  2435,  2452,  2472,  2492,  2513,  2537,  2561,
    2571,  2586,  2603,  2634,  2647,  2661,  2676,  2682,  2689,  2706,
    2719,  2723,  2728,  2742,  2748,  2755,  2762,  2770,  2782,  2791,
    2801,  2812,  2822,  2833,  2845,  2856,  2868,  2881,  2893,  2906,
    2920,  2926,  2933,  2940,  2948,  2956,  2965,  2970,  2974,  2977,
    2981,  2986,  2992,  2999,  3005,  3010,  3016,  3022,  3029,  3037,
    3045,  3054,  3062,  3071,  3076,  3080,  3083,  3089,  3096,  3104,
    3113,  3124,  3131,  3139,  3148,  3158,  3170,  3173,  3179,  3187,
    3191,  3194,  3199,  3202,  3208,  3216,  3219,  3225,  3232,  3240,
    3249,  3260,  3272,  3286,  3301,  3308,  3316,  3325,  3335,  3347,
    3360,  3375,  3391,  3395,  3399,  3403,  3407,  3411,  3414,  3420,
    3425,  3429,  3437,  3444,  3449,  3451,  3454,  3459,  3463,  3469,
    3474,  3478,  3482,  3488,  3493,  3496,  3499,  3502,  3505,  3517,
    3522,  3525,  3528,  3540,  3555,  3568,  3584,  3587,  3595,  3598,
    3604,  3611,  3616,  3623,  3631,  3639,  3646,  3655,  3665,  3673,
    3683,  3694,  3700,  3706,  3714,  3723,  3726,  3731,  3735,  3739,
    3742,  3746,  3749,  3754,  3758,  3761,  3766,  3772,  3775,  3779,
    3784,  3790,  3796,  3802,  3809,  3816,  3823,  3831,  3839,  3848,
    3851,  3855,  3858,  3863,  3870,  3877,  3886,  3889,  3892,  3895,
    3900,  3904,  3910,  3912,  3915,  3918,  3923,  3928,  3933,  3938,
    3943,  3946,  3950,  3953,  3957,  3961,  3965,  3970,  3976,  3983,
    3991,  3996,  4000,  4004,  4008,  4013,  4019,  4022,  4027,  4031,
    4035,  4039,  4043,  4047,  4051,  4055,  4059,  4063,  4067,  4071,
    4076,  4081,  4085,  4089,  4094,  4099,  4104,  4109,  4114,  4119,
    4123,  4127,  4131,  4136,  4140,  4145,  4150,  4155,  4160,  4165,
    4168,  4172,  4176,  4180,  4184,  4188,  4192,  4196,  4199,  4202,
    4206,  4209,  4213,  4216,  4220,  4225,  4229,  4233,  4238,  4243,
    4249,  4255,  4262,  4266,  4271,  4275,  4279,  4283,  4287,  4291,
    4294,  4298,  4302,  4306,  4310,  4314,  4319,  4324,  4329,  4334,
    4339,  4344,  4349,  4353,  4356,  4359,  4363,  4367,  4370,  4374,
    4379,  4385,  4392,  4400,  4403,  4407,  4412,  4416,  4420,  4424,
    4428,  4431,  4434,  4438,  4443,  4449,  4453,  4458,  4462,  4466,
    4470,  4474,  4480,  4486,  4491,  4496,  4503,  4510,  4514,  4518,
    4523,  4527,  4531,  4535,  4539,  4543,  4547,  4550,  4554,  4559,
    4561,  4564,  4568,  4573,  4579,  4581,  4584,  4588,  4591,  4595,
    4600,  4603,  4606,  4610,  4615,  4621,  4626,  4631,  4636,  4639,
    4642,  4645,  4647,  4650,  4655,  4662,  4668,  4670,  4673,  4678,
    4682,  4685,  4689,  4692,  4695,  4697,  4699,  4702,  4706,  4710,
    4714,  4718,  4723,  4728,  4734,  4742,  4748,  4756,  4761,  4767,
    4771,  4776,  4781,  4787,  4794,  4802,  4811,  4815,  4820,  4827,
    4830,  4834,  4839,  4842,  4845,  4856,  4869,  4873,  4878,  4881,
    4886,  4894,  4904,  4914,  4923,  4935,  4948,  4951,  4961,  4972,
    4985,  4995,  5006,  5019,  5029,  5040,  5053,  5063,  5071,  5077,
    5087,  5095,  5101,  5111,  5119,  5125,  5134,  5145,  5154,  5165,
    5173,  5181,  5190,  5199,  5209,  5220,  5231,  5243,  5256,  5281,
    5289,  5298,  5308,  5319,  5331,  5337,  5343,  5346,  5351,  5357,
    5358,  5361,  5362,  5365,  5370,  5377,  5386,  5389,  5393,  5396,
    5399,  5402,  5405,  5408,  5412,  5415,  5419,  5422,  5425,  5427,
    5430,  5434,  5437,  5441,  5444,  5448,  5451,  5454,  5458,  5461,
    5464,  5468,  5473,  5476,  5479,  5482,  5485,  5489,  5494,  5497,
    5500,  5503,  5506,  5509,  5513,  5518,  5521,  5524,  5527,  5530,
    5534,  5539,  5542,  5545,  5548,  5551,  5554,  5558,  5561,  5565,
    5568,  5572,  5574,  5576,  5578,  5580,  5582
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     463,     0,    -1,   464,    91,    -1,   465,    -1,   464,   465,
      -1,   637,    -1,   547,    -1,   603,    -1,   604,    -1,   614,
      -1,   618,    -1,   626,    -1,   646,    -1,   648,    -1,   645,
      -1,   651,    -1,   652,    -1,   655,    -1,   653,    -1,   654,
      -1,   624,    -1,   659,    -1,   656,    -1,   657,    -1,   658,
      -1,   689,    -1,   596,    -1,   595,    -1,   598,    -1,   599,
      -1,   597,    -1,   600,    -1,   601,    -1,   602,    -1,   500,
      -1,   501,    -1,   503,    -1,   502,    -1,   508,    -1,   514,
      -1,   509,    -1,   632,    -1,   519,    -1,   506,    -1,   679,
      -1,   495,    -1,   483,    -1,   484,    -1,   493,    -1,   480,
      -1,   481,    -1,   482,    -1,   608,    -1,   610,    -1,   612,
      -1,   530,    -1,   531,    -1,   496,    -1,   533,    -1,   535,
      -1,   536,    -1,   532,    -1,   497,    -1,   498,    -1,   499,
      -1,   682,    -1,   683,    -1,   476,    -1,   681,    -1,   522,
      -1,   524,    -1,   525,    -1,   526,    -1,   528,    -1,   529,
      -1,   511,    -1,   510,    -1,   555,    -1,   556,    -1,   557,
      -1,   549,    -1,   552,    -1,   558,    -1,   512,    -1,   663,
      -1,   677,    -1,   667,    -1,   673,    -1,   675,    -1,   565,
      -1,   583,    -1,   584,    -1,   670,    -1,   559,    -1,   562,
      -1,   568,    -1,   570,    -1,   572,    -1,   574,    -1,   622,
      -1,   546,    -1,   542,    -1,   543,    -1,   544,    -1,   585,
      -1,   587,    -1,   589,    -1,   590,    -1,   591,    -1,   593,
      -1,   475,    -1,   578,    -1,   579,    -1,   580,    -1,   581,
      -1,   582,    -1,   477,    -1,   478,    -1,   479,    -1,   684,
      -1,   527,    -1,   466,    -1,   467,    -1,   468,    -1,   469,
      -1,   470,    -1,   685,    -1,   686,    -1,   627,    -1,   628,
      -1,   630,    -1,   629,    -1,   631,    -1,   634,    -1,   635,
      -1,   548,    -1,   650,    -1,   538,    -1,   636,    -1,   665,
      -1,   690,    -1,   692,    -1,   693,    -1,   694,    -1,   695,
      -1,   698,    -1,   504,    -1,   224,   217,   700,   700,   217,
      -1,   169,   217,   700,   217,    -1,   134,   217,    -1,   468,
     135,   700,   700,   217,    -1,   468,   135,   700,   700,   700,
     217,    -1,   468,   135,   322,   700,   700,   217,    -1,   283,
     217,    -1,   469,   700,   284,   701,   701,   217,    -1,   373,
     217,   374,   701,   217,    -1,   373,   217,   116,   701,   701,
     700,   217,    -1,   373,   217,   117,   701,   701,   700,   217,
      -1,   373,   217,   115,   701,   701,   700,   700,   217,    -1,
     373,   217,   118,   701,   701,   700,   281,   701,   700,   700,
     700,   700,   217,    -1,   373,   217,   118,   701,   701,   700,
     281,   217,    -1,   373,   700,   700,   217,   374,   701,   217,
      -1,   373,   700,   217,   116,   701,   701,   700,   217,    -1,
     373,   700,   217,   117,   701,   701,   700,   217,    -1,   373,
     700,   217,   115,   701,   701,   700,   700,   217,    -1,   373,
     700,   217,   118,   701,   701,   700,   281,   701,   700,   700,
     700,   700,   217,    -1,   373,   700,   217,   118,   701,   701,
     700,   281,   217,    -1,   373,   217,   473,    -1,   470,   474,
      -1,   470,   541,    -1,   470,   471,    -1,   470,   472,    -1,
     259,   701,   217,    -1,   262,   217,    -1,   262,   701,   701,
     217,    -1,   115,   701,   217,    -1,   473,   701,   700,   217,
      -1,   280,   281,   700,   217,    -1,   280,   281,   700,   700,
     700,   217,    -1,   359,   217,    -1,   475,    39,   316,   217,
      -1,   475,    39,   316,   107,   217,    -1,   475,    40,   316,
     217,    -1,   475,    40,   316,   107,   217,    -1,   475,   360,
     107,   217,    -1,   475,   361,   107,   217,    -1,   475,   362,
     107,   217,    -1,   475,   363,   107,   217,    -1,   475,   364,
     107,   217,    -1,     9,   700,   217,    -1,    79,   217,    -1,
     477,    80,   107,   217,    -1,   477,   403,   700,   217,    -1,
     477,   407,   217,    -1,   477,   406,   217,    -1,   477,   404,
     217,    -1,   477,   405,   217,    -1,   477,    78,   217,    -1,
     477,   348,   217,    -1,   477,   119,   217,    -1,   411,   217,
      -1,   478,   700,   701,   217,    -1,   112,   217,    -1,   479,
      17,   700,   217,    -1,   479,   208,   700,   217,    -1,   479,
     139,   700,   217,    -1,   479,   114,   700,   217,    -1,   203,
     217,    -1,   203,   700,   217,    -1,   480,   701,   701,   217,
      -1,   143,   217,    -1,   143,   700,   217,    -1,   481,   701,
     701,   217,    -1,   183,   217,    -1,   183,   700,   217,    -1,
     482,   700,   701,   217,    -1,   482,   700,   203,   700,   217,
      -1,   482,   700,   143,   700,   217,    -1,    48,   217,    -1,
     483,   485,    -1,   483,   487,    -1,   483,   488,    -1,   483,
     489,    -1,   483,   490,    -1,    45,   217,    -1,   484,   647,
      -1,    44,   700,   217,   486,    -1,    44,   700,   701,   701,
     701,   701,   701,   701,   217,   486,    -1,   700,   700,   701,
     217,    -1,   486,   700,   700,   701,   217,    -1,   175,   700,
     217,    -1,   487,   491,    -1,   176,   700,   217,    -1,   488,
     491,    -1,   177,   700,   217,    -1,   489,   492,    -1,   178,
     700,   217,    -1,   490,   492,    -1,   700,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   217,    -1,   700,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     217,    -1,   700,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   217,    -1,   700,   700,   701,
     701,   217,    -1,   179,   217,    -1,   493,   494,    -1,   700,
     701,   701,   701,   701,   701,   217,    -1,   700,   701,   701,
     701,   701,   701,   701,   701,   217,    -1,   700,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   217,    -1,   700,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     217,    -1,    72,   217,    -1,   495,   700,   700,   701,   217,
      -1,   495,   700,   700,   700,   701,   217,    -1,   130,   217,
      -1,   496,   701,   701,   701,   217,    -1,   279,   217,    -1,
     497,   107,   107,   107,   217,    -1,   497,   107,   107,   217,
      -1,   497,   107,   700,   217,    -1,   182,   107,   217,    -1,
     349,   107,   217,    -1,   293,   217,   605,    -1,     3,   217,
     605,    -1,   351,   217,   605,    -1,   350,   217,   605,    -1,
     454,   217,    -1,   504,   505,   217,    -1,   456,   457,   700,
      -1,   505,   310,    -1,   231,   217,    -1,   232,   217,    -1,
     231,   700,   217,    -1,   232,   700,   217,    -1,   231,   172,
     217,    -1,   232,   172,   217,    -1,   231,   172,   700,   217,
      -1,   232,   172,   700,   217,    -1,   506,   507,   217,    -1,
     308,   107,   700,    -1,   308,   700,   700,   107,   700,    -1,
     308,   107,   700,   700,    -1,   308,   107,   700,   135,   700,
      -1,   308,   700,   700,   107,   700,   700,    -1,   308,   700,
     700,   107,   700,   135,   700,    -1,   455,   107,   700,    -1,
     455,   700,   700,   107,   700,    -1,   345,   107,   700,    -1,
     345,   700,   700,   107,   700,    -1,   345,   107,   700,   700,
      -1,   345,   107,   700,   135,   700,    -1,   345,   700,   700,
     107,   700,   700,    -1,   345,   700,   700,   107,   700,   135,
     700,    -1,   507,   223,   700,    -1,   507,   310,    -1,   507,
     701,   701,    -1,   507,    15,    -1,   507,    61,    -1,   507,
      61,   700,    -1,   507,    12,    -1,   507,   287,    -1,   507,
     381,    -1,   507,   288,   316,    -1,   507,   214,    -1,   507,
     214,   700,    -1,   507,   233,    -1,   507,   190,    -1,   461,
     107,    -1,   459,   107,    -1,   460,   107,    -1,   513,    -1,
      88,   217,    -1,    88,   700,   217,    -1,   458,   316,   217,
      -1,   215,   700,   217,    -1,   309,   217,    -1,   309,   700,
     701,   701,   217,    -1,   225,   700,   217,    -1,   328,   701,
     217,    -1,   330,   701,   217,    -1,    94,   217,    -1,   295,
     701,   217,    -1,    20,   217,    -1,    20,   107,   217,    -1,
      20,   107,   700,   217,    -1,    20,   701,   701,   217,    -1,
      20,   701,   701,   700,   217,    -1,   113,   316,   217,    -1,
     171,   700,   217,    -1,   184,   217,    -1,   508,   191,   700,
     217,    -1,   346,   217,    -1,    31,   700,   217,    -1,   413,
     217,    -1,   414,   701,   217,    -1,   188,   217,    -1,    49,
     217,   701,   700,   217,    -1,    49,   217,    -1,   332,   217,
      -1,    77,   217,    -1,   513,   516,    -1,   513,   537,    -1,
     513,   541,    -1,   513,   194,   217,    -1,   513,   317,   700,
     217,    -1,   513,   317,   316,   217,    -1,   513,   317,   700,
     701,   701,   700,   700,   217,    -1,   513,   155,   316,   217,
      -1,   513,   358,   316,   217,    -1,   513,   228,   217,    -1,
     513,   325,   316,   217,    -1,   513,    42,   217,    -1,   513,
      42,   701,   701,   217,    -1,    68,   217,   515,   217,    -1,
     701,   701,   701,   701,   701,    -1,   216,   217,    -1,   208,
     517,   217,    -1,    17,   517,   217,    -1,   139,   518,   217,
      -1,   701,   701,    -1,   701,   701,   701,   701,    -1,   701,
      -1,   701,   701,   701,    -1,   701,    -1,   234,   217,    -1,
     519,   520,    -1,   519,   521,    -1,   519,   194,   217,    -1,
     208,   701,   701,   700,   217,    -1,   208,   701,   701,   700,
     701,   217,    -1,   139,   701,   701,   701,   700,   217,    -1,
       4,     6,   217,    -1,     4,   217,     6,   217,    -1,     4,
     217,     6,   701,   701,   217,    -1,     4,   217,     6,   701,
     217,    -1,   522,    46,   217,    -1,   522,   363,   107,   217,
      -1,   522,   523,   217,    -1,    32,    -1,   523,   700,    -1,
       5,   217,     6,   701,   701,   217,    -1,   336,   217,    -1,
     335,   217,    -1,   366,   217,    -1,   150,   217,    -1,   415,
     217,    -1,   334,   217,   701,   217,    -1,   131,   217,   701,
     701,   217,    -1,   131,   217,   701,   217,    -1,   131,   217,
      -1,   131,   217,   701,   701,   700,   217,    -1,   210,   700,
     217,    -1,   210,   217,    -1,   271,   700,   217,    -1,   271,
     217,    -1,   271,   217,   534,   217,    -1,   700,    -1,   534,
     700,    -1,   151,   217,    -1,   416,   217,    -1,   327,   701,
     701,   701,   217,    -1,   236,   217,   700,   700,   539,    -1,
     236,   217,   700,   700,   700,   539,    -1,   240,   217,   700,
     700,   700,   700,   539,    -1,   217,   540,   539,    -1,   217,
      -1,   239,    -1,   241,   700,    -1,   242,    -1,   243,    -1,
     244,    -1,   245,   701,    -1,   246,    -1,   247,    -1,   247,
     701,    -1,   248,    -1,   270,   701,   701,   217,    -1,   270,
     701,   701,   211,   217,    -1,   209,   217,   606,    -1,   138,
     217,   623,    -1,   165,   217,   701,   701,   701,   217,    -1,
     165,   217,   701,   701,   217,    -1,   166,   217,   700,   701,
     701,   217,    -1,   166,   700,   217,    -1,   544,   545,    -1,
     700,   701,   701,   701,   217,    -1,   167,   217,   701,   701,
     701,   217,    -1,    73,   217,    -1,   547,   641,    -1,   547,
     700,   347,   700,   700,   701,   217,    -1,   547,   700,   347,
     700,   319,   700,   700,   701,   217,    -1,   547,   322,   641,
      -1,    64,   217,    -1,   548,   700,    66,   700,   700,   217,
      -1,   417,   217,   550,    -1,   551,    -1,   550,   551,    -1,
     700,   701,   217,    -1,   700,   217,    -1,   419,   217,   553,
      -1,   554,    -1,   553,   554,    -1,   700,   700,   564,   217,
      -1,   326,   217,    -1,   555,   700,   701,   217,    -1,   555,
     700,   347,   700,   701,   217,    -1,   555,   700,   347,   700,
     319,   700,   701,   217,    -1,   555,   322,   700,   701,   217,
      -1,   108,   217,    -1,   108,   700,   217,    -1,   556,   700,
     701,   217,    -1,   556,   322,   700,   701,   217,    -1,    47,
     217,    -1,    47,   700,   217,    -1,   557,   700,   701,   701,
     701,   217,    -1,   269,   217,    -1,   269,   700,   217,    -1,
     558,   700,   701,   701,   701,   217,    -1,   142,   217,   560,
      -1,   561,    -1,   560,   561,    -1,   700,   700,   701,   217,
      -1,   700,   700,   700,   701,   217,    -1,   700,   700,   700,
     700,   701,   217,    -1,    92,   217,   563,    -1,   562,   563,
      -1,   700,   700,   564,   217,    -1,   700,    -1,   564,   700,
      -1,   313,   217,   566,    -1,   567,    -1,   566,   567,    -1,
     700,   700,   701,   217,    -1,   700,   700,   700,   701,   217,
      -1,   700,   700,   700,   700,   701,   217,    -1,    90,   217,
     569,    -1,   568,   569,    -1,   700,   700,   564,   217,    -1,
     145,   217,   571,    -1,   570,   571,    -1,   700,   700,   564,
     217,    -1,   149,   217,   573,    -1,   572,   573,    -1,   700,
     700,   564,   217,    -1,   148,   217,   575,    -1,   574,   575,
      -1,   700,   700,   564,   217,    -1,   700,   701,   217,    -1,
     576,    -1,   577,   576,    -1,    21,   217,   577,    -1,    22,
     217,   577,    -1,    22,   700,   217,   577,    -1,    16,   701,
     217,    -1,   580,   563,    -1,    18,   701,   217,    -1,   581,
     571,    -1,    19,   701,   701,   701,   217,    -1,   582,   575,
      -1,   104,   700,   217,    -1,   104,   700,   700,   217,    -1,
     105,   700,   217,   671,    -1,    23,   217,   701,   701,   217,
     701,   701,   701,   217,    -1,   585,   586,    -1,   700,   217,
      -1,    24,   217,   701,   701,   217,   701,   701,   701,   217,
      -1,   587,   588,    -1,   700,   700,   217,    -1,   700,   700,
     700,   217,    -1,    25,   700,   217,    -1,    26,   700,   217,
      -1,    27,   217,   700,   217,   701,   701,   217,    -1,   591,
     592,    -1,   700,   700,   217,    -1,   700,   700,   700,   217,
      -1,    28,   217,    -1,   593,   594,    -1,   700,   701,   701,
     700,   217,    -1,   700,   701,   701,   700,   701,   701,   217,
      -1,   700,   701,   701,   700,   701,   701,   701,   701,   701,
     217,    -1,   273,   107,   217,    -1,   273,   107,   700,   217,
      -1,   273,   107,   107,   217,    -1,   273,   107,   107,   700,
     217,    -1,   595,   453,   316,   217,    -1,   157,   217,    -1,
     157,   358,   217,    -1,   157,   217,   605,    -1,   596,   194,
     217,   606,    -1,   158,   701,   217,    -1,   158,   217,    -1,
     597,   700,   701,   701,   701,   701,   701,   701,   217,    -1,
     597,   700,   701,   701,   701,   217,    -1,   128,   217,    -1,
      36,   217,    -1,   237,   700,   217,    -1,   598,   700,   701,
     701,   701,   701,   701,   701,   217,    -1,   238,   700,   217,
      -1,   599,   700,   701,   701,   701,   701,   701,   701,   217,
      -1,   163,   217,    -1,   163,   217,   605,    -1,   600,   194,
     217,   606,    -1,   161,   217,   607,    -1,    93,   217,   607,
      -1,   109,   217,    -1,   109,   700,   217,    -1,   603,   641,
      -1,   603,   700,   347,   700,   700,   701,   217,    -1,   603,
     700,   347,   700,   319,   700,   700,   701,   217,    -1,   603,
     700,   347,   700,   700,   701,   211,   217,    -1,   603,   700,
     347,   700,   319,   700,   700,   701,   211,   217,    -1,   603,
     322,   641,    -1,   109,   217,   194,   217,   606,    -1,   109,
     700,   217,   194,   217,   606,    -1,   641,    -1,   605,   641,
      -1,   700,   347,   700,   700,   701,   217,    -1,   605,   700,
     347,   700,   700,   701,   217,    -1,   700,   347,   700,   319,
     700,   700,   701,   217,    -1,   605,   700,   347,   700,   319,
     700,   700,   701,   217,    -1,   642,    -1,   606,   642,    -1,
     643,    -1,   607,   643,    -1,   357,   217,    -1,   357,   217,
     609,    -1,    53,   700,   217,   701,   701,   217,    -1,   609,
     701,   701,   217,    -1,   609,    53,   700,   217,   701,   701,
     217,    -1,   337,   217,    -1,   337,   217,   611,    -1,    53,
     700,   217,   701,   701,   217,    -1,   611,   701,   701,   217,
      -1,   611,    53,   700,   217,   701,   701,   217,    -1,   292,
     217,    -1,   292,   217,   613,    -1,    53,   700,   217,   701,
     701,   217,    -1,   613,   701,   701,   217,    -1,   613,    53,
     700,   217,   701,   701,   217,    -1,   181,   217,    -1,   181,
     217,   615,    -1,   616,   617,    -1,   615,   617,    -1,   615,
     616,   617,    -1,   700,   217,    -1,   700,   701,   217,    -1,
     700,   701,   366,   700,   217,    -1,   700,   701,   666,   217,
      -1,   700,   701,   366,   700,   666,   217,    -1,   700,   700,
     701,   217,    -1,   152,   217,    -1,   152,   217,   619,    -1,
     620,   621,    -1,   619,   621,    -1,   619,   620,   621,    -1,
     700,    56,   701,   701,   217,    -1,   700,   701,   217,    -1,
     700,   217,    -1,   700,   700,   701,   701,   217,    -1,   700,
     700,   701,   217,    -1,   141,   217,   623,    -1,   141,   700,
     217,   623,    -1,   644,    -1,   623,   644,    -1,   189,   217,
     625,    -1,   624,   625,    -1,   700,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     217,    -1,   700,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   270,   701,   701,
     217,    -1,   700,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   291,   701,   701,
     217,    -1,   700,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   217,    -1,   700,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   270,   701,   701,   217,    -1,   700,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   291,   701,   701,   217,
      -1,   700,   701,   701,   701,   701,   701,   701,   701,   217,
      -1,   700,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   217,    -1,   700,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,    69,
     701,   217,    -1,   700,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   700,   701,   701,   700,   700,   701,
     701,   701,   217,    -1,   700,     8,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   217,    -1,   700,     8,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   217,
      -1,   700,     8,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   217,    -1,   700,     8,   701,   701,
     217,    -1,   700,     8,   701,   701,   701,   217,    -1,   700,
      98,   700,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   700,   700,   700,   217,    -1,   700,   344,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   217,    -1,   700,
      62,   217,    -1,   700,    62,   666,   217,    -1,   700,   188,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     217,    -1,   700,    62,   188,   701,   217,    -1,   700,    62,
     666,   188,   701,   217,    -1,   700,    62,   188,   701,   701,
     217,    -1,   700,    62,   666,   188,   701,   701,   217,    -1,
     700,    62,   700,   701,   701,   701,   701,   701,   701,   700,
     217,    -1,   700,    62,    97,   701,   701,   701,   701,   217,
      -1,   700,    62,    97,   701,   701,   701,   701,   701,   217,
      -1,   700,    62,    97,   701,   701,   701,   701,   701,   701,
     217,    -1,   700,    62,   666,    97,   701,   701,   701,   701,
     217,    -1,   700,    62,   666,    97,   701,   701,   701,   701,
     701,   217,    -1,   700,    62,   666,    97,   701,   701,   701,
     701,   701,   701,   217,    -1,   700,    62,    97,   701,   701,
     701,   701,   323,   701,   217,    -1,   700,    62,    97,   701,
     701,   701,   701,   701,   323,   701,   217,    -1,   700,    62,
      97,   701,   701,   701,   701,   701,   701,   323,   701,   217,
      -1,   700,    62,   666,    97,   701,   701,   701,   701,   323,
     701,   217,    -1,   700,    62,   666,    97,   701,   701,   701,
     701,   701,   323,   701,   217,    -1,   700,    62,   666,    97,
     701,   701,   701,   701,   701,   701,   323,   701,   217,    -1,
     700,    62,   323,   701,   217,    -1,   700,    62,   666,   323,
     701,   217,    -1,   700,    62,   323,   701,   701,   217,    -1,
     700,    62,   666,   323,   701,   701,   217,    -1,   700,    62,
     323,   701,   701,   701,   217,    -1,   700,    62,   666,   323,
     701,   701,   701,   217,    -1,   700,   323,   701,   217,    -1,
     333,   217,   639,    -1,   626,   639,    -1,   395,   700,   217,
      -1,   395,   700,   286,   217,    -1,   395,   700,   321,   701,
     217,    -1,   395,   700,   321,   701,   286,   217,    -1,   627,
     700,   700,   640,   217,    -1,   396,   700,   700,   217,    -1,
     396,   700,   700,   400,   217,    -1,   396,   700,   700,   401,
     217,    -1,   396,   700,   700,   399,   701,   217,    -1,   396,
     700,   700,   399,   701,   701,   217,    -1,   396,   700,   700,
     400,   399,   701,   217,    -1,   396,   700,   700,   400,   399,
     701,   701,   217,    -1,   396,   700,   700,   401,   399,   701,
     217,    -1,   396,   700,   700,   401,   399,   701,   701,   217,
      -1,   402,   700,   700,   217,    -1,   402,   700,   217,    -1,
     340,   217,    -1,   630,   700,   700,   700,   217,    -1,   630,
     700,   700,   700,   700,   217,    -1,   630,   700,   700,   700,
     700,   701,   217,    -1,   630,   700,   700,   700,   700,   701,
     701,   217,    -1,   630,   700,   700,   700,   700,   701,   701,
     700,   701,   217,    -1,   630,   700,   700,   700,   666,   217,
      -1,   630,   700,   700,   700,   700,   666,   217,    -1,   630,
     700,   700,   700,   700,   701,   666,   217,    -1,   630,   700,
     700,   700,   700,   701,   701,   666,   217,    -1,   630,   700,
     700,   700,   700,   701,   701,   700,   701,   666,   217,    -1,
     120,   217,    -1,   631,   700,   700,   700,   217,    -1,   631,
     700,   700,   700,   701,   701,   217,    -1,   418,   217,   633,
      -1,   632,   633,    -1,   700,   700,   564,   217,    -1,   420,
     217,    -1,   634,   700,   700,   700,   217,    -1,   634,   700,
     700,   700,   701,   701,   217,    -1,    58,   217,    -1,   635,
     700,   700,   700,   217,    -1,   635,   700,   700,   700,   700,
     217,    -1,   635,   700,   700,   700,   700,   701,   217,    -1,
     635,   700,   700,   700,   700,   701,   701,   217,    -1,   635,
     700,   700,   700,   700,   701,   701,   700,   701,   217,    -1,
     635,   700,   700,   700,   700,   701,   701,   700,   701,   701,
     217,    -1,   635,   700,   700,   700,   700,   701,   701,   700,
     701,   701,   701,   701,   217,    -1,   635,   700,   700,   700,
     700,   701,   701,   700,   701,   701,   701,   701,   701,   217,
      -1,   635,   700,   700,   700,   666,   217,    -1,   635,   700,
     700,   700,   700,   666,   217,    -1,   635,   700,   700,   700,
     700,   701,   666,   217,    -1,   635,   700,   700,   700,   700,
     701,   701,   666,   217,    -1,   635,   700,   700,   700,   700,
     701,   701,   700,   701,   666,   217,    -1,   635,   700,   700,
     700,   700,   701,   701,   700,   701,   701,   666,   217,    -1,
     635,   700,   700,   700,   700,   701,   701,   700,   701,   701,
     701,   701,   666,   217,    -1,   635,   700,   700,   700,   700,
     701,   701,   700,   701,   701,   701,   701,   701,   666,   217,
      -1,    30,   700,   217,    -1,   125,   700,   217,    -1,   397,
     701,   217,    -1,   398,   700,   217,    -1,   223,   217,   638,
      -1,   637,   638,    -1,   700,   701,   701,   701,   217,    -1,
     700,   701,   701,   217,    -1,   700,   701,   217,    -1,   700,
     701,   701,   701,   700,   700,   217,    -1,   700,   701,   701,
     701,   700,   217,    -1,   700,   700,   640,   217,    -1,   700,
      -1,   640,   700,    -1,   700,   700,   701,   217,    -1,   700,
     700,   217,    -1,   700,   700,   701,   211,   217,    -1,   700,
     700,   211,   217,    -1,   700,   701,   217,    -1,   700,   701,
     217,    -1,   700,   700,   701,   701,   217,    -1,   700,   700,
     701,   217,    -1,    65,   217,    -1,   645,   647,    -1,    89,
     217,    -1,   646,   647,    -1,   700,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   217,    -1,   700,   343,   700,
     217,    -1,   229,   217,    -1,   648,   649,    -1,   700,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   217,    -1,
     700,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   217,    -1,   700,   127,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   217,    -1,   700,   127,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   217,    -1,    35,   217,    -1,   650,   700,   700,
     701,   701,   701,   217,    -1,    11,   217,    -1,   651,   700,
     342,   701,   217,    -1,   651,   700,   342,   701,    95,   217,
      -1,   651,   700,   700,   217,    -1,   651,   700,   700,   342,
     701,   217,    -1,   651,   700,   700,   342,   275,   701,   217,
      -1,   651,   700,   700,   342,   701,    95,   217,    -1,   651,
     700,   700,   700,   700,   217,    -1,   651,   700,   700,   700,
     700,   342,   701,   217,    -1,   651,   700,   700,   700,   700,
     342,   701,    95,   217,    -1,   651,   700,   700,   700,   341,
     701,   217,    -1,   651,   700,   700,   700,   341,   701,   342,
     701,   217,    -1,   651,   700,   700,   700,   341,   701,   342,
     701,    95,   217,    -1,   651,   700,   700,   156,   217,    -1,
     651,   700,   700,   700,   217,    -1,   651,   700,   700,   700,
     700,   700,   217,    -1,   651,   700,   700,   700,   700,   341,
     701,   217,    -1,    87,   217,    -1,   652,   700,   701,   217,
      -1,   652,   275,   217,    -1,   652,    95,   217,    -1,   447,
     217,    -1,   653,   701,   217,    -1,   448,   217,    -1,   654,
     700,   700,   217,    -1,   654,   701,   217,    -1,   446,   217,
      -1,   655,   700,   700,   217,    -1,   655,   700,   700,   700,
     217,    -1,   252,   217,    -1,   252,   700,   217,    -1,   656,
     700,   701,   217,    -1,   656,   700,   700,   701,   217,    -1,
     656,   322,   700,   701,   217,    -1,   656,   700,   701,   316,
     217,    -1,   656,   700,   700,   701,   316,   217,    -1,   656,
     322,   700,   701,   316,   217,    -1,   656,   700,    99,   700,
     701,   217,    -1,   656,   700,   700,    99,   700,   701,   217,
      -1,   656,   700,    99,   700,   701,   316,   217,    -1,   656,
     700,   700,    99,   700,   701,   316,   217,    -1,   187,   217,
      -1,   187,   700,   217,    -1,   251,   217,    -1,   658,   700,
     701,   217,    -1,   658,   700,   347,   700,   701,   217,    -1,
     658,   700,   701,   701,   701,   217,    -1,   658,   700,   347,
     700,   701,   701,   701,   217,    -1,   307,   217,    -1,   659,
     662,    -1,   659,   661,    -1,   659,    63,   660,   217,    -1,
     659,   268,   217,    -1,   659,   268,   701,   701,   217,    -1,
     700,    -1,   660,   700,    -1,   162,   217,    -1,   661,   249,
     700,   217,    -1,   661,   191,   700,   217,    -1,   661,   331,
     701,   217,    -1,   661,   192,   700,   217,    -1,   661,   318,
     700,   217,    -1,    74,   217,    -1,    74,   700,   217,    -1,
     294,   217,    -1,   294,   255,   217,    -1,   294,   353,   217,
      -1,   162,   700,   217,    -1,   162,   700,   701,   217,    -1,
     162,   700,   701,   700,   217,    -1,   162,   700,   701,   700,
     700,   217,    -1,   162,   700,   701,   700,   700,   700,   217,
      -1,   101,   700,   701,   217,    -1,   101,   700,   217,    -1,
     101,    76,   217,    -1,   101,   375,   217,    -1,   101,   700,
     102,   217,    -1,   101,   700,   701,   700,   217,    -1,   101,
     217,    -1,    34,   217,   294,   217,    -1,   305,   700,   217,
      -1,   306,   700,   217,    -1,   296,   701,   217,    -1,   298,
     700,   217,    -1,   299,   700,   217,    -1,   297,   700,   217,
      -1,   300,   701,   217,    -1,   301,   700,   217,    -1,   303,
     316,   217,    -1,   302,   700,   217,    -1,   304,   700,   217,
      -1,   206,   700,   700,   217,    -1,   207,   700,   701,   217,
      -1,   136,   701,   217,    -1,   137,   316,   217,    -1,   662,
     191,   700,   217,    -1,    83,   700,   700,   217,    -1,    82,
     700,   701,   217,    -1,   662,   249,   103,   217,    -1,   662,
     249,   187,   217,    -1,   662,   249,   700,   217,    -1,   256,
     700,   217,    -1,   256,   257,   217,    -1,   329,   701,   217,
      -1,   329,   701,   701,   217,    -1,   314,   701,   217,    -1,
     314,   701,   701,   217,    -1,   265,   701,   701,   217,    -1,
     267,   700,   700,   217,    -1,   266,   701,   701,   217,    -1,
     662,   192,   700,   217,    -1,   222,   217,    -1,   254,   700,
     217,    -1,   289,   700,   217,    -1,   289,   290,   217,    -1,
     199,   700,   217,    -1,   199,   290,   217,    -1,   121,   700,
     217,    -1,   121,   290,   217,    -1,   200,   217,    -1,   122,
     217,    -1,   124,   700,   217,    -1,   123,   217,    -1,    57,
     701,   217,    -1,   356,   217,    -1,    51,    52,   217,    -1,
      51,    52,   700,   217,    -1,    51,    14,   217,    -1,    13,
      14,   217,    -1,    13,    14,   272,   217,    -1,    13,   370,
     700,   217,    -1,    13,   370,   371,   700,   217,    -1,    13,
     370,   700,   376,   217,    -1,    13,   370,   371,   700,   376,
     217,    -1,   372,   701,   217,    -1,   372,   701,   701,   217,
      -1,   129,   701,   217,    -1,   133,   701,   217,    -1,    54,
     701,   217,    -1,   282,   316,   217,    -1,   352,   700,   217,
      -1,   253,   217,    -1,   253,   107,   217,    -1,   185,   294,
     217,    -1,    43,   294,   217,    -1,    29,   294,   217,    -1,
      55,   294,   217,    -1,   185,   294,   255,   217,    -1,   185,
     294,   315,   217,    -1,    43,   294,   255,   217,    -1,    29,
     294,   255,   217,    -1,    43,   294,   315,   217,    -1,    55,
     294,   255,   217,    -1,    55,   294,   315,   217,    -1,   355,
     700,   217,    -1,   132,   217,    -1,   412,   217,    -1,   412,
     316,   217,    -1,   258,   700,   217,    -1,   221,   217,    -1,
     221,   700,   217,    -1,   140,   700,   701,   217,    -1,   140,
     700,   701,   700,   217,    -1,   140,   700,   701,   700,   701,
     217,    -1,   140,   700,   701,   700,   701,   700,   217,    -1,
     140,   217,    -1,   227,   700,   217,    -1,   227,   700,   701,
     217,    -1,   338,   701,   217,    -1,   312,   700,   217,    -1,
     173,   700,   217,    -1,   173,   168,   217,    -1,   160,   217,
      -1,   311,   217,    -1,   339,   700,   217,    -1,   339,   700,
     701,   217,    -1,   339,   700,   701,   701,   217,    -1,   369,
     162,   217,    -1,   369,   162,   153,   217,    -1,   197,   700,
     217,    -1,   197,   198,   217,    -1,   195,   700,   217,    -1,
     195,   196,   217,    -1,   195,   700,   202,   700,   217,    -1,
     195,   196,   202,   700,   217,    -1,   195,   700,   700,   217,
      -1,   195,   196,   201,   217,    -1,   195,   700,   700,   202,
     700,   217,    -1,   195,   196,   201,   202,   700,   217,    -1,
     204,   700,   217,    -1,   260,   701,   217,    -1,   186,   700,
     701,   217,    -1,    38,   316,   217,    -1,    81,   316,   217,
      -1,    59,   316,   217,    -1,    60,   316,   217,    -1,   205,
     700,   217,    -1,   100,   217,   664,    -1,   701,   217,    -1,
      84,   666,   217,    -1,    84,   217,   666,   217,    -1,    74,
      -1,    74,   701,    -1,    74,   701,   701,    -1,    74,   701,
     701,   701,    -1,    74,   701,   701,   701,   701,    -1,    85,
      -1,    86,   701,    -1,    85,    86,   701,    -1,    33,   701,
      -1,   666,   154,   700,    -1,   666,   154,   700,   701,    -1,
     144,   217,    -1,   111,   217,    -1,    37,   700,   217,    -1,
      37,   700,   701,   217,    -1,    37,   700,   701,   701,   217,
      -1,   263,   700,   217,   668,    -1,   264,   700,   217,   668,
      -1,   164,   700,   217,   668,    -1,   667,   541,    -1,   667,
     474,    -1,   667,   471,    -1,   669,    -1,   668,   669,    -1,
     701,   701,   701,   217,    -1,   174,   217,   701,   701,   701,
     217,    -1,   670,   701,   701,   701,   217,    -1,   672,    -1,
     671,   672,    -1,   701,   701,   701,   217,    -1,   146,   217,
     674,    -1,   701,   217,    -1,   147,   217,   676,    -1,   701,
     217,    -1,    71,   217,    -1,     7,    -1,   700,    -1,   219,
     217,    -1,   679,   218,   217,    -1,   679,    10,   217,    -1,
     679,   220,   217,    -1,   679,   381,   217,    -1,   679,   191,
     700,   217,    -1,   679,   226,   701,   217,    -1,   679,   226,
     701,   701,   217,    -1,   679,   226,   701,   701,   701,   701,
     217,    -1,   679,    75,   701,   701,   217,    -1,   679,    75,
     701,   701,   700,   700,   217,    -1,   679,   106,   700,   217,
      -1,   679,   106,   700,   700,   217,    -1,   679,   353,   217,
      -1,   679,   180,   701,   217,    -1,   679,   186,   678,   217,
      -1,   679,   186,   678,   700,   217,    -1,   679,   186,   678,
     700,   701,   217,    -1,   679,   186,   678,   700,   701,   701,
     217,    -1,   679,   186,   678,   700,   701,   701,   700,   217,
      -1,   679,   126,   217,    -1,   679,   126,   701,   217,    -1,
     679,    86,   700,   701,   701,   217,    -1,   679,   680,    -1,
     274,   700,   217,    -1,   274,   700,   700,   217,    -1,   278,
     217,    -1,    50,   217,    -1,    50,   217,   107,   217,   700,
     217,   107,   217,   107,   217,    -1,    50,   217,   107,   217,
     700,   217,   107,   217,   107,   217,   107,   217,    -1,   230,
     107,   217,    -1,   683,   217,   107,   217,    -1,   365,   217,
      -1,   365,   366,   700,   217,    -1,   684,   700,   700,   701,
     701,   701,   217,    -1,   684,   700,   700,   701,   701,   701,
     368,   701,   217,    -1,   684,   700,   700,   701,   701,   701,
     366,   700,   217,    -1,   684,   700,   700,   701,   701,   701,
     666,   217,    -1,   684,   700,   700,   701,   701,   701,   368,
     701,   366,   700,   217,    -1,   684,   700,   700,   701,   701,
     701,   368,   701,   366,   700,   666,   217,    -1,   377,   217,
      -1,   685,   700,   379,   701,   701,   701,   701,   701,   217,
      -1,   685,   700,   379,   701,   701,   701,   701,   701,   701,
     217,    -1,   685,   700,   379,   701,   701,   701,   701,   701,
     701,   701,   701,   217,    -1,   685,   700,   380,   701,   701,
     701,   701,   701,   217,    -1,   685,   700,   380,   701,   701,
     701,   701,   701,   701,   217,    -1,   685,   700,   380,   701,
     701,   701,   701,   701,   701,   701,   701,   217,    -1,   685,
     700,   393,   701,   701,   701,   701,   701,   217,    -1,   685,
     700,   393,   701,   701,   701,   701,   701,   701,   217,    -1,
     685,   700,   393,   701,   701,   701,   701,   701,   701,   701,
     701,   217,    -1,   685,   700,   381,   701,   701,   701,   701,
     701,   217,    -1,   685,   700,   381,   701,   701,   701,   217,
      -1,   685,   700,   381,   701,   217,    -1,   685,   700,   382,
     701,   701,   701,   701,   701,   217,    -1,   685,   700,   382,
     701,   701,   701,   217,    -1,   685,   700,   382,   701,   217,
      -1,   685,   700,   392,   701,   701,   701,   701,   701,   217,
      -1,   685,   700,   392,   701,   701,   701,   217,    -1,   685,
     700,   392,   701,   217,    -1,   685,   700,   383,   701,   701,
     701,   701,   217,    -1,   685,   700,   383,   701,   701,   701,
     701,   701,   701,   217,    -1,   685,   700,   394,   701,   701,
     701,   701,   217,    -1,   685,   700,   394,   701,   701,   701,
     701,   701,   701,   217,    -1,   685,   700,   386,   701,   701,
     701,   217,    -1,   685,   700,   387,   701,   701,   701,   217,
      -1,   685,   700,   387,   701,   701,   701,   701,   217,    -1,
     685,   700,   391,   701,   701,   701,   701,   217,    -1,   685,
     700,   391,   701,   701,   701,   701,   701,   217,    -1,   685,
     700,   388,   701,   701,   701,   701,   701,   701,   217,    -1,
     685,   700,   389,   701,   701,   701,   701,   701,   701,   217,
      -1,   685,   700,   389,   701,   701,   701,   701,   701,   701,
     701,   217,    -1,   685,   700,   389,   701,   701,   701,   701,
     701,   701,   701,   701,   217,    -1,   685,   700,   385,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   217,
      -1,   685,   700,   385,   701,   701,   701,   217,    -1,   685,
     700,   385,   701,   701,   701,   701,   217,    -1,   685,   700,
     385,   701,   701,   701,   701,   701,   217,    -1,   685,   700,
     385,   701,   701,   701,   701,   701,   701,   217,    -1,   685,
     700,   385,   701,   701,   701,   701,   701,   701,   701,   217,
      -1,   685,   384,   107,   107,   217,    -1,   685,   700,   107,
     687,   217,    -1,   378,   217,    -1,   686,   700,   700,   217,
      -1,   686,   700,   700,   700,   217,    -1,    -1,   687,   701,
      -1,    -1,   688,   107,    -1,   276,   217,   277,   217,    -1,
     276,   217,   277,   217,   277,   217,    -1,   276,   217,   277,
     217,   277,   217,   277,   217,    -1,   426,   217,    -1,   690,
     691,   217,    -1,   421,   688,    -1,   431,   700,    -1,   434,
     700,    -1,   434,   316,    -1,   436,   687,    -1,   436,   316,
     687,    -1,   405,   700,    -1,   405,   700,   700,    -1,   437,
     688,    -1,    41,   700,    -1,   697,    -1,   450,   217,    -1,
     692,   696,   217,    -1,   451,   217,    -1,   693,   696,   217,
      -1,   429,   217,    -1,   694,   696,   217,    -1,   694,   697,
      -1,   430,   217,    -1,   695,   696,   217,    -1,   422,   107,
      -1,   423,   107,    -1,   423,   107,   107,    -1,   423,   107,
     107,   107,    -1,   433,   701,    -1,   405,   700,    -1,   424,
     700,    -1,   431,   700,    -1,   431,   700,   700,    -1,   431,
     700,   700,   700,    -1,   453,   316,    -1,   212,   316,    -1,
     275,   316,    -1,    95,   316,    -1,   444,   107,    -1,   444,
     107,   107,    -1,   444,   107,   107,   107,    -1,   289,   316,
      -1,   408,   409,    -1,   410,   701,    -1,   449,   107,    -1,
     449,   107,   700,    -1,   449,   107,   700,   700,    -1,   442,
     316,    -1,   443,   316,    -1,   445,   316,    -1,   452,   316,
      -1,   441,   217,    -1,   697,   515,   217,    -1,   427,   217,
      -1,   698,   699,   217,    -1,   428,   107,    -1,   428,   107,
     700,    -1,   159,    -1,   213,    -1,   159,    -1,    70,    -1,
     170,    -1,    96,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   168,   168,   174,   175,   178,   179,   181,   183,   184,
     185,   186,   187,   188,   189,   190,   192,   193,   194,   195,
     196,   197,   198,   199,   201,   203,   204,   205,   206,   207,
     208,   210,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     230,   232,   233,   234,   235,   236,   237,   238,   239,   240,
     241,   242,   243,   244,   245,   246,   247,   248,   249,   250,
     251,   252,   253,   254,   255,   256,   257,   258,   260,   262,
     264,   266,   268,   270,   271,   272,   273,   274,   275,   276,
     277,   278,   279,   280,   281,   282,   283,   284,   285,   286,
     288,   289,   291,   292,   293,   294,   295,   296,   297,   298,
     299,   300,   302,   304,   305,   306,   307,   308,   309,   310,
     312,   313,   314,   315,   316,   318,   319,   320,   321,   323,
     325,   327,   329,   331,   333,   335,   337,   338,   339,   340,
     341,   342,   343,   344,   345,   346,   347,   350,   357,   363,
     364,   369,   381,   387,   388,   392,   396,   398,   400,   402,
     422,   448,   452,   454,   456,   458,   478,   504,   505,   506,
     507,   508,   511,   515,   517,   522,   525,   529,   568,   617,
     618,   620,   628,   630,   638,   640,   642,   644,   646,   650,
     658,   660,   662,   664,   666,   668,   670,   672,   674,   676,
     680,   682,   686,   688,   690,   692,   694,   697,   699,   701,
     705,   707,   709,   713,   715,   717,   719,   721,   725,   726,
     727,   728,   729,   730,   733,   734,   738,   740,   750,   752,
     756,   758,   762,   764,   768,   770,   774,   776,   780,   787,
     794,   803,   807,   808,   812,   818,   823,   828,   834,   835,
     837,   841,   842,   846,   847,   851,   854,   859,   864,   868,
     873,   879,   886,   893,   898,   902,   905,   909,   911,   913,
     915,   917,   919,   921,   923,   925,   929,   931,   933,   935,
     938,   940,   942,   944,   946,   948,   950,   952,   954,   956,
     959,   961,   963,   965,   967,   969,   971,   973,   975,   977,
     979,   981,   983,   985,   987,   989,   991,   995,   996,   999,
    1002,  1004,  1006,  1008,  1010,  1012,  1014,  1016,  1018,  1020,
    1022,  1025,  1029,  1032,  1035,  1037,  1039,  1042,  1044,  1046,
    1050,  1052,  1056,  1060,  1063,  1067,  1071,  1073,  1074,  1075,
    1076,  1078,  1080,  1082,  1089,  1091,  1093,  1095,  1097,  1099,
    1105,  1110,  1125,  1127,  1128,  1130,  1133,  1135,  1137,  1139,
    1145,  1156,  1159,  1160,  1161,  1165,  1167,  1171,  1181,  1184,
    1192,  1208,  1213,  1215,  1217,  1221,  1223,  1227,  1231,  1235,
    1239,  1243,  1247,  1251,  1255,  1258,  1261,  1264,  1270,  1272,
    1276,  1278,  1280,  1283,  1291,  1301,  1305,  1309,  1313,  1319,
    1324,  1333,  1334,  1337,  1339,  1341,  1343,  1345,  1347,  1349,
    1351,  1353,  1355,  1359,  1361,  1364,  1369,  1373,  1379,  1387,
    1393,  1398,  1401,  1407,  1415,  1417,  1419,  1421,  1423,  1437,
    1439,  1449,  1453,  1455,  1459,  1461,  1465,  1506,  1510,  1516,
    1522,  1524,  1527,  1529,  1531,  1540,  1542,  1544,  1547,  1558,
    1560,  1562,  1567,  1569,  1571,  1576,  1579,  1580,  1583,  1585,
    1587,  1591,  1592,  1595,  1601,  1603,  1608,  1611,  1612,  1615,
    1618,  1621,  1626,  1627,  1630,  1635,  1636,  1639,  1643,  1644,
    1647,  1653,  1654,  1657,  1661,  1665,  1667,  1671,  1675,  1677,
    1681,  1683,  1686,  1688,  1691,  1695,  1698,  1700,  1704,  1710,
    1715,  1719,  1723,  1728,  1732,  1734,  1738,  1746,  1753,  1758,
    1761,  1763,  1767,  1769,  1773,  1775,  1777,  1781,  1784,  1788,
    1792,  1797,  1801,  1803,  1805,  1808,  1814,  1816,  1818,  1829,
    1840,  1843,  1848,  1850,  1863,  1865,  1877,  1879,  1882,  1888,
    1893,  1899,  1901,  1903,  1905,  1907,  1909,  1911,  1913,  1922,
    1925,  1930,  1932,  1934,  1936,  1938,  1940,  1944,  1946,  1950,
    1952,  1956,  1957,  1960,  1962,  1964,  1968,  1969,  1972,  1974,
    1976,  1980,  1981,  1984,  1986,  1988,  1992,  1993,  1996,  2000,
    2002,  2008,  2011,  2014,  2018,  2023,  2031,  2043,  2044,  2047,
    2049,  2051,  2055,  2057,  2059,  2063,  2071,  2081,  2083,  2088,
    2090,  2094,  2095,  2098,  2105,  2113,  2121,  2129,  2138,  2147,
    2153,  2160,  2167,  2180,  2196,  2213,  2230,  2238,  2246,  2264,
    2271,  2278,  2289,  2305,  2312,  2323,  2331,  2343,  2358,  2370,
    2383,  2397,  2413,  2430,  2448,  2461,  2475,  2490,  2507,  2525,
    2544,  2552,  2564,  2573,  2586,  2596,  2610,  2619,  2620,  2623,
    2629,  2635,  2642,  2650,  2661,  2663,  2665,  2669,  2671,  2673,
    2675,  2677,  2681,  2687,  2689,  2697,  2699,  2705,  2712,  2719,
    2726,  2734,  2741,  2749,  2757,  2765,  2776,  2778,  2784,  2792,
    2795,  2798,  2804,  2806,  2812,  2822,  2824,  2831,  2838,  2845,
    2852,  2860,  2869,  2878,  2887,  2895,  2903,  2911,  2919,  2928,
    2938,  2948,  2960,  2962,  2964,  2966,  2969,  2971,  2975,  2977,
    2979,  2981,  2985,  2990,  2995,  2997,  3002,  3004,  3006,  3008,
    3012,  3016,  3020,  3022,  3026,  3027,  3031,  3032,  3036,  3041,
    3046,  3047,  3051,  3058,  3065,  3072,  3081,  3082,  3089,  3092,
    3095,  3100,  3102,  3106,  3111,  3117,  3119,  3123,  3129,  3131,
    3135,  3141,  3146,  3151,  3156,  3163,  3165,  3167,  3169,  3173,
    3175,  3179,  3180,  3183,  3187,  3188,  3191,  3197,  3199,  3201,
    3205,  3211,  3216,  3220,  3226,  3232,  3237,  3244,  3249,  3258,
    3263,  3277,  3279,  3281,  3286,  3289,  3297,  3299,  3300,  3301,
    3302,  3319,  3340,  3342,  3346,  3349,  3351,  3353,  3355,  3357,
    3361,  3363,  3365,  3367,  3371,  3374,  3376,  3378,  3380,  3382,
    3384,  3389,  3392,  3396,  3401,  3406,  3411,  3413,  3418,  3420,
    3422,  3424,  3426,  3433,  3435,  3442,  3444,  3446,  3448,  3450,
    3452,  3454,  3456,  3458,  3460,  3462,  3468,  3470,  3472,  3479,
    3486,  3493,  3495,  3498,  3500,  3503,  3506,  3509,  3512,  3514,
    3516,  3536,  3539,  3541,  3544,  3546,  3549,  3551,  3553,  3555,
    3557,  3559,  3561,  3563,  3565,  3569,  3578,  3599,  3617,  3622,
    3628,  3634,  3641,  3643,  3646,  3648,  3650,  3652,  3654,  3656,
    3658,  3661,  3663,  3665,  3667,  3669,  3673,  3676,  3680,  3684,
    3687,  3691,  3695,  3721,  3723,  3725,  3727,  3733,  3738,  3743,
    3756,  3770,  3787,  3803,  3813,  3815,  3818,  3820,  3822,  3824,
    3826,  3828,  3830,  3835,  3841,  3848,  3850,  3853,  3855,  3857,
    3859,  3861,  3864,  3867,  3870,  3873,  3877,  3881,  3883,  3885,
    3888,  3890,  3892,  3894,  3896,  3900,  3903,  3914,  3920,  3928,
    3935,  3942,  3950,  3959,  3969,  3975,  3981,  3987,  3993,  3996,
    4001,  4006,  4016,  4021,  4027,  4034,  4039,  4044,  4049,  4050,
    4051,  4054,  4055,  4058,  4062,  4065,  4070,  4071,  4074,  4078,
    4081,  4092,  4095,  4106,  4112,  4114,  4117,  4128,  4130,  4133,
    4138,  4140,  4142,  4144,  4147,  4152,  4155,  4160,  4163,  4166,
    4168,  4170,  4172,  4175,  4181,  4187,  4194,  4196,  4199,  4203,
    4206,  4211,  4248,  4252,  4253,  4255,  4265,  4271,  4279,  4280,
    4284,  4287,  4290,  4292,  4295,  4297,  4301,  4302,  4307,  4312,
    4317,  4322,  4327,  4332,  4337,  4342,  4347,  4352,  4357,  4362,
    4367,  4372,  4377,  4382,  4387,  4392,  4397,  4402,  4407,  4412,
    4418,  4424,  4430,  4436,  4442,  4448,  4454,  4460,  4466,  4471,
    4476,  4481,  4486,  4491,  4496,  4500,  4506,  4507,  4509,  4516,
    4517,  4528,  4529,  4540,  4544,  4547,  4554,  4558,  4562,  4573,
    4575,  4577,  4579,  4581,  4583,  4585,  4588,  4590,  4592,  4596,
    4600,  4604,  4608,  4612,  4616,  4617,  4621,  4625,  4629,  4631,
    4633,  4636,  4640,  4642,  4644,  4646,  4648,  4651,  4655,  4657,
    4659,  4661,  4663,  4665,  4668,  4672,  4674,  4676,  4678,  4680,
    4683,  4687,  4689,  4691,  4693,  4698,  4699,  4704,  4708,  4712,
    4715,  4722,  4724,  4729,  4731,  4733,  4735
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ACTUATORS", "AERO", "AEROH", "AEROTYPE",
  "ALPROC", "AMAT", "ANALYSIS", "ARCLENGTH", "ATTRIBUTES",
  "ANGULAROUTTYPE", "AUGMENT", "AUGMENTTYPE", "AVERAGED", "ATDARB", "ACOU",
  "ATDDNB", "ATDROB", "ARPACK", "ATDDIR", "ATDNEU", "AXIHDIR", "AXIHNEU",
  "AXINUMMODES", "AXINUMSLICES", "AXIHSOMMER", "AXIMPC", "AUXCOARSESOLVER",
  "ACMECNTL", "ADDEDMASS", "AEROEMBED", "AUGMENTED", "BLOCKDIAG",
  "BOFFSET", "BUCKLE", "BGTL", "BMPC", "BINARYINPUT", "BINARYOUTPUT",
  "BLOCKSIZE", "CHECKTOKEN", "COARSESOLVER", "COEF", "CFRAMES",
  "COLLOCATEDTYPE", "CONVECTION", "COMPOSITE", "CONDITION", "CONTROL",
  "CORNER", "CORNERTYPE", "CURVE", "CCTTOL", "CCTSOLVER", "CRHS",
  "COUPLEDSCALE", "CONTACTSURFACES", "CMPC", "CNORM", "COMPLEXOUTTYPE",
  "CONSTRMAT", "CASES", "CONSTRAINEDSURFACES", "CSFRAMES", "CSTYPE",
  "CONSTANT", "CONWEP", "DAMPING", "DblConstant", "DEM", "DIMASS", "DISP",
  "DIRECT", "DLAMBDA", "DP", "DYNAM", "DETER", "DECOMPOSE", "DECOMPFILE",
  "DMPC", "DEBUGCNTL", "DEBUGICNTL", "CONSTRAINTS", "MULTIPLIERS",
  "PENALTY", "ELLUMP", "EIGEN", "EFRAMES", "ELSCATTERER", "END",
  "ELHSOMMERFELD", "ETEMP", "EXPLICIT", "EXTFOL", "EPSILON",
  "ELEMENTARYFUNCTIONTYPE", "FABMAT", "FACE", "FACOUSTICS", "FETI",
  "FETI2TYPE", "FETIPREC", "FFP", "FFPDIR", "FITALG", "FNAME", "FLUX",
  "FORCE", "FRONTAL", "FETIH", "FIELDWEIGHTLIST", "FILTEREIG", "FLUID",
  "FREQSWEEP", "FREQSWEEP1", "FREQSWEEP2", "FREQSWEEPA", "FSGL",
  "FSINTERFACE", "FSISCALING", "FSIELEMENT", "NOLOCALFSISPLITING",
  "FSICORNER", "FFIDEBUG", "FAILSAFE", "FRAMETYPE", "GEPS", "GLOBALTOL",
  "GRAVITY", "GRBM", "GTGSOLVER", "GLOBALCRBMTOL", "GROUP", "GROUPTYPE",
  "GOLDFARBTOL", "GOLDFARBCHECK", "HDIRICHLET", "HEAT", "HFETI", "HNEUMAN",
  "HSOMMERFELD", "HFTT", "HELMHOLTZ", "HNBO", "HELMMF", "HELMSO", "HSCBO",
  "HWIBO", "HZEM", "HZEMFILTER", "HLMPC", "HERMITIAN", "HESSIAN", "IACC",
  "IDENTITY", "IDIS", "IDIS6", "IntConstant", "INTERFACELUMPED", "ITEMP",
  "ITERTYPE", "IVEL", "INCIDENCE", "IHDIRICHLET", "IHDSWEEP", "IHNEUMANN",
  "ISOLVERTYPE", "INPC", "INFINTY", "JACOBI", "KEYLETTER", "KRYLOVTYPE",
  "KIRLOC", "LAYC", "LAYN", "LAYD", "LAYO", "LAYMAT", "LFACTOR", "LMPC",
  "LOAD", "LOADCASE", "LOBPCG", "LOCALSOLVER", "LINESEARCH", "LUMPED",
  "MASS", "MATERIALS", "MATLAB", "MAXITR", "MAXORTHO", "MAXVEC", "MODAL",
  "MPCPRECNO", "MPCPRECNOID", "MPCTYPE", "MPCTYPEID", "MPCSCALING",
  "MPCELEMENT", "MPCBLOCKID", "MPCBLK_OVERLAP", "MFTT", "MRHS", "MPCCHECK",
  "MUMPSICNTL", "MUMPSCNTL", "MECH", "MODDAMP", "MODEFILTER", "MOMENTTYPE",
  "MPROJECT", "MAXIMUM", "NDTYPE", "NEIGPA", "NEWMARK", "NewLine",
  "NEWTON", "NL", "NLMAT", "NLPREC", "NOCOARSE", "NODETOKEN", "NONINPC",
  "NSBSPV", "NLTOL", "NUMCGM", "NOSECONDARY", "NFRAMES", "OPTIMIZATION",
  "OUTPUT", "OUTPUT6", "OUTPUTFRAME", "QSTATIC", "QLOAD", "PITA",
  "PITADISP6", "PITAVEL6", "NOFORCE", "MDPITA", "GLOBALBASES",
  "LOCALBASES", "TIMEREVERSIBLE", "REMOTECOARSE", "ORTHOPROJTOL",
  "READINITSEED", "JUMPCVG", "JUMPOUTPUT", "PRECNO", "PRECONDITIONER",
  "PRELOAD", "PRESSURE", "PRINTMATLAB", "PROJ", "PIVOT", "PRECTYPE",
  "PRECTYPEID", "PICKANYCORNER", "PADEPIVOT", "PROPORTIONING", "PLOAD",
  "PADEPOLES", "POINTSOURCE", "PLANEWAVE", "PTOL", "PLANTOL", "PMAXIT",
  "PIECEWISE", "RADIATION", "RAYDAMP", "RBMFILTER", "RBMSET", "READMODE",
  "REBUILD", "REDFOL", "RENUM", "RENUMBERID", "REORTHO", "RESTART",
  "RECONS", "RECONSALG", "REBUILDCCT", "RANDOM", "RPROP", "RNORM",
  "REVERSENORMALS", "ROTVECOUTTYPE", "RESCALING", "SCALING", "SCALINGTYPE",
  "STRDAMP", "SDETAFT", "SENSORS", "SOLVERTYPE", "SHIFT", "SPOOLESTAU",
  "SPOOLESSEED", "SPOOLESMAXSIZE", "SPOOLESMAXDOMAINSIZE",
  "SPOOLESMAXZEROS", "SPOOLESMSGLVL", "SPOOLESSCALE", "SPOOLESPIVOT",
  "SPOOLESRENUM", "SPARSEMAXSUP", "SPARSEDEFBLK", "STATS", "STRESSID",
  "SUBSPACE", "SURFACE", "SAVEMEMCOARSE", "SPACEDIMENSION", "SCATTERER",
  "STAGTOL", "SCALED", "SWITCH", "STABLE", "SUBTYPE", "STEP", "SOWER",
  "SHELLTHICKNESS", "SURF", "SPRINGMAT", "TANGENT", "TDENFORCE", "TEMP",
  "TIME", "TOLEIG", "TOLFETI", "TOLJAC", "TOLPCG", "TOPFILE", "TOPOLOGY",
  "TRBM", "THERMOE", "THERMOH", "TETT", "TOLCGM", "TURKEL", "TIEDSURFACES",
  "THETA", "HRC", "THIRDNODE", "THERMMAT", "TDENFORC", "TESTULRICH",
  "THRU", "TRIVIAL", "USE", "USERDEFINEDISP", "USERDEFINEFORCE", "UPROJ",
  "UNSYMMETRIC", "USING", "VERSION", "WETCORNERS", "YMTT", "ZERO",
  "BINARY", "GEOMETRY", "DECOMPOSITION", "GLOBAL", "MATCHER", "CPUMAP",
  "NODALCONTACT", "MODE", "FRIC", "GAP", "OUTERLOOP", "EDGEWS", "WAVETYPE",
  "ORTHOTOL", "IMPE", "FREQ", "DPH", "WAVEMETHOD", "MATSPEC", "MATUSAGE",
  "BILINEARPLASTIC", "FINITESTRAINPLASTIC", "LINEARELASTIC",
  "STVENANTKIRCHHOFF", "LINPLSTRESS", "READ", "OPTCTV",
  "ISOTROPICLINEARELASTIC", "NEOHOOKEAN",
  "ISOTROPICLINEARELASTICJ2PLASTIC",
  "ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS", "HYPERELASTIC",
  "MOONEYRIVLIN", "HENCKY", "LOGSTRAINPLASTIC", "SVKPLSTRESS",
  "SURFACETOPOLOGY", "MORTARTIED", "MORTARSCALING",
  "MORTARINTEGRATIONRULE", "SEARCHTOL", "STDMORTAR", "DUALMORTAR",
  "WETINTERFACE", "NSUBS", "EXITAFTERDEC", "SKIP", "OUTPUTMEMORY",
  "OUTPUTWEIGHT", "SOLVER", "SPNNLSSOLVERTYPE", "MAXSIZE", "WEIGHTLIST",
  "GMRESRESIDUAL", "SLOSH", "SLGRAV", "SLZEM", "SLZEMFILTER", "PDIR",
  "HEFSB", "HEFRS", "HEINTERFACE", "SNAPFI", "PODROB", "TRNVCT", "OFFSET",
  "ORTHOG", "SVDTOKEN", "CONVERSIONTOKEN", "CONVFI", "SAMPLING",
  "SNAPSHOTPROJECT", "PODSIZEMAX", "REFSUBSTRACT", "TOLER",
  "NORMALIZETOKEN", "FNUMBER", "SNAPWEIGHT", "ROBFI", "STAVCT", "VELVCT",
  "ACCVCT", "CONWEPCFG", "PSEUDOGNAT", "PSEUDOGNATELEM", "VECTORNORM",
  "REBUILDFORCE", "SAMPNODESLOT", "REDUCEDSTIFFNESS", "UDEIMBASIS",
  "FORCEROB", "DEIMINDICES", "UDEIMINDICES", "SVDFORCESNAP",
  "USEMASSNORMALIZEDBASIS", "OPTSENSITIVITY", "SENSITIVITYID",
  "SENSITIVITYTYPE", "SENSITIVITYMETHOD", "QRFACTORIZATION", "QMATRIX",
  "RMATRIX", "XMATRIX", "$accept", "FinalizedData", "All", "Component",
  "Noninpc", "Inpc", "Group", "Random", "Impe", "PadePivotInfo",
  "PadePolesInfo", "FreqSweep", "ReconsInfo", "BinarySpec", "AnalysisInfo",
  "Decompose", "WeightList", "FieldWeightList", "MFTTInfo", "HFTTInfo",
  "LoadCase", "Composites", "Cframes", "CoefInfo", "CoefList", "LaycInfo",
  "LaynInfo", "LaydInfo", "LayoInfo", "LayData", "LayoData", "LayMat",
  "LayMatData", "DiscrMasses", "Gravity", "Restart", "LoadCInfo",
  "UseCInfo", "SensorLocations", "ActuatorLocations", "UsdfLocations",
  "UsddLocations", "OptSensitivity", "SenInfo", "Output", "OutInfo",
  "DynInfo", "SloshInfo", "MassInfo", "CondInfo", "TopInfo", "DynamInfo",
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
  "AtdRobinScatterer", "FarFieldPattern", "FarFieldPatternDirs", "AxiHDir",
  "AxiHD", "AxiHNeu", "AxiHN", "AxiNumModes", "AxiNumSlices", "AxiHSommer",
  "AxiHSData", "AxiMPC", "AxiLmpc", "Mode", "IDisp", "IDisp6",
  "IDisp6Pita", "IVel6Pita", "IVel", "ITemp", "ETemp", "NeumanBC",
  "ModalNeumanBC", "BCDataList", "ModalValList", "TBCDataList", "YMTTable",
  "YMTTList", "TETTable", "TETTList", "SDETAFTable", "SDETAFList",
  "LMPConstrain", "MPCList", "MPCHeader", "MPCLine", "ComplexLMPConstrain",
  "ComplexMPCList", "ComplexMPCHeader", "ComplexMPCLine",
  "ComplexNeumanBC", "ComplexBCDataList", "Materials", "MatData",
  "ElemSet", "FaceSet", "MortarCondition", "WetInterface", "TiedSurfaces",
  "FSInterface", "HEVibInfo", "HEVInterfaceElement", "HEInterface",
  "ContactSurfaces", "AcmeControls", "NodeSet", "Node", "Element",
  "NodeNums", "BC_Data", "ModalVal", "TBC_Data", "ComplexBC_Data",
  "ConstrainedSurfaceFrameDList", "FrameDList", "Frame", "NodalFrameDList",
  "NodalFrame", "BoffsetList", "Attributes", "Ellump", "ReducedStiffness",
  "UDeimBasis", "SampNodeSlot", "Pressure", "Lumped", "Preload", "Statics",
  "CasesList", "IterSolver", "Solver", "OldHelmInfo", "FAcousticData",
  "Constraints", "ConstraintOptionsData", "HelmInfo", "IncidenceList",
  "IncidenceVector", "KirchhoffLocations", "FFPDirList", "FFPDirVector",
  "HelmMFInfo", "FAcousticDataMF", "HelmSOInfo", "FAcousticDataSO",
  "DEMInfo", "AlProc", "NLInfo", "NewtonInfo", "OrthoInfo", "Control",
  "Optimization", "NodalContact", "MatSpec", "MatUsage", "FloatList",
  "StringList", "Renumbering", "SvdToken", "SvdOption", "DeimIndices",
  "UDeimIndices", "Sampling", "SnapshotProject", "SamplingOption",
  "ConwepConfig", "ConversionToken", "ConversionOption", "Integer",
  "Float", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
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
     715,   716
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   462,   463,   464,   464,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   465,   465,   465,
     465,   465,   465,   465,   465,   465,   465,   466,   467,   468,
     468,   468,   468,   469,   469,   470,   470,   470,   470,   470,
     470,   470,   470,   470,   470,   470,   470,   470,   470,   470,
     470,   470,   471,   472,   472,   473,   473,   474,   474,   475,
     475,   475,   475,   475,   475,   475,   475,   475,   475,   476,
     477,   477,   477,   477,   477,   477,   477,   477,   477,   477,
     478,   478,   479,   479,   479,   479,   479,   480,   480,   480,
     481,   481,   481,   482,   482,   482,   482,   482,   483,   483,
     483,   483,   483,   483,   484,   484,   485,   485,   486,   486,
     487,   487,   488,   488,   489,   489,   490,   490,   491,   491,
     491,   492,   493,   493,   494,   494,   494,   494,   495,   495,
     495,   496,   496,   497,   497,   497,   497,   498,   499,   500,
     501,   502,   503,   504,   504,   505,   505,   506,   506,   506,
     506,   506,   506,   506,   506,   506,   507,   507,   507,   507,
     507,   507,   507,   507,   507,   507,   507,   507,   507,   507,
     507,   507,   507,   507,   507,   507,   507,   507,   507,   507,
     507,   507,   507,   507,   507,   507,   507,   508,   508,   508,
     508,   508,   508,   508,   508,   508,   508,   508,   508,   508,
     508,   508,   508,   508,   508,   508,   508,   508,   508,   508,
     509,   509,   510,   511,   511,   512,   513,   513,   513,   513,
     513,   513,   513,   513,   513,   513,   513,   513,   513,   513,
     514,   515,   516,   516,   516,   516,   517,   517,   517,   517,
     518,   519,   519,   519,   519,   520,   520,   521,   522,   522,
     522,   522,   522,   522,   522,   523,   523,   524,   525,   526,
     527,   528,   529,   530,   531,   531,   531,   531,   532,   532,
     533,   533,   533,   534,   534,   535,   536,   537,   538,   538,
     538,   539,   539,   540,   540,   540,   540,   540,   540,   540,
     540,   540,   540,   541,   541,   541,   542,   543,   543,   544,
     544,   544,   545,   546,   547,   547,   547,   547,   547,   548,
     548,   549,   550,   550,   551,   551,   552,   553,   553,   554,
     555,   555,   555,   555,   555,   556,   556,   556,   556,   557,
     557,   557,   558,   558,   558,   559,   560,   560,   561,   561,
     561,   562,   562,   563,   564,   564,   565,   566,   566,   567,
     567,   567,   568,   568,   569,   570,   570,   571,   572,   572,
     573,   574,   574,   575,   576,   577,   577,   578,   579,   579,
     580,   580,   581,   581,   582,   582,   583,   583,   584,   585,
     585,   586,   587,   587,   588,   588,   589,   590,   591,   591,
     592,   592,   593,   593,   594,   594,   594,   595,   595,   595,
     595,   595,   596,   596,   596,   596,   597,   597,   597,   597,
     597,   597,   598,   598,   599,   599,   600,   600,   600,   601,
     602,   603,   603,   603,   603,   603,   603,   603,   603,   604,
     604,   605,   605,   605,   605,   605,   605,   606,   606,   607,
     607,   608,   608,   609,   609,   609,   610,   610,   611,   611,
     611,   612,   612,   613,   613,   613,   614,   614,   615,   615,
     615,   616,   616,   616,   616,   616,   617,   618,   618,   619,
     619,   619,   620,   620,   620,   621,   621,   622,   622,   623,
     623,   624,   624,   625,   625,   625,   625,   625,   625,   625,
     625,   625,   625,   625,   625,   625,   625,   625,   625,   625,
     625,   625,   625,   625,   625,   625,   625,   625,   625,   625,
     625,   625,   625,   625,   625,   625,   625,   625,   625,   625,
     625,   625,   625,   625,   625,   625,   625,   626,   626,   627,
     627,   627,   627,   627,   628,   628,   628,   628,   628,   628,
     628,   628,   628,   629,   629,   630,   630,   630,   630,   630,
     630,   630,   630,   630,   630,   630,   631,   631,   631,   632,
     632,   633,   634,   634,   634,   635,   635,   635,   635,   635,
     635,   635,   635,   635,   635,   635,   635,   635,   635,   635,
     635,   635,   636,   636,   636,   636,   637,   637,   638,   638,
     638,   638,   638,   639,   640,   640,   641,   641,   641,   641,
     642,   643,   644,   644,   645,   645,   646,   646,   647,   647,
     648,   648,   649,   649,   649,   649,   650,   650,   651,   651,
     651,   651,   651,   651,   651,   651,   651,   651,   651,   651,
     651,   651,   651,   651,   651,   652,   652,   652,   652,   653,
     653,   654,   654,   654,   655,   655,   655,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   657,
     657,   658,   658,   658,   658,   658,   659,   659,   659,   659,
     659,   659,   660,   660,   661,   661,   661,   661,   661,   661,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   662,   662,   662,   662,   662,
     662,   662,   662,   662,   662,   663,   664,   665,   665,   666,
     666,   666,   666,   666,   666,   666,   666,   666,   666,   666,
     667,   667,   667,   667,   667,   667,   667,   667,   667,   667,
     667,   668,   668,   669,   670,   670,   671,   671,   672,   673,
     674,   675,   676,   677,   678,   678,   679,   679,   679,   679,
     679,   679,   679,   679,   679,   679,   679,   679,   679,   679,
     679,   679,   679,   679,   679,   679,   679,   679,   679,   679,
     680,   680,   681,   682,   682,   682,   683,   683,   684,   684,
     684,   684,   684,   684,   684,   684,   685,   685,   685,   685,
     685,   685,   685,   685,   685,   685,   685,   685,   685,   685,
     685,   685,   685,   685,   685,   685,   685,   685,   685,   685,
     685,   685,   685,   685,   685,   685,   685,   685,   685,   685,
     685,   685,   685,   685,   685,   685,   686,   686,   686,   687,
     687,   688,   688,   689,   689,   689,   690,   690,   691,   691,
     691,   691,   691,   691,   691,   691,   691,   691,   691,   692,
     692,   693,   693,   694,   694,   694,   695,   695,   696,   696,
     696,   696,   696,   696,   696,   696,   696,   696,   696,   696,
     696,   696,   696,   696,   696,   696,   696,   696,   696,   696,
     696,   696,   696,   696,   696,   697,   697,   698,   698,   699,
     699,   700,   700,   701,   701,   701,   701
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
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
       1,     1,     1,     1,     1,     1,     1,     5,     4,     2,
       5,     6,     6,     2,     6,     5,     7,     7,     8,    13,
       8,     7,     8,     8,     9,    14,     9,     3,     2,     2,
       2,     2,     3,     2,     4,     3,     4,     4,     6,     2,
       4,     5,     4,     5,     4,     4,     4,     4,     4,     3,
       2,     4,     4,     3,     3,     3,     3,     3,     3,     3,
       2,     4,     2,     4,     4,     4,     4,     2,     3,     4,
       2,     3,     4,     2,     3,     4,     5,     5,     2,     2,
       2,     2,     2,     2,     2,     2,     4,    10,     4,     5,
       3,     2,     3,     2,     3,     2,     3,     2,    11,    13,
      14,     5,     2,     2,     7,     9,    11,    12,     2,     5,
       6,     2,     5,     2,     5,     4,     4,     3,     3,     3,
       3,     3,     3,     2,     3,     3,     2,     2,     2,     3,
       3,     3,     3,     4,     4,     3,     3,     5,     4,     5,
       6,     7,     3,     5,     3,     5,     4,     5,     6,     7,
       3,     2,     3,     2,     2,     3,     2,     2,     2,     3,
       2,     3,     2,     2,     2,     2,     2,     1,     2,     3,
       3,     3,     2,     5,     3,     3,     3,     2,     3,     2,
       3,     4,     4,     5,     3,     3,     2,     4,     2,     3,
       2,     3,     2,     5,     2,     2,     2,     2,     2,     2,
       3,     4,     4,     8,     4,     4,     3,     4,     3,     5,
       4,     5,     2,     3,     3,     3,     2,     4,     1,     3,
       1,     2,     2,     2,     3,     5,     6,     6,     3,     4,
       6,     5,     3,     4,     3,     1,     2,     6,     2,     2,
       2,     2,     2,     4,     5,     4,     2,     6,     3,     2,
       3,     2,     4,     1,     2,     2,     2,     5,     5,     6,
       7,     3,     1,     1,     2,     1,     1,     1,     2,     1,
       1,     2,     1,     4,     5,     3,     3,     6,     5,     6,
       3,     2,     5,     6,     2,     2,     7,     9,     3,     2,
       6,     3,     1,     2,     3,     2,     3,     1,     2,     4,
       2,     4,     6,     8,     5,     2,     3,     4,     5,     2,
       3,     6,     2,     3,     6,     3,     1,     2,     4,     5,
       6,     3,     2,     4,     1,     2,     3,     1,     2,     4,
       5,     6,     3,     2,     4,     3,     2,     4,     3,     2,
       4,     3,     2,     4,     3,     1,     2,     3,     3,     4,
       3,     2,     3,     2,     5,     2,     3,     4,     4,     9,
       2,     2,     9,     2,     3,     4,     3,     3,     7,     2,
       3,     4,     2,     2,     5,     7,    10,     3,     4,     4,
       5,     4,     2,     3,     3,     4,     3,     2,     9,     6,
       2,     2,     3,     9,     3,     9,     2,     3,     4,     3,
       3,     2,     3,     2,     7,     9,     8,    10,     3,     5,
       6,     1,     2,     6,     7,     8,     9,     1,     2,     1,
       2,     2,     3,     6,     4,     7,     2,     3,     6,     4,
       7,     2,     3,     6,     4,     7,     2,     3,     2,     2,
       3,     2,     3,     5,     4,     6,     4,     2,     3,     2,
       2,     3,     5,     3,     2,     5,     4,     3,     4,     1,
       2,     3,     2,    16,    19,    19,    20,    23,    23,     9,
      14,    16,    30,    12,    13,    14,     5,     6,    16,    12,
       3,     4,    13,     5,     6,     6,     7,    11,     8,     9,
      10,     9,    10,    11,    10,    11,    12,    11,    12,    13,
       5,     6,     6,     7,     7,     8,     4,     3,     2,     3,
       4,     5,     6,     5,     4,     5,     5,     6,     7,     7,
       8,     7,     8,     4,     3,     2,     5,     6,     7,     8,
      10,     6,     7,     8,     9,    11,     2,     5,     7,     3,
       2,     4,     2,     5,     7,     2,     5,     6,     7,     8,
      10,    11,    13,    14,     6,     7,     8,     9,    11,    12,
      14,    15,     3,     3,     3,     3,     3,     2,     5,     4,
       3,     7,     6,     4,     1,     2,     4,     3,     5,     4,
       3,     3,     5,     4,     2,     2,     2,     2,    11,     4,
       2,     2,    11,    14,    12,    15,     2,     7,     2,     5,
       6,     4,     6,     7,     7,     6,     8,     9,     7,     9,
      10,     5,     5,     7,     8,     2,     4,     3,     3,     2,
       3,     2,     4,     3,     2,     4,     5,     2,     3,     4,
       5,     5,     5,     6,     6,     6,     7,     7,     8,     2,
       3,     2,     4,     6,     6,     8,     2,     2,     2,     4,
       3,     5,     1,     2,     2,     4,     4,     4,     4,     4,
       2,     3,     2,     3,     3,     3,     4,     5,     6,     7,
       4,     3,     3,     3,     4,     5,     2,     4,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     4,
       4,     3,     3,     4,     4,     4,     4,     4,     4,     3,
       3,     3,     4,     3,     4,     4,     4,     4,     4,     2,
       3,     3,     3,     3,     3,     3,     3,     2,     2,     3,
       2,     3,     2,     3,     4,     3,     3,     4,     4,     5,
       5,     6,     3,     4,     3,     3,     3,     3,     3,     2,
       3,     3,     3,     3,     3,     4,     4,     4,     4,     4,
       4,     4,     3,     2,     2,     3,     3,     2,     3,     4,
       5,     6,     7,     2,     3,     4,     3,     3,     3,     3,
       2,     2,     3,     4,     5,     3,     4,     3,     3,     3,
       3,     5,     5,     4,     4,     6,     6,     3,     3,     4,
       3,     3,     3,     3,     3,     3,     2,     3,     4,     1,
       2,     3,     4,     5,     1,     2,     3,     2,     3,     4,
       2,     2,     3,     4,     5,     4,     4,     4,     2,     2,
       2,     1,     2,     4,     6,     5,     1,     2,     4,     3,
       2,     3,     2,     2,     1,     1,     2,     3,     3,     3,
       3,     4,     4,     5,     7,     5,     7,     4,     5,     3,
       4,     4,     5,     6,     7,     8,     3,     4,     6,     2,
       3,     4,     2,     2,    10,    12,     3,     4,     2,     4,
       7,     9,     9,     8,    11,    12,     2,     9,    10,    12,
       9,    10,    12,     9,    10,    12,     9,     7,     5,     9,
       7,     5,     9,     7,     5,     8,    10,     8,    10,     7,
       7,     8,     8,     9,    10,    10,    11,    12,    24,     7,
       8,     9,    10,    11,     5,     5,     2,     4,     5,     0,
       2,     0,     2,     4,     6,     8,     2,     3,     2,     2,
       2,     2,     2,     3,     2,     3,     2,     2,     1,     2,
       3,     2,     3,     2,     3,     2,     2,     3,     2,     2,
       3,     4,     2,     2,     2,     2,     3,     4,     2,     2,
       2,     2,     2,     3,     4,     2,     2,     2,     2,     3,
       4,     2,     2,     2,     2,     2,     3,     2,     3,     2,
       3,     1,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
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
       0,     0,     3,   121,   122,   123,   124,   125,   110,    67,
     116,   117,   118,    49,    50,    51,    46,    47,    48,    45,
      57,    62,    63,    64,    34,    35,    37,    36,   146,    43,
      38,    40,    76,    75,    83,   307,    39,    42,    69,    70,
      71,    72,   120,    73,    74,    55,    56,    61,    58,    59,
      60,   137,   101,   102,   103,   100,     6,   135,    80,    81,
      77,    78,    79,    82,    93,    94,    89,    95,    96,    97,
      98,   111,   112,   113,   114,   115,    90,    91,   104,   105,
     106,   107,   108,   109,    27,    26,    30,    28,    29,    31,
      32,    33,     7,     8,    52,    53,    54,     9,    10,    99,
      20,    11,   128,   129,   131,   130,   132,    41,   133,   134,
     138,     5,    14,    12,    13,   136,    15,    16,    18,    19,
      17,    22,    23,    24,    21,    84,   139,    86,    92,    87,
      88,    85,    44,    68,    65,    66,   119,   126,   127,    25,
     140,   141,   142,   143,   144,   145,     0,     0,     0,     0,
    1121,  1122,     0,   738,  1124,  1126,  1123,  1125,     0,     0,
       0,     0,   319,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   512,     0,     0,   736,   531,     0,   224,   449,
       0,   218,   334,  1003,   685,   429,   724,     0,   973,   248,
     424,   336,   190,     0,   939,   944,     0,     0,     0,   755,
     308,     0,   726,     0,     0,     0,   317,     0,     0,     0,
     445,     0,   541,     0,   951,   202,     0,   676,     0,   530,
     251,   386,   149,     0,     0,     0,     0,   210,     0,   950,
       0,     0,     0,     0,     0,   381,   395,   587,   522,     0,
     527,     0,     0,   536,     0,     0,     0,     0,     0,     0,
       0,     0,   242,   576,     0,   213,     0,   326,   779,     0,
     332,     0,   207,     0,   389,     0,     0,   976,     0,     0,
       0,   730,     0,     0,   267,     0,     0,   268,     0,   361,
       0,     0,     0,     0,   781,   767,     0,     0,     0,   452,
       0,   391,     0,     0,     0,  1002,   253,   153,   571,     0,
       0,   786,   312,     0,     0,   440,     0,     0,   335,     0,
       0,   379,   378,   566,   665,   328,     0,     0,     0,   561,
     179,  1008,     0,   380,     0,     0,  1016,  1056,     0,     0,
       0,     0,     0,   200,   330,     0,   382,   396,     0,     0,
       0,   682,  1066,  1117,  1083,  1086,   764,   759,   761,  1079,
    1081,   263,     0,     1,     2,     4,     0,     0,     0,     0,
       0,     0,     0,   170,   171,   168,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   219,   220,   221,   222,
     223,   225,     0,   243,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     337,   338,   339,     0,     0,     0,   362,   363,   375,     0,
       0,     0,   421,     0,     0,   425,     0,     0,     0,     0,
       0,     0,     0,     0,   462,     0,   473,     0,   476,     0,
     479,     0,   482,     0,   491,   493,   495,   500,     0,   503,
       0,   509,     0,   513,     0,     0,     0,     0,     0,     0,
       0,     0,   543,     0,   602,     0,   648,     0,     0,     0,
       0,   680,     0,     0,     0,   707,     0,   725,   727,   731,
       0,     0,     0,     0,     0,     0,     0,  1123,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   788,   787,   960,   959,   958,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   999,     0,     0,     0,     0,     0,     0,     0,
    1061,     0,     0,  1059,  1061,     0,     0,  1078,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1085,     0,     0,     0,   260,   551,     0,   368,     0,     0,
     189,   490,   492,     0,   320,     0,     0,   485,   487,     0,
     488,     0,     0,     0,   506,   507,     0,   702,   329,   952,
       0,   450,     0,     0,     0,     0,   947,   940,     0,   945,
       0,     0,   937,   309,   472,   461,   540,   559,     0,   935,
       0,   496,     0,     0,   446,     0,   542,   324,   703,     0,
     416,   599,     0,   597,     0,   455,   456,     0,   211,   475,
     969,     0,   971,     0,   481,   478,   588,     0,     0,   524,
     523,   526,   539,   537,     0,     0,     0,   420,     0,     0,
     325,     0,   577,     0,     0,   257,   214,   780,   601,   208,
     388,   311,   706,     0,   314,  1006,   271,     0,   269,   272,
       0,   270,     0,   532,   534,     0,   768,     0,     0,   453,
       0,   393,   390,     0,   517,     0,     0,     0,   572,   259,
     318,     0,   466,   467,     0,   315,   316,   647,     0,     0,
     567,   258,   262,   261,     0,   562,     0,     0,     0,     0,
       0,     0,   167,     0,     0,   649,     0,     0,     0,   704,
     705,   664,     0,   331,   431,   432,     0,   679,   436,   437,
       0,   310,     0,     0,     0,     0,     0,   173,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   197,     0,
     199,   198,     0,   195,   196,   194,   193,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   231,     0,   233,   235,     0,   237,     0,     0,
       0,     0,     0,     0,     0,     0,   264,   266,     0,     0,
       0,     0,     0,     0,   305,   306,   304,   296,   293,   294,
     303,   300,   275,     0,   302,   297,     0,   291,   298,     0,
       0,     0,   358,   348,     0,     0,   360,     0,   340,     0,
     352,   346,     0,     0,     0,     0,     0,     0,   364,     0,
     372,     0,   374,   376,     0,   428,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   501,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   548,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   758,   757,     0,   760,     0,   763,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   792,   800,     0,     0,     0,     0,     0,   816,     0,
       0,     0,     0,   858,   860,     0,     0,   893,     0,     0,
       0,   903,     0,   910,   794,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   857,     0,     0,     0,
       0,   897,     0,   849,     0,     0,   879,     0,     0,     0,
       0,     0,     0,     0,     0,   790,     0,     0,     0,     0,
     802,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   911,     0,     0,     0,     0,     0,
       0,     0,   862,     0,     0,   894,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   978,     0,     0,     0,
     996,     0,     0,   974,     0,   975,     0,   977,   979,     0,
       0,   989,   980,     0,     0,     0,  1059,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1077,  1074,  1068,  1069,  1071,  1070,  1059,  1072,
    1076,  1115,  1067,     0,  1101,  1099,  1100,  1105,  1093,  1106,
    1107,  1088,  1089,  1094,  1095,  1092,  1111,  1112,  1102,  1113,
    1108,  1114,  1098,  1080,  1082,  1084,  1087,  1119,  1118,   552,
       0,     0,   369,     0,     0,     0,   321,   322,     0,   486,
       0,   489,     0,     0,     0,   953,     0,     0,     0,   350,
       0,   941,   946,   938,   948,   560,     0,   936,   497,   498,
     966,     0,     0,     0,   385,     0,   600,     0,   598,   457,
       0,   970,   972,     0,   590,     0,   589,     0,     0,   594,
       0,   957,   961,     0,     0,     0,     0,   148,     0,     0,
     579,     0,   578,     0,   581,     0,     0,   273,   274,     0,
       0,   955,   956,   392,   394,   519,     0,   518,  1063,     0,
       0,     0,     0,   468,     0,   383,     0,     0,     0,     0,
       0,     0,  1009,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   650,     0,   654,     0,     0,     0,
     663,   433,   435,     0,   438,     0,     0,     0,     0,   415,
     557,     0,   172,     0,     0,     0,     0,   180,     0,   182,
     184,   185,   186,   187,   188,   191,   192,   201,   203,   206,
     205,   204,   209,   212,     0,     0,   215,     0,     0,   230,
     232,   234,   236,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   255,   256,   265,   276,     0,   284,     0,   282,
       0,   295,   301,   290,   299,   292,   327,   354,   356,     0,
     355,   344,   353,   342,   341,     0,   347,     0,   345,     0,
       0,   373,     0,     0,     0,   717,     0,     0,     0,     0,
     441,     0,   447,     0,     0,     0,   464,     0,     0,     0,
       0,   504,     0,   510,     0,     0,   521,   525,     0,     0,
       0,   538,     0,     0,     0,     0,   620,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   714,     0,     0,     0,
       0,     0,     0,   710,     0,     0,     0,     0,     0,     0,
     741,     0,     0,   756,   762,   765,     0,     0,     0,     0,
       0,   769,     0,     0,   782,     0,   866,     0,     0,     0,
     883,     0,     0,   930,   882,     0,     0,   865,   863,     0,
     876,   884,     0,     0,   861,   932,   933,   789,   793,   801,
     931,     0,     0,   812,   813,     0,   811,     0,   856,   855,
     859,   874,   875,   831,   832,     0,   805,     0,   909,   908,
     881,     0,     0,     0,     0,     0,   920,     0,   919,     0,
     918,   917,   854,   853,   927,   934,     0,     0,   898,   904,
       0,   880,   850,   840,   839,   896,   928,     0,     0,     0,
       0,   877,   852,   851,   803,   804,   820,   823,   821,   822,
     824,   825,   827,   826,   828,   818,   819,   907,   843,     0,
     841,     0,   906,   912,     0,   878,   892,     0,   915,   872,
       0,   895,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   987,     0,   997,   990,   991,
       0,   981,   982,     0,  1000,     0,  1007,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1057,     0,  1075,  1062,  1073,  1060,
    1116,  1090,  1096,  1103,  1109,  1120,     0,     0,   371,     0,
       0,   494,   323,   484,     0,     0,     0,   954,   333,     0,
       0,   942,   949,   721,   967,     0,   549,     0,   384,     0,
       0,     0,     0,   591,     0,     0,   593,   962,     0,   418,
       0,     0,     0,     0,   580,     0,   582,     0,     0,   147,
     402,   398,     0,     0,   520,     0,     0,     0,     0,   313,
       0,     0,     0,     0,     0,     0,     0,     0,   175,     0,
       0,     0,     0,   155,     0,     0,     0,     0,     0,     0,
     651,     0,     0,   655,     0,   656,     0,   434,     0,     0,
     150,     0,     0,   558,     0,   174,     0,   413,   177,     0,
     181,   183,   217,   216,   226,     0,     0,     0,     0,   729,
       0,     0,     0,   249,   252,   254,     0,   278,     0,     0,
     286,     0,     0,   359,   349,     0,     0,     0,     0,     0,
       0,     0,   719,     0,   716,     0,   444,     0,     0,   448,
       0,     0,   463,   465,   474,   477,   480,   483,   505,   511,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   621,     0,     0,     0,     0,   646,     0,     0,
     713,   715,   653,   666,     0,     0,   677,     0,   681,   683,
       0,   686,     0,     0,   709,     0,     0,     0,     0,     0,
     739,   751,     0,     0,   752,     0,     0,   766,   771,     0,
       0,     0,   770,     0,   772,     0,     0,   867,     0,   868,
       0,   888,   817,   887,   889,   864,   890,   891,   835,   834,
     814,   810,     0,   899,     0,   806,     0,   885,   886,   929,
       0,   924,     0,     0,     0,   923,   829,   830,   905,   845,
     847,   846,   791,   844,   842,   913,     0,   916,   873,   796,
     798,   795,   799,   797,   833,   848,   836,   837,   838,   965,
     985,     0,     0,   988,   992,     0,   983,     0,  1001,     0,
    1054,  1055,     0,     0,  1028,     0,  1031,     0,     0,     0,
       0,     0,     0,     0,     0,  1034,     0,     0,     0,  1058,
    1091,  1097,  1104,  1110,     0,     0,     0,   370,   377,     0,
       0,     0,     0,     0,   943,     0,   550,   387,   723,     0,
       0,     0,   458,     0,     0,     0,   417,   419,   423,   964,
       0,     0,   584,   403,     0,   405,   406,   407,     0,   409,
     410,   412,     0,   399,     0,  1064,     0,     0,   574,     0,
       0,   469,     0,     0,   569,     0,     0,   564,     0,     0,
       0,     0,   176,     0,     0,     0,     0,     0,   652,   657,
       0,     0,     0,   439,   152,   151,   154,   720,   414,     0,
       0,     0,     0,     0,     0,     0,     0,   250,   279,   277,
     287,   285,   283,   357,     0,   397,     0,   365,     0,   422,
       0,     0,   718,   430,     0,   442,   451,   454,   514,     0,
     529,     0,     0,     0,     0,     0,   616,     0,     0,   623,
       0,   640,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   671,   667,     0,     0,     0,     0,   694,   687,     0,
       0,   708,     0,     0,     0,     0,   740,     0,     0,   742,
       0,   745,     0,     0,     0,   774,   775,     0,     0,   773,
     783,     0,   784,   869,     0,   870,   815,   900,     0,   807,
       0,     0,   922,   921,     0,   914,     0,   998,   993,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   508,     0,   351,   968,   722,     0,   459,   596,
       0,   592,   963,   586,   583,     0,   404,   408,   411,   401,
     400,     0,     0,     0,     0,   470,     0,     0,     0,     0,
       0,   156,   157,     0,     0,     0,     0,     0,   161,   658,
     659,     0,   661,     0,   178,     0,     0,     0,     0,   241,
       0,     0,     0,   280,     0,   288,     0,   367,   366,     0,
     426,     0,     0,     0,     0,     0,     0,     0,   544,   617,
       0,     0,   625,   642,     0,     0,   624,     0,   641,     0,
       0,     0,     0,     0,     0,   672,   668,     0,     0,   678,
     684,   695,   688,     0,     0,   712,     0,     0,     0,   737,
     743,   744,   748,     0,     0,     0,   753,   777,   776,     0,
       0,   871,   901,     0,   808,     0,   926,   925,   986,   994,
       0,   984,  1010,     0,     0,     0,     0,     0,  1027,     0,
    1030,     0,     0,  1049,     0,  1039,  1040,     0,     0,     0,
       0,  1033,     0,     0,     0,     0,     0,     0,   553,     0,
       0,     0,   460,   595,   585,  1065,   573,     0,   471,   568,
       0,   563,     0,   158,   160,     0,     0,   162,   163,     0,
     660,   662,     0,   228,     0,     0,     0,   244,     0,   281,
     289,   343,     0,   443,   515,     0,     0,     0,     0,     0,
     546,     0,     0,   644,     0,   626,   643,     0,     0,     0,
       0,     0,     0,   673,   669,     0,     0,   696,   689,     0,
       0,   711,     0,     0,     0,   754,     0,   746,   778,   785,
     902,   809,   995,     0,     0,  1013,     0,     0,     0,     0,
    1035,     0,  1050,     0,  1041,     0,     0,  1042,     0,     0,
       0,  1037,     0,     0,   554,     0,   499,   502,     0,   575,
     570,   565,     0,   164,   166,     0,   229,     0,     0,     0,
       0,   427,     0,   528,   533,   535,     0,   545,     0,   628,
       0,     0,     0,   645,     0,     0,     0,     0,     0,   674,
       0,   697,     0,     0,     0,     0,   749,   747,  1012,  1011,
       0,  1017,     0,  1020,     0,  1026,  1029,     0,  1051,     0,
       0,     0,  1043,  1032,  1023,     0,     0,     0,   555,  1004,
       0,     0,     0,     0,     0,   245,     0,     0,   547,     0,
       0,   629,     0,     0,   631,     0,     0,     0,     0,     0,
       0,   609,     0,   670,     0,   690,     0,     0,     0,     0,
     750,     0,  1018,     0,  1021,     0,  1036,  1052,     0,  1044,
    1045,     0,  1024,     0,  1038,   556,     0,     0,     0,   227,
       0,     0,     0,   516,     0,   634,     0,   630,     0,     0,
     632,     0,     0,     0,     0,     0,     0,     0,   675,   698,
     691,     0,     0,     0,     0,  1014,     0,     0,     0,  1053,
       0,  1046,     0,     0,  1005,     0,     0,     0,   728,   246,
       0,     0,   635,     0,   637,     0,   633,     0,   627,     0,
       0,     0,     0,   699,     0,     0,   732,     0,  1015,  1019,
    1022,     0,  1047,  1025,   159,     0,   238,     0,   247,   613,
       0,   636,   638,     0,     0,     0,   619,     0,   692,     0,
       0,   734,     0,     0,     0,   165,     0,   614,     0,   639,
       0,   622,     0,   700,   693,     0,     0,     0,     0,   239,
       0,   615,     0,     0,   610,     0,   701,     0,   733,     0,
     240,     0,     0,     0,   735,     0,   618,   611,   603,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   604,   605,     0,     0,   606,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   607,   608,
       0,  1048,     0,     0,     0,     0,     0,     0,   612
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   170,   171,   172,   173,   174,   175,   176,   177,   533,
     534,   972,   535,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   566,  1864,   567,   568,   569,   570,  1032,
    1035,   188,   573,   189,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   579,   199,   586,   200,   201,   202,   203,
     204,   205,   206,  1313,   600,  1071,  1075,   207,   606,   607,
     208,   611,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   940,   219,   220,   601,   221,  1811,  2092,   536,
     222,   223,   224,   612,   225,   226,   227,   228,   984,   985,
     229,   988,   989,   230,   231,   232,   233,   234,   885,   886,
     235,   624,  1545,   236,   952,   953,   237,   626,   238,   628,
     239,   630,   240,   632,   837,   838,   241,   242,   243,   244,
     245,   246,   247,   248,   637,   249,   639,   250,   251,   252,
     641,   253,   643,   254,   255,   256,   257,   258,   259,   260,
     261,   262,   263,   824,  1459,   866,   264,   965,   265,   960,
     266,   948,   267,   912,   913,  1400,   268,   896,   897,  1384,
     269,   880,   270,   654,   271,   272,   273,   274,   275,   276,
     277,   661,   278,   279,   280,   281,   665,   656,  1575,   825,
    1460,   867,   881,   282,   283,   571,   284,   669,   285,   286,
     287,   288,   289,   290,   291,   292,   293,   294,  1170,   762,
     763,   295,   869,   296,   368,   297,  1391,  1392,   298,  1369,
    1370,   299,   890,   300,   892,   301,  1274,   302,   782,   303,
     304,   305,   306,   307,   308,  1309,  1304,   309,   310,   796,
     311,   312,   313,   314,   817,   797,   315,   823,  1461,  1393
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -1787
static const yytype_int16 yypact[] =
{
    5848,  -162,    37,   -55,   410,   -50,   649,   649,   649,   699,
      50,   760,    58,   147,   410,   410,   193,   203,   410,   410,
     224,   301,   410,   304,   890,   312,   319,   326,   331,   340,
     374,   377,   387,   404,   411,   488,   492,   188,   542,   927,
     553,   556,   567,   582,   601,   609,   410,   410,   990,  1011,
     611,   625,  -199,   644,   410,   667,   671,   695,   717,   736,
    1021,   739,  1059,   743,   755,   757,   758,   761,   764,   777,
     782,   783,  -173,   713,   799,   802,   410,   805,  1115,   815,
     820,   410,   821,   827,   829,    23,  1118,   834,  1149,   838,
     857,  1167,  1181,   410,   858,   870,   872,   410,   878,   116,
     766,   852,   885,   887,   410,   410,   891,   893,  1234,   410,
     410,  1244,  1246,   143,   899,   922,   924,   925,   929,   937,
     649,   938,  1267,   940,   945,   649,   649,   959,   961,   962,
     965,   966,   968,   970,   971,   804,   975,   978,   979,   981,
    -184,   985,  1272,   988,   995,   410,   410,   649,   410,   410,
     998,  1004,   649,  1005,  1006,  1010,  1012,  1014,  1015,  1016,
    1025,  1028,  1032,  1035,  1046,  1047,  1049,  1050,  1052,  -188,
    1275,  5396, -1787, -1787, -1787,  1165,   410,   530,    39, -1787,
     -15,   410,    33,   649,   649,   410,   398,   410,   410,   410,
     649,  1188, -1787, -1787, -1787, -1787, -1787, -1787,   848,   165,
    1110, -1787, -1787, -1787, -1787,    20, -1787,   471,    -8, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787,   410, -1787,    14,   410, -1787, -1787,
      72,    97,   410,   410, -1787,   410, -1787,   410,   410,   410,
     410, -1787, -1787,   410,   410,   410, -1787, -1787,   410,   410,
   -1787, -1787,   410,   410,   854,  1112,   410,   410,   410,  1117,
   -1787, -1787,   104, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
     410,   410,   410, -1787, -1787,   410,   410,   410,   410,   410,
   -1787,   410,   410,   410,   410,   410,   410,   358,   649,  1906,
     410,   109, -1787,   410,  4907, -1787, -1787,   727,   649, -1787,
   -1787, -1787,     6, -1787, -1787,  1099,   410,  -123,   410, -1787,
     -23,   433,   433,    27,   433,   901,   410,  1119,  1321,  1328,
   -1787, -1787,  1122, -1787, -1787, -1787, -1787, -1787,  1124,  1128,
     649,  1280, -1787,   649,   410,   410,  1134,   649,   649,  1135,
    1144,   410, -1787,  1147,  1161, -1787, -1787,   720, -1787, -1787,
    1166, -1787,   649,  1301, -1787, -1787, -1787,   649, -1787, -1787,
   -1787, -1787, -1787,   649,   649,  1315,   649,  1584,   360, -1787,
   -1787,  1187, -1787,   410,   410,   410, -1787,   649,  1291,  1194,
   -1787,  1195,  1220,  1199, -1787, -1787,  1202, -1787,  1206, -1787,
   -1787,   649, -1787,   410,   410,  1215,   410, -1787,  1225, -1787,
     410,   649,   649,   410,   410, -1787, -1787,   410,   410,  1228,
   -1787,  1231,   410,   410,  1237,   649,   410,  1241,   649,   410,
    1257,   649, -1787,   410,  1258, -1787,  1261, -1787, -1787,  1262,
   -1787,   410, -1787,  1277, -1787,  1278,  1281, -1787,   410,   410,
    1282, -1787,  1283,  1294, -1787,  1285,  1297, -1787,  1299, -1787,
     410,  1300,  1305,   410, -1787, -1787,  1307,  1317,  1319, -1787,
    1325,   410,  1330,   322,  1192, -1787, -1787, -1787,  1384,   410,
    1331, -1787, -1787,   649,   410, -1787,  1338,  1339, -1787,   410,
     649, -1787, -1787,  1465, -1787, -1787,  1345,   410,   410,  1533,
   -1787, -1787,   410, -1787,   -87,  1371, -1787, -1787,   444,   410,
    1372,  1374,  1462, -1787, -1787,  1379, -1787, -1787,   410,   410,
     410, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787,  1381, -1787, -1787, -1787,   130,  1334,  1392,   649,
     734,   649,  1341, -1787, -1787, -1787, -1787,  1312,  1313,  1516,
    1524,  1530,  1534,  1544,  1402,  1547,  1423,  1430,   410,  1432,
    1438,  1439,  1447,   649,   410,   410,   410,   410,   649,   649,
    1979,   410,   410,   410,   410,   410, -1787,   410,   410,   410,
     410, -1787,   -22, -1787,   649,   410,   649,   347,  1209,   -59,
     368,   426,   427,  1560,  1561,  1566,    74,   410,   649,   750,
     649,  1362,  1466,   649,  1469,  1471,    90,  1378,   649,  1380,
   -1787, -1787, -1787,   649,  1478,   649, -1787, -1787, -1787,  1483,
    1600,  1467, -1787,   649,   410, -1787,  -110,  1642,   410,   -36,
     410,   649,   649,   649, -1787,   410, -1787,   410, -1787,   410,
   -1787,   410, -1787,   410, -1787, -1787, -1787, -1787,  1497, -1787,
     410, -1787,   410, -1787,   649,  1409,  1509,   649,   649,   649,
    1536,   410, -1787,  -102, -1787,    31, -1787,   410,   410,   410,
     410, -1787,   410,   410,   410, -1787,   649, -1787, -1787, -1787,
    1947,   410,   -44,  1541,  1543,   649,  1545,  1067,   410,  1554,
     410,   410,  1900,   -28,     3,  1457,  1555,  1448,  1479,   497,
     649,  1484,   649,  1464,  1470,   410,  1468,  1481,   410,   410,
     -49,     4,  1570,  1573,   410,   649,  1574,   649,   649,  1482,
    1476,  1577,  1498,   776,  1502,   410,   947,   898,   300,  1585,
     410,   410,   410,   410,  1500,  1586,   410,   163,   410,   443,
     410,   649,   649,   649,   410,   751,  1491,   303,   -74,   649,
     410,   410,   410,   649,   410,   410,  1499,   410,   410,   410,
    1596,   410,   649,   649,   649,   410,   410,   410,  1601,  1659,
     649,  -121,   349,   756, -1787, -1787, -1787,   649,  1608,   649,
     410,   410,   826,   649,    45,   410,  1618,  1622,   649,   410,
    1624,  1625, -1787,  1736,   410,  1739,   103,   410,   410,   410,
   -1787,   410,   204,  1531, -1787,  1631,  1632,   649,  1557,  1563,
    1565,  1568,   410,  1477,   649,  1750,  1789,   410,   410,   649,
    1571,  1592,  1809,  1609,  1819,  1616,  1617,  1718,  1721,  1727,
     649,  1728,  1845,  1744,   410, -1787,  -100, -1787,   828,   649,
   -1787, -1787, -1787,   649, -1787,  1747,  1506, -1787,   410,   649,
     410,   410,   649,   649, -1787, -1787,  1748, -1787, -1787, -1787,
     831, -1787,   410,  1752,  1755,   649, -1787,   649,   649, -1787,
     473,   410, -1787, -1787, -1787, -1787,   410, -1787,   649, -1787,
    1756, -1787,  1757,   649, -1787,  1760,  1786, -1787, -1787,   836,
     410, -1787,   410,   410,   410,   410, -1787,   410, -1787, -1787,
   -1787,  1766, -1787,  1774, -1787, -1787,   410,   410,   598,   410,
   -1787, -1787,   410,   410,   649,   649,   649, -1787,   649,  1776,
   -1787,   649,   410,   410,   843, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787,   410, -1787, -1787, -1787,  1783, -1787, -1787,
    1787, -1787,   410, -1787, -1787,   410, -1787,   649,   649, -1787,
    1646, -1787, -1787,  1653, -1787,  1791,  1793,   410,   880,   410,
   -1787,   649,   410, -1787,   410, -1787, -1787, -1787,  1796,   410,
    1031, -1787,   410,   410,   410,  1306,  1799,   649,   649,   649,
     649,   649,   649,  1252,  1803, -1787,  1805,   649,   -51, -1787,
   -1787, -1787,  1807, -1787,   410, -1787,   844, -1787,   410, -1787,
     410, -1787,   410,   410,   649,   410,  1811, -1787,   649,   649,
     410,   221,   276,  1813,  1822,  1829,  1831,  1839, -1787,  1844,
   -1787, -1787,  1851, -1787, -1787, -1787, -1787,  1852,  1854,  1855,
    1856,  1863,  1864,  1872,   410,   410,  1875,   845,  1877,  1878,
    1884,  1886, -1787,   649, -1787, -1787,   410, -1787,   410,   649,
     649,  1906,   649,   314,  1887,   410, -1787, -1787,   410,   410,
     410,   410,   410,   410, -1787, -1787, -1787, -1787, -1787,   410,
   -1787,   410, -1787,   410, -1787, -1787,  1687, -1787, -1787,   649,
    1890,  1891,   649, -1787,   649,  1894, -1787,  1895, -1787,  1903,
   -1787, -1787,  1907,   853,  1911,   649,  1916,   649, -1787,   649,
   -1787,  1917, -1787, -1787,   649, -1787,   410,   410,   680,   410,
     649,   410,  1930,   649,  1934,   649,   649,   410,   410,   410,
     410,   410, -1787,  1661,  1675,   649,  1937,   410,   649,   649,
     649,   410, -1787,   410,   649,    21,   410,   649,   649,   649,
     649,   410,   410,   410,   410,   410,   410,   410,   859,   649,
     649,   649,   649,    -2, -1787, -1787,  1946, -1787,  1959, -1787,
    1678,   649,   410,  1896,    52,   410,   861,  -134,  -119,   498,
    1870,  1963,   420,  1966,  1762,  1968,   436,  1969,  1973,  1976,
    1768, -1787, -1787,  1984,  1986,   649,   410,  1990, -1787,  1992,
     711,  1994,  1998, -1787, -1787,  2003,  2004, -1787,  2005,  2007,
    2008, -1787,   649, -1787, -1787,   888,  2015,  2023,   481,   649,
     808,   547,  2025,  2032,  2034,  2041, -1787,  2044,  2045,   410,
     649, -1787,  2046, -1787,   912,  2047, -1787,  2048,  2049,  2051,
    2052,  2053,   649,   649,   410, -1787,   649,  2055,  2057,  2058,
   -1787,  2069,  2070,  2071,  2073,  2076,  2079,  2080,  2081,  2082,
    2083,  2084,  2085,  2087, -1787,  2088,   921,   956,  2094,   958,
    2095,  2096, -1787,   384,   960, -1787,  2099,   410,   410,   410,
     410,   649,   410,   410,   345,   649, -1787,   649,   649,  1769,
   -1787,  2100,  2101, -1787,  1777, -1787,  2103, -1787, -1787,   989,
    1792, -1787, -1787,  2104,   649,  2147, -1787,   649,   649,   649,
     649,   649,   649,   649,   649,   649,   649,   649,   649,   649,
     649,  1794, -1787,   410,  2215, -1787, -1787, -1787, -1787,   649,
    2215, -1787, -1787,  2112, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787,  2223, -1787,   410, -1787, -1787, -1787,  2224, -1787,
     410, -1787, -1787, -1787, -1787, -1787, -1787,   410, -1787, -1787,
     -61,   410, -1787,   994,   649,  2116, -1787, -1787,  2117, -1787,
    2118,   410,  2119,  2120,   649, -1787,  2121,  2122,   410, -1787,
     649,   649, -1787, -1787,   649, -1787,  2123, -1787, -1787,   649,
   -1787,   649,   410,  2126, -1787,  1866, -1787,   649,   410, -1787,
    1906, -1787, -1787,   410, -1787,   449, -1787,   410,   649, -1787,
    2127,   649, -1787,   649,   997,   649,   649, -1787,   649,   410,
   -1787,   561, -1787,   410, -1787,   -10,  2128, -1787, -1787,  1885,
     410,   649,   649, -1787, -1787, -1787,  2129, -1787,  2072,  2131,
     410,   649,  2133, -1787,  1906, -1787,  2134,   410,   649,  2135,
     410,   649, -1787,  1002,   649,   649,   649,  2136,   410,   649,
     649,   649,   649,  1980, -1787,   121, -1787,   649,  -164,  -159,
   -1787, -1787, -1787,  2138, -1787,   410,   410,  1901,   649,   410,
   -1787,   649, -1787,  2139,   -62,  1918,  2140, -1787,  2142, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787,  2143,  2144, -1787,   410,   649, -1787,
   -1787, -1787, -1787,   649,   649,  2145,   649,   649,   649,  2146,
    2148,  2149, -1787, -1787, -1787,   581,  2257,   621,  2260, -1787,
    2262, -1787, -1787, -1787, -1787, -1787, -1787, -1787,   649,  2155,
   -1787, -1787, -1787, -1787, -1787,   649, -1787,   649, -1787,   649,
     410, -1787,   649,   -89,  2157, -1787,   506,   410,  2160,     1,
   -1787,  2162, -1787,   649,   649,  1919, -1787,  1926,  1927,  1928,
    1929, -1787,  2163, -1787,  2164,   410, -1787,   410,   649,   649,
     649,   410,   -20,   649,   649,   649, -1787,   649,   -52,   649,
     649,   649,  2168,   649,   649,  1962, -1787,  1964,   210,  1088,
    1978,  1095,   227, -1787,  1098,   649,   649,   649,   125,  2172,
   -1787,    30,   -26, -1787, -1787, -1787,  2173,   199,   649,   410,
     222, -1787,  2175,   649, -1787,   649, -1787,  2176,   410,  -170,
   -1787,  2177,  2178, -1787, -1787,  2179,  2180, -1787, -1787,  2181,
   -1787, -1787,  2182,  2183, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787,  2184,  2185, -1787, -1787,  2186, -1787,  1991, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787,  1993, -1787,  1999, -1787, -1787,
   -1787,  2187,  2188,  2189,    76,   410, -1787,   410, -1787,   226,
   -1787, -1787, -1787, -1787, -1787, -1787,  2190,  2191, -1787, -1787,
    2192, -1787, -1787, -1787, -1787, -1787, -1787,  2193,  2194,  2195,
    2197, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,  2198,
   -1787,  2199, -1787, -1787,  1101, -1787, -1787,  2200, -1787, -1787,
    2204, -1787,  2206,  2207,  2208,  2209,  2210,  2211,  2213,  2214,
    2216,  2217,  2220,  2000,   649, -1787,  2222, -1787, -1787, -1787,
    1139, -1787, -1787,  1140, -1787,  2230, -1787,   649,  2231,  1143,
     649,   649,  1155,  1160,   649,   649,   649,   649,   649,   649,
     649,  1174,   649,   649, -1787,  2232, -1787, -1787,   649, -1787,
   -1787,  2325,   410,  2333,   410, -1787,   410,    53, -1787,  2234,
    2235, -1787, -1787, -1787,   649,   649,   649, -1787, -1787,  2236,
     649,   649, -1787, -1787, -1787,   649,   410,   410, -1787,  2238,
    1183,  1906,  2241, -1787,   649,   649, -1787, -1787,   649, -1787,
    2242,  2243,  2244,  2245, -1787,   649, -1787,   410,   480, -1787,
    2037, -1787,  2246,   410, -1787,  2247,   649,  2248,  2249, -1787,
    1906,  2250,   649,  2251,  2252,   649,  2255,  2256, -1787,   410,
     410,   410,   410, -1787,  2259,   649,   649,   649,   649,   649,
   -1787,  2261,  1226, -1787,   649, -1787,   649, -1787,  2001,  2264,
   -1787,  2265,  2266, -1787,  2268, -1787,  2271, -1787, -1787,   410,
   -1787, -1787, -1787, -1787,   410,   410,   649,   649,   649, -1787,
     649,   649,  2272, -1787, -1787, -1787,   410, -1787,   410,   410,
   -1787,   410,   410,   649, -1787,   410,  2273,   410,  1227,  2274,
     410,   649, -1787,  2275, -1787,  2276, -1787,   410,  2277, -1787,
    2278,  2279, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
    1247,  1250,   649,   649,   410,   649,  1251,   649,  1254,  1279,
     649,   649, -1787,   649,   649,   649,   649, -1787,   649,   649,
   -1787, -1787, -1787, -1787,   490,   476, -1787,   649, -1787, -1787,
     649, -1787,   496,   493, -1787,  2006,   649,   649,   649,  2280,
   -1787, -1787,   649,   160, -1787,   649,   122, -1787, -1787,  2281,
     248,   649, -1787,  2282, -1787,  1303,  2283, -1787,  -151, -1787,
    2284, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787,  2288, -1787,  1311, -1787,  2030, -1787, -1787, -1787,
     410, -1787,  2290,  2291,   410, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787,  2292, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787,   410,  2293, -1787, -1787,  1370, -1787,   649, -1787,   649,
   -1787, -1787,   649,   649, -1787,   649, -1787,   649,   649,   649,
     649,   649,   649,   649,   649, -1787,   649,   649,   649, -1787,
   -1787, -1787, -1787, -1787,   211,   410,   649, -1787, -1787,   649,
     649,  2294,  2334,   649, -1787,  2296,   410, -1787, -1787,  2300,
     649,  2306, -1787,  1376,  2307,  2308, -1787, -1787, -1787, -1787,
    2312,   244, -1787, -1787,   410, -1787, -1787, -1787,   649, -1787,
     649, -1787,  2246, -1787,  2246,  2253,   649,   649, -1787,   649,
    2314, -1787,   649,   649, -1787,   649,   649, -1787,   410,  2315,
    2316,  2254, -1787,   410,   410,   410,   410,  2317, -1787, -1787,
    2320,  1390,  1407, -1787, -1787, -1787, -1787, -1787, -1787,  2322,
     410,   649,   649,   649,  2323,   649,   649, -1787, -1787,   694,
   -1787,   730, -1787, -1787,   410, -1787,  2324, -1787,  2327, -1787,
     410,  2329, -1787, -1787,   649, -1787, -1787, -1787, -1787,   649,
   -1787,   649,   649,   649,   410,   507, -1787,  1413,   649, -1787,
    2330, -1787,  1431,   649,  1435,  1436,   649,   649,   649,   649,
     649, -1787, -1787,   535,   526,  2331,  2335, -1787, -1787,   541,
     562, -1787,  2033,   649,   649,  2336, -1787,  2337,  2338, -1787,
    -137, -1787,   649,   649,  2339, -1787, -1787,  2342,   290, -1787,
   -1787,   649, -1787, -1787,  2344, -1787, -1787, -1787,  2035, -1787,
    2040,  2345, -1787, -1787,  2347, -1787,  2348, -1787, -1787,  2043,
    2350,    -1,   649,   649,  1443,  1445,   649,  1455,  2351,  1473,
     649,   649,   649,  1474,   649,   649,   410,   649,   410,  2354,
     649,   649, -1787,  2358, -1787, -1787, -1787,  2360, -1787, -1787,
    2361, -1787, -1787, -1787, -1787,   608, -1787, -1787, -1787, -1787,
   -1787,  2362,  2363,   649,  2367, -1787,  2368,   649,  2369,   649,
    2371, -1787, -1787,  1475,   410,  2373,  2374,  2269, -1787, -1787,
   -1787,  2376, -1787,  2377, -1787,   649,  2379,   649,   649, -1787,
     649,  1480,   410, -1787,   410, -1787,  2381, -1787, -1787,   649,
   -1787,  2382,  1504,   649,   649,   649,   649,  2383, -1787, -1787,
     649,   649, -1787, -1787,  2384,   649, -1787,  2385, -1787,  1512,
     649,   649,   649,   649,   649, -1787, -1787,   628,   392, -1787,
   -1787, -1787, -1787,   631,   470, -1787,  2386,   649,   649, -1787,
   -1787, -1787, -1787,   649,  2388,   162, -1787, -1787, -1787,  2390,
    2392, -1787, -1787,  2393, -1787,  2394, -1787, -1787, -1787, -1787,
    2396, -1787, -1787,   410,   649,   651,   649,   649, -1787,   649,
   -1787,   649,  1540, -1787,  1542, -1787, -1787,  2397,   649,   649,
    1546, -1787,   649,   649,  1550,   410,  2399,   649, -1787,  2400,
    2405,  2435, -1787, -1787, -1787, -1787, -1787,  2406, -1787, -1787,
    2407, -1787,  2408, -1787, -1787,   410,  2412, -1787, -1787,  1591,
   -1787, -1787,  2413, -1787,   649,   649,   649, -1787,   649, -1787,
   -1787, -1787,  2414, -1787, -1787,   649,  2415,  2416,  2418,   735,
   -1787,   649,   -19, -1787,   649, -1787, -1787,  2419,   649,   649,
     649,   649,   649, -1787, -1787,   657,   649, -1787, -1787,   663,
     649, -1787,   649,   649,   230, -1787,  2420, -1787, -1787, -1787,
   -1787, -1787, -1787,  2421,  -148, -1787,  1606,  1607,  2422,  2423,
   -1787,   649, -1787,  1634, -1787,   649,   649, -1787,  2425,  2426,
    1636, -1787,   649,   649, -1787,  2427, -1787, -1787,  2429, -1787,
   -1787, -1787,   410, -1787, -1787,   410, -1787,  2430,   649,   649,
    1652, -1787,   649, -1787, -1787, -1787,  2431, -1787,   649, -1787,
     649,   -14,    -9, -1787,   649,   649,   649,   649,  1658, -1787,
     348, -1787,   585,   649,   649,  2432, -1787, -1787, -1787, -1787,
     410, -1787,  1685, -1787,  1686, -1787, -1787,  2433, -1787,  1713,
    2438,  1730, -1787, -1787, -1787,  1731,  2440,  2443, -1787,  2444,
     410,   410,   410,   649,   649, -1787,   649,  2445, -1787,   649,
    2446, -1787,   649,  -129, -1787,   649,    29,   410,   649,   649,
     649, -1787,   649, -1787,   670, -1787,   677,   618,   649,   649,
   -1787,   371, -1787,   649, -1787,   649, -1787, -1787,  1780, -1787,
   -1787,  1784, -1787,   649, -1787, -1787,  2447,   410,   410,   410,
     649,  2449,  1797, -1787,   649, -1787,  2450, -1787,   649,  2451,
   -1787,   649,   219,  2452,   649,   649,   649,   649, -1787, -1787,
   -1787,   678,   649,   649,  1798, -1787,   686,  2453,  2455, -1787,
     649, -1787,  2456,  2458, -1787,  2459,   410,  1801, -1787, -1787,
    2460,  1828, -1787,  2461, -1787,  2462, -1787,   649, -1787,   649,
     649,  2463,   649, -1787,   701,  1888, -1787,   649, -1787, -1787,
   -1787,   649, -1787, -1787, -1787,  2465, -1787,   649, -1787, -1787,
    1892, -1787, -1787,  2466,   410,  2467, -1787,   649, -1787,   700,
     413, -1787,   649,   649,   649, -1787,  1893, -1787,  2468, -1787,
     410, -1787,   642, -1787, -1787,   709,   649,  2472,   649, -1787,
    2474, -1787,   410,   649, -1787,   649, -1787,  2475, -1787,   649,
   -1787,  2476,  2477,    83, -1787,   649, -1787, -1787, -1787,   649,
     649,   649,   649,   649,   649,   649,   649,  2478,  2480,   649,
     649, -1787, -1787,   106,   649, -1787,   649,   649,   649,   649,
     649,   649,   410,   649,  2481,  2482,   649,  2484, -1787, -1787,
     649, -1787,   410,   410,   649,   649,   649,  2486, -1787
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1787, -1787, -1787,  2303, -1787, -1787, -1787, -1787, -1787,  2355,
   -1787, -1787,  2409, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787,   142, -1787, -1787, -1787, -1787,  2137,
    2141, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787,  2353, -1787,  2115, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1786, -1787,   -84,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,  1732,
   -1787, -1787,  1724, -1787, -1787, -1787, -1787, -1787, -1787,  1830,
   -1787,  -152, -1089, -1787, -1787,  1770, -1787,  2341, -1787,  -128,
   -1787,  2328, -1787,  -220,  -829,  -329, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787,   712, -1113,  2309, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787,  1806,  -899, -1787, -1787,  1824,  -882,
   -1787,  -382, -1787,  2304, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787,  2227, -1787, -1787, -1787, -1787,  2301,  2267,  1599,  -221,
   -1452,  -794,  -870, -1787, -1787,   390, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787,  -332, -1787,    32,  -260, -1787, -1787,
    1382, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787, -1787, -1112,  1954, -1787, -1787, -1787,
   -1787, -1787, -1787, -1787,   728,  2436, -1787, -1787,  2751,    -6
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1122
static const yytype_int16 yytable[] =
{
     328,   329,   330,   333,  1557,   615,   840,  1853,  1561,  1349,
    1376,  1349,   883,  1376,  1402,  1386,   768,  1157,   788,  1547,
    1548,  1549,  1550,   363,   608,   636,  2093,  1177,   967,   968,
     969,   970,   363,   491,   324,   860,   320,   588,   609,  1124,
     320,   652,   324,   317,   408,  1920,  1580,  1969,   324,   320,
     554,   324,  1273,  1843,   363,   316,   324,   320,  1845,   320,
     325,   324,   589,   544,   364,   545,  2213,   411,   325,  2539,
     320,   324,  1365,   364,   325,   365,   366,   325,   537,   538,
    2352,   769,   325,  1606,   365,   366,  1057,   325,  2617,  1058,
     321,   634,   770,  1125,   321,   364,  1255,   325,   320,   324,
     324,   324,   861,   321,   546,  1853,   365,   366,  1365,  1853,
     320,   321,   771,   321,   470,   320,   635,   386,  1564,   476,
     477,   602,   798,   326,   321,   325,   325,   325,   522,  1126,
     424,   326,   772,   320,   327,  1059,  1921,   326,  1607,   320,
     326,   500,   327,  1230,   324,   326,   505,   555,   327,  1856,
     326,   327,   321,   324,  1589,  1857,   327,   320,  1046,   590,
     326,   327,   319,   320,   321,  1922,  1446,   323,  1178,   321,
     325,   327,   556,   320,  1739,   591,   324,   558,   559,   325,
     320,  1231,   492,   894,   576,   409,   773,   321,   326,   326,
     326,  1954,   774,   321,  2618,  1256,  1758,   775,  2519,   327,
     327,   327,   325,  2571,   320,  2353,  1970,  1806,  2574,  1565,
    1286,   321,   320,   766,   592,  1590,  2372,   321,  2540,  1127,
    1949,   363,   865,   442,   776,  2214,   777,   321,   593,   528,
    1890,   320,   778,   326,   321,  1844,   594,  1097,  1566,   799,
    1846,   557,   326,   363,   327,  1123,  2620,  1341,   595,   320,
     463,  1047,  1608,   327,   318,  2198,   320,  2466,   321,  1786,
     363,   785,   364,   320,  1060,   326,   321,   334,   320,  1601,
    1215,  1923,   889,   365,   366,   337,   327,   363,  1990,  1232,
     779,   320,   676,   679,   364,   321,  1766,   971,  1061,   320,
     531,  1062,   767,  1991,  1181,   365,   366,  1063,  1142,  1914,
    2728,   364,   800,   321,  2520,  1952,  2269,  1064,  2270,  2572,
     321,  1101,   365,   366,  2575,  1955,   801,   321,   364,  1155,
    1897,  1038,   321,  2745,   833,  2535,  1179,   836,  1466,   365,
     366,   842,   843,   547,  1853,   321,   614,   596,  1840,  2201,
    1591,   850,  1950,   321,  1567,   597,   852,   598,  1447,  1448,
    1449,   855,  2621,  2729,  1128,   610,  1807,   856,   857,   780,
     859,  1065,  1066,   320,   338,  2373,  1848,  2374,  1602,   320,
     320,   870,  2055,  1158,  2730,  1129,  2746,  2199,   599,  2467,
    1216,   363,   789,  1468,  1067,   879,   320,   781,   548,   549,
     550,   551,   552,  1095,   618,   891,   893,  2747,   790,   539,
     540,   541,   542,   543,   363,   367,  1082,  1841,   791,   905,
     341,   792,   908,   793,   794,   911,  1958,   321,   795,   620,
     342,  1501,   364,   321,   321,   363,   651,  1933,  1994,   943,
    1122,   681,   802,   365,   366,   803,  2656,   804,  1467,  1962,
     321,   345,   561,  1995,  1941,   364,   363,  2536,  1719,   805,
     806,   807,   992,   673,  1043,  1068,   365,   366,   808,   320,
     809,  2264,   320,  2202,  2203,  2206,   364,   951,   795,   810,
     811,   812,   813,   580,   958,  1048,   814,   365,   366,   815,
     816,   320,  1287,  1288,  1289,  1290,  1291,   364,  1292,  1293,
    1294,  1295,  1296,  1469,  1297,  1298,  1299,  1300,   365,   366,
    1804,  1793,  1378,   363,   320,  1388,   320,  2358,  1376,   363,
     581,  1163,  1351,   321,   861,  1959,   321,   320,   346,   324,
    1306,   348,  1349,   996,   998,   999,   363,   320,   798,   351,
    2246,  1502,  1720,  1050,  1052,   321,   352,  1707,  1963,   944,
    1257,  1258,  2657,   353,   364,   325,   324,  1017,   354,  1164,
     364,   320,  1022,  1023,  1026,   365,   366,   355,   321,   363,
     321,   365,   366,   324,  2207,  2583,  1039,   364,  1040,   320,
    1042,   321,   325,   562,   563,   564,   565,   862,   365,   366,
    1069,   321,  1072,  1074,  1076,   320,   320,  1072,  2635,   325,
    1204,   356,  1085,  1228,   357,   363,   324,  1087,  1259,  1089,
     364,  1708,   320,  1339,   358,   321,  2359,  1094,   677,  2454,
     603,   365,   366,  1102,  1853,  1104,  1105,  1106,   363,   327,
     582,   359,   325,   321,   583,   584,   585,   861,   360,   320,
    2704,   324,   324,   674,   861,   326,   364,  1614,  1115,   321,
     321,  1118,  1119,  1120,   861,   799,   327,   365,   366,  1130,
     861,   363,   326,  1621,  1388,   324,   321,   325,   325,   364,
    1138,   975,   321,   327,  1140,   604,  1389,  1260,   324,  1146,
     365,   366,   667,   668,  2066,  1615,  1154,  1156,  1339,   605,
    1261,   325,  1339,   321,  1165,   326,  1167,  2458,   324,   861,
    1363,  1622,   364,  2182,   325,   861,   327,  2082,  1650,  1186,
    1218,  1188,  1189,   365,   366,   361,   320,  2181,   800,   362,
    2188,  2713,   324,  2187,   325,  1610,  1876,  1893,  2317,   324,
     677,   326,   801,  1894,  2318,  1221,  1222,  1223,  1339,  1226,
     976,   327,   327,  1233,   363,  1616,  1651,  1237,   325,   528,
     320,  1339,  1339,  2336,   326,   325,  1246,  1247,  1248,  1657,
     324,  1623,  2335,  1611,  1254,   327,  1879,   326,  2341,   369,
     321,  1265,   861,  1267,  1658,   977,  1271,  1272,   327,   324,
     372,   324,  1279,   373,   321,   364,   325,   326,  1404,  2342,
     320,   324,   861,   324,   374,   861,   365,   366,   327,   529,
     324,   855,   530,  1568,   321,   325,  1652,   325,  1320,   375,
     531,   326,  2585,  1325,   324,   861,   331,   325,   326,   325,
     532,   861,   327,  1635,   855,  1389,   325,   861,   376,   327,
     324,   324,  1343,  1344,   861,  2404,   377,  1345,   384,  2302,
     325,   861,   861,  1350,   321,  2630,  1352,  1353,   802,   326,
     861,   803,   385,   804,  1356,  2453,   325,   325,  2457,  1360,
     327,  1361,  1362,   320,   861,   805,   806,   807,   326,  2714,
     326,   387,  1366,   861,   808,  2304,   809,  1371,  2475,   327,
     326,   327,   326,  1375,  2529,   810,   811,   812,   813,   326,
    2531,   327,   814,   327,   389,   815,   816,  2628,   390,   320,
     327,  1534,  1390,   326,  2629,  2663,   324,  1535,   324,  1394,
    1395,   324,  1396,  2668,   327,  1398,   324,   321,  1405,   326,
     326,   486,   391,   324,   324,   324,   332,  2703,  2688,   320,
     327,   327,   325,   324,   325,   320,  2716,   325,  1636,   324,
     410,   324,   325,  1420,   392,   320,   528,   849,   443,   325,
     325,   325,  1421,   321,  1196,  1422,  2516,  1262,  1263,   325,
     324,   997,  2517,   393,  1428,   325,   396,   325,   324,  1431,
     399,  1433,  1434,  1435,  1436,  1437,  1438,  1073,  1225,  1411,
    1412,  1445,   400,   321,   401,   402,   325,   335,   403,   321,
    1453,   404,   324,   444,   325,   326,   529,   326,  1458,   321,
     326,   324,  1463,  1464,   405,   326,   327,   531,   327,   406,
     407,   327,   326,   326,   326,  1264,   327,   532,   325,  1654,
    1655,   320,   326,   327,   327,   327,   412,   325,   326,   413,
     326,  1488,   415,   327,   446,  1656,   324,  1493,   324,   327,
     324,   327,   418,  1496,  1497,  1499,  1500,   419,   421,   326,
     818,   819,   821,  1270,   422,  1342,   423,   326,  1355,   320,
     327,   427,   325,  1374,   325,   430,   325,   320,   327,   324,
    1404,  1452,  1487,  1515,   324,   321,  1518,   324,  1519,   447,
    1524,   326,   324,  1808,   431,   437,  1583,  1525,  1604,  1527,
     326,  1529,   327,  1530,  1427,   325,   320,   438,  1532,   439,
     325,   327,  1536,   325,  1538,   441,  1202,  1541,   325,  1543,
    1544,   324,   449,   321,   450,  1646,   320,   349,   453,  1555,
     454,   321,  1558,  1559,  1560,   326,   464,   326,  1563,   326,
     899,  1571,  1572,  1573,  1574,   903,   327,   325,   327,  1669,
     327,  1797,  1584,  1585,  1586,  1587,  1588, -1121,  1698,   465,
     321,   466,   467,  1200,   370,  1597,   468,  1600,   326,   320,
    1605,  1797,  1797,   326,   469,   471,   326,   474,   324,   327,
     321,   326,   475, -1121,   327,   324, -1121,   327,   324,  1631,
     320,   324,   327,  1700,  1637,  1703,   478,  1709,   479,   480,
     320,   949,   481,   482,   325,   483,  1645,   484,   485,  1647,
     326,   325,   487,  1653,   325,   488,   489,   325,   490,   962,
     963,   327,   493,   321,  1667,   496,  1732,   380,  1670,   324,
     324,  1768,   497,   324,  1799,   503,  1677,  1678,   320,  1828,
    1680,   504,   506,   507,   321,   324, -1121,   508,   382,   509,
     324,   510,   511,   512,   321,   325,   325, -1121,   394,   325,
    1699,  1701,   513,  1704,   324,   514,  1934,   326,  1710,   515,
    1942,   325,   516,   324,   326,  1716,   325,   326,   327,  1722,
     326,  1723,  1724,   517,   518,   327,   519,   520,   327,   521,
     325,   327,   321,  1733,   320,   523,   397,   320,  1737,   325,
   -1121,  1740,  1741,  1742,  1743,  1744,  1745,  1746,  1747,  1748,
    1749,  1750,  1751,  1752,  1753,   577,   324,   324,   326,   326,
     526,   587,   326,  1759,   578,  1936,   646,   645,   320,   327,
     327,   650,  1939,   327,   326,  1944,   783,   324,  2005,   326,
     324,   324,   325,   325,   324,   327,   320,   828,   321,   822,
     327,   321,   416,   326,   829,   425,   827,  1769,  1770,   830,
     320,   831,   326,   325,   327,   832,   325,   325,  1776,   324,
     325,   841,   844,   327,  1780,  1781,  2024,  2026,  1782,  1430,
    2031,   845,   321,  1371,   847,  1785,   428,  1439,  1440,  1441,
    1442,  1790,  2034,   324,  1792,   325,   324,  2036,   848,  1390,
     321,   324,  1795,   851,   432,   326,   326,  1798,  1800,  1801,
    1802,  2045,  1803,   320,   321,  1405,   327,   327,   434,   325,
    2068,   858,   325,   320,   863,   320,   326,   325,   853,   326,
     326,   873,   874,   326,   875,  1818,   876,   327,  1821,   877,
     327,   327,  1824,   878,   327,  1827,   320,  1829,  1830,  1831,
    1832,   320,   884,  1835,  1836,  1837,  1838,   947,   326,   320,
     324,  1842,   888,  2119,  2147,   900,   324,   321,   901,   327,
     320,   455,  1852,   320,   904,  1854,   320,   321,   907,   321,
     324,   459,   326,   461,  2158,   326,   325,  2160,  2166,   946,
     326,  2169,   325,   327,   910,   915,   327,   324,   916,   917,
     321,   327,  1866,   324,   472,   321,   325,  1867,  1868,   494,
    1870,  1871,  1872,   321,   919,   920,  2171,   834,   921,   924,
     925,   324,   928,   325,   321,   324,   324,   321,   871,   325,
     321,   926,  1883,   324,   929,   324,   931,   933,   959,  1885,
    2210,  1886,   934,  1887,   936,   324,  1889,   325,  2217,   326,
     320,   325,   325,  1898,   937,   326,   938,  1900,  1901,   325,
     327,   325,   939,   324,   324,   324,   327,   942,   950,   326,
     324,   325,  1911,  1912,  1913,   955,   956,  1916,  1917,  1918,
     327,  1919,   961,  1924,  1925,  1926,   326,  1928,  1929,   325,
     325,   325,   326,  1937,   324,  1940,   325,   327,  1945,  1946,
    1947,  1948,   324,   327,   321,  1953,   964,  2228,   973,   979,
     326,   980,  1960,  2259,   326,   326,   983,  1965,   991,  1966,
     325,   327,   326,  2183,   326,   327,   327,  2290,   325,   995,
     324,  2189,   324,   327,   326,   327,   324,   363,   994,  1008,
     324,   320,  1000,  1003,  2292,   327,   320,   320,  1001,  1002,
    2319,  1004,   326,   326,   326,   320,   325,  1005,   325,   326,
    1010,  1006,   325,   327,   327,   327,   325,  1011,  2323,  1013,
     327,  1007,  2326,  2328,  1009,  1014,  1015,   320,   364,   320,
    2378,   324,  2380,   326,  1016,   320,  1045,  1054,  1055,   365,
     366,   326,  2383,  1056,   327,   321,   324,   324,  1077,   981,
     321,   321,   327,  1078,  1092,  1172,  1080,   325,  1081,   321,
    2386,  2391,  2414,  1191,  1084,  1088,  1086,  2427,  2006,   326,
    1090,   326,   325,   325,   324,   326,   324,  1091,  1099,   326,
     327,   321,   327,   321,  1112,  1194,   327,  1211,  2022,   321,
     327,  2434,   324,  1347,  2025,  1116,  1117,  2027,   324,  2446,
     325,  2029,   325,  1759,  2032,  2033,  2035,  2037,  2038,  2039,
    2040,  2041,  2042,  2043,  2044,  2046,  2047,  2048,   325,  2265,
     326,  1159,  1759,  1121,   325,   324,   324,  2480,  1144,  2482,
    1145,   327,  1147,  2487,  1161,   326,   326,  2491,  2059,  2060,
    2061,  1149,  1160,  1162,  2063,  2064,   327,   327,  1166,  2065,
    1168,   325,   325,   324,  2069,  2071,  1169,  1183,  2073,  2074,
    1184,  1187,  2075,   326,  1193,   326,  1198,  1174,  1190,  2080,
     324,   324,  1206,  1213,   327,   320,   327,  1227,  2504,   325,
    2096,   326,   320,  1244,  2100,  1240,  2102,   326,  1252,  2105,
     320,  1253,   327,  2541,  2543,  1266,   325,   325,   327,  2113,
    2114,  2115,  2116,  2117,   320,  1277,  2120,   320,  2121,  1278,
    2122,  1281,  1282,  1283,   326,   326,  1285,  1308,  1311,  1312,
     324,  2548,  2337,  2554,   324,   327,   327,  1321,  2343,   321,
    2132,  2133,  2134,  1413,  2135,  2136,   321,   324,   324,  2565,
    1415,   324,   326,  1314,   321,  2581,   325,  2143,  1551,  1315,
     325,  1316,  2148,   327,  1317,  2151,  1319,  1326,   321,   326,
     326,   321,  1553,   325,   325,  1595,  1322,   325,   324,  2375,
     327,   327,  2592,  2594,  2159,  2161,  2162,  2163,  1327,  2165,
    2167,  2168,  2170,  2172,  2173,  2174,  1328,  2175,  2176,  2177,
    2178,   320,  2179,  2180,   325,  1329,  1330,   320,   320,  2184,
    2597,  2185,  1331,  1332,  2186,  1333,   320,  2190,  1334,   326,
    2193,  2194,  2195,   326,  1335,  1336,  2197,  2600,  2602,  2200,
     327,   320,  1337,   320,   327,  2208,   326,   326,   324,  2211,
     326,  1338,   324,   324,  1346,  1354,   324,   327,   327,  1358,
     324,   327,  1359,  1367,  1368,   321,   324,  1372,  2218,  1618,
    1373,   321,   321,  1381,   325,  1627,  1725,   326,   325,   325,
     321,  1382,   325,  1397,  1729,  1599,   325,  2639,   327,  1152,
    1407,  2641,   325,  1514,  1408,   321,  2455,   321,  1417,  1734,
    1418,  1754,  2459,  1425,  2649,  2666,  1432,   324,  2676,  2229,
    1443,  2230,  1444,  2231,  1450,   320,  2232,  2233,  1462,  2234,
    1470,  2235,  2236,  2237,  2238,  2239,  2240,  2241,  2242,  1471,
    2243,  2244,  2245,   325,   320,  2679,  1472,   326,  1473,   324,
    2249,   326,   326,  2250,  2251,   326,  1474,  2254,   327,   677,
     320,  1475,   327,   327,  2257,   677,   327,  2260,  1476,  1477,
     327,  1478,  1479,  1480,  1139,   325,   327,   320,   320,   321,
    1481,  1482,  2267,  1788,  2268,   320,   320,   320,   320,  1483,
    2272,  2273,  1486,  2274,  1489,  1490,  2276,  2277,   321,  2278,
    2279,  1491,  1810,  1492,  1503,  2691,   326,  1516,  1517,  2697,
    2709,  1520,  1521,   321,   321,  2291,  2293,   327,  1850,   321,
    1522,   320,  1024,   320,  1523,  2296,  2297,  2298,  1526,  2300,
    2301,   321,   321,  1528,  1531,  1858,  1902,   320,   326,   321,
     321,   321,   321,  1904,  1905,  1906,  1907,  1540,  2311,   327,
     320,  1542,   320,  2312,  1556,  2313,  2314,  2315,   320,   320,
     320,  2320,  2321,  1593,  1612,   320,  2324,  2325,  2327,  2329,
    2330,  2331,  2332,  2333,  2334,   321,  1594,   321,  2338,  1930,
    1613,  1932,  1025,  1617,  2344,  1620,  1624,  2347,  2348,   320,
    1625,   321,   320,  1626,   320,  1938,  2354,  2355,  2584,   320,
    2586,  1629,   320,  1630,   321,  2360,   321,  1633,  1981,  1634,
    1983,  1638,   321,   321,   321,  1639,  1985,  2020,  2123,   321,
    1640,  1641,  1642,  2191,  1643,  1644,  2376,  2377,  2379,  2381,
    2382,  2384,  1648,  2387,  2388,  2389,  2390,  2392,  2393,  2394,
    1649,  2396,  1660,   321,  2399,  2400,   321,  2219,   321,  1661,
    2345,  1662,  2362,   321,  1738,  2631,   321,  2364,  1663,  2636,
    2369,  1664,  1665,  1668,  1671,  1672,  1673,  2407,  1674,  1675,
    1676,  2410,  1681,  2412,  1682,  1683,  2083,  2415,  2084,  2085,
    2086,  2087,  2088,  2089,  2090,  2091,  1684,  1685,  1686,  2422,
    1687,  2424,  2425,  1688,  2426,  2428,  1689,  1690,  1691,  1692,
    1693,  1694,  1695,  2432,  1696,  1697,  2435,  2436,  2437,  2438,
    2439,  1702,  1705,  1706,  2441,  2442,  1711,  1727,  1728,  2444,
    1731,  1736,  1757,  2447,  2448,  2449,  2450,  2451,  2452,  1760,
    1761,  1763,  2689,  1771,  1772,  1773,  1774,  1775,  1777,  1778,
    1783,  2462,  2463,  1787,  1796,  1809,  1814,  2464,  1816,  1815,
    1819,  1822,  1825,  1833,  1839,  1847,  1855,  1860,  2705,  1861,
    1862,  1863,  1869,  1873,  1878,  1874,  1875,  1881,  2474,  1882,
    2476,  2477,  1884,  2478,  1892,  2479,  2481,  1896,  2483,  1899,
    1908,  1909,  2485,  2486,  2488,  1927,  2489,  2490,  2492,  1951,
    1957,  2495,  1964,  1967,  1971,  1972,  1973,  1974,  1975,  1976,
    1977,  1978,  1979,  1980,  1987,  1988,  1989,  1996,  1997,  1998,
    1999,  2000,  2001,  2505,  2002,  2003,  2004,  2007,  2507,  2508,
    2509,  2008,  2510,  2009,  2010,  2011,  2012,  2013,  2014,  2512,
    2015,  2016,  2050,  2017,  2018,  2518,  2521,  2019,  2522,  2023,
    2052,  2253,  2524,  2525,  2526,  2527,  2528,  2028,  2030,  2049,
    2530,  2057,  2058,  2062,  2532,  2067,  2533,  2534,  2072,  2076,
    2077,  2078,  2079,  1810,  2095,  2097,  2098,  2101,  2103,  2104,
    2542,  2544,  2106,  2107,   525,  2547,  2112,  2549,  2118,  2550,
    2551,  2124,  2125,  2126,  2555,  2127,  2556,  2557,  2128,  2137,
    2145,  2149,  2152,  2153,  2155,  2156,  2157,  2196,  2205,  2209,
    2212,  2215,  2563,  2564,  2566,  2216,  2567,  2222,  2223,  2225,
    2227,  2252,  2569,  2255,  2570,  2573,  2576,  2256,  2577,  2578,
    2579,  2580,  2582,  2258,  2261,  2262,  2587,  2588,  2589,  2263,
    2271,  2275,  2281,  2282,  2288,  2283,  2593,  2289,  2595,  2294,
    2299,  2307,  2498,  2598,  2308,  2601,  2310,  2322,  2339,  2603,
    2419,  2606,  2340,  2349,  2350,  2351,  2356,  2610,  2611,  2357,
    2612,  2361,  2366,  2614,  2367,  2368,  2616,  2371,  2385,  2619,
    2622,  2398,  2624,  2625,  2626,  2401,  2627,  2402,  2403,  2405,
    2406,  2632,  2633,  2634,  2408,  2409,  2411,  2637,  2413,  2638,
    2417,  2418,  2640,  2420,  2421,  2642,  2423,  2643,  2431,  2433,
    2440,  2443,  2445,  2461,  2647,  2465,  2650,  2468,  2651,  2469,
    2470,  2471,  2653,  2472,  2484,  2655,  2494,  2496,  2659,  2660,
    2661,  2662,  2497,  2499,  2500,  2501,  2664,  2665,  2667,  2503,
    2506,  2511,  2513,  2514,  2671,  2515,  2523,  2537,  2538,  2545,
    2546,  2677,  2552,  2553,  2558,  2680,  2559,  2562,  2568,  2590,
    2596,  2683,   764,  2684,  2685,  2599,  2687,  2604,  2690,  2692,
    2605,  2693,  2613,  2615,  2644,  2694,  2648,  2652,  2654,  2658,
    2669,  2696,  2670,  2672,  2698,  2673,  2674,  2678,  2681,  2682,
    2686,  2702,  2695,  2699,  2701,  2711,  2706,  2707,  2708,  2718,
    2710,  2720,  2724,  2726,  2727,  2741,  2715,  2742,  2758,  2759,
    2717,  2761,  2719,  2768,  2609,  1034,   765,  2722,  1079,  2723,
     854,  1037,  1454,  2725,   864,  1379,  1451,  2731,  1399,  2732,
    1383,   902,  1423,  2733,  2734,  2735,  2736,  2737,  2738,  2739,
    2740,  1577,   895,  2743,  2744,   918,   987,  2748,  2749,   922,
    2750,  2751,  2752,  2753,  2754,  2755,   957,  2757,  1310,   820,
    2760,  1784,     0,     0,  2762,   322,     0,     0,  2765,  2766,
    2767,     0,   336,     0,     0,   339,   340,     0,     0,   343,
     344,     0,     0,   347,     0,   350,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     371,     0,     0,     0,     0,     0,     0,   378,   379,   381,
     383,     0,     0,     0,     0,   388,     0,     0,     0,     0,
       0,   395,     0,   398,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   414,     0,   417,
       0,     0,   420,     0,     0,     0,     0,   426,     0,   429,
       0,     0,   433,   435,   436,     0,     0,     0,   440,     0,
       0,   445,   448,     0,     0,   451,   452,     0,     0,   456,
     457,   458,   460,   462,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   473,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   495,     0,     0,   498,   499,     0,   501,
     502,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   527,     0,     0,
       0,     0,   553,     0,     0,     0,   560,     0,   572,   574,
     575,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   613,     0,   616,   617,     0,
       0,   619,   621,   622,   623,     0,   625,     0,   627,   629,
     631,   633,     0,     0,   625,   629,   633,     0,     0,   638,
     640,     0,     0,   642,   644,     0,     0,   647,   648,   649,
       0,     0,     0,   653,     0,     0,     0,     0,     0,     0,
       0,   655,   657,   658,     0,     0,   659,   660,   662,   663,
     664,     0,   666,   572,   572,   670,   671,   672,   675,     0,
     678,   680,   682,     0,   683,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   784,   786,   787,
       0,     0,     0,     0,     0,     0,     0,   826,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   835,     0,     0,   839,   839,     0,     0,     0,
       0,     0,   846,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   627,   625,   868,     0,     0,   872,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   882,   882,     0,   887,     0,     0,
       0,   629,     0,     0,   633,   631,     0,     0,   898,   826,
       0,     0,     0,   868,   826,     0,     0,   906,     0,     0,
     909,     0,     0,     0,   914,     0,     0,     0,     0,     0,
       0,     0,   655,     0,     0,     0,     0,     0,     0,   666,
     923,     0,     0,     0,   927,     0,     0,   930,     0,     0,
       0,   932,     0,     0,   935,     0,     0,     0,     0,     0,
       0,     0,   941,     0,   945,     0,     0,     0,     0,     0,
     826,     0,     0,     0,     0,   954,     0,     0,     0,     0,
     657,     0,     0,     0,     0,     0,     0,     0,   826,   826,
       0,     0,     0,   966,     0,     0,   974,     0,     0,     0,
     978,     0,     0,   982,     0,     0,     0,     0,     0,   986,
     662,   990,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   993,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1012,
       0,     0,     0,     0,     0,  1018,  1019,  1020,  1021,     0,
       0,     0,  1027,  1028,  1029,  1030,  1031,     0,  1033,  1033,
    1036,  1036,     0,     0,     0,     0,  1041,     0,  1044,     0,
       0,  1049,  1051,  1053,     0,     0,     0,     0,  1070,     0,
       0,     0,     0,     0,     0,     0,     0,  1083,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1093,     0,     0,  1096,     0,  1098,     0,  1100,
       0,  1103,     0,     0,     0,     0,  1107,     0,  1108,     0,
    1109,     0,  1110,     0,  1111,     0,     0,     0,     0,     0,
       0,  1113,     0,  1114,     0,     0,     0,     0,     0,     0,
       0,     0,  1096,     0,  1098,     0,     0,     0,  1131,  1132,
    1133,  1134,     0,  1135,  1136,  1137,     0,     0,     0,     0,
       0,     0,  1141,  1143,     0,     0,     0,     0,     0,  1148,
       0,  1150,  1151,  1153,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1171,  1173,     0,  1175,
    1176,  1180,  1182,     0,     0,  1185,     0,     0,     0,     0,
       0,  1192,     0,  1195,  1197,     0,  1199,  1201,  1203,  1205,
       0,  1207,  1208,  1209,  1210,  1212,     0,  1214,     0,  1217,
    1219,  1220,     0,     0,     0,  1224,     0,     0,  1229,     0,
       0,  1234,  1235,  1236,     0,  1238,  1239,     0,  1241,  1242,
    1243,     0,  1245,     0,     0,     0,  1249,  1250,  1251,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1268,  1269,     0,     0,  1275,  1276,     0,     0,     0,
    1280,     0,     0,     0,     0,  1284,     0,     0,  1301,  1302,
    1303,     0,  1305,  1307,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1318,     0,     0,     0,     0,  1323,  1324,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1340,     0,  1098,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1348,     0,   839,
       0,   839,   839,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1357,     0,     0,     0,     0,     0,     0,
       0,     0,  1364,     0,     0,     0,     0,   868,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   882,     0,  1377,   882,   882,   887,     0,  1380,     0,
       0,     0,     0,     0,     0,     0,     0,  1385,  1387,     0,
    1340,     0,     0,   868,  1340,     0,     0,     0,     0,     0,
       0,     0,     0,  1401,  1403,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1406,     0,     0,     0,     0,     0,
       0,     0,     0,  1409,     0,     0,  1410,     0,     0,     0,
       0,  1414,     0,     0,  1416,     0,     0,     0,  1419,     0,
    1340,     0,     0,   954,     0,  1424,     0,     0,     0,     0,
    1426,     0,     0,  1340,  1340,  1429,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   986,     0,     0,     0,   990,
       0,  1455,     0,  1456,  1457,     0,     0,     0,     0,     0,
       0,  1465,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1484,  1485,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1494,     0,  1495,
       0,     0,  1498,     0,     0,     0,  1504,     0,     0,  1505,
    1506,  1507,  1508,  1509,  1510,     0,     0,     0,     0,     0,
    1511,     0,  1512,     0,  1513,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1098,  1533,     0,
    1537,     0,  1539,     0,     0,     0,     0,     0,  1546,  1546,
    1546,  1546,  1546,     0,  1552,  1554,     0,     0,     0,     0,
       0,     0,     0,     0,  1562,     0,  1569,  1570,     0,     0,
       0,     0,  1576,  1576,  1578,  1579,  1546,  1581,  1582,     0,
       0,     0,     0,     0,  1592,     0,     0,     0,     0,     0,
       0,  1596,     0,  1598,     0,     0,  1603,     0,     0,  1609,
       0,     0,     0,     0,     0,  1619,     0,     0,     0,     0,
       0,  1628,     0,     0,     0,     0,     0,  1632,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1659,     0,     0,     0,     0,     0,     0,     0,
    1666,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1679,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1712,  1713,
    1714,  1715,     0,  1717,  1718,  1721,     0,     0,     0,     0,
    1726,     0,     0,     0,     0,  1730,     0,     0,     0,     0,
       0,  1735,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1755,     0,  1756,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1762,     0,     0,     0,     0,
       0,  1764,     0,     0,     0,     0,     0,     0,  1765,     0,
       0,  1098,  1767,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   839,     0,     0,     0,     0,     0,     0,  1779,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1789,     0,     0,   882,
       0,  1791,     0,     0,  1387,     0,  1794,     0,  1794,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1403,     0,  1805,     0,  1805,     0,     0,     0,     0,     0,
    1812,  1813,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1817,     0,     0,     0,  1820,     0,     0,  1823,     0,
       0,  1826,     0,     0,     0,     0,     0,     0,     0,  1834,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1546,  1849,  1851,     0,
       0,     0,     0,     0,     0,     0,  1859,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1865,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1877,     0,  1880,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1888,     0,     0,  1891,     0,     0,     0,  1895,     0,
       0,     0,     0,     0,     0,     0,  1903,     0,  1903,  1903,
    1903,  1903,     0,     0,     0,     0,  1910,     0,     0,     0,
       0,     0,     0,  1915,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1931,     0,  1931,  1935,
       0,  1903,     0,  1943,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1956,     0,     0,     0,     0,     0,     0,
    1961,     0,     0,     0,     0,     0,     0,     0,     0,  1968,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1982,     0,
       0,     0,     0,     0,     0,     0,  1984,     0,  1986,     0,
       0,     0,     0,     0,     0,     0,  1992,     0,  1993,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2021,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2051,     0,  2053,     0,  2054,  2056,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2070,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2081,     0,
       0,     0,     0,     0,  2094,     0,     0,     0,     0,     0,
       0,  2099,     0,     0,     0,     0,     0,     0,     0,     0,
    2108,  2109,  2110,  2111,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1903,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2129,     0,     0,     0,     0,  2130,  2131,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2138,     0,  2139,
    2140,     0,  2141,  2142,     0,     0,  2144,     0,  2146,     0,
       0,  2150,     0,     0,     0,     0,     0,     0,  2154,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2164,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2192,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2204,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2220,     0,     0,
       0,  2221,     0,     0,     0,  2224,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2226,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2247,  2248,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2266,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2280,
       0,     0,     0,     0,  2284,  2285,  2286,  2287,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2295,     0,     0,     0,     0,     0,     0,     0,     0,
    2303,     0,  2305,     0,     0,  2306,     0,     0,     0,     0,
       0,  2309,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2316,     0,     0,     0,     0,
     684,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   685,     0,     0,     0,
       0,   686,     0,  2346,     0,   687,     0,     0,     0,     0,
     688,     0,     0,     0,     0,     0,     0,     0,   689,     0,
       0,   690,   691,     0,   692,     0,   693,   694,     0,  2363,
     695,  2365,     0,     0,     0,     0,     0,     0,     0,     0,
    2370,   696,     0,     0,     0,     0,     0,     0,   697,   698,
     699,     0,     0,     0,     0,     0,     0,  2395,     0,  2397,
       0,     0,     0,     0,     0,     0,     0,     0,   700,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   701,   702,
     703,   704,     0,     0,     0,  2416,   705,     0,     0,   706,
     707,     0,     0,   708,   709,     0,     0,   710,     0,     0,
       0,     0,     0,  2429,     0,  2430,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   711,     0,   712,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     713,     0,     0,     0,     0,     0,     0,     0,     0,  2456,
       0,     0,   714,   715,     0,  2460,     0,     0,     0,     0,
       0,     0,   716,     0,   717,     0,   718,   719,     0,     0,
       0,   720,   721,   722,   723,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2473,     0,     0,     0,   724,   725,
       0,     0,     0,     0,   726,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2493,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     727,   728,     0,   729,     0,   730,  2502,   731,     0,     0,
       0,     0,   732,   733,   734,   735,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   736,
       0,     0,     0,     0,     0,     0,   737,     0,     0,     0,
       0,   738,     0,   739,   740,   741,   742,   743,   744,   745,
     746,   747,   748,   749,     0,     0,     0,     0,   750,   751,
       0,   752,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   753,     0,     0,     0,
       0,     0,     0,     0,     0,   754,   755,     0,     0,     0,
       0,     0,     0,  2560,     0,     0,  2561,     0,     0,   756,
       0,     0,   757,   758,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   759,     0,     0,   760,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2591,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2607,  2608,  1865,     0,     0,     0,     0,     0,   761,
       0,     0,     0,     0,     0,     0,     0,     0,  2623,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2645,  2646,
    2130,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2675,     0,     1,
       2,     3,     0,     0,     0,     4,     0,     5,     0,     0,
       0,     0,     6,     0,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,     0,    18,    19,     0,     0,
       0,    20,    21,    22,     0,  2700,     0,     0,     0,     0,
       0,    23,     0,    24,    25,    26,    27,     0,     0,     0,
       0,  2712,     0,     0,    28,     0,     0,     0,     0,     0,
      29,    30,     0,  2721,    31,     0,     0,    32,    33,    34,
       0,     0,     0,    35,     0,    36,     0,     0,     0,     0,
      37,     0,     0,    38,    39,    40,    41,   524,    42,    43,
      44,     0,     0,     0,     0,     0,    45,     0,     0,     0,
      46,    47,     0,  2756,    48,    49,     0,    50,    51,    52,
       0,     0,     0,  2763,  2764,     0,    53,     0,     0,     0,
       0,    54,     0,     0,    55,     0,    56,    57,     0,     0,
      58,     0,     0,     0,    59,     0,     0,    60,    61,    62,
      63,    64,    65,    66,    67,    68,    69,    70,    71,     0,
       0,     0,     0,    72,    73,     0,     0,    74,     0,    75,
      76,    77,    78,    79,     0,    80,     0,    81,     0,     0,
      82,     0,     0,     0,     0,    83,     0,    84,    85,    86,
      87,     0,     0,    88,    89,    90,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    91,
       0,     0,     0,     0,     0,     0,    92,     0,     0,     0,
       0,    93,     0,     0,     0,    94,     0,     0,     0,    95,
      96,    97,     0,     0,     0,    98,    99,   100,   101,     0,
     102,     0,   103,   104,   105,     0,   106,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   107,   108,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   109,
     110,     0,     0,     0,     0,   111,     0,   112,     0,   113,
       0,     0,   114,     0,   115,   116,     0,     0,     0,   117,
       0,     0,     0,     0,     0,     0,     0,     0,   118,   119,
       0,   120,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   121,     0,   122,     0,     0,     0,   123,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   124,     0,   125,     0,   126,     0,   127,   128,
     129,   130,   131,   132,     0,     0,   133,     0,     0,     0,
       0,     0,   134,     0,     0,   135,   136,   137,     0,     0,
       0,     0,     0,   138,     0,   139,     0,     0,     0,     0,
       0,   140,   141,     0,     0,     0,     0,     0,     0,   142,
       0,     0,     0,   143,   144,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   145,   146,   147,   148,     0,     0,     0,   149,     0,
       0,     0,     0,     0,     0,     0,     0,   150,     0,   151,
     152,   153,   154,   155,   156,   157,   158,     0,     0,     0,
       0,     0,   159,   160,     0,   161,   162,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   163,   164,   165,     0,   166,   167,     0,     0,
     168,     1,     2,     3,   169,     0,     0,     4,     0,     5,
       0,     0,     0,     0,     6,     0,     7,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,     0,    18,    19,
       0,     0,     0,    20,    21,    22,     0,     0,     0,     0,
       0,     0,     0,    23,     0,    24,    25,    26,    27,     0,
       0,     0,     0,     0,     0,     0,    28,     0,     0,     0,
       0,     0,    29,    30,     0,     0,    31,     0,     0,    32,
      33,    34,     0,     0,     0,    35,     0,    36,     0,     0,
       0,     0,    37,     0,     0,    38,    39,    40,    41,     0,
      42,    43,    44,     0,     0,     0,     0,     0,    45,     0,
       0,     0,    46,    47,     0,     0,    48,    49,     0,    50,
      51,    52,     0,     0,     0,     0,     0,     0,    53,     0,
       0,     0,     0,    54,     0,     0,    55,     0,    56,    57,
       0,     0,    58,     0,     0,     0,    59,     0,     0,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    70,
      71,     0,     0,     0,     0,    72,    73,     0,     0,    74,
       0,    75,    76,    77,    78,    79,     0,    80,     0,    81,
       0,     0,    82,     0,     0,     0,     0,    83,     0,    84,
      85,    86,    87,     0,     0,    88,    89,    90,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    91,     0,     0,     0,     0,     0,     0,    92,     0,
       0,     0,     0,    93,     0,     0,     0,    94,     0,     0,
       0,    95,    96,    97,     0,     0,     0,    98,    99,   100,
     101,     0,   102,     0,   103,   104,   105,     0,   106,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   107,
     108,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   109,   110,     0,     0,     0,     0,   111,     0,   112,
       0,   113,     0,     0,   114,     0,   115,   116,     0,     0,
       0,   117,     0,     0,     0,     0,     0,     0,     0,     0,
     118,   119,     0,   120,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   121,     0,   122,     0,     0,
       0,   123,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   124,     0,   125,     0,   126,     0,
     127,   128,   129,   130,   131,   132,     0,     0,   133,     0,
       0,     0,     0,     0,   134,     0,     0,   135,   136,   137,
       0,     0,     0,     0,     0,   138,     0,   139,     0,     0,
       0,     0,     0,   140,   141,     0,     0,     0,     0,     0,
       0,   142,     0,     0,     0,   143,   144,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   145,   146,   147,   148,     0,     0,     0,
     149,     0,     0,     0,     0,     0,     0,     0,     0,   150,
       0,   151,   152,   153,   154,   155,   156,   157,   158,     0,
       0,     0,     0,     0,   159,   160,     0,   161,   162,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   163,   164,   165,     0,   166,   167,
       0,     0,   168,     0,     0,     0,   169
};

static const yytype_int16 yycheck[] =
{
       6,     7,     8,     9,  1117,   226,   335,  1459,  1121,   838,
     880,   840,   394,   883,   913,   897,    10,    14,    41,  1108,
    1109,  1110,  1111,    33,    32,   245,  1812,    76,   115,   116,
     117,   118,    33,   217,    70,   367,   159,    17,    46,     8,
     159,   262,    70,     6,   217,    97,  1135,   217,    70,   159,
      17,    70,     7,   217,    33,   217,    70,   159,   217,   159,
      96,    70,    42,    78,    74,    80,   217,    73,    96,   217,
     159,    70,   866,    74,    96,    85,    86,    96,    39,    40,
     217,    75,    96,   217,    85,    86,    12,    96,   217,    15,
     213,   243,    86,    62,   213,    74,   217,    96,   159,    70,
      70,    70,   154,   213,   119,  1557,    85,    86,   902,  1561,
     159,   213,   106,   213,   120,   159,   244,   316,    97,   125,
     126,   205,    95,   159,   213,    96,    96,    96,   316,    98,
     107,   159,   126,   159,   170,    61,   188,   159,   272,   159,
     159,   147,   170,   217,    70,   159,   152,   114,   170,   211,
     159,   170,   213,    70,   156,   217,   170,   159,   217,   139,
     159,   170,   217,   159,   213,   217,   217,   217,   217,   213,
      96,   170,   139,   159,  1286,   155,    70,   183,   184,    96,
     159,   255,   366,   403,   190,   358,   180,   213,   159,   159,
     159,   217,   186,   213,   323,   316,  1308,   191,   217,   170,
     170,   170,    96,   217,   159,   342,   376,   217,   217,   188,
     107,   213,   159,   297,   194,   217,   217,   213,   366,   188,
      95,    33,   374,   107,   218,   376,   220,   213,   208,   209,
     319,   159,   226,   159,   213,   399,   216,   347,   217,   212,
     399,   208,   159,    33,   170,   347,   217,   347,   228,   159,
     107,   310,   371,   170,   217,    95,   159,    95,   213,  1372,
      33,   384,    74,   159,   190,   159,   213,   217,   159,   217,
     107,   323,   400,    85,    86,   217,   170,    33,   202,   353,
     274,   159,   288,   289,    74,   213,   347,   374,   214,   159,
     270,   217,   298,   217,   290,    85,    86,   223,   342,   319,
     217,    74,   275,   213,   323,   275,  2092,   233,  2094,   323,
     213,   347,    85,    86,   323,   341,   289,   213,    74,   347,
     319,   343,   213,   217,   330,    95,   375,   333,   107,    85,
      86,   337,   338,   348,  1786,   213,   322,   317,   217,   217,
     342,   347,   217,   213,   323,   325,   352,   327,   399,   400,
     401,   357,   323,   270,   323,   363,   366,   363,   364,   353,
     366,   287,   288,   159,   217,   366,  1455,   368,   316,   159,
     159,   377,   319,   370,   291,   344,   270,   217,   358,   217,
     217,    33,   405,   107,   310,   391,   159,   381,   403,   404,
     405,   406,   407,   614,   322,   401,   402,   291,   421,   360,
     361,   362,   363,   364,    33,   217,   316,   286,   431,   415,
     217,   434,   418,   436,   437,   421,   217,   213,   441,   322,
     217,   107,    74,   213,   213,    33,   322,   217,   202,   107,
     651,   322,   405,    85,    86,   408,   217,   410,   217,   217,
     213,   217,    44,   217,   217,    74,    33,   217,   103,   422,
     423,   424,   322,    95,   107,   381,    85,    86,   431,   159,
     433,   217,   159,   341,   342,   217,    74,   473,   441,   442,
     443,   444,   445,   308,   480,   107,   449,    85,    86,   452,
     453,   159,   379,   380,   381,   382,   383,    74,   385,   386,
     387,   388,   389,   217,   391,   392,   393,   394,    85,    86,
    1399,  1383,   884,    33,   159,    56,   159,   217,  1378,    33,
     345,    14,   841,   213,   154,   316,   213,   159,   217,    70,
     316,   217,  1351,   529,   530,   531,    33,   159,    95,   217,
     319,   217,   187,   107,   107,   213,   217,   153,   316,   217,
     191,   192,   323,   217,    74,    96,    70,   553,   217,    52,
      74,   159,   558,   559,   560,    85,    86,   217,   213,    33,
     213,    85,    86,    70,   316,   217,   572,    74,   574,   159,
     576,   213,    96,   175,   176,   177,   178,   217,    85,    86,
     586,   213,   588,   589,   590,   159,   159,   593,   217,    96,
     290,   217,   598,   290,   217,    33,    70,   603,   249,   605,
      74,   217,   159,   824,   217,   213,   316,   613,   159,   217,
     139,    85,    86,   619,  2066,   621,   622,   623,    33,   170,
     455,   217,    96,   213,   459,   460,   461,   154,   217,   159,
     217,    70,    70,   275,   154,   159,    74,   217,   644,   213,
     213,   647,   648,   649,   154,   212,   170,    85,    86,   655,
     154,    33,   159,   217,    56,    70,   213,    96,    96,    74,
     666,   217,   213,   170,   670,   194,   217,   318,    70,   675,
      85,    86,   282,   283,  1787,   255,   682,   683,   899,   208,
     331,    96,   903,   213,   690,   159,   692,   217,    70,   154,
     217,   255,    74,   217,    96,   154,   170,   217,   217,   705,
     257,   707,   708,    85,    86,   217,   159,   217,   275,   217,
     217,    69,    70,   217,    96,   217,   135,   211,   211,    70,
     159,   159,   289,   217,   217,   731,   732,   733,   949,   735,
     286,   170,   170,   739,    33,   315,   255,   743,    96,   209,
     159,   962,   963,   217,   159,    96,   752,   753,   754,   202,
      70,   315,   217,   255,   760,   170,   135,   159,   217,   217,
     213,   767,   154,   769,   217,   321,   772,   773,   170,    70,
     217,    70,   778,   217,   213,    74,    96,   159,   217,   217,
     159,    70,   154,    70,   217,   154,    85,    86,   170,   259,
      70,   797,   262,  1125,   213,    96,   315,    96,   804,   217,
     270,   159,   217,   809,    70,   154,   107,    96,   159,    96,
     280,   154,   170,   102,   820,   217,    96,   154,   217,   170,
      70,    70,   828,   829,   154,   217,   217,   833,   217,   135,
      96,   154,   154,   839,   213,   217,   842,   843,   405,   159,
     154,   408,   217,   410,   850,   217,    96,    96,   217,   855,
     170,   857,   858,   159,   154,   422,   423,   424,   159,   217,
     159,   217,   868,   154,   431,   135,   433,   873,   217,   170,
     159,   170,   159,   879,   217,   442,   443,   444,   445,   159,
     217,   170,   449,   170,   217,   452,   453,   217,   217,   159,
     170,   211,   898,   159,   217,   217,    70,   217,    70,   905,
     906,    70,   908,   217,   170,   911,    70,   213,   914,   159,
     159,   107,   217,    70,    70,    70,   217,   217,   217,   159,
     170,   170,    96,    70,    96,   159,   217,    96,   217,    70,
     217,    70,    96,    53,   217,   159,   209,   217,   172,    96,
      96,    96,   948,   213,   168,   951,   211,   191,   192,    96,
      70,   217,   217,   217,   960,    96,   217,    96,    70,   965,
     217,   967,   968,   969,   970,   971,   972,   217,   217,   937,
     938,   977,   217,   213,   217,   217,    96,   217,   217,   213,
     986,   217,    70,   217,    96,   159,   259,   159,   994,   213,
     159,    70,   998,   999,   217,   159,   170,   270,   170,   217,
     217,   170,   159,   159,   159,   249,   170,   280,    96,   201,
     202,   159,   159,   170,   170,   170,   217,    96,   159,   217,
     159,  1027,   217,   170,   172,   217,    70,  1033,    70,   170,
      70,   170,   217,  1039,  1040,  1041,  1042,   217,   217,   159,
     312,   313,   314,   217,   217,   217,   217,   159,   217,   159,
     170,   217,    96,   217,    96,   217,    96,   159,   170,    70,
     217,   217,   217,  1069,    70,   213,  1072,    70,  1074,   217,
     217,   159,    70,  1405,   217,   217,   217,  1083,   217,  1085,
     159,  1087,   170,  1089,    53,    96,   159,   217,  1094,   217,
      96,   170,  1098,    96,  1100,   217,   198,  1103,    96,  1105,
    1106,    70,   217,   213,   217,   217,   159,   217,   217,  1115,
     217,   213,  1118,  1119,  1120,   159,   217,   159,  1124,   159,
     408,  1127,  1128,  1129,  1130,   413,   170,    96,   170,   217,
     170,  1391,  1138,  1139,  1140,  1141,  1142,    70,   217,   217,
     213,   217,   217,   196,   217,  1151,   217,  1153,   159,   159,
    1156,  1411,  1412,   159,   217,   217,   159,   217,    70,   170,
     213,   159,   217,    96,   170,    70,    99,   170,    70,  1175,
     159,    70,   170,   217,  1180,   217,   217,   217,   217,   217,
     159,   469,   217,   217,    96,   217,  1192,   217,   217,  1195,
     159,    96,   217,  1199,    96,   217,   217,    96,   217,   487,
     488,   170,   217,   213,  1210,   217,   217,   217,  1214,    70,
      70,   217,   217,    70,   217,   217,  1222,  1223,   159,   217,
    1226,   217,   217,   217,   213,    70,   159,   217,   217,   217,
      70,   217,   217,   217,   213,    96,    96,   170,   217,    96,
    1246,  1247,   217,  1249,    70,   217,  1578,   159,  1254,   217,
    1582,    96,   217,    70,   159,  1261,    96,   159,   170,  1265,
     159,  1267,  1268,   217,   217,   170,   217,   217,   170,   217,
      96,   170,   213,  1279,   159,     0,   217,   159,  1284,    96,
     213,  1287,  1288,  1289,  1290,  1291,  1292,  1293,  1294,  1295,
    1296,  1297,  1298,  1299,  1300,   107,    70,    70,   159,   159,
     135,   191,   159,  1309,   456,   217,   194,   453,   159,   170,
     170,   194,   217,   170,   159,   217,   217,    70,   217,   159,
      70,    70,    96,    96,    70,   170,   159,     6,   213,   428,
     170,   213,   217,   159,     6,   217,   217,  1343,  1344,   217,
     159,   217,   159,    96,   170,   217,    96,    96,  1354,    70,
      96,   217,   217,   170,  1360,  1361,   217,   217,  1364,    53,
     217,   217,   213,  1369,   217,  1371,   217,   115,   116,   117,
     118,  1377,   217,    70,  1380,    96,    70,   217,   217,  1385,
     213,    70,  1388,   217,   217,   159,   159,  1393,  1394,  1395,
    1396,   217,  1398,   159,   213,  1401,   170,   170,   217,    96,
     217,    86,    96,   159,   217,   159,   159,    96,   107,   159,
     159,   217,   217,   159,   194,  1421,   217,   170,  1424,   217,
     170,   170,  1428,   217,   170,  1431,   159,  1433,  1434,  1435,
    1436,   159,   217,  1439,  1440,  1441,  1442,    53,   159,   159,
      70,  1447,   217,   217,   217,   217,    70,   213,   217,   170,
     159,   217,  1458,   159,   217,  1461,   159,   213,   217,   213,
      70,   217,   159,   217,   217,   159,    96,   217,   217,   277,
     159,   217,    96,   170,   217,   217,   170,    70,   217,   217,
     213,   170,  1488,    70,   217,   213,    96,  1493,  1494,   217,
    1496,  1497,  1498,   213,   217,   217,   217,   217,   217,   217,
     217,    70,   217,    96,   213,    70,    70,   213,   217,    96,
     213,   217,  1518,    70,   217,    70,   217,   217,    53,  1525,
     217,  1527,   217,  1529,   217,    70,  1532,    96,   217,   159,
     159,    96,    96,  1539,   217,   159,   217,  1543,  1544,    96,
     170,    96,   217,    70,    70,    70,   170,   217,   217,   159,
      70,    96,  1558,  1559,  1560,   217,   217,  1563,  1564,  1565,
     170,  1567,   217,  1569,  1570,  1571,   159,  1573,  1574,    96,
      96,    96,   159,  1579,    70,  1581,    96,   170,  1584,  1585,
    1586,  1587,    70,   170,   213,  1591,    53,   217,   217,   217,
     159,   217,  1598,   217,   159,   159,   217,  1603,   217,  1605,
      96,   170,   159,  1935,   159,   170,   170,   217,    96,   217,
      70,  1943,    70,   170,   159,   170,    70,    33,   284,   217,
      70,   159,   281,   107,   217,   170,   159,   159,   316,   316,
     217,   107,   159,   159,   159,   159,    96,   107,    96,   159,
     217,   107,    96,   170,   170,   170,    96,   217,   217,   217,
     170,   107,   217,   217,   107,   217,   217,   159,    74,   159,
     217,    70,   217,   159,   217,   159,   457,   107,   107,    85,
      86,   159,   217,   107,   170,   213,    70,    70,   316,   217,
     213,   213,   170,   217,   217,   217,   217,    96,   217,   213,
     217,   217,   217,   217,   316,   217,   316,   217,  1704,   159,
     217,   159,    96,    96,    70,   159,    70,   107,    66,   159,
     170,   213,   170,   213,   217,   217,   170,   217,  1724,   213,
     170,   217,    70,   217,  1730,   316,   217,  1733,    70,   217,
      96,  1737,    96,  1739,  1740,  1741,  1742,  1743,  1744,  1745,
    1746,  1747,  1748,  1749,  1750,  1751,  1752,  1753,    96,  2081,
     159,   294,  1758,   217,    96,    70,    70,   217,   217,   217,
     217,   170,   217,   217,   316,   159,   159,   217,  1774,  1775,
    1776,   217,   217,   294,  1780,  1781,   170,   170,   294,  1785,
     316,    96,    96,    70,  1790,  1791,   316,   217,  1794,  1795,
     217,   217,  1798,   159,   217,   159,   294,   316,   316,  1805,
      70,    70,   217,   217,   170,   159,   170,   316,   217,    96,
    1816,   159,   159,   217,  1820,   316,  1822,   159,   217,  1825,
     159,   162,   170,   217,   217,   217,    96,    96,   170,  1835,
    1836,  1837,  1838,  1839,   159,   217,  1842,   159,  1844,   217,
    1846,   217,   217,   107,   159,   159,   107,   316,   217,   217,
      70,   217,  2184,   217,    70,   170,   170,   107,  2190,   213,
    1866,  1867,  1868,   217,  1870,  1871,   213,    70,    70,   217,
     217,    70,   159,   316,   213,   217,    96,  1883,   217,   316,
      96,   316,  1888,   170,   316,  1891,   409,   316,   213,   159,
     159,   213,   217,    96,    96,   217,   107,    96,    70,  2231,
     170,   170,   217,   217,  1910,  1911,  1912,  1913,   316,  1915,
    1916,  1917,  1918,  1919,  1920,  1921,   107,  1923,  1924,  1925,
    1926,   159,  1928,  1929,    96,   316,   107,   159,   159,  1935,
     217,  1937,   316,   316,  1940,   217,   159,  1943,   217,   159,
    1946,  1947,  1948,   159,   217,   217,  1952,   217,   217,  1955,
     170,   159,   107,   159,   170,  1961,   159,   159,    70,  1965,
     159,   217,    70,    70,   217,   217,    70,   170,   170,   217,
      70,   170,   217,   217,   217,   213,    70,   217,  1984,   217,
     194,   213,   213,   217,    96,   217,   217,   159,    96,    96,
     213,   217,    96,   217,   217,    99,    96,   217,   170,    99,
     217,   217,    96,   316,   217,   213,  2338,   213,   217,   217,
     217,   217,  2344,   217,   217,   217,   217,    70,   217,  2025,
     217,  2027,   217,  2029,   217,   159,  2032,  2033,   217,  2035,
     217,  2037,  2038,  2039,  2040,  2041,  2042,  2043,  2044,   217,
    2046,  2047,  2048,    96,   159,   217,   217,   159,   217,    70,
    2056,   159,   159,  2059,  2060,   159,   217,  2063,   170,   159,
     159,   217,   170,   170,  2070,   159,   170,  2073,   217,   217,
     170,   217,   217,   217,   127,    96,   170,   159,   159,   213,
     217,   217,  2088,   217,  2090,   159,   159,   159,   159,   217,
    2096,  2097,   217,  2099,   217,   217,  2102,  2103,   213,  2105,
    2106,   217,   217,   217,   217,   217,   159,   217,   217,   217,
     217,   217,   217,   213,   213,  2121,  2122,   170,   217,   213,
     217,   159,   143,   159,   217,  2131,  2132,  2133,   217,  2135,
    2136,   213,   213,   217,   217,   217,   217,   159,   159,   213,
     213,   213,   213,   217,   217,   217,   217,   217,  2154,   170,
     159,   217,   159,  2159,   217,  2161,  2162,  2163,   159,   159,
     159,  2167,  2168,   217,   294,   159,  2172,  2173,  2174,  2175,
    2176,  2177,  2178,  2179,  2180,   213,   217,   213,  2184,   217,
     217,   217,   203,   217,  2190,   217,   217,  2193,  2194,   159,
     217,   213,   159,   217,   159,   217,  2202,  2203,  2530,   159,
    2532,   217,   159,   217,   213,  2211,   213,   217,   217,   217,
     217,   217,   213,   213,   213,   217,   217,   217,   217,   213,
     217,   217,   217,   217,   217,   217,  2232,  2233,  2234,  2235,
    2236,  2237,   217,  2239,  2240,  2241,  2242,  2243,  2244,  2245,
     217,  2247,   217,   213,  2250,  2251,   213,   217,   213,   217,
     217,   217,   217,   213,   107,  2587,   213,   217,   217,  2591,
     217,   217,   217,   217,   217,   217,   217,  2273,   217,   217,
     217,  2277,   217,  2279,   217,   217,   239,  2283,   241,   242,
     243,   244,   245,   246,   247,   248,   217,   217,   217,  2295,
     217,  2297,  2298,   217,  2300,  2301,   217,   217,   217,   217,
     217,   217,   217,  2309,   217,   217,  2312,  2313,  2314,  2315,
    2316,   217,   217,   217,  2320,  2321,   217,   217,   217,  2325,
     217,   217,   107,  2329,  2330,  2331,  2332,  2333,  2334,   217,
     107,   107,  2664,   217,   217,   217,   217,   217,   217,   217,
     217,  2347,  2348,   217,   217,   217,   217,  2353,   217,   277,
     217,   217,   217,   217,   374,   217,   217,   217,  2690,   217,
     217,   217,   217,   217,   107,   217,   217,   107,  2374,   107,
    2376,  2377,   217,  2379,   217,  2381,  2382,   217,  2384,   217,
     217,   217,  2388,  2389,  2390,   217,  2392,  2393,  2394,   217,
     217,  2397,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   217,  2419,   217,   217,   217,   217,  2424,  2425,
    2426,   217,  2428,   217,   217,   217,   217,   217,   217,  2435,
     217,   217,   107,   217,   217,  2441,  2442,   217,  2444,   217,
     107,   107,  2448,  2449,  2450,  2451,  2452,   217,   217,   217,
    2456,   217,   217,   217,  2460,   217,  2462,  2463,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
    2476,  2477,   217,   217,   171,  2481,   217,  2483,   217,  2485,
    2486,   217,   217,   217,  2490,   217,  2492,  2493,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,  2508,  2509,  2510,   217,  2512,   217,   217,   217,
     217,   217,  2518,   217,  2520,  2521,  2522,   217,  2524,  2525,
    2526,  2527,  2528,   217,   217,   217,  2532,  2533,  2534,   217,
     277,   217,   217,   217,   217,   281,  2542,   217,  2544,   217,
     217,   217,   107,  2549,   217,  2551,   217,   217,   217,  2555,
     281,   107,   217,   217,   217,   217,   217,  2563,  2564,   217,
    2566,   217,   217,  2569,   217,   217,  2572,   217,   217,  2575,
    2576,   217,  2578,  2579,  2580,   217,  2582,   217,   217,   217,
     217,  2587,  2588,  2589,   217,   217,   217,  2593,   217,  2595,
     217,   217,  2598,   217,   217,  2601,   217,  2603,   217,   217,
     217,   217,   217,   217,  2610,   217,  2612,   217,  2614,   217,
     217,   217,  2618,   217,   217,  2621,   217,   217,  2624,  2625,
    2626,  2627,   217,   217,   217,   217,  2632,  2633,  2634,   217,
     217,   217,   217,   217,  2640,   217,   217,   217,   217,   217,
     217,  2647,   217,   217,   217,  2651,   217,   217,   217,   217,
     217,  2657,   297,  2659,  2660,   217,  2662,   217,  2664,  2665,
     217,  2667,   217,   217,   217,  2671,   217,   217,   217,   217,
     217,  2677,   217,   217,  2680,   217,   217,   217,   217,   217,
     217,  2687,   217,   217,   217,   217,  2692,  2693,  2694,   217,
    2696,   217,   217,   217,   217,   217,  2702,   217,   217,   217,
    2706,   217,  2708,   217,  2562,   568,   297,  2713,   593,  2715,
     357,   570,   988,  2719,   373,   885,   984,  2723,   912,  2725,
     896,   412,   952,  2729,  2730,  2731,  2732,  2733,  2734,  2735,
    2736,  1132,   404,  2739,  2740,   431,   509,  2743,  2744,   438,
    2746,  2747,  2748,  2749,  2750,  2751,   479,  2753,   794,   313,
    2756,  1369,    -1,    -1,  2760,     4,    -1,    -1,  2764,  2765,
    2766,    -1,    11,    -1,    -1,    14,    15,    -1,    -1,    18,
      19,    -1,    -1,    22,    -1,    24,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      39,    -1,    -1,    -1,    -1,    -1,    -1,    46,    47,    48,
      49,    -1,    -1,    -1,    -1,    54,    -1,    -1,    -1,    -1,
      -1,    60,    -1,    62,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    76,    -1,    78,
      -1,    -1,    81,    -1,    -1,    -1,    -1,    86,    -1,    88,
      -1,    -1,    91,    92,    93,    -1,    -1,    -1,    97,    -1,
      -1,   100,   101,    -1,    -1,   104,   105,    -1,    -1,   108,
     109,   110,   111,   112,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   122,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   142,    -1,    -1,   145,   146,    -1,   148,
     149,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   176,    -1,    -1,
      -1,    -1,   181,    -1,    -1,    -1,   185,    -1,   187,   188,
     189,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   224,    -1,   226,   227,    -1,
      -1,   230,   231,   232,   233,    -1,   235,    -1,   237,   238,
     239,   240,    -1,    -1,   243,   244,   245,    -1,    -1,   248,
     249,    -1,    -1,   252,   253,    -1,    -1,   256,   257,   258,
      -1,    -1,    -1,   262,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   270,   271,   272,    -1,    -1,   275,   276,   277,   278,
     279,    -1,   281,   282,   283,   284,   285,   286,   287,    -1,
     289,   290,   291,    -1,   293,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   306,   307,   308,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   316,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   331,    -1,    -1,   334,   335,    -1,    -1,    -1,
      -1,    -1,   341,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   373,   374,   375,    -1,    -1,   378,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   393,   394,    -1,   396,    -1,    -1,
      -1,   400,    -1,    -1,   403,   404,    -1,    -1,   407,   408,
      -1,    -1,    -1,   412,   413,    -1,    -1,   416,    -1,    -1,
     419,    -1,    -1,    -1,   423,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   431,    -1,    -1,    -1,    -1,    -1,    -1,   438,
     439,    -1,    -1,    -1,   443,    -1,    -1,   446,    -1,    -1,
      -1,   450,    -1,    -1,   453,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   461,    -1,   463,    -1,    -1,    -1,    -1,    -1,
     469,    -1,    -1,    -1,    -1,   474,    -1,    -1,    -1,    -1,
     479,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   487,   488,
      -1,    -1,    -1,   492,    -1,    -1,   495,    -1,    -1,    -1,
     499,    -1,    -1,   502,    -1,    -1,    -1,    -1,    -1,   508,
     509,   510,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   526,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   548,
      -1,    -1,    -1,    -1,    -1,   554,   555,   556,   557,    -1,
      -1,    -1,   561,   562,   563,   564,   565,    -1,   567,   568,
     569,   570,    -1,    -1,    -1,    -1,   575,    -1,   577,    -1,
      -1,   580,   581,   582,    -1,    -1,    -1,    -1,   587,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   596,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   611,    -1,    -1,   614,    -1,   616,    -1,   618,
      -1,   620,    -1,    -1,    -1,    -1,   625,    -1,   627,    -1,
     629,    -1,   631,    -1,   633,    -1,    -1,    -1,    -1,    -1,
      -1,   640,    -1,   642,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   651,    -1,   653,    -1,    -1,    -1,   657,   658,
     659,   660,    -1,   662,   663,   664,    -1,    -1,    -1,    -1,
      -1,    -1,   671,   672,    -1,    -1,    -1,    -1,    -1,   678,
      -1,   680,   681,   682,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   695,   696,    -1,   698,
     699,   700,   701,    -1,    -1,   704,    -1,    -1,    -1,    -1,
      -1,   710,    -1,   712,   713,    -1,   715,   716,   717,   718,
      -1,   720,   721,   722,   723,   724,    -1,   726,    -1,   728,
     729,   730,    -1,    -1,    -1,   734,    -1,    -1,   737,    -1,
      -1,   740,   741,   742,    -1,   744,   745,    -1,   747,   748,
     749,    -1,   751,    -1,    -1,    -1,   755,   756,   757,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   770,   771,    -1,    -1,   774,   775,    -1,    -1,    -1,
     779,    -1,    -1,    -1,    -1,   784,    -1,    -1,   787,   788,
     789,    -1,   791,   792,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   802,    -1,    -1,    -1,    -1,   807,   808,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   824,    -1,   826,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   836,    -1,   838,
      -1,   840,   841,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   852,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   861,    -1,    -1,    -1,    -1,   866,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   880,    -1,   882,   883,   884,   885,    -1,   887,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   896,   897,    -1,
     899,    -1,    -1,   902,   903,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   912,   913,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   923,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   932,    -1,    -1,   935,    -1,    -1,    -1,
      -1,   940,    -1,    -1,   943,    -1,    -1,    -1,   947,    -1,
     949,    -1,    -1,   952,    -1,   954,    -1,    -1,    -1,    -1,
     959,    -1,    -1,   962,   963,   964,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   984,    -1,    -1,    -1,   988,
      -1,   990,    -1,   992,   993,    -1,    -1,    -1,    -1,    -1,
      -1,  1000,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1024,  1025,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1036,    -1,  1038,
      -1,    -1,  1041,    -1,    -1,    -1,  1045,    -1,    -1,  1048,
    1049,  1050,  1051,  1052,  1053,    -1,    -1,    -1,    -1,    -1,
    1059,    -1,  1061,    -1,  1063,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1096,  1097,    -1,
    1099,    -1,  1101,    -1,    -1,    -1,    -1,    -1,  1107,  1108,
    1109,  1110,  1111,    -1,  1113,  1114,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1123,    -1,  1125,  1126,    -1,    -1,
      -1,    -1,  1131,  1132,  1133,  1134,  1135,  1136,  1137,    -1,
      -1,    -1,    -1,    -1,  1143,    -1,    -1,    -1,    -1,    -1,
      -1,  1150,    -1,  1152,    -1,    -1,  1155,    -1,    -1,  1158,
      -1,    -1,    -1,    -1,    -1,  1164,    -1,    -1,    -1,    -1,
      -1,  1170,    -1,    -1,    -1,    -1,    -1,  1176,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1201,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1209,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1224,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1257,  1258,
    1259,  1260,    -1,  1262,  1263,  1264,    -1,    -1,    -1,    -1,
    1269,    -1,    -1,    -1,    -1,  1274,    -1,    -1,    -1,    -1,
      -1,  1280,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1301,    -1,  1303,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1324,    -1,    -1,    -1,    -1,
      -1,  1330,    -1,    -1,    -1,    -1,    -1,    -1,  1337,    -1,
      -1,  1340,  1341,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1351,    -1,    -1,    -1,    -1,    -1,    -1,  1358,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1375,    -1,    -1,  1378,
      -1,  1380,    -1,    -1,  1383,    -1,  1385,    -1,  1387,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1399,    -1,  1401,    -1,  1403,    -1,    -1,    -1,    -1,    -1,
    1409,  1410,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1420,    -1,    -1,    -1,  1424,    -1,    -1,  1427,    -1,
      -1,  1430,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1438,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1455,  1456,  1457,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1465,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1487,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1505,    -1,  1507,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1530,    -1,    -1,  1533,    -1,    -1,    -1,  1537,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1545,    -1,  1547,  1548,
    1549,  1550,    -1,    -1,    -1,    -1,  1555,    -1,    -1,    -1,
      -1,    -1,    -1,  1562,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1575,    -1,  1577,  1578,
      -1,  1580,    -1,  1582,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1592,    -1,    -1,    -1,    -1,    -1,    -1,
    1599,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1608,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1637,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1645,    -1,  1647,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1655,    -1,  1657,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1723,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1762,    -1,  1764,    -1,  1766,  1767,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1791,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1807,    -1,
      -1,    -1,    -1,    -1,  1813,    -1,    -1,    -1,    -1,    -1,
      -1,  1820,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1829,  1830,  1831,  1832,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1848,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1859,    -1,    -1,    -1,    -1,  1864,  1865,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1876,    -1,  1878,
    1879,    -1,  1881,  1882,    -1,    -1,  1885,    -1,  1887,    -1,
      -1,  1890,    -1,    -1,    -1,    -1,    -1,    -1,  1897,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1914,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1945,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1956,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1986,    -1,    -1,
      -1,  1990,    -1,    -1,    -1,  1994,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2021,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2054,  2055,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2084,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2108,
      -1,    -1,    -1,    -1,  2113,  2114,  2115,  2116,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2130,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2139,    -1,  2141,    -1,    -1,  2144,    -1,    -1,    -1,    -1,
      -1,  2150,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2164,    -1,    -1,    -1,    -1,
      13,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    29,    -1,    -1,    -1,
      -1,    34,    -1,  2192,    -1,    38,    -1,    -1,    -1,    -1,
      43,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    51,    -1,
      -1,    54,    55,    -1,    57,    -1,    59,    60,    -1,  2218,
      63,  2220,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2229,    74,    -1,    -1,    -1,    -1,    -1,    -1,    81,    82,
      83,    -1,    -1,    -1,    -1,    -1,    -1,  2246,    -1,  2248,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   101,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   121,   122,
     123,   124,    -1,    -1,    -1,  2284,   129,    -1,    -1,   132,
     133,    -1,    -1,   136,   137,    -1,    -1,   140,    -1,    -1,
      -1,    -1,    -1,  2302,    -1,  2304,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   160,    -1,   162,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     173,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2338,
      -1,    -1,   185,   186,    -1,  2344,    -1,    -1,    -1,    -1,
      -1,    -1,   195,    -1,   197,    -1,   199,   200,    -1,    -1,
      -1,   204,   205,   206,   207,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2373,    -1,    -1,    -1,   221,   222,
      -1,    -1,    -1,    -1,   227,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2395,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     253,   254,    -1,   256,    -1,   258,  2415,   260,    -1,    -1,
      -1,    -1,   265,   266,   267,   268,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   282,
      -1,    -1,    -1,    -1,    -1,    -1,   289,    -1,    -1,    -1,
      -1,   294,    -1,   296,   297,   298,   299,   300,   301,   302,
     303,   304,   305,   306,    -1,    -1,    -1,    -1,   311,   312,
      -1,   314,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   329,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   338,   339,    -1,    -1,    -1,
      -1,    -1,    -1,  2502,    -1,    -1,  2505,    -1,    -1,   352,
      -1,    -1,   355,   356,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   369,    -1,    -1,   372,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2540,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2560,  2561,  2562,    -1,    -1,    -1,    -1,    -1,   412,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2577,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2607,  2608,
    2609,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2646,    -1,     3,
       4,     5,    -1,    -1,    -1,     9,    -1,    11,    -1,    -1,
      -1,    -1,    16,    -1,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    -1,    30,    31,    -1,    -1,
      -1,    35,    36,    37,    -1,  2684,    -1,    -1,    -1,    -1,
      -1,    45,    -1,    47,    48,    49,    50,    -1,    -1,    -1,
      -1,  2700,    -1,    -1,    58,    -1,    -1,    -1,    -1,    -1,
      64,    65,    -1,  2712,    68,    -1,    -1,    71,    72,    73,
      -1,    -1,    -1,    77,    -1,    79,    -1,    -1,    -1,    -1,
      84,    -1,    -1,    87,    88,    89,    90,    91,    92,    93,
      94,    -1,    -1,    -1,    -1,    -1,   100,    -1,    -1,    -1,
     104,   105,    -1,  2752,   108,   109,    -1,   111,   112,   113,
      -1,    -1,    -1,  2762,  2763,    -1,   120,    -1,    -1,    -1,
      -1,   125,    -1,    -1,   128,    -1,   130,   131,    -1,    -1,
     134,    -1,    -1,    -1,   138,    -1,    -1,   141,   142,   143,
     144,   145,   146,   147,   148,   149,   150,   151,   152,    -1,
      -1,    -1,    -1,   157,   158,    -1,    -1,   161,    -1,   163,
     164,   165,   166,   167,    -1,   169,    -1,   171,    -1,    -1,
     174,    -1,    -1,    -1,    -1,   179,    -1,   181,   182,   183,
     184,    -1,    -1,   187,   188,   189,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   203,
      -1,    -1,    -1,    -1,    -1,    -1,   210,    -1,    -1,    -1,
      -1,   215,    -1,    -1,    -1,   219,    -1,    -1,    -1,   223,
     224,   225,    -1,    -1,    -1,   229,   230,   231,   232,    -1,
     234,    -1,   236,   237,   238,    -1,   240,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   251,   252,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   263,
     264,    -1,    -1,    -1,    -1,   269,    -1,   271,    -1,   273,
      -1,    -1,   276,    -1,   278,   279,    -1,    -1,    -1,   283,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   292,   293,
      -1,   295,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   307,    -1,   309,    -1,    -1,    -1,   313,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   326,    -1,   328,    -1,   330,    -1,   332,   333,
     334,   335,   336,   337,    -1,    -1,   340,    -1,    -1,    -1,
      -1,    -1,   346,    -1,    -1,   349,   350,   351,    -1,    -1,
      -1,    -1,    -1,   357,    -1,   359,    -1,    -1,    -1,    -1,
      -1,   365,   366,    -1,    -1,    -1,    -1,    -1,    -1,   373,
      -1,    -1,    -1,   377,   378,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   395,   396,   397,   398,    -1,    -1,    -1,   402,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   411,    -1,   413,
     414,   415,   416,   417,   418,   419,   420,    -1,    -1,    -1,
      -1,    -1,   426,   427,    -1,   429,   430,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   446,   447,   448,    -1,   450,   451,    -1,    -1,
     454,     3,     4,     5,   458,    -1,    -1,     9,    -1,    11,
      -1,    -1,    -1,    -1,    16,    -1,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    -1,    30,    31,
      -1,    -1,    -1,    35,    36,    37,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    45,    -1,    47,    48,    49,    50,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    58,    -1,    -1,    -1,
      -1,    -1,    64,    65,    -1,    -1,    68,    -1,    -1,    71,
      72,    73,    -1,    -1,    -1,    77,    -1,    79,    -1,    -1,
      -1,    -1,    84,    -1,    -1,    87,    88,    89,    90,    -1,
      92,    93,    94,    -1,    -1,    -1,    -1,    -1,   100,    -1,
      -1,    -1,   104,   105,    -1,    -1,   108,   109,    -1,   111,
     112,   113,    -1,    -1,    -1,    -1,    -1,    -1,   120,    -1,
      -1,    -1,    -1,   125,    -1,    -1,   128,    -1,   130,   131,
      -1,    -1,   134,    -1,    -1,    -1,   138,    -1,    -1,   141,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   151,
     152,    -1,    -1,    -1,    -1,   157,   158,    -1,    -1,   161,
      -1,   163,   164,   165,   166,   167,    -1,   169,    -1,   171,
      -1,    -1,   174,    -1,    -1,    -1,    -1,   179,    -1,   181,
     182,   183,   184,    -1,    -1,   187,   188,   189,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   203,    -1,    -1,    -1,    -1,    -1,    -1,   210,    -1,
      -1,    -1,    -1,   215,    -1,    -1,    -1,   219,    -1,    -1,
      -1,   223,   224,   225,    -1,    -1,    -1,   229,   230,   231,
     232,    -1,   234,    -1,   236,   237,   238,    -1,   240,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   251,
     252,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   263,   264,    -1,    -1,    -1,    -1,   269,    -1,   271,
      -1,   273,    -1,    -1,   276,    -1,   278,   279,    -1,    -1,
      -1,   283,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     292,   293,    -1,   295,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   307,    -1,   309,    -1,    -1,
      -1,   313,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   326,    -1,   328,    -1,   330,    -1,
     332,   333,   334,   335,   336,   337,    -1,    -1,   340,    -1,
      -1,    -1,    -1,    -1,   346,    -1,    -1,   349,   350,   351,
      -1,    -1,    -1,    -1,    -1,   357,    -1,   359,    -1,    -1,
      -1,    -1,    -1,   365,   366,    -1,    -1,    -1,    -1,    -1,
      -1,   373,    -1,    -1,    -1,   377,   378,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   395,   396,   397,   398,    -1,    -1,    -1,
     402,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   411,
      -1,   413,   414,   415,   416,   417,   418,   419,   420,    -1,
      -1,    -1,    -1,    -1,   426,   427,    -1,   429,   430,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   446,   447,   448,    -1,   450,   451,
      -1,    -1,   454,    -1,    -1,    -1,   458
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     3,     4,     5,     9,    11,    16,    18,    19,    20,
      21,    22,    23,    24,    25,    26,    27,    28,    30,    31,
      35,    36,    37,    45,    47,    48,    49,    50,    58,    64,
      65,    68,    71,    72,    73,    77,    79,    84,    87,    88,
      89,    90,    92,    93,    94,   100,   104,   105,   108,   109,
     111,   112,   113,   120,   125,   128,   130,   131,   134,   138,
     141,   142,   143,   144,   145,   146,   147,   148,   149,   150,
     151,   152,   157,   158,   161,   163,   164,   165,   166,   167,
     169,   171,   174,   179,   181,   182,   183,   184,   187,   188,
     189,   203,   210,   215,   219,   223,   224,   225,   229,   230,
     231,   232,   234,   236,   237,   238,   240,   251,   252,   263,
     264,   269,   271,   273,   276,   278,   279,   283,   292,   293,
     295,   307,   309,   313,   326,   328,   330,   332,   333,   334,
     335,   336,   337,   340,   346,   349,   350,   351,   357,   359,
     365,   366,   373,   377,   378,   395,   396,   397,   398,   402,
     411,   413,   414,   415,   416,   417,   418,   419,   420,   426,
     427,   429,   430,   446,   447,   448,   450,   451,   454,   458,
     463,   464,   465,   466,   467,   468,   469,   470,   475,   476,
     477,   478,   479,   480,   481,   482,   483,   484,   493,   495,
     496,   497,   498,   499,   500,   501,   502,   503,   504,   506,
     508,   509,   510,   511,   512,   513,   514,   519,   522,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   535,
     536,   538,   542,   543,   544,   546,   547,   548,   549,   552,
     555,   556,   557,   558,   559,   562,   565,   568,   570,   572,
     574,   578,   579,   580,   581,   582,   583,   584,   585,   587,
     589,   590,   591,   593,   595,   596,   597,   598,   599,   600,
     601,   602,   603,   604,   608,   610,   612,   614,   618,   622,
     624,   626,   627,   628,   629,   630,   631,   632,   634,   635,
     636,   637,   645,   646,   648,   650,   651,   652,   653,   654,
     655,   656,   657,   658,   659,   663,   665,   667,   670,   673,
     675,   677,   679,   681,   682,   683,   684,   685,   686,   689,
     690,   692,   693,   694,   695,   698,   217,     6,   217,   217,
     159,   213,   700,   217,    70,    96,   159,   170,   701,   701,
     701,   107,   217,   701,   217,   217,   700,   217,   217,   700,
     700,   217,   217,   700,   700,   217,   217,   700,   217,   217,
     700,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   217,    33,    74,    85,    86,   217,   666,   217,
     217,   700,   217,   217,   217,   217,   217,   217,   700,   700,
     217,   700,   217,   700,   217,   217,   316,   217,   700,   217,
     217,   217,   217,   217,   217,   700,   217,   217,   700,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   358,
     217,   701,   217,   217,   700,   217,   217,   700,   217,   217,
     700,   217,   217,   217,   107,   217,   700,   217,   217,   700,
     217,   217,   217,   700,   217,   700,   700,   217,   217,   217,
     700,   217,   107,   172,   217,   700,   172,   217,   700,   217,
     217,   700,   700,   217,   217,   217,   700,   700,   700,   217,
     700,   217,   700,   107,   217,   217,   217,   217,   217,   217,
     701,   217,   217,   700,   217,   217,   701,   701,   217,   217,
     217,   217,   217,   217,   217,   217,   107,   217,   217,   217,
     217,   217,   366,   217,   217,   700,   217,   217,   700,   700,
     701,   700,   700,   217,   217,   701,   217,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   316,     0,    91,   465,   135,   700,   209,   259,
     262,   270,   280,   471,   472,   474,   541,    39,    40,   360,
     361,   362,   363,   364,    78,    80,   119,   348,   403,   404,
     405,   406,   407,   700,    17,   114,   139,   208,   701,   701,
     700,    44,   175,   176,   177,   178,   485,   487,   488,   489,
     490,   647,   700,   494,   700,   700,   701,   107,   456,   505,
     308,   345,   455,   459,   460,   461,   507,   191,    17,    42,
     139,   155,   194,   208,   216,   228,   317,   325,   327,   358,
     516,   537,   541,   139,   194,   208,   520,   521,    32,    46,
     363,   523,   545,   700,   322,   641,   700,   700,   322,   700,
     322,   700,   700,   700,   563,   700,   569,   700,   571,   700,
     573,   700,   575,   700,   563,   571,   575,   586,   700,   588,
     700,   592,   700,   594,   700,   453,   194,   700,   700,   700,
     194,   322,   641,   700,   625,   700,   639,   700,   700,   700,
     700,   633,   700,   700,   700,   638,   700,   647,   647,   649,
     700,   700,   700,    95,   275,   700,   701,   159,   700,   701,
     700,   322,   700,   700,    13,    29,    34,    38,    43,    51,
      54,    55,    57,    59,    60,    63,    74,    81,    82,    83,
     101,   121,   122,   123,   124,   129,   132,   133,   136,   137,
     140,   160,   162,   173,   185,   186,   195,   197,   199,   200,
     204,   205,   206,   207,   221,   222,   227,   253,   254,   256,
     258,   260,   265,   266,   267,   268,   282,   289,   294,   296,
     297,   298,   299,   300,   301,   302,   303,   304,   305,   306,
     311,   312,   314,   329,   338,   339,   352,   355,   356,   369,
     372,   412,   661,   662,   471,   474,   541,   701,    10,    75,
      86,   106,   126,   180,   186,   191,   218,   220,   226,   274,
     353,   381,   680,   217,   700,   384,   700,   700,    41,   405,
     421,   431,   434,   436,   437,   441,   691,   697,    95,   212,
     275,   289,   405,   408,   410,   422,   423,   424,   431,   433,
     442,   443,   444,   445,   449,   452,   453,   696,   696,   696,
     697,   696,   428,   699,   605,   641,   700,   217,     6,     6,
     217,   217,   217,   701,   217,   700,   701,   576,   577,   700,
     577,   217,   701,   701,   217,   217,   700,   217,   217,   217,
     701,   217,   701,   107,   515,   701,   701,   701,    86,   701,
     666,   154,   217,   217,   569,   563,   607,   643,   700,   664,
     701,   217,   700,   217,   217,   194,   217,   217,   217,   701,
     623,   644,   700,   623,   217,   560,   561,   700,   217,   571,
     674,   701,   676,   701,   575,   573,   619,   620,   700,   605,
     217,   217,   607,   605,   217,   701,   700,   217,   701,   700,
     217,   701,   615,   616,   700,   217,   217,   217,   625,   217,
     217,   217,   638,   700,   217,   217,   217,   700,   217,   217,
     700,   217,   700,   217,   217,   700,   217,   217,   217,   217,
     534,   700,   217,   107,   217,   700,   277,    53,   613,   605,
     217,   701,   566,   567,   700,   217,   217,   639,   701,    53,
     611,   217,   605,   605,    53,   609,   700,   115,   116,   117,
     118,   374,   473,   217,   700,   217,   286,   321,   700,   217,
     217,   217,   700,   217,   550,   551,   700,   633,   553,   554,
     700,   217,   322,   700,   284,   217,   701,   217,   701,   701,
     281,   316,   316,   107,   107,   107,   107,   107,   217,   107,
     217,   217,   700,   217,   217,   217,   217,   701,   700,   700,
     700,   700,   701,   701,   143,   203,   701,   700,   700,   700,
     700,   700,   491,   700,   491,   492,   700,   492,   343,   701,
     701,   700,   701,   107,   700,   457,   217,   310,   107,   700,
     107,   700,   107,   700,   107,   107,   107,    12,    15,    61,
     190,   214,   217,   223,   233,   287,   288,   310,   381,   701,
     700,   517,   701,   217,   701,   518,   701,   316,   217,   517,
     217,   217,   316,   700,   316,   701,   316,   701,   217,   701,
     217,   107,   217,   700,   701,   641,   700,   347,   700,    66,
     700,   347,   701,   700,   701,   701,   701,   700,   700,   700,
     700,   700,   217,   700,   700,   701,   316,   217,   701,   701,
     701,   217,   641,   347,     8,    62,    98,   188,   323,   344,
     701,   700,   700,   700,   700,   700,   700,   700,   701,   127,
     701,   700,   342,   700,   217,   217,   701,   217,   700,   217,
     700,   700,    99,   700,   701,   347,   701,    14,   370,   294,
     217,   316,   294,    14,    52,   701,   294,   701,   316,   316,
     660,   700,   217,   700,   316,   700,   700,    76,   217,   375,
     700,   290,   700,   217,   217,   700,   701,   217,   701,   701,
     316,   217,   700,   217,   217,   700,   168,   700,   294,   700,
     196,   700,   198,   700,   290,   700,   217,   700,   700,   700,
     700,   217,   700,   217,   700,   107,   217,   700,   257,   700,
     700,   701,   701,   701,   700,   217,   701,   316,   290,   700,
     217,   255,   353,   701,   700,   700,   700,   701,   700,   700,
     316,   700,   700,   700,   217,   700,   701,   701,   701,   700,
     700,   700,   217,   162,   701,   217,   316,   191,   192,   249,
     318,   331,   191,   192,   249,   701,   217,   701,   700,   700,
     217,   701,   701,     7,   678,   700,   700,   217,   217,   701,
     700,   217,   217,   107,   700,   107,   107,   379,   380,   381,
     382,   383,   385,   386,   387,   388,   389,   391,   392,   393,
     394,   700,   700,   700,   688,   700,   316,   700,   316,   687,
     688,   217,   217,   515,   316,   316,   316,   316,   700,   409,
     701,   107,   107,   700,   700,   701,   316,   316,   107,   316,
     107,   316,   316,   217,   217,   217,   217,   107,   217,   641,
     700,   347,   217,   701,   701,   701,   217,   217,   700,   576,
     701,   577,   701,   701,   217,   217,   701,   700,   217,   217,
     701,   701,   701,   217,   700,   643,   701,   217,   217,   671,
     672,   701,   217,   194,   217,   701,   644,   700,   623,   561,
     700,   217,   217,   620,   621,   700,   621,   700,    56,   217,
     701,   668,   669,   701,   701,   701,   701,   217,   701,   616,
     617,   700,   617,   700,   217,   701,   700,   217,   217,   700,
     700,   668,   668,   217,   700,   217,   700,   217,   217,   700,
      53,   701,   701,   567,   700,   217,   700,    53,   701,   700,
      53,   701,   217,   701,   701,   701,   701,   701,   701,   115,
     116,   117,   118,   217,   217,   701,   217,   399,   400,   401,
     217,   551,   217,   701,   554,   700,   700,   700,   701,   606,
     642,   700,   217,   701,   701,   700,   107,   217,   107,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   217,   217,   700,   700,   217,   217,   701,   217,
     217,   217,   217,   701,   700,   700,   701,   701,   700,   701,
     701,   107,   217,   217,   700,   700,   700,   700,   700,   700,
     700,   700,   700,   700,   316,   701,   217,   217,   701,   701,
     217,   217,   217,   217,   217,   701,   217,   701,   217,   701,
     701,   217,   701,   700,   211,   217,   701,   700,   701,   700,
     217,   701,   217,   701,   701,   564,   700,   564,   564,   564,
     564,   217,   700,   217,   700,   701,   217,   606,   701,   701,
     701,   606,   700,   701,    97,   188,   217,   323,   666,   700,
     700,   701,   701,   701,   701,   640,   700,   640,   700,   700,
     564,   700,   700,   217,   701,   701,   701,   701,   701,   156,
     217,   342,   700,   217,   217,   217,   700,   701,   700,    99,
     701,   217,   316,   700,   217,   701,   217,   272,   371,   700,
     217,   255,   294,   217,   217,   255,   315,   217,   217,   700,
     217,   217,   255,   315,   217,   217,   217,   217,   700,   217,
     217,   701,   700,   217,   217,   102,   217,   701,   217,   217,
     217,   217,   217,   217,   217,   701,   217,   701,   217,   217,
     217,   255,   315,   701,   201,   202,   217,   202,   217,   700,
     217,   217,   217,   217,   217,   217,   700,   701,   217,   217,
     701,   217,   217,   217,   217,   217,   217,   701,   701,   700,
     701,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   701,
     217,   701,   217,   217,   701,   217,   217,   153,   217,   217,
     701,   217,   700,   700,   700,   700,   701,   700,   700,   103,
     187,   700,   701,   701,   701,   217,   700,   217,   217,   217,
     700,   217,   217,   701,   217,   700,   217,   701,   107,   687,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   217,   700,   700,   107,   687,   701,
     217,   107,   700,   107,   700,   700,   347,   700,   217,   701,
     701,   217,   217,   217,   217,   217,   701,   217,   217,   700,
     701,   701,   701,   217,   672,   701,   606,   217,   217,   700,
     701,   700,   701,   621,   700,   701,   217,   669,   701,   217,
     701,   701,   701,   701,   617,   700,   217,   366,   666,   217,
     217,   539,   700,   700,   217,   277,   217,   700,   701,   217,
     700,   701,   217,   700,   701,   217,   700,   701,   217,   701,
     701,   701,   701,   217,   700,   701,   701,   701,   701,   374,
     217,   286,   701,   217,   399,   217,   399,   217,   564,   700,
     217,   700,   701,   642,   701,   217,   211,   217,   217,   700,
     217,   217,   217,   217,   486,   700,   701,   701,   701,   217,
     701,   701,   701,   217,   217,   217,   135,   700,   107,   135,
     700,   107,   107,   701,   217,   701,   701,   701,   700,   701,
     319,   700,   217,   211,   217,   700,   217,   319,   701,   217,
     701,   701,   217,   700,   217,   217,   217,   217,   217,   217,
     700,   701,   701,   701,   319,   700,   701,   701,   701,   701,
      97,   188,   217,   323,   701,   701,   701,   217,   701,   701,
     217,   700,   217,   217,   666,   700,   217,   701,   217,   217,
     701,   217,   666,   700,   217,   701,   701,   701,   701,    95,
     217,   217,   275,   701,   217,   341,   700,   217,   217,   316,
     701,   700,   217,   316,   217,   701,   701,   217,   700,   217,
     376,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   217,   700,   217,   700,   217,   700,   217,   217,   217,
     202,   217,   700,   700,   202,   217,   217,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   701,   217,   217,   217,
     217,   217,   217,   217,   217,   217,   217,   217,   217,   217,
     217,   700,   701,   217,   217,   701,   217,   701,   217,   701,
     217,   217,   701,   701,   217,   701,   217,   701,   701,   701,
     701,   701,   701,   701,   701,   217,   701,   701,   701,   217,
     107,   700,   107,   700,   700,   319,   700,   217,   217,   701,
     701,   701,   217,   701,   701,   701,   606,   217,   217,   701,
     700,   701,   217,   701,   701,   701,   217,   217,   217,   217,
     701,   700,   217,   239,   241,   242,   243,   244,   245,   246,
     247,   248,   540,   539,   700,   217,   701,   217,   217,   700,
     701,   217,   701,   217,   217,   701,   217,   217,   700,   700,
     700,   700,   217,   701,   701,   701,   701,   701,   217,   217,
     701,   701,   701,   217,   217,   217,   217,   217,   217,   700,
     700,   700,   701,   701,   701,   701,   701,   217,   700,   700,
     700,   700,   700,   701,   700,   217,   700,   217,   701,   217,
     700,   701,   217,   217,   700,   217,   217,   217,   217,   701,
     217,   701,   701,   701,   700,   701,   217,   701,   701,   217,
     701,   217,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   217,   217,   666,   701,   701,   701,   217,   217,   666,
     701,   217,   700,   701,   701,   701,   217,   701,    95,   217,
     701,   217,   341,   342,   700,   217,   217,   316,   701,   217,
     217,   701,   217,   217,   376,   217,   217,   217,   701,   217,
     700,   700,   217,   217,   700,   217,   700,   217,   217,   701,
     701,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   701,   701,   701,   701,   701,   319,   700,   700,   701,
     701,   701,   217,   107,   701,   217,   217,   701,   217,   217,
     701,   217,   217,   217,   217,   666,   700,   701,   701,   539,
     539,   277,   701,   701,   701,   217,   701,   701,   701,   701,
     700,   217,   217,   281,   700,   700,   700,   700,   217,   217,
     217,   701,   217,   701,   217,   700,   701,   701,   701,   217,
     701,   701,   135,   700,   135,   700,   700,   217,   217,   700,
     217,   701,   701,   701,   701,   701,   700,   211,   217,   217,
     701,   701,   217,   217,   701,   701,   217,   701,   217,   701,
     701,   701,   701,   701,   701,   217,   217,   666,   701,   217,
     217,   217,   217,   666,   701,   217,   700,   701,   701,   217,
     217,   217,   217,   342,   701,   701,   217,   217,   217,   316,
     701,   217,   217,   700,   217,   700,   217,   217,   217,   217,
     700,   217,   217,   366,   368,   666,   701,   701,   217,   701,
     217,   701,   701,   217,   701,   217,   217,   701,   701,   701,
     701,   217,   701,   701,   701,   700,   701,   700,   217,   701,
     701,   217,   217,   217,   217,   217,   217,   701,   217,   217,
     701,   217,   701,   217,   217,   701,   700,   217,   217,   281,
     217,   217,   701,   217,   701,   701,   701,   217,   701,   700,
     700,   217,   701,   217,   217,   701,   701,   701,   701,   701,
     217,   701,   701,   217,   701,   217,   217,   701,   701,   701,
     701,   701,   701,   217,   217,   666,   700,   217,   217,   666,
     700,   217,   701,   701,   701,   217,    95,   217,   217,   217,
     217,   217,   217,   700,   701,   217,   701,   701,   701,   701,
     217,   701,   217,   701,   217,   701,   701,   217,   701,   701,
     701,   217,   701,   700,   217,   701,   217,   217,   107,   217,
     217,   217,   700,   217,   217,   701,   217,   701,   701,   701,
     701,   217,   701,   217,   217,   217,   211,   217,   701,   217,
     323,   701,   701,   217,   701,   701,   701,   701,   701,   217,
     701,   217,   701,   701,   701,    95,   217,   217,   217,   217,
     366,   217,   701,   217,   701,   217,   217,   701,   217,   701,
     701,   701,   217,   217,   217,   701,   701,   701,   217,   217,
     700,   700,   217,   701,   701,   217,   701,   701,   217,   701,
     701,   217,   323,   701,   217,   323,   701,   701,   701,   701,
     701,   217,   701,   217,   666,   217,   666,   701,   701,   701,
     217,   700,   217,   701,   217,   701,   217,   217,   701,   217,
     217,   701,   217,   701,   217,   217,   107,   700,   700,   486,
     701,   701,   701,   217,   701,   217,   701,   217,   323,   701,
     217,   323,   701,   700,   701,   701,   701,   701,   217,   217,
     217,   666,   701,   701,   701,   217,   666,   701,   701,   217,
     701,   217,   701,   701,   217,   700,   700,   701,   217,   217,
     701,   701,   217,   701,   217,   701,   217,   323,   217,   701,
     701,   701,   701,   217,   701,   701,   217,   701,   217,   217,
     217,   701,   217,   217,   217,   700,   217,   701,   217,   217,
     701,   217,   217,   701,   701,   701,   217,   701,   217,   666,
     701,   217,   701,   701,   701,   217,   701,   217,   701,   217,
     700,   217,   701,   217,   217,   666,   701,   701,   701,   217,
     701,   217,   700,    69,   217,   701,   217,   701,   217,   701,
     217,   700,   701,   701,   217,   701,   217,   217,   217,   270,
     291,   701,   701,   701,   701,   701,   701,   701,   701,   701,
     701,   217,   217,   701,   701,   217,   270,   291,   701,   701,
     701,   701,   701,   701,   701,   701,   700,   701,   217,   217,
     701,   217,   701,   700,   700,   701,   701,   701,   217
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

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
#ifndef	YYINITDEPTH
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
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
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
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

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
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
  int yytoken;
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

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

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
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
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
      if (yyn == 0 || yyn == YYTABLE_NINF)
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
  *++yyvsp = yylval;

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
     `$$ = $1'.

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

/* Line 1455 of yacc.c  */
#line 169 "p.y"
    { 
          return 0;
         ;}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 180 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; delete (yyvsp[(1) - (1)].bclist); ;}
    break;

  case 7:

/* Line 1455 of yacc.c  */
#line 182 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 191 "p.y"
    {;}
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 200 "p.y"
    {;}
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 202 "p.y"
    {;}
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 209 "p.y"
    {;}
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 229 "p.y"
    { domain->setMFTT((yyvsp[(1) - (1)].mftval).first, (yyvsp[(1) - (1)].mftval).second); ;}
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 231 "p.y"
    { domain->setHFTT((yyvsp[(1) - (1)].hftval).first, (yyvsp[(1) - (1)].hftval).second); ;}
    break;

  case 77:

/* Line 1455 of yacc.c  */
#line 259 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 78:

/* Line 1455 of yacc.c  */
#line 261 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 79:

/* Line 1455 of yacc.c  */
#line 263 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 265 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 267 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 269 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 287 "p.y"
    { if(domain->setComplexNeuman((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; ;}
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 290 "p.y"
    { if(domain->setComplexDirichlet((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; ;}
    break;

  case 111:

/* Line 1455 of yacc.c  */
#line 301 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 112:

/* Line 1455 of yacc.c  */
#line 303 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 311 "p.y"
    {;}
    break;

  case 128:

/* Line 1455 of yacc.c  */
#line 322 "p.y"
    {;}
    break;

  case 129:

/* Line 1455 of yacc.c  */
#line 324 "p.y"
    {;}
    break;

  case 130:

/* Line 1455 of yacc.c  */
#line 326 "p.y"
    {;}
    break;

  case 131:

/* Line 1455 of yacc.c  */
#line 328 "p.y"
    {;}
    break;

  case 132:

/* Line 1455 of yacc.c  */
#line 330 "p.y"
    {;}
    break;

  case 133:

/* Line 1455 of yacc.c  */
#line 332 "p.y"
    {;}
    break;

  case 134:

/* Line 1455 of yacc.c  */
#line 334 "p.y"
    {;}
    break;

  case 135:

/* Line 1455 of yacc.c  */
#line 336 "p.y"
    {;}
    break;

  case 147:

/* Line 1455 of yacc.c  */
#line 351 "p.y"
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[(3) - (5)].ival)); 
            domain->solInfo().nsample = (yyvsp[(4) - (5)].ival);
          ;}
    break;

  case 148:

/* Line 1455 of yacc.c  */
#line 358 "p.y"
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[(3) - (4)].ival));
          ;}
    break;

  case 150:

/* Line 1455 of yacc.c  */
#line 365 "p.y"
    { if ((yyvsp[(2) - (5)].ival) == OutputInfo::Attribute)  geoSource->setGroupAttribute((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival)-1);
          else if ((yyvsp[(2) - (5)].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[(2) - (5)].ival));  exit(-1); }
        ;}
    break;

  case 151:

/* Line 1455 of yacc.c  */
#line 370 "p.y"
    { int i;
          if ((yyvsp[(2) - (6)].ival) == OutputInfo::Attribute)  {
            for(i=(yyvsp[(3) - (6)].ival); i<(yyvsp[(4) - (6)].ival)+1; ++i)
              geoSource->setGroupAttribute(i-1,(yyvsp[(5) - (6)].ival)-1);
          }
          else if ((yyvsp[(2) - (6)].ival) == OutputInfo::Nodal)  {
            for(i=(yyvsp[(3) - (6)].ival); i<(yyvsp[(4) - (6)].ival)+1; ++i)
              geoSource->setNodeGroup(i-1, (yyvsp[(5) - (6)].ival));
          }
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[(2) - (6)].ival));  exit(-1); }
        ;}
    break;

  case 152:

/* Line 1455 of yacc.c  */
#line 382 "p.y"
    { if ((yyvsp[(2) - (6)].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[(2) - (6)].ival));  exit(-1); }
        ;}
    break;

  case 154:

/* Line 1455 of yacc.c  */
#line 389 "p.y"
    { geoSource->setGroupRandomProperty((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].rprop),(yyvsp[(4) - (6)].fval),(yyvsp[(5) - (6)].fval)); ;}
    break;

  case 155:

/* Line 1455 of yacc.c  */
#line 393 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (5)].fval)); ;}
    break;

  case 156:

/* Line 1455 of yacc.c  */
#line 397 "p.y"
    { domain->solInfo().curSweepParam = 0; domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies1(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); ;}
    break;

  case 157:

/* Line 1455 of yacc.c  */
#line 399 "p.y"
    { domain->solInfo().curSweepParam = 0; domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies2(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); ;}
    break;

  case 158:

/* Line 1455 of yacc.c  */
#line 401 "p.y"
    { domain->solInfo().curSweepParam = 0; domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (8)].fval)); domain->addFrequencies(2.0*PI*(yyvsp[(4) - (8)].fval), 2.0*PI*(yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].ival), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 159:

/* Line 1455 of yacc.c  */
#line 403 "p.y"
    {
          domain->solInfo().curSweepParam = 0;
          domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (13)].fval));
          domain->addFrequencies(2.0*PI*(yyvsp[(4) - (13)].fval), 2.0*PI*(yyvsp[(5) - (13)].fval), 2,(yyvsp[(6) - (13)].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = (yyvsp[(9) - (13)].ival);
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[(6) - (13)].ival);
          if ((yyvsp[(7) - (13)].ival) == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = false; 
          else 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = true;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[(4) - (13)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[(5) - (13)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.atol = (yyvsp[(8) - (13)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = (yyvsp[(10) - (13)].ival);
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = (yyvsp[(11) - (13)].ival);
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = (yyvsp[(12) - (13)].ival);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[(11) - (13)].ival);
        ;}
    break;

  case 160:

/* Line 1455 of yacc.c  */
#line 423 "p.y"
    {
          domain->solInfo().curSweepParam = 0;
          domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (8)].fval));
          domain->addFrequencies(2.0*PI*(yyvsp[(4) - (8)].fval), 2.0*PI*(yyvsp[(5) - (8)].fval), 2,(yyvsp[(6) - (8)].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = 6;
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[(6) - (8)].ival);
          if ((yyvsp[(5) - (8)].fval) == SweepParams::KrylovGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = false; 
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = true;
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 16;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[(4) - (8)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[(5) - (8)].fval);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
        ;}
    break;

  case 161:

/* Line 1455 of yacc.c  */
#line 449 "p.y"
    { domain->solInfo().curSweepParam = (yyvsp[(3) - (7)].ival); if ((yyvsp[(3) - (7)].ival) == 0) geoSource->setImpe((yyvsp[(6) - (7)].fval)); ;}
    break;

  case 162:

/* Line 1455 of yacc.c  */
#line 453 "p.y"
    { domain->setFrequencySet((yyvsp[(2) - (8)].ival)); domain->solInfo().curSweepParam = (yyvsp[(2) - (8)].ival); if ((yyvsp[(2) - (8)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (8)].fval)); domain->addFrequencies1(2.0*PI*(yyvsp[(5) - (8)].fval), 2.0*PI*(yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 163:

/* Line 1455 of yacc.c  */
#line 455 "p.y"
    { domain->setFrequencySet((yyvsp[(2) - (8)].ival)); domain->solInfo().curSweepParam = (yyvsp[(2) - (8)].ival); if ((yyvsp[(2) - (8)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (8)].fval)); domain->addFrequencies2(2.0*PI*(yyvsp[(5) - (8)].fval), 2.0*PI*(yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 164:

/* Line 1455 of yacc.c  */
#line 457 "p.y"
    { domain->setFrequencySet((yyvsp[(2) - (9)].ival)); domain->solInfo().curSweepParam = (yyvsp[(2) - (9)].ival); if ((yyvsp[(2) - (9)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (9)].fval)); domain->addFrequencies(2.0*PI*(yyvsp[(5) - (9)].fval), 2.0*PI*(yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].ival), (yyvsp[(8) - (9)].ival)); ;}
    break;

  case 165:

/* Line 1455 of yacc.c  */
#line 459 "p.y"
    {
          domain->setFrequencySet((yyvsp[(2) - (14)].ival));  domain->solInfo().curSweepParam = (yyvsp[(2) - (14)].ival);
          if ((yyvsp[(2) - (14)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (14)].fval));
          domain->addFrequencies(2.0*PI*(yyvsp[(5) - (14)].fval), 2.0*PI*(yyvsp[(6) - (14)].fval), 2,(yyvsp[(7) - (14)].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = (yyvsp[(10) - (14)].ival);
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[(7) - (14)].ival);
          if ((yyvsp[(8) - (14)].ival) == SweepParams::KrylovGalProjection) 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = false; 
          else 
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = true;
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[(5) - (14)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[(6) - (14)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.atol = (yyvsp[(9) - (14)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.minRHS = (yyvsp[(11) - (14)].ival);
          domain->solInfo().getSweepParams()->adaptSweep.maxRHS = (yyvsp[(12) - (14)].ival);
          domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = (yyvsp[(13) - (14)].ival);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[(12) - (14)].ival);
        ;}
    break;

  case 166:

/* Line 1455 of yacc.c  */
#line 479 "p.y"
    {
          domain->setFrequencySet((yyvsp[(2) - (9)].ival));  domain->solInfo().curSweepParam = (yyvsp[(2) - (9)].ival);
          if ((yyvsp[(2) - (9)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (9)].fval));
          domain->addFrequencies(2.0*PI*(yyvsp[(5) - (9)].fval), 2.0*PI*(yyvsp[(6) - (9)].fval), 2,(yyvsp[(7) - (9)].ival));
          domain->solInfo().getSweepParams()->isAdaptSweep = true;
          domain->solInfo().getSweepParams()->adaptSweep.maxP = 6;
          domain->solInfo().getSweepParams()->adaptSweep.numS = (yyvsp[(7) - (9)].ival);
          if ((yyvsp[(8) - (9)].ival) == SweepParams::KrylovGalProjection) {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = false; 
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 48;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          else {
             domain->solInfo().getSweepParams()->adaptSweep.dgp_flag = true;
             domain->solInfo().getSweepParams()->adaptSweep.atol = 1e-2;
             domain->solInfo().getSweepParams()->adaptSweep.minRHS = 8;
             domain->solInfo().getSweepParams()->adaptSweep.maxRHS = 16;
             domain->solInfo().getSweepParams()->adaptSweep.deltaRHS = 4;
          }
          domain->solInfo().getSweepParams()->adaptSweep.w1 = 2.0*PI*(yyvsp[(5) - (9)].fval);
          domain->solInfo().getSweepParams()->adaptSweep.w2 = 2.0*PI*(yyvsp[(6) - (9)].fval);
          domain->solInfo().getSweepParams()->nFreqSweepRHS = domain->solInfo().getSweepParams()->adaptSweep.maxRHS;
        ;}
    break;

  case 172:

/* Line 1455 of yacc.c  */
#line 512 "p.y"
    { domain->solInfo().getSweepParams()->pade_pivot = true; domain->solInfo().getSweepParams()->pade_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 173:

/* Line 1455 of yacc.c  */
#line 516 "p.y"
    { domain->solInfo().getSweepParams()->pade_poles = true; ;}
    break;

  case 174:

/* Line 1455 of yacc.c  */
#line 518 "p.y"
    { domain->solInfo().getSweepParams()->pade_poles = true; 
          domain->solInfo().getSweepParams()->pade_poles_sigmaL = (yyvsp[(2) - (4)].fval); domain->solInfo().getSweepParams()->pade_poles_sigmaU = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 175:

/* Line 1455 of yacc.c  */
#line 523 "p.y"
    { geoSource->setImpe((yyvsp[(2) - (3)].fval)); domain->addCoarseFrequency(2.0*PI*(yyvsp[(2) - (3)].fval)); ;}
    break;

  case 176:

/* Line 1455 of yacc.c  */
#line 526 "p.y"
    { domain->addFrequencies(2.0*PI*(yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].ival)); ;}
    break;

  case 177:

/* Line 1455 of yacc.c  */
#line 530 "p.y"
    { domain->solInfo().getSweepParams()->freqSweepMethod = (yyvsp[(2) - (4)].ival); 
          int &l = domain->solInfo().getSweepParams()->padeL, &m = domain->solInfo().getSweepParams()->padeM, &n = domain->solInfo().getSweepParams()->padeN;
          switch((yyvsp[(2) - (4)].ival)) {
            case SweepParams::Taylor:
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[(3) - (4)].ival)+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = (yyvsp[(3) - (4)].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = (yyvsp[(3) - (4)].ival);
              if(m%n != 0) m = m/n*(n+1)-m%n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = (yyvsp[(3) - (4)].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = (yyvsp[(3) - (4)].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::QRGalProjection:
              n = (yyvsp[(3) - (4)].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        ;}
    break;

  case 178:

/* Line 1455 of yacc.c  */
#line 569 "p.y"
    { domain->solInfo().getSweepParams()->freqSweepMethod = (yyvsp[(2) - (6)].ival);
          int &l = domain->solInfo().getSweepParams()->padeL, &m = domain->solInfo().getSweepParams()->padeM, &n = domain->solInfo().getSweepParams()->padeN;
          switch((yyvsp[(2) - (6)].ival)) {
            case SweepParams::Taylor:
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (yyvsp[(3) - (6)].ival)+1; // taylor
              break;
            case SweepParams::Pade1:
              n = 1;
              l = (yyvsp[(4) - (6)].ival);
              m = (yyvsp[(5) - (6)].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+m+1;
              break;
            case SweepParams::Pade:
            case SweepParams::Fourier:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival); 
              m = (yyvsp[(5) - (6)].ival);
              domain->solInfo().getSweepParams()->nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SweepParams::PadeLanczos:
              n = (yyvsp[(3) - (6)].ival);
              m = (yyvsp[(5) - (6)].ival);
              if(m%n != 0) m = m/n*(n+1)-m%n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = m/n;
              break;
            case SweepParams::GalProjection:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::KrylovGalProjection:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
            case SweepParams::QRGalProjection:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival);
              m = 1;
              domain->solInfo().getSweepParams()->nFreqSweepRHS = l+1;
              break;
          }
        ;}
    break;

  case 180:

/* Line 1455 of yacc.c  */
#line 619 "p.y"
    { geoSource->binaryInput = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 181:

/* Line 1455 of yacc.c  */
#line 621 "p.y"
    { geoSource->binaryInput = bool((yyvsp[(3) - (5)].ival));
            std::string prefix = (yyvsp[(4) - (5)].strval);
            clusterData_ = prefix + ".msh";
            decomposition_ = prefix + ".dec";
            connectivity_ = prefix + ".con";
            subdomains_ = prefix + ".sub";
          ;}
    break;

  case 182:

/* Line 1455 of yacc.c  */
#line 629 "p.y"
    { geoSource->binaryOutput = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 183:

/* Line 1455 of yacc.c  */
#line 631 "p.y"
    { geoSource->binaryOutput = bool((yyvsp[(3) - (5)].ival));
            int len = strlen((yyvsp[(4) - (5)].strval));
            char *file = new char[len+5];
            strcpy(file, (yyvsp[(4) - (5)].strval));
            strcat(file,".con");
            geoSource->setGlob(file);
          ;}
    break;

  case 184:

/* Line 1455 of yacc.c  */
#line 639 "p.y"
    { geoSource->setGeo((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 185:

/* Line 1455 of yacc.c  */
#line 641 "p.y"
    { geoSource->setDecomp((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 186:

/* Line 1455 of yacc.c  */
#line 643 "p.y"
    { geoSource->setGlob((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 187:

/* Line 1455 of yacc.c  */
#line 645 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 188:

/* Line 1455 of yacc.c  */
#line 647 "p.y"
    { geoSource->setCpuMap((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 189:

/* Line 1455 of yacc.c  */
#line 651 "p.y"
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[(2) - (3)].ival)); 
#endif
	;}
    break;

  case 190:

/* Line 1455 of yacc.c  */
#line 659 "p.y"
    {if(decInit==0) decInit = new DecInit(); ;}
    break;

  case 191:

/* Line 1455 of yacc.c  */
#line 661 "p.y"
    {decInit->file = strdup((yyvsp[(3) - (4)].strval));;}
    break;

  case 192:

/* Line 1455 of yacc.c  */
#line 663 "p.y"
    {decInit->nsubs = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 193:

/* Line 1455 of yacc.c  */
#line 665 "p.y"
    {decInit->weight = true; ;}
    break;

  case 194:

/* Line 1455 of yacc.c  */
#line 667 "p.y"
    {decInit->memory = true; ;}
    break;

  case 195:

/* Line 1455 of yacc.c  */
#line 669 "p.y"
    {decInit->exitAfterDec = true;;}
    break;

  case 196:

/* Line 1455 of yacc.c  */
#line 671 "p.y"
    {decInit->skip = true;;}
    break;

  case 197:

/* Line 1455 of yacc.c  */
#line 673 "p.y"
    {decInit->nosa = true; ;}
    break;

  case 198:

/* Line 1455 of yacc.c  */
#line 675 "p.y"
    {decInit->trivial = true; ;}
    break;

  case 199:

/* Line 1455 of yacc.c  */
#line 677 "p.y"
    {decInit->fsgl = true; ;}
    break;

  case 200:

/* Line 1455 of yacc.c  */
#line 681 "p.y"
    {;}
    break;

  case 201:

/* Line 1455 of yacc.c  */
#line 683 "p.y"
    { weightList[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 202:

/* Line 1455 of yacc.c  */
#line 687 "p.y"
    {;}
    break;

  case 203:

/* Line 1455 of yacc.c  */
#line 689 "p.y"
    { fieldWeightList[(int)Element::Acoustic] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 204:

/* Line 1455 of yacc.c  */
#line 691 "p.y"
    { fieldWeightList[(int)Element::Structural] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 205:

/* Line 1455 of yacc.c  */
#line 693 "p.y"
    { fieldWeightList[(int)Element::Thermal] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 206:

/* Line 1455 of yacc.c  */
#line 695 "p.y"
    { fieldWeightList[(int)Element::Fluid] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 207:

/* Line 1455 of yacc.c  */
#line 698 "p.y"
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = 0; ;}
    break;

  case 208:

/* Line 1455 of yacc.c  */
#line 700 "p.y"
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 209:

/* Line 1455 of yacc.c  */
#line 702 "p.y"
    { (yyval.mftval).first->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 210:

/* Line 1455 of yacc.c  */
#line 706 "p.y"
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = 0; ;}
    break;

  case 211:

/* Line 1455 of yacc.c  */
#line 708 "p.y"
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 212:

/* Line 1455 of yacc.c  */
#line 710 "p.y"
    { (yyval.hftval).first->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 213:

/* Line 1455 of yacc.c  */
#line 714 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 214:

/* Line 1455 of yacc.c  */
#line 716 "p.y"
    { (yyval.ival) = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 215:

/* Line 1455 of yacc.c  */
#line 718 "p.y"
    { domain->setLoadFactor((yyval.ival), (yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 216:

/* Line 1455 of yacc.c  */
#line 720 "p.y"
    { domain->setLoadFactorMFTT((yyval.ival), (yyvsp[(2) - (5)].ival), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 217:

/* Line 1455 of yacc.c  */
#line 722 "p.y"
    { domain->setLoadFactorHFTT((yyval.ival), (yyvsp[(2) - (5)].ival), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 225:

/* Line 1455 of yacc.c  */
#line 735 "p.y"
    { geoSource->addCFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 226:

/* Line 1455 of yacc.c  */
#line 739 "p.y"
    { geoSource->addCoefInfo((yyvsp[(2) - (4)].ival)-1,(yyvsp[(4) - (4)].coefdata)); ;}
    break;

  case 227:

/* Line 1455 of yacc.c  */
#line 741 "p.y"
    { (yyvsp[(10) - (10)].coefdata).c[6][0] = (yyvsp[(3) - (10)].fval);
          (yyvsp[(10) - (10)].coefdata).c[6][1] = (yyvsp[(4) - (10)].fval);
          (yyvsp[(10) - (10)].coefdata).c[6][2] = (yyvsp[(5) - (10)].fval);
          (yyvsp[(10) - (10)].coefdata).c[6][3] = (yyvsp[(6) - (10)].fval);
          (yyvsp[(10) - (10)].coefdata).c[6][4] = (yyvsp[(7) - (10)].fval);
          (yyvsp[(10) - (10)].coefdata).c[6][5] = (yyvsp[(8) - (10)].fval);
          geoSource->addCoefInfo((yyvsp[(2) - (10)].ival)-1,(yyvsp[(10) - (10)].coefdata)); ;}
    break;

  case 228:

/* Line 1455 of yacc.c  */
#line 751 "p.y"
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[(1) - (4)].ival)-1,(yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 229:

/* Line 1455 of yacc.c  */
#line 753 "p.y"
    { (yyval.coefdata).setCoef((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 230:

/* Line 1455 of yacc.c  */
#line 757 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 231:

/* Line 1455 of yacc.c  */
#line 759 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 232:

/* Line 1455 of yacc.c  */
#line 763 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 233:

/* Line 1455 of yacc.c  */
#line 765 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 234:

/* Line 1455 of yacc.c  */
#line 769 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 235:

/* Line 1455 of yacc.c  */
#line 771 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 236:

/* Line 1455 of yacc.c  */
#line 775 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 237:

/* Line 1455 of yacc.c  */
#line 777 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 238:

/* Line 1455 of yacc.c  */
#line 781 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (11)].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
	  (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval);
	  (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval);
          (yyval.ldata).d[9] = 0;  (yyval.ldata).d[10] = 0; (yyval.ldata).d[11] = 0; ;}
    break;

  case 239:

/* Line 1455 of yacc.c  */
#line 788 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (13)].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (13)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (13)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (13)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (13)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (13)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (13)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (13)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (13)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (13)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (13)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (13)].fval);(yyval.ldata).d[11] = 0; ;}
    break;

  case 240:

/* Line 1455 of yacc.c  */
#line 795 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (14)].ival)-1;
          (yyval.ldata).matid = -1; // this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (14)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (14)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (14)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (14)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (14)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (14)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (14)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (14)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (14)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (14)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (14)].fval); (yyval.ldata).d[11] = (yyvsp[(13) - (14)].fval); ;}
    break;

  case 241:

/* Line 1455 of yacc.c  */
#line 804 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (5)].ival)-1;  (yyval.ldata).matid = (yyvsp[(2) - (5)].ival)-1; (yyval.ldata).d[7] = (yyvsp[(3) - (5)].fval); (yyval.ldata).d[8] = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 243:

/* Line 1455 of yacc.c  */
#line 809 "p.y"
    { geoSource->addLayMat((yyvsp[(2) - (2)].ldata).matid, (yyvsp[(2) - (2)].ldata).d); ;}
    break;

  case 244:

/* Line 1455 of yacc.c  */
#line 815 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (7)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (7)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (7)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (7)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (7)].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[(6) - (7)].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; ;}
    break;

  case 245:

/* Line 1455 of yacc.c  */
#line 820 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (9)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (9)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (9)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (9)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (9)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (9)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (9)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (9)].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; ;}
    break;

  case 246:

/* Line 1455 of yacc.c  */
#line 825 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (11)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval);
          (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); (yyval.ldata).d[9] = 0; ;}
    break;

  case 247:

/* Line 1455 of yacc.c  */
#line 829 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (12)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (12)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (12)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (12)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (12)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (12)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (12)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (12)].fval); 
          (yyval.ldata).d[7] = (yyvsp[(9) - (12)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (12)].fval); (yyval.ldata).d[9] = (yyvsp[(11) - (12)].fval); ;}
    break;

  case 249:

/* Line 1455 of yacc.c  */
#line 836 "p.y"
    { domain->addDMass((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 250:

/* Line 1455 of yacc.c  */
#line 838 "p.y"
    { domain->addDMass((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(5) - (6)].fval),(yyvsp[(4) - (6)].ival)-1); ;}
    break;

  case 252:

/* Line 1455 of yacc.c  */
#line 843 "p.y"
    { domain->setGravity((yyvsp[(2) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 254:

/* Line 1455 of yacc.c  */
#line 848 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (5)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (5)].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[(4) - (5)].strval); ;}
    break;

  case 255:

/* Line 1455 of yacc.c  */
#line 852 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (4)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (4)].strval);;}
    break;

  case 256:

/* Line 1455 of yacc.c  */
#line 855 "p.y"
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[(2) - (4)].strval);
          domain->solInfo().nRestart = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 257:

/* Line 1455 of yacc.c  */
#line 860 "p.y"
    { geoSource->setControlFile((yyvsp[(2) - (3)].strval));
         geoSource->setControlRoutine((char *) "controlObj");;}
    break;

  case 258:

/* Line 1455 of yacc.c  */
#line 865 "p.y"
    { geoSource->setControlRoutine((yyvsp[(2) - (3)].strval)); ;}
    break;

  case 259:

/* Line 1455 of yacc.c  */
#line 869 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 260:

/* Line 1455 of yacc.c  */
#line 874 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Actuators; }
          if(geoSource->setActuatorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)            return -1; ;}
    break;

  case 261:

/* Line 1455 of yacc.c  */
#line 880 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdf; }
          if(geoSource->setUsdfLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)       return -1; ;}
    break;

  case 262:

/* Line 1455 of yacc.c  */
#line 887 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdd;
          if(geoSource->setUsddLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)    return -1; ;}
    break;

  case 263:

/* Line 1455 of yacc.c  */
#line 894 "p.y"
    { 
    domain->solInfo().sensitivity = true;
    domain->senInfo = new SensitivityInfo[50];  // maximum number of sensitivities are fixed to 50
  ;}
    break;

  case 264:

/* Line 1455 of yacc.c  */
#line 899 "p.y"
    { domain->addSensitivity((yyvsp[(2) - (3)].sinfo)); ;}
    break;

  case 265:

/* Line 1455 of yacc.c  */
#line 903 "p.y"
    { (yyval.sinfo).initialize(); (yyval.sinfo).type = (SensitivityInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.sinfo).method = (SensitivityInfo::Method) (yyvsp[(2) - (3)].ival); (yyval.sinfo).numParam = (yyvsp[(3) - (3)].ival);  
  ;}
    break;

  case 266:

/* Line 1455 of yacc.c  */
#line 906 "p.y"
    { (yyval.sinfo).surface = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 267:

/* Line 1455 of yacc.c  */
#line 910 "p.y"
    { numColumns = 3; ;}
    break;

  case 268:

/* Line 1455 of yacc.c  */
#line 912 "p.y"
    { numColumns = 6; ;}
    break;

  case 269:

/* Line 1455 of yacc.c  */
#line 914 "p.y"
    { numColumns = 3; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 270:

/* Line 1455 of yacc.c  */
#line 916 "p.y"
    { numColumns = 6; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 271:

/* Line 1455 of yacc.c  */
#line 918 "p.y"
    { numColumns = 3; domain->outFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 272:

/* Line 1455 of yacc.c  */
#line 920 "p.y"
    { numColumns = 6; domain->outFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 273:

/* Line 1455 of yacc.c  */
#line 922 "p.y"
    { numColumns = 3; domain->outFlag = (yyvsp[(2) - (4)].ival); geoSource->setOutLimit((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 274:

/* Line 1455 of yacc.c  */
#line 924 "p.y"
    { numColumns = 6; domain->outFlag = (yyvsp[(2) - (4)].ival); geoSource->setOutLimit((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 275:

/* Line 1455 of yacc.c  */
#line 926 "p.y"
    { (yyvsp[(2) - (3)].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[(2) - (3)].oinfo)); ;}
    break;

  case 276:

/* Line 1455 of yacc.c  */
#line 930 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 277:

/* Line 1455 of yacc.c  */
#line 932 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); ;}
    break;

  case 278:

/* Line 1455 of yacc.c  */
#line 934 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; ;}
    break;

  case 279:

/* Line 1455 of yacc.c  */
#line 936 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); 
          if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1;;}
    break;

  case 280:

/* Line 1455 of yacc.c  */
#line 939 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; ;}
    break;

  case 281:

/* Line 1455 of yacc.c  */
#line 941 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; ;}
    break;

  case 282:

/* Line 1455 of yacc.c  */
#line 943 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); (yyval.oinfo).sentype = 1; ;}
    break;

  case 283:

/* Line 1455 of yacc.c  */
#line 945 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); (yyval.oinfo).sentype = 1; ;}
    break;

  case 284:

/* Line 1455 of yacc.c  */
#line 947 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 285:

/* Line 1455 of yacc.c  */
#line 949 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); ;}
    break;

  case 286:

/* Line 1455 of yacc.c  */
#line 951 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; ;}
    break;

  case 287:

/* Line 1455 of yacc.c  */
#line 953 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1; ;}
    break;

  case 288:

/* Line 1455 of yacc.c  */
#line 955 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; ;}
    break;

  case 289:

/* Line 1455 of yacc.c  */
#line 957 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; ;}
    break;

  case 290:

/* Line 1455 of yacc.c  */
#line 960 "p.y"
    { (yyval.oinfo).nodeNumber = (yyvsp[(3) - (3)].ival)-1; ;}
    break;

  case 291:

/* Line 1455 of yacc.c  */
#line 962 "p.y"
    { (yyval.oinfo).surface = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 292:

/* Line 1455 of yacc.c  */
#line 964 "p.y"
    { (yyval.oinfo).ylayer = (yyvsp[(2) - (3)].fval); (yyval.oinfo).zlayer = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 293:

/* Line 1455 of yacc.c  */
#line 966 "p.y"
    { (yyval.oinfo).averageFlg = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 294:

/* Line 1455 of yacc.c  */
#line 968 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 295:

/* Line 1455 of yacc.c  */
#line 970 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (3)].ival); (yyval.oinfo).ncomplexout = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 296:

/* Line 1455 of yacc.c  */
#line 972 "p.y"
    { (yyval.oinfo).angularouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 297:

/* Line 1455 of yacc.c  */
#line 974 "p.y"
    { (yyval.oinfo).rotvecouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 298:

/* Line 1455 of yacc.c  */
#line 976 "p.y"
    { (yyval.oinfo).rotvecouttype = OutputInfo::Linear; ;}
    break;

  case 299:

/* Line 1455 of yacc.c  */
#line 978 "p.y"
    { (yyval.oinfo).rescaling = bool((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 300:

/* Line 1455 of yacc.c  */
#line 980 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 301:

/* Line 1455 of yacc.c  */
#line 982 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (3)].ival); sfem->setnsamp_out((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 302:

/* Line 1455 of yacc.c  */
#line 984 "p.y"
    { (yyval.oinfo).oframe = (OutputInfo::FrameType) (yyvsp[(2) - (2)].ival); ;}
    break;

  case 303:

/* Line 1455 of yacc.c  */
#line 986 "p.y"
    { (yyval.oinfo).matlab = true; ;}
    break;

  case 304:

/* Line 1455 of yacc.c  */
#line 988 "p.y"
    { domain->solInfo().xmatrixname = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 305:

/* Line 1455 of yacc.c  */
#line 990 "p.y"
    { domain->solInfo().qmatrixname = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 306:

/* Line 1455 of yacc.c  */
#line 992 "p.y"
    { domain->solInfo().rmatrixname = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 308:

/* Line 1455 of yacc.c  */
#line 997 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; ;}
    break;

  case 309:

/* Line 1455 of yacc.c  */
#line 1000 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[(2) - (3)].ival);;}
    break;

  case 310:

/* Line 1455 of yacc.c  */
#line 1003 "p.y"
    { domain->solInfo().qrfactorization = (yyvsp[(2) - (3)].ival);;}
    break;

  case 311:

/* Line 1455 of yacc.c  */
#line 1005 "p.y"
    { domain->solInfo().nEig = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 312:

/* Line 1455 of yacc.c  */
#line 1007 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;;}
    break;

  case 313:

/* Line 1455 of yacc.c  */
#line 1009 "p.y"
    { domain->solInfo().setSubSpaceInfo((yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 314:

/* Line 1455 of yacc.c  */
#line 1011 "p.y"
    { domain->solInfo().subspaceSize = (yyvsp[(2) - (3)].ival);;}
    break;

  case 315:

/* Line 1455 of yacc.c  */
#line 1013 "p.y"
    { domain->solInfo().tolEig = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 316:

/* Line 1455 of yacc.c  */
#line 1015 "p.y"
    { domain->solInfo().tolJac = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 317:

/* Line 1455 of yacc.c  */
#line 1017 "p.y"
    { domain->solInfo().explicitK = true; ;}
    break;

  case 318:

/* Line 1455 of yacc.c  */
#line 1019 "p.y"
    { geoSource->setShift((yyvsp[(2) - (3)].fval)); ;}
    break;

  case 319:

/* Line 1455 of yacc.c  */
#line 1021 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; ;}
    break;

  case 320:

/* Line 1455 of yacc.c  */
#line 1023 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (3)].strval); ;}
    break;

  case 321:

/* Line 1455 of yacc.c  */
#line 1026 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (4)].strval); 
          domain->solInfo().arpack_mode = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 322:

/* Line 1455 of yacc.c  */
#line 1030 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[(2) - (4)].fval), int((yyvsp[(3) - (4)].fval))); ;}
    break;

  case 323:

/* Line 1455 of yacc.c  */
#line 1033 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival));;}
    break;

  case 324:

/* Line 1455 of yacc.c  */
#line 1036 "p.y"
    { domain->solInfo().filtereig = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 325:

/* Line 1455 of yacc.c  */
#line 1038 "p.y"
    { domain->solInfo().eigenSolverSubType = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 326:

/* Line 1455 of yacc.c  */
#line 1040 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;;}
    break;

  case 327:

/* Line 1455 of yacc.c  */
#line 1043 "p.y"
    { domain->solInfo().maxitEig = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 328:

/* Line 1455 of yacc.c  */
#line 1045 "p.y"
    { domain->solInfo().test_ulrich = true; ;}
    break;

  case 329:

/* Line 1455 of yacc.c  */
#line 1047 "p.y"
    { domain->solInfo().addedMass = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 330:

/* Line 1455 of yacc.c  */
#line 1051 "p.y"
    { domain->solInfo().sloshing = 1; ;}
    break;

  case 331:

/* Line 1455 of yacc.c  */
#line 1053 "p.y"
    { domain->setGravitySloshing((yyvsp[(2) - (3)].fval)); ;}
    break;

  case 332:

/* Line 1455 of yacc.c  */
#line 1057 "p.y"
    { domain->solInfo().massFlag = 1; ;}
    break;

  case 333:

/* Line 1455 of yacc.c  */
#line 1061 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 334:

/* Line 1455 of yacc.c  */
#line 1064 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);;}
    break;

  case 335:

/* Line 1455 of yacc.c  */
#line 1068 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Top); ;}
    break;

  case 336:

/* Line 1455 of yacc.c  */
#line 1072 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); ;}
    break;

  case 340:

/* Line 1455 of yacc.c  */
#line 1077 "p.y"
    { domain->solInfo().modal = true; ;}
    break;

  case 341:

/* Line 1455 of yacc.c  */
#line 1079 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 342:

/* Line 1455 of yacc.c  */
#line 1081 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 343:

/* Line 1455 of yacc.c  */
#line 1083 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (8)].ival);
          domain->solInfo().stable_cfl = (yyvsp[(4) - (8)].fval);
          domain->solInfo().stable_tol = (yyvsp[(5) - (8)].fval);
          domain->solInfo().stable_maxit = (yyvsp[(6) - (8)].ival);
          domain->solInfo().stable_freq = (yyvsp[(7) - (8)].ival);
        ;}
    break;

  case 344:

/* Line 1455 of yacc.c  */
#line 1090 "p.y"
    { domain->solInfo().iacc_switch = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 345:

/* Line 1455 of yacc.c  */
#line 1092 "p.y"
    { domain->solInfo().zeroRot = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 346:

/* Line 1455 of yacc.c  */
#line 1094 "p.y"
    { domain->solInfo().no_secondary = true; ;}
    break;

  case 347:

/* Line 1455 of yacc.c  */
#line 1096 "p.y"
    { domain->solInfo().tdenforceFlag = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 348:

/* Line 1455 of yacc.c  */
#line 1098 "p.y"
    { domain->solInfo().check_energy_balance = true; ;}
    break;

  case 349:

/* Line 1455 of yacc.c  */
#line 1100 "p.y"
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[(3) - (5)].fval); 
          domain->solInfo().epsilon2 = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 350:

/* Line 1455 of yacc.c  */
#line 1106 "p.y"
    { domain->solInfo().ConwepOnOff = true;
          BlastLoading::InputFileData = (yyvsp[(3) - (4)].blastData); ;}
    break;

  case 351:

/* Line 1455 of yacc.c  */
#line 1111 "p.y"
    { // Note: chargeWeight must be entered in the units of mass of the problem, not units of force.
          (yyval.blastData).ExplosivePosition[0] = (yyvsp[(1) - (5)].fval);
          (yyval.blastData).ExplosivePosition[1] = (yyvsp[(2) - (5)].fval);
          (yyval.blastData).ExplosivePosition[2] = (yyvsp[(3) - (5)].fval);
          (yyval.blastData).ExplosiveDetonationTime = (yyvsp[(5) - (5)].fval);
          (yyval.blastData).BlastType = BlastLoading::BlastData::AirBurst; // ($5 == 0 ? BlastLoading::BlastData::SurfaceBurst : BlastLoading::BlastData::AirBurst);
          (yyval.blastData).ScaleLength = 1.0;
          (yyval.blastData).ScaleTime = 1.0;
          (yyval.blastData).ScaleMass = 1.0;
          (yyval.blastData).ExplosiveWeight = (yyvsp[(4) - (5)].fval)*2.2; // The 2.2 factor is to convert from kilograms to pounds force.
          (yyval.blastData).ExplosiveWeightCubeRoot = pow((yyval.blastData).ExplosiveWeight,1.0/3.0);
        ;}
    break;

  case 352:

/* Line 1455 of yacc.c  */
#line 1126 "p.y"
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; ;}
    break;

  case 354:

/* Line 1455 of yacc.c  */
#line 1129 "p.y"
    { domain->solInfo().acoustic = true; ;}
    break;

  case 356:

/* Line 1455 of yacc.c  */
#line 1134 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (2)].fval),(yyvsp[(2) - (2)].fval)); ;}
    break;

  case 357:

/* Line 1455 of yacc.c  */
#line 1136 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (4)].fval),(yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval),(yyvsp[(4) - (4)].fval)); ;}
    break;

  case 358:

/* Line 1455 of yacc.c  */
#line 1138 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[(1) - (1)].fval)); ;}
    break;

  case 359:

/* Line 1455 of yacc.c  */
#line 1140 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (3)].fval),(yyvsp[(2) - (3)].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 360:

/* Line 1455 of yacc.c  */
#line 1146 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
            domain->solInfo().order = 1;
          }
          else 
            domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo((yyvsp[(1) - (1)].fval)); 
        ;}
    break;

  case 361:

/* Line 1455 of yacc.c  */
#line 1157 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; ;}
    break;

  case 364:

/* Line 1455 of yacc.c  */
#line 1162 "p.y"
    { domain->solInfo().modal = true; ;}
    break;

  case 365:

/* Line 1455 of yacc.c  */
#line 1166 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (5)].fval), 0, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 366:

/* Line 1455 of yacc.c  */
#line 1168 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), 0, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].ival), (yyvsp[(5) - (6)].fval)); ;}
    break;

  case 367:

/* Line 1455 of yacc.c  */
#line 1172 "p.y"
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].ival)); ;}
    break;

  case 368:

/* Line 1455 of yacc.c  */
#line 1182 "p.y"
    { domain->solInfo().setAero((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCollocated = 0; ;}
    break;

  case 369:

/* Line 1455 of yacc.c  */
#line 1185 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().isCollocated = 0;
          if((yyvsp[(3) - (4)].ival) == 20) { // set default alphas for C0
            domain->solInfo().alphas[0] = 0.5+0.375;
            domain->solInfo().alphas[1] = -0.375;
          }
        ;}
    break;

  case 370:

/* Line 1455 of yacc.c  */
#line 1193 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (6)].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[(3) - (6)].ival) == 8) {
            // MPP uses only the first of the two inputted alphas
            domain->solInfo().mppFactor = (yyvsp[(4) - (6)].fval);
          }
          else {
            // These alphas are used in FlExchanger::sendDisplacements and DistFlExchanger::sendDisplacements
            // As of 4/14/2014 the following schemes can use the displacement predictor on the structure side:
            // A0, A4, A5, A6, A7 and C0. Furthermore, we now apply a separate "anti-predictor" for A6 and A7
            // to compensate for the legacy predictor on the fluid side (in MatchNodeSet::getDisplacement)
            domain->solInfo().alphas[0] = (yyvsp[(4) - (6)].fval)+(yyvsp[(5) - (6)].fval);
            domain->solInfo().alphas[1] = -(yyvsp[(5) - (6)].fval);
          }
        ;}
    break;

  case 371:

/* Line 1455 of yacc.c  */
#line 1209 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (5)].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[(4) - (5)].fval);
        ;}
    break;

  case 372:

/* Line 1455 of yacc.c  */
#line 1214 "p.y"
    { domain->solInfo().isCollocated = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 373:

/* Line 1455 of yacc.c  */
#line 1216 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 374:

/* Line 1455 of yacc.c  */
#line 1218 "p.y"
    {;}
    break;

  case 375:

/* Line 1455 of yacc.c  */
#line 1222 "p.y"
    {;}
    break;

  case 376:

/* Line 1455 of yacc.c  */
#line 1224 "p.y"
    { domain->AddAeroEmbedSurfaceId((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 377:

/* Line 1455 of yacc.c  */
#line 1228 "p.y"
    { domain->solInfo().setAeroHeat((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); ;}
    break;

  case 378:

/* Line 1455 of yacc.c  */
#line 1232 "p.y"
    { domain->solInfo().setThermoh(1); ;}
    break;

  case 379:

/* Line 1455 of yacc.c  */
#line 1236 "p.y"
    { domain->solInfo().setThermoe(1); ;}
    break;

  case 380:

/* Line 1455 of yacc.c  */
#line 1240 "p.y"
    { domain->solInfo().setModeDecomp(1); ;}
    break;

  case 381:

/* Line 1455 of yacc.c  */
#line 1244 "p.y"
    { domain->solInfo().hzemFlag=1; ;}
    break;

  case 382:

/* Line 1455 of yacc.c  */
#line 1248 "p.y"
    { domain->solInfo().slzemFlag=1; ;}
    break;

  case 383:

/* Line 1455 of yacc.c  */
#line 1252 "p.y"
    { domain->solInfo().setTrbm((yyvsp[(3) - (4)].fval)); ;}
    break;

  case 384:

/* Line 1455 of yacc.c  */
#line 1256 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); 
          filePrint(stderr," ... Using Geometric RBM Method     ...\n"); ;}
    break;

  case 385:

/* Line 1455 of yacc.c  */
#line 1259 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (4)].fval)); 
          filePrint(stderr," ... Using Geometric RBM Method     ...\n"); ;}
    break;

  case 386:

/* Line 1455 of yacc.c  */
#line 1262 "p.y"
    { domain->solInfo().setGrbm();
          filePrint(stderr," ... Using Geometric RBM Method     ...\n"); ;}
    break;

  case 387:

/* Line 1455 of yacc.c  */
#line 1265 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (6)].fval),(yyvsp[(4) - (6)].fval));
          domain->solInfo().grbm_use_lmpc = bool((yyvsp[(5) - (6)].ival));
          filePrint(stderr," ... Using Geometric RBM Method     ...\n"); ;}
    break;

  case 388:

/* Line 1455 of yacc.c  */
#line 1271 "p.y"
    { domain->solInfo().modeFilterFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 389:

/* Line 1455 of yacc.c  */
#line 1273 "p.y"
    { domain->solInfo().modeFilterFlag = 1; ;}
    break;

  case 390:

/* Line 1455 of yacc.c  */
#line 1277 "p.y"
    { domain->solInfo().useRbmFilter((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 391:

/* Line 1455 of yacc.c  */
#line 1279 "p.y"
    { domain->solInfo().useRbmFilter(1); ;}
    break;

  case 393:

/* Line 1455 of yacc.c  */
#line 1284 "p.y"
    { if((yyvsp[(1) - (1)].ival) < 1 || (yyvsp[(1) - (1)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(1) - (1)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(1) - (1)].ival)-1] = 1;
    ;}
    break;

  case 394:

/* Line 1455 of yacc.c  */
#line 1292 "p.y"
    { if((yyvsp[(2) - (2)].ival) < 1 || (yyvsp[(2) - (2)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(2) - (2)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(2) - (2)].ival)-1] = 1;
    ;}
    break;

  case 395:

/* Line 1455 of yacc.c  */
#line 1302 "p.y"
    { domain->solInfo().hzemFilterFlag=1; ;}
    break;

  case 396:

/* Line 1455 of yacc.c  */
#line 1306 "p.y"
    { domain->solInfo().slzemFilterFlag=1; ;}
    break;

  case 397:

/* Line 1455 of yacc.c  */
#line 1310 "p.y"
    { domain->solInfo().setTimes((yyvsp[(4) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(2) - (5)].fval)); ;}
    break;

  case 398:

/* Line 1455 of yacc.c  */
#line 1314 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (5)].ival),(yyvsp[(4) - (5)].ival),1);
        ;}
    break;

  case 399:

/* Line 1455 of yacc.c  */
#line 1320 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].ival));
        ;}
    break;

  case 400:

/* Line 1455 of yacc.c  */
#line 1325 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (7)].ival),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        ;}
    break;

  case 403:

/* Line 1455 of yacc.c  */
#line 1338 "p.y"
    { domain->solInfo().pitaNoForce = true; ;}
    break;

  case 404:

/* Line 1455 of yacc.c  */
#line 1340 "p.y"
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 405:

/* Line 1455 of yacc.c  */
#line 1342 "p.y"
    { domain->solInfo().pitaLocalBasisImprovement = 1; ;}
    break;

  case 406:

/* Line 1455 of yacc.c  */
#line 1344 "p.y"
    { domain->solInfo().pitaTimeReversible = true; ;}
    break;

  case 407:

/* Line 1455 of yacc.c  */
#line 1346 "p.y"
    { domain->solInfo().pitaRemoteCoarse = true; ;}
    break;

  case 408:

/* Line 1455 of yacc.c  */
#line 1348 "p.y"
    { domain->solInfo().pitaProjTol = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 409:

/* Line 1455 of yacc.c  */
#line 1350 "p.y"
    { domain->solInfo().pitaReadInitSeed = true; ;}
    break;

  case 410:

/* Line 1455 of yacc.c  */
#line 1352 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = 0.0; ;}
    break;

  case 411:

/* Line 1455 of yacc.c  */
#line 1354 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 412:

/* Line 1455 of yacc.c  */
#line 1356 "p.y"
    { domain->solInfo().pitaJumpMagnOutput = true; ;}
    break;

  case 413:

/* Line 1455 of yacc.c  */
#line 1360 "p.y"
    { domain->solInfo().setDamping((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 414:

/* Line 1455 of yacc.c  */
#line 1362 "p.y"
    { domain->solInfo().setDamping((yyvsp[(2) - (5)].fval),(yyvsp[(3) - (5)].fval));
          domain->solInfo().mtypeDamp = (int)(yyvsp[(4) - (5)].ival); ;}
    break;

  case 415:

/* Line 1455 of yacc.c  */
#line 1365 "p.y"
    { if(geoSource->setModalDamping((yyvsp[(3) - (3)].bclist)->n, (yyvsp[(3) - (3)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 416:

/* Line 1455 of yacc.c  */
#line 1370 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); ;}
    break;

  case 417:

/* Line 1455 of yacc.c  */
#line 1374 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 418:

/* Line 1455 of yacc.c  */
#line 1380 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval), 0.0);
        ;}
    break;

  case 419:

/* Line 1455 of yacc.c  */
#line 1388 "p.y"
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 420:

/* Line 1455 of yacc.c  */
#line 1394 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[(2) - (3)].ival),0.0,0.0,0.0);
        ;}
    break;

  case 422:

/* Line 1455 of yacc.c  */
#line 1402 "p.y"
    {
           domain->setWaveDirections((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        ;}
    break;

  case 423:

/* Line 1455 of yacc.c  */
#line 1408 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 424:

/* Line 1455 of yacc.c  */
#line 1416 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 425:

/* Line 1455 of yacc.c  */
#line 1418 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 426:

/* Line 1455 of yacc.c  */
#line 1420 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } ;}
    break;

  case 427:

/* Line 1455 of yacc.c  */
#line 1422 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } ;}
    break;

  case 428:

/* Line 1455 of yacc.c  */
#line 1424 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; ;}
    break;

  case 430:

/* Line 1455 of yacc.c  */
#line 1440 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(2) - (6)].ival)-1;
          surf_bc[0].type = (BCond::BCType) (yyvsp[(3) - (6)].ival); //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = (yyvsp[(4) - (6)].ival)-1;
          surf_bc[0].val = (yyvsp[(5) - (6)].ival)-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
          if(geoSource->getNumSurfaceConstraint() > 1) delete [] surf_bc;
        ;}
    break;

  case 431:

/* Line 1455 of yacc.c  */
#line 1450 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 432:

/* Line 1455 of yacc.c  */
#line 1454 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 433:

/* Line 1455 of yacc.c  */
#line 1456 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 434:

/* Line 1455 of yacc.c  */
#line 1460 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 435:

/* Line 1455 of yacc.c  */
#line 1462 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (2)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; ;}
    break;

  case 436:

/* Line 1455 of yacc.c  */
#line 1466 "p.y"
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; 

          int* allPDirNodes = new int[(yyvsp[(3) - (3)].bclist)->n];

          for (int ii=0; ii < (yyvsp[(3) - (3)].bclist)->n; ii++)
            allPDirNodes[ii]=((yyvsp[(3) - (3)].bclist)->d[ii]).nnum;
          std::sort(allPDirNodes,allPDirNodes + (yyvsp[(3) - (3)].bclist)->n);

          int maxFSNodes = 32;
          int* allHEVFSNodes = new int[maxFSNodes];
          allHEVFSNodes[0] = allPDirNodes[0];

          int numHEVFSNodes = 1;
          for (int ii = 1; ii < (yyvsp[(3) - (3)].bclist)->n; ++ii)  {
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
        ;}
    break;

  case 437:

/* Line 1455 of yacc.c  */
#line 1507 "p.y"
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[(1) - (1)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(1) - (1)].bclist)->d)[ii]); ;}
    break;

  case 438:

/* Line 1455 of yacc.c  */
#line 1511 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); 
          for (int ii = 0; ii < (yyvsp[(2) - (2)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(2) - (2)].bclist)->d)[ii]); ;}
    break;

  case 439:

/* Line 1455 of yacc.c  */
#line 1517 "p.y"
    { (yyval.bclist) = new BCList;
          for(int i=0; i<(yyvsp[(3) - (4)].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[(3) - (4)].nl).nd[i],10,0.0); } ;}
    break;

  case 440:

/* Line 1455 of yacc.c  */
#line 1523 "p.y"
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;;}
    break;

  case 441:

/* Line 1455 of yacc.c  */
#line 1525 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); ;}
    break;

  case 442:

/* Line 1455 of yacc.c  */
#line 1528 "p.y"
    { for(int i=(yyvsp[(2) - (6)].ival); i<=(yyvsp[(4) - (6)].ival); ++i) { BCond bc; bc.setData(i-1, 6, (yyvsp[(5) - (6)].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } ;}
    break;

  case 443:

/* Line 1455 of yacc.c  */
#line 1530 "p.y"
    { for(int i=(yyvsp[(2) - (8)].ival); i<=(yyvsp[(4) - (8)].ival); i+=(yyvsp[(6) - (8)].ival)) { BCond bc; bc.setData(i-1, 6, (yyvsp[(7) - (8)].fval), BCond::Temperatures); (yyval.bclist)->add(bc); } ;}
    break;

  case 444:

/* Line 1455 of yacc.c  */
#line 1532 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); ;}
    break;

  case 445:

/* Line 1455 of yacc.c  */
#line 1541 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 446:

/* Line 1455 of yacc.c  */
#line 1543 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 447:

/* Line 1455 of yacc.c  */
#line 1545 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Flux; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); ;}
    break;

  case 448:

/* Line 1455 of yacc.c  */
#line 1548 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; ;}
    break;

  case 449:

/* Line 1455 of yacc.c  */
#line 1559 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 450:

/* Line 1455 of yacc.c  */
#line 1561 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 451:

/* Line 1455 of yacc.c  */
#line 1563 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Convection; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); ;}
    break;

  case 452:

/* Line 1455 of yacc.c  */
#line 1568 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 453:

/* Line 1455 of yacc.c  */
#line 1570 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 454:

/* Line 1455 of yacc.c  */
#line 1572 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Radiation; (yyval.bclist)->add(bc); ;}
    break;

  case 458:

/* Line 1455 of yacc.c  */
#line 1584 "p.y"
    { domain->addSommer(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); ;}
    break;

  case 459:

/* Line 1455 of yacc.c  */
#line 1586 "p.y"
    { domain->addSommer(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); ;}
    break;

  case 460:

/* Line 1455 of yacc.c  */
#line 1588 "p.y"
    { domain->addSommer(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); ;}
    break;

  case 463:

/* Line 1455 of yacc.c  */
#line 1596 "p.y"
    { domain->addSommerElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        ;}
    break;

  case 464:

/* Line 1455 of yacc.c  */
#line 1602 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1; ;}
    break;

  case 465:

/* Line 1455 of yacc.c  */
#line 1604 "p.y"
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++; ;}
    break;

  case 469:

/* Line 1455 of yacc.c  */
#line 1616 "p.y"
    { domain->addScatter(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); ;}
    break;

  case 470:

/* Line 1455 of yacc.c  */
#line 1619 "p.y"
    { domain->addScatter(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); ;}
    break;

  case 471:

/* Line 1455 of yacc.c  */
#line 1622 "p.y"
    { domain->addScatter(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); ;}
    break;

  case 474:

/* Line 1455 of yacc.c  */
#line 1631 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 477:

/* Line 1455 of yacc.c  */
#line 1640 "p.y"
    { domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 480:

/* Line 1455 of yacc.c  */
#line 1648 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; ;}
    break;

  case 483:

/* Line 1455 of yacc.c  */
#line 1658 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);;}
    break;

  case 484:

/* Line 1455 of yacc.c  */
#line 1662 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 485:

/* Line 1455 of yacc.c  */
#line 1666 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 486:

/* Line 1455 of yacc.c  */
#line 1668 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 487:

/* Line 1455 of yacc.c  */
#line 1672 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 488:

/* Line 1455 of yacc.c  */
#line 1676 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 489:

/* Line 1455 of yacc.c  */
#line 1678 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) { (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[(4) - (4)].bclist); ;}
    break;

  case 490:

/* Line 1455 of yacc.c  */
#line 1682 "p.y"
    { domain->solInfo().ATDARBFlag = (yyvsp[(2) - (3)].fval);;}
    break;

  case 492:

/* Line 1455 of yacc.c  */
#line 1687 "p.y"
    { domain->solInfo().ATDDNBVal = (yyvsp[(2) - (3)].fval);;}
    break;

  case 494:

/* Line 1455 of yacc.c  */
#line 1692 "p.y"
    { domain->solInfo().ATDROBVal = (yyvsp[(2) - (5)].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[(3) - (5)].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[(4) - (5)].fval);;}
    break;

  case 496:

/* Line 1455 of yacc.c  */
#line 1699 "p.y"
    { domain->setFFP((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 497:

/* Line 1455 of yacc.c  */
#line 1701 "p.y"
    { domain->setFFP((yyvsp[(2) - (4)].ival),(yyvsp[(3) - (4)].ival)); ;}
    break;

  case 498:

/* Line 1455 of yacc.c  */
#line 1705 "p.y"
    {
           domain->setFFP((yyvsp[(2) - (4)].ival));
        ;}
    break;

  case 499:

/* Line 1455 of yacc.c  */
#line 1711 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        ;}
    break;

  case 500:

/* Line 1455 of yacc.c  */
#line 1716 "p.y"
    { fourHelmBC->addDirichlet((yyvsp[(2) - (2)].complexFDBC)); ;}
    break;

  case 501:

/* Line 1455 of yacc.c  */
#line 1720 "p.y"
    { (yyval.complexFDBC) = FDBC((yyvsp[(1) - (2)].ival)-1); ;}
    break;

  case 502:

/* Line 1455 of yacc.c  */
#line 1724 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        ;}
    break;

  case 503:

/* Line 1455 of yacc.c  */
#line 1729 "p.y"
    { fourHelmBC->addNeuman((yyvsp[(2) - (2)].complexFNBC)); ;}
    break;

  case 504:

/* Line 1455 of yacc.c  */
#line 1733 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1); ;}
    break;

  case 505:

/* Line 1455 of yacc.c  */
#line 1735 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 506:

/* Line 1455 of yacc.c  */
#line 1739 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setModes((yyvsp[(2) - (3)].ival));
          domain->solInfo().setProbType(SolverInfo::AxiHelm);
        ;}
    break;

  case 507:

/* Line 1455 of yacc.c  */
#line 1747 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSlices((yyvsp[(2) - (3)].ival));
        ;}
    break;

  case 508:

/* Line 1455 of yacc.c  */
#line 1754 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSomType((yyvsp[(3) - (7)].ival));
          fourHelmBC->setSurf((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval));
        ;}
    break;

  case 510:

/* Line 1455 of yacc.c  */
#line 1762 "p.y"
    { fourHelmBC->addSommer(new LineAxiSommer((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1)); ;}
    break;

  case 511:

/* Line 1455 of yacc.c  */
#line 1764 "p.y"
    { fourHelmBC->addSommer(new Line2AxiSommer((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1)); ;}
    break;

  case 512:

/* Line 1455 of yacc.c  */
#line 1768 "p.y"
    { if( globalMPCs== NULL) globalMPCs = new MPCData(); ;}
    break;

  case 513:

/* Line 1455 of yacc.c  */
#line 1770 "p.y"
    { globalMPCs->addMPC((yyvsp[(2) - (2)].axiMPC)); ;}
    break;

  case 514:

/* Line 1455 of yacc.c  */
#line 1774 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (5)].ival)-1, (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival), DComplex(1.0,0.0), 0.0, 0.0, 0.0); ;}
    break;

  case 515:

/* Line 1455 of yacc.c  */
#line 1776 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (7)].ival)-1, (yyvsp[(2) - (7)].fval), (yyvsp[(3) - (7)].fval), (yyvsp[(4) - (7)].ival), DComplex((yyvsp[(5) - (7)].fval),(yyvsp[(6) - (7)].fval)) , 0.0, 0.0, 0.0); ;}
    break;

  case 516:

/* Line 1455 of yacc.c  */
#line 1778 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (10)].ival)-1, (yyvsp[(2) - (10)].fval), (yyvsp[(3) - (10)].fval), (yyvsp[(4) - (10)].ival), DComplex((yyvsp[(5) - (10)].fval),(yyvsp[(6) - (10)].fval)) , (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)); ;}
    break;

  case 517:

/* Line 1455 of yacc.c  */
#line 1782 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (3)].strval);
	  domain->solInfo().readmodeCalled = true; ;}
    break;

  case 518:

/* Line 1455 of yacc.c  */
#line 1785 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (4)].strval);
          domain->solInfo().readmodeCalled = true; 
 	  domain->solInfo().maxSizePodRom = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 519:

/* Line 1455 of yacc.c  */
#line 1789 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (4)].strval);
          domain->solInfo().readInModes = (yyvsp[(3) - (4)].strval);
          domain->solInfo().readmodeCalled = true; ;}
    break;

  case 520:

/* Line 1455 of yacc.c  */
#line 1793 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (5)].strval);
          domain->solInfo().readInModes = (yyvsp[(3) - (5)].strval);
          domain->solInfo().readmodeCalled = true;
          domain->solInfo().maxSizePodRom = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 521:

/* Line 1455 of yacc.c  */
#line 1798 "p.y"
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 522:

/* Line 1455 of yacc.c  */
#line 1802 "p.y"
    { ;}
    break;

  case 523:

/* Line 1455 of yacc.c  */
#line 1804 "p.y"
    { domain->solInfo().zeroInitialDisp = 1; ;}
    break;

  case 524:

/* Line 1455 of yacc.c  */
#line 1806 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 525:

/* Line 1455 of yacc.c  */
#line 1809 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 526:

/* Line 1455 of yacc.c  */
#line 1815 "p.y"
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[(2) - (3)].fval);  ;}
    break;

  case 527:

/* Line 1455 of yacc.c  */
#line 1817 "p.y"
    { (yyval.bclist) = new BCList; amplitude = 1.0; ;}
    break;

  case 528:

/* Line 1455 of yacc.c  */
#line 1819 "p.y"
    { BCond bc; /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (9)].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[(3) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[(4) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[(5) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = amplitude*(yyvsp[(6) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = amplitude*(yyvsp[(7) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = amplitude*(yyvsp[(8) - (9)].fval); (yyval.bclist)->add(bc);
          for(int i=0; i<(yyval.bclist)->n; ++i) (yyval.bclist)->d[i].type = BCond::Idisp6;
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        ;}
    break;

  case 529:

/* Line 1455 of yacc.c  */
#line 1830 "p.y"
    { BCond bc; /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[(3) - (6)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[(4) - (6)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[(5) - (6)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = 0.0         ; (yyval.bclist)->add(bc);
          for(int i=0; i<(yyval.bclist)->n; ++i) (yyval.bclist)->d[i].type = BCond::Idisp6;
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        ;}
    break;

  case 530:

/* Line 1455 of yacc.c  */
#line 1841 "p.y"
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); ;}
    break;

  case 531:

/* Line 1455 of yacc.c  */
#line 1844 "p.y"
    { domain->solInfo().buckling = 1; ;}
    break;

  case 532:

/* Line 1455 of yacc.c  */
#line 1849 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 533:

/* Line 1455 of yacc.c  */
#line 1851 "p.y"
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (9)].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[(3) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[(4) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[(5) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[(6) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[(7) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[(8) - (9)].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIDis6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        ;}
    break;

  case 534:

/* Line 1455 of yacc.c  */
#line 1864 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 535:

/* Line 1455 of yacc.c  */
#line 1866 "p.y"
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (9)].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[(3) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[(4) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[(5) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[(6) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[(7) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[(8) - (9)].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIVel6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        ;}
    break;

  case 536:

/* Line 1455 of yacc.c  */
#line 1878 "p.y"
    { ;}
    break;

  case 537:

/* Line 1455 of yacc.c  */
#line 1880 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 538:

/* Line 1455 of yacc.c  */
#line 1883 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 539:

/* Line 1455 of yacc.c  */
#line 1889 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 540:

/* Line 1455 of yacc.c  */
#line 1894 "p.y"
    { domain->solInfo().setGEPS();
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Etemperatures;
          if(geoSource->setIDis6((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 541:

/* Line 1455 of yacc.c  */
#line 1900 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 542:

/* Line 1455 of yacc.c  */
#line 1902 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 543:

/* Line 1455 of yacc.c  */
#line 1904 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Forces; (yyvsp[(2) - (2)].bcval).loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 544:

/* Line 1455 of yacc.c  */
#line 1906 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } ;}
    break;

  case 545:

/* Line 1455 of yacc.c  */
#line 1908 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } ;}
    break;

  case 546:

/* Line 1455 of yacc.c  */
#line 1910 "p.y"
    { for(int i=(yyvsp[(2) - (8)].ival); i<=(yyvsp[(4) - (8)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (8)].ival)-1, (yyvsp[(6) - (8)].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[(7) - (8)].ival)); (yyval.bclist)->add(bc); } ;}
    break;

  case 547:

/* Line 1455 of yacc.c  */
#line 1912 "p.y"
    { for(int i=(yyvsp[(2) - (10)].ival); i<=(yyvsp[(4) - (10)].ival); i+=(yyvsp[(6) - (10)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (10)].ival)-1, (yyvsp[(8) - (10)].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[(7) - (10)].ival)); (yyval.bclist)->add(bc); } ;}
    break;

  case 548:

/* Line 1455 of yacc.c  */
#line 1914 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; ;}
    break;

  case 549:

/* Line 1455 of yacc.c  */
#line 1923 "p.y"
    { for(int i=0; i<(yyvsp[(5) - (5)].bclist)->n; ++i) (yyvsp[(5) - (5)].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[(5) - (5)].bclist)->n, (yyvsp[(5) - (5)].bclist)->d) < 0) return -1; ;}
    break;

  case 550:

/* Line 1455 of yacc.c  */
#line 1926 "p.y"
    { for(int i=0; i<(yyvsp[(6) - (6)].bclist)->n; ++i) { (yyvsp[(6) - (6)].bclist)->d[i].type = BCond::Forces; (yyvsp[(6) - (6)].bclist)->d[i].loadsetid = (yyvsp[(2) - (6)].ival); }
          if(geoSource->setNeumanModal((yyvsp[(6) - (6)].bclist)->n, (yyvsp[(6) - (6)].bclist)->d) < 0) return -1; ;}
    break;

  case 551:

/* Line 1455 of yacc.c  */
#line 1931 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 552:

/* Line 1455 of yacc.c  */
#line 1933 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 553:

/* Line 1455 of yacc.c  */
#line 1935 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (6)].ival); i<=(yyvsp[(3) - (6)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 554:

/* Line 1455 of yacc.c  */
#line 1937 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (7)].bclist); for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 555:

/* Line 1455 of yacc.c  */
#line 1939 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (8)].ival); i<=(yyvsp[(3) - (8)].ival); i+=(yyvsp[(5) - (8)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(6) - (8)].ival)-1, (yyvsp[(7) - (8)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 556:

/* Line 1455 of yacc.c  */
#line 1941 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (9)].bclist); for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 557:

/* Line 1455 of yacc.c  */
#line 1945 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 558:

/* Line 1455 of yacc.c  */
#line 1947 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 559:

/* Line 1455 of yacc.c  */
#line 1951 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 560:

/* Line 1455 of yacc.c  */
#line 1953 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 563:

/* Line 1455 of yacc.c  */
#line 1961 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ymtt)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addYMTT((yyval.ymtt));;}
    break;

  case 564:

/* Line 1455 of yacc.c  */
#line 1963 "p.y"
    { (yyval.ymtt)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 565:

/* Line 1455 of yacc.c  */
#line 1965 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ymtt)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addYMTT((yyval.ymtt));;}
    break;

  case 568:

/* Line 1455 of yacc.c  */
#line 1973 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ctett)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addCTETT((yyval.ctett));;}
    break;

  case 569:

/* Line 1455 of yacc.c  */
#line 1975 "p.y"
    { (yyval.ctett)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 570:

/* Line 1455 of yacc.c  */
#line 1977 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ctett)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addCTETT((yyval.ctett));;}
    break;

  case 573:

/* Line 1455 of yacc.c  */
#line 1985 "p.y"
    { (yyval.sdetaft) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.sdetaft)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addSDETAFT((yyval.sdetaft));;}
    break;

  case 574:

/* Line 1455 of yacc.c  */
#line 1987 "p.y"
    { (yyval.sdetaft)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 575:

/* Line 1455 of yacc.c  */
#line 1989 "p.y"
    { (yyval.sdetaft) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.sdetaft)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addSDETAFT((yyval.sdetaft));;}
    break;

  case 578:

/* Line 1455 of yacc.c  */
#line 1997 "p.y"
    { (yyval.lmpcons) = (yyvsp[(1) - (2)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 579:

/* Line 1455 of yacc.c  */
#line 2001 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); ;}
    break;

  case 580:

/* Line 1455 of yacc.c  */
#line 2003 "p.y"
    { (yyval.lmpcons) = (yyvsp[(2) - (3)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(3) - (3)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 581:

/* Line 1455 of yacc.c  */
#line 2009 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 582:

/* Line 1455 of yacc.c  */
#line 2012 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (3)].ival), (yyvsp[(2) - (3)].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 583:

/* Line 1455 of yacc.c  */
#line 2015 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval));
          (yyval.lmpcons)->type = (yyvsp[(4) - (5)].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 584:

/* Line 1455 of yacc.c  */
#line 2019 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (4)].ival), (yyvsp[(2) - (4)].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[(3) - (4)].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 585:

/* Line 1455 of yacc.c  */
#line 2024 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (6)].ival), (yyvsp[(2) - (6)].fval));
          (yyval.lmpcons)->type = (yyvsp[(4) - (6)].ival);
          (yyval.lmpcons)->lagrangeMult = (yyvsp[(5) - (6)].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[(5) - (6)].copt).penalty;
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 586:

/* Line 1455 of yacc.c  */
#line 2032 "p.y"
    { if((yyvsp[(3) - (4)].fval) == 0.0) {
            fprintf(stderr," *** WARNING: zero coefficient in LMPC\n");
            fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival));
          }
          (yyval.mpcterm) = new LMPCTerm();
          (yyval.mpcterm)->nnum = (yyvsp[(1) - (4)].ival)-1;
          (yyval.mpcterm)->dofnum = (yyvsp[(2) - (4)].ival)-1;
          (yyval.mpcterm)->coef.r_value = (yyvsp[(3) - (4)].fval);
        ;}
    break;

  case 589:

/* Line 1455 of yacc.c  */
#line 2048 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].cxbcval).nnum,(yyvsp[(1) - (2)].cxbcval).reval,(yyvsp[(1) - (2)].cxbcval).imval,(yyvsp[(2) - (2)].mpcterm)); domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 590:

/* Line 1455 of yacc.c  */
#line 2050 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); ;}
    break;

  case 591:

/* Line 1455 of yacc.c  */
#line 2052 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(2) - (3)].cxbcval).nnum,(yyvsp[(2) - (3)].cxbcval).reval,(yyvsp[(2) - (3)].cxbcval).imval,(yyvsp[(3) - (3)].mpcterm)); domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 592:

/* Line 1455 of yacc.c  */
#line 2056 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (5)].ival); (yyval.cxbcval).reval=(yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval=(yyvsp[(4) - (5)].fval); ;}
    break;

  case 593:

/* Line 1455 of yacc.c  */
#line 2058 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (3)].ival); (yyval.cxbcval).reval=(yyvsp[(2) - (3)].fval); (yyval.cxbcval).imval=0.0; ;}
    break;

  case 594:

/* Line 1455 of yacc.c  */
#line 2060 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (2)].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; ;}
    break;

  case 595:

/* Line 1455 of yacc.c  */
#line 2064 "p.y"
    { if(((yyvsp[(3) - (5)].fval)==0.0) && ((yyvsp[(4) - (5)].fval)==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (5)].ival),(yyvsp[(2) - (5)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (5)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (5)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval)); }
        ;}
    break;

  case 596:

/* Line 1455 of yacc.c  */
#line 2072 "p.y"
    { if((yyvsp[(3) - (4)].fval)==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (4)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (4)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (4)].fval),0.0); }
        ;}
    break;

  case 597:

/* Line 1455 of yacc.c  */
#line 2082 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); ;}
    break;

  case 598:

/* Line 1455 of yacc.c  */
#line 2084 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].cxbclist)->n; ++i) (yyvsp[(4) - (4)].cxbclist)->d[i].loadsetid = (yyvsp[(2) - (4)].ival);
          (yyval.cxbclist) = (yyvsp[(4) - (4)].cxbclist); ;}
    break;

  case 599:

/* Line 1455 of yacc.c  */
#line 2089 "p.y"
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[(1) - (1)].cxbcval)); ;}
    break;

  case 600:

/* Line 1455 of yacc.c  */
#line 2091 "p.y"
    { (yyval.cxbclist) = (yyvsp[(1) - (2)].cxbclist); (yyval.cxbclist)->add((yyvsp[(2) - (2)].cxbcval)); ;}
    break;

  case 603:

/* Line 1455 of yacc.c  */
#line 2099 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval); 
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.Iyy = (yyvsp[(14) - (16)].fval); sp.Izz = (yyvsp[(15) - (16)].fval);
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 604:

/* Line 1455 of yacc.c  */
#line 2106 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.betaDamp = (yyvsp[(17) - (19)].fval); sp.alphaDamp = (yyvsp[(18) - (19)].fval);
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        ;}
    break;

  case 605:

/* Line 1455 of yacc.c  */
#line 2114 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.etaDamp = (yyvsp[(17) - (19)].fval); sp.betaDamp = (yyvsp[(18) - (19)].fval);
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        ;}
    break;

  case 606:

/* Line 1455 of yacc.c  */
#line 2122 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (20)].fval);  sp.E = (yyvsp[(3) - (20)].fval);  sp.nu  = (yyvsp[(4) - (20)].fval);  sp.rho = (yyvsp[(5) - (20)].fval);
          sp.c = (yyvsp[(6) - (20)].fval);  sp.k = (yyvsp[(7) - (20)].fval);  sp.eh  = (yyvsp[(8) - (20)].fval);  sp.P   = (yyvsp[(9) - (20)].fval);  sp.Ta  = (yyvsp[(10) - (20)].fval); 
          sp.Q = (yyvsp[(11) - (20)].fval); sp.W = (yyvsp[(12) - (20)].fval); sp.Ixx = (yyvsp[(13) - (20)].fval); sp.Iyy = (yyvsp[(14) - (20)].fval); sp.Izz = (yyvsp[(15) - (20)].fval);
	  sp.ymin = (yyvsp[(16) - (20)].fval); sp.ymax = (yyvsp[(17) - (20)].fval); sp.zmin = (yyvsp[(18) - (20)].fval); sp.zmax = (yyvsp[(19) - (20)].fval);
          geoSource->addMat( (yyvsp[(1) - (20)].ival)-1, sp );
        ;}
    break;

  case 607:

/* Line 1455 of yacc.c  */
#line 2130 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.betaDamp = (yyvsp[(21) - (23)].fval); sp.alphaDamp = (yyvsp[(22) - (23)].fval);
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        ;}
    break;

  case 608:

/* Line 1455 of yacc.c  */
#line 2139 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.etaDamp = (yyvsp[(21) - (23)].fval); sp.betaDamp = (yyvsp[(22) - (23)].fval);
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        ;}
    break;

  case 609:

/* Line 1455 of yacc.c  */
#line 2148 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (9)].fval); sp.E = (yyvsp[(3) - (9)].fval); sp.nu = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(5) - (9)].fval);
          sp.c = (yyvsp[(6) - (9)].fval); sp.k = (yyvsp[(7) - (9)].fval); sp.eh = (yyvsp[(8) - (9)].fval);
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp ); 
        ;}
    break;

  case 610:

/* Line 1455 of yacc.c  */
#line 2154 "p.y"
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (14)].fval);  sp.E = (yyvsp[(3) - (14)].fval);  sp.nu  = (yyvsp[(4) - (14)].fval);  sp.rho = (yyvsp[(5) - (14)].fval);
          sp.c = (yyvsp[(6) - (14)].fval);  sp.k = (yyvsp[(7) - (14)].fval);  sp.eh  = (yyvsp[(8) - (14)].fval);  sp.P   = (yyvsp[(9) - (14)].fval);  sp.Ta  = (yyvsp[(10) - (14)].fval);
          sp.Q = (yyvsp[(11) - (14)].fval); sp.W = (yyvsp[(12) - (14)].fval); sp.Ixx = (yyvsp[(13) - (14)].fval);  
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
        ;}
    break;

  case 611:

/* Line 1455 of yacc.c  */
#line 2161 "p.y"
    { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval);
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.betaDamp = (yyvsp[(15) - (16)].fval);
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 612:

/* Line 1455 of yacc.c  */
#line 2169 "p.y"
    { StructProp sp; // this is used for the reduced mesh file output in Rom.d/MeshOutput.C
                         // all properties relevant to structural nonlinear dynamics should be included
          sp.A = (yyvsp[(2) - (30)].fval);  sp.E = (yyvsp[(3) - (30)].fval);  sp.nu  = (yyvsp[(4) - (30)].fval);  sp.rho = (yyvsp[(5) - (30)].fval);
          sp.c = (yyvsp[(6) - (30)].fval);  sp.k = (yyvsp[(7) - (30)].fval);  sp.eh  = (yyvsp[(8) - (30)].fval);  sp.P   = (yyvsp[(9) - (30)].fval);  sp.Ta  = (yyvsp[(10) - (30)].fval);
          sp.Q = (yyvsp[(11) - (30)].fval); sp.W = (yyvsp[(12) - (30)].fval); sp.Ixx = (yyvsp[(13) - (30)].fval); sp.Iyy = (yyvsp[(14) - (30)].fval); sp.Izz = (yyvsp[(15) - (30)].fval);
          sp.ymin = (yyvsp[(16) - (30)].fval); sp.ymax = (yyvsp[(17) - (30)].fval); sp.zmin = (yyvsp[(18) - (30)].fval); sp.zmax = (yyvsp[(19) - (30)].fval);
          sp.betaDamp = (yyvsp[(20) - (30)].fval); sp.alphaDamp = (yyvsp[(21) - (30)].fval); 
          sp.lagrangeMult = bool((yyvsp[(22) - (30)].ival)); sp.penalty = (yyvsp[(23) - (30)].fval); sp.initialPenalty = (yyvsp[(24) - (30)].fval);
          sp.funtype = (yyvsp[(25) - (30)].ival); sp.type = StructProp::PropType((yyvsp[(26) - (30)].ival)); sp.k1 = (yyvsp[(27) - (30)].fval); sp.k2 = (yyvsp[(28) - (30)].fval); sp.k3 = (yyvsp[(29) - (30)].fval);
          geoSource->addMat( (yyvsp[(1) - (30)].ival)-1, sp );
        ;}
    break;

  case 613:

/* Line 1455 of yacc.c  */
#line 2181 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (12)].fval),0.0);
          sp.fp.PMLtype = int((yyvsp[(4) - (12)].fval));
          sp.fp.gamma = (yyvsp[(5) - (12)].fval);
          sp.fp.Rx = (yyvsp[(6) - (12)].fval);
          sp.fp.Sx = (yyvsp[(7) - (12)].fval);
          sp.fp.Ry = (yyvsp[(8) - (12)].fval);
          sp.fp.Sy = (yyvsp[(9) - (12)].fval);
          sp.fp.Rz = (yyvsp[(10) - (12)].fval);
          sp.fp.Sz = (yyvsp[(11) - (12)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 614:

/* Line 1455 of yacc.c  */
#line 2197 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (13)].fval),0.0);
          sp.rho = (yyvsp[(4) - (13)].fval);
          sp.fp.PMLtype = int((yyvsp[(5) - (13)].fval));
          sp.fp.gamma = (yyvsp[(6) - (13)].fval);
          sp.fp.Rx = (yyvsp[(7) - (13)].fval);
          sp.fp.Sx = (yyvsp[(8) - (13)].fval);
          sp.fp.Ry = (yyvsp[(9) - (13)].fval);
          sp.fp.Sy = (yyvsp[(10) - (13)].fval);
          sp.fp.Rz = (yyvsp[(11) - (13)].fval);
          sp.fp.Sz = (yyvsp[(12) - (13)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (13)].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 615:

/* Line 1455 of yacc.c  */
#line 2214 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (14)].fval),(yyvsp[(4) - (14)].fval));
          sp.rho = (yyvsp[(5) - (14)].fval);
          sp.fp.PMLtype = int((yyvsp[(6) - (14)].fval));
          sp.fp.gamma = (yyvsp[(7) - (14)].fval);
          sp.fp.Rx = (yyvsp[(8) - (14)].fval);
          sp.fp.Sx = (yyvsp[(9) - (14)].fval);
          sp.fp.Ry = (yyvsp[(10) - (14)].fval);
          sp.fp.Sy = (yyvsp[(11) - (14)].fval);
          sp.fp.Rz = (yyvsp[(12) - (14)].fval);
          sp.fp.Sz = (yyvsp[(13) - (14)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 616:

/* Line 1455 of yacc.c  */
#line 2231 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (5)].fval),0.0);
          sp.rho = (yyvsp[(4) - (5)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 617:

/* Line 1455 of yacc.c  */
#line 2239 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (6)].fval),(yyvsp[(4) - (6)].fval));
          sp.rho = (yyvsp[(5) - (6)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 618:

/* Line 1455 of yacc.c  */
#line 2247 "p.y"
    { StructProp sp;
          sp.F_op = (yyvsp[(3) - (16)].ival);
          sp.E = (yyvsp[(4) - (16)].fval);
          sp.rho = (yyvsp[(5) - (16)].fval);
	  sp.A = (yyvsp[(6) - (16)].fval);
          sp.F_Uc = (yyvsp[(7) - (16)].fval);
          sp.F_Uf = (yyvsp[(8) - (16)].fval);
          sp.lambda = (yyvsp[(9) - (16)].fval);
          sp.F_h = (yyvsp[(10) - (16)].fval);
          sp.F_d = (yyvsp[(11) - (16)].fval);
          sp.F_dlambda = (yyvsp[(12) - (16)].fval);
          sp.F_np = (yyvsp[(13) - (16)].ival);
          sp.F_Nf = (yyvsp[(14) - (16)].ival);
	  sp.Seed = (yyvsp[(15) - (16)].ival);
          sp.type = StructProp::Fabric;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 619:

/* Line 1455 of yacc.c  */
#line 2265 "p.y"
    { StructProp sp; 
          sp.A = (yyvsp[(3) - (12)].fval); sp.rho = (yyvsp[(4) - (12)].fval); sp.Q = (yyvsp[(5) - (12)].fval); sp.c = (yyvsp[(6) - (12)].fval); 
          sp.sigma = (yyvsp[(7) - (12)].fval); sp.k = (yyvsp[(8) - (12)].fval); sp.eh = (yyvsp[(9) - (12)].fval); sp.P = (yyvsp[(10) - (12)].fval); sp.Ta = (yyvsp[(11) - (12)].fval);
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 620:

/* Line 1455 of yacc.c  */
#line 2272 "p.y"
    { // rigid element or joint with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (3)].ival)-1, sp );
        ;}
    break;

  case 621:

/* Line 1455 of yacc.c  */
#line 2279 "p.y"
    { // rigid element or joint
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (4)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (4)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (4)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (4)].ival)-1, sp );
        ;}
    break;

  case 622:

/* Line 1455 of yacc.c  */
#line 2290 "p.y"
    { //           m     Ixx   Iyy   Izz   Ixy   Iyz   Ixz   cx    cy    cz
          // discrete mass with offset
          StructProp sp;
          sp.rho = (yyvsp[(3) - (13)].fval);
          sp.Ixx = (yyvsp[(4) - (13)].fval);
          sp.Iyy = (yyvsp[(5) - (13)].fval);
          sp.Izz = (yyvsp[(6) - (13)].fval);
          sp.Ixy = (yyvsp[(7) - (13)].fval);
          sp.Iyz = (yyvsp[(8) - (13)].fval);
          sp.Ixz = (yyvsp[(9) - (13)].fval);
          sp.cx  = (yyvsp[(10) - (13)].fval);
          sp.cy  = (yyvsp[(11) - (13)].fval);
          sp.cz  = (yyvsp[(12) - (13)].fval);
          geoSource->addMat( (yyvsp[(1) - (13)].ival)-1, sp );
        ;}
    break;

  case 623:

/* Line 1455 of yacc.c  */
#line 2306 "p.y"
    { // rigid 8-node brick element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[(4) - (5)].fval);
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        ;}
    break;

  case 624:

/* Line 1455 of yacc.c  */
#line 2313 "p.y"
    { // rigid 8-node brick element with mass
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (6)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (6)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (6)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (6)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = (yyvsp[(5) - (6)].fval);
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 625:

/* Line 1455 of yacc.c  */
#line 2324 "p.y"
    { // rigid beam or shell element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[(4) - (6)].fval);
          sp.A = sp.eh = (yyvsp[(5) - (6)].fval);
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 626:

/* Line 1455 of yacc.c  */
#line 2332 "p.y"
    { // rigid beam or shell element with mass
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (7)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (7)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (7)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (7)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.rho = (yyvsp[(5) - (7)].fval);
          sp.A = sp.eh = (yyvsp[(6) - (7)].fval);
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 627:

/* Line 1455 of yacc.c  */
#line 2344 "p.y"
    { // constraint function element XXX
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (11)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (11)].fval);
          sp.amplitude = (yyvsp[(5) - (11)].fval);
          sp.omega = (yyvsp[(6) - (11)].fval);
          sp.phase = (yyvsp[(7) - (11)].fval);
          sp.B = (yyvsp[(8) - (11)].fval);
          sp.C = (yyvsp[(9) - (11)].fval);
          sp.relop = (yyvsp[(10) - (11)].ival);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 628:

/* Line 1455 of yacc.c  */
#line 2359 "p.y"
    { // joint-with-driver, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (8)].ival);
          sp.amplitude = (yyvsp[(4) - (8)].fval);
          sp.offset = (yyvsp[(5) - (8)].fval);
          sp.c1 = (yyvsp[(6) - (8)].fval);
          sp.c2 = (yyvsp[(7) - (8)].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        ;}
    break;

  case 629:

/* Line 1455 of yacc.c  */
#line 2371 "p.y"
    { // joint-with-driver, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (9)].ival);
          sp.amplitude = (yyvsp[(4) - (9)].fval);
          sp.offset = (yyvsp[(5) - (9)].fval);
          sp.c1 = (yyvsp[(6) - (9)].fval);
          sp.c2 = (yyvsp[(7) - (9)].fval);
          sp.c3 = (yyvsp[(8) - (9)].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        ;}
    break;

  case 630:

/* Line 1455 of yacc.c  */
#line 2384 "p.y"
    { // joint-with-driver, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (10)].ival);
          sp.amplitude = (yyvsp[(4) - (10)].fval);
          sp.offset = (yyvsp[(5) - (10)].fval);
          sp.c1 = (yyvsp[(6) - (10)].fval);
          sp.c2 = (yyvsp[(7) - (10)].fval);
          sp.c3 = (yyvsp[(8) - (10)].fval);
          sp.c4 = (yyvsp[(9) - (10)].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 631:

/* Line 1455 of yacc.c  */
#line 2398 "p.y"
    { // joint-with-driver, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (9)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (9)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (9)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (9)].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[(4) - (9)].ival);
          sp.amplitude = (yyvsp[(5) - (9)].fval);
          sp.offset = (yyvsp[(6) - (9)].fval);
          sp.c1 = (yyvsp[(7) - (9)].fval);
          sp.c2 = (yyvsp[(8) - (9)].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        ;}
    break;

  case 632:

/* Line 1455 of yacc.c  */
#line 2414 "p.y"
    { // joint-with-driver, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (10)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (10)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (10)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (10)].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[(4) - (10)].ival);
          sp.amplitude = (yyvsp[(5) - (10)].fval);
          sp.offset = (yyvsp[(6) - (10)].fval);
          sp.c1 = (yyvsp[(7) - (10)].fval);
          sp.c2 = (yyvsp[(8) - (10)].fval);
          sp.c3 = (yyvsp[(9) - (10)].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 633:

/* Line 1455 of yacc.c  */
#line 2431 "p.y"
    { // joint-with-driver, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (11)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (11)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (11)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (11)].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[(4) - (11)].ival);
          sp.amplitude = (yyvsp[(5) - (11)].fval);
          sp.offset = (yyvsp[(6) - (11)].fval);
          sp.c1 = (yyvsp[(7) - (11)].fval);
          sp.c2 = (yyvsp[(8) - (11)].fval);
          sp.c3 = (yyvsp[(9) - (11)].fval);
          sp.c4 = (yyvsp[(10) - (11)].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 634:

/* Line 1455 of yacc.c  */
#line 2449 "p.y"
    { // actuated joint, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (10)].ival);
          sp.amplitude = (yyvsp[(4) - (10)].fval);
          sp.offset = (yyvsp[(5) - (10)].fval);
          sp.c1 = (yyvsp[(6) - (10)].fval);
          sp.c2 = (yyvsp[(7) - (10)].fval);
          sp.k1 = (yyvsp[(9) - (10)].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 635:

/* Line 1455 of yacc.c  */
#line 2462 "p.y"
    { // actuated joint, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (11)].ival);
          sp.amplitude = (yyvsp[(4) - (11)].fval);
          sp.offset = (yyvsp[(5) - (11)].fval);
          sp.c1 = (yyvsp[(6) - (11)].fval);
          sp.c2 = (yyvsp[(7) - (11)].fval);
          sp.c3 = (yyvsp[(8) - (11)].fval);
          sp.k1 = (yyvsp[(10) - (11)].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 636:

/* Line 1455 of yacc.c  */
#line 2476 "p.y"
    { // actuated joints, 4-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (12)].ival);
          sp.amplitude = (yyvsp[(4) - (12)].fval);
          sp.offset = (yyvsp[(5) - (12)].fval);
          sp.c1 = (yyvsp[(6) - (12)].fval);
          sp.c2 = (yyvsp[(7) - (12)].fval);
          sp.c3 = (yyvsp[(8) - (12)].fval);
          sp.c4 = (yyvsp[(9) - (12)].fval);
          sp.k1 = (yyvsp[(11) - (12)].fval);
          sp.type = StructProp::Undefined;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 637:

/* Line 1455 of yacc.c  */
#line 2491 "p.y"
    { // actuated joint, 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (11)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (11)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (11)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (11)].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[(4) - (11)].ival);
          sp.amplitude = (yyvsp[(5) - (11)].fval);
          sp.offset = (yyvsp[(6) - (11)].fval);
          sp.c1 = (yyvsp[(7) - (11)].fval);
          sp.c2 = (yyvsp[(8) - (11)].fval);
          sp.k1 = (yyvsp[(10) - (11)].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 638:

/* Line 1455 of yacc.c  */
#line 2508 "p.y"
    { // actuated joint, 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (12)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (12)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (12)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (12)].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[(4) - (12)].ival);
          sp.amplitude = (yyvsp[(5) - (12)].fval);
          sp.offset = (yyvsp[(6) - (12)].fval);
          sp.c1 = (yyvsp[(7) - (12)].fval);
          sp.c2 = (yyvsp[(8) - (12)].fval);
          sp.c3 = (yyvsp[(9) - (12)].fval);
          sp.k1 = (yyvsp[(11) - (12)].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 639:

/* Line 1455 of yacc.c  */
#line 2526 "p.y"
    { // actuated joints, 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (13)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (13)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (13)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (13)].copt).constraint_hess_eps;
          sp.funtype = (yyvsp[(4) - (13)].ival);
          sp.amplitude = (yyvsp[(5) - (13)].fval);
          sp.offset = (yyvsp[(6) - (13)].fval);
          sp.c1 = (yyvsp[(7) - (13)].fval);
          sp.c2 = (yyvsp[(8) - (13)].fval);
          sp.c3 = (yyvsp[(9) - (13)].fval);
          sp.c4 = (yyvsp[(10) - (13)].fval);
          sp.k1 = (yyvsp[(12) - (13)].fval);
          sp.type = StructProp::Constraint;
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (13)].ival)-1, sp );
        ;}
    break;

  case 640:

/* Line 1455 of yacc.c  */
#line 2545 "p.y"
    { // RevoluteJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[(4) - (5)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        ;}
    break;

  case 641:

/* Line 1455 of yacc.c  */
#line 2553 "p.y"
    { // RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (6)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (6)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (6)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (6)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(5) - (6)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 642:

/* Line 1455 of yacc.c  */
#line 2565 "p.y"
    { // UniversalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[(4) - (6)].fval);
          sp.k2 = (yyvsp[(5) - (6)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 643:

/* Line 1455 of yacc.c  */
#line 2574 "p.y"
    { // UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (7)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (7)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (7)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (7)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(5) - (7)].fval);
          sp.k2 = (yyvsp[(6) - (7)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 644:

/* Line 1455 of yacc.c  */
#line 2587 "p.y"
    { // SphericalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[(4) - (7)].fval);
          sp.k2 = (yyvsp[(5) - (7)].fval);
          sp.k3 = (yyvsp[(6) - (7)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 645:

/* Line 1455 of yacc.c  */
#line 2597 "p.y"
    { // SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (8)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (8)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (8)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (8)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(5) - (8)].fval);
          sp.k2 = (yyvsp[(6) - (8)].fval);
          sp.k3 = (yyvsp[(7) - (8)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        ;}
    break;

  case 646:

/* Line 1455 of yacc.c  */
#line 2611 "p.y"
    { // TorsionalSpringType1 or TranslationalSpring
          StructProp sp;
          sp.k1 = (yyvsp[(3) - (4)].fval);
          sp.rho = 0;
          geoSource->addMat( (yyvsp[(1) - (4)].ival)-1, sp );
        ;}
    break;

  case 649:

/* Line 1455 of yacc.c  */
#line 2624 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (3)].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 650:

/* Line 1455 of yacc.c  */
#line 2630 "p.y"
    { if((yyvsp[(2) - (4)].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (4)].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 651:

/* Line 1455 of yacc.c  */
#line 2636 "p.y"
    { if((yyvsp[(2) - (5)].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (5)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (5)].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 652:

/* Line 1455 of yacc.c  */
#line 2643 "p.y"
    { if((yyvsp[(2) - (6)].ival) == 0) { std::cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (6)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (6)].fval));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 653:

/* Line 1455 of yacc.c  */
#line 2651 "p.y"
    { if((yyval.SurfObj)->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[(yyvsp[(4) - (5)].nl).num];
            for(int i=0; i<(yyvsp[(4) - (5)].nl).num; ++i) nodes[(yyvsp[(4) - (5)].nl).num-1-i] = (yyvsp[(4) - (5)].nl).nd[i];
            (yyval.SurfObj)->AddFaceElement((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].nl).num, nodes);
            delete [] nodes;
          }
          else (yyval.SurfObj)->AddFaceElement((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].nl).num, (yyvsp[(4) - (5)].nl).nd);
        ;}
    break;

  case 654:

/* Line 1455 of yacc.c  */
#line 2662 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 655:

/* Line 1455 of yacc.c  */
#line 2664 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 656:

/* Line 1455 of yacc.c  */
#line 2666 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 657:

/* Line 1455 of yacc.c  */
#line 2670 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (6)].ival), (yyvsp[(3) - (6)].ival), (yyvsp[(5) - (6)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 658:

/* Line 1455 of yacc.c  */
#line 2672 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 659:

/* Line 1455 of yacc.c  */
#line 2674 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 660:

/* Line 1455 of yacc.c  */
#line 2676 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 661:

/* Line 1455 of yacc.c  */
#line 2678 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 662:

/* Line 1455 of yacc.c  */
#line 2682 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 663:

/* Line 1455 of yacc.c  */
#line 2688 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->solInfo().isCoupled = true; ;}
    break;

  case 664:

/* Line 1455 of yacc.c  */
#line 2690 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (3)].ival), (yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; ;}
    break;

  case 665:

/* Line 1455 of yacc.c  */
#line 2698 "p.y"
    { ;}
    break;

  case 666:

/* Line 1455 of yacc.c  */
#line 2700 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 667:

/* Line 1455 of yacc.c  */
#line 2706 "p.y"
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 668:

/* Line 1455 of yacc.c  */
#line 2713 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 669:

/* Line 1455 of yacc.c  */
#line 2720 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 670:

/* Line 1455 of yacc.c  */
#line 2727 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 671:

/* Line 1455 of yacc.c  */
#line 2735 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 672:

/* Line 1455 of yacc.c  */
#line 2742 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(6) - (7)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 673:

/* Line 1455 of yacc.c  */
#line 2750 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(7) - (8)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 674:

/* Line 1455 of yacc.c  */
#line 2758 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 675:

/* Line 1455 of yacc.c  */
#line 2766 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(10) - (11)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 676:

/* Line 1455 of yacc.c  */
#line 2777 "p.y"
    { ;}
    break;

  case 677:

/* Line 1455 of yacc.c  */
#line 2779 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 678:

/* Line 1455 of yacc.c  */
#line 2785 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 679:

/* Line 1455 of yacc.c  */
#line 2793 "p.y"
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; ;}
    break;

  case 681:

/* Line 1455 of yacc.c  */
#line 2799 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; ;}
    break;

  case 682:

/* Line 1455 of yacc.c  */
#line 2805 "p.y"
    { ;}
    break;

  case 683:

/* Line 1455 of yacc.c  */
#line 2807 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 684:

/* Line 1455 of yacc.c  */
#line 2813 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 685:

/* Line 1455 of yacc.c  */
#line 2823 "p.y"
    { ;}
    break;

  case 686:

/* Line 1455 of yacc.c  */
#line 2825 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 687:

/* Line 1455 of yacc.c  */
#line 2832 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 688:

/* Line 1455 of yacc.c  */
#line 2839 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 689:

/* Line 1455 of yacc.c  */
#line 2846 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 690:

/* Line 1455 of yacc.c  */
#line 2853 "p.y"
    { /* this one is for frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 691:

/* Line 1455 of yacc.c  */
#line 2861 "p.y"
    { /* this one is for constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 692:

/* Line 1455 of yacc.c  */
#line 2870 "p.y"
    { /* this one is for velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (13)].ival), (yyvsp[(4) - (13)].ival), (yyvsp[(6) - (13)].fval), (yyvsp[(7) - (13)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (13)].ival), (yyvsp[(9) - (13)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (13)].fval), (yyvsp[(11) - (13)].fval), (yyvsp[(12) - (13)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (13)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 693:

/* Line 1455 of yacc.c  */
#line 2879 "p.y"
    { /* this one is for pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (14)].ival), (yyvsp[(4) - (14)].ival), (yyvsp[(6) - (14)].fval), (yyvsp[(7) - (14)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (14)].ival), (yyvsp[(9) - (14)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (14)].fval), (yyvsp[(11) - (14)].fval), (yyvsp[(12) - (14)].fval), (yyvsp[(13) - (14)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (14)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 694:

/* Line 1455 of yacc.c  */
#line 2888 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 695:

/* Line 1455 of yacc.c  */
#line 2896 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(6) - (7)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 696:

/* Line 1455 of yacc.c  */
#line 2904 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(7) - (8)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 697:

/* Line 1455 of yacc.c  */
#line 2912 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival)); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 698:

/* Line 1455 of yacc.c  */
#line 2920 "p.y"
    { /* this one is for frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(10) - (11)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 699:

/* Line 1455 of yacc.c  */
#line 2929 "p.y"
    { /* this one is for constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (12)].ival), (yyvsp[(4) - (12)].ival), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (12)].ival), (yyvsp[(9) - (12)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (12)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (12)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(11) - (12)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 700:

/* Line 1455 of yacc.c  */
#line 2939 "p.y"
    { /* this one is for velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (14)].ival), (yyvsp[(4) - (14)].ival), (yyvsp[(6) - (14)].fval), (yyvsp[(7) - (14)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (14)].ival), (yyvsp[(9) - (14)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (14)].fval), (yyvsp[(11) - (14)].fval), (yyvsp[(12) - (14)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (14)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(13) - (14)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 701:

/* Line 1455 of yacc.c  */
#line 2949 "p.y"
    { /* this one is for pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (15)].ival), (yyvsp[(4) - (15)].ival), (yyvsp[(6) - (15)].fval), (yyvsp[(7) - (15)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (15)].ival), (yyvsp[(9) - (15)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (15)].fval), (yyvsp[(11) - (15)].fval), (yyvsp[(12) - (15)].fval), (yyvsp[(13) - (15)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (15)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(14) - (15)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 702:

/* Line 1455 of yacc.c  */
#line 2961 "p.y"
    { domain->solInfo().dist_acme = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 703:

/* Line 1455 of yacc.c  */
#line 2963 "p.y"
    { domain->solInfo().ffi_debug = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 704:

/* Line 1455 of yacc.c  */
#line 2965 "p.y"
    { domain->solInfo().mortar_scaling = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 705:

/* Line 1455 of yacc.c  */
#line 2967 "p.y"
    { domain->solInfo().mortar_integration_rule = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 706:

/* Line 1455 of yacc.c  */
#line 2970 "p.y"
    { geoSource->addNode((yyvsp[(3) - (3)].nval).num, (yyvsp[(3) - (3)].nval).xyz, (yyvsp[(3) - (3)].nval).cp, (yyvsp[(3) - (3)].nval).cd); ;}
    break;

  case 707:

/* Line 1455 of yacc.c  */
#line 2972 "p.y"
    { geoSource->addNode((yyvsp[(2) - (2)].nval).num, (yyvsp[(2) - (2)].nval).xyz, (yyvsp[(2) - (2)].nval).cp, (yyvsp[(2) - (2)].nval).cd); ;}
    break;

  case 708:

/* Line 1455 of yacc.c  */
#line 2976 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (5)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (5)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (5)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (5)].fval);  (yyval.nval).cp = 0;  (yyval.nval).cd = 0; ;}
    break;

  case 709:

/* Line 1455 of yacc.c  */
#line 2978 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (4)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (4)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (4)].fval);  (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; ;}
    break;

  case 710:

/* Line 1455 of yacc.c  */
#line 2980 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (3)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (3)].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; ;}
    break;

  case 711:

/* Line 1455 of yacc.c  */
#line 2982 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (7)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (7)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (7)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (7)].fval);  (yyval.nval).cp = (yyvsp[(5) - (7)].ival); (yyval.nval).cd = (yyvsp[(6) - (7)].ival);
          if((yyvsp[(5) - (7)].ival) != 0) domain->solInfo().basicPosCoords = false;
          if((yyvsp[(6) - (7)].ival) != 0) domain->solInfo().basicDofCoords = false; ;}
    break;

  case 712:

/* Line 1455 of yacc.c  */
#line 2986 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (6)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (6)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (6)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (6)].fval);  (yyval.nval).cp = (yyvsp[(5) - (6)].ival); (yyval.nval).cd = (yyvsp[(5) - (6)].ival);
          if((yyvsp[(5) - (6)].ival) != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } ;}
    break;

  case 713:

/* Line 1455 of yacc.c  */
#line 2991 "p.y"
    { /* Define each Element */
          geoSource->addElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);;}
    break;

  case 714:

/* Line 1455 of yacc.c  */
#line 2996 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1;;}
    break;

  case 715:

/* Line 1455 of yacc.c  */
#line 2998 "p.y"
    { if((yyval.nl).num == 125) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++;;}
    break;

  case 716:

/* Line 1455 of yacc.c  */
#line 3003 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.bcval).val = (yyvsp[(3) - (4)].fval); (yyval.bcval).mtype = BCond::Axial; ;}
    break;

  case 717:

/* Line 1455 of yacc.c  */
#line 3005 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (3)].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = BCond::Axial; ;}
    break;

  case 718:

/* Line 1455 of yacc.c  */
#line 3007 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (5)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (5)].ival)-1; (yyval.bcval).val = (yyvsp[(3) - (5)].fval); (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[(4) - (5)].ival); ;}
    break;

  case 719:

/* Line 1455 of yacc.c  */
#line 3009 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 720:

/* Line 1455 of yacc.c  */
#line 3013 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 721:

/* Line 1455 of yacc.c  */
#line 3017 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 722:

/* Line 1455 of yacc.c  */
#line 3021 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (5)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (5)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval = (yyvsp[(4) - (5)].fval);  ;}
    break;

  case 723:

/* Line 1455 of yacc.c  */
#line 3023 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (4)].fval); (yyval.cxbcval).imval = 0.0; ;}
    break;

  case 725:

/* Line 1455 of yacc.c  */
#line 3028 "p.y"
    { geoSource->setCSFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 727:

/* Line 1455 of yacc.c  */
#line 3033 "p.y"
    { geoSource->setFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 728:

/* Line 1455 of yacc.c  */
#line 3037 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (11)].ival)-1; 
          (yyval.frame).d[0] = (yyvsp[(2) - (11)].fval); (yyval.frame).d[1] = (yyvsp[(3) - (11)].fval); (yyval.frame).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.frame).d[3] = (yyvsp[(5) - (11)].fval); (yyval.frame).d[4] = (yyvsp[(6) - (11)].fval); (yyval.frame).d[5] = (yyvsp[(7) - (11)].fval);
          (yyval.frame).d[6] = (yyvsp[(8) - (11)].fval); (yyval.frame).d[7] = (yyvsp[(9) - (11)].fval); (yyval.frame).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 729:

/* Line 1455 of yacc.c  */
#line 3042 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (4)].ival)-1;
          geoSource->makeEframe((yyvsp[(1) - (4)].ival)-1, (yyvsp[(3) - (4)].ival), (yyval.frame).d); ;}
    break;

  case 731:

/* Line 1455 of yacc.c  */
#line 3048 "p.y"
    { geoSource->setNodalFrame((yyvsp[(2) - (2)].nframe).id,(yyvsp[(2) - (2)].nframe).o,(yyvsp[(2) - (2)].nframe).d,(yyvsp[(2) - (2)].nframe).type); ;}
    break;

  case 732:

/* Line 1455 of yacc.c  */
#line 3052 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (11)].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;  (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[(2) - (11)].fval); (yyval.nframe).d[1] = (yyvsp[(3) - (11)].fval); (yyval.nframe).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.nframe).d[3] = (yyvsp[(5) - (11)].fval); (yyval.nframe).d[4] = (yyvsp[(6) - (11)].fval); (yyval.nframe).d[5] = (yyvsp[(7) - (11)].fval);
          (yyval.nframe).d[6] = (yyvsp[(8) - (11)].fval); (yyval.nframe).d[7] = (yyvsp[(9) - (11)].fval); (yyval.nframe).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 733:

/* Line 1455 of yacc.c  */
#line 3059 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (14)].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = (yyvsp[(2) - (14)].fval);  (yyval.nframe).o[1] = (yyvsp[(3) - (14)].fval);  (yyval.nframe).o[2] = (yyvsp[(4) - (14)].fval);
          (yyval.nframe).d[0] = (yyvsp[(5) - (14)].fval);  (yyval.nframe).d[1] = (yyvsp[(6) - (14)].fval);  (yyval.nframe).d[2] = (yyvsp[(7) - (14)].fval);
          (yyval.nframe).d[3] = (yyvsp[(8) - (14)].fval);  (yyval.nframe).d[4] = (yyvsp[(9) - (14)].fval);  (yyval.nframe).d[5] = (yyvsp[(10) - (14)].fval);
          (yyval.nframe).d[6] = (yyvsp[(11) - (14)].fval); (yyval.nframe).d[7] = (yyvsp[(12) - (14)].fval); (yyval.nframe).d[8] = (yyvsp[(13) - (14)].fval); ;}
    break;

  case 734:

/* Line 1455 of yacc.c  */
#line 3066 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (12)].ival);
          (yyval.nframe).type = (yyvsp[(2) - (12)].ival);
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;   (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[(3) - (12)].fval); (yyval.nframe).d[1] = (yyvsp[(4) - (12)].fval);  (yyval.nframe).d[2] = (yyvsp[(5) - (12)].fval);
          (yyval.nframe).d[3] = (yyvsp[(6) - (12)].fval); (yyval.nframe).d[4] = (yyvsp[(7) - (12)].fval);  (yyval.nframe).d[5] = (yyvsp[(8) - (12)].fval);
          (yyval.nframe).d[6] = (yyvsp[(9) - (12)].fval); (yyval.nframe).d[7] = (yyvsp[(10) - (12)].fval); (yyval.nframe).d[8] = (yyvsp[(11) - (12)].fval); ;}
    break;

  case 735:

/* Line 1455 of yacc.c  */
#line 3073 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (15)].ival);
          (yyval.nframe).type = (yyvsp[(2) - (15)].ival);
          (yyval.nframe).o[0] = (yyvsp[(3) - (15)].fval);  (yyval.nframe).o[1] = (yyvsp[(4) - (15)].fval);  (yyval.nframe).o[2] = (yyvsp[(5) - (15)].fval);
          (yyval.nframe).d[0] = (yyvsp[(6) - (15)].fval);  (yyval.nframe).d[1] = (yyvsp[(7) - (15)].fval);  (yyval.nframe).d[2] = (yyvsp[(8) - (15)].fval);
          (yyval.nframe).d[3] = (yyvsp[(9) - (15)].fval);  (yyval.nframe).d[4] = (yyvsp[(10) - (15)].fval); (yyval.nframe).d[5] = (yyvsp[(11) - (15)].fval);
          (yyval.nframe).d[6] = (yyvsp[(12) - (15)].fval); (yyval.nframe).d[7] = (yyvsp[(13) - (15)].fval); (yyval.nframe).d[8] = (yyvsp[(14) - (15)].fval); ;}
    break;

  case 737:

/* Line 1455 of yacc.c  */
#line 3083 "p.y"
    { OffsetData od;
	  od.first = (yyvsp[(2) - (7)].ival)-1; od.last = (yyvsp[(3) - (7)].ival)-1;
	  od.o[0] = (yyvsp[(4) - (7)].fval); od.o[1] = (yyvsp[(5) - (7)].fval); od.o[2] = (yyvsp[(6) - (7)].fval); 
	  geoSource->addOffset(od); ;}
    break;

  case 738:

/* Line 1455 of yacc.c  */
#line 3090 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 739:

/* Line 1455 of yacc.c  */
#line 3093 "p.y"
    { geoSource->setElementLumpingWeight((yyvsp[(2) - (5)].ival)-1,(yyvsp[(4) - (5)].fval));
          domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 740:

/* Line 1455 of yacc.c  */
#line 3096 "p.y"
    { geoSource->setElementLumpingWeight((yyvsp[(2) - (6)].ival)-1,(yyvsp[(4) - (6)].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; ;}
    break;

  case 741:

/* Line 1455 of yacc.c  */
#line 3101 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 742:

/* Line 1455 of yacc.c  */
#line 3103 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1); 
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (6)].ival)-1,(yyvsp[(5) - (6)].fval));
	  domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 743:

/* Line 1455 of yacc.c  */
#line 3107 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[(2) - (7)].ival)-1,(yyvsp[(6) - (7)].fval));
          domain->solInfo().elemLumpPodRom = true; 
          domain->solInfo().reduceFollower = true; ;}
    break;

  case 744:

/* Line 1455 of yacc.c  */
#line 3112 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[(2) - (7)].ival)-1,(yyvsp[(5) - (7)].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; ;}
    break;

  case 745:

/* Line 1455 of yacc.c  */
#line 3118 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(4) - (6)].ival)-1,(yyvsp[(5) - (6)].ival)-1); ;}
    break;

  case 746:

/* Line 1455 of yacc.c  */
#line 3120 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (8)].ival)-1,(yyvsp[(3) - (8)].ival)-1,(yyvsp[(4) - (8)].ival)-1,(yyvsp[(5) - (8)].ival)-1);
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (8)].ival)-1,(yyvsp[(7) - (8)].fval)); 
	  domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 747:

/* Line 1455 of yacc.c  */
#line 3124 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (9)].ival)-1,(yyvsp[(3) - (9)].ival)-1,(yyvsp[(4) - (9)].ival)-1,(yyvsp[(5) - (9)].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[(2) - (9)].ival)-1,(yyvsp[(7) - (9)].fval));   
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; ;}
    break;

  case 748:

/* Line 1455 of yacc.c  */
#line 3130 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1,(yyvsp[(4) - (7)].ival)-1,-1,(yyvsp[(6) - (7)].fval)); ;}
    break;

  case 749:

/* Line 1455 of yacc.c  */
#line 3132 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (9)].ival)-1,(yyvsp[(3) - (9)].ival)-1,(yyvsp[(4) - (9)].ival)-1,-1,(yyvsp[(6) - (9)].fval));
          geoSource->setElementLumpingWeight((yyvsp[(2) - (9)].ival)-1,(yyvsp[(8) - (9)].fval));
          domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 750:

/* Line 1455 of yacc.c  */
#line 3136 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (10)].ival)-1,(yyvsp[(3) - (10)].ival)-1,(yyvsp[(4) - (10)].ival)-1,-1,(yyvsp[(6) - (10)].fval)); 
          geoSource->setElementLumpingWeight((yyvsp[(2) - (10)].ival)-1,(yyvsp[(8) - (10)].fval));
          domain->solInfo().elemLumpPodRom = true;
          domain->solInfo().reduceFollower = true; ;}
    break;

  case 751:

/* Line 1455 of yacc.c  */
#line 3142 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        ;}
    break;

  case 752:

/* Line 1455 of yacc.c  */
#line 3147 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
 	    geoSource->setAttrib(i-1,(yyvsp[(4) - (5)].ival)-1);
	;}
    break;

  case 753:

/* Line 1455 of yacc.c  */
#line 3152 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (7)].ival); i<(yyvsp[(3) - (7)].ival)+1; ++i)
	    geoSource->setAttrib(i-1, (yyvsp[(4) - (7)].ival)-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].ival)-1);
	;}
    break;

  case 754:

/* Line 1455 of yacc.c  */
#line 3157 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (8)].ival); i<(yyvsp[(3) - (8)].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[(4) - (8)].ival)-1, (yyvsp[(5) - (8)].ival)-1, -1, (yyvsp[(7) - (8)].fval));
        ;}
    break;

  case 755:

/* Line 1455 of yacc.c  */
#line 3164 "p.y"
    { domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 756:

/* Line 1455 of yacc.c  */
#line 3166 "p.y"
    { geoSource->setElementLumpingWeight((yyvsp[(2) - (4)].ival) - 1, (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 757:

/* Line 1455 of yacc.c  */
#line 3168 "p.y"
    { domain->solInfo().reduceFollower = true;;}
    break;

  case 758:

/* Line 1455 of yacc.c  */
#line 3170 "p.y"
    { domain->solInfo().reduceFollower = true;;}
    break;

  case 759:

/* Line 1455 of yacc.c  */
#line 3174 "p.y"
    { domain->solInfo().ReducedStiffness = true;;}
    break;

  case 760:

/* Line 1455 of yacc.c  */
#line 3176 "p.y"
    { geoSource->pushBackStiffVec((yyvsp[(2) - (3)].fval));;}
    break;

  case 762:

/* Line 1455 of yacc.c  */
#line 3181 "p.y"
    { domain->solInfo().forcePodSize = (yyvsp[(2) - (4)].ival);
          domain->solInfo().maxDeimBasisSize = (yyvsp[(3) - (4)].ival);;}
    break;

  case 763:

/* Line 1455 of yacc.c  */
#line 3184 "p.y"
    { geoSource->pushBackUDEIMVec((yyvsp[(2) - (3)].fval));;}
    break;

  case 765:

/* Line 1455 of yacc.c  */
#line 3189 "p.y"
    { domain->solInfo().DEIMPodRom = true;
          geoSource->setSampleNodesAndSlots((yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].ival));;}
    break;

  case 766:

/* Line 1455 of yacc.c  */
#line 3192 "p.y"
    { geoSource->setSampleElemsAndDOFs((yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].ival));
          domain->solInfo().UDEIMPodRom = true;;}
    break;

  case 767:

/* Line 1455 of yacc.c  */
#line 3198 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 768:

/* Line 1455 of yacc.c  */
#line 3200 "p.y"
    { (yyval.ival) = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 769:

/* Line 1455 of yacc.c  */
#line 3202 "p.y"
    { PressureBCond pbc;
          pbc.setData((yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval), (yyval.ival), true);
          geoSource->setElementPressure(pbc); ;}
    break;

  case 770:

/* Line 1455 of yacc.c  */
#line 3206 "p.y"
    { for(int i = (yyvsp[(2) - (5)].ival); i < ((yyvsp[(3) - (5)].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[(4) - (5)].fval), (yyval.ival), true);
            geoSource->setElementPressure(pbc);
          } ;}
    break;

  case 771:

/* Line 1455 of yacc.c  */
#line 3212 "p.y"
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].fval), (yyval.ival), true);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; ;}
    break;

  case 772:

/* Line 1455 of yacc.c  */
#line 3217 "p.y"
    { PressureBCond pbc;
          pbc.setData((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].fval), (yyval.ival), (yyvsp[(4) - (5)].ival));
          geoSource->setElementPressure(pbc); ;}
    break;

  case 773:

/* Line 1455 of yacc.c  */
#line 3221 "p.y"
    { for(int i = (yyvsp[(2) - (6)].ival); i < ((yyvsp[(3) - (6)].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[(4) - (6)].fval), (yyval.ival), (yyvsp[(5) - (6)].ival));
            geoSource->setElementPressure(pbc);
          } ;}
    break;

  case 774:

/* Line 1455 of yacc.c  */
#line 3227 "p.y"
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].fval), (yyval.ival), (yyvsp[(5) - (6)].ival));
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; ;}
    break;

  case 775:

/* Line 1455 of yacc.c  */
#line 3233 "p.y"
    { PressureBCond pbc;
          pbc.setData((yyvsp[(2) - (6)].ival)-1, (yyvsp[(5) - (6)].fval), (yyval.ival), true);
          pbc.face = (yyvsp[(4) - (6)].ival)-1;
          geoSource->setElementPressure(pbc); ;}
    break;

  case 776:

/* Line 1455 of yacc.c  */
#line 3238 "p.y"
    { for(int i = (yyvsp[(2) - (7)].ival); i < ((yyvsp[(3) - (7)].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[(6) - (7)].fval), (yyval.ival), true);
            pbc.face = (yyvsp[(5) - (7)].ival)-1;
            geoSource->setElementPressure(pbc);
          } ;}
    break;

  case 777:

/* Line 1455 of yacc.c  */
#line 3245 "p.y"
    { PressureBCond pbc;
          pbc.setData((yyvsp[(2) - (7)].ival)-1, (yyvsp[(5) - (7)].fval), (yyval.ival), (yyvsp[(6) - (7)].ival));
          pbc.face = (yyvsp[(4) - (7)].ival)-1;
          geoSource->setElementPressure(pbc); ;}
    break;

  case 778:

/* Line 1455 of yacc.c  */
#line 3250 "p.y"
    { for(int i = (yyvsp[(2) - (8)].ival); i < ((yyvsp[(3) - (8)].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[(6) - (8)].fval), (yyval.ival), (yyvsp[(7) - (8)].ival));
            pbc.face = (yyvsp[(5) - (8)].ival)-1;
            geoSource->setElementPressure(pbc);
          } ;}
    break;

  case 779:

/* Line 1455 of yacc.c  */
#line 3259 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        ;}
    break;

  case 780:

/* Line 1455 of yacc.c  */
#line 3264 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[(2) - (3)].ival));
          geoSource->setConsistentPFlag(false);
        ;}
    break;

  case 781:

/* Line 1455 of yacc.c  */
#line 3278 "p.y"
    { ;}
    break;

  case 782:

/* Line 1455 of yacc.c  */
#line 3280 "p.y"
    { geoSource->setElementPreLoad( (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval) ); ;}
    break;

  case 783:

/* Line 1455 of yacc.c  */
#line 3282 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (6)].ival); i<((yyvsp[(4) - (6)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[(5) - (6)].fval) );
        ;}
    break;

  case 784:

/* Line 1455 of yacc.c  */
#line 3287 "p.y"
    { double load[3] = { (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval) };
          geoSource->setElementPreLoad( (yyvsp[(2) - (6)].ival)-1, load ); ;}
    break;

  case 785:

/* Line 1455 of yacc.c  */
#line 3290 "p.y"
    { double load[3] = { (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
          int i;
          for(i=(yyvsp[(2) - (8)].ival); i<((yyvsp[(4) - (8)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        ;}
    break;

  case 786:

/* Line 1455 of yacc.c  */
#line 3298 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 790:

/* Line 1455 of yacc.c  */
#line 3303 "p.y"
    { // activate piecewise constant configuration dependent external forces for a linear dynamic analysis
          if(!domain->solInfo().isNonLin()) { 
            if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
              domain->solInfo().probType = SolverInfo::NonLinStatic;
            else if(domain->solInfo().probType == SolverInfo::Dynamic)
              domain->solInfo().probType = SolverInfo::NonLinDynam;
            else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
              domain->solInfo().order = 1;
              domain->solInfo().probType = SolverInfo::NonLinDynam;
            }
            domain->solInfo().setNewton(std::numeric_limits<int>::max());
            domain->solInfo().getNLInfo().stepUpdateK = std::numeric_limits<int>::max();
            domain->solInfo().getNLInfo().linearelastic = true;
            domain->solInfo().getNLInfo().maxiter = 1;
          }
        ;}
    break;

  case 791:

/* Line 1455 of yacc.c  */
#line 3320 "p.y"
    { // activate piecewise constant configuration dependent external forces for a linear static analysis
          if(!domain->solInfo().isNonLin()) {  
            if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
              domain->solInfo().probType = SolverInfo::NonLinStatic;
            else if(domain->solInfo().probType == SolverInfo::Dynamic)
              domain->solInfo().probType = SolverInfo::NonLinDynam;
            else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
              domain->solInfo().order = 1;
              domain->solInfo().probType = SolverInfo::NonLinDynam;
            }
            domain->solInfo().setNewton(std::numeric_limits<int>::max());
            domain->solInfo().getNLInfo().stepUpdateK = std::numeric_limits<int>::max();
            domain->solInfo().getNLInfo().linearelastic = true;
            domain->solInfo().getNLInfo().maxiter = 1;
            domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (5)].fval);
            domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (5)].fval);
          }
        ;}
    break;

  case 792:

/* Line 1455 of yacc.c  */
#line 3341 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(1) - (1)].ival)); ;}
    break;

  case 793:

/* Line 1455 of yacc.c  */
#line 3343 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 794:

/* Line 1455 of yacc.c  */
#line 3347 "p.y"
    { domain->solInfo().type = 1;
          domain->solInfo().iterType = (yyvsp[(1) - (2)].ival); ;}
    break;

  case 795:

/* Line 1455 of yacc.c  */
#line 3350 "p.y"
    { domain->solInfo().precond = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 796:

/* Line 1455 of yacc.c  */
#line 3352 "p.y"
    { domain->solInfo().maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 797:

/* Line 1455 of yacc.c  */
#line 3354 "p.y"
    { domain->solInfo().tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 798:

/* Line 1455 of yacc.c  */
#line 3356 "p.y"
    { domain->solInfo().maxvecsize = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 799:

/* Line 1455 of yacc.c  */
#line 3358 "p.y"
    { domain->solInfo().iterSubtype = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 800:

/* Line 1455 of yacc.c  */
#line 3362 "p.y"
    { domain->solInfo().setSolver(0); ;}
    break;

  case 801:

/* Line 1455 of yacc.c  */
#line 3364 "p.y"
    { domain->solInfo().setSolver((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 802:

/* Line 1455 of yacc.c  */
#line 3366 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (2)].ival)); ;}
    break;

  case 803:

/* Line 1455 of yacc.c  */
#line 3368 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (3)].ival));
          if((yyvsp[(1) - (3)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 804:

/* Line 1455 of yacc.c  */
#line 3372 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (3)].ival));
          domain->solInfo().getNLInfo().unsymmetric = true; ;}
    break;

  case 805:

/* Line 1455 of yacc.c  */
#line 3375 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (3)].ival),(yyvsp[(2) - (3)].ival)); ;}
    break;

  case 806:

/* Line 1455 of yacc.c  */
#line 3377 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 807:

/* Line 1455 of yacc.c  */
#line 3379 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (5)].ival),(yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].ival)); ;}
    break;

  case 808:

/* Line 1455 of yacc.c  */
#line 3381 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (6)].ival),(yyvsp[(2) - (6)].ival),(yyvsp[(3) - (6)].fval),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].ival)); ;}
    break;

  case 809:

/* Line 1455 of yacc.c  */
#line 3383 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (7)].ival),(yyvsp[(2) - (7)].ival),(yyvsp[(3) - (7)].fval),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].ival),(yyvsp[(6) - (7)].ival)); ;}
    break;

  case 810:

/* Line 1455 of yacc.c  */
#line 3385 "p.y"
    { domain->solInfo().fetiInfo.maxit    = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(3) - (4)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(2) - (4)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 811:

/* Line 1455 of yacc.c  */
#line 3390 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (3)].ival)-1); ;}
    break;

  case 812:

/* Line 1455 of yacc.c  */
#line 3393 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; ;}
    break;

  case 813:

/* Line 1455 of yacc.c  */
#line 3397 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; ;}
    break;

  case 814:

/* Line 1455 of yacc.c  */
#line 3402 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (4)].ival)-1); 
          domain->solInfo().fetiInfo.feti2version 
                  = (FetiInfo::Feti2Version) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 815:

/* Line 1455 of yacc.c  */
#line 3407 "p.y"
    { domain->solInfo().fetiInfo.maxit    = (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(3) - (5)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(4) - (5)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 816:

/* Line 1455 of yacc.c  */
#line 3412 "p.y"
    { domain->solInfo().type =(2); ;}
    break;

  case 817:

/* Line 1455 of yacc.c  */
#line 3414 "p.y"
    { domain->solInfo().type = 3;
          domain->solInfo().subtype = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getFetiInfo().solvertype = (FetiInfo::Solvertype)((yyvsp[(3) - (4)].ival));
	;}
    break;

  case 818:

/* Line 1455 of yacc.c  */
#line 3419 "p.y"
    { domain->solInfo().sparse_maxsup  = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 819:

/* Line 1455 of yacc.c  */
#line 3421 "p.y"
    { domain->solInfo().sparse_defblk  = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 820:

/* Line 1455 of yacc.c  */
#line 3423 "p.y"
    { domain->solInfo().spooles_tau  = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 821:

/* Line 1455 of yacc.c  */
#line 3425 "p.y"
    { domain->solInfo().spooles_maxsize = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 822:

/* Line 1455 of yacc.c  */
#line 3427 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 0) {
            (yyvsp[(2) - (3)].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          domain->solInfo().spooles_maxdomainsize = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 823:

/* Line 1455 of yacc.c  */
#line 3434 "p.y"
    { domain->solInfo().spooles_seed = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 824:

/* Line 1455 of yacc.c  */
#line 3436 "p.y"
    { if(((yyvsp[(2) - (3)].fval) < 0.0) || ((yyvsp[(2) - (3)].fval) > 1.0)) {
            (yyvsp[(2) - (3)].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          domain->solInfo().spooles_maxzeros = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 825:

/* Line 1455 of yacc.c  */
#line 3443 "p.y"
    { domain->solInfo().spooles_msglvl = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 826:

/* Line 1455 of yacc.c  */
#line 3445 "p.y"
    { domain->solInfo().pivot = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 827:

/* Line 1455 of yacc.c  */
#line 3447 "p.y"
    { domain->solInfo().spooles_scale = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 828:

/* Line 1455 of yacc.c  */
#line 3449 "p.y"
    { domain->solInfo().spooles_renum = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 829:

/* Line 1455 of yacc.c  */
#line 3451 "p.y"
    { domain->solInfo().mumps_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 830:

/* Line 1455 of yacc.c  */
#line 3453 "p.y"
    { domain->solInfo().mumps_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 831:

/* Line 1455 of yacc.c  */
#line 3455 "p.y"
    { domain->solInfo().goldfarb_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 832:

/* Line 1455 of yacc.c  */
#line 3457 "p.y"
    { domain->solInfo().goldfarb_check = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 833:

/* Line 1455 of yacc.c  */
#line 3459 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 834:

/* Line 1455 of yacc.c  */
#line 3461 "p.y"
    { domain->solInfo().debug_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 835:

/* Line 1455 of yacc.c  */
#line 3463 "p.y"
    { domain->solInfo().debug_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 836:

/* Line 1455 of yacc.c  */
#line 3469 "p.y"
    { domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 837:

/* Line 1455 of yacc.c  */
#line 3471 "p.y"
    { domain->solInfo().fetiInfo.precno = FetiInfo::lumped; ;}
    break;

  case 838:

/* Line 1455 of yacc.c  */
#line 3473 "p.y"
    { if(((yyvsp[(3) - (4)].ival) < 0) || ((yyvsp[(3) - (4)].ival) > 3)) { 
            (yyvsp[(3) - (4)].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival);
	;}
    break;

  case 839:

/* Line 1455 of yacc.c  */
#line 3480 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        ;}
    break;

  case 840:

/* Line 1455 of yacc.c  */
#line 3487 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        ;}
    break;

  case 841:

/* Line 1455 of yacc.c  */
#line 3494 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 842:

/* Line 1455 of yacc.c  */
#line 3496 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.absolute_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 843:

/* Line 1455 of yacc.c  */
#line 3499 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 844:

/* Line 1455 of yacc.c  */
#line 3501 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.absolute_stagnation_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 845:

/* Line 1455 of yacc.c  */
#line 3504 "p.y"
    { domain->solInfo().fetiInfo.primal_proj_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_proj_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 846:

/* Line 1455 of yacc.c  */
#line 3507 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.dual_plan_maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 847:

/* Line 1455 of yacc.c  */
#line 3510 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_plan_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 848:

/* Line 1455 of yacc.c  */
#line 3513 "p.y"
    { domain->solInfo().fetiInfo.maxortho = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 849:

/* Line 1455 of yacc.c  */
#line 3515 "p.y"
    { domain->solInfo().fetiInfo.noCoarse = 1; ;}
    break;

  case 850:

/* Line 1455 of yacc.c  */
#line 3517 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 1) 
            domain->solInfo().fetiInfo.nonLocalQ = 0;
          else if((yyvsp[(2) - (3)].ival) == 2) {
            domain->solInfo().fetiInfo.nonLocalQ = 1;
            if(domain->solInfo().fetiInfo.version == FetiInfo::feti2) {
              domain->solInfo().fetiInfo.nonLocalQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used"
                             " with FETI 2\n");
            }
          } else if((yyvsp[(2) - (3)].ival) == 3) {
            domain->solInfo().fetiInfo.nonLocalQ = 1;
            domain->solInfo().fetiInfo.nQ = 3;
          } else if((yyvsp[(2) - (3)].ival) == 4) {
            domain->solInfo().fetiInfo.nonLocalQ = 1;
            domain->solInfo().fetiInfo.nQ = 4;
          } else
            fprintf(stderr," *** WARNING: This projector does not exist,"
                           " using basic projector\n");
        ;}
    break;

  case 851:

/* Line 1455 of yacc.c  */
#line 3537 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 1; 
          domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 852:

/* Line 1455 of yacc.c  */
#line 3540 "p.y"
    { domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 853:

/* Line 1455 of yacc.c  */
#line 3542 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 854:

/* Line 1455 of yacc.c  */
#line 3545 "p.y"
    { domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 855:

/* Line 1455 of yacc.c  */
#line 3547 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 856:

/* Line 1455 of yacc.c  */
#line 3550 "p.y"
    { domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 857:

/* Line 1455 of yacc.c  */
#line 3552 "p.y"
    { domain->solInfo().fetiInfo.mpc_element = true; ;}
    break;

  case 858:

/* Line 1455 of yacc.c  */
#line 3554 "p.y"
    { domain->solInfo().fetiInfo.fsi_element = true; ;}
    break;

  case 859:

/* Line 1455 of yacc.c  */
#line 3556 "p.y"
    { domain->solInfo().fetiInfo.fsi_corner = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 860:

/* Line 1455 of yacc.c  */
#line 3558 "p.y"
    { domain->solInfo().fetiInfo.splitLocalFsi = false; ;}
    break;

  case 861:

/* Line 1455 of yacc.c  */
#line 3560 "p.y"
    { domain->solInfo().coupled_scale = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 862:

/* Line 1455 of yacc.c  */
#line 3562 "p.y"
    { domain->solInfo().fetiInfo.wetcorners = true; ;}
    break;

  case 863:

/* Line 1455 of yacc.c  */
#line 3564 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 864:

/* Line 1455 of yacc.c  */
#line 3566 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.pick_unsafe_corners = bool((yyvsp[(3) - (4)].ival));
        ;}
    break;

  case 865:

/* Line 1455 of yacc.c  */
#line 3570 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) {
            domain->solInfo().fetiInfo.corners = FetiInfo::noCorners;
            domain->solInfo().fetiInfo.pickAnyCorner = 0; 
            domain->solInfo().fetiInfo.bmpc = true;
            domain->solInfo().fetiInfo.pick_unsafe_corners = false;
            domain->solInfo().fetiInfo.augment = FetiInfo::none;
          }
        ;}
    break;

  case 866:

/* Line 1455 of yacc.c  */
#line 3579 "p.y"
    {
          if(domain->solInfo().fetiInfo.dph_flag && ((yyvsp[(2) - (3)].ival) == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          }
          else domain->solInfo().fetiInfo.augment = (FetiInfo::AugmentType) (yyvsp[(2) - (3)].ival);

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
        ;}
    break;

  case 867:

/* Line 1455 of yacc.c  */
#line 3600 "p.y"
    {
          if(domain->solInfo().fetiInfo.dph_flag && ((yyvsp[(2) - (4)].ival) == 1)) {
            std::cerr << "WARNING: Selected augment type is unsupported for FETI-DPH, set to EdgeGs \n";
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          }
          else domain->solInfo().fetiInfo.augment = (FetiInfo::AugmentType) (yyvsp[(2) - (4)].ival);
          if(domain->solInfo().fetiInfo.dph_flag && ((yyvsp[(3) - (4)].ival) > 2) && ((yyvsp[(3) - (4)].ival) < 6)) {
            std::cerr << "WARNING: Selected rbm type is unsupported for FETI-DPH, set to translation \n";
            domain->solInfo().fetiInfo.rbmType = FetiInfo::translation;
          }
          else domain->solInfo().fetiInfo.rbmType = (FetiInfo::RbmType) (yyvsp[(3) - (4)].ival);

          if(domain->solInfo().fetiInfo.rbmType == FetiInfo::all)
            domain->solInfo().fetiInfo.nGs = 6;
          else if(domain->solInfo().fetiInfo.rbmType != FetiInfo::None)
            domain->solInfo().fetiInfo.nGs = 3;
        ;}
    break;

  case 868:

/* Line 1455 of yacc.c  */
#line 3618 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (4)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 869:

/* Line 1455 of yacc.c  */
#line 3623 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (5)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 870:

/* Line 1455 of yacc.c  */
#line 3629 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(4) - (5)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 871:

/* Line 1455 of yacc.c  */
#line 3635 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (6)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(5) - (6)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (6)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 872:

/* Line 1455 of yacc.c  */
#line 3642 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 873:

/* Line 1455 of yacc.c  */
#line 3644 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.orthotol2 = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 874:

/* Line 1455 of yacc.c  */
#line 3647 "p.y"
    { domain->solInfo().fetiInfo.grbm_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 875:

/* Line 1455 of yacc.c  */
#line 3649 "p.y"
    { domain->solInfo().fetiInfo.crbm_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 876:

/* Line 1455 of yacc.c  */
#line 3651 "p.y"
    { domain->solInfo().fetiInfo.cct_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 877:

/* Line 1455 of yacc.c  */
#line 3653 "p.y"
    { domain->solInfo().fetiInfo.rebuildcct = int((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 878:

/* Line 1455 of yacc.c  */
#line 3655 "p.y"
    { domain->solInfo().fetiInfo.uproj = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 879:

/* Line 1455 of yacc.c  */
#line 3657 "p.y"
    { domain->solInfo().fetiInfo.printMatLab = 1; ;}
    break;

  case 880:

/* Line 1455 of yacc.c  */
#line 3659 "p.y"
    { domain->solInfo().printMatLab = 1;
          domain->solInfo().printMatLabFile = (yyvsp[(2) - (3)].strval); ;}
    break;

  case 881:

/* Line 1455 of yacc.c  */
#line 3662 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 882:

/* Line 1455 of yacc.c  */
#line 3664 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 883:

/* Line 1455 of yacc.c  */
#line 3666 "p.y"
    {  domain->solInfo().fetiInfo.auxCoarseSolver = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 884:

/* Line 1455 of yacc.c  */
#line 3668 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 885:

/* Line 1455 of yacc.c  */
#line 3670 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 886:

/* Line 1455 of yacc.c  */
#line 3674 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().localScaled = true; ;}
    break;

  case 887:

/* Line 1455 of yacc.c  */
#line 3677 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 888:

/* Line 1455 of yacc.c  */
#line 3681 "p.y"
    { domain->solInfo().fetiInfo.auxCoarseSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 889:

/* Line 1455 of yacc.c  */
#line 3685 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().coarseScaled = true; ;}
    break;

  case 890:

/* Line 1455 of yacc.c  */
#line 3688 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 891:

/* Line 1455 of yacc.c  */
#line 3692 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival)!=0) fprintf(stderr," *** WARNING: Scaling not supported for this CCt solver \n");
          else domain->solInfo().fetiInfo.cctScaled = true; ;}
    break;

  case 892:

/* Line 1455 of yacc.c  */
#line 3696 "p.y"
    { 
          if((yyvsp[(2) - (3)].ival) == 1)
            domain->solInfo().fetiInfo.version = FetiInfo::feti1; 
          else if((yyvsp[(2) - (3)].ival) == 2) {
            domain->solInfo().fetiInfo.version = FetiInfo::feti2;
            if(domain->solInfo().fetiInfo.nonLocalQ == 1) {
              domain->solInfo().fetiInfo.nonLocalQ = 0;
              domain->solInfo().fetiInfo.nQ = 0;
              fprintf(stderr," *** WARNING: Basic projector is used "
                             "with FETI 2\n");
            }
          } else if((yyvsp[(2) - (3)].ival) == 3) {
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
	;}
    break;

  case 893:

/* Line 1455 of yacc.c  */
#line 3722 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(1) - (2)].ival); ;}
    break;

  case 894:

/* Line 1455 of yacc.c  */
#line 3724 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = true; ;}
    break;

  case 895:

/* Line 1455 of yacc.c  */
#line 3726 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 896:

/* Line 1455 of yacc.c  */
#line 3728 "p.y"
    { domain->solInfo().fetiInfo.pickAnyCorner = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 897:

/* Line 1455 of yacc.c  */
#line 3734 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().fetiInfo.nlPrecFlg = 1; 
          domain->solInfo().setKrylov(); 
        ;}
    break;

  case 898:

/* Line 1455 of yacc.c  */
#line 3739 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
	  domain->solInfo().fetiInfo.nlPrecFlg = (yyvsp[(2) - (3)].ival);
	  domain->solInfo().setKrylov();
	;}
    break;

  case 899:

/* Line 1455 of yacc.c  */
#line 3744 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.tol = (yyvsp[(3) - (4)].fval);
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        ;}
    break;

  case 900:

/* Line 1455 of yacc.c  */
#line 3757 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.tol = (yyvsp[(3) - (5)].fval);
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.numcgm = (yyvsp[(4) - (5)].ival);
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        ;}
    break;

  case 901:

/* Line 1455 of yacc.c  */
#line 3771 "p.y"
    {
          domain->solInfo().fetiInfo.maxit = (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.tol = (yyvsp[(3) - (6)].fval);
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.numcgm = (yyvsp[(4) - (6)].ival);
          domain->solInfo().fetiInfo.tolcgm = (yyvsp[(5) - (6)].fval);
          domain->solInfo().fetiInfo.spaceDimension = 2;
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        ;}
    break;

  case 902:

/* Line 1455 of yacc.c  */
#line 3788 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(2) - (7)].ival);
          domain->solInfo().fetiInfo.tol = (yyvsp[(3) - (7)].fval);
          domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.numcgm = (yyvsp[(4) - (7)].ival);
          domain->solInfo().fetiInfo.tolcgm = (yyvsp[(5) - (7)].fval);
          domain->solInfo().fetiInfo.spaceDimension = (yyvsp[(6) - (7)].ival);
          domain->solInfo().fetiInfo.krylovtype = 1;
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling;
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        ;}
    break;

  case 903:

/* Line 1455 of yacc.c  */
#line 3804 "p.y"
    { domain->solInfo().type =(2); 
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        ;}
    break;

  case 904:

/* Line 1455 of yacc.c  */
#line 3814 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 905:

/* Line 1455 of yacc.c  */
#line 3816 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.numcgm2 = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 906:

/* Line 1455 of yacc.c  */
#line 3819 "p.y"
    { domain->solInfo().fetiInfo.tolcgm = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 907:

/* Line 1455 of yacc.c  */
#line 3821 "p.y"
    { domain->solInfo().fetiInfo.spaceDimension = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 908:

/* Line 1455 of yacc.c  */
#line 3823 "p.y"
    { domain->solInfo().fetiInfo.krylovtype = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 909:

/* Line 1455 of yacc.c  */
#line 3825 "p.y"
    { domain->solInfo().fetiInfo.krylovtype =  (yyvsp[(2) - (3)].ival); ;}
    break;

  case 910:

/* Line 1455 of yacc.c  */
#line 3827 "p.y"
    { domain->solInfo().fetiInfo.lumpedinterface = 1; ;}
    break;

  case 911:

/* Line 1455 of yacc.c  */
#line 3829 "p.y"
    { domain->solInfo().fetiInfo.saveMemCoarse = 1; ;}
    break;

  case 912:

/* Line 1455 of yacc.c  */
#line 3831 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        ;}
    break;

  case 913:

/* Line 1455 of yacc.c  */
#line 3836 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        ;}
    break;

  case 914:

/* Line 1455 of yacc.c  */
#line 3842 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        ;}
    break;

  case 915:

/* Line 1455 of yacc.c  */
#line 3849 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 916:

/* Line 1455 of yacc.c  */
#line 3851 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.complex_hermitian = true; ;}
    break;

  case 917:

/* Line 1455 of yacc.c  */
#line 3854 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 918:

/* Line 1455 of yacc.c  */
#line 3856 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 919:

/* Line 1455 of yacc.c  */
#line 3858 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 920:

/* Line 1455 of yacc.c  */
#line 3860 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 921:

/* Line 1455 of yacc.c  */
#line 3862 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 922:

/* Line 1455 of yacc.c  */
#line 3865 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 923:

/* Line 1455 of yacc.c  */
#line 3868 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 924:

/* Line 1455 of yacc.c  */
#line 3871 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 925:

/* Line 1455 of yacc.c  */
#line 3874 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); ;}
    break;

  case 926:

/* Line 1455 of yacc.c  */
#line 3878 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); ;}
    break;

  case 927:

/* Line 1455 of yacc.c  */
#line 3882 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 1) domain->solInfo().fetiInfo.useMRHS = false; ;}
    break;

  case 928:

/* Line 1455 of yacc.c  */
#line 3884 "p.y"
    { domain->solInfo().fetiInfo.gamma = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 929:

/* Line 1455 of yacc.c  */
#line 3886 "p.y"
    { domain->solInfo().fetiInfo.linesearch_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.linesearch_tau = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 930:

/* Line 1455 of yacc.c  */
#line 3889 "p.y"
    { domain->solInfo().fetiInfo.bmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 931:

/* Line 1455 of yacc.c  */
#line 3891 "p.y"
    { domain->solInfo().fetiInfo.dmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 932:

/* Line 1455 of yacc.c  */
#line 3893 "p.y"
    { domain->solInfo().fetiInfo.cmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 933:

/* Line 1455 of yacc.c  */
#line 3895 "p.y"
    { domain->solInfo().fetiInfo.c_normalize = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 934:

/* Line 1455 of yacc.c  */
#line 3897 "p.y"
    { domain->solInfo().dbccheck = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 936:

/* Line 1455 of yacc.c  */
#line 3904 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp; 
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        ;}
    break;

  case 937:

/* Line 1455 of yacc.c  */
#line 3915 "p.y"
    { if(!(yyvsp[(2) - (3)].copt).lagrangeMult && (yyvsp[(2) - (3)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(2) - (3)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(2) - (3)].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[(2) - (3)].copt).constraint_hess; 
          domain->solInfo().constraint_hess_eps = (yyvsp[(2) - (3)].copt).constraint_hess_eps; ;}
    break;

  case 938:

/* Line 1455 of yacc.c  */
#line 3921 "p.y"
    { if(!(yyvsp[(3) - (4)].copt).lagrangeMult && (yyvsp[(3) - (4)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(3) - (4)].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[(3) - (4)].copt).constraint_hess;
          domain->solInfo().constraint_hess_eps = (yyvsp[(3) - (4)].copt).constraint_hess_eps; ;}
    break;

  case 939:

/* Line 1455 of yacc.c  */
#line 3929 "p.y"
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        ;}
    break;

  case 940:

/* Line 1455 of yacc.c  */
#line 3936 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 941:

/* Line 1455 of yacc.c  */
#line 3943 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (3)].fval);
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 942:

/* Line 1455 of yacc.c  */
#line 3951 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (4)].fval);
          domain->solInfo().rhsZeroTol = (yyvsp[(4) - (4)].fval); ;}
    break;

  case 943:

/* Line 1455 of yacc.c  */
#line 3960 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (5)].fval);
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (5)].fval); 
          domain->solInfo().rhsZeroTol = (yyvsp[(4) - (5)].fval);
          domain->solInfo().inconsistentTol = (yyvsp[(5) - (5)].fval); ;}
    break;

  case 944:

/* Line 1455 of yacc.c  */
#line 3970 "p.y"
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 945:

/* Line 1455 of yacc.c  */
#line 3976 "p.y"
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[(2) - (2)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 946:

/* Line 1455 of yacc.c  */
#line 3982 "p.y"
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[(3) - (3)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 947:

/* Line 1455 of yacc.c  */
#line 3988 "p.y"
    { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[(2) - (2)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 948:

/* Line 1455 of yacc.c  */
#line 3994 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (3)].ival);
          (yyval.copt).constraint_hess_eps = 0; ;}
    break;

  case 949:

/* Line 1455 of yacc.c  */
#line 3997 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (4)].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[(4) - (4)].fval); ;}
    break;

  case 950:

/* Line 1455 of yacc.c  */
#line 4002 "p.y"
    { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        ;}
    break;

  case 951:

/* Line 1455 of yacc.c  */
#line 4007 "p.y"
    { domain->solInfo().type = (2); 
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        ;}
    break;

  case 952:

/* Line 1455 of yacc.c  */
#line 4017 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        ;}
    break;

  case 953:

/* Line 1455 of yacc.c  */
#line 4022 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        ;}
    break;

  case 954:

/* Line 1455 of yacc.c  */
#line 4028 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        ;}
    break;

  case 955:

/* Line 1455 of yacc.c  */
#line 4035 "p.y"
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        ;}
    break;

  case 956:

/* Line 1455 of yacc.c  */
#line 4040 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        ;}
    break;

  case 957:

/* Line 1455 of yacc.c  */
#line 4045 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        ;}
    break;

  case 963:

/* Line 1455 of yacc.c  */
#line 4059 "p.y"
    { domain->setWaveDirections(0, (yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 964:

/* Line 1455 of yacc.c  */
#line 4062 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 965:

/* Line 1455 of yacc.c  */
#line 4065 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        ;}
    break;

  case 968:

/* Line 1455 of yacc.c  */
#line 4075 "p.y"
    { domain->setFFPDirections((yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 970:

/* Line 1455 of yacc.c  */
#line 4082 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        ;}
    break;

  case 972:

/* Line 1455 of yacc.c  */
#line 4096 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        ;}
    break;

  case 973:

/* Line 1455 of yacc.c  */
#line 4107 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        ;}
    break;

  case 974:

/* Line 1455 of yacc.c  */
#line 4113 "p.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 975:

/* Line 1455 of yacc.c  */
#line 4115 "p.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 976:

/* Line 1455 of yacc.c  */
#line 4118 "p.y"
    { if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
            domain->solInfo().probType = SolverInfo::NonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::Dynamic)
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
            domain->solInfo().order = 1;
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          }
          domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().getNLInfo().setDefaults(); /* just in case PIECEWISE is used under statics */ ;}
    break;

  case 977:

/* Line 1455 of yacc.c  */
#line 4129 "p.y"
    {;}
    break;

  case 978:

/* Line 1455 of yacc.c  */
#line 4131 "p.y"
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength; ;}
    break;

  case 979:

/* Line 1455 of yacc.c  */
#line 4134 "p.y"
    { if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam; ;}
    break;

  case 980:

/* Line 1455 of yacc.c  */
#line 4139 "p.y"
    { domain->solInfo().getNLInfo().linearelastic = true; ;}
    break;

  case 981:

/* Line 1455 of yacc.c  */
#line 4141 "p.y"
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 982:

/* Line 1455 of yacc.c  */
#line 4143 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 983:

/* Line 1455 of yacc.c  */
#line 4145 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 984:

/* Line 1455 of yacc.c  */
#line 4148 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (7)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[(5) - (7)].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[(6) - (7)].fval); ;}
    break;

  case 985:

/* Line 1455 of yacc.c  */
#line 4153 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 986:

/* Line 1455 of yacc.c  */
#line 4156 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (7)].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[(5) - (7)].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[(6) - (7)].ival); ;}
    break;

  case 987:

/* Line 1455 of yacc.c  */
#line 4161 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 988:

/* Line 1455 of yacc.c  */
#line 4164 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (5)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 989:

/* Line 1455 of yacc.c  */
#line 4167 "p.y"
    { domain->solInfo().getNLInfo().unsymmetric = true; ;}
    break;

  case 990:

/* Line 1455 of yacc.c  */
#line 4169 "p.y"
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 991:

/* Line 1455 of yacc.c  */
#line 4171 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 992:

/* Line 1455 of yacc.c  */
#line 4173 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (5)].ival); 
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 993:

/* Line 1455 of yacc.c  */
#line 4176 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (6)].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (6)].ival);
          // note: currently we use either c1 or c2, but never both
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[(5) - (6)].fval);
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[(5) - (6)].fval); ;}
    break;

  case 994:

/* Line 1455 of yacc.c  */
#line 4182 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (7)].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (7)].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[(5) - (7)].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[(5) - (7)].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[(6) - (7)].fval); ;}
    break;

  case 995:

/* Line 1455 of yacc.c  */
#line 4188 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (8)].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (8)].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[(5) - (8)].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[(5) - (8)].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[(6) - (8)].fval);
          domain->solInfo().getNLInfo().linesearch.verbose = bool((yyvsp[(7) - (8)].ival)); ;}
    break;

  case 996:

/* Line 1455 of yacc.c  */
#line 4195 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true; ;}
    break;

  case 997:

/* Line 1455 of yacc.c  */
#line 4197 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 998:

/* Line 1455 of yacc.c  */
#line 4200 "p.y"
    { domain->solInfo().num_penalty_its = (yyvsp[(3) - (6)].ival); 
          domain->solInfo().penalty_tol = (yyvsp[(4) - (6)].fval);
          domain->solInfo().penalty_beta = (yyvsp[(5) - (6)].fval); ;}
    break;

  case 1000:

/* Line 1455 of yacc.c  */
#line 4207 "p.y"
    { 
          domain->solInfo().setNewton((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear; 
        ;}
    break;

  case 1001:

/* Line 1455 of yacc.c  */
#line 4212 "p.y"
    {
          domain->solInfo().setNewton((yyvsp[(2) - (4)].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[(3) - (4)].ival);
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        ;}
    break;

  case 1002:

/* Line 1455 of yacc.c  */
#line 4249 "p.y"
    { domain->solInfo().setReOrtho(); ;}
    break;

  case 1004:

/* Line 1455 of yacc.c  */
#line 4254 "p.y"
    { geoSource->setControl((yyvsp[(3) - (10)].strval),(yyvsp[(7) - (10)].strval),(yyvsp[(9) - (10)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (10)].ival); ;}
    break;

  case 1005:

/* Line 1455 of yacc.c  */
#line 4256 "p.y"
    { geoSource->setControl((yyvsp[(3) - (12)].strval),(yyvsp[(7) - (12)].strval),(yyvsp[(9) - (12)].strval),(yyvsp[(11) - (12)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (12)].ival); ;}
    break;

  case 1006:

/* Line 1455 of yacc.c  */
#line 4266 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(2) - (3)].strval);
#endif
        ;}
    break;

  case 1007:

/* Line 1455 of yacc.c  */
#line 4272 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(3) - (4)].strval);
#endif
 ;}
    break;

  case 1009:

/* Line 1455 of yacc.c  */
#line 4281 "p.y"
    { domain->solInfo().contact_mode = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 1010:

/* Line 1455 of yacc.c  */
#line 4285 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (7)].ival)-1, (yyvsp[(3) - (7)].ival)-1, (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); ;}
    break;

  case 1011:

/* Line 1455 of yacc.c  */
#line 4288 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(8) - (9)].fval));;}
    break;

  case 1012:

/* Line 1455 of yacc.c  */
#line 4291 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), 0.0, (yyvsp[(8) - (9)].ival));;}
    break;

  case 1013:

/* Line 1455 of yacc.c  */
#line 4293 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (8)].ival)-1, (yyvsp[(3) - (8)].ival)-1, (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), 0.0, -1, (yyvsp[(7) - (8)].copt).lagrangeMult, (yyvsp[(7) - (8)].copt).penalty);;}
    break;

  case 1014:

/* Line 1455 of yacc.c  */
#line 4296 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (11)].ival)-1, (yyvsp[(3) - (11)].ival)-1, (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(10) - (11)].ival));;}
    break;

  case 1015:

/* Line 1455 of yacc.c  */
#line 4298 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (12)].ival)-1, (yyvsp[(3) - (12)].ival)-1, (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(10) - (12)].ival), (yyvsp[(11) - (12)].copt).lagrangeMult, (yyvsp[(11) - (12)].copt).penalty);;}
    break;

  case 1017:

/* Line 1455 of yacc.c  */
#line 4303 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 1018:

/* Line 1455 of yacc.c  */
#line 4308 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new BilinPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 1019:

/* Line 1455 of yacc.c  */
#line 4313 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
             new BilinPlasKinHardMat((yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval)) );
         ;}
    break;

  case 1020:

/* Line 1455 of yacc.c  */
#line 4318 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 1021:

/* Line 1455 of yacc.c  */
#line 4323 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 1022:

/* Line 1455 of yacc.c  */
#line 4328 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval)) );
         ;}
    break;

  case 1023:

/* Line 1455 of yacc.c  */
#line 4333 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 1024:

/* Line 1455 of yacc.c  */
#line 4338 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 1025:

/* Line 1455 of yacc.c  */
#line 4343 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval)) );
         ;}
    break;

  case 1026:

/* Line 1455 of yacc.c  */
#line 4348 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new ElaLinIsoMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)));
         ;}
    break;

  case 1027:

/* Line 1455 of yacc.c  */
#line 4353 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1, 
             new ElaLinIsoMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), 0, 0));
	 ;}
    break;

  case 1028:

/* Line 1455 of yacc.c  */
#line 4358 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1,
             new ElaLinIsoMat((yyvsp[(4) - (5)].fval), 0, 0, 0, 0));
         ;}
    break;

  case 1029:

/* Line 1455 of yacc.c  */
#line 4363 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new StVenantKirchhoffMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)));
         ;}
    break;

  case 1030:

/* Line 1455 of yacc.c  */
#line 4368 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new StVenantKirchhoffMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), 0, 0));
         ;}
    break;

  case 1031:

/* Line 1455 of yacc.c  */
#line 4373 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1,
             new StVenantKirchhoffMat((yyvsp[(4) - (5)].fval), 0, 0, 0, 0));
         ;}
    break;

  case 1032:

/* Line 1455 of yacc.c  */
#line 4378 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new HenckyMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)));
         ;}
    break;

  case 1033:

/* Line 1455 of yacc.c  */
#line 4383 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new HenckyMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), 0, 0));
         ;}
    break;

  case 1034:

/* Line 1455 of yacc.c  */
#line 4388 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1,
             new HenckyMat((yyvsp[(4) - (5)].fval), 0, 0, 0, 0));
         ;}
    break;

  case 1035:

/* Line 1455 of yacc.c  */
#line 4393 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ElaLinIsoMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), 0, 0));
         ;}
    break;

  case 1036:

/* Line 1455 of yacc.c  */
#line 4398 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new ElaLinIsoMat2D((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)));
         ;}
    break;

  case 1037:

/* Line 1455 of yacc.c  */
#line 4403 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), 0, 0));
         ;}
    break;

  case 1038:

/* Line 1455 of yacc.c  */
#line 4408 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)));
         ;}
    break;

  case 1039:

/* Line 1455 of yacc.c  */
#line 4413 "p.y"
    {
            double params[3] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          ;}
    break;

  case 1040:

/* Line 1455 of yacc.c  */
#line 4419 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          ;}
    break;

  case 1041:

/* Line 1455 of yacc.c  */
#line 4425 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          ;}
    break;

  case 1042:

/* Line 1455 of yacc.c  */
#line 4431 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          ;}
    break;

  case 1043:

/* Line 1455 of yacc.c  */
#line 4437 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          ;}
    break;

  case 1044:

/* Line 1455 of yacc.c  */
#line 4443 "p.y"
    {
            double params[6] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          ;}
    break;

  case 1045:

/* Line 1455 of yacc.c  */
#line 4449 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 1.0e-6, -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 1046:

/* Line 1455 of yacc.c  */
#line 4455 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 1047:

/* Line 1455 of yacc.c  */
#line 4461 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 1048:

/* Line 1455 of yacc.c  */
#line 4467 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (24)].ival)-1,
             new ExpMat((yyvsp[(3) - (24)].ival), (yyvsp[(4) - (24)].fval), (yyvsp[(5) - (24)].fval), (yyvsp[(6) - (24)].fval), (yyvsp[(7) - (24)].fval), (yyvsp[(8) - (24)].fval), (yyvsp[(9) - (24)].fval), (yyvsp[(10) - (24)].fval), (yyvsp[(11) - (24)].fval), (yyvsp[(12) - (24)].fval), (yyvsp[(13) - (24)].fval), (yyvsp[(14) - (24)].fval), (yyvsp[(15) - (24)].fval), (yyvsp[(16) - (24)].fval), (yyvsp[(17) - (24)].fval), (yyvsp[(18) - (24)].fval), (yyvsp[(19) - (24)].fval), (yyvsp[(20) - (24)].fval), (yyvsp[(21) - (24)].fval), (yyvsp[(22) - (24)].fval), (yyvsp[(23) - (24)].fval)));
         ;}
    break;

  case 1049:

/* Line 1455 of yacc.c  */
#line 4472 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new ExpMat((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1050:

/* Line 1455 of yacc.c  */
#line 4477 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ExpMat((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1051:

/* Line 1455 of yacc.c  */
#line 4482 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new ExpMat((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1052:

/* Line 1455 of yacc.c  */
#line 4487 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new ExpMat((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1053:

/* Line 1455 of yacc.c  */
#line 4492 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
             new ExpMat((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1054:

/* Line 1455 of yacc.c  */
#line 4497 "p.y"
    {
	   geoSource->loadMaterial((yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].strval));
	 ;}
    break;

  case 1055:

/* Line 1455 of yacc.c  */
#line 4501 "p.y"
    {
	   geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].dlist));
	 ;}
    break;

  case 1057:

/* Line 1455 of yacc.c  */
#line 4508 "p.y"
    { geoSource->setMatUsage((yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 1058:

/* Line 1455 of yacc.c  */
#line 4510 "p.y"
    {
            for(int i = (yyvsp[(2) - (5)].ival)-1; i < (yyvsp[(3) - (5)].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[(4) - (5)].ival)-1);
	  ;}
    break;

  case 1059:

/* Line 1455 of yacc.c  */
#line 4516 "p.y"
    { (yyval.dlist).nval = 0; ;}
    break;

  case 1060:

/* Line 1455 of yacc.c  */
#line 4518 "p.y"
    { 
          if((yyvsp[(1) - (2)].dlist).nval == 32) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[(1) - (2)].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[(2) - (2)].fval);
 	;}
    break;

  case 1061:

/* Line 1455 of yacc.c  */
#line 4528 "p.y"
    { (yyval.slist).nval = 0; ;}
    break;

  case 1062:

/* Line 1455 of yacc.c  */
#line 4530 "p.y"
    { 
          if((yyvsp[(1) - (2)].slist).nval == 32) {
             fprintf(stderr, "Too many files!\n");
	     exit(-1);
          }
          (yyval.slist) = (yyvsp[(1) - (2)].slist);
          (yyval.slist).v[(yyval.slist).nval++] = (yyvsp[(2) - (2)].strval);
 	;}
    break;

  case 1063:

/* Line 1455 of yacc.c  */
#line 4541 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (4)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 1064:

/* Line 1455 of yacc.c  */
#line 4545 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (6)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (6)].ival)); ;}
    break;

  case 1065:

/* Line 1455 of yacc.c  */
#line 4548 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (8)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (8)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(7) - (8)].ival)); ;}
    break;

  case 1066:

/* Line 1455 of yacc.c  */
#line 4555 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().svdPodRom = true;;}
    break;

  case 1068:

/* Line 1455 of yacc.c  */
#line 4563 "p.y"
    { for(int i=0; i<(yyvsp[(2) - (2)].slist).nval; ++i) domain->solInfo().snapfiPodRom.push_back(std::string((yyvsp[(2) - (2)].slist).v[i])); ;}
    break;

  case 1069:

/* Line 1455 of yacc.c  */
#line 4574 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1070:

/* Line 1455 of yacc.c  */
#line 4576 "p.y"
    { domain->solInfo().normalize = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1071:

/* Line 1455 of yacc.c  */
#line 4578 "p.y"
    { domain->solInfo().normalize = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1072:

/* Line 1455 of yacc.c  */
#line 4580 "p.y"
    { for(int i=0; i<(yyvsp[(2) - (2)].dlist).nval; ++i) domain->solInfo().snapshotWeights.push_back((yyvsp[(2) - (2)].dlist).v[i]); ;}
    break;

  case 1073:

/* Line 1455 of yacc.c  */
#line 4582 "p.y"
    { if((yyvsp[(2) - (3)].ival)) for(int i=0; i<(yyvsp[(3) - (3)].dlist).nval; ++i) domain->solInfo().snapshotWeights.push_back((yyvsp[(3) - (3)].dlist).v[i]); ;}
    break;

  case 1074:

/* Line 1455 of yacc.c  */
#line 4584 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1075:

/* Line 1455 of yacc.c  */
#line 4586 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (3)].ival);
    domain->solInfo().skipOffSet = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 1076:

/* Line 1455 of yacc.c  */
#line 4589 "p.y"
    { for(int i=0; i<(yyvsp[(2) - (2)].slist).nval; ++i) domain->solInfo().robfi.push_back(std::string((yyvsp[(2) - (2)].slist).v[i])); ;}
    break;

  case 1077:

/* Line 1455 of yacc.c  */
#line 4591 "p.y"
    { domain->solInfo().svdBlockSize = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1079:

/* Line 1455 of yacc.c  */
#line 4597 "p.y"
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().DEIMBasisPod = true; ;}
    break;

  case 1081:

/* Line 1455 of yacc.c  */
#line 4605 "p.y"
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().UDEIMBasisPod = true; ;}
    break;

  case 1083:

/* Line 1455 of yacc.c  */
#line 4613 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().samplingPodRom = true; ;}
    break;

  case 1086:

/* Line 1455 of yacc.c  */
#line 4622 "p.y"
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().snapProjPodRom = true; ;}
    break;

  case 1088:

/* Line 1455 of yacc.c  */
#line 4630 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 1089:

/* Line 1455 of yacc.c  */
#line 4632 "p.y"
    { domain->solInfo().statePodRomFile.push_back((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 1090:

/* Line 1455 of yacc.c  */
#line 4634 "p.y"
    { domain->solInfo().statePodRomFile.push_back((yyvsp[(2) - (3)].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[(3) - (3)].strval)); ;}
    break;

  case 1091:

/* Line 1455 of yacc.c  */
#line 4637 "p.y"
    { domain->solInfo().statePodRomFile.push_back((yyvsp[(2) - (4)].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[(3) - (4)].strval));
    domain->solInfo().accelPodRomFile.push_back((yyvsp[(4) - (4)].strval)); ;}
    break;

  case 1092:

/* Line 1455 of yacc.c  */
#line 4641 "p.y"
    { domain->solInfo().tolPodRom = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 1093:

/* Line 1455 of yacc.c  */
#line 4643 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1094:

/* Line 1455 of yacc.c  */
#line 4645 "p.y"
    { domain->solInfo().skipOffSet = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1095:

/* Line 1455 of yacc.c  */
#line 4647 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1096:

/* Line 1455 of yacc.c  */
#line 4649 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (3)].ival); 
    domain->solInfo().forcePodSize = (yyvsp[(3) - (3)].ival);;}
    break;

  case 1097:

/* Line 1455 of yacc.c  */
#line 4652 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (4)].ival); 
    domain->solInfo().forcePodSize = (yyvsp[(3) - (4)].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[(4) - (4)].ival); ;}
    break;

  case 1098:

/* Line 1455 of yacc.c  */
#line 4656 "p.y"
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1099:

/* Line 1455 of yacc.c  */
#line 4658 "p.y"
    { domain->solInfo().useMassOrthogonalProjection = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1100:

/* Line 1455 of yacc.c  */
#line 4660 "p.y"
    { domain->solInfo().reduceFollower = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1101:

/* Line 1455 of yacc.c  */
#line 4662 "p.y"
    { domain->solInfo().reduceFollower = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1102:

/* Line 1455 of yacc.c  */
#line 4664 "p.y"
    { domain->solInfo().PODerrornorm.push_back((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 1103:

/* Line 1455 of yacc.c  */
#line 4666 "p.y"
    { domain->solInfo().PODerrornorm.push_back((yyvsp[(2) - (3)].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[(3) - (3)].strval)); ;}
    break;

  case 1104:

/* Line 1455 of yacc.c  */
#line 4669 "p.y"
    { domain->solInfo().PODerrornorm.push_back((yyvsp[(2) - (4)].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[(3) - (4)].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[(4) - (4)].strval)); ;}
    break;

  case 1105:

/* Line 1455 of yacc.c  */
#line 4673 "p.y"
    { domain->solInfo().useScalingSpnnls = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1106:

/* Line 1455 of yacc.c  */
#line 4675 "p.y"
    { domain->solInfo().solverTypeSpnnls = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1107:

/* Line 1455 of yacc.c  */
#line 4677 "p.y"
    { domain->solInfo().maxSizeSpnnls = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 1108:

/* Line 1455 of yacc.c  */
#line 4679 "p.y"
    { domain->solInfo().forcePodRomFile = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 1109:

/* Line 1455 of yacc.c  */
#line 4681 "p.y"
    { domain->solInfo().forcePodRomFile = (yyvsp[(2) - (3)].strval);
    domain->solInfo().forcePodSize = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 1110:

/* Line 1455 of yacc.c  */
#line 4684 "p.y"
    { domain->solInfo().forcePodRomFile = (yyvsp[(2) - (4)].strval); 
    domain->solInfo().forcePodSize = (yyvsp[(3) - (4)].ival); 
    domain->solInfo().maxDeimBasisSize = (yyvsp[(4) - (4)].ival); ;}
    break;

  case 1111:

/* Line 1455 of yacc.c  */
#line 4688 "p.y"
    { domain->solInfo().selectFullNode = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1112:

/* Line 1455 of yacc.c  */
#line 4690 "p.y"
    { domain->solInfo().selectFullElem = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1113:

/* Line 1455 of yacc.c  */
#line 4692 "p.y"
    { domain->solInfo().computeForceSnap = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1114:

/* Line 1455 of yacc.c  */
#line 4694 "p.y"
    { domain->solInfo().orthogForceSnap = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1116:

/* Line 1455 of yacc.c  */
#line 4700 "p.y"
    { domain->solInfo().conwepConfigurations.push_back((yyvsp[(2) - (3)].blastData)); ;}
    break;

  case 1117:

/* Line 1455 of yacc.c  */
#line 4705 "p.y"
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().ROMPostProcess = true; ;}
    break;

  case 1119:

/* Line 1455 of yacc.c  */
#line 4713 "p.y"
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[(2) - (2)].strval)); 
    domain->solInfo().numRODFile += 1; ;}
    break;

  case 1120:

/* Line 1455 of yacc.c  */
#line 4716 "p.y"
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[(2) - (3)].strval));
    domain->solInfo().numRODFile += 1; 
    domain->solInfo().skipPodRom = (yyvsp[(3) - (3)].ival);;}
    break;

  case 1121:

/* Line 1455 of yacc.c  */
#line 4723 "p.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 1122:

/* Line 1455 of yacc.c  */
#line 4725 "p.y"
    { (yyval.ival) = std::numeric_limits<int>::max(); ;}
    break;

  case 1123:

/* Line 1455 of yacc.c  */
#line 4730 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 1124:

/* Line 1455 of yacc.c  */
#line 4732 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].fval); ;}
    break;

  case 1125:

/* Line 1455 of yacc.c  */
#line 4734 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::infinity(); ;}
    break;

  case 1126:

/* Line 1455 of yacc.c  */
#line 4736 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); ;}
    break;



/* Line 1455 of yacc.c  */
#line 13504 "/home/tac688/codes/AEROS/Parser.d/parser.cpp"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
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

  /* Do not reclaim the symbols of the rule which action triggered
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
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
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

  *++yyvsp = yylval;


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

#if !defined(yyoverflow) || YYERROR_VERBOSE
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
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
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
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 4738 "p.y"


