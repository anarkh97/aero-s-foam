/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2011 Free Software Foundation, Inc.
   
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
#define YYBISON_VERSION "2.5"

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

/* Line 268 of yacc.c  */
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
#ifdef STRUCTOPT
#include <Structopt.d/Driver_opt.d/Domain_opt.h>
#endif

 int numColumns = 3;
 double amplitude = 1.0;
 int PitaTS = 1;         //CD: Pita


/* Line 268 of yacc.c  */
#line 93 "/home/tac688/Research/FEM1/FEM3/Parser.d/parser.cpp"

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
     ANALYSIS = 262,
     ARCLENGTH = 263,
     ATTRIBUTES = 264,
     AUGMENT = 265,
     AUGMENTTYPE = 266,
     AVERAGED = 267,
     ATDARB = 268,
     ACOU = 269,
     ATDDNB = 270,
     ATDROB = 271,
     ARPACK = 272,
     ATDDIR = 273,
     ATDNEU = 274,
     AXIHDIR = 275,
     AXIHNEU = 276,
     AXINUMMODES = 277,
     AXINUMSLICES = 278,
     AXIHSOMMER = 279,
     AXIMPC = 280,
     AUXCOARSESOLVER = 281,
     ACMECNTL = 282,
     ADDEDMASS = 283,
     AEROEMBED = 284,
     BLOCKDIAG = 285,
     BOFFSET = 286,
     BUCKLE = 287,
     BGTL = 288,
     BMPC = 289,
     BINARYINPUT = 290,
     BINARYOUTPUT = 291,
     CHECKTOKEN = 292,
     COARSESOLVER = 293,
     COEF = 294,
     CFRAMES = 295,
     COLLOCATEDTYPE = 296,
     CONVECTION = 297,
     COMPOSITE = 298,
     CONDITION = 299,
     CONTROL = 300,
     CORNER = 301,
     CORNERTYPE = 302,
     CURVE = 303,
     CCTTOL = 304,
     CCTSOLVER = 305,
     CRHS = 306,
     COUPLEDSCALE = 307,
     CONTACTSURFACES = 308,
     CMPC = 309,
     CNORM = 310,
     COMPLEXOUTTYPE = 311,
     CONSTRMAT = 312,
     CASES = 313,
     CONSTRAINEDSURFACES = 314,
     CSFRAMES = 315,
     CSTYPE = 316,
     DAMPING = 317,
     DblConstant = 318,
     DEM = 319,
     DIMASS = 320,
     DISP = 321,
     DIRECT = 322,
     DLAMBDA = 323,
     DP = 324,
     DYNAM = 325,
     DETER = 326,
     DECOMPOSE = 327,
     DECOMPFILE = 328,
     DMPC = 329,
     DEBUGCNTL = 330,
     DEBUGICNTL = 331,
     CONSTRAINTS = 332,
     MULTIPLIERS = 333,
     PENALTY = 334,
     EIGEN = 335,
     EFRAMES = 336,
     ELSCATTERER = 337,
     END = 338,
     ELHSOMMERFELD = 339,
     EXPLICIT = 340,
     EPSILON = 341,
     ELEMENTARYFUNCTIONTYPE = 342,
     FABMAT = 343,
     FACOUSTICS = 344,
     FETI = 345,
     FETI2TYPE = 346,
     FETIPREC = 347,
     FFP = 348,
     FFPDIR = 349,
     FITALG = 350,
     FLUMAT = 351,
     FNAME = 352,
     FLUX = 353,
     FORCE = 354,
     FRONTAL = 355,
     FETIH = 356,
     FILTEREIG = 357,
     FREQSWEEP = 358,
     FREQSWEEP1 = 359,
     FREQSWEEP2 = 360,
     FSINTERFACE = 361,
     FSISCALING = 362,
     FSIELEMENT = 363,
     NOLOCALFSISPLITING = 364,
     FSICORNER = 365,
     FFIDEBUG = 366,
     FAILSAFE = 367,
     GEPS = 368,
     GLOBALTOL = 369,
     GRAVITY = 370,
     GRBM = 371,
     GTGSOLVER = 372,
     GLOBALCRBMTOL = 373,
     GROUP = 374,
     GROUPTYPE = 375,
     GOLDFARBTOL = 376,
     GOLDFARBCHECK = 377,
     HDIRICHLET = 378,
     HEAT = 379,
     HFETI = 380,
     HNEUMAN = 381,
     HSOMMERFELD = 382,
     HFTT = 383,
     HELMHOLTZ = 384,
     HNBO = 385,
     HELMMF = 386,
     HELMSO = 387,
     HSCBO = 388,
     HWIBO = 389,
     HZEM = 390,
     HZEMFILTER = 391,
     HLMPC = 392,
     HELMSWEEP = 393,
     HELMSWEEP1 = 394,
     HELMSWEEP2 = 395,
     HERMITIAN = 396,
     HESSIAN = 397,
     IACC = 398,
     IDENTITY = 399,
     IDIS = 400,
     IDIS6 = 401,
     IntConstant = 402,
     INTERFACELUMPED = 403,
     ITEMP = 404,
     ITERTYPE = 405,
     IVEL = 406,
     INCIDENCE = 407,
     IHDIRICHLET = 408,
     IHDSWEEP = 409,
     IHNEUMANN = 410,
     ISOLVERTYPE = 411,
     INPC = 412,
     INFINTY = 413,
     JACOBI = 414,
     KRYLOVTYPE = 415,
     KIRLOC = 416,
     LAYC = 417,
     LAYN = 418,
     LAYD = 419,
     LAYO = 420,
     LAYMAT = 421,
     LFACTOR = 422,
     LMPC = 423,
     LOAD = 424,
     LOBPCG = 425,
     LOCALSOLVER = 426,
     LINESEARCH = 427,
     LUMPED = 428,
     MASS = 429,
     MATERIALS = 430,
     MATLAB = 431,
     MAXITR = 432,
     MAXORTHO = 433,
     MAXVEC = 434,
     MODAL = 435,
     MPCPRECNO = 436,
     MPCPRECNOID = 437,
     MPCTYPE = 438,
     MPCTYPEID = 439,
     MPCSCALING = 440,
     MPCELEMENT = 441,
     MPCBLOCKID = 442,
     MPCBLK_OVERLAP = 443,
     MFTT = 444,
     MPTT = 445,
     MRHS = 446,
     MPCCHECK = 447,
     MUMPSICNTL = 448,
     MUMPSCNTL = 449,
     MECH = 450,
     MODEFILTER = 451,
     NDTYPE = 452,
     NEIGPA = 453,
     NEWMARK = 454,
     NewLine = 455,
     NL = 456,
     NLMAT = 457,
     NLPREC = 458,
     NOCOARSE = 459,
     NODETOKEN = 460,
     NONINPC = 461,
     NSBSPV = 462,
     NLTOL = 463,
     NUMCGM = 464,
     NOSECONDARY = 465,
     OPTIMIZATION = 466,
     OUTPUT = 467,
     OUTPUT6 = 468,
     QSTATIC = 469,
     QLOAD = 470,
     PITA = 471,
     PITADISP6 = 472,
     PITAVEL6 = 473,
     NOFORCE = 474,
     MDPITA = 475,
     GLOBALBASES = 476,
     LOCALBASES = 477,
     TIMEREVERSIBLE = 478,
     REMOTECOARSE = 479,
     ORTHOPROJTOL = 480,
     READINITSEED = 481,
     JUMPCVG = 482,
     JUMPOUTPUT = 483,
     PRECNO = 484,
     PRECONDITIONER = 485,
     PRELOAD = 486,
     PRESSURE = 487,
     PRINTMATLAB = 488,
     PROJ = 489,
     PIVOT = 490,
     PRECTYPE = 491,
     PRECTYPEID = 492,
     PICKANYCORNER = 493,
     PADEPIVOT = 494,
     PROPORTIONING = 495,
     PLOAD = 496,
     PADEPOLES = 497,
     POINTSOURCE = 498,
     PLANEWAVE = 499,
     PTOL = 500,
     PLANTOL = 501,
     PMAXIT = 502,
     RADIATION = 503,
     RBMFILTER = 504,
     RBMSET = 505,
     READMODE = 506,
     REBUILD = 507,
     RENUM = 508,
     RENUMBERID = 509,
     REORTHO = 510,
     RESTART = 511,
     RECONS = 512,
     RECONSALG = 513,
     REBUILDCCT = 514,
     RANDOM = 515,
     RPROP = 516,
     RNORM = 517,
     REVERSENORMALS = 518,
     RIGID = 519,
     SCALING = 520,
     SCALINGTYPE = 521,
     SENSORS = 522,
     SOLVERTYPE = 523,
     SHIFT = 524,
     SPOOLESTAU = 525,
     SPOOLESSEED = 526,
     SPOOLESMAXSIZE = 527,
     SPOOLESMAXDOMAINSIZE = 528,
     SPOOLESMAXZEROS = 529,
     SPOOLESMSGLVL = 530,
     SPOOLESSCALE = 531,
     SPOOLESPIVOT = 532,
     SPOOLESRENUM = 533,
     SPARSEMAXSUP = 534,
     SPARSEDEFBLK = 535,
     STATS = 536,
     STRESSID = 537,
     SUBSPACE = 538,
     SURFACE = 539,
     SAVEMEMCOARSE = 540,
     SPACEDIMENSION = 541,
     SCATTERER = 542,
     STAGTOL = 543,
     SCALED = 544,
     SWITCH = 545,
     STABLE = 546,
     SUBTYPE = 547,
     STEP = 548,
     SOWER = 549,
     SHELLTHICKNESS = 550,
     SURF = 551,
     SPRINGMAT = 552,
     TANGENT = 553,
     TEMP = 554,
     TIME = 555,
     TOLEIG = 556,
     TOLFETI = 557,
     TOLJAC = 558,
     TOLPCG = 559,
     TOPFILE = 560,
     TOPOLOGY = 561,
     TRBM = 562,
     THERMOE = 563,
     THERMOH = 564,
     TETT = 565,
     TOLCGM = 566,
     TURKEL = 567,
     TIEDSURFACES = 568,
     THETA = 569,
     HRC = 570,
     THIRDNODE = 571,
     THERMMAT = 572,
     TDENFORC = 573,
     TESTULRICH = 574,
     THRU = 575,
     TOPFLAG = 576,
     USE = 577,
     USERDEFINEDISP = 578,
     USERDEFINEFORCE = 579,
     UPROJ = 580,
     UNSYMMETRIC = 581,
     USING = 582,
     VERSION = 583,
     WAVENUMBER = 584,
     WETCORNERS = 585,
     XPOST = 586,
     YMTT = 587,
     ZERO = 588,
     BINARY = 589,
     GEOMETRY = 590,
     DECOMPOSITION = 591,
     GLOBAL = 592,
     MATCHER = 593,
     CPUMAP = 594,
     NODALCONTACT = 595,
     MODE = 596,
     FRIC = 597,
     GAP = 598,
     OUTERLOOP = 599,
     EDGEWS = 600,
     WAVETYPE = 601,
     ORTHOTOL = 602,
     IMPE = 603,
     FREQ = 604,
     DPH = 605,
     WAVEMETHOD = 606,
     MATSPEC = 607,
     MATUSAGE = 608,
     BILINEARPLASTIC = 609,
     FINITESTRAINPLASTIC = 610,
     LINEARELASTIC = 611,
     STVENANTKIRCHHOFF = 612,
     LINPLSTRESS = 613,
     READ = 614,
     OPTCTV = 615,
     ISOTROPICLINEARELASTIC = 616,
     NEOHOOKEAN = 617,
     ISOTROPICLINEARELASTICJ2PLASTIC = 618,
     ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 619,
     HYPERELASTIC = 620,
     MOONEYRIVLIN = 621,
     HENCKY = 622,
     LOGSTRAINPLASTIC = 623,
     SVKPLSTRESS = 624,
     SURFACETOPOLOGY = 625,
     MORTARTIED = 626,
     MORTARSCALING = 627,
     MORTARINTEGRATIONRULE = 628,
     SEARCHTOL = 629,
     STDMORTAR = 630,
     DUALMORTAR = 631,
     WETINTERFACE = 632,
     NSUBS = 633,
     EXITAFTERDEC = 634,
     SKIP = 635,
     OUTPUTMEMORY = 636,
     OUTPUTWEIGHT = 637,
     WEIGHTLIST = 638,
     GMRESRESIDUAL = 639,
     SLOSH = 640,
     SLGRAV = 641,
     SLZEM = 642,
     SLZEMFILTER = 643,
     PDIR = 644,
     HEFSB = 645,
     HEFRS = 646,
     HEINTERFACE = 647,
     SNAPFI = 648,
     PODROB = 649,
     TRNVCT = 650,
     ORTHOG = 651,
     SVDTOKEN = 652,
     SAMPLING = 653,
     PODSIZEMAX = 654,
     REFSUBSTRACT = 655,
     TOLER = 656
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 293 of yacc.c  */
#line 23 "p.y"

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



/* Line 293 of yacc.c  */
#line 565 "/home/tac688/Research/FEM1/FEM3/Parser.d/parser.cpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 343 of yacc.c  */
#line 577 "/home/tac688/Research/FEM1/FEM3/Parser.d/parser.cpp"

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
# if defined YYENABLE_NLS && YYENABLE_NLS
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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
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
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
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

# define YYCOPY_NEEDED 1

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

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
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
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  640
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   5232

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  402
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  220
/* YYNRULES -- Number of rules.  */
#define YYNRULES  930
/* YYNRULES -- Number of states.  */
#define YYNSTATES  2290

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   656

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
     395,   396,   397,   398,   399,   400,   401
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
     261,   263,   265,   267,   269,   275,   280,   283,   289,   296,
     303,   306,   313,   319,   327,   335,   344,   348,   351,   354,
     357,   360,   364,   367,   372,   376,   381,   386,   393,   396,
     400,   404,   407,   412,   417,   422,   427,   432,   436,   439,
     444,   449,   453,   457,   461,   465,   469,   472,   477,   480,
     485,   488,   493,   496,   501,   504,   507,   510,   513,   516,
     519,   522,   525,   530,   535,   541,   545,   548,   552,   555,
     559,   562,   566,   569,   581,   595,   610,   616,   619,   622,
     630,   640,   652,   665,   668,   674,   681,   684,   690,   693,
     699,   704,   709,   713,   717,   721,   725,   729,   733,   736,
     739,   743,   747,   751,   755,   761,   766,   772,   779,   787,
     791,   797,   802,   808,   815,   823,   827,   830,   834,   837,
     840,   844,   847,   851,   854,   858,   860,   863,   867,   871,
     874,   880,   884,   888,   892,   895,   899,   902,   906,   911,
     916,   922,   926,   930,   933,   938,   941,   945,   948,   952,
     955,   961,   964,   967,   970,   973,   976,   979,   983,   988,
     993,  1002,  1007,  1012,  1016,  1020,  1026,  1029,  1033,  1037,
    1041,  1044,  1049,  1051,  1055,  1057,  1060,  1063,  1066,  1070,
    1076,  1083,  1090,  1094,  1099,  1106,  1112,  1116,  1121,  1125,
    1127,  1130,  1137,  1140,  1143,  1146,  1149,  1152,  1157,  1163,
    1168,  1171,  1175,  1178,  1182,  1185,  1190,  1192,  1195,  1198,
    1201,  1207,  1213,  1220,  1228,  1232,  1234,  1236,  1239,  1241,
    1243,  1245,  1248,  1250,  1252,  1255,  1257,  1262,  1267,  1271,
    1278,  1284,  1291,  1295,  1298,  1304,  1311,  1314,  1317,  1325,
    1335,  1339,  1342,  1349,  1353,  1355,  1358,  1362,  1365,  1369,
    1371,  1374,  1379,  1382,  1387,  1393,  1396,  1401,  1407,  1410,
    1417,  1420,  1427,  1431,  1433,  1436,  1441,  1447,  1454,  1458,
    1461,  1466,  1468,  1471,  1475,  1477,  1480,  1485,  1491,  1498,
    1502,  1505,  1510,  1514,  1517,  1522,  1526,  1529,  1534,  1538,
    1541,  1546,  1550,  1552,  1555,  1559,  1563,  1567,  1570,  1574,
    1577,  1583,  1586,  1590,  1595,  1600,  1610,  1613,  1616,  1626,
    1629,  1633,  1638,  1642,  1646,  1654,  1657,  1661,  1666,  1669,
    1672,  1678,  1686,  1697,  1701,  1706,  1709,  1713,  1717,  1722,
    1726,  1729,  1739,  1746,  1749,  1752,  1756,  1766,  1770,  1780,
    1783,  1787,  1792,  1796,  1799,  1803,  1806,  1814,  1824,  1828,
    1830,  1833,  1840,  1848,  1857,  1867,  1869,  1872,  1874,  1877,
    1880,  1884,  1891,  1896,  1904,  1907,  1911,  1918,  1923,  1931,
    1934,  1938,  1941,  1944,  1948,  1951,  1955,  1961,  1966,  1971,
    1974,  1978,  1981,  1984,  1988,  1994,  1998,  2001,  2007,  2012,
    2016,  2021,  2023,  2026,  2030,  2033,  2050,  2070,  2090,  2111,
    2135,  2159,  2169,  2184,  2201,  2214,  2228,  2243,  2249,  2256,
    2273,  2286,  2292,  2300,  2309,  2319,  2331,  2342,  2354,  2367,
    2375,  2384,  2394,  2399,  2403,  2406,  2410,  2415,  2421,  2428,
    2434,  2439,  2445,  2451,  2458,  2466,  2474,  2483,  2491,  2500,
    2505,  2509,  2512,  2518,  2525,  2533,  2542,  2553,  2565,  2572,
    2582,  2585,  2591,  2599,  2603,  2606,  2611,  2614,  2620,  2628,
    2631,  2637,  2644,  2652,  2661,  2672,  2684,  2698,  2713,  2720,
    2730,  2734,  2738,  2742,  2746,  2750,  2753,  2759,  2764,  2768,
    2773,  2775,  2778,  2783,  2787,  2791,  2795,  2801,  2806,  2809,
    2812,  2815,  2818,  2830,  2835,  2838,  2846,  2849,  2854,  2861,
    2868,  2877,  2885,  2891,  2897,  2905,  2914,  2917,  2922,  2928,
    2932,  2935,  2939,  2942,  2947,  2954,  2961,  2970,  2972,  2974,
    2979,  2981,  2984,  2989,  2994,  2999,  3004,  3009,  3014,  3019,
    3025,  3030,  3036,  3042,  3048,  3055,  3063,  3072,  3082,  3089,
    3095,  3101,  3107,  3114,  3122,  3128,  3133,  3136,  3140,  3144,
    3148,  3153,  3158,  3162,  3166,  3170,  3174,  3178,  3182,  3186,
    3190,  3194,  3198,  3202,  3207,  3212,  3216,  3220,  3225,  3230,
    3235,  3240,  3245,  3250,  3254,  3258,  3262,  3267,  3271,  3276,
    3281,  3286,  3291,  3296,  3299,  3303,  3307,  3311,  3315,  3319,
    3323,  3327,  3330,  3333,  3337,  3340,  3344,  3347,  3351,  3356,
    3360,  3364,  3369,  3374,  3380,  3386,  3393,  3397,  3402,  3406,
    3410,  3414,  3418,  3422,  3425,  3429,  3433,  3437,  3441,  3446,
    3451,  3456,  3461,  3466,  3471,  3476,  3480,  3483,  3486,  3490,
    3494,  3497,  3501,  3506,  3512,  3519,  3527,  3530,  3534,  3539,
    3543,  3547,  3551,  3555,  3558,  3561,  3565,  3570,  3576,  3580,
    3585,  3589,  3593,  3597,  3601,  3607,  3613,  3618,  3623,  3630,
    3637,  3641,  3645,  3650,  3654,  3658,  3662,  3666,  3670,  3674,
    3677,  3681,  3686,  3688,  3691,  3695,  3700,  3706,  3708,  3711,
    3715,  3719,  3724,  3727,  3730,  3734,  3740,  3744,  3749,  3755,
    3760,  3765,  3770,  3778,  3786,  3795,  3799,  3809,  3819,  3830,
    3833,  3836,  3839,  3843,  3849,  3854,  3856,  3859,  3864,  3871,
    3877,  3879,  3882,  3887,  3891,  3894,  3898,  3901,  3904,  3907,
    3911,  3915,  3920,  3925,  3931,  3939,  3945,  3953,  3958,  3964,
    3968,  3973,  3977,  3982,  3985,  3989,  3994,  4000,  4009,  4012,
    4015,  4026,  4039,  4043,  4048,  4051,  4056,  4064,  4074,  4084,
    4093,  4105,  4118,  4121,  4131,  4142,  4152,  4163,  4173,  4184,
    4192,  4200,  4208,  4217,  4226,  4234,  4242,  4251,  4260,  4270,
    4281,  4292,  4304,  4317,  4342,  4348,  4354,  4357,  4362,  4368,
    4369,  4372,  4377,  4384,  4393,  4396,  4400,  4403,  4407,  4410,
    4413,  4417,  4420,  4423,  4426,  4429,  4432,  4434,  4436,  4438,
    4440
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     403,     0,    -1,   404,    83,    -1,   405,    -1,   404,   405,
      -1,   570,    -1,   484,    -1,   539,    -1,   547,    -1,   551,
      -1,   559,    -1,   579,    -1,   578,    -1,   582,    -1,   557,
      -1,   586,    -1,   583,    -1,   584,    -1,   585,    -1,   615,
      -1,   533,    -1,   532,    -1,   535,    -1,   536,    -1,   534,
      -1,   537,    -1,   538,    -1,   440,    -1,   441,    -1,   443,
      -1,   442,    -1,   447,    -1,   448,    -1,   565,    -1,   456,
      -1,   444,    -1,   446,    -1,   606,    -1,   435,    -1,   423,
      -1,   424,    -1,   433,    -1,   420,    -1,   421,    -1,   422,
      -1,   543,    -1,   545,    -1,   467,    -1,   468,    -1,   436,
      -1,   470,    -1,   472,    -1,   473,    -1,   469,    -1,   437,
      -1,   438,    -1,   439,    -1,   609,    -1,   610,    -1,   417,
      -1,   608,    -1,   459,    -1,   461,    -1,   462,    -1,   463,
      -1,   465,    -1,   466,    -1,   450,    -1,   449,    -1,   492,
      -1,   493,    -1,   494,    -1,   486,    -1,   489,    -1,   495,
      -1,   451,    -1,   590,    -1,   605,    -1,   594,    -1,   601,
      -1,   603,    -1,   502,    -1,   520,    -1,   521,    -1,   598,
      -1,   496,    -1,   499,    -1,   505,    -1,   507,    -1,   509,
      -1,   511,    -1,   555,    -1,   483,    -1,   479,    -1,   480,
      -1,   481,    -1,   522,    -1,   524,    -1,   526,    -1,   527,
      -1,   528,    -1,   530,    -1,   415,    -1,   416,    -1,   515,
      -1,   516,    -1,   517,    -1,   518,    -1,   519,    -1,   418,
      -1,   419,    -1,   611,    -1,   464,    -1,   406,    -1,   407,
      -1,   408,    -1,   409,    -1,   410,    -1,   612,    -1,   613,
      -1,   560,    -1,   561,    -1,   563,    -1,   562,    -1,   564,
      -1,   567,    -1,   568,    -1,   485,    -1,   581,    -1,   475,
      -1,   569,    -1,   592,    -1,   616,    -1,   618,    -1,   206,
     200,   620,   620,   200,    -1,   157,   200,   620,   200,    -1,
     119,   200,    -1,   408,   120,   620,   620,   200,    -1,   408,
     120,   620,   620,   620,   200,    -1,   408,   120,   296,   620,
     620,   200,    -1,   260,   200,    -1,   409,   620,   261,   621,
     621,   200,    -1,   348,   200,   349,   621,   200,    -1,   348,
     200,   104,   621,   621,   620,   200,    -1,   348,   200,   105,
     621,   621,   620,   200,    -1,   348,   200,   103,   621,   621,
     620,   620,   200,    -1,   348,   200,   413,    -1,   410,   414,
      -1,   410,   478,    -1,   410,   411,    -1,   410,   412,    -1,
     239,   621,   200,    -1,   242,   200,    -1,   242,   621,   621,
     200,    -1,   103,   621,   200,    -1,   413,   621,   620,   200,
      -1,   257,   258,   620,   200,    -1,   257,   258,   620,   620,
     620,   200,    -1,   294,   200,    -1,    35,   290,   200,    -1,
      36,   290,   200,    -1,   334,   200,    -1,   416,   335,    97,
     200,    -1,   416,   336,    97,   200,    -1,   416,   337,    97,
     200,    -1,   416,   338,    97,   200,    -1,   416,   339,    97,
     200,    -1,     7,   620,   200,    -1,    72,   200,    -1,   418,
      73,    97,   200,    -1,   418,   378,   620,   200,    -1,   418,
     382,   200,    -1,   418,   381,   200,    -1,   418,   379,   200,
      -1,   418,   380,   200,    -1,   418,    71,   200,    -1,   383,
     200,    -1,   419,   620,   621,   200,    -1,   189,   200,    -1,
     420,   621,   621,   200,    -1,   190,   200,    -1,   421,   621,
     621,   200,    -1,   128,   200,    -1,   422,   621,   621,   200,
      -1,    43,   200,    -1,   423,   425,    -1,   423,   427,    -1,
     423,   428,    -1,   423,   429,    -1,   423,   430,    -1,    40,
     200,    -1,   424,   580,    -1,    39,   620,   200,   426,    -1,
     620,   620,   621,   200,    -1,   426,   620,   620,   621,   200,
      -1,   162,   620,   200,    -1,   427,   431,    -1,   163,   620,
     200,    -1,   428,   431,    -1,   164,   620,   200,    -1,   429,
     432,    -1,   165,   620,   200,    -1,   430,   432,    -1,   620,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   200,
      -1,   620,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   200,    -1,   620,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   200,    -1,
     620,   620,   621,   621,   200,    -1,   166,   200,    -1,   433,
     434,    -1,   620,   621,   621,   621,   621,   621,   200,    -1,
     620,   621,   621,   621,   621,   621,   621,   621,   200,    -1,
     620,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     200,    -1,   620,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   200,    -1,    65,   200,    -1,   435,   620,
     620,   621,   200,    -1,   435,   620,   620,   620,   621,   200,
      -1,   115,   200,    -1,   436,   621,   621,   621,   200,    -1,
     256,   200,    -1,   437,    97,    97,    97,   200,    -1,   437,
      97,    97,   200,    -1,   437,    97,   620,   200,    -1,   169,
      97,   200,    -1,   322,    97,   200,    -1,   267,   200,   540,
      -1,     3,   200,   540,    -1,   324,   200,   540,    -1,   323,
     200,   540,    -1,   212,   200,    -1,   213,   200,    -1,   212,
     620,   200,    -1,   213,   620,   200,    -1,   444,   445,   200,
      -1,   282,    97,   620,    -1,   282,   620,   620,    97,   620,
      -1,   282,    97,   620,   620,    -1,   282,    97,   620,   120,
     620,    -1,   282,   620,   620,    97,   620,   620,    -1,   282,
     620,   620,    97,   620,   120,   620,    -1,   318,    97,   620,
      -1,   318,   620,   620,    97,   620,    -1,   318,    97,   620,
     620,    -1,   318,    97,   620,   120,   620,    -1,   318,   620,
     620,    97,   620,   620,    -1,   318,   620,   620,    97,   620,
     120,   620,    -1,   445,   205,   620,    -1,   445,   284,    -1,
     445,   621,   621,    -1,   445,    12,    -1,   445,    56,    -1,
     445,    56,   620,    -1,   445,   197,    -1,   445,   197,   620,
      -1,   445,   176,    -1,   331,   321,   200,    -1,   452,    -1,
      80,   200,    -1,    80,   620,   200,    -1,   198,   620,   200,
      -1,   283,   200,    -1,   283,   620,   621,   621,   200,    -1,
     207,   620,   200,    -1,   301,   621,   200,    -1,   303,   621,
     200,    -1,    85,   200,    -1,   269,   621,   200,    -1,    17,
     200,    -1,    17,    97,   200,    -1,    17,    97,   620,   200,
      -1,    17,   621,   621,   200,    -1,    17,   621,   621,   620,
     200,    -1,   102,   290,   200,    -1,   159,   620,   200,    -1,
     170,   200,    -1,   447,   177,   620,   200,    -1,   319,   200,
      -1,    28,   620,   200,    -1,   385,   200,    -1,   386,   621,
     200,    -1,   174,   200,    -1,    44,   200,   621,   620,   200,
      -1,    44,   200,    -1,   305,   200,    -1,    70,   200,    -1,
     452,   453,    -1,   452,   474,    -1,   452,   478,    -1,   452,
     180,   200,    -1,   452,   291,   620,   200,    -1,   452,   291,
     290,   200,    -1,   452,   291,   620,   621,   621,   620,   620,
     200,    -1,   452,   143,   290,   200,    -1,   452,   333,   290,
     200,    -1,   452,   210,   200,    -1,   452,    37,   200,    -1,
     452,    37,   621,   621,   200,    -1,   199,   200,    -1,   195,
     454,   200,    -1,    14,   454,   200,    -1,   124,   455,   200,
      -1,   621,   621,    -1,   621,   621,   621,   621,    -1,   621,
      -1,   621,   621,   621,    -1,   621,    -1,   214,   200,    -1,
     456,   457,    -1,   456,   458,    -1,   456,   180,   200,    -1,
     195,   621,   621,   620,   200,    -1,   195,   621,   621,   620,
     621,   200,    -1,   124,   621,   621,   621,   620,   200,    -1,
       4,     6,   200,    -1,     4,   200,     6,   200,    -1,     4,
     200,     6,   621,   621,   200,    -1,     4,   200,     6,   621,
     200,    -1,   459,    41,   200,    -1,   459,   338,    97,   200,
      -1,   459,   460,   200,    -1,    29,    -1,   460,   620,    -1,
       5,   200,     6,   621,   621,   200,    -1,   309,   200,    -1,
     308,   200,    -1,   341,   200,    -1,   135,   200,    -1,   387,
     200,    -1,   307,   200,   621,   200,    -1,   116,   200,   621,
     621,   200,    -1,   116,   200,   621,   200,    -1,   116,   200,
      -1,   196,   620,   200,    -1,   196,   200,    -1,   249,   620,
     200,    -1,   249,   200,    -1,   249,   200,   471,   200,    -1,
     620,    -1,   471,   620,    -1,   136,   200,    -1,   388,   200,
      -1,   300,   621,   621,   621,   200,    -1,   216,   200,   620,
     620,   476,    -1,   216,   200,   620,   620,   620,   476,    -1,
     220,   200,   620,   620,   620,   620,   476,    -1,   200,   477,
     476,    -1,   200,    -1,   219,    -1,   221,   620,    -1,   222,
      -1,   223,    -1,   224,    -1,   225,   621,    -1,   226,    -1,
     227,    -1,   227,   621,    -1,   228,    -1,    62,   621,   621,
     200,    -1,    62,   180,   200,   541,    -1,   123,   200,   556,
      -1,   153,   200,   621,   621,   621,   200,    -1,   153,   200,
     621,   621,   200,    -1,   154,   200,   620,   621,   621,   200,
      -1,   154,   620,   200,    -1,   481,   482,    -1,   620,   621,
     621,   621,   200,    -1,   155,   200,   621,   621,   621,   200,
      -1,    66,   200,    -1,   484,   574,    -1,   484,   620,   320,
     620,   620,   621,   200,    -1,   484,   620,   320,   620,   293,
     620,   620,   621,   200,    -1,   484,   296,   574,    -1,    59,
     200,    -1,   485,   620,    61,   620,   620,   200,    -1,   389,
     200,   487,    -1,   488,    -1,   487,   488,    -1,   620,   621,
     200,    -1,   620,   200,    -1,   391,   200,   490,    -1,   491,
      -1,   490,   491,    -1,   620,   620,   501,   200,    -1,   299,
     200,    -1,   492,   620,   621,   200,    -1,   492,   296,   620,
     621,   200,    -1,    98,   200,    -1,   493,   620,   621,   200,
      -1,   493,   296,   620,   621,   200,    -1,    42,   200,    -1,
     494,   620,   621,   621,   621,   200,    -1,   248,   200,    -1,
     495,   620,   621,   621,   621,   200,    -1,   127,   200,   497,
      -1,   498,    -1,   497,   498,    -1,   620,   620,   621,   200,
      -1,   620,   620,   620,   621,   200,    -1,   620,   620,   620,
     620,   621,   200,    -1,    84,   200,   500,    -1,   499,   500,
      -1,   620,   620,   501,   200,    -1,   620,    -1,   501,   620,
      -1,   287,   200,   503,    -1,   504,    -1,   503,   504,    -1,
     620,   620,   621,   200,    -1,   620,   620,   620,   621,   200,
      -1,   620,   620,   620,   620,   621,   200,    -1,    82,   200,
     506,    -1,   505,   506,    -1,   620,   620,   501,   200,    -1,
     130,   200,   508,    -1,   507,   508,    -1,   620,   620,   501,
     200,    -1,   134,   200,   510,    -1,   509,   510,    -1,   620,
     620,   501,   200,    -1,   133,   200,   512,    -1,   511,   512,
      -1,   620,   620,   501,   200,    -1,   620,   621,   200,    -1,
     513,    -1,   514,   513,    -1,    18,   200,   514,    -1,    19,
     200,   514,    -1,    13,   621,   200,    -1,   517,   500,    -1,
      15,   621,   200,    -1,   518,   508,    -1,    16,   621,   621,
     621,   200,    -1,   519,   512,    -1,    93,   620,   200,    -1,
      93,   620,   620,   200,    -1,    94,   620,   200,   599,    -1,
      20,   200,   621,   621,   200,   621,   621,   621,   200,    -1,
     522,   523,    -1,   620,   200,    -1,    21,   200,   621,   621,
     200,   621,   621,   621,   200,    -1,   524,   525,    -1,   620,
     620,   200,    -1,   620,   620,   620,   200,    -1,    22,   620,
     200,    -1,    23,   620,   200,    -1,    24,   200,   620,   200,
     621,   621,   200,    -1,   528,   529,    -1,   620,   620,   200,
      -1,   620,   620,   620,   200,    -1,    25,   200,    -1,   530,
     531,    -1,   620,   621,   621,   620,   200,    -1,   620,   621,
     621,   620,   621,   621,   200,    -1,   620,   621,   621,   620,
     621,   621,   621,   621,   621,   200,    -1,   251,    97,   200,
      -1,   251,    97,   620,   200,    -1,   145,   200,    -1,   145,
     333,   200,    -1,   145,   200,   540,    -1,   533,   180,   200,
     541,    -1,   146,   621,   200,    -1,   146,   200,    -1,   534,
     620,   621,   621,   621,   621,   621,   621,   200,    -1,   534,
     620,   621,   621,   621,   200,    -1,   113,   200,    -1,    32,
     200,    -1,   217,   620,   200,    -1,   535,   620,   621,   621,
     621,   621,   621,   621,   200,    -1,   218,   620,   200,    -1,
     536,   620,   621,   621,   621,   621,   621,   621,   200,    -1,
     151,   200,    -1,   151,   200,   540,    -1,   537,   180,   200,
     541,    -1,   149,   200,   542,    -1,    99,   200,    -1,    99,
     620,   200,    -1,   539,   574,    -1,   539,   620,   320,   620,
     620,   621,   200,    -1,   539,   620,   320,   620,   293,   620,
     620,   621,   200,    -1,   539,   296,   574,    -1,   574,    -1,
     540,   574,    -1,   620,   320,   620,   620,   621,   200,    -1,
     540,   620,   320,   620,   620,   621,   200,    -1,   620,   320,
     620,   293,   620,   620,   621,   200,    -1,   540,   620,   320,
     620,   293,   620,   620,   621,   200,    -1,   575,    -1,   541,
     575,    -1,   576,    -1,   542,   576,    -1,   332,   200,    -1,
     332,   200,   544,    -1,    48,   620,   200,   621,   621,   200,
      -1,   544,   621,   621,   200,    -1,   544,    48,   620,   200,
     621,   621,   200,    -1,   310,   200,    -1,   310,   200,   546,
      -1,    48,   620,   200,   621,   621,   200,    -1,   546,   621,
     621,   200,    -1,   546,    48,   620,   200,   621,   621,   200,
      -1,   168,   200,    -1,   168,   200,   548,    -1,   549,   550,
      -1,   548,   550,    -1,   548,   549,   550,    -1,   620,   200,
      -1,   620,   621,   200,    -1,   620,   621,   341,   620,   200,
      -1,   620,   621,   593,   200,    -1,   620,   620,   621,   200,
      -1,   137,   200,    -1,   137,   200,   552,    -1,   553,   554,
      -1,   552,   554,    -1,   552,   553,   554,    -1,   620,    51,
     621,   621,   200,    -1,   620,   621,   200,    -1,   620,   200,
      -1,   620,   620,   621,   621,   200,    -1,   620,   620,   621,
     200,    -1,   126,   200,   556,    -1,   126,   620,   200,   556,
      -1,   577,    -1,   556,   577,    -1,   175,   200,   558,    -1,
     557,   558,    -1,   620,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   200,    -1,
     620,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,    62,   621,   621,   200,    -1,
     620,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   264,   620,   621,   200,    -1,
     620,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   200,
      -1,   620,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
      62,   621,   621,   200,    -1,   620,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   264,   620,   621,   200,    -1,   620,
     621,   621,   621,   621,   621,   621,   621,   200,    -1,   620,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   200,    -1,   620,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,    62,   621,   200,
      -1,   620,    96,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   200,    -1,   620,    96,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   200,    -1,   620,    96,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   200,    -1,   620,    96,   621,   621,   200,    -1,   620,
      96,   621,   621,   621,   200,    -1,   620,    88,   620,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   620,   620,
     620,   200,    -1,   620,   317,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   200,    -1,   620,    57,   620,   621,
     200,    -1,   620,    57,   620,   621,   621,   621,   200,    -1,
     620,    57,   620,   621,   621,   621,   621,   200,    -1,   620,
      57,   620,   621,   621,   621,   621,   621,   200,    -1,   620,
      57,   620,   621,   621,   621,   621,   621,   621,   620,   200,
      -1,   620,    57,   620,   621,    87,   621,   621,   621,   621,
     200,    -1,   620,    57,   620,   621,    87,   621,   621,   621,
     621,   621,   200,    -1,   620,    57,   620,   621,    87,   621,
     621,   621,   621,   621,   621,   200,    -1,   620,    57,   620,
     621,   297,   621,   200,    -1,   620,    57,   620,   621,   297,
     621,   621,   200,    -1,   620,    57,   620,   621,   297,   621,
     621,   621,   200,    -1,   620,   297,   621,   200,    -1,   306,
     200,   572,    -1,   559,   572,    -1,   370,   620,   200,    -1,
     370,   620,   263,   200,    -1,   370,   620,   295,   621,   200,
      -1,   370,   620,   295,   621,   263,   200,    -1,   560,   620,
     620,   573,   200,    -1,   371,   620,   620,   200,    -1,   371,
     620,   620,   375,   200,    -1,   371,   620,   620,   376,   200,
      -1,   371,   620,   620,   374,   621,   200,    -1,   371,   620,
     620,   374,   621,   621,   200,    -1,   371,   620,   620,   375,
     374,   621,   200,    -1,   371,   620,   620,   375,   374,   621,
     621,   200,    -1,   371,   620,   620,   376,   374,   621,   200,
      -1,   371,   620,   620,   376,   374,   621,   621,   200,    -1,
     377,   620,   620,   200,    -1,   377,   620,   200,    -1,   313,
     200,    -1,   563,   620,   620,   620,   200,    -1,   563,   620,
     620,   620,   620,   200,    -1,   563,   620,   620,   620,   620,
     621,   200,    -1,   563,   620,   620,   620,   620,   621,   621,
     200,    -1,   563,   620,   620,   620,   620,   621,   621,   620,
     621,   200,    -1,   563,   620,   620,   620,   620,   621,   621,
     620,   621,   621,   200,    -1,   563,   620,   620,   620,   593,
     200,    -1,   563,   620,   620,   620,   620,   621,   621,   593,
     200,    -1,   106,   200,    -1,   564,   620,   620,   620,   200,
      -1,   564,   620,   620,   620,   621,   621,   200,    -1,   390,
     200,   566,    -1,   565,   566,    -1,   620,   620,   501,   200,
      -1,   392,   200,    -1,   567,   620,   620,   620,   200,    -1,
     567,   620,   620,   620,   621,   621,   200,    -1,    53,   200,
      -1,   568,   620,   620,   620,   200,    -1,   568,   620,   620,
     620,   620,   200,    -1,   568,   620,   620,   620,   620,   621,
     200,    -1,   568,   620,   620,   620,   620,   621,   621,   200,
      -1,   568,   620,   620,   620,   620,   621,   621,   620,   621,
     200,    -1,   568,   620,   620,   620,   620,   621,   621,   620,
     621,   621,   200,    -1,   568,   620,   620,   620,   620,   621,
     621,   620,   621,   621,   621,   621,   200,    -1,   568,   620,
     620,   620,   620,   621,   621,   620,   621,   621,   621,   621,
     621,   200,    -1,   568,   620,   620,   620,   593,   200,    -1,
     568,   620,   620,   620,   620,   621,   621,   593,   200,    -1,
      27,   620,   200,    -1,   111,   620,   200,    -1,   372,   621,
     200,    -1,   373,   620,   200,    -1,   205,   200,   571,    -1,
     570,   571,    -1,   620,   621,   621,   621,   200,    -1,   620,
     621,   621,   200,    -1,   620,   621,   200,    -1,   620,   620,
     573,   200,    -1,   620,    -1,   573,   620,    -1,   620,   620,
     621,   200,    -1,   620,   620,   200,    -1,   620,   621,   200,
      -1,   620,   621,   200,    -1,   620,   620,   621,   621,   200,
      -1,   620,   620,   621,   200,    -1,    60,   200,    -1,   578,
     580,    -1,    81,   200,    -1,   579,   580,    -1,   620,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   200,    -1,
     620,   316,   620,   200,    -1,    31,   200,    -1,   581,   620,
     620,   621,   621,   621,   200,    -1,     9,   200,    -1,   582,
     620,   620,   200,    -1,   582,   620,   620,   315,   621,   200,
      -1,   582,   620,   620,   620,   620,   200,    -1,   582,   620,
     620,   620,   620,   315,   621,   200,    -1,   582,   620,   620,
     620,   314,   621,   200,    -1,   582,   620,   620,   144,   200,
      -1,   582,   620,   620,   620,   200,    -1,   582,   620,   620,
     620,   620,   620,   200,    -1,   582,   620,   620,   620,   620,
     314,   621,   200,    -1,   232,   200,    -1,   583,   620,   621,
     200,    -1,   583,   620,   620,   621,   200,    -1,   583,   296,
     513,    -1,   173,   200,    -1,   173,   620,   200,    -1,   231,
     200,    -1,   585,   620,   621,   200,    -1,   585,   620,   320,
     620,   621,   200,    -1,   585,   620,   621,   621,   621,   200,
      -1,   585,   620,   320,   620,   621,   621,   621,   200,    -1,
     589,    -1,   588,    -1,   586,    58,   587,   200,    -1,   620,
      -1,   587,   620,    -1,   281,   200,   150,   200,    -1,   588,
     229,   620,   200,    -1,   588,   177,   620,   200,    -1,   588,
     304,   621,   200,    -1,   588,   178,   620,   200,    -1,   588,
     292,   620,   200,    -1,   281,   200,    67,   200,    -1,   281,
     200,    67,   620,   200,    -1,   281,   200,   268,   200,    -1,
     281,   200,   268,   235,   200,    -1,   281,   200,   268,   326,
     200,    -1,   281,   200,   150,   620,   200,    -1,   281,   200,
     150,   620,   621,   200,    -1,   281,   200,   150,   620,   621,
     620,   200,    -1,   281,   200,   150,   620,   621,   620,   620,
     200,    -1,   281,   200,   150,   620,   621,   620,   620,   620,
     200,    -1,   281,   200,    90,   620,   621,   200,    -1,   281,
     200,    90,   620,   200,    -1,   281,   200,    90,    69,   200,
      -1,   281,   200,    90,   350,   200,    -1,   281,   200,    90,
     620,    91,   200,    -1,   281,   200,    90,   620,   621,   620,
     200,    -1,    90,   620,   621,   620,   200,    -1,   281,   200,
      90,   200,    -1,    90,   200,    -1,    90,   620,   200,    -1,
      90,    69,   200,    -1,    90,   350,   200,    -1,    90,   620,
      91,   200,    -1,    30,   200,   268,   200,    -1,   279,   620,
     200,    -1,   280,   620,   200,    -1,   270,   621,   200,    -1,
     272,   620,   200,    -1,   273,   620,   200,    -1,   271,   620,
     200,    -1,   274,   621,   200,    -1,   275,   620,   200,    -1,
     277,   290,   200,    -1,   276,   620,   200,    -1,   278,   620,
     200,    -1,   193,   620,   620,   200,    -1,   194,   620,   621,
     200,    -1,   121,   621,   200,    -1,   122,   290,   200,    -1,
     589,   177,   620,   200,    -1,    76,   620,   620,   200,    -1,
      75,   620,   621,   200,    -1,   589,   229,    92,   200,    -1,
     589,   229,   173,   200,    -1,   589,   229,   620,   200,    -1,
     236,   620,   200,    -1,   236,   237,   200,    -1,   302,   621,
     200,    -1,   302,   621,   621,   200,    -1,   288,   621,   200,
      -1,   288,   621,   621,   200,    -1,   245,   621,   621,   200,
      -1,   247,   620,   620,   200,    -1,   246,   621,   621,   200,
      -1,   589,   178,   620,   200,    -1,   204,   200,    -1,   234,
     620,   200,    -1,   265,   620,   200,    -1,   265,   266,   200,
      -1,   185,   620,   200,    -1,   185,   266,   200,    -1,   107,
     620,   200,    -1,   107,   266,   200,    -1,   186,   200,    -1,
     108,   200,    -1,   110,   620,   200,    -1,   109,   200,    -1,
      52,   621,   200,    -1,   330,   200,    -1,    46,    47,   200,
      -1,    46,    47,   620,   200,    -1,    46,    11,   200,    -1,
      10,    11,   200,    -1,    10,    11,   250,   200,    -1,    10,
     345,   620,   200,    -1,    10,   345,   346,   620,   200,    -1,
      10,   345,   620,   351,   200,    -1,    10,   345,   346,   620,
     351,   200,    -1,   347,   621,   200,    -1,   347,   621,   621,
     200,    -1,   114,   621,   200,    -1,   118,   621,   200,    -1,
      49,   621,   200,    -1,   259,   290,   200,    -1,   325,   620,
     200,    -1,   233,   200,    -1,   171,   268,   200,    -1,    38,
     268,   200,    -1,    26,   268,   200,    -1,    50,   268,   200,
      -1,   171,   268,   235,   200,    -1,   171,   268,   289,   200,
      -1,    38,   268,   235,   200,    -1,    26,   268,   235,   200,
      -1,    38,   268,   289,   200,    -1,    50,   268,   235,   200,
      -1,    50,   268,   289,   200,    -1,   328,   620,   200,    -1,
     117,   200,    -1,   384,   200,    -1,   384,   290,   200,    -1,
     238,   620,   200,    -1,   203,   200,    -1,   203,   620,   200,
      -1,   125,   620,   621,   200,    -1,   125,   620,   621,   620,
     200,    -1,   125,   620,   621,   620,   621,   200,    -1,   125,
     620,   621,   620,   621,   620,   200,    -1,   125,   200,    -1,
     209,   620,   200,    -1,   209,   620,   621,   200,    -1,   311,
     621,   200,    -1,   286,   620,   200,    -1,   160,   620,   200,
      -1,   160,   156,   200,    -1,   148,   200,    -1,   285,   200,
      -1,   312,   620,   200,    -1,   312,   620,   621,   200,    -1,
     312,   620,   621,   621,   200,    -1,   344,   150,   200,    -1,
     344,   150,   141,   200,    -1,   183,   620,   200,    -1,   183,
     184,   200,    -1,   181,   620,   200,    -1,   181,   182,   200,
      -1,   181,   620,   188,   620,   200,    -1,   181,   182,   188,
     620,   200,    -1,   181,   620,   620,   200,    -1,   181,   182,
     187,   200,    -1,   181,   620,   620,   188,   620,   200,    -1,
     181,   182,   187,   188,   620,   200,    -1,   191,   620,   200,
      -1,   240,   621,   200,    -1,   172,   620,   621,   200,    -1,
      34,   290,   200,    -1,    74,   290,   200,    -1,    54,   290,
     200,    -1,    55,   290,   200,    -1,   192,   620,   200,    -1,
      89,   200,   591,    -1,   621,   200,    -1,    77,   593,   200,
      -1,    77,   200,   593,   200,    -1,    67,    -1,    67,   621,
      -1,    67,   621,   621,    -1,    67,   621,   621,   621,    -1,
      67,   621,   621,   621,   621,    -1,    78,    -1,    79,   621,
      -1,    78,    79,   621,    -1,   593,   142,   620,    -1,   593,
     142,   620,   621,    -1,   129,   200,    -1,   101,   200,    -1,
     329,   621,   200,    -1,   329,   621,   621,   621,   200,    -1,
      33,   620,   200,    -1,    33,   620,   621,   200,    -1,    33,
     620,   621,   621,   200,    -1,   243,   620,   200,   596,    -1,
     244,   620,   200,   596,    -1,   152,   620,   200,   596,    -1,
     129,   200,   139,   621,   621,   620,   200,    -1,   129,   200,
     140,   621,   621,   620,   200,    -1,   129,   200,   138,   621,
     621,   620,   620,   200,    -1,   129,   200,   595,    -1,   129,
     200,   139,   621,   621,   620,   621,   621,   200,    -1,   129,
     200,   140,   621,   621,   620,   621,   621,   200,    -1,   129,
     200,   138,   621,   621,   620,   620,   621,   621,   200,    -1,
     594,   478,    -1,   594,   414,    -1,   594,   411,    -1,   138,
     621,   200,    -1,   138,   621,   621,   621,   200,    -1,   595,
     621,   620,   200,    -1,   597,    -1,   596,   597,    -1,   621,
     621,   621,   200,    -1,   161,   200,   621,   621,   621,   200,
      -1,   598,   621,   621,   621,   200,    -1,   600,    -1,   599,
     600,    -1,   621,   621,   621,   200,    -1,   131,   200,   602,
      -1,   621,   200,    -1,   132,   200,   604,    -1,   621,   200,
      -1,    64,   200,    -1,   201,   200,    -1,   606,     8,   200,
      -1,   606,   202,   200,    -1,   606,   177,   620,   200,    -1,
     606,   208,   621,   200,    -1,   606,   208,   621,   621,   200,
      -1,   606,   208,   621,   621,   621,   621,   200,    -1,   606,
      68,   621,   621,   200,    -1,   606,    68,   621,   621,   620,
     620,   200,    -1,   606,    95,   620,   200,    -1,   606,    95,
     620,   620,   200,    -1,   606,   326,   200,    -1,   606,   167,
     621,   200,    -1,   606,   112,   200,    -1,   606,   112,   621,
     200,    -1,   606,   607,    -1,   252,   620,   200,    -1,   252,
     620,   620,   200,    -1,   252,   200,   298,   620,   200,    -1,
     252,   200,   298,   620,   200,   230,   620,   200,    -1,   255,
     200,    -1,    45,   200,    -1,    45,   200,    97,   200,   620,
     200,    97,   200,    97,   200,    -1,    45,   200,    97,   200,
     620,   200,    97,   200,    97,   200,    97,   200,    -1,   211,
      97,   200,    -1,   610,   200,    97,   200,    -1,   340,   200,
      -1,   340,   341,   620,   200,    -1,   611,   620,   620,   621,
     621,   621,   200,    -1,   611,   620,   620,   621,   621,   621,
     343,   621,   200,    -1,   611,   620,   620,   621,   621,   621,
     341,   620,   200,    -1,   611,   620,   620,   621,   621,   621,
     593,   200,    -1,   611,   620,   620,   621,   621,   621,   343,
     621,   341,   620,   200,    -1,   611,   620,   620,   621,   621,
     621,   343,   621,   341,   620,   593,   200,    -1,   352,   200,
      -1,   612,   620,   354,   621,   621,   621,   621,   621,   200,
      -1,   612,   620,   354,   621,   621,   621,   621,   621,   621,
     200,    -1,   612,   620,   355,   621,   621,   621,   621,   621,
     200,    -1,   612,   620,   355,   621,   621,   621,   621,   621,
     621,   200,    -1,   612,   620,   368,   621,   621,   621,   621,
     621,   200,    -1,   612,   620,   368,   621,   621,   621,   621,
     621,   621,   200,    -1,   612,   620,   356,   621,   621,   621,
     200,    -1,   612,   620,   357,   621,   621,   621,   200,    -1,
     612,   620,   367,   621,   621,   621,   200,    -1,   612,   620,
     358,   621,   621,   621,   621,   200,    -1,   612,   620,   369,
     621,   621,   621,   621,   200,    -1,   612,   620,   361,   621,
     621,   621,   200,    -1,   612,   620,   362,   621,   621,   621,
     200,    -1,   612,   620,   362,   621,   621,   621,   621,   200,
      -1,   612,   620,   366,   621,   621,   621,   621,   200,    -1,
     612,   620,   366,   621,   621,   621,   621,   621,   200,    -1,
     612,   620,   363,   621,   621,   621,   621,   621,   621,   200,
      -1,   612,   620,   364,   621,   621,   621,   621,   621,   621,
     200,    -1,   612,   620,   364,   621,   621,   621,   621,   621,
     621,   621,   200,    -1,   612,   620,   364,   621,   621,   621,
     621,   621,   621,   621,   621,   200,    -1,   612,   620,   360,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     200,    -1,   612,   359,    97,    97,   200,    -1,   612,   620,
      97,   614,   200,    -1,   353,   200,    -1,   613,   620,   620,
     200,    -1,   613,   620,   620,   620,   200,    -1,    -1,   614,
     621,    -1,   253,   200,   254,   200,    -1,   253,   200,   254,
     200,   254,   200,    -1,   253,   200,   254,   200,   254,   200,
     254,   200,    -1,   397,   200,    -1,   616,   617,   200,    -1,
     393,    97,    -1,   393,    97,   620,    -1,   399,   620,    -1,
     398,   200,    -1,   618,   619,   200,    -1,   394,    97,    -1,
     395,    97,    -1,   401,   621,    -1,   380,   620,    -1,   399,
     620,    -1,   147,    -1,   147,    -1,    63,    -1,   158,    -1,
      86,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   149,   149,   155,   156,   159,   160,   162,   164,   165,
     166,   167,   168,   169,   171,   172,   173,   174,   176,   178,
     179,   180,   181,   182,   183,   185,   186,   187,   188,   189,
     190,   191,   192,   193,   194,   195,   196,   197,   198,   199,
     200,   201,   202,   204,   206,   208,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   220,   221,   222,
     223,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     235,   237,   239,   241,   243,   245,   246,   247,   248,   249,
     250,   251,   252,   253,   254,   255,   256,   257,   258,   259,
     260,   261,   263,   264,   266,   267,   268,   269,   270,   271,
     272,   273,   274,   275,   276,   278,   280,   281,   282,   283,
     284,   285,   287,   288,   289,   290,   291,   293,   294,   295,
     296,   298,   300,   302,   304,   306,   308,   310,   312,   313,
     314,   315,   316,   317,   320,   327,   333,   334,   339,   351,
     357,   358,   362,   366,   368,   370,   372,   373,   374,   375,
     376,   379,   383,   385,   390,   393,   397,   436,   485,   487,
     489,   493,   494,   496,   498,   500,   502,   506,   514,   516,
     518,   520,   522,   524,   526,   528,   532,   534,   542,   544,
     548,   550,   554,   556,   560,   561,   562,   563,   564,   565,
     568,   569,   573,   577,   579,   583,   585,   589,   591,   595,
     597,   601,   603,   607,   613,   620,   629,   633,   634,   638,
     645,   651,   656,   662,   663,   665,   669,   670,   674,   675,
     679,   682,   687,   692,   696,   701,   707,   714,   721,   723,
     725,   727,   729,   733,   735,   737,   739,   742,   744,   746,
     748,   750,   752,   754,   756,   759,   761,   763,   765,   767,
     769,   771,   773,   775,   779,   783,   784,   787,   790,   792,
     794,   796,   798,   800,   802,   804,   806,   808,   811,   815,
     818,   821,   823,   825,   828,   830,   832,   836,   838,   842,
     846,   849,   853,   857,   859,   860,   861,   862,   864,   866,
     868,   875,   877,   879,   881,   883,   889,   891,   892,   894,
     897,   899,   901,   903,   909,   920,   923,   924,   925,   929,
     931,   935,   945,   948,   951,   959,   964,   966,   968,   972,
     974,   978,   982,   986,   990,   994,   998,  1002,  1010,  1013,
    1016,  1021,  1023,  1027,  1029,  1031,  1034,  1042,  1052,  1056,
    1060,  1064,  1070,  1075,  1084,  1085,  1088,  1090,  1092,  1094,
    1096,  1098,  1100,  1102,  1104,  1106,  1110,  1112,  1117,  1121,
    1127,  1135,  1141,  1146,  1149,  1155,  1163,  1165,  1167,  1169,
    1171,  1184,  1186,  1195,  1199,  1201,  1205,  1207,  1211,  1252,
    1256,  1262,  1270,  1272,  1275,  1284,  1286,  1289,  1298,  1300,
    1305,  1307,  1312,  1315,  1316,  1319,  1321,  1323,  1327,  1328,
    1331,  1337,  1339,  1344,  1347,  1348,  1351,  1354,  1357,  1362,
    1363,  1366,  1371,  1372,  1375,  1379,  1380,  1383,  1390,  1391,
    1394,  1398,  1402,  1404,  1408,  1412,  1416,  1418,  1421,  1423,
    1426,  1430,  1433,  1435,  1439,  1445,  1450,  1454,  1458,  1463,
    1467,  1469,  1473,  1481,  1488,  1493,  1496,  1498,  1502,  1504,
    1508,  1510,  1512,  1516,  1519,  1525,  1527,  1529,  1532,  1538,
    1540,  1542,  1553,  1564,  1567,  1572,  1574,  1587,  1589,  1601,
    1603,  1606,  1612,  1617,  1619,  1621,  1623,  1625,  1627,  1635,
    1637,  1639,  1641,  1643,  1645,  1649,  1651,  1655,  1657,  1661,
    1662,  1665,  1667,  1669,  1673,  1674,  1677,  1679,  1681,  1685,
    1686,  1689,  1693,  1695,  1701,  1704,  1707,  1711,  1718,  1730,
    1731,  1734,  1736,  1738,  1742,  1744,  1746,  1750,  1758,  1768,
    1770,  1775,  1777,  1781,  1782,  1785,  1793,  1802,  1813,  1822,
    1832,  1844,  1851,  1859,  1867,  1883,  1900,  1918,  1926,  1934,
    1953,  1961,  1968,  1977,  1987,  1998,  2011,  2024,  2038,  2053,
    2062,  2072,  2083,  2091,  2092,  2095,  2101,  2107,  2114,  2122,
    2146,  2148,  2150,  2154,  2156,  2158,  2160,  2162,  2166,  2172,
    2174,  2182,  2184,  2190,  2197,  2204,  2211,  2219,  2228,  2235,
    2245,  2247,  2253,  2261,  2264,  2267,  2273,  2275,  2281,  2291,
    2293,  2300,  2307,  2314,  2321,  2329,  2338,  2347,  2356,  2364,
    2374,  2376,  2378,  2380,  2383,  2385,  2389,  2391,  2393,  2397,
    2402,  2404,  2409,  2411,  2415,  2419,  2423,  2425,  2429,  2430,
    2434,  2439,  2443,  2448,  2452,  2453,  2460,  2462,  2464,  2468,
    2470,  2474,  2476,  2481,  2486,  2491,  2498,  2499,  2501,  2506,
    2512,  2517,  2524,  2526,  2528,  2533,  2536,  2544,  2545,  2546,
    2549,  2551,  2555,  2559,  2561,  2563,  2565,  2567,  2571,  2574,
    2577,  2580,  2585,  2589,  2592,  2595,  2598,  2601,  2609,  2615,
    2619,  2624,  2630,  2636,  2642,  2648,  2651,  2653,  2656,  2660,
    2665,  2670,  2676,  2678,  2680,  2682,  2684,  2691,  2693,  2700,
    2702,  2704,  2706,  2708,  2710,  2712,  2714,  2716,  2718,  2720,
    2726,  2728,  2730,  2737,  2744,  2751,  2753,  2756,  2758,  2761,
    2764,  2767,  2770,  2772,  2774,  2794,  2797,  2799,  2802,  2804,
    2807,  2809,  2811,  2813,  2815,  2817,  2819,  2821,  2823,  2827,
    2836,  2857,  2875,  2880,  2886,  2892,  2899,  2901,  2904,  2906,
    2908,  2910,  2912,  2914,  2916,  2918,  2920,  2922,  2924,  2928,
    2931,  2935,  2939,  2942,  2946,  2950,  2976,  2978,  2980,  2982,
    2988,  2993,  2998,  3011,  3025,  3042,  3058,  3068,  3070,  3073,
    3075,  3077,  3079,  3081,  3083,  3085,  3090,  3096,  3103,  3105,
    3108,  3110,  3112,  3114,  3116,  3119,  3122,  3125,  3128,  3132,
    3136,  3138,  3140,  3143,  3145,  3147,  3149,  3151,  3155,  3158,
    3169,  3173,  3179,  3186,  3193,  3201,  3210,  3220,  3226,  3232,
    3238,  3241,  3246,  3249,  3259,  3266,  3274,  3279,  3285,  3292,
    3297,  3302,  3307,  3315,  3323,  3331,  3332,  3341,  3350,  3359,
    3360,  3361,  3364,  3373,  3382,  3387,  3388,  3391,  3395,  3398,
    3403,  3404,  3407,  3411,  3414,  3425,  3428,  3439,  3445,  3458,
    3463,  3470,  3472,  3474,  3477,  3482,  3485,  3490,  3493,  3496,
    3498,  3504,  3506,  3509,  3512,  3517,  3527,  3534,  3546,  3550,
    3551,  3553,  3563,  3569,  3577,  3578,  3582,  3585,  3588,  3590,
    3593,  3595,  3599,  3600,  3605,  3610,  3615,  3620,  3625,  3630,
    3635,  3640,  3645,  3650,  3655,  3661,  3667,  3673,  3679,  3685,
    3691,  3697,  3703,  3709,  3714,  3718,  3724,  3725,  3727,  3734,
    3735,  3746,  3750,  3753,  3760,  3764,  3768,  3770,  3777,  3782,
    3786,  3790,  3792,  3794,  3796,  3798,  3803,  3808,  3810,  3812,
    3814
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ACTUATORS", "AERO", "AEROH", "AEROTYPE",
  "ANALYSIS", "ARCLENGTH", "ATTRIBUTES", "AUGMENT", "AUGMENTTYPE",
  "AVERAGED", "ATDARB", "ACOU", "ATDDNB", "ATDROB", "ARPACK", "ATDDIR",
  "ATDNEU", "AXIHDIR", "AXIHNEU", "AXINUMMODES", "AXINUMSLICES",
  "AXIHSOMMER", "AXIMPC", "AUXCOARSESOLVER", "ACMECNTL", "ADDEDMASS",
  "AEROEMBED", "BLOCKDIAG", "BOFFSET", "BUCKLE", "BGTL", "BMPC",
  "BINARYINPUT", "BINARYOUTPUT", "CHECKTOKEN", "COARSESOLVER", "COEF",
  "CFRAMES", "COLLOCATEDTYPE", "CONVECTION", "COMPOSITE", "CONDITION",
  "CONTROL", "CORNER", "CORNERTYPE", "CURVE", "CCTTOL", "CCTSOLVER",
  "CRHS", "COUPLEDSCALE", "CONTACTSURFACES", "CMPC", "CNORM",
  "COMPLEXOUTTYPE", "CONSTRMAT", "CASES", "CONSTRAINEDSURFACES",
  "CSFRAMES", "CSTYPE", "DAMPING", "DblConstant", "DEM", "DIMASS", "DISP",
  "DIRECT", "DLAMBDA", "DP", "DYNAM", "DETER", "DECOMPOSE", "DECOMPFILE",
  "DMPC", "DEBUGCNTL", "DEBUGICNTL", "CONSTRAINTS", "MULTIPLIERS",
  "PENALTY", "EIGEN", "EFRAMES", "ELSCATTERER", "END", "ELHSOMMERFELD",
  "EXPLICIT", "EPSILON", "ELEMENTARYFUNCTIONTYPE", "FABMAT", "FACOUSTICS",
  "FETI", "FETI2TYPE", "FETIPREC", "FFP", "FFPDIR", "FITALG", "FLUMAT",
  "FNAME", "FLUX", "FORCE", "FRONTAL", "FETIH", "FILTEREIG", "FREQSWEEP",
  "FREQSWEEP1", "FREQSWEEP2", "FSINTERFACE", "FSISCALING", "FSIELEMENT",
  "NOLOCALFSISPLITING", "FSICORNER", "FFIDEBUG", "FAILSAFE", "GEPS",
  "GLOBALTOL", "GRAVITY", "GRBM", "GTGSOLVER", "GLOBALCRBMTOL", "GROUP",
  "GROUPTYPE", "GOLDFARBTOL", "GOLDFARBCHECK", "HDIRICHLET", "HEAT",
  "HFETI", "HNEUMAN", "HSOMMERFELD", "HFTT", "HELMHOLTZ", "HNBO", "HELMMF",
  "HELMSO", "HSCBO", "HWIBO", "HZEM", "HZEMFILTER", "HLMPC", "HELMSWEEP",
  "HELMSWEEP1", "HELMSWEEP2", "HERMITIAN", "HESSIAN", "IACC", "IDENTITY",
  "IDIS", "IDIS6", "IntConstant", "INTERFACELUMPED", "ITEMP", "ITERTYPE",
  "IVEL", "INCIDENCE", "IHDIRICHLET", "IHDSWEEP", "IHNEUMANN",
  "ISOLVERTYPE", "INPC", "INFINTY", "JACOBI", "KRYLOVTYPE", "KIRLOC",
  "LAYC", "LAYN", "LAYD", "LAYO", "LAYMAT", "LFACTOR", "LMPC", "LOAD",
  "LOBPCG", "LOCALSOLVER", "LINESEARCH", "LUMPED", "MASS", "MATERIALS",
  "MATLAB", "MAXITR", "MAXORTHO", "MAXVEC", "MODAL", "MPCPRECNO",
  "MPCPRECNOID", "MPCTYPE", "MPCTYPEID", "MPCSCALING", "MPCELEMENT",
  "MPCBLOCKID", "MPCBLK_OVERLAP", "MFTT", "MPTT", "MRHS", "MPCCHECK",
  "MUMPSICNTL", "MUMPSCNTL", "MECH", "MODEFILTER", "NDTYPE", "NEIGPA",
  "NEWMARK", "NewLine", "NL", "NLMAT", "NLPREC", "NOCOARSE", "NODETOKEN",
  "NONINPC", "NSBSPV", "NLTOL", "NUMCGM", "NOSECONDARY", "OPTIMIZATION",
  "OUTPUT", "OUTPUT6", "QSTATIC", "QLOAD", "PITA", "PITADISP6", "PITAVEL6",
  "NOFORCE", "MDPITA", "GLOBALBASES", "LOCALBASES", "TIMEREVERSIBLE",
  "REMOTECOARSE", "ORTHOPROJTOL", "READINITSEED", "JUMPCVG", "JUMPOUTPUT",
  "PRECNO", "PRECONDITIONER", "PRELOAD", "PRESSURE", "PRINTMATLAB", "PROJ",
  "PIVOT", "PRECTYPE", "PRECTYPEID", "PICKANYCORNER", "PADEPIVOT",
  "PROPORTIONING", "PLOAD", "PADEPOLES", "POINTSOURCE", "PLANEWAVE",
  "PTOL", "PLANTOL", "PMAXIT", "RADIATION", "RBMFILTER", "RBMSET",
  "READMODE", "REBUILD", "RENUM", "RENUMBERID", "REORTHO", "RESTART",
  "RECONS", "RECONSALG", "REBUILDCCT", "RANDOM", "RPROP", "RNORM",
  "REVERSENORMALS", "RIGID", "SCALING", "SCALINGTYPE", "SENSORS",
  "SOLVERTYPE", "SHIFT", "SPOOLESTAU", "SPOOLESSEED", "SPOOLESMAXSIZE",
  "SPOOLESMAXDOMAINSIZE", "SPOOLESMAXZEROS", "SPOOLESMSGLVL",
  "SPOOLESSCALE", "SPOOLESPIVOT", "SPOOLESRENUM", "SPARSEMAXSUP",
  "SPARSEDEFBLK", "STATS", "STRESSID", "SUBSPACE", "SURFACE",
  "SAVEMEMCOARSE", "SPACEDIMENSION", "SCATTERER", "STAGTOL", "SCALED",
  "SWITCH", "STABLE", "SUBTYPE", "STEP", "SOWER", "SHELLTHICKNESS", "SURF",
  "SPRINGMAT", "TANGENT", "TEMP", "TIME", "TOLEIG", "TOLFETI", "TOLJAC",
  "TOLPCG", "TOPFILE", "TOPOLOGY", "TRBM", "THERMOE", "THERMOH", "TETT",
  "TOLCGM", "TURKEL", "TIEDSURFACES", "THETA", "HRC", "THIRDNODE",
  "THERMMAT", "TDENFORC", "TESTULRICH", "THRU", "TOPFLAG", "USE",
  "USERDEFINEDISP", "USERDEFINEFORCE", "UPROJ", "UNSYMMETRIC", "USING",
  "VERSION", "WAVENUMBER", "WETCORNERS", "XPOST", "YMTT", "ZERO", "BINARY",
  "GEOMETRY", "DECOMPOSITION", "GLOBAL", "MATCHER", "CPUMAP",
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
  "OUTPUTWEIGHT", "WEIGHTLIST", "GMRESRESIDUAL", "SLOSH", "SLGRAV",
  "SLZEM", "SLZEMFILTER", "PDIR", "HEFSB", "HEFRS", "HEINTERFACE",
  "SNAPFI", "PODROB", "TRNVCT", "ORTHOG", "SVDTOKEN", "SAMPLING",
  "PODSIZEMAX", "REFSUBSTRACT", "TOLER", "$accept", "FinalizedData", "All",
  "Component", "Noninpc", "Inpc", "Group", "Random", "Impe",
  "PadePivotInfo", "PadePolesInfo", "FreqSweep", "ReconsInfo", "Sower",
  "BinarySpec", "AnalysisInfo", "Decompose", "WeightList", "MFTTInfo",
  "MPTTInfo", "HFTTInfo", "Composites", "Cframes", "CoefInfo", "CoefList",
  "LaycInfo", "LaynInfo", "LaydInfo", "LayoInfo", "LayData", "LayoData",
  "LayMat", "LayMatData", "DiscrMasses", "Gravity", "Restart", "LoadCInfo",
  "UseCInfo", "SensorLocations", "ActuatorLocations", "UsdfLocations",
  "UsddLocations", "Output", "OutInfo", "RenumberOutput", "DynInfo",
  "SloshInfo", "MassInfo", "CondInfo", "TopInfo", "DynamInfo",
  "TimeIntegration", "NewmarkSecondOrder", "NewmarkFirstOrder",
  "QstaticInfo", "QMechInfo", "QHeatInfo", "AeroInfo",
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
  "IDisp6Pita", "IVel6Pita", "IVel", "ITemp", "NeumanBC", "BCDataList",
  "ModalValList", "TBCDataList", "YMTTable", "YMTTList", "TETTable",
  "TETTList", "LMPConstrain", "MPCList", "MPCHeader", "MPCLine",
  "ComplexLMPConstrain", "ComplexMPCList", "ComplexMPCHeader",
  "ComplexMPCLine", "ComplexNeumanBC", "ComplexBCDataList", "Materials",
  "MatData", "ElemSet", "FaceSet", "MortarCondition", "WetInterface",
  "TiedSurfaces", "FSInterface", "HEVibInfo", "HEVInterfaceElement",
  "HEInterface", "ContactSurfaces", "AcmeControls", "NodeSet", "Node",
  "Element", "NodeNums", "BC_Data", "ModalVal", "TBC_Data",
  "ComplexBC_Data", "ConstrainedSurfaceFrameDList", "FrameDList", "Frame",
  "BoffsetList", "Attributes", "Pressure", "Lumped", "Preload", "Statics",
  "CasesList", "IterSolver", "Solver", "OldHelmInfo", "FAcousticData",
  "Constraints", "ConstraintOptionsData", "HelmInfo", "HelmSweep",
  "IncidenceList", "IncidenceVector", "KirchhoffLocations", "FFPDirList",
  "FFPDirVector", "HelmMFInfo", "FAcousticDataMF", "HelmSOInfo",
  "FAcousticDataSO", "DEMInfo", "NLInfo", "NewtonInfo", "OrthoInfo",
  "Control", "Optimization", "NodalContact", "MatSpec", "MatUsage",
  "FloatList", "Renumbering", "SvdToken", "SvdOption", "Sampling",
  "SamplingOption", "Integer", "Float", 0
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
     655,   656
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   402,   403,   404,   404,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   405,   405,   405,   405,   405,   405,
     405,   405,   405,   405,   406,   407,   408,   408,   408,   408,
     409,   409,   410,   410,   410,   410,   410,   410,   410,   410,
     410,   411,   412,   412,   413,   413,   414,   414,   415,   415,
     415,   416,   416,   416,   416,   416,   416,   417,   418,   418,
     418,   418,   418,   418,   418,   418,   419,   419,   420,   420,
     421,   421,   422,   422,   423,   423,   423,   423,   423,   423,
     424,   424,   425,   426,   426,   427,   427,   428,   428,   429,
     429,   430,   430,   431,   431,   431,   432,   433,   433,   434,
     434,   434,   434,   435,   435,   435,   436,   436,   437,   437,
     437,   437,   438,   439,   440,   441,   442,   443,   444,   444,
     444,   444,   444,   445,   445,   445,   445,   445,   445,   445,
     445,   445,   445,   445,   445,   445,   445,   445,   445,   445,
     445,   445,   445,   445,   446,   447,   447,   447,   447,   447,
     447,   447,   447,   447,   447,   447,   447,   447,   447,   447,
     447,   447,   447,   447,   447,   447,   447,   448,   448,   449,
     450,   450,   451,   452,   452,   452,   452,   452,   452,   452,
     452,   452,   452,   452,   452,   452,   453,   453,   453,   453,
     454,   454,   454,   454,   455,   456,   456,   456,   456,   457,
     457,   458,   459,   459,   459,   459,   459,   459,   459,   460,
     460,   461,   462,   463,   464,   465,   466,   467,   468,   468,
     468,   469,   469,   470,   470,   470,   471,   471,   472,   473,
     474,   475,   475,   475,   476,   476,   477,   477,   477,   477,
     477,   477,   477,   477,   477,   477,   478,   478,   479,   480,
     480,   481,   481,   481,   482,   483,   484,   484,   484,   484,
     484,   485,   485,   486,   487,   487,   488,   488,   489,   490,
     490,   491,   492,   492,   492,   493,   493,   493,   494,   494,
     495,   495,   496,   497,   497,   498,   498,   498,   499,   499,
     500,   501,   501,   502,   503,   503,   504,   504,   504,   505,
     505,   506,   507,   507,   508,   509,   509,   510,   511,   511,
     512,   513,   514,   514,   515,   516,   517,   517,   518,   518,
     519,   519,   520,   520,   521,   522,   522,   523,   524,   524,
     525,   525,   526,   527,   528,   528,   529,   529,   530,   530,
     531,   531,   531,   532,   532,   533,   533,   533,   533,   534,
     534,   534,   534,   534,   534,   535,   535,   536,   536,   537,
     537,   537,   538,   539,   539,   539,   539,   539,   539,   540,
     540,   540,   540,   540,   540,   541,   541,   542,   542,   543,
     543,   544,   544,   544,   545,   545,   546,   546,   546,   547,
     547,   548,   548,   548,   549,   549,   549,   549,   550,   551,
     551,   552,   552,   552,   553,   553,   553,   554,   554,   555,
     555,   556,   556,   557,   557,   558,   558,   558,   558,   558,
     558,   558,   558,   558,   558,   558,   558,   558,   558,   558,
     558,   558,   558,   558,   558,   558,   558,   558,   558,   558,
     558,   558,   558,   559,   559,   560,   560,   560,   560,   560,
     561,   561,   561,   561,   561,   561,   561,   561,   561,   562,
     562,   563,   563,   563,   563,   563,   563,   563,   563,   563,
     564,   564,   564,   565,   565,   566,   567,   567,   567,   568,
     568,   568,   568,   568,   568,   568,   568,   568,   568,   568,
     569,   569,   569,   569,   570,   570,   571,   571,   571,   572,
     573,   573,   574,   574,   575,   576,   577,   577,   578,   578,
     579,   579,   580,   580,   581,   581,   582,   582,   582,   582,
     582,   582,   582,   582,   582,   582,   583,   583,   583,   583,
     584,   584,   585,   585,   585,   585,   585,   586,   586,   586,
     587,   587,   588,   588,   588,   588,   588,   588,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   589,   589,
     589,   589,   589,   589,   589,   589,   589,   589,   590,   591,
     592,   592,   593,   593,   593,   593,   593,   593,   593,   593,
     593,   593,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   594,   594,   594,   594,   594,   594,   594,   594,
     594,   594,   595,   595,   595,   596,   596,   597,   598,   598,
     599,   599,   600,   601,   602,   603,   604,   605,   606,   606,
     606,   606,   606,   606,   606,   606,   606,   606,   606,   606,
     606,   606,   606,   606,   607,   607,   607,   607,   608,   609,
     609,   609,   610,   610,   611,   611,   611,   611,   611,   611,
     611,   611,   612,   612,   612,   612,   612,   612,   612,   612,
     612,   612,   612,   612,   612,   612,   612,   612,   612,   612,
     612,   612,   612,   612,   612,   612,   613,   613,   613,   614,
     614,   615,   615,   615,   616,   616,   617,   617,   617,   618,
     618,   619,   619,   619,   619,   619,   620,   621,   621,   621,
     621
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
       1,     1,     1,     1,     5,     4,     2,     5,     6,     6,
       2,     6,     5,     7,     7,     8,     3,     2,     2,     2,
       2,     3,     2,     4,     3,     4,     4,     6,     2,     3,
       3,     2,     4,     4,     4,     4,     4,     3,     2,     4,
       4,     3,     3,     3,     3,     3,     2,     4,     2,     4,
       2,     4,     2,     4,     2,     2,     2,     2,     2,     2,
       2,     2,     4,     4,     5,     3,     2,     3,     2,     3,
       2,     3,     2,    11,    13,    14,     5,     2,     2,     7,
       9,    11,    12,     2,     5,     6,     2,     5,     2,     5,
       4,     4,     3,     3,     3,     3,     3,     3,     2,     2,
       3,     3,     3,     3,     5,     4,     5,     6,     7,     3,
       5,     4,     5,     6,     7,     3,     2,     3,     2,     2,
       3,     2,     3,     2,     3,     1,     2,     3,     3,     2,
       5,     3,     3,     3,     2,     3,     2,     3,     4,     4,
       5,     3,     3,     2,     4,     2,     3,     2,     3,     2,
       5,     2,     2,     2,     2,     2,     2,     3,     4,     4,
       8,     4,     4,     3,     3,     5,     2,     3,     3,     3,
       2,     4,     1,     3,     1,     2,     2,     2,     3,     5,
       6,     6,     3,     4,     6,     5,     3,     4,     3,     1,
       2,     6,     2,     2,     2,     2,     2,     4,     5,     4,
       2,     3,     2,     3,     2,     4,     1,     2,     2,     2,
       5,     5,     6,     7,     3,     1,     1,     2,     1,     1,
       1,     2,     1,     1,     2,     1,     4,     4,     3,     6,
       5,     6,     3,     2,     5,     6,     2,     2,     7,     9,
       3,     2,     6,     3,     1,     2,     3,     2,     3,     1,
       2,     4,     2,     4,     5,     2,     4,     5,     2,     6,
       2,     6,     3,     1,     2,     4,     5,     6,     3,     2,
       4,     1,     2,     3,     1,     2,     4,     5,     6,     3,
       2,     4,     3,     2,     4,     3,     2,     4,     3,     2,
       4,     3,     1,     2,     3,     3,     3,     2,     3,     2,
       5,     2,     3,     4,     4,     9,     2,     2,     9,     2,
       3,     4,     3,     3,     7,     2,     3,     4,     2,     2,
       5,     7,    10,     3,     4,     2,     3,     3,     4,     3,
       2,     9,     6,     2,     2,     3,     9,     3,     9,     2,
       3,     4,     3,     2,     3,     2,     7,     9,     3,     1,
       2,     6,     7,     8,     9,     1,     2,     1,     2,     2,
       3,     6,     4,     7,     2,     3,     6,     4,     7,     2,
       3,     2,     2,     3,     2,     3,     5,     4,     4,     2,
       3,     2,     2,     3,     5,     3,     2,     5,     4,     3,
       4,     1,     2,     3,     2,    16,    19,    19,    20,    23,
      23,     9,    14,    16,    12,    13,    14,     5,     6,    16,
      12,     5,     7,     8,     9,    11,    10,    11,    12,     7,
       8,     9,     4,     3,     2,     3,     4,     5,     6,     5,
       4,     5,     5,     6,     7,     7,     8,     7,     8,     4,
       3,     2,     5,     6,     7,     8,    10,    11,     6,     9,
       2,     5,     7,     3,     2,     4,     2,     5,     7,     2,
       5,     6,     7,     8,    10,    11,    13,    14,     6,     9,
       3,     3,     3,     3,     3,     2,     5,     4,     3,     4,
       1,     2,     4,     3,     3,     3,     5,     4,     2,     2,
       2,     2,    11,     4,     2,     7,     2,     4,     6,     6,
       8,     7,     5,     5,     7,     8,     2,     4,     5,     3,
       2,     3,     2,     4,     6,     6,     8,     1,     1,     4,
       1,     2,     4,     4,     4,     4,     4,     4,     4,     5,
       4,     5,     5,     5,     6,     7,     8,     9,     6,     5,
       5,     5,     6,     7,     5,     4,     2,     3,     3,     3,
       4,     4,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     4,     4,     3,     3,     4,     4,     4,
       4,     4,     4,     3,     3,     3,     4,     3,     4,     4,
       4,     4,     4,     2,     3,     3,     3,     3,     3,     3,
       3,     2,     2,     3,     2,     3,     2,     3,     4,     3,
       3,     4,     4,     5,     5,     6,     3,     4,     3,     3,
       3,     3,     3,     2,     3,     3,     3,     3,     4,     4,
       4,     4,     4,     4,     4,     3,     2,     2,     3,     3,
       2,     3,     4,     5,     6,     7,     2,     3,     4,     3,
       3,     3,     3,     2,     2,     3,     4,     5,     3,     4,
       3,     3,     3,     3,     5,     5,     4,     4,     6,     6,
       3,     3,     4,     3,     3,     3,     3,     3,     3,     2,
       3,     4,     1,     2,     3,     4,     5,     1,     2,     3,
       3,     4,     2,     2,     3,     5,     3,     4,     5,     4,
       4,     4,     7,     7,     8,     3,     9,     9,    10,     2,
       2,     2,     3,     5,     4,     1,     2,     4,     6,     5,
       1,     2,     4,     3,     2,     3,     2,     2,     2,     3,
       3,     4,     4,     5,     7,     5,     7,     4,     5,     3,
       4,     3,     4,     2,     3,     4,     5,     8,     2,     2,
      10,    12,     3,     4,     2,     4,     7,     9,     9,     8,
      11,    12,     2,     9,    10,     9,    10,     9,    10,     7,
       7,     7,     8,     8,     7,     7,     8,     8,     9,    10,
      10,    11,    12,    24,     5,     5,     2,     4,     5,     0,
       2,     4,     6,     8,     2,     3,     2,     3,     2,     2,
       3,     2,     2,     2,     2,     2,     1,     1,     1,     1,
       1
};

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
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
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     3,   113,   114,   115,   116,
     117,   102,   103,    59,   109,   110,    42,    43,    44,    39,
      40,    41,    38,    49,    54,    55,    56,    27,    28,    30,
      29,    35,    36,    31,    32,    68,    67,    75,   255,    34,
      61,    62,    63,    64,   112,    65,    66,    47,    48,    53,
      50,    51,    52,   129,    93,    94,    95,    92,     6,   127,
      72,    73,    69,    70,    71,    74,    85,    86,    81,    87,
      88,    89,    90,   104,   105,   106,   107,   108,    82,    83,
      96,    97,    98,    99,   100,   101,    21,    20,    24,    22,
      23,    25,    26,     7,    45,    46,     8,     9,    91,    14,
      10,   120,   121,   123,   122,   124,    33,   125,   126,   130,
       5,    12,    11,   128,    13,    16,    17,    18,    15,   648,
     647,    76,   131,    78,    84,    79,    80,    77,    37,    60,
      57,    58,   111,   118,   119,    19,   132,   133,     0,     0,
       0,     0,   926,     0,   626,     0,     0,   928,   930,   927,
     929,     0,     0,     0,     0,   266,     0,     0,     0,     0,
       0,     0,     0,     0,   448,     0,     0,     0,     0,   624,
     464,     0,     0,     0,     0,     0,   190,   388,   184,   281,
     869,     0,     0,     0,     0,     0,   589,     0,     0,   371,
     618,   847,   213,   366,   283,   168,     0,     0,     0,   802,
     807,     0,     0,     0,   256,     0,   620,     0,     0,   264,
       0,     0,   676,     0,     0,     0,     0,   385,   473,     0,
     813,     0,   580,     0,     0,   722,   724,     0,     0,   463,
       0,   216,   330,   756,     0,   136,     0,     0,     0,   766,
       0,     0,     0,     0,   182,   812,     0,     0,     0,     0,
       0,   325,   338,   509,   455,     0,   460,     0,   773,     0,
     469,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   207,   499,     0,   273,     0,     0,   640,     0,   279,
       0,     0,     0,     0,     0,     0,     0,   721,   178,   180,
       0,     0,     0,     0,   332,     0,     0,   848,   760,     0,
     713,     0,     0,     0,     0,     0,   228,     0,   229,     0,
     305,     0,     0,     0,     0,   642,   636,   743,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   390,   334,
       0,     0,     0,   868,   218,     0,   140,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   259,     0,   774,     0,     0,     0,   158,
     382,     0,     0,     0,   282,     0,     0,   323,   322,   494,
       0,     0,   571,   275,     0,     0,     0,     0,     0,     0,
     726,     0,   489,   161,   874,     0,   324,     0,     0,     0,
     882,   906,     0,     0,     0,     0,     0,   176,   757,     0,
     277,     0,   326,   339,     0,     0,     0,   586,   914,   919,
       1,     2,     4,     0,     0,     0,     0,     0,     0,   149,
     150,   147,   148,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   185,   186,   187,   188,   189,   191,
       0,   208,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   284,   285,   286,     0,     0,     0,   306,   307,   319,
       0,     0,     0,   363,     0,     0,   367,     0,     0,     0,
       0,     0,     0,     0,     0,   399,     0,   410,     0,   413,
       0,   416,     0,   419,     0,   427,   429,   431,   436,     0,
     439,     0,   445,     0,   449,     0,     0,     0,     0,     0,
       0,     0,   475,     0,   524,     0,   554,     0,     0,     0,
       0,   584,     0,     0,     0,   605,     0,   619,   621,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   831,   830,   829,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   863,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   225,   479,     0,   312,     0,     0,   167,   730,
       0,     0,     0,   426,   428,     0,   267,     0,     0,   422,
     424,     0,   425,     0,     0,   442,   443,     0,   746,     0,
     600,   276,     0,   816,     0,   793,   159,   160,   745,     0,
       0,     0,     0,   729,   727,     0,   740,   747,     0,     0,
     725,   795,   796,   794,     0,     0,   803,     0,   808,     0,
       0,   800,   257,   409,   398,   798,     0,   678,   679,     0,
     677,     0,   432,     0,     0,   474,   271,   720,   719,   723,
     601,   738,     0,   739,   695,   696,   358,   521,     0,     0,
     519,     0,   392,   393,     0,     0,     0,     0,   825,   412,
     843,     0,   845,     0,   418,   415,   510,     0,     0,   457,
     456,   459,   472,   487,     0,   470,     0,     0,     0,   362,
       0,     0,   272,   772,   771,     0,   500,     0,     0,   222,
     744,     0,     0,     0,   641,   523,     0,     0,   783,     0,
     782,     0,   781,   780,   718,   717,   790,   797,     0,     0,
     331,   258,   761,   604,     0,   261,   767,     0,   872,   230,
     231,     0,   465,   467,     0,   714,   704,   703,   759,   791,
       0,     0,     0,     0,     0,     0,   336,   333,   453,     0,
       0,   741,   716,   715,   224,   265,   684,   687,   685,   686,
     688,   689,   691,   690,   692,   682,   683,     0,     0,     0,
       0,     0,   770,   403,   404,     0,   707,     0,   262,   705,
       0,   263,   553,     0,     0,   495,   769,   775,     0,   223,
     227,   226,   742,   755,   814,     0,   254,     0,   490,     0,
       0,   778,   736,     0,     0,     0,     0,     0,   146,   555,
       0,     0,     0,   602,   603,   570,     0,   758,   278,   373,
     374,     0,   583,   378,   379,     0,     0,     0,     0,     0,
       0,     0,   152,     0,     0,     0,     0,     0,     0,     0,
     175,     0,     0,   173,   174,   172,   171,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   196,     0,   198,   200,
       0,   202,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   248,   249,   253,   251,   232,     0,   246,
       0,     0,     0,   302,   294,     0,     0,   304,     0,   287,
       0,   296,   293,     0,     0,     0,     0,     0,   308,     0,
     316,     0,   318,   320,     0,   370,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   437,     0,     0,     0,     0,     0,     0,     0,     0,
     478,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   639,   927,
       0,     0,     0,     0,     0,   650,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   849,     0,     0,
     861,     0,     0,     0,   850,     0,     0,     0,   859,     0,
       0,     0,   909,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   916,   918,
     915,   924,   921,   922,   925,   923,   920,   480,     0,     0,
     313,     0,     0,   731,     0,   732,     0,     0,   268,   269,
       0,   423,     0,     0,     0,     0,   751,   681,   817,     0,
     750,   752,     0,     0,   728,   753,   754,   699,   698,   804,
     809,   801,   810,   799,   680,     0,   433,   434,   840,     0,
     329,     0,   522,     0,   762,     0,   520,   394,     0,     0,
       0,     0,     0,   844,   846,     0,   512,     0,   511,     0,
       0,   516,     0,   488,     0,   821,   835,     0,     0,     0,
       0,   135,     0,     0,   502,     0,   501,     0,   504,     0,
     748,   749,   792,     0,   787,     0,     0,     0,   786,   693,
     694,     0,   768,     0,     0,   819,   820,   709,   711,   710,
     335,   337,   454,   911,   658,     0,     0,   675,     0,     0,
     652,     0,   660,     0,     0,     0,   405,     0,   708,   706,
     327,     0,     0,     0,   776,     0,     0,     0,     0,     0,
     875,   779,   737,     0,     0,     0,     0,     0,   556,     0,
     560,     0,     0,     0,   569,   375,   377,     0,   380,     0,
       0,     0,     0,     0,     0,   151,     0,     0,   162,   163,
     164,   165,   166,   169,   170,   177,   179,   181,   183,     0,
     195,   197,   199,   201,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   220,   221,   233,     0,   239,     0,   250,
     252,   245,   247,   274,   298,   300,     0,   299,   291,   297,
     289,   288,     0,     0,   292,     0,     0,   317,     0,     0,
     613,     0,     0,     0,   383,     0,   386,     0,     0,     0,
     401,     0,     0,     0,     0,   440,     0,   446,     0,     0,
     458,   485,     0,     0,     0,     0,   471,     0,     0,     0,
       0,     0,     0,     0,     0,   610,     0,     0,     0,     0,
       0,     0,   608,     0,     0,     0,   627,     0,     0,     0,
     637,     0,   643,     0,   649,   651,   654,   656,   653,   657,
     655,   697,   712,   700,   701,   702,     0,     0,   857,     0,
     862,   860,   851,   852,     0,     0,   864,     0,   873,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   907,     0,   917,     0,
       0,   315,     0,     0,   733,     0,   734,   430,   270,   421,
       0,     0,     0,   818,   280,     0,   805,   811,   674,   841,
       0,   328,     0,   763,     0,     0,     0,   832,     0,     0,
       0,     0,   513,     0,     0,   515,   615,   836,     0,   360,
       0,     0,     0,     0,   503,     0,   505,     0,     0,     0,
     785,   784,     0,   134,   345,   341,     0,     0,     0,   659,
     670,   671,     0,   669,     0,   663,     0,   661,   662,   260,
       0,     0,     0,     0,     0,   777,   815,     0,     0,     0,
     154,     0,     0,     0,   142,     0,   557,     0,     0,   561,
       0,   562,     0,   376,     0,     0,   137,     0,     0,   357,
     356,   153,   156,     0,   192,     0,     0,     0,   623,     0,
       0,     0,   214,   217,   219,     0,   235,     0,     0,   241,
       0,   303,   295,     0,     0,     0,     0,     0,     0,     0,
     612,     0,   384,   387,     0,     0,   400,   402,   411,   414,
     417,   420,   441,   447,     0,   486,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   552,     0,     0,   609,   611,
     559,   572,     0,     0,   581,     0,   585,   587,     0,   590,
       0,     0,   607,     0,     0,   632,     0,   633,     0,     0,
     638,     0,     0,   839,   855,     0,   858,   853,     0,     0,
     865,     0,   904,   905,   910,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   908,
       0,     0,     0,   314,   321,   735,     0,     0,     0,     0,
     806,     0,   617,     0,   764,     0,     0,     0,   395,     0,
       0,     0,     0,   834,     0,     0,     0,   359,   361,   365,
     838,     0,     0,   507,   789,   788,   346,     0,   348,   349,
     350,     0,   352,   353,   355,     0,   342,     0,   912,   672,
     668,     0,   664,     0,     0,     0,   406,     0,     0,   497,
       0,     0,   492,     0,     0,     0,   155,   558,   563,     0,
       0,     0,   381,   139,   138,   141,     0,     0,     0,     0,
       0,     0,     0,   215,   236,   234,   242,   240,   301,     0,
     340,     0,   309,     0,   364,     0,     0,   372,   389,   391,
     450,     0,   614,   462,     0,     0,     0,     0,     0,     0,
     541,     0,     0,     0,   537,     0,     0,     0,   578,   573,
       0,     0,     0,   598,   591,     0,   606,     0,   628,     0,
     629,     0,     0,     0,   644,     0,   645,     0,     0,   866,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   444,     0,   842,   616,   765,     0,   396,     0,   833,
     822,     0,   823,     0,   518,     0,   514,   837,   508,   506,
     347,   351,   354,   344,   343,     0,   673,   665,     0,     0,
     407,     0,     0,     0,     0,     0,   143,   144,   564,   565,
       0,   567,     0,   157,     0,     0,     0,   206,     0,     0,
       0,   237,     0,   243,     0,   311,   310,     0,   368,     0,
       0,     0,     0,     0,   476,     0,     0,     0,     0,   538,
       0,     0,     0,   574,     0,   582,   588,   592,     0,   625,
     631,     0,     0,   634,     0,   856,   854,     0,   876,     0,
       0,     0,     0,     0,   889,   890,     0,     0,   894,   895,
       0,     0,     0,     0,   891,     0,     0,     0,     0,     0,
     481,     0,     0,     0,   397,   824,     0,     0,     0,   517,
     913,   666,     0,   408,   496,     0,   491,     0,   145,   566,
     568,     0,   193,     0,     0,   209,     0,   238,   244,   290,
       0,   451,     0,     0,     0,     0,     0,     0,   549,     0,
     542,     0,     0,     0,     0,     0,   575,     0,     0,   593,
       0,     0,   635,   630,   646,     0,     0,     0,   879,     0,
       0,   892,     0,   896,     0,     0,   897,     0,     0,   893,
       0,   482,     0,   435,   438,     0,     0,   826,   827,   667,
     498,   493,   194,     0,     0,     0,   369,     0,   461,   466,
     468,   477,     0,   550,     0,   543,     0,     0,     0,     0,
       0,   579,     0,   599,     0,   867,   878,   877,     0,   883,
       0,   885,     0,     0,     0,     0,   898,   887,     0,     0,
     483,   870,   828,     0,     0,   210,     0,     0,     0,   551,
     544,     0,     0,     0,     0,   531,     0,   576,     0,   594,
       0,     0,   884,   886,     0,   899,   900,     0,   888,   484,
       0,     0,     0,     0,   452,   546,     0,     0,     0,     0,
       0,     0,   577,   595,     0,   880,     0,     0,   901,     0,
     871,     0,   622,   211,     0,   547,     0,   545,     0,     0,
       0,     0,     0,   881,     0,   902,   203,     0,   212,   548,
       0,   534,     0,   540,     0,   596,     0,     0,     0,     0,
     535,     0,     0,   597,     0,   204,     0,     0,   536,     0,
     532,     0,     0,   205,     0,     0,     0,     0,   539,   533,
       0,   525,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   526,   527,     0,     0,     0,   528,     0,
       0,     0,     0,     0,     0,     0,     0,   529,   530,   903
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   233,   234,   235,   236,   237,   238,   239,   240,   649,
     650,  1048,   651,   241,   242,   243,   244,   245,   246,   247,
     248,   249,   250,   674,  1664,   675,   676,   677,   678,  1096,
    1099,   251,   681,   252,   253,   254,   255,   256,   257,   258,
     259,   260,   261,   688,   262,   263,   264,   265,   266,   267,
     268,   701,  1122,  1126,   269,   707,   708,   270,   712,   271,
     272,   273,   274,   275,   276,   277,   278,   279,   280,   985,
     281,   282,   702,   283,  1615,  1815,   652,   284,   285,   286,
     713,   287,   288,   289,   290,  1059,  1060,   291,  1063,  1064,
     292,   293,   294,   295,   296,   902,   903,   297,   725,  1469,
     298,  1013,  1014,   299,   727,   300,   729,   301,   731,   302,
     733,   829,   830,   303,   304,   305,   306,   307,   308,   309,
     310,   738,   311,   740,   312,   313,   314,   742,   315,   744,
     316,   317,   318,   319,   320,   321,   322,   323,   812,  1480,
     922,   324,  1038,   325,  1025,   326,   936,   937,  1324,   327,
     916,   917,  1306,   328,   896,   329,   754,   330,   331,   332,
     333,   334,   335,   336,   761,   337,   338,   339,   340,   765,
     756,  1494,   813,  1481,   923,   897,   341,   342,   679,   343,
     344,   345,   346,   347,   348,  1194,   349,   350,   351,   875,
     352,   433,   353,   908,  1315,  1316,   354,  1287,  1288,   355,
     910,   356,   912,   357,   358,   797,   359,   360,   361,   362,
     363,   364,  1541,   365,   366,   805,   367,   811,  1470,  1317
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -1587
static const yytype_int16 yypact[] =
{
    4834,  -106,    32,     8,    14,    50,    -3,   263,   263,   263,
     243,   134,   148,   152,   167,    14,    14,   171,   175,   -62,
      14,    14,   181,   242,   264,    14,  -157,   187,   190,   229,
     305,   323,   363,   371,   382,   124,   263,   277,   263,   388,
     332,   340,   431,   456,   462,   463,   467,   486,   487,   404,
      14,    14,   325,   -31,   490,   493,   495,   512,   526,   -56,
      14,    14,   533,   189,   536,   465,   557,    26,   562,   573,
      14,    14,   577,   263,   578,   598,   617,   263,   620,   263,
     538,   636,   301,   315,   643,   646,   647,   650,   653,   656,
     663,   664,   678,   680,   683,  -115,   268,   687,   688,   690,
      14,   691,   356,   705,   715,    14,   206,   749,   751,   752,
     856,   756,   692,    14,   372,   759,   762,   -54,   -22,    27,
     769,   777,   778,    14,    14,    14,    14,   374,    14,   785,
     384,   788,   791,   797,    14,    14,   928,   399,   400,   798,
     801,    14,    14,   814,   817,   834,   843,    14,   -47,    14,
     263,    14,    14,   263,   263,    14,   845,   426,   961,   869,
     874,   877,   794,   892,    37,   894,   263,   263,    14,    14,
      14,   263,    14,    14,   812,    14,    14,    14,   900,   433,
     906,    14,   907,   263,   908,   909,   263,   263,   263,   914,
     915,   925,   931,   934,   936,   263,    14,   939,   954,  1059,
     965,   966,    14,    14,   263,   972,   852,   975,   977,   -69,
     980,  1032,   263,   984,   988,   990,    14,    14,   263,    14,
      14,   991,   -97,   995,   263,  1000,  1002,  1004,  1006,  1009,
    1010,  1011,  1030,  1231,  4438, -1587, -1587, -1587,  1113,    14,
     -19, -1587,  1001, -1587,   -59,    14,   263,   263,   263,   345,
      14,    14,    14,   263,  1137, -1587, -1587, -1587, -1587, -1587,
   -1587,    28, -1587,  1058, -1587, -1587, -1587, -1587,     5,    12,
     -14, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587,    14, -1587,   -84,    14,
   -1587, -1587,   -45,   -42,    14,    14, -1587,    14, -1587,    14,
      14,    14,    14, -1587, -1587,    14,    14,    14, -1587, -1587,
      14,    14, -1587, -1587,    14,    14, -1587,  1061,    14,    14,
      14,  1063, -1587,   -39, -1587, -1587, -1587, -1587, -1587,    14,
      14,    14, -1587, -1587,    14,    14,    14,    14,    14, -1587,
      14,    14,    14,    14,    14,   -29, -1587,    14,  1194,   140,
     494, -1587, -1587,     6,   263, -1587, -1587, -1587,     1, -1587,
   -1587,  1054,    14,  -125,    14, -1587,   -27,   870,    14,  1060,
    1253,  1257, -1587,  1067, -1587,   143,  -113, -1587, -1587, -1587,
   -1587,  1068,  1077,   263,   444, -1587,   263,    14,    14,   263,
     263,  1078,  1079,    14, -1587,   -96,  1081,  1083,  1017, -1587,
   -1587,   272,  1089,  1096,  1097,   -79, -1587, -1587, -1587,   263,
    1201,  1101,   454,  1102,   -68,  1104, -1587,  1105,  1106, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587,  1107,   263,    14,   263,
    1230,   263,  1071,    -5, -1587,  1114, -1587,    14,    14, -1587,
     263,  1116, -1587,  1117,   253,   470,  1119, -1587, -1587,  1120,
   -1587,  1121, -1587,  1124,  1125, -1587, -1587,  1127,  1128, -1587,
    1130, -1587,   263, -1587,  1132, -1587,  1141,  1143,    14, -1587,
     263,    14,  1144,    14, -1587,   507,    14,   263,   263,    14,
      14, -1587, -1587,    14,    14,  1145, -1587,  1147, -1587,    14,
      14,  1148,   263,    14,  1149,   263,    14,  1150,  1151,  1153,
     263, -1587,    14,  1155, -1587,   125,   263, -1587,  1156, -1587,
      14,   544,   -30,  1157,  1158,  1159,  1161, -1587, -1587, -1587,
    1162,  1163,    14,   263, -1587,  1165,  1167, -1587, -1587,  1168,
   -1587,    14,    14,  1169,   278,  1170, -1587,  1171, -1587,  1173,
   -1587,    14,  1174,  1179,    14, -1587, -1587, -1587,  1180,  1181,
    1182,  1193,  1196,  1197,  1198,   263,   263,    14, -1587,    14,
    1199,   474,  1110, -1587, -1587,  1200, -1587,  1202,  1203,    14,
    1204,  1205,  1206,  1207,  1208,  1209,  1210,  1211,  1212,  1213,
    1214,  1215,   -26, -1587,   263, -1587,  1216,    14,   282, -1587,
   -1587,  1219,   287,  1223, -1587,    14,   263, -1587, -1587,  1327,
    1224,   292, -1587, -1587,  1225,    14,    14,  1226,  1227,   294,
   -1587,  1228,  1346, -1587, -1587,    14, -1587,   130,   300,   -58,
   -1587, -1587,  -108,    14,  1229,  1232,   482, -1587, -1587,  1233,
   -1587,  1234, -1587, -1587,    14,    14,    14, -1587, -1587, -1587,
   -1587, -1587, -1587,   -28,  1140,   542,   263,   302,  1172, -1587,
   -1587, -1587, -1587,  1334,  1338,  1339,  1340,  1342,  1240,  1344,
      14,  1242,  1243,  1244,  1247,   263,   263,   263,   263,    14,
      14,    14,    14,    14, -1587,    14,    14,    14,    14, -1587,
     -38, -1587,   263,    14,   263,    30,   377,   489,    25,    14,
     263,   375,   263,  1160,  1249,   263,  1252,  1254,   -21,   263,
    1166, -1587, -1587, -1587,   263,  1255,   263, -1587, -1587, -1587,
    1258,  1356,   518, -1587,   263,    14, -1587,  -124,  1396,    14,
     263,    14,   263,   263,   263, -1587,    14, -1587,    14, -1587,
      14, -1587,    14, -1587,    14, -1587, -1587, -1587, -1587,  1259,
   -1587,    14, -1587,    14, -1587,   263,  1262,   263,   263,   263,
    1263,    14, -1587,  -121, -1587,   -13, -1587,    14,    14,    14,
      14, -1587,    14,    14,    14, -1587,   263, -1587, -1587,    14,
      14,    14,  1175,    -9,    14,    14,    14,    14,    14,   263,
      14,    14,   240, -1587, -1587, -1587,   263,  1264,   263,    14,
     383,   263,    14,  1265,   263,   561,  1266, -1587,  1370,    14,
    1371,   -81,    14,  1372,    14,  1270,    14,  1374,  1375,    14,
     263,  1273,    14, -1587,   -73, -1587,   393,   263, -1587, -1587,
    1274,    14,  -138, -1587, -1587,   263, -1587,  1279,   574, -1587,
      14,   263,    14,   263,   263, -1587, -1587,  1280, -1587,  1284,
   -1587, -1587,  1287, -1587,   403, -1587, -1587, -1587, -1587,  1288,
    1289,    14,  1290, -1587, -1587,  1292, -1587, -1587,  1294,  1298,
   -1587, -1587, -1587, -1587,  1299,  1301,   263,   263, -1587,   191,
      14, -1587, -1587, -1587, -1587, -1587,  1302, -1587, -1587,  1303,
   -1587,    14, -1587,  1305,   263, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587,   410, -1587, -1587, -1587,    14, -1587,    14,   582,
      14,    14,    14, -1587,    14,   263,   263,   263,   263, -1587,
   -1587,  1307, -1587,  1308, -1587, -1587,    14,    14,   112,    14,
   -1587, -1587,    14, -1587,   263,    14,   263,   263,   263, -1587,
     263,  1309, -1587, -1587, -1587,   263,    14,    14,   412, -1587,
   -1587,  1310,  1311,  1312, -1587, -1587,   366,    14, -1587,    14,
   -1587,   661, -1587, -1587, -1587, -1587, -1587, -1587,  1313,  1314,
   -1587, -1587, -1587, -1587,    14, -1587, -1587,  1315, -1587, -1587,
   -1587,    14, -1587, -1587,    14, -1587, -1587, -1587, -1587, -1587,
     263,   263,  1316,  1317,  1318,   602, -1587, -1587, -1587,  1319,
    1320, -1587, -1587, -1587,    14, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587,   607,   -49,   622,
     -11,   263, -1587,    14, -1587,    14, -1587,  1321, -1587, -1587,
    1322, -1587, -1587,  1323,    14,  1184, -1587, -1587,   418, -1587,
      14,    14, -1587, -1587, -1587,   263, -1587,    14,  1188,  1324,
    1325, -1587, -1587,  1326,   263,   263,   263,   263,   263, -1587,
    1328,   263,  -139, -1587, -1587, -1587,  1329, -1587, -1587,    14,
   -1587,   420, -1587,    14, -1587,    14,    14,    14,   263,  1330,
     263,  1331, -1587,   263,    14,  1333,  1349,  1350,  1351,  1352,
   -1587,  1353,  1354, -1587, -1587, -1587, -1587,  1355,  1357,  1358,
    1359,  1360,  1361,  1362,  1366,  1367, -1587,   263, -1587, -1587,
      14, -1587,    14,   263,   263,  1175,   263,   -18,  1368,    14,
      14,    14,    14, -1587,    14, -1587,    14, -1587,    14, -1587,
     263,  1376,  1377,   263, -1587,   263,  1379, -1587,  1380, -1587,
    1382, -1587, -1587,  1383,   432,   263,  1384,   263, -1587,   263,
   -1587,  1385, -1587, -1587,   263, -1587,    14,    14,   448,    14,
     263,  1388,   263,  1390,   263,   263,    14,    14,    14,    14,
      14, -1587,   669,   676,   263,    14,   263,   263,   263,    14,
   -1587,    14,    14,    14,   263,   263,   263,   263,    14,    14,
      14,    14,    14,    14,    14,   449,   263,   -60, -1587,  1176,
     263,  1392,    14,   450,   737, -1587,  1394,  1395,  1397,  1399,
    1400,  1401,  1402,  1403,  1404,  1405,   263, -1587,   263,   775,
   -1587,  1406,  1407,  1408, -1587,   451,  1272,   779, -1587,  1409,
     263,  1459, -1587,   263,   263,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   263,   263,   263,   780,    14, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,   -71,    14,
   -1587,   457,   263, -1587,  -135, -1587,  1410,  1411, -1587, -1587,
    1412, -1587,  1413,  1414,  1415,   263, -1587, -1587, -1587,  1416,
   -1587, -1587,  1417,    14, -1587, -1587, -1587, -1587, -1587,   263,
   -1587, -1587,   263, -1587, -1587,  1418, -1587,   263, -1587,   263,
   -1587,  1419, -1587,   263, -1587,   466,    14, -1587,  1175,   469,
     263,   263,    14, -1587, -1587,    14, -1587,   251, -1587,    14,
     263, -1587,  1420, -1587,  1421,   263, -1587,   263,   476,   263,
     263, -1587,   263,    14, -1587,   506, -1587,    14, -1587,   -43,
   -1587, -1587, -1587,    14, -1587,  1422,  1424,    14, -1587, -1587,
   -1587,  1426, -1587,   795,    14,   263,   263, -1587, -1587, -1587,
   -1587, -1587, -1587,  1373, -1587,  1428,  1429, -1587,  1431,   265,
   -1587,   551, -1587,  1432,  1433,  1434, -1587,  1175, -1587, -1587,
   -1587,  1435,    14,   263, -1587,  1436,  1437,  1438,    14,   263,
   -1587, -1587, -1587,   556,   263,   263,  1439,    14, -1587,  -110,
   -1587,   263,  -160,  -147, -1587, -1587, -1587,  1440, -1587,    14,
      14,   799,   263,    14,  1442, -1587,  1444,   803, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,    14,
   -1587, -1587, -1587, -1587,   263,   263,  1446,   263,   263,   263,
    1447,  1448,  1449, -1587, -1587,   -25,  1474,   406,  1475, -1587,
   -1587, -1587, -1587, -1587, -1587,   263,  1450, -1587, -1587, -1587,
   -1587, -1587,   263,   263, -1587,   263,    14, -1587,   263,   -67,
   -1587,  1451,    14,  1452, -1587,  1453, -1587,   263,   263,   824,
   -1587,   829,   849,   873,   885, -1587,  1454, -1587,  1455,    14,
      14, -1587,   263,   263,   263,   263,    14,   -37,   263,   263,
     263,  1456,   263,   263,   888, -1587,   889,   697,   592,   895,
     605,   734, -1587,   633,   263,  1457, -1587,   263,   -48,  1458,
   -1587,   263, -1587,   263, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587,  1461,   905, -1587,  1464,
   -1587, -1587, -1587, -1587,   638,    14, -1587,  1465, -1587,   263,
    1466,   639,   263,   263,   263,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   263,   263, -1587,  1467, -1587,    14,
     -33, -1587,  1468,  1469, -1587,  1470, -1587, -1587, -1587, -1587,
     263,   263,   263, -1587, -1587,  1471,   263, -1587, -1587, -1587,
     263, -1587,   641, -1587,   910,  1175,  1472, -1587,  1175,    14,
      14,  1473, -1587,   263,   263, -1587, -1587, -1587,   263, -1587,
    1476,  1477,  1478,  1480, -1587,   263, -1587,    14,   228,  1481,
   -1587, -1587,  1483, -1587,  1164, -1587,  1484,    14,  1485, -1587,
   -1587, -1587,  1486, -1587,   917, -1587,   924, -1587, -1587, -1587,
    1175,  1487,   263,  1488,  1489, -1587, -1587,   263,  1490,  1491,
   -1587,    14,    14,    14, -1587,  1492, -1587,  1493,   642, -1587,
     263, -1587,   263, -1587,   933,  1494, -1587,  1495,  1496,    14,
   -1587, -1587, -1587,    14,    14,    14,   263,   263, -1587,   263,
     263,  1498, -1587, -1587, -1587,    14, -1587,    14,    14, -1587,
      14,   263, -1587,    14,  1499,    14,   648,  1503,    14,   263,
   -1587,  1508, -1587, -1587,  1511,  1512, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587,   657, -1587,  1513,   667,   263,   263,
      14,   263,   -35,   263,   672, -1587,   263,   263, -1587, -1587,
   -1587, -1587,   254,   674, -1587,   263, -1587, -1587,   263, -1587,
     257,   682, -1587,  1514,   263, -1587,  1515, -1587,   263,   -50,
   -1587,   704,  1517, -1587, -1587,    14, -1587, -1587,   263,  1519,
   -1587,   263, -1587, -1587, -1587,   263,   263,   263,   263,   263,
     263,   263,   263,   263,   263,   263,   263,   263,   263, -1587,
     -32,    14,   263, -1587, -1587, -1587,   263,   263,  1520,  1578,
   -1587,  1522, -1587,  1523, -1587,  1525,   263,  1526, -1587,    14,
    1528,   707,   708, -1587,   745,  1529,  1530, -1587, -1587, -1587,
   -1587,  1532,  1533, -1587, -1587, -1587, -1587,    14, -1587, -1587,
   -1587,   263, -1587,   263, -1587,  1484, -1587,  1484,  1482, -1587,
   -1587,  1535, -1587,   946,   263,  1537, -1587,   263,   263, -1587,
     263,   263, -1587,    14,  1538,  1539, -1587, -1587, -1587,  1540,
     748,   772, -1587, -1587, -1587, -1587,  1542,    14,   263,   263,
    1543,   263,   263, -1587, -1587,   440, -1587,   505, -1587,    14,
   -1587,  1545, -1587,  1546, -1587,    14,  1547, -1587, -1587, -1587,
   -1587,   263, -1587, -1587,   263,   263,   263,    14,  1562,   263,
   -1587,   263,   263,   263, -1587,   782,   263,   263, -1587, -1587,
     789,  1563,  1564, -1587, -1587,   807, -1587,  1567, -1587,  1568,
   -1587,   263,   263,  1571, -1587,   263, -1587,  1572,  1573,  1544,
     -46,   263,   263,  1575,  1576,   263,   263,  1577,   808,   263,
     263,   263,  1580,   263,   263,    14,   263,    14,  1581,   263,
     263, -1587,  1582, -1587, -1587, -1587,  1583, -1587,   810, -1587,
   -1587,   263, -1587,   263, -1587,  1586, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587,  1588, -1587, -1587,   974,  1589,
   -1587,  1590,   263,  1591,   263,  1592, -1587, -1587, -1587, -1587,
    1593, -1587,  1594, -1587,   263,  1595,   263, -1587,   263,   816,
      14, -1587,    14, -1587,  1596, -1587, -1587,   263, -1587,   823,
     263,   263,   263,   263, -1587,   263,   826,   847,   263, -1587,
     263,   263,   263, -1587,   846, -1587, -1587, -1587,   865, -1587,
   -1587,  1597,  1598, -1587,  1599, -1587, -1587,    14, -1587,    14,
     263,   276,   263,   263, -1587, -1587,  1600,   263, -1587, -1587,
    1601,   263,   263,   853, -1587,   263,  1602,    14,  1603,   263,
   -1587,  1605,  1607,  1681, -1587, -1587,   263,  1608,  1609, -1587,
   -1587, -1587,  1610, -1587, -1587,  1611, -1587,  1612, -1587, -1587,
   -1587,  1613, -1587,   263,   263, -1587,   263, -1587, -1587, -1587,
    1614, -1587,   263,  1615,  1616,  1618,  1619,   263, -1587,   855,
   -1587,   872,   263,   263,   263,   263, -1587,   317,   263, -1587,
     342,   263, -1587, -1587, -1587,  1622,  1625,   -34, -1587,   875,
     901, -1587,   263, -1587,   263,   263, -1587,  1626,   904, -1587,
     263, -1587,  1627, -1587, -1587,  1628,  1629, -1587, -1587, -1587,
   -1587, -1587, -1587,   263,   263,   920, -1587,   263, -1587, -1587,
   -1587, -1587,   263, -1587,  1630, -1587,   923,   263,   263,   263,
     929, -1587,   964, -1587,   993, -1587, -1587, -1587,    14, -1587,
    1631, -1587,  1632,   263,  1635,   997, -1587, -1587,  1636,  1637,
   -1587,  1741, -1587,   263,   263, -1587,   263,  1639,  1005, -1587,
   -1587,    14,   263,   263,   263, -1587,   263, -1587,  1640, -1587,
    1040,   338, -1587, -1587,   263, -1587, -1587,  1049, -1587, -1587,
    1643,   263,  1646,  1056, -1587, -1587,  1090,  1647,   263,   263,
     263,   263, -1587, -1587,   263, -1587,   344,   263, -1587,  1648,
   -1587,  1095, -1587, -1587,  1649, -1587,  1650, -1587,   263,  1099,
    1651,   263,  1108, -1587,   263, -1587, -1587,   263, -1587, -1587,
      14, -1587,  1126, -1587,   263, -1587,  1652,   263,  1129,    14,
   -1587,  1653,   227, -1587,   263, -1587,  1654,    14, -1587,   263,
   -1587,   263,   263, -1587,  1655,  1656,    44,   263, -1587, -1587,
     263, -1587,    14,   263,   263,   263,   263,   263,   263,  1657,
    1658,   263,   263, -1587, -1587,    -6,   263,   263, -1587,    14,
     263,   263,   263,   263,  1659,  1660,  1661, -1587, -1587, -1587
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1587, -1587, -1587,  1662, -1587, -1587, -1587, -1587, -1587,  1221,
   -1587, -1587,  1509, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,  1187,
     918, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587,  1178, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1586, -1587,  -167, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587,   806, -1587, -1587,   818,
   -1587, -1587, -1587, -1587, -1587, -1587,   968, -1587,  -219, -1100,
   -1587, -1587,   858, -1587,  1445, -1587,  -228, -1587,  1386, -1587,
    -236,  -471,  1497, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,    74, -1151,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587,   941,  -926, -1587,
   -1587,   970,  -910, -1587,  -466, -1587,  1381, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587,  1267, -1587, -1587, -1587, -1587,  1369,
    1295,   710,  -284, -1431,   962,  -890, -1587, -1587,  -200, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587,  -415, -1587, -1587,  -490,  -938, -1587, -1587,   600, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,
   -1587, -1587, -1587, -1587, -1587, -1587, -1587, -1587,  2247,    -7
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -927
static const yytype_int16 yytable[] =
{
     381,   382,   383,   386,   716,   900,  1292,  1308,   375,   787,
    1292,  1326,   658,   441,   659,   709,  1222,   869,  1486,   690,
    1356,   429,   372,   372,   429,   377,   372,   710,   377,   413,
    1816,   415,   430,   431,   372,   430,   431,  1113,   369,   752,
    1649,  1007,   691,   645,  1172,  1044,  1045,  1046,   378,  1705,
     377,   378,  1879,  1651,   377,  1705,  2277,  1471,  1472,  1473,
    1474,  1390,  1255,   372,  1008,  1564,   460,   645,   645,   788,
     464,   737,   466,   378,   372,  1173,   372,   378,   736,  1432,
     372,  1114,  1499,  1174,  1505,   484,   735,   372,   377,   487,
    1646,   372,  1049,   372,   368,  1675,   789,   372,   372,   372,
     372,   703,   372,   628,   838,   372,  2260,   377,   372,   379,
     372,   378,   379,   790,   372,   372,   372,   372,   372,   372,
     380,   848,   372,   380,  1009,   372,   372,  1107,   511,   692,
     378,   614,   857,   402,   379,   411,   704,   870,   379,   839,
    1506,   767,   768,   552,   442,   380,   555,   556,   693,   380,
    1900,  1357,  1737,  1647,  2018,  1050,   849,  1606,   949,   570,
     571,   372,   513,  1310,   575,  1880,  2147,   858,   791,   434,
     950,   412,   379,   372,   372,   377,   588,   372,   792,   591,
     592,   593,  1433,   380,   372,   694,   785,  1051,   600,  1362,
     549,   379,   705,   629,  2278,   871,  1147,   609,   378,  1171,
     695,  1115,   380,   793,   696,   618,   395,   706,   371,   794,
     850,   624,   715,  1256,  1650,   697,  1565,   631,   485,   874,
     646,   859,  1116,   647,  1363,  1117,  1688,  1652,  1705,  1953,
    1118,  1954,   370,   821,   800,  1391,  1392,  1393,   648,   666,
     667,   668,  1010,   914,  2261,   646,   684,  1249,   909,  1559,
     374,   719,  1659,   795,   721,  1507,  1710,   751,  2279,   379,
    1771,  1925,  1881,   648,  1901,  1902,  1738,   771,  1066,  1133,
     380,  1040,   615,  1223,  1224,  1225,  1226,  1227,  1102,  1228,
    1229,  1230,  1231,  1232,  1175,  1233,  1234,  1235,  1236,  2249,
     377,  1047,   453,   515,   443,  2019,   698,  2020,  1607,  1654,
    1188,  1358,  1310,   567,  1176,   699,   377,  2148,  2262,  1119,
     686,  1192,  1311,   378,   377,  1364,   377,   775,   776,   660,
     661,   662,   663,   664,   711,   940,   377,   796,   377,   378,
    1041,   377,  1203,   870,   387,   377,   372,   378,   700,   378,
     384,   377,   376,   819,   879,   377,   687,   786,   388,   378,
     377,   378,   389,   372,   378,   377,  1622,   377,   378,  1261,
     941,  1261,   498,   377,   378,   377,   803,   390,   378,   777,
     870,   393,   804,   378,   379,   394,   825,  1597,   378,   828,
     378,   398,   833,   834,   669,   380,   378,   372,   378,   448,
     379,  1281,   429,   820,   844,  1592,   870,  1604,  1189,   870,
     379,   380,   851,   430,   431,   429,  1292,  1597,  1597,   380,
     379,   380,   379,  1204,   942,   379,   430,   431,   870,   379,
     864,   380,   866,   380,   868,   379,   380,  2250,  1803,   379,
     380,  1145,   778,   876,   379,  1296,   380,   881,   377,   379,
     380,   379,   399,   385,   779,   380,   377,   379,   372,   379,
     380,  1311,   380,   880,  1888,   892,   377,  1893,   380,   870,
     380,   378,   372,   899,   400,  1623,   377,  1170,   486,   378,
     911,   913,   843,   377,  1109,   377,  2098,   403,   966,   378,
     404,   377,  1016,   377,   870,   927,   870,  1019,   930,   378,
    1345,  1346,  1027,   935,  1034,   377,   378,   405,   378,   943,
    1042,   469,  1072,   372,   378,   406,   378,   670,   671,   672,
     673,   377,   377,   377,   377,   471,   959,  2141,   378,   372,
     377,   372,   379,   407,   372,   432,  1678,   967,  1247,   377,
     379,   372,   377,   380,   378,   378,   378,   378,  2205,   377,
     379,   380,  2143,   378,  2223,   414,   372,   372,   982,   983,
     379,   380,   378,   372,  1333,   378,   493,   379,   919,   379,
    1980,   380,   378,   408,   925,   379,  1334,   379,   380,   377,
     380,   409,   507,   372,   524,  1124,   380,  1011,   380,   379,
     372,  1017,   410,  1210,   528,  1020,  1111,   372,   416,  1023,
     380,   372,   378,  1250,  1028,   379,   379,   379,   379,   536,
     538,   372,  1035,  1268,   379,   377,   380,   380,   380,   380,
    1290,  1043,  1328,   379,   377,   380,   379,   372,  1374,   377,
    1396,   372,   417,   379,   380,  1982,   559,   380,   378,   372,
     418,   419,  1451,   583,   380,  1247,   372,   378,  1070,  1071,
    1073,  1247,   378,   994,   826,   905,   906,   907,  1460,  1502,
    1512,  1533,   372,  1189,   854,   377,   420,  1561,  1087,  1088,
    1089,  1090,   421,   422,   380,   372,  1583,   423,   377,  1587,
     882,   780,   781,  1103,   988,  1104,  1599,  1106,   378,  1030,
    1031,  1120,  1055,  1123,  1125,  1127,   424,   425,  1123,   379,
     436,   378,  1135,   437,   426,   438,   377,  1137,   379,  1139,
     380,   377,   377,   379,   377,   377,  1328,  1144,   372,   380,
    1247,   377,   439,  1151,   380,  1153,  1154,  1155,  1142,   378,
     377,   372,  1069,   782,   378,   378,   440,   378,   378,   372,
     377,   946,   947,   447,   378,   377,   450,   377,  1164,   379,
    1166,  1167,  1168,   378,   948,   377,  1247,  1247,  1177,   372,
     380,  1625,   379,   378,   372,   451,  1640,   452,   378,  1185,
     378,  1216,   455,   380,   429,  1191,  1193,   377,   378,   372,
     377,   377,  1200,   456,  1259,   430,   431,   459,   461,  1206,
     379,  1208,  1294,  1211,  1212,   379,   379,  1215,   379,   379,
     378,   380,  1724,   378,   378,   379,   380,   380,   462,   380,
     380,   429,  1350,  1245,   379,  1727,   380,  1354,   377,  1251,
    1252,   377,   430,   431,   379,   380,   372,   463,  1257,   379,
     465,   379,  1360,   372,  1262,   380,  1263,  1264,   467,   379,
     380,   378,   380,  1732,   378,   377,   468,  1269,  1747,  1753,
     380,  1782,  1838,   473,   372,   377,   474,   475,  1862,  1337,
     476,   379,   377,   477,   379,   379,   478,  1870,   378,  1279,
    1280,  1338,   380,   479,   480,   380,   380,  1873,   378,  1475,
     377,   377,  1884,   377,  1889,   378,  1477,  1289,   481,   377,
     482,   372,  1894,   483,   372,  1291,   377,   488,   489,   377,
     490,   492,   379,   378,   378,   379,   378,  1721,  1299,  1300,
    1301,  1302,   378,   380,  1904,   495,   380,  1940,  1942,   378,
     377,  1312,   378,   429,  1608,   496,   377,  1314,   377,   379,
    1318,  1319,   372,  1320,   430,   431,   372,   372,  1322,   379,
     380,  1329,   429,   378,  1729,   377,   379,  1514,   377,   378,
     380,   378,   372,   430,   431,  1944,   372,   380,  1969,   500,
     372,   501,   502,   503,   379,   379,   504,   379,   378,   509,
     505,   378,   510,   379,   377,   380,   380,   377,   380,   517,
     379,   372,  1971,   379,   380,  1528,   372,   518,   519,  1536,
    1556,   380,  1999,   377,   380,   527,   377,   378,   530,  2003,
     378,   531,   377,   372,   379,  1614,   372,   532,   540,  1656,
     379,   541,   379,  1662,  1365,   380,   378,  2007,  2029,   378,
    2045,   380,   372,   380,   544,   378,  2065,   545,  1373,   379,
     372,  1375,   379,  2071,  1696,   535,  2078,   377,  1376,  1698,
     380,  1379,   372,   380,   546,   372,   372,  1383,  1384,  1385,
    1386,  1387,   372,   547,  1389,   558,  2086,  2080,   379,  1699,
     378,   379,   372,  2106,  1397,  2133,   377,   372,   561,   380,
     377,  1402,   380,  1404,   372,  2089,  1406,   379,   377,   562,
     379,   372,  2135,  1700,   563,  2149,   379,   564,   380,   378,
     372,   380,  1722,   378,   565,  1701,  1730,   380,  1718,  1720,
    1424,   378,   566,   372,   569,  1726,  1427,  1428,  1430,  1431,
     582,  2151,   578,   377,  2157,  1744,   585,   587,   589,   590,
    1784,   379,   377,  1442,   594,   595,  1445,  1820,  1446,   377,
    2165,   372,   380,  2170,  1822,   596,   378,  1452,  1453,  2175,
    1455,   597,  1456,  1842,   598,   378,   599,  1458,   429,   602,
     379,  1461,   378,  1463,   379,  1465,  1957,  1467,  1468,   430,
     431,   380,   379,   377,   603,   380,   604,  1479,   377,  1483,
    1484,  1485,   377,   380,  2177,   605,   606,  1490,  1491,  1492,
    1493,   377,   610,   611,  2051,   612,   378,   613,  1503,  1504,
     616,   378,   617,  1509,   619,   378,  1513,   379,   620,   377,
     621,   627,   377,  2179,   378,   630,   379,  2186,   380,  1526,
     632,  1527,   633,   379,   634,  2195,   635,   380,  1534,   636,
     637,   638,   378,  1539,   380,   378,  1542,  1543,  1544,  1545,
    1546,  1547,  1548,  1549,  1550,  1551,  1552,  1553,  1554,  1555,
     639,   640,  1372,   643,   685,   689,  1378,   379,   377,  -926,
    2203,   746,   379,   750,  1562,  1563,   379,   377,   380,  2208,
     806,   377,   774,   380,   798,   379,  2213,   380,  1572,   816,
     815,   378,  -926,   817,   807,   808,   380,   818,   823,   809,
     378,   810,  1576,   379,   378,  1577,   379,   824,   835,   836,
    1289,   840,  1580,   841,   380,   842,  1582,   380,  1584,   845,
    2215,  1586,  1588,  1589,  1590,  2226,   846,   847,   852,  2231,
    1312,   853,   856,  1594,   860,   861,   862,   863,  2235,   867,
    1598,  1600,  1601,  1602,   872,  1603,   877,   878,  1329,   884,
     885,   886,  1189,  -926,   887,   888,  2240,   889,   890,  2245,
     891,   379,   893,   380,  -926,   379,   653,   654,   655,   656,
     657,   894,   380,   895,   901,   920,   380,   921,   926,   929,
     932,   933,  1624,   934,  1626,   939,   944,   952,   953,   954,
    1631,   955,   956,   957,   990,   960,  1634,   961,   962,   965,
     968,   969,  1639,   970,   972,  1024,  1641,  1642,  1643,   973,
     975,   976,   977,  1806,  1648,  1807,  1808,  1809,  1810,  1811,
    1812,  1813,  1814,   978,  1037,  1658,   979,   980,   981,   987,
     991,  1068,   992,   993,   995,   996,   997,   998,   999,  1000,
    1001,  1002,  1003,  1004,  1005,  1006,  1012,  1666,  1667,  1018,
    1669,  1670,  1671,  1021,  1026,  1029,  1032,  1033,  1036,  1053,
    1074,  1075,  1054,  1057,  1058,  1076,  1077,  1078,  1681,  1079,
    1080,  1081,  1083,  1084,  1085,  1683,  1684,  1086,  1685,  1129,
    1128,  1687,  1131,  1141,  1132,  1138,  1136,  1149,  1140,  1161,
    1694,  1695,  1165,  1169,  1207,  1214,  1218,  1219,  1221,  1238,
    1240,  1242,  1243,  1246,  1253,  1706,  1707,  1708,  1709,  1258,
    1265,  1712,  1713,  1714,  1266,  1716,  1717,  1267,  1270,  1271,
    1273,  1725,  1274,  1728,  1275,  2021,  1733,  1734,  1276,  1277,
    1736,  1278,  1283,  1284,  1741,  1286,  1742,  1303,  1304,  1321,
    1330,  1331,  1332,  1339,  1340,  1342,  1347,  1348,  1349,  1352,
    1353,  1368,  1369,  1370,  1380,  1381,  1382,  1748,  1388,  1394,
    1403,  1405,  1751,  1408,  1754,  1755,  1756,  1757,  1758,  1759,
    1760,  1761,  1762,  1763,  1764,  1765,  1766,  1767,  1768,  1409,
    1410,  1411,  1412,  1413,  1414,  1415,  1540,  1416,  1417,  1418,
    1419,  1420,  1421,  1776,  1777,  1778,  1422,  1423,  1434,  1780,
    1535,  1677,  1680,  1781,   783,  1783,  1443,  1444,  1787,  1447,
    1448,  1790,  1449,  1450,  1454,  1457,  1794,  1795,  1464,  2087,
    1466,  1796,  1510,  2090,  1516,  1517,  1101,  1518,  1801,  1519,
    1520,  1521,  1522,  1523,  1524,  1525,  1530,  1531,  1532,  1538,
    1566,  1567,  1568,  1569,  1570,  1571,  1573,  1574,  1578,  1581,
    1595,  1596,  1610,  1825,  1611,  1827,  1613,  1618,  1619,  1620,
    1830,  1621,  1627,  1628,  1629,  1632,  1635,  1636,  1637,  1644,
    1653,  1839,  1660,  1840,  1661,  1841,  1668,  1672,  1673,  1674,
    1682,  1690,  1692,  1693,  1702,  1703,  1715,  1735,  1740,  1849,
    1850,  1743,  1851,  1852,  1746,  1750,  1752,  1769,  1773,  1774,
    1775,  1779,  1788,  1793,  1858,  1932,  1797,  1798,  1799,  1863,
    1800,  1804,  1866,  1805,  1614,  1818,  1819,  1826,  1828,  1829,
    1831,  1832,  1836,  1837,  1843,  1844,  1845,  1871,  1853,  1860,
    1874,  1875,  1876,  1864,  1878,  1882,  1883,  1885,  1867,  1886,
    1887,  1868,  1869,  1872,  1896,  1898,  1890,  1906,  1891,  1909,
    1931,  1892,  1933,  1934,  1895,  1935,  1937,  1897,  1939,  1946,
    1947,  1899,  1948,  1949,  1905,  1956,  1955,  1960,  1966,  1967,
    1968,  1908,  1973,  1977,  1910,  1985,  1986,  1988,  1911,  1912,
    1913,  1914,  1915,  1916,  1917,  1918,  1919,  1920,  1921,  1922,
    1923,  1924,  1994,  2005,  2006,  1928,  2206,  2009,  2010,  1929,
    1930,  2013,  2015,  2016,  2017,  2024,  2025,  2028,  2115,  1936,
    2034,  2040,  2043,  2044,  1941,  1943,  2049,  1945,  2050,  2053,
    2054,  2056,  2058,  2059,  2060,  2062,  2069,  2092,  2093,  2094,
    2101,  2103,  2109,  2111,  1951,  2113,  1952,  2114,  2117,  2118,
    2119,  2120,  2121,  2122,  2126,  2128,  2129,  1959,  2130,  2131,
    1961,  1962,  2145,  1963,  1964,  2146,  2156,  2160,  2161,  2162,
    2169,  2182,  2183,  1970,  1972,  2185,  2188,  2189,  2190,  2194,
    2202,  1975,  1976,  2210,  1978,  1979,  2212,  2217,  2225,  2228,
    2229,  2233,  2243,  2248,  2253,  2258,  2259,  2273,  2274,  2287,
    2288,  2289,   784,  1098,  1989,  1395,   915,  1990,  1991,  1992,
    1297,  1366,  1995,  1130,  1996,  1997,  1998,  1323,  2000,  2001,
    2002,  1398,   873,  2004,  1313,   832,  1305,  1579,  2008,  1496,
    1022,   945,     0,     0,  2011,  2012,   642,     0,  2014,     0,
     963,     0,  1062,     0,  2022,  2023,     0,     0,  2026,  2027,
       0,  2030,  2031,  2032,  2033,     0,  2035,  2036,     0,  2038,
       0,     0,  2041,  2042,     0,     0,     0,     0,     0,     0,
       0,  2046,     0,     0,  2047,     0,  2048,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2055,     0,  2057,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2061,     0,  2063,
       0,  2064,  2066,     0,     0,     0,     0,     0,     0,     0,
    2070,     0,  2072,  2073,  2074,  2075,  2076,     0,  2077,  2079,
    2081,  2082,     0,  2083,  2084,  2085,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2097,     0,  2099,  2100,     0,     0,     0,
    2102,     0,     0,     0,  2104,  2105,  2107,     0,  2108,     0,
       0,     0,  2112,     0,     0,     0,     0,     0,     0,  2116,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2123,  2124,     0,  2125,
       0,     0,     0,     0,     0,  2127,     0,     0,     0,     0,
    2132,     0,  2134,     0,  2136,  2137,  2138,  2139,  2140,     0,
       0,  2142,     0,     0,  2144,     0,     0,     0,     0,     0,
       0,     0,  2150,  2152,     0,  2153,     0,  2154,  2155,     0,
       0,  2158,     0,  2159,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2163,  2164,  2166,     0,
    2167,     0,     0,     0,     0,  2168,     0,     0,     0,  2171,
    2172,  2173,  2174,  2176,     0,  2178,     0,  2180,     0,     0,
       0,     0,     0,     0,     0,     0,  2184,     0,  2187,     0,
       0,     0,     0,     0,     0,     0,  2191,  2192,     0,  2193,
       0,  2196,     0,     0,     0,  2198,  2199,  2200,     0,  2201,
       0,     0,     0,  2204,     0,     0,     0,  2207,     0,     0,
    2209,     0,     0,     0,  2211,     0,  2214,     0,     0,  2216,
       0,  2218,  2219,  2220,  2221,     0,     0,  2222,     0,     0,
    2224,     0,     0,     0,  2227,     0,     0,     0,     0,     0,
       0,  2230,  2232,     0,  2234,  2236,     0,  2237,     0,     0,
    2238,     0,     0,     0,     0,  2241,     0,  2242,     0,     0,
    2244,  2246,     0,     0,     0,  2251,     0,  2252,     0,     0,
       0,     0,  2255,     0,  2256,  2257,     0,     0,     0,  2263,
    2264,   373,     0,  2265,     0,     0,  2267,  2268,  2269,  2270,
    2271,  2272,   391,   392,  2275,  2276,     0,   396,   397,  2280,
    2281,     0,   401,  2283,  2284,  2285,  2286,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   427,   428,     0,
     435,     0,     0,     0,     0,     0,   444,   445,   446,     0,
     449,     0,     0,     0,   454,     0,     0,   457,   458,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   470,
     472,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   491,     0,   494,
       0,     0,   497,   499,     0,     0,     0,     0,     0,     0,
     506,   508,     0,     0,   512,   514,   516,     0,     0,     0,
     520,   521,   522,   523,   525,   526,     0,   529,     0,     0,
       0,   533,   534,     0,   537,   539,     0,     0,   542,   543,
       0,     0,     0,     0,   548,   550,   551,     0,   553,   554,
       0,     0,   557,     0,   560,     0,     0,     0,     0,     0,
       0,   568,     0,     0,     0,   572,   573,   574,     0,   576,
     577,     0,   579,   580,   581,     0,   584,     0,   586,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   601,     0,     0,     0,     0,     0,   607,
     608,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   622,   623,     0,   625,   626,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   644,     0,     0,     0,
       0,     0,   665,     0,     0,     0,     0,   680,   682,   683,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   714,     0,   717,   718,     0,     0,   720,
     722,   723,   724,     0,   726,     0,   728,   730,   732,   734,
       0,     0,   726,   730,   734,     0,     0,   739,   741,     0,
       0,   743,   745,     0,     0,   747,   748,   749,     0,     0,
     753,     0,     0,     0,     0,     0,   755,   757,   758,     0,
       0,   759,   760,   762,   763,   764,     0,   766,   680,   680,
     769,   770,   772,     0,   773,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   799,
     801,   802,     0,     0,     0,   814,     0,     0,     0,     0,
       0,     0,     0,   822,     0,     0,     0,     0,     0,     0,
       0,   827,     0,     0,   831,   831,     0,     0,     0,     0,
     837,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   855,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   865,     0,     0,     0,     0,
       0,     0,     0,     0,   728,   726,     0,     0,     0,     0,
       0,     0,   883,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   898,     0,     0,   898,     0,
     904,     0,     0,   730,     0,     0,   734,   732,     0,     0,
     918,   814,     0,     0,     0,     0,   924,   814,     0,     0,
     928,     0,     0,   931,     0,     0,     0,     0,     0,   938,
       0,     0,     0,     0,     0,     0,     0,   755,     0,   951,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   958,
       0,     0,     0,     0,     0,     0,     0,     0,   766,   964,
       0,     0,     0,     0,     0,     0,     0,     0,   971,     0,
       0,   974,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   984,     0,   986,     0,   989,     0,
       0,     0,     0,     0,     0,     0,   814,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1015,     0,     0,     0,     0,     0,
       0,     0,   757,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   814,   814,     0,     0,     0,     0,     0,     0,
       0,     0,  1039,     0,     0,     0,     0,     0,     0,     0,
    1052,     0,     0,  1056,     0,     0,     0,     0,     0,     0,
       0,  1061,   762,  1065,     0,     0,     0,     0,     0,     0,
    1067,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1082,     0,     0,
       0,     0,     0,     0,     0,     0,  1091,  1092,  1093,  1094,
    1095,     0,  1097,  1097,  1100,  1100,     0,     0,     0,     0,
    1105,     0,  1108,  1110,  1112,     0,  1121,     0,     0,     0,
       0,     0,     0,     0,     0,  1134,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1143,
       0,     0,  1146,     0,  1148,     0,  1150,     0,  1152,     0,
       0,     0,     0,  1156,     0,  1157,     0,  1158,     0,  1159,
       0,  1160,     0,     0,     0,     0,     0,     0,  1162,     0,
    1163,     0,     0,     0,     0,     0,     0,     0,  1146,     0,
    1148,     0,     0,     0,  1178,  1179,  1180,  1181,     0,  1182,
    1183,  1184,     0,     0,     0,     0,  1186,  1187,   831,  1190,
       0,  1195,  1196,  1197,  1198,  1199,     0,  1201,  1202,  1205,
       0,     0,     0,     0,     0,     0,  1209,     0,     0,  1213,
       0,     0,  1217,     0,     0,     0,  1220,     0,     0,  1237,
       0,  1239,     0,  1241,     0,     0,  1244,     0,     0,  1248,
       0,  1148,     0,     0,     0,     0,     0,     0,  1254,     0,
       0,     0,     0,     0,     0,  1260,     0,   831,     0,   831,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1272,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1282,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1285,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   898,     0,  1293,  1295,   898,   898,   904,
       0,  1298,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1307,  1309,     0,  1248,     0,     0,   924,
       0,     0,  1248,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1325,  1327,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1335,     0,  1336,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1341,     0,     0,     0,     0,     0,     0,  1343,     0,
       0,  1344,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1351,     0,     0,     0,     0,     0,     0,     0,
       0,  1248,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1355,  1359,  1361,     0,     0,     0,
    1015,     0,  1367,     0,     0,     0,     0,     0,     0,     0,
       0,  1371,     0,     0,     0,     0,     0,  1248,  1248,     0,
       0,     0,     0,     0,  1377,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1061,     0,     0,     0,
    1065,     0,  1399,  1400,  1401,     0,     0,     0,     0,     0,
       0,  1407,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1425,     0,  1426,
       0,     0,  1429,     0,     0,     0,  1435,  1436,  1437,  1438,
       0,  1439,     0,  1440,     0,  1441,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1148,  1459,     0,  1462,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1476,
    1478,     0,  1482,     0,     0,     0,  1482,     0,  1487,  1488,
    1489,     0,     0,     0,     0,  1495,  1495,  1497,  1498,     0,
    1500,  1501,     0,     0,  1508,     0,     0,     0,     0,  1511,
       0,  1515,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1529,     0,     0,     0,
       0,     0,     0,     0,  1537,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1557,  1558,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1148,  1560,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1575,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   898,     0,  1585,     0,     0,     0,  1591,
       0,     0,  1309,     0,  1593,     0,  1593,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1327,     0,  1605,     0,  1605,     0,     0,     0,     0,     0,
    1609,     0,     0,     0,  1612,     0,     0,     0,     0,     0,
    1616,  1617,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1630,     0,     0,     0,     0,  1633,
       0,     0,     0,     0,     0,  1638,     0,     0,     0,     0,
       0,     0,     0,     0,  1645,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1655,  1657,     0,
    1482,     0,     0,     0,  1663,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1665,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1676,     0,  1679,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1686,     0,     0,  1689,     0,     0,  1691,
       0,     0,     0,     0,     0,     0,  1697,     0,  1697,  1697,
    1697,  1697,     0,     0,     0,     0,  1704,  1482,     0,     0,
       0,     0,     0,  1482,  1711,     0,     0,     0,     0,     0,
       0,  1719,     0,  1719,  1723,     0,  1697,     0,  1731,     0,
       0,     0,     0,     0,     0,  1739,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1745,     0,     0,     0,     0,     0,
       0,     0,  1749,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1770,  1772,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1785,  1786,     0,     0,  1789,  1791,  1792,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1802,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1817,     0,     0,     0,     0,     0,
       0,  1821,     0,  1823,     0,     0,     0,  1824,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1833,  1834,
    1835,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1697,     0,     0,     0,     0,  1482,     0,     0,     0,
    1846,  1847,  1848,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1854,     0,  1855,  1856,     0,  1857,     0,     0,
    1859,     0,  1861,     0,     0,  1865,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1877,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1903,     0,     0,     0,
       0,     0,  1907,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1926,  1927,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1938,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1950,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1958,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1965,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1974,     0,     0,     0,     0,     0,
       0,     0,  1981,     0,  1983,     0,  1984,     0,     0,     0,
       0,     0,  1987,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1993,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2037,     0,  2039,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2052,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2067,     0,  2068,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2088,     0,     0,     0,  2091,     0,     0,     0,     0,
       0,     0,     0,     0,  2095,     0,  2096,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2110,     0,     0,     0,     0,     0,
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
       0,     0,     0,     0,     0,  2181,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2197,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     1,     2,     3,     0,     4,     0,     5,     6,     0,
       0,     7,     0,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,     0,    22,    23,
      24,    25,    26,    27,    28,     0,    29,  2239,    30,     0,
      31,    32,    33,    34,    35,     0,  2247,    36,    37,     0,
      38,    39,    40,    41,  2254,     0,     0,    42,    43,     0,
       0,     0,    44,    45,    46,     0,     0,     0,    47,  2266,
      48,     0,    49,    50,    51,    52,     0,     0,    53,    54,
      55,   641,    56,    57,     0,     0,  2282,    58,    59,     0,
       0,    60,    61,     0,     0,     0,    62,    63,     0,    64,
      65,     0,     0,     0,    66,    67,    68,    69,    70,    71,
       0,    72,    73,    74,    75,    76,    77,    78,     0,    79,
      80,    81,     0,    82,    83,    84,    85,    86,    87,    88,
      89,    90,    91,    92,    93,    94,     0,     0,     0,     0,
       0,     0,     0,    95,    96,     0,    97,    98,     0,    99,
     100,   101,   102,   103,     0,   104,     0,   105,   106,   107,
       0,     0,     0,     0,   108,     0,   109,   110,   111,   112,
     113,   114,   115,   116,     0,     0,     0,     0,     0,   117,
       0,   118,     0,   119,   120,     0,     0,   121,   122,   123,
     124,   125,   126,     0,   127,     0,   128,     0,     0,   129,
       0,   130,   131,   132,   133,   134,     0,   135,     0,   136,
     137,   138,   139,     0,   140,   141,   142,     0,   143,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   144,
     145,   146,   147,     0,   148,     0,   149,     0,   150,     0,
       0,   151,   152,   153,   154,   155,   156,   157,     0,   158,
       0,   159,     0,   160,   161,     0,     0,   162,   163,     0,
       0,     0,     0,   164,     0,   165,     0,   166,   167,   168,
     169,   170,   171,   172,   173,   174,   175,   176,   177,   178,
       0,   179,     0,   180,   181,   182,   183,     0,     0,     0,
       0,     0,   184,     0,     0,     0,     0,   185,     0,   186,
     187,   188,     0,   189,   190,   191,   192,   193,   194,   195,
     196,   197,     0,     0,     0,     0,     0,   198,     0,     0,
     199,   200,   201,   202,     0,     0,   203,   204,   205,   206,
     207,     0,   208,     0,     0,     0,     0,     0,   209,   210,
       0,     0,   211,     0,     0,   212,   213,     0,     0,     0,
     214,   215,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   216,   217,
     218,   219,     0,     0,     0,   220,     0,     0,     0,     0,
       0,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,     0,     0,     0,     0,   231,   232,     1,     2,     3,
       0,     4,     0,     5,     6,     0,     0,     7,     0,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,     0,    22,    23,    24,    25,    26,    27,
      28,     0,    29,     0,    30,     0,    31,    32,    33,    34,
      35,     0,     0,    36,    37,     0,    38,    39,    40,    41,
       0,     0,     0,    42,    43,     0,     0,     0,    44,    45,
      46,     0,     0,     0,    47,     0,    48,     0,    49,    50,
      51,    52,     0,     0,    53,    54,    55,     0,    56,    57,
       0,     0,     0,    58,    59,     0,     0,    60,    61,     0,
       0,     0,    62,    63,     0,    64,    65,     0,     0,     0,
      66,    67,    68,    69,    70,    71,     0,    72,    73,    74,
      75,    76,    77,    78,     0,    79,    80,    81,     0,    82,
      83,    84,    85,    86,    87,    88,    89,    90,    91,    92,
      93,    94,     0,     0,     0,     0,     0,     0,     0,    95,
      96,     0,    97,    98,     0,    99,   100,   101,   102,   103,
       0,   104,     0,   105,   106,   107,     0,     0,     0,     0,
     108,     0,   109,   110,   111,   112,   113,   114,   115,   116,
       0,     0,     0,     0,     0,   117,     0,   118,     0,   119,
     120,     0,     0,   121,   122,   123,   124,   125,   126,     0,
     127,     0,   128,     0,     0,   129,     0,   130,   131,   132,
     133,   134,     0,   135,     0,   136,   137,   138,   139,     0,
     140,   141,   142,     0,   143,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   144,   145,   146,   147,     0,
     148,     0,   149,     0,   150,     0,     0,   151,   152,   153,
     154,   155,   156,   157,     0,   158,     0,   159,     0,   160,
     161,     0,     0,   162,   163,     0,     0,     0,     0,   164,
       0,   165,     0,   166,   167,   168,   169,   170,   171,   172,
     173,   174,   175,   176,   177,   178,     0,   179,     0,   180,
     181,   182,   183,     0,     0,     0,     0,     0,   184,     0,
       0,     0,     0,   185,     0,   186,   187,   188,     0,   189,
     190,   191,   192,   193,   194,   195,   196,   197,     0,     0,
       0,     0,     0,   198,     0,     0,   199,   200,   201,   202,
       0,     0,   203,   204,   205,   206,   207,     0,   208,     0,
       0,     0,     0,     0,   209,   210,     0,     0,   211,     0,
       0,   212,   213,     0,     0,     0,   214,   215,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   216,   217,   218,   219,     0,     0,
       0,   220,     0,     0,     0,     0,     0,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,     0,     0,     0,
       0,   231,   232
};

#define yypact_value_is_default(yystate) \
  ((yystate) == (-1587))

#define yytable_value_is_error(yytable_value) \
  YYID (0)

static const yytype_int16 yycheck[] =
{
       7,     8,     9,    10,   288,   471,   896,   917,    11,     8,
     900,   937,    71,    69,    73,    29,    97,   432,  1169,    14,
      69,    67,   147,   147,    67,    63,   147,    41,    63,    36,
    1616,    38,    78,    79,   147,    78,    79,    12,     6,   323,
     200,    67,    37,    62,    57,   103,   104,   105,    86,  1480,
      63,    86,    87,   200,    63,  1486,    62,  1157,  1158,  1159,
    1160,   200,   200,   147,    90,   200,    73,    62,    62,    68,
      77,   307,    79,    86,   147,    88,   147,    86,   306,    97,
     147,    56,  1182,    96,   144,   200,   305,   147,    63,    96,
     200,   147,   200,   147,   200,   120,    95,   147,   147,   147,
     147,   268,   147,   200,   200,   147,    62,    63,   147,   147,
     147,    86,   147,   112,   147,   147,   147,   147,   147,   147,
     158,   200,   147,   158,   150,   147,   147,    97,   182,   124,
      86,   200,   200,   290,   147,    11,   124,   142,   147,   235,
     200,   341,   342,   150,   200,   158,   153,   154,   143,   158,
     200,   200,   200,   263,   200,   263,   235,   200,   188,   166,
     167,   147,   184,    51,   171,   200,   200,   235,   167,   200,
     200,    47,   147,   147,   147,    63,   183,   147,   177,   186,
     187,   188,   200,   158,   147,   180,   353,   295,   195,   200,
     237,   147,   180,   290,   200,   200,   320,   204,    86,   320,
     195,   176,   158,   202,   199,   212,   268,   195,   200,   208,
     289,   218,   296,   351,   374,   210,   351,   224,   333,   438,
     239,   289,   197,   242,   235,   200,   293,   374,  1659,  1815,
     205,  1817,   200,   346,   359,   374,   375,   376,   257,   246,
     247,   248,   268,   479,   200,   239,   253,   320,   476,   320,
     200,   296,  1403,   252,   296,   315,   293,   296,   264,   147,
     293,   293,   297,   257,   314,   315,   314,   296,   296,   290,
     158,   141,   341,   354,   355,   356,   357,   358,   316,   360,
     361,   362,   363,   364,   297,   366,   367,   368,   369,    62,
      63,   349,   266,   266,   350,   341,   291,   343,   341,  1399,
     771,   350,    51,   266,   317,   300,    63,   341,   264,   284,
     282,   320,   200,    86,    63,   326,    63,   177,   178,   378,
     379,   380,   381,   382,   338,   200,    63,   326,    63,    86,
     200,    63,    92,   142,   200,    63,   147,    86,   333,    86,
      97,    63,   345,   200,    91,    63,   318,   354,   200,    86,
      63,    86,   200,   147,    86,    63,    91,    63,    86,   830,
     235,   832,   156,    63,    86,    63,   393,   200,    86,   229,
     142,   200,   399,    86,   147,   200,   383,  1315,    86,   386,
      86,   200,   389,   390,    39,   158,    86,   147,    86,   200,
     147,   200,    67,   250,   401,  1305,   142,  1323,   147,   142,
     147,   158,   409,    78,    79,    67,  1296,  1345,  1346,   158,
     147,   158,   147,   173,   289,   147,    78,    79,   142,   147,
     427,   158,   429,   158,   431,   147,   158,   200,   200,   147,
     158,   715,   292,   440,   147,   901,   158,   444,    63,   147,
     158,   147,   200,   200,   304,   158,    63,   147,   147,   147,
     158,   200,   158,   200,   200,   462,    63,   200,   158,   142,
     158,    86,   147,   470,   200,   200,    63,   751,   200,    86,
     477,   478,   200,    63,    97,    63,   200,   290,   200,    86,
     290,    63,   200,    63,   142,   492,   142,   200,   495,    86,
     980,   981,   200,   500,   200,    63,    86,   268,    86,   506,
     200,   200,   200,   147,    86,   200,    86,   162,   163,   164,
     165,    63,    63,    63,    63,   200,   523,   200,    86,   147,
      63,   147,   147,   200,   147,   200,   120,   534,   812,    63,
     147,   147,    63,   158,    86,    86,    86,    86,   200,    63,
     147,   158,   200,    86,   200,   268,   147,   147,   555,   556,
     147,   158,    86,   147,   188,    86,   200,   147,   484,   147,
     120,   158,    86,   200,   490,   147,   200,   147,   158,    63,
     158,   200,   200,   147,   200,   200,   158,   584,   158,   147,
     147,   588,   200,   200,   200,   592,    97,   147,   200,   596,
     158,   147,    86,   200,   601,   147,   147,   147,   147,   200,
     200,   147,   609,   200,   147,    63,   158,   158,   158,   158,
     200,   618,   200,   147,    63,   158,   147,   147,   200,    63,
     200,   147,   290,   147,   158,   120,   200,   158,    86,   147,
     290,   200,   200,   200,   158,   919,   147,    86,   645,   646,
     647,   925,    86,   569,   200,   138,   139,   140,   200,   200,
     200,   200,   147,   147,   200,    63,   200,   200,   665,   666,
     667,   668,   200,   200,   158,   147,   200,   200,    63,   200,
     200,   177,   178,   680,   200,   682,   200,   684,    86,   605,
     606,   688,   200,   690,   691,   692,   200,   200,   695,   147,
     200,    86,   699,   200,   290,   200,    63,   704,   147,   706,
     158,    63,    63,   147,    63,    63,   200,   714,   147,   158,
     994,    63,   200,   720,   158,   722,   723,   724,   200,    86,
      63,   147,   180,   229,    86,    86,   200,    86,    86,   147,
      63,   187,   188,   200,    86,    63,   200,    63,   745,   147,
     747,   748,   749,    86,   200,    63,  1030,  1031,   755,   147,
     158,   200,   147,    86,   147,   290,   200,   200,    86,   766,
      86,   200,   200,   158,    67,   772,   773,    63,    86,   147,
      63,    63,   779,   200,   200,    78,    79,   200,   200,   786,
     147,   788,   200,   790,   791,   147,   147,   794,   147,   147,
      86,   158,   200,    86,    86,   147,   158,   158,   200,   158,
     158,    67,   200,   810,   147,   200,   158,   200,    63,   816,
     817,    63,    78,    79,   147,   158,   147,   200,   825,   147,
     200,   147,   200,   147,   831,   158,   833,   834,   290,   147,
     158,    86,   158,   200,    86,    63,   200,   844,   200,   200,
     158,   200,   200,   200,   147,    63,   200,   200,   200,   188,
     200,   147,    63,   200,   147,   147,   200,   200,    86,   866,
     867,   200,   158,   200,   200,   158,   158,   200,    86,   200,
      63,    63,   200,    63,   200,    86,   200,   884,   200,    63,
     200,   147,   200,   200,   147,   892,    63,   200,   200,    63,
     200,   200,   147,    86,    86,   147,    86,   200,   905,   906,
     907,   908,    86,   158,   200,   200,   158,   200,   200,    86,
      63,   918,    86,    67,  1329,   200,    63,   924,    63,   147,
     927,   928,   147,   930,    78,    79,   147,   147,   935,   147,
     158,   938,    67,    86,   200,    63,   147,   200,    63,    86,
     158,    86,   147,    78,    79,   200,   147,   158,   200,   200,
     147,   200,   200,    97,   147,   147,   200,   147,    86,   200,
     268,    86,   200,   147,    63,   158,   158,    63,   158,   200,
     147,   147,   200,   147,   158,   200,   147,   200,   200,   200,
     200,   158,   200,    63,   158,   200,    63,    86,   200,   200,
      86,   200,    63,   147,   147,   200,   147,   200,   200,   200,
     147,   200,   147,   200,  1011,   158,    86,   200,   200,    86,
     200,   158,   147,   158,   200,    86,   200,   200,  1025,   147,
     147,  1028,   147,   200,   200,    97,   200,    63,  1035,   200,
     158,  1038,   147,   158,   200,   147,   147,  1044,  1045,  1046,
    1047,  1048,   147,   200,  1051,   200,   200,   200,   147,   200,
      86,   147,   147,   200,  1061,   200,    63,   147,    97,   158,
      63,  1068,   158,  1070,   147,   200,  1073,   147,    63,   200,
     147,   147,   200,   200,   200,   200,   147,   200,   158,    86,
     147,   158,  1497,    86,   290,   200,  1501,   158,   200,   200,
    1097,    86,   200,   147,   200,   200,  1103,  1104,  1105,  1106,
     200,   200,   290,    63,   200,   200,   200,   200,   200,   200,
     200,   147,    63,  1120,   200,   200,  1123,   200,  1125,    63,
     200,   147,   158,   200,   200,   200,    86,  1134,  1135,   200,
    1137,   200,  1139,   200,   200,    86,   200,  1144,    67,   200,
     147,  1148,    86,  1150,   147,  1152,   200,  1154,  1155,    78,
      79,   158,   147,    63,   200,   158,    97,  1164,    63,  1166,
    1167,  1168,    63,   158,   200,   200,   200,  1174,  1175,  1176,
    1177,    63,   200,   321,   200,   200,    86,   200,  1185,  1186,
     200,    86,   150,  1190,   200,    86,  1193,   147,   200,    63,
     200,   200,    63,   200,    86,   200,   147,   200,   158,  1206,
     200,  1208,   200,   147,   200,   200,   200,   158,  1215,   200,
     200,   200,    86,  1220,   158,    86,  1223,  1224,  1225,  1226,
    1227,  1228,  1229,  1230,  1231,  1232,  1233,  1234,  1235,  1236,
     200,     0,    48,   120,    97,   177,    48,   147,    63,    63,
     200,   180,   147,   180,  1251,  1252,   147,    63,   158,   200,
     380,    63,    58,   158,   200,   147,   200,   158,  1265,     6,
     200,    86,    86,     6,   394,   395,   158,   200,   200,   399,
      86,   401,  1279,   147,    86,  1282,   147,   200,   200,   200,
    1287,   200,  1289,   200,   158,   268,  1293,   158,  1295,   200,
     200,  1298,  1299,  1300,  1301,   200,   200,   200,    97,   200,
    1307,   200,   200,  1310,   200,   200,   200,   200,   200,    79,
    1317,  1318,  1319,  1320,   200,  1322,   200,   200,  1325,   200,
     200,   200,   147,   147,   200,   200,   200,   200,   200,   200,
     200,   147,   200,   158,   158,   147,   335,   336,   337,   338,
     339,   200,   158,   200,   200,   200,   158,   200,   200,   200,
     200,   200,  1359,   200,  1361,   200,   200,   200,   200,   200,
    1367,   200,   200,   200,   254,   200,  1373,   200,   200,   200,
     200,   200,  1379,   200,   200,    48,  1383,  1384,  1385,   200,
     200,   200,   200,   219,  1391,   221,   222,   223,   224,   225,
     226,   227,   228,   200,    48,  1402,   200,   200,   200,   200,
     200,   261,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,  1424,  1425,   200,
    1427,  1428,  1429,   200,   200,   200,   200,   200,   200,   200,
     258,    97,   200,   200,   200,    97,    97,    97,  1445,    97,
     200,    97,   200,   200,   200,  1452,  1453,   200,  1455,   200,
     290,  1458,   200,    97,   200,   200,   290,    61,   200,   200,
    1467,  1468,   200,   200,   200,   200,   200,    97,    97,    97,
     200,    97,    97,   200,   200,  1482,  1483,  1484,  1485,   200,
     200,  1488,  1489,  1490,   200,  1492,  1493,   200,   200,   200,
     200,  1498,   200,  1500,   200,  1910,  1503,  1504,   200,   200,
    1507,   200,   200,   200,  1511,   200,  1513,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,  1534,   200,   200,
     200,   200,  1539,   200,  1541,  1542,  1543,  1544,  1545,  1546,
    1547,  1548,  1549,  1550,  1551,  1552,  1553,  1554,  1555,   200,
     200,   200,   200,   200,   200,   200,    97,   200,   200,   200,
     200,   200,   200,  1570,  1571,  1572,   200,   200,   200,  1576,
     298,    97,    97,  1580,   353,  1582,   200,   200,  1585,   200,
     200,  1588,   200,   200,   200,   200,  1593,  1594,   200,  2004,
     200,  1598,   200,  2008,   200,   200,   678,   200,  1605,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,  1630,   200,  1632,   200,   254,   200,   200,
    1637,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,  1648,   200,  1650,   200,  1652,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,  1666,
    1667,   200,  1669,  1670,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,  1681,    97,   200,   200,   200,  1686,
     200,   200,  1689,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,  1704,   200,   200,
    1707,  1708,  1709,   200,  1711,  1712,  1713,  1714,   200,  1716,
    1717,   200,   200,   200,   200,   200,  1723,   200,  1725,   200,
     200,  1728,   200,   200,  1731,   200,   200,  1734,   200,   200,
     200,  1738,   200,   200,  1741,   200,   254,   200,   200,   200,
     200,  1748,   200,   200,  1751,   200,   200,   200,  1755,  1756,
    1757,  1758,  1759,  1760,  1761,  1762,  1763,  1764,  1765,  1766,
    1767,  1768,   200,   200,   200,  1772,  2181,   200,   200,  1776,
    1777,   200,   200,   200,   230,   200,   200,   200,    97,  1786,
     200,   200,   200,   200,  1791,  1792,   200,  1794,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,  1811,   200,  1813,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,  1824,   200,   200,
    1827,  1828,   200,  1830,  1831,   200,   200,   200,   200,   200,
     200,   200,   200,  1840,  1841,   200,   200,   200,    97,   200,
     200,  1848,  1849,   200,  1851,  1852,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   353,   676,  1871,  1059,   480,  1874,  1875,  1876,
     902,  1013,  1879,   695,  1881,  1882,  1883,   936,  1885,  1886,
    1887,  1063,   437,  1890,   922,   388,   916,  1287,  1895,  1179,
     595,   510,    -1,    -1,  1901,  1902,   234,    -1,  1905,    -1,
     531,    -1,   635,    -1,  1911,  1912,    -1,    -1,  1915,  1916,
      -1,  1918,  1919,  1920,  1921,    -1,  1923,  1924,    -1,  1926,
      -1,    -1,  1929,  1930,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1938,    -1,    -1,  1941,    -1,  1943,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1962,    -1,  1964,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1974,    -1,  1976,
      -1,  1978,  1979,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1987,    -1,  1989,  1990,  1991,  1992,  1993,    -1,  1995,  1996,
    1997,  1998,    -1,  2000,  2001,  2002,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2020,    -1,  2022,  2023,    -1,    -1,    -1,
    2027,    -1,    -1,    -1,  2031,  2032,  2033,    -1,  2035,    -1,
      -1,    -1,  2039,    -1,    -1,    -1,    -1,    -1,    -1,  2046,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2063,  2064,    -1,  2066,
      -1,    -1,    -1,    -1,    -1,  2072,    -1,    -1,    -1,    -1,
    2077,    -1,  2079,    -1,  2081,  2082,  2083,  2084,  2085,    -1,
      -1,  2088,    -1,    -1,  2091,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2099,  2100,    -1,  2102,    -1,  2104,  2105,    -1,
      -1,  2108,    -1,  2110,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2123,  2124,  2125,    -1,
    2127,    -1,    -1,    -1,    -1,  2132,    -1,    -1,    -1,  2136,
    2137,  2138,  2139,  2140,    -1,  2142,    -1,  2144,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2153,    -1,  2155,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2163,  2164,    -1,  2166,
      -1,  2168,    -1,    -1,    -1,  2172,  2173,  2174,    -1,  2176,
      -1,    -1,    -1,  2180,    -1,    -1,    -1,  2184,    -1,    -1,
    2187,    -1,    -1,    -1,  2191,    -1,  2193,    -1,    -1,  2196,
      -1,  2198,  2199,  2200,  2201,    -1,    -1,  2204,    -1,    -1,
    2207,    -1,    -1,    -1,  2211,    -1,    -1,    -1,    -1,    -1,
      -1,  2218,  2219,    -1,  2221,  2222,    -1,  2224,    -1,    -1,
    2227,    -1,    -1,    -1,    -1,  2232,    -1,  2234,    -1,    -1,
    2237,  2238,    -1,    -1,    -1,  2242,    -1,  2244,    -1,    -1,
      -1,    -1,  2249,    -1,  2251,  2252,    -1,    -1,    -1,  2256,
    2257,     4,    -1,  2260,    -1,    -1,  2263,  2264,  2265,  2266,
    2267,  2268,    15,    16,  2271,  2272,    -1,    20,    21,  2276,
    2277,    -1,    25,  2280,  2281,  2282,  2283,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    50,    51,    -1,
      53,    -1,    -1,    -1,    -1,    -1,    59,    60,    61,    -1,
      63,    -1,    -1,    -1,    67,    -1,    -1,    70,    71,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    82,
      83,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   100,    -1,   102,
      -1,    -1,   105,   106,    -1,    -1,    -1,    -1,    -1,    -1,
     113,   114,    -1,    -1,   117,   118,   119,    -1,    -1,    -1,
     123,   124,   125,   126,   127,   128,    -1,   130,    -1,    -1,
      -1,   134,   135,    -1,   137,   138,    -1,    -1,   141,   142,
      -1,    -1,    -1,    -1,   147,   148,   149,    -1,   151,   152,
      -1,    -1,   155,    -1,   157,    -1,    -1,    -1,    -1,    -1,
      -1,   164,    -1,    -1,    -1,   168,   169,   170,    -1,   172,
     173,    -1,   175,   176,   177,    -1,   179,    -1,   181,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   196,    -1,    -1,    -1,    -1,    -1,   202,
     203,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   216,   217,    -1,   219,   220,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   239,    -1,    -1,    -1,
      -1,    -1,   245,    -1,    -1,    -1,    -1,   250,   251,   252,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   286,    -1,   288,   289,    -1,    -1,   292,
     293,   294,   295,    -1,   297,    -1,   299,   300,   301,   302,
      -1,    -1,   305,   306,   307,    -1,    -1,   310,   311,    -1,
      -1,   314,   315,    -1,    -1,   318,   319,   320,    -1,    -1,
     323,    -1,    -1,    -1,    -1,    -1,   329,   330,   331,    -1,
      -1,   334,   335,   336,   337,   338,    -1,   340,   341,   342,
     343,   344,   345,    -1,   347,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   362,
     363,   364,    -1,    -1,    -1,   368,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   376,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   384,    -1,    -1,   387,   388,    -1,    -1,    -1,    -1,
     393,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   412,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   428,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   437,   438,    -1,    -1,    -1,    -1,
      -1,    -1,   445,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   468,    -1,    -1,   471,    -1,
     473,    -1,    -1,   476,    -1,    -1,   479,   480,    -1,    -1,
     483,   484,    -1,    -1,    -1,    -1,   489,   490,    -1,    -1,
     493,    -1,    -1,   496,    -1,    -1,    -1,    -1,    -1,   502,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   510,    -1,   512,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   522,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   531,   532,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   541,    -1,
      -1,   544,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   557,    -1,   559,    -1,   561,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   569,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   587,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   595,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   605,   606,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   615,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     623,    -1,    -1,   626,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   634,   635,   636,    -1,    -1,    -1,    -1,    -1,    -1,
     643,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   660,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   669,   670,   671,   672,
     673,    -1,   675,   676,   677,   678,    -1,    -1,    -1,    -1,
     683,    -1,   685,   686,   687,    -1,   689,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   698,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   712,
      -1,    -1,   715,    -1,   717,    -1,   719,    -1,   721,    -1,
      -1,    -1,    -1,   726,    -1,   728,    -1,   730,    -1,   732,
      -1,   734,    -1,    -1,    -1,    -1,    -1,    -1,   741,    -1,
     743,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   751,    -1,
     753,    -1,    -1,    -1,   757,   758,   759,   760,    -1,   762,
     763,   764,    -1,    -1,    -1,    -1,   769,   770,   771,   772,
      -1,   774,   775,   776,   777,   778,    -1,   780,   781,   782,
      -1,    -1,    -1,    -1,    -1,    -1,   789,    -1,    -1,   792,
      -1,    -1,   795,    -1,    -1,    -1,   799,    -1,    -1,   802,
      -1,   804,    -1,   806,    -1,    -1,   809,    -1,    -1,   812,
      -1,   814,    -1,    -1,    -1,    -1,    -1,    -1,   821,    -1,
      -1,    -1,    -1,    -1,    -1,   828,    -1,   830,    -1,   832,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   851,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   870,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   881,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   896,    -1,   898,   899,   900,   901,   902,
      -1,   904,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   916,   917,    -1,   919,    -1,    -1,   922,
      -1,    -1,   925,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   936,   937,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   947,    -1,   949,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   964,    -1,    -1,    -1,    -1,    -1,    -1,   971,    -1,
      -1,   974,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   985,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   994,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1007,  1008,  1009,    -1,    -1,    -1,
    1013,    -1,  1015,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1024,    -1,    -1,    -1,    -1,    -1,  1030,  1031,    -1,
      -1,    -1,    -1,    -1,  1037,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1059,    -1,    -1,    -1,
    1063,    -1,  1065,  1066,  1067,    -1,    -1,    -1,    -1,    -1,
      -1,  1074,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1100,    -1,  1102,
      -1,    -1,  1105,    -1,    -1,    -1,  1109,  1110,  1111,  1112,
      -1,  1114,    -1,  1116,    -1,  1118,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1146,  1147,    -1,  1149,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1162,
    1163,    -1,  1165,    -1,    -1,    -1,  1169,    -1,  1171,  1172,
    1173,    -1,    -1,    -1,    -1,  1178,  1179,  1180,  1181,    -1,
    1183,  1184,    -1,    -1,  1187,    -1,    -1,    -1,    -1,  1192,
      -1,  1194,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1209,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1217,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1237,  1238,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1248,  1249,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1273,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1296,    -1,  1298,    -1,    -1,    -1,  1302,
      -1,    -1,  1305,    -1,  1307,    -1,  1309,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1323,    -1,  1325,    -1,  1327,    -1,    -1,    -1,    -1,    -1,
    1333,    -1,    -1,    -1,  1337,    -1,    -1,    -1,    -1,    -1,
    1343,  1344,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1367,    -1,    -1,    -1,    -1,  1372,
      -1,    -1,    -1,    -1,    -1,  1378,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1387,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1400,  1401,    -1,
    1403,    -1,    -1,    -1,  1407,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1419,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1435,    -1,  1437,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1456,    -1,    -1,  1459,    -1,    -1,  1462,
      -1,    -1,    -1,    -1,    -1,    -1,  1469,    -1,  1471,  1472,
    1473,  1474,    -1,    -1,    -1,    -1,  1479,  1480,    -1,    -1,
      -1,    -1,    -1,  1486,  1487,    -1,    -1,    -1,    -1,    -1,
      -1,  1494,    -1,  1496,  1497,    -1,  1499,    -1,  1501,    -1,
      -1,    -1,    -1,    -1,    -1,  1508,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1527,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1535,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1559,  1560,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1584,  1585,    -1,    -1,  1588,  1589,  1590,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1607,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1617,    -1,    -1,    -1,    -1,    -1,
      -1,  1624,    -1,  1626,    -1,    -1,    -1,  1630,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1641,  1642,
    1643,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1654,    -1,    -1,    -1,    -1,  1659,    -1,    -1,    -1,
    1663,  1664,  1665,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1675,    -1,  1677,  1678,    -1,  1680,    -1,    -1,
    1683,    -1,  1685,    -1,    -1,  1688,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1710,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1739,    -1,    -1,    -1,
      -1,    -1,  1745,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1770,  1771,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1789,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1807,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1823,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1833,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1847,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1855,    -1,  1857,    -1,  1859,    -1,    -1,    -1,
      -1,    -1,  1865,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1877,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1925,    -1,  1927,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1958,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1980,    -1,  1982,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2004,    -1,    -1,    -1,  2008,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2017,    -1,  2019,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2037,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2148,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2171,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,     3,     4,     5,    -1,     7,    -1,     9,    10,    -1,
      -1,    13,    -1,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    -1,    30,    31,
      32,    33,    34,    35,    36,    -1,    38,  2230,    40,    -1,
      42,    43,    44,    45,    46,    -1,  2239,    49,    50,    -1,
      52,    53,    54,    55,  2247,    -1,    -1,    59,    60,    -1,
      -1,    -1,    64,    65,    66,    -1,    -1,    -1,    70,  2262,
      72,    -1,    74,    75,    76,    77,    -1,    -1,    80,    81,
      82,    83,    84,    85,    -1,    -1,  2279,    89,    90,    -1,
      -1,    93,    94,    -1,    -1,    -1,    98,    99,    -1,   101,
     102,    -1,    -1,    -1,   106,   107,   108,   109,   110,   111,
      -1,   113,   114,   115,   116,   117,   118,   119,    -1,   121,
     122,   123,    -1,   125,   126,   127,   128,   129,   130,   131,
     132,   133,   134,   135,   136,   137,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   145,   146,    -1,   148,   149,    -1,   151,
     152,   153,   154,   155,    -1,   157,    -1,   159,   160,   161,
      -1,    -1,    -1,    -1,   166,    -1,   168,   169,   170,   171,
     172,   173,   174,   175,    -1,    -1,    -1,    -1,    -1,   181,
      -1,   183,    -1,   185,   186,    -1,    -1,   189,   190,   191,
     192,   193,   194,    -1,   196,    -1,   198,    -1,    -1,   201,
      -1,   203,   204,   205,   206,   207,    -1,   209,    -1,   211,
     212,   213,   214,    -1,   216,   217,   218,    -1,   220,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   231,
     232,   233,   234,    -1,   236,    -1,   238,    -1,   240,    -1,
      -1,   243,   244,   245,   246,   247,   248,   249,    -1,   251,
      -1,   253,    -1,   255,   256,    -1,    -1,   259,   260,    -1,
      -1,    -1,    -1,   265,    -1,   267,    -1,   269,   270,   271,
     272,   273,   274,   275,   276,   277,   278,   279,   280,   281,
      -1,   283,    -1,   285,   286,   287,   288,    -1,    -1,    -1,
      -1,    -1,   294,    -1,    -1,    -1,    -1,   299,    -1,   301,
     302,   303,    -1,   305,   306,   307,   308,   309,   310,   311,
     312,   313,    -1,    -1,    -1,    -1,    -1,   319,    -1,    -1,
     322,   323,   324,   325,    -1,    -1,   328,   329,   330,   331,
     332,    -1,   334,    -1,    -1,    -1,    -1,    -1,   340,   341,
      -1,    -1,   344,    -1,    -1,   347,   348,    -1,    -1,    -1,
     352,   353,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   370,   371,
     372,   373,    -1,    -1,    -1,   377,    -1,    -1,    -1,    -1,
      -1,   383,   384,   385,   386,   387,   388,   389,   390,   391,
     392,    -1,    -1,    -1,    -1,   397,   398,     3,     4,     5,
      -1,     7,    -1,     9,    10,    -1,    -1,    13,    -1,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    -1,    30,    31,    32,    33,    34,    35,
      36,    -1,    38,    -1,    40,    -1,    42,    43,    44,    45,
      46,    -1,    -1,    49,    50,    -1,    52,    53,    54,    55,
      -1,    -1,    -1,    59,    60,    -1,    -1,    -1,    64,    65,
      66,    -1,    -1,    -1,    70,    -1,    72,    -1,    74,    75,
      76,    77,    -1,    -1,    80,    81,    82,    -1,    84,    85,
      -1,    -1,    -1,    89,    90,    -1,    -1,    93,    94,    -1,
      -1,    -1,    98,    99,    -1,   101,   102,    -1,    -1,    -1,
     106,   107,   108,   109,   110,   111,    -1,   113,   114,   115,
     116,   117,   118,   119,    -1,   121,   122,   123,    -1,   125,
     126,   127,   128,   129,   130,   131,   132,   133,   134,   135,
     136,   137,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   145,
     146,    -1,   148,   149,    -1,   151,   152,   153,   154,   155,
      -1,   157,    -1,   159,   160,   161,    -1,    -1,    -1,    -1,
     166,    -1,   168,   169,   170,   171,   172,   173,   174,   175,
      -1,    -1,    -1,    -1,    -1,   181,    -1,   183,    -1,   185,
     186,    -1,    -1,   189,   190,   191,   192,   193,   194,    -1,
     196,    -1,   198,    -1,    -1,   201,    -1,   203,   204,   205,
     206,   207,    -1,   209,    -1,   211,   212,   213,   214,    -1,
     216,   217,   218,    -1,   220,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   231,   232,   233,   234,    -1,
     236,    -1,   238,    -1,   240,    -1,    -1,   243,   244,   245,
     246,   247,   248,   249,    -1,   251,    -1,   253,    -1,   255,
     256,    -1,    -1,   259,   260,    -1,    -1,    -1,    -1,   265,
      -1,   267,    -1,   269,   270,   271,   272,   273,   274,   275,
     276,   277,   278,   279,   280,   281,    -1,   283,    -1,   285,
     286,   287,   288,    -1,    -1,    -1,    -1,    -1,   294,    -1,
      -1,    -1,    -1,   299,    -1,   301,   302,   303,    -1,   305,
     306,   307,   308,   309,   310,   311,   312,   313,    -1,    -1,
      -1,    -1,    -1,   319,    -1,    -1,   322,   323,   324,   325,
      -1,    -1,   328,   329,   330,   331,   332,    -1,   334,    -1,
      -1,    -1,    -1,    -1,   340,   341,    -1,    -1,   344,    -1,
      -1,   347,   348,    -1,    -1,    -1,   352,   353,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   370,   371,   372,   373,    -1,    -1,
      -1,   377,    -1,    -1,    -1,    -1,    -1,   383,   384,   385,
     386,   387,   388,   389,   390,   391,   392,    -1,    -1,    -1,
      -1,   397,   398
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     3,     4,     5,     7,     9,    10,    13,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    30,    31,    32,    33,    34,    35,    36,    38,
      40,    42,    43,    44,    45,    46,    49,    50,    52,    53,
      54,    55,    59,    60,    64,    65,    66,    70,    72,    74,
      75,    76,    77,    80,    81,    82,    84,    85,    89,    90,
      93,    94,    98,    99,   101,   102,   106,   107,   108,   109,
     110,   111,   113,   114,   115,   116,   117,   118,   119,   121,
     122,   123,   125,   126,   127,   128,   129,   130,   131,   132,
     133,   134,   135,   136,   137,   145,   146,   148,   149,   151,
     152,   153,   154,   155,   157,   159,   160,   161,   166,   168,
     169,   170,   171,   172,   173,   174,   175,   181,   183,   185,
     186,   189,   190,   191,   192,   193,   194,   196,   198,   201,
     203,   204,   205,   206,   207,   209,   211,   212,   213,   214,
     216,   217,   218,   220,   231,   232,   233,   234,   236,   238,
     240,   243,   244,   245,   246,   247,   248,   249,   251,   253,
     255,   256,   259,   260,   265,   267,   269,   270,   271,   272,
     273,   274,   275,   276,   277,   278,   279,   280,   281,   283,
     285,   286,   287,   288,   294,   299,   301,   302,   303,   305,
     306,   307,   308,   309,   310,   311,   312,   313,   319,   322,
     323,   324,   325,   328,   329,   330,   331,   332,   334,   340,
     341,   344,   347,   348,   352,   353,   370,   371,   372,   373,
     377,   383,   384,   385,   386,   387,   388,   389,   390,   391,
     392,   397,   398,   403,   404,   405,   406,   407,   408,   409,
     410,   415,   416,   417,   418,   419,   420,   421,   422,   423,
     424,   433,   435,   436,   437,   438,   439,   440,   441,   442,
     443,   444,   446,   447,   448,   449,   450,   451,   452,   456,
     459,   461,   462,   463,   464,   465,   466,   467,   468,   469,
     470,   472,   473,   475,   479,   480,   481,   483,   484,   485,
     486,   489,   492,   493,   494,   495,   496,   499,   502,   505,
     507,   509,   511,   515,   516,   517,   518,   519,   520,   521,
     522,   524,   526,   527,   528,   530,   532,   533,   534,   535,
     536,   537,   538,   539,   543,   545,   547,   551,   555,   557,
     559,   560,   561,   562,   563,   564,   565,   567,   568,   569,
     570,   578,   579,   581,   582,   583,   584,   585,   586,   588,
     589,   590,   592,   594,   598,   601,   603,   605,   606,   608,
     609,   610,   611,   612,   613,   615,   616,   618,   200,     6,
     200,   200,   147,   620,   200,    11,   345,    63,    86,   147,
     158,   621,   621,   621,    97,   200,   621,   200,   200,   200,
     200,   620,   620,   200,   200,   268,   620,   620,   200,   200,
     200,   620,   290,   290,   290,   268,   200,   200,   200,   200,
     200,    11,    47,   621,   268,   621,   200,   290,   290,   200,
     200,   200,   200,   200,   200,   200,   290,   620,   620,    67,
      78,    79,   200,   593,   200,   620,   200,   200,   200,   200,
     200,    69,   200,   350,   620,   620,   620,   200,   200,   620,
     200,   290,   200,   266,   620,   200,   200,   620,   620,   200,
     621,   200,   200,   200,   621,   200,   621,   290,   200,   200,
     620,   200,   620,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   333,   200,   621,   200,   200,
     200,   620,   200,   200,   620,   200,   200,   620,   156,   620,
     200,   200,   200,    97,   200,   268,   620,   200,   620,   200,
     200,   182,   620,   184,   620,   266,   620,   200,   200,   200,
     620,   620,   620,   620,   200,   620,   620,   200,   200,   620,
     200,   200,   200,   620,   620,    97,   200,   620,   200,   620,
     200,   200,   620,   620,   200,   200,   200,   200,   620,   237,
     620,   620,   621,   620,   620,   621,   621,   620,   200,   200,
     620,    97,   200,   200,   200,   290,   200,   266,   620,   200,
     621,   621,   620,   620,   620,   621,   620,   620,   290,   620,
     620,   620,   200,   200,   620,   200,   620,   200,   621,   200,
     200,   621,   621,   621,   200,   200,   200,   200,   200,   200,
     621,   620,   200,   200,    97,   200,   200,   620,   620,   621,
     200,   321,   200,   200,   200,   341,   200,   150,   621,   200,
     200,   200,   620,   620,   621,   620,   620,   200,   200,   290,
     200,   621,   200,   200,   200,   200,   200,   200,   200,   200,
       0,    83,   405,   120,   620,    62,   239,   242,   257,   411,
     412,   414,   478,   335,   336,   337,   338,   339,    71,    73,
     378,   379,   380,   381,   382,   620,   621,   621,   621,    39,
     162,   163,   164,   165,   425,   427,   428,   429,   430,   580,
     620,   434,   620,   620,   621,    97,   282,   318,   445,   177,
      14,    37,   124,   143,   180,   195,   199,   210,   291,   300,
     333,   453,   474,   478,   124,   180,   195,   457,   458,    29,
      41,   338,   460,   482,   620,   296,   574,   620,   620,   296,
     620,   296,   620,   620,   620,   500,   620,   506,   620,   508,
     620,   510,   620,   512,   620,   500,   508,   512,   523,   620,
     525,   620,   529,   620,   531,   620,   180,   620,   620,   620,
     180,   296,   574,   620,   558,   620,   572,   620,   620,   620,
     620,   566,   620,   620,   620,   571,   620,   580,   580,   620,
     620,   296,   620,   620,    58,   177,   178,   229,   292,   304,
     177,   178,   229,   411,   414,   478,   621,     8,    68,    95,
     112,   167,   177,   202,   208,   252,   326,   607,   200,   620,
     359,   620,   620,   393,   399,   617,   380,   394,   395,   399,
     401,   619,   540,   574,   620,   200,     6,     6,   200,   200,
     250,   346,   620,   200,   200,   621,   200,   620,   621,   513,
     514,   620,   514,   621,   621,   200,   200,   620,   200,   235,
     200,   200,   268,   200,   621,   200,   200,   200,   200,   235,
     289,   621,    97,   200,   200,   620,   200,   200,   235,   289,
     200,   200,   200,   200,   621,   620,   621,    79,   621,   593,
     142,   200,   200,   506,   500,   591,   621,   200,   200,    91,
     200,   621,   200,   620,   200,   200,   200,   200,   200,   200,
     200,   200,   621,   200,   200,   200,   556,   577,   620,   621,
     556,   200,   497,   498,   620,   138,   139,   140,   595,   508,
     602,   621,   604,   621,   512,   510,   552,   553,   620,   540,
     200,   200,   542,   576,   620,   540,   200,   621,   620,   200,
     621,   620,   200,   200,   200,   621,   548,   549,   620,   200,
     200,   235,   289,   621,   200,   558,   187,   188,   200,   188,
     200,   620,   200,   200,   200,   200,   200,   200,   620,   621,
     200,   200,   200,   571,   620,   200,   200,   621,   200,   200,
     200,   620,   200,   200,   620,   200,   200,   200,   200,   200,
     200,   200,   621,   621,   620,   471,   620,   200,   200,   620,
     254,   200,   200,   200,   540,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   200,    67,    90,   150,
     268,   621,   200,   503,   504,   620,   200,   621,   200,   200,
     621,   200,   572,   621,    48,   546,   200,   200,   621,   200,
     540,   540,   200,   200,   200,   621,   200,    48,   544,   620,
     141,   200,   200,   621,   103,   104,   105,   349,   413,   200,
     263,   295,   620,   200,   200,   200,   620,   200,   200,   487,
     488,   620,   566,   490,   491,   620,   296,   620,   261,   180,
     621,   621,   200,   621,   258,    97,    97,    97,    97,    97,
     200,    97,   620,   200,   200,   200,   200,   621,   621,   621,
     621,   620,   620,   620,   620,   620,   431,   620,   431,   432,
     620,   432,   316,   621,   621,   620,   621,    97,   620,    97,
     620,    97,   620,    12,    56,   176,   197,   200,   205,   284,
     621,   620,   454,   621,   200,   621,   455,   621,   290,   200,
     454,   200,   200,   290,   620,   621,   290,   621,   200,   621,
     200,    97,   200,   620,   621,   574,   620,   320,   620,    61,
     620,   621,   620,   621,   621,   621,   620,   620,   620,   620,
     620,   200,   620,   620,   621,   200,   621,   621,   621,   200,
     574,   320,    57,    88,    96,   297,   317,   621,   620,   620,
     620,   620,   620,   620,   620,   621,   620,   620,   513,   147,
     620,   621,   320,   621,   587,   620,   620,   620,   620,   620,
     621,   620,   620,    92,   173,   620,   621,   200,   621,   620,
     200,   621,   621,   620,   200,   621,   200,   620,   200,    97,
     620,    97,    97,   354,   355,   356,   357,   358,   360,   361,
     362,   363,   364,   366,   367,   368,   369,   620,    97,   620,
     200,   620,    97,    97,   620,   621,   200,   574,   620,   320,
     200,   621,   621,   200,   620,   200,   351,   621,   200,   200,
     620,   513,   621,   621,   621,   200,   200,   200,   200,   621,
     200,   200,   620,   200,   200,   200,   200,   200,   200,   621,
     621,   200,   620,   200,   200,   620,   200,   599,   600,   621,
     200,   621,   577,   620,   200,   620,   556,   498,   620,   621,
     621,   621,   621,   200,   200,   553,   554,   620,   554,   620,
      51,   200,   621,   576,   621,   596,   597,   621,   621,   621,
     621,   200,   621,   549,   550,   620,   550,   620,   200,   621,
     200,   200,   200,   188,   200,   620,   620,   188,   200,   200,
     200,   620,   200,   620,   620,   596,   596,   200,   200,   200,
     200,   620,   200,   200,   200,   620,    69,   200,   350,   620,
     200,   620,   200,   235,   326,   621,   504,   620,   200,   200,
     200,   620,    48,   621,   200,   621,   621,   620,    48,   621,
     200,   200,   200,   621,   621,   621,   621,   621,   200,   621,
     200,   374,   375,   376,   200,   488,   200,   621,   491,   620,
     620,   620,   621,   200,   621,   200,   621,   620,   200,   200,
     200,   200,   200,   200,   200,   200,   200,   200,   200,   200,
     200,   200,   200,   200,   621,   620,   620,   621,   621,   620,
     621,   621,    97,   200,   200,   620,   620,   620,   620,   620,
     620,   620,   621,   200,   200,   621,   621,   200,   200,   200,
     200,   200,   621,   621,   200,   621,   621,   200,   621,   620,
     200,   621,   620,   621,   200,   621,   200,   621,   621,   501,
     620,   501,   501,   501,   501,   200,   620,   200,   620,   621,
     541,   575,   620,   621,   621,   621,   541,   620,   620,   620,
     621,   621,   621,   621,   573,   620,   573,   620,   620,   501,
     620,   620,   200,   621,   621,   144,   200,   315,   620,   621,
     200,   620,   200,   621,   200,   620,   200,   200,   200,   200,
     200,   200,   200,   200,   200,   200,   621,   621,   200,   620,
     200,   200,   200,   200,   621,   298,   200,   620,   200,   621,
      97,   614,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   200,   620,   620,   320,
     620,   200,   621,   621,   200,   351,   200,   200,   200,   200,
     200,   200,   621,   200,   200,   620,   621,   621,   200,   600,
     621,   200,   621,   200,   621,   620,   621,   200,   621,   621,
     621,   620,   554,   620,   621,   200,   200,   597,   621,   200,
     621,   621,   621,   621,   550,   620,   200,   341,   593,   620,
     200,   200,   620,   200,   200,   476,   620,   620,   254,   200,
     200,   200,    91,   200,   621,   200,   621,   200,   200,   200,
     620,   621,   200,   620,   621,   200,   200,   200,   620,   621,
     200,   621,   621,   621,   200,   620,   200,   263,   621,   200,
     374,   200,   374,   200,   501,   620,   200,   620,   621,   541,
     200,   200,   200,   620,   426,   620,   621,   621,   200,   621,
     621,   621,   200,   200,   200,   120,   620,    97,   120,   620,
      97,   621,   200,   621,   621,   621,   620,   621,   293,   620,
     200,   620,   200,   200,   621,   621,   200,   620,   200,   200,
     200,   200,   200,   200,   620,   575,   621,   621,   621,   621,
     293,   620,   621,   621,   621,   200,   621,   621,   200,   620,
     200,   200,   593,   620,   200,   621,   200,   200,   621,   200,
     593,   620,   200,   621,   621,   200,   621,   200,   314,   620,
     200,   621,   621,   200,   200,   620,   200,   200,   621,   620,
     200,   621,   200,   200,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   200,
     620,   293,   620,   200,   200,   200,   621,   621,   621,   200,
     621,   621,   200,   621,   200,   620,   620,   621,   200,   620,
     621,   620,   620,   200,   621,   621,   621,   200,   200,   200,
     200,   621,   620,   200,   200,   200,   219,   221,   222,   223,
     224,   225,   226,   227,   228,   477,   476,   620,   200,   200,
     200,   620,   200,   620,   620,   621,   200,   621,   200,   200,
     621,   200,   200,   620,   620,   620,   200,   200,   200,   621,
     621,   621,   200,   200,   200,   200,   620,   620,   620,   621,
     621,   621,   621,   200,   620,   620,   620,   620,   621,   620,
     200,   620,   200,   621,   200,   620,   621,   200,   200,   200,
     200,   621,   200,   200,   621,   621,   621,   620,   621,    87,
     200,   297,   621,   621,   200,   621,   621,   621,   200,   200,
     621,   621,   621,   200,   200,   621,   200,   621,   200,   621,
     200,   314,   315,   620,   200,   621,   200,   620,   621,   200,
     621,   621,   621,   621,   621,   621,   621,   621,   621,   621,
     621,   621,   621,   621,   621,   293,   620,   620,   621,   621,
     621,   200,    97,   200,   200,   200,   621,   200,   620,   200,
     200,   621,   200,   621,   200,   621,   200,   200,   200,   200,
     620,   621,   621,   476,   476,   254,   200,   200,   620,   621,
     200,   621,   621,   621,   621,   620,   200,   200,   200,   200,
     621,   200,   621,   200,   620,   621,   621,   200,   621,   621,
     120,   620,   120,   620,   620,   200,   200,   620,   200,   621,
     621,   621,   621,   620,   200,   621,   621,   621,   621,   200,
     621,   621,   621,   200,   621,   200,   200,   200,   621,   200,
     200,   621,   621,   200,   621,   200,   200,   230,   200,   341,
     343,   593,   621,   621,   200,   200,   621,   621,   200,   200,
     621,   621,   621,   621,   200,   621,   621,   620,   621,   620,
     200,   621,   621,   200,   200,   200,   621,   621,   621,   200,
     200,   200,   620,   200,   200,   621,   200,   621,   200,   200,
     200,   621,   200,   621,   621,   200,   621,   620,   620,   200,
     621,   200,   621,   621,   621,   621,   621,   621,   200,   621,
     200,   621,   621,   621,   621,   621,   200,   593,   620,   200,
     593,   620,   200,   200,   200,   620,   620,   621,   200,   621,
     621,   200,   621,   200,   621,   621,   200,   621,   621,   200,
     620,   200,   621,   200,   200,    97,   621,   200,   200,   200,
     200,   200,   200,   621,   621,   621,   200,   621,   200,   200,
     200,   200,   621,   200,   621,   200,   621,   621,   621,   621,
     621,   200,   621,   200,   621,   200,   200,   200,   341,   200,
     621,   200,   621,   621,   621,   621,   200,   200,   621,   621,
     200,   200,   200,   621,   621,   200,   621,   621,   621,   200,
     200,   621,   621,   621,   621,   200,   621,   200,   621,   200,
     621,   620,   200,   200,   621,   200,   200,   621,   200,   200,
      97,   621,   621,   621,   200,   200,   621,   620,   621,   621,
     621,   621,   200,   200,   621,   200,   593,   621,   200,   621,
     200,   621,   200,   200,   621,   200,   621,   200,   621,   621,
     621,   621,   621,   200,   621,   200,   200,   621,   200,   200,
     621,   200,   621,   200,   621,   200,   621,   621,   621,   620,
     200,   621,   621,   200,   621,   200,   621,   620,   200,    62,
     200,   621,   621,   200,   620,   621,   621,   621,   200,   200,
      62,   200,   264,   621,   621,   621,   620,   621,   621,   621,
     621,   621,   621,   200,   200,   621,   621,    62,   200,   264,
     621,   621,   620,   621,   621,   621,   621,   200,   200,   200
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
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
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


/* This macro is provided for backward compatibility. */

#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
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
  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  YYSIZE_T yysize1;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = 0;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
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
                yysize1 = yysize + yytnamerr (0, yytname[yyx]);
                if (! (yysize <= yysize1
                       && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                  return 2;
                yysize = yysize1;
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

  yysize1 = yysize + yystrlen (yyformat);
  if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
    return 2;
  yysize = yysize1;

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


/*----------.
| yyparse.  |
`----------*/

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
  if (yypact_value_is_default (yyn))
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

/* Line 1806 of yacc.c  */
#line 150 "p.y"
    { 
          return 0;
         }
    break;

  case 6:

/* Line 1806 of yacc.c  */
#line 161 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; delete (yyvsp[(1) - (1)].bclist); }
    break;

  case 7:

/* Line 1806 of yacc.c  */
#line 163 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 13:

/* Line 1806 of yacc.c  */
#line 170 "p.y"
    {}
    break;

  case 17:

/* Line 1806 of yacc.c  */
#line 175 "p.y"
    {}
    break;

  case 18:

/* Line 1806 of yacc.c  */
#line 177 "p.y"
    {}
    break;

  case 24:

/* Line 1806 of yacc.c  */
#line 184 "p.y"
    {}
    break;

  case 42:

/* Line 1806 of yacc.c  */
#line 203 "p.y"
    { domain->setMFTT((yyvsp[(1) - (1)].mftval)); }
    break;

  case 43:

/* Line 1806 of yacc.c  */
#line 205 "p.y"
    { domain->setMPTT((yyvsp[(1) - (1)].mptval)); }
    break;

  case 44:

/* Line 1806 of yacc.c  */
#line 207 "p.y"
    { domain->setHFTT((yyvsp[(1) - (1)].hftval)); }
    break;

  case 69:

/* Line 1806 of yacc.c  */
#line 234 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 70:

/* Line 1806 of yacc.c  */
#line 236 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 71:

/* Line 1806 of yacc.c  */
#line 238 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 72:

/* Line 1806 of yacc.c  */
#line 240 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 73:

/* Line 1806 of yacc.c  */
#line 242 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 74:

/* Line 1806 of yacc.c  */
#line 244 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 91:

/* Line 1806 of yacc.c  */
#line 262 "p.y"
    { if(domain->setComplexNeuman((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; }
    break;

  case 93:

/* Line 1806 of yacc.c  */
#line 265 "p.y"
    { if(domain->setComplexDirichlet((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; }
    break;

  case 104:

/* Line 1806 of yacc.c  */
#line 277 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 105:

/* Line 1806 of yacc.c  */
#line 279 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; }
    break;

  case 111:

/* Line 1806 of yacc.c  */
#line 286 "p.y"
    {}
    break;

  case 120:

/* Line 1806 of yacc.c  */
#line 297 "p.y"
    {}
    break;

  case 121:

/* Line 1806 of yacc.c  */
#line 299 "p.y"
    {}
    break;

  case 122:

/* Line 1806 of yacc.c  */
#line 301 "p.y"
    {}
    break;

  case 123:

/* Line 1806 of yacc.c  */
#line 303 "p.y"
    {}
    break;

  case 124:

/* Line 1806 of yacc.c  */
#line 305 "p.y"
    {}
    break;

  case 125:

/* Line 1806 of yacc.c  */
#line 307 "p.y"
    {}
    break;

  case 126:

/* Line 1806 of yacc.c  */
#line 309 "p.y"
    {}
    break;

  case 127:

/* Line 1806 of yacc.c  */
#line 311 "p.y"
    {}
    break;

  case 134:

/* Line 1806 of yacc.c  */
#line 321 "p.y"
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[(3) - (5)].ival)); 
            domain->solInfo().nsample = (yyvsp[(4) - (5)].ival);
          }
    break;

  case 135:

/* Line 1806 of yacc.c  */
#line 328 "p.y"
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[(3) - (4)].ival));
          }
    break;

  case 137:

/* Line 1806 of yacc.c  */
#line 335 "p.y"
    { if ((yyvsp[(2) - (5)].ival) == OutputInfo::Attribute)  geoSource->setGroupAttribute((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival)-1);
          else if ((yyvsp[(2) - (5)].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[(2) - (5)].ival));  exit(-1); }
        }
    break;

  case 138:

/* Line 1806 of yacc.c  */
#line 340 "p.y"
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
        }
    break;

  case 139:

/* Line 1806 of yacc.c  */
#line 352 "p.y"
    { if ((yyvsp[(2) - (6)].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[(2) - (6)].ival));  exit(-1); }
        }
    break;

  case 141:

/* Line 1806 of yacc.c  */
#line 359 "p.y"
    { geoSource->setGroupRandomProperty((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].rprop),(yyvsp[(4) - (6)].fval),(yyvsp[(5) - (6)].fval)); }
    break;

  case 142:

/* Line 1806 of yacc.c  */
#line 363 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (5)].fval)); }
    break;

  case 143:

/* Line 1806 of yacc.c  */
#line 367 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies1(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); }
    break;

  case 144:

/* Line 1806 of yacc.c  */
#line 369 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies2(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); }
    break;

  case 145:

/* Line 1806 of yacc.c  */
#line 371 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (8)].fval)); domain->addFrequencies(2.0*PI*(yyvsp[(4) - (8)].fval), 2.0*PI*(yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].ival), (yyvsp[(7) - (8)].ival)); }
    break;

  case 151:

/* Line 1806 of yacc.c  */
#line 380 "p.y"
    { domain->solInfo().pade_pivot = true; domain->solInfo().pade_tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 152:

/* Line 1806 of yacc.c  */
#line 384 "p.y"
    { domain->solInfo().pade_poles = true; }
    break;

  case 153:

/* Line 1806 of yacc.c  */
#line 386 "p.y"
    { domain->solInfo().pade_poles = true; 
          domain->solInfo().pade_poles_sigmaL = (yyvsp[(2) - (4)].fval); domain->solInfo().pade_poles_sigmaU = (yyvsp[(3) - (4)].fval); }
    break;

  case 154:

/* Line 1806 of yacc.c  */
#line 391 "p.y"
    { geoSource->setImpe((yyvsp[(2) - (3)].fval)); domain->addCoarseFrequency(2.0*PI*(yyvsp[(2) - (3)].fval)); }
    break;

  case 155:

/* Line 1806 of yacc.c  */
#line 394 "p.y"
    { domain->addFrequencies(2.0*PI*(yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].ival)); }
    break;

  case 156:

/* Line 1806 of yacc.c  */
#line 398 "p.y"
    { domain->solInfo().freqSweepMethod = (yyvsp[(2) - (4)].ival); 
          int &l = domain->solInfo().padeL, &m = domain->solInfo().padeM, &n = domain->solInfo().padeN;
          switch((yyvsp[(2) - (4)].ival)) {
            case SolverInfo::Taylor:
              domain->solInfo().nFreqSweepRHS = (yyvsp[(3) - (4)].ival)+1; // taylor
              break;
            case SolverInfo::Pade1:
              n = 1;
              domain->solInfo().nFreqSweepRHS = l+m+1;
              break;
            case SolverInfo::Pade:
            case SolverInfo::Fourier:
              n = (yyvsp[(3) - (4)].ival);
              domain->solInfo().nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SolverInfo::PadeLanczos:
              n = (yyvsp[(3) - (4)].ival);
              if(m%n != 0) m = m/n*(n+1)-m%n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().nFreqSweepRHS = m/n;
              break;
            case SolverInfo::GalProjection:
              n = (yyvsp[(3) - (4)].ival);
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::KrylovGalProjection:
              n = (yyvsp[(3) - (4)].ival);
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::QRGalProjection:
              n = (yyvsp[(3) - (4)].ival);
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
          }
        }
    break;

  case 157:

/* Line 1806 of yacc.c  */
#line 437 "p.y"
    { domain->solInfo().freqSweepMethod = (yyvsp[(2) - (6)].ival);
          int &l = domain->solInfo().padeL, &m = domain->solInfo().padeM, &n = domain->solInfo().padeN;
          switch((yyvsp[(2) - (6)].ival)) {
            case SolverInfo::Taylor:
              domain->solInfo().nFreqSweepRHS = (yyvsp[(3) - (6)].ival)+1; // taylor
              break;
            case SolverInfo::Pade1:
              n = 1;
              l = (yyvsp[(4) - (6)].ival);
              m = (yyvsp[(5) - (6)].ival);
              domain->solInfo().nFreqSweepRHS = l+m+1;
              break;
            case SolverInfo::Pade:
            case SolverInfo::Fourier:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival); 
              m = (yyvsp[(5) - (6)].ival);
              domain->solInfo().nFreqSweepRHS = (int) ceil(float(l+m+1)/float(n));
              break;
            case SolverInfo::PadeLanczos:
              n = (yyvsp[(3) - (6)].ival);
              m = (yyvsp[(5) - (6)].ival);
              if(m%n != 0) m = m/n*(n+1)-m%n; // round m up to the nearest multiple of n
              l = m-1;
              domain->solInfo().nFreqSweepRHS = m/n;
              break;
            case SolverInfo::GalProjection:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival);
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::KrylovGalProjection:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival);
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
            case SolverInfo::QRGalProjection:
              n = (yyvsp[(3) - (6)].ival);
              l = (yyvsp[(4) - (6)].ival);
              m = 1;
              domain->solInfo().nFreqSweepRHS = l+1;
              break;
          }
        }
    break;

  case 158:

/* Line 1806 of yacc.c  */
#line 486 "p.y"
    { geoSource->binaryInput = true; geoSource->binaryOutput = true; }
    break;

  case 159:

/* Line 1806 of yacc.c  */
#line 488 "p.y"
    { geoSource->binaryInput = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 160:

/* Line 1806 of yacc.c  */
#line 490 "p.y"
    { geoSource->binaryOutput = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 162:

/* Line 1806 of yacc.c  */
#line 495 "p.y"
    { geoSource->setGeo((yyvsp[(3) - (4)].strval)); }
    break;

  case 163:

/* Line 1806 of yacc.c  */
#line 497 "p.y"
    { geoSource->setDecomp((yyvsp[(3) - (4)].strval)); }
    break;

  case 164:

/* Line 1806 of yacc.c  */
#line 499 "p.y"
    { geoSource->setGlob((yyvsp[(3) - (4)].strval)); }
    break;

  case 165:

/* Line 1806 of yacc.c  */
#line 501 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); }
    break;

  case 166:

/* Line 1806 of yacc.c  */
#line 503 "p.y"
    { geoSource->setCpuMap((yyvsp[(3) - (4)].strval)); }
    break;

  case 167:

/* Line 1806 of yacc.c  */
#line 507 "p.y"
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[(2) - (3)].ival)); 
#endif
	}
    break;

  case 168:

/* Line 1806 of yacc.c  */
#line 515 "p.y"
    {if(decInit==0) decInit = new DecInit(); }
    break;

  case 169:

/* Line 1806 of yacc.c  */
#line 517 "p.y"
    {decInit->file = strdup((yyvsp[(3) - (4)].strval));}
    break;

  case 170:

/* Line 1806 of yacc.c  */
#line 519 "p.y"
    {decInit->nsubs = (yyvsp[(3) - (4)].ival); }
    break;

  case 171:

/* Line 1806 of yacc.c  */
#line 521 "p.y"
    {decInit->weight = true; }
    break;

  case 172:

/* Line 1806 of yacc.c  */
#line 523 "p.y"
    {decInit->memory = true; }
    break;

  case 173:

/* Line 1806 of yacc.c  */
#line 525 "p.y"
    {decInit->exitAfterDec = true;}
    break;

  case 174:

/* Line 1806 of yacc.c  */
#line 527 "p.y"
    {decInit->skip = true;}
    break;

  case 175:

/* Line 1806 of yacc.c  */
#line 529 "p.y"
    {decInit->nosa = true; }
    break;

  case 176:

/* Line 1806 of yacc.c  */
#line 533 "p.y"
    {}
    break;

  case 177:

/* Line 1806 of yacc.c  */
#line 535 "p.y"
    {
	   // map<int,double >::iterator it = weightList.find($2);
	   //if(it == weightList.end())
	     weightList[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval);
	 }
    break;

  case 178:

/* Line 1806 of yacc.c  */
#line 543 "p.y"
    { (yyval.mftval) = new MFTTData; }
    break;

  case 179:

/* Line 1806 of yacc.c  */
#line 545 "p.y"
    { (yyval.mftval) = (yyvsp[(1) - (4)].mftval); (yyval.mftval)->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); }
    break;

  case 180:

/* Line 1806 of yacc.c  */
#line 549 "p.y"
    { (yyval.mptval) = new MFTTData; }
    break;

  case 181:

/* Line 1806 of yacc.c  */
#line 551 "p.y"
    { (yyval.mptval) = (yyvsp[(1) - (4)].mptval); (yyval.mptval)->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); }
    break;

  case 182:

/* Line 1806 of yacc.c  */
#line 555 "p.y"
    { (yyval.hftval) = new MFTTData; }
    break;

  case 183:

/* Line 1806 of yacc.c  */
#line 557 "p.y"
    { (yyval.hftval) = (yyvsp[(1) - (4)].hftval); (yyval.hftval)->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); }
    break;

  case 191:

/* Line 1806 of yacc.c  */
#line 570 "p.y"
    { geoSource->addCFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); }
    break;

  case 192:

/* Line 1806 of yacc.c  */
#line 574 "p.y"
    { geoSource->addCoefInfo((yyvsp[(2) - (4)].ival)-1,(yyvsp[(4) - (4)].coefdata)); }
    break;

  case 193:

/* Line 1806 of yacc.c  */
#line 578 "p.y"
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[(1) - (4)].ival)-1,(yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].fval)); }
    break;

  case 194:

/* Line 1806 of yacc.c  */
#line 580 "p.y"
    { (yyval.coefdata).setCoef((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); }
    break;

  case 195:

/* Line 1806 of yacc.c  */
#line 584 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); }
    break;

  case 196:

/* Line 1806 of yacc.c  */
#line 586 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); }
    break;

  case 197:

/* Line 1806 of yacc.c  */
#line 590 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); }
    break;

  case 198:

/* Line 1806 of yacc.c  */
#line 592 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); }
    break;

  case 199:

/* Line 1806 of yacc.c  */
#line 596 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); }
    break;

  case 200:

/* Line 1806 of yacc.c  */
#line 598 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); }
    break;

  case 201:

/* Line 1806 of yacc.c  */
#line 602 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); }
    break;

  case 202:

/* Line 1806 of yacc.c  */
#line 604 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); }
    break;

  case 203:

/* Line 1806 of yacc.c  */
#line 608 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (11)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
	  (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval);
	  (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); }
    break;

  case 204:

/* Line 1806 of yacc.c  */
#line 614 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (13)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (13)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (13)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (13)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (13)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (13)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (13)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (13)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (13)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (13)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (13)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (13)].fval); }
    break;

  case 205:

/* Line 1806 of yacc.c  */
#line 621 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (14)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (14)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (14)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (14)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (14)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (14)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (14)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (14)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (14)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (14)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (14)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (14)].fval); (yyval.ldata).d[11] = (yyvsp[(13) - (14)].fval); }
    break;

  case 206:

/* Line 1806 of yacc.c  */
#line 630 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (5)].ival)-1;  (yyval.ldata).matid = (yyvsp[(2) - (5)].ival)-1; (yyval.ldata).d[7] = (yyvsp[(3) - (5)].fval); (yyval.ldata).d[8] = (yyvsp[(4) - (5)].fval); }
    break;

  case 208:

/* Line 1806 of yacc.c  */
#line 635 "p.y"
    { geoSource->addLayMat((yyvsp[(2) - (2)].ldata).matid, (yyvsp[(2) - (2)].ldata).d); }
    break;

  case 209:

/* Line 1806 of yacc.c  */
#line 642 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (7)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (7)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (7)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (7)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (7)].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[(6) - (7)].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
    break;

  case 210:

/* Line 1806 of yacc.c  */
#line 648 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (9)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (9)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (9)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (9)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (9)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (9)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (9)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (9)].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; }
    break;

  case 211:

/* Line 1806 of yacc.c  */
#line 653 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (11)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval);
          (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); (yyval.ldata).d[9] = 0; }
    break;

  case 212:

/* Line 1806 of yacc.c  */
#line 657 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (12)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (12)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (12)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (12)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (12)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (12)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (12)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (12)].fval); 
          (yyval.ldata).d[7] = (yyvsp[(9) - (12)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (12)].fval); (yyval.ldata).d[9] = (yyvsp[(11) - (12)].fval); }
    break;

  case 214:

/* Line 1806 of yacc.c  */
#line 664 "p.y"
    { domain->addDMass((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); }
    break;

  case 215:

/* Line 1806 of yacc.c  */
#line 666 "p.y"
    { domain->addDMass((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(5) - (6)].fval),(yyvsp[(4) - (6)].ival)-1); }
    break;

  case 217:

/* Line 1806 of yacc.c  */
#line 671 "p.y"
    { domain->setGravity((yyvsp[(2) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); }
    break;

  case 219:

/* Line 1806 of yacc.c  */
#line 676 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (5)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (5)].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[(4) - (5)].strval); }
    break;

  case 220:

/* Line 1806 of yacc.c  */
#line 680 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (4)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (4)].strval);}
    break;

  case 221:

/* Line 1806 of yacc.c  */
#line 683 "p.y"
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[(2) - (4)].strval);
          domain->solInfo().nRestart = (yyvsp[(3) - (4)].ival); }
    break;

  case 222:

/* Line 1806 of yacc.c  */
#line 688 "p.y"
    { geoSource->setControlFile((yyvsp[(2) - (3)].strval));
         geoSource->setControlRoutine((char *) "controlObj");}
    break;

  case 223:

/* Line 1806 of yacc.c  */
#line 693 "p.y"
    { geoSource->setControlRoutine((yyvsp[(2) - (3)].strval)); }
    break;

  case 224:

/* Line 1806 of yacc.c  */
#line 697 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; }
    break;

  case 225:

/* Line 1806 of yacc.c  */
#line 702 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Actuators;
          if(geoSource->setActuatorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)            return -1; }
    break;

  case 226:

/* Line 1806 of yacc.c  */
#line 708 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdf;
          if(geoSource->setUsdfLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)       return -1; }
    break;

  case 227:

/* Line 1806 of yacc.c  */
#line 715 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdd;
          if(geoSource->setUsddLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)    return -1; }
    break;

  case 228:

/* Line 1806 of yacc.c  */
#line 722 "p.y"
    { numColumns = 3; }
    break;

  case 229:

/* Line 1806 of yacc.c  */
#line 724 "p.y"
    { numColumns = 6; }
    break;

  case 230:

/* Line 1806 of yacc.c  */
#line 726 "p.y"
    { numColumns = 3; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); }
    break;

  case 231:

/* Line 1806 of yacc.c  */
#line 728 "p.y"
    { numColumns = 6; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); }
    break;

  case 232:

/* Line 1806 of yacc.c  */
#line 730 "p.y"
    { (yyvsp[(2) - (3)].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[(2) - (3)].oinfo)); }
    break;

  case 233:

/* Line 1806 of yacc.c  */
#line 734 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); }
    break;

  case 234:

/* Line 1806 of yacc.c  */
#line 736 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); }
    break;

  case 235:

/* Line 1806 of yacc.c  */
#line 738 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; }
    break;

  case 236:

/* Line 1806 of yacc.c  */
#line 740 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); 
          if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1;}
    break;

  case 237:

/* Line 1806 of yacc.c  */
#line 743 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; }
    break;

  case 238:

/* Line 1806 of yacc.c  */
#line 745 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; }
    break;

  case 239:

/* Line 1806 of yacc.c  */
#line 747 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); }
    break;

  case 240:

/* Line 1806 of yacc.c  */
#line 749 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); }
    break;

  case 241:

/* Line 1806 of yacc.c  */
#line 751 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; }
    break;

  case 242:

/* Line 1806 of yacc.c  */
#line 753 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1; }
    break;

  case 243:

/* Line 1806 of yacc.c  */
#line 755 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; }
    break;

  case 244:

/* Line 1806 of yacc.c  */
#line 757 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; }
    break;

  case 245:

/* Line 1806 of yacc.c  */
#line 760 "p.y"
    { (yyval.oinfo).nodeNumber = (yyvsp[(3) - (3)].ival)-1; }
    break;

  case 246:

/* Line 1806 of yacc.c  */
#line 762 "p.y"
    { (yyval.oinfo).surface = (yyvsp[(2) - (2)].ival); }
    break;

  case 247:

/* Line 1806 of yacc.c  */
#line 764 "p.y"
    { (yyval.oinfo).ylayer = (yyvsp[(2) - (3)].fval); (yyval.oinfo).zlayer = (yyvsp[(3) - (3)].fval); }
    break;

  case 248:

/* Line 1806 of yacc.c  */
#line 766 "p.y"
    { (yyval.oinfo).averageFlg = (yyvsp[(2) - (2)].ival); }
    break;

  case 249:

/* Line 1806 of yacc.c  */
#line 768 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (2)].ival); }
    break;

  case 250:

/* Line 1806 of yacc.c  */
#line 770 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (3)].ival); (yyval.oinfo).ncomplexout = (yyvsp[(3) - (3)].ival); }
    break;

  case 251:

/* Line 1806 of yacc.c  */
#line 772 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (2)].ival); }
    break;

  case 252:

/* Line 1806 of yacc.c  */
#line 774 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (3)].ival); sfem->setnsamp_out((yyvsp[(3) - (3)].ival)); }
    break;

  case 253:

/* Line 1806 of yacc.c  */
#line 776 "p.y"
    { (yyval.oinfo).matlab = true; }
    break;

  case 254:

/* Line 1806 of yacc.c  */
#line 780 "p.y"
    { domain->outFlag = (yyvsp[(2) - (3)].ival); }
    break;

  case 256:

/* Line 1806 of yacc.c  */
#line 785 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; }
    break;

  case 257:

/* Line 1806 of yacc.c  */
#line 788 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[(2) - (3)].ival);}
    break;

  case 258:

/* Line 1806 of yacc.c  */
#line 791 "p.y"
    { domain->solInfo().nEig = (yyvsp[(2) - (3)].ival); }
    break;

  case 259:

/* Line 1806 of yacc.c  */
#line 793 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;}
    break;

  case 260:

/* Line 1806 of yacc.c  */
#line 795 "p.y"
    { domain->solInfo().setSubSpaceInfo((yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); }
    break;

  case 261:

/* Line 1806 of yacc.c  */
#line 797 "p.y"
    { domain->solInfo().subspaceSize = (yyvsp[(2) - (3)].ival);}
    break;

  case 262:

/* Line 1806 of yacc.c  */
#line 799 "p.y"
    { domain->solInfo().tolEig = (yyvsp[(2) - (3)].fval); }
    break;

  case 263:

/* Line 1806 of yacc.c  */
#line 801 "p.y"
    { domain->solInfo().tolJac = (yyvsp[(2) - (3)].fval); }
    break;

  case 264:

/* Line 1806 of yacc.c  */
#line 803 "p.y"
    { domain->solInfo().explicitK = true; }
    break;

  case 265:

/* Line 1806 of yacc.c  */
#line 805 "p.y"
    { geoSource->setShift((yyvsp[(2) - (3)].fval)); }
    break;

  case 266:

/* Line 1806 of yacc.c  */
#line 807 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; }
    break;

  case 267:

/* Line 1806 of yacc.c  */
#line 809 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (3)].strval); }
    break;

  case 268:

/* Line 1806 of yacc.c  */
#line 812 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (4)].strval); 
          domain->solInfo().arpack_mode = (yyvsp[(3) - (4)].ival); }
    break;

  case 269:

/* Line 1806 of yacc.c  */
#line 816 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[(2) - (4)].fval), int((yyvsp[(3) - (4)].fval))); }
    break;

  case 270:

/* Line 1806 of yacc.c  */
#line 819 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival));}
    break;

  case 271:

/* Line 1806 of yacc.c  */
#line 822 "p.y"
    { domain->solInfo().filtereig = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 272:

/* Line 1806 of yacc.c  */
#line 824 "p.y"
    { domain->solInfo().eigenSolverSubType = (yyvsp[(2) - (3)].ival); }
    break;

  case 273:

/* Line 1806 of yacc.c  */
#line 826 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;}
    break;

  case 274:

/* Line 1806 of yacc.c  */
#line 829 "p.y"
    { domain->solInfo().maxitEig = (yyvsp[(3) - (4)].ival); }
    break;

  case 275:

/* Line 1806 of yacc.c  */
#line 831 "p.y"
    { domain->solInfo().test_ulrich = true; }
    break;

  case 276:

/* Line 1806 of yacc.c  */
#line 833 "p.y"
    { domain->solInfo().addedMass = (yyvsp[(2) - (3)].ival); }
    break;

  case 277:

/* Line 1806 of yacc.c  */
#line 837 "p.y"
    { domain->solInfo().sloshing = 1; }
    break;

  case 278:

/* Line 1806 of yacc.c  */
#line 839 "p.y"
    { domain->setGravitySloshing((yyvsp[(2) - (3)].fval)); }
    break;

  case 279:

/* Line 1806 of yacc.c  */
#line 843 "p.y"
    { domain->solInfo().massFlag = 1; }
    break;

  case 280:

/* Line 1806 of yacc.c  */
#line 847 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); }
    break;

  case 281:

/* Line 1806 of yacc.c  */
#line 850 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);}
    break;

  case 282:

/* Line 1806 of yacc.c  */
#line 854 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Top); }
    break;

  case 283:

/* Line 1806 of yacc.c  */
#line 858 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); }
    break;

  case 287:

/* Line 1806 of yacc.c  */
#line 863 "p.y"
    { domain->solInfo().modal = true; }
    break;

  case 288:

/* Line 1806 of yacc.c  */
#line 865 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); }
    break;

  case 289:

/* Line 1806 of yacc.c  */
#line 867 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); }
    break;

  case 290:

/* Line 1806 of yacc.c  */
#line 869 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (8)].ival);
          domain->solInfo().stable_cfl = (yyvsp[(4) - (8)].fval);
          domain->solInfo().stable_tol = (yyvsp[(5) - (8)].fval);
          domain->solInfo().stable_maxit = (yyvsp[(6) - (8)].ival);
          domain->solInfo().stable_freq = (yyvsp[(7) - (8)].ival);
        }
    break;

  case 291:

/* Line 1806 of yacc.c  */
#line 876 "p.y"
    { domain->solInfo().iacc_switch = bool((yyvsp[(3) - (4)].ival)); }
    break;

  case 292:

/* Line 1806 of yacc.c  */
#line 878 "p.y"
    { domain->solInfo().zeroRot = bool((yyvsp[(3) - (4)].ival)); }
    break;

  case 293:

/* Line 1806 of yacc.c  */
#line 880 "p.y"
    { domain->solInfo().no_secondary = true; }
    break;

  case 294:

/* Line 1806 of yacc.c  */
#line 882 "p.y"
    { domain->solInfo().check_energy_balance = true; }
    break;

  case 295:

/* Line 1806 of yacc.c  */
#line 884 "p.y"
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[(3) - (5)].fval); 
          domain->solInfo().epsilon2 = (yyvsp[(4) - (5)].fval); }
    break;

  case 296:

/* Line 1806 of yacc.c  */
#line 890 "p.y"
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; }
    break;

  case 298:

/* Line 1806 of yacc.c  */
#line 893 "p.y"
    { domain->solInfo().acoustic = true; }
    break;

  case 300:

/* Line 1806 of yacc.c  */
#line 898 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (2)].fval),(yyvsp[(2) - (2)].fval)); }
    break;

  case 301:

/* Line 1806 of yacc.c  */
#line 900 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (4)].fval),(yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval),(yyvsp[(4) - (4)].fval)); }
    break;

  case 302:

/* Line 1806 of yacc.c  */
#line 902 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[(1) - (1)].fval)); }
    break;

  case 303:

/* Line 1806 of yacc.c  */
#line 904 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (3)].fval),(yyvsp[(2) - (3)].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[(3) - (3)].fval); }
    break;

  case 304:

/* Line 1806 of yacc.c  */
#line 910 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
          domain->solInfo().order = 1;
         }
          else 
          domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo((yyvsp[(1) - (1)].fval)); 
        }
    break;

  case 305:

/* Line 1806 of yacc.c  */
#line 921 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; }
    break;

  case 308:

/* Line 1806 of yacc.c  */
#line 926 "p.y"
    { domain->solInfo().modal = true; }
    break;

  case 309:

/* Line 1806 of yacc.c  */
#line 930 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (5)].fval), 0, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); }
    break;

  case 310:

/* Line 1806 of yacc.c  */
#line 932 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), 0, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].ival), (yyvsp[(5) - (6)].fval)); }
    break;

  case 311:

/* Line 1806 of yacc.c  */
#line 936 "p.y"
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].ival)); }
    break;

  case 312:

/* Line 1806 of yacc.c  */
#line 946 "p.y"
    { domain->solInfo().setAero((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCollocated = 0; }
    break;

  case 313:

/* Line 1806 of yacc.c  */
#line 949 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().isCollocated = 0; }
    break;

  case 314:

/* Line 1806 of yacc.c  */
#line 952 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (6)].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[(3) - (6)].ival) < 6 || (yyvsp[(3) - (6)].ival) == 20) {
              domain->solInfo().alphas[0] = (yyvsp[(4) - (6)].fval)+(yyvsp[(5) - (6)].fval);
              domain->solInfo().alphas[1] = -(yyvsp[(5) - (6)].fval);
          }
        }
    break;

  case 315:

/* Line 1806 of yacc.c  */
#line 960 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (5)].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[(4) - (5)].fval);
        }
    break;

  case 316:

/* Line 1806 of yacc.c  */
#line 965 "p.y"
    { domain->solInfo().isCollocated = (yyvsp[(2) - (3)].ival); }
    break;

  case 317:

/* Line 1806 of yacc.c  */
#line 967 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); }
    break;

  case 318:

/* Line 1806 of yacc.c  */
#line 969 "p.y"
    {}
    break;

  case 319:

/* Line 1806 of yacc.c  */
#line 973 "p.y"
    {}
    break;

  case 320:

/* Line 1806 of yacc.c  */
#line 975 "p.y"
    { domain->AddAeroEmbedSurfaceId((yyvsp[(2) - (2)].ival)); }
    break;

  case 321:

/* Line 1806 of yacc.c  */
#line 979 "p.y"
    { domain->solInfo().setAeroHeat((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); }
    break;

  case 322:

/* Line 1806 of yacc.c  */
#line 983 "p.y"
    { domain->solInfo().setThermoh(1); }
    break;

  case 323:

/* Line 1806 of yacc.c  */
#line 987 "p.y"
    { domain->solInfo().setThermoe(1); }
    break;

  case 324:

/* Line 1806 of yacc.c  */
#line 991 "p.y"
    { domain->solInfo().setModeDecomp(1); }
    break;

  case 325:

/* Line 1806 of yacc.c  */
#line 995 "p.y"
    { domain->solInfo().hzemFlag=1; }
    break;

  case 326:

/* Line 1806 of yacc.c  */
#line 999 "p.y"
    { domain->solInfo().slzemFlag=1; }
    break;

  case 327:

/* Line 1806 of yacc.c  */
#line 1003 "p.y"
    { domain->solInfo().setTrbm((yyvsp[(3) - (4)].fval)); }
    break;

  case 328:

/* Line 1806 of yacc.c  */
#line 1011 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");}
    break;

  case 329:

/* Line 1806 of yacc.c  */
#line 1014 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (4)].fval)); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");}
    break;

  case 330:

/* Line 1806 of yacc.c  */
#line 1017 "p.y"
    { domain->solInfo().setGrbm();
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");}
    break;

  case 331:

/* Line 1806 of yacc.c  */
#line 1022 "p.y"
    { domain->solInfo().modeFilterFlag = (yyvsp[(2) - (3)].ival); }
    break;

  case 332:

/* Line 1806 of yacc.c  */
#line 1024 "p.y"
    { domain->solInfo().modeFilterFlag = 1; }
    break;

  case 333:

/* Line 1806 of yacc.c  */
#line 1028 "p.y"
    { domain->solInfo().useRbmFilter((yyvsp[(2) - (3)].ival)); }
    break;

  case 334:

/* Line 1806 of yacc.c  */
#line 1030 "p.y"
    { domain->solInfo().useRbmFilter(1); }
    break;

  case 336:

/* Line 1806 of yacc.c  */
#line 1035 "p.y"
    { if((yyvsp[(1) - (1)].ival) < 1 || (yyvsp[(1) - (1)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(1) - (1)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(1) - (1)].ival)-1] = 1;
    }
    break;

  case 337:

/* Line 1806 of yacc.c  */
#line 1043 "p.y"
    { if((yyvsp[(2) - (2)].ival) < 1 || (yyvsp[(2) - (2)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(2) - (2)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(2) - (2)].ival)-1] = 1;
    }
    break;

  case 338:

/* Line 1806 of yacc.c  */
#line 1053 "p.y"
    { domain->solInfo().hzemFilterFlag=1; }
    break;

  case 339:

/* Line 1806 of yacc.c  */
#line 1057 "p.y"
    { domain->solInfo().slzemFilterFlag=1; }
    break;

  case 340:

/* Line 1806 of yacc.c  */
#line 1061 "p.y"
    { domain->solInfo().setTimes((yyvsp[(4) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(2) - (5)].fval)); }
    break;

  case 341:

/* Line 1806 of yacc.c  */
#line 1065 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (5)].ival),(yyvsp[(4) - (5)].ival),1);
        }
    break;

  case 342:

/* Line 1806 of yacc.c  */
#line 1071 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].ival));
        }
    break;

  case 343:

/* Line 1806 of yacc.c  */
#line 1076 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (7)].ival),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        }
    break;

  case 346:

/* Line 1806 of yacc.c  */
#line 1089 "p.y"
    { domain->solInfo().pitaNoForce = true; }
    break;

  case 347:

/* Line 1806 of yacc.c  */
#line 1091 "p.y"
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[(2) - (2)].ival); }
    break;

  case 348:

/* Line 1806 of yacc.c  */
#line 1093 "p.y"
    { domain->solInfo().pitaLocalBasisImprovement = 1; }
    break;

  case 349:

/* Line 1806 of yacc.c  */
#line 1095 "p.y"
    { domain->solInfo().pitaTimeReversible = true; }
    break;

  case 350:

/* Line 1806 of yacc.c  */
#line 1097 "p.y"
    { domain->solInfo().pitaRemoteCoarse = true; }
    break;

  case 351:

/* Line 1806 of yacc.c  */
#line 1099 "p.y"
    { domain->solInfo().pitaProjTol = (yyvsp[(2) - (2)].fval); }
    break;

  case 352:

/* Line 1806 of yacc.c  */
#line 1101 "p.y"
    { domain->solInfo().pitaReadInitSeed = true; }
    break;

  case 353:

/* Line 1806 of yacc.c  */
#line 1103 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = 0.0; }
    break;

  case 354:

/* Line 1806 of yacc.c  */
#line 1105 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[(2) - (2)].fval); }
    break;

  case 355:

/* Line 1806 of yacc.c  */
#line 1107 "p.y"
    { domain->solInfo().pitaJumpMagnOutput = true; }
    break;

  case 356:

/* Line 1806 of yacc.c  */
#line 1111 "p.y"
    { domain->solInfo().setDamping((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); }
    break;

  case 357:

/* Line 1806 of yacc.c  */
#line 1113 "p.y"
    { if(geoSource->setModalDamping((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
    break;

  case 358:

/* Line 1806 of yacc.c  */
#line 1118 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); }
    break;

  case 359:

/* Line 1806 of yacc.c  */
#line 1122 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        }
    break;

  case 360:

/* Line 1806 of yacc.c  */
#line 1128 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval), 0.0);
        }
    break;

  case 361:

/* Line 1806 of yacc.c  */
#line 1136 "p.y"
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        }
    break;

  case 362:

/* Line 1806 of yacc.c  */
#line 1142 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[(2) - (3)].ival),0.0,0.0,0.0);
        }
    break;

  case 364:

/* Line 1806 of yacc.c  */
#line 1150 "p.y"
    {
           domain->setWaveDirections((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        }
    break;

  case 365:

/* Line 1806 of yacc.c  */
#line 1156 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        }
    break;

  case 366:

/* Line 1806 of yacc.c  */
#line 1164 "p.y"
    { (yyval.bclist) = new BCList; }
    break;

  case 367:

/* Line 1806 of yacc.c  */
#line 1166 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 368:

/* Line 1806 of yacc.c  */
#line 1168 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
    break;

  case 369:

/* Line 1806 of yacc.c  */
#line 1170 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } }
    break;

  case 370:

/* Line 1806 of yacc.c  */
#line 1172 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
    break;

  case 372:

/* Line 1806 of yacc.c  */
#line 1187 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(2) - (6)].ival)-1;
          surf_bc[0].type = (BCond::BCType) (yyvsp[(3) - (6)].ival); //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = (yyvsp[(4) - (6)].ival)-1;
          surf_bc[0].val = (yyvsp[(5) - (6)].ival)-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
        }
    break;

  case 373:

/* Line 1806 of yacc.c  */
#line 1196 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); }
    break;

  case 374:

/* Line 1806 of yacc.c  */
#line 1200 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); }
    break;

  case 375:

/* Line 1806 of yacc.c  */
#line 1202 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 376:

/* Line 1806 of yacc.c  */
#line 1206 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); }
    break;

  case 377:

/* Line 1806 of yacc.c  */
#line 1208 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (2)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; }
    break;

  case 378:

/* Line 1806 of yacc.c  */
#line 1212 "p.y"
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; 

          int* allPDirNodes = new int[(yyvsp[(3) - (3)].bclist)->n];

          for (int ii=0; ii < (yyvsp[(3) - (3)].bclist)->n; ii++)
            allPDirNodes[ii]=((yyvsp[(3) - (3)].bclist)->d[ii]).nnum;
          sort(allPDirNodes,allPDirNodes + (yyvsp[(3) - (3)].bclist)->n);

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
        }
    break;

  case 379:

/* Line 1806 of yacc.c  */
#line 1253 "p.y"
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[(1) - (1)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(1) - (1)].bclist)->d)[ii]); }
    break;

  case 380:

/* Line 1806 of yacc.c  */
#line 1257 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); 
          for (int ii = 0; ii < (yyvsp[(2) - (2)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(2) - (2)].bclist)->d)[ii]); }
    break;

  case 381:

/* Line 1806 of yacc.c  */
#line 1263 "p.y"
    { (yyval.bclist) = new BCList;
          //BCond *bclist = new BCond[$3.num];
          for(int i=0; i<(yyvsp[(3) - (4)].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[(3) - (4)].nl).nd[i],10,0.0); } }
    break;

  case 382:

/* Line 1806 of yacc.c  */
#line 1271 "p.y"
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;}
    break;

  case 383:

/* Line 1806 of yacc.c  */
#line 1273 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); }
    break;

  case 384:

/* Line 1806 of yacc.c  */
#line 1276 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); }
    break;

  case 385:

/* Line 1806 of yacc.c  */
#line 1285 "p.y"
    { (yyval.bclist) = new BCList; }
    break;

  case 386:

/* Line 1806 of yacc.c  */
#line 1287 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Flux; (yyval.bclist)->add(bc); }
    break;

  case 387:

/* Line 1806 of yacc.c  */
#line 1290 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].type = BCond::Flux;
          geoSource->addSurfaceNeuman(1,surf_bc); }
    break;

  case 388:

/* Line 1806 of yacc.c  */
#line 1299 "p.y"
    { (yyval.bclist) = new BCList; }
    break;

  case 389:

/* Line 1806 of yacc.c  */
#line 1301 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Convection; (yyval.bclist)->add(bc); }
    break;

  case 390:

/* Line 1806 of yacc.c  */
#line 1306 "p.y"
    { (yyval.bclist) = new BCList; }
    break;

  case 391:

/* Line 1806 of yacc.c  */
#line 1308 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Radiation; (yyval.bclist)->add(bc); }
    break;

  case 395:

/* Line 1806 of yacc.c  */
#line 1320 "p.y"
    { domain->addSommer(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); }
    break;

  case 396:

/* Line 1806 of yacc.c  */
#line 1322 "p.y"
    { domain->addSommer(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); }
    break;

  case 397:

/* Line 1806 of yacc.c  */
#line 1324 "p.y"
    { domain->addSommer(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); }
    break;

  case 400:

/* Line 1806 of yacc.c  */
#line 1332 "p.y"
    { domain->addSommerElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        }
    break;

  case 401:

/* Line 1806 of yacc.c  */
#line 1338 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1; }
    break;

  case 402:

/* Line 1806 of yacc.c  */
#line 1340 "p.y"
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++; }
    break;

  case 406:

/* Line 1806 of yacc.c  */
#line 1352 "p.y"
    { domain->addScatter(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); }
    break;

  case 407:

/* Line 1806 of yacc.c  */
#line 1355 "p.y"
    { domain->addScatter(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); }
    break;

  case 408:

/* Line 1806 of yacc.c  */
#line 1358 "p.y"
    { domain->addScatter(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); }
    break;

  case 411:

/* Line 1806 of yacc.c  */
#line 1367 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); }
    break;

  case 414:

/* Line 1806 of yacc.c  */
#line 1376 "p.y"
    { domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); }
    break;

  case 417:

/* Line 1806 of yacc.c  */
#line 1384 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; }
    break;

  case 420:

/* Line 1806 of yacc.c  */
#line 1395 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);}
    break;

  case 421:

/* Line 1806 of yacc.c  */
#line 1399 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); }
    break;

  case 422:

/* Line 1806 of yacc.c  */
#line 1403 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); }
    break;

  case 423:

/* Line 1806 of yacc.c  */
#line 1405 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 424:

/* Line 1806 of yacc.c  */
#line 1409 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); }
    break;

  case 425:

/* Line 1806 of yacc.c  */
#line 1413 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atdneu; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); }
    break;

  case 426:

/* Line 1806 of yacc.c  */
#line 1417 "p.y"
    { domain->solInfo().ATDARBFlag = (yyvsp[(2) - (3)].fval);}
    break;

  case 428:

/* Line 1806 of yacc.c  */
#line 1422 "p.y"
    { domain->solInfo().ATDDNBVal = (yyvsp[(2) - (3)].fval);}
    break;

  case 430:

/* Line 1806 of yacc.c  */
#line 1427 "p.y"
    { domain->solInfo().ATDROBVal = (yyvsp[(2) - (5)].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[(3) - (5)].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[(4) - (5)].fval);}
    break;

  case 432:

/* Line 1806 of yacc.c  */
#line 1434 "p.y"
    { domain->setFFP((yyvsp[(2) - (3)].ival)); }
    break;

  case 433:

/* Line 1806 of yacc.c  */
#line 1436 "p.y"
    { domain->setFFP((yyvsp[(2) - (4)].ival),(yyvsp[(3) - (4)].ival)); }
    break;

  case 434:

/* Line 1806 of yacc.c  */
#line 1440 "p.y"
    {
           domain->setFFP((yyvsp[(2) - (4)].ival));
        }
    break;

  case 435:

/* Line 1806 of yacc.c  */
#line 1446 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        }
    break;

  case 436:

/* Line 1806 of yacc.c  */
#line 1451 "p.y"
    { fourHelmBC->addDirichlet((yyvsp[(2) - (2)].complexFDBC)); }
    break;

  case 437:

/* Line 1806 of yacc.c  */
#line 1455 "p.y"
    { (yyval.complexFDBC) = FDBC((yyvsp[(1) - (2)].ival)-1); }
    break;

  case 438:

/* Line 1806 of yacc.c  */
#line 1459 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        }
    break;

  case 439:

/* Line 1806 of yacc.c  */
#line 1464 "p.y"
    { fourHelmBC->addNeuman((yyvsp[(2) - (2)].complexFNBC)); }
    break;

  case 440:

/* Line 1806 of yacc.c  */
#line 1468 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1); }
    break;

  case 441:

/* Line 1806 of yacc.c  */
#line 1470 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); }
    break;

  case 442:

/* Line 1806 of yacc.c  */
#line 1474 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setModes((yyvsp[(2) - (3)].ival));
          domain->solInfo().setProbType(SolverInfo::AxiHelm);
        }
    break;

  case 443:

/* Line 1806 of yacc.c  */
#line 1482 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSlices((yyvsp[(2) - (3)].ival));
        }
    break;

  case 444:

/* Line 1806 of yacc.c  */
#line 1489 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSomType((yyvsp[(3) - (7)].ival));
          fourHelmBC->setSurf((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval));
        }
    break;

  case 446:

/* Line 1806 of yacc.c  */
#line 1497 "p.y"
    { fourHelmBC->addSommer(new LineAxiSommer((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1)); }
    break;

  case 447:

/* Line 1806 of yacc.c  */
#line 1499 "p.y"
    { fourHelmBC->addSommer(new Line2AxiSommer((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1)); }
    break;

  case 448:

/* Line 1806 of yacc.c  */
#line 1503 "p.y"
    { if( globalMPCs== NULL) globalMPCs = new MPCData(); }
    break;

  case 449:

/* Line 1806 of yacc.c  */
#line 1505 "p.y"
    { globalMPCs->addMPC((yyvsp[(2) - (2)].axiMPC)); }
    break;

  case 450:

/* Line 1806 of yacc.c  */
#line 1509 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (5)].ival)-1, (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival), DComplex(1.0,0.0), 0.0, 0.0, 0.0); }
    break;

  case 451:

/* Line 1806 of yacc.c  */
#line 1511 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (7)].ival)-1, (yyvsp[(2) - (7)].fval), (yyvsp[(3) - (7)].fval), (yyvsp[(4) - (7)].ival), DComplex((yyvsp[(5) - (7)].fval),(yyvsp[(6) - (7)].fval)) , 0.0, 0.0, 0.0); }
    break;

  case 452:

/* Line 1806 of yacc.c  */
#line 1513 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (10)].ival)-1, (yyvsp[(2) - (10)].fval), (yyvsp[(3) - (10)].fval), (yyvsp[(4) - (10)].ival), DComplex((yyvsp[(5) - (10)].fval),(yyvsp[(6) - (10)].fval)) , (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)); }
    break;

  case 453:

/* Line 1806 of yacc.c  */
#line 1517 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (3)].strval);
	  domain->solInfo().readmodeCalled = true; }
    break;

  case 454:

/* Line 1806 of yacc.c  */
#line 1520 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (4)].strval);
          domain->solInfo().readmodeCalled = true; 
 	  domain->solInfo().maxSizePodRom = (yyvsp[(3) - (4)].ival); }
    break;

  case 455:

/* Line 1806 of yacc.c  */
#line 1526 "p.y"
    { }
    break;

  case 456:

/* Line 1806 of yacc.c  */
#line 1528 "p.y"
    { domain->solInfo().zeroInitialDisp = 1; }
    break;

  case 457:

/* Line 1806 of yacc.c  */
#line 1530 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; }
    break;

  case 458:

/* Line 1806 of yacc.c  */
#line 1533 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
    break;

  case 459:

/* Line 1806 of yacc.c  */
#line 1539 "p.y"
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[(2) - (3)].fval);  }
    break;

  case 460:

/* Line 1806 of yacc.c  */
#line 1541 "p.y"
    { (yyval.bclist) = new BCList; amplitude = 1.0; }
    break;

  case 461:

/* Line 1806 of yacc.c  */
#line 1543 "p.y"
    { BCond bc; /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (9)].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[(3) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[(4) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[(5) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = amplitude*(yyvsp[(6) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = amplitude*(yyvsp[(7) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = amplitude*(yyvsp[(8) - (9)].fval); (yyval.bclist)->add(bc);
          for(int i=0; i<(yyval.bclist)->n; ++i) (yyval.bclist)->d[i].type = BCond::Idisp6;
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        }
    break;

  case 462:

/* Line 1806 of yacc.c  */
#line 1554 "p.y"
    { BCond bc; /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 0; bc.val = amplitude*(yyvsp[(3) - (6)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = amplitude*(yyvsp[(4) - (6)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = amplitude*(yyvsp[(5) - (6)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = 0.0         ; (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = 0.0         ; (yyval.bclist)->add(bc);
          for(int i=0; i<(yyval.bclist)->n; ++i) (yyval.bclist)->d[i].type = BCond::Idisp6;
          geoSource->setIDis6((yyval.bclist)->n, (yyval.bclist)->d);
        }
    break;

  case 463:

/* Line 1806 of yacc.c  */
#line 1565 "p.y"
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); }
    break;

  case 464:

/* Line 1806 of yacc.c  */
#line 1568 "p.y"
    { domain->solInfo().buckling = 1; }
    break;

  case 465:

/* Line 1806 of yacc.c  */
#line 1573 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); }
    break;

  case 466:

/* Line 1806 of yacc.c  */
#line 1575 "p.y"
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (9)].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[(3) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[(4) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[(5) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[(6) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[(7) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[(8) - (9)].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIDis6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        }
    break;

  case 467:

/* Line 1806 of yacc.c  */
#line 1588 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); }
    break;

  case 468:

/* Line 1806 of yacc.c  */
#line 1590 "p.y"
    { BCond bc;                          /* add 6 boundary conditions */
          bc.nnum = (yyvsp[(2) - (9)].ival)-1; bc.dofnum = 0; bc.val = (yyvsp[(3) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 1; bc.val = (yyvsp[(4) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 2; bc.val = (yyvsp[(5) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 3; bc.val = (yyvsp[(6) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 4; bc.val = (yyvsp[(7) - (9)].fval); (yyval.bclist)->add(bc);
                          bc.dofnum = 5; bc.val = (yyvsp[(8) - (9)].fval); (yyval.bclist)->add(bc);
          geoSource->setPitaIVel6((yyval.bclist)->n, (yyval.bclist)->d, PitaTS);
        }
    break;

  case 469:

/* Line 1806 of yacc.c  */
#line 1602 "p.y"
    { }
    break;

  case 470:

/* Line 1806 of yacc.c  */
#line 1604 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; }
    break;

  case 471:

/* Line 1806 of yacc.c  */
#line 1607 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; }
    break;

  case 472:

/* Line 1806 of yacc.c  */
#line 1613 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; }
    break;

  case 473:

/* Line 1806 of yacc.c  */
#line 1618 "p.y"
    { (yyval.bclist) = new BCList; }
    break;

  case 474:

/* Line 1806 of yacc.c  */
#line 1620 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); }
    break;

  case 475:

/* Line 1806 of yacc.c  */
#line 1622 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Forces; (yyvsp[(2) - (2)].bcval).caseid = (yyval.bclist)->caseid; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 476:

/* Line 1806 of yacc.c  */
#line 1624 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Forces, (yyval.bclist)->caseid); (yyval.bclist)->add(bc); } }
    break;

  case 477:

/* Line 1806 of yacc.c  */
#line 1626 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Forces, (yyval.bclist)->caseid); (yyval.bclist)->add(bc); } }
    break;

  case 478:

/* Line 1806 of yacc.c  */
#line 1628 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].caseid = (yyval.bclist)->caseid;
          geoSource->addSurfaceNeuman(1,surf_bc); }
    break;

  case 479:

/* Line 1806 of yacc.c  */
#line 1636 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); }
    break;

  case 480:

/* Line 1806 of yacc.c  */
#line 1638 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 481:

/* Line 1806 of yacc.c  */
#line 1640 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (6)].ival); i<=(yyvsp[(3) - (6)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].fval)); (yyval.bclist)->add(bc); } }
    break;

  case 482:

/* Line 1806 of yacc.c  */
#line 1642 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (7)].bclist); for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval)); (yyval.bclist)->add(bc); } }
    break;

  case 483:

/* Line 1806 of yacc.c  */
#line 1644 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (8)].ival); i<=(yyvsp[(3) - (8)].ival); i+=(yyvsp[(5) - (8)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(6) - (8)].ival)-1, (yyvsp[(7) - (8)].fval)); (yyval.bclist)->add(bc); } }
    break;

  case 484:

/* Line 1806 of yacc.c  */
#line 1646 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (9)].bclist); for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval)); (yyval.bclist)->add(bc); } }
    break;

  case 485:

/* Line 1806 of yacc.c  */
#line 1650 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); }
    break;

  case 486:

/* Line 1806 of yacc.c  */
#line 1652 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 487:

/* Line 1806 of yacc.c  */
#line 1656 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); }
    break;

  case 488:

/* Line 1806 of yacc.c  */
#line 1658 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); }
    break;

  case 491:

/* Line 1806 of yacc.c  */
#line 1666 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ymtt)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addYMTT((yyval.ymtt));}
    break;

  case 492:

/* Line 1806 of yacc.c  */
#line 1668 "p.y"
    { (yyval.ymtt)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); }
    break;

  case 493:

/* Line 1806 of yacc.c  */
#line 1670 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ymtt)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addYMTT((yyval.ymtt));}
    break;

  case 496:

/* Line 1806 of yacc.c  */
#line 1678 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ctett)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addCTETT((yyval.ctett));}
    break;

  case 497:

/* Line 1806 of yacc.c  */
#line 1680 "p.y"
    { (yyval.ctett)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); }
    break;

  case 498:

/* Line 1806 of yacc.c  */
#line 1682 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ctett)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addCTETT((yyval.ctett));}
    break;

  case 501:

/* Line 1806 of yacc.c  */
#line 1690 "p.y"
    { (yyval.lmpcons) = (yyvsp[(1) - (2)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
    break;

  case 502:

/* Line 1806 of yacc.c  */
#line 1694 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); }
    break;

  case 503:

/* Line 1806 of yacc.c  */
#line 1696 "p.y"
    { (yyval.lmpcons) = (yyvsp[(2) - (3)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(3) - (3)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); }
    break;

  case 504:

/* Line 1806 of yacc.c  */
#line 1702 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
    break;

  case 505:

/* Line 1806 of yacc.c  */
#line 1705 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (3)].ival), (yyvsp[(2) - (3)].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
    break;

  case 506:

/* Line 1806 of yacc.c  */
#line 1708 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval));
          (yyval.lmpcons)->type = (yyvsp[(4) - (5)].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
    break;

  case 507:

/* Line 1806 of yacc.c  */
#line 1712 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (4)].ival), (yyvsp[(2) - (4)].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[(3) - (4)].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); }
    break;

  case 508:

/* Line 1806 of yacc.c  */
#line 1719 "p.y"
    { if((yyvsp[(3) - (4)].fval) == 0.0) {
            fprintf(stderr," *** WARNING: zero coefficient in LMPC\n");
            fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival));
          }
          (yyval.mpcterm) = new LMPCTerm();
          (yyval.mpcterm)->nnum = (yyvsp[(1) - (4)].ival)-1;
          (yyval.mpcterm)->dofnum = (yyvsp[(2) - (4)].ival)-1;
          (yyval.mpcterm)->coef.r_value = (yyvsp[(3) - (4)].fval);
        }
    break;

  case 511:

/* Line 1806 of yacc.c  */
#line 1735 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].cxbcval).nnum,(yyvsp[(1) - (2)].cxbcval).reval,(yyvsp[(1) - (2)].cxbcval).imval,(yyvsp[(2) - (2)].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
    break;

  case 512:

/* Line 1806 of yacc.c  */
#line 1737 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); }
    break;

  case 513:

/* Line 1806 of yacc.c  */
#line 1739 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(2) - (3)].cxbcval).nnum,(yyvsp[(2) - (3)].cxbcval).reval,(yyvsp[(2) - (3)].cxbcval).imval,(yyvsp[(3) - (3)].mpcterm)); domain->addLMPC((yyval.lmpcons)); }
    break;

  case 514:

/* Line 1806 of yacc.c  */
#line 1743 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (5)].ival); (yyval.cxbcval).reval=(yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval=(yyvsp[(4) - (5)].fval); }
    break;

  case 515:

/* Line 1806 of yacc.c  */
#line 1745 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (3)].ival); (yyval.cxbcval).reval=(yyvsp[(2) - (3)].fval); (yyval.cxbcval).imval=0.0; }
    break;

  case 516:

/* Line 1806 of yacc.c  */
#line 1747 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (2)].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; }
    break;

  case 517:

/* Line 1806 of yacc.c  */
#line 1751 "p.y"
    { if(((yyvsp[(3) - (5)].fval)==0.0) && ((yyvsp[(4) - (5)].fval)==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (5)].ival),(yyvsp[(2) - (5)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (5)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (5)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval)); }
        }
    break;

  case 518:

/* Line 1806 of yacc.c  */
#line 1759 "p.y"
    { if((yyvsp[(3) - (4)].fval)==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (4)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (4)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (4)].fval),0.0); }
        }
    break;

  case 519:

/* Line 1806 of yacc.c  */
#line 1769 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); }
    break;

  case 520:

/* Line 1806 of yacc.c  */
#line 1771 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].cxbclist)->n; ++i) (yyvsp[(4) - (4)].cxbclist)->d[i].caseid = (yyvsp[(2) - (4)].ival);
          (yyval.cxbclist) = (yyvsp[(4) - (4)].cxbclist); }
    break;

  case 521:

/* Line 1806 of yacc.c  */
#line 1776 "p.y"
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[(1) - (1)].cxbcval)); }
    break;

  case 522:

/* Line 1806 of yacc.c  */
#line 1778 "p.y"
    { (yyval.cxbclist) = (yyvsp[(1) - (2)].cxbclist); (yyval.cxbclist)->add((yyvsp[(2) - (2)].cxbcval)); }
    break;

  case 525:

/* Line 1806 of yacc.c  */
#line 1786 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval); 
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.Iyy = (yyvsp[(14) - (16)].fval); sp.Izz = (yyvsp[(15) - (16)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        }
    break;

  case 526:

/* Line 1806 of yacc.c  */
#line 1794 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.betaDamp = (yyvsp[(17) - (19)].fval); sp.alphaDamp = (yyvsp[(18) - (19)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        }
    break;

  case 527:

/* Line 1806 of yacc.c  */
#line 1803 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.lagrangeMult = bool((yyvsp[(17) - (19)].ival));
          sp.penalty = (yyvsp[(18) - (19)].fval);
          sp.type = StructProp::Constraint;
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        }
    break;

  case 528:

/* Line 1806 of yacc.c  */
#line 1814 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (20)].fval);  sp.E = (yyvsp[(3) - (20)].fval);  sp.nu  = (yyvsp[(4) - (20)].fval);  sp.rho = (yyvsp[(5) - (20)].fval);
          sp.c = (yyvsp[(6) - (20)].fval);  sp.k = (yyvsp[(7) - (20)].fval);  sp.eh  = (yyvsp[(8) - (20)].fval);  sp.P   = (yyvsp[(9) - (20)].fval);  sp.Ta  = (yyvsp[(10) - (20)].fval); 
          sp.Q = (yyvsp[(11) - (20)].fval); sp.W = (yyvsp[(12) - (20)].fval); sp.Ixx = (yyvsp[(13) - (20)].fval); sp.Iyy = (yyvsp[(14) - (20)].fval); sp.Izz = (yyvsp[(15) - (20)].fval);
	  sp.ymin = (yyvsp[(16) - (20)].fval); sp.ymax = (yyvsp[(17) - (20)].fval); sp.zmin = (yyvsp[(18) - (20)].fval); sp.zmax = (yyvsp[(19) - (20)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (20)].ival)-1, sp );
        }
    break;

  case 529:

/* Line 1806 of yacc.c  */
#line 1823 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.betaDamp = (yyvsp[(21) - (23)].fval); sp.alphaDamp = (yyvsp[(22) - (23)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        }
    break;

  case 530:

/* Line 1806 of yacc.c  */
#line 1833 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.lagrangeMult = bool((yyvsp[(21) - (23)].ival));
          sp.penalty = (yyvsp[(22) - (23)].fval);
          sp.type = StructProp::Constraint;
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        }
    break;

  case 531:

/* Line 1806 of yacc.c  */
#line 1845 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (9)].fval); sp.E = (yyvsp[(3) - (9)].fval); sp.nu = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(5) - (9)].fval);
          sp.c = (yyvsp[(6) - (9)].fval); sp.k = (yyvsp[(7) - (9)].fval); sp.eh = (yyvsp[(8) - (9)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp ); 
        }
    break;

  case 532:

/* Line 1806 of yacc.c  */
#line 1852 "p.y"
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (14)].fval);  sp.E = (yyvsp[(3) - (14)].fval);  sp.nu  = (yyvsp[(4) - (14)].fval);  sp.rho = (yyvsp[(5) - (14)].fval);
          sp.c = (yyvsp[(6) - (14)].fval);  sp.k = (yyvsp[(7) - (14)].fval);  sp.eh  = (yyvsp[(8) - (14)].fval);  sp.P   = (yyvsp[(9) - (14)].fval);  sp.Ta  = (yyvsp[(10) - (14)].fval);
          sp.Q = (yyvsp[(11) - (14)].fval); sp.W = (yyvsp[(12) - (14)].fval); sp.Ixx = (yyvsp[(13) - (14)].fval);  
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
        }
    break;

  case 533:

/* Line 1806 of yacc.c  */
#line 1860 "p.y"
    { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval);
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.betaDamp = (yyvsp[(15) - (16)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        }
    break;

  case 534:

/* Line 1806 of yacc.c  */
#line 1868 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (12)].fval);
          sp.fp.PMLtype = int((yyvsp[(4) - (12)].fval));
          sp.fp.gamma = (yyvsp[(5) - (12)].fval);
          sp.fp.Rx = (yyvsp[(6) - (12)].fval);
          sp.fp.Sx = (yyvsp[(7) - (12)].fval);
          sp.fp.Ry = (yyvsp[(8) - (12)].fval);
          sp.fp.Sy = (yyvsp[(9) - (12)].fval);
          sp.fp.Rz = (yyvsp[(10) - (12)].fval);
          sp.fp.Sz = (yyvsp[(11) - (12)].fval);
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
          domain->PMLFlag = 1;
        }
    break;

  case 535:

/* Line 1806 of yacc.c  */
#line 1884 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (13)].fval);
          sp.rho = (yyvsp[(4) - (13)].fval);
          sp.fp.PMLtype = int((yyvsp[(5) - (13)].fval));
          sp.fp.gamma = (yyvsp[(6) - (13)].fval);
          sp.fp.Rx = (yyvsp[(7) - (13)].fval);
          sp.fp.Sx = (yyvsp[(8) - (13)].fval);
          sp.fp.Ry = (yyvsp[(9) - (13)].fval);
          sp.fp.Sy = (yyvsp[(10) - (13)].fval);
          sp.fp.Rz = (yyvsp[(11) - (13)].fval);
          sp.fp.Sz = (yyvsp[(12) - (13)].fval);
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (13)].ival)-1, sp );
          domain->PMLFlag = 1;
        }
    break;

  case 536:

/* Line 1806 of yacc.c  */
#line 1901 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (14)].fval);
          sp.kappaHelmImag = (yyvsp[(4) - (14)].fval);
          sp.rho = (yyvsp[(5) - (14)].fval);
          sp.fp.PMLtype = int((yyvsp[(6) - (14)].fval));
          sp.fp.gamma = (yyvsp[(7) - (14)].fval);
          sp.fp.Rx = (yyvsp[(8) - (14)].fval);
          sp.fp.Sx = (yyvsp[(9) - (14)].fval);
          sp.fp.Ry = (yyvsp[(10) - (14)].fval);
          sp.fp.Sy = (yyvsp[(11) - (14)].fval);
          sp.fp.Rz = (yyvsp[(12) - (14)].fval);
          sp.fp.Sz = (yyvsp[(13) - (14)].fval);
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
          domain->PMLFlag = 1;
        }
    break;

  case 537:

/* Line 1806 of yacc.c  */
#line 1919 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (5)].fval);
          sp.rho = (yyvsp[(4) - (5)].fval);
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        }
    break;

  case 538:

/* Line 1806 of yacc.c  */
#line 1927 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (6)].fval);
          sp.kappaHelmImag = (yyvsp[(4) - (6)].fval);
          sp.rho = (yyvsp[(5) - (6)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        }
    break;

  case 539:

/* Line 1806 of yacc.c  */
#line 1935 "p.y"
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
	  sp.isReal = true;
          sp.type = StructProp::Fabric;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        }
    break;

  case 540:

/* Line 1806 of yacc.c  */
#line 1954 "p.y"
    { StructProp sp; 
          sp.A = (yyvsp[(3) - (12)].fval);  sp.rho = (yyvsp[(4) - (12)].fval); sp.Q = (yyvsp[(5) - (12)].fval); sp.c = (yyvsp[(6) - (12)].fval); 
          sp.sigma = (yyvsp[(7) - (12)].fval);  sp.k = (yyvsp[(8) - (12)].fval);  sp.eh  = (yyvsp[(9) - (12)].fval);  sp.P   = (yyvsp[(10) - (12)].fval);  sp.Ta  = (yyvsp[(11) - (12)].fval);
          sp.isReal = true;
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        }
    break;

  case 541:

/* Line 1806 of yacc.c  */
#line 1962 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (5)].ival));
          sp.penalty = (yyvsp[(4) - (5)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        }
    break;

  case 542:

/* Line 1806 of yacc.c  */
#line 1969 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (7)].ival));
          sp.penalty = (yyvsp[(4) - (7)].fval);
          sp.amplitude = (yyvsp[(5) - (7)].fval);
          sp.omega = (yyvsp[(6) - (7)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        }
    break;

  case 543:

/* Line 1806 of yacc.c  */
#line 1978 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (8)].ival));
          sp.penalty = (yyvsp[(4) - (8)].fval);
          sp.amplitude = (yyvsp[(5) - (8)].fval);
          sp.omega = (yyvsp[(6) - (8)].fval);
          sp.phase = (yyvsp[(7) - (8)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        }
    break;

  case 544:

/* Line 1806 of yacc.c  */
#line 1988 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (9)].ival));
          sp.penalty = (yyvsp[(4) - (9)].fval);
          sp.amplitude = (yyvsp[(5) - (9)].fval);
          sp.omega = (yyvsp[(6) - (9)].fval);
          sp.phase = (yyvsp[(7) - (9)].fval);
          sp.offset = (yyvsp[(8) - (9)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        }
    break;

  case 545:

/* Line 1806 of yacc.c  */
#line 1999 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (11)].ival));
          sp.penalty = (yyvsp[(4) - (11)].fval);
          sp.amplitude = (yyvsp[(5) - (11)].fval);
          sp.omega = (yyvsp[(6) - (11)].fval);
          sp.phase = (yyvsp[(7) - (11)].fval);
          sp.B = (yyvsp[(8) - (11)].fval);
          sp.C = (yyvsp[(9) - (11)].fval);
          sp.relop = (yyvsp[(10) - (11)].ival);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        }
    break;

  case 546:

/* Line 1806 of yacc.c  */
#line 2012 "p.y"
    { // new style for joints with prescribed motion by 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (10)].ival));
          sp.penalty = (yyvsp[(4) - (10)].fval);
          sp.funtype = (yyvsp[(5) - (10)].ival);
          sp.amplitude = (yyvsp[(6) - (10)].fval);
          sp.offset = (yyvsp[(7) - (10)].fval);
          sp.c1 = (yyvsp[(8) - (10)].fval);
          sp.c2 = (yyvsp[(9) - (10)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        }
    break;

  case 547:

/* Line 1806 of yacc.c  */
#line 2025 "p.y"
    { // new style for joints with prescribed motion by 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (11)].ival));
          sp.penalty = (yyvsp[(4) - (11)].fval);
          sp.funtype = (yyvsp[(5) - (11)].ival);
          sp.amplitude = (yyvsp[(6) - (11)].fval);
          sp.offset = (yyvsp[(7) - (11)].fval);
          sp.c1 = (yyvsp[(8) - (11)].fval);
          sp.c2 = (yyvsp[(9) - (11)].fval);
          sp.c3 = (yyvsp[(10) - (11)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        }
    break;

  case 548:

/* Line 1806 of yacc.c  */
#line 2039 "p.y"
    { // new style for joints with prescribed motion by 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (12)].ival));
          sp.penalty = (yyvsp[(4) - (12)].fval);
          sp.funtype = (yyvsp[(5) - (12)].ival);
          sp.amplitude = (yyvsp[(6) - (12)].fval);
          sp.offset = (yyvsp[(7) - (12)].fval);
          sp.c1 = (yyvsp[(8) - (12)].fval);
          sp.c2 = (yyvsp[(9) - (12)].fval);
          sp.c3 = (yyvsp[(10) - (12)].fval);
          sp.c4 = (yyvsp[(11) - (12)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        }
    break;

  case 549:

/* Line 1806 of yacc.c  */
#line 2054 "p.y"
    { // use for RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (7)].ival));
          sp.penalty = (yyvsp[(4) - (7)].fval);
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(6) - (7)].fval);
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        }
    break;

  case 550:

/* Line 1806 of yacc.c  */
#line 2063 "p.y"
    { // use for UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (8)].ival));
          sp.penalty = (yyvsp[(4) - (8)].fval);
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(6) - (8)].fval);
          sp.k2 = (yyvsp[(7) - (8)].fval);
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        }
    break;

  case 551:

/* Line 1806 of yacc.c  */
#line 2073 "p.y"
    { // use for SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (9)].ival));
          sp.penalty = (yyvsp[(4) - (9)].fval);
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(6) - (9)].fval);
          sp.k2 = (yyvsp[(7) - (9)].fval);
          sp.k3 = (yyvsp[(8) - (9)].fval);
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        }
    break;

  case 552:

/* Line 1806 of yacc.c  */
#line 2084 "p.y"
    { // use for TorsionalSpringType1 or TranslationalSpring
          StructProp sp;
          sp.k1 = (yyvsp[(3) - (4)].fval);
          geoSource->addMat( (yyvsp[(1) - (4)].ival)-1, sp );
        }
    break;

  case 555:

/* Line 1806 of yacc.c  */
#line 2096 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (3)].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
    break;

  case 556:

/* Line 1806 of yacc.c  */
#line 2102 "p.y"
    { if((yyvsp[(2) - (4)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (4)].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
    break;

  case 557:

/* Line 1806 of yacc.c  */
#line 2108 "p.y"
    { if((yyvsp[(2) - (5)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (5)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (5)].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
    break;

  case 558:

/* Line 1806 of yacc.c  */
#line 2115 "p.y"
    { if((yyvsp[(2) - (6)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (6)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (6)].fval));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        }
    break;

  case 559:

/* Line 1806 of yacc.c  */
#line 2123 "p.y"
    { if((yyval.SurfObj)->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[(yyvsp[(4) - (5)].nl).num];
            for(int i=0; i<(yyvsp[(4) - (5)].nl).num; ++i) nodes[(yyvsp[(4) - (5)].nl).num-1-i] = (yyvsp[(4) - (5)].nl).nd[i];
            (yyval.SurfObj)->AddFaceElement((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].nl).num, nodes);
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
          else (yyval.SurfObj)->AddFaceElement((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].nl).num, (yyvsp[(4) - (5)].nl).nd);
        }
    break;

  case 560:

/* Line 1806 of yacc.c  */
#line 2147 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
    break;

  case 561:

/* Line 1806 of yacc.c  */
#line 2149 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); }
    break;

  case 562:

/* Line 1806 of yacc.c  */
#line 2151 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
    break;

  case 563:

/* Line 1806 of yacc.c  */
#line 2155 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (6)].ival), (yyvsp[(3) - (6)].ival), (yyvsp[(5) - (6)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
    break;

  case 564:

/* Line 1806 of yacc.c  */
#line 2157 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
    break;

  case 565:

/* Line 1806 of yacc.c  */
#line 2159 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
    break;

  case 566:

/* Line 1806 of yacc.c  */
#line 2161 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); }
    break;

  case 567:

/* Line 1806 of yacc.c  */
#line 2163 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
    break;

  case 568:

/* Line 1806 of yacc.c  */
#line 2167 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
    break;

  case 569:

/* Line 1806 of yacc.c  */
#line 2173 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->solInfo().isCoupled = true; }
    break;

  case 570:

/* Line 1806 of yacc.c  */
#line 2175 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (3)].ival), (yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; }
    break;

  case 571:

/* Line 1806 of yacc.c  */
#line 2183 "p.y"
    { }
    break;

  case 572:

/* Line 1806 of yacc.c  */
#line 2185 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
    break;

  case 573:

/* Line 1806 of yacc.c  */
#line 2191 "p.y"
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        }
    break;

  case 574:

/* Line 1806 of yacc.c  */
#line 2198 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 575:

/* Line 1806 of yacc.c  */
#line 2205 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 576:

/* Line 1806 of yacc.c  */
#line 2212 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 577:

/* Line 1806 of yacc.c  */
#line 2220 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 578:

/* Line 1806 of yacc.c  */
#line 2229 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 579:

/* Line 1806 of yacc.c  */
#line 2236 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 580:

/* Line 1806 of yacc.c  */
#line 2246 "p.y"
    { }
    break;

  case 581:

/* Line 1806 of yacc.c  */
#line 2248 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        }
    break;

  case 582:

/* Line 1806 of yacc.c  */
#line 2254 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        }
    break;

  case 583:

/* Line 1806 of yacc.c  */
#line 2262 "p.y"
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
    break;

  case 585:

/* Line 1806 of yacc.c  */
#line 2268 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; }
    break;

  case 586:

/* Line 1806 of yacc.c  */
#line 2274 "p.y"
    { }
    break;

  case 587:

/* Line 1806 of yacc.c  */
#line 2276 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        }
    break;

  case 588:

/* Line 1806 of yacc.c  */
#line 2282 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        }
    break;

  case 589:

/* Line 1806 of yacc.c  */
#line 2292 "p.y"
    { }
    break;

  case 590:

/* Line 1806 of yacc.c  */
#line 2294 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 591:

/* Line 1806 of yacc.c  */
#line 2301 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 592:

/* Line 1806 of yacc.c  */
#line 2308 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 593:

/* Line 1806 of yacc.c  */
#line 2315 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 594:

/* Line 1806 of yacc.c  */
#line 2322 "p.y"
    { /* this one is for frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 595:

/* Line 1806 of yacc.c  */
#line 2330 "p.y"
    { /* this one is for constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 596:

/* Line 1806 of yacc.c  */
#line 2339 "p.y"
    { /* this one is for velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (13)].ival), (yyvsp[(4) - (13)].ival), (yyvsp[(6) - (13)].fval), (yyvsp[(7) - (13)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (13)].ival), (yyvsp[(9) - (13)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (13)].fval), (yyvsp[(11) - (13)].fval), (yyvsp[(12) - (13)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (13)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 597:

/* Line 1806 of yacc.c  */
#line 2348 "p.y"
    { /* this one is for pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (14)].ival), (yyvsp[(4) - (14)].ival), (yyvsp[(6) - (14)].fval), (yyvsp[(7) - (14)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (14)].ival), (yyvsp[(9) - (14)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (14)].fval), (yyvsp[(11) - (14)].fval), (yyvsp[(12) - (14)].fval), (yyvsp[(13) - (14)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (14)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 598:

/* Line 1806 of yacc.c  */
#line 2357 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 599:

/* Line 1806 of yacc.c  */
#line 2365 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival)); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        }
    break;

  case 600:

/* Line 1806 of yacc.c  */
#line 2375 "p.y"
    { domain->solInfo().dist_acme = (yyvsp[(2) - (3)].ival); }
    break;

  case 601:

/* Line 1806 of yacc.c  */
#line 2377 "p.y"
    { domain->solInfo().ffi_debug = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 602:

/* Line 1806 of yacc.c  */
#line 2379 "p.y"
    { domain->solInfo().mortar_scaling = (yyvsp[(2) - (3)].fval); }
    break;

  case 603:

/* Line 1806 of yacc.c  */
#line 2381 "p.y"
    { domain->solInfo().mortar_integration_rule = (yyvsp[(2) - (3)].ival); }
    break;

  case 604:

/* Line 1806 of yacc.c  */
#line 2384 "p.y"
    { geoSource->addNode((yyvsp[(3) - (3)].nval).num, (yyvsp[(3) - (3)].nval).xyz); }
    break;

  case 605:

/* Line 1806 of yacc.c  */
#line 2386 "p.y"
    { geoSource->addNode((yyvsp[(2) - (2)].nval).num, (yyvsp[(2) - (2)].nval).xyz); }
    break;

  case 606:

/* Line 1806 of yacc.c  */
#line 2390 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (5)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (5)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (5)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (5)].fval); }
    break;

  case 607:

/* Line 1806 of yacc.c  */
#line 2392 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (4)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (4)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (4)].fval);  (yyval.nval).xyz[2] = 0.0; }
    break;

  case 608:

/* Line 1806 of yacc.c  */
#line 2394 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (3)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (3)].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; }
    break;

  case 609:

/* Line 1806 of yacc.c  */
#line 2398 "p.y"
    { /* Define each Element */
          geoSource->addElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); }
    break;

  case 610:

/* Line 1806 of yacc.c  */
#line 2403 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1; }
    break;

  case 611:

/* Line 1806 of yacc.c  */
#line 2405 "p.y"
    { if((yyval.nl).num == 125) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++; }
    break;

  case 612:

/* Line 1806 of yacc.c  */
#line 2410 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.bcval).val = (yyvsp[(3) - (4)].fval); }
    break;

  case 613:

/* Line 1806 of yacc.c  */
#line 2412 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (3)].ival)-1; (yyval.bcval).val = 0.0; }
    break;

  case 614:

/* Line 1806 of yacc.c  */
#line 2416 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[(2) - (3)].fval); }
    break;

  case 615:

/* Line 1806 of yacc.c  */
#line 2420 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); }
    break;

  case 616:

/* Line 1806 of yacc.c  */
#line 2424 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (5)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (5)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval = (yyvsp[(4) - (5)].fval);  }
    break;

  case 617:

/* Line 1806 of yacc.c  */
#line 2426 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (4)].fval); (yyval.cxbcval).imval = 0.0; }
    break;

  case 619:

/* Line 1806 of yacc.c  */
#line 2431 "p.y"
    { geoSource->setCSFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); }
    break;

  case 621:

/* Line 1806 of yacc.c  */
#line 2440 "p.y"
    { geoSource->setFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); }
    break;

  case 622:

/* Line 1806 of yacc.c  */
#line 2444 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (11)].ival)-1; 
            (yyval.frame).d[0] = (yyvsp[(2) - (11)].fval); (yyval.frame).d[1] = (yyvsp[(3) - (11)].fval); (yyval.frame).d[2] = (yyvsp[(4) - (11)].fval);
            (yyval.frame).d[3] = (yyvsp[(5) - (11)].fval); (yyval.frame).d[4] = (yyvsp[(6) - (11)].fval); (yyval.frame).d[5] = (yyvsp[(7) - (11)].fval);
            (yyval.frame).d[6] = (yyvsp[(8) - (11)].fval); (yyval.frame).d[7] = (yyvsp[(9) - (11)].fval); (yyval.frame).d[8] = (yyvsp[(10) - (11)].fval); }
    break;

  case 623:

/* Line 1806 of yacc.c  */
#line 2449 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (4)].ival)-1; geoSource->makeEframe((yyvsp[(1) - (4)].ival)-1, (yyvsp[(3) - (4)].ival), (yyval.frame).d); }
    break;

  case 625:

/* Line 1806 of yacc.c  */
#line 2454 "p.y"
    { OffsetData od;
	  od.first = (yyvsp[(2) - (7)].ival)-1; od.last = (yyvsp[(3) - (7)].ival)-1;
	  od.o[0] = (yyvsp[(4) - (7)].fval); od.o[1] = (yyvsp[(5) - (7)].fval); od.o[2] = (yyvsp[(6) - (7)].fval); 
	  geoSource->addOffset(od); }
    break;

  case 626:

/* Line 1806 of yacc.c  */
#line 2461 "p.y"
    { (yyval.ival) = 0; }
    break;

  case 627:

/* Line 1806 of yacc.c  */
#line 2463 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].ival)-1); }
    break;

  case 628:

/* Line 1806 of yacc.c  */
#line 2465 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1); 
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (6)].ival) - 1, (yyvsp[(5) - (6)].fval));
	  domain->solInfo().elemLumpPodRom = true; }
    break;

  case 629:

/* Line 1806 of yacc.c  */
#line 2469 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(4) - (6)].ival)-1,(yyvsp[(5) - (6)].ival)-1); }
    break;

  case 630:

/* Line 1806 of yacc.c  */
#line 2471 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (8)].ival)-1,(yyvsp[(3) - (8)].ival)-1,(yyvsp[(4) - (8)].ival)-1,(yyvsp[(5) - (8)].ival)-1);
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (8)].ival) - 1, (yyvsp[(7) - (8)].fval)); 
	  domain->solInfo().elemLumpPodRom = true; }
    break;

  case 631:

/* Line 1806 of yacc.c  */
#line 2475 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1,(yyvsp[(4) - (7)].ival)-1,-1,(yyvsp[(6) - (7)].fval)); }
    break;

  case 632:

/* Line 1806 of yacc.c  */
#line 2477 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        }
    break;

  case 633:

/* Line 1806 of yacc.c  */
#line 2482 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
 	    geoSource->setAttrib(i-1,(yyvsp[(4) - (5)].ival)-1);
	}
    break;

  case 634:

/* Line 1806 of yacc.c  */
#line 2487 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (7)].ival); i<(yyvsp[(3) - (7)].ival)+1; ++i)
	    geoSource->setAttrib(i-1, (yyvsp[(4) - (7)].ival)-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].ival)-1);
	}
    break;

  case 635:

/* Line 1806 of yacc.c  */
#line 2492 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (8)].ival); i<(yyvsp[(3) - (8)].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[(4) - (8)].ival)-1, (yyvsp[(5) - (8)].ival)-1, -1, (yyvsp[(7) - (8)].fval));
        }
    break;

  case 637:

/* Line 1806 of yacc.c  */
#line 2500 "p.y"
    { geoSource->setElementPressure( (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval) ); }
    break;

  case 638:

/* Line 1806 of yacc.c  */
#line 2502 "p.y"
    { int i; 
          for(i=(yyvsp[(2) - (5)].ival); i<((yyvsp[(3) - (5)].ival)+1); ++i)
            geoSource->setElementPressure( i-1, (yyvsp[(4) - (5)].fval) );  
        }
    break;

  case 639:

/* Line 1806 of yacc.c  */
#line 2507 "p.y"
    { BCond *surf_pbc = new BCond[1];
          surf_pbc[0] = (yyvsp[(3) - (3)].bcval);
          geoSource->addSurfacePressure(1,surf_pbc); }
    break;

  case 640:

/* Line 1806 of yacc.c  */
#line 2513 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        }
    break;

  case 641:

/* Line 1806 of yacc.c  */
#line 2518 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[(2) - (3)].ival));
          geoSource->setConsistentPFlag(false);
        }
    break;

  case 642:

/* Line 1806 of yacc.c  */
#line 2525 "p.y"
    { }
    break;

  case 643:

/* Line 1806 of yacc.c  */
#line 2527 "p.y"
    { geoSource->setElementPreLoad( (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval) ); }
    break;

  case 644:

/* Line 1806 of yacc.c  */
#line 2529 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (6)].ival); i<((yyvsp[(4) - (6)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[(5) - (6)].fval) );
        }
    break;

  case 645:

/* Line 1806 of yacc.c  */
#line 2534 "p.y"
    { double load[3] = { (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval) };
          geoSource->setElementPreLoad( (yyvsp[(2) - (6)].ival)-1, load ); }
    break;

  case 646:

/* Line 1806 of yacc.c  */
#line 2537 "p.y"
    { double load[3] = { (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
          int i;
          for(i=(yyvsp[(2) - (8)].ival); i<((yyvsp[(4) - (8)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        }
    break;

  case 650:

/* Line 1806 of yacc.c  */
#line 2550 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(1) - (1)].ival)); }
    break;

  case 651:

/* Line 1806 of yacc.c  */
#line 2552 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(2) - (2)].ival)); }
    break;

  case 652:

/* Line 1806 of yacc.c  */
#line 2556 "p.y"
    { domain->solInfo().type = 1;
          domain->solInfo().iterType = (yyvsp[(3) - (4)].ival);
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 653:

/* Line 1806 of yacc.c  */
#line 2560 "p.y"
    { domain->solInfo().precond = (yyvsp[(3) - (4)].ival); }
    break;

  case 654:

/* Line 1806 of yacc.c  */
#line 2562 "p.y"
    { domain->solInfo().maxit = (yyvsp[(3) - (4)].ival); }
    break;

  case 655:

/* Line 1806 of yacc.c  */
#line 2564 "p.y"
    { domain->solInfo().tol = (yyvsp[(3) - (4)].fval); }
    break;

  case 656:

/* Line 1806 of yacc.c  */
#line 2566 "p.y"
    { domain->solInfo().maxvecsize = (yyvsp[(3) - (4)].ival); }
    break;

  case 657:

/* Line 1806 of yacc.c  */
#line 2568 "p.y"
    { domain->solInfo().iterSubtype = (yyvsp[(3) - (4)].ival); }
    break;

  case 658:

/* Line 1806 of yacc.c  */
#line 2572 "p.y"
    { domain->solInfo().setSolver(0); 
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 659:

/* Line 1806 of yacc.c  */
#line 2575 "p.y"
    { domain->solInfo().setSolver((yyvsp[(4) - (5)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 660:

/* Line 1806 of yacc.c  */
#line 2578 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 661:

/* Line 1806 of yacc.c  */
#line 2581 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (5)].ival));
          domain->solInfo().setProbType(SolverInfo::Static);
          if((yyvsp[(3) - (5)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
    break;

  case 662:

/* Line 1806 of yacc.c  */
#line 2586 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (5)].ival));
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().getNLInfo().unsymmetric = true; }
    break;

  case 663:

/* Line 1806 of yacc.c  */
#line 2590 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (5)].ival),(yyvsp[(4) - (5)].ival));    
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 664:

/* Line 1806 of yacc.c  */
#line 2593 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].fval));    
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 665:

/* Line 1806 of yacc.c  */
#line 2596 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (7)].ival),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].fval),(yyvsp[(6) - (7)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 666:

/* Line 1806 of yacc.c  */
#line 2599 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (8)].ival),(yyvsp[(4) - (8)].ival),(yyvsp[(5) - (8)].fval),(yyvsp[(6) - (8)].ival),(yyvsp[(7) - (8)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 667:

/* Line 1806 of yacc.c  */
#line 2602 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (9)].ival),(yyvsp[(4) - (9)].ival),(yyvsp[(5) - (9)].fval),(yyvsp[(6) - (9)].ival),(yyvsp[(7) - (9)].ival),(yyvsp[(8) - (9)].ival));
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 668:

/* Line 1806 of yacc.c  */
#line 2610 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = (yyvsp[(4) - (6)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(5) - (6)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(4) - (6)].ival);
          domain->solInfo().type =(2); }
    break;

  case 669:

/* Line 1806 of yacc.c  */
#line 2616 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(4) - (5)].ival)-1); }
    break;

  case 670:

/* Line 1806 of yacc.c  */
#line 2620 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; }
    break;

  case 671:

/* Line 1806 of yacc.c  */
#line 2625 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; }
    break;

  case 672:

/* Line 1806 of yacc.c  */
#line 2631 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(4) - (6)].ival)-1); 
          domain->solInfo().fetiInfo.feti2version 
                  = (FetiInfo::Feti2Version) (yyvsp[(5) - (6)].ival); }
    break;

  case 673:

/* Line 1806 of yacc.c  */
#line 2637 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = (yyvsp[(4) - (7)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(5) - (7)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(6) - (7)].ival);
          domain->solInfo().type =(2); }
    break;

  case 674:

/* Line 1806 of yacc.c  */
#line 2643 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(3) - (5)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(4) - (5)].ival);
          domain->solInfo().type =(2); }
    break;

  case 675:

/* Line 1806 of yacc.c  */
#line 2649 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static); }
    break;

  case 676:

/* Line 1806 of yacc.c  */
#line 2652 "p.y"
    { domain->solInfo().type =(2);}
    break;

  case 677:

/* Line 1806 of yacc.c  */
#line 2654 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (3)].ival)-1);}
    break;

  case 678:

/* Line 1806 of yacc.c  */
#line 2657 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; }
    break;

  case 679:

/* Line 1806 of yacc.c  */
#line 2661 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; }
    break;

  case 680:

/* Line 1806 of yacc.c  */
#line 2666 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (4)].ival)-1);
          domain->solInfo().fetiInfo.feti2version = (FetiInfo::Feti2Version) (yyvsp[(3) - (4)].ival); 
        }
    break;

  case 681:

/* Line 1806 of yacc.c  */
#line 2671 "p.y"
    {
	  domain->solInfo().type = 3;
          domain->solInfo().subtype = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getFetiInfo().solvertype = (FetiInfo::Solvertype)((yyvsp[(3) - (4)].ival));
	}
    break;

  case 682:

/* Line 1806 of yacc.c  */
#line 2677 "p.y"
    { domain->solInfo().sparse_maxsup  = (yyvsp[(2) - (3)].ival); }
    break;

  case 683:

/* Line 1806 of yacc.c  */
#line 2679 "p.y"
    { domain->solInfo().sparse_defblk  = (yyvsp[(2) - (3)].ival); }
    break;

  case 684:

/* Line 1806 of yacc.c  */
#line 2681 "p.y"
    { domain->solInfo().spooles_tau  = (yyvsp[(2) - (3)].fval); }
    break;

  case 685:

/* Line 1806 of yacc.c  */
#line 2683 "p.y"
    { domain->solInfo().spooles_maxsize = (yyvsp[(2) - (3)].ival); }
    break;

  case 686:

/* Line 1806 of yacc.c  */
#line 2685 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 0) {
            (yyvsp[(2) - (3)].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          domain->solInfo().spooles_maxdomainsize = (yyvsp[(2) - (3)].ival); }
    break;

  case 687:

/* Line 1806 of yacc.c  */
#line 2692 "p.y"
    { domain->solInfo().spooles_seed = (yyvsp[(2) - (3)].ival); }
    break;

  case 688:

/* Line 1806 of yacc.c  */
#line 2694 "p.y"
    { if(((yyvsp[(2) - (3)].fval) < 0.0) || ((yyvsp[(2) - (3)].fval) > 1.0)) {
            (yyvsp[(2) - (3)].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          domain->solInfo().spooles_maxzeros = (yyvsp[(2) - (3)].fval); }
    break;

  case 689:

/* Line 1806 of yacc.c  */
#line 2701 "p.y"
    { domain->solInfo().spooles_msglvl = (yyvsp[(2) - (3)].ival); }
    break;

  case 690:

/* Line 1806 of yacc.c  */
#line 2703 "p.y"
    { domain->solInfo().pivot = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 691:

/* Line 1806 of yacc.c  */
#line 2705 "p.y"
    { domain->solInfo().spooles_scale = (yyvsp[(2) - (3)].ival); }
    break;

  case 692:

/* Line 1806 of yacc.c  */
#line 2707 "p.y"
    { domain->solInfo().spooles_renum = (yyvsp[(2) - (3)].ival); }
    break;

  case 693:

/* Line 1806 of yacc.c  */
#line 2709 "p.y"
    { domain->solInfo().mumps_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); }
    break;

  case 694:

/* Line 1806 of yacc.c  */
#line 2711 "p.y"
    { domain->solInfo().mumps_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); }
    break;

  case 695:

/* Line 1806 of yacc.c  */
#line 2713 "p.y"
    { domain->solInfo().goldfarb_tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 696:

/* Line 1806 of yacc.c  */
#line 2715 "p.y"
    { domain->solInfo().goldfarb_check = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 697:

/* Line 1806 of yacc.c  */
#line 2717 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(3) - (4)].ival); }
    break;

  case 698:

/* Line 1806 of yacc.c  */
#line 2719 "p.y"
    { domain->solInfo().debug_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); }
    break;

  case 699:

/* Line 1806 of yacc.c  */
#line 2721 "p.y"
    { domain->solInfo().debug_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); }
    break;

  case 700:

/* Line 1806 of yacc.c  */
#line 2727 "p.y"
    { domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival); }
    break;

  case 701:

/* Line 1806 of yacc.c  */
#line 2729 "p.y"
    { domain->solInfo().fetiInfo.precno = FetiInfo::lumped; }
    break;

  case 702:

/* Line 1806 of yacc.c  */
#line 2731 "p.y"
    { if(((yyvsp[(3) - (4)].ival) < 0) || ((yyvsp[(3) - (4)].ival) > 3)) { 
            (yyvsp[(3) - (4)].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival);
	}
    break;

  case 703:

/* Line 1806 of yacc.c  */
#line 2738 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        }
    break;

  case 704:

/* Line 1806 of yacc.c  */
#line 2745 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        }
    break;

  case 705:

/* Line 1806 of yacc.c  */
#line 2752 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 706:

/* Line 1806 of yacc.c  */
#line 2754 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.absolute_tol = (yyvsp[(3) - (4)].fval); }
    break;

  case 707:

/* Line 1806 of yacc.c  */
#line 2757 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 708:

/* Line 1806 of yacc.c  */
#line 2759 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.absolute_stagnation_tol = (yyvsp[(3) - (4)].fval); }
    break;

  case 709:

/* Line 1806 of yacc.c  */
#line 2762 "p.y"
    { domain->solInfo().fetiInfo.primal_proj_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_proj_tol = (yyvsp[(3) - (4)].fval); }
    break;

  case 710:

/* Line 1806 of yacc.c  */
#line 2765 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.dual_plan_maxit = (yyvsp[(3) - (4)].ival); }
    break;

  case 711:

/* Line 1806 of yacc.c  */
#line 2768 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_plan_tol = (yyvsp[(3) - (4)].fval); }
    break;

  case 712:

/* Line 1806 of yacc.c  */
#line 2771 "p.y"
    { domain->solInfo().fetiInfo.maxortho = (yyvsp[(3) - (4)].ival); }
    break;

  case 713:

/* Line 1806 of yacc.c  */
#line 2773 "p.y"
    { domain->solInfo().fetiInfo.noCoarse = 1; }
    break;

  case 714:

/* Line 1806 of yacc.c  */
#line 2775 "p.y"
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
        }
    break;

  case 715:

/* Line 1806 of yacc.c  */
#line 2795 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 1; 
          domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); }
    break;

  case 716:

/* Line 1806 of yacc.c  */
#line 2798 "p.y"
    { domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); }
    break;

  case 717:

/* Line 1806 of yacc.c  */
#line 2800 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); }
    break;

  case 718:

/* Line 1806 of yacc.c  */
#line 2803 "p.y"
    { domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); }
    break;

  case 719:

/* Line 1806 of yacc.c  */
#line 2805 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); }
    break;

  case 720:

/* Line 1806 of yacc.c  */
#line 2808 "p.y"
    { domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); }
    break;

  case 721:

/* Line 1806 of yacc.c  */
#line 2810 "p.y"
    { domain->solInfo().fetiInfo.mpc_element = true; }
    break;

  case 722:

/* Line 1806 of yacc.c  */
#line 2812 "p.y"
    { domain->solInfo().fetiInfo.fsi_element = true; }
    break;

  case 723:

/* Line 1806 of yacc.c  */
#line 2814 "p.y"
    { domain->solInfo().fetiInfo.fsi_corner = (yyvsp[(2) - (3)].ival); }
    break;

  case 724:

/* Line 1806 of yacc.c  */
#line 2816 "p.y"
    { domain->solInfo().fetiInfo.splitLocalFsi = false; }
    break;

  case 725:

/* Line 1806 of yacc.c  */
#line 2818 "p.y"
    { domain->solInfo().coupled_scale = (yyvsp[(2) - (3)].fval); }
    break;

  case 726:

/* Line 1806 of yacc.c  */
#line 2820 "p.y"
    { domain->solInfo().fetiInfo.wetcorners = true; }
    break;

  case 727:

/* Line 1806 of yacc.c  */
#line 2822 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (3)].ival); }
    break;

  case 728:

/* Line 1806 of yacc.c  */
#line 2824 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.pick_unsafe_corners = bool((yyvsp[(3) - (4)].ival));
        }
    break;

  case 729:

/* Line 1806 of yacc.c  */
#line 2828 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) {
            domain->solInfo().fetiInfo.corners = FetiInfo::noCorners;
            domain->solInfo().fetiInfo.pickAnyCorner = 0; 
            domain->solInfo().fetiInfo.bmpc = true;
            domain->solInfo().fetiInfo.pick_unsafe_corners = false;
            domain->solInfo().fetiInfo.augment = FetiInfo::none;
          }
        }
    break;

  case 730:

/* Line 1806 of yacc.c  */
#line 2837 "p.y"
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
        }
    break;

  case 731:

/* Line 1806 of yacc.c  */
#line 2858 "p.y"
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
        }
    break;

  case 732:

/* Line 1806 of yacc.c  */
#line 2876 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (4)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
    break;

  case 733:

/* Line 1806 of yacc.c  */
#line 2881 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (5)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
    break;

  case 734:

/* Line 1806 of yacc.c  */
#line 2887 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(4) - (5)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
    break;

  case 735:

/* Line 1806 of yacc.c  */
#line 2893 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (6)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(5) - (6)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (6)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  }
    break;

  case 736:

/* Line 1806 of yacc.c  */
#line 2900 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (3)].fval); }
    break;

  case 737:

/* Line 1806 of yacc.c  */
#line 2902 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.orthotol2 = (yyvsp[(3) - (4)].fval); }
    break;

  case 738:

/* Line 1806 of yacc.c  */
#line 2905 "p.y"
    { domain->solInfo().fetiInfo.grbm_tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 739:

/* Line 1806 of yacc.c  */
#line 2907 "p.y"
    { domain->solInfo().fetiInfo.crbm_tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 740:

/* Line 1806 of yacc.c  */
#line 2909 "p.y"
    { domain->solInfo().fetiInfo.cct_tol = (yyvsp[(2) - (3)].fval); }
    break;

  case 741:

/* Line 1806 of yacc.c  */
#line 2911 "p.y"
    { domain->solInfo().fetiInfo.rebuildcct = int((yyvsp[(2) - (3)].ival)); }
    break;

  case 742:

/* Line 1806 of yacc.c  */
#line 2913 "p.y"
    { domain->solInfo().fetiInfo.uproj = (yyvsp[(2) - (3)].ival); }
    break;

  case 743:

/* Line 1806 of yacc.c  */
#line 2915 "p.y"
    { domain->solInfo().fetiInfo.printMatLab = 1; }
    break;

  case 744:

/* Line 1806 of yacc.c  */
#line 2917 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); }
    break;

  case 745:

/* Line 1806 of yacc.c  */
#line 2919 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); }
    break;

  case 746:

/* Line 1806 of yacc.c  */
#line 2921 "p.y"
    {  domain->solInfo().fetiInfo.auxCoarseSolver = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); }
    break;

  case 747:

/* Line 1806 of yacc.c  */
#line 2923 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); }
    break;

  case 748:

/* Line 1806 of yacc.c  */
#line 2925 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
    break;

  case 749:

/* Line 1806 of yacc.c  */
#line 2929 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().localScaled = true; }
    break;

  case 750:

/* Line 1806 of yacc.c  */
#line 2932 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
    break;

  case 751:

/* Line 1806 of yacc.c  */
#line 2936 "p.y"
    { domain->solInfo().fetiInfo.auxCoarseSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
    break;

  case 752:

/* Line 1806 of yacc.c  */
#line 2940 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().coarseScaled = true; }
    break;

  case 753:

/* Line 1806 of yacc.c  */
#line 2943 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; }
    break;

  case 754:

/* Line 1806 of yacc.c  */
#line 2947 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival)!=0) fprintf(stderr," *** WARNING: Scaling not supported for this CCt solver \n");
          else domain->solInfo().fetiInfo.cctScaled = true; }
    break;

  case 755:

/* Line 1806 of yacc.c  */
#line 2951 "p.y"
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
	}
    break;

  case 756:

/* Line 1806 of yacc.c  */
#line 2977 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(1) - (2)].ival); }
    break;

  case 757:

/* Line 1806 of yacc.c  */
#line 2979 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = true; }
    break;

  case 758:

/* Line 1806 of yacc.c  */
#line 2981 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 759:

/* Line 1806 of yacc.c  */
#line 2983 "p.y"
    { domain->solInfo().fetiInfo.pickAnyCorner = (yyvsp[(2) - (3)].ival); }
    break;

  case 760:

/* Line 1806 of yacc.c  */
#line 2989 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().fetiInfo.nlPrecFlg = 1; 
          domain->solInfo().setKrylov(); 
        }
    break;

  case 761:

/* Line 1806 of yacc.c  */
#line 2994 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
	  domain->solInfo().fetiInfo.nlPrecFlg = (yyvsp[(2) - (3)].ival);
	  domain->solInfo().setKrylov();
	}
    break;

  case 762:

/* Line 1806 of yacc.c  */
#line 2999 "p.y"
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
        }
    break;

  case 763:

/* Line 1806 of yacc.c  */
#line 3012 "p.y"
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
        }
    break;

  case 764:

/* Line 1806 of yacc.c  */
#line 3026 "p.y"
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
        }
    break;

  case 765:

/* Line 1806 of yacc.c  */
#line 3043 "p.y"
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
        }
    break;

  case 766:

/* Line 1806 of yacc.c  */
#line 3059 "p.y"
    { domain->solInfo().type =(2); 
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
    break;

  case 767:

/* Line 1806 of yacc.c  */
#line 3069 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (3)].ival); }
    break;

  case 768:

/* Line 1806 of yacc.c  */
#line 3071 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.numcgm2 = (yyvsp[(3) - (4)].fval); }
    break;

  case 769:

/* Line 1806 of yacc.c  */
#line 3074 "p.y"
    { domain->solInfo().fetiInfo.tolcgm = (yyvsp[(2) - (3)].fval); }
    break;

  case 770:

/* Line 1806 of yacc.c  */
#line 3076 "p.y"
    { domain->solInfo().fetiInfo.spaceDimension = (yyvsp[(2) - (3)].ival); }
    break;

  case 771:

/* Line 1806 of yacc.c  */
#line 3078 "p.y"
    { domain->solInfo().fetiInfo.krylovtype = (yyvsp[(2) - (3)].ival); }
    break;

  case 772:

/* Line 1806 of yacc.c  */
#line 3080 "p.y"
    { domain->solInfo().fetiInfo.krylovtype =  (yyvsp[(2) - (3)].ival); }
    break;

  case 773:

/* Line 1806 of yacc.c  */
#line 3082 "p.y"
    { domain->solInfo().fetiInfo.lumpedinterface = 1; }
    break;

  case 774:

/* Line 1806 of yacc.c  */
#line 3084 "p.y"
    { domain->solInfo().fetiInfo.saveMemCoarse = 1; }
    break;

  case 775:

/* Line 1806 of yacc.c  */
#line 3086 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        }
    break;

  case 776:

/* Line 1806 of yacc.c  */
#line 3091 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        }
    break;

  case 777:

/* Line 1806 of yacc.c  */
#line 3097 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        }
    break;

  case 778:

/* Line 1806 of yacc.c  */
#line 3104 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (3)].ival); }
    break;

  case 779:

/* Line 1806 of yacc.c  */
#line 3106 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.complex_hermitian = true; }
    break;

  case 780:

/* Line 1806 of yacc.c  */
#line 3109 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); }
    break;

  case 781:

/* Line 1806 of yacc.c  */
#line 3111 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); }
    break;

  case 782:

/* Line 1806 of yacc.c  */
#line 3113 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); }
    break;

  case 783:

/* Line 1806 of yacc.c  */
#line 3115 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); }
    break;

  case 784:

/* Line 1806 of yacc.c  */
#line 3117 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); }
    break;

  case 785:

/* Line 1806 of yacc.c  */
#line 3120 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); }
    break;

  case 786:

/* Line 1806 of yacc.c  */
#line 3123 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); }
    break;

  case 787:

/* Line 1806 of yacc.c  */
#line 3126 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); }
    break;

  case 788:

/* Line 1806 of yacc.c  */
#line 3129 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); }
    break;

  case 789:

/* Line 1806 of yacc.c  */
#line 3133 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); }
    break;

  case 790:

/* Line 1806 of yacc.c  */
#line 3137 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 1) domain->solInfo().fetiInfo.useMRHS = false; }
    break;

  case 791:

/* Line 1806 of yacc.c  */
#line 3139 "p.y"
    { domain->solInfo().fetiInfo.gamma = (yyvsp[(2) - (3)].fval); }
    break;

  case 792:

/* Line 1806 of yacc.c  */
#line 3141 "p.y"
    { domain->solInfo().fetiInfo.linesearch_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.linesearch_tau = (yyvsp[(3) - (4)].fval); }
    break;

  case 793:

/* Line 1806 of yacc.c  */
#line 3144 "p.y"
    { domain->solInfo().fetiInfo.bmpc = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 794:

/* Line 1806 of yacc.c  */
#line 3146 "p.y"
    { domain->solInfo().fetiInfo.dmpc = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 795:

/* Line 1806 of yacc.c  */
#line 3148 "p.y"
    { domain->solInfo().fetiInfo.cmpc = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 796:

/* Line 1806 of yacc.c  */
#line 3150 "p.y"
    { domain->solInfo().fetiInfo.c_normalize = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 797:

/* Line 1806 of yacc.c  */
#line 3152 "p.y"
    { domain->solInfo().dbccheck = bool((yyvsp[(2) - (3)].ival)); }
    break;

  case 799:

/* Line 1806 of yacc.c  */
#line 3159 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp; 
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        }
    break;

  case 800:

/* Line 1806 of yacc.c  */
#line 3170 "p.y"
    { if(!(yyvsp[(2) - (3)].copt).lagrangeMult && (yyvsp[(2) - (3)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(2) - (3)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(2) - (3)].copt).penalty; }
    break;

  case 801:

/* Line 1806 of yacc.c  */
#line 3174 "p.y"
    { if(!(yyvsp[(3) - (4)].copt).lagrangeMult && (yyvsp[(3) - (4)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(3) - (4)].copt).penalty; }
    break;

  case 802:

/* Line 1806 of yacc.c  */
#line 3180 "p.y"
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        }
    break;

  case 803:

/* Line 1806 of yacc.c  */
#line 3187 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (2)].fval); }
    break;

  case 804:

/* Line 1806 of yacc.c  */
#line 3194 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (3)].fval);
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (3)].fval); }
    break;

  case 805:

/* Line 1806 of yacc.c  */
#line 3202 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (4)].fval);
          domain->solInfo().rhsZeroTol = (yyvsp[(4) - (4)].fval); }
    break;

  case 806:

/* Line 1806 of yacc.c  */
#line 3211 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (5)].fval);
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (5)].fval); 
          domain->solInfo().rhsZeroTol = (yyvsp[(4) - (5)].fval);
          domain->solInfo().inconsistentTol = (yyvsp[(5) - (5)].fval); }
    break;

  case 807:

/* Line 1806 of yacc.c  */
#line 3221 "p.y"
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
    break;

  case 808:

/* Line 1806 of yacc.c  */
#line 3227 "p.y"
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[(2) - (2)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
    break;

  case 809:

/* Line 1806 of yacc.c  */
#line 3233 "p.y"
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[(3) - (3)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; }
    break;

  case 810:

/* Line 1806 of yacc.c  */
#line 3239 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (3)].ival);
          (yyval.copt).constraint_hess_eps = 0; }
    break;

  case 811:

/* Line 1806 of yacc.c  */
#line 3242 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (4)].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[(4) - (4)].fval); }
    break;

  case 812:

/* Line 1806 of yacc.c  */
#line 3247 "p.y"
    { // hack??
	  domain->solInfo().acoustic = true; }
    break;

  case 813:

/* Line 1806 of yacc.c  */
#line 3250 "p.y"
    { domain->solInfo().type = (2); 
          domain->solInfo().fetiInfo.scaling = FetiInfo::tscaling; 
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners3;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true;
          domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          domain->solInfo().fetiInfo.rbmType = FetiInfo::None;
          domain->solInfo().fetiInfo.nGs = 0;
        }
    break;

  case 814:

/* Line 1806 of yacc.c  */
#line 3260 "p.y"
    {
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(2) - (3)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (3)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
        }
    break;

  case 815:

/* Line 1806 of yacc.c  */
#line 3267 "p.y"
    {
          // this is for coupled problem: need fluid density and angular frequency
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
          domain->fluidCelerity = (yyvsp[(3) - (5)].fval)/(yyvsp[(2) - (5)].fval);
          geoSource->setOmega((yyvsp[(3) - (5)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (5)].fval); sp.rho = (yyvsp[(4) - (5)].fval); geoSource->addMat(0,sp); // default flumat
        }
    break;

  case 816:

/* Line 1806 of yacc.c  */
#line 3275 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        }
    break;

  case 817:

/* Line 1806 of yacc.c  */
#line 3280 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        }
    break;

  case 818:

/* Line 1806 of yacc.c  */
#line 3286 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        }
    break;

  case 819:

/* Line 1806 of yacc.c  */
#line 3293 "p.y"
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        }
    break;

  case 820:

/* Line 1806 of yacc.c  */
#line 3298 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
    break;

  case 821:

/* Line 1806 of yacc.c  */
#line 3303 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        }
    break;

  case 822:

/* Line 1806 of yacc.c  */
#line 3308 "p.y"
    { 
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(4) - (7)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (7)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies1((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); 
        }
    break;

  case 823:

/* Line 1806 of yacc.c  */
#line 3316 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(4) - (7)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (7)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies2((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival));
        }
    break;

  case 824:

/* Line 1806 of yacc.c  */
#line 3324 "p.y"
    { 
           domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
           domain->fluidCelerity = 1.0; // defines the ratio omega/k
           geoSource->setOmega((yyvsp[(4) - (8)].fval)*domain->fluidCelerity);
           StructProp sp; sp.kappaHelm = (yyvsp[(4) - (8)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
           domain->addFrequencies((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].ival), (yyvsp[(7) - (8)].ival)); 
        }
    break;

  case 826:

/* Line 1806 of yacc.c  */
#line 3333 "p.y"
    {
          // coupled sweep:            k0    delta_k n       omega_0 rho
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(7) - (9)].fval)/(yyvsp[(4) - (9)].fval);
          geoSource->setOmega((yyvsp[(7) - (9)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(8) - (9)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies1((yyvsp[(4) - (9)].fval)*domain->fluidCelerity, (yyvsp[(5) - (9)].fval)*domain->fluidCelerity, (yyvsp[(6) - (9)].ival));
        }
    break;

  case 827:

/* Line 1806 of yacc.c  */
#line 3343 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(7) - (9)].fval)/(yyvsp[(4) - (9)].fval);
          geoSource->setOmega((yyvsp[(7) - (9)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(8) - (9)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies2((yyvsp[(4) - (9)].fval)*domain->fluidCelerity, (yyvsp[(5) - (9)].fval)*domain->fluidCelerity, (yyvsp[(6) - (9)].ival));
        }
    break;

  case 828:

/* Line 1806 of yacc.c  */
#line 3352 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(8) - (10)].fval)/(yyvsp[(4) - (10)].fval);
          geoSource->setOmega((yyvsp[(8) - (10)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (10)].fval); sp.rho = (yyvsp[(9) - (10)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies((yyvsp[(4) - (10)].fval)*domain->fluidCelerity, (yyvsp[(5) - (10)].fval)*domain->fluidCelerity, (yyvsp[(6) - (10)].ival), (yyvsp[(7) - (10)].ival));
        }
    break;

  case 832:

/* Line 1806 of yacc.c  */
#line 3366 "p.y"
    { 
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(2) - (3)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (3)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
          domain->addCoarseFrequency((yyvsp[(2) - (3)].fval)*domain->fluidCelerity); 
        }
    break;

  case 833:

/* Line 1806 of yacc.c  */
#line 3375 "p.y"
    { 
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(3) - (5)].fval)/(yyvsp[(2) - (5)].fval);
          geoSource->setOmega((yyvsp[(3) - (5)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (5)].fval); sp.rho = (yyvsp[(4) - (5)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addCoarseFrequency((yyvsp[(2) - (5)].fval)*domain->fluidCelerity);
        }
    break;

  case 834:

/* Line 1806 of yacc.c  */
#line 3384 "p.y"
    { domain->addFrequencies((yyvsp[(2) - (4)].fval)*domain->fluidCelerity, (yyvsp[(3) - (4)].ival)); }
    break;

  case 837:

/* Line 1806 of yacc.c  */
#line 3392 "p.y"
    { domain->setWaveDirections(0, (yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); }
    break;

  case 838:

/* Line 1806 of yacc.c  */
#line 3395 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        }
    break;

  case 839:

/* Line 1806 of yacc.c  */
#line 3398 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        }
    break;

  case 842:

/* Line 1806 of yacc.c  */
#line 3408 "p.y"
    { domain->setFFPDirections((yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); }
    break;

  case 844:

/* Line 1806 of yacc.c  */
#line 3415 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        }
    break;

  case 846:

/* Line 1806 of yacc.c  */
#line 3429 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        }
    break;

  case 847:

/* Line 1806 of yacc.c  */
#line 3440 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        }
    break;

  case 848:

/* Line 1806 of yacc.c  */
#line 3446 "p.y"
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
    break;

  case 849:

/* Line 1806 of yacc.c  */
#line 3459 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength;
        }
    break;

  case 850:

/* Line 1806 of yacc.c  */
#line 3464 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam;
        }
    break;

  case 851:

/* Line 1806 of yacc.c  */
#line 3471 "p.y"
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[(3) - (4)].ival); }
    break;

  case 852:

/* Line 1806 of yacc.c  */
#line 3473 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (4)].fval); }
    break;

  case 853:

/* Line 1806 of yacc.c  */
#line 3475 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (5)].fval); }
    break;

  case 854:

/* Line 1806 of yacc.c  */
#line 3478 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (7)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[(5) - (7)].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[(6) - (7)].fval); }
    break;

  case 855:

/* Line 1806 of yacc.c  */
#line 3483 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (5)].fval); }
    break;

  case 856:

/* Line 1806 of yacc.c  */
#line 3486 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (7)].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[(5) - (7)].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[(6) - (7)].ival); }
    break;

  case 857:

/* Line 1806 of yacc.c  */
#line 3491 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(3) - (4)].ival); }
    break;

  case 858:

/* Line 1806 of yacc.c  */
#line 3494 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (5)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(4) - (5)].ival); }
    break;

  case 859:

/* Line 1806 of yacc.c  */
#line 3497 "p.y"
    { domain->solInfo().getNLInfo().unsymmetric = true; }
    break;

  case 860:

/* Line 1806 of yacc.c  */
#line 3499 "p.y"
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[(3) - (4)].fval); }
    break;

  case 861:

/* Line 1806 of yacc.c  */
#line 3505 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true; }
    break;

  case 862:

/* Line 1806 of yacc.c  */
#line 3507 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[(3) - (4)].fval); }
    break;

  case 864:

/* Line 1806 of yacc.c  */
#line 3513 "p.y"
    { 
          domain->solInfo().setNewton((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear; 
        }
    break;

  case 865:

/* Line 1806 of yacc.c  */
#line 3518 "p.y"
    { 
          domain->solInfo().setNewton((yyvsp[(2) - (4)].ival)); 
          int rebuildK    = (yyvsp[(2) - (4)].ival); 
          int rebuildPrec = (yyvsp[(3) - (4)].ival);
          if(rebuildK > 1) rebuildPrec = rebuildK;
          domain->solInfo().fetiInfo.nTang = rebuildK;
          domain->solInfo().fetiInfo.nPrec = rebuildPrec;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        }
    break;

  case 866:

/* Line 1806 of yacc.c  */
#line 3528 "p.y"
    {
          domain->solInfo().setNewton((yyvsp[(4) - (5)].ival));
          domain->solInfo().fetiInfo.nPrec = (yyvsp[(4) - (5)].ival);
          domain->solInfo().fetiInfo.nTang = (yyvsp[(4) - (5)].ival);
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        }
    break;

  case 867:

/* Line 1806 of yacc.c  */
#line 3535 "p.y"
    {
	  domain->solInfo().setNewton((yyvsp[(4) - (8)].ival)); 
          int rebuildK    = (yyvsp[(4) - (8)].ival); 
          int rebuildPrec = (yyvsp[(7) - (8)].ival);
          if(rebuildK > 1) rebuildPrec = rebuildK;
          domain->solInfo().fetiInfo.nTang = rebuildK;
          domain->solInfo().fetiInfo.nPrec = rebuildPrec;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
	}
    break;

  case 868:

/* Line 1806 of yacc.c  */
#line 3547 "p.y"
    { domain->solInfo().setReOrtho(); }
    break;

  case 870:

/* Line 1806 of yacc.c  */
#line 3552 "p.y"
    { geoSource->setControl((yyvsp[(3) - (10)].strval),(yyvsp[(7) - (10)].strval),(yyvsp[(9) - (10)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (10)].ival); }
    break;

  case 871:

/* Line 1806 of yacc.c  */
#line 3554 "p.y"
    { geoSource->setControl((yyvsp[(3) - (12)].strval),(yyvsp[(7) - (12)].strval),(yyvsp[(9) - (12)].strval),(yyvsp[(11) - (12)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (12)].ival); }
    break;

  case 872:

/* Line 1806 of yacc.c  */
#line 3564 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(2) - (3)].strval);
#endif
        }
    break;

  case 873:

/* Line 1806 of yacc.c  */
#line 3570 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(3) - (4)].strval);
#endif
 }
    break;

  case 875:

/* Line 1806 of yacc.c  */
#line 3579 "p.y"
    { domain->solInfo().contact_mode = (yyvsp[(3) - (4)].ival); }
    break;

  case 876:

/* Line 1806 of yacc.c  */
#line 3583 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (7)].ival)-1, (yyvsp[(3) - (7)].ival)-1, (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); }
    break;

  case 877:

/* Line 1806 of yacc.c  */
#line 3586 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(8) - (9)].fval));}
    break;

  case 878:

/* Line 1806 of yacc.c  */
#line 3589 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), 0.0, (yyvsp[(8) - (9)].ival));}
    break;

  case 879:

/* Line 1806 of yacc.c  */
#line 3591 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (8)].ival)-1, (yyvsp[(3) - (8)].ival)-1, (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), 0.0, -1, (yyvsp[(7) - (8)].copt).lagrangeMult, (yyvsp[(7) - (8)].copt).penalty);}
    break;

  case 880:

/* Line 1806 of yacc.c  */
#line 3594 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (11)].ival)-1, (yyvsp[(3) - (11)].ival)-1, (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(10) - (11)].ival));}
    break;

  case 881:

/* Line 1806 of yacc.c  */
#line 3596 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (12)].ival)-1, (yyvsp[(3) - (12)].ival)-1, (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(10) - (12)].ival), (yyvsp[(11) - (12)].copt).lagrangeMult, (yyvsp[(11) - (12)].copt).penalty);}
    break;

  case 883:

/* Line 1806 of yacc.c  */
#line 3601 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         }
    break;

  case 884:

/* Line 1806 of yacc.c  */
#line 3606 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new BilinPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         }
    break;

  case 885:

/* Line 1806 of yacc.c  */
#line 3611 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         }
    break;

  case 886:

/* Line 1806 of yacc.c  */
#line 3616 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         }
    break;

  case 887:

/* Line 1806 of yacc.c  */
#line 3621 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         }
    break;

  case 888:

/* Line 1806 of yacc.c  */
#line 3626 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         }
    break;

  case 889:

/* Line 1806 of yacc.c  */
#line 3631 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1, 
             new ElaLinIsoMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
	 }
    break;

  case 890:

/* Line 1806 of yacc.c  */
#line 3636 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new StVenantKirchhoffMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
         }
    break;

  case 891:

/* Line 1806 of yacc.c  */
#line 3641 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new HenckyMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
         }
    break;

  case 892:

/* Line 1806 of yacc.c  */
#line 3646 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ElaLinIsoMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)));
         }
    break;

  case 893:

/* Line 1806 of yacc.c  */
#line 3651 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)));
         }
    break;

  case 894:

/* Line 1806 of yacc.c  */
#line 3656 "p.y"
    {
            double params[3] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          }
    break;

  case 895:

/* Line 1806 of yacc.c  */
#line 3662 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          }
    break;

  case 896:

/* Line 1806 of yacc.c  */
#line 3668 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          }
    break;

  case 897:

/* Line 1806 of yacc.c  */
#line 3674 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          }
    break;

  case 898:

/* Line 1806 of yacc.c  */
#line 3680 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          }
    break;

  case 899:

/* Line 1806 of yacc.c  */
#line 3686 "p.y"
    {
            double params[6] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          }
    break;

  case 900:

/* Line 1806 of yacc.c  */
#line 3692 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 1.0e-6, -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
    break;

  case 901:

/* Line 1806 of yacc.c  */
#line 3698 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
    break;

  case 902:

/* Line 1806 of yacc.c  */
#line 3704 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          }
    break;

  case 903:

/* Line 1806 of yacc.c  */
#line 3710 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (24)].ival)-1,
             new ExpMat((yyvsp[(3) - (24)].ival), (yyvsp[(4) - (24)].fval), (yyvsp[(5) - (24)].fval), (yyvsp[(6) - (24)].fval), (yyvsp[(7) - (24)].fval), (yyvsp[(8) - (24)].fval), (yyvsp[(9) - (24)].fval), (yyvsp[(10) - (24)].fval), (yyvsp[(11) - (24)].fval), (yyvsp[(12) - (24)].fval), (yyvsp[(13) - (24)].fval), (yyvsp[(14) - (24)].fval), (yyvsp[(15) - (24)].fval), (yyvsp[(16) - (24)].fval), (yyvsp[(17) - (24)].fval), (yyvsp[(18) - (24)].fval), (yyvsp[(19) - (24)].fval), (yyvsp[(20) - (24)].fval), (yyvsp[(21) - (24)].fval), (yyvsp[(22) - (24)].fval), (yyvsp[(23) - (24)].fval)));
         }
    break;

  case 904:

/* Line 1806 of yacc.c  */
#line 3715 "p.y"
    {
	   geoSource->loadMaterial((yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].strval));
	 }
    break;

  case 905:

/* Line 1806 of yacc.c  */
#line 3719 "p.y"
    {
	   geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].dlist));
	 }
    break;

  case 907:

/* Line 1806 of yacc.c  */
#line 3726 "p.y"
    { geoSource->setMatUsage((yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); }
    break;

  case 908:

/* Line 1806 of yacc.c  */
#line 3728 "p.y"
    {
            for(int i = (yyvsp[(2) - (5)].ival)-1; i < (yyvsp[(3) - (5)].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[(4) - (5)].ival)-1);
	  }
    break;

  case 909:

/* Line 1806 of yacc.c  */
#line 3734 "p.y"
    { (yyval.dlist).nval = 0; }
    break;

  case 910:

/* Line 1806 of yacc.c  */
#line 3736 "p.y"
    { 
          if((yyvsp[(1) - (2)].dlist).nval == 32) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[(1) - (2)].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[(2) - (2)].fval);
 	}
    break;

  case 911:

/* Line 1806 of yacc.c  */
#line 3747 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (4)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(3) - (4)].ival)); }
    break;

  case 912:

/* Line 1806 of yacc.c  */
#line 3751 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (6)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (6)].ival)); }
    break;

  case 913:

/* Line 1806 of yacc.c  */
#line 3754 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (8)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (8)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(7) - (8)].ival)); }
    break;

  case 914:

/* Line 1806 of yacc.c  */
#line 3761 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().svdPodRom = true;}
    break;

  case 916:

/* Line 1806 of yacc.c  */
#line 3769 "p.y"
    { domain->solInfo().snapfiPodRom = (yyvsp[(2) - (2)].strval); }
    break;

  case 917:

/* Line 1806 of yacc.c  */
#line 3771 "p.y"
    { domain->solInfo().snapfiPodRom = (yyvsp[(2) - (3)].strval);
    if ((yyvsp[(3) - (3)].ival) == 1) domain->solInfo().statevectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 2) domain->solInfo().residvectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 3) domain->solInfo().jacobvectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 4) domain->solInfo().forcevectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 5) domain->solInfo().accelvectPodRom = true;}
    break;

  case 918:

/* Line 1806 of yacc.c  */
#line 3778 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); }
    break;

  case 919:

/* Line 1806 of yacc.c  */
#line 3783 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().samplingPodRom = true; }
    break;

  case 921:

/* Line 1806 of yacc.c  */
#line 3791 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (2)].strval); }
    break;

  case 922:

/* Line 1806 of yacc.c  */
#line 3793 "p.y"
    { domain->solInfo().statePodRomFile = (yyvsp[(2) - (2)].strval); }
    break;

  case 923:

/* Line 1806 of yacc.c  */
#line 3795 "p.y"
    { domain->solInfo().tolPodRom = (yyvsp[(2) - (2)].fval); }
    break;

  case 924:

/* Line 1806 of yacc.c  */
#line 3797 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (2)].ival); }
    break;

  case 925:

/* Line 1806 of yacc.c  */
#line 3799 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); }
    break;

  case 926:

/* Line 1806 of yacc.c  */
#line 3804 "p.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); }
    break;

  case 927:

/* Line 1806 of yacc.c  */
#line 3809 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].ival); }
    break;

  case 928:

/* Line 1806 of yacc.c  */
#line 3811 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].fval); }
    break;

  case 929:

/* Line 1806 of yacc.c  */
#line 3813 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::infinity(); }
    break;

  case 930:

/* Line 1806 of yacc.c  */
#line 3815 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); }
    break;



/* Line 1806 of yacc.c  */
#line 11218 "/home/tac688/Research/FEM1/FEM3/Parser.d/parser.cpp"
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
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
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



/* Line 2067 of yacc.c  */
#line 3817 "p.y"


