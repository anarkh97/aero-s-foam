/* A Bison parser, made by GNU Bison 2.5.  */

/* Bison interface for Yacc-like parsers in C
   
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

/* Line 2068 of yacc.c  */
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



/* Line 2068 of yacc.c  */
#line 486 "/home/tac688/Research/FEM1/FEM3/Parser.d/parser.hpp"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


