/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

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
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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
     AMAT = 262,
     ANALYSIS = 263,
     ARCLENGTH = 264,
     ATTRIBUTES = 265,
     ANGULAROUTTYPE = 266,
     AUGMENT = 267,
     AUGMENTTYPE = 268,
     AVERAGED = 269,
     ATDARB = 270,
     ACOU = 271,
     ATDDNB = 272,
     ATDROB = 273,
     ARPACK = 274,
     ATDDIR = 275,
     ATDNEU = 276,
     AXIHDIR = 277,
     AXIHNEU = 278,
     AXINUMMODES = 279,
     AXINUMSLICES = 280,
     AXIHSOMMER = 281,
     AXIMPC = 282,
     AUXCOARSESOLVER = 283,
     ACMECNTL = 284,
     ADDEDMASS = 285,
     AEROEMBED = 286,
     AUGMENTED = 287,
     BLOCKDIAG = 288,
     BOFFSET = 289,
     BUCKLE = 290,
     BGTL = 291,
     BMPC = 292,
     BINARYINPUT = 293,
     BINARYOUTPUT = 294,
     BLOCKSIZE = 295,
     CHECKTOKEN = 296,
     COARSESOLVER = 297,
     COEF = 298,
     CFRAMES = 299,
     COLLOCATEDTYPE = 300,
     CONVECTION = 301,
     COMPOSITE = 302,
     CONDITION = 303,
     CONTROL = 304,
     CORNER = 305,
     CORNERTYPE = 306,
     CURVE = 307,
     CCTTOL = 308,
     CCTSOLVER = 309,
     CRHS = 310,
     COUPLEDSCALE = 311,
     CONTACTSURFACES = 312,
     CMPC = 313,
     CNORM = 314,
     COMPLEXOUTTYPE = 315,
     CONSTRMAT = 316,
     CASES = 317,
     CONSTRAINEDSURFACES = 318,
     CSFRAMES = 319,
     CSTYPE = 320,
     CONSTANT = 321,
     CONWEP = 322,
     DAMPING = 323,
     DblConstant = 324,
     DEM = 325,
     DIMASS = 326,
     DISP = 327,
     DIRECT = 328,
     DLAMBDA = 329,
     DP = 330,
     DYNAM = 331,
     DETER = 332,
     DECOMPOSE = 333,
     DECOMPFILE = 334,
     DMPC = 335,
     DEBUGCNTL = 336,
     DEBUGICNTL = 337,
     CONSTRAINTS = 338,
     MULTIPLIERS = 339,
     PENALTY = 340,
     ELLUMP = 341,
     EIGEN = 342,
     EFRAMES = 343,
     ELSCATTERER = 344,
     END = 345,
     ELHSOMMERFELD = 346,
     ETEMP = 347,
     EXPLICIT = 348,
     EPSILON = 349,
     ELEMENTARYFUNCTIONTYPE = 350,
     FABMAT = 351,
     FACOUSTICS = 352,
     FETI = 353,
     FETI2TYPE = 354,
     FETIPREC = 355,
     FFP = 356,
     FFPDIR = 357,
     FITALG = 358,
     FNAME = 359,
     FLUX = 360,
     FORCE = 361,
     FRONTAL = 362,
     FETIH = 363,
     FIELDWEIGHTLIST = 364,
     FILTEREIG = 365,
     FLUID = 366,
     FREQSWEEP = 367,
     FREQSWEEP1 = 368,
     FREQSWEEP2 = 369,
     FREQSWEEPA = 370,
     FSGL = 371,
     FSINTERFACE = 372,
     FSISCALING = 373,
     FSIELEMENT = 374,
     NOLOCALFSISPLITING = 375,
     FSICORNER = 376,
     FFIDEBUG = 377,
     FAILSAFE = 378,
     FRAMETYPE = 379,
     GEPS = 380,
     GLOBALTOL = 381,
     GRAVITY = 382,
     GRBM = 383,
     GTGSOLVER = 384,
     GLOBALCRBMTOL = 385,
     GROUP = 386,
     GROUPTYPE = 387,
     GOLDFARBTOL = 388,
     GOLDFARBCHECK = 389,
     HDIRICHLET = 390,
     HEAT = 391,
     HFETI = 392,
     HNEUMAN = 393,
     HSOMMERFELD = 394,
     HFTT = 395,
     HELMHOLTZ = 396,
     HNBO = 397,
     HELMMF = 398,
     HELMSO = 399,
     HSCBO = 400,
     HWIBO = 401,
     HZEM = 402,
     HZEMFILTER = 403,
     HLMPC = 404,
     HERMITIAN = 405,
     HESSIAN = 406,
     IACC = 407,
     IDENTITY = 408,
     IDIS = 409,
     IDIS6 = 410,
     IntConstant = 411,
     INTERFACELUMPED = 412,
     ITEMP = 413,
     ITERTYPE = 414,
     IVEL = 415,
     INCIDENCE = 416,
     IHDIRICHLET = 417,
     IHDSWEEP = 418,
     IHNEUMANN = 419,
     ISOLVERTYPE = 420,
     INPC = 421,
     INFINTY = 422,
     JACOBI = 423,
     KEYLETTER = 424,
     KRYLOVTYPE = 425,
     KIRLOC = 426,
     LAYC = 427,
     LAYN = 428,
     LAYD = 429,
     LAYO = 430,
     LAYMAT = 431,
     LFACTOR = 432,
     LMPC = 433,
     LOAD = 434,
     LOADCASE = 435,
     LOBPCG = 436,
     LOCALSOLVER = 437,
     LINESEARCH = 438,
     LUMPED = 439,
     MASS = 440,
     MATERIALS = 441,
     MATLAB = 442,
     MAXITR = 443,
     MAXORTHO = 444,
     MAXVEC = 445,
     MODAL = 446,
     MPCPRECNO = 447,
     MPCPRECNOID = 448,
     MPCTYPE = 449,
     MPCTYPEID = 450,
     MPCSCALING = 451,
     MPCELEMENT = 452,
     MPCBLOCKID = 453,
     MPCBLK_OVERLAP = 454,
     MFTT = 455,
     MRHS = 456,
     MPCCHECK = 457,
     MUMPSICNTL = 458,
     MUMPSCNTL = 459,
     MECH = 460,
     MODDAMP = 461,
     MODEFILTER = 462,
     MOMENTTYPE = 463,
     MPROJECT = 464,
     MAXIMUM = 465,
     NDTYPE = 466,
     NEIGPA = 467,
     NEWMARK = 468,
     NewLine = 469,
     NL = 470,
     NLMAT = 471,
     NLPREC = 472,
     NOCOARSE = 473,
     NODETOKEN = 474,
     NONINPC = 475,
     NSBSPV = 476,
     NLTOL = 477,
     NUMCGM = 478,
     NOSECONDARY = 479,
     NFRAMES = 480,
     OPTIMIZATION = 481,
     OUTPUT = 482,
     OUTPUT6 = 483,
     OUTPUTFRAME = 484,
     QSTATIC = 485,
     QLOAD = 486,
     PITA = 487,
     PITADISP6 = 488,
     PITAVEL6 = 489,
     NOFORCE = 490,
     MDPITA = 491,
     GLOBALBASES = 492,
     LOCALBASES = 493,
     TIMEREVERSIBLE = 494,
     REMOTECOARSE = 495,
     ORTHOPROJTOL = 496,
     READINITSEED = 497,
     JUMPCVG = 498,
     JUMPOUTPUT = 499,
     PRECNO = 500,
     PRECONDITIONER = 501,
     PRELOAD = 502,
     PRESSURE = 503,
     PRINTMATLAB = 504,
     PROJ = 505,
     PIVOT = 506,
     PRECTYPE = 507,
     PRECTYPEID = 508,
     PICKANYCORNER = 509,
     PADEPIVOT = 510,
     PROPORTIONING = 511,
     PLOAD = 512,
     PADEPOLES = 513,
     POINTSOURCE = 514,
     PLANEWAVE = 515,
     PTOL = 516,
     PLANTOL = 517,
     PMAXIT = 518,
     PIECEWISE = 519,
     RADIATION = 520,
     RAYDAMP = 521,
     RBMFILTER = 522,
     RBMSET = 523,
     READMODE = 524,
     REBUILD = 525,
     RENUM = 526,
     RENUMBERID = 527,
     REORTHO = 528,
     RESTART = 529,
     RECONS = 530,
     RECONSALG = 531,
     REBUILDCCT = 532,
     RANDOM = 533,
     RPROP = 534,
     RNORM = 535,
     REVERSENORMALS = 536,
     ROTVECOUTTYPE = 537,
     RESCALING = 538,
     SCALING = 539,
     SCALINGTYPE = 540,
     STRDAMP = 541,
     SDETAFT = 542,
     SENSORS = 543,
     SOLVERTYPE = 544,
     SHIFT = 545,
     SPOOLESTAU = 546,
     SPOOLESSEED = 547,
     SPOOLESMAXSIZE = 548,
     SPOOLESMAXDOMAINSIZE = 549,
     SPOOLESMAXZEROS = 550,
     SPOOLESMSGLVL = 551,
     SPOOLESSCALE = 552,
     SPOOLESPIVOT = 553,
     SPOOLESRENUM = 554,
     SPARSEMAXSUP = 555,
     SPARSEDEFBLK = 556,
     STATS = 557,
     STRESSID = 558,
     SUBSPACE = 559,
     SURFACE = 560,
     SAVEMEMCOARSE = 561,
     SPACEDIMENSION = 562,
     SCATTERER = 563,
     STAGTOL = 564,
     SCALED = 565,
     SWITCH = 566,
     STABLE = 567,
     SUBTYPE = 568,
     STEP = 569,
     SOWER = 570,
     SHELLTHICKNESS = 571,
     SURF = 572,
     SPRINGMAT = 573,
     TANGENT = 574,
     TEMP = 575,
     TIME = 576,
     TOLEIG = 577,
     TOLFETI = 578,
     TOLJAC = 579,
     TOLPCG = 580,
     TOPFILE = 581,
     TOPOLOGY = 582,
     TRBM = 583,
     THERMOE = 584,
     THERMOH = 585,
     TETT = 586,
     TOLCGM = 587,
     TURKEL = 588,
     TIEDSURFACES = 589,
     THETA = 590,
     REDFOL = 591,
     HRC = 592,
     THIRDNODE = 593,
     THERMMAT = 594,
     TDENFORC = 595,
     TESTULRICH = 596,
     THRU = 597,
     TRIVIAL = 598,
     USE = 599,
     USERDEFINEDISP = 600,
     USERDEFINEFORCE = 601,
     UPROJ = 602,
     UNSYMMETRIC = 603,
     USING = 604,
     VERSION = 605,
     WETCORNERS = 606,
     YMTT = 607,
     ZERO = 608,
     BINARY = 609,
     GEOMETRY = 610,
     DECOMPOSITION = 611,
     GLOBAL = 612,
     MATCHER = 613,
     CPUMAP = 614,
     NODALCONTACT = 615,
     MODE = 616,
     FRIC = 617,
     GAP = 618,
     OUTERLOOP = 619,
     EDGEWS = 620,
     WAVETYPE = 621,
     ORTHOTOL = 622,
     IMPE = 623,
     FREQ = 624,
     DPH = 625,
     WAVEMETHOD = 626,
     MATSPEC = 627,
     MATUSAGE = 628,
     BILINEARPLASTIC = 629,
     FINITESTRAINPLASTIC = 630,
     LINEARELASTIC = 631,
     STVENANTKIRCHHOFF = 632,
     LINPLSTRESS = 633,
     READ = 634,
     OPTCTV = 635,
     ISOTROPICLINEARELASTIC = 636,
     NEOHOOKEAN = 637,
     ISOTROPICLINEARELASTICJ2PLASTIC = 638,
     ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 639,
     HYPERELASTIC = 640,
     MOONEYRIVLIN = 641,
     HENCKY = 642,
     LOGSTRAINPLASTIC = 643,
     SVKPLSTRESS = 644,
     SURFACETOPOLOGY = 645,
     MORTARTIED = 646,
     MORTARSCALING = 647,
     MORTARINTEGRATIONRULE = 648,
     SEARCHTOL = 649,
     STDMORTAR = 650,
     DUALMORTAR = 651,
     WETINTERFACE = 652,
     NSUBS = 653,
     EXITAFTERDEC = 654,
     SKIP = 655,
     OUTPUTMEMORY = 656,
     OUTPUTWEIGHT = 657,
     WEIGHTLIST = 658,
     GMRESRESIDUAL = 659,
     SLOSH = 660,
     SLGRAV = 661,
     SLZEM = 662,
     SLZEMFILTER = 663,
     PDIR = 664,
     HEFSB = 665,
     HEFRS = 666,
     HEINTERFACE = 667,
     SNAPFI = 668,
     PODROB = 669,
     TRNVCT = 670,
     OFFSET = 671,
     ORTHOG = 672,
     SVDTOKEN = 673,
     CONVERSIONTOKEN = 674,
     CONVFI = 675,
     SAMPLING = 676,
     SNAPSHOTPROJECT = 677,
     PODSIZEMAX = 678,
     REFSUBSTRACT = 679,
     TOLER = 680,
     OUTOFCORE = 681,
     NORMALIZETOKEN = 682,
     FNUMBER = 683,
     SNAPWEIGHT = 684,
     ROBFI = 685,
     STAVCT = 686,
     VELVCT = 687,
     ACCVCT = 688,
     CONWEPCFG = 689,
     PSEUDOGNAT = 690,
     PSEUDOGNATELEM = 691,
     VECTORNORM = 692,
     LOCALTOLERANCE = 693,
     REBUILDFORCE = 694,
     SAMPNODESLOT = 695,
     REDUCEDSTIFFNESS = 696,
     UDEIMBASIS = 697,
     FORCEROB = 698,
     DEIMINDICES = 699,
     UDEIMINDICES = 700,
     SVDFORCESNAP = 701,
     USEMASSNORMALIZEDBASIS = 702,
     OPTSENSITIVITY = 703,
     SENSITIVITYID = 704,
     SENSITIVITYTYPE = 705,
     SENSITIVITYMETHOD = 706,
     QRFACTORIZATION = 707,
     QMATRIX = 708,
     RMATRIX = 709,
     XMATRIX = 710
   };
#endif
/* Tokens.  */
#define ACTUATORS 258
#define AERO 259
#define AEROH 260
#define AEROTYPE 261
#define AMAT 262
#define ANALYSIS 263
#define ARCLENGTH 264
#define ATTRIBUTES 265
#define ANGULAROUTTYPE 266
#define AUGMENT 267
#define AUGMENTTYPE 268
#define AVERAGED 269
#define ATDARB 270
#define ACOU 271
#define ATDDNB 272
#define ATDROB 273
#define ARPACK 274
#define ATDDIR 275
#define ATDNEU 276
#define AXIHDIR 277
#define AXIHNEU 278
#define AXINUMMODES 279
#define AXINUMSLICES 280
#define AXIHSOMMER 281
#define AXIMPC 282
#define AUXCOARSESOLVER 283
#define ACMECNTL 284
#define ADDEDMASS 285
#define AEROEMBED 286
#define AUGMENTED 287
#define BLOCKDIAG 288
#define BOFFSET 289
#define BUCKLE 290
#define BGTL 291
#define BMPC 292
#define BINARYINPUT 293
#define BINARYOUTPUT 294
#define BLOCKSIZE 295
#define CHECKTOKEN 296
#define COARSESOLVER 297
#define COEF 298
#define CFRAMES 299
#define COLLOCATEDTYPE 300
#define CONVECTION 301
#define COMPOSITE 302
#define CONDITION 303
#define CONTROL 304
#define CORNER 305
#define CORNERTYPE 306
#define CURVE 307
#define CCTTOL 308
#define CCTSOLVER 309
#define CRHS 310
#define COUPLEDSCALE 311
#define CONTACTSURFACES 312
#define CMPC 313
#define CNORM 314
#define COMPLEXOUTTYPE 315
#define CONSTRMAT 316
#define CASES 317
#define CONSTRAINEDSURFACES 318
#define CSFRAMES 319
#define CSTYPE 320
#define CONSTANT 321
#define CONWEP 322
#define DAMPING 323
#define DblConstant 324
#define DEM 325
#define DIMASS 326
#define DISP 327
#define DIRECT 328
#define DLAMBDA 329
#define DP 330
#define DYNAM 331
#define DETER 332
#define DECOMPOSE 333
#define DECOMPFILE 334
#define DMPC 335
#define DEBUGCNTL 336
#define DEBUGICNTL 337
#define CONSTRAINTS 338
#define MULTIPLIERS 339
#define PENALTY 340
#define ELLUMP 341
#define EIGEN 342
#define EFRAMES 343
#define ELSCATTERER 344
#define END 345
#define ELHSOMMERFELD 346
#define ETEMP 347
#define EXPLICIT 348
#define EPSILON 349
#define ELEMENTARYFUNCTIONTYPE 350
#define FABMAT 351
#define FACOUSTICS 352
#define FETI 353
#define FETI2TYPE 354
#define FETIPREC 355
#define FFP 356
#define FFPDIR 357
#define FITALG 358
#define FNAME 359
#define FLUX 360
#define FORCE 361
#define FRONTAL 362
#define FETIH 363
#define FIELDWEIGHTLIST 364
#define FILTEREIG 365
#define FLUID 366
#define FREQSWEEP 367
#define FREQSWEEP1 368
#define FREQSWEEP2 369
#define FREQSWEEPA 370
#define FSGL 371
#define FSINTERFACE 372
#define FSISCALING 373
#define FSIELEMENT 374
#define NOLOCALFSISPLITING 375
#define FSICORNER 376
#define FFIDEBUG 377
#define FAILSAFE 378
#define FRAMETYPE 379
#define GEPS 380
#define GLOBALTOL 381
#define GRAVITY 382
#define GRBM 383
#define GTGSOLVER 384
#define GLOBALCRBMTOL 385
#define GROUP 386
#define GROUPTYPE 387
#define GOLDFARBTOL 388
#define GOLDFARBCHECK 389
#define HDIRICHLET 390
#define HEAT 391
#define HFETI 392
#define HNEUMAN 393
#define HSOMMERFELD 394
#define HFTT 395
#define HELMHOLTZ 396
#define HNBO 397
#define HELMMF 398
#define HELMSO 399
#define HSCBO 400
#define HWIBO 401
#define HZEM 402
#define HZEMFILTER 403
#define HLMPC 404
#define HERMITIAN 405
#define HESSIAN 406
#define IACC 407
#define IDENTITY 408
#define IDIS 409
#define IDIS6 410
#define IntConstant 411
#define INTERFACELUMPED 412
#define ITEMP 413
#define ITERTYPE 414
#define IVEL 415
#define INCIDENCE 416
#define IHDIRICHLET 417
#define IHDSWEEP 418
#define IHNEUMANN 419
#define ISOLVERTYPE 420
#define INPC 421
#define INFINTY 422
#define JACOBI 423
#define KEYLETTER 424
#define KRYLOVTYPE 425
#define KIRLOC 426
#define LAYC 427
#define LAYN 428
#define LAYD 429
#define LAYO 430
#define LAYMAT 431
#define LFACTOR 432
#define LMPC 433
#define LOAD 434
#define LOADCASE 435
#define LOBPCG 436
#define LOCALSOLVER 437
#define LINESEARCH 438
#define LUMPED 439
#define MASS 440
#define MATERIALS 441
#define MATLAB 442
#define MAXITR 443
#define MAXORTHO 444
#define MAXVEC 445
#define MODAL 446
#define MPCPRECNO 447
#define MPCPRECNOID 448
#define MPCTYPE 449
#define MPCTYPEID 450
#define MPCSCALING 451
#define MPCELEMENT 452
#define MPCBLOCKID 453
#define MPCBLK_OVERLAP 454
#define MFTT 455
#define MRHS 456
#define MPCCHECK 457
#define MUMPSICNTL 458
#define MUMPSCNTL 459
#define MECH 460
#define MODDAMP 461
#define MODEFILTER 462
#define MOMENTTYPE 463
#define MPROJECT 464
#define MAXIMUM 465
#define NDTYPE 466
#define NEIGPA 467
#define NEWMARK 468
#define NewLine 469
#define NL 470
#define NLMAT 471
#define NLPREC 472
#define NOCOARSE 473
#define NODETOKEN 474
#define NONINPC 475
#define NSBSPV 476
#define NLTOL 477
#define NUMCGM 478
#define NOSECONDARY 479
#define NFRAMES 480
#define OPTIMIZATION 481
#define OUTPUT 482
#define OUTPUT6 483
#define OUTPUTFRAME 484
#define QSTATIC 485
#define QLOAD 486
#define PITA 487
#define PITADISP6 488
#define PITAVEL6 489
#define NOFORCE 490
#define MDPITA 491
#define GLOBALBASES 492
#define LOCALBASES 493
#define TIMEREVERSIBLE 494
#define REMOTECOARSE 495
#define ORTHOPROJTOL 496
#define READINITSEED 497
#define JUMPCVG 498
#define JUMPOUTPUT 499
#define PRECNO 500
#define PRECONDITIONER 501
#define PRELOAD 502
#define PRESSURE 503
#define PRINTMATLAB 504
#define PROJ 505
#define PIVOT 506
#define PRECTYPE 507
#define PRECTYPEID 508
#define PICKANYCORNER 509
#define PADEPIVOT 510
#define PROPORTIONING 511
#define PLOAD 512
#define PADEPOLES 513
#define POINTSOURCE 514
#define PLANEWAVE 515
#define PTOL 516
#define PLANTOL 517
#define PMAXIT 518
#define PIECEWISE 519
#define RADIATION 520
#define RAYDAMP 521
#define RBMFILTER 522
#define RBMSET 523
#define READMODE 524
#define REBUILD 525
#define RENUM 526
#define RENUMBERID 527
#define REORTHO 528
#define RESTART 529
#define RECONS 530
#define RECONSALG 531
#define REBUILDCCT 532
#define RANDOM 533
#define RPROP 534
#define RNORM 535
#define REVERSENORMALS 536
#define ROTVECOUTTYPE 537
#define RESCALING 538
#define SCALING 539
#define SCALINGTYPE 540
#define STRDAMP 541
#define SDETAFT 542
#define SENSORS 543
#define SOLVERTYPE 544
#define SHIFT 545
#define SPOOLESTAU 546
#define SPOOLESSEED 547
#define SPOOLESMAXSIZE 548
#define SPOOLESMAXDOMAINSIZE 549
#define SPOOLESMAXZEROS 550
#define SPOOLESMSGLVL 551
#define SPOOLESSCALE 552
#define SPOOLESPIVOT 553
#define SPOOLESRENUM 554
#define SPARSEMAXSUP 555
#define SPARSEDEFBLK 556
#define STATS 557
#define STRESSID 558
#define SUBSPACE 559
#define SURFACE 560
#define SAVEMEMCOARSE 561
#define SPACEDIMENSION 562
#define SCATTERER 563
#define STAGTOL 564
#define SCALED 565
#define SWITCH 566
#define STABLE 567
#define SUBTYPE 568
#define STEP 569
#define SOWER 570
#define SHELLTHICKNESS 571
#define SURF 572
#define SPRINGMAT 573
#define TANGENT 574
#define TEMP 575
#define TIME 576
#define TOLEIG 577
#define TOLFETI 578
#define TOLJAC 579
#define TOLPCG 580
#define TOPFILE 581
#define TOPOLOGY 582
#define TRBM 583
#define THERMOE 584
#define THERMOH 585
#define TETT 586
#define TOLCGM 587
#define TURKEL 588
#define TIEDSURFACES 589
#define THETA 590
#define REDFOL 591
#define HRC 592
#define THIRDNODE 593
#define THERMMAT 594
#define TDENFORC 595
#define TESTULRICH 596
#define THRU 597
#define TRIVIAL 598
#define USE 599
#define USERDEFINEDISP 600
#define USERDEFINEFORCE 601
#define UPROJ 602
#define UNSYMMETRIC 603
#define USING 604
#define VERSION 605
#define WETCORNERS 606
#define YMTT 607
#define ZERO 608
#define BINARY 609
#define GEOMETRY 610
#define DECOMPOSITION 611
#define GLOBAL 612
#define MATCHER 613
#define CPUMAP 614
#define NODALCONTACT 615
#define MODE 616
#define FRIC 617
#define GAP 618
#define OUTERLOOP 619
#define EDGEWS 620
#define WAVETYPE 621
#define ORTHOTOL 622
#define IMPE 623
#define FREQ 624
#define DPH 625
#define WAVEMETHOD 626
#define MATSPEC 627
#define MATUSAGE 628
#define BILINEARPLASTIC 629
#define FINITESTRAINPLASTIC 630
#define LINEARELASTIC 631
#define STVENANTKIRCHHOFF 632
#define LINPLSTRESS 633
#define READ 634
#define OPTCTV 635
#define ISOTROPICLINEARELASTIC 636
#define NEOHOOKEAN 637
#define ISOTROPICLINEARELASTICJ2PLASTIC 638
#define ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS 639
#define HYPERELASTIC 640
#define MOONEYRIVLIN 641
#define HENCKY 642
#define LOGSTRAINPLASTIC 643
#define SVKPLSTRESS 644
#define SURFACETOPOLOGY 645
#define MORTARTIED 646
#define MORTARSCALING 647
#define MORTARINTEGRATIONRULE 648
#define SEARCHTOL 649
#define STDMORTAR 650
#define DUALMORTAR 651
#define WETINTERFACE 652
#define NSUBS 653
#define EXITAFTERDEC 654
#define SKIP 655
#define OUTPUTMEMORY 656
#define OUTPUTWEIGHT 657
#define WEIGHTLIST 658
#define GMRESRESIDUAL 659
#define SLOSH 660
#define SLGRAV 661
#define SLZEM 662
#define SLZEMFILTER 663
#define PDIR 664
#define HEFSB 665
#define HEFRS 666
#define HEINTERFACE 667
#define SNAPFI 668
#define PODROB 669
#define TRNVCT 670
#define OFFSET 671
#define ORTHOG 672
#define SVDTOKEN 673
#define CONVERSIONTOKEN 674
#define CONVFI 675
#define SAMPLING 676
#define SNAPSHOTPROJECT 677
#define PODSIZEMAX 678
#define REFSUBSTRACT 679
#define TOLER 680
#define OUTOFCORE 681
#define NORMALIZETOKEN 682
#define FNUMBER 683
#define SNAPWEIGHT 684
#define ROBFI 685
#define STAVCT 686
#define VELVCT 687
#define ACCVCT 688
#define CONWEPCFG 689
#define PSEUDOGNAT 690
#define PSEUDOGNATELEM 691
#define VECTORNORM 692
#define LOCALTOLERANCE 693
#define REBUILDFORCE 694
#define SAMPNODESLOT 695
#define REDUCEDSTIFFNESS 696
#define UDEIMBASIS 697
#define FORCEROB 698
#define DEIMINDICES 699
#define UDEIMINDICES 700
#define SVDFORCESNAP 701
#define USEMASSNORMALIZEDBASIS 702
#define OPTSENSITIVITY 703
#define SENSITIVITYID 704
#define SENSITIVITYTYPE 705
#define SENSITIVITYMETHOD 706
#define QRFACTORIZATION 707
#define QMATRIX 708
#define RMATRIX 709
#define XMATRIX 710




/* Copy the first part of user declarations.  */
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

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 24 "p.y"
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
}
/* Line 193 of yacc.c.  */
#line 1065 "/lustre/home/tac688/newCodes/Empirical_GNAT/Parser.d/parser.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 1078 "/lustre/home/tac688/newCodes/Empirical_GNAT/Parser.d/parser.cpp"

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
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
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
  yytype_int16 yyss;
  YYSTYPE yyvs;
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
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  523
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   6038

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  456
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  239
/* YYNRULES -- Number of rules.  */
#define YYNRULES  1080
/* YYNRULES -- Number of states.  */
#define YYNSTATES  2649

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   710

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
     455
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
     506,   511,   517,   522,   527,   532,   537,   542,   547,   551,
     554,   559,   564,   568,   572,   576,   580,   584,   588,   592,
     595,   600,   603,   608,   613,   618,   623,   626,   630,   635,
     638,   642,   647,   650,   654,   659,   665,   671,   674,   677,
     680,   683,   686,   689,   692,   695,   700,   705,   711,   715,
     718,   722,   725,   729,   732,   736,   739,   751,   765,   780,
     786,   789,   792,   800,   810,   822,   835,   838,   844,   851,
     854,   860,   863,   869,   874,   879,   883,   887,   891,   895,
     899,   903,   906,   910,   914,   917,   920,   924,   928,   932,
     936,   941,   946,   950,   954,   960,   965,   971,   978,   986,
     990,   996,  1000,  1006,  1011,  1017,  1024,  1032,  1036,  1039,
    1043,  1046,  1049,  1053,  1056,  1059,  1062,  1066,  1069,  1073,
    1076,  1079,  1082,  1085,  1088,  1090,  1093,  1097,  1101,  1105,
    1108,  1114,  1118,  1122,  1126,  1129,  1133,  1136,  1140,  1145,
    1150,  1156,  1160,  1164,  1167,  1172,  1175,  1179,  1182,  1186,
    1189,  1195,  1198,  1201,  1204,  1207,  1210,  1213,  1217,  1222,
    1227,  1236,  1241,  1246,  1250,  1254,  1260,  1265,  1271,  1274,
    1278,  1282,  1286,  1289,  1294,  1296,  1300,  1302,  1305,  1308,
    1311,  1315,  1321,  1328,  1335,  1339,  1344,  1351,  1357,  1361,
    1366,  1370,  1372,  1375,  1382,  1385,  1388,  1391,  1394,  1397,
    1402,  1408,  1413,  1416,  1420,  1423,  1427,  1430,  1435,  1437,
    1440,  1443,  1446,  1452,  1458,  1465,  1473,  1477,  1479,  1481,
    1484,  1486,  1488,  1490,  1493,  1495,  1497,  1500,  1502,  1507,
    1513,  1517,  1521,  1528,  1534,  1541,  1545,  1548,  1554,  1561,
    1564,  1567,  1575,  1585,  1589,  1592,  1599,  1603,  1605,  1608,
    1612,  1615,  1619,  1621,  1624,  1629,  1632,  1637,  1643,  1646,
    1650,  1655,  1661,  1664,  1668,  1675,  1678,  1682,  1689,  1693,
    1695,  1698,  1703,  1709,  1716,  1720,  1723,  1728,  1730,  1733,
    1737,  1739,  1742,  1747,  1753,  1760,  1764,  1767,  1772,  1776,
    1779,  1784,  1788,  1791,  1796,  1800,  1803,  1808,  1812,  1814,
    1817,  1821,  1825,  1830,  1834,  1837,  1841,  1844,  1850,  1853,
    1857,  1862,  1867,  1877,  1880,  1883,  1893,  1896,  1900,  1905,
    1909,  1913,  1921,  1924,  1928,  1933,  1936,  1939,  1945,  1953,
    1964,  1968,  1973,  1978,  1984,  1989,  1992,  1996,  2000,  2005,
    2009,  2012,  2022,  2029,  2032,  2035,  2039,  2049,  2053,  2063,
    2066,  2070,  2075,  2079,  2083,  2086,  2090,  2093,  2101,  2111,
    2120,  2131,  2135,  2141,  2143,  2146,  2153,  2161,  2170,  2180,
    2182,  2185,  2187,  2190,  2193,  2197,  2204,  2209,  2217,  2220,
    2224,  2231,  2236,  2244,  2247,  2251,  2258,  2263,  2271,  2274,
    2278,  2281,  2284,  2288,  2291,  2295,  2301,  2306,  2313,  2318,
    2321,  2325,  2328,  2331,  2335,  2341,  2345,  2348,  2354,  2359,
    2363,  2368,  2370,  2373,  2377,  2380,  2397,  2417,  2437,  2458,
    2482,  2506,  2516,  2531,  2548,  2579,  2592,  2606,  2621,  2627,
    2634,  2651,  2664,  2668,  2673,  2679,  2686,  2693,  2701,  2713,
    2722,  2732,  2743,  2753,  2764,  2776,  2787,  2799,  2812,  2824,
    2837,  2851,  2857,  2864,  2871,  2879,  2887,  2896,  2901,  2905,
    2908,  2912,  2917,  2923,  2930,  2936,  2941,  2947,  2953,  2960,
    2968,  2976,  2985,  2993,  3002,  3007,  3011,  3014,  3020,  3027,
    3035,  3044,  3055,  3067,  3074,  3084,  3087,  3093,  3101,  3105,
    3108,  3113,  3116,  3122,  3130,  3133,  3139,  3146,  3154,  3163,
    3174,  3186,  3200,  3215,  3222,  3232,  3236,  3240,  3244,  3248,
    3252,  3255,  3261,  3266,  3270,  3278,  3285,  3290,  3292,  3295,
    3300,  3304,  3310,  3315,  3319,  3323,  3329,  3334,  3337,  3340,
    3343,  3346,  3358,  3363,  3366,  3369,  3381,  3396,  3409,  3425,
    3428,  3436,  3439,  3444,  3451,  3459,  3466,  3475,  3483,  3489,
    3495,  3503,  3512,  3515,  3520,  3524,  3527,  3531,  3534,  3539,
    3543,  3546,  3551,  3557,  3560,  3564,  3569,  3575,  3581,  3587,
    3594,  3601,  3604,  3608,  3613,  3616,  3621,  3628,  3635,  3644,
    3647,  3650,  3653,  3658,  3662,  3668,  3670,  3673,  3676,  3681,
    3686,  3691,  3696,  3701,  3704,  3708,  3711,  3715,  3719,  3723,
    3728,  3734,  3741,  3749,  3754,  3758,  3762,  3766,  3771,  3777,
    3780,  3785,  3789,  3793,  3797,  3801,  3805,  3809,  3813,  3817,
    3821,  3825,  3829,  3834,  3839,  3843,  3847,  3852,  3857,  3862,
    3867,  3872,  3877,  3881,  3885,  3889,  3894,  3898,  3903,  3908,
    3913,  3918,  3923,  3926,  3930,  3934,  3938,  3942,  3946,  3950,
    3954,  3957,  3960,  3964,  3967,  3971,  3974,  3978,  3983,  3987,
    3991,  3996,  4001,  4007,  4013,  4020,  4024,  4029,  4033,  4037,
    4041,  4045,  4049,  4052,  4056,  4060,  4064,  4068,  4072,  4077,
    4082,  4087,  4092,  4097,  4102,  4107,  4111,  4114,  4117,  4121,
    4125,  4128,  4132,  4137,  4143,  4150,  4158,  4161,  4165,  4170,
    4174,  4178,  4182,  4186,  4189,  4192,  4196,  4201,  4207,  4211,
    4216,  4220,  4224,  4228,  4232,  4238,  4244,  4249,  4254,  4261,
    4268,  4272,  4276,  4281,  4285,  4289,  4293,  4297,  4301,  4305,
    4308,  4312,  4317,  4319,  4322,  4326,  4331,  4337,  4339,  4342,
    4346,  4349,  4353,  4358,  4361,  4364,  4368,  4373,  4379,  4384,
    4389,  4394,  4397,  4400,  4403,  4405,  4408,  4413,  4420,  4426,
    4428,  4431,  4436,  4440,  4443,  4447,  4450,  4453,  4456,  4460,
    4464,  4468,  4473,  4478,  4484,  4492,  4498,  4506,  4511,  4517,
    4521,  4526,  4531,  4537,  4544,  4552,  4561,  4565,  4570,  4577,
    4580,  4584,  4589,  4592,  4595,  4606,  4619,  4623,  4628,  4631,
    4636,  4644,  4654,  4664,  4673,  4685,  4698,  4701,  4711,  4722,
    4732,  4743,  4753,  4764,  4772,  4780,  4788,  4797,  4806,  4814,
    4822,  4831,  4840,  4850,  4861,  4872,  4884,  4897,  4922,  4930,
    4939,  4949,  4960,  4972,  4978,  4984,  4987,  4992,  4998,  4999,
    5002,  5003,  5006,  5011,  5018,  5027,  5030,  5034,  5037,  5040,
    5043,  5046,  5049,  5052,  5055,  5057,  5060,  5064,  5067,  5071,
    5074,  5078,  5081,  5084,  5088,  5091,  5094,  5098,  5103,  5106,
    5109,  5112,  5115,  5119,  5124,  5127,  5130,  5133,  5136,  5139,
    5143,  5148,  5151,  5154,  5158,  5163,  5166,  5169,  5172,  5175,
    5178,  5182,  5185,  5189,  5192,  5196,  5198,  5200,  5202,  5204,
    5206
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     457,     0,    -1,   458,    90,    -1,   459,    -1,   458,   459,
      -1,   631,    -1,   541,    -1,   597,    -1,   598,    -1,   608,
      -1,   612,    -1,   620,    -1,   640,    -1,   642,    -1,   639,
      -1,   645,    -1,   646,    -1,   649,    -1,   647,    -1,   648,
      -1,   618,    -1,   653,    -1,   650,    -1,   651,    -1,   652,
      -1,   682,    -1,   590,    -1,   589,    -1,   592,    -1,   593,
      -1,   591,    -1,   594,    -1,   595,    -1,   596,    -1,   494,
      -1,   495,    -1,   497,    -1,   496,    -1,   502,    -1,   508,
      -1,   503,    -1,   626,    -1,   513,    -1,   500,    -1,   672,
      -1,   489,    -1,   477,    -1,   478,    -1,   487,    -1,   474,
      -1,   475,    -1,   476,    -1,   602,    -1,   604,    -1,   606,
      -1,   524,    -1,   525,    -1,   490,    -1,   527,    -1,   529,
      -1,   530,    -1,   526,    -1,   491,    -1,   492,    -1,   493,
      -1,   675,    -1,   676,    -1,   470,    -1,   674,    -1,   516,
      -1,   518,    -1,   519,    -1,   520,    -1,   522,    -1,   523,
      -1,   505,    -1,   504,    -1,   549,    -1,   550,    -1,   551,
      -1,   543,    -1,   546,    -1,   552,    -1,   506,    -1,   657,
      -1,   671,    -1,   661,    -1,   667,    -1,   669,    -1,   559,
      -1,   577,    -1,   578,    -1,   664,    -1,   553,    -1,   556,
      -1,   562,    -1,   564,    -1,   566,    -1,   568,    -1,   616,
      -1,   540,    -1,   536,    -1,   537,    -1,   538,    -1,   579,
      -1,   581,    -1,   583,    -1,   584,    -1,   585,    -1,   587,
      -1,   469,    -1,   572,    -1,   573,    -1,   574,    -1,   575,
      -1,   576,    -1,   471,    -1,   472,    -1,   473,    -1,   677,
      -1,   521,    -1,   460,    -1,   461,    -1,   462,    -1,   463,
      -1,   464,    -1,   678,    -1,   679,    -1,   621,    -1,   622,
      -1,   624,    -1,   623,    -1,   625,    -1,   628,    -1,   629,
      -1,   542,    -1,   644,    -1,   532,    -1,   630,    -1,   659,
      -1,   683,    -1,   685,    -1,   686,    -1,   687,    -1,   688,
      -1,   691,    -1,   498,    -1,   220,   214,   693,   693,   214,
      -1,   166,   214,   693,   214,    -1,   131,   214,    -1,   462,
     132,   693,   693,   214,    -1,   462,   132,   693,   693,   693,
     214,    -1,   462,   132,   317,   693,   693,   214,    -1,   278,
     214,    -1,   463,   693,   279,   694,   694,   214,    -1,   368,
     214,   369,   694,   214,    -1,   368,   214,   113,   694,   694,
     693,   214,    -1,   368,   214,   114,   694,   694,   693,   214,
      -1,   368,   214,   112,   694,   694,   693,   693,   214,    -1,
     368,   214,   115,   694,   694,   693,   276,   694,   693,   693,
     693,   693,   214,    -1,   368,   214,   115,   694,   694,   693,
     276,   214,    -1,   368,   693,   693,   214,   369,   694,   214,
      -1,   368,   693,   214,   113,   694,   694,   693,   214,    -1,
     368,   693,   214,   114,   694,   694,   693,   214,    -1,   368,
     693,   214,   112,   694,   694,   693,   693,   214,    -1,   368,
     693,   214,   115,   694,   694,   693,   276,   694,   693,   693,
     693,   693,   214,    -1,   368,   693,   214,   115,   694,   694,
     693,   276,   214,    -1,   368,   214,   467,    -1,   464,   468,
      -1,   464,   535,    -1,   464,   465,    -1,   464,   466,    -1,
     255,   694,   214,    -1,   258,   214,    -1,   258,   694,   694,
     214,    -1,   112,   694,   214,    -1,   467,   694,   693,   214,
      -1,   275,   276,   693,   214,    -1,   275,   276,   693,   693,
     693,   214,    -1,   354,   214,    -1,   469,    38,   311,   214,
      -1,   469,    38,   311,   104,   214,    -1,   469,    39,   311,
     214,    -1,   469,   355,   104,   214,    -1,   469,   356,   104,
     214,    -1,   469,   357,   104,   214,    -1,   469,   358,   104,
     214,    -1,   469,   359,   104,   214,    -1,     8,   693,   214,
      -1,    78,   214,    -1,   471,    79,   104,   214,    -1,   471,
     398,   693,   214,    -1,   471,   402,   214,    -1,   471,   401,
     214,    -1,   471,   399,   214,    -1,   471,   400,   214,    -1,
     471,    77,   214,    -1,   471,   343,   214,    -1,   471,   116,
     214,    -1,   403,   214,    -1,   472,   693,   694,   214,    -1,
     109,   214,    -1,   473,    16,   693,   214,    -1,   473,   205,
     693,   214,    -1,   473,   136,   693,   214,    -1,   473,   111,
     693,   214,    -1,   200,   214,    -1,   200,   693,   214,    -1,
     474,   694,   694,   214,    -1,   140,   214,    -1,   140,   693,
     214,    -1,   475,   694,   694,   214,    -1,   180,   214,    -1,
     180,   693,   214,    -1,   476,   693,   694,   214,    -1,   476,
     693,   200,   693,   214,    -1,   476,   693,   140,   693,   214,
      -1,    47,   214,    -1,   477,   479,    -1,   477,   481,    -1,
     477,   482,    -1,   477,   483,    -1,   477,   484,    -1,    44,
     214,    -1,   478,   641,    -1,    43,   693,   214,   480,    -1,
     693,   693,   694,   214,    -1,   480,   693,   693,   694,   214,
      -1,   172,   693,   214,    -1,   481,   485,    -1,   173,   693,
     214,    -1,   482,   485,    -1,   174,   693,   214,    -1,   483,
     486,    -1,   175,   693,   214,    -1,   484,   486,    -1,   693,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   214,
      -1,   693,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   214,    -1,   693,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   214,    -1,
     693,   693,   694,   694,   214,    -1,   176,   214,    -1,   487,
     488,    -1,   693,   694,   694,   694,   694,   694,   214,    -1,
     693,   694,   694,   694,   694,   694,   694,   694,   214,    -1,
     693,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     214,    -1,   693,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   214,    -1,    71,   214,    -1,   489,   693,
     693,   694,   214,    -1,   489,   693,   693,   693,   694,   214,
      -1,   127,   214,    -1,   490,   694,   694,   694,   214,    -1,
     274,   214,    -1,   491,   104,   104,   104,   214,    -1,   491,
     104,   104,   214,    -1,   491,   104,   693,   214,    -1,   179,
     104,   214,    -1,   344,   104,   214,    -1,   288,   214,   599,
      -1,     3,   214,   599,    -1,   346,   214,   599,    -1,   345,
     214,   599,    -1,   448,   214,    -1,   498,   499,   214,    -1,
     450,   451,   693,    -1,   227,   214,    -1,   228,   214,    -1,
     227,   693,   214,    -1,   228,   693,   214,    -1,   227,   169,
     214,    -1,   228,   169,   214,    -1,   227,   169,   693,   214,
      -1,   228,   169,   693,   214,    -1,   500,   501,   214,    -1,
     303,   104,   693,    -1,   303,   693,   693,   104,   693,    -1,
     303,   104,   693,   693,    -1,   303,   104,   693,   132,   693,
      -1,   303,   693,   693,   104,   693,   693,    -1,   303,   693,
     693,   104,   693,   132,   693,    -1,   449,   104,   693,    -1,
     449,   693,   693,   104,   693,    -1,   340,   104,   693,    -1,
     340,   693,   693,   104,   693,    -1,   340,   104,   693,   693,
      -1,   340,   104,   693,   132,   693,    -1,   340,   693,   693,
     104,   693,   693,    -1,   340,   693,   693,   104,   693,   132,
     693,    -1,   501,   219,   693,    -1,   501,   305,    -1,   501,
     694,   694,    -1,   501,    14,    -1,   501,    60,    -1,   501,
      60,   693,    -1,   501,    11,    -1,   501,   282,    -1,   501,
     376,    -1,   501,   283,   311,    -1,   501,   211,    -1,   501,
     211,   693,    -1,   501,   229,    -1,   501,   187,    -1,   455,
     104,    -1,   453,   104,    -1,   454,   104,    -1,   507,    -1,
      87,   214,    -1,    87,   693,   214,    -1,   452,   311,   214,
      -1,   212,   693,   214,    -1,   304,   214,    -1,   304,   693,
     694,   694,   214,    -1,   221,   693,   214,    -1,   322,   694,
     214,    -1,   324,   694,   214,    -1,    93,   214,    -1,   290,
     694,   214,    -1,    19,   214,    -1,    19,   104,   214,    -1,
      19,   104,   693,   214,    -1,    19,   694,   694,   214,    -1,
      19,   694,   694,   693,   214,    -1,   110,   311,   214,    -1,
     168,   693,   214,    -1,   181,   214,    -1,   502,   188,   693,
     214,    -1,   341,   214,    -1,    30,   693,   214,    -1,   405,
     214,    -1,   406,   694,   214,    -1,   185,   214,    -1,    48,
     214,   694,   693,   214,    -1,    48,   214,    -1,   326,   214,
      -1,    76,   214,    -1,   507,   510,    -1,   507,   531,    -1,
     507,   535,    -1,   507,   191,   214,    -1,   507,   312,   693,
     214,    -1,   507,   312,   311,   214,    -1,   507,   312,   693,
     694,   694,   693,   693,   214,    -1,   507,   152,   311,   214,
      -1,   507,   353,   311,   214,    -1,   507,   224,   214,    -1,
     507,    41,   214,    -1,   507,    41,   694,   694,   214,    -1,
      67,   214,   509,   214,    -1,   694,   694,   694,   694,   694,
      -1,   213,   214,    -1,   205,   511,   214,    -1,    16,   511,
     214,    -1,   136,   512,   214,    -1,   694,   694,    -1,   694,
     694,   694,   694,    -1,   694,    -1,   694,   694,   694,    -1,
     694,    -1,   230,   214,    -1,   513,   514,    -1,   513,   515,
      -1,   513,   191,   214,    -1,   205,   694,   694,   693,   214,
      -1,   205,   694,   694,   693,   694,   214,    -1,   136,   694,
     694,   694,   693,   214,    -1,     4,     6,   214,    -1,     4,
     214,     6,   214,    -1,     4,   214,     6,   694,   694,   214,
      -1,     4,   214,     6,   694,   214,    -1,   516,    45,   214,
      -1,   516,   358,   104,   214,    -1,   516,   517,   214,    -1,
      31,    -1,   517,   693,    -1,     5,   214,     6,   694,   694,
     214,    -1,   330,   214,    -1,   329,   214,    -1,   361,   214,
      -1,   147,   214,    -1,   407,   214,    -1,   328,   214,   694,
     214,    -1,   128,   214,   694,   694,   214,    -1,   128,   214,
     694,   214,    -1,   128,   214,    -1,   207,   693,   214,    -1,
     207,   214,    -1,   267,   693,   214,    -1,   267,   214,    -1,
     267,   214,   528,   214,    -1,   693,    -1,   528,   693,    -1,
     148,   214,    -1,   408,   214,    -1,   321,   694,   694,   694,
     214,    -1,   232,   214,   693,   693,   533,    -1,   232,   214,
     693,   693,   693,   533,    -1,   236,   214,   693,   693,   693,
     693,   533,    -1,   214,   534,   533,    -1,   214,    -1,   235,
      -1,   237,   693,    -1,   238,    -1,   239,    -1,   240,    -1,
     241,   694,    -1,   242,    -1,   243,    -1,   243,   694,    -1,
     244,    -1,   266,   694,   694,   214,    -1,   266,   694,   694,
     208,   214,    -1,   206,   214,   600,    -1,   135,   214,   617,
      -1,   162,   214,   694,   694,   694,   214,    -1,   162,   214,
     694,   694,   214,    -1,   163,   214,   693,   694,   694,   214,
      -1,   163,   693,   214,    -1,   538,   539,    -1,   693,   694,
     694,   694,   214,    -1,   164,   214,   694,   694,   694,   214,
      -1,    72,   214,    -1,   541,   635,    -1,   541,   693,   342,
     693,   693,   694,   214,    -1,   541,   693,   342,   693,   314,
     693,   693,   694,   214,    -1,   541,   317,   635,    -1,    63,
     214,    -1,   542,   693,    65,   693,   693,   214,    -1,   409,
     214,   544,    -1,   545,    -1,   544,   545,    -1,   693,   694,
     214,    -1,   693,   214,    -1,   411,   214,   547,    -1,   548,
      -1,   547,   548,    -1,   693,   693,   558,   214,    -1,   320,
     214,    -1,   549,   693,   694,   214,    -1,   549,   317,   693,
     694,   214,    -1,   105,   214,    -1,   105,   693,   214,    -1,
     550,   693,   694,   214,    -1,   550,   317,   693,   694,   214,
      -1,    46,   214,    -1,    46,   693,   214,    -1,   551,   693,
     694,   694,   694,   214,    -1,   265,   214,    -1,   265,   693,
     214,    -1,   552,   693,   694,   694,   694,   214,    -1,   139,
     214,   554,    -1,   555,    -1,   554,   555,    -1,   693,   693,
     694,   214,    -1,   693,   693,   693,   694,   214,    -1,   693,
     693,   693,   693,   694,   214,    -1,    91,   214,   557,    -1,
     556,   557,    -1,   693,   693,   558,   214,    -1,   693,    -1,
     558,   693,    -1,   308,   214,   560,    -1,   561,    -1,   560,
     561,    -1,   693,   693,   694,   214,    -1,   693,   693,   693,
     694,   214,    -1,   693,   693,   693,   693,   694,   214,    -1,
      89,   214,   563,    -1,   562,   563,    -1,   693,   693,   558,
     214,    -1,   142,   214,   565,    -1,   564,   565,    -1,   693,
     693,   558,   214,    -1,   146,   214,   567,    -1,   566,   567,
      -1,   693,   693,   558,   214,    -1,   145,   214,   569,    -1,
     568,   569,    -1,   693,   693,   558,   214,    -1,   693,   694,
     214,    -1,   570,    -1,   571,   570,    -1,    20,   214,   571,
      -1,    21,   214,   571,    -1,    21,   693,   214,   571,    -1,
      15,   694,   214,    -1,   574,   557,    -1,    17,   694,   214,
      -1,   575,   565,    -1,    18,   694,   694,   694,   214,    -1,
     576,   569,    -1,   101,   693,   214,    -1,   101,   693,   693,
     214,    -1,   102,   693,   214,   665,    -1,    22,   214,   694,
     694,   214,   694,   694,   694,   214,    -1,   579,   580,    -1,
     693,   214,    -1,    23,   214,   694,   694,   214,   694,   694,
     694,   214,    -1,   581,   582,    -1,   693,   693,   214,    -1,
     693,   693,   693,   214,    -1,    24,   693,   214,    -1,    25,
     693,   214,    -1,    26,   214,   693,   214,   694,   694,   214,
      -1,   585,   586,    -1,   693,   693,   214,    -1,   693,   693,
     693,   214,    -1,    27,   214,    -1,   587,   588,    -1,   693,
     694,   694,   693,   214,    -1,   693,   694,   694,   693,   694,
     694,   214,    -1,   693,   694,   694,   693,   694,   694,   694,
     694,   694,   214,    -1,   269,   104,   214,    -1,   269,   104,
     693,   214,    -1,   269,   104,   104,   214,    -1,   269,   104,
     104,   693,   214,    -1,   589,   447,   311,   214,    -1,   154,
     214,    -1,   154,   353,   214,    -1,   154,   214,   599,    -1,
     590,   191,   214,   600,    -1,   155,   694,   214,    -1,   155,
     214,    -1,   591,   693,   694,   694,   694,   694,   694,   694,
     214,    -1,   591,   693,   694,   694,   694,   214,    -1,   125,
     214,    -1,    35,   214,    -1,   233,   693,   214,    -1,   592,
     693,   694,   694,   694,   694,   694,   694,   214,    -1,   234,
     693,   214,    -1,   593,   693,   694,   694,   694,   694,   694,
     694,   214,    -1,   160,   214,    -1,   160,   214,   599,    -1,
     594,   191,   214,   600,    -1,   158,   214,   601,    -1,    92,
     214,   601,    -1,   106,   214,    -1,   106,   693,   214,    -1,
     597,   635,    -1,   597,   693,   342,   693,   693,   694,   214,
      -1,   597,   693,   342,   693,   314,   693,   693,   694,   214,
      -1,   597,   693,   342,   693,   693,   694,   208,   214,    -1,
     597,   693,   342,   693,   314,   693,   693,   694,   208,   214,
      -1,   597,   317,   635,    -1,   106,   214,   191,   214,   600,
      -1,   635,    -1,   599,   635,    -1,   693,   342,   693,   693,
     694,   214,    -1,   599,   693,   342,   693,   693,   694,   214,
      -1,   693,   342,   693,   314,   693,   693,   694,   214,    -1,
     599,   693,   342,   693,   314,   693,   693,   694,   214,    -1,
     636,    -1,   600,   636,    -1,   637,    -1,   601,   637,    -1,
     352,   214,    -1,   352,   214,   603,    -1,    52,   693,   214,
     694,   694,   214,    -1,   603,   694,   694,   214,    -1,   603,
      52,   693,   214,   694,   694,   214,    -1,   331,   214,    -1,
     331,   214,   605,    -1,    52,   693,   214,   694,   694,   214,
      -1,   605,   694,   694,   214,    -1,   605,    52,   693,   214,
     694,   694,   214,    -1,   287,   214,    -1,   287,   214,   607,
      -1,    52,   693,   214,   694,   694,   214,    -1,   607,   694,
     694,   214,    -1,   607,    52,   693,   214,   694,   694,   214,
      -1,   178,   214,    -1,   178,   214,   609,    -1,   610,   611,
      -1,   609,   611,    -1,   609,   610,   611,    -1,   693,   214,
      -1,   693,   694,   214,    -1,   693,   694,   361,   693,   214,
      -1,   693,   694,   660,   214,    -1,   693,   694,   361,   693,
     660,   214,    -1,   693,   693,   694,   214,    -1,   149,   214,
      -1,   149,   214,   613,    -1,   614,   615,    -1,   613,   615,
      -1,   613,   614,   615,    -1,   693,    55,   694,   694,   214,
      -1,   693,   694,   214,    -1,   693,   214,    -1,   693,   693,
     694,   694,   214,    -1,   693,   693,   694,   214,    -1,   138,
     214,   617,    -1,   138,   693,   214,   617,    -1,   638,    -1,
     617,   638,    -1,   186,   214,   619,    -1,   618,   619,    -1,
     693,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   214,    -1,   693,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   266,   694,   694,   214,    -1,   693,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   286,   694,   694,   214,    -1,   693,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   214,    -1,   693,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   266,   694,   694,
     214,    -1,   693,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   286,   694,   694,   214,    -1,   693,   694,   694,   694,
     694,   694,   694,   694,   214,    -1,   693,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   214,
      -1,   693,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,    68,   694,   214,    -1,   693,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   693,
     694,   694,   693,   693,   694,   694,   694,   214,    -1,   693,
       7,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     214,    -1,   693,     7,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   214,    -1,   693,     7,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   214,
      -1,   693,     7,   694,   694,   214,    -1,   693,     7,   694,
     694,   694,   214,    -1,   693,    96,   693,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   693,   693,   693,   214,
      -1,   693,   339,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   214,    -1,   693,    61,   214,    -1,   693,    61,
     660,   214,    -1,   693,    61,   185,   694,   214,    -1,   693,
      61,   660,   185,   694,   214,    -1,   693,    61,   185,   694,
     694,   214,    -1,   693,    61,   660,   185,   694,   694,   214,
      -1,   693,    61,   693,   694,   694,   694,   694,   694,   694,
     693,   214,    -1,   693,    61,    95,   694,   694,   694,   694,
     214,    -1,   693,    61,    95,   694,   694,   694,   694,   694,
     214,    -1,   693,    61,    95,   694,   694,   694,   694,   694,
     694,   214,    -1,   693,    61,   660,    95,   694,   694,   694,
     694,   214,    -1,   693,    61,   660,    95,   694,   694,   694,
     694,   694,   214,    -1,   693,    61,   660,    95,   694,   694,
     694,   694,   694,   694,   214,    -1,   693,    61,    95,   694,
     694,   694,   694,   318,   694,   214,    -1,   693,    61,    95,
     694,   694,   694,   694,   694,   318,   694,   214,    -1,   693,
      61,    95,   694,   694,   694,   694,   694,   694,   318,   694,
     214,    -1,   693,    61,   660,    95,   694,   694,   694,   694,
     318,   694,   214,    -1,   693,    61,   660,    95,   694,   694,
     694,   694,   694,   318,   694,   214,    -1,   693,    61,   660,
      95,   694,   694,   694,   694,   694,   694,   318,   694,   214,
      -1,   693,    61,   318,   694,   214,    -1,   693,    61,   660,
     318,   694,   214,    -1,   693,    61,   318,   694,   694,   214,
      -1,   693,    61,   660,   318,   694,   694,   214,    -1,   693,
      61,   318,   694,   694,   694,   214,    -1,   693,    61,   660,
     318,   694,   694,   694,   214,    -1,   693,   318,   694,   214,
      -1,   327,   214,   633,    -1,   620,   633,    -1,   390,   693,
     214,    -1,   390,   693,   281,   214,    -1,   390,   693,   316,
     694,   214,    -1,   390,   693,   316,   694,   281,   214,    -1,
     621,   693,   693,   634,   214,    -1,   391,   693,   693,   214,
      -1,   391,   693,   693,   395,   214,    -1,   391,   693,   693,
     396,   214,    -1,   391,   693,   693,   394,   694,   214,    -1,
     391,   693,   693,   394,   694,   694,   214,    -1,   391,   693,
     693,   395,   394,   694,   214,    -1,   391,   693,   693,   395,
     394,   694,   694,   214,    -1,   391,   693,   693,   396,   394,
     694,   214,    -1,   391,   693,   693,   396,   394,   694,   694,
     214,    -1,   397,   693,   693,   214,    -1,   397,   693,   214,
      -1,   334,   214,    -1,   624,   693,   693,   693,   214,    -1,
     624,   693,   693,   693,   693,   214,    -1,   624,   693,   693,
     693,   693,   694,   214,    -1,   624,   693,   693,   693,   693,
     694,   694,   214,    -1,   624,   693,   693,   693,   693,   694,
     694,   693,   694,   214,    -1,   624,   693,   693,   693,   693,
     694,   694,   693,   694,   694,   214,    -1,   624,   693,   693,
     693,   660,   214,    -1,   624,   693,   693,   693,   693,   694,
     694,   660,   214,    -1,   117,   214,    -1,   625,   693,   693,
     693,   214,    -1,   625,   693,   693,   693,   694,   694,   214,
      -1,   410,   214,   627,    -1,   626,   627,    -1,   693,   693,
     558,   214,    -1,   412,   214,    -1,   628,   693,   693,   693,
     214,    -1,   628,   693,   693,   693,   694,   694,   214,    -1,
      57,   214,    -1,   629,   693,   693,   693,   214,    -1,   629,
     693,   693,   693,   693,   214,    -1,   629,   693,   693,   693,
     693,   694,   214,    -1,   629,   693,   693,   693,   693,   694,
     694,   214,    -1,   629,   693,   693,   693,   693,   694,   694,
     693,   694,   214,    -1,   629,   693,   693,   693,   693,   694,
     694,   693,   694,   694,   214,    -1,   629,   693,   693,   693,
     693,   694,   694,   693,   694,   694,   694,   694,   214,    -1,
     629,   693,   693,   693,   693,   694,   694,   693,   694,   694,
     694,   694,   694,   214,    -1,   629,   693,   693,   693,   660,
     214,    -1,   629,   693,   693,   693,   693,   694,   694,   660,
     214,    -1,    29,   693,   214,    -1,   122,   693,   214,    -1,
     392,   694,   214,    -1,   393,   693,   214,    -1,   219,   214,
     632,    -1,   631,   632,    -1,   693,   694,   694,   694,   214,
      -1,   693,   694,   694,   214,    -1,   693,   694,   214,    -1,
     693,   694,   694,   694,   693,   693,   214,    -1,   693,   694,
     694,   694,   693,   214,    -1,   693,   693,   634,   214,    -1,
     693,    -1,   634,   693,    -1,   693,   693,   694,   214,    -1,
     693,   693,   214,    -1,   693,   693,   694,   208,   214,    -1,
     693,   693,   208,   214,    -1,   693,   694,   214,    -1,   693,
     694,   214,    -1,   693,   693,   694,   694,   214,    -1,   693,
     693,   694,   214,    -1,    64,   214,    -1,   639,   641,    -1,
      88,   214,    -1,   640,   641,    -1,   693,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   214,    -1,   693,   338,
     693,   214,    -1,   225,   214,    -1,   642,   643,    -1,   693,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   214,
      -1,   693,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   214,    -1,   693,   124,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   214,    -1,   693,
     124,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   214,    -1,    34,   214,    -1,   644,   693,
     693,   694,   694,   694,   214,    -1,    10,   214,    -1,   645,
     693,   693,   214,    -1,   645,   693,   693,   337,   694,   214,
      -1,   645,   693,   693,   337,   336,   694,   214,    -1,   645,
     693,   693,   693,   693,   214,    -1,   645,   693,   693,   693,
     693,   337,   694,   214,    -1,   645,   693,   693,   693,   335,
     694,   214,    -1,   645,   693,   693,   153,   214,    -1,   645,
     693,   693,   693,   214,    -1,   645,   693,   693,   693,   693,
     693,   214,    -1,   645,   693,   693,   693,   693,   335,   694,
     214,    -1,    86,   214,    -1,   646,   693,   694,   214,    -1,
     646,   336,   214,    -1,   441,   214,    -1,   647,   694,   214,
      -1,   442,   214,    -1,   648,   693,   693,   214,    -1,   648,
     694,   214,    -1,   440,   214,    -1,   649,   693,   693,   214,
      -1,   649,   693,   693,   693,   214,    -1,   248,   214,    -1,
     248,   693,   214,    -1,   650,   693,   694,   214,    -1,   650,
     693,   693,   694,   214,    -1,   650,   317,   693,   694,   214,
      -1,   650,   693,   694,   311,   214,    -1,   650,   693,   693,
     694,   311,   214,    -1,   650,   317,   693,   694,   311,   214,
      -1,   184,   214,    -1,   184,   693,   214,    -1,   184,   693,
     693,   214,    -1,   247,   214,    -1,   652,   693,   694,   214,
      -1,   652,   693,   342,   693,   694,   214,    -1,   652,   693,
     694,   694,   694,   214,    -1,   652,   693,   342,   693,   694,
     694,   694,   214,    -1,   302,   214,    -1,   653,   656,    -1,
     653,   655,    -1,   653,    62,   654,   214,    -1,   653,   264,
     214,    -1,   653,   264,   694,   694,   214,    -1,   693,    -1,
     654,   693,    -1,   159,   214,    -1,   655,   245,   693,   214,
      -1,   655,   188,   693,   214,    -1,   655,   325,   694,   214,
      -1,   655,   189,   693,   214,    -1,   655,   313,   693,   214,
      -1,    73,   214,    -1,    73,   693,   214,    -1,   289,   214,
      -1,   289,   251,   214,    -1,   289,   348,   214,    -1,   159,
     693,   214,    -1,   159,   693,   694,   214,    -1,   159,   693,
     694,   693,   214,    -1,   159,   693,   694,   693,   693,   214,
      -1,   159,   693,   694,   693,   693,   693,   214,    -1,    98,
     693,   694,   214,    -1,    98,   693,   214,    -1,    98,    75,
     214,    -1,    98,   370,   214,    -1,    98,   693,    99,   214,
      -1,    98,   693,   694,   693,   214,    -1,    98,   214,    -1,
      33,   214,   289,   214,    -1,   300,   693,   214,    -1,   301,
     693,   214,    -1,   291,   694,   214,    -1,   293,   693,   214,
      -1,   294,   693,   214,    -1,   292,   693,   214,    -1,   295,
     694,   214,    -1,   296,   693,   214,    -1,   298,   311,   214,
      -1,   297,   693,   214,    -1,   299,   693,   214,    -1,   203,
     693,   693,   214,    -1,   204,   693,   694,   214,    -1,   133,
     694,   214,    -1,   134,   311,   214,    -1,   656,   188,   693,
     214,    -1,    82,   693,   693,   214,    -1,    81,   693,   694,
     214,    -1,   656,   245,   100,   214,    -1,   656,   245,   184,
     214,    -1,   656,   245,   693,   214,    -1,   252,   693,   214,
      -1,   252,   253,   214,    -1,   323,   694,   214,    -1,   323,
     694,   694,   214,    -1,   309,   694,   214,    -1,   309,   694,
     694,   214,    -1,   261,   694,   694,   214,    -1,   263,   693,
     693,   214,    -1,   262,   694,   694,   214,    -1,   656,   189,
     693,   214,    -1,   218,   214,    -1,   250,   693,   214,    -1,
     284,   693,   214,    -1,   284,   285,   214,    -1,   196,   693,
     214,    -1,   196,   285,   214,    -1,   118,   693,   214,    -1,
     118,   285,   214,    -1,   197,   214,    -1,   119,   214,    -1,
     121,   693,   214,    -1,   120,   214,    -1,    56,   694,   214,
      -1,   351,   214,    -1,    50,    51,   214,    -1,    50,    51,
     693,   214,    -1,    50,    13,   214,    -1,    12,    13,   214,
      -1,    12,    13,   268,   214,    -1,    12,   365,   693,   214,
      -1,    12,   365,   366,   693,   214,    -1,    12,   365,   693,
     371,   214,    -1,    12,   365,   366,   693,   371,   214,    -1,
     367,   694,   214,    -1,   367,   694,   694,   214,    -1,   126,
     694,   214,    -1,   130,   694,   214,    -1,    53,   694,   214,
      -1,   277,   311,   214,    -1,   347,   693,   214,    -1,   249,
     214,    -1,   249,   104,   214,    -1,   182,   289,   214,    -1,
      42,   289,   214,    -1,    28,   289,   214,    -1,    54,   289,
     214,    -1,   182,   289,   251,   214,    -1,   182,   289,   310,
     214,    -1,    42,   289,   251,   214,    -1,    28,   289,   251,
     214,    -1,    42,   289,   310,   214,    -1,    54,   289,   251,
     214,    -1,    54,   289,   310,   214,    -1,   350,   693,   214,
      -1,   129,   214,    -1,   404,   214,    -1,   404,   311,   214,
      -1,   254,   693,   214,    -1,   217,   214,    -1,   217,   693,
     214,    -1,   137,   693,   694,   214,    -1,   137,   693,   694,
     693,   214,    -1,   137,   693,   694,   693,   694,   214,    -1,
     137,   693,   694,   693,   694,   693,   214,    -1,   137,   214,
      -1,   223,   693,   214,    -1,   223,   693,   694,   214,    -1,
     332,   694,   214,    -1,   307,   693,   214,    -1,   170,   693,
     214,    -1,   170,   165,   214,    -1,   157,   214,    -1,   306,
     214,    -1,   333,   693,   214,    -1,   333,   693,   694,   214,
      -1,   333,   693,   694,   694,   214,    -1,   364,   159,   214,
      -1,   364,   159,   150,   214,    -1,   194,   693,   214,    -1,
     194,   195,   214,    -1,   192,   693,   214,    -1,   192,   193,
     214,    -1,   192,   693,   199,   693,   214,    -1,   192,   193,
     199,   693,   214,    -1,   192,   693,   693,   214,    -1,   192,
     193,   198,   214,    -1,   192,   693,   693,   199,   693,   214,
      -1,   192,   193,   198,   199,   693,   214,    -1,   201,   693,
     214,    -1,   256,   694,   214,    -1,   183,   693,   694,   214,
      -1,    37,   311,   214,    -1,    80,   311,   214,    -1,    58,
     311,   214,    -1,    59,   311,   214,    -1,   202,   693,   214,
      -1,    97,   214,   658,    -1,   694,   214,    -1,    83,   660,
     214,    -1,    83,   214,   660,   214,    -1,    73,    -1,    73,
     694,    -1,    73,   694,   694,    -1,    73,   694,   694,   694,
      -1,    73,   694,   694,   694,   694,    -1,    84,    -1,    85,
     694,    -1,    84,    85,   694,    -1,    32,   694,    -1,   660,
     151,   693,    -1,   660,   151,   693,   694,    -1,   141,   214,
      -1,   108,   214,    -1,    36,   693,   214,    -1,    36,   693,
     694,   214,    -1,    36,   693,   694,   694,   214,    -1,   259,
     693,   214,   662,    -1,   260,   693,   214,   662,    -1,   161,
     693,   214,   662,    -1,   661,   535,    -1,   661,   468,    -1,
     661,   465,    -1,   663,    -1,   662,   663,    -1,   694,   694,
     694,   214,    -1,   171,   214,   694,   694,   694,   214,    -1,
     664,   694,   694,   694,   214,    -1,   666,    -1,   665,   666,
      -1,   694,   694,   694,   214,    -1,   143,   214,   668,    -1,
     694,   214,    -1,   144,   214,   670,    -1,   694,   214,    -1,
      70,   214,    -1,   215,   214,    -1,   672,     9,   214,    -1,
     672,   216,   214,    -1,   672,   376,   214,    -1,   672,   188,
     693,   214,    -1,   672,   222,   694,   214,    -1,   672,   222,
     694,   694,   214,    -1,   672,   222,   694,   694,   694,   694,
     214,    -1,   672,    74,   694,   694,   214,    -1,   672,    74,
     694,   694,   693,   693,   214,    -1,   672,   103,   693,   214,
      -1,   672,   103,   693,   693,   214,    -1,   672,   348,   214,
      -1,   672,   177,   694,   214,    -1,   672,   183,   693,   214,
      -1,   672,   183,   693,   693,   214,    -1,   672,   183,   693,
     693,   694,   214,    -1,   672,   183,   693,   693,   694,   694,
     214,    -1,   672,   183,   693,   693,   694,   694,   693,   214,
      -1,   672,   123,   214,    -1,   672,   123,   694,   214,    -1,
     672,    85,   693,   694,   694,   214,    -1,   672,   673,    -1,
     270,   693,   214,    -1,   270,   693,   693,   214,    -1,   273,
     214,    -1,    49,   214,    -1,    49,   214,   104,   214,   693,
     214,   104,   214,   104,   214,    -1,    49,   214,   104,   214,
     693,   214,   104,   214,   104,   214,   104,   214,    -1,   226,
     104,   214,    -1,   676,   214,   104,   214,    -1,   360,   214,
      -1,   360,   361,   693,   214,    -1,   677,   693,   693,   694,
     694,   694,   214,    -1,   677,   693,   693,   694,   694,   694,
     363,   694,   214,    -1,   677,   693,   693,   694,   694,   694,
     361,   693,   214,    -1,   677,   693,   693,   694,   694,   694,
     660,   214,    -1,   677,   693,   693,   694,   694,   694,   363,
     694,   361,   693,   214,    -1,   677,   693,   693,   694,   694,
     694,   363,   694,   361,   693,   660,   214,    -1,   372,   214,
      -1,   678,   693,   374,   694,   694,   694,   694,   694,   214,
      -1,   678,   693,   374,   694,   694,   694,   694,   694,   694,
     214,    -1,   678,   693,   375,   694,   694,   694,   694,   694,
     214,    -1,   678,   693,   375,   694,   694,   694,   694,   694,
     694,   214,    -1,   678,   693,   388,   694,   694,   694,   694,
     694,   214,    -1,   678,   693,   388,   694,   694,   694,   694,
     694,   694,   214,    -1,   678,   693,   376,   694,   694,   694,
     214,    -1,   678,   693,   377,   694,   694,   694,   214,    -1,
     678,   693,   387,   694,   694,   694,   214,    -1,   678,   693,
     378,   694,   694,   694,   694,   214,    -1,   678,   693,   389,
     694,   694,   694,   694,   214,    -1,   678,   693,   381,   694,
     694,   694,   214,    -1,   678,   693,   382,   694,   694,   694,
     214,    -1,   678,   693,   382,   694,   694,   694,   694,   214,
      -1,   678,   693,   386,   694,   694,   694,   694,   214,    -1,
     678,   693,   386,   694,   694,   694,   694,   694,   214,    -1,
     678,   693,   383,   694,   694,   694,   694,   694,   694,   214,
      -1,   678,   693,   384,   694,   694,   694,   694,   694,   694,
     214,    -1,   678,   693,   384,   694,   694,   694,   694,   694,
     694,   694,   214,    -1,   678,   693,   384,   694,   694,   694,
     694,   694,   694,   694,   694,   214,    -1,   678,   693,   380,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     214,    -1,   678,   693,   380,   694,   694,   694,   214,    -1,
     678,   693,   380,   694,   694,   694,   694,   214,    -1,   678,
     693,   380,   694,   694,   694,   694,   694,   214,    -1,   678,
     693,   380,   694,   694,   694,   694,   694,   694,   214,    -1,
     678,   693,   380,   694,   694,   694,   694,   694,   694,   694,
     214,    -1,   678,   379,   104,   104,   214,    -1,   678,   693,
     104,   680,   214,    -1,   373,   214,    -1,   679,   693,   693,
     214,    -1,   679,   693,   693,   693,   214,    -1,    -1,   680,
     694,    -1,    -1,   681,   104,    -1,   271,   214,   272,   214,
      -1,   271,   214,   272,   214,   272,   214,    -1,   271,   214,
     272,   214,   272,   214,   272,   214,    -1,   418,   214,    -1,
     683,   684,   214,    -1,   413,   681,    -1,   423,   693,    -1,
     427,   693,    -1,   429,   680,    -1,   400,   693,    -1,   430,
     681,    -1,    40,   693,    -1,   690,    -1,   444,   214,    -1,
     685,   689,   214,    -1,   445,   214,    -1,   686,   689,   214,
      -1,   421,   214,    -1,   687,   689,   214,    -1,   687,   690,
      -1,   422,   214,    -1,   688,   689,   214,    -1,   414,   104,
      -1,   415,   104,    -1,   415,   104,   104,    -1,   415,   104,
     104,   104,    -1,   425,   694,    -1,   400,   693,    -1,   416,
     693,    -1,   423,   693,    -1,   423,   693,   693,    -1,   423,
     693,   693,   693,    -1,   426,   311,    -1,   447,   311,    -1,
     209,   311,    -1,   336,   311,    -1,   437,   104,    -1,   437,
     104,   104,    -1,   437,   104,   104,   104,    -1,   438,   311,
      -1,   443,   104,    -1,   443,   104,   693,    -1,   443,   104,
     693,   693,    -1,   435,   311,    -1,   436,   311,    -1,   439,
     311,    -1,   446,   311,    -1,   434,   214,    -1,   690,   509,
     214,    -1,   419,   214,    -1,   691,   692,   214,    -1,   420,
     104,    -1,   420,   104,   693,    -1,   156,    -1,   210,    -1,
     156,    -1,    69,    -1,   167,    -1,    94,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   164,   164,   170,   171,   174,   175,   177,   179,   180,
     181,   182,   183,   184,   185,   186,   188,   189,   190,   191,
     192,   193,   194,   195,   197,   199,   200,   201,   202,   203,
     204,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     226,   228,   229,   230,   231,   232,   233,   234,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,   246,
     247,   248,   249,   250,   251,   252,   253,   254,   256,   258,
     260,   262,   264,   266,   267,   268,   269,   270,   271,   272,
     273,   274,   275,   276,   277,   278,   279,   280,   281,   282,
     284,   285,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   298,   300,   301,   302,   303,   304,   305,   306,
     308,   309,   310,   311,   312,   314,   315,   316,   317,   319,
     321,   323,   325,   327,   329,   331,   333,   334,   335,   336,
     337,   338,   339,   340,   341,   342,   343,   346,   353,   359,
     360,   365,   377,   383,   384,   388,   392,   394,   396,   398,
     418,   444,   448,   450,   452,   454,   474,   500,   501,   502,
     503,   504,   507,   511,   513,   518,   521,   525,   564,   613,
     614,   616,   628,   630,   632,   634,   636,   638,   642,   650,
     652,   654,   656,   658,   660,   662,   664,   666,   668,   672,
     674,   678,   680,   682,   684,   686,   689,   691,   693,   697,
     699,   701,   705,   707,   709,   711,   713,   717,   718,   719,
     720,   721,   722,   725,   726,   730,   734,   736,   740,   742,
     746,   748,   752,   754,   758,   760,   764,   770,   777,   786,
     790,   791,   795,   802,   808,   813,   819,   820,   822,   826,
     827,   831,   832,   836,   839,   844,   849,   853,   858,   864,
     871,   878,   883,   887,   892,   894,   896,   898,   900,   902,
     904,   906,   908,   912,   914,   916,   918,   921,   923,   925,
     927,   929,   931,   933,   935,   937,   939,   942,   944,   946,
     948,   950,   952,   954,   956,   958,   960,   962,   964,   966,
     968,   970,   972,   974,   978,   979,   982,   985,   987,   989,
     991,   993,   995,   997,   999,  1001,  1003,  1005,  1008,  1012,
    1015,  1018,  1020,  1022,  1025,  1027,  1029,  1033,  1035,  1039,
    1043,  1046,  1050,  1054,  1056,  1057,  1058,  1059,  1061,  1063,
    1065,  1072,  1074,  1076,  1078,  1080,  1086,  1091,  1106,  1108,
    1109,  1111,  1114,  1116,  1118,  1120,  1126,  1137,  1140,  1141,
    1142,  1146,  1148,  1152,  1162,  1165,  1168,  1176,  1181,  1183,
    1185,  1189,  1191,  1195,  1199,  1203,  1207,  1211,  1215,  1219,
    1227,  1230,  1233,  1238,  1240,  1244,  1246,  1248,  1251,  1259,
    1269,  1273,  1277,  1281,  1287,  1292,  1301,  1302,  1305,  1307,
    1309,  1311,  1313,  1315,  1317,  1319,  1321,  1323,  1327,  1329,
    1332,  1337,  1341,  1347,  1355,  1361,  1366,  1369,  1375,  1383,
    1385,  1387,  1389,  1391,  1405,  1407,  1417,  1421,  1423,  1427,
    1429,  1433,  1474,  1478,  1484,  1490,  1492,  1495,  1504,  1506,
    1508,  1511,  1522,  1524,  1526,  1531,  1533,  1535,  1540,  1543,
    1544,  1547,  1549,  1551,  1555,  1556,  1559,  1565,  1567,  1572,
    1575,  1576,  1579,  1582,  1585,  1590,  1591,  1594,  1599,  1600,
    1603,  1607,  1608,  1611,  1617,  1618,  1621,  1625,  1629,  1631,
    1635,  1639,  1641,  1645,  1647,  1650,  1652,  1655,  1659,  1662,
    1664,  1668,  1674,  1679,  1683,  1687,  1692,  1696,  1698,  1702,
    1710,  1717,  1722,  1725,  1727,  1731,  1733,  1737,  1739,  1741,
    1745,  1748,  1752,  1756,  1761,  1765,  1767,  1769,  1772,  1778,
    1780,  1782,  1793,  1804,  1807,  1812,  1814,  1827,  1829,  1841,
    1843,  1846,  1852,  1857,  1863,  1865,  1867,  1869,  1871,  1873,
    1875,  1877,  1886,  1891,  1893,  1895,  1897,  1899,  1901,  1905,
    1907,  1911,  1913,  1917,  1918,  1921,  1923,  1925,  1929,  1930,
    1933,  1935,  1937,  1941,  1942,  1945,  1947,  1949,  1953,  1954,
    1957,  1961,  1963,  1969,  1972,  1975,  1979,  1984,  1992,  2004,
    2005,  2008,  2010,  2012,  2016,  2018,  2020,  2024,  2032,  2042,
    2044,  2049,  2051,  2055,  2056,  2059,  2067,  2076,  2085,  2094,
    2104,  2114,  2121,  2129,  2137,  2151,  2168,  2186,  2204,  2213,
    2221,  2240,  2248,  2254,  2264,  2271,  2282,  2290,  2302,  2316,
    2327,  2339,  2352,  2367,  2383,  2400,  2412,  2425,  2439,  2455,
    2472,  2490,  2497,  2508,  2516,  2528,  2537,  2550,  2558,  2559,
    2562,  2568,  2574,  2581,  2589,  2600,  2602,  2604,  2608,  2610,
    2612,  2614,  2616,  2620,  2626,  2628,  2636,  2638,  2644,  2651,
    2658,  2665,  2673,  2682,  2689,  2699,  2701,  2707,  2715,  2718,
    2721,  2727,  2729,  2735,  2745,  2747,  2754,  2761,  2768,  2775,
    2783,  2792,  2801,  2810,  2818,  2828,  2830,  2832,  2834,  2837,
    2839,  2843,  2845,  2847,  2849,  2853,  2858,  2863,  2865,  2870,
    2872,  2874,  2876,  2880,  2884,  2888,  2890,  2894,  2895,  2899,
    2900,  2904,  2909,  2914,  2915,  2919,  2926,  2933,  2940,  2949,
    2950,  2957,  2959,  2961,  2965,  2970,  2972,  2976,  2978,  2983,
    2988,  2993,  3000,  3002,  3004,  3009,  3011,  3015,  3016,  3019,
    3023,  3024,  3027,  3033,  3035,  3037,  3041,  3047,  3052,  3056,
    3062,  3069,  3074,  3079,  3087,  3089,  3091,  3096,  3099,  3107,
    3109,  3110,  3111,  3112,  3129,  3150,  3152,  3156,  3159,  3161,
    3163,  3165,  3167,  3171,  3173,  3175,  3177,  3181,  3184,  3186,
    3188,  3190,  3192,  3194,  3199,  3202,  3206,  3211,  3216,  3221,
    3223,  3228,  3230,  3232,  3234,  3236,  3243,  3245,  3252,  3254,
    3256,  3258,  3260,  3262,  3264,  3266,  3268,  3270,  3272,  3278,
    3280,  3282,  3289,  3296,  3303,  3305,  3308,  3310,  3313,  3316,
    3319,  3322,  3324,  3326,  3346,  3349,  3351,  3354,  3356,  3359,
    3361,  3363,  3365,  3367,  3369,  3371,  3373,  3375,  3379,  3388,
    3409,  3427,  3432,  3438,  3444,  3451,  3453,  3456,  3458,  3460,
    3462,  3464,  3466,  3468,  3471,  3473,  3475,  3477,  3479,  3483,
    3486,  3490,  3494,  3497,  3501,  3505,  3531,  3533,  3535,  3537,
    3543,  3548,  3553,  3566,  3580,  3597,  3613,  3623,  3625,  3628,
    3630,  3632,  3634,  3636,  3638,  3640,  3645,  3651,  3658,  3660,
    3663,  3665,  3667,  3669,  3671,  3674,  3677,  3680,  3683,  3687,
    3691,  3693,  3695,  3698,  3700,  3702,  3704,  3706,  3710,  3713,
    3724,  3730,  3738,  3745,  3752,  3760,  3769,  3779,  3785,  3791,
    3797,  3803,  3806,  3811,  3816,  3826,  3831,  3837,  3844,  3849,
    3854,  3859,  3860,  3861,  3864,  3865,  3868,  3872,  3875,  3880,
    3881,  3884,  3888,  3891,  3902,  3905,  3916,  3922,  3935,  3940,
    3947,  3949,  3951,  3953,  3956,  3961,  3964,  3969,  3972,  3975,
    3977,  3979,  3981,  3984,  3990,  3996,  4003,  4005,  4008,  4012,
    4015,  4020,  4057,  4061,  4062,  4064,  4074,  4080,  4088,  4089,
    4093,  4096,  4099,  4101,  4104,  4106,  4110,  4111,  4116,  4121,
    4126,  4131,  4136,  4141,  4146,  4151,  4156,  4161,  4166,  4172,
    4178,  4184,  4190,  4196,  4202,  4208,  4214,  4220,  4225,  4230,
    4235,  4240,  4245,  4250,  4254,  4260,  4261,  4263,  4270,  4271,
    4282,  4283,  4294,  4298,  4301,  4308,  4312,  4316,  4327,  4329,
    4331,  4333,  4335,  4337,  4339,  4343,  4347,  4351,  4355,  4359,
    4363,  4364,  4368,  4372,  4376,  4378,  4380,  4383,  4387,  4389,
    4391,  4393,  4395,  4398,  4402,  4404,  4406,  4408,  4410,  4412,
    4415,  4419,  4421,  4423,  4426,  4430,  4432,  4434,  4436,  4441,
    4442,  4447,  4451,  4455,  4458,  4465,  4467,  4472,  4474,  4476,
    4478
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ACTUATORS", "AERO", "AEROH", "AEROTYPE",
  "AMAT", "ANALYSIS", "ARCLENGTH", "ATTRIBUTES", "ANGULAROUTTYPE",
  "AUGMENT", "AUGMENTTYPE", "AVERAGED", "ATDARB", "ACOU", "ATDDNB",
  "ATDROB", "ARPACK", "ATDDIR", "ATDNEU", "AXIHDIR", "AXIHNEU",
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
  "ELHSOMMERFELD", "ETEMP", "EXPLICIT", "EPSILON",
  "ELEMENTARYFUNCTIONTYPE", "FABMAT", "FACOUSTICS", "FETI", "FETI2TYPE",
  "FETIPREC", "FFP", "FFPDIR", "FITALG", "FNAME", "FLUX", "FORCE",
  "FRONTAL", "FETIH", "FIELDWEIGHTLIST", "FILTEREIG", "FLUID", "FREQSWEEP",
  "FREQSWEEP1", "FREQSWEEP2", "FREQSWEEPA", "FSGL", "FSINTERFACE",
  "FSISCALING", "FSIELEMENT", "NOLOCALFSISPLITING", "FSICORNER",
  "FFIDEBUG", "FAILSAFE", "FRAMETYPE", "GEPS", "GLOBALTOL", "GRAVITY",
  "GRBM", "GTGSOLVER", "GLOBALCRBMTOL", "GROUP", "GROUPTYPE",
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
  "MPROJECT", "MAXIMUM", "NDTYPE", "NEIGPA", "NEWMARK", "NewLine", "NL",
  "NLMAT", "NLPREC", "NOCOARSE", "NODETOKEN", "NONINPC", "NSBSPV", "NLTOL",
  "NUMCGM", "NOSECONDARY", "NFRAMES", "OPTIMIZATION", "OUTPUT", "OUTPUT6",
  "OUTPUTFRAME", "QSTATIC", "QLOAD", "PITA", "PITADISP6", "PITAVEL6",
  "NOFORCE", "MDPITA", "GLOBALBASES", "LOCALBASES", "TIMEREVERSIBLE",
  "REMOTECOARSE", "ORTHOPROJTOL", "READINITSEED", "JUMPCVG", "JUMPOUTPUT",
  "PRECNO", "PRECONDITIONER", "PRELOAD", "PRESSURE", "PRINTMATLAB", "PROJ",
  "PIVOT", "PRECTYPE", "PRECTYPEID", "PICKANYCORNER", "PADEPIVOT",
  "PROPORTIONING", "PLOAD", "PADEPOLES", "POINTSOURCE", "PLANEWAVE",
  "PTOL", "PLANTOL", "PMAXIT", "PIECEWISE", "RADIATION", "RAYDAMP",
  "RBMFILTER", "RBMSET", "READMODE", "REBUILD", "RENUM", "RENUMBERID",
  "REORTHO", "RESTART", "RECONS", "RECONSALG", "REBUILDCCT", "RANDOM",
  "RPROP", "RNORM", "REVERSENORMALS", "ROTVECOUTTYPE", "RESCALING",
  "SCALING", "SCALINGTYPE", "STRDAMP", "SDETAFT", "SENSORS", "SOLVERTYPE",
  "SHIFT", "SPOOLESTAU", "SPOOLESSEED", "SPOOLESMAXSIZE",
  "SPOOLESMAXDOMAINSIZE", "SPOOLESMAXZEROS", "SPOOLESMSGLVL",
  "SPOOLESSCALE", "SPOOLESPIVOT", "SPOOLESRENUM", "SPARSEMAXSUP",
  "SPARSEDEFBLK", "STATS", "STRESSID", "SUBSPACE", "SURFACE",
  "SAVEMEMCOARSE", "SPACEDIMENSION", "SCATTERER", "STAGTOL", "SCALED",
  "SWITCH", "STABLE", "SUBTYPE", "STEP", "SOWER", "SHELLTHICKNESS", "SURF",
  "SPRINGMAT", "TANGENT", "TEMP", "TIME", "TOLEIG", "TOLFETI", "TOLJAC",
  "TOLPCG", "TOPFILE", "TOPOLOGY", "TRBM", "THERMOE", "THERMOH", "TETT",
  "TOLCGM", "TURKEL", "TIEDSURFACES", "THETA", "REDFOL", "HRC",
  "THIRDNODE", "THERMMAT", "TDENFORC", "TESTULRICH", "THRU", "TRIVIAL",
  "USE", "USERDEFINEDISP", "USERDEFINEFORCE", "UPROJ", "UNSYMMETRIC",
  "USING", "VERSION", "WETCORNERS", "YMTT", "ZERO", "BINARY", "GEOMETRY",
  "DECOMPOSITION", "GLOBAL", "MATCHER", "CPUMAP", "NODALCONTACT", "MODE",
  "FRIC", "GAP", "OUTERLOOP", "EDGEWS", "WAVETYPE", "ORTHOTOL", "IMPE",
  "FREQ", "DPH", "WAVEMETHOD", "MATSPEC", "MATUSAGE", "BILINEARPLASTIC",
  "FINITESTRAINPLASTIC", "LINEARELASTIC", "STVENANTKIRCHHOFF",
  "LINPLSTRESS", "READ", "OPTCTV", "ISOTROPICLINEARELASTIC", "NEOHOOKEAN",
  "ISOTROPICLINEARELASTICJ2PLASTIC",
  "ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS", "HYPERELASTIC",
  "MOONEYRIVLIN", "HENCKY", "LOGSTRAINPLASTIC", "SVKPLSTRESS",
  "SURFACETOPOLOGY", "MORTARTIED", "MORTARSCALING",
  "MORTARINTEGRATIONRULE", "SEARCHTOL", "STDMORTAR", "DUALMORTAR",
  "WETINTERFACE", "NSUBS", "EXITAFTERDEC", "SKIP", "OUTPUTMEMORY",
  "OUTPUTWEIGHT", "WEIGHTLIST", "GMRESRESIDUAL", "SLOSH", "SLGRAV",
  "SLZEM", "SLZEMFILTER", "PDIR", "HEFSB", "HEFRS", "HEINTERFACE",
  "SNAPFI", "PODROB", "TRNVCT", "OFFSET", "ORTHOG", "SVDTOKEN",
  "CONVERSIONTOKEN", "CONVFI", "SAMPLING", "SNAPSHOTPROJECT", "PODSIZEMAX",
  "REFSUBSTRACT", "TOLER", "OUTOFCORE", "NORMALIZETOKEN", "FNUMBER",
  "SNAPWEIGHT", "ROBFI", "STAVCT", "VELVCT", "ACCVCT", "CONWEPCFG",
  "PSEUDOGNAT", "PSEUDOGNATELEM", "VECTORNORM", "LOCALTOLERANCE",
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
  "DEMInfo", "NLInfo", "NewtonInfo", "OrthoInfo", "Control",
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
     705,   706,   707,   708,   709,   710
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   456,   457,   458,   458,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   459,   459,   459,
     459,   459,   459,   459,   459,   459,   459,   460,   461,   462,
     462,   462,   462,   463,   463,   464,   464,   464,   464,   464,
     464,   464,   464,   464,   464,   464,   464,   464,   464,   464,
     464,   464,   465,   466,   466,   467,   467,   468,   468,   469,
     469,   469,   469,   469,   469,   469,   469,   469,   470,   471,
     471,   471,   471,   471,   471,   471,   471,   471,   471,   472,
     472,   473,   473,   473,   473,   473,   474,   474,   474,   475,
     475,   475,   476,   476,   476,   476,   476,   477,   477,   477,
     477,   477,   477,   478,   478,   479,   480,   480,   481,   481,
     482,   482,   483,   483,   484,   484,   485,   485,   485,   486,
     487,   487,   488,   488,   488,   488,   489,   489,   489,   490,
     490,   491,   491,   491,   491,   492,   493,   494,   495,   496,
     497,   498,   498,   499,   500,   500,   500,   500,   500,   500,
     500,   500,   500,   501,   501,   501,   501,   501,   501,   501,
     501,   501,   501,   501,   501,   501,   501,   501,   501,   501,
     501,   501,   501,   501,   501,   501,   501,   501,   501,   501,
     501,   501,   501,   501,   502,   502,   502,   502,   502,   502,
     502,   502,   502,   502,   502,   502,   502,   502,   502,   502,
     502,   502,   502,   502,   502,   502,   502,   503,   503,   504,
     505,   505,   506,   507,   507,   507,   507,   507,   507,   507,
     507,   507,   507,   507,   507,   507,   508,   509,   510,   510,
     510,   510,   511,   511,   511,   511,   512,   513,   513,   513,
     513,   514,   514,   515,   516,   516,   516,   516,   516,   516,
     516,   517,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   525,   525,   526,   526,   527,   527,   527,   528,   528,
     529,   530,   531,   532,   532,   532,   533,   533,   534,   534,
     534,   534,   534,   534,   534,   534,   534,   534,   535,   535,
     535,   536,   537,   537,   538,   538,   538,   539,   540,   541,
     541,   541,   541,   541,   542,   542,   543,   544,   544,   545,
     545,   546,   547,   547,   548,   549,   549,   549,   550,   550,
     550,   550,   551,   551,   551,   552,   552,   552,   553,   554,
     554,   555,   555,   555,   556,   556,   557,   558,   558,   559,
     560,   560,   561,   561,   561,   562,   562,   563,   564,   564,
     565,   566,   566,   567,   568,   568,   569,   570,   571,   571,
     572,   573,   573,   574,   574,   575,   575,   576,   576,   577,
     577,   578,   579,   579,   580,   581,   581,   582,   582,   583,
     584,   585,   585,   586,   586,   587,   587,   588,   588,   588,
     589,   589,   589,   589,   589,   590,   590,   590,   590,   591,
     591,   591,   591,   591,   591,   592,   592,   593,   593,   594,
     594,   594,   595,   596,   597,   597,   597,   597,   597,   597,
     597,   597,   598,   599,   599,   599,   599,   599,   599,   600,
     600,   601,   601,   602,   602,   603,   603,   603,   604,   604,
     605,   605,   605,   606,   606,   607,   607,   607,   608,   608,
     609,   609,   609,   610,   610,   610,   610,   610,   611,   612,
     612,   613,   613,   613,   614,   614,   614,   615,   615,   616,
     616,   617,   617,   618,   618,   619,   619,   619,   619,   619,
     619,   619,   619,   619,   619,   619,   619,   619,   619,   619,
     619,   619,   619,   619,   619,   619,   619,   619,   619,   619,
     619,   619,   619,   619,   619,   619,   619,   619,   619,   619,
     619,   619,   619,   619,   619,   619,   619,   619,   620,   620,
     621,   621,   621,   621,   621,   622,   622,   622,   622,   622,
     622,   622,   622,   622,   623,   623,   624,   624,   624,   624,
     624,   624,   624,   624,   624,   625,   625,   625,   626,   626,
     627,   628,   628,   628,   629,   629,   629,   629,   629,   629,
     629,   629,   629,   629,   629,   630,   630,   630,   630,   631,
     631,   632,   632,   632,   632,   632,   633,   634,   634,   635,
     635,   635,   635,   636,   637,   638,   638,   639,   639,   640,
     640,   641,   641,   642,   642,   643,   643,   643,   643,   644,
     644,   645,   645,   645,   645,   645,   645,   645,   645,   645,
     645,   645,   646,   646,   646,   647,   647,   648,   648,   648,
     649,   649,   649,   650,   650,   650,   650,   650,   650,   650,
     650,   651,   651,   651,   652,   652,   652,   652,   652,   653,
     653,   653,   653,   653,   653,   654,   654,   655,   655,   655,
     655,   655,   655,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   656,   656,
     656,   656,   656,   656,   656,   656,   656,   656,   657,   658,
     659,   659,   660,   660,   660,   660,   660,   660,   660,   660,
     660,   660,   660,   661,   661,   661,   661,   661,   661,   661,
     661,   661,   661,   661,   662,   662,   663,   664,   664,   665,
     665,   666,   667,   668,   669,   670,   671,   672,   672,   672,
     672,   672,   672,   672,   672,   672,   672,   672,   672,   672,
     672,   672,   672,   672,   672,   672,   672,   672,   672,   672,
     673,   673,   674,   675,   675,   675,   676,   676,   677,   677,
     677,   677,   677,   677,   677,   677,   678,   678,   678,   678,
     678,   678,   678,   678,   678,   678,   678,   678,   678,   678,
     678,   678,   678,   678,   678,   678,   678,   678,   678,   678,
     678,   678,   678,   678,   678,   679,   679,   679,   680,   680,
     681,   681,   682,   682,   682,   683,   683,   684,   684,   684,
     684,   684,   684,   684,   684,   685,   685,   686,   686,   687,
     687,   687,   688,   688,   689,   689,   689,   689,   689,   689,
     689,   689,   689,   689,   689,   689,   689,   689,   689,   689,
     689,   689,   689,   689,   689,   689,   689,   689,   689,   690,
     690,   691,   691,   692,   692,   693,   693,   694,   694,   694,
     694
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
       4,     5,     4,     4,     4,     4,     4,     4,     3,     2,
       4,     4,     3,     3,     3,     3,     3,     3,     3,     2,
       4,     2,     4,     4,     4,     4,     2,     3,     4,     2,
       3,     4,     2,     3,     4,     5,     5,     2,     2,     2,
       2,     2,     2,     2,     2,     4,     4,     5,     3,     2,
       3,     2,     3,     2,     3,     2,    11,    13,    14,     5,
       2,     2,     7,     9,    11,    12,     2,     5,     6,     2,
       5,     2,     5,     4,     4,     3,     3,     3,     3,     3,
       3,     2,     3,     3,     2,     2,     3,     3,     3,     3,
       4,     4,     3,     3,     5,     4,     5,     6,     7,     3,
       5,     3,     5,     4,     5,     6,     7,     3,     2,     3,
       2,     2,     3,     2,     2,     2,     3,     2,     3,     2,
       2,     2,     2,     2,     1,     2,     3,     3,     3,     2,
       5,     3,     3,     3,     2,     3,     2,     3,     4,     4,
       5,     3,     3,     2,     4,     2,     3,     2,     3,     2,
       5,     2,     2,     2,     2,     2,     2,     3,     4,     4,
       8,     4,     4,     3,     3,     5,     4,     5,     2,     3,
       3,     3,     2,     4,     1,     3,     1,     2,     2,     2,
       3,     5,     6,     6,     3,     4,     6,     5,     3,     4,
       3,     1,     2,     6,     2,     2,     2,     2,     2,     4,
       5,     4,     2,     3,     2,     3,     2,     4,     1,     2,
       2,     2,     5,     5,     6,     7,     3,     1,     1,     2,
       1,     1,     1,     2,     1,     1,     2,     1,     4,     5,
       3,     3,     6,     5,     6,     3,     2,     5,     6,     2,
       2,     7,     9,     3,     2,     6,     3,     1,     2,     3,
       2,     3,     1,     2,     4,     2,     4,     5,     2,     3,
       4,     5,     2,     3,     6,     2,     3,     6,     3,     1,
       2,     4,     5,     6,     3,     2,     4,     1,     2,     3,
       1,     2,     4,     5,     6,     3,     2,     4,     3,     2,
       4,     3,     2,     4,     3,     2,     4,     3,     1,     2,
       3,     3,     4,     3,     2,     3,     2,     5,     2,     3,
       4,     4,     9,     2,     2,     9,     2,     3,     4,     3,
       3,     7,     2,     3,     4,     2,     2,     5,     7,    10,
       3,     4,     4,     5,     4,     2,     3,     3,     4,     3,
       2,     9,     6,     2,     2,     3,     9,     3,     9,     2,
       3,     4,     3,     3,     2,     3,     2,     7,     9,     8,
      10,     3,     5,     1,     2,     6,     7,     8,     9,     1,
       2,     1,     2,     2,     3,     6,     4,     7,     2,     3,
       6,     4,     7,     2,     3,     6,     4,     7,     2,     3,
       2,     2,     3,     2,     3,     5,     4,     6,     4,     2,
       3,     2,     2,     3,     5,     3,     2,     5,     4,     3,
       4,     1,     2,     3,     2,    16,    19,    19,    20,    23,
      23,     9,    14,    16,    30,    12,    13,    14,     5,     6,
      16,    12,     3,     4,     5,     6,     6,     7,    11,     8,
       9,    10,     9,    10,    11,    10,    11,    12,    11,    12,
      13,     5,     6,     6,     7,     7,     8,     4,     3,     2,
       3,     4,     5,     6,     5,     4,     5,     5,     6,     7,
       7,     8,     7,     8,     4,     3,     2,     5,     6,     7,
       8,    10,    11,     6,     9,     2,     5,     7,     3,     2,
       4,     2,     5,     7,     2,     5,     6,     7,     8,    10,
      11,    13,    14,     6,     9,     3,     3,     3,     3,     3,
       2,     5,     4,     3,     7,     6,     4,     1,     2,     4,
       3,     5,     4,     3,     3,     5,     4,     2,     2,     2,
       2,    11,     4,     2,     2,    11,    14,    12,    15,     2,
       7,     2,     4,     6,     7,     6,     8,     7,     5,     5,
       7,     8,     2,     4,     3,     2,     3,     2,     4,     3,
       2,     4,     5,     2,     3,     4,     5,     5,     5,     6,
       6,     2,     3,     4,     2,     4,     6,     6,     8,     2,
       2,     2,     4,     3,     5,     1,     2,     2,     4,     4,
       4,     4,     4,     2,     3,     2,     3,     3,     3,     4,
       5,     6,     7,     4,     3,     3,     3,     4,     5,     2,
       4,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     4,     4,     3,     3,     4,     4,     4,     4,
       4,     4,     3,     3,     3,     4,     3,     4,     4,     4,
       4,     4,     2,     3,     3,     3,     3,     3,     3,     3,
       2,     2,     3,     2,     3,     2,     3,     4,     3,     3,
       4,     4,     5,     5,     6,     3,     4,     3,     3,     3,
       3,     3,     2,     3,     3,     3,     3,     3,     4,     4,
       4,     4,     4,     4,     4,     3,     2,     2,     3,     3,
       2,     3,     4,     5,     6,     7,     2,     3,     4,     3,
       3,     3,     3,     2,     2,     3,     4,     5,     3,     4,
       3,     3,     3,     3,     5,     5,     4,     4,     6,     6,
       3,     3,     4,     3,     3,     3,     3,     3,     3,     2,
       3,     4,     1,     2,     3,     4,     5,     1,     2,     3,
       2,     3,     4,     2,     2,     3,     4,     5,     4,     4,
       4,     2,     2,     2,     1,     2,     4,     6,     5,     1,
       2,     4,     3,     2,     3,     2,     2,     2,     3,     3,
       3,     4,     4,     5,     7,     5,     7,     4,     5,     3,
       4,     4,     5,     6,     7,     8,     3,     4,     6,     2,
       3,     4,     2,     2,    10,    12,     3,     4,     2,     4,
       7,     9,     9,     8,    11,    12,     2,     9,    10,     9,
      10,     9,    10,     7,     7,     7,     8,     8,     7,     7,
       8,     8,     9,    10,    10,    11,    12,    24,     7,     8,
       9,    10,    11,     5,     5,     2,     4,     5,     0,     2,
       0,     2,     4,     6,     8,     2,     3,     2,     2,     2,
       2,     2,     2,     2,     1,     2,     3,     2,     3,     2,
       3,     2,     2,     3,     2,     2,     3,     4,     2,     2,
       2,     2,     3,     4,     2,     2,     2,     2,     2,     3,
       4,     2,     2,     3,     4,     2,     2,     2,     2,     2,
       3,     2,     3,     2,     3,     1,     1,     1,     1,     1,
       1
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
      38,    40,    76,    75,    83,   304,    39,    42,    69,    70,
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
    1075,  1076,     0,   721,  1078,  1080,  1077,  1079,     0,     0,
       0,     0,   316,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   505,     0,     0,   719,   524,     0,   223,   442,
       0,   217,   331,   973,   674,   424,   707,     0,   946,   246,
     419,   333,   189,     0,   912,   917,     0,     0,     0,   732,
     305,     0,   709,     0,     0,     0,   314,     0,     0,     0,
     438,     0,   534,     0,   924,   201,     0,   665,     0,   523,
     249,   382,   149,     0,     0,     0,     0,   209,     0,   923,
       0,     0,     0,     0,     0,   377,   390,   579,   515,     0,
     520,     0,     0,   529,     0,     0,     0,     0,     0,     0,
       0,     0,   240,   568,     0,   212,     0,   323,   751,     0,
     329,     0,   206,     0,   384,     0,     0,   947,     0,     0,
       0,   713,     0,     0,   264,     0,     0,   265,     0,   357,
       0,     0,     0,     0,   754,   743,     0,     0,     0,   445,
       0,   386,     0,     0,     0,   972,   251,   153,   563,     0,
       0,   759,   309,     0,     0,   435,     0,     0,   332,     0,
       0,   375,   374,   558,   656,   325,     0,     0,     0,   553,
     179,   978,     0,   376,     0,     0,   986,  1015,     0,     0,
       0,     0,     0,   199,   327,     0,   378,   391,     0,     0,
       0,   671,  1025,  1071,  1039,  1042,   740,   735,   737,  1035,
    1037,   261,     0,     1,     2,     4,     0,     0,     0,     0,
       0,     0,     0,   170,   171,   168,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   218,   219,   220,   221,
     222,   224,     0,   241,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   334,
     335,   336,     0,     0,     0,   358,   359,   371,     0,     0,
       0,   416,     0,     0,   420,     0,     0,     0,     0,     0,
       0,     0,     0,   455,     0,   466,     0,   469,     0,   472,
       0,   475,     0,   484,   486,   488,   493,     0,   496,     0,
     502,     0,   506,     0,     0,     0,     0,     0,     0,     0,
       0,   536,     0,   594,     0,   639,     0,     0,     0,     0,
     669,     0,     0,     0,   690,     0,   708,   710,   714,     0,
       0,     0,     0,     0,     0,  1077,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     761,   760,   933,   932,   931,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   969,
       0,     0,     0,     0,     0,     0,     0,  1020,     0,     0,
    1018,  1020,     0,     0,  1034,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1041,     0,     0,     0,   258,
     543,     0,   364,     0,     0,   188,   483,   485,     0,   317,
       0,     0,   478,   480,     0,   481,     0,     0,     0,   499,
     500,     0,   685,   326,   925,     0,   443,     0,     0,     0,
       0,   920,   913,     0,   918,     0,     0,   910,   306,   465,
     454,   533,   551,     0,   908,     0,   489,     0,     0,   439,
       0,   535,   321,   686,     0,   411,   591,     0,   589,     0,
     448,   449,     0,   210,   468,   942,     0,   944,     0,   474,
     471,   580,     0,     0,   517,   516,   519,   532,   530,     0,
       0,     0,   415,     0,     0,   322,     0,   569,     0,     0,
     255,   213,   752,     0,   593,   207,   383,   308,   689,     0,
     311,   976,   268,     0,   266,   269,     0,   267,     0,   525,
     527,     0,   744,     0,     0,   446,     0,   388,   385,     0,
     510,     0,     0,     0,   564,   257,   315,     0,   459,   460,
       0,   312,   313,   638,     0,     0,   559,   256,   260,   259,
       0,   554,     0,     0,     0,     0,     0,     0,   167,     0,
       0,   640,     0,     0,     0,   687,   688,   655,     0,   328,
     426,   427,     0,   668,   431,   432,     0,   307,     0,     0,
       0,     0,     0,   173,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   196,     0,   198,   197,     0,   194,
     195,   193,   192,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   229,     0,
     231,   233,     0,   235,     0,     0,     0,     0,     0,     0,
       0,     0,   262,     0,     0,     0,     0,     0,     0,   302,
     303,   301,   293,   290,   291,   300,   297,   272,     0,   299,
     294,     0,   288,   295,     0,     0,     0,   354,   344,     0,
       0,   356,     0,   337,     0,   348,   343,     0,     0,     0,
       0,     0,   360,     0,   368,     0,   370,   372,     0,   423,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   494,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   541,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   734,     0,   736,     0,   739,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     765,   773,     0,     0,     0,     0,     0,   789,     0,     0,
       0,     0,   831,   833,     0,     0,   866,     0,     0,     0,
     876,     0,   883,   767,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   830,     0,     0,     0,     0,
     870,     0,   822,     0,     0,   852,     0,     0,     0,     0,
       0,     0,     0,     0,   763,     0,     0,     0,     0,   775,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   884,     0,     0,     0,     0,     0,     0,
       0,   835,     0,     0,   867,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   948,     0,     0,     0,   966,
       0,     0,     0,     0,   949,     0,     0,   959,   950,     0,
       0,     0,  1018,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1033,  1031,
    1027,  1028,  1029,  1030,  1032,  1069,  1026,     0,  1056,  1057,
    1049,  1044,  1045,  1050,  1051,  1048,  1054,  1065,  1066,  1058,
    1061,  1067,  1062,  1068,  1055,  1036,  1038,  1040,  1043,  1073,
    1072,   544,     0,     0,   365,     0,     0,     0,   318,   319,
       0,   479,     0,   482,     0,     0,     0,   926,     0,     0,
       0,   346,     0,   914,   919,   911,   921,   552,     0,   909,
     490,   491,   939,     0,     0,   381,     0,   592,     0,   590,
     450,     0,   943,   945,     0,   582,     0,   581,     0,     0,
     586,     0,   930,   934,     0,     0,     0,     0,   148,     0,
       0,   571,     0,   570,     0,   573,     0,   753,     0,   270,
     271,     0,     0,   928,   929,   387,   389,   512,     0,   511,
    1022,     0,     0,     0,     0,   461,     0,   379,     0,     0,
       0,     0,     0,     0,   979,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   641,     0,   645,     0,
       0,     0,   654,   428,   430,     0,   433,     0,     0,     0,
       0,   410,   549,     0,   172,     0,     0,     0,     0,   180,
     182,   183,   184,   185,   186,   187,   190,   191,   200,   202,
     205,   204,   203,   208,   211,     0,     0,   214,     0,   228,
     230,   232,   234,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   253,   254,   263,   273,     0,   281,     0,   279,
       0,   292,   298,   287,   296,   289,   324,   350,   352,     0,
     351,   341,   349,   339,   338,     0,     0,   342,     0,     0,
     369,     0,     0,     0,   700,     0,     0,     0,   436,     0,
     440,     0,     0,     0,   457,     0,     0,     0,     0,   497,
       0,   503,     0,     0,   514,   518,     0,     0,     0,   531,
       0,     0,     0,     0,   612,     0,     0,     0,     0,     0,
       0,     0,     0,   697,     0,     0,     0,     0,     0,     0,
     693,     0,     0,     0,     0,     0,   722,     0,     0,   733,
     738,   741,     0,     0,     0,   745,     0,     0,   755,     0,
     839,     0,     0,     0,   856,     0,     0,   903,   855,     0,
       0,   838,   836,     0,   849,   857,     0,     0,   834,   905,
     906,   762,   766,   774,   904,     0,     0,   785,   786,     0,
     784,     0,   829,   828,   832,   847,   848,   804,   805,     0,
     778,     0,   882,   881,   854,     0,     0,     0,     0,     0,
     893,     0,   892,     0,   891,   890,   827,   826,   900,   907,
       0,     0,   871,   877,     0,   853,   823,   813,   812,   869,
     901,     0,     0,     0,     0,   850,   825,   824,   776,   777,
     793,   796,   794,   795,   797,   798,   800,   799,   801,   791,
     792,   880,   816,     0,   814,     0,   879,   885,     0,   851,
     865,     0,   888,   845,     0,   868,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   957,
       0,   967,   960,   961,     0,   951,   952,     0,   970,     0,
     977,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1016,     0,
    1021,  1019,  1070,  1046,  1052,  1059,  1063,  1074,     0,     0,
     367,     0,     0,   487,   320,   477,     0,     0,     0,   927,
     330,     0,     0,   915,   922,   704,   940,     0,   542,   380,
       0,     0,     0,   583,     0,     0,   585,   935,     0,   413,
       0,     0,     0,     0,   572,     0,   574,     0,     0,   147,
     397,   393,     0,     0,   513,     0,     0,     0,     0,   310,
       0,     0,     0,     0,     0,     0,     0,     0,   175,     0,
       0,     0,     0,   155,     0,     0,     0,     0,     0,     0,
     642,     0,     0,   646,     0,   647,     0,   429,     0,     0,
     150,     0,     0,   550,     0,   174,     0,   408,   177,     0,
     181,   216,   215,   225,     0,     0,     0,   712,     0,     0,
       0,   247,   250,   252,     0,   275,     0,     0,   283,     0,
       0,   355,   345,     0,     0,     0,     0,     0,     0,     0,
     702,     0,   699,     0,   437,   441,     0,     0,   456,   458,
     467,   470,   473,   476,   498,   504,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   613,     0,
       0,     0,   637,     0,     0,   696,   698,   644,   657,     0,
       0,   666,     0,   670,   672,     0,   675,     0,     0,   692,
       0,     0,     0,     0,   728,     0,     0,   729,     0,     0,
     742,   747,     0,   746,     0,   748,     0,     0,   840,     0,
     841,     0,   861,   790,   860,   862,   837,   863,   864,   808,
     807,   787,   783,     0,   872,     0,   779,     0,   858,   859,
     902,     0,   897,     0,     0,     0,   896,   802,   803,   878,
     818,   820,   819,   764,   817,   815,   886,     0,   889,   846,
     769,   771,   768,   772,   770,   806,   821,   809,   810,   811,
     938,   955,     0,     0,   958,   962,     0,   953,     0,   971,
       0,  1013,  1014,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1017,  1047,  1053,
    1060,  1064,     0,     0,     0,   366,   373,     0,     0,     0,
       0,     0,   916,     0,   706,     0,     0,     0,   451,     0,
       0,     0,   412,   414,   418,   937,     0,     0,   576,   398,
       0,   400,   401,   402,     0,   404,   405,   407,     0,   394,
       0,  1023,     0,     0,   566,     0,     0,   462,     0,     0,
     561,     0,     0,   556,     0,     0,     0,     0,   176,     0,
       0,     0,     0,     0,   643,   648,     0,     0,     0,   434,
     152,   151,   154,   703,   409,     0,     0,     0,     0,     0,
       0,     0,   248,   276,   274,   284,   282,   280,   353,     0,
     392,     0,   361,     0,   417,     0,     0,   701,   425,   444,
     447,   507,     0,   522,     0,     0,     0,     0,     0,   608,
       0,     0,   614,     0,   631,     0,     0,     0,     0,     0,
       0,     0,     0,   663,   658,     0,     0,     0,   683,   676,
       0,   691,     0,     0,     0,     0,     0,   723,     0,   725,
       0,     0,     0,   750,   749,   756,     0,   757,   842,     0,
     843,   788,   873,     0,   780,     0,     0,   895,   894,     0,
     887,     0,   968,   963,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   501,     0,   347,
     941,   705,     0,   452,   588,     0,   584,   936,   578,   575,
       0,   399,   403,   406,   396,   395,     0,     0,     0,     0,
     463,     0,     0,     0,     0,     0,   156,   157,     0,     0,
       0,     0,     0,   161,   649,   650,     0,   652,     0,   178,
       0,     0,     0,   239,     0,     0,     0,   277,     0,   285,
       0,   363,   362,     0,   421,     0,     0,     0,     0,     0,
       0,   537,   609,     0,     0,   616,   633,     0,     0,   615,
       0,   632,     0,     0,     0,     0,     0,   659,     0,   667,
     673,   677,     0,   695,     0,     0,     0,   720,   724,   727,
       0,     0,   730,     0,   844,   874,     0,   781,     0,   899,
     898,   956,   964,     0,   954,   980,     0,     0,     0,     0,
       0,   993,   994,     0,  1008,     0,   998,   999,     0,     0,
       0,     0,   995,     0,     0,     0,     0,     0,   545,     0,
       0,     0,   453,   587,   577,  1024,   565,     0,   464,   560,
       0,   555,     0,   158,   160,     0,     0,   162,   163,     0,
     651,   653,     0,   226,     0,     0,   242,     0,   278,   286,
     340,     0,   508,     0,     0,     0,     0,     0,   539,     0,
       0,   635,     0,   617,   634,     0,     0,     0,     0,     0,
     660,     0,     0,   678,     0,     0,   694,     0,     0,   731,
     726,   758,   875,   782,   965,     0,     0,   983,     0,     0,
     996,  1009,     0,  1000,     0,     0,  1001,     0,     0,   997,
       0,   546,     0,   492,   495,     0,   567,   562,   557,     0,
     164,   166,     0,   227,     0,     0,     0,   422,     0,   521,
     526,   528,     0,   538,     0,   619,     0,     0,     0,   636,
       0,     0,     0,     0,   664,     0,   684,     0,     0,     0,
     982,   981,     0,   987,     0,   989,     0,  1010,     0,     0,
       0,  1002,   991,     0,     0,   547,   974,     0,     0,     0,
       0,   243,     0,     0,   540,     0,     0,   620,     0,     0,
     622,     0,     0,     0,     0,     0,   601,     0,   661,     0,
     679,     0,     0,     0,     0,   988,   990,  1011,     0,  1003,
    1004,     0,   992,   548,     0,     0,     0,     0,     0,     0,
     509,     0,   625,     0,   621,     0,     0,   623,     0,     0,
       0,     0,     0,     0,   662,   680,     0,     0,     0,   984,
       0,  1012,     0,  1005,     0,   975,     0,     0,     0,   711,
     244,     0,     0,   626,     0,   628,     0,   624,     0,   618,
       0,     0,     0,     0,     0,   715,     0,   985,     0,  1006,
     159,     0,   236,     0,   245,   605,     0,   627,   629,     0,
       0,   611,     0,   681,     0,   717,     0,     0,     0,   165,
       0,   606,     0,   630,     0,     0,   682,     0,     0,     0,
     237,     0,   607,     0,     0,   602,     0,     0,   716,     0,
     238,     0,     0,     0,   718,     0,   610,   603,   595,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   596,   597,     0,     0,   598,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   599,   600,
       0,  1007,     0,     0,     0,     0,     0,     0,   604
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   170,   171,   172,   173,   174,   175,   176,   177,   533,
     534,   968,   535,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   566,  1833,   567,   568,   569,   570,  1028,
    1031,   188,   573,   189,   190,   191,   192,   193,   194,   195,
     196,   197,   198,   579,   199,   586,   200,   201,   202,   203,
     204,   205,   206,  1297,   599,  1066,  1070,   207,   605,   606,
     208,   610,   209,   210,   211,   212,   213,   214,   215,   216,
     217,   218,   936,   219,   220,   600,   221,  1781,  2048,   536,
     222,   223,   224,   611,   225,   226,   227,   228,   980,   981,
     229,   984,   985,   230,   231,   232,   233,   234,   880,   881,
     235,   623,  1523,   236,   948,   949,   237,   625,   238,   627,
     239,   629,   240,   631,   832,   833,   241,   242,   243,   244,
     245,   246,   247,   248,   636,   249,   638,   250,   251,   252,
     640,   253,   642,   254,   255,   256,   257,   258,   259,   260,
     261,   262,   263,   819,  1441,   861,   264,   961,   265,   956,
     266,   944,   267,   907,   908,  1381,   268,   891,   892,  1365,
     269,   875,   270,   653,   271,   272,   273,   274,   275,   276,
     277,   660,   278,   279,   280,   281,   664,   655,  1552,   820,
    1442,   862,   876,   282,   283,   571,   284,   668,   285,   286,
     287,   288,   289,   290,   291,   292,   293,   294,  1159,   760,
     761,   295,   864,   296,   368,   297,  1372,  1373,   298,  1351,
    1352,   299,   885,   300,   887,   301,   302,   779,   303,   304,
     305,   306,   307,   308,  1293,  1290,   309,   310,   793,   311,
     312,   313,   314,   812,   794,   315,   818,  1443,  1374
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -1763
static const yytype_int16 yypact[] =
{
    5586,  -149,    37,   -97,   447,   -60,  1223,  1223,  1223,   437,
     -51,   602,    75,   116,   447,   447,   139,   150,   447,   447,
     158,   210,   447,   219,   663,   297,   315,   325,   341,   349,
     350,   380,   391,   409,   429,   440,   442,    21,   455,   705,
     465,   476,   488,   491,   498,   507,   447,   447,   710,   771,
     515,   522,  -189,   527,   447,   537,   554,   572,   592,   593,
     835,   616,   870,   619,   626,   627,   629,   640,   642,   655,
     657,   662,  -135,   451,   664,   672,   447,   675,   938,   679,
     684,   447,   687,   697,   707,    66,   990,   712,  1029,   723,
     730,  1043,  1152,   447,   734,   740,   742,   447,   755,    87,
     795,  1017,   760,   761,   447,   447,   766,   769,  1182,   447,
     447,  1243,  1249,   888,   785,   787,   798,   799,   806,   808,
    1223,   814,  1268,   819,   821,  1223,  1223,   824,   833,   834,
     836,   838,   843,   850,   865,   979,   871,   882,   890,   896,
     -76,   899,  1270,   906,   919,   447,   447,  1223,   447,   447,
     920,   924,  1223,   930,   931,   936,   937,   941,   956,   957,
     960,   963,   969,   976,   977,   988,   996,   998,  1000,   907,
     753,  5140, -1763, -1763, -1763,   985,   447,  1119,    62, -1763,
     -23,   447,    29,  1223,  1223,   447,   536,   447,   447,   447,
    1223,  1113, -1763, -1763, -1763, -1763, -1763, -1763,   772,   172,
    1037, -1763, -1763, -1763, -1763,    33, -1763,   777,     5, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763,   447, -1763,   -75,   447, -1763, -1763,
     -44,    20,   447,   447, -1763,   447, -1763,   447,   447,   447,
     447, -1763, -1763,   447,   447,   447, -1763, -1763,   447,   447,
   -1763, -1763,   447,   447,   779,  1049,   447,   447,   447,  1055,
   -1763, -1763,    26, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
     447,   447,   447, -1763, -1763,   447,   447,   447,   447,   447,
   -1763,   447,   447,   447,   447,   447,   447,   -81,  1223,  1288,
     447,    50, -1763,   447,  4265, -1763, -1763,   609,  1223, -1763,
   -1763, -1763,     4, -1763, -1763,  1020,   447,  -123,   447, -1763,
     -26,   356,   356,    34,   356,   827,   447,  1035,  1246,  1254,
   -1763, -1763,  1041, -1763, -1763, -1763, -1763, -1763,  1044,  1048,
    1223,  1368, -1763,  1223,   447,   447,  1052,  1223,  1223,  1070,
    1072,   447, -1763,  1075,  1079, -1763, -1763,   464, -1763, -1763,
    1082, -1763,  1223,  1159, -1763, -1763, -1763,  1223, -1763, -1763,
   -1763, -1763, -1763,  1223,  1223,  1213,  1223,  1096,    98, -1763,
   -1763,  1085, -1763,   447,   447,   447, -1763,  1223,  1420,  1088,
   -1763,  1091,  1118,  1101, -1763, -1763,  1110, -1763,  1117, -1763,
   -1763,  1223, -1763,   447,   447,  1121,   447, -1763,  1125, -1763,
     447,  1223,  1223,   447,   447, -1763, -1763,   447,   447,  1130,
   -1763,  1135,   447,   447,  1136,  1223,   447,  1137,  1223,   447,
    1142,  1223, -1763,   447,  1145, -1763,  1150, -1763, -1763,  1425,
   -1763,   447, -1763,  1166, -1763,  1173,  1179, -1763,   447,   447,
    1184, -1763,  1200,  1434, -1763,  1218,  1442, -1763,  1221, -1763,
     447,  1224,  1225,   447, -1763, -1763,  1227,  1229,  1234, -1763,
    1240,   447,  1244,   387,  1060, -1763, -1763, -1763,  1281,   447,
    1247, -1763, -1763,  1223,   447, -1763,  1248,  1250, -1763,   447,
    1223, -1763, -1763,  1408, -1763, -1763,  1251,   447,   447,  1417,
   -1763, -1763,   447, -1763,    36,  1443, -1763, -1763,    40,   447,
    1259,  1260,  1449, -1763, -1763,  1263, -1763, -1763,   447,   447,
     447, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763,  1265, -1763, -1763, -1763,    65,  1204,  1271,  1223,
     466,  1223,  1215, -1763, -1763, -1763, -1763,  1178,  1183,  1391,
    1399,  1402,  1406,  1407,  1303,  1422,  1313,  1315,   447,  1326,
    1355,  1363,  1365,  1223,   447,   447,   447,   447,  1223,  1223,
    1491,   447,   447,   447,   447,   447, -1763,   447,   447,   447,
     447, -1763,   -28, -1763,  1223,   447,  1223,   251,  1129,  1378,
     340,   434,   463,  1496,  1498,  1500,     8,   447,  1223,   505,
    1223,  1286,  1389,  1223,  1392,  1393,   -66,  1223,  1297, -1763,
   -1763, -1763,  1223,  1395,  1223, -1763, -1763, -1763,  1398,  1510,
    1454, -1763,  1223,   447, -1763,  -101,  1552,   447,  1223,   447,
    1223,  1223,  1223, -1763,   447, -1763,   447, -1763,   447, -1763,
     447, -1763,   447, -1763, -1763, -1763, -1763,  1404, -1763,   447,
   -1763,   447, -1763,  1223,  1308,  1410,  1223,  1223,  1223,  1411,
     447, -1763,   -94, -1763,    30, -1763,   447,   447,   447,   447,
   -1763,   447,   447,   447, -1763,  1223, -1763, -1763, -1763,  1499,
     447,   447,  1423,  1223,  1424,  1493,   447,  1426,   447,   447,
    1288,   -31,     3,  1331,  1429,  1311,  1385,    35,  1223,  1433,
    1223,  1343,  1373,   447,  1455,  1377,   447,   447,   -52,   173,
    1431,  1479,   447,  1223,  1482,  1223,  1223,  1386,  1457,  1485,
    1471,   695,  1437,   447,   574,   776,   287,  1492,   447,   447,
     447,   447,  1472,  1513,   447,   -46,   447,   242,   447,  1223,
    1223,  1223,   447,   540,  1421,   300,   -61,  1223,   447,   447,
     447,  1223,   447,   447,  1432,   447,   447,   447,  1516,   447,
    1223,  1223,  1223,   447,   447,   447,  1530,  1586,  1223,  -132,
     333,   461, -1763, -1763, -1763,  1223,  1534,  1223,   447,   447,
     559,  1223,   447,   447,  1536,  1223,   447,  1538,  1543, -1763,
    1657,   447,  1663,   -80,   447,   447,   447, -1763,   447,   447,
   -1763, -1763,  1554,  1556,  1223,  1460,  1461,   447,  1670,  1671,
     447,   447,  1223,  1467,  1476,  1481,  1692,  1501,  1502,  1705,
    1503,  1505,  1562,  1597,  1611,  1223,  1612,  1730,  1624,   447,
   -1763,   -92, -1763,   561,  1223, -1763, -1763, -1763,  1223, -1763,
    1627,  1473, -1763,   447,  1223,   447,   447,  1223,  1223, -1763,
   -1763,  1630, -1763, -1763, -1763,   566, -1763,   447,  1633,  1634,
    1223, -1763,  1223,  1223, -1763,   286,   447, -1763, -1763, -1763,
   -1763,   447, -1763,  1223, -1763,  1635, -1763,  1637,  1223, -1763,
    1640, -1763, -1763, -1763,   576,   447, -1763,   447,   447,   447,
     447, -1763,   447, -1763, -1763, -1763,  1641, -1763,  1650, -1763,
   -1763,   447,   447,   337,   447, -1763, -1763,   447,   447,  1223,
    1223,  1223, -1763,  1223,  1652, -1763,  1223,   447,   447,   583,
   -1763, -1763, -1763,  1653, -1763, -1763, -1763, -1763, -1763,   447,
   -1763, -1763, -1763,  1676, -1763, -1763,  1678, -1763,   447, -1763,
   -1763,   447, -1763,  1223,  1223, -1763,  1480, -1763, -1763,  1490,
   -1763,  1683,  1684,   447,  1494,   447, -1763,  1223,   447, -1763,
     447, -1763, -1763, -1763,  1686,   447,  1495, -1763,   447,   447,
     447,  1522,  1687,  1223,  1223,  1223,  1223,  1223,  1223,   903,
    1689, -1763,  1690,  1223,  -179, -1763, -1763, -1763,  1694, -1763,
     447, -1763,   599, -1763,   447, -1763,   447, -1763,   447,   447,
    1223,   447,  1696, -1763,  1223,  1223,   447,   -12,  1699,  1702,
    1710,  1711,  1715,  1716, -1763,  1717, -1763, -1763,  1719, -1763,
   -1763, -1763, -1763,  1720,  1721,  1722,  1723,  1734,  1735,  1736,
     447,   447,  1737,  1738,  1739,  1740,  1741,  1742, -1763,  1223,
   -1763, -1763,   447, -1763,   447,  1223,  1223,  1288,  1223,    -7,
    1743,   447, -1763,   447,   447,   447,   447,   447,   447, -1763,
   -1763, -1763, -1763, -1763,   447, -1763,   447, -1763,   447, -1763,
   -1763,  1647, -1763, -1763,  1223,  1745,  1746,  1223, -1763,  1223,
    1747, -1763,  1748, -1763,  1749, -1763, -1763,  1750,   607,  1223,
    1751,  1223, -1763,  1223, -1763,  1752, -1763, -1763,  1223, -1763,
     447,   447,   322,   447,  1223,  1753,  1223,  1754,  1223,  1223,
     447,   447,   447,   447,   447, -1763,  1514,  1519,  1223,  1755,
     447,  1223,  1223,  1223,   447, -1763,   447,  1223,    48,   447,
    1223,  1223,  1223,   447,   447,   447,   447,   447,   447,   447,
     620,  1223,  1223,  1223,   -58, -1763,  1756, -1763,  1757, -1763,
    1521,  1223,  1223,   208,   447,   622,   546,  -122,  -104,  1685,
    1758,   326,  1759,  1524,  1761,   388,  1762,  1763,  1764,  1539,
   -1763, -1763,  1765,  1767,  1223,   447,  1769, -1763,  1771,   342,
    1772,  1787, -1763, -1763,  1788,  1789, -1763,  1790,  1791,  1792,
   -1763,  1223, -1763, -1763,   644,  1793,  1795,   467,  1223,   639,
     649,  1796,  1799,  1800,  1802, -1763,  1803,  1804,   447,  1223,
   -1763,  1805, -1763,   668,  1807, -1763,  1808,  1810,  1811,  1812,
    1813,  1223,  1223,   447, -1763,  1223,  1814,  1815,  1816, -1763,
    1817,  1818,  1819,  1820,  1821,  1822,  1823,  1825,  1827,  1828,
    1829,  1830,  1831, -1763,  1834,   714,   716,  1836,   729,  1837,
    1840, -1763,   130,   732, -1763,  1843,   447,   447,   447,   447,
    1223,   447,   447,   510,  1223, -1763,  1223,  1223,  1546, -1763,
    1844,  1845,  1567,  1846, -1763,   735,  1569, -1763, -1763,  1847,
    1223,  1958, -1763,  1223,  1223,  1223,  1223,  1223,  1223,  1223,
    1223,  1223,  1223,  1223,  1223,  1223,  1223,  1580, -1763, -1763,
    1959, -1763, -1763,  1223,  1959, -1763, -1763,  1850, -1763, -1763,
   -1763, -1763,  1961, -1763,   447, -1763, -1763, -1763, -1763,  1962,
   -1763, -1763,   447, -1763, -1763, -1763, -1763, -1763, -1763,   447,
   -1763, -1763,   -73,   447, -1763,   751,  1223,  1853, -1763, -1763,
    1854, -1763,  1855,   447,  1856,  1859,  1223, -1763,  1860,  1861,
     447, -1763,  1223,  1223, -1763, -1763,  1223, -1763,  1862, -1763,
   -1763,  1223, -1763,  1223,   447, -1763,  1863, -1763,  1223,   447,
   -1763,  1288, -1763, -1763,   447, -1763,   359, -1763,   447,  1223,
   -1763,  1864,  1223, -1763,  1223,   756,  1223,  1223, -1763,  1223,
     447, -1763,   361, -1763,   447, -1763,     0, -1763,  1865, -1763,
   -1763,  1581,   447,  1223,  1223, -1763, -1763, -1763,  1866, -1763,
    1824,  1869,   447,  1223,  1872, -1763,  1288, -1763,  1873,   447,
    1223,  1874,   447,  1223, -1763,   758,  1223,  1223,  1223,  1875,
     447,  1223,  1223,  1223,  1223,  1550, -1763,    52, -1763,  1223,
    -183,  -154, -1763, -1763, -1763,  1876, -1763,   447,   447,  1583,
    1223,   447, -1763,  1223, -1763,  1877,   186,  1607,  1878, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763,  1879,  1880, -1763,   447, -1763,
   -1763, -1763, -1763,  1223,  1223,  1881,  1223,  1223,  1223,  1883,
    1884,  1885, -1763, -1763, -1763,   376,  1996,   530,  1997, -1763,
    1998, -1763, -1763, -1763, -1763, -1763, -1763, -1763,  1223,  1889,
   -1763, -1763, -1763, -1763, -1763,  1223,  1223, -1763,  1223,   447,
   -1763,  1223,   -43,  1890, -1763,   424,   447,  1891, -1763,  1893,
   -1763,  1223,  1223,  1608, -1763,  1609,  1610,  1617,  1626, -1763,
    1897, -1763,  1898,   447, -1763,   447,  1223,  1223,  1223,   447,
      67,  1223,  1223,  1223, -1763,  1223,   -55,  1223,  1223,  1899,
    1223,  1223,  1629, -1763,  1632,   128,   773,  1642,   780,   232,
   -1763,   786,  1223,  1223,  1223,  1902, -1763,   -22,    15, -1763,
   -1763, -1763,  1903,   233,   250, -1763,  1904,  1223, -1763,  1223,
   -1763,  1913,   447,  -175, -1763,  1914,  1916, -1763, -1763,  1917,
    1918, -1763, -1763,  1921, -1763, -1763,  1922,  1925, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763,  1926,  1927, -1763, -1763,  1928,
   -1763,  1648, -1763, -1763, -1763, -1763, -1763, -1763, -1763,  1649,
   -1763,  1651, -1763, -1763, -1763,  1929,  1932,  1933,  -138,   447,
   -1763,   447, -1763,   -91, -1763, -1763, -1763, -1763, -1763, -1763,
    1934,  1935, -1763, -1763,  1937, -1763, -1763, -1763, -1763, -1763,
   -1763,  1938,  1939,  1940,  1941, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763,  1942, -1763,  1943, -1763, -1763,   796, -1763,
   -1763,  1944, -1763, -1763,  1945, -1763,  1946,  1947,  1948,  1949,
    1950,  1951,  1952,  1953,  1954,  1955,  1956,  1672,  1223, -1763,
    1960, -1763, -1763, -1763,   810, -1763, -1763,   823, -1763,  1963,
   -1763,  1223,  1968,   837,  1223,  1223,  1223,  1223,  1223,  1223,
    1223,  1223,  1223,  1223,  1223,  1223,  1223,  1223, -1763,  1971,
   -1763, -1763, -1763,  2069,   447,  2083,   447, -1763,   447,   112,
   -1763,  1974,  1977, -1763, -1763, -1763,  1223,  1223,  1223, -1763,
   -1763,  1978,  1223,  1223, -1763, -1763, -1763,  1223,   447, -1763,
     839,  1288,  1979, -1763,  1223,  1223, -1763, -1763,  1223, -1763,
    1980,  1981,  1982,  1983, -1763,  1223, -1763,   447,   289, -1763,
    1703, -1763,  1984,   447, -1763,  1985,  1223,  1986,  1987, -1763,
    1288,  1988,  1223,  1989,  1990,  1223,  1991,  1992, -1763,   447,
     447,   447,   447, -1763,  1993,  1223,  1223,  1223,  1223,  1223,
   -1763,  1994,   840, -1763,  1223, -1763,  1223, -1763,  1679,  1995,
   -1763,  1999,  2000, -1763,  2001, -1763,  2003, -1763, -1763,   447,
   -1763, -1763, -1763,   447,   447,  1223,  1223, -1763,  1223,  1223,
    2005, -1763, -1763, -1763,   447, -1763,   447,   447, -1763,   447,
     447,  1223, -1763,   447,  2006,   447,   841,  2007,   447,  1223,
   -1763,  2009, -1763,  2010, -1763, -1763,  2011,  2012, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763,   904,   909,  1223,  1223,
     447,  1223,   933,  1223,   942,   945,  1223,  1223, -1763,  1223,
    1223,  1223, -1763,  1223,  1223, -1763, -1763, -1763, -1763,   328,
     952, -1763,  1223, -1763, -1763,  1223, -1763,   332,   965, -1763,
    1681,  1223,  1223,  1223, -1763,  1223,  2013, -1763,  1223,   179,
   -1763, -1763,  2014, -1763,  2015, -1763,   973,  2016, -1763,  -168,
   -1763,  2017, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763,  2018, -1763,   974, -1763,  1701, -1763, -1763,
   -1763,   447, -1763,  2019,  2021,   447, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763,  2023, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763,   447,  2026, -1763, -1763,   975, -1763,  1223, -1763,
    1223, -1763, -1763,  1223,  1223,  1223,  1223,  1223,  1223,  1223,
    1223,  1223,  1223,  1223,  1223,  1223,  1223, -1763, -1763, -1763,
   -1763, -1763,   180,   447,  1223, -1763, -1763,  1223,  1223,  2027,
    2106,  1223, -1763,  2028, -1763,  2029,  1223,  2030, -1763,   987,
    2031,  2032, -1763, -1763, -1763, -1763,  2034,   125, -1763, -1763,
     447, -1763, -1763, -1763,  1223, -1763,  1223, -1763,  1984, -1763,
    1984,  2002,  1223,  1223, -1763,  1223,  2040, -1763,  1223,  1223,
   -1763,  1223,  1223, -1763,   447,  2041,  2042,  2008, -1763,   447,
     447,   447,   447,  2045, -1763, -1763,  2046,   993,   997, -1763,
   -1763, -1763, -1763, -1763, -1763,  2047,   447,  1223,  1223,  2049,
    1223,  1223, -1763, -1763,   532, -1763,   539, -1763, -1763,   447,
   -1763,  2050, -1763,  2051, -1763,   447,  2057, -1763, -1763, -1763,
   -1763, -1763,  1223, -1763,  1223,  1223,  1223,   447,   708, -1763,
    1001,  1223, -1763,  2058, -1763,  1005,  1223,  1009,  1028,  1223,
    1223,  1223,  1223, -1763, -1763,  1068,  2059,  2061, -1763, -1763,
    1098, -1763,  1704,  1223,  1223,  2062,  2063, -1763,  2064, -1763,
    1223,  1223,  2067, -1763, -1763, -1763,  1223, -1763, -1763,  2068,
   -1763, -1763, -1763,  1712, -1763,  1713,  2071, -1763, -1763,  2072,
   -1763,  2073, -1763, -1763,  1718,  2074,   -14,  1223,  1223,  2075,
    2076,  1223,  1127,  2077,  1134,  1223,  1223,  1223,  2078,  1223,
    1223,   447,  1223,   447,  2079,  1223,  1223, -1763,  2080, -1763,
   -1763, -1763,  2081, -1763, -1763,  2082, -1763, -1763, -1763, -1763,
     335, -1763, -1763, -1763, -1763, -1763,  2084,  2085,  1223,  2086,
   -1763,  2088,  1223,  2091,  1223,  2092, -1763, -1763,  1139,   447,
    2093,  2094,  2035, -1763, -1763, -1763,  2096, -1763,  2098, -1763,
    1223,  2102,  1223, -1763,  1223,  1144,   447, -1763,   447, -1763,
    2104, -1763, -1763,  1223, -1763,  1147,  1223,  1223,  1223,  1223,
    2105, -1763, -1763,  1223,  1223, -1763, -1763,  2108,  1223, -1763,
    2109, -1763,  1151,  1223,  1223,  1223,  1223, -1763,   378, -1763,
   -1763, -1763,   403, -1763,  2110,  1223,  1223, -1763, -1763, -1763,
    2111,  2112, -1763,  2113, -1763, -1763,  2114, -1763,  2115, -1763,
   -1763, -1763, -1763,  2116, -1763, -1763,   447,  1223,   348,  1223,
    1223, -1763, -1763,  2117, -1763,  1167, -1763, -1763,  2118,  1223,
    1223,  1187, -1763,  1223,  2119,   447,  2120,  1223, -1763,  2121,
    2122,  2107, -1763, -1763, -1763, -1763, -1763,  2123, -1763, -1763,
    2124, -1763,  2125, -1763, -1763,   447,  2126, -1763, -1763,  1216,
   -1763, -1763,  2127, -1763,  1223,  1223, -1763,  1223, -1763, -1763,
   -1763,  2128, -1763,  1223,  2130,  2131,  2132,   731, -1763,  1223,
     -25, -1763,  1223, -1763, -1763,  2133,  1223,  1223,  1223,  1223,
   -1763,   362,  1223, -1763,   367,  1223, -1763,  1223,  1223, -1763,
   -1763, -1763, -1763, -1763, -1763,  2136,    38, -1763,  1219,  1222,
   -1763, -1763,  1228, -1763,  1223,  1223, -1763,  2138,  1235, -1763,
    1223, -1763,  2139, -1763, -1763,  2140, -1763, -1763, -1763,   447,
   -1763, -1763,   447, -1763,  1223,  1223,  1252, -1763,  1223, -1763,
   -1763, -1763,  2141, -1763,  1223, -1763,  1223,   105,   111, -1763,
    1223,  1223,  1223,  1273, -1763,  1319, -1763,  1334,  1223,  1223,
   -1763, -1763,   447, -1763,  2142, -1763,  2144, -1763,  1337,  2145,
    1351, -1763, -1763,  2146,  2147, -1763,  2179,   447,   447,  1223,
    1223, -1763,  1223,  2148, -1763,  1223,  2151, -1763,  1223,  -157,
   -1763,  1223,   114,   447,  1223,  1223, -1763,  1223, -1763,  2153,
   -1763,  1352,  1223,  1223,   381, -1763, -1763, -1763,  1353, -1763,
   -1763,  1356, -1763, -1763,  2154,   447,   447,  1223,  2155,  1358,
   -1763,  1223, -1763,  2160, -1763,  1223,  2161, -1763,  1223,   -42,
    2163,  1223,  1223,  1223, -1763, -1763,  1223,  1223,  1382, -1763,
     397, -1763,  1223, -1763,  2164, -1763,  2166,   447,  1387, -1763,
   -1763,  2169,  1419, -1763,  2170, -1763,  2171, -1763,  1223, -1763,
    1223,  2172,  1223,  1427,  1428, -1763,  1223, -1763,  1223, -1763,
   -1763,  2173, -1763,  1223, -1763, -1763,  1459, -1763, -1763,  2174,
     447, -1763,  1223, -1763,  2175, -1763,  1223,  1223,  1223, -1763,
    1465, -1763,  2176, -1763,   447,   433, -1763,  1223,  2177,  1223,
   -1763,  2180, -1763,   447,  1223, -1763,  1223,  2181, -1763,  1223,
   -1763,  2183,  2186,   159, -1763,  1223, -1763, -1763, -1763,  1223,
    1223,  1223,  1223,  1223,  1223,  1223,  1223,  2187,  2189,  1223,
    1223, -1763, -1763,   271,  1223, -1763,  1223,  1223,  1223,  1223,
    1223,  1223,   447,  1223,  2191,  2192,  1223,  2193, -1763, -1763,
    1223, -1763,   447,   447,  1223,  1223,  1223,  2194, -1763
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1763, -1763, -1763,  2149, -1763, -1763, -1763, -1763, -1763,  2129,
   -1763, -1763,  2143, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,  1729,
    1839, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763,  2053, -1763,  1832, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1762, -1763,   -50,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,  1435,
   -1763, -1763,  1430, -1763, -1763, -1763, -1763, -1763, -1763,  1531,
   -1763,  -201, -1075, -1763, -1763,  1464, -1763,  2043, -1763,  -141,
   -1763,  2020, -1763,  -215,  -826,  -331, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763,   216, -1089,  2033, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763,  1506,  -898, -1763, -1763,  1526,  -880,
   -1763,  -389, -1763,  2022, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763,  1912, -1763, -1763, -1763, -1763,  2024,  1965,  1299,  -211,
   -1424,  -538,  -867, -1763, -1763,   -89, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763,  -361, -1763,  -702,  -353, -1763, -1763,
    1076, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763, -1763,  1157,  1655, -1763, -1763, -1763, -1763,
   -1763, -1763, -1763,  1014,  2134, -1763, -1763,  2631,    -6
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1076
static const yytype_int16 yytable[] =
{
     328,   329,   330,   333,   835,   878,   855,  1331,  1357,  1331,
    1383,  1357,  1367,   766,   785,   614,  1146,  1823,   363,  1052,
    2049,  1535,  1053,  1166,  1272,  1539,  1525,  1526,  1527,  1528,
     635,  1813,   363,   320,   320,  1428,   607,  1117,   324,  1930,
    1886,   324,   633,   317,   324,   554,  2158,   324,  1152,   588,
     608,   651,  1557,   363,   544,   320,   545,  2514,  1204,   364,
    1815,  1951,   320,   325,   320,   316,   325,   411,  1054,   325,
     365,   366,   325,   364,   589,   320,  1952,   324,   767,   408,
     363,   320,  1244,   320,   365,   366,  1153,   321,   321,   768,
     320,  1118,  1448,   546,   364,  1565,   856,  1481,   320,   324,
     537,   538,   325,   634,   320,   365,   366,   769,  1955,   321,
    1584,  1823,   320,   320,   470,  1823,   321,   319,   321,   476,
     477,   364,   386,  1956,   325,   326,  1119,   770,   326,   321,
    1887,   326,   365,   366,   326,   321,   327,   321,   491,   327,
     555,   500,   327,  1542,   321,   327,   505,  1585,   963,   964,
     965,   966,   321,  1219,   323,   601,  1566,   363,   321,  1888,
     363,  2515,  1167,   334,   326,   556,   321,   321,  1205,   590,
     424,   320,  2547,   860,   324,   327,   320,   558,   559,  1245,
     324,   771,   320,   324,   576,   591,   326,   772,   889,  2435,
    1220,   442,   773,   666,   667,  1055,  1931,   327,   364,   325,
    2305,   364,  1449,  2159,   320,   325,   320,  1482,   325,   365,
     366,  1814,   365,   366,  1776,  1429,  1430,  1431,   409,  1056,
     774,   320,  1057,   320,   592,   321,   775,  1058,   324,  1917,
     321,  1393,  1394,  1543,   557,   367,   321,  1059,   593,   528,
    1816,  1091,   613,   795,  1582,  1077,   594,   764,  1116,   856,
    1323,   318,  2451,   325,   971,   672,   782,   595,   321,   884,
     321,   326,  1544,  1889,   363,  1758,  1810,   326,   320,  1738,
     326,  1858,   327,   617,   776,   321,  2548,   321,   327,  1567,
    1681,   327,   674,   677,   320,   492,  2214,  1221,  2215,   337,
    1060,  1061,   765,  2436,  1273,  1274,  1275,  1276,  1277,   531,
    1278,  1279,  1280,  1281,  1282,   364,  1283,  1284,  1285,  1286,
    1034,  1144,   857,  1062,  1915,   326,   365,   366,  1168,  2477,
     547,   972,   321,  1347,   828,  2480,   327,   831,  2517,   320,
     338,   837,   838,  1811,  1823,   320,   320,   619,   321,  2209,
     324,   845,  1898,   650,  1682,   596,   847,  2306,  1120,  2307,
    1918,   850,   777,   341,   597,  1039,   973,   851,   852,  1347,
     854,  1777,  1818,   609,   342,   325,  1545,   679,  1147,  1121,
     796,   865,   345,  2608,   786,   548,   549,   550,   551,   552,
     778,  1880,   988,   321,  1063,   874,   598,   787,   320,   321,
     321,   324,  1369,  2149,  1826,   886,   888,   788,   320,  2452,
    1827,   789,  1089,   790,   791,   967,   324,   320,   792,   900,
     363,   324,   903,   363,  1369,   906,   325,   539,   540,   541,
     542,   543,  1575,  2478,   346,  2609,  2013,   326,   324,  2481,
     324,   325,  2518,   348,   797,   363,   325,   856,   327,  1115,
     856,  1609,   321,   320,  1043,  2610,  1906,  1921,   798,   799,
     800,   364,   321,   325,   364,   325,   320,   801,  1170,   802,
     803,   321,   365,   366,  1923,   365,   366,   947,   792,   804,
     805,   806,   807,   808,   954,   580,   364,   809,   326,   856,
     810,   811,  1774,   856,  1763,  2625,   856,   365,   366,   327,
    1359,   939,  1357,   326,  2191,  1207,   320,   321,   326,   856,
    1345,  2594,   324,  2038,   327,  1333,   324,  1331,  1844,   327,
     321,   351,   581,   856,  2150,   675,  2151,   675,   856,  1576,
     324,  1246,  1247,   992,   994,   995,   327,   325,   327,   352,
    1513,   325,   320,   324,   320,   324,  1514,  2626,  1045,   353,
    1588,   331,  2133,   320,  1922,   325,  2138,  1013,   856,  2334,
     321,  1370,  1018,  1019,  1022,   354,  1610,  2627,   325,   320,
     325,  1924,  2397,   355,   356,   795,  1035,  1047,  1036,   321,
    1038,   321,  1193,  1370,   324,  1385,  2444,  1589,  1248,   561,
    1064,  2446,  1067,  1069,  1071,  1217,   321,  1067,   321,   326,
     320,  1079,  2380,   326,   357,  2529,  1081,   321,  1083,   325,
     327,   940,  1595,   320,   327,   358,  1088,   326,  1321,   324,
    1693,  2557,  1095,   321,  1097,  1098,  1099,  2383,   327,   320,
     326,   582,   326,   359,   894,   583,   584,   585,   324,   898,
     324,   327,  1861,   327,   325,   324,  1590,  1108,  1862,  1596,
    1111,  1112,  1113,   360,   321,   324,  1249,  2595,  1122,  1251,
    1252,   332,   324,   325,   361,   325,   362,   321,  1250,  1130,
     325,   326,  1847,  1132,  2246,   410,   320,  1136,   324,   369,
     325,  2248,   327,   321,  1143,  1145,   324,   325,   844,   372,
     993,  1624,  1154,  1321,  1156,   945,   320,  1321,   320,   324,
     373,   324,   796,   325,  1694,   320,   326,  1175,  1597,  1177,
    1178,   325,   374,   958,   959,   375,  1253,   327,   562,   563,
     564,   565,   376,   324,   325,   326,   325,   326,  1625,  1068,
     321,   377,   326,  1210,  1211,  1212,   327,  1215,   327,   384,
     320,  1222,   326,   327,  1321,  1226,   385,   324,   325,   326,
     321,   387,   321,   327,  1235,  1236,  1237,  1321,  1321,   321,
     327,   389,  1243,   523,  1214,   326,   797,  1546,   320,  1254,
    1580,  1256,   325,   326,  1260,  1261,   327,  1189,   390,  1265,
     798,   799,   800,  1259,   327,  1324,   326,  1626,   326,   801,
    1337,   802,   803,   324,   321,   324,   391,   327,   850,   327,
    1355,   804,   805,   806,   807,   808,  1305,  1385,   324,   809,
     326,   324,   810,   811,   324,   320,   392,   393,   325,   850,
     325,   327,   321,  1434,  1581,   528,   335,  1325,  1326,   320,
     324,  1504,  1327,   325,   326,   324,   325,   324,  1332,   325,
     396,  1334,  1335,   399,  1560,   327,  1578,  1628,  1629,  1338,
     400,   401,   324,   402,  1342,   325,  1343,  1344,  1631,   324,
     325,   320,   325,  1630,   403,   324,   404,  1348,  1620,   321,
    1185,   320,  1353,  1632,   529,   324,   320,   325,  1356,   405,
     326,   406,   326,   321,   325,   531,   407,   349,   412,   324,
     325,   327,  1643,   327,   532,   326,   413,  1371,   326,   415,
     325,   326,   324,   418,  1375,  1376,   327,  1377,   419,   327,
    1379,   421,   327,  1386,   325,   321,   324,   326,   324,   324,
     324,   422,   326,   602,   326,   321,  2260,   325,   327,   370,
     321,   423,  2261,   327,   380,   327,   427,   320,  1672,   326,
    1674,   325,   320,   325,   325,   325,   326,   430,  1403,  2432,
     327,  1404,   326,  1677,   431,  2433,  1683,   327,   437,  1706,
    1410,   320,   326,   327,   438,  1413,   439,  1415,  1416,  1417,
    1418,  1419,  1420,   327,   443,  1740,   326,  1427,   603,   441,
    1769,  1191,  1798,   324,   449,   450,  1435,   327,   324,   326,
     453,   321,   604,   454,  1440,   382,   321,  1901,  1445,  1446,
     327,   320,   463,   326,  1904,   326,   326,   326,   325,   464,
    1909,   465,   324,   325,   327,   321,   327,   327,   327,   444,
    1966,   324,   466,   467,   324,  1421,  1422,  1423,  1424,  1767,
     468,   324,   469,  1473,  1985,  1778,   320,   325,   471,  1476,
    1477,  1479,  1480,   474,   324,   475,   325,  1987,   478,   325,
    1767,  1767,   324,   324,   324,   321,   325,   479,   480,   394,
     481,  1992,   482,  2024,  2075,  2102,   324,   483,  1495,   325,
     326,  1498,   324,  1499,   484,   326,   324,   325,   325,   325,
     324,   327,  1505,  1506,   324,  1508,   327,  1509,   324,   485,
     321,   325,  1511,   486,   397,   487,  1515,   325,  1517,   326,
    1519,   325,  1521,  1522,   320,   325,   488,   324,   326,   325,
     327,   326,  1533,   325,   489,  1536,  1537,  1538,   326,   327,
     490,  1541,   327,   493,  1549,  1550,  1551,   526,  2111,   327,
     496,   326,   325,  2113,  1561,  1562,  1563,  1564,   363,   326,
     326,   326,   327,   497,   503,  1573,  1574,   324,   504,  1579,
     327,   327,   327,   326,   506,   507,   320,  2119,   321,   326,
     508,   509,   416,   326,   327,   510,  2122,   326,  1605,  2124,
     327,   326,   325,  1611,   327,   326,  2134,   324,   327,   364,
     511,   512,   327,   320,   513,  1619,   327,   514,  1621,  2139,
     365,   366,  1627,   515,   326,   320,   446,  2155,  2162,  2173,
     516,   517,   325,  1641,  1899,   327,   324,  1644,  1907,   320,
     321,  2204,   518,   324,   425,  1651,  1652,  2235,   324,  1654,
     519,  2237,   520,   324,   521,  2262,   324,   577,   522,  2266,
     324,   325,   578,  2269,   326,   587,   644,   321,   325,  1673,
    1675,   447,  1678,   325,   780,   327,   324,  1684,   325,   321,
     645,   325,  2271,   428,  1690,   325,   649,   817,  1696,   822,
    1697,  1698,   823,   321,   326,   825,   324,   432,   826,  1707,
     824,   325,   827,   848,  1711,   327,   836,  1714,  1715,  1716,
    1717,  1718,  1719,  1720,  1721,  1722,  1723,  1724,  1725,  1726,
    1727,   325,  2277,   326,   839,   324,   840,  1731,   324,   842,
     326,   324,   324,   843,   327,   326,   846,   324,   853,   858,
     326,   327,   868,   326,   324,   869,   327,   326,   320,   870,
     325,   327,  2281,   325,   327,   871,   325,   325,   327,  1741,
    1742,   324,   325,   326,   872,   528,   813,   814,   816,   325,
    1748,   873,   942,   943,   327,   879,  1752,  1753,   320,   883,
    1754,  2314,   324,   326,   895,  1353,   325,  1757,  2317,   896,
     899,   902,  1760,  2344,   327,  1762,   905,   324,  2356,   910,
    1371,  2362,   321,  1765,   911,  2374,   434,   325,  1768,  1770,
    1771,  1772,   326,  1773,   529,   326,  1386,   530,   326,   326,
     915,  2401,   325,   327,   326,   531,   327,   916,   324,   327,
     327,   326,   321,   917,   532,   327,   455,  1788,   920,   320,
    1791,  2406,   327,   324,  1794,   320,   324,  1797,   326,  1799,
    1800,  1801,  1802,   325,   921,  1805,  1806,  1807,  1808,   327,
     324,   324,   324,  1812,   320,   324,   320,   324,   325,   326,
    2421,   325,   924,  2453,  1822,   927,  2455,  1824,   929,   930,
     327,   932,  2457,   933,   675,   325,   325,   325,   934,  2462,
     325,   324,   325,   321,   935,   327,   324,   459,   938,   321,
     955,   946,   951,   461,   952,   957,  2471,  1835,  1836,   960,
    1838,  1839,  1840,   975,   976,   326,   325,   979,   321,   987,
     321,   325,   472,   990,   494,   991,   327,  2486,   324,   997,
     326,   996,  1851,   326,   998,   999,   324,   324,   321,  1853,
    1854,   327,  1855,  1000,   327,  1857,  1001,   326,   326,   326,
    1002,  1003,   326,   325,   326,  1866,  1867,  1004,   327,   327,
     327,   325,   325,   327,   320,   327,  1005,  1006,   324,  1007,
    1877,  1878,  1879,  2488,   324,  1882,  1883,  1884,   326,  1885,
    1009,  1890,  1891,   326,  1893,  1894,  1402,  1409,  2490,   327,
    1902,  2497,  1905,   325,   327,  1910,  1911,  1912,  1913,   325,
     324,  1916, -1075,   324,   324,  2500,  2525,  2531,   324,  1010,
    2533,  1926,  2540,  1927,  1412,   326,   320,  1011,   321,  1012,
    1041,   320,   829,   326,   326,   325,   327, -1075,   325,   325,
     320,   324,  1042,   325,   327,   327,  2555,  1072,   320,   320,
    1049,  2562,  1050,  1073,  1051,   320,  1075,  1076,  1080,  1082,
     320,   320,  1084,   320,  1085,   326,   325,  1093,  1105,  1109,
    1148,   326,  1150,  1131,  1110,  1114,   327,   320,   320,   320,
     321,  1020,   327,  2565,   866,   321,   320,  1135,  1137,   912,
    1139,  2573,  2575,  1149,   321,  1172,   320,   326,   922, -1075,
     326,   326,   321,   321,  1157,   326,   925,   969,   327,   321,
   -1075,   327,   327,   977,   321,   321,   327,   321,  1086,  1161,
     320,  1180,  1967,  2581,  1151,   320,  2210,   320,   326,  2590,
     320,   321,   321,   321,  1158,  1183,  1200,  1329,  1163,   327,
     321,  1021,  1983,  1173,  1395,   320,  1176,  1179,  1986,  1182,
     321,  1988,   320, -1075,  1397,  1990,  1195,  1731,  1993,  1994,
    1995,  1996,  1997,  1998,  1999,  2000,  2001,  2002,  2003,  2004,
    2005,  2006,  1155,   320,   321,   320,  1187,  1202,  1529,   321,
    1233,   321,  1216,  1531,   321,  1571,   320,   320,  1592,   320,
    2017,  2018,  2019,  1229,  1241,  1242,  2021,  2022,  1255,   321,
    1264,  2023,  1267,  1601,  2025,  2027,   321,  1268,  2029,  2030,
    1699,  1269,  2031,   320,   320,   320,   320,  1271,  1295,  2036,
    1296,  1298,  1299,   320,  1301,  1302,  1315,   321,  1306,   321,
    2052,  1703,   320,  1708,  2056,   320,  2058,  1307,   320,  2061,
     321,   321,  1308,   321,  1728,  1780,  1309,  1820,   320,  2069,
    2070,  2071,  2072,  2073,   320,   320,  2076,   320,  2077,  1312,
    2078,  1316,  1310,  1311,  1313,  2308,  1314,   321,   321,   321,
     321,  1828,  1868,  1870,  1871,  1317,  1318,   321,   320,  2088,
    2089,  1872,  2090,  2091,  1319,   320,   321,   320,  1320,   321,
    1873,  1328,   321,  1895,  1336,  2098,  1897,  1340,  1341,  1349,
    2103,  1350,   321,  2106,  1354,  1362,  1903,   320,   321,   321,
     320,   321,  1942,  1944,  1363,  1946,  1378,  1387,   320,   320,
    2112,  2114,  2115,  2116,   320,  2118,  2120,  2121,  2123,  2125,
    2126,  2127,   321,  2128,  2129,  2130,  1981,  2131,  2132,   321,
    1389,   321,  1390,  2079,  2135,  2141,  2136,  1399,  1400,  2137,
    1407,  1414,  2140,  1425,  1426,  2143,  2144,  2145,  1432,  2146,
    1444,   321,  2148,  1450,   321,  2164,  1451,  2381,  2283,  1809,
    2156,  2384,   321,   321,  1452,  1453,  2295,  2297,   321,  1454,
    1455,  1456,  2302,  1457,  1458,  1459,  1460,  1461,  2039,  2163,
    2040,  2041,  2042,  2043,  2044,  2045,  2046,  2047,  1462,  1463,
    1464,  1467,  1468,  1469,  1470,  1471,  1472,  1483,  1494,  1496,
    1497,  1500,  1501,  1502,  1503,  1507,  1510,  1518,  1520,  1534,
    1569,  1570,  1587,  1591,  1586,  1594,  1598,  1599,  1600,  1603,
    2174,  1604,  2175,  1607,  2176,  1608,  1612,  2177,  2178,  2179,
    2180,  2181,  2182,  2183,  2184,  2185,  2186,  2187,  2188,  2189,
    2190,  1613,  1614,  1615,  1616,  1617,  1618,  1622,  2194,  1623,
    1634,  2195,  2196,  1635,  1636,  2199,  1637,  1638,  1639,  1642,
    2202,  1645,  1646,  2205,  1647,  1648,  1649,  1650,  1655,  1656,
    1657,  1658,  1659,  1660,  1661,  1662,  1663,  1664,  2212,  1665,
    2213,  1666,  1667,  1668,  1669,  1670,  2217,  2218,  1671,  2219,
    1676,  1679,  2221,  2222,  1680,  2223,  2224,  1685,  1701,  1702,
    1705,  1710,  1712,  1730,  1732,  1733,  1735,  1743,  1744,  1745,
    1746,  2236,  2238,  1747,  1749,  1750,  1755,  1759,  1766,  1779,
    1784,  2241,  2242,  1786,  2244,  2245,  1789,  1792,  1795,  1803,
    1817,  1825,  1830,  1831,  1832,  1837,  1785,  1841,  1842,  1843,
    1846,  1849,  1850,  1852,  1860,  1864,  2255,  1865,  2256,  2257,
    2258,  1874,  1875,  1892,  2263,  2264,  1914,  1920,  1925,  2267,
    2268,  2270,  2272,  2273,  2274,  2275,  2276,  1928,  1932,  2278,
    1933,  1934,  1935,  2530,  2282,  1936,  1937,  2285,  2286,  1938,
    1939,  1940,  1941,  1948,  2290,  2291,  1949,  1950,  1957,  1958,
    2293,  1959,  1960,  1961,  1962,  1963,  1964,  1965,  1968,  1969,
    1970,  1971,  1972,  1973,  1974,  1975,  1976,  1977,  1978,  1979,
    1980,  2309,  2310,  2008,  1984,  2313,  2315,  1989,  2318,  2319,
    2320,  2321,  1991,  2323,  2324,  2007,  2326,  2010,  2015,  2329,
    2330,  2016,  2020,  2028,  2032,  2033,  2034,  2035,  1780,  2051,
    2053,  2054,  2057,  2059,  2060,  2062,  2063,  2068,  2074,  2080,
    2198,  2415,  2337,  2081,  2082,  2083,  2340,  2084,  2342,  2092,
    2100,  2104,  2345,  2107,  2108,  2109,  2110,  2147,  2153,  2154,
    2157,  2160,  2161,  2167,  2352,  2168,  2354,  2170,  2355,  2357,
    2172,  2197,  2200,  2201,  2203,  2206,  2207,  2361,  2208,  2363,
    2364,  2365,  2366,  2367,  2220,  2226,  2227,  2369,  2370,  2233,
    2234,  2239,  2372,  2243,  2251,  2252,  2375,  2376,  2377,  2378,
    2379,  2254,  2265,  2279,  2216,  2280,  2287,  2288,  2289,  2387,
    2388,  2292,  2294,  2504,  2228,  2299,  2300,  2301,  2304,  2311,
    2312,  2316,  2322,  2328,  2331,  2332,  2333,  1030,  2335,  2336,
    2338,  2396,  2339,  2398,  2399,  2341,  2343,  2347,  2348,  2402,
    2350,  2349,  2351,  2404,  2405,  2407,  2353,  2408,  2360,  2368,
     525,  2412,  2371,  2373,  2386,  2389,  2390,  2391,  2392,  2393,
    2394,  2400,  2403,  2409,  2411,  2413,  2414,  2416,  2417,  2418,
    2420,  2423,  2427,  2422,  2429,  2430,  2431,  2439,  2424,  2425,
    2450,  2426,  2461,  2465,  2466,  2474,  2495,  2428,  2496,  2499,
    2502,  2503,  2510,  2434,  2437,  2512,  2438,  2524,  2535,  2539,
    2440,  2441,  2442,  2443,  2543,  2545,  2445,  2549,  2559,  2447,
    2560,  2448,  2449,  2564,  2567,  2568,  2571,  2579,  2583,  2586,
    2592,  2598,  2454,  2456,  2600,  2604,  2458,  2606,  2459,  2460,
    2607,  2621,  2463,  2622,  2464,  2638,  2639,  2641,  2648,  1033,
     849,  1360,  1405,  1380,  1436,  1433,   859,  1364,  2469,  2470,
    2472,   983,  2473,  1554,   890,  1074,   762,  1756,  2475,  1713,
    2476,  2479,  2482,     0,  2483,  2484,  2485,  2487,     0,  2489,
     763,  2491,  2492,  2493,   953,   897,  1294,   815,     0,     0,
       0,     0,  2498,   914,  2501,     0,     0,     0,     0,     0,
       0,     0,   918,  2507,  2508,     0,  2509,     0,     0,  2511,
       0,     0,  2513,     0,     0,  2516,  2519,     0,  2521,  2522,
       0,  2523,     0,     0,     0,  2526,  2527,  2528,     0,     0,
       0,     0,  2532,     0,     0,  2534,     0,     0,     0,     0,
       0,  2538,     0,  2541,     0,  2542,     0,     0,     0,  2544,
       0,     0,  2546,     0,     0,  2550,  2551,  2552,     0,     0,
    2553,  2554,  2556,     0,     0,     0,  2558,     0,     0,     0,
       0,     0,  2563,     0,     0,     0,  2566,     0,     0,     0,
       0,     0,  2569,     0,  2570,     0,  2572,  2574,  2576,     0,
    2577,     0,  2578,     0,     0,     0,     0,  2580,     0,     0,
    2582,     0,     0,     0,     0,     0,  2585,     0,     0,     0,
    2587,  2588,  2589,     0,  2591,     0,     0,     0,     0,  2596,
       0,  2597,     0,  2599,     0,     0,     0,     0,  2602,     0,
    2603,     0,     0,  2605,     0,     0,     0,  2611,     0,  2612,
       0,     0,     0,  2613,  2614,  2615,  2616,  2617,  2618,  2619,
    2620,     0,     0,  2623,  2624,     0,     0,  2628,  2629,     0,
    2630,  2631,  2632,  2633,  2634,  2635,     0,  2637,     0,     0,
    2640,     0,     0,     0,  2642,   322,     0,     0,  2645,  2646,
    2647,     0,   336,     0,     0,   339,   340,     0,     0,   343,
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
       0,     0,     0,     0,     0,   612,     0,   615,   616,     0,
       0,   618,   620,   621,   622,     0,   624,     0,   626,   628,
     630,   632,     0,     0,   624,   628,   632,     0,     0,   637,
     639,     0,     0,   641,   643,     0,     0,   646,   647,   648,
       0,     0,     0,   652,     0,     0,     0,     0,     0,     0,
       0,   654,   656,   657,     0,     0,   658,   659,   661,   662,
     663,     0,   665,   572,   572,   669,   670,   671,   673,     0,
     676,   678,   680,     0,   681,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   781,   783,   784,
       0,     0,     0,     0,     0,     0,     0,   821,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   830,     0,     0,   834,   834,     0,     0,     0,
       0,     0,   841,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   626,   624,   863,     0,     0,   867,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   877,   877,     0,   882,     0,     0,
       0,   628,     0,     0,   632,   630,     0,     0,   893,   821,
       0,     0,     0,   863,   821,     0,     0,   901,     0,     0,
     904,     0,     0,     0,   909,     0,     0,     0,     0,     0,
     913,     0,   654,     0,     0,     0,     0,     0,     0,   665,
     919,     0,     0,     0,   923,     0,     0,   926,     0,     0,
       0,   928,     0,     0,   931,     0,     0,     0,     0,     0,
       0,     0,   937,     0,   941,     0,     0,     0,     0,     0,
     821,     0,     0,     0,     0,   950,     0,     0,     0,     0,
     656,     0,     0,     0,     0,     0,     0,     0,   821,   821,
       0,     0,     0,   962,     0,     0,   970,     0,     0,     0,
     974,     0,     0,   978,     0,     0,     0,     0,     0,   982,
     661,   986,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   989,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1008,
       0,     0,     0,     0,     0,  1014,  1015,  1016,  1017,     0,
       0,     0,  1023,  1024,  1025,  1026,  1027,     0,  1029,  1029,
    1032,  1032,     0,     0,     0,     0,  1037,     0,  1040,     0,
       0,  1044,  1046,  1048,     0,     0,     0,     0,  1065,     0,
       0,     0,     0,     0,     0,     0,     0,  1078,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1087,     0,     0,  1090,     0,  1092,     0,  1094,     0,
    1096,     0,     0,     0,     0,  1100,     0,  1101,     0,  1102,
       0,  1103,     0,  1104,     0,     0,     0,     0,     0,     0,
    1106,     0,  1107,     0,     0,     0,     0,     0,     0,     0,
       0,  1090,     0,  1092,     0,     0,     0,  1123,  1124,  1125,
    1126,     0,  1127,  1128,  1129,     0,     0,     0,     0,     0,
       0,  1133,  1134,     0,     0,     0,     0,  1138,     0,  1140,
    1141,  1142,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1160,  1162,     0,  1164,  1165,  1169,
    1171,     0,     0,  1174,     0,     0,     0,     0,     0,  1181,
       0,  1184,  1186,     0,  1188,  1190,  1192,  1194,     0,  1196,
    1197,  1198,  1199,  1201,     0,  1203,     0,  1206,  1208,  1209,
       0,     0,     0,  1213,     0,     0,  1218,     0,     0,  1223,
    1224,  1225,     0,  1227,  1228,     0,  1230,  1231,  1232,     0,
    1234,     0,     0,     0,  1238,  1239,  1240,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1257,
    1258,     0,     0,  1262,  1263,     0,     0,  1266,     0,     0,
       0,     0,  1270,     0,     0,  1287,  1288,  1289,     0,  1291,
    1292,     0,     0,     0,     0,     0,     0,     0,  1300,     0,
       0,  1303,  1304,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1322,     0,  1092,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1330,     0,   834,     0,   834,   834,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1339,     0,
       0,     0,     0,     0,     0,     0,     0,  1346,     0,     0,
       0,     0,   863,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   877,     0,  1358,   877,
     877,   882,     0,  1361,     0,     0,     0,     0,     0,     0,
       0,     0,  1366,  1368,     0,  1322,     0,     0,   863,  1322,
       0,     0,     0,     0,     0,     0,     0,     0,  1382,  1384,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1388,     0,     0,     0,     0,     0,     0,     0,     0,  1391,
       0,     0,  1392,     0,     0,     0,     0,  1396,     0,     0,
    1398,     0,     0,     0,  1401,     0,  1322,     0,     0,   950,
       0,  1406,     0,     0,     0,     0,  1408,     0,     0,  1322,
    1322,  1411,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   982,     0,     0,     0,   986,     0,  1437,     0,  1438,
    1439,     0,     0,     0,     0,     0,     0,  1447,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1465,  1466,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1474,     0,  1475,     0,     0,  1478,     0,
       0,     0,  1484,     0,  1485,  1486,  1487,  1488,  1489,  1490,
       0,     0,     0,     0,     0,  1491,     0,  1492,     0,  1493,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1092,  1512,     0,  1516,     0,     0,     0,     0,     0,
       0,  1524,  1524,  1524,  1524,  1524,     0,  1530,  1532,     0,
       0,     0,     0,     0,     0,     0,     0,  1540,     0,  1547,
    1548,     0,     0,     0,  1553,  1553,  1555,  1556,  1524,  1558,
    1559,     0,     0,     0,     0,  1568,     0,     0,     0,     0,
       0,  1572,     0,     0,     0,  1577,     0,     0,  1583,     0,
       0,     0,     0,     0,  1593,     0,     0,     0,     0,     0,
    1602,     0,     0,     0,     0,     0,  1606,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1633,     0,     0,     0,     0,     0,     0,     0,  1640,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1653,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1686,  1687,  1688,
    1689,     0,  1691,  1692,  1695,     0,     0,     0,     0,  1700,
       0,     0,     0,  1704,     0,     0,     0,  1709,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1729,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1734,     0,     0,     0,     0,
       0,     0,     0,  1736,     0,     0,     0,     0,     0,     0,
    1737,     0,     0,  1092,  1739,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   834,     0,     0,     0,     0,     0,
       0,  1751,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     877,     0,  1761,     0,     0,  1368,     0,  1764,     0,  1764,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1384,     0,  1775,     0,  1775,     0,     0,     0,     0,
       0,     0,  1782,  1783,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1787,     0,     0,     0,  1790,     0,     0,
    1793,     0,     0,  1796,     0,     0,     0,     0,     0,     0,
       0,  1804,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1524,  1819,
    1821,     0,     0,     0,     0,     0,     0,     0,  1829,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1834,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1845,     0,  1848,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1856,     0,     0,  1859,     0,     0,     0,  1863,     0,     0,
       0,     0,     0,     0,  1869,     0,  1869,  1869,  1869,  1869,
       0,     0,     0,     0,  1876,     0,     0,     0,     0,     0,
       0,  1881,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1896,     0,  1896,  1900,     0,  1869,     0,
    1908,     0,     0,     0,     0,     0,     0,     0,     0,  1919,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1929,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1943,     0,     0,     0,     0,     0,     0,     0,
    1945,     0,  1947,     0,     0,     0,     0,     0,     0,     0,
    1953,     0,  1954,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   682,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   683,     0,     0,     0,     0,   684,     0,
       0,     0,   685,     0,     0,     0,     0,   686,     0,     0,
       0,     0,     0,     0,     0,   687,     0,     0,   688,   689,
       0,   690,     0,   691,   692,     0,     0,   693,  1982,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   694,     0,
       0,     0,     0,     0,     0,   695,   696,   697,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   698,     0,  2009,     0,  2011,     0,  2012,
    2014,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   699,   700,   701,   702,     0,     0,     0,
       0,   703,  2026,     0,   704,   705,     0,     0,   706,   707,
       0,     0,   708,     0,     0,     0,     0,     0,  2037,     0,
       0,     0,     0,     0,  2050,     0,     0,     0,     0,     0,
       0,  2055,   709,     0,   710,     0,     0,     0,     0,     0,
    2064,  2065,  2066,  2067,     0,   711,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   712,   713,  1869,
       0,     0,     0,     0,     0,     0,     0,   714,     0,   715,
    2085,   716,   717,     0,  2086,  2087,   718,   719,   720,   721,
       0,     0,     0,     0,     0,  2093,     0,  2094,  2095,     0,
    2096,  2097,   722,   723,  2099,     0,  2101,     0,   724,  2105,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2117,     0,     0,   725,   726,     0,   727,     0,   728,
       0,   729,     0,     0,     0,     0,   730,   731,   732,   733,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2142,   734,     0,     0,     0,     0,     0,     0,   735,
    2152,     0,     0,     0,   736,     0,   737,   738,   739,   740,
     741,   742,   743,   744,   745,   746,   747,     0,     0,     0,
       0,   748,   749,     0,   750,     0,     0,     0,  2165,     0,
       0,     0,  2166,     0,     0,     0,  2169,     0,   751,     0,
       0,     0,     0,     0,     0,     0,     0,   752,   753,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   754,  2171,     0,   755,   756,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   757,
       0,     0,   758,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2192,  2193,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   759,
       0,  2211,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2225,     0,     0,     0,     0,
    2229,  2230,  2231,  2232,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2240,     0,     0,
       0,     0,     0,     0,     0,  2247,     0,  2249,     0,     0,
    2250,     0,     0,     0,     0,     0,  2253,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2259,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2284,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2296,     0,  2298,     0,     0,     0,
       0,     0,     0,     0,     0,  2303,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2325,     0,  2327,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2346,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2358,     0,  2359,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2382,
       0,     0,     0,  2385,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  2395,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2410,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2419,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2467,     0,     0,  2468,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  2494,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2505,  2506,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2520,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2536,  2537,     0,     0,
       0,     0,     0,     1,     2,     3,     0,     0,     4,     0,
       5,     0,     0,     0,     0,     6,     0,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,  2561,    18,
      19,     0,     0,     0,    20,    21,    22,     0,     0,     0,
       0,     0,     0,     0,    23,     0,    24,    25,    26,    27,
       0,     0,     0,     0,     0,     0,     0,    28,     0,     0,
       0,  2584,     0,    29,    30,     0,     0,    31,     0,     0,
      32,    33,    34,     0,     0,  2593,    35,     0,    36,     0,
       0,     0,     0,    37,  2601,     0,    38,    39,    40,    41,
     524,    42,    43,    44,     0,     0,     0,    45,     0,     0,
       0,    46,    47,     0,     0,    48,    49,     0,    50,    51,
      52,     0,     0,     0,     0,     0,     0,    53,     0,     0,
       0,     0,    54,  2636,     0,    55,     0,    56,    57,     0,
       0,    58,     0,  2643,  2644,    59,     0,     0,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    70,    71,
       0,     0,     0,     0,    72,    73,     0,     0,    74,     0,
      75,    76,    77,    78,    79,     0,    80,     0,    81,     0,
       0,    82,     0,     0,     0,     0,    83,     0,    84,    85,
      86,    87,     0,     0,    88,    89,    90,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      91,     0,     0,     0,     0,     0,     0,    92,     0,     0,
       0,     0,    93,     0,     0,    94,     0,     0,     0,    95,
      96,    97,     0,     0,     0,    98,    99,   100,   101,     0,
     102,     0,   103,   104,   105,     0,   106,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   107,   108,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   109,
     110,     0,     0,     0,     0,   111,     0,   112,     0,   113,
       0,   114,     0,   115,   116,     0,     0,     0,   117,     0,
       0,     0,     0,     0,     0,     0,     0,   118,   119,     0,
     120,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   121,     0,   122,     0,     0,     0,   123,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     124,     0,   125,     0,   126,     0,   127,   128,   129,   130,
     131,   132,     0,     0,   133,     0,     0,     0,     0,     0,
       0,   134,     0,     0,   135,   136,   137,     0,     0,     0,
       0,     0,   138,     0,   139,     0,     0,     0,     0,     0,
     140,   141,     0,     0,     0,     0,     0,     0,   142,     0,
       0,     0,   143,   144,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     145,   146,   147,   148,     0,     0,     0,   149,     0,     0,
       0,     0,     0,   150,     0,   151,   152,   153,   154,   155,
     156,   157,   158,     0,     0,     0,     0,     0,   159,   160,
       0,   161,   162,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     163,   164,   165,     0,   166,   167,     0,     0,   168,     1,
       2,     3,   169,     0,     4,     0,     5,     0,     0,     0,
       0,     6,     0,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,     0,    18,    19,     0,     0,     0,
      20,    21,    22,     0,     0,     0,     0,     0,     0,     0,
      23,     0,    24,    25,    26,    27,     0,     0,     0,     0,
       0,     0,     0,    28,     0,     0,     0,     0,     0,    29,
      30,     0,     0,    31,     0,     0,    32,    33,    34,     0,
       0,     0,    35,     0,    36,     0,     0,     0,     0,    37,
       0,     0,    38,    39,    40,    41,     0,    42,    43,    44,
       0,     0,     0,    45,     0,     0,     0,    46,    47,     0,
       0,    48,    49,     0,    50,    51,    52,     0,     0,     0,
       0,     0,     0,    53,     0,     0,     0,     0,    54,     0,
       0,    55,     0,    56,    57,     0,     0,    58,     0,     0,
       0,    59,     0,     0,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    70,    71,     0,     0,     0,     0,
      72,    73,     0,     0,    74,     0,    75,    76,    77,    78,
      79,     0,    80,     0,    81,     0,     0,    82,     0,     0,
       0,     0,    83,     0,    84,    85,    86,    87,     0,     0,
      88,    89,    90,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    91,     0,     0,     0,
       0,     0,     0,    92,     0,     0,     0,     0,    93,     0,
       0,    94,     0,     0,     0,    95,    96,    97,     0,     0,
       0,    98,    99,   100,   101,     0,   102,     0,   103,   104,
     105,     0,   106,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   107,   108,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   109,   110,     0,     0,     0,
       0,   111,     0,   112,     0,   113,     0,   114,     0,   115,
     116,     0,     0,     0,   117,     0,     0,     0,     0,     0,
       0,     0,     0,   118,   119,     0,   120,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   121,     0,
     122,     0,     0,     0,   123,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   124,     0,   125,     0,
     126,     0,   127,   128,   129,   130,   131,   132,     0,     0,
     133,     0,     0,     0,     0,     0,     0,   134,     0,     0,
     135,   136,   137,     0,     0,     0,     0,     0,   138,     0,
     139,     0,     0,     0,     0,     0,   140,   141,     0,     0,
       0,     0,     0,     0,   142,     0,     0,     0,   143,   144,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   145,   146,   147,   148,
       0,     0,     0,   149,     0,     0,     0,     0,     0,   150,
       0,   151,   152,   153,   154,   155,   156,   157,   158,     0,
       0,     0,     0,     0,   159,   160,     0,   161,   162,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   163,   164,   165,     0,
     166,   167,     0,     0,   168,     0,     0,     0,   169
};

static const yytype_int16 yycheck[] =
{
       6,     7,     8,     9,   335,   394,   367,   833,   875,   835,
     908,   878,   892,     9,    40,   226,    13,  1441,    32,    11,
    1782,  1110,    14,    75,   104,  1114,  1101,  1102,  1103,  1104,
     245,   214,    32,   156,   156,   214,    31,     7,    69,   214,
      95,    69,   243,     6,    69,    16,   214,    69,    13,    16,
      45,   262,  1127,    32,    77,   156,    79,   214,   104,    73,
     214,   199,   156,    94,   156,   214,    94,    73,    60,    94,
      84,    85,    94,    73,    41,   156,   214,    69,    74,   214,
      32,   156,   214,   156,    84,    85,    51,   210,   210,    85,
     156,    61,   104,   116,    73,   153,   151,   104,   156,    69,
      38,    39,    94,   244,   156,    84,    85,   103,   199,   210,
     214,  1535,   156,   156,   120,  1539,   210,   214,   210,   125,
     126,    73,   311,   214,    94,   156,    96,   123,   156,   210,
     185,   156,    84,    85,   156,   210,   167,   210,   214,   167,
     111,   147,   167,    95,   210,   167,   152,   251,   112,   113,
     114,   115,   210,   214,   214,   205,   214,    32,   210,   214,
      32,   318,   214,   214,   156,   136,   210,   210,   214,   136,
     104,   156,   214,   374,    69,   167,   156,   183,   184,   311,
      69,   177,   156,    69,   190,   152,   156,   183,   403,   214,
     251,   104,   188,   282,   283,   187,   371,   167,    73,    94,
     214,    73,   214,   371,   156,    94,   156,   214,    94,    84,
      85,   394,    84,    85,   214,   394,   395,   396,   353,   211,
     216,   156,   214,   156,   191,   210,   222,   219,    69,   214,
     210,   933,   934,   185,   205,   214,   210,   229,   205,   206,
     394,   342,   317,   209,   366,   311,   213,   297,   342,   151,
     342,   214,   214,    94,   214,   336,   379,   224,   210,   400,
     210,   156,   214,   318,    32,  1354,   214,   156,   156,   342,
     156,   314,   167,   317,   270,   210,   318,   210,   167,   337,
     150,   167,   288,   289,   156,   361,  2048,   348,  2050,   214,
     282,   283,   298,   318,   374,   375,   376,   377,   378,   266,
     380,   381,   382,   383,   384,    73,   386,   387,   388,   389,
     338,   342,   214,   305,   336,   156,    84,    85,   370,   214,
     343,   281,   210,   861,   330,   214,   167,   333,   214,   156,
     214,   337,   338,   281,  1758,   156,   156,   317,   210,   214,
      69,   347,   214,   317,   214,   312,   352,   361,   318,   363,
     335,   357,   348,   214,   321,   104,   316,   363,   364,   897,
     366,   361,  1437,   358,   214,    94,   318,   317,   365,   339,
     336,   377,   214,   214,   400,   398,   399,   400,   401,   402,
     376,   314,   317,   210,   376,   391,   353,   413,   156,   210,
     210,    69,    55,   214,   208,   401,   402,   423,   156,   361,
     214,   427,   613,   429,   430,   369,    69,   156,   434,   415,
      32,    69,   418,    32,    55,   421,    94,   355,   356,   357,
     358,   359,   214,   318,   214,   266,   314,   156,    69,   318,
      69,    94,   318,   214,   400,    32,    94,   151,   167,   650,
     151,    99,   210,   156,   104,   286,   214,   214,   414,   415,
     416,    73,   210,    94,    73,    94,   156,   423,   285,   425,
     426,   210,    84,    85,   214,    84,    85,   473,   434,   435,
     436,   437,   438,   439,   480,   303,    73,   443,   156,   151,
     446,   447,  1380,   151,  1364,   214,   151,    84,    85,   167,
     879,   104,  1359,   156,   314,   253,   156,   210,   156,   151,
     214,    68,    69,   214,   167,   836,    69,  1333,   132,   167,
     210,   214,   340,   151,   335,   156,   337,   156,   151,   311,
      69,   188,   189,   529,   530,   531,   167,    94,   167,   214,
     208,    94,   156,    69,   156,    69,   214,   266,   104,   214,
     214,   104,   214,   156,   311,    94,   214,   553,   151,   214,
     210,   214,   558,   559,   560,   214,   214,   286,    94,   156,
      94,   311,   214,   214,   214,   209,   572,   104,   574,   210,
     576,   210,   285,   214,    69,   214,   214,   251,   245,    43,
     586,   214,   588,   589,   590,   285,   210,   593,   210,   156,
     156,   597,   214,   156,   214,   214,   602,   210,   604,    94,
     167,   214,   214,   156,   167,   214,   612,   156,   819,    69,
     100,   214,   618,   210,   620,   621,   622,   214,   167,   156,
     156,   449,   156,   214,   408,   453,   454,   455,    69,   413,
      69,   167,   208,   167,    94,    69,   310,   643,   214,   251,
     646,   647,   648,   214,   210,    69,   313,   214,   654,   188,
     189,   214,    69,    94,   214,    94,   214,   210,   325,   665,
      94,   156,   132,   669,   132,   214,   156,   673,    69,   214,
      94,   132,   167,   210,   680,   681,    69,    94,   214,   214,
     214,   214,   688,   894,   690,   469,   156,   898,   156,    69,
     214,    69,   336,    94,   184,   156,   156,   703,   310,   705,
     706,    94,   214,   487,   488,   214,   245,   167,   172,   173,
     174,   175,   214,    69,    94,   156,    94,   156,   251,   214,
     210,   214,   156,   729,   730,   731,   167,   733,   167,   214,
     156,   737,   156,   167,   945,   741,   214,    69,    94,   156,
     210,   214,   210,   167,   750,   751,   752,   958,   959,   210,
     167,   214,   758,     0,   214,   156,   400,  1118,   156,   765,
     214,   767,    94,   156,   770,   771,   167,   193,   214,   775,
     414,   415,   416,   214,   167,   214,   156,   310,   156,   423,
     214,   425,   426,    69,   210,    69,   214,   167,   794,   167,
     214,   435,   436,   437,   438,   439,   802,   214,    69,   443,
     156,    69,   446,   447,    69,   156,   214,   214,    94,   815,
      94,   167,   210,   214,   268,   206,   214,   823,   824,   156,
      69,   214,   828,    94,   156,    69,    94,    69,   834,    94,
     214,   837,   838,   214,   214,   167,   214,   198,   199,   845,
     214,   214,    69,   214,   850,    94,   852,   853,   199,    69,
      94,   156,    94,   214,   214,    69,   214,   863,   214,   210,
     165,   156,   868,   214,   255,    69,   156,    94,   874,   214,
     156,   214,   156,   210,    94,   266,   214,   214,   214,    69,
      94,   167,   214,   167,   275,   156,   214,   893,   156,   214,
      94,   156,    69,   214,   900,   901,   167,   903,   214,   167,
     906,   214,   167,   909,    94,   210,    69,   156,    69,    69,
      69,   214,   156,   136,   156,   210,   208,    94,   167,   214,
     210,   214,   214,   167,   214,   167,   214,   156,   214,   156,
     214,    94,   156,    94,    94,    94,   156,   214,   944,   208,
     167,   947,   156,   214,   214,   214,   214,   167,   214,   214,
     956,   156,   156,   167,   214,   961,   214,   963,   964,   965,
     966,   967,   968,   167,   169,   214,   156,   973,   191,   214,
     214,   195,   214,    69,   214,   214,   982,   167,    69,   156,
     214,   210,   205,   214,   990,   214,   210,   214,   994,   995,
     167,   156,   104,   156,   214,   156,   156,   156,    94,   214,
     214,   214,    69,    94,   167,   210,   167,   167,   167,   214,
     214,    69,   214,   214,    69,   112,   113,   114,   115,  1372,
     214,    69,   214,  1029,   214,  1386,   156,    94,   214,  1035,
    1036,  1037,  1038,   214,    69,   214,    94,   214,   214,    94,
    1393,  1394,    69,    69,    69,   210,    94,   214,   214,   214,
     214,   214,   214,   214,   214,   214,    69,   214,  1064,    94,
     156,  1067,    69,  1069,   214,   156,    69,    94,    94,    94,
      69,   167,  1078,  1079,    69,  1081,   167,  1083,    69,   214,
     210,    94,  1088,   104,   214,   214,  1092,    94,  1094,   156,
    1096,    94,  1098,  1099,   156,    94,   214,    69,   156,    94,
     167,   156,  1108,    94,   214,  1111,  1112,  1113,   156,   167,
     214,  1117,   167,   214,  1120,  1121,  1122,   132,   214,   167,
     214,   156,    94,   214,  1130,  1131,  1132,  1133,    32,   156,
     156,   156,   167,   214,   214,  1141,  1142,    69,   214,  1145,
     167,   167,   167,   156,   214,   214,   156,   214,   210,   156,
     214,   214,   214,   156,   167,   214,   214,   156,  1164,   214,
     167,   156,    94,  1169,   167,   156,   214,    69,   167,    73,
     214,   214,   167,   156,   214,  1181,   167,   214,  1184,   214,
      84,    85,  1188,   214,   156,   156,   169,   214,   214,   214,
     214,   214,    94,  1199,  1555,   167,    69,  1203,  1559,   156,
     210,   214,   214,    69,   214,  1211,  1212,   214,    69,  1215,
     214,   214,   214,    69,   214,   214,    69,   104,   311,   214,
      69,    94,   450,   214,   156,   188,   447,   210,    94,  1235,
    1236,   214,  1238,    94,   214,   167,    69,  1243,    94,   210,
     191,    94,   214,   214,  1250,    94,   191,   420,  1254,   214,
    1256,  1257,     6,   210,   156,   214,    69,   214,   214,  1265,
       6,    94,   214,   104,  1270,   167,   214,  1273,  1274,  1275,
    1276,  1277,  1278,  1279,  1280,  1281,  1282,  1283,  1284,  1285,
    1286,    94,   214,   156,   214,    69,   214,  1293,    69,   214,
     156,    69,    69,   214,   167,   156,   214,    69,    85,   214,
     156,   167,   214,   156,    69,   214,   167,   156,   156,   191,
      94,   167,   214,    94,   167,   214,    94,    94,   167,  1325,
    1326,    69,    94,   156,   214,   206,   312,   313,   314,    94,
    1336,   214,   272,    52,   167,   214,  1342,  1343,   156,   214,
    1346,   214,    69,   156,   214,  1351,    94,  1353,   214,   214,
     214,   214,  1358,   214,   167,  1361,   214,    69,   214,   214,
    1366,   214,   210,  1369,   214,   214,   214,    94,  1374,  1375,
    1376,  1377,   156,  1379,   255,   156,  1382,   258,   156,   156,
     214,   214,    94,   167,   156,   266,   167,   214,    69,   167,
     167,   156,   210,   214,   275,   167,   214,  1403,   214,   156,
    1406,   214,   167,    69,  1410,   156,    69,  1413,   156,  1415,
    1416,  1417,  1418,    94,   214,  1421,  1422,  1423,  1424,   167,
      69,    69,    69,  1429,   156,    69,   156,    69,    94,   156,
     214,    94,   214,   214,  1440,   214,   214,  1443,   214,   214,
     167,   214,   214,   214,   156,    94,    94,    94,   214,   214,
      94,    69,    94,   210,   214,   167,    69,   214,   214,   210,
      52,   214,   214,   214,   214,   214,   214,  1473,  1474,    52,
    1476,  1477,  1478,   214,   214,   156,    94,   214,   210,   214,
     210,    94,   214,   279,   214,   214,   167,   214,    69,   311,
     156,   276,  1498,   156,   311,   104,    69,    69,   210,  1505,
    1506,   167,  1508,   104,   167,  1511,   104,   156,   156,   156,
     104,   104,   156,    94,   156,  1521,  1522,   214,   167,   167,
     167,    94,    94,   167,   156,   167,   104,   214,    69,   214,
    1536,  1537,  1538,   214,    69,  1541,  1542,  1543,   156,  1545,
     214,  1547,  1548,   156,  1550,  1551,    52,    52,   214,   167,
    1556,   214,  1558,    94,   167,  1561,  1562,  1563,  1564,    94,
      69,  1567,    69,    69,    69,   214,   214,   214,    69,   214,
     214,  1577,   214,  1579,    52,   156,   156,   214,   210,   214,
     451,   156,   214,   156,   156,    94,   167,    94,    94,    94,
     156,    69,   214,    94,   167,   167,   214,   311,   156,   156,
     104,   214,   104,   214,   104,   156,   214,   214,   311,   214,
     156,   156,   214,   156,   104,   156,    94,    65,   214,   311,
     289,   156,   311,   124,   214,   214,   167,   156,   156,   156,
     210,   140,   167,   214,   214,   210,   156,   214,   214,   214,
     214,   214,   214,   214,   210,   214,   156,   156,   214,   156,
     156,   156,   210,   210,   311,   156,   214,   214,   167,   210,
     167,   167,   167,   214,   210,   210,   167,   210,   214,   214,
     156,   214,  1678,   214,   289,   156,  2037,   156,   156,   214,
     156,   210,   210,   210,   311,   214,   214,   214,   311,   167,
     210,   200,  1698,   214,   214,   156,   214,   311,  1704,   214,
     210,  1707,   156,   210,   214,  1711,   214,  1713,  1714,  1715,
    1716,  1717,  1718,  1719,  1720,  1721,  1722,  1723,  1724,  1725,
    1726,  1727,   289,   156,   210,   156,   289,   214,   214,   210,
     214,   210,   311,   214,   210,   214,   156,   156,   214,   156,
    1746,  1747,  1748,   311,   214,   159,  1752,  1753,   214,   210,
     214,  1757,   214,   214,  1760,  1761,   210,   214,  1764,  1765,
     214,   104,  1768,   156,   156,   156,   156,   104,   214,  1775,
     214,   311,   311,   156,   104,   104,   214,   210,   311,   210,
    1786,   214,   156,   214,  1790,   156,  1792,   311,   156,  1795,
     210,   210,   311,   210,   214,   214,   104,   214,   156,  1805,
    1806,  1807,  1808,  1809,   156,   156,  1812,   156,  1814,   104,
    1816,   214,   311,   311,   311,  2176,   311,   210,   210,   210,
     210,   214,   214,   214,   214,   214,   214,   210,   156,  1835,
    1836,   214,  1838,  1839,   104,   156,   210,   156,   214,   210,
     214,   214,   210,   214,   214,  1851,   214,   214,   214,   214,
    1856,   214,   210,  1859,   214,   214,   214,   156,   210,   210,
     156,   210,   214,   214,   214,   214,   214,   214,   156,   156,
    1876,  1877,  1878,  1879,   156,  1881,  1882,  1883,  1884,  1885,
    1886,  1887,   210,  1889,  1890,  1891,   214,  1893,  1894,   210,
     214,   210,   214,   214,  1900,   214,  1902,   214,   214,  1905,
     214,   214,  1908,   214,   214,  1911,  1912,  1913,   214,  1915,
     214,   210,  1918,   214,   210,   214,   214,  2278,   214,   369,
    1926,  2282,   210,   210,   214,   214,   214,   214,   210,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   235,  1945,
     237,   238,   239,   240,   241,   242,   243,   244,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   311,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   289,   214,   214,   214,   214,   214,
    1986,   214,  1988,   214,  1990,   214,   214,  1993,  1994,  1995,
    1996,  1997,  1998,  1999,  2000,  2001,  2002,  2003,  2004,  2005,
    2006,   214,   214,   214,   214,   214,   214,   214,  2014,   214,
     214,  2017,  2018,   214,   214,  2021,   214,   214,   214,   214,
    2026,   214,   214,  2029,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,  2044,   214,
    2046,   214,   214,   214,   214,   214,  2052,  2053,   214,  2055,
     214,   214,  2058,  2059,   214,  2061,  2062,   214,   214,   214,
     214,   214,   104,   104,   214,   104,   104,   214,   214,   214,
     214,  2077,  2078,   214,   214,   214,   214,   214,   214,   214,
     214,  2087,  2088,   214,  2090,  2091,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   272,   214,   214,   214,
     104,   104,   104,   214,   214,   214,  2112,   214,  2114,  2115,
    2116,   214,   214,   214,  2120,  2121,   214,   214,   214,  2125,
    2126,  2127,  2128,  2129,  2130,  2131,  2132,   214,   214,  2135,
     214,   214,   214,  2494,  2140,   214,   214,  2143,  2144,   214,
     214,   214,   214,   214,  2150,  2151,   214,   214,   214,   214,
    2156,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,  2177,  2178,   104,   214,  2181,  2182,   214,  2184,  2185,
    2186,  2187,   214,  2189,  2190,   214,  2192,   104,   214,  2195,
    2196,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     104,   104,  2218,   214,   214,   214,  2222,   214,  2224,   214,
     214,   214,  2228,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,  2240,   214,  2242,   214,  2244,  2245,
     214,   214,   214,   214,   214,   214,   214,  2253,   214,  2255,
    2256,  2257,  2258,  2259,   214,   214,   214,  2263,  2264,   214,
     214,   214,  2268,   214,   214,   214,  2272,  2273,  2274,  2275,
    2276,   214,   214,   214,   272,   214,   214,   214,   214,  2285,
    2286,   214,   214,   104,   276,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   568,   214,   214,
     214,  2307,   214,  2309,  2310,   214,   214,   214,   214,  2315,
     214,   276,   214,  2319,  2320,  2321,   214,  2323,   214,   214,
     171,  2327,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,  2349,   214,   214,   214,   214,  2354,  2355,
     214,  2357,   214,   214,   214,   214,   214,  2363,   214,   214,
     214,   214,   214,  2369,  2370,   214,  2372,   214,   214,   214,
    2376,  2377,  2378,  2379,   214,   214,  2382,   214,   214,  2385,
     214,  2387,  2388,   214,   214,   214,   214,   214,   214,   214,
     214,   214,  2398,  2399,   214,   214,  2402,   214,  2404,  2405,
     214,   214,  2408,   214,  2410,   214,   214,   214,   214,   570,
     357,   880,   948,   907,   984,   980,   373,   891,  2424,  2425,
    2426,   509,  2428,  1124,   404,   593,   297,  1351,  2434,  1272,
    2436,  2437,  2438,    -1,  2440,  2441,  2442,  2443,    -1,  2445,
     297,  2447,  2448,  2449,   479,   412,   791,   313,    -1,    -1,
      -1,    -1,  2458,   431,  2460,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   438,  2469,  2470,    -1,  2472,    -1,    -1,  2475,
      -1,    -1,  2478,    -1,    -1,  2481,  2482,    -1,  2484,  2485,
      -1,  2487,    -1,    -1,    -1,  2491,  2492,  2493,    -1,    -1,
      -1,    -1,  2498,    -1,    -1,  2501,    -1,    -1,    -1,    -1,
      -1,  2507,    -1,  2509,    -1,  2511,    -1,    -1,    -1,  2515,
      -1,    -1,  2518,    -1,    -1,  2521,  2522,  2523,    -1,    -1,
    2526,  2527,  2528,    -1,    -1,    -1,  2532,    -1,    -1,    -1,
      -1,    -1,  2538,    -1,    -1,    -1,  2542,    -1,    -1,    -1,
      -1,    -1,  2548,    -1,  2550,    -1,  2552,  2553,  2554,    -1,
    2556,    -1,  2558,    -1,    -1,    -1,    -1,  2563,    -1,    -1,
    2566,    -1,    -1,    -1,    -1,    -1,  2572,    -1,    -1,    -1,
    2576,  2577,  2578,    -1,  2580,    -1,    -1,    -1,    -1,  2585,
      -1,  2587,    -1,  2589,    -1,    -1,    -1,    -1,  2594,    -1,
    2596,    -1,    -1,  2599,    -1,    -1,    -1,  2603,    -1,  2605,
      -1,    -1,    -1,  2609,  2610,  2611,  2612,  2613,  2614,  2615,
    2616,    -1,    -1,  2619,  2620,    -1,    -1,  2623,  2624,    -1,
    2626,  2627,  2628,  2629,  2630,  2631,    -1,  2633,    -1,    -1,
    2636,    -1,    -1,    -1,  2640,     4,    -1,    -1,  2644,  2645,
    2646,    -1,    11,    -1,    -1,    14,    15,    -1,    -1,    18,
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
     429,    -1,   431,    -1,    -1,    -1,    -1,    -1,    -1,   438,
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
      -1,   610,    -1,    -1,   613,    -1,   615,    -1,   617,    -1,
     619,    -1,    -1,    -1,    -1,   624,    -1,   626,    -1,   628,
      -1,   630,    -1,   632,    -1,    -1,    -1,    -1,    -1,    -1,
     639,    -1,   641,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   650,    -1,   652,    -1,    -1,    -1,   656,   657,   658,
     659,    -1,   661,   662,   663,    -1,    -1,    -1,    -1,    -1,
      -1,   670,   671,    -1,    -1,    -1,    -1,   676,    -1,   678,
     679,   680,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   693,   694,    -1,   696,   697,   698,
     699,    -1,    -1,   702,    -1,    -1,    -1,    -1,    -1,   708,
      -1,   710,   711,    -1,   713,   714,   715,   716,    -1,   718,
     719,   720,   721,   722,    -1,   724,    -1,   726,   727,   728,
      -1,    -1,    -1,   732,    -1,    -1,   735,    -1,    -1,   738,
     739,   740,    -1,   742,   743,    -1,   745,   746,   747,    -1,
     749,    -1,    -1,    -1,   753,   754,   755,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   768,
     769,    -1,    -1,   772,   773,    -1,    -1,   776,    -1,    -1,
      -1,    -1,   781,    -1,    -1,   784,   785,   786,    -1,   788,
     789,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   797,    -1,
      -1,   800,   801,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     819,    -1,   821,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   831,    -1,   833,    -1,   835,   836,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   847,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   856,    -1,    -1,
      -1,    -1,   861,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   875,    -1,   877,   878,
     879,   880,    -1,   882,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   891,   892,    -1,   894,    -1,    -1,   897,   898,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   907,   908,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     919,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   928,
      -1,    -1,   931,    -1,    -1,    -1,    -1,   936,    -1,    -1,
     939,    -1,    -1,    -1,   943,    -1,   945,    -1,    -1,   948,
      -1,   950,    -1,    -1,    -1,    -1,   955,    -1,    -1,   958,
     959,   960,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   980,    -1,    -1,    -1,   984,    -1,   986,    -1,   988,
     989,    -1,    -1,    -1,    -1,    -1,    -1,   996,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1020,  1021,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1032,    -1,  1034,    -1,    -1,  1037,    -1,
      -1,    -1,  1041,    -1,  1043,  1044,  1045,  1046,  1047,  1048,
      -1,    -1,    -1,    -1,    -1,  1054,    -1,  1056,    -1,  1058,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1090,  1091,    -1,  1093,    -1,    -1,    -1,    -1,    -1,
      -1,  1100,  1101,  1102,  1103,  1104,    -1,  1106,  1107,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1116,    -1,  1118,
    1119,    -1,    -1,    -1,  1123,  1124,  1125,  1126,  1127,  1128,
    1129,    -1,    -1,    -1,    -1,  1134,    -1,    -1,    -1,    -1,
      -1,  1140,    -1,    -1,    -1,  1144,    -1,    -1,  1147,    -1,
      -1,    -1,    -1,    -1,  1153,    -1,    -1,    -1,    -1,    -1,
    1159,    -1,    -1,    -1,    -1,    -1,  1165,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1190,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1198,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1213,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1246,  1247,  1248,
    1249,    -1,  1251,  1252,  1253,    -1,    -1,    -1,    -1,  1258,
      -1,    -1,    -1,  1262,    -1,    -1,    -1,  1266,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1287,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1304,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1312,    -1,    -1,    -1,    -1,    -1,    -1,
    1319,    -1,    -1,  1322,  1323,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1333,    -1,    -1,    -1,    -1,    -1,
      -1,  1340,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1359,    -1,  1361,    -1,    -1,  1364,    -1,  1366,    -1,  1368,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1380,    -1,  1382,    -1,  1384,    -1,    -1,    -1,    -1,
      -1,    -1,  1391,  1392,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1402,    -1,    -1,    -1,  1406,    -1,    -1,
    1409,    -1,    -1,  1412,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1420,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1437,  1438,
    1439,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1447,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1468,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1485,    -1,  1487,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1509,    -1,    -1,  1512,    -1,    -1,    -1,  1516,    -1,    -1,
      -1,    -1,    -1,    -1,  1523,    -1,  1525,  1526,  1527,  1528,
      -1,    -1,    -1,    -1,  1533,    -1,    -1,    -1,    -1,    -1,
      -1,  1540,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1552,    -1,  1554,  1555,    -1,  1557,    -1,
    1559,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1568,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1582,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1611,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1619,    -1,  1621,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1629,    -1,  1631,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    12,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    28,    -1,    -1,    -1,    -1,    33,    -1,
      -1,    -1,    37,    -1,    -1,    -1,    -1,    42,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    50,    -1,    -1,    53,    54,
      -1,    56,    -1,    58,    59,    -1,    -1,    62,  1697,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    73,    -1,
      -1,    -1,    -1,    -1,    -1,    80,    81,    82,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    98,    -1,  1734,    -1,  1736,    -1,  1738,
    1739,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   118,   119,   120,   121,    -1,    -1,    -1,
      -1,   126,  1761,    -1,   129,   130,    -1,    -1,   133,   134,
      -1,    -1,   137,    -1,    -1,    -1,    -1,    -1,  1777,    -1,
      -1,    -1,    -1,    -1,  1783,    -1,    -1,    -1,    -1,    -1,
      -1,  1790,   157,    -1,   159,    -1,    -1,    -1,    -1,    -1,
    1799,  1800,  1801,  1802,    -1,   170,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   182,   183,  1818,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   192,    -1,   194,
    1829,   196,   197,    -1,  1833,  1834,   201,   202,   203,   204,
      -1,    -1,    -1,    -1,    -1,  1844,    -1,  1846,  1847,    -1,
    1849,  1850,   217,   218,  1853,    -1,  1855,    -1,   223,  1858,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1880,    -1,    -1,   249,   250,    -1,   252,    -1,   254,
      -1,   256,    -1,    -1,    -1,    -1,   261,   262,   263,   264,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1910,   277,    -1,    -1,    -1,    -1,    -1,    -1,   284,
    1919,    -1,    -1,    -1,   289,    -1,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,    -1,    -1,    -1,
      -1,   306,   307,    -1,   309,    -1,    -1,    -1,  1947,    -1,
      -1,    -1,  1951,    -1,    -1,    -1,  1955,    -1,   323,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   332,   333,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   347,  1982,    -1,   350,   351,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   364,
      -1,    -1,   367,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2012,  2013,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   404,
      -1,  2040,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2064,    -1,    -1,    -1,    -1,
    2069,  2070,  2071,  2072,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2086,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2094,    -1,  2096,    -1,    -1,
    2099,    -1,    -1,    -1,    -1,    -1,  2105,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2117,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2142,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2163,    -1,  2165,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2174,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2191,    -1,  2193,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2229,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2246,    -1,  2248,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2278,
      -1,    -1,    -1,  2282,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  2306,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2325,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2345,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2419,    -1,    -1,  2422,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2452,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2467,  2468,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2483,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  2505,  2506,    -1,    -1,
      -1,    -1,    -1,     3,     4,     5,    -1,    -1,     8,    -1,
      10,    -1,    -1,    -1,    -1,    15,    -1,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,  2537,    29,
      30,    -1,    -1,    -1,    34,    35,    36,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    44,    -1,    46,    47,    48,    49,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    57,    -1,    -1,
      -1,  2570,    -1,    63,    64,    -1,    -1,    67,    -1,    -1,
      70,    71,    72,    -1,    -1,  2584,    76,    -1,    78,    -1,
      -1,    -1,    -1,    83,  2593,    -1,    86,    87,    88,    89,
      90,    91,    92,    93,    -1,    -1,    -1,    97,    -1,    -1,
      -1,   101,   102,    -1,    -1,   105,   106,    -1,   108,   109,
     110,    -1,    -1,    -1,    -1,    -1,    -1,   117,    -1,    -1,
      -1,    -1,   122,  2632,    -1,   125,    -1,   127,   128,    -1,
      -1,   131,    -1,  2642,  2643,   135,    -1,    -1,   138,   139,
     140,   141,   142,   143,   144,   145,   146,   147,   148,   149,
      -1,    -1,    -1,    -1,   154,   155,    -1,    -1,   158,    -1,
     160,   161,   162,   163,   164,    -1,   166,    -1,   168,    -1,
      -1,   171,    -1,    -1,    -1,    -1,   176,    -1,   178,   179,
     180,   181,    -1,    -1,   184,   185,   186,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     200,    -1,    -1,    -1,    -1,    -1,    -1,   207,    -1,    -1,
      -1,    -1,   212,    -1,    -1,   215,    -1,    -1,    -1,   219,
     220,   221,    -1,    -1,    -1,   225,   226,   227,   228,    -1,
     230,    -1,   232,   233,   234,    -1,   236,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   247,   248,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   259,
     260,    -1,    -1,    -1,    -1,   265,    -1,   267,    -1,   269,
      -1,   271,    -1,   273,   274,    -1,    -1,    -1,   278,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   287,   288,    -1,
     290,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   302,    -1,   304,    -1,    -1,    -1,   308,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     320,    -1,   322,    -1,   324,    -1,   326,   327,   328,   329,
     330,   331,    -1,    -1,   334,    -1,    -1,    -1,    -1,    -1,
      -1,   341,    -1,    -1,   344,   345,   346,    -1,    -1,    -1,
      -1,    -1,   352,    -1,   354,    -1,    -1,    -1,    -1,    -1,
     360,   361,    -1,    -1,    -1,    -1,    -1,    -1,   368,    -1,
      -1,    -1,   372,   373,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     390,   391,   392,   393,    -1,    -1,    -1,   397,    -1,    -1,
      -1,    -1,    -1,   403,    -1,   405,   406,   407,   408,   409,
     410,   411,   412,    -1,    -1,    -1,    -1,    -1,   418,   419,
      -1,   421,   422,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     440,   441,   442,    -1,   444,   445,    -1,    -1,   448,     3,
       4,     5,   452,    -1,     8,    -1,    10,    -1,    -1,    -1,
      -1,    15,    -1,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    -1,    29,    30,    -1,    -1,    -1,
      34,    35,    36,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      44,    -1,    46,    47,    48,    49,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    57,    -1,    -1,    -1,    -1,    -1,    63,
      64,    -1,    -1,    67,    -1,    -1,    70,    71,    72,    -1,
      -1,    -1,    76,    -1,    78,    -1,    -1,    -1,    -1,    83,
      -1,    -1,    86,    87,    88,    89,    -1,    91,    92,    93,
      -1,    -1,    -1,    97,    -1,    -1,    -1,   101,   102,    -1,
      -1,   105,   106,    -1,   108,   109,   110,    -1,    -1,    -1,
      -1,    -1,    -1,   117,    -1,    -1,    -1,    -1,   122,    -1,
      -1,   125,    -1,   127,   128,    -1,    -1,   131,    -1,    -1,
      -1,   135,    -1,    -1,   138,   139,   140,   141,   142,   143,
     144,   145,   146,   147,   148,   149,    -1,    -1,    -1,    -1,
     154,   155,    -1,    -1,   158,    -1,   160,   161,   162,   163,
     164,    -1,   166,    -1,   168,    -1,    -1,   171,    -1,    -1,
      -1,    -1,   176,    -1,   178,   179,   180,   181,    -1,    -1,
     184,   185,   186,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   200,    -1,    -1,    -1,
      -1,    -1,    -1,   207,    -1,    -1,    -1,    -1,   212,    -1,
      -1,   215,    -1,    -1,    -1,   219,   220,   221,    -1,    -1,
      -1,   225,   226,   227,   228,    -1,   230,    -1,   232,   233,
     234,    -1,   236,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   247,   248,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   259,   260,    -1,    -1,    -1,
      -1,   265,    -1,   267,    -1,   269,    -1,   271,    -1,   273,
     274,    -1,    -1,    -1,   278,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   287,   288,    -1,   290,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   302,    -1,
     304,    -1,    -1,    -1,   308,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   320,    -1,   322,    -1,
     324,    -1,   326,   327,   328,   329,   330,   331,    -1,    -1,
     334,    -1,    -1,    -1,    -1,    -1,    -1,   341,    -1,    -1,
     344,   345,   346,    -1,    -1,    -1,    -1,    -1,   352,    -1,
     354,    -1,    -1,    -1,    -1,    -1,   360,   361,    -1,    -1,
      -1,    -1,    -1,    -1,   368,    -1,    -1,    -1,   372,   373,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   390,   391,   392,   393,
      -1,    -1,    -1,   397,    -1,    -1,    -1,    -1,    -1,   403,
      -1,   405,   406,   407,   408,   409,   410,   411,   412,    -1,
      -1,    -1,    -1,    -1,   418,   419,    -1,   421,   422,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   440,   441,   442,    -1,
     444,   445,    -1,    -1,   448,    -1,    -1,    -1,   452
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint16 yystos[] =
{
       0,     3,     4,     5,     8,    10,    15,    17,    18,    19,
      20,    21,    22,    23,    24,    25,    26,    27,    29,    30,
      34,    35,    36,    44,    46,    47,    48,    49,    57,    63,
      64,    67,    70,    71,    72,    76,    78,    83,    86,    87,
      88,    89,    91,    92,    93,    97,   101,   102,   105,   106,
     108,   109,   110,   117,   122,   125,   127,   128,   131,   135,
     138,   139,   140,   141,   142,   143,   144,   145,   146,   147,
     148,   149,   154,   155,   158,   160,   161,   162,   163,   164,
     166,   168,   171,   176,   178,   179,   180,   181,   184,   185,
     186,   200,   207,   212,   215,   219,   220,   221,   225,   226,
     227,   228,   230,   232,   233,   234,   236,   247,   248,   259,
     260,   265,   267,   269,   271,   273,   274,   278,   287,   288,
     290,   302,   304,   308,   320,   322,   324,   326,   327,   328,
     329,   330,   331,   334,   341,   344,   345,   346,   352,   354,
     360,   361,   368,   372,   373,   390,   391,   392,   393,   397,
     403,   405,   406,   407,   408,   409,   410,   411,   412,   418,
     419,   421,   422,   440,   441,   442,   444,   445,   448,   452,
     457,   458,   459,   460,   461,   462,   463,   464,   469,   470,
     471,   472,   473,   474,   475,   476,   477,   478,   487,   489,
     490,   491,   492,   493,   494,   495,   496,   497,   498,   500,
     502,   503,   504,   505,   506,   507,   508,   513,   516,   518,
     519,   520,   521,   522,   523,   524,   525,   526,   527,   529,
     530,   532,   536,   537,   538,   540,   541,   542,   543,   546,
     549,   550,   551,   552,   553,   556,   559,   562,   564,   566,
     568,   572,   573,   574,   575,   576,   577,   578,   579,   581,
     583,   584,   585,   587,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   602,   604,   606,   608,   612,   616,
     618,   620,   621,   622,   623,   624,   625,   626,   628,   629,
     630,   631,   639,   640,   642,   644,   645,   646,   647,   648,
     649,   650,   651,   652,   653,   657,   659,   661,   664,   667,
     669,   671,   672,   674,   675,   676,   677,   678,   679,   682,
     683,   685,   686,   687,   688,   691,   214,     6,   214,   214,
     156,   210,   693,   214,    69,    94,   156,   167,   694,   694,
     694,   104,   214,   694,   214,   214,   693,   214,   214,   693,
     693,   214,   214,   693,   693,   214,   214,   693,   214,   214,
     693,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,    32,    73,    84,    85,   214,   660,   214,
     214,   693,   214,   214,   214,   214,   214,   214,   693,   693,
     214,   693,   214,   693,   214,   214,   311,   214,   693,   214,
     214,   214,   214,   214,   214,   693,   214,   214,   693,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   353,
     214,   694,   214,   214,   693,   214,   214,   693,   214,   214,
     693,   214,   214,   214,   104,   214,   693,   214,   214,   693,
     214,   214,   214,   693,   214,   693,   693,   214,   214,   214,
     693,   214,   104,   169,   214,   693,   169,   214,   693,   214,
     214,   693,   693,   214,   214,   214,   693,   693,   693,   214,
     693,   214,   693,   104,   214,   214,   214,   214,   214,   214,
     694,   214,   214,   693,   214,   214,   694,   694,   214,   214,
     214,   214,   214,   214,   214,   214,   104,   214,   214,   214,
     214,   214,   361,   214,   214,   693,   214,   214,   693,   693,
     694,   693,   693,   214,   214,   694,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   311,     0,    90,   459,   132,   693,   206,   255,
     258,   266,   275,   465,   466,   468,   535,    38,    39,   355,
     356,   357,   358,   359,    77,    79,   116,   343,   398,   399,
     400,   401,   402,   693,    16,   111,   136,   205,   694,   694,
     693,    43,   172,   173,   174,   175,   479,   481,   482,   483,
     484,   641,   693,   488,   693,   693,   694,   104,   450,   499,
     303,   340,   449,   453,   454,   455,   501,   188,    16,    41,
     136,   152,   191,   205,   213,   224,   312,   321,   353,   510,
     531,   535,   136,   191,   205,   514,   515,    31,    45,   358,
     517,   539,   693,   317,   635,   693,   693,   317,   693,   317,
     693,   693,   693,   557,   693,   563,   693,   565,   693,   567,
     693,   569,   693,   557,   565,   569,   580,   693,   582,   693,
     586,   693,   588,   693,   447,   191,   693,   693,   693,   191,
     317,   635,   693,   619,   693,   633,   693,   693,   693,   693,
     627,   693,   693,   693,   632,   693,   641,   641,   643,   693,
     693,   693,   336,   693,   694,   156,   693,   694,   693,   317,
     693,   693,    12,    28,    33,    37,    42,    50,    53,    54,
      56,    58,    59,    62,    73,    80,    81,    82,    98,   118,
     119,   120,   121,   126,   129,   130,   133,   134,   137,   157,
     159,   170,   182,   183,   192,   194,   196,   197,   201,   202,
     203,   204,   217,   218,   223,   249,   250,   252,   254,   256,
     261,   262,   263,   264,   277,   284,   289,   291,   292,   293,
     294,   295,   296,   297,   298,   299,   300,   301,   306,   307,
     309,   323,   332,   333,   347,   350,   351,   364,   367,   404,
     655,   656,   465,   468,   535,   694,     9,    74,    85,   103,
     123,   177,   183,   188,   216,   222,   270,   348,   376,   673,
     214,   693,   379,   693,   693,    40,   400,   413,   423,   427,
     429,   430,   434,   684,   690,   209,   336,   400,   414,   415,
     416,   423,   425,   426,   435,   436,   437,   438,   439,   443,
     446,   447,   689,   689,   689,   690,   689,   420,   692,   599,
     635,   693,   214,     6,     6,   214,   214,   214,   694,   214,
     693,   694,   570,   571,   693,   571,   214,   694,   694,   214,
     214,   693,   214,   214,   214,   694,   214,   694,   104,   509,
     694,   694,   694,    85,   694,   660,   151,   214,   214,   563,
     557,   601,   637,   693,   658,   694,   214,   693,   214,   214,
     191,   214,   214,   214,   694,   617,   638,   693,   617,   214,
     554,   555,   693,   214,   565,   668,   694,   670,   694,   569,
     567,   613,   614,   693,   599,   214,   214,   601,   599,   214,
     694,   693,   214,   694,   693,   214,   694,   609,   610,   693,
     214,   214,   214,   693,   619,   214,   214,   214,   632,   693,
     214,   214,   214,   693,   214,   214,   693,   214,   693,   214,
     214,   693,   214,   214,   214,   214,   528,   693,   214,   104,
     214,   693,   272,    52,   607,   599,   214,   694,   560,   561,
     693,   214,   214,   633,   694,    52,   605,   214,   599,   599,
      52,   603,   693,   112,   113,   114,   115,   369,   467,   214,
     693,   214,   281,   316,   693,   214,   214,   214,   693,   214,
     544,   545,   693,   627,   547,   548,   693,   214,   317,   693,
     279,   214,   694,   214,   694,   694,   276,   311,   311,   104,
     104,   104,   104,   104,   214,   104,   214,   214,   693,   214,
     214,   214,   214,   694,   693,   693,   693,   693,   694,   694,
     140,   200,   694,   693,   693,   693,   693,   693,   485,   693,
     485,   486,   693,   486,   338,   694,   694,   693,   694,   104,
     693,   451,   214,   104,   693,   104,   693,   104,   693,   104,
     104,   104,    11,    14,    60,   187,   211,   214,   219,   229,
     282,   283,   305,   376,   694,   693,   511,   694,   214,   694,
     512,   694,   311,   214,   511,   214,   214,   311,   693,   694,
     311,   694,   214,   694,   214,   104,   214,   693,   694,   635,
     693,   342,   693,    65,   693,   694,   693,   694,   694,   694,
     693,   693,   693,   693,   693,   214,   693,   693,   694,   311,
     214,   694,   694,   694,   214,   635,   342,     7,    61,    96,
     318,   339,   694,   693,   693,   693,   693,   693,   693,   693,
     694,   124,   694,   693,   693,   214,   694,   214,   693,   214,
     693,   693,   693,   694,   342,   694,    13,   365,   289,   214,
     311,   289,    13,    51,   694,   289,   694,   311,   311,   654,
     693,   214,   693,   311,   693,   693,    75,   214,   370,   693,
     285,   693,   214,   214,   693,   694,   214,   694,   694,   311,
     214,   693,   214,   214,   693,   165,   693,   289,   693,   193,
     693,   195,   693,   285,   693,   214,   693,   693,   693,   693,
     214,   693,   214,   693,   104,   214,   693,   253,   693,   693,
     694,   694,   694,   693,   214,   694,   311,   285,   693,   214,
     251,   348,   694,   693,   693,   693,   694,   693,   693,   311,
     693,   693,   693,   214,   693,   694,   694,   694,   693,   693,
     693,   214,   159,   694,   214,   311,   188,   189,   245,   313,
     325,   188,   189,   245,   694,   214,   694,   693,   693,   214,
     694,   694,   693,   693,   214,   694,   693,   214,   214,   104,
     693,   104,   104,   374,   375,   376,   377,   378,   380,   381,
     382,   383,   384,   386,   387,   388,   389,   693,   693,   693,
     681,   693,   693,   680,   681,   214,   214,   509,   311,   311,
     693,   104,   104,   693,   693,   694,   311,   311,   311,   104,
     311,   311,   104,   311,   311,   214,   214,   214,   214,   104,
     214,   635,   693,   342,   214,   694,   694,   694,   214,   214,
     693,   570,   694,   571,   694,   694,   214,   214,   694,   693,
     214,   214,   694,   694,   694,   214,   693,   637,   694,   214,
     214,   665,   666,   694,   214,   214,   694,   638,   693,   617,
     555,   693,   214,   214,   614,   615,   693,   615,   693,    55,
     214,   694,   662,   663,   694,   694,   694,   694,   214,   694,
     610,   611,   693,   611,   693,   214,   694,   214,   693,   214,
     214,   693,   693,   662,   662,   214,   693,   214,   693,   214,
     214,   693,    52,   694,   694,   561,   693,   214,   693,    52,
     694,   693,    52,   694,   214,   694,   694,   694,   694,   694,
     694,   112,   113,   114,   115,   214,   214,   694,   214,   394,
     395,   396,   214,   545,   214,   694,   548,   693,   693,   693,
     694,   600,   636,   693,   214,   694,   694,   693,   104,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   693,   693,   214,   214,   214,
     214,   214,   214,   694,   693,   693,   694,   694,   693,   694,
     694,   104,   214,   214,   693,   693,   693,   693,   693,   693,
     693,   693,   693,   693,   311,   694,   214,   214,   694,   694,
     214,   214,   214,   214,   214,   694,   694,   214,   694,   694,
     214,   694,   693,   208,   214,   694,   693,   694,   214,   694,
     214,   694,   694,   558,   693,   558,   558,   558,   558,   214,
     693,   214,   693,   694,   214,   600,   694,   694,   694,   600,
     693,   694,    95,   185,   214,   318,   660,   693,   693,   694,
     694,   694,   634,   693,   634,   693,   693,   558,   693,   693,
     214,   694,   694,   694,   694,   153,   214,   337,   693,   214,
     214,   214,   693,   694,   694,   214,   311,   693,   214,   694,
     214,   268,   366,   693,   214,   251,   289,   214,   214,   251,
     310,   214,   214,   693,   214,   214,   251,   310,   214,   214,
     214,   214,   693,   214,   214,   694,   693,   214,   214,    99,
     214,   694,   214,   214,   214,   214,   214,   214,   214,   694,
     214,   694,   214,   214,   214,   251,   310,   694,   198,   199,
     214,   199,   214,   693,   214,   214,   214,   214,   214,   214,
     693,   694,   214,   214,   694,   214,   214,   214,   214,   214,
     214,   694,   694,   693,   694,   214,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   694,   214,   694,   214,   214,   694,   214,
     214,   150,   214,   214,   694,   214,   693,   693,   693,   693,
     694,   693,   693,   100,   184,   693,   694,   694,   694,   214,
     693,   214,   214,   214,   693,   214,   214,   694,   214,   693,
     214,   694,   104,   680,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   214,   693,
     104,   694,   214,   104,   693,   104,   693,   693,   342,   693,
     214,   694,   694,   214,   214,   214,   214,   214,   694,   214,
     214,   693,   694,   694,   694,   214,   666,   694,   600,   214,
     694,   693,   694,   615,   693,   694,   214,   663,   694,   214,
     694,   694,   694,   694,   611,   693,   214,   361,   660,   214,
     214,   533,   693,   693,   214,   272,   214,   693,   694,   214,
     693,   694,   214,   693,   694,   214,   693,   694,   214,   694,
     694,   694,   694,   214,   693,   694,   694,   694,   694,   369,
     214,   281,   694,   214,   394,   214,   394,   214,   558,   693,
     214,   693,   694,   636,   694,   214,   208,   214,   214,   693,
     214,   214,   214,   480,   693,   694,   694,   214,   694,   694,
     694,   214,   214,   214,   132,   693,   104,   132,   693,   104,
     104,   694,   214,   694,   694,   694,   693,   694,   314,   693,
     214,   208,   214,   693,   214,   214,   694,   694,   214,   693,
     214,   214,   214,   214,   214,   214,   693,   694,   694,   694,
     314,   693,   694,   694,   694,   694,    95,   185,   214,   318,
     694,   694,   214,   694,   694,   214,   693,   214,   214,   660,
     693,   214,   694,   214,   214,   694,   214,   660,   693,   214,
     694,   694,   694,   694,   214,   336,   694,   214,   335,   693,
     214,   214,   311,   214,   311,   214,   694,   694,   214,   693,
     214,   371,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   214,   693,   214,   693,   214,   693,   214,   214,
     214,   199,   214,   693,   693,   199,   214,   214,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   694,   214,   214,
     214,   214,   214,   214,   214,   214,   214,   214,   214,   214,
     214,   214,   693,   694,   214,   214,   694,   214,   694,   214,
     694,   214,   214,   694,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   214,   104,   693,
     104,   693,   693,   314,   693,   214,   214,   694,   694,   694,
     214,   694,   694,   694,   214,   694,   693,   694,   214,   694,
     694,   694,   214,   214,   214,   214,   694,   693,   214,   235,
     237,   238,   239,   240,   241,   242,   243,   244,   534,   533,
     693,   214,   694,   214,   214,   693,   694,   214,   694,   214,
     214,   694,   214,   214,   693,   693,   693,   693,   214,   694,
     694,   694,   694,   694,   214,   214,   694,   694,   694,   214,
     214,   214,   214,   214,   214,   693,   693,   693,   694,   694,
     694,   694,   214,   693,   693,   693,   693,   693,   694,   693,
     214,   693,   214,   694,   214,   693,   694,   214,   214,   214,
     214,   214,   694,   214,   694,   694,   694,   693,   694,   214,
     694,   694,   214,   694,   214,   694,   694,   694,   694,   694,
     694,   694,   694,   214,   214,   694,   694,   694,   214,   214,
     694,   214,   693,   694,   694,   694,   694,   214,   694,   214,
     335,   337,   693,   214,   214,   214,   694,   214,   214,   371,
     214,   214,   214,   694,   214,   693,   693,   214,   214,   693,
     214,   693,   214,   214,   694,   694,   694,   694,   694,   694,
     694,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   314,   693,   693,   694,   694,   694,   214,   104,   694,
     214,   214,   694,   214,   214,   694,   214,   214,   214,   214,
     660,   693,   694,   694,   533,   533,   272,   694,   694,   694,
     214,   694,   694,   694,   694,   693,   214,   214,   276,   693,
     693,   693,   693,   214,   214,   214,   694,   214,   694,   214,
     693,   694,   694,   214,   694,   694,   132,   693,   132,   693,
     693,   214,   214,   693,   214,   694,   694,   694,   694,   693,
     208,   214,   214,   694,   694,   214,   214,   694,   694,   214,
     694,   214,   694,   694,   694,   694,   694,   214,   694,   214,
     214,   214,   694,   214,   693,   694,   694,   214,   214,   214,
     694,   694,   214,   694,   214,   214,   693,   214,   693,   214,
     214,   214,   214,   693,   214,   214,   361,   363,   660,   694,
     694,   214,   214,   694,   214,   694,   214,   214,   694,   694,
     694,   694,   214,   694,   694,   693,   694,   693,   214,   694,
     694,   214,   214,   214,   214,   214,   214,   694,   214,   214,
     694,   214,   694,   214,   214,   694,   693,   214,   214,   276,
     214,   214,   694,   214,   694,   694,   214,   694,   693,   693,
     214,   694,   214,   694,   694,   694,   694,   694,   214,   694,
     694,   214,   694,   214,   214,   694,   694,   694,   694,   694,
     214,   660,   693,   214,   660,   693,   214,   694,   694,   214,
     214,   214,   214,   214,   214,   693,   694,   214,   694,   694,
     214,   214,   694,   214,   694,   694,   214,   694,   694,   214,
     693,   214,   694,   214,   214,   104,   214,   214,   214,   693,
     214,   214,   694,   214,   694,   694,   694,   214,   694,   214,
     214,   214,   208,   214,   694,   214,   318,   694,   694,   214,
     694,   694,   694,   694,   214,   694,   214,   694,   694,   694,
     214,   214,   361,   214,   694,   214,   694,   214,   694,   694,
     694,   214,   214,   694,   694,   214,   214,   693,   693,   694,
     694,   214,   694,   694,   214,   694,   694,   214,   318,   694,
     214,   318,   694,   694,   694,   694,   214,   694,   214,   694,
     214,   694,   694,   694,   693,   214,   214,   214,   694,   214,
     214,   694,   214,   214,   104,   693,   693,   694,   694,   694,
     214,   694,   214,   694,   214,   318,   694,   214,   318,   694,
     693,   694,   694,   694,   214,   214,   694,   694,   694,   214,
     660,   214,   694,   214,   694,   214,   693,   693,   694,   214,
     214,   694,   694,   214,   694,   214,   694,   214,   318,   214,
     694,   694,   694,   694,   694,   214,   694,   214,   694,   214,
     214,   693,   214,   694,   214,   214,   694,   214,   214,   694,
     694,   214,   694,   214,   694,   214,   694,   694,   694,   214,
     694,   214,   694,   214,   693,   694,   214,   694,   694,   694,
     214,   694,   214,   693,    68,   214,   694,   694,   214,   694,
     214,   693,   694,   694,   214,   694,   214,   214,   214,   266,
     286,   694,   694,   694,   694,   694,   694,   694,   694,   694,
     694,   214,   214,   694,   694,   214,   266,   286,   694,   694,
     694,   694,   694,   694,   694,   694,   693,   694,   214,   214,
     694,   214,   694,   693,   693,   694,   694,   694,   214
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
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
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
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
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



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
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
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

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
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

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

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
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

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
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
#line 165 "p.y"
    { 
          return 0;
         ;}
    break;

  case 6:
#line 176 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; delete (yyvsp[(1) - (1)].bclist); ;}
    break;

  case 7:
#line 178 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 15:
#line 187 "p.y"
    {;}
    break;

  case 23:
#line 196 "p.y"
    {;}
    break;

  case 24:
#line 198 "p.y"
    {;}
    break;

  case 30:
#line 205 "p.y"
    {;}
    break;

  case 49:
#line 225 "p.y"
    { domain->setMFTT((yyvsp[(1) - (1)].mftval).first, (yyvsp[(1) - (1)].mftval).second); ;}
    break;

  case 50:
#line 227 "p.y"
    { domain->setHFTT((yyvsp[(1) - (1)].hftval).first, (yyvsp[(1) - (1)].hftval).second); ;}
    break;

  case 77:
#line 255 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 78:
#line 257 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 79:
#line 259 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 80:
#line 261 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 81:
#line 263 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 82:
#line 265 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 99:
#line 283 "p.y"
    { if(domain->setComplexNeuman((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; ;}
    break;

  case 101:
#line 286 "p.y"
    { if(domain->setComplexDirichlet((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; ;}
    break;

  case 111:
#line 297 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 112:
#line 299 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 119:
#line 307 "p.y"
    {;}
    break;

  case 128:
#line 318 "p.y"
    {;}
    break;

  case 129:
#line 320 "p.y"
    {;}
    break;

  case 130:
#line 322 "p.y"
    {;}
    break;

  case 131:
#line 324 "p.y"
    {;}
    break;

  case 132:
#line 326 "p.y"
    {;}
    break;

  case 133:
#line 328 "p.y"
    {;}
    break;

  case 134:
#line 330 "p.y"
    {;}
    break;

  case 135:
#line 332 "p.y"
    {;}
    break;

  case 147:
#line 347 "p.y"
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[(3) - (5)].ival)); 
            domain->solInfo().nsample = (yyvsp[(4) - (5)].ival);
          ;}
    break;

  case 148:
#line 354 "p.y"
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[(3) - (4)].ival));
          ;}
    break;

  case 150:
#line 361 "p.y"
    { if ((yyvsp[(2) - (5)].ival) == OutputInfo::Attribute)  geoSource->setGroupAttribute((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival)-1);
          else if ((yyvsp[(2) - (5)].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[(2) - (5)].ival));  exit(-1); }
        ;}
    break;

  case 151:
#line 366 "p.y"
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
#line 378 "p.y"
    { if ((yyvsp[(2) - (6)].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[(2) - (6)].ival));  exit(-1); }
        ;}
    break;

  case 154:
#line 385 "p.y"
    { geoSource->setGroupRandomProperty((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].rprop),(yyvsp[(4) - (6)].fval),(yyvsp[(5) - (6)].fval)); ;}
    break;

  case 155:
#line 389 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (5)].fval)); ;}
    break;

  case 156:
#line 393 "p.y"
    { domain->solInfo().curSweepParam = 0; domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies1(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); ;}
    break;

  case 157:
#line 395 "p.y"
    { domain->solInfo().curSweepParam = 0; domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies2(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); ;}
    break;

  case 158:
#line 397 "p.y"
    { domain->solInfo().curSweepParam = 0; domain->setFrequencySet(0); geoSource->setImpe((yyvsp[(4) - (8)].fval)); domain->addFrequencies(2.0*PI*(yyvsp[(4) - (8)].fval), 2.0*PI*(yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].ival), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 159:
#line 399 "p.y"
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
#line 419 "p.y"
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
#line 445 "p.y"
    { domain->solInfo().curSweepParam = (yyvsp[(3) - (7)].ival); if ((yyvsp[(3) - (7)].ival) == 0) geoSource->setImpe((yyvsp[(6) - (7)].fval)); ;}
    break;

  case 162:
#line 449 "p.y"
    { domain->setFrequencySet((yyvsp[(2) - (8)].ival)); domain->solInfo().curSweepParam = (yyvsp[(2) - (8)].ival); if ((yyvsp[(2) - (8)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (8)].fval)); domain->addFrequencies1(2.0*PI*(yyvsp[(5) - (8)].fval), 2.0*PI*(yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 163:
#line 451 "p.y"
    { domain->setFrequencySet((yyvsp[(2) - (8)].ival)); domain->solInfo().curSweepParam = (yyvsp[(2) - (8)].ival); if ((yyvsp[(2) - (8)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (8)].fval)); domain->addFrequencies2(2.0*PI*(yyvsp[(5) - (8)].fval), 2.0*PI*(yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 164:
#line 453 "p.y"
    { domain->setFrequencySet((yyvsp[(2) - (9)].ival)); domain->solInfo().curSweepParam = (yyvsp[(2) - (9)].ival); if ((yyvsp[(2) - (9)].ival) == 0) geoSource->setImpe((yyvsp[(5) - (9)].fval)); domain->addFrequencies(2.0*PI*(yyvsp[(5) - (9)].fval), 2.0*PI*(yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].ival), (yyvsp[(8) - (9)].ival)); ;}
    break;

  case 165:
#line 455 "p.y"
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
#line 475 "p.y"
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
#line 508 "p.y"
    { domain->solInfo().getSweepParams()->pade_pivot = true; domain->solInfo().getSweepParams()->pade_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 173:
#line 512 "p.y"
    { domain->solInfo().getSweepParams()->pade_poles = true; ;}
    break;

  case 174:
#line 514 "p.y"
    { domain->solInfo().getSweepParams()->pade_poles = true; 
          domain->solInfo().getSweepParams()->pade_poles_sigmaL = (yyvsp[(2) - (4)].fval); domain->solInfo().getSweepParams()->pade_poles_sigmaU = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 175:
#line 519 "p.y"
    { geoSource->setImpe((yyvsp[(2) - (3)].fval)); domain->addCoarseFrequency(2.0*PI*(yyvsp[(2) - (3)].fval)); ;}
    break;

  case 176:
#line 522 "p.y"
    { domain->addFrequencies(2.0*PI*(yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].ival)); ;}
    break;

  case 177:
#line 526 "p.y"
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
#line 565 "p.y"
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
#line 615 "p.y"
    { geoSource->binaryInput = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 181:
#line 617 "p.y"
    { geoSource->binaryInput = bool((yyvsp[(3) - (5)].ival));
            std::string prefix = (yyvsp[(4) - (5)].strval);
            clusterData_ = prefix + ".msh";
            decomposition_ = prefix + ".dec";
            connectivity_ = prefix + ".con";
            subdomains_ = prefix + ".sub";
            //fprintf(stderr, "clusterData_ = %s\n", clusterData_.c_str());
            //fprintf(stderr, "decomposition_ = %s\n", decomposition_.c_str());
            //fprintf(stderr, "connectivity_ = %s\n", connectivity_.c_str());
            //fprintf(stderr, "subdomains_ = %s\n", subdomains_.c_str());
          ;}
    break;

  case 182:
#line 629 "p.y"
    { geoSource->binaryOutput = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 183:
#line 631 "p.y"
    { geoSource->setGeo((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 184:
#line 633 "p.y"
    { geoSource->setDecomp((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 185:
#line 635 "p.y"
    { geoSource->setGlob((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 186:
#line 637 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 187:
#line 639 "p.y"
    { geoSource->setCpuMap((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 188:
#line 643 "p.y"
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[(2) - (3)].ival)); 
#endif
	;}
    break;

  case 189:
#line 651 "p.y"
    {if(decInit==0) decInit = new DecInit(); ;}
    break;

  case 190:
#line 653 "p.y"
    {decInit->file = strdup((yyvsp[(3) - (4)].strval));;}
    break;

  case 191:
#line 655 "p.y"
    {decInit->nsubs = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 192:
#line 657 "p.y"
    {decInit->weight = true; ;}
    break;

  case 193:
#line 659 "p.y"
    {decInit->memory = true; ;}
    break;

  case 194:
#line 661 "p.y"
    {decInit->exitAfterDec = true;;}
    break;

  case 195:
#line 663 "p.y"
    {decInit->skip = true;;}
    break;

  case 196:
#line 665 "p.y"
    {decInit->nosa = true; ;}
    break;

  case 197:
#line 667 "p.y"
    {decInit->trivial = true; ;}
    break;

  case 198:
#line 669 "p.y"
    {decInit->fsgl = true; ;}
    break;

  case 199:
#line 673 "p.y"
    {;}
    break;

  case 200:
#line 675 "p.y"
    { weightList[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 201:
#line 679 "p.y"
    {;}
    break;

  case 202:
#line 681 "p.y"
    { fieldWeightList[(int)Element::Acoustic] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 203:
#line 683 "p.y"
    { fieldWeightList[(int)Element::Structural] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 204:
#line 685 "p.y"
    { fieldWeightList[(int)Element::Thermal] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 205:
#line 687 "p.y"
    { fieldWeightList[(int)Element::Fluid] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 206:
#line 690 "p.y"
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = 0; ;}
    break;

  case 207:
#line 692 "p.y"
    { (yyval.mftval).first = new MFTTData; (yyval.mftval).second = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 208:
#line 694 "p.y"
    { (yyval.mftval).first->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 209:
#line 698 "p.y"
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = 0; ;}
    break;

  case 210:
#line 700 "p.y"
    { (yyval.hftval).first = new MFTTData; (yyval.hftval).second = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 211:
#line 702 "p.y"
    { (yyval.hftval).first->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 212:
#line 706 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 213:
#line 708 "p.y"
    { (yyval.ival) = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 214:
#line 710 "p.y"
    { domain->setLoadFactor((yyval.ival), (yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 215:
#line 712 "p.y"
    { domain->setLoadFactorMFTT((yyval.ival), (yyvsp[(2) - (5)].ival), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 216:
#line 714 "p.y"
    { domain->setLoadFactorHFTT((yyval.ival), (yyvsp[(2) - (5)].ival), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 224:
#line 727 "p.y"
    { geoSource->addCFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 225:
#line 731 "p.y"
    { geoSource->addCoefInfo((yyvsp[(2) - (4)].ival)-1,(yyvsp[(4) - (4)].coefdata)); ;}
    break;

  case 226:
#line 735 "p.y"
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[(1) - (4)].ival)-1,(yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 227:
#line 737 "p.y"
    { (yyval.coefdata).setCoef((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 228:
#line 741 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 229:
#line 743 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 230:
#line 747 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 231:
#line 749 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 232:
#line 753 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 233:
#line 755 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 234:
#line 759 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 235:
#line 761 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 236:
#line 765 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (11)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
	  (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval);
	  (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 237:
#line 771 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (13)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (13)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (13)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (13)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (13)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (13)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (13)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (13)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (13)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (13)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (13)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (13)].fval); ;}
    break;

  case 238:
#line 778 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (14)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (14)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (14)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (14)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (14)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (14)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (14)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (14)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (14)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (14)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (14)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (14)].fval); (yyval.ldata).d[11] = (yyvsp[(13) - (14)].fval); ;}
    break;

  case 239:
#line 787 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (5)].ival)-1;  (yyval.ldata).matid = (yyvsp[(2) - (5)].ival)-1; (yyval.ldata).d[7] = (yyvsp[(3) - (5)].fval); (yyval.ldata).d[8] = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 241:
#line 792 "p.y"
    { geoSource->addLayMat((yyvsp[(2) - (2)].ldata).matid, (yyvsp[(2) - (2)].ldata).d); ;}
    break;

  case 242:
#line 799 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (7)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (7)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (7)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (7)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (7)].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[(6) - (7)].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; ;}
    break;

  case 243:
#line 805 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (9)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (9)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (9)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (9)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (9)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (9)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (9)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (9)].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; ;}
    break;

  case 244:
#line 810 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (11)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval);
          (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); (yyval.ldata).d[9] = 0; ;}
    break;

  case 245:
#line 814 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (12)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (12)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (12)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (12)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (12)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (12)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (12)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (12)].fval); 
          (yyval.ldata).d[7] = (yyvsp[(9) - (12)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (12)].fval); (yyval.ldata).d[9] = (yyvsp[(11) - (12)].fval); ;}
    break;

  case 247:
#line 821 "p.y"
    { domain->addDMass((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 248:
#line 823 "p.y"
    { domain->addDMass((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(5) - (6)].fval),(yyvsp[(4) - (6)].ival)-1); ;}
    break;

  case 250:
#line 828 "p.y"
    { domain->setGravity((yyvsp[(2) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 252:
#line 833 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (5)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (5)].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[(4) - (5)].strval); ;}
    break;

  case 253:
#line 837 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (4)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (4)].strval);;}
    break;

  case 254:
#line 840 "p.y"
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[(2) - (4)].strval);
          domain->solInfo().nRestart = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 255:
#line 845 "p.y"
    { geoSource->setControlFile((yyvsp[(2) - (3)].strval));
         geoSource->setControlRoutine((char *) "controlObj");;}
    break;

  case 256:
#line 850 "p.y"
    { geoSource->setControlRoutine((yyvsp[(2) - (3)].strval)); ;}
    break;

  case 257:
#line 854 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 258:
#line 859 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Actuators; }
          if(geoSource->setActuatorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)            return -1; ;}
    break;

  case 259:
#line 865 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdf; }
          if(geoSource->setUsdfLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)       return -1; ;}
    break;

  case 260:
#line 872 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdd;
          if(geoSource->setUsddLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)    return -1; ;}
    break;

  case 261:
#line 879 "p.y"
    { 
    domain->solInfo().sensitivity = true;
    domain->senInfo = new SensitivityInfo[50];  // maximum number of sensitivities are fixed to 50
  ;}
    break;

  case 262:
#line 884 "p.y"
    { domain->addSensitivity((yyvsp[(2) - (3)].sinfo)); ;}
    break;

  case 263:
#line 888 "p.y"
    { (yyval.sinfo).type = (SensitivityInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.sinfo).method = (SensitivityInfo::Method) (yyvsp[(2) - (3)].ival); (yyval.sinfo).numParam = (yyvsp[(3) - (3)].ival);  
  ;}
    break;

  case 264:
#line 893 "p.y"
    { numColumns = 3; ;}
    break;

  case 265:
#line 895 "p.y"
    { numColumns = 6; ;}
    break;

  case 266:
#line 897 "p.y"
    { numColumns = 3; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 267:
#line 899 "p.y"
    { numColumns = 6; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 268:
#line 901 "p.y"
    { numColumns = 3; domain->outFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 269:
#line 903 "p.y"
    { numColumns = 6; domain->outFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 270:
#line 905 "p.y"
    { numColumns = 3; domain->outFlag = (yyvsp[(2) - (4)].ival); geoSource->setOutLimit((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 271:
#line 907 "p.y"
    { numColumns = 6; domain->outFlag = (yyvsp[(2) - (4)].ival); geoSource->setOutLimit((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 272:
#line 909 "p.y"
    { (yyvsp[(2) - (3)].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[(2) - (3)].oinfo)); ;}
    break;

  case 273:
#line 913 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 274:
#line 915 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); ;}
    break;

  case 275:
#line 917 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; ;}
    break;

  case 276:
#line 919 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); 
          if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1;;}
    break;

  case 277:
#line 922 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; ;}
    break;

  case 278:
#line 924 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; ;}
    break;

  case 279:
#line 926 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); (yyval.oinfo).sentype = 1; ;}
    break;

  case 280:
#line 928 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); (yyval.oinfo).sentype = 1; ;}
    break;

  case 281:
#line 930 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 282:
#line 932 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); ;}
    break;

  case 283:
#line 934 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; ;}
    break;

  case 284:
#line 936 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1; ;}
    break;

  case 285:
#line 938 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; ;}
    break;

  case 286:
#line 940 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; ;}
    break;

  case 287:
#line 943 "p.y"
    { (yyval.oinfo).nodeNumber = (yyvsp[(3) - (3)].ival)-1; ;}
    break;

  case 288:
#line 945 "p.y"
    { (yyval.oinfo).surface = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 289:
#line 947 "p.y"
    { (yyval.oinfo).ylayer = (yyvsp[(2) - (3)].fval); (yyval.oinfo).zlayer = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 290:
#line 949 "p.y"
    { (yyval.oinfo).averageFlg = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 291:
#line 951 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 292:
#line 953 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (3)].ival); (yyval.oinfo).ncomplexout = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 293:
#line 955 "p.y"
    { (yyval.oinfo).angularouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 294:
#line 957 "p.y"
    { (yyval.oinfo).rotvecouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 295:
#line 959 "p.y"
    { (yyval.oinfo).rotvecouttype = OutputInfo::Linear; ;}
    break;

  case 296:
#line 961 "p.y"
    { (yyval.oinfo).rescaling = bool((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 297:
#line 963 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 298:
#line 965 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (3)].ival); sfem->setnsamp_out((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 299:
#line 967 "p.y"
    { (yyval.oinfo).oframe = (OutputInfo::FrameType) (yyvsp[(2) - (2)].ival); ;}
    break;

  case 300:
#line 969 "p.y"
    { (yyval.oinfo).matlab = true; ;}
    break;

  case 301:
#line 971 "p.y"
    { domain->solInfo().xmatrixname = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 302:
#line 973 "p.y"
    { domain->solInfo().qmatrixname = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 303:
#line 975 "p.y"
    { domain->solInfo().rmatrixname = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 305:
#line 980 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; ;}
    break;

  case 306:
#line 983 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[(2) - (3)].ival);;}
    break;

  case 307:
#line 986 "p.y"
    { domain->solInfo().qrfactorization = (yyvsp[(2) - (3)].ival);;}
    break;

  case 308:
#line 988 "p.y"
    { domain->solInfo().nEig = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 309:
#line 990 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;;}
    break;

  case 310:
#line 992 "p.y"
    { domain->solInfo().setSubSpaceInfo((yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 311:
#line 994 "p.y"
    { domain->solInfo().subspaceSize = (yyvsp[(2) - (3)].ival);;}
    break;

  case 312:
#line 996 "p.y"
    { domain->solInfo().tolEig = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 313:
#line 998 "p.y"
    { domain->solInfo().tolJac = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 314:
#line 1000 "p.y"
    { domain->solInfo().explicitK = true; ;}
    break;

  case 315:
#line 1002 "p.y"
    { geoSource->setShift((yyvsp[(2) - (3)].fval)); ;}
    break;

  case 316:
#line 1004 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; ;}
    break;

  case 317:
#line 1006 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (3)].strval); ;}
    break;

  case 318:
#line 1009 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (4)].strval); 
          domain->solInfo().arpack_mode = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 319:
#line 1013 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[(2) - (4)].fval), int((yyvsp[(3) - (4)].fval))); ;}
    break;

  case 320:
#line 1016 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival));;}
    break;

  case 321:
#line 1019 "p.y"
    { domain->solInfo().filtereig = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 322:
#line 1021 "p.y"
    { domain->solInfo().eigenSolverSubType = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 323:
#line 1023 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;;}
    break;

  case 324:
#line 1026 "p.y"
    { domain->solInfo().maxitEig = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 325:
#line 1028 "p.y"
    { domain->solInfo().test_ulrich = true; ;}
    break;

  case 326:
#line 1030 "p.y"
    { domain->solInfo().addedMass = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 327:
#line 1034 "p.y"
    { domain->solInfo().sloshing = 1; ;}
    break;

  case 328:
#line 1036 "p.y"
    { domain->setGravitySloshing((yyvsp[(2) - (3)].fval)); ;}
    break;

  case 329:
#line 1040 "p.y"
    { domain->solInfo().massFlag = 1; ;}
    break;

  case 330:
#line 1044 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 331:
#line 1047 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);;}
    break;

  case 332:
#line 1051 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Top); ;}
    break;

  case 333:
#line 1055 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); ;}
    break;

  case 337:
#line 1060 "p.y"
    { domain->solInfo().modal = true; ;}
    break;

  case 338:
#line 1062 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 339:
#line 1064 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 340:
#line 1066 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (8)].ival);
          domain->solInfo().stable_cfl = (yyvsp[(4) - (8)].fval);
          domain->solInfo().stable_tol = (yyvsp[(5) - (8)].fval);
          domain->solInfo().stable_maxit = (yyvsp[(6) - (8)].ival);
          domain->solInfo().stable_freq = (yyvsp[(7) - (8)].ival);
        ;}
    break;

  case 341:
#line 1073 "p.y"
    { domain->solInfo().iacc_switch = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 342:
#line 1075 "p.y"
    { domain->solInfo().zeroRot = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 343:
#line 1077 "p.y"
    { domain->solInfo().no_secondary = true; ;}
    break;

  case 344:
#line 1079 "p.y"
    { domain->solInfo().check_energy_balance = true; ;}
    break;

  case 345:
#line 1081 "p.y"
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[(3) - (5)].fval); 
          domain->solInfo().epsilon2 = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 346:
#line 1087 "p.y"
    { domain->solInfo().ConwepOnOff = true;
          BlastLoading::InputFileData = (yyvsp[(3) - (4)].blastData); ;}
    break;

  case 347:
#line 1092 "p.y"
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

  case 348:
#line 1107 "p.y"
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; ;}
    break;

  case 350:
#line 1110 "p.y"
    { domain->solInfo().acoustic = true; ;}
    break;

  case 352:
#line 1115 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (2)].fval),(yyvsp[(2) - (2)].fval)); ;}
    break;

  case 353:
#line 1117 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (4)].fval),(yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval),(yyvsp[(4) - (4)].fval)); ;}
    break;

  case 354:
#line 1119 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[(1) - (1)].fval)); ;}
    break;

  case 355:
#line 1121 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (3)].fval),(yyvsp[(2) - (3)].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 356:
#line 1127 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
            domain->solInfo().order = 1;
          }
          else 
            domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo((yyvsp[(1) - (1)].fval)); 
        ;}
    break;

  case 357:
#line 1138 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; ;}
    break;

  case 360:
#line 1143 "p.y"
    { domain->solInfo().modal = true; ;}
    break;

  case 361:
#line 1147 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (5)].fval), 0, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 362:
#line 1149 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), 0, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].ival), (yyvsp[(5) - (6)].fval)); ;}
    break;

  case 363:
#line 1153 "p.y"
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].ival)); ;}
    break;

  case 364:
#line 1163 "p.y"
    { domain->solInfo().setAero((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCollocated = 0; ;}
    break;

  case 365:
#line 1166 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().isCollocated = 0; ;}
    break;

  case 366:
#line 1169 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (6)].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[(3) - (6)].ival) < 6 || (yyvsp[(3) - (6)].ival) == 20) {
              domain->solInfo().alphas[0] = (yyvsp[(4) - (6)].fval)+(yyvsp[(5) - (6)].fval);
              domain->solInfo().alphas[1] = -(yyvsp[(5) - (6)].fval);
          }
        ;}
    break;

  case 367:
#line 1177 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (5)].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[(4) - (5)].fval);
        ;}
    break;

  case 368:
#line 1182 "p.y"
    { domain->solInfo().isCollocated = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 369:
#line 1184 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 370:
#line 1186 "p.y"
    {;}
    break;

  case 371:
#line 1190 "p.y"
    {;}
    break;

  case 372:
#line 1192 "p.y"
    { domain->AddAeroEmbedSurfaceId((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 373:
#line 1196 "p.y"
    { domain->solInfo().setAeroHeat((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); ;}
    break;

  case 374:
#line 1200 "p.y"
    { domain->solInfo().setThermoh(1); ;}
    break;

  case 375:
#line 1204 "p.y"
    { domain->solInfo().setThermoe(1); ;}
    break;

  case 376:
#line 1208 "p.y"
    { domain->solInfo().setModeDecomp(1); ;}
    break;

  case 377:
#line 1212 "p.y"
    { domain->solInfo().hzemFlag=1; ;}
    break;

  case 378:
#line 1216 "p.y"
    { domain->solInfo().slzemFlag=1; ;}
    break;

  case 379:
#line 1220 "p.y"
    { domain->solInfo().setTrbm((yyvsp[(3) - (4)].fval)); ;}
    break;

  case 380:
#line 1228 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");;}
    break;

  case 381:
#line 1231 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (4)].fval)); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");;}
    break;

  case 382:
#line 1234 "p.y"
    { domain->solInfo().setGrbm();
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");;}
    break;

  case 383:
#line 1239 "p.y"
    { domain->solInfo().modeFilterFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 384:
#line 1241 "p.y"
    { domain->solInfo().modeFilterFlag = 1; ;}
    break;

  case 385:
#line 1245 "p.y"
    { domain->solInfo().useRbmFilter((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 386:
#line 1247 "p.y"
    { domain->solInfo().useRbmFilter(1); ;}
    break;

  case 388:
#line 1252 "p.y"
    { if((yyvsp[(1) - (1)].ival) < 1 || (yyvsp[(1) - (1)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(1) - (1)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(1) - (1)].ival)-1] = 1;
    ;}
    break;

  case 389:
#line 1260 "p.y"
    { if((yyvsp[(2) - (2)].ival) < 1 || (yyvsp[(2) - (2)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(2) - (2)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(2) - (2)].ival)-1] = 1;
    ;}
    break;

  case 390:
#line 1270 "p.y"
    { domain->solInfo().hzemFilterFlag=1; ;}
    break;

  case 391:
#line 1274 "p.y"
    { domain->solInfo().slzemFilterFlag=1; ;}
    break;

  case 392:
#line 1278 "p.y"
    { domain->solInfo().setTimes((yyvsp[(4) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(2) - (5)].fval)); ;}
    break;

  case 393:
#line 1282 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (5)].ival),(yyvsp[(4) - (5)].ival),1);
        ;}
    break;

  case 394:
#line 1288 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].ival));
        ;}
    break;

  case 395:
#line 1293 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (7)].ival),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        ;}
    break;

  case 398:
#line 1306 "p.y"
    { domain->solInfo().pitaNoForce = true; ;}
    break;

  case 399:
#line 1308 "p.y"
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 400:
#line 1310 "p.y"
    { domain->solInfo().pitaLocalBasisImprovement = 1; ;}
    break;

  case 401:
#line 1312 "p.y"
    { domain->solInfo().pitaTimeReversible = true; ;}
    break;

  case 402:
#line 1314 "p.y"
    { domain->solInfo().pitaRemoteCoarse = true; ;}
    break;

  case 403:
#line 1316 "p.y"
    { domain->solInfo().pitaProjTol = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 404:
#line 1318 "p.y"
    { domain->solInfo().pitaReadInitSeed = true; ;}
    break;

  case 405:
#line 1320 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = 0.0; ;}
    break;

  case 406:
#line 1322 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 407:
#line 1324 "p.y"
    { domain->solInfo().pitaJumpMagnOutput = true; ;}
    break;

  case 408:
#line 1328 "p.y"
    { domain->solInfo().setDamping((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 409:
#line 1330 "p.y"
    { domain->solInfo().setDamping((yyvsp[(2) - (5)].fval),(yyvsp[(3) - (5)].fval));
          domain->solInfo().mtypeDamp = (int)(yyvsp[(4) - (5)].ival); ;}
    break;

  case 410:
#line 1333 "p.y"
    { if(geoSource->setModalDamping((yyvsp[(3) - (3)].bclist)->n, (yyvsp[(3) - (3)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 411:
#line 1338 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); ;}
    break;

  case 412:
#line 1342 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 413:
#line 1348 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval), 0.0);
        ;}
    break;

  case 414:
#line 1356 "p.y"
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 415:
#line 1362 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[(2) - (3)].ival),0.0,0.0,0.0);
        ;}
    break;

  case 417:
#line 1370 "p.y"
    {
           domain->setWaveDirections((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        ;}
    break;

  case 418:
#line 1376 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 419:
#line 1384 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 420:
#line 1386 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 421:
#line 1388 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } ;}
    break;

  case 422:
#line 1390 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } ;}
    break;

  case 423:
#line 1392 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc);
          if(geoSource->getNumSurfaceDirichlet() > 1) delete [] surf_bc; ;}
    break;

  case 425:
#line 1408 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(2) - (6)].ival)-1;
          surf_bc[0].type = (BCond::BCType) (yyvsp[(3) - (6)].ival); //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = (yyvsp[(4) - (6)].ival)-1;
          surf_bc[0].val = (yyvsp[(5) - (6)].ival)-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
          if(geoSource->getNumSurfaceConstraint() > 1) delete [] surf_bc;
        ;}
    break;

  case 426:
#line 1418 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 427:
#line 1422 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 428:
#line 1424 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 429:
#line 1428 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 430:
#line 1430 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (2)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; ;}
    break;

  case 431:
#line 1434 "p.y"
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
        ;}
    break;

  case 432:
#line 1475 "p.y"
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[(1) - (1)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(1) - (1)].bclist)->d)[ii]); ;}
    break;

  case 433:
#line 1479 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); 
          for (int ii = 0; ii < (yyvsp[(2) - (2)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(2) - (2)].bclist)->d)[ii]); ;}
    break;

  case 434:
#line 1485 "p.y"
    { (yyval.bclist) = new BCList;
          for(int i=0; i<(yyvsp[(3) - (4)].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[(3) - (4)].nl).nd[i],10,0.0); } ;}
    break;

  case 435:
#line 1491 "p.y"
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;;}
    break;

  case 436:
#line 1493 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); ;}
    break;

  case 437:
#line 1496 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); ;}
    break;

  case 438:
#line 1505 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 439:
#line 1507 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 440:
#line 1509 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Flux; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); ;}
    break;

  case 441:
#line 1512 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; ;}
    break;

  case 442:
#line 1523 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 443:
#line 1525 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 444:
#line 1527 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Convection; bc.loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add(bc); ;}
    break;

  case 445:
#line 1532 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 446:
#line 1534 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 447:
#line 1536 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Radiation; (yyval.bclist)->add(bc); ;}
    break;

  case 451:
#line 1548 "p.y"
    { domain->addSommer(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); ;}
    break;

  case 452:
#line 1550 "p.y"
    { domain->addSommer(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); ;}
    break;

  case 453:
#line 1552 "p.y"
    { domain->addSommer(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); ;}
    break;

  case 456:
#line 1560 "p.y"
    { domain->addSommerElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        ;}
    break;

  case 457:
#line 1566 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1; ;}
    break;

  case 458:
#line 1568 "p.y"
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++; ;}
    break;

  case 462:
#line 1580 "p.y"
    { domain->addScatter(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); ;}
    break;

  case 463:
#line 1583 "p.y"
    { domain->addScatter(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); ;}
    break;

  case 464:
#line 1586 "p.y"
    { domain->addScatter(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); ;}
    break;

  case 467:
#line 1595 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 470:
#line 1604 "p.y"
    { domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 473:
#line 1612 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; ;}
    break;

  case 476:
#line 1622 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);;}
    break;

  case 477:
#line 1626 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 478:
#line 1630 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 479:
#line 1632 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 480:
#line 1636 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 481:
#line 1640 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 482:
#line 1642 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) { (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Atdneu; } (yyval.bclist) = (yyvsp[(4) - (4)].bclist); ;}
    break;

  case 483:
#line 1646 "p.y"
    { domain->solInfo().ATDARBFlag = (yyvsp[(2) - (3)].fval);;}
    break;

  case 485:
#line 1651 "p.y"
    { domain->solInfo().ATDDNBVal = (yyvsp[(2) - (3)].fval);;}
    break;

  case 487:
#line 1656 "p.y"
    { domain->solInfo().ATDROBVal = (yyvsp[(2) - (5)].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[(3) - (5)].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[(4) - (5)].fval);;}
    break;

  case 489:
#line 1663 "p.y"
    { domain->setFFP((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 490:
#line 1665 "p.y"
    { domain->setFFP((yyvsp[(2) - (4)].ival),(yyvsp[(3) - (4)].ival)); ;}
    break;

  case 491:
#line 1669 "p.y"
    {
           domain->setFFP((yyvsp[(2) - (4)].ival));
        ;}
    break;

  case 492:
#line 1675 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        ;}
    break;

  case 493:
#line 1680 "p.y"
    { fourHelmBC->addDirichlet((yyvsp[(2) - (2)].complexFDBC)); ;}
    break;

  case 494:
#line 1684 "p.y"
    { (yyval.complexFDBC) = FDBC((yyvsp[(1) - (2)].ival)-1); ;}
    break;

  case 495:
#line 1688 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        ;}
    break;

  case 496:
#line 1693 "p.y"
    { fourHelmBC->addNeuman((yyvsp[(2) - (2)].complexFNBC)); ;}
    break;

  case 497:
#line 1697 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1); ;}
    break;

  case 498:
#line 1699 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 499:
#line 1703 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setModes((yyvsp[(2) - (3)].ival));
          domain->solInfo().setProbType(SolverInfo::AxiHelm);
        ;}
    break;

  case 500:
#line 1711 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSlices((yyvsp[(2) - (3)].ival));
        ;}
    break;

  case 501:
#line 1718 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSomType((yyvsp[(3) - (7)].ival));
          fourHelmBC->setSurf((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval));
        ;}
    break;

  case 503:
#line 1726 "p.y"
    { fourHelmBC->addSommer(new LineAxiSommer((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1)); ;}
    break;

  case 504:
#line 1728 "p.y"
    { fourHelmBC->addSommer(new Line2AxiSommer((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1)); ;}
    break;

  case 505:
#line 1732 "p.y"
    { if( globalMPCs== NULL) globalMPCs = new MPCData(); ;}
    break;

  case 506:
#line 1734 "p.y"
    { globalMPCs->addMPC((yyvsp[(2) - (2)].axiMPC)); ;}
    break;

  case 507:
#line 1738 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (5)].ival)-1, (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival), DComplex(1.0,0.0), 0.0, 0.0, 0.0); ;}
    break;

  case 508:
#line 1740 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (7)].ival)-1, (yyvsp[(2) - (7)].fval), (yyvsp[(3) - (7)].fval), (yyvsp[(4) - (7)].ival), DComplex((yyvsp[(5) - (7)].fval),(yyvsp[(6) - (7)].fval)) , 0.0, 0.0, 0.0); ;}
    break;

  case 509:
#line 1742 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (10)].ival)-1, (yyvsp[(2) - (10)].fval), (yyvsp[(3) - (10)].fval), (yyvsp[(4) - (10)].ival), DComplex((yyvsp[(5) - (10)].fval),(yyvsp[(6) - (10)].fval)) , (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)); ;}
    break;

  case 510:
#line 1746 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (3)].strval);
	  domain->solInfo().readmodeCalled = true; ;}
    break;

  case 511:
#line 1749 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (4)].strval);
          domain->solInfo().readmodeCalled = true; 
 	  domain->solInfo().maxSizePodRom = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 512:
#line 1753 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (4)].strval);
          domain->solInfo().readInModes = (yyvsp[(3) - (4)].strval);
          domain->solInfo().readmodeCalled = true; ;}
    break;

  case 513:
#line 1757 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (5)].strval);
          domain->solInfo().readInModes = (yyvsp[(3) - (5)].strval);
          domain->solInfo().readmodeCalled = true;
          domain->solInfo().maxSizePodRom = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 514:
#line 1762 "p.y"
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 515:
#line 1766 "p.y"
    { ;}
    break;

  case 516:
#line 1768 "p.y"
    { domain->solInfo().zeroInitialDisp = 1; ;}
    break;

  case 517:
#line 1770 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 518:
#line 1773 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 519:
#line 1779 "p.y"
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[(2) - (3)].fval);  ;}
    break;

  case 520:
#line 1781 "p.y"
    { (yyval.bclist) = new BCList; amplitude = 1.0; ;}
    break;

  case 521:
#line 1783 "p.y"
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

  case 522:
#line 1794 "p.y"
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

  case 523:
#line 1805 "p.y"
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); ;}
    break;

  case 524:
#line 1808 "p.y"
    { domain->solInfo().buckling = 1; ;}
    break;

  case 525:
#line 1813 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 526:
#line 1815 "p.y"
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

  case 527:
#line 1828 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 528:
#line 1830 "p.y"
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

  case 529:
#line 1842 "p.y"
    { ;}
    break;

  case 530:
#line 1844 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 531:
#line 1847 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 532:
#line 1853 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 533:
#line 1858 "p.y"
    { domain->solInfo().setGEPS();
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Etemperatures;
          if(geoSource->setIDis6((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 534:
#line 1864 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 535:
#line 1866 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 536:
#line 1868 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Forces; (yyvsp[(2) - (2)].bcval).loadsetid = (yyval.bclist)->loadsetid; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 537:
#line 1870 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } ;}
    break;

  case 538:
#line 1872 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Forces, (yyval.bclist)->loadsetid); (yyval.bclist)->add(bc); } ;}
    break;

  case 539:
#line 1874 "p.y"
    { for(int i=(yyvsp[(2) - (8)].ival); i<=(yyvsp[(4) - (8)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (8)].ival)-1, (yyvsp[(6) - (8)].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[(7) - (8)].ival)); (yyval.bclist)->add(bc); } ;}
    break;

  case 540:
#line 1876 "p.y"
    { for(int i=(yyvsp[(2) - (10)].ival); i<=(yyvsp[(4) - (10)].ival); i+=(yyvsp[(6) - (10)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (10)].ival)-1, (yyvsp[(8) - (10)].fval), BCond::Forces, (yyval.bclist)->loadsetid, (BCond::MomentType) (yyvsp[(7) - (10)].ival)); (yyval.bclist)->add(bc); } ;}
    break;

  case 541:
#line 1878 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].loadsetid = (yyval.bclist)->loadsetid;
          geoSource->addSurfaceNeuman(1,surf_bc);
          if(geoSource->getNumSurfaceNeuman() > 1) delete [] surf_bc; ;}
    break;

  case 542:
#line 1887 "p.y"
    { for(int i=0; i<(yyvsp[(5) - (5)].bclist)->n; ++i) (yyvsp[(5) - (5)].bclist)->d[i].type = BCond::Forces;
          if(geoSource->setNeumanModal((yyvsp[(5) - (5)].bclist)->n, (yyvsp[(5) - (5)].bclist)->d) < 0) return -1; ;}
    break;

  case 543:
#line 1892 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 544:
#line 1894 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 545:
#line 1896 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (6)].ival); i<=(yyvsp[(3) - (6)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 546:
#line 1898 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (7)].bclist); for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 547:
#line 1900 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (8)].ival); i<=(yyvsp[(3) - (8)].ival); i+=(yyvsp[(5) - (8)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(6) - (8)].ival)-1, (yyvsp[(7) - (8)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 548:
#line 1902 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (9)].bclist); for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 549:
#line 1906 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 550:
#line 1908 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 551:
#line 1912 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 552:
#line 1914 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 555:
#line 1922 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ymtt)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addYMTT((yyval.ymtt));;}
    break;

  case 556:
#line 1924 "p.y"
    { (yyval.ymtt)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 557:
#line 1926 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ymtt)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addYMTT((yyval.ymtt));;}
    break;

  case 560:
#line 1934 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ctett)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addCTETT((yyval.ctett));;}
    break;

  case 561:
#line 1936 "p.y"
    { (yyval.ctett)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 562:
#line 1938 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ctett)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addCTETT((yyval.ctett));;}
    break;

  case 565:
#line 1946 "p.y"
    { (yyval.sdetaft) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.sdetaft)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addSDETAFT((yyval.sdetaft));;}
    break;

  case 566:
#line 1948 "p.y"
    { (yyval.sdetaft)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 567:
#line 1950 "p.y"
    { (yyval.sdetaft) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.sdetaft)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addSDETAFT((yyval.sdetaft));;}
    break;

  case 570:
#line 1958 "p.y"
    { (yyval.lmpcons) = (yyvsp[(1) - (2)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 571:
#line 1962 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); ;}
    break;

  case 572:
#line 1964 "p.y"
    { (yyval.lmpcons) = (yyvsp[(2) - (3)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(3) - (3)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 573:
#line 1970 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 574:
#line 1973 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (3)].ival), (yyvsp[(2) - (3)].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 575:
#line 1976 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval));
          (yyval.lmpcons)->type = (yyvsp[(4) - (5)].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 576:
#line 1980 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (4)].ival), (yyvsp[(2) - (4)].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[(3) - (4)].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 577:
#line 1985 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (6)].ival), (yyvsp[(2) - (6)].fval));
          (yyval.lmpcons)->type = (yyvsp[(4) - (6)].ival);
          (yyval.lmpcons)->lagrangeMult = (yyvsp[(5) - (6)].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[(5) - (6)].copt).penalty;
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 578:
#line 1993 "p.y"
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

  case 581:
#line 2009 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].cxbcval).nnum,(yyvsp[(1) - (2)].cxbcval).reval,(yyvsp[(1) - (2)].cxbcval).imval,(yyvsp[(2) - (2)].mpcterm)); domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 582:
#line 2011 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); ;}
    break;

  case 583:
#line 2013 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(2) - (3)].cxbcval).nnum,(yyvsp[(2) - (3)].cxbcval).reval,(yyvsp[(2) - (3)].cxbcval).imval,(yyvsp[(3) - (3)].mpcterm)); domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 584:
#line 2017 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (5)].ival); (yyval.cxbcval).reval=(yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval=(yyvsp[(4) - (5)].fval); ;}
    break;

  case 585:
#line 2019 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (3)].ival); (yyval.cxbcval).reval=(yyvsp[(2) - (3)].fval); (yyval.cxbcval).imval=0.0; ;}
    break;

  case 586:
#line 2021 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (2)].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; ;}
    break;

  case 587:
#line 2025 "p.y"
    { if(((yyvsp[(3) - (5)].fval)==0.0) && ((yyvsp[(4) - (5)].fval)==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (5)].ival),(yyvsp[(2) - (5)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (5)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (5)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval)); }
        ;}
    break;

  case 588:
#line 2033 "p.y"
    { if((yyvsp[(3) - (4)].fval)==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (4)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (4)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (4)].fval),0.0); }
        ;}
    break;

  case 589:
#line 2043 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); ;}
    break;

  case 590:
#line 2045 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].cxbclist)->n; ++i) (yyvsp[(4) - (4)].cxbclist)->d[i].loadsetid = (yyvsp[(2) - (4)].ival);
          (yyval.cxbclist) = (yyvsp[(4) - (4)].cxbclist); ;}
    break;

  case 591:
#line 2050 "p.y"
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[(1) - (1)].cxbcval)); ;}
    break;

  case 592:
#line 2052 "p.y"
    { (yyval.cxbclist) = (yyvsp[(1) - (2)].cxbclist); (yyval.cxbclist)->add((yyvsp[(2) - (2)].cxbcval)); ;}
    break;

  case 595:
#line 2060 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval); 
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.Iyy = (yyvsp[(14) - (16)].fval); sp.Izz = (yyvsp[(15) - (16)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 596:
#line 2068 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.betaDamp = (yyvsp[(17) - (19)].fval); sp.alphaDamp = (yyvsp[(18) - (19)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        ;}
    break;

  case 597:
#line 2077 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.etaDamp = (yyvsp[(17) - (19)].fval); sp.betaDamp = (yyvsp[(18) - (19)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        ;}
    break;

  case 598:
#line 2086 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (20)].fval);  sp.E = (yyvsp[(3) - (20)].fval);  sp.nu  = (yyvsp[(4) - (20)].fval);  sp.rho = (yyvsp[(5) - (20)].fval);
          sp.c = (yyvsp[(6) - (20)].fval);  sp.k = (yyvsp[(7) - (20)].fval);  sp.eh  = (yyvsp[(8) - (20)].fval);  sp.P   = (yyvsp[(9) - (20)].fval);  sp.Ta  = (yyvsp[(10) - (20)].fval); 
          sp.Q = (yyvsp[(11) - (20)].fval); sp.W = (yyvsp[(12) - (20)].fval); sp.Ixx = (yyvsp[(13) - (20)].fval); sp.Iyy = (yyvsp[(14) - (20)].fval); sp.Izz = (yyvsp[(15) - (20)].fval);
	  sp.ymin = (yyvsp[(16) - (20)].fval); sp.ymax = (yyvsp[(17) - (20)].fval); sp.zmin = (yyvsp[(18) - (20)].fval); sp.zmax = (yyvsp[(19) - (20)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (20)].ival)-1, sp );
        ;}
    break;

  case 599:
#line 2095 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.betaDamp = (yyvsp[(21) - (23)].fval); sp.alphaDamp = (yyvsp[(22) - (23)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        ;}
    break;

  case 600:
#line 2105 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.etaDamp = (yyvsp[(21) - (23)].fval); sp.betaDamp = (yyvsp[(22) - (23)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        ;}
    break;

  case 601:
#line 2115 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (9)].fval); sp.E = (yyvsp[(3) - (9)].fval); sp.nu = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(5) - (9)].fval);
          sp.c = (yyvsp[(6) - (9)].fval); sp.k = (yyvsp[(7) - (9)].fval); sp.eh = (yyvsp[(8) - (9)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp ); 
        ;}
    break;

  case 602:
#line 2122 "p.y"
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (14)].fval);  sp.E = (yyvsp[(3) - (14)].fval);  sp.nu  = (yyvsp[(4) - (14)].fval);  sp.rho = (yyvsp[(5) - (14)].fval);
          sp.c = (yyvsp[(6) - (14)].fval);  sp.k = (yyvsp[(7) - (14)].fval);  sp.eh  = (yyvsp[(8) - (14)].fval);  sp.P   = (yyvsp[(9) - (14)].fval);  sp.Ta  = (yyvsp[(10) - (14)].fval);
          sp.Q = (yyvsp[(11) - (14)].fval); sp.W = (yyvsp[(12) - (14)].fval); sp.Ixx = (yyvsp[(13) - (14)].fval);  
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
        ;}
    break;

  case 603:
#line 2130 "p.y"
    { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval);
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.betaDamp = (yyvsp[(15) - (16)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 604:
#line 2139 "p.y"
    { StructProp sp; // this is used for the reduced mesh file output in Rom.d/MeshOutput.C
                         // all properties relevant to structural nonlinear dynamics should be included
          sp.A = (yyvsp[(2) - (30)].fval);  sp.E = (yyvsp[(3) - (30)].fval);  sp.nu  = (yyvsp[(4) - (30)].fval);  sp.rho = (yyvsp[(5) - (30)].fval);
          sp.c = (yyvsp[(6) - (30)].fval);  sp.k = (yyvsp[(7) - (30)].fval);  sp.eh  = (yyvsp[(8) - (30)].fval);  sp.P   = (yyvsp[(9) - (30)].fval);  sp.Ta  = (yyvsp[(10) - (30)].fval);
          sp.Q = (yyvsp[(11) - (30)].fval); sp.W = (yyvsp[(12) - (30)].fval); sp.Ixx = (yyvsp[(13) - (30)].fval); sp.Iyy = (yyvsp[(14) - (30)].fval); sp.Izz = (yyvsp[(15) - (30)].fval);
          sp.ymin = (yyvsp[(16) - (30)].fval); sp.ymax = (yyvsp[(17) - (30)].fval); sp.zmin = (yyvsp[(18) - (30)].fval); sp.zmax = (yyvsp[(19) - (30)].fval);
          sp.betaDamp = (yyvsp[(20) - (30)].fval); sp.alphaDamp = (yyvsp[(21) - (30)].fval); 
          sp.lagrangeMult = bool((yyvsp[(22) - (30)].ival)); sp.penalty = (yyvsp[(23) - (30)].fval); sp.initialPenalty = (yyvsp[(24) - (30)].fval);
          sp.funtype = (yyvsp[(25) - (30)].ival); sp.type = StructProp::PropType((yyvsp[(26) - (30)].ival)); sp.k1 = (yyvsp[(27) - (30)].fval); sp.k2 = (yyvsp[(28) - (30)].fval); sp.k3 = (yyvsp[(29) - (30)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (30)].ival)-1, sp );
        ;}
    break;

  case 605:
#line 2152 "p.y"
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
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 606:
#line 2169 "p.y"
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
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (13)].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 607:
#line 2187 "p.y"
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
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
          domain->PMLFlag = 1;
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 608:
#line 2205 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (5)].fval),0.0);
          sp.rho = (yyvsp[(4) - (5)].fval);
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 609:
#line 2214 "p.y"
    { StructProp sp;
          sp.soundSpeed = complex<double>((yyvsp[(3) - (6)].fval),(yyvsp[(4) - (6)].fval));
          sp.rho = (yyvsp[(5) - (6)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
          domain->solInfo().acoustic = true;
        ;}
    break;

  case 610:
#line 2222 "p.y"
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
        ;}
    break;

  case 611:
#line 2241 "p.y"
    { StructProp sp; 
          sp.A = (yyvsp[(3) - (12)].fval); sp.rho = (yyvsp[(4) - (12)].fval); sp.Q = (yyvsp[(5) - (12)].fval); sp.c = (yyvsp[(6) - (12)].fval); 
          sp.sigma = (yyvsp[(7) - (12)].fval); sp.k = (yyvsp[(8) - (12)].fval); sp.eh = (yyvsp[(9) - (12)].fval); sp.P = (yyvsp[(10) - (12)].fval); sp.Ta = (yyvsp[(11) - (12)].fval);
          sp.isReal = true;
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 612:
#line 2249 "p.y"
    { // rigid element or joint with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          geoSource->addMat( (yyvsp[(1) - (3)].ival)-1, sp );
        ;}
    break;

  case 613:
#line 2255 "p.y"
    { // rigid element or joint
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (4)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (4)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (4)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (4)].ival)-1, sp );
        ;}
    break;

  case 614:
#line 2265 "p.y"
    { // rigid solid element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[(4) - (5)].fval);
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        ;}
    break;

  case 615:
#line 2272 "p.y"
    { // rigid solid element with mass
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

  case 616:
#line 2283 "p.y"
    { // rigid beam or shell element with mass, and default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.rho = (yyvsp[(4) - (6)].fval);
          sp.A = sp.eh = (yyvsp[(5) - (6)].fval);
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 617:
#line 2291 "p.y"
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

  case 618:
#line 2303 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 619:
#line 2317 "p.y"
    { // joint-with-driver, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (8)].ival);
          sp.amplitude = (yyvsp[(4) - (8)].fval);
          sp.offset = (yyvsp[(5) - (8)].fval);
          sp.c1 = (yyvsp[(6) - (8)].fval);
          sp.c2 = (yyvsp[(7) - (8)].fval);
          sp.type = StructProp::Undefined;
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        ;}
    break;

  case 620:
#line 2328 "p.y"
    { // joint-with-driver, 3-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (9)].ival);
          sp.amplitude = (yyvsp[(4) - (9)].fval);
          sp.offset = (yyvsp[(5) - (9)].fval);
          sp.c1 = (yyvsp[(6) - (9)].fval);
          sp.c2 = (yyvsp[(7) - (9)].fval);
          sp.c3 = (yyvsp[(8) - (9)].fval);
          sp.type = StructProp::Undefined;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        ;}
    break;

  case 621:
#line 2340 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 622:
#line 2353 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        ;}
    break;

  case 623:
#line 2368 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 624:
#line 2384 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 625:
#line 2401 "p.y"
    { // actuated joint, 2-parameter elementary function, and default constraint options
          StructProp sp;
          sp.funtype = (yyvsp[(3) - (10)].ival);
          sp.amplitude = (yyvsp[(4) - (10)].fval);
          sp.offset = (yyvsp[(5) - (10)].fval);
          sp.c1 = (yyvsp[(6) - (10)].fval);
          sp.c2 = (yyvsp[(7) - (10)].fval);
          sp.k1 = (yyvsp[(9) - (10)].fval);
          sp.type = StructProp::Undefined;
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 626:
#line 2413 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 627:
#line 2426 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 628:
#line 2440 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 629:
#line 2456 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 630:
#line 2473 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (13)].ival)-1, sp );
        ;}
    break;

  case 631:
#line 2491 "p.y"
    { // RevoluteJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[(4) - (5)].fval);
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        ;}
    break;

  case 632:
#line 2498 "p.y"
    { // RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (6)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (6)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (6)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (6)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(5) - (6)].fval);
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 633:
#line 2509 "p.y"
    { // UniversalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[(4) - (6)].fval);
          sp.k2 = (yyvsp[(5) - (6)].fval);
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 634:
#line 2517 "p.y"
    { // UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = (yyvsp[(3) - (7)].copt).lagrangeMult;
          sp.initialPenalty = sp.penalty = (yyvsp[(3) - (7)].copt).penalty;
          sp.constraint_hess = (yyvsp[(3) - (7)].copt).constraint_hess;
          sp.constraint_hess_eps = (yyvsp[(3) - (7)].copt).constraint_hess_eps;
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(5) - (7)].fval);
          sp.k2 = (yyvsp[(6) - (7)].fval);
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 635:
#line 2529 "p.y"
    { // SphericalJointSpringCombo with default constraint options
          StructProp sp;
          sp.type = StructProp::Undefined;
          sp.k1 = (yyvsp[(4) - (7)].fval);
          sp.k2 = (yyvsp[(5) - (7)].fval);
          sp.k3 = (yyvsp[(6) - (7)].fval);
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 636:
#line 2538 "p.y"
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
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        ;}
    break;

  case 637:
#line 2551 "p.y"
    { // TorsionalSpringType1 or TranslationalSpring
          StructProp sp;
          sp.k1 = (yyvsp[(3) - (4)].fval);
          geoSource->addMat( (yyvsp[(1) - (4)].ival)-1, sp );
        ;}
    break;

  case 640:
#line 2563 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (3)].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 641:
#line 2569 "p.y"
    { if((yyvsp[(2) - (4)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (4)].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 642:
#line 2575 "p.y"
    { if((yyvsp[(2) - (5)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (5)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (5)].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 643:
#line 2582 "p.y"
    { if((yyvsp[(2) - (6)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (6)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (6)].fval));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 644:
#line 2590 "p.y"
    { if((yyval.SurfObj)->GetReverseNormals()) { // reverse the node numbering
            int *nodes = new int[(yyvsp[(4) - (5)].nl).num];
            for(int i=0; i<(yyvsp[(4) - (5)].nl).num; ++i) nodes[(yyvsp[(4) - (5)].nl).num-1-i] = (yyvsp[(4) - (5)].nl).nd[i];
            (yyval.SurfObj)->AddFaceElement((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].nl).num, nodes);
            delete [] nodes;
          }
          else (yyval.SurfObj)->AddFaceElement((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].nl).num, (yyvsp[(4) - (5)].nl).nd);
        ;}
    break;

  case 645:
#line 2601 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 646:
#line 2603 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 647:
#line 2605 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 648:
#line 2609 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (6)].ival), (yyvsp[(3) - (6)].ival), (yyvsp[(5) - (6)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 649:
#line 2611 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 650:
#line 2613 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 651:
#line 2615 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 652:
#line 2617 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 653:
#line 2621 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 654:
#line 2627 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->solInfo().isCoupled = true; ;}
    break;

  case 655:
#line 2629 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (3)].ival), (yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; ;}
    break;

  case 656:
#line 2637 "p.y"
    { ;}
    break;

  case 657:
#line 2639 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 658:
#line 2645 "p.y"
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 659:
#line 2652 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 660:
#line 2659 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 661:
#line 2666 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 662:
#line 2674 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 663:
#line 2683 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 664:
#line 2690 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 665:
#line 2700 "p.y"
    { ;}
    break;

  case 666:
#line 2702 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 667:
#line 2708 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 668:
#line 2716 "p.y"
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; ;}
    break;

  case 670:
#line 2722 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; ;}
    break;

  case 671:
#line 2728 "p.y"
    { ;}
    break;

  case 672:
#line 2730 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 673:
#line 2736 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 674:
#line 2746 "p.y"
    { ;}
    break;

  case 675:
#line 2748 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 676:
#line 2755 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 677:
#line 2762 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 678:
#line 2769 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 679:
#line 2776 "p.y"
    { /* this one is for frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 680:
#line 2784 "p.y"
    { /* this one is for constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 681:
#line 2793 "p.y"
    { /* this one is for velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (13)].ival), (yyvsp[(4) - (13)].ival), (yyvsp[(6) - (13)].fval), (yyvsp[(7) - (13)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (13)].ival), (yyvsp[(9) - (13)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (13)].fval), (yyvsp[(11) - (13)].fval), (yyvsp[(12) - (13)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (13)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 682:
#line 2802 "p.y"
    { /* this one is for pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (14)].ival), (yyvsp[(4) - (14)].ival), (yyvsp[(6) - (14)].fval), (yyvsp[(7) - (14)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (14)].ival), (yyvsp[(9) - (14)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (14)].fval), (yyvsp[(11) - (14)].fval), (yyvsp[(12) - (14)].fval), (yyvsp[(13) - (14)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (14)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 683:
#line 2811 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 684:
#line 2819 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival)); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 685:
#line 2829 "p.y"
    { domain->solInfo().dist_acme = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 686:
#line 2831 "p.y"
    { domain->solInfo().ffi_debug = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 687:
#line 2833 "p.y"
    { domain->solInfo().mortar_scaling = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 688:
#line 2835 "p.y"
    { domain->solInfo().mortar_integration_rule = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 689:
#line 2838 "p.y"
    { geoSource->addNode((yyvsp[(3) - (3)].nval).num, (yyvsp[(3) - (3)].nval).xyz, (yyvsp[(3) - (3)].nval).cp, (yyvsp[(3) - (3)].nval).cd); ;}
    break;

  case 690:
#line 2840 "p.y"
    { geoSource->addNode((yyvsp[(2) - (2)].nval).num, (yyvsp[(2) - (2)].nval).xyz, (yyvsp[(2) - (2)].nval).cp, (yyvsp[(2) - (2)].nval).cd); ;}
    break;

  case 691:
#line 2844 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (5)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (5)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (5)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (5)].fval);  (yyval.nval).cp = 0;  (yyval.nval).cd = 0; ;}
    break;

  case 692:
#line 2846 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (4)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (4)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (4)].fval);  (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; ;}
    break;

  case 693:
#line 2848 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (3)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (3)].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; (yyval.nval).cp = 0;  (yyval.nval).cd = 0; ;}
    break;

  case 694:
#line 2850 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (7)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (7)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (7)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (7)].fval);  (yyval.nval).cp = (yyvsp[(5) - (7)].ival); (yyval.nval).cd = (yyvsp[(6) - (7)].ival);
          if((yyvsp[(5) - (7)].ival) != 0) domain->solInfo().basicPosCoords = false;
          if((yyvsp[(6) - (7)].ival) != 0) domain->solInfo().basicDofCoords = false; ;}
    break;

  case 695:
#line 2854 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (6)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (6)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (6)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (6)].fval);  (yyval.nval).cp = (yyvsp[(5) - (6)].ival); (yyval.nval).cd = (yyvsp[(5) - (6)].ival);
          if((yyvsp[(5) - (6)].ival) != 0) { domain->solInfo().basicPosCoords = false; domain->solInfo().basicDofCoords = false; } ;}
    break;

  case 696:
#line 2859 "p.y"
    { /* Define each Element */
          geoSource->addElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);;}
    break;

  case 697:
#line 2864 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1;;}
    break;

  case 698:
#line 2866 "p.y"
    { if((yyval.nl).num == 125) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++;;}
    break;

  case 699:
#line 2871 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.bcval).val = (yyvsp[(3) - (4)].fval); (yyval.bcval).mtype = BCond::Axial; ;}
    break;

  case 700:
#line 2873 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (3)].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = BCond::Axial; ;}
    break;

  case 701:
#line 2875 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (5)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (5)].ival)-1; (yyval.bcval).val = (yyvsp[(3) - (5)].fval); (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[(4) - (5)].ival); ;}
    break;

  case 702:
#line 2877 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.bcval).val = 0.0; (yyval.bcval).mtype = (BCond::MomentType) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 703:
#line 2881 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 704:
#line 2885 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 705:
#line 2889 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (5)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (5)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval = (yyvsp[(4) - (5)].fval);  ;}
    break;

  case 706:
#line 2891 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (4)].fval); (yyval.cxbcval).imval = 0.0; ;}
    break;

  case 708:
#line 2896 "p.y"
    { geoSource->setCSFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 710:
#line 2901 "p.y"
    { geoSource->setFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 711:
#line 2905 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (11)].ival)-1; 
          (yyval.frame).d[0] = (yyvsp[(2) - (11)].fval); (yyval.frame).d[1] = (yyvsp[(3) - (11)].fval); (yyval.frame).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.frame).d[3] = (yyvsp[(5) - (11)].fval); (yyval.frame).d[4] = (yyvsp[(6) - (11)].fval); (yyval.frame).d[5] = (yyvsp[(7) - (11)].fval);
          (yyval.frame).d[6] = (yyvsp[(8) - (11)].fval); (yyval.frame).d[7] = (yyvsp[(9) - (11)].fval); (yyval.frame).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 712:
#line 2910 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (4)].ival)-1;
          geoSource->makeEframe((yyvsp[(1) - (4)].ival)-1, (yyvsp[(3) - (4)].ival), (yyval.frame).d); ;}
    break;

  case 714:
#line 2916 "p.y"
    { geoSource->setNodalFrame((yyvsp[(2) - (2)].nframe).id,(yyvsp[(2) - (2)].nframe).o,(yyvsp[(2) - (2)].nframe).d,(yyvsp[(2) - (2)].nframe).type); ;}
    break;

  case 715:
#line 2920 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (11)].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;  (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[(2) - (11)].fval); (yyval.nframe).d[1] = (yyvsp[(3) - (11)].fval); (yyval.nframe).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.nframe).d[3] = (yyvsp[(5) - (11)].fval); (yyval.nframe).d[4] = (yyvsp[(6) - (11)].fval); (yyval.nframe).d[5] = (yyvsp[(7) - (11)].fval);
          (yyval.nframe).d[6] = (yyvsp[(8) - (11)].fval); (yyval.nframe).d[7] = (yyvsp[(9) - (11)].fval); (yyval.nframe).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 716:
#line 2927 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (14)].ival);
          (yyval.nframe).type = NFrameData::Rectangular;
          (yyval.nframe).o[0] = (yyvsp[(2) - (14)].fval);  (yyval.nframe).o[1] = (yyvsp[(3) - (14)].fval);  (yyval.nframe).o[2] = (yyvsp[(4) - (14)].fval);
          (yyval.nframe).d[0] = (yyvsp[(5) - (14)].fval);  (yyval.nframe).d[1] = (yyvsp[(6) - (14)].fval);  (yyval.nframe).d[2] = (yyvsp[(7) - (14)].fval);
          (yyval.nframe).d[3] = (yyvsp[(8) - (14)].fval);  (yyval.nframe).d[4] = (yyvsp[(9) - (14)].fval);  (yyval.nframe).d[5] = (yyvsp[(10) - (14)].fval);
          (yyval.nframe).d[6] = (yyvsp[(11) - (14)].fval); (yyval.nframe).d[7] = (yyvsp[(12) - (14)].fval); (yyval.nframe).d[8] = (yyvsp[(13) - (14)].fval); ;}
    break;

  case 717:
#line 2934 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (12)].ival);
          (yyval.nframe).type = (yyvsp[(2) - (12)].ival);
          (yyval.nframe).o[0] = 0;  (yyval.nframe).o[1] = 0;   (yyval.nframe).o[2] = 0;
          (yyval.nframe).d[0] = (yyvsp[(3) - (12)].fval); (yyval.nframe).d[1] = (yyvsp[(4) - (12)].fval);  (yyval.nframe).d[2] = (yyvsp[(5) - (12)].fval);
          (yyval.nframe).d[3] = (yyvsp[(6) - (12)].fval); (yyval.nframe).d[4] = (yyvsp[(7) - (12)].fval);  (yyval.nframe).d[5] = (yyvsp[(8) - (12)].fval);
          (yyval.nframe).d[6] = (yyvsp[(9) - (12)].fval); (yyval.nframe).d[7] = (yyvsp[(10) - (12)].fval); (yyval.nframe).d[8] = (yyvsp[(11) - (12)].fval); ;}
    break;

  case 718:
#line 2941 "p.y"
    { (yyval.nframe).id = (yyvsp[(1) - (15)].ival);
          (yyval.nframe).type = (yyvsp[(2) - (15)].ival);
          (yyval.nframe).o[0] = (yyvsp[(3) - (15)].fval);  (yyval.nframe).o[1] = (yyvsp[(4) - (15)].fval);  (yyval.nframe).o[2] = (yyvsp[(5) - (15)].fval);
          (yyval.nframe).d[0] = (yyvsp[(6) - (15)].fval);  (yyval.nframe).d[1] = (yyvsp[(7) - (15)].fval);  (yyval.nframe).d[2] = (yyvsp[(8) - (15)].fval);
          (yyval.nframe).d[3] = (yyvsp[(9) - (15)].fval);  (yyval.nframe).d[4] = (yyvsp[(10) - (15)].fval); (yyval.nframe).d[5] = (yyvsp[(11) - (15)].fval);
          (yyval.nframe).d[6] = (yyvsp[(12) - (15)].fval); (yyval.nframe).d[7] = (yyvsp[(13) - (15)].fval); (yyval.nframe).d[8] = (yyvsp[(14) - (15)].fval); ;}
    break;

  case 720:
#line 2951 "p.y"
    { OffsetData od;
	  od.first = (yyvsp[(2) - (7)].ival)-1; od.last = (yyvsp[(3) - (7)].ival)-1;
	  od.o[0] = (yyvsp[(4) - (7)].fval); od.o[1] = (yyvsp[(5) - (7)].fval); od.o[2] = (yyvsp[(6) - (7)].fval); 
	  geoSource->addOffset(od); ;}
    break;

  case 721:
#line 2958 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 722:
#line 2960 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 723:
#line 2962 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1); 
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (6)].ival) - 1, (yyvsp[(5) - (6)].fval));
	  domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 724:
#line 2966 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1);
          geoSource->setElementLumpingWeight((yyvsp[(2) - (7)].ival) - 1, (yyvsp[(6) - (7)].fval));
          domain->solInfo().elemLumpPodRom = true; 
          domain->solInfo().reduceFollower = true;;}
    break;

  case 725:
#line 2971 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(4) - (6)].ival)-1,(yyvsp[(5) - (6)].ival)-1); ;}
    break;

  case 726:
#line 2973 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (8)].ival)-1,(yyvsp[(3) - (8)].ival)-1,(yyvsp[(4) - (8)].ival)-1,(yyvsp[(5) - (8)].ival)-1);
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (8)].ival) - 1, (yyvsp[(7) - (8)].fval)); 
	  domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 727:
#line 2977 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1,(yyvsp[(4) - (7)].ival)-1,-1,(yyvsp[(6) - (7)].fval)); ;}
    break;

  case 728:
#line 2979 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        ;}
    break;

  case 729:
#line 2984 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
 	    geoSource->setAttrib(i-1,(yyvsp[(4) - (5)].ival)-1);
	;}
    break;

  case 730:
#line 2989 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (7)].ival); i<(yyvsp[(3) - (7)].ival)+1; ++i)
	    geoSource->setAttrib(i-1, (yyvsp[(4) - (7)].ival)-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].ival)-1);
	;}
    break;

  case 731:
#line 2994 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (8)].ival); i<(yyvsp[(3) - (8)].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[(4) - (8)].ival)-1, (yyvsp[(5) - (8)].ival)-1, -1, (yyvsp[(7) - (8)].fval));
        ;}
    break;

  case 732:
#line 3001 "p.y"
    { domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 733:
#line 3003 "p.y"
    { geoSource->setElementLumpingWeight((yyvsp[(2) - (4)].ival) - 1, (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 734:
#line 3005 "p.y"
    { domain->solInfo().reduceFollower = true;;}
    break;

  case 735:
#line 3010 "p.y"
    { domain->solInfo().ReducedStiffness = true;;}
    break;

  case 736:
#line 3012 "p.y"
    { geoSource->pushBackStiffVec((yyvsp[(2) - (3)].fval));;}
    break;

  case 738:
#line 3017 "p.y"
    { domain->solInfo().forcePodSize = (yyvsp[(2) - (4)].ival);
          domain->solInfo().maxDeimBasisSize = (yyvsp[(3) - (4)].ival);;}
    break;

  case 739:
#line 3020 "p.y"
    { geoSource->pushBackUDEIMVec((yyvsp[(2) - (3)].fval));;}
    break;

  case 741:
#line 3025 "p.y"
    { domain->solInfo().DEIMPodRom = true;
          geoSource->setSampleNodesAndSlots((yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].ival));;}
    break;

  case 742:
#line 3028 "p.y"
    { geoSource->setSampleElemsAndDOFs((yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].ival));
          domain->solInfo().UDEIMPodRom = true;;}
    break;

  case 743:
#line 3034 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 744:
#line 3036 "p.y"
    { (yyval.ival) = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 745:
#line 3038 "p.y"
    { PressureBCond pbc;
          pbc.setData((yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval), (yyval.ival), true);
          geoSource->setElementPressure(pbc); ;}
    break;

  case 746:
#line 3042 "p.y"
    { for(int i = (yyvsp[(2) - (5)].ival); i < ((yyvsp[(3) - (5)].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[(4) - (5)].fval), (yyval.ival), true);
            geoSource->setElementPressure(pbc);
          } ;}
    break;

  case 747:
#line 3048 "p.y"
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].fval), (yyval.ival), true);
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; ;}
    break;

  case 748:
#line 3053 "p.y"
    { PressureBCond pbc;
          pbc.setData((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].fval), (yyval.ival), (yyvsp[(4) - (5)].ival));
          geoSource->setElementPressure(pbc); ;}
    break;

  case 749:
#line 3057 "p.y"
    { for(int i = (yyvsp[(2) - (6)].ival); i < ((yyvsp[(3) - (6)].ival)+1); ++i) {
            PressureBCond pbc;
            pbc.setData(i-1, (yyvsp[(4) - (6)].fval), (yyval.ival), (yyvsp[(5) - (6)].ival));
            geoSource->setElementPressure(pbc);
          } ;}
    break;

  case 750:
#line 3063 "p.y"
    { PressureBCond *pbc = new PressureBCond[1];
          pbc[0].setData((yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].fval), (yyval.ival), (yyvsp[(5) - (6)].ival));
          geoSource->addSurfacePressure(1, pbc);
          if(geoSource->getNumSurfacePressure() > 1) delete [] pbc; ;}
    break;

  case 751:
#line 3070 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        ;}
    break;

  case 752:
#line 3075 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[(2) - (3)].ival));
          geoSource->setConsistentPFlag(false);
        ;}
    break;

  case 753:
#line 3080 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[(2) - (4)].ival));
          geoSource->setConsistentPFlag(false);
          domain->solInfo().inertiaLumping = (yyvsp[(3) - (4)].ival);
        ;}
    break;

  case 754:
#line 3088 "p.y"
    { ;}
    break;

  case 755:
#line 3090 "p.y"
    { geoSource->setElementPreLoad( (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval) ); ;}
    break;

  case 756:
#line 3092 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (6)].ival); i<((yyvsp[(4) - (6)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[(5) - (6)].fval) );
        ;}
    break;

  case 757:
#line 3097 "p.y"
    { double load[3] = { (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval) };
          geoSource->setElementPreLoad( (yyvsp[(2) - (6)].ival)-1, load ); ;}
    break;

  case 758:
#line 3100 "p.y"
    { double load[3] = { (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
          int i;
          for(i=(yyvsp[(2) - (8)].ival); i<((yyvsp[(4) - (8)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        ;}
    break;

  case 759:
#line 3108 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 763:
#line 3113 "p.y"
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

  case 764:
#line 3130 "p.y"
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

  case 765:
#line 3151 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(1) - (1)].ival)); ;}
    break;

  case 766:
#line 3153 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 767:
#line 3157 "p.y"
    { domain->solInfo().type = 1;
          domain->solInfo().iterType = (yyvsp[(1) - (2)].ival); ;}
    break;

  case 768:
#line 3160 "p.y"
    { domain->solInfo().precond = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 769:
#line 3162 "p.y"
    { domain->solInfo().maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 770:
#line 3164 "p.y"
    { domain->solInfo().tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 771:
#line 3166 "p.y"
    { domain->solInfo().maxvecsize = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 772:
#line 3168 "p.y"
    { domain->solInfo().iterSubtype = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 773:
#line 3172 "p.y"
    { domain->solInfo().setSolver(0); ;}
    break;

  case 774:
#line 3174 "p.y"
    { domain->solInfo().setSolver((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 775:
#line 3176 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (2)].ival)); ;}
    break;

  case 776:
#line 3178 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (3)].ival));
          if((yyvsp[(1) - (3)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 777:
#line 3182 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (3)].ival));
          domain->solInfo().getNLInfo().unsymmetric = true; ;}
    break;

  case 778:
#line 3185 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (3)].ival),(yyvsp[(2) - (3)].ival)); ;}
    break;

  case 779:
#line 3187 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 780:
#line 3189 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (5)].ival),(yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].ival)); ;}
    break;

  case 781:
#line 3191 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (6)].ival),(yyvsp[(2) - (6)].ival),(yyvsp[(3) - (6)].fval),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].ival)); ;}
    break;

  case 782:
#line 3193 "p.y"
    { domain->solInfo().setSolver((yyvsp[(1) - (7)].ival),(yyvsp[(2) - (7)].ival),(yyvsp[(3) - (7)].fval),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].ival),(yyvsp[(6) - (7)].ival)); ;}
    break;

  case 783:
#line 3195 "p.y"
    { domain->solInfo().fetiInfo.maxit    = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(3) - (4)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(2) - (4)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 784:
#line 3200 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (3)].ival)-1); ;}
    break;

  case 785:
#line 3203 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; ;}
    break;

  case 786:
#line 3207 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; ;}
    break;

  case 787:
#line 3212 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (4)].ival)-1); 
          domain->solInfo().fetiInfo.feti2version 
                  = (FetiInfo::Feti2Version) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 788:
#line 3217 "p.y"
    { domain->solInfo().fetiInfo.maxit    = (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(3) - (5)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(4) - (5)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 789:
#line 3222 "p.y"
    { domain->solInfo().type =(2); ;}
    break;

  case 790:
#line 3224 "p.y"
    { domain->solInfo().type = 3;
          domain->solInfo().subtype = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getFetiInfo().solvertype = (FetiInfo::Solvertype)((yyvsp[(3) - (4)].ival));
	;}
    break;

  case 791:
#line 3229 "p.y"
    { domain->solInfo().sparse_maxsup  = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 792:
#line 3231 "p.y"
    { domain->solInfo().sparse_defblk  = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 793:
#line 3233 "p.y"
    { domain->solInfo().spooles_tau  = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 794:
#line 3235 "p.y"
    { domain->solInfo().spooles_maxsize = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 795:
#line 3237 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 0) {
            (yyvsp[(2) - (3)].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          domain->solInfo().spooles_maxdomainsize = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 796:
#line 3244 "p.y"
    { domain->solInfo().spooles_seed = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 797:
#line 3246 "p.y"
    { if(((yyvsp[(2) - (3)].fval) < 0.0) || ((yyvsp[(2) - (3)].fval) > 1.0)) {
            (yyvsp[(2) - (3)].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          domain->solInfo().spooles_maxzeros = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 798:
#line 3253 "p.y"
    { domain->solInfo().spooles_msglvl = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 799:
#line 3255 "p.y"
    { domain->solInfo().pivot = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 800:
#line 3257 "p.y"
    { domain->solInfo().spooles_scale = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 801:
#line 3259 "p.y"
    { domain->solInfo().spooles_renum = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 802:
#line 3261 "p.y"
    { domain->solInfo().mumps_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 803:
#line 3263 "p.y"
    { domain->solInfo().mumps_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 804:
#line 3265 "p.y"
    { domain->solInfo().goldfarb_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 805:
#line 3267 "p.y"
    { domain->solInfo().goldfarb_check = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 806:
#line 3269 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 807:
#line 3271 "p.y"
    { domain->solInfo().debug_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 808:
#line 3273 "p.y"
    { domain->solInfo().debug_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 809:
#line 3279 "p.y"
    { domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 810:
#line 3281 "p.y"
    { domain->solInfo().fetiInfo.precno = FetiInfo::lumped; ;}
    break;

  case 811:
#line 3283 "p.y"
    { if(((yyvsp[(3) - (4)].ival) < 0) || ((yyvsp[(3) - (4)].ival) > 3)) { 
            (yyvsp[(3) - (4)].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival);
	;}
    break;

  case 812:
#line 3290 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        ;}
    break;

  case 813:
#line 3297 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        ;}
    break;

  case 814:
#line 3304 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 815:
#line 3306 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.absolute_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 816:
#line 3309 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 817:
#line 3311 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.absolute_stagnation_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 818:
#line 3314 "p.y"
    { domain->solInfo().fetiInfo.primal_proj_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_proj_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 819:
#line 3317 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.dual_plan_maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 820:
#line 3320 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_plan_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 821:
#line 3323 "p.y"
    { domain->solInfo().fetiInfo.maxortho = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 822:
#line 3325 "p.y"
    { domain->solInfo().fetiInfo.noCoarse = 1; ;}
    break;

  case 823:
#line 3327 "p.y"
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

  case 824:
#line 3347 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 1; 
          domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 825:
#line 3350 "p.y"
    { domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 826:
#line 3352 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 827:
#line 3355 "p.y"
    { domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 828:
#line 3357 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 829:
#line 3360 "p.y"
    { domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 830:
#line 3362 "p.y"
    { domain->solInfo().fetiInfo.mpc_element = true; ;}
    break;

  case 831:
#line 3364 "p.y"
    { domain->solInfo().fetiInfo.fsi_element = true; ;}
    break;

  case 832:
#line 3366 "p.y"
    { domain->solInfo().fetiInfo.fsi_corner = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 833:
#line 3368 "p.y"
    { domain->solInfo().fetiInfo.splitLocalFsi = false; ;}
    break;

  case 834:
#line 3370 "p.y"
    { domain->solInfo().coupled_scale = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 835:
#line 3372 "p.y"
    { domain->solInfo().fetiInfo.wetcorners = true; ;}
    break;

  case 836:
#line 3374 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 837:
#line 3376 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.pick_unsafe_corners = bool((yyvsp[(3) - (4)].ival));
        ;}
    break;

  case 838:
#line 3380 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) {
            domain->solInfo().fetiInfo.corners = FetiInfo::noCorners;
            domain->solInfo().fetiInfo.pickAnyCorner = 0; 
            domain->solInfo().fetiInfo.bmpc = true;
            domain->solInfo().fetiInfo.pick_unsafe_corners = false;
            domain->solInfo().fetiInfo.augment = FetiInfo::none;
          }
        ;}
    break;

  case 839:
#line 3389 "p.y"
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

  case 840:
#line 3410 "p.y"
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

  case 841:
#line 3428 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (4)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 842:
#line 3433 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (5)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 843:
#line 3439 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(4) - (5)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 844:
#line 3445 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (6)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(5) - (6)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (6)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 845:
#line 3452 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 846:
#line 3454 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.orthotol2 = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 847:
#line 3457 "p.y"
    { domain->solInfo().fetiInfo.grbm_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 848:
#line 3459 "p.y"
    { domain->solInfo().fetiInfo.crbm_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 849:
#line 3461 "p.y"
    { domain->solInfo().fetiInfo.cct_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 850:
#line 3463 "p.y"
    { domain->solInfo().fetiInfo.rebuildcct = int((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 851:
#line 3465 "p.y"
    { domain->solInfo().fetiInfo.uproj = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 852:
#line 3467 "p.y"
    { domain->solInfo().fetiInfo.printMatLab = 1; ;}
    break;

  case 853:
#line 3469 "p.y"
    { domain->solInfo().printMatLab = 1;
          domain->solInfo().printMatLabFile = (yyvsp[(2) - (3)].strval); ;}
    break;

  case 854:
#line 3472 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 855:
#line 3474 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 856:
#line 3476 "p.y"
    {  domain->solInfo().fetiInfo.auxCoarseSolver = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 857:
#line 3478 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 858:
#line 3480 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 859:
#line 3484 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().localScaled = true; ;}
    break;

  case 860:
#line 3487 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 861:
#line 3491 "p.y"
    { domain->solInfo().fetiInfo.auxCoarseSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 862:
#line 3495 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().coarseScaled = true; ;}
    break;

  case 863:
#line 3498 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 864:
#line 3502 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival)!=0) fprintf(stderr," *** WARNING: Scaling not supported for this CCt solver \n");
          else domain->solInfo().fetiInfo.cctScaled = true; ;}
    break;

  case 865:
#line 3506 "p.y"
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

  case 866:
#line 3532 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(1) - (2)].ival); ;}
    break;

  case 867:
#line 3534 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = true; ;}
    break;

  case 868:
#line 3536 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 869:
#line 3538 "p.y"
    { domain->solInfo().fetiInfo.pickAnyCorner = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 870:
#line 3544 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().fetiInfo.nlPrecFlg = 1; 
          domain->solInfo().setKrylov(); 
        ;}
    break;

  case 871:
#line 3549 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
	  domain->solInfo().fetiInfo.nlPrecFlg = (yyvsp[(2) - (3)].ival);
	  domain->solInfo().setKrylov();
	;}
    break;

  case 872:
#line 3554 "p.y"
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

  case 873:
#line 3567 "p.y"
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

  case 874:
#line 3581 "p.y"
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

  case 875:
#line 3598 "p.y"
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

  case 876:
#line 3614 "p.y"
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

  case 877:
#line 3624 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 878:
#line 3626 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.numcgm2 = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 879:
#line 3629 "p.y"
    { domain->solInfo().fetiInfo.tolcgm = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 880:
#line 3631 "p.y"
    { domain->solInfo().fetiInfo.spaceDimension = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 881:
#line 3633 "p.y"
    { domain->solInfo().fetiInfo.krylovtype = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 882:
#line 3635 "p.y"
    { domain->solInfo().fetiInfo.krylovtype =  (yyvsp[(2) - (3)].ival); ;}
    break;

  case 883:
#line 3637 "p.y"
    { domain->solInfo().fetiInfo.lumpedinterface = 1; ;}
    break;

  case 884:
#line 3639 "p.y"
    { domain->solInfo().fetiInfo.saveMemCoarse = 1; ;}
    break;

  case 885:
#line 3641 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        ;}
    break;

  case 886:
#line 3646 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        ;}
    break;

  case 887:
#line 3652 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        ;}
    break;

  case 888:
#line 3659 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 889:
#line 3661 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.complex_hermitian = true; ;}
    break;

  case 890:
#line 3664 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 891:
#line 3666 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 892:
#line 3668 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 893:
#line 3670 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 894:
#line 3672 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 895:
#line 3675 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 896:
#line 3678 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 897:
#line 3681 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 898:
#line 3684 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); ;}
    break;

  case 899:
#line 3688 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); ;}
    break;

  case 900:
#line 3692 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 1) domain->solInfo().fetiInfo.useMRHS = false; ;}
    break;

  case 901:
#line 3694 "p.y"
    { domain->solInfo().fetiInfo.gamma = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 902:
#line 3696 "p.y"
    { domain->solInfo().fetiInfo.linesearch_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.linesearch_tau = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 903:
#line 3699 "p.y"
    { domain->solInfo().fetiInfo.bmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 904:
#line 3701 "p.y"
    { domain->solInfo().fetiInfo.dmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 905:
#line 3703 "p.y"
    { domain->solInfo().fetiInfo.cmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 906:
#line 3705 "p.y"
    { domain->solInfo().fetiInfo.c_normalize = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 907:
#line 3707 "p.y"
    { domain->solInfo().dbccheck = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 909:
#line 3714 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp; 
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        ;}
    break;

  case 910:
#line 3725 "p.y"
    { if(!(yyvsp[(2) - (3)].copt).lagrangeMult && (yyvsp[(2) - (3)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(2) - (3)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(2) - (3)].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[(2) - (3)].copt).constraint_hess; 
          domain->solInfo().constraint_hess_eps = (yyvsp[(2) - (3)].copt).constraint_hess_eps; ;}
    break;

  case 911:
#line 3731 "p.y"
    { if(!(yyvsp[(3) - (4)].copt).lagrangeMult && (yyvsp[(3) - (4)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(3) - (4)].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[(3) - (4)].copt).constraint_hess;
          domain->solInfo().constraint_hess_eps = (yyvsp[(3) - (4)].copt).constraint_hess_eps; ;}
    break;

  case 912:
#line 3739 "p.y"
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        ;}
    break;

  case 913:
#line 3746 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 914:
#line 3753 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (3)].fval);
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 915:
#line 3761 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (4)].fval);
          domain->solInfo().rhsZeroTol = (yyvsp[(4) - (4)].fval); ;}
    break;

  case 916:
#line 3770 "p.y"
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

  case 917:
#line 3780 "p.y"
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 918:
#line 3786 "p.y"
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[(2) - (2)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 919:
#line 3792 "p.y"
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[(3) - (3)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 920:
#line 3798 "p.y"
    { // Alternative input syntax for treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[(2) - (2)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 921:
#line 3804 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (3)].ival);
          (yyval.copt).constraint_hess_eps = 0; ;}
    break;

  case 922:
#line 3807 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (4)].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[(4) - (4)].fval); ;}
    break;

  case 923:
#line 3812 "p.y"
    { // hack??
	  domain->solInfo().acoustic = true;
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
        ;}
    break;

  case 924:
#line 3817 "p.y"
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

  case 925:
#line 3827 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        ;}
    break;

  case 926:
#line 3832 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        ;}
    break;

  case 927:
#line 3838 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        ;}
    break;

  case 928:
#line 3845 "p.y"
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        ;}
    break;

  case 929:
#line 3850 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        ;}
    break;

  case 930:
#line 3855 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        ;}
    break;

  case 936:
#line 3869 "p.y"
    { domain->setWaveDirections(0, (yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 937:
#line 3872 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 938:
#line 3875 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        ;}
    break;

  case 941:
#line 3885 "p.y"
    { domain->setFFPDirections((yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 943:
#line 3892 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        ;}
    break;

  case 945:
#line 3906 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        ;}
    break;

  case 946:
#line 3917 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        ;}
    break;

  case 947:
#line 3923 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::Static || domain->solInfo().probType == SolverInfo::None)
            domain->solInfo().probType = SolverInfo::NonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::Dynamic)
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          else if(domain->solInfo().probType == SolverInfo::TempDynamic) {
            domain->solInfo().order = 1;
            domain->solInfo().probType = SolverInfo::NonLinDynam;
          }
          domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().getNLInfo().setDefaults(); // just in case PIECEWISE is used under statics
        ;}
    break;

  case 948:
#line 3936 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength;
        ;}
    break;

  case 949:
#line 3941 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam;
        ;}
    break;

  case 950:
#line 3948 "p.y"
    { domain->solInfo().getNLInfo().linearelastic = true; ;}
    break;

  case 951:
#line 3950 "p.y"
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 952:
#line 3952 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 953:
#line 3954 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 954:
#line 3957 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (7)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[(5) - (7)].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[(6) - (7)].fval); ;}
    break;

  case 955:
#line 3962 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 956:
#line 3965 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (7)].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[(5) - (7)].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[(6) - (7)].ival); ;}
    break;

  case 957:
#line 3970 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 958:
#line 3973 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (5)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 959:
#line 3976 "p.y"
    { domain->solInfo().getNLInfo().unsymmetric = true; ;}
    break;

  case 960:
#line 3978 "p.y"
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 961:
#line 3980 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 962:
#line 3982 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (5)].ival); 
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 963:
#line 3985 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (6)].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (6)].ival);
          // note: currently we use either c1 or c2, but never both
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[(5) - (6)].fval);
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[(5) - (6)].fval); ;}
    break;

  case 964:
#line 3991 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (7)].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (7)].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[(5) - (7)].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[(5) - (7)].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[(6) - (7)].fval); ;}
    break;

  case 965:
#line 3997 "p.y"
    { domain->solInfo().getNLInfo().linesearch.type = (yyvsp[(3) - (8)].ival);
          domain->solInfo().getNLInfo().linesearch.maxit = (yyvsp[(4) - (8)].ival);
          domain->solInfo().getNLInfo().linesearch.c1 = (yyvsp[(5) - (8)].fval); 
          domain->solInfo().getNLInfo().linesearch.c2 = (yyvsp[(5) - (8)].fval);
          domain->solInfo().getNLInfo().linesearch.tau = (yyvsp[(6) - (8)].fval);
          domain->solInfo().getNLInfo().linesearch.verbose = bool((yyvsp[(7) - (8)].ival)); ;}
    break;

  case 966:
#line 4004 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true; ;}
    break;

  case 967:
#line 4006 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 968:
#line 4009 "p.y"
    { domain->solInfo().num_penalty_its = (yyvsp[(3) - (6)].ival); 
          domain->solInfo().penalty_tol = (yyvsp[(4) - (6)].fval);
          domain->solInfo().penalty_beta = (yyvsp[(5) - (6)].fval); ;}
    break;

  case 970:
#line 4016 "p.y"
    { 
          domain->solInfo().setNewton((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear; 
        ;}
    break;

  case 971:
#line 4021 "p.y"
    {
          domain->solInfo().setNewton((yyvsp[(2) - (4)].ival));
          domain->solInfo().getNLInfo().stepUpdateK = (yyvsp[(3) - (4)].ival);
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        ;}
    break;

  case 972:
#line 4058 "p.y"
    { domain->solInfo().setReOrtho(); ;}
    break;

  case 974:
#line 4063 "p.y"
    { geoSource->setControl((yyvsp[(3) - (10)].strval),(yyvsp[(7) - (10)].strval),(yyvsp[(9) - (10)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (10)].ival); ;}
    break;

  case 975:
#line 4065 "p.y"
    { geoSource->setControl((yyvsp[(3) - (12)].strval),(yyvsp[(7) - (12)].strval),(yyvsp[(9) - (12)].strval),(yyvsp[(11) - (12)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (12)].ival); ;}
    break;

  case 976:
#line 4075 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(2) - (3)].strval);
#endif
        ;}
    break;

  case 977:
#line 4081 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(3) - (4)].strval);
#endif
 ;}
    break;

  case 979:
#line 4090 "p.y"
    { domain->solInfo().contact_mode = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 980:
#line 4094 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (7)].ival)-1, (yyvsp[(3) - (7)].ival)-1, (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); ;}
    break;

  case 981:
#line 4097 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(8) - (9)].fval));;}
    break;

  case 982:
#line 4100 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), 0.0, (yyvsp[(8) - (9)].ival));;}
    break;

  case 983:
#line 4102 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (8)].ival)-1, (yyvsp[(3) - (8)].ival)-1, (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), 0.0, -1, (yyvsp[(7) - (8)].copt).lagrangeMult, (yyvsp[(7) - (8)].copt).penalty);;}
    break;

  case 984:
#line 4105 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (11)].ival)-1, (yyvsp[(3) - (11)].ival)-1, (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(10) - (11)].ival));;}
    break;

  case 985:
#line 4107 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (12)].ival)-1, (yyvsp[(3) - (12)].ival)-1, (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(10) - (12)].ival), (yyvsp[(11) - (12)].copt).lagrangeMult, (yyvsp[(11) - (12)].copt).penalty);;}
    break;

  case 987:
#line 4112 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 988:
#line 4117 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new BilinPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 989:
#line 4122 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 990:
#line 4127 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 991:
#line 4132 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 992:
#line 4137 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 993:
#line 4142 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1, 
             new ElaLinIsoMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
	 ;}
    break;

  case 994:
#line 4147 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new StVenantKirchhoffMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
         ;}
    break;

  case 995:
#line 4152 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new HenckyMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
         ;}
    break;

  case 996:
#line 4157 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ElaLinIsoMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)));
         ;}
    break;

  case 997:
#line 4162 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)));
         ;}
    break;

  case 998:
#line 4167 "p.y"
    {
            double params[3] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          ;}
    break;

  case 999:
#line 4173 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          ;}
    break;

  case 1000:
#line 4179 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          ;}
    break;

  case 1001:
#line 4185 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          ;}
    break;

  case 1002:
#line 4191 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          ;}
    break;

  case 1003:
#line 4197 "p.y"
    {
            double params[6] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          ;}
    break;

  case 1004:
#line 4203 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 1.0e-6, -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 1005:
#line 4209 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 1006:
#line 4215 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 1007:
#line 4221 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (24)].ival)-1,
             new ExpMat((yyvsp[(3) - (24)].ival), (yyvsp[(4) - (24)].fval), (yyvsp[(5) - (24)].fval), (yyvsp[(6) - (24)].fval), (yyvsp[(7) - (24)].fval), (yyvsp[(8) - (24)].fval), (yyvsp[(9) - (24)].fval), (yyvsp[(10) - (24)].fval), (yyvsp[(11) - (24)].fval), (yyvsp[(12) - (24)].fval), (yyvsp[(13) - (24)].fval), (yyvsp[(14) - (24)].fval), (yyvsp[(15) - (24)].fval), (yyvsp[(16) - (24)].fval), (yyvsp[(17) - (24)].fval), (yyvsp[(18) - (24)].fval), (yyvsp[(19) - (24)].fval), (yyvsp[(20) - (24)].fval), (yyvsp[(21) - (24)].fval), (yyvsp[(22) - (24)].fval), (yyvsp[(23) - (24)].fval)));
         ;}
    break;

  case 1008:
#line 4226 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new ExpMat((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1009:
#line 4231 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ExpMat((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1010:
#line 4236 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new ExpMat((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1011:
#line 4241 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new ExpMat((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1012:
#line 4246 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
             new ExpMat((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 1013:
#line 4251 "p.y"
    {
	   geoSource->loadMaterial((yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].strval));
	 ;}
    break;

  case 1014:
#line 4255 "p.y"
    {
	   geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].dlist));
	 ;}
    break;

  case 1016:
#line 4262 "p.y"
    { geoSource->setMatUsage((yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 1017:
#line 4264 "p.y"
    {
            for(int i = (yyvsp[(2) - (5)].ival)-1; i < (yyvsp[(3) - (5)].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[(4) - (5)].ival)-1);
	  ;}
    break;

  case 1018:
#line 4270 "p.y"
    { (yyval.dlist).nval = 0; ;}
    break;

  case 1019:
#line 4272 "p.y"
    { 
          if((yyvsp[(1) - (2)].dlist).nval == 32) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[(1) - (2)].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[(2) - (2)].fval);
 	;}
    break;

  case 1020:
#line 4282 "p.y"
    { (yyval.slist).nval = 0; ;}
    break;

  case 1021:
#line 4284 "p.y"
    { 
          if((yyvsp[(1) - (2)].slist).nval == 32) {
             fprintf(stderr, "Too many files!\n");
	     exit(-1);
          }
          (yyval.slist) = (yyvsp[(1) - (2)].slist);
          (yyval.slist).v[(yyval.slist).nval++] = (yyvsp[(2) - (2)].strval);
 	;}
    break;

  case 1022:
#line 4295 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (4)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 1023:
#line 4299 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (6)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (6)].ival)); ;}
    break;

  case 1024:
#line 4302 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (8)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (8)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(7) - (8)].ival)); ;}
    break;

  case 1025:
#line 4309 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().svdPodRom = true;;}
    break;

  case 1027:
#line 4317 "p.y"
    { for(int i=0; i<(yyvsp[(2) - (2)].slist).nval; ++i) domain->solInfo().snapfiPodRom.push_back(std::string((yyvsp[(2) - (2)].slist).v[i])); ;}
    break;

  case 1028:
#line 4328 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1029:
#line 4330 "p.y"
    { domain->solInfo().normalize = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1030:
#line 4332 "p.y"
    { for(int i=0; i<(yyvsp[(2) - (2)].dlist).nval; ++i) domain->solInfo().snapshotWeights.push_back((yyvsp[(2) - (2)].dlist).v[i]); ;}
    break;

  case 1031:
#line 4334 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1032:
#line 4336 "p.y"
    { for(int i=0; i<(yyvsp[(2) - (2)].slist).nval; ++i) domain->solInfo().robfi.push_back(std::string((yyvsp[(2) - (2)].slist).v[i])); ;}
    break;

  case 1033:
#line 4338 "p.y"
    { domain->solInfo().svdBlockSize = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1035:
#line 4344 "p.y"
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().DEIMBasisPod = true; ;}
    break;

  case 1037:
#line 4352 "p.y"
    { domain->solInfo().activatePodRom = true;
     domain->solInfo().setProbType(SolverInfo::PodRomOffline);
     domain->solInfo().UDEIMBasisPod = true; ;}
    break;

  case 1039:
#line 4360 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().samplingPodRom = true; ;}
    break;

  case 1042:
#line 4369 "p.y"
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().snapProjPodRom = true; ;}
    break;

  case 1044:
#line 4377 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 1045:
#line 4379 "p.y"
    { domain->solInfo().statePodRomFile.push_back((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 1046:
#line 4381 "p.y"
    { domain->solInfo().statePodRomFile.push_back((yyvsp[(2) - (3)].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[(3) - (3)].strval)); ;}
    break;

  case 1047:
#line 4384 "p.y"
    { domain->solInfo().statePodRomFile.push_back((yyvsp[(2) - (4)].strval));
    domain->solInfo().velocPodRomFile.push_back((yyvsp[(3) - (4)].strval));
    domain->solInfo().accelPodRomFile.push_back((yyvsp[(4) - (4)].strval)); ;}
    break;

  case 1048:
#line 4388 "p.y"
    { domain->solInfo().tolPodRom = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 1049:
#line 4390 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1050:
#line 4392 "p.y"
    { domain->solInfo().skipOffSet = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1051:
#line 4394 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 1052:
#line 4396 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (3)].ival); 
    domain->solInfo().forcePodSize = (yyvsp[(3) - (3)].ival);;}
    break;

  case 1053:
#line 4399 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (4)].ival); 
    domain->solInfo().forcePodSize = (yyvsp[(3) - (4)].ival);
    domain->solInfo().maxDeimBasisSize = (yyvsp[(4) - (4)].ival); ;}
    break;

  case 1054:
#line 4403 "p.y"
    { domain->solInfo().oocPodRom = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1055:
#line 4405 "p.y"
    { domain->solInfo().useMassNormalizedBasis = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1056:
#line 4407 "p.y"
    { domain->solInfo().useMassOrthogonalProjection = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1057:
#line 4409 "p.y"
    { domain->solInfo().reduceFollower = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1058:
#line 4411 "p.y"
    { domain->solInfo().PODerrornorm.push_back((yyvsp[(2) - (2)].strval)); ;}
    break;

  case 1059:
#line 4413 "p.y"
    { domain->solInfo().PODerrornorm.push_back((yyvsp[(2) - (3)].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[(3) - (3)].strval)); ;}
    break;

  case 1060:
#line 4416 "p.y"
    { domain->solInfo().PODerrornorm.push_back((yyvsp[(2) - (4)].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[(3) - (4)].strval));
    domain->solInfo().PODerrornorm.push_back((yyvsp[(4) - (4)].strval)); ;}
    break;

  case 1061:
#line 4420 "p.y"
    { domain->solInfo().localTol = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1062:
#line 4422 "p.y"
    { domain->solInfo().forcePodRomFile = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 1063:
#line 4424 "p.y"
    { domain->solInfo().forcePodRomFile = (yyvsp[(2) - (3)].strval);
    domain->solInfo().forcePodSize = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 1064:
#line 4427 "p.y"
    { domain->solInfo().forcePodRomFile = (yyvsp[(2) - (4)].strval); 
    domain->solInfo().forcePodSize = (yyvsp[(3) - (4)].ival); 
    domain->solInfo().maxDeimBasisSize = (yyvsp[(4) - (4)].ival); ;}
    break;

  case 1065:
#line 4431 "p.y"
    { domain->solInfo().selectFullNode = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1066:
#line 4433 "p.y"
    { domain->solInfo().selectFullElem = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1067:
#line 4435 "p.y"
    { domain->solInfo().computeForceSnap = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1068:
#line 4437 "p.y"
    { domain->solInfo().orthogForceSnap = bool((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 1070:
#line 4443 "p.y"
    { domain->solInfo().conwepConfigurations.push_back((yyvsp[(2) - (3)].blastData)); ;}
    break;

  case 1071:
#line 4448 "p.y"
    { domain->solInfo().activatePodRom = true;
    domain->solInfo().setProbType(SolverInfo::PodRomOffline);
    domain->solInfo().ROMPostProcess = true; ;}
    break;

  case 1073:
#line 4456 "p.y"
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[(2) - (2)].strval)); 
    domain->solInfo().numRODFile += 1; ;}
    break;

  case 1074:
#line 4459 "p.y"
    { domain->solInfo().RODConversionFiles.push_back((yyvsp[(2) - (3)].strval));
    domain->solInfo().numRODFile += 1; 
    domain->solInfo().skipPodRom = (yyvsp[(3) - (3)].ival);;}
    break;

  case 1075:
#line 4466 "p.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 1076:
#line 4468 "p.y"
    { (yyval.ival) = std::numeric_limits<int>::max(); ;}
    break;

  case 1077:
#line 4473 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 1078:
#line 4475 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].fval); ;}
    break;

  case 1079:
#line 4477 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::infinity(); ;}
    break;

  case 1080:
#line 4479 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 11590 "/lustre/home/tac688/newCodes/Empirical_GNAT/Parser.d/parser.cpp"
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
      /* If just tried and failed to reuse look-ahead token after an
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

  /* Else will try to reuse look-ahead token after shifting the error
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

  if (yyn == YYFINAL)
    YYACCEPT;

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

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
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


#line 4481 "p.y"


