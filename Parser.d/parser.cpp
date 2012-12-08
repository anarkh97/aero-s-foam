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
     LOCAL = 429,
     MASS = 430,
     MATERIALS = 431,
     MATLAB = 432,
     MAXITR = 433,
     MAXORTHO = 434,
     MAXVEC = 435,
     MODAL = 436,
     MPCPRECNO = 437,
     MPCPRECNOID = 438,
     MPCTYPE = 439,
     MPCTYPEID = 440,
     MPCSCALING = 441,
     MPCELEMENT = 442,
     MPCBLOCKID = 443,
     MPCBLK_OVERLAP = 444,
     MFTT = 445,
     MPTT = 446,
     MRHS = 447,
     MPCCHECK = 448,
     MUMPSICNTL = 449,
     MUMPSCNTL = 450,
     MECH = 451,
     MODEFILTER = 452,
     MOMENT = 453,
     NDTYPE = 454,
     NEIGPA = 455,
     NEWMARK = 456,
     NewLine = 457,
     NL = 458,
     NLMAT = 459,
     NLPREC = 460,
     NOCOARSE = 461,
     NODETOKEN = 462,
     NONINPC = 463,
     NSBSPV = 464,
     NLTOL = 465,
     NUMCGM = 466,
     NOSECONDARY = 467,
     NFRAMES = 468,
     OPTIMIZATION = 469,
     OUTPUT = 470,
     OUTPUT6 = 471,
     QSTATIC = 472,
     QLOAD = 473,
     PITA = 474,
     PITADISP6 = 475,
     PITAVEL6 = 476,
     NOFORCE = 477,
     MDPITA = 478,
     GLOBALBASES = 479,
     LOCALBASES = 480,
     TIMEREVERSIBLE = 481,
     REMOTECOARSE = 482,
     ORTHOPROJTOL = 483,
     READINITSEED = 484,
     JUMPCVG = 485,
     JUMPOUTPUT = 486,
     PRECNO = 487,
     PRECONDITIONER = 488,
     PRELOAD = 489,
     PRESSURE = 490,
     PRINTMATLAB = 491,
     PROJ = 492,
     PIVOT = 493,
     PRECTYPE = 494,
     PRECTYPEID = 495,
     PICKANYCORNER = 496,
     PADEPIVOT = 497,
     PROPORTIONING = 498,
     PLOAD = 499,
     PADEPOLES = 500,
     POINTSOURCE = 501,
     PLANEWAVE = 502,
     PTOL = 503,
     PLANTOL = 504,
     PMAXIT = 505,
     RADIATION = 506,
     RBMFILTER = 507,
     RBMSET = 508,
     READMODE = 509,
     REBUILD = 510,
     RENUM = 511,
     RENUMBERID = 512,
     REORTHO = 513,
     RESTART = 514,
     RECONS = 515,
     RECONSALG = 516,
     REBUILDCCT = 517,
     RANDOM = 518,
     RPROP = 519,
     RNORM = 520,
     REVERSENORMALS = 521,
     RIGID = 522,
     SCALING = 523,
     SCALINGTYPE = 524,
     SENSORS = 525,
     SOLVERTYPE = 526,
     SHIFT = 527,
     SPOOLESTAU = 528,
     SPOOLESSEED = 529,
     SPOOLESMAXSIZE = 530,
     SPOOLESMAXDOMAINSIZE = 531,
     SPOOLESMAXZEROS = 532,
     SPOOLESMSGLVL = 533,
     SPOOLESSCALE = 534,
     SPOOLESPIVOT = 535,
     SPOOLESRENUM = 536,
     SPARSEMAXSUP = 537,
     SPARSEDEFBLK = 538,
     STATS = 539,
     STRESSID = 540,
     SUBSPACE = 541,
     SURFACE = 542,
     SAVEMEMCOARSE = 543,
     SPACEDIMENSION = 544,
     SCATTERER = 545,
     STAGTOL = 546,
     SCALED = 547,
     SWITCH = 548,
     STABLE = 549,
     SUBTYPE = 550,
     STEP = 551,
     SOWER = 552,
     SHELLTHICKNESS = 553,
     SURF = 554,
     SPRINGMAT = 555,
     TANGENT = 556,
     TEMP = 557,
     TIME = 558,
     TOLEIG = 559,
     TOLFETI = 560,
     TOLJAC = 561,
     TOLPCG = 562,
     TOPFILE = 563,
     TOPOLOGY = 564,
     TRBM = 565,
     THERMOE = 566,
     THERMOH = 567,
     TETT = 568,
     TOLCGM = 569,
     TURKEL = 570,
     TIEDSURFACES = 571,
     THETA = 572,
     HRC = 573,
     THIRDNODE = 574,
     THERMMAT = 575,
     TDENFORC = 576,
     TESTULRICH = 577,
     THRU = 578,
     TOPFLAG = 579,
     USE = 580,
     USERDEFINEDISP = 581,
     USERDEFINEFORCE = 582,
     UPROJ = 583,
     UNSYMMETRIC = 584,
     USING = 585,
     VERSION = 586,
     WAVENUMBER = 587,
     WETCORNERS = 588,
     XPOST = 589,
     YMTT = 590,
     ZERO = 591,
     BINARY = 592,
     GEOMETRY = 593,
     DECOMPOSITION = 594,
     GLOBAL = 595,
     MATCHER = 596,
     CPUMAP = 597,
     NODALCONTACT = 598,
     MODE = 599,
     FRIC = 600,
     GAP = 601,
     OUTERLOOP = 602,
     EDGEWS = 603,
     WAVETYPE = 604,
     ORTHOTOL = 605,
     IMPE = 606,
     FREQ = 607,
     DPH = 608,
     WAVEMETHOD = 609,
     MATSPEC = 610,
     MATUSAGE = 611,
     BILINEARPLASTIC = 612,
     FINITESTRAINPLASTIC = 613,
     LINEARELASTIC = 614,
     STVENANTKIRCHHOFF = 615,
     LINPLSTRESS = 616,
     READ = 617,
     OPTCTV = 618,
     ISOTROPICLINEARELASTIC = 619,
     NEOHOOKEAN = 620,
     ISOTROPICLINEARELASTICJ2PLASTIC = 621,
     ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS = 622,
     HYPERELASTIC = 623,
     MOONEYRIVLIN = 624,
     HENCKY = 625,
     LOGSTRAINPLASTIC = 626,
     SVKPLSTRESS = 627,
     SURFACETOPOLOGY = 628,
     MORTARTIED = 629,
     MORTARSCALING = 630,
     MORTARINTEGRATIONRULE = 631,
     SEARCHTOL = 632,
     STDMORTAR = 633,
     DUALMORTAR = 634,
     WETINTERFACE = 635,
     NSUBS = 636,
     EXITAFTERDEC = 637,
     SKIP = 638,
     OUTPUTMEMORY = 639,
     OUTPUTWEIGHT = 640,
     WEIGHTLIST = 641,
     GMRESRESIDUAL = 642,
     SLOSH = 643,
     SLGRAV = 644,
     SLZEM = 645,
     SLZEMFILTER = 646,
     PDIR = 647,
     HEFSB = 648,
     HEFRS = 649,
     HEINTERFACE = 650,
     SNAPFI = 651,
     PODROB = 652,
     TRNVCT = 653,
     OFFSET = 654,
     ORTHOG = 655,
     SVDTOKEN = 656,
     SAMPLING = 657,
     PODSIZEMAX = 658,
     REFSUBSTRACT = 659,
     TOLER = 660
   };
#endif
/* Tokens.  */
#define ACTUATORS 258
#define AERO 259
#define AEROH 260
#define AEROTYPE 261
#define ANALYSIS 262
#define ARCLENGTH 263
#define ATTRIBUTES 264
#define AUGMENT 265
#define AUGMENTTYPE 266
#define AVERAGED 267
#define ATDARB 268
#define ACOU 269
#define ATDDNB 270
#define ATDROB 271
#define ARPACK 272
#define ATDDIR 273
#define ATDNEU 274
#define AXIHDIR 275
#define AXIHNEU 276
#define AXINUMMODES 277
#define AXINUMSLICES 278
#define AXIHSOMMER 279
#define AXIMPC 280
#define AUXCOARSESOLVER 281
#define ACMECNTL 282
#define ADDEDMASS 283
#define AEROEMBED 284
#define BLOCKDIAG 285
#define BOFFSET 286
#define BUCKLE 287
#define BGTL 288
#define BMPC 289
#define BINARYINPUT 290
#define BINARYOUTPUT 291
#define CHECKTOKEN 292
#define COARSESOLVER 293
#define COEF 294
#define CFRAMES 295
#define COLLOCATEDTYPE 296
#define CONVECTION 297
#define COMPOSITE 298
#define CONDITION 299
#define CONTROL 300
#define CORNER 301
#define CORNERTYPE 302
#define CURVE 303
#define CCTTOL 304
#define CCTSOLVER 305
#define CRHS 306
#define COUPLEDSCALE 307
#define CONTACTSURFACES 308
#define CMPC 309
#define CNORM 310
#define COMPLEXOUTTYPE 311
#define CONSTRMAT 312
#define CASES 313
#define CONSTRAINEDSURFACES 314
#define CSFRAMES 315
#define CSTYPE 316
#define DAMPING 317
#define DblConstant 318
#define DEM 319
#define DIMASS 320
#define DISP 321
#define DIRECT 322
#define DLAMBDA 323
#define DP 324
#define DYNAM 325
#define DETER 326
#define DECOMPOSE 327
#define DECOMPFILE 328
#define DMPC 329
#define DEBUGCNTL 330
#define DEBUGICNTL 331
#define CONSTRAINTS 332
#define MULTIPLIERS 333
#define PENALTY 334
#define EIGEN 335
#define EFRAMES 336
#define ELSCATTERER 337
#define END 338
#define ELHSOMMERFELD 339
#define EXPLICIT 340
#define EPSILON 341
#define ELEMENTARYFUNCTIONTYPE 342
#define FABMAT 343
#define FACOUSTICS 344
#define FETI 345
#define FETI2TYPE 346
#define FETIPREC 347
#define FFP 348
#define FFPDIR 349
#define FITALG 350
#define FLUMAT 351
#define FNAME 352
#define FLUX 353
#define FORCE 354
#define FRONTAL 355
#define FETIH 356
#define FILTEREIG 357
#define FREQSWEEP 358
#define FREQSWEEP1 359
#define FREQSWEEP2 360
#define FSINTERFACE 361
#define FSISCALING 362
#define FSIELEMENT 363
#define NOLOCALFSISPLITING 364
#define FSICORNER 365
#define FFIDEBUG 366
#define FAILSAFE 367
#define GEPS 368
#define GLOBALTOL 369
#define GRAVITY 370
#define GRBM 371
#define GTGSOLVER 372
#define GLOBALCRBMTOL 373
#define GROUP 374
#define GROUPTYPE 375
#define GOLDFARBTOL 376
#define GOLDFARBCHECK 377
#define HDIRICHLET 378
#define HEAT 379
#define HFETI 380
#define HNEUMAN 381
#define HSOMMERFELD 382
#define HFTT 383
#define HELMHOLTZ 384
#define HNBO 385
#define HELMMF 386
#define HELMSO 387
#define HSCBO 388
#define HWIBO 389
#define HZEM 390
#define HZEMFILTER 391
#define HLMPC 392
#define HELMSWEEP 393
#define HELMSWEEP1 394
#define HELMSWEEP2 395
#define HERMITIAN 396
#define HESSIAN 397
#define IACC 398
#define IDENTITY 399
#define IDIS 400
#define IDIS6 401
#define IntConstant 402
#define INTERFACELUMPED 403
#define ITEMP 404
#define ITERTYPE 405
#define IVEL 406
#define INCIDENCE 407
#define IHDIRICHLET 408
#define IHDSWEEP 409
#define IHNEUMANN 410
#define ISOLVERTYPE 411
#define INPC 412
#define INFINTY 413
#define JACOBI 414
#define KRYLOVTYPE 415
#define KIRLOC 416
#define LAYC 417
#define LAYN 418
#define LAYD 419
#define LAYO 420
#define LAYMAT 421
#define LFACTOR 422
#define LMPC 423
#define LOAD 424
#define LOBPCG 425
#define LOCALSOLVER 426
#define LINESEARCH 427
#define LUMPED 428
#define LOCAL 429
#define MASS 430
#define MATERIALS 431
#define MATLAB 432
#define MAXITR 433
#define MAXORTHO 434
#define MAXVEC 435
#define MODAL 436
#define MPCPRECNO 437
#define MPCPRECNOID 438
#define MPCTYPE 439
#define MPCTYPEID 440
#define MPCSCALING 441
#define MPCELEMENT 442
#define MPCBLOCKID 443
#define MPCBLK_OVERLAP 444
#define MFTT 445
#define MPTT 446
#define MRHS 447
#define MPCCHECK 448
#define MUMPSICNTL 449
#define MUMPSCNTL 450
#define MECH 451
#define MODEFILTER 452
#define MOMENT 453
#define NDTYPE 454
#define NEIGPA 455
#define NEWMARK 456
#define NewLine 457
#define NL 458
#define NLMAT 459
#define NLPREC 460
#define NOCOARSE 461
#define NODETOKEN 462
#define NONINPC 463
#define NSBSPV 464
#define NLTOL 465
#define NUMCGM 466
#define NOSECONDARY 467
#define NFRAMES 468
#define OPTIMIZATION 469
#define OUTPUT 470
#define OUTPUT6 471
#define QSTATIC 472
#define QLOAD 473
#define PITA 474
#define PITADISP6 475
#define PITAVEL6 476
#define NOFORCE 477
#define MDPITA 478
#define GLOBALBASES 479
#define LOCALBASES 480
#define TIMEREVERSIBLE 481
#define REMOTECOARSE 482
#define ORTHOPROJTOL 483
#define READINITSEED 484
#define JUMPCVG 485
#define JUMPOUTPUT 486
#define PRECNO 487
#define PRECONDITIONER 488
#define PRELOAD 489
#define PRESSURE 490
#define PRINTMATLAB 491
#define PROJ 492
#define PIVOT 493
#define PRECTYPE 494
#define PRECTYPEID 495
#define PICKANYCORNER 496
#define PADEPIVOT 497
#define PROPORTIONING 498
#define PLOAD 499
#define PADEPOLES 500
#define POINTSOURCE 501
#define PLANEWAVE 502
#define PTOL 503
#define PLANTOL 504
#define PMAXIT 505
#define RADIATION 506
#define RBMFILTER 507
#define RBMSET 508
#define READMODE 509
#define REBUILD 510
#define RENUM 511
#define RENUMBERID 512
#define REORTHO 513
#define RESTART 514
#define RECONS 515
#define RECONSALG 516
#define REBUILDCCT 517
#define RANDOM 518
#define RPROP 519
#define RNORM 520
#define REVERSENORMALS 521
#define RIGID 522
#define SCALING 523
#define SCALINGTYPE 524
#define SENSORS 525
#define SOLVERTYPE 526
#define SHIFT 527
#define SPOOLESTAU 528
#define SPOOLESSEED 529
#define SPOOLESMAXSIZE 530
#define SPOOLESMAXDOMAINSIZE 531
#define SPOOLESMAXZEROS 532
#define SPOOLESMSGLVL 533
#define SPOOLESSCALE 534
#define SPOOLESPIVOT 535
#define SPOOLESRENUM 536
#define SPARSEMAXSUP 537
#define SPARSEDEFBLK 538
#define STATS 539
#define STRESSID 540
#define SUBSPACE 541
#define SURFACE 542
#define SAVEMEMCOARSE 543
#define SPACEDIMENSION 544
#define SCATTERER 545
#define STAGTOL 546
#define SCALED 547
#define SWITCH 548
#define STABLE 549
#define SUBTYPE 550
#define STEP 551
#define SOWER 552
#define SHELLTHICKNESS 553
#define SURF 554
#define SPRINGMAT 555
#define TANGENT 556
#define TEMP 557
#define TIME 558
#define TOLEIG 559
#define TOLFETI 560
#define TOLJAC 561
#define TOLPCG 562
#define TOPFILE 563
#define TOPOLOGY 564
#define TRBM 565
#define THERMOE 566
#define THERMOH 567
#define TETT 568
#define TOLCGM 569
#define TURKEL 570
#define TIEDSURFACES 571
#define THETA 572
#define HRC 573
#define THIRDNODE 574
#define THERMMAT 575
#define TDENFORC 576
#define TESTULRICH 577
#define THRU 578
#define TOPFLAG 579
#define USE 580
#define USERDEFINEDISP 581
#define USERDEFINEFORCE 582
#define UPROJ 583
#define UNSYMMETRIC 584
#define USING 585
#define VERSION 586
#define WAVENUMBER 587
#define WETCORNERS 588
#define XPOST 589
#define YMTT 590
#define ZERO 591
#define BINARY 592
#define GEOMETRY 593
#define DECOMPOSITION 594
#define GLOBAL 595
#define MATCHER 596
#define CPUMAP 597
#define NODALCONTACT 598
#define MODE 599
#define FRIC 600
#define GAP 601
#define OUTERLOOP 602
#define EDGEWS 603
#define WAVETYPE 604
#define ORTHOTOL 605
#define IMPE 606
#define FREQ 607
#define DPH 608
#define WAVEMETHOD 609
#define MATSPEC 610
#define MATUSAGE 611
#define BILINEARPLASTIC 612
#define FINITESTRAINPLASTIC 613
#define LINEARELASTIC 614
#define STVENANTKIRCHHOFF 615
#define LINPLSTRESS 616
#define READ 617
#define OPTCTV 618
#define ISOTROPICLINEARELASTIC 619
#define NEOHOOKEAN 620
#define ISOTROPICLINEARELASTICJ2PLASTIC 621
#define ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS 622
#define HYPERELASTIC 623
#define MOONEYRIVLIN 624
#define HENCKY 625
#define LOGSTRAINPLASTIC 626
#define SVKPLSTRESS 627
#define SURFACETOPOLOGY 628
#define MORTARTIED 629
#define MORTARSCALING 630
#define MORTARINTEGRATIONRULE 631
#define SEARCHTOL 632
#define STDMORTAR 633
#define DUALMORTAR 634
#define WETINTERFACE 635
#define NSUBS 636
#define EXITAFTERDEC 637
#define SKIP 638
#define OUTPUTMEMORY 639
#define OUTPUTWEIGHT 640
#define WEIGHTLIST 641
#define GMRESRESIDUAL 642
#define SLOSH 643
#define SLGRAV 644
#define SLZEM 645
#define SLZEMFILTER 646
#define PDIR 647
#define HEFSB 648
#define HEFRS 649
#define HEINTERFACE 650
#define SNAPFI 651
#define PODROB 652
#define TRNVCT 653
#define OFFSET 654
#define ORTHOG 655
#define SVDTOKEN 656
#define SAMPLING 657
#define PODSIZEMAX 658
#define REFSUBSTRACT 659
#define TOLER 660




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
#line 23 "p.y"
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
/* Line 193 of yacc.c.  */
#line 960 "/lustre/home/mpotts/FEM/Parser.d/parser.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 973 "/lustre/home/mpotts/FEM/Parser.d/parser.cpp"

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
#define YYFINAL  646
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   5194

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  406
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  221
/* YYNRULES -- Number of rules.  */
#define YYNRULES  948
/* YYNRULES -- Number of states.  */
#define YYNSTATES  2321

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   660

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
     405
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
     261,   263,   265,   267,   269,   271,   277,   282,   285,   291,
     298,   305,   308,   315,   321,   329,   337,   346,   350,   353,
     356,   359,   362,   366,   369,   374,   378,   383,   388,   395,
     398,   402,   406,   409,   414,   419,   424,   429,   434,   438,
     441,   446,   451,   455,   459,   463,   467,   471,   474,   479,
     482,   487,   490,   495,   498,   503,   506,   509,   512,   515,
     518,   521,   524,   527,   532,   537,   543,   547,   550,   554,
     557,   561,   564,   568,   571,   583,   597,   612,   618,   621,
     624,   632,   642,   654,   667,   670,   676,   683,   686,   692,
     695,   701,   706,   711,   715,   719,   723,   727,   731,   735,
     738,   741,   745,   749,   753,   757,   763,   768,   774,   781,
     789,   793,   799,   804,   810,   817,   825,   829,   832,   836,
     839,   842,   846,   849,   853,   856,   859,   862,   866,   868,
     871,   875,   879,   882,   888,   892,   896,   900,   903,   907,
     910,   914,   919,   924,   930,   934,   938,   941,   946,   949,
     953,   956,   960,   963,   969,   972,   975,   978,   981,   984,
     987,   991,   996,  1001,  1010,  1015,  1020,  1024,  1028,  1034,
    1037,  1041,  1045,  1049,  1052,  1057,  1059,  1063,  1065,  1068,
    1071,  1074,  1078,  1084,  1091,  1098,  1102,  1107,  1114,  1120,
    1124,  1129,  1133,  1135,  1138,  1145,  1148,  1151,  1154,  1157,
    1160,  1165,  1171,  1176,  1179,  1183,  1186,  1190,  1193,  1198,
    1200,  1203,  1206,  1209,  1215,  1221,  1228,  1236,  1240,  1242,
    1244,  1247,  1249,  1251,  1253,  1256,  1258,  1260,  1263,  1265,
    1270,  1275,  1279,  1286,  1292,  1299,  1303,  1306,  1312,  1319,
    1322,  1325,  1333,  1343,  1347,  1350,  1357,  1361,  1363,  1366,
    1370,  1373,  1377,  1379,  1382,  1387,  1390,  1395,  1401,  1404,
    1408,  1413,  1419,  1422,  1426,  1433,  1436,  1440,  1447,  1451,
    1453,  1456,  1461,  1467,  1474,  1478,  1481,  1486,  1488,  1491,
    1495,  1497,  1500,  1505,  1511,  1518,  1522,  1525,  1530,  1534,
    1537,  1542,  1546,  1549,  1554,  1558,  1561,  1566,  1570,  1572,
    1575,  1579,  1583,  1587,  1590,  1594,  1597,  1603,  1606,  1610,
    1615,  1620,  1630,  1633,  1636,  1646,  1649,  1653,  1658,  1662,
    1666,  1674,  1677,  1681,  1686,  1689,  1692,  1698,  1706,  1717,
    1721,  1726,  1729,  1733,  1737,  1742,  1746,  1749,  1759,  1766,
    1769,  1772,  1776,  1786,  1790,  1800,  1803,  1807,  1812,  1816,
    1819,  1823,  1826,  1834,  1844,  1848,  1850,  1853,  1860,  1868,
    1877,  1887,  1889,  1892,  1894,  1897,  1900,  1904,  1911,  1916,
    1924,  1927,  1931,  1938,  1943,  1951,  1954,  1958,  1961,  1964,
    1968,  1971,  1975,  1981,  1986,  1991,  1994,  1998,  2001,  2004,
    2008,  2014,  2018,  2021,  2027,  2032,  2036,  2041,  2043,  2046,
    2050,  2053,  2070,  2090,  2110,  2131,  2155,  2179,  2189,  2204,
    2221,  2234,  2248,  2263,  2269,  2276,  2293,  2306,  2312,  2320,
    2329,  2339,  2351,  2362,  2374,  2387,  2395,  2404,  2414,  2419,
    2423,  2426,  2430,  2435,  2441,  2448,  2454,  2459,  2465,  2471,
    2478,  2486,  2494,  2503,  2511,  2520,  2525,  2529,  2532,  2538,
    2545,  2553,  2562,  2573,  2585,  2592,  2602,  2605,  2611,  2619,
    2623,  2626,  2631,  2634,  2640,  2648,  2651,  2657,  2664,  2672,
    2681,  2692,  2704,  2718,  2733,  2740,  2750,  2754,  2758,  2762,
    2766,  2770,  2773,  2779,  2784,  2788,  2796,  2803,  2808,  2810,
    2813,  2818,  2822,  2826,  2830,  2836,  2841,  2844,  2847,  2850,
    2853,  2856,  2859,  2871,  2876,  2879,  2887,  2890,  2895,  2902,
    2909,  2918,  2926,  2932,  2938,  2946,  2955,  2958,  2963,  2969,
    2973,  2976,  2980,  2983,  2988,  2995,  3002,  3011,  3013,  3015,
    3020,  3022,  3025,  3030,  3035,  3040,  3045,  3050,  3055,  3060,
    3066,  3071,  3077,  3083,  3089,  3096,  3104,  3113,  3123,  3130,
    3136,  3142,  3148,  3155,  3163,  3169,  3174,  3177,  3181,  3185,
    3189,  3194,  3199,  3203,  3207,  3211,  3215,  3219,  3223,  3227,
    3231,  3235,  3239,  3243,  3248,  3253,  3257,  3261,  3266,  3271,
    3276,  3281,  3286,  3291,  3295,  3299,  3303,  3308,  3312,  3317,
    3322,  3327,  3332,  3337,  3340,  3344,  3348,  3352,  3356,  3360,
    3364,  3368,  3371,  3374,  3378,  3381,  3385,  3388,  3392,  3397,
    3401,  3405,  3410,  3415,  3421,  3427,  3434,  3438,  3443,  3447,
    3451,  3455,  3459,  3463,  3466,  3470,  3474,  3478,  3482,  3487,
    3492,  3497,  3502,  3507,  3512,  3517,  3521,  3524,  3527,  3531,
    3535,  3538,  3542,  3547,  3553,  3560,  3568,  3571,  3575,  3580,
    3584,  3588,  3592,  3596,  3599,  3602,  3606,  3611,  3617,  3621,
    3626,  3630,  3634,  3638,  3642,  3648,  3654,  3659,  3664,  3671,
    3678,  3682,  3686,  3691,  3695,  3699,  3703,  3707,  3711,  3715,
    3718,  3722,  3727,  3729,  3732,  3736,  3741,  3747,  3749,  3752,
    3756,  3760,  3765,  3768,  3771,  3775,  3781,  3785,  3790,  3796,
    3801,  3806,  3811,  3819,  3827,  3836,  3840,  3850,  3860,  3871,
    3874,  3877,  3880,  3884,  3890,  3895,  3897,  3900,  3905,  3912,
    3918,  3920,  3923,  3928,  3932,  3935,  3939,  3942,  3945,  3948,
    3952,  3956,  3961,  3966,  3972,  3980,  3986,  3994,  3999,  4005,
    4009,  4014,  4018,  4023,  4028,  4035,  4038,  4042,  4047,  4053,
    4062,  4065,  4068,  4079,  4092,  4096,  4101,  4104,  4109,  4117,
    4127,  4137,  4146,  4158,  4171,  4174,  4184,  4195,  4205,  4216,
    4226,  4237,  4245,  4253,  4261,  4270,  4279,  4287,  4295,  4304,
    4313,  4323,  4334,  4345,  4357,  4370,  4395,  4403,  4412,  4422,
    4433,  4445,  4451,  4457,  4460,  4465,  4471,  4472,  4475,  4480,
    4487,  4496,  4499,  4503,  4506,  4510,  4513,  4516,  4520,  4523,
    4526,  4529,  4532,  4535,  4538,  4540,  4542,  4544,  4546
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     407,     0,    -1,   408,    83,    -1,   409,    -1,   408,   409,
      -1,   574,    -1,   488,    -1,   543,    -1,   551,    -1,   555,
      -1,   563,    -1,   584,    -1,   583,    -1,   582,    -1,   587,
      -1,   561,    -1,   591,    -1,   588,    -1,   589,    -1,   590,
      -1,   620,    -1,   537,    -1,   536,    -1,   539,    -1,   540,
      -1,   538,    -1,   541,    -1,   542,    -1,   444,    -1,   445,
      -1,   447,    -1,   446,    -1,   451,    -1,   452,    -1,   569,
      -1,   460,    -1,   448,    -1,   450,    -1,   611,    -1,   439,
      -1,   427,    -1,   428,    -1,   437,    -1,   424,    -1,   425,
      -1,   426,    -1,   547,    -1,   549,    -1,   471,    -1,   472,
      -1,   440,    -1,   474,    -1,   476,    -1,   477,    -1,   473,
      -1,   441,    -1,   442,    -1,   443,    -1,   614,    -1,   615,
      -1,   421,    -1,   613,    -1,   463,    -1,   465,    -1,   466,
      -1,   467,    -1,   469,    -1,   470,    -1,   454,    -1,   453,
      -1,   496,    -1,   497,    -1,   498,    -1,   490,    -1,   493,
      -1,   499,    -1,   455,    -1,   595,    -1,   610,    -1,   599,
      -1,   606,    -1,   608,    -1,   506,    -1,   524,    -1,   525,
      -1,   603,    -1,   500,    -1,   503,    -1,   509,    -1,   511,
      -1,   513,    -1,   515,    -1,   559,    -1,   487,    -1,   483,
      -1,   484,    -1,   485,    -1,   526,    -1,   528,    -1,   530,
      -1,   531,    -1,   532,    -1,   534,    -1,   419,    -1,   420,
      -1,   519,    -1,   520,    -1,   521,    -1,   522,    -1,   523,
      -1,   422,    -1,   423,    -1,   616,    -1,   468,    -1,   410,
      -1,   411,    -1,   412,    -1,   413,    -1,   414,    -1,   617,
      -1,   618,    -1,   564,    -1,   565,    -1,   567,    -1,   566,
      -1,   568,    -1,   571,    -1,   572,    -1,   489,    -1,   586,
      -1,   479,    -1,   573,    -1,   597,    -1,   621,    -1,   623,
      -1,   208,   202,   625,   625,   202,    -1,   157,   202,   625,
     202,    -1,   119,   202,    -1,   412,   120,   625,   625,   202,
      -1,   412,   120,   625,   625,   625,   202,    -1,   412,   120,
     299,   625,   625,   202,    -1,   263,   202,    -1,   413,   625,
     264,   626,   626,   202,    -1,   351,   202,   352,   626,   202,
      -1,   351,   202,   104,   626,   626,   625,   202,    -1,   351,
     202,   105,   626,   626,   625,   202,    -1,   351,   202,   103,
     626,   626,   625,   625,   202,    -1,   351,   202,   417,    -1,
     414,   418,    -1,   414,   482,    -1,   414,   415,    -1,   414,
     416,    -1,   242,   626,   202,    -1,   245,   202,    -1,   245,
     626,   626,   202,    -1,   103,   626,   202,    -1,   417,   626,
     625,   202,    -1,   260,   261,   625,   202,    -1,   260,   261,
     625,   625,   625,   202,    -1,   297,   202,    -1,    35,   293,
     202,    -1,    36,   293,   202,    -1,   337,   202,    -1,   420,
     338,    97,   202,    -1,   420,   339,    97,   202,    -1,   420,
     340,    97,   202,    -1,   420,   341,    97,   202,    -1,   420,
     342,    97,   202,    -1,     7,   625,   202,    -1,    72,   202,
      -1,   422,    73,    97,   202,    -1,   422,   381,   625,   202,
      -1,   422,   385,   202,    -1,   422,   384,   202,    -1,   422,
     382,   202,    -1,   422,   383,   202,    -1,   422,    71,   202,
      -1,   386,   202,    -1,   423,   625,   626,   202,    -1,   190,
     202,    -1,   424,   626,   626,   202,    -1,   191,   202,    -1,
     425,   626,   626,   202,    -1,   128,   202,    -1,   426,   626,
     626,   202,    -1,    43,   202,    -1,   427,   429,    -1,   427,
     431,    -1,   427,   432,    -1,   427,   433,    -1,   427,   434,
      -1,    40,   202,    -1,   428,   585,    -1,    39,   625,   202,
     430,    -1,   625,   625,   626,   202,    -1,   430,   625,   625,
     626,   202,    -1,   162,   625,   202,    -1,   431,   435,    -1,
     163,   625,   202,    -1,   432,   435,    -1,   164,   625,   202,
      -1,   433,   436,    -1,   165,   625,   202,    -1,   434,   436,
      -1,   625,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   202,    -1,   625,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   202,    -1,   625,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     202,    -1,   625,   625,   626,   626,   202,    -1,   166,   202,
      -1,   437,   438,    -1,   625,   626,   626,   626,   626,   626,
     202,    -1,   625,   626,   626,   626,   626,   626,   626,   626,
     202,    -1,   625,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   202,    -1,   625,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   202,    -1,    65,   202,    -1,
     439,   625,   625,   626,   202,    -1,   439,   625,   625,   625,
     626,   202,    -1,   115,   202,    -1,   440,   626,   626,   626,
     202,    -1,   259,   202,    -1,   441,    97,    97,    97,   202,
      -1,   441,    97,    97,   202,    -1,   441,    97,   625,   202,
      -1,   169,    97,   202,    -1,   325,    97,   202,    -1,   270,
     202,   544,    -1,     3,   202,   544,    -1,   327,   202,   544,
      -1,   326,   202,   544,    -1,   215,   202,    -1,   216,   202,
      -1,   215,   625,   202,    -1,   216,   625,   202,    -1,   448,
     449,   202,    -1,   285,    97,   625,    -1,   285,   625,   625,
      97,   625,    -1,   285,    97,   625,   625,    -1,   285,    97,
     625,   120,   625,    -1,   285,   625,   625,    97,   625,   625,
      -1,   285,   625,   625,    97,   625,   120,   625,    -1,   321,
      97,   625,    -1,   321,   625,   625,    97,   625,    -1,   321,
      97,   625,   625,    -1,   321,    97,   625,   120,   625,    -1,
     321,   625,   625,    97,   625,   625,    -1,   321,   625,   625,
      97,   625,   120,   625,    -1,   449,   207,   625,    -1,   449,
     287,    -1,   449,   626,   626,    -1,   449,    12,    -1,   449,
      56,    -1,   449,    56,   625,    -1,   449,   199,    -1,   449,
     199,   625,    -1,   449,   340,    -1,   449,   174,    -1,   449,
     177,    -1,   334,   324,   202,    -1,   456,    -1,    80,   202,
      -1,    80,   625,   202,    -1,   200,   625,   202,    -1,   286,
     202,    -1,   286,   625,   626,   626,   202,    -1,   209,   625,
     202,    -1,   304,   626,   202,    -1,   306,   626,   202,    -1,
      85,   202,    -1,   272,   626,   202,    -1,    17,   202,    -1,
      17,    97,   202,    -1,    17,    97,   625,   202,    -1,    17,
     626,   626,   202,    -1,    17,   626,   626,   625,   202,    -1,
     102,   293,   202,    -1,   159,   625,   202,    -1,   170,   202,
      -1,   451,   178,   625,   202,    -1,   322,   202,    -1,    28,
     625,   202,    -1,   388,   202,    -1,   389,   626,   202,    -1,
     175,   202,    -1,    44,   202,   626,   625,   202,    -1,    44,
     202,    -1,   308,   202,    -1,    70,   202,    -1,   456,   457,
      -1,   456,   478,    -1,   456,   482,    -1,   456,   181,   202,
      -1,   456,   294,   625,   202,    -1,   456,   294,   293,   202,
      -1,   456,   294,   625,   626,   626,   625,   625,   202,    -1,
     456,   143,   293,   202,    -1,   456,   336,   293,   202,    -1,
     456,   212,   202,    -1,   456,    37,   202,    -1,   456,    37,
     626,   626,   202,    -1,   201,   202,    -1,   196,   458,   202,
      -1,    14,   458,   202,    -1,   124,   459,   202,    -1,   626,
     626,    -1,   626,   626,   626,   626,    -1,   626,    -1,   626,
     626,   626,    -1,   626,    -1,   217,   202,    -1,   460,   461,
      -1,   460,   462,    -1,   460,   181,   202,    -1,   196,   626,
     626,   625,   202,    -1,   196,   626,   626,   625,   626,   202,
      -1,   124,   626,   626,   626,   625,   202,    -1,     4,     6,
     202,    -1,     4,   202,     6,   202,    -1,     4,   202,     6,
     626,   626,   202,    -1,     4,   202,     6,   626,   202,    -1,
     463,    41,   202,    -1,   463,   341,    97,   202,    -1,   463,
     464,   202,    -1,    29,    -1,   464,   625,    -1,     5,   202,
       6,   626,   626,   202,    -1,   312,   202,    -1,   311,   202,
      -1,   344,   202,    -1,   135,   202,    -1,   390,   202,    -1,
     310,   202,   626,   202,    -1,   116,   202,   626,   626,   202,
      -1,   116,   202,   626,   202,    -1,   116,   202,    -1,   197,
     625,   202,    -1,   197,   202,    -1,   252,   625,   202,    -1,
     252,   202,    -1,   252,   202,   475,   202,    -1,   625,    -1,
     475,   625,    -1,   136,   202,    -1,   391,   202,    -1,   303,
     626,   626,   626,   202,    -1,   219,   202,   625,   625,   480,
      -1,   219,   202,   625,   625,   625,   480,    -1,   223,   202,
     625,   625,   625,   625,   480,    -1,   202,   481,   480,    -1,
     202,    -1,   222,    -1,   224,   625,    -1,   225,    -1,   226,
      -1,   227,    -1,   228,   626,    -1,   229,    -1,   230,    -1,
     230,   626,    -1,   231,    -1,    62,   626,   626,   202,    -1,
      62,   181,   202,   545,    -1,   123,   202,   560,    -1,   153,
     202,   626,   626,   626,   202,    -1,   153,   202,   626,   626,
     202,    -1,   154,   202,   625,   626,   626,   202,    -1,   154,
     625,   202,    -1,   485,   486,    -1,   625,   626,   626,   626,
     202,    -1,   155,   202,   626,   626,   626,   202,    -1,    66,
     202,    -1,   488,   578,    -1,   488,   625,   323,   625,   625,
     626,   202,    -1,   488,   625,   323,   625,   296,   625,   625,
     626,   202,    -1,   488,   299,   578,    -1,    59,   202,    -1,
     489,   625,    61,   625,   625,   202,    -1,   392,   202,   491,
      -1,   492,    -1,   491,   492,    -1,   625,   626,   202,    -1,
     625,   202,    -1,   394,   202,   494,    -1,   495,    -1,   494,
     495,    -1,   625,   625,   505,   202,    -1,   302,   202,    -1,
     496,   625,   626,   202,    -1,   496,   299,   625,   626,   202,
      -1,    98,   202,    -1,    98,   625,   202,    -1,   497,   625,
     626,   202,    -1,   497,   299,   625,   626,   202,    -1,    42,
     202,    -1,    42,   625,   202,    -1,   498,   625,   626,   626,
     626,   202,    -1,   251,   202,    -1,   251,   625,   202,    -1,
     499,   625,   626,   626,   626,   202,    -1,   127,   202,   501,
      -1,   502,    -1,   501,   502,    -1,   625,   625,   626,   202,
      -1,   625,   625,   625,   626,   202,    -1,   625,   625,   625,
     625,   626,   202,    -1,    84,   202,   504,    -1,   503,   504,
      -1,   625,   625,   505,   202,    -1,   625,    -1,   505,   625,
      -1,   290,   202,   507,    -1,   508,    -1,   507,   508,    -1,
     625,   625,   626,   202,    -1,   625,   625,   625,   626,   202,
      -1,   625,   625,   625,   625,   626,   202,    -1,    82,   202,
     510,    -1,   509,   510,    -1,   625,   625,   505,   202,    -1,
     130,   202,   512,    -1,   511,   512,    -1,   625,   625,   505,
     202,    -1,   134,   202,   514,    -1,   513,   514,    -1,   625,
     625,   505,   202,    -1,   133,   202,   516,    -1,   515,   516,
      -1,   625,   625,   505,   202,    -1,   625,   626,   202,    -1,
     517,    -1,   518,   517,    -1,    18,   202,   518,    -1,    19,
     202,   518,    -1,    13,   626,   202,    -1,   521,   504,    -1,
      15,   626,   202,    -1,   522,   512,    -1,    16,   626,   626,
     626,   202,    -1,   523,   516,    -1,    93,   625,   202,    -1,
      93,   625,   625,   202,    -1,    94,   625,   202,   604,    -1,
      20,   202,   626,   626,   202,   626,   626,   626,   202,    -1,
     526,   527,    -1,   625,   202,    -1,    21,   202,   626,   626,
     202,   626,   626,   626,   202,    -1,   528,   529,    -1,   625,
     625,   202,    -1,   625,   625,   625,   202,    -1,    22,   625,
     202,    -1,    23,   625,   202,    -1,    24,   202,   625,   202,
     626,   626,   202,    -1,   532,   533,    -1,   625,   625,   202,
      -1,   625,   625,   625,   202,    -1,    25,   202,    -1,   534,
     535,    -1,   625,   626,   626,   625,   202,    -1,   625,   626,
     626,   625,   626,   626,   202,    -1,   625,   626,   626,   625,
     626,   626,   626,   626,   626,   202,    -1,   254,    97,   202,
      -1,   254,    97,   625,   202,    -1,   145,   202,    -1,   145,
     336,   202,    -1,   145,   202,   544,    -1,   537,   181,   202,
     545,    -1,   146,   626,   202,    -1,   146,   202,    -1,   538,
     625,   626,   626,   626,   626,   626,   626,   202,    -1,   538,
     625,   626,   626,   626,   202,    -1,   113,   202,    -1,    32,
     202,    -1,   220,   625,   202,    -1,   539,   625,   626,   626,
     626,   626,   626,   626,   202,    -1,   221,   625,   202,    -1,
     540,   625,   626,   626,   626,   626,   626,   626,   202,    -1,
     151,   202,    -1,   151,   202,   544,    -1,   541,   181,   202,
     545,    -1,   149,   202,   546,    -1,    99,   202,    -1,    99,
     625,   202,    -1,   543,   578,    -1,   543,   625,   323,   625,
     625,   626,   202,    -1,   543,   625,   323,   625,   296,   625,
     625,   626,   202,    -1,   543,   299,   578,    -1,   578,    -1,
     544,   578,    -1,   625,   323,   625,   625,   626,   202,    -1,
     544,   625,   323,   625,   625,   626,   202,    -1,   625,   323,
     625,   296,   625,   625,   626,   202,    -1,   544,   625,   323,
     625,   296,   625,   625,   626,   202,    -1,   579,    -1,   545,
     579,    -1,   580,    -1,   546,   580,    -1,   335,   202,    -1,
     335,   202,   548,    -1,    48,   625,   202,   626,   626,   202,
      -1,   548,   626,   626,   202,    -1,   548,    48,   625,   202,
     626,   626,   202,    -1,   313,   202,    -1,   313,   202,   550,
      -1,    48,   625,   202,   626,   626,   202,    -1,   550,   626,
     626,   202,    -1,   550,    48,   625,   202,   626,   626,   202,
      -1,   168,   202,    -1,   168,   202,   552,    -1,   553,   554,
      -1,   552,   554,    -1,   552,   553,   554,    -1,   625,   202,
      -1,   625,   626,   202,    -1,   625,   626,   344,   625,   202,
      -1,   625,   626,   598,   202,    -1,   625,   625,   626,   202,
      -1,   137,   202,    -1,   137,   202,   556,    -1,   557,   558,
      -1,   556,   558,    -1,   556,   557,   558,    -1,   625,    51,
     626,   626,   202,    -1,   625,   626,   202,    -1,   625,   202,
      -1,   625,   625,   626,   626,   202,    -1,   625,   625,   626,
     202,    -1,   126,   202,   560,    -1,   126,   625,   202,   560,
      -1,   581,    -1,   560,   581,    -1,   176,   202,   562,    -1,
     561,   562,    -1,   625,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   202,    -1,
     625,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,    62,   626,   626,   202,    -1,
     625,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   267,   625,   626,   202,    -1,
     625,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   202,
      -1,   625,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
      62,   626,   626,   202,    -1,   625,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   267,   625,   626,   202,    -1,   625,
     626,   626,   626,   626,   626,   626,   626,   202,    -1,   625,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   202,    -1,   625,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,    62,   626,   202,
      -1,   625,    96,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   202,    -1,   625,    96,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   202,    -1,   625,    96,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   202,    -1,   625,    96,   626,   626,   202,    -1,   625,
      96,   626,   626,   626,   202,    -1,   625,    88,   625,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   625,   625,
     625,   202,    -1,   625,   320,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   202,    -1,   625,    57,   625,   626,
     202,    -1,   625,    57,   625,   626,   626,   626,   202,    -1,
     625,    57,   625,   626,   626,   626,   626,   202,    -1,   625,
      57,   625,   626,   626,   626,   626,   626,   202,    -1,   625,
      57,   625,   626,   626,   626,   626,   626,   626,   625,   202,
      -1,   625,    57,   625,   626,    87,   626,   626,   626,   626,
     202,    -1,   625,    57,   625,   626,    87,   626,   626,   626,
     626,   626,   202,    -1,   625,    57,   625,   626,    87,   626,
     626,   626,   626,   626,   626,   202,    -1,   625,    57,   625,
     626,   300,   626,   202,    -1,   625,    57,   625,   626,   300,
     626,   626,   202,    -1,   625,    57,   625,   626,   300,   626,
     626,   626,   202,    -1,   625,   300,   626,   202,    -1,   309,
     202,   576,    -1,   563,   576,    -1,   373,   625,   202,    -1,
     373,   625,   266,   202,    -1,   373,   625,   298,   626,   202,
      -1,   373,   625,   298,   626,   266,   202,    -1,   564,   625,
     625,   577,   202,    -1,   374,   625,   625,   202,    -1,   374,
     625,   625,   378,   202,    -1,   374,   625,   625,   379,   202,
      -1,   374,   625,   625,   377,   626,   202,    -1,   374,   625,
     625,   377,   626,   626,   202,    -1,   374,   625,   625,   378,
     377,   626,   202,    -1,   374,   625,   625,   378,   377,   626,
     626,   202,    -1,   374,   625,   625,   379,   377,   626,   202,
      -1,   374,   625,   625,   379,   377,   626,   626,   202,    -1,
     380,   625,   625,   202,    -1,   380,   625,   202,    -1,   316,
     202,    -1,   567,   625,   625,   625,   202,    -1,   567,   625,
     625,   625,   625,   202,    -1,   567,   625,   625,   625,   625,
     626,   202,    -1,   567,   625,   625,   625,   625,   626,   626,
     202,    -1,   567,   625,   625,   625,   625,   626,   626,   625,
     626,   202,    -1,   567,   625,   625,   625,   625,   626,   626,
     625,   626,   626,   202,    -1,   567,   625,   625,   625,   598,
     202,    -1,   567,   625,   625,   625,   625,   626,   626,   598,
     202,    -1,   106,   202,    -1,   568,   625,   625,   625,   202,
      -1,   568,   625,   625,   625,   626,   626,   202,    -1,   393,
     202,   570,    -1,   569,   570,    -1,   625,   625,   505,   202,
      -1,   395,   202,    -1,   571,   625,   625,   625,   202,    -1,
     571,   625,   625,   625,   626,   626,   202,    -1,    53,   202,
      -1,   572,   625,   625,   625,   202,    -1,   572,   625,   625,
     625,   625,   202,    -1,   572,   625,   625,   625,   625,   626,
     202,    -1,   572,   625,   625,   625,   625,   626,   626,   202,
      -1,   572,   625,   625,   625,   625,   626,   626,   625,   626,
     202,    -1,   572,   625,   625,   625,   625,   626,   626,   625,
     626,   626,   202,    -1,   572,   625,   625,   625,   625,   626,
     626,   625,   626,   626,   626,   626,   202,    -1,   572,   625,
     625,   625,   625,   626,   626,   625,   626,   626,   626,   626,
     626,   202,    -1,   572,   625,   625,   625,   598,   202,    -1,
     572,   625,   625,   625,   625,   626,   626,   598,   202,    -1,
      27,   625,   202,    -1,   111,   625,   202,    -1,   375,   626,
     202,    -1,   376,   625,   202,    -1,   207,   202,   575,    -1,
     574,   575,    -1,   625,   626,   626,   626,   202,    -1,   625,
     626,   626,   202,    -1,   625,   626,   202,    -1,   625,   626,
     626,   626,   625,   625,   202,    -1,   625,   626,   626,   626,
     625,   202,    -1,   625,   625,   577,   202,    -1,   625,    -1,
     577,   625,    -1,   625,   625,   626,   202,    -1,   625,   625,
     202,    -1,   625,   626,   202,    -1,   625,   626,   202,    -1,
     625,   625,   626,   626,   202,    -1,   625,   625,   626,   202,
      -1,    60,   202,    -1,   582,   585,    -1,   213,   202,    -1,
     583,   585,    -1,    81,   202,    -1,   584,   585,    -1,   625,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   202,
      -1,   625,   319,   625,   202,    -1,    31,   202,    -1,   586,
     625,   625,   626,   626,   626,   202,    -1,     9,   202,    -1,
     587,   625,   625,   202,    -1,   587,   625,   625,   318,   626,
     202,    -1,   587,   625,   625,   625,   625,   202,    -1,   587,
     625,   625,   625,   625,   318,   626,   202,    -1,   587,   625,
     625,   625,   317,   626,   202,    -1,   587,   625,   625,   144,
     202,    -1,   587,   625,   625,   625,   202,    -1,   587,   625,
     625,   625,   625,   625,   202,    -1,   587,   625,   625,   625,
     625,   317,   626,   202,    -1,   235,   202,    -1,   588,   625,
     626,   202,    -1,   588,   625,   625,   626,   202,    -1,   588,
     299,   517,    -1,   173,   202,    -1,   173,   625,   202,    -1,
     234,   202,    -1,   590,   625,   626,   202,    -1,   590,   625,
     323,   625,   626,   202,    -1,   590,   625,   626,   626,   626,
     202,    -1,   590,   625,   323,   625,   626,   626,   626,   202,
      -1,   594,    -1,   593,    -1,   591,    58,   592,   202,    -1,
     625,    -1,   592,   625,    -1,   284,   202,   150,   202,    -1,
     593,   232,   625,   202,    -1,   593,   178,   625,   202,    -1,
     593,   307,   626,   202,    -1,   593,   179,   625,   202,    -1,
     593,   295,   625,   202,    -1,   284,   202,    67,   202,    -1,
     284,   202,    67,   625,   202,    -1,   284,   202,   271,   202,
      -1,   284,   202,   271,   238,   202,    -1,   284,   202,   271,
     329,   202,    -1,   284,   202,   150,   625,   202,    -1,   284,
     202,   150,   625,   626,   202,    -1,   284,   202,   150,   625,
     626,   625,   202,    -1,   284,   202,   150,   625,   626,   625,
     625,   202,    -1,   284,   202,   150,   625,   626,   625,   625,
     625,   202,    -1,   284,   202,    90,   625,   626,   202,    -1,
     284,   202,    90,   625,   202,    -1,   284,   202,    90,    69,
     202,    -1,   284,   202,    90,   353,   202,    -1,   284,   202,
      90,   625,    91,   202,    -1,   284,   202,    90,   625,   626,
     625,   202,    -1,    90,   625,   626,   625,   202,    -1,   284,
     202,    90,   202,    -1,    90,   202,    -1,    90,   625,   202,
      -1,    90,    69,   202,    -1,    90,   353,   202,    -1,    90,
     625,    91,   202,    -1,    30,   202,   271,   202,    -1,   282,
     625,   202,    -1,   283,   625,   202,    -1,   273,   626,   202,
      -1,   275,   625,   202,    -1,   276,   625,   202,    -1,   274,
     625,   202,    -1,   277,   626,   202,    -1,   278,   625,   202,
      -1,   280,   293,   202,    -1,   279,   625,   202,    -1,   281,
     625,   202,    -1,   194,   625,   625,   202,    -1,   195,   625,
     626,   202,    -1,   121,   626,   202,    -1,   122,   293,   202,
      -1,   594,   178,   625,   202,    -1,    76,   625,   625,   202,
      -1,    75,   625,   626,   202,    -1,   594,   232,    92,   202,
      -1,   594,   232,   173,   202,    -1,   594,   232,   625,   202,
      -1,   239,   625,   202,    -1,   239,   240,   202,    -1,   305,
     626,   202,    -1,   305,   626,   626,   202,    -1,   291,   626,
     202,    -1,   291,   626,   626,   202,    -1,   248,   626,   626,
     202,    -1,   250,   625,   625,   202,    -1,   249,   626,   626,
     202,    -1,   594,   179,   625,   202,    -1,   206,   202,    -1,
     237,   625,   202,    -1,   268,   625,   202,    -1,   268,   269,
     202,    -1,   186,   625,   202,    -1,   186,   269,   202,    -1,
     107,   625,   202,    -1,   107,   269,   202,    -1,   187,   202,
      -1,   108,   202,    -1,   110,   625,   202,    -1,   109,   202,
      -1,    52,   626,   202,    -1,   333,   202,    -1,    46,    47,
     202,    -1,    46,    47,   625,   202,    -1,    46,    11,   202,
      -1,    10,    11,   202,    -1,    10,    11,   253,   202,    -1,
      10,   348,   625,   202,    -1,    10,   348,   349,   625,   202,
      -1,    10,   348,   625,   354,   202,    -1,    10,   348,   349,
     625,   354,   202,    -1,   350,   626,   202,    -1,   350,   626,
     626,   202,    -1,   114,   626,   202,    -1,   118,   626,   202,
      -1,    49,   626,   202,    -1,   262,   293,   202,    -1,   328,
     625,   202,    -1,   236,   202,    -1,   171,   271,   202,    -1,
      38,   271,   202,    -1,    26,   271,   202,    -1,    50,   271,
     202,    -1,   171,   271,   238,   202,    -1,   171,   271,   292,
     202,    -1,    38,   271,   238,   202,    -1,    26,   271,   238,
     202,    -1,    38,   271,   292,   202,    -1,    50,   271,   238,
     202,    -1,    50,   271,   292,   202,    -1,   331,   625,   202,
      -1,   117,   202,    -1,   387,   202,    -1,   387,   293,   202,
      -1,   241,   625,   202,    -1,   205,   202,    -1,   205,   625,
     202,    -1,   125,   625,   626,   202,    -1,   125,   625,   626,
     625,   202,    -1,   125,   625,   626,   625,   626,   202,    -1,
     125,   625,   626,   625,   626,   625,   202,    -1,   125,   202,
      -1,   211,   625,   202,    -1,   211,   625,   626,   202,    -1,
     314,   626,   202,    -1,   289,   625,   202,    -1,   160,   625,
     202,    -1,   160,   156,   202,    -1,   148,   202,    -1,   288,
     202,    -1,   315,   625,   202,    -1,   315,   625,   626,   202,
      -1,   315,   625,   626,   626,   202,    -1,   347,   150,   202,
      -1,   347,   150,   141,   202,    -1,   184,   625,   202,    -1,
     184,   185,   202,    -1,   182,   625,   202,    -1,   182,   183,
     202,    -1,   182,   625,   189,   625,   202,    -1,   182,   183,
     189,   625,   202,    -1,   182,   625,   625,   202,    -1,   182,
     183,   188,   202,    -1,   182,   625,   625,   189,   625,   202,
      -1,   182,   183,   188,   189,   625,   202,    -1,   192,   625,
     202,    -1,   243,   626,   202,    -1,   172,   625,   626,   202,
      -1,    34,   293,   202,    -1,    74,   293,   202,    -1,    54,
     293,   202,    -1,    55,   293,   202,    -1,   193,   625,   202,
      -1,    89,   202,   596,    -1,   626,   202,    -1,    77,   598,
     202,    -1,    77,   202,   598,   202,    -1,    67,    -1,    67,
     626,    -1,    67,   626,   626,    -1,    67,   626,   626,   626,
      -1,    67,   626,   626,   626,   626,    -1,    78,    -1,    79,
     626,    -1,    78,    79,   626,    -1,   598,   142,   625,    -1,
     598,   142,   625,   626,    -1,   129,   202,    -1,   101,   202,
      -1,   332,   626,   202,    -1,   332,   626,   626,   626,   202,
      -1,    33,   625,   202,    -1,    33,   625,   626,   202,    -1,
      33,   625,   626,   626,   202,    -1,   246,   625,   202,   601,
      -1,   247,   625,   202,   601,    -1,   152,   625,   202,   601,
      -1,   129,   202,   139,   626,   626,   625,   202,    -1,   129,
     202,   140,   626,   626,   625,   202,    -1,   129,   202,   138,
     626,   626,   625,   625,   202,    -1,   129,   202,   600,    -1,
     129,   202,   139,   626,   626,   625,   626,   626,   202,    -1,
     129,   202,   140,   626,   626,   625,   626,   626,   202,    -1,
     129,   202,   138,   626,   626,   625,   625,   626,   626,   202,
      -1,   599,   482,    -1,   599,   418,    -1,   599,   415,    -1,
     138,   626,   202,    -1,   138,   626,   626,   626,   202,    -1,
     600,   626,   625,   202,    -1,   602,    -1,   601,   602,    -1,
     626,   626,   626,   202,    -1,   161,   202,   626,   626,   626,
     202,    -1,   603,   626,   626,   626,   202,    -1,   605,    -1,
     604,   605,    -1,   626,   626,   626,   202,    -1,   131,   202,
     607,    -1,   626,   202,    -1,   132,   202,   609,    -1,   626,
     202,    -1,    64,   202,    -1,   203,   202,    -1,   611,     8,
     202,    -1,   611,   204,   202,    -1,   611,   178,   625,   202,
      -1,   611,   210,   626,   202,    -1,   611,   210,   626,   626,
     202,    -1,   611,   210,   626,   626,   626,   626,   202,    -1,
     611,    68,   626,   626,   202,    -1,   611,    68,   626,   626,
     625,   625,   202,    -1,   611,    95,   625,   202,    -1,   611,
      95,   625,   625,   202,    -1,   611,   329,   202,    -1,   611,
     167,   626,   202,    -1,   611,   112,   202,    -1,   611,   112,
     626,   202,    -1,   611,   198,   625,   202,    -1,   611,    79,
     625,   626,   626,   202,    -1,   611,   612,    -1,   255,   625,
     202,    -1,   255,   625,   625,   202,    -1,   255,   202,   301,
     625,   202,    -1,   255,   202,   301,   625,   202,   233,   625,
     202,    -1,   258,   202,    -1,    45,   202,    -1,    45,   202,
      97,   202,   625,   202,    97,   202,    97,   202,    -1,    45,
     202,    97,   202,   625,   202,    97,   202,    97,   202,    97,
     202,    -1,   214,    97,   202,    -1,   615,   202,    97,   202,
      -1,   343,   202,    -1,   343,   344,   625,   202,    -1,   616,
     625,   625,   626,   626,   626,   202,    -1,   616,   625,   625,
     626,   626,   626,   346,   626,   202,    -1,   616,   625,   625,
     626,   626,   626,   344,   625,   202,    -1,   616,   625,   625,
     626,   626,   626,   598,   202,    -1,   616,   625,   625,   626,
     626,   626,   346,   626,   344,   625,   202,    -1,   616,   625,
     625,   626,   626,   626,   346,   626,   344,   625,   598,   202,
      -1,   355,   202,    -1,   617,   625,   357,   626,   626,   626,
     626,   626,   202,    -1,   617,   625,   357,   626,   626,   626,
     626,   626,   626,   202,    -1,   617,   625,   358,   626,   626,
     626,   626,   626,   202,    -1,   617,   625,   358,   626,   626,
     626,   626,   626,   626,   202,    -1,   617,   625,   371,   626,
     626,   626,   626,   626,   202,    -1,   617,   625,   371,   626,
     626,   626,   626,   626,   626,   202,    -1,   617,   625,   359,
     626,   626,   626,   202,    -1,   617,   625,   360,   626,   626,
     626,   202,    -1,   617,   625,   370,   626,   626,   626,   202,
      -1,   617,   625,   361,   626,   626,   626,   626,   202,    -1,
     617,   625,   372,   626,   626,   626,   626,   202,    -1,   617,
     625,   364,   626,   626,   626,   202,    -1,   617,   625,   365,
     626,   626,   626,   202,    -1,   617,   625,   365,   626,   626,
     626,   626,   202,    -1,   617,   625,   369,   626,   626,   626,
     626,   202,    -1,   617,   625,   369,   626,   626,   626,   626,
     626,   202,    -1,   617,   625,   366,   626,   626,   626,   626,
     626,   626,   202,    -1,   617,   625,   367,   626,   626,   626,
     626,   626,   626,   202,    -1,   617,   625,   367,   626,   626,
     626,   626,   626,   626,   626,   202,    -1,   617,   625,   367,
     626,   626,   626,   626,   626,   626,   626,   626,   202,    -1,
     617,   625,   363,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   202,    -1,   617,   625,   363,   626,   626,
     626,   202,    -1,   617,   625,   363,   626,   626,   626,   626,
     202,    -1,   617,   625,   363,   626,   626,   626,   626,   626,
     202,    -1,   617,   625,   363,   626,   626,   626,   626,   626,
     626,   202,    -1,   617,   625,   363,   626,   626,   626,   626,
     626,   626,   626,   202,    -1,   617,   362,    97,    97,   202,
      -1,   617,   625,    97,   619,   202,    -1,   356,   202,    -1,
     618,   625,   625,   202,    -1,   618,   625,   625,   625,   202,
      -1,    -1,   619,   626,    -1,   256,   202,   257,   202,    -1,
     256,   202,   257,   202,   257,   202,    -1,   256,   202,   257,
     202,   257,   202,   257,   202,    -1,   401,   202,    -1,   621,
     622,   202,    -1,   396,    97,    -1,   396,    97,   625,    -1,
     403,   625,    -1,   402,   202,    -1,   623,   624,   202,    -1,
     397,    97,    -1,   398,    97,    -1,   405,   626,    -1,   383,
     625,    -1,   399,   625,    -1,   403,   625,    -1,   147,    -1,
     147,    -1,    63,    -1,   158,    -1,    86,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   149,   149,   155,   156,   159,   160,   162,   164,   165,
     166,   167,   168,   169,   170,   172,   173,   174,   175,   177,
     179,   180,   181,   182,   183,   184,   186,   187,   188,   189,
     190,   191,   192,   193,   194,   195,   196,   197,   198,   199,
     200,   201,   202,   203,   205,   207,   209,   210,   211,   212,
     213,   214,   215,   216,   217,   218,   219,   220,   221,   222,
     223,   224,   226,   227,   228,   229,   230,   231,   232,   233,
     234,   236,   238,   240,   242,   244,   246,   247,   248,   249,
     250,   251,   252,   253,   254,   255,   256,   257,   258,   259,
     260,   261,   262,   264,   265,   267,   268,   269,   270,   271,
     272,   273,   274,   275,   276,   277,   279,   281,   282,   283,
     284,   285,   286,   288,   289,   290,   291,   292,   294,   295,
     296,   297,   299,   301,   303,   305,   307,   309,   311,   313,
     314,   315,   316,   317,   318,   321,   328,   334,   335,   340,
     352,   358,   359,   363,   367,   369,   371,   373,   374,   375,
     376,   377,   380,   384,   386,   391,   394,   398,   437,   486,
     488,   490,   494,   495,   497,   499,   501,   503,   507,   515,
     517,   519,   521,   523,   525,   527,   529,   533,   535,   543,
     545,   549,   551,   555,   557,   561,   562,   563,   564,   565,
     566,   569,   570,   574,   578,   580,   584,   586,   590,   592,
     596,   598,   602,   604,   608,   614,   621,   630,   634,   635,
     639,   646,   652,   657,   663,   664,   666,   670,   671,   675,
     676,   680,   683,   688,   693,   697,   702,   708,   715,   722,
     724,   726,   728,   730,   734,   736,   738,   740,   743,   745,
     747,   749,   751,   753,   755,   757,   760,   762,   764,   766,
     768,   770,   772,   774,   776,   778,   780,   784,   788,   789,
     792,   795,   797,   799,   801,   803,   805,   807,   809,   811,
     813,   816,   820,   823,   826,   828,   830,   833,   835,   837,
     841,   843,   847,   851,   854,   858,   862,   864,   865,   866,
     867,   869,   871,   873,   880,   882,   884,   886,   888,   894,
     896,   897,   899,   902,   904,   906,   908,   914,   925,   928,
     929,   930,   934,   936,   940,   950,   953,   956,   964,   969,
     971,   973,   977,   979,   983,   987,   991,   995,   999,  1003,
    1007,  1015,  1018,  1021,  1026,  1028,  1032,  1034,  1036,  1039,
    1047,  1057,  1061,  1065,  1069,  1075,  1080,  1089,  1090,  1093,
    1095,  1097,  1099,  1101,  1103,  1105,  1107,  1109,  1111,  1115,
    1117,  1122,  1126,  1132,  1140,  1146,  1151,  1154,  1160,  1168,
    1170,  1172,  1174,  1176,  1189,  1191,  1200,  1204,  1206,  1210,
    1212,  1216,  1257,  1261,  1267,  1275,  1277,  1280,  1289,  1291,
    1293,  1296,  1306,  1308,  1310,  1315,  1317,  1319,  1324,  1327,
    1328,  1331,  1333,  1335,  1339,  1340,  1343,  1349,  1351,  1356,
    1359,  1360,  1363,  1366,  1369,  1374,  1375,  1378,  1383,  1384,
    1387,  1391,  1392,  1395,  1402,  1403,  1406,  1410,  1414,  1416,
    1420,  1424,  1428,  1430,  1433,  1435,  1438,  1442,  1445,  1447,
    1451,  1457,  1462,  1466,  1470,  1475,  1479,  1481,  1485,  1493,
    1500,  1505,  1508,  1510,  1514,  1516,  1520,  1522,  1524,  1528,
    1531,  1537,  1539,  1541,  1544,  1550,  1552,  1554,  1565,  1576,
    1579,  1584,  1586,  1599,  1601,  1613,  1615,  1618,  1624,  1629,
    1631,  1633,  1635,  1637,  1639,  1647,  1649,  1651,  1653,  1655,
    1657,  1661,  1663,  1667,  1669,  1673,  1674,  1677,  1679,  1681,
    1685,  1686,  1689,  1691,  1693,  1697,  1698,  1701,  1705,  1707,
    1713,  1716,  1719,  1723,  1730,  1742,  1743,  1746,  1748,  1750,
    1754,  1756,  1758,  1762,  1770,  1780,  1782,  1787,  1789,  1793,
    1794,  1797,  1805,  1814,  1825,  1834,  1844,  1856,  1863,  1871,
    1879,  1895,  1912,  1930,  1938,  1946,  1965,  1973,  1980,  1989,
    1999,  2010,  2023,  2036,  2050,  2065,  2074,  2084,  2095,  2103,
    2104,  2107,  2113,  2119,  2126,  2134,  2158,  2160,  2162,  2166,
    2168,  2170,  2172,  2174,  2178,  2184,  2186,  2194,  2196,  2202,
    2209,  2216,  2223,  2231,  2240,  2247,  2257,  2259,  2265,  2273,
    2276,  2279,  2285,  2287,  2293,  2303,  2305,  2312,  2319,  2326,
    2333,  2341,  2350,  2359,  2368,  2376,  2386,  2388,  2390,  2392,
    2395,  2397,  2401,  2403,  2405,  2407,  2409,  2413,  2418,  2420,
    2425,  2427,  2431,  2435,  2439,  2441,  2445,  2446,  2450,  2451,
    2455,  2456,  2460,  2465,  2469,  2470,  2477,  2479,  2481,  2485,
    2487,  2491,  2493,  2498,  2503,  2508,  2515,  2516,  2518,  2523,
    2529,  2534,  2541,  2543,  2545,  2550,  2553,  2561,  2562,  2563,
    2566,  2568,  2572,  2576,  2578,  2580,  2582,  2584,  2588,  2591,
    2594,  2597,  2602,  2606,  2609,  2612,  2615,  2618,  2626,  2632,
    2636,  2641,  2647,  2653,  2659,  2665,  2668,  2670,  2673,  2677,
    2682,  2687,  2693,  2695,  2697,  2699,  2701,  2708,  2710,  2717,
    2719,  2721,  2723,  2725,  2727,  2729,  2731,  2733,  2735,  2737,
    2743,  2745,  2747,  2754,  2761,  2768,  2770,  2773,  2775,  2778,
    2781,  2784,  2787,  2789,  2791,  2811,  2814,  2816,  2819,  2821,
    2824,  2826,  2828,  2830,  2832,  2834,  2836,  2838,  2840,  2844,
    2853,  2874,  2892,  2897,  2903,  2909,  2916,  2918,  2921,  2923,
    2925,  2927,  2929,  2931,  2933,  2935,  2937,  2939,  2941,  2945,
    2948,  2952,  2956,  2959,  2963,  2967,  2993,  2995,  2997,  2999,
    3005,  3010,  3015,  3028,  3042,  3059,  3075,  3085,  3087,  3090,
    3092,  3094,  3096,  3098,  3100,  3102,  3107,  3113,  3120,  3122,
    3125,  3127,  3129,  3131,  3133,  3136,  3139,  3142,  3145,  3149,
    3153,  3155,  3157,  3160,  3162,  3164,  3166,  3168,  3172,  3175,
    3186,  3192,  3200,  3207,  3214,  3222,  3231,  3241,  3247,  3253,
    3259,  3262,  3267,  3270,  3280,  3287,  3295,  3300,  3306,  3313,
    3318,  3323,  3328,  3336,  3344,  3352,  3353,  3362,  3371,  3380,
    3381,  3382,  3385,  3394,  3403,  3408,  3409,  3412,  3416,  3419,
    3424,  3425,  3428,  3432,  3435,  3446,  3449,  3460,  3466,  3479,
    3484,  3491,  3493,  3495,  3498,  3503,  3506,  3511,  3514,  3517,
    3519,  3525,  3527,  3530,  3532,  3536,  3539,  3544,  3554,  3561,
    3573,  3577,  3578,  3580,  3590,  3596,  3604,  3605,  3609,  3612,
    3615,  3617,  3620,  3622,  3626,  3627,  3632,  3637,  3642,  3647,
    3652,  3657,  3662,  3667,  3672,  3677,  3682,  3688,  3694,  3700,
    3706,  3712,  3718,  3724,  3730,  3736,  3741,  3746,  3751,  3756,
    3761,  3766,  3770,  3776,  3777,  3779,  3786,  3787,  3798,  3802,
    3805,  3812,  3816,  3820,  3822,  3829,  3834,  3838,  3842,  3844,
    3846,  3848,  3850,  3852,  3857,  3862,  3864,  3866,  3868
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
  "LOBPCG", "LOCALSOLVER", "LINESEARCH", "LUMPED", "LOCAL", "MASS",
  "MATERIALS", "MATLAB", "MAXITR", "MAXORTHO", "MAXVEC", "MODAL",
  "MPCPRECNO", "MPCPRECNOID", "MPCTYPE", "MPCTYPEID", "MPCSCALING",
  "MPCELEMENT", "MPCBLOCKID", "MPCBLK_OVERLAP", "MFTT", "MPTT", "MRHS",
  "MPCCHECK", "MUMPSICNTL", "MUMPSCNTL", "MECH", "MODEFILTER", "MOMENT",
  "NDTYPE", "NEIGPA", "NEWMARK", "NewLine", "NL", "NLMAT", "NLPREC",
  "NOCOARSE", "NODETOKEN", "NONINPC", "NSBSPV", "NLTOL", "NUMCGM",
  "NOSECONDARY", "NFRAMES", "OPTIMIZATION", "OUTPUT", "OUTPUT6", "QSTATIC",
  "QLOAD", "PITA", "PITADISP6", "PITAVEL6", "NOFORCE", "MDPITA",
  "GLOBALBASES", "LOCALBASES", "TIMEREVERSIBLE", "REMOTECOARSE",
  "ORTHOPROJTOL", "READINITSEED", "JUMPCVG", "JUMPOUTPUT", "PRECNO",
  "PRECONDITIONER", "PRELOAD", "PRESSURE", "PRINTMATLAB", "PROJ", "PIVOT",
  "PRECTYPE", "PRECTYPEID", "PICKANYCORNER", "PADEPIVOT", "PROPORTIONING",
  "PLOAD", "PADEPOLES", "POINTSOURCE", "PLANEWAVE", "PTOL", "PLANTOL",
  "PMAXIT", "RADIATION", "RBMFILTER", "RBMSET", "READMODE", "REBUILD",
  "RENUM", "RENUMBERID", "REORTHO", "RESTART", "RECONS", "RECONSALG",
  "REBUILDCCT", "RANDOM", "RPROP", "RNORM", "REVERSENORMALS", "RIGID",
  "SCALING", "SCALINGTYPE", "SENSORS", "SOLVERTYPE", "SHIFT", "SPOOLESTAU",
  "SPOOLESSEED", "SPOOLESMAXSIZE", "SPOOLESMAXDOMAINSIZE",
  "SPOOLESMAXZEROS", "SPOOLESMSGLVL", "SPOOLESSCALE", "SPOOLESPIVOT",
  "SPOOLESRENUM", "SPARSEMAXSUP", "SPARSEDEFBLK", "STATS", "STRESSID",
  "SUBSPACE", "SURFACE", "SAVEMEMCOARSE", "SPACEDIMENSION", "SCATTERER",
  "STAGTOL", "SCALED", "SWITCH", "STABLE", "SUBTYPE", "STEP", "SOWER",
  "SHELLTHICKNESS", "SURF", "SPRINGMAT", "TANGENT", "TEMP", "TIME",
  "TOLEIG", "TOLFETI", "TOLJAC", "TOLPCG", "TOPFILE", "TOPOLOGY", "TRBM",
  "THERMOE", "THERMOH", "TETT", "TOLCGM", "TURKEL", "TIEDSURFACES",
  "THETA", "HRC", "THIRDNODE", "THERMMAT", "TDENFORC", "TESTULRICH",
  "THRU", "TOPFLAG", "USE", "USERDEFINEDISP", "USERDEFINEFORCE", "UPROJ",
  "UNSYMMETRIC", "USING", "VERSION", "WAVENUMBER", "WETCORNERS", "XPOST",
  "YMTT", "ZERO", "BINARY", "GEOMETRY", "DECOMPOSITION", "GLOBAL",
  "MATCHER", "CPUMAP", "NODALCONTACT", "MODE", "FRIC", "GAP", "OUTERLOOP",
  "EDGEWS", "WAVETYPE", "ORTHOTOL", "IMPE", "FREQ", "DPH", "WAVEMETHOD",
  "MATSPEC", "MATUSAGE", "BILINEARPLASTIC", "FINITESTRAINPLASTIC",
  "LINEARELASTIC", "STVENANTKIRCHHOFF", "LINPLSTRESS", "READ", "OPTCTV",
  "ISOTROPICLINEARELASTIC", "NEOHOOKEAN",
  "ISOTROPICLINEARELASTICJ2PLASTIC",
  "ISOTROPICLINEARELASTICJ2PLASTICPLANESTRESS", "HYPERELASTIC",
  "MOONEYRIVLIN", "HENCKY", "LOGSTRAINPLASTIC", "SVKPLSTRESS",
  "SURFACETOPOLOGY", "MORTARTIED", "MORTARSCALING",
  "MORTARINTEGRATIONRULE", "SEARCHTOL", "STDMORTAR", "DUALMORTAR",
  "WETINTERFACE", "NSUBS", "EXITAFTERDEC", "SKIP", "OUTPUTMEMORY",
  "OUTPUTWEIGHT", "WEIGHTLIST", "GMRESRESIDUAL", "SLOSH", "SLGRAV",
  "SLZEM", "SLZEMFILTER", "PDIR", "HEFSB", "HEFRS", "HEINTERFACE",
  "SNAPFI", "PODROB", "TRNVCT", "OFFSET", "ORTHOG", "SVDTOKEN", "SAMPLING",
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
  "ComplexBC_Data", "ConstrainedSurfaceFrameDList", "NodeFrameDList",
  "FrameDList", "Frame", "BoffsetList", "Attributes", "Pressure", "Lumped",
  "Preload", "Statics", "CasesList", "IterSolver", "Solver", "OldHelmInfo",
  "FAcousticData", "Constraints", "ConstraintOptionsData", "HelmInfo",
  "HelmSweep", "IncidenceList", "IncidenceVector", "KirchhoffLocations",
  "FFPDirList", "FFPDirVector", "HelmMFInfo", "FAcousticDataMF",
  "HelmSOInfo", "FAcousticDataSO", "DEMInfo", "NLInfo", "NewtonInfo",
  "OrthoInfo", "Control", "Optimization", "NodalContact", "MatSpec",
  "MatUsage", "FloatList", "Renumbering", "SvdToken", "SvdOption",
  "Sampling", "SamplingOption", "Integer", "Float", 0
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
     655,   656,   657,   658,   659,   660
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
{
       0,   406,   407,   408,   408,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   409,   409,   409,   409,   409,
     409,   409,   409,   409,   409,   410,   411,   412,   412,   412,
     412,   413,   413,   414,   414,   414,   414,   414,   414,   414,
     414,   414,   415,   416,   416,   417,   417,   418,   418,   419,
     419,   419,   420,   420,   420,   420,   420,   420,   421,   422,
     422,   422,   422,   422,   422,   422,   422,   423,   423,   424,
     424,   425,   425,   426,   426,   427,   427,   427,   427,   427,
     427,   428,   428,   429,   430,   430,   431,   431,   432,   432,
     433,   433,   434,   434,   435,   435,   435,   436,   437,   437,
     438,   438,   438,   438,   439,   439,   439,   440,   440,   441,
     441,   441,   441,   442,   443,   444,   445,   446,   447,   448,
     448,   448,   448,   448,   449,   449,   449,   449,   449,   449,
     449,   449,   449,   449,   449,   449,   449,   449,   449,   449,
     449,   449,   449,   449,   449,   449,   449,   450,   451,   451,
     451,   451,   451,   451,   451,   451,   451,   451,   451,   451,
     451,   451,   451,   451,   451,   451,   451,   451,   451,   451,
     452,   452,   453,   454,   454,   455,   456,   456,   456,   456,
     456,   456,   456,   456,   456,   456,   456,   456,   456,   457,
     457,   457,   457,   458,   458,   458,   458,   459,   460,   460,
     460,   460,   461,   461,   462,   463,   463,   463,   463,   463,
     463,   463,   464,   464,   465,   466,   467,   468,   469,   470,
     471,   472,   472,   472,   473,   473,   474,   474,   474,   475,
     475,   476,   477,   478,   479,   479,   479,   480,   480,   481,
     481,   481,   481,   481,   481,   481,   481,   481,   481,   482,
     482,   483,   484,   484,   485,   485,   485,   486,   487,   488,
     488,   488,   488,   488,   489,   489,   490,   491,   491,   492,
     492,   493,   494,   494,   495,   496,   496,   496,   497,   497,
     497,   497,   498,   498,   498,   499,   499,   499,   500,   501,
     501,   502,   502,   502,   503,   503,   504,   505,   505,   506,
     507,   507,   508,   508,   508,   509,   509,   510,   511,   511,
     512,   513,   513,   514,   515,   515,   516,   517,   518,   518,
     519,   520,   521,   521,   522,   522,   523,   523,   524,   524,
     525,   526,   526,   527,   528,   528,   529,   529,   530,   531,
     532,   532,   533,   533,   534,   534,   535,   535,   535,   536,
     536,   537,   537,   537,   537,   538,   538,   538,   538,   538,
     538,   539,   539,   540,   540,   541,   541,   541,   542,   543,
     543,   543,   543,   543,   543,   544,   544,   544,   544,   544,
     544,   545,   545,   546,   546,   547,   547,   548,   548,   548,
     549,   549,   550,   550,   550,   551,   551,   552,   552,   552,
     553,   553,   553,   553,   554,   555,   555,   556,   556,   556,
     557,   557,   557,   558,   558,   559,   559,   560,   560,   561,
     561,   562,   562,   562,   562,   562,   562,   562,   562,   562,
     562,   562,   562,   562,   562,   562,   562,   562,   562,   562,
     562,   562,   562,   562,   562,   562,   562,   562,   562,   563,
     563,   564,   564,   564,   564,   564,   565,   565,   565,   565,
     565,   565,   565,   565,   565,   566,   566,   567,   567,   567,
     567,   567,   567,   567,   567,   567,   568,   568,   568,   569,
     569,   570,   571,   571,   571,   572,   572,   572,   572,   572,
     572,   572,   572,   572,   572,   572,   573,   573,   573,   573,
     574,   574,   575,   575,   575,   575,   575,   576,   577,   577,
     578,   578,   579,   580,   581,   581,   582,   582,   583,   583,
     584,   584,   585,   585,   586,   586,   587,   587,   587,   587,
     587,   587,   587,   587,   587,   587,   588,   588,   588,   588,
     589,   589,   590,   590,   590,   590,   590,   591,   591,   591,
     592,   592,   593,   593,   593,   593,   593,   593,   594,   594,
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
     594,   594,   594,   594,   594,   594,   594,   594,   595,   596,
     597,   597,   598,   598,   598,   598,   598,   598,   598,   598,
     598,   598,   599,   599,   599,   599,   599,   599,   599,   599,
     599,   599,   599,   599,   599,   599,   599,   599,   599,   599,
     599,   599,   600,   600,   600,   601,   601,   602,   603,   603,
     604,   604,   605,   606,   607,   608,   609,   610,   611,   611,
     611,   611,   611,   611,   611,   611,   611,   611,   611,   611,
     611,   611,   611,   611,   611,   611,   612,   612,   612,   612,
     613,   614,   614,   614,   615,   615,   616,   616,   616,   616,
     616,   616,   616,   616,   617,   617,   617,   617,   617,   617,
     617,   617,   617,   617,   617,   617,   617,   617,   617,   617,
     617,   617,   617,   617,   617,   617,   617,   617,   617,   617,
     617,   617,   617,   618,   618,   618,   619,   619,   620,   620,
     620,   621,   621,   622,   622,   622,   623,   623,   624,   624,
     624,   624,   624,   624,   625,   626,   626,   626,   626
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
       1,     1,     1,     1,     1,     5,     4,     2,     5,     6,
       6,     2,     6,     5,     7,     7,     8,     3,     2,     2,
       2,     2,     3,     2,     4,     3,     4,     4,     6,     2,
       3,     3,     2,     4,     4,     4,     4,     4,     3,     2,
       4,     4,     3,     3,     3,     3,     3,     2,     4,     2,
       4,     2,     4,     2,     4,     2,     2,     2,     2,     2,
       2,     2,     2,     4,     4,     5,     3,     2,     3,     2,
       3,     2,     3,     2,    11,    13,    14,     5,     2,     2,
       7,     9,    11,    12,     2,     5,     6,     2,     5,     2,
       5,     4,     4,     3,     3,     3,     3,     3,     3,     2,
       2,     3,     3,     3,     3,     5,     4,     5,     6,     7,
       3,     5,     4,     5,     6,     7,     3,     2,     3,     2,
       2,     3,     2,     3,     2,     2,     2,     3,     1,     2,
       3,     3,     2,     5,     3,     3,     3,     2,     3,     2,
       3,     4,     4,     5,     3,     3,     2,     4,     2,     3,
       2,     3,     2,     5,     2,     2,     2,     2,     2,     2,
       3,     4,     4,     8,     4,     4,     3,     3,     5,     2,
       3,     3,     3,     2,     4,     1,     3,     1,     2,     2,
       2,     3,     5,     6,     6,     3,     4,     6,     5,     3,
       4,     3,     1,     2,     6,     2,     2,     2,     2,     2,
       4,     5,     4,     2,     3,     2,     3,     2,     4,     1,
       2,     2,     2,     5,     5,     6,     7,     3,     1,     1,
       2,     1,     1,     1,     2,     1,     1,     2,     1,     4,
       4,     3,     6,     5,     6,     3,     2,     5,     6,     2,
       2,     7,     9,     3,     2,     6,     3,     1,     2,     3,
       2,     3,     1,     2,     4,     2,     4,     5,     2,     3,
       4,     5,     2,     3,     6,     2,     3,     6,     3,     1,
       2,     4,     5,     6,     3,     2,     4,     1,     2,     3,
       1,     2,     4,     5,     6,     3,     2,     4,     3,     2,
       4,     3,     2,     4,     3,     2,     4,     3,     1,     2,
       3,     3,     3,     2,     3,     2,     5,     2,     3,     4,
       4,     9,     2,     2,     9,     2,     3,     4,     3,     3,
       7,     2,     3,     4,     2,     2,     5,     7,    10,     3,
       4,     2,     3,     3,     4,     3,     2,     9,     6,     2,
       2,     3,     9,     3,     9,     2,     3,     4,     3,     2,
       3,     2,     7,     9,     3,     1,     2,     6,     7,     8,
       9,     1,     2,     1,     2,     2,     3,     6,     4,     7,
       2,     3,     6,     4,     7,     2,     3,     2,     2,     3,
       2,     3,     5,     4,     4,     2,     3,     2,     2,     3,
       5,     3,     2,     5,     4,     3,     4,     1,     2,     3,
       2,    16,    19,    19,    20,    23,    23,     9,    14,    16,
      12,    13,    14,     5,     6,    16,    12,     5,     7,     8,
       9,    11,    10,    11,    12,     7,     8,     9,     4,     3,
       2,     3,     4,     5,     6,     5,     4,     5,     5,     6,
       7,     7,     8,     7,     8,     4,     3,     2,     5,     6,
       7,     8,    10,    11,     6,     9,     2,     5,     7,     3,
       2,     4,     2,     5,     7,     2,     5,     6,     7,     8,
      10,    11,    13,    14,     6,     9,     3,     3,     3,     3,
       3,     2,     5,     4,     3,     7,     6,     4,     1,     2,
       4,     3,     3,     3,     5,     4,     2,     2,     2,     2,
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
       4,     3,     4,     4,     6,     2,     3,     4,     5,     8,
       2,     2,    10,    12,     3,     4,     2,     4,     7,     9,
       9,     8,    11,    12,     2,     9,    10,     9,    10,     9,
      10,     7,     7,     7,     8,     8,     7,     7,     8,     8,
       9,    10,    10,    11,    12,    24,     7,     8,     9,    10,
      11,     5,     5,     2,     4,     5,     0,     2,     4,     6,
       8,     2,     3,     2,     3,     2,     2,     3,     2,     2,
       2,     2,     2,     2,     1,     1,     1,     1,     1
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
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     3,   114,   115,   116,
     117,   118,   103,   104,    60,   110,   111,    43,    44,    45,
      40,    41,    42,    39,    50,    55,    56,    57,    28,    29,
      31,    30,    36,    37,    32,    33,    69,    68,    76,   258,
      35,    62,    63,    64,    65,   113,    66,    67,    48,    49,
      54,    51,    52,    53,   130,    94,    95,    96,    93,     6,
     128,    73,    74,    70,    71,    72,    75,    86,    87,    82,
      88,    89,    90,    91,   105,   106,   107,   108,   109,    83,
      84,    97,    98,    99,   100,   101,   102,    22,    21,    25,
      23,    24,    26,    27,     7,    46,    47,     8,     9,    92,
      15,    10,   121,   122,   124,   123,   125,    34,   126,   127,
     131,     5,    13,    12,    11,   129,    14,    17,    18,    19,
      16,   658,   657,    77,   132,    79,    85,    80,    81,    78,
      38,    61,    58,    59,   112,   119,   120,    20,   133,   134,
       0,     0,     0,     0,   944,     0,   636,     0,     0,   946,
     948,   945,   947,     0,     0,     0,     0,   269,     0,     0,
       0,     0,     0,     0,     0,     0,   454,     0,     0,     0,
       0,   634,   470,     0,     0,     0,     0,     0,   191,   392,
       0,   185,   284,   881,     0,     0,     0,     0,     0,   595,
       0,     0,   374,   626,   857,   214,   369,   286,   169,     0,
       0,     0,   812,   817,     0,     0,     0,   259,     0,   630,
       0,     0,   267,     0,     0,   686,     0,     0,     0,     0,
     388,     0,   479,     0,   823,     0,   586,     0,     0,   732,
     734,     0,     0,   469,     0,   217,   333,   766,     0,   137,
       0,     0,     0,   776,     0,     0,     0,     0,   183,   822,
       0,     0,     0,     0,     0,   328,   341,   515,   461,     0,
     466,     0,   783,     0,   475,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   208,   505,     0,   276,     0,
       0,   650,     0,   282,     0,     0,     0,     0,     0,     0,
       0,   731,   179,   181,     0,     0,     0,     0,   335,     0,
       0,   858,   770,     0,   723,     0,     0,     0,     0,   628,
       0,   229,     0,   230,     0,   308,     0,     0,     0,     0,
     652,   646,   753,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   395,     0,   337,     0,     0,     0,   880,
     219,     0,   141,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   262,
       0,   784,     0,     0,     0,   159,   385,     0,     0,     0,
     285,     0,     0,   326,   325,   500,     0,     0,   577,   278,
       0,     0,     0,     0,     0,     0,   736,     0,   495,   162,
     886,     0,   327,     0,     0,     0,   894,   923,     0,     0,
       0,     0,     0,   177,   767,     0,   280,     0,   329,   342,
       0,     0,     0,   592,   931,   936,     1,     2,     4,     0,
       0,     0,     0,     0,     0,   150,   151,   148,   149,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     186,   187,   188,   189,   190,   192,     0,   209,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   287,   288,   289,
       0,     0,     0,   309,   310,   322,     0,     0,     0,   366,
       0,     0,   370,     0,     0,     0,     0,     0,     0,     0,
       0,   405,     0,   416,     0,   419,     0,   422,     0,   425,
       0,   433,   435,   437,   442,     0,   445,     0,   451,     0,
     455,     0,     0,     0,     0,     0,     0,     0,   481,     0,
     530,     0,   560,     0,     0,     0,     0,   590,     0,     0,
       0,   611,     0,   627,   629,   631,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     841,   840,   839,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   875,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   226,   485,     0,   315,     0,     0,   168,   740,
       0,     0,     0,   432,   434,     0,   270,     0,     0,   428,
     430,     0,   431,     0,     0,   448,   449,     0,   756,     0,
     606,   279,     0,   826,     0,   803,   160,   161,   755,     0,
       0,   393,     0,     0,   739,   737,     0,   750,   757,     0,
       0,   735,   805,   806,   804,     0,     0,   813,     0,   818,
       0,     0,   810,   260,   415,   404,   808,     0,   688,   689,
       0,   687,     0,   438,     0,     0,   389,   480,   274,   730,
     729,   733,   607,   748,     0,   749,   705,   706,   361,   527,
       0,     0,   525,     0,   398,   399,     0,     0,     0,     0,
     835,   418,   853,     0,   855,     0,   424,   421,   516,     0,
       0,   463,   462,   465,   478,   493,     0,   476,     0,     0,
       0,   365,     0,     0,   275,   782,   781,     0,   506,     0,
       0,   223,   754,     0,     0,     0,   651,   529,     0,     0,
     793,     0,   792,     0,   791,   790,   728,   727,   800,   807,
       0,     0,   334,   261,   771,   610,     0,   264,   777,     0,
     884,   231,   232,     0,   471,   473,     0,   724,   714,   713,
     769,   801,     0,     0,     0,     0,     0,   396,     0,   339,
     336,   459,     0,     0,   751,   726,   725,   225,   268,   694,
     697,   695,   696,   698,   699,   701,   700,   702,   692,   693,
       0,     0,     0,     0,     0,   780,   409,   410,     0,   717,
       0,   265,   715,     0,   266,   559,     0,     0,   501,   779,
     785,     0,   224,   228,   227,   752,   765,   824,     0,   257,
       0,   496,     0,     0,   788,   746,     0,     0,     0,     0,
       0,   147,   561,     0,     0,     0,   608,   609,   576,     0,
     768,   281,   376,   377,     0,   589,   381,   382,     0,     0,
       0,     0,     0,     0,     0,   153,     0,     0,     0,     0,
       0,     0,     0,   176,     0,     0,   174,   175,   173,   172,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   197,
       0,   199,   201,     0,   203,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   249,   250,   255,   256,
     252,   233,     0,   247,   254,     0,     0,     0,   305,   297,
       0,     0,   307,     0,   290,     0,   299,   296,     0,     0,
       0,     0,     0,   311,     0,   319,     0,   321,   323,     0,
     373,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   443,     0,     0,     0,
       0,     0,     0,     0,     0,   484,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   649,   945,     0,     0,     0,     0,     0,
     660,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   859,     0,     0,     0,   871,     0,     0,     0,
       0,   860,     0,     0,     0,   869,     0,     0,     0,   926,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   933,   935,   932,   941,   938,
     939,   942,   943,   940,   937,   486,     0,     0,   316,     0,
       0,   741,     0,   742,     0,     0,   271,   272,     0,   429,
       0,     0,     0,     0,   761,   691,   827,     0,   760,   762,
       0,     0,   738,   763,   764,   709,   708,   814,   819,   811,
     820,   809,   690,     0,   439,   440,   850,     0,   332,     0,
     528,     0,   772,     0,   526,   400,     0,     0,     0,     0,
       0,   854,   856,     0,   518,     0,   517,     0,     0,   522,
       0,   494,     0,   831,   845,     0,     0,     0,     0,   136,
       0,     0,   508,     0,   507,     0,   510,     0,   758,   759,
     802,     0,   797,     0,     0,     0,   796,   703,   704,     0,
     778,     0,     0,   829,   830,   719,   721,   720,   338,   340,
     460,   928,   668,     0,     0,   685,     0,     0,   662,     0,
     670,     0,     0,     0,   411,     0,   718,   716,   330,     0,
       0,     0,   786,     0,     0,     0,     0,     0,   887,   789,
     747,     0,     0,     0,     0,     0,   562,     0,   566,     0,
       0,     0,   575,   378,   380,     0,   383,     0,     0,     0,
       0,     0,     0,   152,     0,     0,   163,   164,   165,   166,
     167,   170,   171,   178,   180,   182,   184,     0,   196,   198,
     200,   202,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   221,   222,   234,     0,   240,     0,   251,   253,   246,
     248,   277,   301,   303,     0,   302,   294,   300,   292,   291,
       0,     0,   295,     0,     0,   320,     0,     0,   621,     0,
       0,     0,   386,     0,   390,     0,     0,     0,   407,     0,
       0,     0,     0,   446,     0,   452,     0,     0,   464,   491,
       0,     0,     0,     0,   477,     0,     0,     0,     0,     0,
       0,     0,     0,   618,     0,     0,     0,     0,     0,     0,
     614,     0,     0,     0,   637,     0,     0,     0,   647,     0,
     653,     0,   659,   661,   664,   666,   663,   667,   665,   707,
     722,   710,   711,   712,     0,     0,     0,   867,     0,   872,
     870,   861,   873,   862,     0,     0,   876,     0,   885,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   924,     0,   934,     0,
       0,   318,     0,     0,   743,     0,   744,   436,   273,   427,
       0,     0,     0,   828,   283,     0,   815,   821,   684,   851,
       0,   331,     0,   773,     0,     0,     0,   842,     0,     0,
       0,     0,   519,     0,     0,   521,   623,   846,     0,   363,
       0,     0,     0,     0,   509,     0,   511,     0,     0,     0,
     795,   794,     0,   135,   348,   344,     0,     0,     0,   669,
     680,   681,     0,   679,     0,   673,     0,   671,   672,   263,
       0,     0,     0,     0,     0,   787,   825,     0,     0,     0,
     155,     0,     0,     0,   143,     0,   563,     0,     0,   567,
       0,   568,     0,   379,     0,     0,   138,     0,     0,   360,
     359,   154,   157,     0,   193,     0,     0,     0,   633,     0,
       0,     0,   215,   218,   220,     0,   236,     0,     0,   242,
       0,   306,   298,     0,     0,     0,     0,     0,     0,     0,
     620,     0,   387,   391,     0,     0,   406,   408,   417,   420,
     423,   426,   447,   453,     0,   492,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   558,     0,     0,   617,   619,
     565,   578,     0,     0,   587,     0,   591,   593,     0,   596,
       0,     0,   613,     0,     0,   642,     0,   643,     0,     0,
     648,     0,     0,   849,   865,     0,     0,   868,   863,     0,
       0,   877,     0,   921,   922,   927,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     925,     0,     0,     0,   317,   324,   745,     0,     0,     0,
       0,   816,     0,   625,     0,   774,     0,     0,     0,   401,
       0,     0,     0,     0,   844,     0,     0,     0,   362,   364,
     368,   848,     0,     0,   513,   799,   798,   349,     0,   351,
     352,   353,     0,   355,   356,   358,     0,   345,     0,   929,
     682,   678,     0,   674,     0,     0,     0,   412,     0,     0,
     503,     0,     0,   498,     0,     0,     0,   156,   564,   569,
       0,     0,     0,   384,   140,   139,   142,     0,     0,     0,
       0,     0,     0,     0,   216,   237,   235,   243,   241,   304,
       0,   343,     0,   312,     0,   367,     0,     0,   375,   394,
     397,   456,     0,   622,   468,     0,     0,     0,     0,     0,
       0,   547,     0,     0,     0,   543,     0,     0,     0,   584,
     579,     0,     0,     0,   604,   597,     0,   612,     0,     0,
     638,     0,   639,     0,     0,     0,   654,     0,   655,     0,
     874,     0,   878,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   450,     0,   852,   624,   775,     0,
     402,     0,   843,   832,     0,   833,     0,   524,     0,   520,
     847,   514,   512,   350,   354,   357,   347,   346,     0,   683,
     675,     0,     0,   413,     0,     0,     0,     0,     0,   144,
     145,   570,   571,     0,   573,     0,   158,     0,     0,     0,
     207,     0,     0,     0,   238,     0,   244,     0,   314,   313,
       0,   371,     0,     0,     0,     0,     0,   482,     0,     0,
       0,     0,   544,     0,     0,     0,   580,     0,   588,   594,
     598,     0,   616,     0,   635,   641,     0,     0,   644,     0,
     866,   864,     0,   888,     0,     0,     0,     0,     0,   901,
     902,     0,   916,     0,   906,   907,     0,     0,     0,     0,
     903,     0,     0,     0,     0,     0,   487,     0,     0,     0,
     403,   834,     0,     0,     0,   523,   930,   676,     0,   414,
     502,     0,   497,     0,   146,   572,   574,     0,   194,     0,
       0,   210,     0,   239,   245,   293,     0,   457,     0,     0,
       0,     0,     0,     0,   555,     0,   548,     0,     0,     0,
       0,     0,   581,     0,     0,   599,     0,     0,   615,   645,
     640,   656,     0,     0,     0,   891,     0,     0,   904,   917,
       0,   908,     0,     0,   909,     0,     0,   905,     0,   488,
       0,   441,   444,     0,     0,   836,   837,   677,   504,   499,
     195,     0,     0,     0,   372,     0,   467,   472,   474,   483,
       0,   556,     0,   549,     0,     0,     0,     0,     0,   585,
       0,   605,     0,   879,   890,   889,     0,   895,     0,   897,
       0,   918,     0,     0,     0,   910,   899,     0,     0,   489,
     882,   838,     0,     0,   211,     0,     0,     0,   557,   550,
       0,     0,     0,     0,   537,     0,   582,     0,   600,     0,
       0,   896,   898,   919,     0,   911,   912,     0,   900,   490,
       0,     0,     0,     0,   458,   552,     0,     0,     0,     0,
       0,     0,   583,   601,     0,   892,     0,   920,     0,   913,
       0,   883,     0,   632,   212,     0,   553,     0,   551,     0,
       0,     0,     0,     0,   893,     0,   914,   204,     0,   213,
     554,     0,   540,     0,   546,     0,   602,     0,     0,     0,
       0,   541,     0,     0,   603,     0,   205,     0,     0,   542,
       0,   538,     0,     0,   206,     0,     0,     0,     0,   545,
     539,     0,   531,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   532,   533,     0,     0,     0,   534,
       0,     0,     0,     0,     0,     0,     0,     0,   535,   536,
     915
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   234,   235,   236,   237,   238,   239,   240,   241,   655,
     656,  1061,   657,   242,   243,   244,   245,   246,   247,   248,
     249,   250,   251,   680,  1684,   681,   682,   683,   684,  1109,
    1112,   252,   687,   253,   254,   255,   256,   257,   258,   259,
     260,   261,   262,   694,   263,   264,   265,   266,   267,   268,
     269,   707,  1137,  1141,   270,   713,   714,   271,   718,   272,
     273,   274,   275,   276,   277,   278,   279,   280,   281,   998,
     282,   283,   708,   284,  1635,  1836,   658,   285,   286,   287,
     719,   288,   289,   290,   291,  1072,  1073,   292,  1076,  1077,
     293,   294,   295,   296,   297,   914,   915,   298,   731,  1487,
     299,  1026,  1027,   300,   733,   301,   735,   302,   737,   303,
     739,   839,   840,   304,   305,   306,   307,   308,   309,   310,
     311,   744,   312,   746,   313,   314,   315,   748,   316,   750,
     317,   318,   319,   320,   321,   322,   323,   324,   822,  1498,
     934,   325,  1051,   326,  1038,   327,   948,   949,  1342,   328,
     928,   929,  1324,   329,   908,   330,   760,   331,   332,   333,
     334,   335,   336,   337,   767,   338,   339,   340,   341,   771,
     762,  1512,   823,  1499,   935,   909,   342,   343,   344,   685,
     345,   346,   347,   348,   349,   350,  1209,   351,   352,   353,
     886,   354,   436,   355,   920,  1333,  1334,   356,  1305,  1306,
     357,   922,   358,   924,   359,   360,   806,   361,   362,   363,
     364,   365,   366,  1561,   367,   368,   814,   369,   821,  1488,
    1335
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -1611
static const yytype_int16 yypact[] =
{
    4792,   -62,    31,   -61,    24,    47,    -4,   842,   842,   842,
     251,    66,    87,   111,   117,    24,    24,   128,   143,   -96,
      24,    24,   160,   190,   212,    24,   158,   191,   201,   267,
     360,   153,   362,   378,   383,   123,   842,   338,   842,   419,
     335,   349,   429,   446,   447,   457,   466,   472,   478,   407,
      24,    24,   294,   206,   491,   539,   551,   556,   580,   -51,
      24,    24,   247,   274,   582,   440,   591,    41,   600,   603,
      24,    24,   612,   842,   614,   615,   630,   842,   642,   842,
     483,   653,   317,   332,   654,   661,   662,   666,   678,   687,
     698,   712,   713,   719,   722,   -69,   232,   723,   729,   734,
      24,   735,   401,   737,   743,    24,   467,   749,   751,   754,
     860,   756,   688,    24,   404,   760,   764,   375,   -21,    43,
     781,   788,   794,    24,    24,    24,    24,   409,    24,   795,
     422,   796,   802,   814,    24,    24,   817,   929,   437,   452,
     831,   835,    24,    24,   837,   853,   859,   861,    24,   -26,
      24,   842,    24,    24,   842,   842,    24,   460,   468,   933,
     866,   876,   888,   739,   900,    77,   902,   842,   842,    24,
      24,    24,   842,    24,    24,   812,    24,    24,    24,   904,
     514,   912,    24,   917,   842,   918,   919,   842,   842,   842,
     927,   936,   942,   944,   947,   951,   842,    24,   965,   976,
    1082,   986,   987,    24,    24,   842,   988,   867,  1004,  1005,
    -135,  1013,  1081,   842,  1019,  1030,  1047,    24,    24,   842,
      24,    24,  1048,  -107,  1049,   842,  1051,  1052,  1054,  1056,
    1058,  1063,  1066,  1069,  1272,  4392, -1611, -1611, -1611,  1154,
      24,     9, -1611,  1091, -1611,   -59,    24,   842,   842,   842,
     830,    24,    24,    24,   842,  1178, -1611, -1611, -1611, -1611,
   -1611, -1611,    62, -1611,  1099, -1611, -1611, -1611, -1611,     3,
     273,   -13, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611,    24, -1611,   -79,
      24, -1611, -1611,   -65,   -54,    24,    24, -1611,    24, -1611,
      24,    24,    24,    24, -1611, -1611,    24,    24,    24, -1611,
   -1611,    24,    24, -1611, -1611,    24,    24, -1611,  1104,    24,
      24,    24,  1105, -1611,   -49, -1611, -1611, -1611, -1611, -1611,
      24,    24,    24, -1611, -1611,    24,    24,    24,    24,    24,
   -1611,    24,    24,    24,    24,    24,    24,   -43, -1611,    24,
    1224,   828,   -25, -1611, -1611,    -5,   842, -1611, -1611, -1611,
       7, -1611, -1611,  1093,    24,  -126,    24, -1611,   137,   982,
      24,  1094,  1291,  1299, -1611,  1106, -1611,   -74,  -106, -1611,
   -1611, -1611, -1611,  1111,  1112,   842,   516, -1611,   842,    24,
      24,   842,   842,  1113,  1114,    24, -1611,   325,  1115,  1117,
    1036, -1611, -1611,   270,  1118,  1120,  1121,   243, -1611, -1611,
    1122, -1611,   842,  1229,  1125,   528,  1130,   306,  1132, -1611,
    1135,  1137, -1611, -1611, -1611, -1611, -1611, -1611, -1611,  1138,
     842,    24,   842,  1262,   842,   574,   136, -1611,  1140, -1611,
      24,    24, -1611,   842,  1142, -1611,  1144,   230,   538,  1145,
   -1611,  1147, -1611,  1149, -1611,  1150, -1611,  1152,  1153, -1611,
   -1611,  1155,  1156, -1611,  1158, -1611,   842, -1611,  1164, -1611,
    1165,  1166,    24, -1611,   842,    24,  1167,    24, -1611,   572,
      24,   842,   842,    24,    24, -1611, -1611,    24,    24,  1169,
   -1611,  1173, -1611,    24,    24,  1174,   842,    24,  1175,   842,
      24,  1180,  1181,  1184,   842, -1611,    24,  1186, -1611,   327,
     842, -1611,  1187, -1611,    24,   -52,   400,  1189,  1190,  1191,
    1195, -1611, -1611, -1611,  1196,  1199,    24,   842, -1611,  1202,
    1204, -1611, -1611,  1205, -1611,    24,    24,  1207,   271, -1611,
    1208, -1611,  1209, -1611,  1210, -1611,    24,  1212,  1213,    24,
   -1611, -1611, -1611,  1214,  1218,  1232,  1235,  1239,  1240,  1241,
     842,   842,    24, -1611,  1242,    24,  1243,   573,  1192, -1611,
   -1611,  1244, -1611,  1245,  1246,    24,  1248,  1249,  1250,  1251,
    1252,  1253,  1255,  1256,  1257,  1258,  1259,  1260,   -32, -1611,
     842, -1611,  1263,    24,   288, -1611, -1611,  1265,   300,  1266,
   -1611,    24,   842, -1611, -1611,  1351,  1268,   303, -1611, -1611,
    1269,    24,    24,  1270,  1271,   305, -1611,  1273,  1426, -1611,
   -1611,    24, -1611,    -6,   324,   -80, -1611, -1611,  -103,    24,
    1274,  1275,   595, -1611, -1611,  1278, -1611,  1279, -1611, -1611,
      24,    24,    24, -1611, -1611, -1611, -1611, -1611, -1611,   -39,
    1219,   511,   842,   353,  1221, -1611, -1611, -1611, -1611,  1387,
    1388,  1389,  1390,  1391,  1287,  1393,    24,  1289,  1290,  1295,
    1296,   842,   842,   842,   842,    24,    24,    24,    24,    24,
   -1611,    24,    24,    24,    24, -1611,   -27, -1611,   842,    24,
     842,    25,    26,    46,    20,    24,   842,   357,   842,  1217,
    1300,   842,  1303,  1305,   -35,   842,  1220, -1611, -1611, -1611,
     842,  1306,   842, -1611, -1611, -1611,  1310,  1419,   621, -1611,
     842,    24, -1611,  -102,  1456,    24,   842,    24,   842,   842,
     842, -1611,    24, -1611,    24, -1611,    24, -1611,    24, -1611,
      24, -1611, -1611, -1611, -1611,  1317, -1611,    24, -1611,    24,
   -1611,   842,  1318,   842,   842,   842,  1319,    24, -1611,  -100,
   -1611,    -9, -1611,    24,    24,    24,    24, -1611,    24,    24,
      24, -1611,   842, -1611, -1611, -1611,    24,    24,    24,   855,
     -33,    24,    24,    24,    24,    24,   842,    24,    24,    19,
   -1611, -1611, -1611,   842,  1321,   842,    24,    24,   370,   842,
      24,    24,  1323,   842,   622,  1324, -1611,  1430,    24,  1431,
     -84,    24,  1432,    24,  1328,    24,  1434,  1435,    24,    24,
     842,  1331,    24, -1611,   -91, -1611,   374,   842, -1611, -1611,
    1332,    24,   -95, -1611, -1611,   842, -1611,  1333,   627, -1611,
      24,   842,    24,   842,   842, -1611, -1611,  1334, -1611,  1335,
   -1611, -1611,  1336, -1611,   392, -1611, -1611, -1611, -1611,  1338,
    1339, -1611,    24,  1340, -1611, -1611,  1341, -1611, -1611,  1342,
    1343, -1611, -1611, -1611, -1611,  1344,  1346,   842,   842, -1611,
     163,    24, -1611, -1611, -1611, -1611, -1611,  1347, -1611, -1611,
    1348, -1611,    24, -1611,  1349,   842, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611,   394, -1611, -1611, -1611,    24, -1611,
      24,   640,    24,    24,    24, -1611,    24,   842,   842,   842,
     842, -1611, -1611,  1367, -1611,  1368, -1611, -1611,    24,    24,
     264,    24, -1611, -1611,    24, -1611,   842,    24,   842,   842,
     842, -1611,   842,  1369, -1611, -1611, -1611,   842,    24,    24,
     399, -1611, -1611,  1370,  1371,  1372, -1611, -1611,   178,    24,
   -1611,    24, -1611,   180, -1611, -1611, -1611, -1611, -1611, -1611,
    1373,  1374, -1611, -1611, -1611, -1611,    24, -1611, -1611,  1375,
   -1611, -1611, -1611,    24, -1611, -1611,    24, -1611, -1611, -1611,
   -1611, -1611,   842,   842,  1376,  1377,  1378, -1611,   652, -1611,
   -1611, -1611,  1379,  1380, -1611, -1611, -1611,    24, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
     693,   -50,   704,  -121,   842, -1611,    24, -1611,    24, -1611,
    1384, -1611, -1611,  1385, -1611, -1611,  1386,    24,   272, -1611,
   -1611,   402, -1611,    24,    24, -1611, -1611, -1611,   842, -1611,
      24,   570,  1392,  1394, -1611, -1611,  1395,   842,   842,   842,
     842,   842, -1611,  1397,   842,  -114, -1611, -1611, -1611,  1400,
   -1611, -1611,    24, -1611,   420, -1611,    24, -1611,    24,    24,
      24,   842,  1401,   842,  1403, -1611,   842,    24,  1406,  1407,
    1408,  1410,  1411, -1611,  1412,  1413, -1611, -1611, -1611, -1611,
    1414,  1415,  1417,  1418,  1420,  1421,  1422,  1423,  1424, -1611,
     842, -1611, -1611,    24, -1611,    24,   842,   842,   855,   842,
       8,  1425,    24,    24,    24,    24, -1611,    24, -1611, -1611,
      24, -1611,    24, -1611, -1611,   842,  1427,  1428,   842, -1611,
     842,  1429, -1611,  1433, -1611,  1436, -1611, -1611,  1437,   423,
     842,  1438,   842, -1611,   842, -1611,  1439, -1611, -1611,   842,
   -1611,    24,    24,   424,    24,   842,  1440,   842,  1442,   842,
     842,    24,    24,    24,    24,    24, -1611,   715,   746,   842,
      24,   842,   842,   842,    24, -1611,    24,    24,    24,   842,
     842,   842,   842,    24,    24,    24,    24,    24,    24,    24,
     428,   842,   -34, -1611,   901,   842,  1444,    24,   430,   803,
   -1611,  1445,  1446,  1447,  1449,  1450,  1451,  1452,  1453,  1454,
    1455,   842, -1611,   842,   842,   816, -1611,  1457,  1458,  1460,
    1462, -1611,   432,  1320,   820, -1611,  1464,   842,  1493, -1611,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   825,    24, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611,   -85,    24, -1611,   480,
     842, -1611,   -93, -1611,  1465,  1466, -1611, -1611,  1467, -1611,
    1468,  1469,  1470,   842, -1611, -1611, -1611,  1471, -1611, -1611,
    1472,    24, -1611, -1611, -1611, -1611, -1611,   842, -1611, -1611,
     842, -1611, -1611,  1473, -1611,   842, -1611,   842, -1611,  1474,
   -1611,   842, -1611,   482,    24, -1611,   855,   530,   842,   842,
      24, -1611, -1611,    24, -1611,   266, -1611,    24,   842, -1611,
    1475, -1611,  1476,   842, -1611,   842,   549,   842,   842, -1611,
     842,    24, -1611,   550, -1611,    24, -1611,   -40, -1611, -1611,
   -1611,    24, -1611,  1479,  1482,    24, -1611, -1611, -1611,  1483,
   -1611,   826,    24,   842,   842, -1611, -1611, -1611, -1611, -1611,
   -1611,  1441, -1611,  1484,  1485, -1611,  1486,   268, -1611,   557,
   -1611,  1487,  1488,  1489, -1611,   855, -1611, -1611, -1611,  1490,
      24,   842, -1611,  1491,  1494,  1495,    24,   842, -1611, -1611,
   -1611,   587,   842,   842,  1498,    24, -1611,  -108, -1611,   842,
    -147,  -142, -1611, -1611, -1611,  1499, -1611,    24,    24,   856,
     842,    24,  1501, -1611,  1502,   882, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611,    24, -1611, -1611,
   -1611, -1611,   842,   842,  1503,   842,   842,   842,  1504,  1505,
    1506, -1611, -1611,    12,  1531,    22,  1535, -1611, -1611, -1611,
   -1611, -1611, -1611,   842,  1507, -1611, -1611, -1611, -1611, -1611,
     842,   842, -1611,   842,    24, -1611,   842,  -105, -1611,  1508,
      24,  1509, -1611,  1510, -1611,   842,   842,   893, -1611,   895,
     896,   899,   930, -1611,  1511, -1611,  1512,    24,    24, -1611,
     842,   842,   842,   842,    24,   -83,   842,   842,   842,  1513,
     842,   842,   935, -1611,   952,   297,   592,   968,   609,   445,
   -1611,   613,   842,  1514, -1611,   842,  -101,  1516, -1611,   842,
   -1611,   842, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611,  1517,   975,   842, -1611,  1521, -1611,
   -1611, -1611, -1611, -1611,   623,    24, -1611,  1526, -1611,   842,
    1529,   643,   842,   842,   842,   842,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842, -1611,  1530, -1611,    24,
     -67, -1611,  1532,  1533, -1611,  1537, -1611, -1611, -1611, -1611,
     842,   842,   842, -1611, -1611,  1538,   842, -1611, -1611, -1611,
     842, -1611,   651, -1611,   979,   855,  1540, -1611,   855,    24,
      24,  1541, -1611,   842,   842, -1611, -1611, -1611,   842, -1611,
    1543,  1544,  1546,  1547, -1611,   842, -1611,    24,   198,  1548,
   -1611, -1611,  1550, -1611,  1197, -1611,  1551,    24,  1553, -1611,
   -1611, -1611,  1554, -1611,   993, -1611,  1018, -1611, -1611, -1611,
     855,  1555,   842,  1556,  1557, -1611, -1611,   842,  1558,  1559,
   -1611,    24,    24,    24, -1611,  1561, -1611,  1562,   663, -1611,
     842, -1611,   842, -1611,  1021,  1564, -1611,  1565,  1566,    24,
   -1611, -1611, -1611,    24,    24,    24,   842,   842, -1611,   842,
     842,  1582, -1611, -1611, -1611,    24, -1611,    24,    24, -1611,
      24,   842, -1611,    24,  1583,    24,   675,  1585,    24,   842,
   -1611,  1586, -1611, -1611,  1587,  1590, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611,   680, -1611,  1591,   692,   842,   842,
      24,   842,    -2,   842,   694, -1611,   842,   842, -1611, -1611,
   -1611, -1611,   199,   699, -1611,   842, -1611, -1611,   842, -1611,
     200,   714, -1611,  1024,   842, -1611,  1592, -1611,   842,   -47,
   -1611,   720,  1593, -1611, -1611,    24,  1594, -1611, -1611,   842,
    1595, -1611,   842, -1611, -1611, -1611,   842,   842,   842,   842,
     842,   842,   842,   842,   842,   842,   842,   842,   842,   842,
   -1611,   -44,    24,   842, -1611, -1611, -1611,   842,   842,  1596,
    1536, -1611,  1597, -1611,  1599, -1611,  1600,   842,  1601, -1611,
      24,  1602,   732,   740, -1611,   744,  1605,  1607, -1611, -1611,
   -1611, -1611,  1608,  1609, -1611, -1611, -1611, -1611,    24, -1611,
   -1611, -1611,   842, -1611,   842, -1611,  1551, -1611,  1551,  1480,
   -1611, -1611,  1610, -1611,  1026,   842,  1611, -1611,   842,   842,
   -1611,   842,   842, -1611,    24,  1612,  1613, -1611, -1611, -1611,
    1614,   745,   772, -1611, -1611, -1611, -1611,  1615,    24,   842,
     842,  1616,   842,   842, -1611, -1611,   223, -1611,   490, -1611,
      24, -1611,  1617, -1611,  1618, -1611,    24,  1619, -1611, -1611,
   -1611, -1611,   842, -1611, -1611,   842,   842,   842,    24,  1620,
     842, -1611,   842,   842,   842, -1611,   780,   842,   842, -1611,
   -1611,   797,  1621,  1622, -1611, -1611,   813, -1611,  1045,  1624,
   -1611,  1626, -1611,   842,   842,  1627, -1611,   842, -1611,  1628,
   -1611,  1629,  1358,   -45,   842,   842,  1630,  1631,   842,   818,
    1632,   821,   842,   842,   842,  1633,   842,   842,    24,   842,
      24,  1634,   842,   842, -1611,  1635, -1611, -1611, -1611,  1637,
   -1611,   822, -1611, -1611,   842, -1611,   842, -1611,  1638, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,  1641, -1611,
   -1611,  1065,  1644, -1611,  1645,   842,  1646,   842,  1647, -1611,
   -1611, -1611, -1611,  1648, -1611,  1649, -1611,   842,  1650,   842,
   -1611,   842,   823,    24, -1611,    24, -1611,  1651, -1611, -1611,
     842, -1611,   863,   842,   842,   842,   842, -1611,   842,   889,
     891,   842, -1611,   842,   842,   842, -1611,   458, -1611, -1611,
   -1611,   657, -1611,  1654, -1611, -1611,  1655,  1656, -1611,  1657,
   -1611, -1611,    24, -1611,    24,   842,   229,   842,   842, -1611,
   -1611,  1658, -1611,   898, -1611, -1611,  1659,   842,   842,   915,
   -1611,   842,  1662,    24,  1665,   842, -1611,  1666,  1667,  1539,
   -1611, -1611,   842,  1668,  1669, -1611, -1611, -1611,  1670, -1611,
   -1611,  1671, -1611,  1672, -1611, -1611, -1611,  1673, -1611,   842,
     842, -1611,   842, -1611, -1611, -1611,  1674, -1611,   842,  1675,
    1676,  1677,  1678,   842, -1611,   922, -1611,   923,   842,   842,
     842,   842, -1611,   296,   842, -1611,   359,   842, -1611, -1611,
   -1611, -1611,  1679,  1680,   -78, -1611,   925,   928, -1611, -1611,
     949, -1611,   842,   842, -1611,  1681,  1001, -1611,   842, -1611,
    1682, -1611, -1611,  1684,  1685, -1611, -1611, -1611, -1611, -1611,
   -1611,   842,   842,  1003, -1611,   842, -1611, -1611, -1611, -1611,
     842, -1611,  1689, -1611,  1008,   842,   842,   842,  1022, -1611,
    1050, -1611,  1053, -1611, -1611, -1611,    24, -1611,  1690, -1611,
    1692, -1611,  1055,  1696,  1071, -1611, -1611,  1700,  1701, -1611,
    1598, -1611,   842,   842, -1611,   842,  1703,  1101, -1611, -1611,
      24,   842,   842,   842, -1611,   842, -1611,  1704, -1611,  1123,
     328, -1611, -1611, -1611,  1133, -1611, -1611,  1136, -1611, -1611,
    1705,   842,  1706,  1141, -1611, -1611,  1198,  1708,   842,   842,
     842,   842, -1611, -1611,   842, -1611,   371, -1611,   842, -1611,
    1709, -1611,  1201, -1611, -1611,  1710, -1611,  1711, -1611,   842,
    1203,  1712,   842,  1206, -1611,   842, -1611, -1611,   842, -1611,
   -1611,    24, -1611,  1215, -1611,   842, -1611,  1713,   842,  1216,
      24, -1611,  1716,   246, -1611,   842, -1611,  1717,    24, -1611,
     842, -1611,   842,   842, -1611,  1719,  1720,    29,   842, -1611,
   -1611,   842, -1611,    24,   842,   842,   842,   842,   842,   842,
    1721,  1722,   842,   842, -1611, -1611,   -19,   842,   842, -1611,
      24,   842,   842,   842,   842,  1723,  1724,  1727, -1611, -1611,
   -1611
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
   -1611, -1611, -1611,  1357, -1611, -1611, -1611, -1611, -1611,  1282,
   -1611, -1611,  1575, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,   871,
     950, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611,   957, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1610, -1611,  -153, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611,   869, -1611, -1611,   857,
   -1611, -1611, -1611, -1611, -1611, -1611,   819, -1611,  -216, -1123,
   -1611, -1611,   921, -1611,  1515, -1611,  -234, -1611,  1459, -1611,
    -230,  -639,  1560, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,   237, -1164,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611,   990,  -938, -1611,
   -1611,  1016,  -924, -1611,  -471, -1611,  1448, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611,  1307, -1611, -1611, -1611, -1611,  1416,
    1352,   755,  -280, -1435,  1027,  -902, -1611, -1611, -1611,   531,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611,  -427, -1611, -1611,  -863,  -662, -1611, -1611,   647,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,
   -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611, -1611,  2278,
      -7
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -945
static const yytype_int16 yytable[] =
{
     383,   384,   385,   388,   912,  1326,  1310,   377,   880,   722,
    1310,  1344,   664,  1239,   665,   794,   715,   696,   444,  1374,
    1504,   374,   432,  1057,  1058,  1059,  1837,   432,   716,   416,
     379,   418,  1126,   433,   434,  1020,   379,   371,   433,   434,
     697,   374,   374,  2308,   758,   374,   374,   374,  1187,  1489,
    1490,  1491,  1492,   380,   379,  1669,   374,   651,  1021,   380,
    1671,   379,   374,  1725,   374,   651,   464,   620,   374,  1725,
     468,   651,   470,   742,  1517,   795,  1127,   380,   743,  1188,
     374,  1380,   374,   379,   380,  1900,   796,  1189,  1408,   491,
     741,  2291,   379,   374,  1666,   634,   374,   374,   374,  1062,
     374,  1757,   797,   374,   374,  1450,   380,  1273,   374,  1584,
    1523,  1218,   374,   374,   381,   380,   709,  1381,  1022,   798,
     381,   374,  1120,  1122,  2175,   382,   374,   698,   829,  1363,
    1364,   382,  1695,   488,   414,  1053,   958,   959,   381,  1203,
     370,   373,  1698,  1124,   557,   381,   699,   560,   561,   382,
     960,   445,  1375,   787,   788,  1922,   382,  2043,  1667,   374,
     576,   577,  1626,  1063,   517,   581,   374,   381,  1524,   374,
     415,   374,   374,   374,   799,   397,   381,   594,   382,   830,
     597,   598,   599,  2309,   700,   800,   635,   382,   374,   606,
     374,  1708,  1219,   374,  1128,  1064,  1054,  1129,   615,   701,
    1901,  1279,   792,  1279,   702,   801,   624,   789,  1382,   621,
    1451,   802,   630,  1730,   554,   703,  1758,   803,   637,  1130,
     721,  1162,  1131,  1186,   374,   885,  1976,  1132,  1977,  1792,
    1670,  2292,  1267,   372,   725,  1672,   809,   652,  1579,  1023,
     672,   673,   674,   831,  1725,   727,   921,   690,  2310,   376,
     757,   652,  1948,   926,   653,   654,   778,  1679,  1148,  1274,
    1079,  1585,   804,  1409,  1410,  1411,  2176,   489,   389,   654,
    1923,  1924,  1060,  1240,  1241,  1242,  1243,  1244,   881,  1245,
    1246,  1247,  1248,  1249,  1525,  1250,  1251,  1252,  1253,   390,
    1207,  1190,  1115,   379,  1674,   379,  2293,   704,  1902,  2044,
     374,  2045,   446,  1376,  1627,   881,   705,  1133,  2280,   379,
     457,  1191,   519,   391,   379,  1328,   380,  1328,   380,   392,
    1390,   890,   666,   667,   668,   669,   670,   379,   717,   379,
     395,   379,   380,   379,   379,   379,   805,   380,   882,   706,
     881,   881,   881,  2003,   378,   396,   573,   692,   386,   793,
     380,   379,   380,   374,   380,   409,   380,   380,   380,  1642,
    1134,   432,   400,   379,   432,  1299,   379,  1351,   379,  1355,
     374,   881,   433,   434,   380,   433,   434,   381,   835,   381,
    1352,   838,  1356,   693,   843,   844,   380,   379,   382,   380,
     382,   380,   401,   381,   374,   432,   854,   710,   381,  1612,
    1824,  1909,  1914,  1624,   382,   862,   433,   434,   437,   382,
     380,   381,  1310,  1204,   402,   381,   379,   381,   381,   381,
     379,   374,   382,   875,   382,   877,   382,   879,   382,   382,
     382,  2125,   891,   379,   490,   381,   887,   379,   881,   380,
     892,  1160,  1314,   380,   374,   858,   382,   381,  2281,   450,
     381,   404,   381,   387,   711,   379,   380,   379,   382,   904,
     380,   382,   379,   382,   374,   379,  1329,   911,  1329,   712,
    1643,   381,   853,   978,   923,   925,   452,  1185,   380,   374,
     380,   859,   382,   379,   405,   380,   379,   379,   380,   939,
    1029,   379,   942,   379,   406,   379,   435,   947,  2169,  1741,
     381,   881,  1032,   955,   381,  1040,   380,  1047,   868,   380,
     380,   382,   432,   881,   380,   382,   380,   381,   380,   473,
     971,   381,   374,   433,   434,   432,  1055,   848,   382,   952,
    2235,   979,   382,   812,   475,   860,   433,   434,   407,   381,
     813,   381,  1265,   379,   869,   379,   381,   374,   374,   381,
     382,   374,   382,   994,   995,  1085,   374,   382,   515,  1139,
     382,  2171,   408,   849,   411,   953,   380,   381,   380,   374,
     381,   381,  1226,  2254,   379,   381,  1268,   381,   382,   381,
     412,   382,   382,  1024,   374,   413,   382,  1030,   382,   961,
     382,  1033,   374,   379,  1286,  1036,  1308,   380,   870,   374,
    1041,  1346,   962,   497,  1392,   374,   511,   374,  1048,   417,
    2005,   528,   379,   379,   374,   374,   380,  1056,  1396,   954,
     379,   419,  1414,   502,   532,  1469,  1478,   381,   420,   381,
    1520,   422,  1530,   379,  1553,   380,   380,   374,   382,   541,
     382,   432,   421,   380,  1083,  1084,  1086,  1749,   423,   424,
     379,  1265,   433,   434,   543,   379,   380,  1265,   381,   425,
    2112,   374,   563,   374,  1100,  1101,  1102,  1103,   426,   382,
     565,  1617,   379,   380,   427,   374,   379,   381,   380,  1116,
     428,  1117,  1581,  1119,  1603,   374,   379,  1135,   382,  1138,
    1140,  1142,  1082,   439,  1138,   380,   381,  1204,  1150,   380,
     429,  1617,  1617,  1152,   381,  1154,   379,   382,   382,   380,
     917,   918,   919,  1159,   379,   382,   589,   381,   836,  1166,
     374,  1168,  1169,  1170,   432,   931,   379,  1265,   382,   380,
     865,   937,  1607,   455,   381,   433,   434,   380,   379,   381,
     893,   440,   374,   379,  1179,   382,  1181,  1182,  1183,   380,
     382,  1619,  1346,   441,  1192,   379,   381,   379,   442,  1645,
     381,   380,   379,  1265,  1265,  1200,   380,   382,   374,   374,
     381,   382,  1206,  1208,   374,  1001,   471,   379,   380,  1215,
     380,   382,   443,   379,   454,   380,  1221,   374,  1223,  1660,
     381,  1227,  1228,   456,  1744,   379,  1232,  1068,   381,   374,
     380,   382,   459,   379,   374,   460,   380,   379,   379,   382,
     381,  1747,  1007,  1263,   463,  1752,   465,   466,   380,  1269,
    1270,   382,   381,  1157,  1233,  1768,   380,   381,  1275,  1277,
     380,   380,   467,   382,  1280,   379,  1281,  1282,   382,   381,
     374,   381,  1312,   379,   469,  1774,   381,  1287,  1043,  1044,
     382,   374,   382,  1803,  1368,   472,   477,   382,   380,  2115,
     379,   381,   374,   478,   479,  1859,   380,   381,   480,   675,
    1297,  1298,   382,   773,   774,   775,   379,  1883,   382,   381,
     481,   379,  1891,   380,   379,   379,   379,   381,  1307,   482,
     382,   381,   381,   374,  1894,  1372,  1905,  1309,   382,   380,
     483,  1910,   382,   382,   380,   379,  1378,   380,   380,   380,
    1317,  1318,  1319,  1320,   484,   485,  1915,  1493,   379,   381,
    1628,   486,  1926,  1330,   487,   492,   379,   381,   380,  1332,
     382,   493,  1336,  1337,  1963,  1338,   494,   496,   382,   499,
    1340,   380,  1965,  1347,   381,   500,  1967,  1992,  1495,   380,
     374,   504,   379,   505,   379,   382,   506,   507,   508,   509,
     381,   379,   513,   374,  -944,   381,   514,   374,   381,   381,
     381,   382,   374,   374,  1994,   380,   382,   380,   379,   382,
     382,   382,  2022,   521,   380,   379,   379,  -944,   379,   381,
     522,   379,   676,   677,   678,   679,   523,   531,   534,  2026,
     382,   380,  1204,   374,   535,  1532,   782,   783,   380,   380,
     381,   380,   379,   382,   380,  2030,   536,  1383,  1547,   539,
    2052,   382,  1556,  2055,  2071,  2091,   540,  1576,  1634,   374,
     567,  1391,   571,   545,  1393,   380,   381,   546,   381,   549,
     374,  1394,   374,   374,  1397,   381,   374,   382,  -944,   382,
    1401,  1402,  1403,  1404,  1405,   550,   382,  1407,  1676,  -944,
     784,   551,   381,   552,   379,  2097,   379,  1415,   568,   381,
     381,   379,   381,   382,  1420,   381,  1422,   374,   569,  1424,
     382,   382,   374,   382,  1682,   379,   382,   380,  1742,   380,
     570,  2104,  1750,  2106,   380,  1716,   381,  1718,  1719,   374,
    2129,  1720,   572,  1442,   575,   584,   588,   382,   380,  1445,
    1446,  1448,  1449,   379,   591,   374,   379,  2134,   379,   593,
     595,   596,   374,   785,  2161,  2163,   374,  2177,  1460,   600,
    2179,  1463,  1721,  1464,   379,   786,   380,  1738,   601,   380,
     374,   380,  1470,  1471,   602,  1473,   603,  1474,   381,   604,
     381,  2181,  1476,   605,  1740,   381,  1479,   380,  1481,   382,
    1483,   382,  1485,  1486,   379,   374,   382,   608,   374,   381,
    1746,   374,  1497,   374,  1501,  1502,  1503,  1764,   609,   610,
     382,  1805,  1508,  1509,  1510,  1511,   379,   380,   611,   612,
     616,   617,   374,  1521,  1522,  1841,   379,   381,  1527,   379,
     381,  1531,   381,  2186,   379,  2194,   618,   619,   382,   380,
    2199,   382,   374,   382,  1544,   622,  1545,  1546,   381,   380,
    1843,   625,   380,  1863,  2204,  1554,  1917,   380,  1980,   382,
    1559,   623,   626,  1562,  1563,  1564,  1565,  1566,  1567,  1568,
    1569,  1570,  1571,  1572,  1573,  1574,  1575,  2032,   381,   627,
     633,   636,  2206,   638,   639,  2208,   640,  2213,   641,   382,
     642,   379,  1582,  1583,   379,   643,   379,  2077,   644,   379,
     381,   645,   646,  2216,   649,   691,  1592,   695,   379,   379,
     381,   382,   781,   381,   380,   752,   756,   380,   381,   380,
    1596,   382,   380,  1597,   382,   807,   825,   826,  1307,   382,
    1600,   380,   380,  2225,  1602,   827,  1604,   852,   828,  1606,
    1608,  1609,  1610,   833,   834,   845,   846,   850,  1330,   851,
     855,  1614,   856,   857,   861,  2233,   863,   864,  1618,  1620,
    1621,  1622,   867,  1623,   871,  2237,  1347,   872,  2239,   873,
     874,   878,   883,  2244,   888,   381,   889,   895,   381,   896,
     381,   897,   898,   381,   899,   900,   382,   901,   902,   382,
     903,   382,   381,   381,   382,   815,   905,   906,   907,   913,
    1644,   932,  1646,   382,   382,   933,   938,   941,  1651,   816,
     817,   818,   944,   945,  1654,   819,   946,   820,   951,   956,
    1659,   964,   965,   966,  1661,  1662,  1663,   967,   968,  1037,
    2246,   969,  1668,  2257,   972,  2262,   973,   974,  2266,   977,
     980,   981,   982,  1678,   984,   985,   987,  2271,  2276,  1827,
     988,  1828,  1829,  1830,  1831,  1832,  1833,  1834,  1835,   659,
     660,   661,   662,   663,   989,  1686,  1687,   990,  1689,  1690,
    1691,   991,   992,   993,   997,  1000,  1004,  1005,  1006,  1003,
    1008,  1009,  1010,  1011,  1012,  1013,  1701,  1014,  1015,  1016,
    1017,  1018,  1019,  1703,  1704,  1025,  1705,  1031,  1034,  1707,
    1039,  1042,  1045,  1046,  1050,  1049,  1066,  1067,  1714,  1715,
    1070,  1071,  1087,  1081,  1088,  1089,  1090,  1091,  1092,  1093,
    1094,  1096,  1097,  1726,  1727,  1728,  1729,  1098,  1099,  1732,
    1733,  1734,  1144,  1736,  1737,  1146,  2046,  1147,  1153,  1745,
    1143,  1748,  1155,  1151,  1753,  1754,  1156,  1164,  1756,  1176,
    1180,  1184,  1761,  1222,  1762,  1231,  1235,  1236,  1238,  1255,
    1257,  1259,  1260,  1264,  1271,  1276,  1283,  1284,  1285,  1766,
    1288,  1289,  1291,  1292,  1293,  1294,  1295,  1769,  1296,  1301,
    1302,  1304,  1772,  1111,  1775,  1776,  1777,  1778,  1779,  1780,
    1781,  1782,  1783,  1784,  1785,  1786,  1787,  1788,  1789,  1321,
    1322,  1339,  1348,  1349,  1350,  1357,  1358,  1360,  1365,  1366,
    1367,  1370,  1371,  1797,  1798,  1799,  1386,  1387,  1388,  1801,
    1560,  2042,   648,  1802,  1398,  1804,  1399,  1400,  1808,  1406,
    2113,  1811,  1412,  1421,  2116,  1423,  1815,  1816,  1426,  1427,
    1428,  1817,  1429,  1430,  1431,  1432,  1433,  1434,  1822,  1435,
    1436,  1555,  1437,  1438,  1439,  1440,  1441,  1452,  1697,  1461,
    1462,  1465,  1700,  1955,  1114,  1466,  2143,   790,  1467,  1468,
    1472,  1475,  1482,  1846,  1484,  1848,  1528,  1534,  1535,  1536,
    1851,  1537,  1538,  1539,  1540,  1541,  1542,  1543,  1145,  1549,
    1550,  1860,  1551,  1861,  1552,  1862,  1558,  1586,  1587,  1588,
    1589,  1590,  1591,  1593,  1594,  1598,  1601,  1615,  1616,  1870,
    1871,  1630,  1872,  1873,  1631,  1633,  1639,  1640,  1641,  1647,
    1648,  1649,  1652,  1655,  1879,  2220,  1656,  1657,  1638,  1884,
    1664,  1673,  1887,  1680,  1681,  1688,  1692,  1693,  1694,  1702,
    1710,  1712,  1713,  1722,  1723,  1735,  1755,  1892,  1760,  1763,
    1895,  1896,  1897,  1767,  1899,  1903,  1904,  1906,  1771,  1907,
    1908,  1773,  1790,  1315,  1794,  1795,  1911,  1978,  1912,  1796,
    1800,  1913,  1809,  1814,  1916,  1818,  1819,  1919,  1820,  1821,
    1825,  1921,  1826,  1634,  1927,  1839,  1840,  1847,  1849,  1850,
    1852,  1853,  1931,  1857,  1858,  1933,  1864,  1865,  1866,  1934,
    1935,  1936,  1937,  1938,  1939,  1940,  1941,  1942,  1943,  1944,
    1945,  1946,  1947,  2236,  1874,  1881,  1951,  1885,  1888,  1889,
    1952,  1953,  1890,  1893,  1920,  1928,  1930,  1932,  1954,  1956,
    1959,  1957,  1958,  1960,  1962,  1964,  1966,  1969,  1968,  1970,
    1971,  1972,  1979,  1983,  1989,  1990,  1991,  1996,  2000,  2008,
    2009,  2011,  2017,  2028,  2029,  1974,  2034,  1975,  2035,  2038,
    2040,  2041,  2049,  2050,  2054,  2060,  2066,  2069,  1982,  2070,
    2075,  1984,  1985,  2076,  1986,  1987,  2079,  2080,  2082,  2084,
    2085,  2086,  2088,  2095,  1993,  1995,  2118,  2119,  2120,  2121,
    2128,  2131,  1998,  1999,  2137,  2001,  2002,  2139,  2141,  2142,
    2145,  2146,  2147,  2148,  2149,  2150,  2154,  2156,  2157,  2158,
    2159,  2173,  2174,  2185,  2189,  2012,  2190,  2191,  2013,  2014,
    2015,  2198,  2211,  2018,  2212,  2019,  2020,  2021,  2215,  2023,
    2024,  2025,  2218,  2219,  2027,  2224,  2232,  2241,  2243,  2031,
    2248,  2256,  2259,  2260,  2264,  2274,  2036,  2037,  2279,  2284,
    2039,  2289,  2290,  2304,  2305,  2318,  2319,  2047,  2048,  2320,
     791,  2051,  2053,  1416,  2056,  2057,  2058,  2059,  1341,  2061,
    2062,  1413,  2064,   927,  1323,  2067,  2068,  1384,  1075,  1514,
     842,   975,  1599,  1035,  2072,   884,     0,  2073,     0,  2074,
       0,  1331,   957,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2081,     0,
    2083,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2087,     0,  2089,     0,  2090,  2092,     0,     0,     0,     0,
       0,     0,     0,  2096,     0,  2098,  2099,  2100,  2101,  2102,
       0,  2103,  2105,  2107,  2108,     0,  2109,  2110,  2111,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  2124,     0,
    2126,  2127,     0,     0,     0,     0,  2130,     0,     0,     0,
    2132,  2133,  2135,     0,  2136,     0,     0,     0,  2140,     0,
       0,     0,     0,     0,     0,  2144,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  2151,  2152,     0,  2153,     0,     0,     0,     0,
       0,  2155,     0,     0,     0,     0,  2160,     0,  2162,     0,
    2164,  2165,  2166,  2167,  2168,     0,     0,  2170,     0,     0,
    2172,     0,     0,     0,     0,     0,     0,     0,     0,  2178,
    2180,     0,     0,  2182,     0,  2183,  2184,     0,     0,  2187,
       0,  2188,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  2192,  2193,  2195,     0,  2196,     0,
       0,     0,     0,  2197,     0,     0,     0,  2200,  2201,  2202,
    2203,  2205,     0,  2207,     0,  2209,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2214,     0,  2217,     0,     0,
       0,     0,     0,     0,     0,  2221,  2222,     0,  2223,     0,
    2226,     0,     0,     0,  2228,  2229,  2230,     0,  2231,     0,
       0,     0,  2234,     0,     0,     0,     0,  2238,     0,     0,
    2240,     0,     0,     0,  2242,     0,  2245,     0,     0,  2247,
       0,  2249,  2250,  2251,  2252,     0,     0,  2253,     0,     0,
       0,  2255,     0,     0,     0,  2258,     0,     0,     0,     0,
       0,     0,  2261,  2263,     0,  2265,  2267,     0,  2268,     0,
       0,  2269,     0,     0,     0,     0,  2272,     0,  2273,     0,
       0,  2275,  2277,     0,     0,     0,  2282,     0,  2283,     0,
       0,     0,     0,  2286,     0,  2287,  2288,     0,     0,     0,
    2294,  2295,   375,     0,  2296,     0,     0,  2298,  2299,  2300,
    2301,  2302,  2303,   393,   394,  2306,  2307,     0,   398,   399,
    2311,  2312,     0,   403,  2314,  2315,  2316,  2317,     0,   410,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   430,   431,
       0,   438,     0,     0,     0,     0,     0,   447,   448,   449,
     451,   453,     0,     0,     0,   458,     0,     0,   461,   462,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     474,   476,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   495,     0,
     498,     0,     0,   501,   503,     0,     0,     0,     0,     0,
       0,   510,   512,     0,     0,   516,   518,   520,     0,     0,
       0,   524,   525,   526,   527,   529,   530,     0,   533,     0,
       0,     0,   537,   538,     0,     0,   542,   544,     0,     0,
     547,   548,     0,     0,     0,     0,   553,   555,   556,     0,
     558,   559,     0,     0,   562,   564,   566,     0,     0,     0,
       0,     0,     0,   574,     0,     0,     0,   578,   579,   580,
       0,   582,   583,     0,   585,   586,   587,     0,   590,     0,
     592,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   607,     0,     0,     0,     0,
       0,   613,   614,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   628,   629,     0,   631,   632,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   650,     0,
       0,     0,     0,     0,   671,     0,     0,     0,     0,   686,
     688,   689,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   720,     0,   723,   724,     0,
       0,   726,   728,   729,   730,     0,   732,     0,   734,   736,
     738,   740,     0,     0,   732,   736,   740,     0,     0,   745,
     747,     0,     0,   749,   751,     0,     0,   753,   754,   755,
       0,     0,   759,     0,     0,     0,     0,     0,   761,   763,
     764,     0,     0,   765,   766,   768,   769,   770,     0,   772,
     686,   686,   686,   776,   777,   779,     0,   780,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   808,   810,   811,     0,     0,     0,   824,     0,
       0,     0,     0,     0,     0,     0,   832,     0,     0,     0,
       0,     0,     0,     0,   837,     0,     0,   841,   841,     0,
       0,     0,     0,   847,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   866,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   876,
       0,     0,     0,     0,     0,     0,     0,     0,   734,   732,
       0,     0,     0,     0,     0,     0,   894,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     910,     0,     0,   910,     0,   916,     0,     0,   736,     0,
       0,   740,   738,     0,     0,   930,   824,     0,     0,     0,
       0,   936,   824,     0,     0,   940,     0,     0,   943,     0,
       0,     0,     0,     0,   950,     0,     0,     0,     0,     0,
       0,     0,   761,     0,   963,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   970,     0,     0,     0,     0,     0,
       0,     0,     0,   772,   976,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   983,     0,     0,   986,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     996,     0,     0,   999,     0,  1002,     0,     0,     0,     0,
       0,     0,     0,   824,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1028,     0,     0,     0,     0,     0,     0,     0,   763,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   824,
     824,     0,     0,     0,     0,     0,     0,     0,     0,  1052,
       0,     0,     0,     0,     0,     0,     0,  1065,     0,     0,
    1069,     0,     0,     0,     0,     0,     0,     0,  1074,   768,
    1078,     0,     0,     0,     0,     0,     0,  1080,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1095,     0,     0,     0,     0,     0,
       0,     0,     0,  1104,  1105,  1106,  1107,  1108,     0,  1110,
    1110,  1113,  1113,     0,     0,     0,     0,  1118,     0,  1121,
    1123,  1125,     0,  1136,     0,     0,     0,     0,     0,     0,
       0,     0,  1149,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1158,     0,     0,  1161,
       0,  1163,     0,  1165,     0,  1167,     0,     0,     0,     0,
    1171,     0,  1172,     0,  1173,     0,  1174,     0,  1175,     0,
       0,     0,     0,     0,     0,  1177,     0,  1178,     0,     0,
       0,     0,     0,     0,     0,  1161,     0,  1163,     0,     0,
       0,  1193,  1194,  1195,  1196,     0,  1197,  1198,  1199,     0,
       0,     0,     0,     0,  1201,  1202,   841,  1205,     0,  1210,
    1211,  1212,  1213,  1214,     0,  1216,  1217,  1220,     0,     0,
       0,     0,     0,     0,  1224,  1225,     0,     0,  1229,  1230,
       0,     0,  1234,     0,     0,     0,  1237,     0,     0,  1254,
       0,  1256,     0,  1258,     0,     0,  1261,  1262,     0,     0,
    1266,     0,  1163,     0,     0,     0,     0,     0,     0,  1272,
       0,     0,     0,     0,     0,     0,  1278,     0,   841,     0,
     841,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1290,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1300,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1303,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   910,     0,  1311,  1313,
     910,   910,   916,     0,  1316,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1325,  1327,     0,  1266,
       0,     0,   936,     0,     0,  1266,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1343,  1345,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1353,     0,  1354,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1359,     0,     0,     0,     0,     0,
       0,  1361,     0,     0,  1362,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1369,     0,     0,     0,
       0,     0,     0,     0,     0,  1266,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1373,  1377,
    1379,     0,     0,     0,  1028,     0,  1385,     0,     0,     0,
       0,     0,     0,     0,     0,  1389,     0,     0,     0,     0,
       0,  1266,  1266,     0,     0,     0,     0,     0,  1395,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    1074,     0,     0,     0,  1078,     0,  1417,  1418,  1419,     0,
       0,     0,     0,     0,     0,  1425,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1443,     0,  1444,     0,     0,  1447,     0,     0,     0,
    1453,  1454,  1455,  1456,     0,  1457,     0,     0,  1458,     0,
    1459,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1163,
    1477,     0,  1480,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1494,  1496,     0,  1500,     0,
       0,     0,  1500,     0,  1505,  1506,  1507,     0,     0,     0,
       0,  1513,  1513,  1515,  1516,     0,  1518,  1519,     0,     0,
    1526,     0,     0,     0,     0,  1529,     0,  1533,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1548,     0,     0,     0,     0,     0,     0,
       0,     0,  1557,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1577,  1578,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1163,  1580,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1595,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   910,     0,  1605,     0,     0,     0,  1611,     0,
       0,  1327,     0,  1613,     0,  1613,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1345,
       0,  1625,     0,  1625,     0,     0,     0,     0,     0,  1629,
       0,     0,     0,  1632,     0,     0,     0,     0,     0,  1636,
    1637,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1650,     0,     0,     0,     0,  1653,     0,
       0,     0,     0,     0,  1658,     0,     0,     0,     0,     0,
       0,     0,     0,  1665,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1675,  1677,     0,  1500,
       0,     0,     0,  1683,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1685,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1696,     0,  1699,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1706,     0,     0,  1709,     0,     0,  1711,     0,
       0,     0,     0,     0,     0,  1717,     0,  1717,  1717,  1717,
    1717,     0,     0,     0,     0,  1724,  1500,     0,     0,     0,
       0,     0,  1500,  1731,     0,     0,     0,     0,     0,     0,
    1739,     0,  1739,  1743,     0,  1717,     0,  1751,     0,     0,
       0,     0,     0,     0,  1759,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1765,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1770,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,  1791,  1793,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1806,  1807,     0,     0,  1810,  1812,  1813,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1823,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  1838,     0,     0,     0,     0,
       0,     0,  1842,     0,  1844,     0,     0,     0,  1845,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1854,
    1855,  1856,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1717,     0,     0,     0,     0,  1500,     0,     0,
       0,  1867,  1868,  1869,     0,     0,     0,     0,     0,     0,
       0,     0,     0,  1875,     0,  1876,  1877,     0,  1878,     0,
       0,  1880,     0,  1882,     0,     0,  1886,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1898,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1918,     0,     0,     0,     0,     0,  1925,     0,     0,
       0,     0,     0,  1929,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1949,
    1950,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1961,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1973,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1981,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1988,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1997,     0,     0,     0,
       0,     0,     0,     0,  2004,     0,  2006,     0,  2007,     0,
       0,     0,     0,     0,  2010,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2016,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2033,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  2063,     0,  2065,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  2078,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2093,     0,  2094,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,  2114,     0,     0,     0,  2117,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    2122,     0,  2123,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  2138,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     1,     2,     3,     0,     4,
       0,     5,     6,     0,     0,     7,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,     0,    22,    23,    24,    25,    26,    27,    28,     0,
      29,     0,    30,     0,    31,    32,    33,    34,    35,     0,
       0,    36,    37,     0,    38,    39,    40,    41,     0,     0,
       0,    42,    43,     0,  2210,     0,    44,    45,    46,     0,
       0,     0,    47,     0,    48,     0,    49,    50,    51,    52,
       0,     0,    53,    54,    55,   647,    56,    57,  2227,     0,
       0,    58,    59,     0,     0,    60,    61,     0,     0,     0,
      62,    63,     0,    64,    65,     0,     0,     0,    66,    67,
      68,    69,    70,    71,     0,    72,    73,    74,    75,    76,
      77,    78,     0,    79,    80,    81,     0,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
       0,     0,     0,     0,     0,     0,     0,    95,    96,  2270,
      97,    98,     0,    99,   100,   101,   102,   103,  2278,   104,
       0,   105,   106,   107,     0,     0,  2285,     0,   108,     0,
     109,   110,   111,   112,   113,   114,     0,   115,   116,     0,
       0,  2297,     0,     0,   117,     0,   118,     0,   119,   120,
       0,     0,   121,   122,   123,   124,   125,   126,  2313,   127,
       0,     0,   128,     0,     0,   129,     0,   130,   131,   132,
     133,   134,     0,   135,     0,   136,   137,   138,   139,   140,
       0,   141,   142,   143,     0,   144,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   145,   146,   147,   148,
       0,   149,     0,   150,     0,   151,     0,     0,   152,   153,
     154,   155,   156,   157,   158,     0,   159,     0,   160,     0,
     161,   162,     0,     0,   163,   164,     0,     0,     0,     0,
     165,     0,   166,     0,   167,   168,   169,   170,   171,   172,
     173,   174,   175,   176,   177,   178,   179,     0,   180,     0,
     181,   182,   183,   184,     0,     0,     0,     0,     0,   185,
       0,     0,     0,     0,   186,     0,   187,   188,   189,     0,
     190,   191,   192,   193,   194,   195,   196,   197,   198,     0,
       0,     0,     0,     0,   199,     0,     0,   200,   201,   202,
     203,     0,     0,   204,   205,   206,   207,   208,     0,   209,
       0,     0,     0,     0,     0,   210,   211,     0,     0,   212,
       0,     0,   213,   214,     0,     0,     0,   215,   216,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   217,   218,   219,   220,     0,
       0,     0,   221,     0,     0,     0,     0,     0,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,     0,     0,
       0,     0,     0,   232,   233,     1,     2,     3,     0,     4,
       0,     5,     6,     0,     0,     7,     0,     8,     9,    10,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,     0,    22,    23,    24,    25,    26,    27,    28,     0,
      29,     0,    30,     0,    31,    32,    33,    34,    35,     0,
       0,    36,    37,     0,    38,    39,    40,    41,     0,     0,
       0,    42,    43,     0,     0,     0,    44,    45,    46,     0,
       0,     0,    47,     0,    48,     0,    49,    50,    51,    52,
       0,     0,    53,    54,    55,     0,    56,    57,     0,     0,
       0,    58,    59,     0,     0,    60,    61,     0,     0,     0,
      62,    63,     0,    64,    65,     0,     0,     0,    66,    67,
      68,    69,    70,    71,     0,    72,    73,    74,    75,    76,
      77,    78,     0,    79,    80,    81,     0,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
       0,     0,     0,     0,     0,     0,     0,    95,    96,     0,
      97,    98,     0,    99,   100,   101,   102,   103,     0,   104,
       0,   105,   106,   107,     0,     0,     0,     0,   108,     0,
     109,   110,   111,   112,   113,   114,     0,   115,   116,     0,
       0,     0,     0,     0,   117,     0,   118,     0,   119,   120,
       0,     0,   121,   122,   123,   124,   125,   126,     0,   127,
       0,     0,   128,     0,     0,   129,     0,   130,   131,   132,
     133,   134,     0,   135,     0,   136,   137,   138,   139,   140,
       0,   141,   142,   143,     0,   144,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   145,   146,   147,   148,
       0,   149,     0,   150,     0,   151,     0,     0,   152,   153,
     154,   155,   156,   157,   158,     0,   159,     0,   160,     0,
     161,   162,     0,     0,   163,   164,     0,     0,     0,     0,
     165,     0,   166,     0,   167,   168,   169,   170,   171,   172,
     173,   174,   175,   176,   177,   178,   179,     0,   180,     0,
     181,   182,   183,   184,     0,     0,     0,     0,     0,   185,
       0,     0,     0,     0,   186,     0,   187,   188,   189,     0,
     190,   191,   192,   193,   194,   195,   196,   197,   198,     0,
       0,     0,     0,     0,   199,     0,     0,   200,   201,   202,
     203,     0,     0,   204,   205,   206,   207,   208,     0,   209,
       0,     0,     0,     0,     0,   210,   211,     0,     0,   212,
       0,     0,   213,   214,     0,     0,     0,   215,   216,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   217,   218,   219,   220,     0,
       0,     0,   221,     0,     0,     0,     0,     0,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,     0,     0,
       0,     0,     0,   232,   233
};

static const yytype_int16 yycheck[] =
{
       7,     8,     9,    10,   475,   929,   908,    11,   435,   289,
     912,   949,    71,    97,    73,     8,    29,    14,    69,    69,
    1184,   147,    67,   103,   104,   105,  1636,    67,    41,    36,
      63,    38,    12,    78,    79,    67,    63,     6,    78,    79,
      37,   147,   147,    62,   324,   147,   147,   147,    57,  1172,
    1173,  1174,  1175,    86,    63,   202,   147,    62,    90,    86,
     202,    63,   147,  1498,   147,    62,    73,   202,   147,  1504,
      77,    62,    79,   307,  1197,    68,    56,    86,   308,    88,
     147,   202,   147,    63,    86,    87,    79,    96,   202,    96,
     306,    62,    63,   147,   202,   202,   147,   147,   147,   202,
     147,   202,    95,   147,   147,    97,    86,   202,   147,   202,
     144,    92,   147,   147,   147,    86,   269,   238,   150,   112,
     147,   147,    97,    97,   202,   158,   147,   124,   202,   992,
     993,   158,   120,   202,    11,   141,   188,   189,   147,   778,
     202,   202,   120,    97,   151,   147,   143,   154,   155,   158,
     202,   202,   202,   178,   179,   202,   158,   202,   266,   147,
     167,   168,   202,   266,   185,   172,   147,   147,   202,   147,
      47,   147,   147,   147,   167,   271,   147,   184,   158,   253,
     187,   188,   189,   202,   181,   178,   293,   158,   147,   196,
     147,   296,   173,   147,   174,   298,   202,   177,   205,   196,
     202,   840,   355,   842,   201,   198,   213,   232,   329,   344,
     202,   204,   219,   296,   240,   212,   317,   210,   225,   199,
     299,   323,   202,   323,   147,   441,  1836,   207,  1838,   296,
     377,   202,   323,   202,   299,   377,   362,   242,   323,   271,
     247,   248,   249,   349,  1679,   299,   480,   254,   267,   202,
     299,   242,   296,   483,   245,   260,   299,  1421,   293,   354,
     299,   354,   255,   377,   378,   379,   344,   336,   202,   260,
     317,   318,   352,   357,   358,   359,   360,   361,   142,   363,
     364,   365,   366,   367,   318,   369,   370,   371,   372,   202,
     323,   300,   319,    63,  1417,    63,   267,   294,   300,   344,
     147,   346,   353,   353,   344,   142,   303,   287,    62,    63,
     269,   320,   269,   202,    63,    51,    86,    51,    86,   202,
      48,    91,   381,   382,   383,   384,   385,    63,   341,    63,
     202,    63,    86,    63,    63,    63,   329,    86,   202,   336,
     142,   142,   142,   120,   348,   202,   269,   285,    97,   356,
      86,    63,    86,   147,    86,   202,    86,    86,    86,    91,
     340,    67,   202,    63,    67,   202,    63,   189,    63,   189,
     147,   142,    78,    79,    86,    78,    79,   147,   385,   147,
     202,   388,   202,   321,   391,   392,    86,    63,   158,    86,
     158,    86,   202,   147,   147,    67,   403,   124,   147,  1323,
     202,   202,   202,  1341,   158,   412,    78,    79,   202,   158,
      86,   147,  1314,   147,   202,   147,    63,   147,   147,   147,
      63,   147,   158,   430,   158,   432,   158,   434,   158,   158,
     158,   202,   202,    63,   202,   147,   443,    63,   142,    86,
     447,   721,   913,    86,   147,   202,   158,   147,   202,   202,
     147,   293,   147,   202,   181,    63,    86,    63,   158,   466,
      86,   158,    63,   158,   147,    63,   202,   474,   202,   196,
     202,   147,   202,   202,   481,   482,   202,   757,    86,   147,
      86,   238,   158,    63,   293,    86,    63,    63,    86,   496,
     202,    63,   499,    63,   293,    63,   202,   504,   202,   202,
     147,   142,   202,   510,   147,   202,    86,   202,   202,    86,
      86,   158,    67,   142,    86,   158,    86,   147,    86,   202,
     527,   147,   147,    78,    79,    67,   202,   202,   158,   202,
     202,   538,   158,   396,   202,   292,    78,    79,   271,   147,
     403,   147,   822,    63,   238,    63,   147,   147,   147,   147,
     158,   147,   158,   560,   561,   202,   147,   158,   183,   202,
     158,   202,   202,   238,   202,   238,    86,   147,    86,   147,
     147,   147,   202,   202,    63,   147,   202,   147,   158,   147,
     202,   158,   158,   590,   147,   202,   158,   594,   158,   189,
     158,   598,   147,    63,   202,   602,   202,    86,   292,   147,
     607,   202,   202,   202,   202,   147,   202,   147,   615,   271,
     120,   202,    63,    63,   147,   147,    86,   624,    48,   292,
      63,   202,   202,   156,   202,   202,   202,   147,   293,   147,
     202,   202,   202,    63,   202,    86,    86,   147,   158,   202,
     158,    67,   293,    86,   651,   652,   653,   202,   202,   202,
      63,   931,    78,    79,   202,    63,    86,   937,   147,   202,
     202,   147,   202,   147,   671,   672,   673,   674,   202,   158,
     202,  1333,    63,    86,   202,   147,    63,   147,    86,   686,
     202,   688,   202,   690,   202,   147,    63,   694,   158,   696,
     697,   698,   181,   202,   701,    86,   147,   147,   705,    86,
     293,  1363,  1364,   710,   147,   712,    63,   158,   158,    86,
     138,   139,   140,   720,    63,   158,   202,   147,   202,   726,
     147,   728,   729,   730,    67,   488,    63,  1007,   158,    86,
     202,   494,   202,   293,   147,    78,    79,    86,    63,   147,
     202,   202,   147,    63,   751,   158,   753,   754,   755,    86,
     158,   202,   202,   202,   761,    63,   147,    63,   202,   202,
     147,    86,    63,  1043,  1044,   772,    86,   158,   147,   147,
     147,   158,   779,   780,   147,   202,   293,    63,    86,   786,
      86,   158,   202,    63,   202,    86,   793,   147,   795,   202,
     147,   798,   799,   202,   202,    63,   803,   202,   147,   147,
      86,   158,   202,    63,   147,   202,    86,    63,    63,   158,
     147,   202,   575,   820,   202,   202,   202,   202,    86,   826,
     827,   158,   147,   202,   202,   202,    86,   147,   835,   202,
      86,    86,   202,   158,   841,    63,   843,   844,   158,   147,
     147,   147,   202,    63,   202,   202,   147,   854,   611,   612,
     158,   147,   158,   202,   202,   202,   202,   158,    86,   202,
      63,   147,   147,   202,   202,   202,    86,   147,   202,    39,
     877,   878,   158,   342,   343,   344,    63,   202,   158,   147,
     202,    63,   202,    86,    63,    63,    63,   147,   895,   202,
     158,   147,   147,   147,   202,   202,   202,   904,   158,    86,
     202,   202,   158,   158,    86,    63,   202,    86,    86,    86,
     917,   918,   919,   920,   202,   202,   202,   202,    63,   147,
    1347,   202,   202,   930,   202,   202,    63,   147,    86,   936,
     158,   202,   939,   940,   202,   942,   202,   202,   158,   202,
     947,    86,   202,   950,   147,   202,   202,   202,   202,    86,
     147,   202,    63,   202,    63,   158,   202,    97,   202,   271,
     147,    63,   202,   147,    63,   147,   202,   147,   147,   147,
     147,   158,   147,   147,   202,    86,   158,    86,    63,   158,
     158,   158,   202,   202,    86,    63,    63,    86,    63,   147,
     202,    63,   162,   163,   164,   165,   202,   202,   202,   202,
     158,    86,   147,   147,   202,   202,   178,   179,    86,    86,
     147,    86,    63,   158,    86,   202,   202,  1024,   202,   202,
     202,   158,   202,   202,   202,   202,    97,   202,   202,   147,
      97,  1038,   293,   202,  1041,    86,   147,   202,   147,   202,
     147,  1048,   147,   147,  1051,   147,   147,   158,   147,   158,
    1057,  1058,  1059,  1060,  1061,   202,   158,  1064,   202,   158,
     232,   202,   147,   202,    63,   202,    63,  1074,   202,   147,
     147,    63,   147,   158,  1081,   147,  1083,   147,   202,  1086,
     158,   158,   147,   158,   202,    63,   158,    86,  1515,    86,
     202,   202,  1519,   202,    86,   202,   147,   202,   202,   147,
     202,   202,   202,  1110,   202,   293,   202,   158,    86,  1116,
    1117,  1118,  1119,    63,   202,   147,    63,   202,    63,   202,
     202,   202,   147,   295,   202,   202,   147,   202,  1135,   202,
     202,  1138,   202,  1140,    63,   307,    86,   202,   202,    86,
     147,    86,  1149,  1150,   202,  1152,   202,  1154,   147,   202,
     147,   202,  1159,   202,   202,   147,  1163,    86,  1165,   158,
    1167,   158,  1169,  1170,    63,   147,   158,   202,   147,   147,
     202,   147,  1179,   147,  1181,  1182,  1183,   202,   202,    97,
     158,   202,  1189,  1190,  1191,  1192,    63,    86,   202,   202,
     202,   324,   147,  1200,  1201,   202,    63,   147,  1205,    63,
     147,  1208,   147,   202,    63,   202,   202,   202,   158,    86,
     202,   158,   147,   158,  1221,   202,  1223,  1224,   147,    86,
     202,   202,    86,   202,   202,  1232,   202,    86,   202,   158,
    1237,   150,   202,  1240,  1241,  1242,  1243,  1244,  1245,  1246,
    1247,  1248,  1249,  1250,  1251,  1252,  1253,   202,   147,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   158,
     202,    63,  1269,  1270,    63,   202,    63,   202,   202,    63,
     147,   202,     0,   202,   120,    97,  1283,   178,    63,    63,
     147,   158,    58,   147,    86,   181,   181,    86,   147,    86,
    1297,   158,    86,  1300,   158,   202,   202,     6,  1305,   158,
    1307,    86,    86,   202,  1311,     6,  1313,   271,   202,  1316,
    1317,  1318,  1319,   202,   202,   202,   202,   202,  1325,   202,
     202,  1328,   202,   202,   202,   202,    97,   202,  1335,  1336,
    1337,  1338,   202,  1340,   202,   202,  1343,   202,   202,   202,
     202,    79,   202,   202,   202,   147,   202,   202,   147,   202,
     147,   202,   202,   147,   202,   202,   158,   202,   202,   158,
     202,   158,   147,   147,   158,   383,   202,   202,   202,   202,
    1377,   202,  1379,   158,   158,   202,   202,   202,  1385,   397,
     398,   399,   202,   202,  1391,   403,   202,   405,   202,   202,
    1397,   202,   202,   202,  1401,  1402,  1403,   202,   202,    48,
     202,   202,  1409,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,  1420,   202,   202,   202,   202,   202,   222,
     202,   224,   225,   226,   227,   228,   229,   230,   231,   338,
     339,   340,   341,   342,   202,  1442,  1443,   202,  1445,  1446,
    1447,   202,   202,   202,   202,   202,   202,   202,   202,   257,
     202,   202,   202,   202,   202,   202,  1463,   202,   202,   202,
     202,   202,   202,  1470,  1471,   202,  1473,   202,   202,  1476,
     202,   202,   202,   202,    48,   202,   202,   202,  1485,  1486,
     202,   202,   261,   264,    97,    97,    97,    97,    97,   202,
      97,   202,   202,  1500,  1501,  1502,  1503,   202,   202,  1506,
    1507,  1508,   202,  1510,  1511,   202,  1933,   202,   202,  1516,
     293,  1518,   202,   293,  1521,  1522,    97,    61,  1525,   202,
     202,   202,  1529,   202,  1531,   202,   202,    97,    97,    97,
     202,    97,    97,   202,   202,   202,   202,   202,   202,  1546,
     202,   202,   202,   202,   202,   202,   202,  1554,   202,   202,
     202,   202,  1559,   682,  1561,  1562,  1563,  1564,  1565,  1566,
    1567,  1568,  1569,  1570,  1571,  1572,  1573,  1574,  1575,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,  1590,  1591,  1592,   202,   202,   202,  1596,
      97,   233,   235,  1600,   202,  1602,   202,   202,  1605,   202,
    2027,  1608,   202,   202,  2031,   202,  1613,  1614,   202,   202,
     202,  1618,   202,   202,   202,   202,   202,   202,  1625,   202,
     202,   301,   202,   202,   202,   202,   202,   202,    97,   202,
     202,   202,    97,    97,   684,   202,    97,   355,   202,   202,
     202,   202,   202,  1650,   202,  1652,   202,   202,   202,   202,
    1657,   202,   202,   202,   202,   202,   202,   202,   701,   202,
     202,  1668,   202,  1670,   202,  1672,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,  1686,
    1687,   202,  1689,  1690,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,  1701,    97,   202,   202,   257,  1706,
     202,   202,  1709,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,  1724,   202,   202,
    1727,  1728,  1729,   202,  1731,  1732,  1733,  1734,   202,  1736,
    1737,   202,   202,   914,   202,   202,  1743,   257,  1745,   202,
     202,  1748,   202,   202,  1751,   202,   202,  1754,   202,   202,
     202,  1758,   202,   202,  1761,   202,   202,   202,   202,   202,
     202,   202,  1769,   202,   202,  1772,   202,   202,   202,  1776,
    1777,  1778,  1779,  1780,  1781,  1782,  1783,  1784,  1785,  1786,
    1787,  1788,  1789,  2210,   202,   202,  1793,   202,   202,   202,
    1797,  1798,   202,   202,   202,   202,   202,   202,   202,   202,
    1807,   202,   202,   202,   202,  1812,  1813,   202,  1815,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,  1832,   202,  1834,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,  1845,   202,
     202,  1848,  1849,   202,  1851,  1852,   202,   202,   202,   202,
     202,   202,   202,   202,  1861,  1862,   202,   202,   202,   202,
     202,   202,  1869,  1870,   202,  1872,  1873,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   202,  1892,   202,   202,  1895,  1896,
    1897,   202,   202,  1900,   202,  1902,  1903,  1904,   202,  1906,
    1907,  1908,   202,   202,  1911,   202,   202,   202,   202,  1916,
     202,   202,   202,   202,   202,   202,  1923,  1924,   202,   202,
    1927,   202,   202,   202,   202,   202,   202,  1934,  1935,   202,
     355,  1938,  1939,  1076,  1941,  1942,  1943,  1944,   948,  1946,
    1947,  1072,  1949,   484,   928,  1952,  1953,  1026,   641,  1194,
     390,   535,  1305,   601,  1961,   440,    -1,  1964,    -1,  1966,
      -1,   934,   514,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1985,    -1,
    1987,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1997,    -1,  1999,    -1,  2001,  2002,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  2010,    -1,  2012,  2013,  2014,  2015,  2016,
      -1,  2018,  2019,  2020,  2021,    -1,  2023,  2024,  2025,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2045,    -1,
    2047,  2048,    -1,    -1,    -1,    -1,  2053,    -1,    -1,    -1,
    2057,  2058,  2059,    -1,  2061,    -1,    -1,    -1,  2065,    -1,
      -1,    -1,    -1,    -1,    -1,  2072,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  2089,  2090,    -1,  2092,    -1,    -1,    -1,    -1,
      -1,  2098,    -1,    -1,    -1,    -1,  2103,    -1,  2105,    -1,
    2107,  2108,  2109,  2110,  2111,    -1,    -1,  2114,    -1,    -1,
    2117,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  2126,
    2127,    -1,    -1,  2130,    -1,  2132,  2133,    -1,    -1,  2136,
      -1,  2138,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  2151,  2152,  2153,    -1,  2155,    -1,
      -1,    -1,    -1,  2160,    -1,    -1,    -1,  2164,  2165,  2166,
    2167,  2168,    -1,  2170,    -1,  2172,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2182,    -1,  2184,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2192,  2193,    -1,  2195,    -1,
    2197,    -1,    -1,    -1,  2201,  2202,  2203,    -1,  2205,    -1,
      -1,    -1,  2209,    -1,    -1,    -1,    -1,  2214,    -1,    -1,
    2217,    -1,    -1,    -1,  2221,    -1,  2223,    -1,    -1,  2226,
      -1,  2228,  2229,  2230,  2231,    -1,    -1,  2234,    -1,    -1,
      -1,  2238,    -1,    -1,    -1,  2242,    -1,    -1,    -1,    -1,
      -1,    -1,  2249,  2250,    -1,  2252,  2253,    -1,  2255,    -1,
      -1,  2258,    -1,    -1,    -1,    -1,  2263,    -1,  2265,    -1,
      -1,  2268,  2269,    -1,    -1,    -1,  2273,    -1,  2275,    -1,
      -1,    -1,    -1,  2280,    -1,  2282,  2283,    -1,    -1,    -1,
    2287,  2288,     4,    -1,  2291,    -1,    -1,  2294,  2295,  2296,
    2297,  2298,  2299,    15,    16,  2302,  2303,    -1,    20,    21,
    2307,  2308,    -1,    25,  2311,  2312,  2313,  2314,    -1,    31,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    50,    51,
      -1,    53,    -1,    -1,    -1,    -1,    -1,    59,    60,    61,
      62,    63,    -1,    -1,    -1,    67,    -1,    -1,    70,    71,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      82,    83,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   100,    -1,
     102,    -1,    -1,   105,   106,    -1,    -1,    -1,    -1,    -1,
      -1,   113,   114,    -1,    -1,   117,   118,   119,    -1,    -1,
      -1,   123,   124,   125,   126,   127,   128,    -1,   130,    -1,
      -1,    -1,   134,   135,    -1,    -1,   138,   139,    -1,    -1,
     142,   143,    -1,    -1,    -1,    -1,   148,   149,   150,    -1,
     152,   153,    -1,    -1,   156,   157,   158,    -1,    -1,    -1,
      -1,    -1,    -1,   165,    -1,    -1,    -1,   169,   170,   171,
      -1,   173,   174,    -1,   176,   177,   178,    -1,   180,    -1,
     182,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   197,    -1,    -1,    -1,    -1,
      -1,   203,   204,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   217,   218,    -1,   220,   221,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   240,    -1,
      -1,    -1,    -1,    -1,   246,    -1,    -1,    -1,    -1,   251,
     252,   253,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   287,    -1,   289,   290,    -1,
      -1,   293,   294,   295,   296,    -1,   298,    -1,   300,   301,
     302,   303,    -1,    -1,   306,   307,   308,    -1,    -1,   311,
     312,    -1,    -1,   315,   316,    -1,    -1,   319,   320,   321,
      -1,    -1,   324,    -1,    -1,    -1,    -1,    -1,   330,   331,
     332,    -1,    -1,   335,   336,   337,   338,   339,    -1,   341,
     342,   343,   344,   345,   346,   347,    -1,   349,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   364,   365,   366,    -1,    -1,    -1,   370,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   378,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   386,    -1,    -1,   389,   390,    -1,
      -1,    -1,    -1,   395,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   415,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   431,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   440,   441,
      -1,    -1,    -1,    -1,    -1,    -1,   448,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     472,    -1,    -1,   475,    -1,   477,    -1,    -1,   480,    -1,
      -1,   483,   484,    -1,    -1,   487,   488,    -1,    -1,    -1,
      -1,   493,   494,    -1,    -1,   497,    -1,    -1,   500,    -1,
      -1,    -1,    -1,    -1,   506,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   526,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   535,   536,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   546,    -1,    -1,   549,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     562,    -1,    -1,   565,    -1,   567,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   575,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   593,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   601,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   611,
     612,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   621,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   629,    -1,    -1,
     632,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   640,   641,
     642,    -1,    -1,    -1,    -1,    -1,    -1,   649,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   666,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   675,   676,   677,   678,   679,    -1,   681,
     682,   683,   684,    -1,    -1,    -1,    -1,   689,    -1,   691,
     692,   693,    -1,   695,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   704,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   718,    -1,    -1,   721,
      -1,   723,    -1,   725,    -1,   727,    -1,    -1,    -1,    -1,
     732,    -1,   734,    -1,   736,    -1,   738,    -1,   740,    -1,
      -1,    -1,    -1,    -1,    -1,   747,    -1,   749,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   757,    -1,   759,    -1,    -1,
      -1,   763,   764,   765,   766,    -1,   768,   769,   770,    -1,
      -1,    -1,    -1,    -1,   776,   777,   778,   779,    -1,   781,
     782,   783,   784,   785,    -1,   787,   788,   789,    -1,    -1,
      -1,    -1,    -1,    -1,   796,   797,    -1,    -1,   800,   801,
      -1,    -1,   804,    -1,    -1,    -1,   808,    -1,    -1,   811,
      -1,   813,    -1,   815,    -1,    -1,   818,   819,    -1,    -1,
     822,    -1,   824,    -1,    -1,    -1,    -1,    -1,    -1,   831,
      -1,    -1,    -1,    -1,    -1,    -1,   838,    -1,   840,    -1,
     842,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     862,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   881,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     892,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   908,    -1,   910,   911,
     912,   913,   914,    -1,   916,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   928,   929,    -1,   931,
      -1,    -1,   934,    -1,    -1,   937,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   948,   949,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   959,    -1,   961,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   976,    -1,    -1,    -1,    -1,    -1,
      -1,   983,    -1,    -1,   986,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   998,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1007,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1020,  1021,
    1022,    -1,    -1,    -1,  1026,    -1,  1028,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1037,    -1,    -1,    -1,    -1,
      -1,  1043,  1044,    -1,    -1,    -1,    -1,    -1,  1050,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    1072,    -1,    -1,    -1,  1076,    -1,  1078,  1079,  1080,    -1,
      -1,    -1,    -1,    -1,    -1,  1087,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1113,    -1,  1115,    -1,    -1,  1118,    -1,    -1,    -1,
    1122,  1123,  1124,  1125,    -1,  1127,    -1,    -1,  1130,    -1,
    1132,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1161,
    1162,    -1,  1164,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1177,  1178,    -1,  1180,    -1,
      -1,    -1,  1184,    -1,  1186,  1187,  1188,    -1,    -1,    -1,
      -1,  1193,  1194,  1195,  1196,    -1,  1198,  1199,    -1,    -1,
    1202,    -1,    -1,    -1,    -1,  1207,    -1,  1209,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1225,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1234,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1254,  1255,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1266,  1267,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1291,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1314,    -1,  1316,    -1,    -1,    -1,  1320,    -1,
      -1,  1323,    -1,  1325,    -1,  1327,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1341,
      -1,  1343,    -1,  1345,    -1,    -1,    -1,    -1,    -1,  1351,
      -1,    -1,    -1,  1355,    -1,    -1,    -1,    -1,    -1,  1361,
    1362,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1385,    -1,    -1,    -1,    -1,  1390,    -1,
      -1,    -1,    -1,    -1,  1396,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1405,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1418,  1419,    -1,  1421,
      -1,    -1,    -1,  1425,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1437,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1453,    -1,  1455,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1474,    -1,    -1,  1477,    -1,    -1,  1480,    -1,
      -1,    -1,    -1,    -1,    -1,  1487,    -1,  1489,  1490,  1491,
    1492,    -1,    -1,    -1,    -1,  1497,  1498,    -1,    -1,    -1,
      -1,    -1,  1504,  1505,    -1,    -1,    -1,    -1,    -1,    -1,
    1512,    -1,  1514,  1515,    -1,  1517,    -1,  1519,    -1,    -1,
      -1,    -1,    -1,    -1,  1526,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1545,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1555,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,  1579,  1580,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1604,  1605,    -1,    -1,  1608,  1609,  1610,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1627,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  1637,    -1,    -1,    -1,    -1,
      -1,    -1,  1644,    -1,  1646,    -1,    -1,    -1,  1650,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1661,
    1662,  1663,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1674,    -1,    -1,    -1,    -1,  1679,    -1,    -1,
      -1,  1683,  1684,  1685,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1695,    -1,  1697,  1698,    -1,  1700,    -1,
      -1,  1703,    -1,  1705,    -1,    -1,  1708,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1730,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1753,    -1,    -1,    -1,    -1,    -1,  1759,    -1,    -1,
      -1,    -1,    -1,  1765,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1791,
    1792,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1810,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1828,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1844,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,  1854,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1868,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1876,    -1,  1878,    -1,  1880,    -1,
      -1,    -1,    -1,    -1,  1886,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1898,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1918,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1948,    -1,  1950,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,  1981,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2003,    -1,  2005,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,  2027,    -1,    -1,    -1,  2031,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    2042,    -1,  2044,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  2063,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,     3,     4,     5,    -1,     7,
      -1,     9,    10,    -1,    -1,    13,    -1,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    -1,    30,    31,    32,    33,    34,    35,    36,    -1,
      38,    -1,    40,    -1,    42,    43,    44,    45,    46,    -1,
      -1,    49,    50,    -1,    52,    53,    54,    55,    -1,    -1,
      -1,    59,    60,    -1,  2176,    -1,    64,    65,    66,    -1,
      -1,    -1,    70,    -1,    72,    -1,    74,    75,    76,    77,
      -1,    -1,    80,    81,    82,    83,    84,    85,  2200,    -1,
      -1,    89,    90,    -1,    -1,    93,    94,    -1,    -1,    -1,
      98,    99,    -1,   101,   102,    -1,    -1,    -1,   106,   107,
     108,   109,   110,   111,    -1,   113,   114,   115,   116,   117,
     118,   119,    -1,   121,   122,   123,    -1,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   145,   146,  2261,
     148,   149,    -1,   151,   152,   153,   154,   155,  2270,   157,
      -1,   159,   160,   161,    -1,    -1,  2278,    -1,   166,    -1,
     168,   169,   170,   171,   172,   173,    -1,   175,   176,    -1,
      -1,  2293,    -1,    -1,   182,    -1,   184,    -1,   186,   187,
      -1,    -1,   190,   191,   192,   193,   194,   195,  2310,   197,
      -1,    -1,   200,    -1,    -1,   203,    -1,   205,   206,   207,
     208,   209,    -1,   211,    -1,   213,   214,   215,   216,   217,
      -1,   219,   220,   221,    -1,   223,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   234,   235,   236,   237,
      -1,   239,    -1,   241,    -1,   243,    -1,    -1,   246,   247,
     248,   249,   250,   251,   252,    -1,   254,    -1,   256,    -1,
     258,   259,    -1,    -1,   262,   263,    -1,    -1,    -1,    -1,
     268,    -1,   270,    -1,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,    -1,   286,    -1,
     288,   289,   290,   291,    -1,    -1,    -1,    -1,    -1,   297,
      -1,    -1,    -1,    -1,   302,    -1,   304,   305,   306,    -1,
     308,   309,   310,   311,   312,   313,   314,   315,   316,    -1,
      -1,    -1,    -1,    -1,   322,    -1,    -1,   325,   326,   327,
     328,    -1,    -1,   331,   332,   333,   334,   335,    -1,   337,
      -1,    -1,    -1,    -1,    -1,   343,   344,    -1,    -1,   347,
      -1,    -1,   350,   351,    -1,    -1,    -1,   355,   356,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   373,   374,   375,   376,    -1,
      -1,    -1,   380,    -1,    -1,    -1,    -1,    -1,   386,   387,
     388,   389,   390,   391,   392,   393,   394,   395,    -1,    -1,
      -1,    -1,    -1,   401,   402,     3,     4,     5,    -1,     7,
      -1,     9,    10,    -1,    -1,    13,    -1,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    -1,    30,    31,    32,    33,    34,    35,    36,    -1,
      38,    -1,    40,    -1,    42,    43,    44,    45,    46,    -1,
      -1,    49,    50,    -1,    52,    53,    54,    55,    -1,    -1,
      -1,    59,    60,    -1,    -1,    -1,    64,    65,    66,    -1,
      -1,    -1,    70,    -1,    72,    -1,    74,    75,    76,    77,
      -1,    -1,    80,    81,    82,    -1,    84,    85,    -1,    -1,
      -1,    89,    90,    -1,    -1,    93,    94,    -1,    -1,    -1,
      98,    99,    -1,   101,   102,    -1,    -1,    -1,   106,   107,
     108,   109,   110,   111,    -1,   113,   114,   115,   116,   117,
     118,   119,    -1,   121,   122,   123,    -1,   125,   126,   127,
     128,   129,   130,   131,   132,   133,   134,   135,   136,   137,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   145,   146,    -1,
     148,   149,    -1,   151,   152,   153,   154,   155,    -1,   157,
      -1,   159,   160,   161,    -1,    -1,    -1,    -1,   166,    -1,
     168,   169,   170,   171,   172,   173,    -1,   175,   176,    -1,
      -1,    -1,    -1,    -1,   182,    -1,   184,    -1,   186,   187,
      -1,    -1,   190,   191,   192,   193,   194,   195,    -1,   197,
      -1,    -1,   200,    -1,    -1,   203,    -1,   205,   206,   207,
     208,   209,    -1,   211,    -1,   213,   214,   215,   216,   217,
      -1,   219,   220,   221,    -1,   223,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   234,   235,   236,   237,
      -1,   239,    -1,   241,    -1,   243,    -1,    -1,   246,   247,
     248,   249,   250,   251,   252,    -1,   254,    -1,   256,    -1,
     258,   259,    -1,    -1,   262,   263,    -1,    -1,    -1,    -1,
     268,    -1,   270,    -1,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   284,    -1,   286,    -1,
     288,   289,   290,   291,    -1,    -1,    -1,    -1,    -1,   297,
      -1,    -1,    -1,    -1,   302,    -1,   304,   305,   306,    -1,
     308,   309,   310,   311,   312,   313,   314,   315,   316,    -1,
      -1,    -1,    -1,    -1,   322,    -1,    -1,   325,   326,   327,
     328,    -1,    -1,   331,   332,   333,   334,   335,    -1,   337,
      -1,    -1,    -1,    -1,    -1,   343,   344,    -1,    -1,   347,
      -1,    -1,   350,   351,    -1,    -1,    -1,   355,   356,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   373,   374,   375,   376,    -1,
      -1,    -1,   380,    -1,    -1,    -1,    -1,    -1,   386,   387,
     388,   389,   390,   391,   392,   393,   394,   395,    -1,    -1,
      -1,    -1,    -1,   401,   402
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
     169,   170,   171,   172,   173,   175,   176,   182,   184,   186,
     187,   190,   191,   192,   193,   194,   195,   197,   200,   203,
     205,   206,   207,   208,   209,   211,   213,   214,   215,   216,
     217,   219,   220,   221,   223,   234,   235,   236,   237,   239,
     241,   243,   246,   247,   248,   249,   250,   251,   252,   254,
     256,   258,   259,   262,   263,   268,   270,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     286,   288,   289,   290,   291,   297,   302,   304,   305,   306,
     308,   309,   310,   311,   312,   313,   314,   315,   316,   322,
     325,   326,   327,   328,   331,   332,   333,   334,   335,   337,
     343,   344,   347,   350,   351,   355,   356,   373,   374,   375,
     376,   380,   386,   387,   388,   389,   390,   391,   392,   393,
     394,   395,   401,   402,   407,   408,   409,   410,   411,   412,
     413,   414,   419,   420,   421,   422,   423,   424,   425,   426,
     427,   428,   437,   439,   440,   441,   442,   443,   444,   445,
     446,   447,   448,   450,   451,   452,   453,   454,   455,   456,
     460,   463,   465,   466,   467,   468,   469,   470,   471,   472,
     473,   474,   476,   477,   479,   483,   484,   485,   487,   488,
     489,   490,   493,   496,   497,   498,   499,   500,   503,   506,
     509,   511,   513,   515,   519,   520,   521,   522,   523,   524,
     525,   526,   528,   530,   531,   532,   534,   536,   537,   538,
     539,   540,   541,   542,   543,   547,   549,   551,   555,   559,
     561,   563,   564,   565,   566,   567,   568,   569,   571,   572,
     573,   574,   582,   583,   584,   586,   587,   588,   589,   590,
     591,   593,   594,   595,   597,   599,   603,   606,   608,   610,
     611,   613,   614,   615,   616,   617,   618,   620,   621,   623,
     202,     6,   202,   202,   147,   625,   202,    11,   348,    63,
      86,   147,   158,   626,   626,   626,    97,   202,   626,   202,
     202,   202,   202,   625,   625,   202,   202,   271,   625,   625,
     202,   202,   202,   625,   293,   293,   293,   271,   202,   202,
     625,   202,   202,   202,    11,    47,   626,   271,   626,   202,
     293,   293,   202,   202,   202,   202,   202,   202,   202,   293,
     625,   625,    67,    78,    79,   202,   598,   202,   625,   202,
     202,   202,   202,   202,    69,   202,   353,   625,   625,   625,
     202,   625,   202,   625,   202,   293,   202,   269,   625,   202,
     202,   625,   625,   202,   626,   202,   202,   202,   626,   202,
     626,   293,   202,   202,   625,   202,   625,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   336,
     202,   626,   202,   202,   202,   625,   202,   202,   625,   202,
     202,   625,   156,   625,   202,   202,   202,    97,   202,   271,
     625,   202,   625,   202,   202,   183,   625,   185,   625,   269,
     625,   202,   202,   202,   625,   625,   625,   625,   202,   625,
     625,   202,   202,   625,   202,   202,   202,   625,   625,   202,
      97,   202,   625,   202,   625,   202,   202,   625,   625,   202,
     202,   202,   202,   625,   240,   625,   625,   626,   625,   625,
     626,   626,   625,   202,   625,   202,   625,    97,   202,   202,
     202,   293,   202,   269,   625,   202,   626,   626,   625,   625,
     625,   626,   625,   625,   293,   625,   625,   625,   202,   202,
     625,   202,   625,   202,   626,   202,   202,   626,   626,   626,
     202,   202,   202,   202,   202,   202,   626,   625,   202,   202,
      97,   202,   202,   625,   625,   626,   202,   324,   202,   202,
     202,   344,   202,   150,   626,   202,   202,   202,   625,   625,
     626,   625,   625,   202,   202,   293,   202,   626,   202,   202,
     202,   202,   202,   202,   202,   202,     0,    83,   409,   120,
     625,    62,   242,   245,   260,   415,   416,   418,   482,   338,
     339,   340,   341,   342,    71,    73,   381,   382,   383,   384,
     385,   625,   626,   626,   626,    39,   162,   163,   164,   165,
     429,   431,   432,   433,   434,   585,   625,   438,   625,   625,
     626,    97,   285,   321,   449,   178,    14,    37,   124,   143,
     181,   196,   201,   212,   294,   303,   336,   457,   478,   482,
     124,   181,   196,   461,   462,    29,    41,   341,   464,   486,
     625,   299,   578,   625,   625,   299,   625,   299,   625,   625,
     625,   504,   625,   510,   625,   512,   625,   514,   625,   516,
     625,   504,   512,   516,   527,   625,   529,   625,   533,   625,
     535,   625,   181,   625,   625,   625,   181,   299,   578,   625,
     562,   625,   576,   625,   625,   625,   625,   570,   625,   625,
     625,   575,   625,   585,   585,   585,   625,   625,   299,   625,
     625,    58,   178,   179,   232,   295,   307,   178,   179,   232,
     415,   418,   482,   626,     8,    68,    79,    95,   112,   167,
     178,   198,   204,   210,   255,   329,   612,   202,   625,   362,
     625,   625,   396,   403,   622,   383,   397,   398,   399,   403,
     405,   624,   544,   578,   625,   202,     6,     6,   202,   202,
     253,   349,   625,   202,   202,   626,   202,   625,   626,   517,
     518,   625,   518,   626,   626,   202,   202,   625,   202,   238,
     202,   202,   271,   202,   626,   202,   202,   202,   202,   238,
     292,   202,   626,    97,   202,   202,   625,   202,   202,   238,
     292,   202,   202,   202,   202,   626,   625,   626,    79,   626,
     598,   142,   202,   202,   510,   504,   596,   626,   202,   202,
      91,   202,   626,   202,   625,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   626,   202,   202,   202,   560,   581,
     625,   626,   560,   202,   501,   502,   625,   138,   139,   140,
     600,   512,   607,   626,   609,   626,   516,   514,   556,   557,
     625,   544,   202,   202,   546,   580,   625,   544,   202,   626,
     625,   202,   626,   625,   202,   202,   202,   626,   552,   553,
     625,   202,   202,   238,   292,   626,   202,   562,   188,   189,
     202,   189,   202,   625,   202,   202,   202,   202,   202,   202,
     625,   626,   202,   202,   202,   575,   625,   202,   202,   626,
     202,   202,   202,   625,   202,   202,   625,   202,   202,   202,
     202,   202,   202,   202,   626,   626,   625,   202,   475,   625,
     202,   202,   625,   257,   202,   202,   202,   544,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
      67,    90,   150,   271,   626,   202,   507,   508,   625,   202,
     626,   202,   202,   626,   202,   576,   626,    48,   550,   202,
     202,   626,   202,   544,   544,   202,   202,   202,   626,   202,
      48,   548,   625,   141,   202,   202,   626,   103,   104,   105,
     352,   417,   202,   266,   298,   625,   202,   202,   202,   625,
     202,   202,   491,   492,   625,   570,   494,   495,   625,   299,
     625,   264,   181,   626,   626,   202,   626,   261,    97,    97,
      97,    97,    97,   202,    97,   625,   202,   202,   202,   202,
     626,   626,   626,   626,   625,   625,   625,   625,   625,   435,
     625,   435,   436,   625,   436,   319,   626,   626,   625,   626,
      97,   625,    97,   625,    97,   625,    12,    56,   174,   177,
     199,   202,   207,   287,   340,   626,   625,   458,   626,   202,
     626,   459,   626,   293,   202,   458,   202,   202,   293,   625,
     626,   293,   626,   202,   626,   202,    97,   202,   625,   626,
     578,   625,   323,   625,    61,   625,   626,   625,   626,   626,
     626,   625,   625,   625,   625,   625,   202,   625,   625,   626,
     202,   626,   626,   626,   202,   578,   323,    57,    88,    96,
     300,   320,   626,   625,   625,   625,   625,   625,   625,   625,
     626,   625,   625,   517,   147,   625,   626,   323,   626,   592,
     625,   625,   625,   625,   625,   626,   625,   625,    92,   173,
     625,   626,   202,   626,   625,   625,   202,   626,   626,   625,
     625,   202,   626,   202,   625,   202,    97,   625,    97,    97,
     357,   358,   359,   360,   361,   363,   364,   365,   366,   367,
     369,   370,   371,   372,   625,    97,   625,   202,   625,    97,
      97,   625,   625,   626,   202,   578,   625,   323,   202,   626,
     626,   202,   625,   202,   354,   626,   202,   202,   625,   517,
     626,   626,   626,   202,   202,   202,   202,   626,   202,   202,
     625,   202,   202,   202,   202,   202,   202,   626,   626,   202,
     625,   202,   202,   625,   202,   604,   605,   626,   202,   626,
     581,   625,   202,   625,   560,   502,   625,   626,   626,   626,
     626,   202,   202,   557,   558,   625,   558,   625,    51,   202,
     626,   580,   626,   601,   602,   626,   626,   626,   626,   202,
     626,   553,   554,   625,   554,   625,   202,   626,   202,   202,
     202,   189,   202,   625,   625,   189,   202,   202,   202,   625,
     202,   625,   625,   601,   601,   202,   202,   202,   202,   625,
     202,   202,   202,   625,    69,   202,   353,   625,   202,   625,
     202,   238,   329,   626,   508,   625,   202,   202,   202,   625,
      48,   626,   202,   626,   626,   625,    48,   626,   202,   202,
     202,   626,   626,   626,   626,   626,   202,   626,   202,   377,
     378,   379,   202,   492,   202,   626,   495,   625,   625,   625,
     626,   202,   626,   202,   626,   625,   202,   202,   202,   202,
     202,   202,   202,   202,   202,   202,   202,   202,   202,   202,
     202,   202,   626,   625,   625,   626,   626,   625,   626,   626,
      97,   202,   202,   625,   625,   625,   625,   625,   625,   625,
     626,   202,   202,   626,   626,   202,   202,   202,   202,   202,
     626,   626,   202,   626,   626,   202,   626,   625,   202,   626,
     625,   626,   202,   626,   202,   626,   626,   505,   625,   505,
     505,   505,   505,   202,   625,   202,   625,   626,   545,   579,
     625,   626,   626,   626,   545,   625,   625,   625,   626,   626,
     626,   626,   577,   625,   577,   625,   625,   505,   625,   625,
     202,   626,   626,   144,   202,   318,   625,   626,   202,   625,
     202,   626,   202,   625,   202,   202,   202,   202,   202,   202,
     202,   202,   202,   202,   626,   626,   626,   202,   625,   202,
     202,   202,   202,   202,   626,   301,   202,   625,   202,   626,
      97,   619,   626,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   202,   625,   625,   323,
     625,   202,   626,   626,   202,   354,   202,   202,   202,   202,
     202,   202,   626,   202,   202,   625,   626,   626,   202,   605,
     626,   202,   626,   202,   626,   625,   626,   202,   626,   626,
     626,   625,   558,   625,   626,   202,   202,   602,   626,   202,
     626,   626,   626,   626,   554,   625,   202,   344,   598,   625,
     202,   202,   625,   202,   202,   480,   625,   625,   257,   202,
     202,   202,    91,   202,   626,   202,   626,   202,   202,   202,
     625,   626,   202,   625,   626,   202,   202,   202,   625,   626,
     202,   626,   626,   626,   202,   625,   202,   266,   626,   202,
     377,   202,   377,   202,   505,   625,   202,   625,   626,   545,
     202,   202,   202,   625,   430,   625,   626,   626,   202,   626,
     626,   626,   202,   202,   202,   120,   625,    97,   120,   625,
      97,   626,   202,   626,   626,   626,   625,   626,   296,   625,
     202,   625,   202,   202,   626,   626,   202,   625,   202,   202,
     202,   202,   202,   202,   625,   579,   626,   626,   626,   626,
     296,   625,   626,   626,   626,   202,   626,   626,   202,   625,
     202,   202,   598,   625,   202,   626,   202,   202,   626,   202,
     598,   625,   202,   626,   626,   202,   626,   202,   317,   625,
     202,   626,   626,   202,   202,   625,   626,   202,   202,   626,
     625,   202,   626,   202,   202,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   626,   626,
     202,   625,   296,   625,   202,   202,   202,   626,   626,   626,
     202,   626,   626,   202,   626,   202,   625,   625,   626,   202,
     625,   626,   625,   625,   202,   626,   626,   626,   202,   202,
     202,   202,   626,   625,   202,   202,   202,   222,   224,   225,
     226,   227,   228,   229,   230,   231,   481,   480,   625,   202,
     202,   202,   625,   202,   625,   625,   626,   202,   626,   202,
     202,   626,   202,   202,   625,   625,   625,   202,   202,   202,
     626,   626,   626,   202,   202,   202,   202,   625,   625,   625,
     626,   626,   626,   626,   202,   625,   625,   625,   625,   626,
     625,   202,   625,   202,   626,   202,   625,   626,   202,   202,
     202,   202,   626,   202,   202,   626,   626,   626,   625,   626,
      87,   202,   300,   626,   626,   202,   626,   626,   626,   202,
     202,   626,   626,   626,   202,   202,   626,   202,   625,   626,
     202,   626,   202,   317,   318,   625,   202,   626,   202,   625,
     202,   626,   202,   626,   626,   626,   626,   626,   626,   626,
     626,   626,   626,   626,   626,   626,   626,   626,   296,   625,
     625,   626,   626,   626,   202,    97,   202,   202,   202,   626,
     202,   625,   202,   202,   626,   202,   626,   202,   626,   202,
     202,   202,   202,   625,   626,   626,   480,   480,   257,   202,
     202,   625,   626,   202,   626,   626,   626,   626,   625,   202,
     202,   202,   202,   626,   202,   626,   202,   625,   626,   626,
     202,   626,   626,   120,   625,   120,   625,   625,   202,   202,
     625,   202,   626,   626,   626,   626,   625,   202,   626,   626,
     626,   626,   202,   626,   626,   626,   202,   626,   202,   202,
     202,   626,   202,   625,   202,   202,   626,   626,   202,   626,
     202,   202,   233,   202,   344,   346,   598,   626,   626,   202,
     202,   626,   202,   626,   202,   202,   626,   626,   626,   626,
     202,   626,   626,   625,   626,   625,   202,   626,   626,   202,
     202,   202,   626,   626,   626,   202,   202,   202,   625,   202,
     202,   626,   202,   626,   202,   202,   202,   626,   202,   626,
     626,   202,   626,   625,   625,   202,   626,   202,   626,   626,
     626,   626,   626,   626,   202,   626,   202,   626,   626,   626,
     626,   626,   202,   598,   625,   202,   598,   625,   202,   202,
     202,   202,   625,   625,   626,   202,   626,   626,   202,   202,
     626,   202,   626,   626,   202,   626,   626,   202,   625,   202,
     626,   202,   202,    97,   626,   202,   202,   202,   202,   202,
     202,   626,   626,   626,   202,   626,   202,   202,   202,   202,
     626,   202,   626,   202,   626,   626,   626,   626,   626,   202,
     626,   202,   626,   202,   202,   202,   344,   202,   626,   202,
     626,   202,   626,   626,   626,   202,   202,   626,   626,   202,
     202,   202,   626,   626,   202,   626,   626,   626,   202,   202,
     626,   626,   626,   626,   202,   626,   202,   626,   202,   626,
     625,   202,   202,   202,   626,   202,   202,   626,   202,   202,
      97,   626,   626,   626,   202,   202,   626,   625,   626,   626,
     626,   626,   202,   202,   626,   202,   598,   202,   626,   202,
     626,   202,   626,   202,   202,   626,   202,   626,   202,   626,
     626,   626,   626,   626,   202,   626,   202,   202,   626,   202,
     202,   626,   202,   626,   202,   626,   202,   626,   626,   626,
     625,   202,   626,   626,   202,   626,   202,   626,   625,   202,
      62,   202,   626,   626,   202,   625,   626,   626,   626,   202,
     202,    62,   202,   267,   626,   626,   626,   625,   626,   626,
     626,   626,   626,   626,   202,   202,   626,   626,    62,   202,
     267,   626,   626,   625,   626,   626,   626,   626,   202,   202,
     202
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
#line 150 "p.y"
    { 
          return 0;
         ;}
    break;

  case 6:
#line 161 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; delete (yyvsp[(1) - (1)].bclist); ;}
    break;

  case 7:
#line 163 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 14:
#line 171 "p.y"
    {;}
    break;

  case 18:
#line 176 "p.y"
    {;}
    break;

  case 19:
#line 178 "p.y"
    {;}
    break;

  case 25:
#line 185 "p.y"
    {;}
    break;

  case 43:
#line 204 "p.y"
    { domain->setMFTT((yyvsp[(1) - (1)].mftval)); ;}
    break;

  case 44:
#line 206 "p.y"
    { domain->setMPTT((yyvsp[(1) - (1)].mptval)); ;}
    break;

  case 45:
#line 208 "p.y"
    { domain->setHFTT((yyvsp[(1) - (1)].hftval)); ;}
    break;

  case 70:
#line 235 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 71:
#line 237 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 72:
#line 239 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 73:
#line 241 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 74:
#line 243 "p.y"
    { if(geoSource->setDirichletFluid((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 75:
#line 245 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 92:
#line 263 "p.y"
    { if(domain->setComplexNeuman((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; ;}
    break;

  case 94:
#line 266 "p.y"
    { if(domain->setComplexDirichlet((yyvsp[(1) - (1)].cxbclist)->n,(yyvsp[(1) - (1)].cxbclist)->d) < 0) return -1; ;}
    break;

  case 105:
#line 278 "p.y"
    { if(geoSource->setDirichlet((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 106:
#line 280 "p.y"
    { if(geoSource->setNeuman((yyvsp[(1) - (1)].bclist)->n,(yyvsp[(1) - (1)].bclist)->d) < 0) return -1; ;}
    break;

  case 112:
#line 287 "p.y"
    {;}
    break;

  case 121:
#line 298 "p.y"
    {;}
    break;

  case 122:
#line 300 "p.y"
    {;}
    break;

  case 123:
#line 302 "p.y"
    {;}
    break;

  case 124:
#line 304 "p.y"
    {;}
    break;

  case 125:
#line 306 "p.y"
    {;}
    break;

  case 126:
#line 308 "p.y"
    {;}
    break;

  case 127:
#line 310 "p.y"
    {;}
    break;

  case 128:
#line 312 "p.y"
    {;}
    break;

  case 135:
#line 322 "p.y"
    { domain->solInfo().noninpc = true;
            sfem->setOrder((yyvsp[(3) - (5)].ival)); 
            domain->solInfo().nsample = (yyvsp[(4) - (5)].ival);
          ;}
    break;

  case 136:
#line 329 "p.y"
    { domain->solInfo().inpc = true;
            sfem->setOrder((yyvsp[(3) - (4)].ival));
          ;}
    break;

  case 138:
#line 336 "p.y"
    { if ((yyvsp[(2) - (5)].ival) == OutputInfo::Attribute)  geoSource->setGroupAttribute((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival)-1);
          else if ((yyvsp[(2) - (5)].ival) == OutputInfo::Nodal)  geoSource->setNodeGroup((yyvsp[(3) - (5)].ival)-1, (yyvsp[(4) - (5)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Group Type: %d\n", (yyvsp[(2) - (5)].ival));  exit(-1); }
        ;}
    break;

  case 139:
#line 341 "p.y"
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

  case 140:
#line 353 "p.y"
    { if ((yyvsp[(2) - (6)].ival) == OutputInfo::Nodal) geoSource->setSurfaceGroup((yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].ival));
          else  {  fprintf(stderr, " ### AS.ERR: Unrecognized Surface Group Type: %d\n", (yyvsp[(2) - (6)].ival));  exit(-1); }
        ;}
    break;

  case 142:
#line 360 "p.y"
    { geoSource->setGroupRandomProperty((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].rprop),(yyvsp[(4) - (6)].fval),(yyvsp[(5) - (6)].fval)); ;}
    break;

  case 143:
#line 364 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (5)].fval)); ;}
    break;

  case 144:
#line 368 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies1(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); ;}
    break;

  case 145:
#line 370 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (7)].fval)); domain->addFrequencies2(2.0*PI*(yyvsp[(4) - (7)].fval), 2.0*PI*(yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); ;}
    break;

  case 146:
#line 372 "p.y"
    { geoSource->setImpe((yyvsp[(4) - (8)].fval)); domain->addFrequencies(2.0*PI*(yyvsp[(4) - (8)].fval), 2.0*PI*(yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].ival), (yyvsp[(7) - (8)].ival)); ;}
    break;

  case 152:
#line 381 "p.y"
    { domain->solInfo().pade_pivot = true; domain->solInfo().pade_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 153:
#line 385 "p.y"
    { domain->solInfo().pade_poles = true; ;}
    break;

  case 154:
#line 387 "p.y"
    { domain->solInfo().pade_poles = true; 
          domain->solInfo().pade_poles_sigmaL = (yyvsp[(2) - (4)].fval); domain->solInfo().pade_poles_sigmaU = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 155:
#line 392 "p.y"
    { geoSource->setImpe((yyvsp[(2) - (3)].fval)); domain->addCoarseFrequency(2.0*PI*(yyvsp[(2) - (3)].fval)); ;}
    break;

  case 156:
#line 395 "p.y"
    { domain->addFrequencies(2.0*PI*(yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].ival)); ;}
    break;

  case 157:
#line 399 "p.y"
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
        ;}
    break;

  case 158:
#line 438 "p.y"
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
        ;}
    break;

  case 159:
#line 487 "p.y"
    { geoSource->binaryInput = true; geoSource->binaryOutput = true; ;}
    break;

  case 160:
#line 489 "p.y"
    { geoSource->binaryInput = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 161:
#line 491 "p.y"
    { geoSource->binaryOutput = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 163:
#line 496 "p.y"
    { geoSource->setGeo((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 164:
#line 498 "p.y"
    { geoSource->setDecomp((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 165:
#line 500 "p.y"
    { geoSource->setGlob((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 166:
#line 502 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 167:
#line 504 "p.y"
    { geoSource->setCpuMap((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 168:
#line 508 "p.y"
    { 
#ifdef STRUCTOPT	  
	  dynamic_cast<Domain_opt*>(domain)->addAnalysis((yyvsp[(2) - (3)].ival)); 
#endif
	;}
    break;

  case 169:
#line 516 "p.y"
    {if(decInit==0) decInit = new DecInit(); ;}
    break;

  case 170:
#line 518 "p.y"
    {decInit->file = strdup((yyvsp[(3) - (4)].strval));;}
    break;

  case 171:
#line 520 "p.y"
    {decInit->nsubs = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 172:
#line 522 "p.y"
    {decInit->weight = true; ;}
    break;

  case 173:
#line 524 "p.y"
    {decInit->memory = true; ;}
    break;

  case 174:
#line 526 "p.y"
    {decInit->exitAfterDec = true;;}
    break;

  case 175:
#line 528 "p.y"
    {decInit->skip = true;;}
    break;

  case 176:
#line 530 "p.y"
    {decInit->nosa = true; ;}
    break;

  case 177:
#line 534 "p.y"
    {;}
    break;

  case 178:
#line 536 "p.y"
    {
	   // map<int,double >::iterator it = weightList.find($2);
	   //if(it == weightList.end())
	     weightList[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval);
	 ;}
    break;

  case 179:
#line 544 "p.y"
    { (yyval.mftval) = new MFTTData; ;}
    break;

  case 180:
#line 546 "p.y"
    { (yyval.mftval) = (yyvsp[(1) - (4)].mftval); (yyval.mftval)->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 181:
#line 550 "p.y"
    { (yyval.mptval) = new MFTTData; ;}
    break;

  case 182:
#line 552 "p.y"
    { (yyval.mptval) = (yyvsp[(1) - (4)].mptval); (yyval.mptval)->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 183:
#line 556 "p.y"
    { (yyval.hftval) = new MFTTData; ;}
    break;

  case 184:
#line 558 "p.y"
    { (yyval.hftval) = (yyvsp[(1) - (4)].hftval); (yyval.hftval)->add((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 192:
#line 571 "p.y"
    { geoSource->addCFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 193:
#line 575 "p.y"
    { geoSource->addCoefInfo((yyvsp[(2) - (4)].ival)-1,(yyvsp[(4) - (4)].coefdata)); ;}
    break;

  case 194:
#line 579 "p.y"
    { (yyval.coefdata).zero(); (yyval.coefdata).setCoef((yyvsp[(1) - (4)].ival)-1,(yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 195:
#line 581 "p.y"
    { (yyval.coefdata).setCoef((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 196:
#line 585 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 197:
#line 587 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 198:
#line 591 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 199:
#line 593 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 200:
#line 597 "p.y"
    { (yyval.linfo) = new LayInfo(0); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 201:
#line 599 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 202:
#line 603 "p.y"
    { (yyval.linfo) = new LayInfo(1); geoSource->addLay((yyvsp[(2) - (3)].ival)-1,(yyval.linfo)); ;}
    break;

  case 203:
#line 605 "p.y"
    { (yyvsp[(1) - (2)].linfo)->add((yyvsp[(2) - (2)].ldata).lnum,(yyvsp[(2) - (2)].ldata).d,(yyvsp[(2) - (2)].ldata).matid); ;}
    break;

  case 204:
#line 609 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (11)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
	  (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval);
	  (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 205:
#line 615 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (13)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (13)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (13)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (13)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (13)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (13)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (13)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (13)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (13)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (13)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (13)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (13)].fval); ;}
    break;

  case 206:
#line 622 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (14)].ival)-1;
          (yyval.ldata).matid = -1; // PJSA 3-30-05: this means elastic constants are defined
          (yyval.ldata).d[0] = (yyvsp[(2) - (14)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (14)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (14)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (14)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (14)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (14)].fval);
          (yyval.ldata).d[6] = (yyvsp[(8) - (14)].fval); (yyval.ldata).d[7] = (yyvsp[(9) - (14)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (14)].fval);
          (yyval.ldata).d[9] = (yyvsp[(11) - (14)].fval);(yyval.ldata).d[10]= (yyvsp[(12) - (14)].fval); (yyval.ldata).d[11] = (yyvsp[(13) - (14)].fval); ;}
    break;

  case 207:
#line 631 "p.y"
    { (yyval.ldata).lnum = (yyvsp[(1) - (5)].ival)-1;  (yyval.ldata).matid = (yyvsp[(2) - (5)].ival)-1; (yyval.ldata).d[7] = (yyvsp[(3) - (5)].fval); (yyval.ldata).d[8] = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 209:
#line 636 "p.y"
    { geoSource->addLayMat((yyvsp[(2) - (2)].ldata).matid, (yyvsp[(2) - (2)].ldata).d); ;}
    break;

  case 210:
#line 643 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (7)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (7)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (7)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (7)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (7)].fval); (yyval.ldata).d[4] = 0.0; (yyval.ldata).d[5] = 0.0; (yyval.ldata).d[6] = (yyvsp[(6) - (7)].fval); 
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; ;}
    break;

  case 211:
#line 649 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (9)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (9)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (9)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (9)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (9)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (9)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (9)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (9)].fval);
          (yyval.ldata).d[7] = 0; (yyval.ldata).d[8] = 0; (yyval.ldata).d[9] = 0; ;}
    break;

  case 212:
#line 654 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (11)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (11)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (11)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (11)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (11)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (11)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (11)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (11)].fval);
          (yyval.ldata).d[7] = (yyvsp[(9) - (11)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (11)].fval); (yyval.ldata).d[9] = 0; ;}
    break;

  case 213:
#line 658 "p.y"
    { (yyval.ldata).matid = (yyvsp[(1) - (12)].ival)-1; (yyval.ldata).d[0] = (yyvsp[(2) - (12)].fval); (yyval.ldata).d[1] = (yyvsp[(3) - (12)].fval); (yyval.ldata).d[2] = (yyvsp[(4) - (12)].fval);
          (yyval.ldata).d[3] = (yyvsp[(5) - (12)].fval); (yyval.ldata).d[4] = (yyvsp[(6) - (12)].fval); (yyval.ldata).d[5] = (yyvsp[(7) - (12)].fval); (yyval.ldata).d[6] = (yyvsp[(8) - (12)].fval); 
          (yyval.ldata).d[7] = (yyvsp[(9) - (12)].fval); (yyval.ldata).d[8] = (yyvsp[(10) - (12)].fval); (yyval.ldata).d[9] = (yyvsp[(11) - (12)].fval); ;}
    break;

  case 215:
#line 665 "p.y"
    { domain->addDMass((yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1,(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 216:
#line 667 "p.y"
    { domain->addDMass((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(5) - (6)].fval),(yyvsp[(4) - (6)].ival)-1); ;}
    break;

  case 218:
#line 672 "p.y"
    { domain->setGravity((yyvsp[(2) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 220:
#line 677 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (5)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (5)].strval);
          geoSource->getCheckFileInfo()->FlagRST = (yyvsp[(4) - (5)].strval); ;}
    break;

  case 221:
#line 681 "p.y"
    { geoSource->getCheckFileInfo()->lastRestartFile = (yyvsp[(2) - (4)].strval);
          geoSource->getCheckFileInfo()->outputExt = (yyvsp[(3) - (4)].strval);;}
    break;

  case 222:
#line 684 "p.y"
    { geoSource->getCheckFileInfo()->currentRestartFile = (yyvsp[(2) - (4)].strval);
          domain->solInfo().nRestart = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 223:
#line 689 "p.y"
    { geoSource->setControlFile((yyvsp[(2) - (3)].strval));
         geoSource->setControlRoutine((char *) "controlObj");;}
    break;

  case 224:
#line 694 "p.y"
    { geoSource->setControlRoutine((yyvsp[(2) - (3)].strval)); ;}
    break;

  case 225:
#line 698 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Sensors;
          if(geoSource->setSensorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 226:
#line 703 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Actuators; (yyvsp[(3) - (3)].bclist)->d[i].caseid = 0; }
          if(geoSource->setActuatorLocations((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; 
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)            return -1; ;}
    break;

  case 227:
#line 709 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdf; (yyvsp[(3) - (3)].bclist)->d[i].caseid = 0; }
          if(geoSource->setUsdfLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setNeuman((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)       return -1; ;}
    break;

  case 228:
#line 716 "p.y"
    { geoSource->binaryInputControlLeft = true;
          for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Usdd;
          if(geoSource->setUsddLocation((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1;
          if(geoSource->setDirichlet((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0)    return -1; ;}
    break;

  case 229:
#line 723 "p.y"
    { numColumns = 3; ;}
    break;

  case 230:
#line 725 "p.y"
    { numColumns = 6; ;}
    break;

  case 231:
#line 727 "p.y"
    { numColumns = 3; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 232:
#line 729 "p.y"
    { numColumns = 6; geoSource->setOutLimit((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 233:
#line 731 "p.y"
    { (yyvsp[(2) - (3)].oinfo).finalize(numColumns); geoSource->addOutput((yyvsp[(2) - (3)].oinfo)); ;}
    break;

  case 234:
#line 735 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 235:
#line 737 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); ;}
    break;

  case 236:
#line 739 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; ;}
    break;

  case 237:
#line 741 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); 
          if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1;;}
    break;

  case 238:
#line 744 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; ;}
    break;

  case 239:
#line 746 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = (OutputInfo::Type) (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; ;}
    break;

  case 240:
#line 748 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (3)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (3)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 241:
#line 750 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).width = (yyvsp[(2) - (5)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (5)].ival); ;}
    break;

  case 242:
#line 752 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (4)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (4)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (4)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(4) - (4)].ival)-1; ;}
    break;

  case 243:
#line 754 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (5)].ival); (yyval.oinfo).filename = (yyvsp[(2) - (5)].strval); (yyval.oinfo).interval = (yyvsp[(3) - (5)].ival); if ((yyvsp[(4) - (5)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(5) - (5)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(5) - (5)].ival)-1; ;}
    break;

  case 244:
#line 756 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (6)].ival); (yyval.oinfo).width = (yyvsp[(2) - (6)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (6)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (6)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (6)].ival); (yyval.oinfo).nodeNumber = (yyvsp[(6) - (6)].ival)-1; ;}
    break;

  case 245:
#line 758 "p.y"
    { (yyval.oinfo).initialize(); (yyval.oinfo).type = OutputInfo::TDEnforcement; (yyval.oinfo).tdenforc_var = (yyvsp[(1) - (7)].ival); (yyval.oinfo).width = (yyvsp[(2) - (7)].ival); (yyval.oinfo).precision = (yyvsp[(3) - (7)].ival); (yyval.oinfo).filename = (yyvsp[(4) - (7)].strval); (yyval.oinfo).interval = (yyvsp[(5) - (7)].ival); if ((yyvsp[(6) - (7)].ival) == OutputInfo::NodeGroup) (yyval.oinfo).groupNumber = (yyvsp[(7) - (7)].ival); else (yyval.oinfo).nodeNumber = (yyvsp[(7) - (7)].ival)-1; ;}
    break;

  case 246:
#line 761 "p.y"
    { (yyval.oinfo).nodeNumber = (yyvsp[(3) - (3)].ival)-1; ;}
    break;

  case 247:
#line 763 "p.y"
    { (yyval.oinfo).surface = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 248:
#line 765 "p.y"
    { (yyval.oinfo).ylayer = (yyvsp[(2) - (3)].fval); (yyval.oinfo).zlayer = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 249:
#line 767 "p.y"
    { (yyval.oinfo).averageFlg = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 250:
#line 769 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 251:
#line 771 "p.y"
    { (yyval.oinfo).complexouttype = (yyvsp[(2) - (3)].ival); (yyval.oinfo).ncomplexout = (yyvsp[(3) - (3)].ival); ;}
    break;

  case 252:
#line 773 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 253:
#line 775 "p.y"
    { (yyval.oinfo).ndtype = (yyvsp[(2) - (3)].ival); sfem->setnsamp_out((yyvsp[(3) - (3)].ival)); ;}
    break;

  case 254:
#line 777 "p.y"
    { /* TODO: $$.oframe = OutputInfo::Global; */ ;}
    break;

  case 255:
#line 779 "p.y"
    { /* TODO: $$.oframe = OutputInfo::Local; */ ;}
    break;

  case 256:
#line 781 "p.y"
    { (yyval.oinfo).matlab = true; ;}
    break;

  case 257:
#line 785 "p.y"
    { domain->outFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 259:
#line 790 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
          domain->solInfo().eigenSolverType = SolverInfo::SubSpace; ;}
    break;

  case 260:
#line 793 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Modal);
	  domain->solInfo().nEig = (yyvsp[(2) - (3)].ival);;}
    break;

  case 261:
#line 796 "p.y"
    { domain->solInfo().nEig = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 262:
#line 798 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::SubSpace;;}
    break;

  case 263:
#line 800 "p.y"
    { domain->solInfo().setSubSpaceInfo((yyvsp[(2) - (5)].ival),(yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); ;}
    break;

  case 264:
#line 802 "p.y"
    { domain->solInfo().subspaceSize = (yyvsp[(2) - (3)].ival);;}
    break;

  case 265:
#line 804 "p.y"
    { domain->solInfo().tolEig = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 266:
#line 806 "p.y"
    { domain->solInfo().tolJac = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 267:
#line 808 "p.y"
    { domain->solInfo().explicitK = true; ;}
    break;

  case 268:
#line 810 "p.y"
    { geoSource->setShift((yyvsp[(2) - (3)].fval)); ;}
    break;

  case 269:
#line 812 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack; ;}
    break;

  case 270:
#line 814 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (3)].strval); ;}
    break;

  case 271:
#line 817 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->solInfo().which = (yyvsp[(2) - (4)].strval); 
          domain->solInfo().arpack_mode = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 272:
#line 821 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValue((yyvsp[(2) - (4)].fval), int((yyvsp[(3) - (4)].fval))); ;}
    break;

  case 273:
#line 824 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::Arpack;
          domain->setEigenValues((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival));;}
    break;

  case 274:
#line 827 "p.y"
    { domain->solInfo().filtereig = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 275:
#line 829 "p.y"
    { domain->solInfo().eigenSolverSubType = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 276:
#line 831 "p.y"
    { domain->solInfo().eigenSolverType = SolverInfo::LobPcg;
          domain->solInfo().explicitK = true;;}
    break;

  case 277:
#line 834 "p.y"
    { domain->solInfo().maxitEig = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 278:
#line 836 "p.y"
    { domain->solInfo().test_ulrich = true; ;}
    break;

  case 279:
#line 838 "p.y"
    { domain->solInfo().addedMass = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 280:
#line 842 "p.y"
    { domain->solInfo().sloshing = 1; ;}
    break;

  case 281:
#line 844 "p.y"
    { domain->setGravitySloshing((yyvsp[(2) - (3)].fval)); ;}
    break;

  case 282:
#line 848 "p.y"
    { domain->solInfo().massFlag = 1; ;}
    break;

  case 283:
#line 852 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber); 
	  domain->solInfo().setCondNumTol((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 284:
#line 855 "p.y"
    { domain->solInfo().setProbType(SolverInfo::ConditionNumber);;}
    break;

  case 285:
#line 859 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Top); ;}
    break;

  case 286:
#line 863 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); ;}
    break;

  case 290:
#line 868 "p.y"
    { domain->solInfo().modal = true; ;}
    break;

  case 291:
#line 870 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 292:
#line 872 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 293:
#line 874 "p.y"
    { domain->solInfo().stable = (yyvsp[(3) - (8)].ival);
          domain->solInfo().stable_cfl = (yyvsp[(4) - (8)].fval);
          domain->solInfo().stable_tol = (yyvsp[(5) - (8)].fval);
          domain->solInfo().stable_maxit = (yyvsp[(6) - (8)].ival);
          domain->solInfo().stable_freq = (yyvsp[(7) - (8)].ival);
        ;}
    break;

  case 294:
#line 881 "p.y"
    { domain->solInfo().iacc_switch = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 295:
#line 883 "p.y"
    { domain->solInfo().zeroRot = bool((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 296:
#line 885 "p.y"
    { domain->solInfo().no_secondary = true; ;}
    break;

  case 297:
#line 887 "p.y"
    { domain->solInfo().check_energy_balance = true; ;}
    break;

  case 298:
#line 889 "p.y"
    { domain->solInfo().check_energy_balance = true;
          domain->solInfo().epsilon1 = (yyvsp[(3) - (5)].fval); 
          domain->solInfo().epsilon2 = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 299:
#line 895 "p.y"
    { domain->solInfo().timeIntegration = SolverInfo::Newmark; ;}
    break;

  case 301:
#line 898 "p.y"
    { domain->solInfo().acoustic = true; ;}
    break;

  case 303:
#line 903 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (2)].fval),(yyvsp[(2) - (2)].fval)); ;}
    break;

  case 304:
#line 905 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (4)].fval),(yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval),(yyvsp[(4) - (4)].fval)); ;}
    break;

  case 305:
#line 907 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo(0.0,0.0,10.0,10.0,(yyvsp[(1) - (1)].fval)); ;}
    break;

  case 306:
#line 909 "p.y"
    { domain->solInfo().setNewmarkSecondOrderInfo((yyvsp[(1) - (3)].fval),(yyvsp[(2) - (3)].fval));
          domain->solInfo().modifiedWaveEquation = true;
          domain->solInfo().modifiedWaveEquationCoef = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 307:
#line 915 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinDynam) {
          domain->solInfo().order = 1;
         }
          else 
          domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setNewmarkFirstOrderInfo((yyvsp[(1) - (1)].fval)); 
        ;}
    break;

  case 308:
#line 926 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Dynamic); 
          domain->solInfo().timeIntegration = SolverInfo::Qstatic; ;}
    break;

  case 311:
#line 931 "p.y"
    { domain->solInfo().modal = true; ;}
    break;

  case 312:
#line 935 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (5)].fval), 0, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival)); ;}
    break;

  case 313:
#line 937 "p.y"
    { domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), 0, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].ival), (yyvsp[(5) - (6)].fval)); ;}
    break;

  case 314:
#line 941 "p.y"
    { domain->solInfo().setProbType(SolverInfo::TempDynamic);
          domain->solInfo().setQuasistaticInfo((yyvsp[(2) - (6)].fval), (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].ival)); ;}
    break;

  case 315:
#line 951 "p.y"
    { domain->solInfo().setAero((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCollocated = 0; ;}
    break;

  case 316:
#line 954 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().isCollocated = 0; ;}
    break;

  case 317:
#line 957 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (6)].ival));
          domain->solInfo().isCollocated = 0;
          if((yyvsp[(3) - (6)].ival) < 6 || (yyvsp[(3) - (6)].ival) == 20) {
              domain->solInfo().alphas[0] = (yyvsp[(4) - (6)].fval)+(yyvsp[(5) - (6)].fval);
              domain->solInfo().alphas[1] = -(yyvsp[(5) - (6)].fval);
          }
        ;}
    break;

  case 318:
#line 965 "p.y"
    { domain->solInfo().setAero((yyvsp[(3) - (5)].ival));
          domain->solInfo().isCollocated = 0;
          domain->solInfo().mppFactor = (yyvsp[(4) - (5)].fval);
        ;}
    break;

  case 319:
#line 970 "p.y"
    { domain->solInfo().isCollocated = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 320:
#line 972 "p.y"
    { geoSource->setMatch((yyvsp[(3) - (4)].strval)); ;}
    break;

  case 321:
#line 974 "p.y"
    {;}
    break;

  case 322:
#line 978 "p.y"
    {;}
    break;

  case 323:
#line 980 "p.y"
    { domain->AddAeroEmbedSurfaceId((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 324:
#line 984 "p.y"
    { domain->solInfo().setAeroHeat((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); ;}
    break;

  case 325:
#line 988 "p.y"
    { domain->solInfo().setThermoh(1); ;}
    break;

  case 326:
#line 992 "p.y"
    { domain->solInfo().setThermoe(1); ;}
    break;

  case 327:
#line 996 "p.y"
    { domain->solInfo().setModeDecomp(1); ;}
    break;

  case 328:
#line 1000 "p.y"
    { domain->solInfo().hzemFlag=1; ;}
    break;

  case 329:
#line 1004 "p.y"
    { domain->solInfo().slzemFlag=1; ;}
    break;

  case 330:
#line 1008 "p.y"
    { domain->solInfo().setTrbm((yyvsp[(3) - (4)].fval)); ;}
    break;

  case 331:
#line 1016 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (5)].fval),(yyvsp[(4) - (5)].fval)); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");;}
    break;

  case 332:
#line 1019 "p.y"
    { domain->solInfo().setGrbm((yyvsp[(3) - (4)].fval)); 
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");;}
    break;

  case 333:
#line 1022 "p.y"
    { domain->solInfo().setGrbm();
         filePrint(stderr," ... Using Geometric RBM Method     ...\n");;}
    break;

  case 334:
#line 1027 "p.y"
    { domain->solInfo().modeFilterFlag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 335:
#line 1029 "p.y"
    { domain->solInfo().modeFilterFlag = 1; ;}
    break;

  case 336:
#line 1033 "p.y"
    { domain->solInfo().useRbmFilter((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 337:
#line 1035 "p.y"
    { domain->solInfo().useRbmFilter(1); ;}
    break;

  case 339:
#line 1040 "p.y"
    { if((yyvsp[(1) - (1)].ival) < 1 || (yyvsp[(1) - (1)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(1) - (1)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(1) - (1)].ival)-1] = 1;
    ;}
    break;

  case 340:
#line 1048 "p.y"
    { if((yyvsp[(2) - (2)].ival) < 1 || (yyvsp[(2) - (2)].ival) > 6){
        fprintf(stderr, " *** ERROR: RBMF specifier must be in the range 1-6, found: %d\n", (yyvsp[(2) - (2)].ival));
        yyerror(NULL);
        exit(-1);
      }
      domain->solInfo().rbmFilters[(yyvsp[(2) - (2)].ival)-1] = 1;
    ;}
    break;

  case 341:
#line 1058 "p.y"
    { domain->solInfo().hzemFilterFlag=1; ;}
    break;

  case 342:
#line 1062 "p.y"
    { domain->solInfo().slzemFilterFlag=1; ;}
    break;

  case 343:
#line 1066 "p.y"
    { domain->solInfo().setTimes((yyvsp[(4) - (5)].fval),(yyvsp[(3) - (5)].fval),(yyvsp[(2) - (5)].fval)); ;}
    break;

  case 344:
#line 1070 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (5)].ival),(yyvsp[(4) - (5)].ival),1);
        ;}
    break;

  case 345:
#line 1076 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].ival));
        ;}
    break;

  case 346:
#line 1081 "p.y"
    {
          domain->solInfo().activatePita = true;
          domain->solInfo().mdPita = true;
          domain->solInfo().setParallelInTime((yyvsp[(3) - (7)].ival),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].ival)); 
          /*domain->solInfo().numSpaceMPIProc = $6;*/
        ;}
    break;

  case 349:
#line 1094 "p.y"
    { domain->solInfo().pitaNoForce = true; ;}
    break;

  case 350:
#line 1096 "p.y"
    { domain->solInfo().pitaGlobalBasisImprovement = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 351:
#line 1098 "p.y"
    { domain->solInfo().pitaLocalBasisImprovement = 1; ;}
    break;

  case 352:
#line 1100 "p.y"
    { domain->solInfo().pitaTimeReversible = true; ;}
    break;

  case 353:
#line 1102 "p.y"
    { domain->solInfo().pitaRemoteCoarse = true; ;}
    break;

  case 354:
#line 1104 "p.y"
    { domain->solInfo().pitaProjTol = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 355:
#line 1106 "p.y"
    { domain->solInfo().pitaReadInitSeed = true; ;}
    break;

  case 356:
#line 1108 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = 0.0; ;}
    break;

  case 357:
#line 1110 "p.y"
    { domain->solInfo().pitaJumpCvgRatio = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 358:
#line 1112 "p.y"
    { domain->solInfo().pitaJumpMagnOutput = true; ;}
    break;

  case 359:
#line 1116 "p.y"
    { domain->solInfo().setDamping((yyvsp[(2) - (4)].fval),(yyvsp[(3) - (4)].fval)); ;}
    break;

  case 360:
#line 1118 "p.y"
    { if(geoSource->setModalDamping((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 361:
#line 1123 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); ;}
    break;

  case 362:
#line 1127 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 363:
#line 1133 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval), 0.0);
        ;}
    break;

  case 364:
#line 1141 "p.y"
    {
           domain->implicitFlag = 1;
           domain->solInfo().setProbType(SolverInfo::HelmholtzDirSweep);
           domain->setWaveDirections((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 365:
#line 1147 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections((yyvsp[(2) - (3)].ival),0.0,0.0,0.0);
        ;}
    break;

  case 367:
#line 1155 "p.y"
    {
           domain->setWaveDirections((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        ;}
    break;

  case 368:
#line 1161 "p.y"
    {
           domain->implicitFlag = 1;
           domain->setWaveDirections(1,0.0,0.0,0.0);
           domain->setWaveDirections(1, (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 369:
#line 1169 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 370:
#line 1171 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Displacements; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 371:
#line 1173 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } ;}
    break;

  case 372:
#line 1175 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Displacements); (yyval.bclist)->add(bc); } ;}
    break;

  case 373:
#line 1177 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Displacements;
          geoSource->addSurfaceDirichlet(1,surf_bc); ;}
    break;

  case 375:
#line 1192 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(2) - (6)].ival)-1;
          surf_bc[0].type = (BCond::BCType) (yyvsp[(3) - (6)].ival); //BCond::PointPlaneDistance;
          surf_bc[0].dofnum = (yyvsp[(4) - (6)].ival)-1;
          surf_bc[0].val = (yyvsp[(5) - (6)].ival)-1;
          geoSource->addSurfaceConstraint(1,surf_bc);
        ;}
    break;

  case 376:
#line 1201 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Pdir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 377:
#line 1205 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 378:
#line 1207 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 379:
#line 1211 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 380:
#line 1213 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (2)].ival)-1; (yyval.bcval).dofnum = 10; (yyval.bcval).val = 0.0; ;}
    break;

  case 381:
#line 1217 "p.y"
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

  case 382:
#line 1258 "p.y"
    { (yyval.bclist) = new BCList; 
          for (int ii = 0; ii < (yyvsp[(1) - (1)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(1) - (1)].bclist)->d)[ii]); ;}
    break;

  case 383:
#line 1262 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); 
          for (int ii = 0; ii < (yyvsp[(2) - (2)].bclist)->n; ii++) 
           (yyval.bclist)->add(((yyvsp[(2) - (2)].bclist)->d)[ii]); ;}
    break;

  case 384:
#line 1268 "p.y"
    { (yyval.bclist) = new BCList;
          //BCond *bclist = new BCond[$3.num];
          for(int i=0; i<(yyvsp[(3) - (4)].nl).num; ++i) 
          { (yyval.bclist)->add((yyvsp[(3) - (4)].nl).nd[i],10,0.0); } ;}
    break;

  case 385:
#line 1276 "p.y"
    { (yyval.bclist) = new BCList; if(domain->solInfo().soltyp != 2) domain->solInfo().thermalLoadFlag = 1;;}
    break;

  case 386:
#line 1278 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Temperatures; (yyval.bclist)->add(bc); ;}
    break;

  case 387:
#line 1281 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].dofnum = 6;
          surf_bc[0].type = BCond::Temperatures;
          geoSource->addSurfaceDirichlet(1,surf_bc); ;}
    break;

  case 388:
#line 1290 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 389:
#line 1292 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 390:
#line 1294 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (4)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (4)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (4)].fval); bc.type = BCond::Flux; bc.caseid = (yyval.bclist)->caseid; (yyval.bclist)->add(bc); ;}
    break;

  case 391:
#line 1297 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0].nnum = (yyvsp[(3) - (5)].ival)-1;
          surf_bc[0].dofnum = 6;
          surf_bc[0].val = (yyvsp[(4) - (5)].fval);
          surf_bc[0].type = BCond::Flux;
          surf_bc[0].caseid = (yyval.bclist)->caseid;
          geoSource->addSurfaceNeuman(1,surf_bc); ;}
    break;

  case 392:
#line 1307 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 393:
#line 1309 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 394:
#line 1311 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = (yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Convection; bc.caseid = (yyval.bclist)->caseid; (yyval.bclist)->add(bc); ;}
    break;

  case 395:
#line 1316 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 396:
#line 1318 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 397:
#line 1320 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (6)].bclist); BCond bc; bc.nnum = (yyvsp[(2) - (6)].ival)-1; bc.dofnum = 6;
          bc.val = 5.670400E-8*(yyvsp[(3) - (6)].fval)*(yyvsp[(4) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval)*(yyvsp[(5) - (6)].fval); bc.type = BCond::Radiation; bc.caseid = (yyval.bclist)->caseid; (yyval.bclist)->add(bc); ;}
    break;

  case 401:
#line 1332 "p.y"
    { domain->addSommer(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); ;}
    break;

  case 402:
#line 1334 "p.y"
    { domain->addSommer(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); ;}
    break;

  case 403:
#line 1336 "p.y"
    { domain->addSommer(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); ;}
    break;

  case 406:
#line 1344 "p.y"
    { domain->addSommerElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          /*geoSource->addElem($1-1, $2, $3.num, $3.nd);include Sommer nodes in PackedEset -JF*/
        ;}
    break;

  case 407:
#line 1350 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1; ;}
    break;

  case 408:
#line 1352 "p.y"
    { if((yyval.nl).num == 64) return -1;
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++; ;}
    break;

  case 412:
#line 1364 "p.y"
    { domain->addScatter(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1));
          domain->addNeum(new LineSommerBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1)); ;}
    break;

  case 413:
#line 1367 "p.y"
    { domain->addScatter(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1));
          domain->addNeum(new TriangleSommerBC((yyvsp[(1) - (5)].ival)-1,(yyvsp[(2) - (5)].ival)-1,(yyvsp[(3) - (5)].ival)-1)); ;}
    break;

  case 414:
#line 1370 "p.y"
    { domain->addScatter(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1));
          domain->addNeum(new QuadSommerBC((yyvsp[(1) - (6)].ival)-1,(yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1, (yyvsp[(4) - (6)].ival)-1)); ;}
    break;

  case 417:
#line 1379 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 420:
#line 1388 "p.y"
    { domain->addNeumElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 423:
#line 1396 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); 
          domain->solInfo().isCoupled = true; 
          domain->solInfo().isMatching = true; ;}
    break;

  case 426:
#line 1407 "p.y"
    { domain->addScatterElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);;}
    break;

  case 427:
#line 1411 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 7; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 428:
#line 1415 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 429:
#line 1417 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 430:
#line 1421 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atddir; (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 431:
#line 1425 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) { (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Atdneu; (yyvsp[(3) - (3)].bclist)->d[i].caseid = 0; } (yyval.bclist) = (yyvsp[(3) - (3)].bclist); ;}
    break;

  case 432:
#line 1429 "p.y"
    { domain->solInfo().ATDARBFlag = (yyvsp[(2) - (3)].fval);;}
    break;

  case 434:
#line 1434 "p.y"
    { domain->solInfo().ATDDNBVal = (yyvsp[(2) - (3)].fval);;}
    break;

  case 436:
#line 1439 "p.y"
    { domain->solInfo().ATDROBVal = (yyvsp[(2) - (5)].fval);
          domain->solInfo().ATDROBalpha = (yyvsp[(3) - (5)].fval);
          domain->solInfo().ATDROBbeta = (yyvsp[(4) - (5)].fval);;}
    break;

  case 438:
#line 1446 "p.y"
    { domain->setFFP((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 439:
#line 1448 "p.y"
    { domain->setFFP((yyvsp[(2) - (4)].ival),(yyvsp[(3) - (4)].ival)); ;}
    break;

  case 440:
#line 1452 "p.y"
    {
           domain->setFFP((yyvsp[(2) - (4)].ival));
        ;}
    break;

  case 441:
#line 1458 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        ;}
    break;

  case 442:
#line 1463 "p.y"
    { fourHelmBC->addDirichlet((yyvsp[(2) - (2)].complexFDBC)); ;}
    break;

  case 443:
#line 1467 "p.y"
    { (yyval.complexFDBC) = FDBC((yyvsp[(1) - (2)].ival)-1); ;}
    break;

  case 444:
#line 1471 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setConst(DComplex((yyvsp[(3) - (9)].fval),(yyvsp[(4) - (9)].fval)));
          fourHelmBC->setDir((yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval));
        ;}
    break;

  case 445:
#line 1476 "p.y"
    { fourHelmBC->addNeuman((yyvsp[(2) - (2)].complexFNBC)); ;}
    break;

  case 446:
#line 1480 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1); ;}
    break;

  case 447:
#line 1482 "p.y"
    { (yyval.complexFNBC) = FNBC((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 448:
#line 1486 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setModes((yyvsp[(2) - (3)].ival));
          domain->solInfo().setProbType(SolverInfo::AxiHelm);
        ;}
    break;

  case 449:
#line 1494 "p.y"
    {
          if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSlices((yyvsp[(2) - (3)].ival));
        ;}
    break;

  case 450:
#line 1501 "p.y"
    { if(fourHelmBC == 0) fourHelmBC = new FourierHelmBCs();
          fourHelmBC->setSomType((yyvsp[(3) - (7)].ival));
          fourHelmBC->setSurf((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval));
        ;}
    break;

  case 452:
#line 1509 "p.y"
    { fourHelmBC->addSommer(new LineAxiSommer((yyvsp[(1) - (3)].ival)-1, (yyvsp[(2) - (3)].ival)-1)); ;}
    break;

  case 453:
#line 1511 "p.y"
    { fourHelmBC->addSommer(new Line2AxiSommer((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1)); ;}
    break;

  case 454:
#line 1515 "p.y"
    { if( globalMPCs== NULL) globalMPCs = new MPCData(); ;}
    break;

  case 455:
#line 1517 "p.y"
    { globalMPCs->addMPC((yyvsp[(2) - (2)].axiMPC)); ;}
    break;

  case 456:
#line 1521 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (5)].ival)-1, (yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].ival), DComplex(1.0,0.0), 0.0, 0.0, 0.0); ;}
    break;

  case 457:
#line 1523 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (7)].ival)-1, (yyvsp[(2) - (7)].fval), (yyvsp[(3) - (7)].fval), (yyvsp[(4) - (7)].ival), DComplex((yyvsp[(5) - (7)].fval),(yyvsp[(6) - (7)].fval)) , 0.0, 0.0, 0.0); ;}
    break;

  case 458:
#line 1525 "p.y"
    { (yyval.axiMPC) = MPC((yyvsp[(1) - (10)].ival)-1, (yyvsp[(2) - (10)].fval), (yyvsp[(3) - (10)].fval), (yyvsp[(4) - (10)].ival), DComplex((yyvsp[(5) - (10)].fval),(yyvsp[(6) - (10)].fval)) , (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)); ;}
    break;

  case 459:
#line 1529 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (3)].strval);
	  domain->solInfo().readmodeCalled = true; ;}
    break;

  case 460:
#line 1532 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (4)].strval);
          domain->solInfo().readmodeCalled = true; 
 	  domain->solInfo().maxSizePodRom = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 461:
#line 1538 "p.y"
    { ;}
    break;

  case 462:
#line 1540 "p.y"
    { domain->solInfo().zeroInitialDisp = 1; ;}
    break;

  case 463:
#line 1542 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 464:
#line 1545 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Idisplacements;
          if(geoSource->setIDisModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 465:
#line 1551 "p.y"
    { (yyval.bclist) = new BCList; amplitude = (yyvsp[(2) - (3)].fval);  ;}
    break;

  case 466:
#line 1553 "p.y"
    { (yyval.bclist) = new BCList; amplitude = 1.0; ;}
    break;

  case 467:
#line 1555 "p.y"
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

  case 468:
#line 1566 "p.y"
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

  case 469:
#line 1577 "p.y"
    { fprintf(stderr," ... Geometric Pre-Stress Effects   ... \n"); 
          domain->solInfo().setGEPS(); ;}
    break;

  case 470:
#line 1580 "p.y"
    { domain->solInfo().buckling = 1; ;}
    break;

  case 471:
#line 1585 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 472:
#line 1587 "p.y"
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

  case 473:
#line 1600 "p.y"
    { (yyval.bclist) = new BCList; PitaTS = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 474:
#line 1602 "p.y"
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

  case 475:
#line 1614 "p.y"
    { ;}
    break;

  case 476:
#line 1616 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVel((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 477:
#line 1619 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].bclist)->n; ++i) (yyvsp[(4) - (4)].bclist)->d[i].type = BCond::Ivelocities;
          if(geoSource->setIVelModal((yyvsp[(4) - (4)].bclist)->n, (yyvsp[(4) - (4)].bclist)->d) < 0) return -1; 
	  domain->solInfo().modalCalled = true; ;}
    break;

  case 478:
#line 1625 "p.y"
    { for(int i=0; i<(yyvsp[(3) - (3)].bclist)->n; ++i) (yyvsp[(3) - (3)].bclist)->d[i].type = BCond::Itemperatures;
          if(geoSource->setIDis((yyvsp[(3) - (3)].bclist)->n,(yyvsp[(3) - (3)].bclist)->d) < 0) return -1; ;}
    break;

  case 479:
#line 1630 "p.y"
    { (yyval.bclist) = new BCList; ;}
    break;

  case 480:
#line 1632 "p.y"
    { (yyval.bclist) = new BCList((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 481:
#line 1634 "p.y"
    { (yyvsp[(2) - (2)].bcval).type = BCond::Forces; (yyvsp[(2) - (2)].bcval).caseid = (yyval.bclist)->caseid; (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 482:
#line 1636 "p.y"
    { for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval), BCond::Forces, (yyval.bclist)->caseid); (yyval.bclist)->add(bc); } ;}
    break;

  case 483:
#line 1638 "p.y"
    { for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval), BCond::Forces, (yyval.bclist)->caseid); (yyval.bclist)->add(bc); } ;}
    break;

  case 484:
#line 1640 "p.y"
    { BCond *surf_bc = new BCond[1];
          surf_bc[0] = (yyvsp[(3) - (3)].bcval);
          surf_bc[0].type = BCond::Forces;
          surf_bc[0].caseid = (yyval.bclist)->caseid;
          geoSource->addSurfaceNeuman(1,surf_bc); ;}
    break;

  case 485:
#line 1648 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 486:
#line 1650 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 487:
#line 1652 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (6)].ival); i<=(yyvsp[(3) - (6)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(4) - (6)].ival)-1, (yyvsp[(5) - (6)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 488:
#line 1654 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (7)].bclist); for(int i=(yyvsp[(2) - (7)].ival); i<=(yyvsp[(4) - (7)].ival); ++i) { BCond bc; bc.setData(i-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 489:
#line 1656 "p.y"
    { (yyval.bclist) = new BCList; for(int i=(yyvsp[(1) - (8)].ival); i<=(yyvsp[(3) - (8)].ival); i+=(yyvsp[(5) - (8)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(6) - (8)].ival)-1, (yyvsp[(7) - (8)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 490:
#line 1658 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (9)].bclist); for(int i=(yyvsp[(2) - (9)].ival); i<=(yyvsp[(4) - (9)].ival); i+=(yyvsp[(6) - (9)].ival)) { BCond bc; bc.setData(i-1, (yyvsp[(7) - (9)].ival)-1, (yyvsp[(8) - (9)].fval)); (yyval.bclist)->add(bc); } ;}
    break;

  case 491:
#line 1662 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 492:
#line 1664 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 493:
#line 1668 "p.y"
    { (yyval.bclist) = new BCList; (yyval.bclist)->add((yyvsp[(1) - (1)].bcval)); ;}
    break;

  case 494:
#line 1670 "p.y"
    { (yyval.bclist) = (yyvsp[(1) - (2)].bclist); (yyval.bclist)->add((yyvsp[(2) - (2)].bcval)); ;}
    break;

  case 497:
#line 1678 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ymtt)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addYMTT((yyval.ymtt));;}
    break;

  case 498:
#line 1680 "p.y"
    { (yyval.ymtt)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 499:
#line 1682 "p.y"
    { (yyval.ymtt) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ymtt)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addYMTT((yyval.ymtt));;}
    break;

  case 502:
#line 1690 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(2) - (6)].ival)); (yyval.ctett)->add((yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval)); domain->addCTETT((yyval.ctett));;}
    break;

  case 503:
#line 1692 "p.y"
    { (yyval.ctett)->add((yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 504:
#line 1694 "p.y"
    { (yyval.ctett) = new MFTTData((yyvsp[(3) - (7)].ival)); (yyval.ctett)->add((yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->addCTETT((yyval.ctett));;}
    break;

  case 507:
#line 1702 "p.y"
    { (yyval.lmpcons) = (yyvsp[(1) - (2)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 508:
#line 1706 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); ;}
    break;

  case 509:
#line 1708 "p.y"
    { (yyval.lmpcons) = (yyvsp[(2) - (3)].lmpcons);
          (yyval.lmpcons)->addterm((yyvsp[(3) - (3)].mpcterm));
          domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 510:
#line 1714 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].ival), 0.0); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 511:
#line 1717 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (3)].ival), (yyvsp[(2) - (3)].fval)); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 512:
#line 1720 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (5)].ival), (yyvsp[(2) - (5)].fval));
          (yyval.lmpcons)->type = (yyvsp[(4) - (5)].ival); 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 513:
#line 1724 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (4)].ival), (yyvsp[(2) - (4)].fval));
          (yyval.lmpcons)->lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          (yyval.lmpcons)->penalty = (yyvsp[(3) - (4)].copt).penalty; 
          (yyval.lmpcons)->setSource(mpc::Lmpc); ;}
    break;

  case 514:
#line 1731 "p.y"
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

  case 517:
#line 1747 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(1) - (2)].cxbcval).nnum,(yyvsp[(1) - (2)].cxbcval).reval,(yyvsp[(1) - (2)].cxbcval).imval,(yyvsp[(2) - (2)].mpcterm)); domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 518:
#line 1749 "p.y"
    { (yyval.lmpcons)->addterm((yyvsp[(2) - (2)].mpcterm)); ;}
    break;

  case 519:
#line 1751 "p.y"
    { (yyval.lmpcons) = new LMPCons((yyvsp[(2) - (3)].cxbcval).nnum,(yyvsp[(2) - (3)].cxbcval).reval,(yyvsp[(2) - (3)].cxbcval).imval,(yyvsp[(3) - (3)].mpcterm)); domain->addLMPC((yyval.lmpcons)); ;}
    break;

  case 520:
#line 1755 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (5)].ival); (yyval.cxbcval).reval=(yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval=(yyvsp[(4) - (5)].fval); ;}
    break;

  case 521:
#line 1757 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (3)].ival); (yyval.cxbcval).reval=(yyvsp[(2) - (3)].fval); (yyval.cxbcval).imval=0.0; ;}
    break;

  case 522:
#line 1759 "p.y"
    { (yyval.cxbcval).nnum=(yyvsp[(1) - (2)].ival); (yyval.cxbcval).reval=0.0; (yyval.cxbcval).imval=0.0; ;}
    break;

  case 523:
#line 1763 "p.y"
    { if(((yyvsp[(3) - (5)].fval)==0.0) && ((yyvsp[(4) - (5)].fval)==0.0)) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (5)].ival),(yyvsp[(2) - (5)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (5)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (5)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval)); }
        ;}
    break;

  case 524:
#line 1771 "p.y"
    { if((yyvsp[(3) - (4)].fval)==0.0) {
          fprintf(stderr," *** ERROR: zero coefficient in LMPC\n");
          fprintf(stderr," ***          node %d dof %d\n",(yyvsp[(1) - (4)].ival),(yyvsp[(2) - (4)].ival));
          return -1;
          }
          else { (yyval.mpcterm) = new LMPCTerm(true); (yyval.mpcterm)->nnum=((yyvsp[(1) - (4)].ival)-1); (yyval.mpcterm)->dofnum=((yyvsp[(2) - (4)].ival)-1); (yyval.mpcterm)->coef.c_value=DComplex((yyvsp[(3) - (4)].fval),0.0); }
        ;}
    break;

  case 525:
#line 1781 "p.y"
    { (yyval.cxbclist) = (yyvsp[(3) - (3)].cxbclist); ;}
    break;

  case 526:
#line 1783 "p.y"
    { for(int i=0; i<(yyvsp[(4) - (4)].cxbclist)->n; ++i) (yyvsp[(4) - (4)].cxbclist)->d[i].caseid = (yyvsp[(2) - (4)].ival);
          (yyval.cxbclist) = (yyvsp[(4) - (4)].cxbclist); ;}
    break;

  case 527:
#line 1788 "p.y"
    { (yyval.cxbclist) = new ComplexBCList; (yyval.cxbclist)->add((yyvsp[(1) - (1)].cxbcval)); ;}
    break;

  case 528:
#line 1790 "p.y"
    { (yyval.cxbclist) = (yyvsp[(1) - (2)].cxbclist); (yyval.cxbclist)->add((yyvsp[(2) - (2)].cxbcval)); ;}
    break;

  case 531:
#line 1798 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval); 
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.Iyy = (yyvsp[(14) - (16)].fval); sp.Izz = (yyvsp[(15) - (16)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 532:
#line 1806 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.betaDamp = (yyvsp[(17) - (19)].fval); sp.alphaDamp = (yyvsp[(18) - (19)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        ;}
    break;

  case 533:
#line 1815 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (19)].fval);  sp.E = (yyvsp[(3) - (19)].fval);  sp.nu  = (yyvsp[(4) - (19)].fval);  sp.rho = (yyvsp[(5) - (19)].fval);
          sp.c = (yyvsp[(6) - (19)].fval);  sp.k = (yyvsp[(7) - (19)].fval);  sp.eh  = (yyvsp[(8) - (19)].fval);  sp.P   = (yyvsp[(9) - (19)].fval);  sp.Ta  = (yyvsp[(10) - (19)].fval);
          sp.Q = (yyvsp[(11) - (19)].fval); sp.W = (yyvsp[(12) - (19)].fval); sp.Ixx = (yyvsp[(13) - (19)].fval); sp.Iyy = (yyvsp[(14) - (19)].fval); sp.Izz = (yyvsp[(15) - (19)].fval);
          sp.lagrangeMult = bool((yyvsp[(17) - (19)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(18) - (19)].fval);
          sp.type = StructProp::Constraint;
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (19)].ival)-1, sp );
        ;}
    break;

  case 534:
#line 1826 "p.y"
    { StructProp sp; 
	  sp.A = (yyvsp[(2) - (20)].fval);  sp.E = (yyvsp[(3) - (20)].fval);  sp.nu  = (yyvsp[(4) - (20)].fval);  sp.rho = (yyvsp[(5) - (20)].fval);
          sp.c = (yyvsp[(6) - (20)].fval);  sp.k = (yyvsp[(7) - (20)].fval);  sp.eh  = (yyvsp[(8) - (20)].fval);  sp.P   = (yyvsp[(9) - (20)].fval);  sp.Ta  = (yyvsp[(10) - (20)].fval); 
          sp.Q = (yyvsp[(11) - (20)].fval); sp.W = (yyvsp[(12) - (20)].fval); sp.Ixx = (yyvsp[(13) - (20)].fval); sp.Iyy = (yyvsp[(14) - (20)].fval); sp.Izz = (yyvsp[(15) - (20)].fval);
	  sp.ymin = (yyvsp[(16) - (20)].fval); sp.ymax = (yyvsp[(17) - (20)].fval); sp.zmin = (yyvsp[(18) - (20)].fval); sp.zmax = (yyvsp[(19) - (20)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (20)].ival)-1, sp );
        ;}
    break;

  case 535:
#line 1835 "p.y"
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

  case 536:
#line 1845 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (23)].fval);  sp.E = (yyvsp[(3) - (23)].fval);  sp.nu  = (yyvsp[(4) - (23)].fval);  sp.rho = (yyvsp[(5) - (23)].fval);
          sp.c = (yyvsp[(6) - (23)].fval);  sp.k = (yyvsp[(7) - (23)].fval);  sp.eh  = (yyvsp[(8) - (23)].fval);  sp.P   = (yyvsp[(9) - (23)].fval);  sp.Ta  = (yyvsp[(10) - (23)].fval);
          sp.Q = (yyvsp[(11) - (23)].fval); sp.W = (yyvsp[(12) - (23)].fval); sp.Ixx = (yyvsp[(13) - (23)].fval); sp.Iyy = (yyvsp[(14) - (23)].fval); sp.Izz = (yyvsp[(15) - (23)].fval);
          sp.ymin = (yyvsp[(16) - (23)].fval); sp.ymax = (yyvsp[(17) - (23)].fval); sp.zmin = (yyvsp[(18) - (23)].fval); sp.zmax = (yyvsp[(19) - (23)].fval);
          sp.lagrangeMult = bool((yyvsp[(21) - (23)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(22) - (23)].fval);
          sp.type = StructProp::Constraint;
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (23)].ival)-1, sp );
        ;}
    break;

  case 537:
#line 1857 "p.y"
    { StructProp sp;
          sp.A = (yyvsp[(2) - (9)].fval); sp.E = (yyvsp[(3) - (9)].fval); sp.nu = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(5) - (9)].fval);
          sp.c = (yyvsp[(6) - (9)].fval); sp.k = (yyvsp[(7) - (9)].fval); sp.eh = (yyvsp[(8) - (9)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp ); 
        ;}
    break;

  case 538:
#line 1864 "p.y"
    { StructProp sp;  // this is for spring: GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (14)].fval);  sp.E = (yyvsp[(3) - (14)].fval);  sp.nu  = (yyvsp[(4) - (14)].fval);  sp.rho = (yyvsp[(5) - (14)].fval);
          sp.c = (yyvsp[(6) - (14)].fval);  sp.k = (yyvsp[(7) - (14)].fval);  sp.eh  = (yyvsp[(8) - (14)].fval);  sp.P   = (yyvsp[(9) - (14)].fval);  sp.Ta  = (yyvsp[(10) - (14)].fval);
          sp.Q = (yyvsp[(11) - (14)].fval); sp.W = (yyvsp[(12) - (14)].fval); sp.Ixx = (yyvsp[(13) - (14)].fval);  
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (14)].ival)-1, sp );
        ;}
    break;

  case 539:
#line 1872 "p.y"
    { StructProp sp;  // this is for spring with stiffness-proportional damping : GID Kx Ky Kz lx1 ...
          sp.A = (yyvsp[(2) - (16)].fval);  sp.E = (yyvsp[(3) - (16)].fval);  sp.nu  = (yyvsp[(4) - (16)].fval);  sp.rho = (yyvsp[(5) - (16)].fval);
          sp.c = (yyvsp[(6) - (16)].fval);  sp.k = (yyvsp[(7) - (16)].fval);  sp.eh  = (yyvsp[(8) - (16)].fval);  sp.P   = (yyvsp[(9) - (16)].fval);  sp.Ta  = (yyvsp[(10) - (16)].fval);
          sp.Q = (yyvsp[(11) - (16)].fval); sp.W = (yyvsp[(12) - (16)].fval); sp.Ixx = (yyvsp[(13) - (16)].fval); sp.betaDamp = (yyvsp[(15) - (16)].fval);
          sp.isReal = true;
          geoSource->addMat( (yyvsp[(1) - (16)].ival)-1, sp );
        ;}
    break;

  case 540:
#line 1880 "p.y"
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
        ;}
    break;

  case 541:
#line 1896 "p.y"
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
        ;}
    break;

  case 542:
#line 1913 "p.y"
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
        ;}
    break;

  case 543:
#line 1931 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (5)].fval);
          sp.rho = (yyvsp[(4) - (5)].fval);
          sp.isReal = true;
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        ;}
    break;

  case 544:
#line 1939 "p.y"
    { StructProp sp;
          sp.kappaHelm = (yyvsp[(3) - (6)].fval);
          sp.kappaHelmImag = (yyvsp[(4) - (6)].fval);
          sp.rho = (yyvsp[(5) - (6)].fval);
          sp.type = StructProp::Fluid;
          geoSource->addMat( (yyvsp[(1) - (6)].ival)-1, sp );
        ;}
    break;

  case 545:
#line 1947 "p.y"
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

  case 546:
#line 1966 "p.y"
    { StructProp sp; 
          sp.A = (yyvsp[(3) - (12)].fval);  sp.rho = (yyvsp[(4) - (12)].fval); sp.Q = (yyvsp[(5) - (12)].fval); sp.c = (yyvsp[(6) - (12)].fval); 
          sp.sigma = (yyvsp[(7) - (12)].fval);  sp.k = (yyvsp[(8) - (12)].fval);  sp.eh  = (yyvsp[(9) - (12)].fval);  sp.P   = (yyvsp[(10) - (12)].fval);  sp.Ta  = (yyvsp[(11) - (12)].fval);
          sp.isReal = true;
          sp.type = StructProp::Thermal;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 547:
#line 1974 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (5)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (5)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (5)].ival)-1, sp );
        ;}
    break;

  case 548:
#line 1981 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (7)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (7)].fval);
          sp.amplitude = (yyvsp[(5) - (7)].fval);
          sp.omega = (yyvsp[(6) - (7)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 549:
#line 1990 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (8)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (8)].fval);
          sp.amplitude = (yyvsp[(5) - (8)].fval);
          sp.omega = (yyvsp[(6) - (8)].fval);
          sp.phase = (yyvsp[(7) - (8)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        ;}
    break;

  case 550:
#line 2000 "p.y"
    { StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (9)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (9)].fval);
          sp.amplitude = (yyvsp[(5) - (9)].fval);
          sp.omega = (yyvsp[(6) - (9)].fval);
          sp.phase = (yyvsp[(7) - (9)].fval);
          sp.offset = (yyvsp[(8) - (9)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        ;}
    break;

  case 551:
#line 2011 "p.y"
    { StructProp sp;
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

  case 552:
#line 2024 "p.y"
    { // new style for joints with prescribed motion by 2-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (10)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (10)].fval);
          sp.funtype = (yyvsp[(5) - (10)].ival);
          sp.amplitude = (yyvsp[(6) - (10)].fval);
          sp.offset = (yyvsp[(7) - (10)].fval);
          sp.c1 = (yyvsp[(8) - (10)].fval);
          sp.c2 = (yyvsp[(9) - (10)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (10)].ival)-1, sp );
        ;}
    break;

  case 553:
#line 2037 "p.y"
    { // new style for joints with prescribed motion by 3-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (11)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (11)].fval);
          sp.funtype = (yyvsp[(5) - (11)].ival);
          sp.amplitude = (yyvsp[(6) - (11)].fval);
          sp.offset = (yyvsp[(7) - (11)].fval);
          sp.c1 = (yyvsp[(8) - (11)].fval);
          sp.c2 = (yyvsp[(9) - (11)].fval);
          sp.c3 = (yyvsp[(10) - (11)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (11)].ival)-1, sp );
        ;}
    break;

  case 554:
#line 2051 "p.y"
    { // new style for joints with prescribed motion by 4-parameter elementary function
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (12)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (12)].fval);
          sp.funtype = (yyvsp[(5) - (12)].ival);
          sp.amplitude = (yyvsp[(6) - (12)].fval);
          sp.offset = (yyvsp[(7) - (12)].fval);
          sp.c1 = (yyvsp[(8) - (12)].fval);
          sp.c2 = (yyvsp[(9) - (12)].fval);
          sp.c3 = (yyvsp[(10) - (12)].fval);
          sp.c4 = (yyvsp[(11) - (12)].fval);
          sp.type = StructProp::Constraint;
          geoSource->addMat( (yyvsp[(1) - (12)].ival)-1, sp );
        ;}
    break;

  case 555:
#line 2066 "p.y"
    { // use for RevoluteJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (7)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (7)].fval);
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(6) - (7)].fval);
          geoSource->addMat( (yyvsp[(1) - (7)].ival)-1, sp );
        ;}
    break;

  case 556:
#line 2075 "p.y"
    { // use for UniversalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (8)].ival));
          sp.penalty = (yyvsp[(4) - (8)].fval);
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(6) - (8)].fval);
          sp.k2 = (yyvsp[(7) - (8)].fval);
          geoSource->addMat( (yyvsp[(1) - (8)].ival)-1, sp );
        ;}
    break;

  case 557:
#line 2085 "p.y"
    { // use for SphericalJointSpringCombo
          StructProp sp;
          sp.lagrangeMult = bool((yyvsp[(3) - (9)].ival));
          sp.initialPenalty = sp.penalty = (yyvsp[(4) - (9)].fval);
          sp.type = StructProp::Constraint;
          sp.k1 = (yyvsp[(6) - (9)].fval);
          sp.k2 = (yyvsp[(7) - (9)].fval);
          sp.k3 = (yyvsp[(8) - (9)].fval);
          geoSource->addMat( (yyvsp[(1) - (9)].ival)-1, sp );
        ;}
    break;

  case 558:
#line 2096 "p.y"
    { // use for TorsionalSpringType1 or TranslationalSpring
          StructProp sp;
          sp.k1 = (yyvsp[(3) - (4)].fval);
          geoSource->addMat( (yyvsp[(1) - (4)].ival)-1, sp );
        ;}
    break;

  case 561:
#line 2108 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (3)].ival));
          (yyval.SurfObj)->SetReverseNormals(false);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 562:
#line 2114 "p.y"
    { if((yyvsp[(2) - (4)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (4)].ival));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 563:
#line 2120 "p.y"
    { if((yyvsp[(2) - (5)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (5)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (5)].fval));
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 564:
#line 2127 "p.y"
    { if((yyvsp[(2) - (6)].ival) == 0) { cerr << " *** ERROR: surface id must be non-zero integer\n"; exit(-1); } // zero reserved for self-contact
          (yyval.SurfObj) = new SurfaceEntity((yyvsp[(2) - (6)].ival));
          (yyval.SurfObj)->SetIsShellFace(true);
          (yyval.SurfObj)->SetShellThickness((yyvsp[(4) - (6)].fval));
          (yyval.SurfObj)->SetReverseNormals(true);
          domain->AddSurfaceEntity((yyval.SurfObj));
        ;}
    break;

  case 565:
#line 2135 "p.y"
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
        ;}
    break;

  case 566:
#line 2159 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 567:
#line 2161 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 568:
#line 2163 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (5)].ival), (yyvsp[(3) - (5)].ival)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 569:
#line 2167 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (6)].ival), (yyvsp[(3) - (6)].ival), (yyvsp[(5) - (6)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 570:
#line 2169 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 571:
#line 2171 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 572:
#line 2173 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); domain->AddMortarCond((yyval.MortarCondObj)); ;}
    break;

  case 573:
#line 2175 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (7)].ival), (yyvsp[(3) - (7)].ival), (yyvsp[(6) - (7)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 574:
#line 2179 "p.y"
    { (yyval.MortarCondObj) = new MortarHandler((yyvsp[(2) - (8)].ival), (yyvsp[(3) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)); (yyval.MortarCondObj)->SetMortarType(MortarHandler::DUAL);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 575:
#line 2185 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].ival)); domain->solInfo().isCoupled = true; ;}
    break;

  case 576:
#line 2187 "p.y"
    { domain->addWetInterface((yyvsp[(2) - (3)].ival), (yyvsp[(2) - (3)].ival)); 
          domain->solInfo().isCoupled  = true; 
          domain->solInfo().isMatching = true; ;}
    break;

  case 577:
#line 2195 "p.y"
    { ;}
    break;

  case 578:
#line 2197 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 579:
#line 2203 "p.y"
    { 
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival));
          domain->AddMortarCond((yyval.MortarCondObj)); 
        ;}
    break;

  case 580:
#line 2210 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 581:
#line 2217 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 582:
#line 2224 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 583:
#line 2232 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival));
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 584:
#line 2241 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 585:
#line 2248 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::TIED);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival));
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 586:
#line 2258 "p.y"
    { ;}
    break;

  case 587:
#line 2260 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().isCoupled = true; 
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 588:
#line 2266 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().isCoupled = true;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 589:
#line 2274 "p.y"
    { domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; ;}
    break;

  case 591:
#line 2280 "p.y"
    { domain->addWetElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), 1.0, (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd);
          domain->solInfo().HEV = 1;
          domain->solInfo().isMatching = true; ;}
    break;

  case 592:
#line 2286 "p.y"
    { ;}
    break;

  case 593:
#line 2288 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (5)].ival) == (yyvsp[(4) - (5)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 594:
#line 2294 "p.y"
    { /* alternative format preferred by charbel. $2 is an id which will later be associated
           with an attribute and properties */
          domain->addWetInterface((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); domain->solInfo().HEV = 1;
          if((yyvsp[(3) - (7)].ival) == (yyvsp[(4) - (7)].ival)) domain->solInfo().isMatching = true;
        ;}
    break;

  case 595:
#line 2304 "p.y"
    { ;}
    break;

  case 596:
#line 2306 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (5)].ival), (yyvsp[(4) - (5)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 597:
#line 2313 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (6)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 598:
#line 2320 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].ival), (yyvsp[(6) - (7)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (7)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 599:
#line 2327 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].ival), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (8)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 600:
#line 2334 "p.y"
    { /* this one is for frictionless */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].ival), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (10)].ival), (yyvsp[(9) - (10)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (10)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 601:
#line 2342 "p.y"
    { /* this one is for constant friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].ival), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (11)].ival), (yyvsp[(9) - (11)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (11)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (11)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 602:
#line 2351 "p.y"
    { /* this one is for velocity dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (13)].ival), (yyvsp[(4) - (13)].ival), (yyvsp[(6) - (13)].fval), (yyvsp[(7) - (13)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (13)].ival), (yyvsp[(9) - (13)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (13)].fval), (yyvsp[(11) - (13)].fval), (yyvsp[(12) - (13)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (13)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 603:
#line 2360 "p.y"
    { /* this one is for pressure dependent friction */
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (14)].ival), (yyvsp[(4) - (14)].ival), (yyvsp[(6) - (14)].fval), (yyvsp[(7) - (14)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetTDEnfParams((yyvsp[(8) - (14)].ival), (yyvsp[(9) - (14)].fval));
          (yyval.MortarCondObj)->SetFrictionCoef((yyvsp[(10) - (14)].fval), (yyvsp[(11) - (14)].fval), (yyvsp[(12) - (14)].fval), (yyvsp[(13) - (14)].fval));
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (14)].ival)); 
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 604:
#line 2369 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (6)].ival), (yyvsp[(4) - (6)].ival)); 
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC); 
          (yyval.MortarCondObj)->SetMortarType(MortarHandler::STD); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(5) - (6)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 605:
#line 2377 "p.y"
    {
          (yyval.MortarCondObj) = new MortarHandler((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].ival), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval));
          (yyval.MortarCondObj)->SetInteractionType(MortarHandler::CTC);
          (yyval.MortarCondObj)->SetMortarType((yyvsp[(5) - (9)].ival)); 
          (yyval.MortarCondObj)->SetConstraintOptions((yyvsp[(8) - (9)].copt));
          domain->AddMortarCond((yyval.MortarCondObj));
        ;}
    break;

  case 606:
#line 2387 "p.y"
    { domain->solInfo().dist_acme = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 607:
#line 2389 "p.y"
    { domain->solInfo().ffi_debug = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 608:
#line 2391 "p.y"
    { domain->solInfo().mortar_scaling = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 609:
#line 2393 "p.y"
    { domain->solInfo().mortar_integration_rule = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 610:
#line 2396 "p.y"
    { geoSource->addNode((yyvsp[(3) - (3)].nval).num, (yyvsp[(3) - (3)].nval).xyz); ;}
    break;

  case 611:
#line 2398 "p.y"
    { geoSource->addNode((yyvsp[(2) - (2)].nval).num, (yyvsp[(2) - (2)].nval).xyz); ;}
    break;

  case 612:
#line 2402 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (5)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (5)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (5)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 613:
#line 2404 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (4)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (4)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (4)].fval);  (yyval.nval).xyz[2] = 0.0; ;}
    break;

  case 614:
#line 2406 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (3)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (3)].fval); (yyval.nval).xyz[1] = 0.0; (yyval.nval).xyz[2] = 0.0; ;}
    break;

  case 615:
#line 2408 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (7)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (7)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (7)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (7)].fval); /* $$.cp = $5; $$.cd = $6; */ ;}
    break;

  case 616:
#line 2410 "p.y"
    { (yyval.nval).num = (yyvsp[(1) - (6)].ival)-1; (yyval.nval).xyz[0] = (yyvsp[(2) - (6)].fval); (yyval.nval).xyz[1] = (yyvsp[(3) - (6)].fval);  (yyval.nval).xyz[2] = (yyvsp[(4) - (6)].fval); /* $$.cp = $$.cd = $5; */ ;}
    break;

  case 617:
#line 2414 "p.y"
    { /* Define each Element */
          geoSource->addElem((yyvsp[(1) - (4)].ival)-1, (yyvsp[(2) - (4)].ival), (yyvsp[(3) - (4)].nl).num, (yyvsp[(3) - (4)].nl).nd); ;}
    break;

  case 618:
#line 2419 "p.y"
    { (yyval.nl).num = 1; (yyval.nl).nd[0] = (yyvsp[(1) - (1)].ival)-1; ;}
    break;

  case 619:
#line 2421 "p.y"
    { if((yyval.nl).num == 125) return -1; 
          (yyval.nl).nd[(yyval.nl).num] = (yyvsp[(2) - (2)].ival)-1; (yyval.nl).num++; ;}
    break;

  case 620:
#line 2426 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.bcval).val = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 621:
#line 2428 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = (yyvsp[(2) - (3)].ival)-1; (yyval.bcval).val = 0.0; ;}
    break;

  case 622:
#line 2432 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1;  (yyval.bcval).dofnum = -1;  (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 623:
#line 2436 "p.y"
    { (yyval.bcval).nnum = (yyvsp[(1) - (3)].ival)-1; (yyval.bcval).dofnum = 6; (yyval.bcval).val = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 624:
#line 2440 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (5)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (5)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (5)].fval); (yyval.cxbcval).imval = (yyvsp[(4) - (5)].fval);  ;}
    break;

  case 625:
#line 2442 "p.y"
    { (yyval.cxbcval).nnum = (yyvsp[(1) - (4)].ival)-1; (yyval.cxbcval).dofnum = (yyvsp[(2) - (4)].ival)-1; (yyval.cxbcval).reval = (yyvsp[(3) - (4)].fval); (yyval.cxbcval).imval = 0.0; ;}
    break;

  case 627:
#line 2447 "p.y"
    { geoSource->setCSFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 629:
#line 2452 "p.y"
    { /*TODO: geoSource->setNodeFrame($2.num,$2.d);*/ ;}
    break;

  case 631:
#line 2457 "p.y"
    { geoSource->setFrame((yyvsp[(2) - (2)].frame).num,(yyvsp[(2) - (2)].frame).d); ;}
    break;

  case 632:
#line 2461 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (11)].ival)-1; 
            (yyval.frame).d[0] = (yyvsp[(2) - (11)].fval); (yyval.frame).d[1] = (yyvsp[(3) - (11)].fval); (yyval.frame).d[2] = (yyvsp[(4) - (11)].fval);
            (yyval.frame).d[3] = (yyvsp[(5) - (11)].fval); (yyval.frame).d[4] = (yyvsp[(6) - (11)].fval); (yyval.frame).d[5] = (yyvsp[(7) - (11)].fval);
            (yyval.frame).d[6] = (yyvsp[(8) - (11)].fval); (yyval.frame).d[7] = (yyvsp[(9) - (11)].fval); (yyval.frame).d[8] = (yyvsp[(10) - (11)].fval); ;}
    break;

  case 633:
#line 2466 "p.y"
    { (yyval.frame).num = (yyvsp[(1) - (4)].ival)-1; geoSource->makeEframe((yyvsp[(1) - (4)].ival)-1, (yyvsp[(3) - (4)].ival), (yyval.frame).d); ;}
    break;

  case 635:
#line 2471 "p.y"
    { OffsetData od;
	  od.first = (yyvsp[(2) - (7)].ival)-1; od.last = (yyvsp[(3) - (7)].ival)-1;
	  od.o[0] = (yyvsp[(4) - (7)].fval); od.o[1] = (yyvsp[(5) - (7)].fval); od.o[2] = (yyvsp[(6) - (7)].fval); 
	  geoSource->addOffset(od); ;}
    break;

  case 636:
#line 2478 "p.y"
    { (yyval.ival) = 0; ;}
    break;

  case 637:
#line 2480 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (4)].ival)-1,(yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 638:
#line 2482 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1); 
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (6)].ival) - 1, (yyvsp[(5) - (6)].fval));
	  domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 639:
#line 2486 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (6)].ival)-1,(yyvsp[(3) - (6)].ival)-1,(yyvsp[(4) - (6)].ival)-1,(yyvsp[(5) - (6)].ival)-1); ;}
    break;

  case 640:
#line 2488 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (8)].ival)-1,(yyvsp[(3) - (8)].ival)-1,(yyvsp[(4) - (8)].ival)-1,(yyvsp[(5) - (8)].ival)-1);
	  geoSource->setElementLumpingWeight((yyvsp[(2) - (8)].ival) - 1, (yyvsp[(7) - (8)].fval)); 
	  domain->solInfo().elemLumpPodRom = true; ;}
    break;

  case 641:
#line 2492 "p.y"
    { geoSource->setAttrib((yyvsp[(2) - (7)].ival)-1,(yyvsp[(3) - (7)].ival)-1,(yyvsp[(4) - (7)].ival)-1,-1,(yyvsp[(6) - (7)].fval)); ;}
    break;

  case 642:
#line 2494 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
            geoSource->setAttrib(i-1,i-1);
        ;}
    break;

  case 643:
#line 2499 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (5)].ival); i<(yyvsp[(3) - (5)].ival)+1; ++i)
 	    geoSource->setAttrib(i-1,(yyvsp[(4) - (5)].ival)-1);
	;}
    break;

  case 644:
#line 2504 "p.y"
    { int i;
	  for(i=(yyvsp[(2) - (7)].ival); i<(yyvsp[(3) - (7)].ival)+1; ++i)
	    geoSource->setAttrib(i-1, (yyvsp[(4) - (7)].ival)-1, (yyvsp[(5) - (7)].ival)-1, (yyvsp[(6) - (7)].ival)-1);
	;}
    break;

  case 645:
#line 2509 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (8)].ival); i<(yyvsp[(3) - (8)].ival)+1; ++i)
            geoSource->setAttrib(i-1, (yyvsp[(4) - (8)].ival)-1, (yyvsp[(5) - (8)].ival)-1, -1, (yyvsp[(7) - (8)].fval));
        ;}
    break;

  case 647:
#line 2517 "p.y"
    { geoSource->setElementPressure( (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval) ); ;}
    break;

  case 648:
#line 2519 "p.y"
    { int i; 
          for(i=(yyvsp[(2) - (5)].ival); i<((yyvsp[(3) - (5)].ival)+1); ++i)
            geoSource->setElementPressure( i-1, (yyvsp[(4) - (5)].fval) );  
        ;}
    break;

  case 649:
#line 2524 "p.y"
    { BCond *surf_pbc = new BCond[1];
          surf_pbc[0] = (yyvsp[(3) - (3)].bcval);
          geoSource->addSurfacePressure(1,surf_pbc); ;}
    break;

  case 650:
#line 2530 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false); 
          geoSource->setConsistentPFlag(false); 
        ;}
    break;

  case 651:
#line 2535 "p.y"
    { geoSource->setMRatio(0.0);
          geoSource->setConsistentQFlag(false, (yyvsp[(2) - (3)].ival));
          geoSource->setConsistentPFlag(false);
        ;}
    break;

  case 652:
#line 2542 "p.y"
    { ;}
    break;

  case 653:
#line 2544 "p.y"
    { geoSource->setElementPreLoad( (yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].fval) ); ;}
    break;

  case 654:
#line 2546 "p.y"
    { int i;
          for(i=(yyvsp[(2) - (6)].ival); i<((yyvsp[(4) - (6)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, (yyvsp[(5) - (6)].fval) );
        ;}
    break;

  case 655:
#line 2551 "p.y"
    { double load[3] = { (yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval) };
          geoSource->setElementPreLoad( (yyvsp[(2) - (6)].ival)-1, load ); ;}
    break;

  case 656:
#line 2554 "p.y"
    { double load[3] = { (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
          int i;
          for(i=(yyvsp[(2) - (8)].ival); i<((yyvsp[(4) - (8)].ival)+1); ++i)
            geoSource->setElementPreLoad( i-1, load );
        ;}
    break;

  case 660:
#line 2567 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(1) - (1)].ival)); ;}
    break;

  case 661:
#line 2569 "p.y"
    { domain->solInfo().loadcases.push_back((yyvsp[(2) - (2)].ival)); ;}
    break;

  case 662:
#line 2573 "p.y"
    { domain->solInfo().type = 1;
          domain->solInfo().iterType = (yyvsp[(3) - (4)].ival);
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 663:
#line 2577 "p.y"
    { domain->solInfo().precond = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 664:
#line 2579 "p.y"
    { domain->solInfo().maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 665:
#line 2581 "p.y"
    { domain->solInfo().tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 666:
#line 2583 "p.y"
    { domain->solInfo().maxvecsize = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 667:
#line 2585 "p.y"
    { domain->solInfo().iterSubtype = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 668:
#line 2589 "p.y"
    { domain->solInfo().setSolver(0); 
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 669:
#line 2592 "p.y"
    { domain->solInfo().setSolver((yyvsp[(4) - (5)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 670:
#line 2595 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 671:
#line 2598 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (5)].ival));
          domain->solInfo().setProbType(SolverInfo::Static);
          if((yyvsp[(3) - (5)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 672:
#line 2603 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (5)].ival));
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().getNLInfo().unsymmetric = true; ;}
    break;

  case 673:
#line 2607 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (5)].ival),(yyvsp[(4) - (5)].ival));    
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 674:
#line 2610 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (6)].ival),(yyvsp[(4) - (6)].ival),(yyvsp[(5) - (6)].fval));    
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 675:
#line 2613 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (7)].ival),(yyvsp[(4) - (7)].ival),(yyvsp[(5) - (7)].fval),(yyvsp[(6) - (7)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 676:
#line 2616 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (8)].ival),(yyvsp[(4) - (8)].ival),(yyvsp[(5) - (8)].fval),(yyvsp[(6) - (8)].ival),(yyvsp[(7) - (8)].ival)); 
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 677:
#line 2619 "p.y"
    { domain->solInfo().setSolver((yyvsp[(3) - (9)].ival),(yyvsp[(4) - (9)].ival),(yyvsp[(5) - (9)].fval),(yyvsp[(6) - (9)].ival),(yyvsp[(7) - (9)].ival),(yyvsp[(8) - (9)].ival));
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 678:
#line 2627 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = (yyvsp[(4) - (6)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(5) - (6)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(4) - (6)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 679:
#line 2633 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(4) - (5)].ival)-1); ;}
    break;

  case 680:
#line 2637 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; ;}
    break;

  case 681:
#line 2642 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; ;}
    break;

  case 682:
#line 2648 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(4) - (6)].ival)-1); 
          domain->solInfo().fetiInfo.feti2version 
                  = (FetiInfo::Feti2Version) (yyvsp[(5) - (6)].ival); ;}
    break;

  case 683:
#line 2654 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = (yyvsp[(4) - (7)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(5) - (7)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(6) - (7)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 684:
#line 2660 "p.y"
    { domain->solInfo().setProbType(SolverInfo::Static);
          domain->solInfo().fetiInfo.maxit    = (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.tol      = (yyvsp[(3) - (5)].fval);
          domain->solInfo().fetiInfo.maxortho = (yyvsp[(4) - (5)].ival);
          domain->solInfo().type =(2); ;}
    break;

  case 685:
#line 2666 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().setProbType(SolverInfo::Static); ;}
    break;

  case 686:
#line 2669 "p.y"
    { domain->solInfo().type =(2);;}
    break;

  case 687:
#line 2671 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (3)].ival)-1);;}
    break;

  case 688:
#line 2674 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp; ;}
    break;

  case 689:
#line 2678 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.corners = FetiInfo::allCorners6;
          domain->solInfo().fetiInfo.version = FetiInfo::fetidp;
          domain->solInfo().fetiInfo.dph_flag = true; ;}
    break;

  case 690:
#line 2683 "p.y"
    { domain->solInfo().type =(2);
          domain->solInfo().fetiInfo.version = (FetiInfo::Version) ((yyvsp[(2) - (4)].ival)-1);
          domain->solInfo().fetiInfo.feti2version = (FetiInfo::Feti2Version) (yyvsp[(3) - (4)].ival); 
        ;}
    break;

  case 691:
#line 2688 "p.y"
    {
	  domain->solInfo().type = 3;
          domain->solInfo().subtype = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getFetiInfo().solvertype = (FetiInfo::Solvertype)((yyvsp[(3) - (4)].ival));
	;}
    break;

  case 692:
#line 2694 "p.y"
    { domain->solInfo().sparse_maxsup  = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 693:
#line 2696 "p.y"
    { domain->solInfo().sparse_defblk  = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 694:
#line 2698 "p.y"
    { domain->solInfo().spooles_tau  = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 695:
#line 2700 "p.y"
    { domain->solInfo().spooles_maxsize = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 696:
#line 2702 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 0) {
            (yyvsp[(2) - (3)].ival) = 24;
            fprintf(stderr," *** WARNING: spooles_maxdomainsize must be > 0,"
                           " using 24\n");
          }
          domain->solInfo().spooles_maxdomainsize = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 697:
#line 2709 "p.y"
    { domain->solInfo().spooles_seed = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 698:
#line 2711 "p.y"
    { if(((yyvsp[(2) - (3)].fval) < 0.0) || ((yyvsp[(2) - (3)].fval) > 1.0)) {
            (yyvsp[(2) - (3)].fval) = 0.04;
            fprintf(stderr," *** WARNING: spooles_maxzeros outside acceptable limits (0..1),"
                           " using 0.04\n");
          }
          domain->solInfo().spooles_maxzeros = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 699:
#line 2718 "p.y"
    { domain->solInfo().spooles_msglvl = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 700:
#line 2720 "p.y"
    { domain->solInfo().pivot = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 701:
#line 2722 "p.y"
    { domain->solInfo().spooles_scale = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 702:
#line 2724 "p.y"
    { domain->solInfo().spooles_renum = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 703:
#line 2726 "p.y"
    { domain->solInfo().mumps_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 704:
#line 2728 "p.y"
    { domain->solInfo().mumps_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 705:
#line 2730 "p.y"
    { domain->solInfo().goldfarb_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 706:
#line 2732 "p.y"
    { domain->solInfo().goldfarb_check = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 707:
#line 2734 "p.y"
    { domain->solInfo().fetiInfo.maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 708:
#line 2736 "p.y"
    { domain->solInfo().debug_icntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 709:
#line 2738 "p.y"
    { domain->solInfo().debug_cntl[(yyvsp[(2) - (4)].ival)] = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 710:
#line 2744 "p.y"
    { domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 711:
#line 2746 "p.y"
    { domain->solInfo().fetiInfo.precno = FetiInfo::lumped; ;}
    break;

  case 712:
#line 2748 "p.y"
    { if(((yyvsp[(3) - (4)].ival) < 0) || ((yyvsp[(3) - (4)].ival) > 3)) { 
            (yyvsp[(3) - (4)].ival) = 1;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner selected, using lumped\n");
          }
          domain->solInfo().fetiInfo.precno = (FetiInfo::Preconditioner) (yyvsp[(3) - (4)].ival);
	;}
    break;

  case 713:
#line 2755 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        ;}
    break;

  case 714:
#line 2762 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 1)) {
            (yyvsp[(2) - (3)].ival) = 0;
            fprintf(stderr," *** WARNING: Incorrect Preconditioner Type selected, using nonshifted\n");
          }
          domain->solInfo().fetiInfo.prectype = (FetiInfo::PreconditionerType) (yyvsp[(2) - (3)].ival);
        ;}
    break;

  case 715:
#line 2769 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 716:
#line 2771 "p.y"
    { domain->solInfo().fetiInfo.tol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.absolute_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 717:
#line 2774 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 718:
#line 2776 "p.y"
    { domain->solInfo().fetiInfo.stagnation_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.absolute_stagnation_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 719:
#line 2779 "p.y"
    { domain->solInfo().fetiInfo.primal_proj_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_proj_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 720:
#line 2782 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.dual_plan_maxit = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 721:
#line 2785 "p.y"
    { domain->solInfo().fetiInfo.primal_plan_tol = (yyvsp[(2) - (4)].fval);
          domain->solInfo().fetiInfo.dual_plan_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 722:
#line 2788 "p.y"
    { domain->solInfo().fetiInfo.maxortho = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 723:
#line 2790 "p.y"
    { domain->solInfo().fetiInfo.noCoarse = 1; ;}
    break;

  case 724:
#line 2792 "p.y"
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

  case 725:
#line 2812 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 1; 
          domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 726:
#line 2815 "p.y"
    { domain->solInfo().fetiInfo.scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 727:
#line 2817 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 728:
#line 2820 "p.y"
    { domain->solInfo().fetiInfo.mpc_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 729:
#line 2822 "p.y"
    { if(((yyvsp[(2) - (3)].ival) < 0) || ((yyvsp[(2) - (3)].ival) > 2)) (yyvsp[(2) - (3)].ival) = 2;
          domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 730:
#line 2825 "p.y"
    { domain->solInfo().fetiInfo.fsi_scaling = (FetiInfo::Scaling) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 731:
#line 2827 "p.y"
    { domain->solInfo().fetiInfo.mpc_element = true; ;}
    break;

  case 732:
#line 2829 "p.y"
    { domain->solInfo().fetiInfo.fsi_element = true; ;}
    break;

  case 733:
#line 2831 "p.y"
    { domain->solInfo().fetiInfo.fsi_corner = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 734:
#line 2833 "p.y"
    { domain->solInfo().fetiInfo.splitLocalFsi = false; ;}
    break;

  case 735:
#line 2835 "p.y"
    { domain->solInfo().coupled_scale = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 736:
#line 2837 "p.y"
    { domain->solInfo().fetiInfo.wetcorners = true; ;}
    break;

  case 737:
#line 2839 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 738:
#line 2841 "p.y"
    { domain->solInfo().fetiInfo.corners = (FetiInfo::CornerType) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.pick_unsafe_corners = bool((yyvsp[(3) - (4)].ival));
        ;}
    break;

  case 739:
#line 2845 "p.y"
    { if((yyvsp[(2) - (3)].ival) == 0) {
            domain->solInfo().fetiInfo.corners = FetiInfo::noCorners;
            domain->solInfo().fetiInfo.pickAnyCorner = 0; 
            domain->solInfo().fetiInfo.bmpc = true;
            domain->solInfo().fetiInfo.pick_unsafe_corners = false;
            domain->solInfo().fetiInfo.augment = FetiInfo::none;
          }
        ;}
    break;

  case 740:
#line 2854 "p.y"
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

  case 741:
#line 2875 "p.y"
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

  case 742:
#line 2893 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (4)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 743:
#line 2898 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (5)].ival); 
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 744:
#line 2904 "p.y"
    { domain->solInfo().fetiInfo.numdir = (yyvsp[(3) - (5)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(4) - (5)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 745:
#line 2910 "p.y"
    { domain->solInfo().fetiInfo.waveType = (FetiInfo::WaveType) (yyvsp[(3) - (6)].ival);
          domain->solInfo().fetiInfo.waveMethod = (FetiInfo::WaveMethod) (yyvsp[(5) - (6)].ival);
          domain->solInfo().fetiInfo.numdir = (yyvsp[(4) - (6)].ival);
          if(domain->solInfo().fetiInfo.augment == FetiInfo::none)
            domain->solInfo().fetiInfo.augment = FetiInfo::Edges;
          /*geoSource->initShift();*/  ;}
    break;

  case 746:
#line 2917 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 747:
#line 2919 "p.y"
    { domain->solInfo().fetiInfo.orthotol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().fetiInfo.orthotol2 = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 748:
#line 2922 "p.y"
    { domain->solInfo().fetiInfo.grbm_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 749:
#line 2924 "p.y"
    { domain->solInfo().fetiInfo.crbm_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 750:
#line 2926 "p.y"
    { domain->solInfo().fetiInfo.cct_tol = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 751:
#line 2928 "p.y"
    { domain->solInfo().fetiInfo.rebuildcct = int((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 752:
#line 2930 "p.y"
    { domain->solInfo().fetiInfo.uproj = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 753:
#line 2932 "p.y"
    { domain->solInfo().fetiInfo.printMatLab = 1; ;}
    break;

  case 754:
#line 2934 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 755:
#line 2936 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 756:
#line 2938 "p.y"
    {  domain->solInfo().fetiInfo.auxCoarseSolver = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 757:
#line 2940 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 758:
#line 2942 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 759:
#line 2946 "p.y"
    { domain->solInfo().fetiInfo.solvertype = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().localScaled = true; ;}
    break;

  case 760:
#line 2949 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 761:
#line 2953 "p.y"
    { domain->solInfo().fetiInfo.auxCoarseSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 762:
#line 2957 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival);
          domain->solInfo().coarseScaled = true; ;}
    break;

  case 763:
#line 2960 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival) < 8) fprintf(stderr," *** WARNING: Pivoting not supported for this solver \n");
          else domain->solInfo().pivot = true; ;}
    break;

  case 764:
#line 2964 "p.y"
    { domain->solInfo().fetiInfo.cctSolver  = (FetiInfo::Solvertype) (yyvsp[(2) - (4)].ival); 
          if((yyvsp[(2) - (4)].ival)!=0) fprintf(stderr," *** WARNING: Scaling not supported for this CCt solver \n");
          else domain->solInfo().fetiInfo.cctScaled = true; ;}
    break;

  case 765:
#line 2968 "p.y"
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

  case 766:
#line 2994 "p.y"
    { domain->solInfo().fetiInfo.gtgSolver  = (FetiInfo::Solvertype) (yyvsp[(1) - (2)].ival); ;}
    break;

  case 767:
#line 2996 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = true; ;}
    break;

  case 768:
#line 2998 "p.y"
    { domain->solInfo().fetiInfo.gmresResidual = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 769:
#line 3000 "p.y"
    { domain->solInfo().fetiInfo.pickAnyCorner = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 770:
#line 3006 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
          domain->solInfo().fetiInfo.nlPrecFlg = 1; 
          domain->solInfo().setKrylov(); 
        ;}
    break;

  case 771:
#line 3011 "p.y"
    { domain->solInfo().fetiInfo.type = FetiInfo::nonlinear;
	  domain->solInfo().fetiInfo.nlPrecFlg = (yyvsp[(2) - (3)].ival);
	  domain->solInfo().setKrylov();
	;}
    break;

  case 772:
#line 3016 "p.y"
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

  case 773:
#line 3029 "p.y"
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

  case 774:
#line 3043 "p.y"
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

  case 775:
#line 3060 "p.y"
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

  case 776:
#line 3076 "p.y"
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

  case 777:
#line 3086 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 778:
#line 3088 "p.y"
    { domain->solInfo().fetiInfo.numcgm = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.numcgm2 = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 779:
#line 3091 "p.y"
    { domain->solInfo().fetiInfo.tolcgm = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 780:
#line 3093 "p.y"
    { domain->solInfo().fetiInfo.spaceDimension = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 781:
#line 3095 "p.y"
    { domain->solInfo().fetiInfo.krylovtype = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 782:
#line 3097 "p.y"
    { domain->solInfo().fetiInfo.krylovtype =  (yyvsp[(2) - (3)].ival); ;}
    break;

  case 783:
#line 3099 "p.y"
    { domain->solInfo().fetiInfo.lumpedinterface = 1; ;}
    break;

  case 784:
#line 3101 "p.y"
    { domain->solInfo().fetiInfo.saveMemCoarse = 1; ;}
    break;

  case 785:
#line 3103 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        ;}
    break;

  case 786:
#line 3108 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        ;}
    break;

  case 787:
#line 3114 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        ;}
    break;

  case 788:
#line 3121 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 789:
#line 3123 "p.y"
    { domain->solInfo().fetiInfo.outerloop = (FetiInfo::OuterloopType) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.complex_hermitian = true; ;}
    break;

  case 790:
#line 3126 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 791:
#line 3128 "p.y"
    { domain->solInfo().fetiInfo.mpcflag = (yyvsp[(2) - (3)].ival); ;}
    break;

  case 792:
#line 3130 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 793:
#line 3132 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (3)].ival); ;}
    break;

  case 794:
#line 3134 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 795:
#line 3137 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (5)].ival);
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 796:
#line 3140 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival); 
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 797:
#line 3143 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (4)].ival); ;}
    break;

  case 798:
#line 3146 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); ;}
    break;

  case 799:
#line 3150 "p.y"
    { domain->solInfo().fetiInfo.mpc_precno = (FetiInfo::MpcPreconditioner) (yyvsp[(2) - (6)].ival);
          domain->solInfo().fetiInfo.mpc_block = (FetiInfo::MpcBlock) (yyvsp[(3) - (6)].ival); 
          domain->solInfo().fetiInfo.mpcBlkOverlap = (yyvsp[(5) - (6)].ival); ;}
    break;

  case 800:
#line 3154 "p.y"
    { if((yyvsp[(2) - (3)].ival) < 1) domain->solInfo().fetiInfo.useMRHS = false; ;}
    break;

  case 801:
#line 3156 "p.y"
    { domain->solInfo().fetiInfo.gamma = (yyvsp[(2) - (3)].fval); ;}
    break;

  case 802:
#line 3158 "p.y"
    { domain->solInfo().fetiInfo.linesearch_maxit = (yyvsp[(2) - (4)].ival);
          domain->solInfo().fetiInfo.linesearch_tau = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 803:
#line 3161 "p.y"
    { domain->solInfo().fetiInfo.bmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 804:
#line 3163 "p.y"
    { domain->solInfo().fetiInfo.dmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 805:
#line 3165 "p.y"
    { domain->solInfo().fetiInfo.cmpc = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 806:
#line 3167 "p.y"
    { domain->solInfo().fetiInfo.c_normalize = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 807:
#line 3169 "p.y"
    { domain->solInfo().dbccheck = bool((yyvsp[(2) - (3)].ival)); ;}
    break;

  case 809:
#line 3176 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp; 
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::Helmholtz);
        ;}
    break;

  case 810:
#line 3187 "p.y"
    { if(!(yyvsp[(2) - (3)].copt).lagrangeMult && (yyvsp[(2) - (3)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(2) - (3)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(2) - (3)].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[(2) - (3)].copt).constraint_hess; 
          domain->solInfo().constraint_hess_eps = (yyvsp[(2) - (3)].copt).constraint_hess_eps; ;}
    break;

  case 811:
#line 3193 "p.y"
    { if(!(yyvsp[(3) - (4)].copt).lagrangeMult && (yyvsp[(3) - (4)].copt).penalty == 0) domain->solInfo().setDirectMPC(true);
          domain->solInfo().lagrangeMult = (yyvsp[(3) - (4)].copt).lagrangeMult;
          domain->solInfo().penalty = (yyvsp[(3) - (4)].copt).penalty;
          domain->solInfo().constraint_hess = (yyvsp[(3) - (4)].copt).constraint_hess;
          domain->solInfo().constraint_hess_eps = (yyvsp[(3) - (4)].copt).constraint_hess_eps; ;}
    break;

  case 812:
#line 3201 "p.y"
    { // Direct elimination of slave dofs
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
        ;}
    break;

  case 813:
#line 3208 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 814:
#line 3215 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (3)].fval);
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (3)].fval); ;}
    break;

  case 815:
#line 3223 "p.y"
    { (yyval.copt).lagrangeMult = false; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 0;
          (yyval.copt).constraint_hess_eps = 0.0;
          domain->solInfo().usePrescribedThreshold = true;
          domain->solInfo().mpcDirectTol = (yyvsp[(2) - (4)].fval); 
          domain->solInfo().coefFilterTol = (yyvsp[(3) - (4)].fval);
          domain->solInfo().rhsZeroTol = (yyvsp[(4) - (4)].fval); ;}
    break;

  case 816:
#line 3232 "p.y"
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

  case 817:
#line 3242 "p.y"
    { // Treatment of constraints through Lagrange multipliers method
          (yyval.copt).lagrangeMult = true; 
          (yyval.copt).penalty = 0.0;
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 818:
#line 3248 "p.y"
    { // Treatment of constraints through penalty method
          (yyval.copt).lagrangeMult = false;
          (yyval.copt).penalty = (yyvsp[(2) - (2)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 819:
#line 3254 "p.y"
    { // Treatment of constraints through augmented Lagrangian method
          (yyval.copt).lagrangeMult = true;
          (yyval.copt).penalty = (yyvsp[(3) - (3)].fval);
          (yyval.copt).constraint_hess = 1;
          (yyval.copt).constraint_hess_eps = 0.0; ;}
    break;

  case 820:
#line 3260 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (3)].ival);
          (yyval.copt).constraint_hess_eps = 0; ;}
    break;

  case 821:
#line 3263 "p.y"
    { (yyval.copt).constraint_hess = (yyvsp[(3) - (4)].ival);
          (yyval.copt).constraint_hess_eps = (yyvsp[(4) - (4)].fval); ;}
    break;

  case 822:
#line 3268 "p.y"
    { // hack??
	  domain->solInfo().acoustic = true; ;}
    break;

  case 823:
#line 3271 "p.y"
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

  case 824:
#line 3281 "p.y"
    {
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(2) - (3)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (3)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
        ;}
    break;

  case 825:
#line 3288 "p.y"
    {
          // this is for coupled problem: need fluid density and angular frequency
          if(domain->solInfo().probType != SolverInfo::HelmholtzDirSweep) domain->solInfo().setProbType(SolverInfo::Helmholtz);
          domain->fluidCelerity = (yyvsp[(3) - (5)].fval)/(yyvsp[(2) - (5)].fval);
          geoSource->setOmega((yyvsp[(3) - (5)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (5)].fval); sp.rho = (yyvsp[(4) - (5)].fval); geoSource->addMat(0,sp); // default flumat
        ;}
    break;

  case 826:
#line 3296 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (3)].ival);
          domain->curvatureFlag = 0;
        ;}
    break;

  case 827:
#line 3301 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (4)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (4)].fval);
          domain->curvatureFlag = 1;
        ;}
    break;

  case 828:
#line 3307 "p.y"
    {
          domain->sommerfeldType = (yyvsp[(2) - (5)].ival);
          domain->curvatureConst1 = (yyvsp[(3) - (5)].fval);
          domain->curvatureConst2 = (yyvsp[(4) - (5)].fval);
          domain->curvatureFlag = 2;
        ;}
    break;

  case 829:
#line 3314 "p.y"
    {
          domain->pointSourceFlag = 1;
          domain->implicitFlag = 1;
        ;}
    break;

  case 830:
#line 3319 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        ;}
    break;

  case 831:
#line 3324 "p.y"
    {
           domain->implicitFlag = 1;
           domain->pointSourceFlag = 0;
        ;}
    break;

  case 832:
#line 3329 "p.y"
    { 
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(4) - (7)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (7)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies1((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival)); 
        ;}
    break;

  case 833:
#line 3337 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(4) - (7)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (7)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies2((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].ival));
        ;}
    break;

  case 834:
#line 3345 "p.y"
    { 
           domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
           domain->fluidCelerity = 1.0; // defines the ratio omega/k
           geoSource->setOmega((yyvsp[(4) - (8)].fval)*domain->fluidCelerity);
           StructProp sp; sp.kappaHelm = (yyvsp[(4) - (8)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
           domain->addFrequencies((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].ival), (yyvsp[(7) - (8)].ival)); 
        ;}
    break;

  case 836:
#line 3354 "p.y"
    {
          // coupled sweep:            k0    delta_k n       omega_0 rho
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(7) - (9)].fval)/(yyvsp[(4) - (9)].fval);
          geoSource->setOmega((yyvsp[(7) - (9)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(8) - (9)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies1((yyvsp[(4) - (9)].fval)*domain->fluidCelerity, (yyvsp[(5) - (9)].fval)*domain->fluidCelerity, (yyvsp[(6) - (9)].ival));
        ;}
    break;

  case 837:
#line 3364 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(7) - (9)].fval)/(yyvsp[(4) - (9)].fval);
          geoSource->setOmega((yyvsp[(7) - (9)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (9)].fval); sp.rho = (yyvsp[(8) - (9)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies2((yyvsp[(4) - (9)].fval)*domain->fluidCelerity, (yyvsp[(5) - (9)].fval)*domain->fluidCelerity, (yyvsp[(6) - (9)].ival));
        ;}
    break;

  case 838:
#line 3373 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(8) - (10)].fval)/(yyvsp[(4) - (10)].fval);
          geoSource->setOmega((yyvsp[(8) - (10)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(4) - (10)].fval); sp.rho = (yyvsp[(9) - (10)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addFrequencies((yyvsp[(4) - (10)].fval)*domain->fluidCelerity, (yyvsp[(5) - (10)].fval)*domain->fluidCelerity, (yyvsp[(6) - (10)].ival), (yyvsp[(7) - (10)].ival));
        ;}
    break;

  case 842:
#line 3387 "p.y"
    { 
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = 1.0; // defines the ratio omega/k
          geoSource->setOmega((yyvsp[(2) - (3)].fval)*domain->fluidCelerity);
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (3)].fval); sp.rho = 1.0; geoSource->addMat(0, sp); // default flumat
          domain->addCoarseFrequency((yyvsp[(2) - (3)].fval)*domain->fluidCelerity); 
        ;}
    break;

  case 843:
#line 3396 "p.y"
    { 
          domain->solInfo().setProbType(SolverInfo::HelmholtzFreqSweep);
          domain->fluidCelerity = (yyvsp[(3) - (5)].fval)/(yyvsp[(2) - (5)].fval);
          geoSource->setOmega((yyvsp[(3) - (5)].fval));
          StructProp sp; sp.kappaHelm = (yyvsp[(2) - (5)].fval); sp.rho = (yyvsp[(4) - (5)].fval); geoSource->addMat(0, sp); // default flumat
          domain->addCoarseFrequency((yyvsp[(2) - (5)].fval)*domain->fluidCelerity);
        ;}
    break;

  case 844:
#line 3405 "p.y"
    { domain->addFrequencies((yyvsp[(2) - (4)].fval)*domain->fluidCelerity, (yyvsp[(3) - (4)].ival)); ;}
    break;

  case 847:
#line 3413 "p.y"
    { domain->setWaveDirections(0, (yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 848:
#line 3416 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(3) - (6)].fval), (yyvsp[(4) - (6)].fval), (yyvsp[(5) - (6)].fval));
        ;}
    break;

  case 849:
#line 3419 "p.y"
    {
          domain->setKirchhoffLocations((yyvsp[(2) - (5)].fval), (yyvsp[(3) - (5)].fval), (yyvsp[(4) - (5)].fval));
        ;}
    break;

  case 852:
#line 3429 "p.y"
    { domain->setFFPDirections((yyvsp[(1) - (4)].fval), (yyvsp[(2) - (4)].fval), (yyvsp[(3) - (4)].fval)); ;}
    break;

  case 854:
#line 3436 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzMF);
        ;}
    break;

  case 856:
#line 3450 "p.y"
    {
          /*domain->omega = $1;*/ geoSource->setOmega((yyvsp[(1) - (2)].fval));
          StructProp sp;
          sp.kappaHelm = (yyvsp[(1) - (2)].fval);
//          domain->setWaveNumber($1);
          geoSource->addMat(0,sp);
          domain->solInfo().setProbType(SolverInfo::HelmholtzSO);
        ;}
    break;

  case 857:
#line 3461 "p.y"
    {
          domain->solInfo().setProbType(SolverInfo::DisEnrM);
        ;}
    break;

  case 858:
#line 3467 "p.y"
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
        ;}
    break;

  case 859:
#line 3480 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::ArcLength;
        ;}
    break;

  case 860:
#line 3485 "p.y"
    { 
          if(domain->solInfo().probType == SolverInfo::NonLinStatic)
            domain->solInfo().probType = SolverInfo::MatNonLinStatic;
          else if(domain->solInfo().probType == SolverInfo::NonLinDynam)
            domain->solInfo().probType = SolverInfo::MatNonLinDynam;
        ;}
    break;

  case 861:
#line 3492 "p.y"
    { domain->solInfo().getNLInfo().maxiter = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 862:
#line 3494 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 863:
#line 3496 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 864:
#line 3499 "p.y"
    { domain->solInfo().getNLInfo().tolRes = (yyvsp[(3) - (7)].fval);
          domain->solInfo().getNLInfo().tolInc = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().absTolRes = (yyvsp[(5) - (7)].fval);
          domain->solInfo().getNLInfo().absTolInc = (yyvsp[(6) - (7)].fval); ;}
    break;

  case 865:
#line 3504 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (5)].fval);
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (5)].fval); ;}
    break;

  case 866:
#line 3507 "p.y"
    { domain->solInfo().getNLInfo().dlambda = (yyvsp[(3) - (7)].fval); 
          domain->solInfo().getNLInfo().maxLambda = (yyvsp[(4) - (7)].fval);
          domain->solInfo().getNLInfo().extMin = (yyvsp[(5) - (7)].ival);
          domain->solInfo().getNLInfo().extMax = (yyvsp[(6) - (7)].ival); ;}
    break;

  case 867:
#line 3512 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (4)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 868:
#line 3515 "p.y"
    { domain->solInfo().getNLInfo().fitAlgShell = (yyvsp[(3) - (5)].ival);
          domain->solInfo().getNLInfo().fitAlgBeam  = (yyvsp[(4) - (5)].ival); ;}
    break;

  case 869:
#line 3518 "p.y"
    { domain->solInfo().getNLInfo().unsymmetric = true; ;}
    break;

  case 870:
#line 3520 "p.y"
    { domain->solInfo().getNLInfo().lfactor = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 871:
#line 3526 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true; ;}
    break;

  case 872:
#line 3528 "p.y"
    { domain->solInfo().getNLInfo().failsafe = true;
          domain->solInfo().getNLInfo().failsafe_tol = (yyvsp[(3) - (4)].fval); ;}
    break;

  case 873:
#line 3531 "p.y"
    { domain->solInfo().momentType = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 874:
#line 3533 "p.y"
    { domain->solInfo().num_penalty_its = (yyvsp[(3) - (6)].ival); 
          domain->solInfo().penalty_tol = (yyvsp[(4) - (6)].fval);
          domain->solInfo().penalty_beta = (yyvsp[(5) - (6)].fval); ;}
    break;

  case 876:
#line 3540 "p.y"
    { 
          domain->solInfo().setNewton((yyvsp[(2) - (3)].ival)); 
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear; 
        ;}
    break;

  case 877:
#line 3545 "p.y"
    { 
          domain->solInfo().setNewton((yyvsp[(2) - (4)].ival)); 
          int rebuildK    = (yyvsp[(2) - (4)].ival); 
          int rebuildPrec = (yyvsp[(3) - (4)].ival);
          if(rebuildK > 1) rebuildPrec = rebuildK;
          domain->solInfo().fetiInfo.nTang = rebuildK;
          domain->solInfo().fetiInfo.nPrec = rebuildPrec;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        ;}
    break;

  case 878:
#line 3555 "p.y"
    {
          domain->solInfo().setNewton((yyvsp[(4) - (5)].ival));
          domain->solInfo().fetiInfo.nPrec = (yyvsp[(4) - (5)].ival);
          domain->solInfo().fetiInfo.nTang = (yyvsp[(4) - (5)].ival);
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
        ;}
    break;

  case 879:
#line 3562 "p.y"
    {
	  domain->solInfo().setNewton((yyvsp[(4) - (8)].ival)); 
          int rebuildK    = (yyvsp[(4) - (8)].ival); 
          int rebuildPrec = (yyvsp[(7) - (8)].ival);
          if(rebuildK > 1) rebuildPrec = rebuildK;
          domain->solInfo().fetiInfo.nTang = rebuildK;
          domain->solInfo().fetiInfo.nPrec = rebuildPrec;
          domain->solInfo().fetiInfo.type  = FetiInfo::nonlinear;
	;}
    break;

  case 880:
#line 3574 "p.y"
    { domain->solInfo().setReOrtho(); ;}
    break;

  case 882:
#line 3579 "p.y"
    { geoSource->setControl((yyvsp[(3) - (10)].strval),(yyvsp[(7) - (10)].strval),(yyvsp[(9) - (10)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (10)].ival); ;}
    break;

  case 883:
#line 3581 "p.y"
    { geoSource->setControl((yyvsp[(3) - (12)].strval),(yyvsp[(7) - (12)].strval),(yyvsp[(9) - (12)].strval),(yyvsp[(11) - (12)].strval)); domain->solInfo().soltyp = (yyvsp[(5) - (12)].ival); ;}
    break;

  case 884:
#line 3591 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(2) - (3)].strval);
#endif
        ;}
    break;

  case 885:
#line 3597 "p.y"
    { 
#ifdef STRUCTOPT
	  dynamic_cast<Domain_opt*>(domain)->setStructoptFlag(1); dynamic_cast<Domain_opt*>(domain)->optinputfile = (yyvsp[(3) - (4)].strval);
#endif
 ;}
    break;

  case 887:
#line 3606 "p.y"
    { domain->solInfo().contact_mode = (yyvsp[(3) - (4)].ival); ;}
    break;

  case 888:
#line 3610 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (7)].ival)-1, (yyvsp[(3) - (7)].ival)-1, (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)); ;}
    break;

  case 889:
#line 3613 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(8) - (9)].fval));;}
    break;

  case 890:
#line 3616 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (9)].ival)-1, (yyvsp[(3) - (9)].ival)-1, (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), 0.0, (yyvsp[(8) - (9)].ival));;}
    break;

  case 891:
#line 3618 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (8)].ival)-1, (yyvsp[(3) - (8)].ival)-1, (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), 0.0, -1, (yyvsp[(7) - (8)].copt).lagrangeMult, (yyvsp[(7) - (8)].copt).penalty);;}
    break;

  case 892:
#line 3621 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (11)].ival)-1, (yyvsp[(3) - (11)].ival)-1, (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(10) - (11)].ival));;}
    break;

  case 893:
#line 3623 "p.y"
    { domain->addNodalCTC((yyvsp[(2) - (12)].ival)-1, (yyvsp[(3) - (12)].ival)-1, (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(10) - (12)].ival), (yyvsp[(11) - (12)].copt).lagrangeMult, (yyvsp[(11) - (12)].copt).penalty);;}
    break;

  case 895:
#line 3628 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1, 
             new BilinPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 896:
#line 3633 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new BilinPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 897:
#line 3638 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 898:
#line 3643 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new FiniteStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 899:
#line 3648 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval)) );
         ;}
    break;

  case 900:
#line 3653 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new LogStrainPlasKinHardMat((yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval)) );
         ;}
    break;

  case 901:
#line 3658 "p.y"
    { 
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1, 
             new ElaLinIsoMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
	 ;}
    break;

  case 902:
#line 3663 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new StVenantKirchhoffMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
         ;}
    break;

  case 903:
#line 3668 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new HenckyMat((yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval)));
         ;}
    break;

  case 904:
#line 3673 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ElaLinIsoMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)));
         ;}
    break;

  case 905:
#line 3678 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new StVenantKirchhoffMat2D((yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval)));
         ;}
    break;

  case 906:
#line 3683 "p.y"
    {
            double params[3] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<IsotropicLinearElastic>(params));
          ;}
    break;

  case 907:
#line 3689 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          ;}
    break;

  case 908:
#line 3695 "p.y"
    {
            double params[4] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<NeoHookean>(params));
          ;}
    break;

  case 909:
#line 3701 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), -1 };
            geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          ;}
    break;

  case 910:
#line 3707 "p.y"
    {
            double params[5] = { (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
              new MaterialWrapper<MooneyRivlin>(params));
          ;}
    break;

  case 911:
#line 3713 "p.y"
    {
            double params[6] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>(params));
          ;}
    break;

  case 912:
#line 3719 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 1.0e-6, -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 913:
#line 3725 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), -std::numeric_limits<double>::infinity() };
            geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 914:
#line 3731 "p.y"
    {
            double params[8] = { (yyvsp[(4) - (12)].fval), (yyvsp[(5) - (12)].fval), (yyvsp[(6) - (12)].fval), (yyvsp[(7) - (12)].fval), (yyvsp[(8) - (12)].fval), (yyvsp[(9) - (12)].fval), (yyvsp[(10) - (12)].fval), (yyvsp[(11) - (12)].fval) };
            geoSource->addMaterial((yyvsp[(2) - (12)].ival)-1,
              new MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(params));
          ;}
    break;

  case 915:
#line 3737 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (24)].ival)-1,
             new ExpMat((yyvsp[(3) - (24)].ival), (yyvsp[(4) - (24)].fval), (yyvsp[(5) - (24)].fval), (yyvsp[(6) - (24)].fval), (yyvsp[(7) - (24)].fval), (yyvsp[(8) - (24)].fval), (yyvsp[(9) - (24)].fval), (yyvsp[(10) - (24)].fval), (yyvsp[(11) - (24)].fval), (yyvsp[(12) - (24)].fval), (yyvsp[(13) - (24)].fval), (yyvsp[(14) - (24)].fval), (yyvsp[(15) - (24)].fval), (yyvsp[(16) - (24)].fval), (yyvsp[(17) - (24)].fval), (yyvsp[(18) - (24)].fval), (yyvsp[(19) - (24)].fval), (yyvsp[(20) - (24)].fval), (yyvsp[(21) - (24)].fval), (yyvsp[(22) - (24)].fval), (yyvsp[(23) - (24)].fval)));
         ;}
    break;

  case 916:
#line 3742 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (7)].ival)-1,
             new ExpMat((yyvsp[(3) - (7)].ival), (yyvsp[(4) - (7)].fval), (yyvsp[(5) - (7)].fval), (yyvsp[(6) - (7)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 917:
#line 3747 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (8)].ival)-1,
             new ExpMat((yyvsp[(3) - (8)].ival), (yyvsp[(4) - (8)].fval), (yyvsp[(5) - (8)].fval), (yyvsp[(6) - (8)].fval), (yyvsp[(7) - (8)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 918:
#line 3752 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (9)].ival)-1,
             new ExpMat((yyvsp[(3) - (9)].ival), (yyvsp[(4) - (9)].fval), (yyvsp[(5) - (9)].fval), (yyvsp[(6) - (9)].fval), (yyvsp[(7) - (9)].fval), (yyvsp[(8) - (9)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 919:
#line 3757 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (10)].ival)-1,
             new ExpMat((yyvsp[(3) - (10)].ival), (yyvsp[(4) - (10)].fval), (yyvsp[(5) - (10)].fval), (yyvsp[(6) - (10)].fval), (yyvsp[(7) - (10)].fval), (yyvsp[(8) - (10)].fval), (yyvsp[(9) - (10)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 920:
#line 3762 "p.y"
    {
           geoSource->addMaterial((yyvsp[(2) - (11)].ival)-1,
             new ExpMat((yyvsp[(3) - (11)].ival), (yyvsp[(4) - (11)].fval), (yyvsp[(5) - (11)].fval), (yyvsp[(6) - (11)].fval), (yyvsp[(7) - (11)].fval), (yyvsp[(8) - (11)].fval), (yyvsp[(9) - (11)].fval), (yyvsp[(10) - (11)].fval), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
         ;}
    break;

  case 921:
#line 3767 "p.y"
    {
	   geoSource->loadMaterial((yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].strval));
	 ;}
    break;

  case 922:
#line 3771 "p.y"
    {
	   geoSource->addMaterial((yyvsp[(2) - (5)].ival)-1, (yyvsp[(3) - (5)].strval), (yyvsp[(4) - (5)].dlist));
	 ;}
    break;

  case 924:
#line 3778 "p.y"
    { geoSource->setMatUsage((yyvsp[(2) - (4)].ival)-1, (yyvsp[(3) - (4)].ival)-1); ;}
    break;

  case 925:
#line 3780 "p.y"
    {
            for(int i = (yyvsp[(2) - (5)].ival)-1; i < (yyvsp[(3) - (5)].ival); ++i)
	      geoSource->setMatUsage(i, (yyvsp[(4) - (5)].ival)-1);
	  ;}
    break;

  case 926:
#line 3786 "p.y"
    { (yyval.dlist).nval = 0; ;}
    break;

  case 927:
#line 3788 "p.y"
    { 
          if((yyvsp[(1) - (2)].dlist).nval == 32) {
             fprintf(stderr, "You'd better invent another material model!\n");
	     exit(-1);
          }
          (yyval.dlist) = (yyvsp[(1) - (2)].dlist);
          (yyval.dlist).v[(yyval.dlist).nval++] = (yyvsp[(2) - (2)].fval);
 	;}
    break;

  case 928:
#line 3799 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (4)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(3) - (4)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(3) - (4)].ival)); ;}
    break;

  case 929:
#line 3803 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (6)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (6)].ival)); ;}
    break;

  case 930:
#line 3806 "p.y"
    { domain->solInfo().setRenum((yyvsp[(3) - (8)].ival));
          domain->solInfo().setSparseRenum((yyvsp[(5) - (8)].ival)); 
          domain->solInfo().setSpoolesRenum((yyvsp[(7) - (8)].ival)); ;}
    break;

  case 931:
#line 3813 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().svdPodRom = true;;}
    break;

  case 933:
#line 3821 "p.y"
    { domain->solInfo().snapfiPodRom = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 934:
#line 3823 "p.y"
    { domain->solInfo().snapfiPodRom = (yyvsp[(2) - (3)].strval);
    if ((yyvsp[(3) - (3)].ival) == 1) domain->solInfo().statevectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 2) domain->solInfo().residvectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 3) domain->solInfo().jacobvectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 4) domain->solInfo().forcevectPodRom = true;
    if ((yyvsp[(3) - (3)].ival) == 5) domain->solInfo().accelvectPodRom = true;;}
    break;

  case 935:
#line 3830 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 936:
#line 3835 "p.y"
    { domain->solInfo().activatePodRom = true; 
    domain->solInfo().probType = SolverInfo::PodRomOffline;
    domain->solInfo().samplingPodRom = true; ;}
    break;

  case 938:
#line 3843 "p.y"
    { domain->solInfo().readInROBorModes = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 939:
#line 3845 "p.y"
    { domain->solInfo().statePodRomFile = (yyvsp[(2) - (2)].strval); ;}
    break;

  case 940:
#line 3847 "p.y"
    { domain->solInfo().tolPodRom = (yyvsp[(2) - (2)].fval); ;}
    break;

  case 941:
#line 3849 "p.y"
    { domain->solInfo().skipPodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 942:
#line 3851 "p.y"
    { domain->solInfo().skipOffSet = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 943:
#line 3853 "p.y"
    { domain->solInfo().maxSizePodRom = (yyvsp[(2) - (2)].ival); ;}
    break;

  case 944:
#line 3858 "p.y"
    { (yyval.ival) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 945:
#line 3863 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].ival); ;}
    break;

  case 946:
#line 3865 "p.y"
    { (yyval.fval) = (yyvsp[(1) - (1)].fval); ;}
    break;

  case 947:
#line 3867 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::infinity(); ;}
    break;

  case 948:
#line 3869 "p.y"
    { (yyval.fval) = std::numeric_limits<double>::epsilon(); ;}
    break;


/* Line 1267 of yacc.c.  */
#line 10235 "/lustre/home/mpotts/FEM/Parser.d/parser.cpp"
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


#line 3871 "p.y"


