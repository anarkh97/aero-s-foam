/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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
/* Line 1529 of yacc.c.  */
#line 996 "/lustre/home/tac688/newCodes/Empirical_GNAT/Parser.d/parser.hpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

