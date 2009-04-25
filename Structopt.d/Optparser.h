/* A Bison parser, made by GNU Bison 1.875c.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

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
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ACCELERATION = 258,
     AERO = 259,
     ANALYSIS = 260,
     ANALYTIC = 261,
     ADJOINT = 262,
     AREA = 263,
     AOATTACK = 264,
     AUTOMATIC = 265,
     AVEFLAG = 266,
     CENTRAL = 267,
     COEF = 268,
     CONSTRAINT = 269,
     CONTINUE = 270,
     CONVECT = 271,
     CTRLCOST = 272,
     CRITERIA = 273,
     CTE = 274,
     CTE1 = 275,
     CTE2 = 276,
     DAMPING = 277,
     DEFOPR = 278,
     DENSITY = 279,
     DENSPROJ = 280,
     DISPLEVEL = 281,
     DISPLACMENT = 282,
     DIRECT = 283,
     DSGVAR = 284,
     DX = 285,
     DY = 286,
     DZ = 287,
     DCMASS = 288,
     ELCPOW = 289,
     ELCRHS = 290,
     ELCVOLT = 291,
     EMOD1 = 292,
     EMOD2 = 293,
     END = 294,
     ENERGY = 295,
     EQUALITY = 296,
     ESTAT = 297,
     ESTATCRIT = 298,
     FAILPROB = 299,
     FORWARD = 300,
     FREQUENCY = 301,
     FUNCRIT = 302,
     FUNVAR = 303,
     FUNOPR = 304,
     GENERATE = 305,
     GMOD = 306,
     GRDTYP = 307,
     GRDMTH = 308,
     GRDRELINC = 309,
     GRDABSINC = 310,
     GRDFILTER = 311,
     GRDSMI = 312,
     INEQUALITY = 313,
     INIT = 314,
     INERTIA = 315,
     INTFORCE = 316,
     KINETIC = 317,
     LAMBDA = 318,
     LAYER = 319,
     LCASE = 320,
     MASS = 321,
     MMAALG = 322,
     MOVINGNODES = 323,
     GCMALG = 324,
     MMAINT = 325,
     MMAFLT = 326,
     MOMAER = 327,
     NLPALG = 328,
     NLPINT = 329,
     NLPFLT = 330,
     NODELE = 331,
     NODDENS = 332,
     NODPOS = 333,
     NUMERIC = 334,
     NORM = 335,
     OBJECTIVE = 336,
     OCMALG = 337,
     OCMFLT = 338,
     OCMINT = 339,
     OPCOSINUS = 340,
     OPKSFUNC = 341,
     OPMULTI = 342,
     OPSINUS = 343,
     OPSUM = 344,
     OPEXP = 345,
     OPLOG = 346,
     OPUDEF = 347,
     PERMTVY = 348,
     PHI = 349,
     PMACRITERIA = 350,
     POISSON = 351,
     PULLIN = 352,
     RX = 353,
     RY = 354,
     RZ = 355,
     RELINDEX = 356,
     RHSDSGNDEP = 357,
     RGDVELOC = 358,
     SALALG = 359,
     SALFLT = 360,
     SALINT = 361,
     SBOOM = 362,
     SDTEMP = 363,
     SLPALG = 364,
     SLPFLT = 365,
     SLPINT = 366,
     SOLVER = 367,
     STCVAR = 368,
     STRESS = 369,
     STRLEVEL = 370,
     STRVOL = 371,
     S1L = 372,
     S1M = 373,
     S1U = 374,
     S2L = 375,
     S2M = 376,
     S2U = 377,
     SXL = 378,
     SXM = 379,
     SXU = 380,
     SYL = 381,
     SYM = 382,
     SYU = 383,
     SZL = 384,
     SZM = 385,
     SZU = 386,
     STCANA = 387,
     STCDYN = 388,
     STCDAM = 389,
     STCTRA = 390,
     STCELEVAR = 391,
     STCVARTYP = 392,
     TEMP = 393,
     THICK_opt = 394,
     THPOWER = 395,
     TIME = 396,
     THERMCOND = 397,
     TRADIATE = 398,
     VELOCITY = 399,
     VML = 400,
     VMM = 401,
     VMU = 402,
     VOLTAGE = 403,
     XMACH = 404,
     YAWANGLE = 405,
     YMODUS = 406,
     NewLine = 407,
     IntConstant = 408,
     DblConstant = 409,
     IPOPTALG = 410,
     IPOPT_INT = 411,
     IPOPT_FLT = 412,
     SNOPTALG = 413,
     SNOPT_INT = 414,
     SNOPT_FLT = 415
   };
#endif
#define ACCELERATION 258
#define AERO 259
#define ANALYSIS 260
#define ANALYTIC 261
#define ADJOINT 262
#define AREA 263
#define AOATTACK 264
#define AUTOMATIC 265
#define AVEFLAG 266
#define CENTRAL 267
#define COEF 268
#define CONSTRAINT 269
#define CONTINUE 270
#define CONVECT 271
#define CTRLCOST 272
#define CRITERIA 273
#define CTE 274
#define CTE1 275
#define CTE2 276
#define DAMPING 277
#define DEFOPR 278
#define DENSITY 279
#define DENSPROJ 280
#define DISPLEVEL 281
#define DISPLACMENT 282
#define DIRECT 283
#define DSGVAR 284
#define DX 285
#define DY 286
#define DZ 287
#define DCMASS 288
#define ELCPOW 289
#define ELCRHS 290
#define ELCVOLT 291
#define EMOD1 292
#define EMOD2 293
#define END 294
#define ENERGY 295
#define EQUALITY 296
#define ESTAT 297
#define ESTATCRIT 298
#define FAILPROB 299
#define FORWARD 300
#define FREQUENCY 301
#define FUNCRIT 302
#define FUNVAR 303
#define FUNOPR 304
#define GENERATE 305
#define GMOD 306
#define GRDTYP 307
#define GRDMTH 308
#define GRDRELINC 309
#define GRDABSINC 310
#define GRDFILTER 311
#define GRDSMI 312
#define INEQUALITY 313
#define INIT 314
#define INERTIA 315
#define INTFORCE 316
#define KINETIC 317
#define LAMBDA 318
#define LAYER 319
#define LCASE 320
#define MASS 321
#define MMAALG 322
#define MOVINGNODES 323
#define GCMALG 324
#define MMAINT 325
#define MMAFLT 326
#define MOMAER 327
#define NLPALG 328
#define NLPINT 329
#define NLPFLT 330
#define NODELE 331
#define NODDENS 332
#define NODPOS 333
#define NUMERIC 334
#define NORM 335
#define OBJECTIVE 336
#define OCMALG 337
#define OCMFLT 338
#define OCMINT 339
#define OPCOSINUS 340
#define OPKSFUNC 341
#define OPMULTI 342
#define OPSINUS 343
#define OPSUM 344
#define OPEXP 345
#define OPLOG 346
#define OPUDEF 347
#define PERMTVY 348
#define PHI 349
#define PMACRITERIA 350
#define POISSON 351
#define PULLIN 352
#define RX 353
#define RY 354
#define RZ 355
#define RELINDEX 356
#define RHSDSGNDEP 357
#define RGDVELOC 358
#define SALALG 359
#define SALFLT 360
#define SALINT 361
#define SBOOM 362
#define SDTEMP 363
#define SLPALG 364
#define SLPFLT 365
#define SLPINT 366
#define SOLVER 367
#define STCVAR 368
#define STRESS 369
#define STRLEVEL 370
#define STRVOL 371
#define S1L 372
#define S1M 373
#define S1U 374
#define S2L 375
#define S2M 376
#define S2U 377
#define SXL 378
#define SXM 379
#define SXU 380
#define SYL 381
#define SYM 382
#define SYU 383
#define SZL 384
#define SZM 385
#define SZU 386
#define STCANA 387
#define STCDYN 388
#define STCDAM 389
#define STCTRA 390
#define STCELEVAR 391
#define STCVARTYP 392
#define TEMP 393
#define THICK_opt 394
#define THPOWER 395
#define TIME 396
#define THERMCOND 397
#define TRADIATE 398
#define VELOCITY 399
#define VML 400
#define VMM 401
#define VMU 402
#define VOLTAGE 403
#define XMACH 404
#define YAWANGLE 405
#define YMODUS 406
#define NewLine 407
#define IntConstant 408
#define DblConstant 409
#define IPOPTALG 410
#define IPOPT_INT 411
#define IPOPT_FLT 412
#define SNOPTALG 413
#define SNOPT_INT 414
#define SNOPT_FLT 415




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 19 "Optparser.y"
typedef union YYSTYPE {
 int        ival;
 double     fval;
 criterion  crit;
 critdata   cdata;
 absvardata absvar;
 solverdata solv;
 nlpdata    nlp;
 graddata   grad;
 funcall    fall;
 funcdata*  func;
 funcdef    fdef;
 sumdata    sum;
 oprdata    opr;
 gendata    gen;
 anadata    anadat;
 stcelvdata sevdata;
 intlist    ilist;
 dllist     dpllist;
} YYSTYPE;
/* Line 1275 of yacc.c.  */
#line 378 "Optparser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yyoptlval;



