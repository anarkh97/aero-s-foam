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
     FORMALG = 258,
     FAILURE = 259,
     RNDVAR = 260,
     RNDTYPE = 261,
     MONTEALG = 262,
     MONTEINT = 263,
     MONTEFLT = 264,
     PMAFLAG = 265,
     STCRND = 266,
     ACCELERATION = 267,
     AERO = 268,
     ANALYSIS = 269,
     ANALYTIC = 270,
     ADJOINT = 271,
     AREA = 272,
     AOATTACK = 273,
     AUTOMATIC = 274,
     AVEFLAG = 275,
     CENTRAL = 276,
     CONSTRAINT = 277,
     CONTINUE = 278,
     CONVECT = 279,
     CTRLCOST = 280,
     CRITERIA = 281,
     CTE = 282,
     CTE1 = 283,
     CTE2 = 284,
     DAMPING = 285,
     DEFOPR = 286,
     DENSITY = 287,
     DISPLEVEL = 288,
     DISPLACMENT = 289,
     DIRECT = 290,
     DSGVAR = 291,
     DX = 292,
     DY = 293,
     DZ = 294,
     ELCPOW = 295,
     ELCRHS = 296,
     ELCVOLT = 297,
     EMOD1 = 298,
     EMOD2 = 299,
     END = 300,
     ENERGY = 301,
     EQUALITY = 302,
     ESTAT = 303,
     ESTATCRIT = 304,
     FAILPROB = 305,
     FORWARD = 306,
     FREQUENCY = 307,
     FUNCRIT = 308,
     FUNVAR = 309,
     FUNOPR = 310,
     GENERATE = 311,
     GMOD = 312,
     GRDTYP = 313,
     GRDMTH = 314,
     GRDRELINC = 315,
     GRDABSINC = 316,
     GRDFILTER = 317,
     GRDSMI = 318,
     INEQUALITY = 319,
     INIT = 320,
     INERTIA = 321,
     INTFORCE = 322,
     KINETIC = 323,
     LAYER = 324,
     LCASE = 325,
     MASS = 326,
     MIXED = 327,
     MMAALG = 328,
     GCMALG = 329,
     MMAINT = 330,
     MMAFLT = 331,
     MOMAER = 332,
     NLPALG = 333,
     NLPINT = 334,
     NLPFLT = 335,
     NODELE = 336,
     NODPOS = 337,
     NUMERIC = 338,
     OBJECTIVE = 339,
     OCMALG = 340,
     OCMFLT = 341,
     OCMINT = 342,
     OPCOSINUS = 343,
     OPKSFUNC = 344,
     OPMULTI = 345,
     OPSINUS = 346,
     OPSUM = 347,
     PERMTVY = 348,
     PHI = 349,
     PMACRITERIA = 350,
     PULLIN = 351,
     RX = 352,
     RY = 353,
     RZ = 354,
     RELINDEX = 355,
     RGDVELOC = 356,
     SALALG = 357,
     SALFLT = 358,
     SALINT = 359,
     SBOOM = 360,
     SDTEMP = 361,
     SLPALG = 362,
     SLPFLT = 363,
     SLPINT = 364,
     SOLVER = 365,
     STCVAR = 366,
     STRESS = 367,
     STRLEVEL = 368,
     STRVOL = 369,
     S1L = 370,
     S1M = 371,
     S1U = 372,
     S2L = 373,
     S2M = 374,
     S2U = 375,
     SXL = 376,
     SXM = 377,
     SXU = 378,
     SYL = 379,
     SYM = 380,
     SYU = 381,
     SZL = 382,
     SZM = 383,
     SZU = 384,
     STCANA = 385,
     STCDYN = 386,
     STCDAM = 387,
     STCTRA = 388,
     STCELEVAR = 389,
     STCVARTYP = 390,
     TEMP = 391,
     THICK_rel = 392,
     THPOWER = 393,
     TIME = 394,
     THERMCOND = 395,
     TRADIATE = 396,
     VELOCITY = 397,
     VML = 398,
     VMM = 399,
     VMU = 400,
     VOLTAGE = 401,
     XMACH = 402,
     YAWANGLE = 403,
     YMODUS = 404,
     NewLine = 405,
     IntConstant = 406,
     DblConstant = 407,
     String = 408
   };
#endif
#define FORMALG 258
#define FAILURE 259
#define RNDVAR 260
#define RNDTYPE 261
#define MONTEALG 262
#define MONTEINT 263
#define MONTEFLT 264
#define PMAFLAG 265
#define STCRND 266
#define ACCELERATION 267
#define AERO 268
#define ANALYSIS 269
#define ANALYTIC 270
#define ADJOINT 271
#define AREA 272
#define AOATTACK 273
#define AUTOMATIC 274
#define AVEFLAG 275
#define CENTRAL 276
#define CONSTRAINT 277
#define CONTINUE 278
#define CONVECT 279
#define CTRLCOST 280
#define CRITERIA 281
#define CTE 282
#define CTE1 283
#define CTE2 284
#define DAMPING 285
#define DEFOPR 286
#define DENSITY 287
#define DISPLEVEL 288
#define DISPLACMENT 289
#define DIRECT 290
#define DSGVAR 291
#define DX 292
#define DY 293
#define DZ 294
#define ELCPOW 295
#define ELCRHS 296
#define ELCVOLT 297
#define EMOD1 298
#define EMOD2 299
#define END 300
#define ENERGY 301
#define EQUALITY 302
#define ESTAT 303
#define ESTATCRIT 304
#define FAILPROB 305
#define FORWARD 306
#define FREQUENCY 307
#define FUNCRIT 308
#define FUNVAR 309
#define FUNOPR 310
#define GENERATE 311
#define GMOD 312
#define GRDTYP 313
#define GRDMTH 314
#define GRDRELINC 315
#define GRDABSINC 316
#define GRDFILTER 317
#define GRDSMI 318
#define INEQUALITY 319
#define INIT 320
#define INERTIA 321
#define INTFORCE 322
#define KINETIC 323
#define LAYER 324
#define LCASE 325
#define MASS 326
#define MIXED 327
#define MMAALG 328
#define GCMALG 329
#define MMAINT 330
#define MMAFLT 331
#define MOMAER 332
#define NLPALG 333
#define NLPINT 334
#define NLPFLT 335
#define NODELE 336
#define NODPOS 337
#define NUMERIC 338
#define OBJECTIVE 339
#define OCMALG 340
#define OCMFLT 341
#define OCMINT 342
#define OPCOSINUS 343
#define OPKSFUNC 344
#define OPMULTI 345
#define OPSINUS 346
#define OPSUM 347
#define PERMTVY 348
#define PHI 349
#define PMACRITERIA 350
#define PULLIN 351
#define RX 352
#define RY 353
#define RZ 354
#define RELINDEX 355
#define RGDVELOC 356
#define SALALG 357
#define SALFLT 358
#define SALINT 359
#define SBOOM 360
#define SDTEMP 361
#define SLPALG 362
#define SLPFLT 363
#define SLPINT 364
#define SOLVER 365
#define STCVAR 366
#define STRESS 367
#define STRLEVEL 368
#define STRVOL 369
#define S1L 370
#define S1M 371
#define S1U 372
#define S2L 373
#define S2M 374
#define S2U 375
#define SXL 376
#define SXM 377
#define SXU 378
#define SYL 379
#define SYM 380
#define SYU 381
#define SZL 382
#define SZM 383
#define SZU 384
#define STCANA 385
#define STCDYN 386
#define STCDAM 387
#define STCTRA 388
#define STCELEVAR 389
#define STCVARTYP 390
#define TEMP 391
#define THICK_rel 392
#define THPOWER 393
#define TIME 394
#define THERMCOND 395
#define TRADIATE 396
#define VELOCITY 397
#define VML 398
#define VMM 399
#define VMU 400
#define VOLTAGE 401
#define XMACH 402
#define YAWANGLE 403
#define YMODUS 404
#define NewLine 405
#define IntConstant 406
#define DblConstant 407
#define String 408




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 17 "Relparser.y"
typedef union YYSTYPE {
 int        ival;
 double     fval;
 char *     strval;
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
#line 365 "Relparser.h"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yyrellval;



