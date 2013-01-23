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

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* If NAME_PREFIX is specified substitute the variables and functions
   names.  */
#define yyparse yyoptparse
#define yylex   yyoptlex
#define yyerror yyopterror
#define yylval  yyoptlval
#define yychar  yyoptchar
#define yydebug yyoptdebug
#define yynerrs yyoptnerrs


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




/* Copy the first part of user declarations.  */
#line 1 "Optparser.y"

#ifdef STRUCTOPT

#include <cstdio>
#include <cstdlib>
#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optsol.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Structopt_sd.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>

extern int yyoptlex(void);
extern void yyopterror(const char  *);
extern Domain * domain;


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
/* Line 191 of yacc.c.  */
#line 443 "Optparser.C"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 455 "Optparser.C"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   define YYSTACK_ALLOC alloca
#  endif
# else
#  if defined (alloca) || defined (_ALLOCA_H)
#   define YYSTACK_ALLOC alloca
#  else
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <cstdlib> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))				\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  register YYSIZE_T yyi;		\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (0)
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
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  43
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   901

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  176
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  58
/* YYNRULES -- Number of rules. */
#define YYNRULES  309
/* YYNRULES -- Number of states. */
#define YYNSTATES  651

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   415

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     173,   175,   157,   155,   174,   156,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   167,     2,
     169,   172,   168,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   170,     2,   171,   158,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   165,     2,   166,     2,     2,     2,     2,
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
     159,   160,   161,   162,   163,   164
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     6,     8,    11,    13,    15,    17,    19,
      21,    23,    25,    27,    29,    31,    33,    37,    40,    44,
      50,    55,    62,    67,    73,    80,    88,    90,    95,    97,
     102,   105,   109,   113,   118,   122,   126,   128,   130,   133,
     136,   138,   142,   144,   151,   161,   169,   177,   186,   195,
     198,   201,   204,   207,   210,   213,   216,   218,   222,   225,
     228,   230,   232,   235,   237,   239,   244,   247,   253,   255,
     257,   259,   261,   263,   265,   267,   269,   271,   273,   275,
     277,   279,   281,   283,   285,   287,   289,   291,   293,   295,
     297,   299,   301,   303,   305,   307,   310,   313,   319,   323,
     331,   337,   344,   349,   351,   353,   357,   360,   368,   375,
     380,   384,   392,   398,   403,   411,   418,   425,   431,   439,
     446,   448,   450,   452,   454,   456,   458,   460,   462,   464,
     466,   468,   470,   472,   474,   476,   478,   480,   484,   488,
     490,   492,   494,   496,   498,   500,   502,   504,   506,   508,
     510,   512,   515,   517,   521,   524,   528,   532,   537,   545,
     547,   549,   551,   553,   555,   560,   564,   569,   574,   579,
     584,   589,   594,   599,   604,   608,   613,   617,   622,   626,
     631,   635,   640,   644,   649,   653,   658,   662,   667,   671,
     676,   680,   685,   689,   694,   695,   699,   704,   708,   713,
     714,   718,   723,   727,   732,   736,   741,   746,   751,   755,
     768,   785,   787,   789,   791,   793,   795,   797,   801,   804,
     809,   815,   821,   828,   836,   843,   851,   860,   862,   864,
     866,   868,   870,   872,   874,   876,   879,   883,   887,   892,
     894,   897,   900,   904,   906,   909,   912,   920,   928,   934,
     940,   946,   950,   956,   962,   966,   970,   974,   976,   985,
     994,  1001,  1008,  1015,  1020,  1027,  1034,  1039,  1044,  1049,
    1052,  1058,  1062,  1067,  1072,  1077,  1084,  1093,  1100,  1102,
    1105,  1108,  1111,  1113,  1116,  1119,  1122,  1127,  1133,  1139,
    1142,  1147,  1153,  1159,  1162,  1164,  1166,  1168,  1170,  1173,
    1175,  1177,  1181,  1185,  1189,  1194,  1200,  1204,  1209,  1215
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const short yyrhs[] =
{
     177,     0,    -1,   178,    39,    -1,   179,    -1,   178,   179,
      -1,   180,    -1,   186,    -1,   187,    -1,   189,    -1,   191,
      -1,   196,    -1,   212,    -1,   200,    -1,   228,    -1,   233,
      -1,   232,    -1,    18,   152,   181,    -1,   180,   181,    -1,
     226,   182,   152,    -1,   226,   182,     5,   226,   152,    -1,
     226,   182,   185,   152,    -1,   226,   182,   185,     5,   226,
     152,    -1,   226,   182,   221,   152,    -1,   226,   182,   185,
     221,   152,    -1,   226,   182,     5,   226,   221,   152,    -1,
     226,   182,   185,     5,   226,   221,   152,    -1,    40,    -1,
      40,   165,   222,   166,    -1,    66,    -1,    66,   165,   222,
     166,    -1,    46,   226,    -1,   114,   226,   183,    -1,    27,
     226,   184,    -1,    27,   226,   184,   226,    -1,   144,   226,
     184,    -1,     3,   226,   184,    -1,    62,    -1,    22,    -1,
       4,   226,    -1,    72,   226,    -1,   107,    -1,    61,   226,
     184,    -1,    17,    -1,   116,   183,   227,   227,   226,    76,
      -1,   116,   183,   227,   227,   226,    76,   165,   222,   166,
      -1,    26,   226,   226,   227,   165,   224,   166,    -1,    26,
      11,   226,   227,   165,   224,   166,    -1,   115,   226,   226,
     227,    76,   165,   225,   166,    -1,   115,    11,   226,   227,
      76,   165,   225,   166,    -1,    44,   226,    -1,    42,    43,
      -1,    42,   226,    -1,   101,   226,    -1,    95,   226,    -1,
      60,   184,    -1,   138,   226,    -1,    97,    -1,    78,   226,
     184,    -1,    36,   226,    -1,    35,   226,    -1,    34,    -1,
     140,    -1,   108,   226,    -1,    63,    -1,    33,    -1,    33,
     165,   222,   166,    -1,    33,   227,    -1,    33,   227,   165,
     222,   166,    -1,   145,    -1,   146,    -1,   147,    -1,   117,
      -1,   118,    -1,   119,    -1,   120,    -1,   121,    -1,   122,
      -1,   123,    -1,   124,    -1,   125,    -1,   126,    -1,   127,
      -1,   128,    -1,   129,    -1,   130,    -1,   131,    -1,   226,
      -1,    30,    -1,    31,    -1,    32,    -1,    98,    -1,    99,
      -1,   100,    -1,    80,    -1,   138,    -1,    65,   226,    -1,
     141,   227,    -1,    81,   152,   227,   157,   214,    -1,    81,
     152,   214,    -1,    14,   152,   226,   188,   227,   157,   214,
      -1,    14,   152,   226,   188,   214,    -1,   187,   226,   188,
     227,   157,   214,    -1,   187,   226,   188,   214,    -1,    41,
      -1,    58,    -1,    29,   152,   190,    -1,   189,   190,    -1,
     226,   227,   227,   227,   227,   227,   152,    -1,   226,   227,
     227,   227,   227,   152,    -1,   226,   227,   227,   152,    -1,
     226,   227,   152,    -1,   226,   227,   227,   227,   227,   221,
     152,    -1,   226,   227,   227,   221,   152,    -1,   226,   227,
     221,   152,    -1,   113,   152,   226,   137,   226,   192,   214,
      -1,   191,   226,   137,   226,   192,   214,    -1,   113,   152,
     226,   137,   194,   214,    -1,   191,   226,   137,   194,   214,
      -1,   113,   152,   226,   137,   226,   195,   214,    -1,   191,
     226,   137,   226,   195,   214,    -1,   226,    -1,     8,    -1,
     151,    -1,    96,    -1,    24,    -1,    16,    -1,    19,    -1,
      93,    -1,   142,    -1,   143,    -1,   139,    -1,    30,    -1,
      31,    -1,    32,    -1,    98,    -1,    99,    -1,   100,    -1,
      64,   226,   193,    -1,    13,   226,   226,    -1,   226,    -1,
      37,    -1,    38,    -1,    51,    -1,    24,    -1,   139,    -1,
      94,    -1,    20,    -1,    21,    -1,   149,    -1,     9,    -1,
     150,    -1,   103,   184,    -1,   148,    -1,   132,   152,   197,
      -1,   196,   197,    -1,   133,   198,   152,    -1,   134,   199,
     152,    -1,   134,   199,   226,   152,    -1,   135,   226,   226,
     227,   227,   226,   152,    -1,   226,    -1,    59,    -1,    15,
      -1,   226,    -1,    10,    -1,   112,   152,   226,   201,    -1,
     200,   226,   201,    -1,    73,   152,   202,   209,    -1,    82,
     152,   203,   209,    -1,   109,   152,   204,   209,    -1,    67,
     152,   205,   209,    -1,   104,   152,   206,   209,    -1,    69,
     152,   205,   209,    -1,   159,   152,   207,   209,    -1,   162,
     152,   208,   209,    -1,    74,   226,   152,    -1,   202,    74,
     226,   152,    -1,    75,   227,   152,    -1,   202,    75,   227,
     152,    -1,    84,   226,   152,    -1,   203,    84,   226,   152,
      -1,    83,   227,   152,    -1,   203,    83,   227,   152,    -1,
     111,   226,   152,    -1,   204,   111,   226,   152,    -1,   110,
     227,   152,    -1,   204,   110,   227,   152,    -1,    70,   226,
     152,    -1,   205,    70,   226,   152,    -1,    71,   227,   152,
      -1,   205,    71,   227,   152,    -1,   106,   226,   152,    -1,
     206,   106,   226,   152,    -1,   105,   227,   152,    -1,   206,
     105,   227,   152,    -1,    -1,   160,   226,   152,    -1,   207,
     160,   226,   152,    -1,   161,   227,   152,    -1,   207,   161,
     227,   152,    -1,    -1,   163,   226,   152,    -1,   208,   163,
     226,   152,    -1,   164,   227,   152,    -1,   208,   164,   227,
     152,    -1,    52,   210,   152,    -1,   209,    53,   211,   152,
      -1,   209,    54,   227,   152,    -1,   209,    55,   227,   152,
      -1,   209,    57,   152,    -1,   209,    56,   226,   226,   227,
     227,   227,   227,   165,   223,   166,   152,    -1,   209,    56,
     226,   226,   227,   227,   227,   227,   165,   223,   166,   226,
     165,   222,   166,   152,    -1,     6,    -1,    79,    -1,    28,
      -1,     7,    -1,    45,    -1,    12,    -1,   136,   152,   213,
      -1,   212,   213,    -1,   226,   226,   137,   152,    -1,   226,
     226,   226,   137,   152,    -1,   215,   165,   216,   166,   152,
      -1,   215,   165,   152,   216,   166,   152,    -1,   215,   165,
     152,   216,   152,   166,   152,    -1,   215,   152,   165,   216,
     166,   152,    -1,   215,   152,   165,   152,   216,   166,   152,
      -1,   215,   152,   165,   152,   216,   152,   166,   152,    -1,
      89,    -1,    86,    -1,    87,    -1,    88,    -1,    85,    -1,
      90,    -1,    92,    -1,    91,    -1,   217,   218,    -1,   217,
     218,   152,    -1,   216,   217,   218,    -1,   216,   217,   218,
     152,    -1,   218,    -1,   218,   152,    -1,   216,   218,    -1,
     216,   218,   152,    -1,   220,    -1,   216,   220,    -1,   226,
     167,    -1,   227,   157,   219,   158,   227,   155,   227,    -1,
     227,   157,   219,   158,   227,   156,   227,    -1,   227,   157,
     219,   158,   227,    -1,   227,   157,   219,   155,   227,    -1,
     227,   157,   219,   156,   227,    -1,   227,   157,   219,    -1,
     219,   158,   227,   155,   227,    -1,   219,   158,   227,   156,
     227,    -1,   219,   158,   227,    -1,   219,   155,   227,    -1,
     219,   156,   227,    -1,   219,    -1,   227,   157,   219,   158,
     227,   155,   227,   221,    -1,   227,   157,   219,   158,   227,
     156,   227,   221,    -1,   227,   157,   219,   158,   227,   221,
      -1,   227,   157,   219,   155,   227,   221,    -1,   227,   157,
     219,   156,   227,   221,    -1,   227,   157,   219,   221,    -1,
     219,   158,   227,   155,   227,   221,    -1,   219,   158,   227,
     156,   227,   221,    -1,   219,   158,   227,   221,    -1,   219,
     155,   227,   221,    -1,   219,   156,   227,   221,    -1,   219,
     221,    -1,   227,   168,   219,   169,   227,    -1,   227,   168,
     219,    -1,    47,   170,   226,   171,    -1,    49,   170,   226,
     171,    -1,    48,   170,   226,   171,    -1,    23,   170,   226,
     171,   172,   214,    -1,    50,   173,   226,   174,   226,   174,
     226,   175,    -1,    50,   173,   226,   174,   226,   175,    -1,
     226,    -1,   152,   226,    -1,   222,   226,    -1,   222,   152,
      -1,   226,    -1,   152,   226,    -1,   223,   226,    -1,   223,
     152,    -1,   226,   184,   227,   227,    -1,   152,   226,   184,
     227,   227,    -1,   224,   226,   184,   227,   227,    -1,   224,
     152,    -1,   226,   183,   227,   227,    -1,   152,   226,   183,
     227,   227,    -1,   225,   226,   183,   227,   227,    -1,   225,
     152,    -1,   153,    -1,   153,    -1,   154,    -1,   229,    -1,
     228,   229,    -1,   230,    -1,   231,    -1,    68,   226,   152,
      -1,   102,   226,   152,    -1,    77,   226,   152,    -1,    77,
     226,   227,   152,    -1,    77,   226,   227,   227,   152,    -1,
      25,   226,   152,    -1,    25,   226,   227,   152,    -1,    25,
     226,   227,   227,   152,    -1,    25,   226,   227,   227,   227,
     152,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   112,   112,   120,   121,   125,   126,   127,   128,   129,
     130,   131,   132,   133,   134,   135,   139,   141,   146,   149,
     152,   155,   158,   161,   164,   167,   173,   175,   177,   179,
     181,   183,   185,   187,   189,   191,   193,   195,   197,   199,
     201,   203,   205,   207,   210,   213,   216,   219,   222,   225,
     227,   229,   231,   233,   235,   237,   239,   241,   243,   245,
     247,   249,   251,   253,   255,   257,   259,   261,   266,   267,
     268,   269,   270,   271,   272,   273,   274,   275,   276,   277,
     278,   279,   280,   281,   282,   283,   287,   288,   289,   290,
     291,   292,   293,   294,   295,   299,   301,   306,   308,   313,
     315,   317,   319,   324,   325,   329,   331,   336,   339,   342,
     345,   348,   351,   354,   360,   363,   366,   369,   372,   375,
     381,   382,   383,   384,   385,   386,   387,   388,   389,   390,
     391,   392,   393,   394,   395,   396,   397,   398,   399,   403,
     404,   405,   406,   407,   408,   409,   410,   411,   415,   416,
     417,   418,   422,   426,   428,   433,   435,   437,   439,   446,
     447,   448,   452,   453,   457,   459,   464,   466,   468,   470,
     472,   474,   476,   478,   483,   486,   488,   491,   496,   499,
     501,   504,   509,   512,   514,   517,   522,   525,   527,   530,
     535,   537,   539,   541,   546,   547,   550,   552,   555,   560,
     561,   564,   566,   569,   575,   577,   579,   581,   583,   585,
     591,   600,   601,   605,   606,   607,   608,   612,   614,   619,
     621,   626,   628,   630,   632,   634,   636,   641,   642,   643,
     644,   645,   646,   647,   648,   652,   663,   674,   690,   706,
     718,   730,   747,   763,   771,   782,   787,   790,   793,   796,
     799,   802,   805,   808,   811,   814,   817,   820,   823,   826,
     829,   832,   835,   838,   841,   844,   847,   850,   853,   856,
     859,   862,   868,   870,   872,   877,   882,   884,   889,   891,
     893,   895,   900,   902,   904,   906,   911,   913,   915,   917,
     922,   924,   926,   928,   933,   938,   940,   945,   946,   950,
     951,   955,   959,   963,   965,   968,   975,   977,   980,   984
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ACCELERATION", "AERO", "ANALYSIS",
  "ANALYTIC", "ADJOINT", "AREA", "AOATTACK", "AUTOMATIC", "AVEFLAG",
  "CENTRAL", "COEF", "CONSTRAINT", "CONTINUE", "CONVECT", "CTRLCOST",
  "CRITERIA", "CTE", "CTE1", "CTE2", "DAMPING", "DEFOPR", "DENSITY",
  "DENSPROJ", "DISPLEVEL", "DISPLACMENT", "DIRECT", "DSGVAR", "DX", "DY",
  "DZ", "DCMASS", "ELCPOW", "ELCRHS", "ELCVOLT", "EMOD1", "EMOD2", "END",
  "ENERGY", "EQUALITY", "ESTAT", "ESTATCRIT", "FAILPROB", "FORWARD",
  "FREQUENCY", "FUNCRIT", "FUNVAR", "FUNOPR", "GENERATE", "GMOD", "GRDTYP",
  "GRDMTH", "GRDRELINC", "GRDABSINC", "GRDFILTER", "GRDSMI", "INEQUALITY",
  "INIT", "INERTIA", "INTFORCE", "KINETIC", "LAMBDA", "LAYER", "LCASE",
  "MASS", "MMAALG", "MOVINGNODES", "GCMALG", "MMAINT", "MMAFLT", "MOMAER",
  "NLPALG", "NLPINT", "NLPFLT", "NODELE", "NODDENS", "NODPOS", "NUMERIC",
  "NORM", "OBJECTIVE", "OCMALG", "OCMFLT", "OCMINT", "OPCOSINUS",
  "OPKSFUNC", "OPMULTI", "OPSINUS", "OPSUM", "OPEXP", "OPLOG", "OPUDEF",
  "PERMTVY", "PHI", "PMACRITERIA", "POISSON", "PULLIN", "RX", "RY", "RZ",
  "RELINDEX", "RHSDSGNDEP", "RGDVELOC", "SALALG", "SALFLT", "SALINT",
  "SBOOM", "SDTEMP", "SLPALG", "SLPFLT", "SLPINT", "SOLVER", "STCVAR",
  "STRESS", "STRLEVEL", "STRVOL", "S1L", "S1M", "S1U", "S2L", "S2M", "S2U",
  "SXL", "SXM", "SXU", "SYL", "SYM", "SYU", "SZL", "SZM", "SZU", "STCANA",
  "STCDYN", "STCDAM", "STCTRA", "STCELEVAR", "STCVARTYP", "TEMP",
  "THICK_opt", "THPOWER", "TIME", "THERMCOND", "TRADIATE", "VELOCITY",
  "VML", "VMM", "VMU", "VOLTAGE", "XMACH", "YAWANGLE", "YMODUS", "NewLine",
  "IntConstant", "DblConstant", "'+'", "'-'", "'*'", "'^'", "IPOPTALG",
  "IPOPT_INT", "IPOPT_FLT", "SNOPTALG", "SNOPT_INT", "SNOPT_FLT", "'{'",
  "'}'", "':'", "'>'", "'<'", "'['", "']'", "'='", "'('", "','", "')'",
  "$accept", "FinalizedData", "All", "Component", "CriteriaSet",
  "Criteria", "CritTyp", "StrTyp", "DispTyp", "TimeCase", "Objective",
  "ConstraintSet", "ConstraintTyp", "DsgVariableSet", "DsgVariable",
  "StcVariableSet", "AccTyp", "CompTyp", "FldAttrTyp", "ElecTyp",
  "StcAnalysisSet", "StcAnalysisData", "StcDynTyp", "StcDampTyp",
  "SolverSet", "SolverTyp", "NlpData", "OcmData", "SlpData", "MmaData",
  "SalData", "IPOPTData", "SNOPTData", "GradData", "GradTyp", "GradMethod",
  "StcEleVarSet", "StcEleVarInfo", "Function", "FuncTyp", "StatementSet",
  "StatementNumber", "Statement", "OprTyp", "FunctionDef", "Generate",
  "IntegerList", "DirectIntegerList", "DispLevelList", "StrLevelList",
  "Integer", "Float", "PerfFlags", "PerfFlag", "MovingNodes", "RhsDsgnDep",
  "NodalDensInfo", "DensProjInfo", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
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
     405,   406,   407,   408,   409,    43,    45,    42,    94,   410,
     411,   412,   413,   414,   415,   123,   125,    58,    62,    60,
      91,    93,    61,    40,    44,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,   176,   177,   178,   178,   179,   179,   179,   179,   179,
     179,   179,   179,   179,   179,   179,   180,   180,   181,   181,
     181,   181,   181,   181,   181,   181,   182,   182,   182,   182,
     182,   182,   182,   182,   182,   182,   182,   182,   182,   182,
     182,   182,   182,   182,   182,   182,   182,   182,   182,   182,
     182,   182,   182,   182,   182,   182,   182,   182,   182,   182,
     182,   182,   182,   182,   182,   182,   182,   182,   183,   183,
     183,   183,   183,   183,   183,   183,   183,   183,   183,   183,
     183,   183,   183,   183,   183,   183,   184,   184,   184,   184,
     184,   184,   184,   184,   184,   185,   185,   186,   186,   187,
     187,   187,   187,   188,   188,   189,   189,   190,   190,   190,
     190,   190,   190,   190,   191,   191,   191,   191,   191,   191,
     192,   192,   192,   192,   192,   192,   192,   192,   192,   192,
     192,   192,   192,   192,   192,   192,   192,   192,   192,   193,
     193,   193,   193,   193,   193,   193,   193,   193,   194,   194,
     194,   194,   195,   196,   196,   197,   197,   197,   197,   198,
     198,   198,   199,   199,   200,   200,   201,   201,   201,   201,
     201,   201,   201,   201,   202,   202,   202,   202,   203,   203,
     203,   203,   204,   204,   204,   204,   205,   205,   205,   205,
     206,   206,   206,   206,   207,   207,   207,   207,   207,   208,
     208,   208,   208,   208,   209,   209,   209,   209,   209,   209,
     209,   210,   210,   211,   211,   211,   211,   212,   212,   213,
     213,   214,   214,   214,   214,   214,   214,   215,   215,   215,
     215,   215,   215,   215,   215,   216,   216,   216,   216,   216,
     216,   216,   216,   216,   216,   217,   218,   218,   218,   218,
     218,   218,   218,   218,   218,   218,   218,   218,   218,   218,
     218,   218,   218,   218,   218,   218,   218,   218,   218,   218,
     218,   218,   219,   219,   219,   220,   221,   221,   222,   222,
     222,   222,   223,   223,   223,   223,   224,   224,   224,   224,
     225,   225,   225,   225,   226,   227,   227,   228,   228,   229,
     229,   230,   231,   232,   232,   232,   233,   233,   233,   233
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     3,     2,     3,     5,
       4,     6,     4,     5,     6,     7,     1,     4,     1,     4,
       2,     3,     3,     4,     3,     3,     1,     1,     2,     2,
       1,     3,     1,     6,     9,     7,     7,     8,     8,     2,
       2,     2,     2,     2,     2,     2,     1,     3,     2,     2,
       1,     1,     2,     1,     1,     4,     2,     5,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     2,     2,     5,     3,     7,
       5,     6,     4,     1,     1,     3,     2,     7,     6,     4,
       3,     7,     5,     4,     7,     6,     6,     5,     7,     6,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     3,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     1,     3,     2,     3,     3,     4,     7,     1,
       1,     1,     1,     1,     4,     3,     4,     4,     4,     4,
       4,     4,     4,     4,     3,     4,     3,     4,     3,     4,
       3,     4,     3,     4,     3,     4,     3,     4,     3,     4,
       3,     4,     3,     4,     0,     3,     4,     3,     4,     0,
       3,     4,     3,     4,     3,     4,     4,     4,     3,    12,
      16,     1,     1,     1,     1,     1,     1,     3,     2,     4,
       5,     5,     6,     7,     6,     7,     8,     1,     1,     1,
       1,     1,     1,     1,     1,     2,     3,     3,     4,     1,
       2,     2,     3,     1,     2,     2,     7,     7,     5,     5,
       5,     3,     5,     5,     3,     3,     3,     1,     8,     8,
       6,     6,     6,     4,     6,     6,     4,     4,     4,     2,
       5,     3,     4,     4,     4,     6,     8,     6,     1,     2,
       2,     2,     1,     2,     2,     2,     4,     5,     5,     2,
       4,     5,     5,     2,     1,     1,     1,     1,     2,     1,
       1,     3,     3,     3,     4,     5,     3,     4,     5,     6
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned short yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     3,     5,     6,     7,     8,
       9,    10,    12,    11,    13,   297,   299,   300,    15,    14,
       0,     0,   294,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     1,     2,     4,    17,     0,     0,   106,
       0,     0,     0,     0,     0,   154,     0,   218,     0,   298,
       0,    16,   306,   295,   296,     0,   105,   301,   303,     0,
     231,   228,   229,   230,   227,   232,   234,   233,    98,     0,
       0,   302,     0,     0,   153,   217,     0,     0,    42,    37,
       0,     0,    64,    60,     0,     0,    26,     0,     0,     0,
       0,     0,    36,    63,    28,     0,     0,     0,    56,     0,
      40,     0,     0,     0,     0,     0,    61,     0,     0,   103,
     104,     0,     0,     0,   161,   160,     0,   159,   163,     0,
     162,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     165,     0,     0,   307,     0,   304,     0,     0,     0,     0,
     164,     0,     0,    38,     0,     0,     0,     0,    66,    59,
      58,     0,    50,    51,    49,    30,    87,    88,    89,    93,
      90,    91,    92,    94,    54,    86,     0,     0,    39,     0,
      53,    52,    62,     0,     0,     0,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    68,    69,    70,     0,    55,     0,     0,     0,     0,
       0,    18,     0,     0,   102,     0,   110,     0,     0,   149,
       0,   148,   150,     0,     0,   155,   156,     0,     0,     0,
       0,     0,     0,     0,     0,   194,   199,     0,     0,   100,
       0,   308,     0,   305,     0,     0,     0,     0,     0,     0,
     295,     0,     0,   239,   257,   243,     0,     0,    97,     0,
       0,    35,     0,     0,    32,     0,     0,   278,     0,     0,
      41,     0,    57,    31,     0,     0,     0,    34,     0,     0,
      95,    96,     0,    20,     0,    22,     0,   113,   109,     0,
       0,   151,   117,   121,     0,   125,   126,   124,   131,   132,
     133,     0,   127,   123,   134,   135,   136,   130,   128,   129,
     152,   122,     0,     0,   120,   157,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   219,
       0,     0,   309,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   241,   244,   235,   240,     0,     0,     0,   269,
     245,     0,     0,   116,     0,     0,     0,     0,    33,   279,
     281,    65,   280,     0,    27,    29,     0,     0,     0,    19,
       0,     0,     0,    23,   101,   112,     0,     0,     0,   115,
     119,     0,     0,     0,     0,     0,     0,   169,   171,     0,
       0,     0,     0,   166,     0,     0,     0,     0,   167,     0,
       0,     0,     0,   170,     0,     0,     0,     0,   168,     0,
       0,     0,     0,   172,     0,     0,     0,     0,   173,   220,
      99,     0,     0,     0,     0,     0,     0,     0,     0,   221,
     237,   242,   236,   255,   256,   254,   251,   271,   114,   118,
       0,     0,    67,     0,     0,     0,    24,     0,    21,     0,
     108,     0,     0,   138,   146,   147,   143,   140,   141,   142,
     145,   144,   137,   139,     0,   186,   188,   211,   212,     0,
       0,     0,     0,     0,     0,     0,     0,   174,   176,     0,
       0,   180,   178,     0,     0,   192,   190,     0,     0,   184,
     182,     0,     0,   195,   197,     0,     0,   200,   202,     0,
       0,     0,     0,   224,     0,   272,   274,   273,     0,   222,
     238,   267,   268,     0,     0,   266,     0,     0,     0,   263,
       0,     0,     0,     0,     0,     0,     0,    43,     0,    25,
     111,   107,   158,   204,   187,   189,   214,   216,   213,   215,
       0,     0,     0,     0,   208,   175,   177,   181,   179,   193,
     191,   185,   183,   196,   198,   201,   203,     0,   225,     0,
     223,   252,   253,   249,   250,   248,   270,     0,   289,    46,
       0,     0,    45,     0,     0,     0,     0,     0,     0,   277,
     205,   206,   207,     0,   226,   275,   264,   265,   261,   262,
       0,     0,   260,     0,     0,     0,     0,   293,    48,     0,
       0,    47,     0,     0,     0,   246,   247,     0,     0,   286,
       0,     0,     0,    44,   276,     0,   258,   259,   287,   288,
       0,     0,   290,     0,   291,   292,     0,     0,     0,     0,
     282,   283,   285,     0,   284,   209,     0,     0,     0,     0,
     210
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,    13,    14,    15,    16,    46,   118,   204,   174,   212,
      17,    18,   121,    19,    49,    20,   312,   472,   223,   313,
      21,    55,   126,   129,    22,   140,   323,   326,   332,   319,
     329,   335,   338,   397,   479,   550,    23,    57,    78,    79,
     251,   252,   253,   254,   255,   213,   266,   639,   532,   584,
     175,   257,    24,    25,    26,    27,    28,    29
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -245
static const short yypact[] =
{
     579,  -119,  -104,  -101,   -90,  -101,  -101,   -80,  -101,   -70,
     -63,   -59,   -53,   101,   765,  -245,  -101,  -245,  -101,  -101,
    -101,   147,  -101,  -101,   -36,  -245,  -245,  -245,  -245,  -245,
    -101,  -101,  -245,   243,  -101,   -35,   348,   467,   -11,  -101,
    -101,   147,  -101,  -245,  -245,  -245,  -245,   715,     5,  -245,
      47,     6,    14,    -9,  -101,  -245,   403,  -245,  -101,  -245,
       5,  -245,  -245,  -245,  -245,   363,  -245,  -245,  -245,   368,
    -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,   -75,
      -2,  -245,   403,    20,  -245,  -245,  -101,  -101,  -245,  -245,
      -7,  -101,    90,  -245,  -101,  -101,   -16,   -22,  -101,  -101,
     -20,  -101,  -245,  -245,    21,  -101,  -101,  -101,  -245,  -101,
    -245,  -101,  -101,     7,   548,  -101,  -245,  -101,     9,  -245,
    -245,   467,    57,    16,  -245,  -245,    10,  -245,  -245,    86,
    -245,  -101,    25,    27,    40,    44,    56,    61,    63,    70,
    -245,   -84,   467,  -245,   390,  -245,    75,    71,   263,   449,
    -245,    16,   -20,  -245,  -101,  -101,   -20,   144,    77,  -245,
    -245,   144,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,
    -245,  -245,  -245,  -245,  -245,  -245,   -20,   144,  -245,   -20,
    -245,  -245,  -245,   548,  -101,  -101,  -245,  -245,  -245,  -245,
    -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,
    -245,  -245,  -245,  -245,    47,  -245,   -20,  -101,    76,  -101,
      47,  -245,     4,   106,  -245,    96,  -245,   112,   206,  -245,
     -20,  -245,  -245,   449,   587,  -245,  -245,   115,    47,   238,
     238,   264,   288,   296,   297,   275,   278,   118,   140,  -245,
     135,  -245,   163,  -245,   303,   157,   159,   162,   166,   252,
     197,   150,   181,   165,    35,  -245,   216,  -127,  -245,   449,
     587,  -245,    47,    47,  -101,  -101,    67,  -245,   144,   119,
    -245,   121,  -245,  -245,    47,    47,    47,  -245,   -15,  -101,
    -245,  -245,  -101,  -245,   168,  -245,   449,  -245,  -245,   224,
      47,  -245,  -245,  -245,  -101,  -245,  -245,  -245,  -245,  -245,
    -245,  -101,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,
    -245,  -245,   449,   449,  -245,  -245,    47,  -101,    47,   -14,
     -14,  -101,    47,   142,    47,  -101,    -8,    47,  -101,    19,
      47,  -101,    78,  -101,    47,   -26,  -101,    47,    12,  -245,
     239,   449,  -245,   252,   246,  -101,  -101,  -101,  -101,   136,
     259,   181,   261,  -245,   268,  -245,    47,    47,    47,  -245,
    -245,   300,   300,  -245,   449,   449,   257,   293,  -245,  -245,
    -245,  -245,  -245,   153,  -245,  -245,   292,   350,  -101,  -245,
     277,   291,    -5,  -245,  -245,  -245,   209,  -101,   429,  -245,
    -245,  -101,   299,   310,    36,  -101,    47,   441,   441,   317,
     327,  -101,    47,   441,   331,   334,    47,  -101,   441,   339,
     354,    47,  -101,   441,   356,   361,    47,  -101,   441,   366,
     372,  -101,    47,   441,   397,   399,  -101,    47,   441,  -245,
    -245,   203,   408,   392,   395,   396,   398,   404,   409,  -245,
     412,  -245,  -245,   526,   526,     8,    79,   410,  -245,  -245,
     321,   321,  -245,   413,   423,   501,  -245,  -101,  -245,   437,
    -245,   438,   439,  -245,  -245,  -245,  -245,  -245,  -245,  -245,
    -245,  -245,  -245,  -245,   440,  -245,  -245,  -245,  -245,   442,
     444,   446,   195,    47,    47,  -101,   447,  -245,  -245,   450,
     453,  -245,  -245,   455,   458,  -245,  -245,   460,   461,  -245,
    -245,   462,   470,  -245,  -245,   471,   472,  -245,  -245,   473,
     476,   435,   480,  -245,   465,  -245,  -245,  -245,   486,  -245,
    -245,  -245,  -245,    47,    47,  -245,    47,    47,    47,  -245,
      47,  -101,   171,   -20,   201,   351,   351,   474,   336,  -245,
    -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,  -245,
     489,   490,   491,  -101,  -245,  -245,  -245,  -245,  -245,  -245,
    -245,  -245,  -245,  -245,  -245,  -245,  -245,   494,  -245,   449,
    -245,   526,   526,   526,   526,    15,  -245,   -20,  -245,  -245,
     -20,    47,  -245,  -101,   221,   548,   228,   144,  -101,  -245,
    -245,  -245,  -245,    47,  -245,  -245,  -245,  -245,  -245,  -245,
      47,    47,  -245,    47,    47,    47,   548,  -245,  -245,   548,
      47,  -245,   237,   475,    47,   526,   526,    47,    47,  -245,
      47,    47,    47,  -245,  -245,    47,  -245,  -245,  -245,  -245,
      47,    47,  -245,    47,  -245,  -245,   483,   374,  -101,   286,
    -245,  -245,  -245,   393,  -245,  -245,   484,   144,   311,   500,
    -245
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
    -245,  -245,  -245,   639,  -245,   623,  -245,  -175,  -109,  -245,
    -245,  -245,   595,  -245,   624,  -245,   401,  -245,   506,   402,
    -245,   618,  -245,  -245,  -245,   581,  -245,  -245,  -245,   434,
    -245,  -245,  -245,   155,  -245,  -245,  -245,   640,   -81,  -245,
    -163,  -244,  -228,   186,  -217,     0,  -155,  -245,   247,   152,
      -3,   103,  -245,   665,  -245,  -245,  -245,  -245
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -295
static const short yytable[] =
{
      33,   128,    35,    36,   154,    38,   269,   351,   273,   282,
     166,   167,   168,    47,   207,    48,    50,    51,   184,    56,
      58,   162,   271,   352,   354,   219,   394,    60,    47,   124,
     361,    50,     5,    30,   353,   208,    82,    83,   394,    58,
     214,   362,   477,   261,   394,   208,   119,   264,    31,   127,
     130,   131,    32,   237,   208,   141,   395,   396,   208,   208,
     169,   239,    34,   120,   394,   208,     8,   270,   258,    32,
     272,   394,    37,   125,   209,   406,   407,   147,   170,   171,
     172,   344,    39,   152,   153,   208,   349,   155,   156,    40,
     148,   159,   160,    41,   163,   164,   165,   277,   176,    42,
     351,    43,   178,   179,   180,   351,   181,   208,   182,   183,
     185,   291,   205,   373,   206,   478,   352,    67,   173,   220,
     224,   352,   217,   440,   411,   412,   227,   353,   228,   208,
     394,    32,   353,    32,   421,   422,    65,   379,   238,    69,
      80,    81,   292,   123,    32,   256,    32,   458,   260,   161,
     210,   262,   263,   122,   267,   149,   283,   151,   267,   245,
      32,   211,   225,   523,   524,   221,   222,    32,   144,    32,
     600,   601,   146,   245,   267,   426,   427,   229,   363,   230,
     431,   274,   275,   246,   247,   248,   177,   351,   416,   417,
     356,   357,   231,   358,   394,   158,   232,   246,   247,   248,
      63,    64,   546,   352,   278,   384,   280,   547,   233,   216,
      63,    64,   284,   234,   353,   235,   401,   402,   289,   370,
      32,   314,   236,   548,   215,   218,   245,   243,   246,   247,
     248,   389,   390,   371,   526,   527,   244,   528,   226,    32,
     549,   256,   268,    63,    64,   240,   256,   242,   256,   279,
     246,   247,   248,   286,   359,   157,   208,   314,   285,   208,
     430,   368,   369,   372,   287,   267,   372,   315,   372,   245,
     339,   370,    32,   370,    32,   245,   381,   340,   380,   382,
      52,    53,    54,   448,   449,   374,   245,   375,   437,   250,
      64,   387,   341,   246,   247,   248,   265,    32,   388,   246,
     247,   248,   438,   250,    64,   370,    32,   276,   317,   318,
     246,   247,   248,   281,   392,   342,   350,   355,   399,   452,
     383,   290,   405,   578,    32,   410,   245,   345,   415,   346,
     419,   316,   347,   424,    63,    64,   348,   579,   321,   322,
     256,   256,   433,   434,   435,   436,   256,   246,   247,   248,
     246,   247,   248,   578,    32,   511,   250,    64,   288,    63,
      64,   460,    63,    64,  -294,   366,   367,   582,   453,   512,
     372,   324,   325,   607,    32,   455,   385,   376,   377,   378,
     607,    32,   459,   360,   463,   473,   461,   608,   474,   370,
      32,   429,   480,   386,   611,    62,    63,    64,   489,   250,
      64,   327,   328,   623,   494,   250,    64,   330,   331,   498,
     610,   439,   432,   441,   502,   249,   250,    64,   505,   391,
     442,   393,   450,   509,   581,   400,   454,   404,   256,   456,
     409,   620,   612,   414,   621,   333,   334,   420,   642,    32,
     425,   336,   337,   521,   522,   525,   529,   533,   533,   464,
     465,   475,   643,   466,   538,   343,   250,    64,   451,   443,
     444,   445,   476,   370,    32,   457,   467,   468,   603,   487,
     132,   604,   133,   531,    32,   398,   134,   649,   403,   488,
     469,   408,   553,   491,   413,   135,   492,   418,   595,   462,
     423,   495,   648,   428,   482,   483,   484,   485,   486,   481,
      68,    63,    64,   583,    32,   490,   496,   136,   499,   493,
     588,   589,   137,   500,   497,   143,    63,    64,   503,   501,
     145,    63,    64,   470,   504,   506,   638,    32,   577,   580,
     510,   580,   585,   585,    70,    71,    72,    73,    74,    75,
      76,    77,   241,    63,    64,   645,    32,   446,   447,   507,
     593,   508,    70,    71,    72,    73,    74,    75,    76,    77,
     513,   519,   138,   514,   520,   139,   515,   516,   471,   517,
     518,   596,   597,   598,   599,   602,   208,   537,   535,   530,
     606,   609,    32,   609,   267,   613,   551,   552,   536,   539,
     540,   541,   542,     1,   543,   293,   544,     2,   545,   554,
     294,   567,   555,   295,     3,   556,   296,   557,     4,   372,
     558,   297,   559,   560,   561,   626,   627,   298,   299,   300,
      63,    64,   562,   563,   564,   565,   571,   572,   566,   573,
     574,   575,   568,   576,   640,   641,   644,   569,   570,   587,
     646,   590,   591,   592,   267,   372,   594,     5,   637,   647,
     624,   301,   650,    45,    61,   142,     6,   259,    66,    84,
       7,   364,   365,   150,   320,   186,   187,   188,   189,   190,
     191,   192,   193,   194,   195,   196,   197,   198,   199,   200,
     302,     8,    85,   303,   605,   304,   305,   306,   586,    59,
       0,     9,    10,   201,   202,   203,   614,     0,   534,     0,
       0,     0,     0,   615,   616,     0,   617,   618,   619,     0,
       0,    11,     0,   622,     0,    12,     0,   625,    86,    87,
     628,   629,     0,   630,   631,   632,   307,     0,   633,   308,
     309,     0,    88,   634,   635,   310,   636,    89,   311,     0,
      32,    90,    91,     0,     0,     0,     0,     0,    92,    93,
      94,    95,     0,     0,     0,    96,     0,    97,     0,    98,
       0,    99,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   100,   101,   102,   103,     1,
       0,   104,     0,     2,     0,     0,     0,   105,     0,     0,
       3,     0,     0,   106,     4,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    44,     0,     0,     0,     0,     0,
     107,     0,   108,     0,     0,     0,   109,     0,     0,     0,
       0,     0,   110,   111,     0,     0,     0,     0,     0,   112,
     113,   114,     0,     5,     0,     0,     0,     0,     0,     0,
       0,     0,     6,     0,     0,     0,     7,     0,     0,     0,
       0,     0,     0,   115,     0,   116,     0,     0,     0,   117,
       0,     0,     0,     0,     0,     0,     0,     8,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     9,    10,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    11,     0,     0,
       0,    12
};

static const short yycheck[] =
{
       3,    10,     5,     6,    11,     8,   161,   251,   183,     5,
      30,    31,    32,    16,     5,    18,    19,    20,    11,    22,
      23,    43,   177,   251,   252,     9,    52,    30,    31,    15,
     157,    34,    68,   152,   251,    50,    39,    40,    52,    42,
     121,   168,     6,   152,    52,    50,    41,   156,   152,    52,
      53,    54,   153,   137,    50,    58,    70,    71,    50,    50,
      80,   142,   152,    58,    52,    50,   102,   176,   149,   153,
     179,    52,   152,    59,    65,    83,    84,   152,    98,    99,
     100,   244,   152,    86,    87,    50,   249,    90,    91,   152,
     165,    94,    95,   152,    97,    98,    99,   206,   101,   152,
     344,     0,   105,   106,   107,   349,   109,    50,   111,   112,
     113,   220,   115,   268,   117,    79,   344,   152,   138,   103,
     123,   349,   122,   351,   105,   106,   129,   344,   131,    50,
      52,   153,   349,   153,   160,   161,    33,   152,   141,    36,
      37,   152,   223,   137,   153,   148,   153,   152,   151,   165,
     141,   154,   155,    50,   157,   157,   152,   137,   161,    23,
     153,   152,   152,   155,   156,   149,   150,   153,    65,   153,
     155,   156,    69,    23,   177,   163,   164,   152,   259,   152,
     343,   184,   185,    47,    48,    49,   165,   431,   110,   111,
     155,   156,   152,   158,    52,    92,   152,    47,    48,    49,
     153,   154,     7,   431,   207,   286,   209,    12,   152,   152,
     153,   154,   212,   152,   431,   152,    74,    75,   218,   152,
     153,   224,   152,    28,   121,   122,    23,   152,    47,    48,
      49,   312,   313,   166,   155,   156,   165,   158,   152,   153,
      45,   244,   165,   153,   154,   142,   249,   144,   251,   173,
      47,    48,    49,   157,   254,   165,    50,   260,   152,    50,
     341,   264,   265,   266,   152,   268,   269,   152,   271,    23,
     152,   152,   153,   152,   153,    23,   279,   137,   278,   282,
     133,   134,   135,   364,   365,   166,    23,   166,   152,   153,
     154,   294,   157,    47,    48,    49,   152,   153,   301,    47,
      48,    49,   166,   153,   154,   152,   153,   204,    70,    71,
      47,    48,    49,   210,   317,   152,   166,   152,   321,   166,
     152,   218,   325,   152,   153,   328,    23,   170,   331,   170,
     333,   228,   170,   336,   153,   154,   170,   166,    74,    75,
     343,   344,   345,   346,   347,   348,   349,    47,    48,    49,
      47,    48,    49,   152,   153,   152,   153,   154,   152,   153,
     154,   152,   153,   154,   167,   262,   263,   166,    76,   166,
     373,    83,    84,   152,   153,   378,   152,   274,   275,   276,
     152,   153,   382,   167,   387,   388,   386,   166,   391,   152,
     153,   152,   395,   290,   166,   152,   153,   154,   401,   153,
     154,   105,   106,   166,   407,   153,   154,   110,   111,   412,
     585,   152,   166,   152,   417,   152,   153,   154,   421,   316,
     152,   318,   165,   426,   533,   322,    76,   324,   431,   152,
     327,   606,   587,   330,   609,   160,   161,   334,   152,   153,
     337,   163,   164,   443,   444,   445,   446,   450,   451,    20,
      21,   152,   166,    24,   457,   152,   153,   154,   165,   356,
     357,   358,   152,   152,   153,   174,    37,    38,   577,   152,
      67,   580,    69,   152,   153,   320,    73,   166,   323,   152,
      51,   326,   485,   152,   329,    82,   152,   332,   569,   386,
     335,   152,   647,   338,    53,    54,    55,    56,    57,   396,
     152,   153,   154,   152,   153,   402,   152,   104,   152,   406,
     174,   175,   109,   152,   411,   152,   153,   154,   152,   416,
     152,   153,   154,    94,   152,   422,   152,   153,   531,   532,
     427,   534,   535,   536,    85,    86,    87,    88,    89,    90,
      91,    92,   152,   153,   154,   152,   153,   361,   362,   152,
     553,   152,    85,    86,    87,    88,    89,    90,    91,    92,
     152,   152,   159,   171,   152,   162,   171,   171,   139,   171,
     166,   571,   572,   573,   574,   575,    50,    76,   165,   169,
     583,   584,   153,   586,   587,   588,   483,   484,   165,   152,
     152,   152,   152,    14,   152,     8,   152,    18,   152,   152,
      13,   166,   152,    16,    25,   152,    19,   152,    29,   612,
     152,    24,   152,   152,   152,   615,   616,    30,    31,    32,
     153,   154,   152,   152,   152,   152,   523,   524,   152,   526,
     527,   528,   152,   530,   637,   638,   639,   172,   152,   165,
     643,   152,   152,   152,   647,   648,   152,    68,   165,   165,
     175,    64,   152,    14,    31,    60,    77,   151,    34,    41,
      81,   260,   260,    82,   230,   117,   118,   119,   120,   121,
     122,   123,   124,   125,   126,   127,   128,   129,   130,   131,
      93,   102,    42,    96,   581,    98,    99,   100,   536,    24,
      -1,   112,   113,   145,   146,   147,   593,    -1,   451,    -1,
      -1,    -1,    -1,   600,   601,    -1,   603,   604,   605,    -1,
      -1,   132,    -1,   610,    -1,   136,    -1,   614,     3,     4,
     617,   618,    -1,   620,   621,   622,   139,    -1,   625,   142,
     143,    -1,    17,   630,   631,   148,   633,    22,   151,    -1,
     153,    26,    27,    -1,    -1,    -1,    -1,    -1,    33,    34,
      35,    36,    -1,    -1,    -1,    40,    -1,    42,    -1,    44,
      -1,    46,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    60,    61,    62,    63,    14,
      -1,    66,    -1,    18,    -1,    -1,    -1,    72,    -1,    -1,
      25,    -1,    -1,    78,    29,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    39,    -1,    -1,    -1,    -1,    -1,
      95,    -1,    97,    -1,    -1,    -1,   101,    -1,    -1,    -1,
      -1,    -1,   107,   108,    -1,    -1,    -1,    -1,    -1,   114,
     115,   116,    -1,    68,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    77,    -1,    -1,    -1,    81,    -1,    -1,    -1,
      -1,    -1,    -1,   138,    -1,   140,    -1,    -1,    -1,   144,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   102,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   112,   113,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   132,    -1,    -1,
      -1,   136
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    14,    18,    25,    29,    68,    77,    81,   102,   112,
     113,   132,   136,   177,   178,   179,   180,   186,   187,   189,
     191,   196,   200,   212,   228,   229,   230,   231,   232,   233,
     152,   152,   153,   226,   152,   226,   226,   152,   226,   152,
     152,   152,   152,     0,    39,   179,   181,   226,   226,   190,
     226,   226,   133,   134,   135,   197,   226,   213,   226,   229,
     226,   181,   152,   153,   154,   227,   190,   152,   152,   227,
      85,    86,    87,    88,    89,    90,    91,    92,   214,   215,
     227,   152,   226,   226,   197,   213,     3,     4,    17,    22,
      26,    27,    33,    34,    35,    36,    40,    42,    44,    46,
      60,    61,    62,    63,    66,    72,    78,    95,    97,   101,
     107,   108,   114,   115,   116,   138,   140,   144,   182,    41,
      58,   188,   227,   137,    15,    59,   198,   226,    10,   199,
     226,   226,    67,    69,    73,    82,   104,   109,   159,   162,
     201,   226,   188,   152,   227,   152,   227,   152,   165,   157,
     201,   137,   226,   226,    11,   226,   226,   165,   227,   226,
     226,   165,    43,   226,   226,   226,    30,    31,    32,    80,
      98,    99,   100,   138,   184,   226,   226,   165,   226,   226,
     226,   226,   226,   226,    11,   226,   117,   118,   119,   120,
     121,   122,   123,   124,   125,   126,   127,   128,   129,   130,
     131,   145,   146,   147,   183,   226,   226,     5,    50,    65,
     141,   152,   185,   221,   214,   227,   152,   221,   227,     9,
     103,   149,   150,   194,   226,   152,   152,   226,   226,   152,
     152,   152,   152,   152,   152,   152,   152,   137,   226,   214,
     227,   152,   227,   152,   165,    23,    47,    48,    49,   152,
     153,   216,   217,   218,   219,   220,   226,   227,   214,   194,
     226,   184,   226,   226,   184,   152,   222,   226,   165,   222,
     184,   222,   184,   183,   226,   226,   227,   184,   226,   173,
     226,   227,     5,   152,   221,   152,   157,   152,   152,   221,
     227,   184,   214,     8,    13,    16,    19,    24,    30,    31,
      32,    64,    93,    96,    98,    99,   100,   139,   142,   143,
     148,   151,   192,   195,   226,   152,   227,    70,    71,   205,
     205,    74,    75,   202,    83,    84,   203,   105,   106,   206,
     110,   111,   204,   160,   161,   207,   163,   164,   208,   152,
     137,   157,   152,   152,   216,   170,   170,   170,   170,   216,
     166,   217,   218,   220,   218,   152,   155,   156,   158,   221,
     167,   157,   168,   214,   192,   195,   227,   227,   226,   226,
     152,   166,   226,   222,   166,   166,   227,   227,   227,   152,
     221,   226,   226,   152,   214,   152,   227,   226,   226,   214,
     214,   227,   226,   227,    52,    70,    71,   209,   209,   226,
     227,    74,    75,   209,   227,   226,    83,    84,   209,   227,
     226,   105,   106,   209,   227,   226,   110,   111,   209,   226,
     227,   160,   161,   209,   226,   227,   163,   164,   209,   152,
     214,   216,   166,   226,   226,   226,   226,   152,   166,   152,
     218,   152,   152,   227,   227,   227,   219,   219,   214,   214,
     165,   165,   166,    76,    76,   226,   152,   174,   152,   221,
     152,   221,   227,   226,    20,    21,    24,    37,    38,    51,
      94,   139,   193,   226,   226,   152,   152,     6,    79,   210,
     226,   227,    53,    54,    55,    56,    57,   152,   152,   226,
     227,   152,   152,   227,   226,   152,   152,   227,   226,   152,
     152,   227,   226,   152,   152,   226,   227,   152,   152,   226,
     227,   152,   166,   152,   171,   171,   171,   171,   166,   152,
     152,   221,   221,   155,   156,   221,   155,   156,   158,   221,
     169,   152,   224,   226,   224,   165,   165,    76,   226,   152,
     152,   152,   152,   152,   152,   152,     7,    12,    28,    45,
     211,   227,   227,   226,   152,   152,   152,   152,   152,   152,
     152,   152,   152,   152,   152,   152,   152,   166,   152,   172,
     152,   227,   227,   227,   227,   227,   227,   226,   152,   166,
     226,   184,   166,   152,   225,   226,   225,   165,   174,   175,
     152,   152,   152,   226,   152,   214,   221,   221,   221,   221,
     155,   156,   221,   184,   184,   227,   226,   152,   166,   226,
     183,   166,   222,   226,   227,   227,   227,   227,   227,   227,
     183,   183,   227,   166,   175,   227,   221,   221,   227,   227,
     227,   227,   227,   227,   227,   227,   227,   165,   152,   223,
     226,   226,   152,   166,   226,   152,   226,   165,   222,   166,
     152
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <cstddef> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

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
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { 								\
      yyerror ("syntax error: cannot back up");\
      YYERROR;							\
    }								\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)		\
   ((Current).first_line   = (Rhs)[1].first_line,	\
    (Current).first_column = (Rhs)[1].first_column,	\
    (Current).last_line    = (Rhs)[N].last_line,	\
    (Current).last_column  = (Rhs)[N].last_column)
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
#  include <cstdio> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (0)

# define YYDSYMPRINT(Args)			\
do {						\
  if (yydebug)					\
    yysymprint Args;				\
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr, 					\
                  Token, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
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
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if defined (YYMAXDEPTH) && YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  register const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  register char *yyd = yydest;
  register const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
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
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  register int yystate;
  register int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short	yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  register short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  register YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

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
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow ("parser stack overflow",
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	short *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyoverflowlab;
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

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

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
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
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

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
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
#line 113 "Optparser.y"
    { 
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->inputcheck();
	  return 0; 
        ;}
    break;

  case 16:
#line 140 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addCriteria(yyvsp[0].crit); ;}
    break;

  case 17:
#line 142 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addCriteria(yyvsp[0].crit); ;}
    break;

  case 18:
#line 147 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-2].ival; yyval.crit.data=yyvsp[-1].cdata; yyval.crit.t=-1; yyval.crit.anaId=0;
          yyval.crit.igen=0; ;}
    break;

  case 19:
#line 150 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-4].ival; yyval.crit.data=yyvsp[-3].cdata; yyval.crit.t=-1; yyval.crit.anaId=yyvsp[-1].ival-1;
          yyval.crit.igen=0; ;}
    break;

  case 20:
#line 153 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-3].ival; yyval.crit.data=yyvsp[-2].cdata; yyval.crit.t=yyvsp[-1].fval; yyval.crit.anaId=0;
          yyval.crit.igen=0; ;}
    break;

  case 21:
#line 156 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-5].ival; yyval.crit.data=yyvsp[-4].cdata; yyval.crit.t=yyvsp[-3].fval; yyval.crit.anaId=yyvsp[-1].ival-1;
          yyval.crit.igen=0; ;}
    break;

  case 22:
#line 159 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-3].ival; yyval.crit.data=yyvsp[-2].cdata; yyval.crit.t=-1; yyval.crit.anaId=0; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 23:
#line 162 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-4].ival; yyval.crit.data=yyvsp[-3].cdata;   yyval.crit.t=yyvsp[-2].fval; yyval.crit.anaId=0; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 24:
#line 165 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-5].ival; yyval.crit.data=yyvsp[-4].cdata;   yyval.crit.t=-1; yyval.crit.anaId=yyvsp[-2].ival-1; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 25:
#line 168 "Optparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-6].ival; yyval.crit.data=yyvsp[-5].cdata;   yyval.crit.t=yyvsp[-4].fval; yyval.crit.anaId=yyvsp[-2].ival-1; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 26:
#line 174 "Optparser.y"
    { yyval.cdata.typ =  0; ;}
    break;

  case 27:
#line 176 "Optparser.y"
    { yyval.cdata.typ =  0; yyval.cdata.ip = yyvsp[-1].ilist.ip; yyval.cdata.ipsize = yyvsp[-1].ilist.ipsize; ;}
    break;

  case 28:
#line 178 "Optparser.y"
    { yyval.cdata.typ =  1; ;}
    break;

  case 29:
#line 180 "Optparser.y"
    { yyval.cdata.typ =  1; yyval.cdata.ip = yyvsp[-1].ilist.ip; yyval.cdata.ipsize = yyvsp[-1].ilist.ipsize; ;}
    break;

  case 30:
#line 182 "Optparser.y"
    { yyval.cdata.typ =  2; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 31:
#line 184 "Optparser.y"
    { yyval.cdata.typ =  3; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 32:
#line 186 "Optparser.y"
    { yyval.cdata.typ =  4; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; yyval.cdata.i[2]=0; ;}
    break;

  case 33:
#line 188 "Optparser.y"
    { yyval.cdata.typ =  4; yyval.cdata.i[0]=yyvsp[-2].ival-1; yyval.cdata.i[1]=yyvsp[-1].ival; yyval.cdata.i[2]=yyvsp[0].ival-1;;}
    break;

  case 34:
#line 190 "Optparser.y"
    { yyval.cdata.typ =  5; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 35:
#line 192 "Optparser.y"
    { yyval.cdata.typ =  6; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 36:
#line 194 "Optparser.y"
    { yyval.cdata.typ =  7; ;}
    break;

  case 37:
#line 196 "Optparser.y"
    { yyval.cdata.typ =  8; ;}
    break;

  case 38:
#line 198 "Optparser.y"
    { yyval.cdata.typ =  9; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 39:
#line 200 "Optparser.y"
    { yyval.cdata.typ =  9; yyval.cdata.i[0]=yyvsp[0].ival+2; ;}
    break;

  case 40:
#line 202 "Optparser.y"
    { yyval.cdata.typ =  9; yyval.cdata.i[0]=6; ;}
    break;

  case 41:
#line 204 "Optparser.y"
    { yyval.cdata.typ = 10; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 42:
#line 206 "Optparser.y"
    { yyval.cdata.typ = 11; ;}
    break;

  case 43:
#line 208 "Optparser.y"
    { yyval.cdata.typ = 12; yyval.cdata.i[0]=yyvsp[-4].ival; yyval.cdata.d[0]=yyvsp[-3].fval; yyval.cdata.d[1]=yyvsp[-2].fval; yyval.cdata.i[1]=yyvsp[-1].ival; yyval.cdata.i[2]=yyvsp[0].ival; 
          yyval.cdata.ipsize=0; ;}
    break;

  case 44:
#line 211 "Optparser.y"
    { yyval.cdata.typ = 12; yyval.cdata.i[0]=yyvsp[-7].ival; yyval.cdata.d[0]=yyvsp[-6].fval; yyval.cdata.d[1]=yyvsp[-5].fval; yyval.cdata.i[1]=yyvsp[-4].ival; yyval.cdata.i[2]=yyvsp[-3].ival;
          yyval.cdata.ip  =yyvsp[-1].ilist.ip; yyval.cdata.ipsize= yyvsp[-1].ilist.ipsize;;}
    break;

  case 45:
#line 214 "Optparser.y"
    { yyval.cdata.typ =13; yyval.cdata.i[0]=yyvsp[-5].ival; yyval.cdata.i[1]=yyvsp[-4].ival; yyval.cdata.d[0]=yyvsp[-3].fval; 
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 46:
#line 217 "Optparser.y"
    { yyval.cdata.typ =13; yyval.cdata.i[0]=1; yyval.cdata.i[1]=yyvsp[-4].ival; yyval.cdata.d[0]=yyvsp[-3].fval;
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 47:
#line 220 "Optparser.y"
    { yyval.cdata.typ =20; yyval.cdata.i[0]=yyvsp[-6].ival; yyval.cdata.i[1]=yyvsp[-5].ival; yyval.cdata.d[0]=yyvsp[-4].fval; yyval.cdata.i[2]=yyvsp[-3].ival;
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 48:
#line 223 "Optparser.y"
    { yyval.cdata.typ =20; yyval.cdata.i[0]=1; yyval.cdata.i[1]=yyvsp[-5].ival; yyval.cdata.d[0]=yyvsp[-4].fval; yyval.cdata.i[2]=yyvsp[-3].ival;
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 49:
#line 226 "Optparser.y"
    { yyval.cdata.typ =  15; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 50:
#line 228 "Optparser.y"
    { yyval.cdata.typ =  16; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 51:
#line 230 "Optparser.y"
    { yyval.cdata.typ =  16; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 52:
#line 232 "Optparser.y"
    { yyval.cdata.typ =  17; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 53:
#line 234 "Optparser.y"
    { yyval.cdata.typ =  18; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 54:
#line 236 "Optparser.y"
    { yyval.cdata.typ =  19; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 55:
#line 238 "Optparser.y"
    { yyval.cdata.typ =  21; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 56:
#line 240 "Optparser.y"
    { yyval.cdata.typ = 22; ;}
    break;

  case 57:
#line 242 "Optparser.y"
    { yyval.cdata.typ = 23; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 58:
#line 244 "Optparser.y"
    { yyval.cdata.typ =  24; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 59:
#line 246 "Optparser.y"
    { yyval.cdata.typ =  25; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 60:
#line 248 "Optparser.y"
    { yyval.cdata.typ =  26; ;}
    break;

  case 61:
#line 250 "Optparser.y"
    { yyval.cdata.typ = 27;  ;}
    break;

  case 62:
#line 252 "Optparser.y"
    { yyval.cdata.typ = 28; yyval.cdata.i[0] = yyvsp[0].ival-1; ;}
    break;

  case 63:
#line 254 "Optparser.y"
    { yyval.cdata.typ = 29;  ;}
    break;

  case 64:
#line 256 "Optparser.y"
    { yyval.cdata.typ = 30; yyval.cdata.d[0]=2.0; ;}
    break;

  case 65:
#line 258 "Optparser.y"
    { yyval.cdata.typ = 30; yyval.cdata.d[0]=2.0; yyval.cdata.ip = yyvsp[-1].ilist.ip; yyval.cdata.ipsize = yyvsp[-1].ilist.ipsize; ;}
    break;

  case 66:
#line 260 "Optparser.y"
    { yyval.cdata.typ = 30; yyval.cdata.d[0]=yyvsp[0].fval; ;}
    break;

  case 67:
#line 262 "Optparser.y"
    { yyval.cdata.typ = 30; yyval.cdata.d[0]=yyvsp[-3].fval; yyval.cdata.ip = yyvsp[-1].ilist.ip; yyval.cdata.ipsize = yyvsp[-1].ilist.ipsize;;}
    break;

  case 68:
#line 266 "Optparser.y"
    { yyval.ival = 0;  ;}
    break;

  case 69:
#line 267 "Optparser.y"
    { yyval.ival = 1;  ;}
    break;

  case 70:
#line 268 "Optparser.y"
    { yyval.ival = 2;  ;}
    break;

  case 71:
#line 269 "Optparser.y"
    { yyval.ival = 10; ;}
    break;

  case 72:
#line 270 "Optparser.y"
    { yyval.ival = 11; ;}
    break;

  case 73:
#line 271 "Optparser.y"
    { yyval.ival = 12; ;}
    break;

  case 74:
#line 272 "Optparser.y"
    { yyval.ival = 20; ;}
    break;

  case 75:
#line 273 "Optparser.y"
    { yyval.ival = 21; ;}
    break;

  case 76:
#line 274 "Optparser.y"
    { yyval.ival = 22; ;}
    break;

  case 77:
#line 275 "Optparser.y"
    { yyval.ival = 30; ;}
    break;

  case 78:
#line 276 "Optparser.y"
    { yyval.ival = 31; ;}
    break;

  case 79:
#line 277 "Optparser.y"
    { yyval.ival = 32; ;}
    break;

  case 80:
#line 278 "Optparser.y"
    { yyval.ival = 40; ;}
    break;

  case 81:
#line 279 "Optparser.y"
    { yyval.ival = 41; ;}
    break;

  case 82:
#line 280 "Optparser.y"
    { yyval.ival = 42; ;}
    break;

  case 83:
#line 281 "Optparser.y"
    { yyval.ival = 50; ;}
    break;

  case 84:
#line 282 "Optparser.y"
    { yyval.ival = 51; ;}
    break;

  case 85:
#line 283 "Optparser.y"
    { yyval.ival = 52; ;}
    break;

  case 86:
#line 287 "Optparser.y"
    { yyval.ival = yyvsp[0].ival-1; ;}
    break;

  case 87:
#line 288 "Optparser.y"
    { yyval.ival = 0;    ;}
    break;

  case 88:
#line 289 "Optparser.y"
    { yyval.ival = 1;    ;}
    break;

  case 89:
#line 290 "Optparser.y"
    { yyval.ival = 2;    ;}
    break;

  case 90:
#line 291 "Optparser.y"
    { yyval.ival = 3;    ;}
    break;

  case 91:
#line 292 "Optparser.y"
    { yyval.ival = 4;    ;}
    break;

  case 92:
#line 293 "Optparser.y"
    { yyval.ival = 5;    ;}
    break;

  case 93:
#line 294 "Optparser.y"
    { yyval.ival = 6;    ;}
    break;

  case 94:
#line 295 "Optparser.y"
    { yyval.ival = 7;    ;}
    break;

  case 95:
#line 300 "Optparser.y"
    { yyval.fval = yyvsp[0].ival-1; ;}
    break;

  case 96:
#line 302 "Optparser.y"
    { yyval.fval = yyvsp[0].fval-1; ;}
    break;

  case 97:
#line 307 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addObjective(yyvsp[-2].fval,yyvsp[0].fall ); ;}
    break;

  case 98:
#line 309 "Optparser.y"
    { double one=1.0; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addObjective(one,yyvsp[0].fall ); ;}
    break;

  case 99:
#line 314 "Optparser.y"
    { int num=yyvsp[-4].ival; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addConstraint(num,yyvsp[-3].ival,yyvsp[-2].fval,yyvsp[0].fall); ;}
    break;

  case 100:
#line 316 "Optparser.y"
    { int num=yyvsp[-2].ival; double one=1.0; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addConstraint(num,yyvsp[-1].ival,one,yyvsp[0].fall); ;}
    break;

  case 101:
#line 318 "Optparser.y"
    { int num=yyvsp[-4].ival; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addConstraint(num,yyvsp[-3].ival,yyvsp[-2].fval,yyvsp[0].fall);;}
    break;

  case 102:
#line 320 "Optparser.y"
    { int num=yyvsp[-2].ival; double one=1.0; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addConstraint(num,yyvsp[-1].ival,one,yyvsp[0].fall);;}
    break;

  case 103:
#line 324 "Optparser.y"
    { yyval.ival=0; ;}
    break;

  case 104:
#line 325 "Optparser.y"
    { yyval.ival=1; ;}
    break;

  case 105:
#line 330 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addDsgvar(yyvsp[0].absvar); ;}
    break;

  case 106:
#line 332 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addDsgvar(yyvsp[0].absvar); ;}
    break;

  case 107:
#line 337 "Optparser.y"
    { yyval.absvar.num=yyvsp[-6].ival; yyval.absvar.val=yyvsp[-5].fval ; yyval.absvar.scl=yyvsp[-4].fval; yyval.absvar.low=yyvsp[-3].fval; yyval.absvar.upp=yyvsp[-2].fval; yyval.absvar.vin=yyvsp[-1].fval;
	  yyval.absvar.igen=0; ;}
    break;

  case 108:
#line 340 "Optparser.y"
    { yyval.absvar.num=yyvsp[-5].ival; yyval.absvar.val=yyvsp[-4].fval ; yyval.absvar.scl=yyvsp[-3].fval; yyval.absvar.low=yyvsp[-2].fval; yyval.absvar.upp=yyvsp[-1].fval; yyval.absvar.vin=0;
	  yyval.absvar.igen=0; ;}
    break;

  case 109:
#line 343 "Optparser.y"
    { yyval.absvar.num=yyvsp[-3].ival; yyval.absvar.val=yyvsp[-2].fval ; yyval.absvar.scl=yyvsp[-1].fval; yyval.absvar.low=-1.0e10; yyval.absvar.upp=1.0e10; 
	  yyval.absvar.vin=0 ; yyval.absvar.igen=0; ;}
    break;

  case 110:
#line 346 "Optparser.y"
    { yyval.absvar.num=yyvsp[-2].ival; yyval.absvar.val=yyvsp[-1].fval ; yyval.absvar.scl=1.0; yyval.absvar.low=-1.0e10; yyval.absvar.upp=1.0e10; 
	  yyval.absvar.vin=0 ; yyval.absvar.igen=0; ;}
    break;

  case 111:
#line 349 "Optparser.y"
    { yyval.absvar.num=yyvsp[-6].ival; yyval.absvar.val=yyvsp[-5].fval ; yyval.absvar.scl=yyvsp[-4].fval; yyval.absvar.low=yyvsp[-3].fval; yyval.absvar.upp=yyvsp[-2].fval; 
	  yyval.absvar.vin=0 ; yyval.absvar.igen=1; yyval.absvar.gena=yyvsp[-1].gen.a; yyval.absvar.genb=yyvsp[-1].gen.b ; yyval.absvar.gens=yyvsp[-1].gen.s; ;}
    break;

  case 112:
#line 352 "Optparser.y"
    { yyval.absvar.num=yyvsp[-4].ival; yyval.absvar.val=yyvsp[-3].fval ; yyval.absvar.scl=yyvsp[-2].fval; yyval.absvar.low=-1.0e10; yyval.absvar.upp=1.0e10; 
	  yyval.absvar.vin=0 ; yyval.absvar.igen=1; yyval.absvar.gena=yyvsp[-1].gen.a; yyval.absvar.genb=yyvsp[-1].gen.b ; yyval.absvar.gens=yyvsp[-1].gen.s; ;}
    break;

  case 113:
#line 355 "Optparser.y"
    { yyval.absvar.num=yyvsp[-3].ival; yyval.absvar.val=yyvsp[-2].fval ; yyval.absvar.scl=1.0; yyval.absvar.low=-1.0e10; yyval.absvar.upp=1.0e10; 
	  yyval.absvar.vin=0 ; yyval.absvar.igen=1; yyval.absvar.gena=yyvsp[-1].gen.a; yyval.absvar.genb=yyvsp[-1].gen.b ; yyval.absvar.gens=yyvsp[-1].gen.s; ;}
    break;

  case 114:
#line 361 "Optparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall); ;}
    break;

  case 115:
#line 364 "Optparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall); ;}
    break;

  case 116:
#line 367 "Optparser.y"
    { int num=yyvsp[-3].ival; int loc1=yyvsp[-1].ival-1; int loc2=0;
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addStcvar(num,yyvsp[-2].ival,loc1,loc2,yyvsp[0].fall); ;}
    break;

  case 117:
#line 370 "Optparser.y"
    { int num=yyvsp[-3].ival; int loc1=yyvsp[-1].ival-1; int loc2=0;
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addStcvar(num,yyvsp[-2].ival,loc1,loc2,yyvsp[0].fall); ;}
    break;

  case 118:
#line 373 "Optparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall); ;}
    break;

  case 119:
#line 376 "Optparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/optpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall); ;}
    break;

  case 120:
#line 381 "Optparser.y"
    { yyval.ival = yyvsp[0].ival-1; ;}
    break;

  case 121:
#line 382 "Optparser.y"
    { yyval.ival = 0; ;}
    break;

  case 122:
#line 383 "Optparser.y"
    { yyval.ival = 1; ;}
    break;

  case 123:
#line 384 "Optparser.y"
    { yyval.ival = 2; ;}
    break;

  case 124:
#line 385 "Optparser.y"
    { yyval.ival = 3; ;}
    break;

  case 125:
#line 386 "Optparser.y"
    { yyval.ival = 4; ;}
    break;

  case 126:
#line 387 "Optparser.y"
    { yyval.ival =10; ;}
    break;

  case 127:
#line 388 "Optparser.y"
    { yyval.ival = 5; ;}
    break;

  case 128:
#line 389 "Optparser.y"
    { yyval.ival = 5; ;}
    break;

  case 129:
#line 390 "Optparser.y"
    { yyval.ival = 2; ;}
    break;

  case 130:
#line 391 "Optparser.y"
    { yyval.ival = 6; ;}
    break;

  case 131:
#line 392 "Optparser.y"
    { yyval.ival = 0; ;}
    break;

  case 132:
#line 393 "Optparser.y"
    { yyval.ival = 1; ;}
    break;

  case 133:
#line 394 "Optparser.y"
    { yyval.ival = 2; ;}
    break;

  case 134:
#line 395 "Optparser.y"
    { yyval.ival = 3; ;}
    break;

  case 135:
#line 396 "Optparser.y"
    { yyval.ival = 4; ;}
    break;

  case 136:
#line 397 "Optparser.y"
    { yyval.ival = 5; ;}
    break;

  case 137:
#line 398 "Optparser.y"
    { yyval.ival = 100000 * yyvsp[-1].ival + yyvsp[0].ival; ;}
    break;

  case 138:
#line 399 "Optparser.y"
    { yyval.ival = -1-(6*(yyvsp[-1].ival-1) + (yyvsp[0].ival-1)); ;}
    break;

  case 139:
#line 403 "Optparser.y"
    { yyval.ival = yyvsp[0].ival-1; ;}
    break;

  case 140:
#line 404 "Optparser.y"
    { yyval.ival = 1; ;}
    break;

  case 141:
#line 405 "Optparser.y"
    { yyval.ival = 2; ;}
    break;

  case 142:
#line 406 "Optparser.y"
    { yyval.ival = 4; ;}
    break;

  case 143:
#line 407 "Optparser.y"
    { yyval.ival = 7; ;}
    break;

  case 144:
#line 408 "Optparser.y"
    { yyval.ival = 8; ;}
    break;

  case 145:
#line 409 "Optparser.y"
    { yyval.ival = 9; ;}
    break;

  case 146:
#line 410 "Optparser.y"
    { yyval.ival = 10; ;}
    break;

  case 147:
#line 411 "Optparser.y"
    { yyval.ival = 11; ;}
    break;

  case 148:
#line 415 "Optparser.y"
    { yyval.ival = 1; ;}
    break;

  case 149:
#line 416 "Optparser.y"
    { yyval.ival = 2; ;}
    break;

  case 150:
#line 417 "Optparser.y"
    { yyval.ival = 3; ;}
    break;

  case 151:
#line 418 "Optparser.y"
    { yyval.ival = 4 + yyvsp[0].ival; ;}
    break;

  case 152:
#line 422 "Optparser.y"
    { yyval.ival = 100; ;}
    break;

  case 153:
#line 427 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structopt->addAnalysisData(yyvsp[0].anadat);  ;}
    break;

  case 154:
#line 429 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structopt->addAnalysisData(yyvsp[0].anadat); ;}
    break;

  case 155:
#line 434 "Optparser.y"
    { yyval.anadat.typ=0; yyval.anadat.ival[0]=yyvsp[-1].ival; ;}
    break;

  case 156:
#line 436 "Optparser.y"
    { yyval.anadat.typ=1; yyval.anadat.ival[0]=yyvsp[-1].ival; ;}
    break;

  case 157:
#line 438 "Optparser.y"
    { yyval.anadat.typ=1; yyval.anadat.ival[0]=yyvsp[-2].ival; yyval.anadat.ival[1]=yyvsp[-1].ival; ;}
    break;

  case 158:
#line 440 "Optparser.y"
    { yyval.anadat.typ=2; yyval.anadat.ival[0]=yyvsp[-5].ival; yyval.anadat.ival[1]=yyvsp[-4].ival; 
	            yyval.anadat.rval[0]=yyvsp[-3].fval; yyval.anadat.rval[1]=yyvsp[-2].fval;
                    yyval.anadat.ival[2]=yyvsp[-1].ival; ;}
    break;

  case 159:
#line 446 "Optparser.y"
    { yyval.ival=yyvsp[0].ival; ;}
    break;

  case 160:
#line 447 "Optparser.y"
    { yyval.ival=0;  ;}
    break;

  case 161:
#line 448 "Optparser.y"
    { yyval.ival=1;  ;}
    break;

  case 162:
#line 452 "Optparser.y"
    { yyval.ival=yyvsp[0].ival; ;}
    break;

  case 163:
#line 453 "Optparser.y"
    { yyval.ival=-1; ;}
    break;

  case 164:
#line 458 "Optparser.y"
    { int num=yyvsp[-1].ival; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->optsol->addSolver(num,yyvsp[0].solv.typ,*(yyvsp[0].solv.param),*(yyvsp[0].solv.pgrad)); ;}
    break;

  case 165:
#line 460 "Optparser.y"
    { int num=yyvsp[-1].ival; /*dynamic_cast<Domain_opt*>(domain)->*/optpro->optsol->addSolver(num,yyvsp[0].solv.typ,*(yyvsp[0].solv.param),*(yyvsp[0].solv.pgrad)); ;}
    break;

  case 166:
#line 465 "Optparser.y"
    { yyval.solv.typ=0; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 167:
#line 467 "Optparser.y"
    { yyval.solv.typ=1; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 168:
#line 469 "Optparser.y"
    { yyval.solv.typ=2; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 169:
#line 471 "Optparser.y"
    { yyval.solv.typ=3; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 170:
#line 473 "Optparser.y"
    { yyval.solv.typ=4; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 171:
#line 475 "Optparser.y"
    { yyval.solv.typ=5; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 172:
#line 477 "Optparser.y"
    { yyval.solv.typ=6; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 173:
#line 479 "Optparser.y"
    { yyval.solv.typ=7; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 174:
#line 484 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 175:
#line 487 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 176:
#line 489 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 177:
#line 492 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 178:
#line 497 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 179:
#line 500 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 180:
#line 502 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 181:
#line 505 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 182:
#line 510 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 183:
#line 513 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 184:
#line 515 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 185:
#line 518 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 186:
#line 523 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 187:
#line 526 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 188:
#line 528 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 189:
#line 531 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 190:
#line 536 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 191:
#line 538 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 192:
#line 540 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 193:
#line 542 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 194:
#line 546 "Optparser.y"
    { yyval.nlp.initialize(); ;}
    break;

  case 195:
#line 548 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 196:
#line 551 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 197:
#line 553 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 198:
#line 556 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 199:
#line 560 "Optparser.y"
    { yyval.nlp.initialize(); ;}
    break;

  case 200:
#line 562 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 201:
#line 565 "Optparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 202:
#line 567 "Optparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 203:
#line 570 "Optparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 204:
#line 576 "Optparser.y"
    { yyval.grad.typ=yyvsp[-1].ival; yyval.grad.filter=0; ;}
    break;

  case 205:
#line 578 "Optparser.y"
    { yyval.grad.mth=yyvsp[-1].ival; ;}
    break;

  case 206:
#line 580 "Optparser.y"
    { yyval.grad.epstyp=0 ; yyval.grad.epsval=yyvsp[-1].fval; ;}
    break;

  case 207:
#line 582 "Optparser.y"
    { yyval.grad.epstyp=1 ; yyval.grad.epsval=yyvsp[-1].fval; ;}
    break;

  case 208:
#line 584 "Optparser.y"
    { dynamic_cast<Domain_opt*>(domain)->setSemiSAflag(1); ;}
    break;

  case 209:
#line 586 "Optparser.y"
    { yyval.grad.filter=1;  yyval.grad.filterTyp=yyvsp[-9].ival; yyval.grad.filterScale=yyvsp[-8].ival;  
          yyval.grad.radius=yyvsp[-7].fval; yyval.grad.maxCount=yyvsp[-6].fval;  yyval.grad.minExp=yyvsp[-5].fval;  yyval.grad.maxExp=yyvsp[-4].fval;
          yyval.grad.numFilCrit=yyvsp[-2].ilist.ipsize;      yyval.grad.filcritList=yyvsp[-2].ilist.ip; 
	  yyval.grad.numGroups=0;              
	  yyval.grad.numFilGrps=0;               yyval.grad.filGrpsList=0; ;}
    break;

  case 210:
#line 592 "Optparser.y"
    { yyval.grad.filter=1;  yyval.grad.filterTyp=yyvsp[-13].ival; yyval.grad.filterScale=yyvsp[-12].ival;  
          yyval.grad.radius=yyvsp[-11].fval; yyval.grad.maxCount=yyvsp[-10].fval;  yyval.grad.minExp=yyvsp[-9].fval;  yyval.grad.maxExp=yyvsp[-8].fval;
          yyval.grad.numFilCrit=yyvsp[-6].ilist.ipsize;      yyval.grad.filcritList=yyvsp[-6].ilist.ip; 
	  yyval.grad.numGroups=yyvsp[-4].ival;              
	  yyval.grad.numFilGrps=yyvsp[-2].ilist.ipsize;      yyval.grad.filGrpsList=yyvsp[-2].ilist.ip; ;}
    break;

  case 211:
#line 600 "Optparser.y"
    { yyval.ival=0; ;}
    break;

  case 212:
#line 601 "Optparser.y"
    { yyval.ival=1; ;}
    break;

  case 213:
#line 605 "Optparser.y"
    { yyval.ival=0; ;}
    break;

  case 214:
#line 606 "Optparser.y"
    { yyval.ival=1; ;}
    break;

  case 215:
#line 607 "Optparser.y"
    { yyval.ival=2; ;}
    break;

  case 216:
#line 608 "Optparser.y"
    { yyval.ival=3; ;}
    break;

  case 217:
#line 613 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structopt->addOptInf(yyvsp[0].sevdata.var,yyvsp[0].sevdata.ele1,yyvsp[0].sevdata.ele2,yyvsp[0].sevdata.typ); ;}
    break;

  case 218:
#line 615 "Optparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structopt->addOptInf(yyvsp[0].sevdata.var,yyvsp[0].sevdata.ele1,yyvsp[0].sevdata.ele2,yyvsp[0].sevdata.typ); ;}
    break;

  case 219:
#line 620 "Optparser.y"
    { yyval.sevdata.var=yyvsp[-3].ival-1; yyval.sevdata.ele1=yyvsp[-2].ival-1; yyval.sevdata.ele2=-1; yyval.sevdata.typ=yyvsp[-1].ival; ;}
    break;

  case 220:
#line 622 "Optparser.y"
    { yyval.sevdata.var=yyvsp[-4].ival-1; yyval.sevdata.ele1=yyvsp[-3].ival-1; yyval.sevdata.ele2=yyvsp[-2].ival-1; yyval.sevdata.typ=yyvsp[-1].ival; ;}
    break;

  case 221:
#line 627 "Optparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-4].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 222:
#line 629 "Optparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-5].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 223:
#line 631 "Optparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-6].ival; yyval.fall.fdata=yyvsp[-3].func; ;}
    break;

  case 224:
#line 633 "Optparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-5].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 225:
#line 635 "Optparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-6].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 226:
#line 637 "Optparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-7].ival; yyval.fall.fdata=yyvsp[-3].func; ;}
    break;

  case 227:
#line 641 "Optparser.y"
    { yyval.ival=0; ;}
    break;

  case 228:
#line 642 "Optparser.y"
    { yyval.ival=1; ;}
    break;

  case 229:
#line 643 "Optparser.y"
    { yyval.ival=2; ;}
    break;

  case 230:
#line 644 "Optparser.y"
    { yyval.ival=3; ;}
    break;

  case 231:
#line 645 "Optparser.y"
    { yyval.ival=4; ;}
    break;

  case 232:
#line 646 "Optparser.y"
    { yyval.ival=5; ;}
    break;

  case 233:
#line 647 "Optparser.y"
    { yyval.ival=6; ;}
    break;

  case 234:
#line 648 "Optparser.y"
    { yyval.ival=7; ;}
    break;

  case 235:
#line 653 "Optparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=1; yyval.func->numgen=0; 
	  yyval.func->numstat[0]=yyvsp[-1].ival-1;  
	  yyval.func->oprtyp[0]=yyvsp[0].sum.oprtyp; yyval.func->oprnum[0]=yyvsp[0].sum.oprnum;
	  yyval.func->a[0]=yyvsp[0].sum.a; yyval.func->p[0]=yyvsp[0].sum.p; yyval.func->b[0]=yyvsp[0].sum.b; 
	  if ( yyvsp[0].sum.igen ) {
	    int ngen=yyval.func->numgen;
	    yyval.func->gena[ngen]=yyvsp[0].sum.gena; yyval.func->genb[ngen]=yyvsp[0].sum.genb;
	    yyval.func->gens[ngen]=yyvsp[0].sum.gens; yyval.func->numgen++;
	  }
	;}
    break;

  case 236:
#line 664 "Optparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=1; yyval.func->numgen=0; 
	  yyval.func->numstat[0]=yyvsp[-2].ival-1;  
	  yyval.func->oprtyp[0]=yyvsp[-1].sum.oprtyp; yyval.func->oprnum[0]=yyvsp[-1].sum.oprnum;
	  yyval.func->a[0]=yyvsp[-1].sum.a; yyval.func->p[0]=yyvsp[-1].sum.p; yyval.func->b[0]=yyvsp[-1].sum.b; 
	  if ( yyvsp[-1].sum.igen ) {
	    int ngen=yyval.func->numgen;
	    yyval.func->gena[ngen]=yyvsp[-1].sum.gena; yyval.func->genb[ngen]=yyvsp[-1].sum.genb;
	    yyval.func->gens[ngen]=yyvsp[-1].sum.gens; yyval.func->numgen++;
	  }
	;}
    break;

  case 237:
#line 675 "Optparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	  yyval.func->numstat[num]=yyvsp[-1].ival-1;  
	  yyval.func->oprtyp[num]=yyvsp[0].sum.oprtyp; yyval.func->oprnum[num]=yyvsp[0].sum.oprnum; 
	  yyval.func->a[num]=yyvsp[0].sum.a; yyval.func->p[num]=yyvsp[0].sum.p; yyval.func->b[num]=yyvsp[0].sum.b;
	  if ( yyvsp[0].sum.igen ) {
	    int ngen=yyval.func->numgen;
	    yyval.func->gena[ngen]=yyvsp[0].sum.gena; yyval.func->genb[ngen]=yyvsp[0].sum.genb;
	    yyval.func->gens[ngen]=yyvsp[0].sum.gens; yyval.func->numgen++;
	  }
	;}
    break;

  case 238:
#line 691 "Optparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	  yyval.func->numstat[num]=yyvsp[-2].ival-1;  
	  yyval.func->oprtyp[num]=yyvsp[-1].sum.oprtyp; yyval.func->oprnum[num]=yyvsp[-1].sum.oprnum; 
	  yyval.func->a[num]=yyvsp[-1].sum.a; yyval.func->p[num]=yyvsp[-1].sum.p; yyval.func->b[num]=yyvsp[-1].sum.b;
	  if ( yyvsp[-1].sum.igen ) {
	    int ngen=yyval.func->numgen;
	    yyval.func->gena[ngen]=yyvsp[-1].sum.gena; yyval.func->genb[ngen]=yyvsp[-1].sum.genb;
	    yyval.func->gens[ngen]=yyvsp[-1].sum.gens; yyval.func->numgen++;
	  }
	;}
    break;

  case 239:
#line 707 "Optparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=1; yyval.func->numgen=0; 
	  yyval.func->numstat[0]=-1;  
	  yyval.func->oprtyp[0]=yyvsp[0].sum.oprtyp; yyval.func->oprnum[0]=yyvsp[0].sum.oprnum;
	  yyval.func->a[0]=yyvsp[0].sum.a; yyval.func->p[0]=yyvsp[0].sum.p; yyval.func->b[0]=yyvsp[0].sum.b; 
	  if ( yyvsp[0].sum.igen ) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"no number for summands assigned\n");
	    fprintf (stderr," *** This is required for generating summands\n");
	    return -1;
	  }
	;}
    break;

  case 240:
#line 719 "Optparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=1; yyval.func->numgen=0; 
	  yyval.func->numstat[0]=-1;  
	  yyval.func->oprtyp[0]=yyvsp[-1].sum.oprtyp; yyval.func->oprnum[0]=yyvsp[-1].sum.oprnum;
	  yyval.func->a[0]=yyvsp[-1].sum.a; yyval.func->p[0]=yyvsp[-1].sum.p; yyval.func->b[0]=yyvsp[-1].sum.b; 
	  if ( yyvsp[-1].sum.igen ) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"no number for summands assigned\n");
	    fprintf (stderr," *** This is required for generating summands!\n");
	    return -1;
	  }
	;}
    break;

  case 241:
#line 731 "Optparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	  yyval.func->numstat[num]=-1;  
	  yyval.func->oprtyp[num]=yyvsp[0].sum.oprtyp; yyval.func->oprnum[num]=yyvsp[0].sum.oprnum; 
	  yyval.func->a[num]=yyvsp[0].sum.a; yyval.func->p[num]=yyvsp[0].sum.p; yyval.func->b[num]=yyvsp[0].sum.b;
	  if ( yyvsp[0].sum.igen ) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"no number for summands assigned\n");
	    fprintf (stderr," *** This is required for generating summands!\n");
	    return -1;
	  }
	;}
    break;

  case 242:
#line 748 "Optparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	  yyval.func->numstat[num]=-1;  
	  yyval.func->oprtyp[num]=yyvsp[-1].sum.oprtyp; yyval.func->oprnum[num]=yyvsp[-1].sum.oprnum; 
	  yyval.func->a[num]=yyvsp[-1].sum.a; yyval.func->p[num]=yyvsp[-1].sum.p; yyval.func->b[num]=yyvsp[-1].sum.b;
	  if ( yyvsp[-1].sum.igen ) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	;}
    break;

  case 243:
#line 764 "Optparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=0; yyval.func->numgen=0; int num=yyvsp[0].fdef.num;
	  if (num > MAXOPR) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of subfunctions > MAXOPR");
	    return -1;
	  }
	  yyval.func->subfunc[num]=yyvsp[0].fdef.falldata;;}
    break;

  case 244:
#line 772 "Optparser.y"
    { int num=yyvsp[0].fdef.num;
	  if (num > MAXOPR) {
	    fprintf (stderr," *** Error reading optimization input file: ");
	    fprintf (stderr,"number of subfunctions > MAXOPR");
	    return -1;
	  }
	  yyval.func->subfunc[num]=yyvsp[0].fdef.falldata; ;}
    break;

  case 245:
#line 783 "Optparser.y"
    { yyval.ival = yyvsp[-1].ival; ;}
    break;

  case 246:
#line 788 "Optparser.y"
    { yyval.sum.a=yyvsp[-6].fval; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 247:
#line 791 "Optparser.y"
    { yyval.sum.a=yyvsp[-6].fval; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 248:
#line 794 "Optparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=yyvsp[0].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 249:
#line 797 "Optparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 250:
#line 800 "Optparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 251:
#line 803 "Optparser.y"
    { yyval.sum.a=yyvsp[-2].fval; yyval.sum.oprtyp=yyvsp[0].opr.oprtyp; yyval.sum.oprnum=yyvsp[0].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 252:
#line 806 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 253:
#line 809 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 254:
#line 812 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=yyvsp[0].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 255:
#line 815 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 256:
#line 818 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 257:
#line 821 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[0].opr.oprtyp; yyval.sum.oprnum=yyvsp[0].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 258:
#line 824 "Optparser.y"
    { yyval.sum.a=yyvsp[-7].fval; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;  ;}
    break;

  case 259:
#line 827 "Optparser.y"
    { yyval.sum.a=yyvsp[-7].fval; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 260:
#line 830 "Optparser.y"
    { yyval.sum.a=yyvsp[-5].fval; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=yyvsp[-1].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 261:
#line 833 "Optparser.y"
    { yyval.sum.a=yyvsp[-5].fval; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 262:
#line 836 "Optparser.y"
    { yyval.sum.a=yyvsp[-5].fval; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 263:
#line 839 "Optparser.y"
    { yyval.sum.a=yyvsp[-3].fval; yyval.sum.oprtyp=yyvsp[-1].opr.oprtyp; yyval.sum.oprnum=yyvsp[-1].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 264:
#line 842 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 265:
#line 845 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 266:
#line 848 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=yyvsp[-1].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 267:
#line 851 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 268:
#line 854 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 269:
#line 857 "Optparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-1].opr.oprtyp; yyval.sum.oprnum=yyvsp[-1].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;  ;}
    break;

  case 270:
#line 860 "Optparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=0.0; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 271:
#line 863 "Optparser.y"
    { yyval.sum.a=yyvsp[-2].fval; yyval.sum.oprtyp=yyvsp[0].opr.oprtyp; yyval.sum.oprnum=yyvsp[0].opr.oprnum; yyval.sum.p=0.0; yyval.sum.b=1;  
	  yyval.sum.igen=0; ;}
    break;

  case 272:
#line 869 "Optparser.y"
    { yyval.opr.oprtyp=0; yyval.opr.oprnum=yyvsp[-1].ival-1; ;}
    break;

  case 273:
#line 871 "Optparser.y"
    { yyval.opr.oprtyp=1; yyval.opr.oprnum=yyvsp[-1].ival-1; ;}
    break;

  case 274:
#line 873 "Optparser.y"
    { yyval.opr.oprtyp=2; yyval.opr.oprnum=yyvsp[-1].ival-1; ;}
    break;

  case 275:
#line 878 "Optparser.y"
    { yyval.fdef.num=yyvsp[-3].ival-1; yyval.fdef.falldata=yyvsp[0].fall; ;}
    break;

  case 276:
#line 883 "Optparser.y"
    { yyval.gen.a=yyvsp[-5].ival-1 ; yyval.gen.b = yyvsp[-3].ival-1 ; yyval.gen.s = yyvsp[-1].ival ; ;}
    break;

  case 277:
#line 885 "Optparser.y"
    { yyval.gen.a=yyvsp[-3].ival-1 ; yyval.gen.b = yyvsp[-1].ival-1 ; yyval.gen.s = 1 ; ;}
    break;

  case 278:
#line 890 "Optparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival-1); ;}
    break;

  case 279:
#line 892 "Optparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival-1); ;}
    break;

  case 280:
#line 894 "Optparser.y"
    { yyval.ilist.add(yyvsp[0].ival-1); ;}
    break;

  case 281:
#line 896 "Optparser.y"
    { ;}
    break;

  case 282:
#line 901 "Optparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival); ;}
    break;

  case 283:
#line 903 "Optparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival); ;}
    break;

  case 284:
#line 905 "Optparser.y"
    { yyval.ilist.add(yyvsp[0].ival); ;}
    break;

  case 285:
#line 907 "Optparser.y"
    { ;}
    break;

  case 286:
#line 912 "Optparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 287:
#line 914 "Optparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 288:
#line 916 "Optparser.y"
    { yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 289:
#line 918 "Optparser.y"
    { ;}
    break;

  case 290:
#line 923 "Optparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 291:
#line 925 "Optparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 292:
#line 927 "Optparser.y"
    { yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 293:
#line 929 "Optparser.y"
    { ;}
    break;

  case 294:
#line 934 "Optparser.y"
    { yyval.ival = yyvsp[0].ival; ;}
    break;

  case 295:
#line 939 "Optparser.y"
    { yyval.fval = yyvsp[0].ival; ;}
    break;

  case 296:
#line 941 "Optparser.y"
    { yyval.fval = yyvsp[0].fval; ;}
    break;

  case 301:
#line 956 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->setMovingNodesFlag(yyvsp[-1].ival != 0); ;}
    break;

  case 302:
#line 960 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->setDsgnDepRHSFlag(yyvsp[-1].ival != 0); ;}
    break;

  case 303:
#line 964 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData().nodDensFunc    = yyvsp[-1].ival; ;}
    break;

  case 304:
#line 966 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData().nodDensFunc    = yyvsp[-2].ival;
         static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData().nodDensParam1  = yyvsp[-1].fval; ;}
    break;

  case 305:
#line 969 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData().nodDensFunc    = yyvsp[-3].ival;
         static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData().nodDensParam1  = yyvsp[-2].fval;
         static_cast<GeoSource_opt*>(geoSource)->getNodalDensityData().nodDensParam2  = yyvsp[-1].fval; ;}
    break;

  case 306:
#line 976 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjFlag = yyvsp[-1].ival; ;}
    break;

  case 307:
#line 978 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjFlag = yyvsp[-2].ival;
         static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjParam1 = yyvsp[-1].fval; ;}
    break;

  case 308:
#line 981 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjFlag = yyvsp[-3].ival;
         static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjParam1 = yyvsp[-2].fval;
         static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjParam2 = yyvsp[-1].fval; ;}
    break;

  case 309:
#line 985 "Optparser.y"
    { static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjFlag = yyvsp[-4].ival;
         static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjParam1 = yyvsp[-3].fval;
         static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjParam2 = yyvsp[-2].fval;
	 static_cast<GeoSource_opt*>(geoSource)->getDensityProjData().densProjParam3 = yyvsp[-1].fval; ;}
    break;


    }

/* Line 1000 of yacc.c.  */
#line 3667 "Optparser.C"

  yyvsp -= yylen;
  yyssp -= yylen;


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
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  YYSIZE_T yysize = 0;
	  int yytype = YYTRANSLATE (yychar);
	  const char* yyprefix;
	  char *yymsg;
	  int yyx;

	  /* Start YYX at -YYN if negative to avoid negative indexes in
	     YYCHECK.  */
	  int yyxbegin = yyn < 0 ? -yyn : 0;

	  /* Stay within bounds of both yycheck and yytname.  */
	  int yychecklim = YYLAST - yyn;
	  int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
	  int yycount = 0;

	  yyprefix = ", expecting ";
	  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	      {
		yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
		yycount += 1;
		if (yycount == 5)
		  {
		    yysize = 0;
		    break;
		  }
	      }
	  yysize += (sizeof ("syntax error, unexpected ")
		     + yystrlen (yytname[yytype]));
	  yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg != 0)
	    {
	      char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
	      yyp = yystpcpy (yyp, yytname[yytype]);

	      if (yycount < 5)
		{
		  yyprefix = ", expecting ";
		  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
		    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
		      {
			yyp = yystpcpy (yyp, yyprefix);
			yyp = yystpcpy (yyp, yytname[yyx]);
			yyprefix = " or ";
		      }
		}
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    yyerror ("syntax error; also virtual memory exhausted");
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
	     then the rest of the stack, then return failure.  */
	  if (yychar == YYEOF)
	     for (;;)
	       {
		 YYPOPSTACK;
		 if (yyssp == yyss)
		   YYABORT;
		 YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
		 yydestruct (yystos[*yyssp], yyvsp);
	       }
        }
      else
	{
	  YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
	  yydestruct (yytoken, &yylval);
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

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

  yyvsp -= yylen;
  yyssp -= yylen;
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

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


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
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 992 "Optparser.y"


#endif


