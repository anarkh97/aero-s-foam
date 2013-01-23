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
#define yyparse yyrelparse
#define yylex   yyrellex
#define yyerror yyrelerror
#define yylval  yyrellval
#define yychar  yyrelchar
#define yydebug yyreldebug
#define yynerrs yyrelnerrs


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




/* Copy the first part of user declarations.  */
#line 1 "Relparser.y"

#ifdef STRUCTOPT

#include <cstdio>
#include <cstdlib>
#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Relsol.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Structopt_sd.h>
extern int yyrellex(void);
extern void yyrelerror(const char  *);
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
/* Line 191 of yacc.c.  */
#line 428 "Relparser.C"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 440 "Relparser.C"

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
#define YYFINAL  25
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   754

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  167
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  49
/* YYNRULES -- Number of rules. */
#define YYNRULES  264
/* YYNRULES -- Number of states. */
#define YYNSTATES  570

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   408

#define YYTRANSLATE(YYX) 						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     164,   166,   156,   154,   165,   155,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   160,     2,
       2,   163,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   161,     2,   162,   157,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   158,     2,   159,     2,     2,     2,     2,
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
     145,   146,   147,   148,   149,   150,   151,   152,   153
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     6,     8,    11,    13,    15,    17,    19,
      21,    23,    25,    29,    32,    36,    42,    47,    54,    59,
      65,    72,    80,    82,    87,    89,    94,    97,   101,   105,
     110,   114,   118,   120,   122,   125,   128,   130,   134,   136,
     143,   153,   161,   169,   178,   187,   190,   193,   196,   199,
     202,   205,   208,   210,   214,   217,   220,   222,   224,   227,
     229,   231,   233,   235,   237,   239,   241,   243,   245,   247,
     249,   251,   253,   255,   257,   259,   261,   263,   265,   267,
     269,   271,   273,   275,   277,   280,   283,   288,   295,   299,
     305,   309,   312,   320,   326,   335,   342,   350,   357,   364,
     370,   378,   385,   387,   389,   391,   393,   395,   397,   399,
     401,   403,   405,   407,   409,   411,   413,   415,   417,   421,
     423,   425,   427,   429,   431,   433,   435,   437,   439,   441,
     443,   445,   448,   450,   454,   457,   461,   465,   470,   478,
     480,   482,   484,   486,   488,   493,   497,   501,   508,   515,
     522,   529,   536,   543,   547,   551,   556,   561,   565,   570,
     574,   579,   583,   588,   592,   597,   601,   606,   610,   615,
     619,   624,   628,   633,   637,   642,   646,   651,   655,   660,
     665,   670,   674,   687,   704,   706,   708,   710,   712,   714,
     716,   720,   723,   728,   734,   740,   747,   755,   762,   770,
     779,   781,   783,   785,   787,   789,   792,   796,   800,   805,
     807,   810,   813,   817,   819,   822,   825,   833,   841,   847,
     853,   859,   863,   869,   875,   879,   883,   887,   889,   898,
     907,   914,   921,   928,   933,   940,   947,   952,   957,   962,
     965,   970,   975,   980,   987,   996,  1003,  1005,  1008,  1011,
    1014,  1016,  1019,  1022,  1025,  1030,  1036,  1042,  1045,  1050,
    1056,  1062,  1065,  1067,  1069
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const short yyrhs[] =
{
     168,     0,    -1,   169,    45,    -1,   170,    -1,   169,   170,
      -1,   171,    -1,   177,    -1,   178,    -1,   180,    -1,   185,
      -1,   200,    -1,   189,    -1,    26,   150,   172,    -1,   171,
     172,    -1,   214,   173,   150,    -1,   214,   173,    14,   214,
     150,    -1,   214,   173,   176,   150,    -1,   214,   173,   176,
      14,   214,   150,    -1,   214,   173,   209,   150,    -1,   214,
     173,   176,   209,   150,    -1,   214,   173,    14,   214,   209,
     150,    -1,   214,   173,   176,    14,   214,   209,   150,    -1,
      46,    -1,    46,   158,   210,   159,    -1,    71,    -1,    71,
     158,   210,   159,    -1,    52,   214,    -1,   112,   214,   174,
      -1,    34,   214,   175,    -1,    34,   214,   175,   214,    -1,
     142,   214,   175,    -1,    12,   214,   175,    -1,    68,    -1,
      30,    -1,    13,   214,    -1,    77,   214,    -1,   105,    -1,
      67,   214,   175,    -1,    25,    -1,   114,   174,   215,   215,
     214,    81,    -1,   114,   174,   215,   215,   214,    81,   158,
     210,   159,    -1,    33,   214,   214,   215,   158,   212,   159,
      -1,    33,    20,   214,   215,   158,   212,   159,    -1,   113,
     214,   214,   215,    81,   158,   213,   159,    -1,   113,    20,
     214,   215,    81,   158,   213,   159,    -1,    50,   214,    -1,
      48,    49,    -1,    48,   214,    -1,   100,   214,    -1,    95,
     214,    -1,    66,   175,    -1,   136,   214,    -1,    96,    -1,
      82,   214,   175,    -1,    42,   214,    -1,    41,   214,    -1,
      40,    -1,   138,    -1,   106,   214,    -1,   143,    -1,   144,
      -1,   145,    -1,   115,    -1,   116,    -1,   117,    -1,   118,
      -1,   119,    -1,   120,    -1,   121,    -1,   122,    -1,   123,
      -1,   124,    -1,   125,    -1,   126,    -1,   127,    -1,   128,
      -1,   129,    -1,   214,    -1,    37,    -1,    38,    -1,    39,
      -1,    97,    -1,    98,    -1,    99,    -1,    70,   214,    -1,
     139,   215,    -1,     4,   150,   214,   202,    -1,     4,   150,
     214,    10,   215,   202,    -1,   177,   214,   202,    -1,   177,
     214,    10,   215,   202,    -1,     5,   150,   179,    -1,   178,
     179,    -1,   214,     6,   215,   215,   215,   215,   150,    -1,
     214,     6,   215,   215,   150,    -1,   214,     6,   215,   215,
     215,   215,   209,   150,    -1,   214,     6,   215,   215,   209,
     150,    -1,    11,   150,   214,   135,   214,   181,   202,    -1,
     180,   214,   135,   214,   181,   202,    -1,    11,   150,   214,
     135,   183,   202,    -1,   180,   214,   135,   183,   202,    -1,
      11,   150,   214,   135,   214,   184,   202,    -1,   180,   214,
     135,   214,   184,   202,    -1,   214,    -1,    17,    -1,   149,
      -1,    32,    -1,    24,    -1,    27,    -1,    93,    -1,   140,
      -1,   141,    -1,   137,    -1,    37,    -1,    38,    -1,    39,
      -1,    97,    -1,    98,    -1,    99,    -1,    69,   214,   182,
      -1,   214,    -1,    43,    -1,    44,    -1,    57,    -1,    32,
      -1,   137,    -1,    94,    -1,    28,    -1,    29,    -1,   147,
      -1,    18,    -1,   148,    -1,   101,   175,    -1,   146,    -1,
     130,   150,   186,    -1,   185,   186,    -1,   131,   187,   150,
      -1,   132,   188,   150,    -1,   132,   188,   214,   150,    -1,
     133,   214,   214,   215,   215,   214,   150,    -1,   214,    -1,
      65,    -1,    23,    -1,   214,    -1,    19,    -1,   110,   150,
     214,   190,    -1,   189,   214,   190,    -1,     7,   150,   191,
      -1,     3,   150,    78,   150,   192,   197,    -1,     3,   150,
      85,   150,   193,   197,    -1,     3,   150,   107,   150,   194,
     197,    -1,     3,   150,    73,   150,   195,   197,    -1,     3,
     150,   102,   150,   196,   197,    -1,     3,   150,    74,   150,
     195,   197,    -1,     8,   214,   150,    -1,     9,   215,   150,
      -1,   191,     8,   214,   150,    -1,   191,     9,   215,   150,
      -1,    79,   214,   150,    -1,   192,    79,   214,   150,    -1,
      80,   215,   150,    -1,   192,    80,   215,   150,    -1,    87,
     214,   150,    -1,   193,    87,   214,   150,    -1,    86,   215,
     150,    -1,   193,    86,   215,   150,    -1,   109,   214,   150,
      -1,   194,   109,   214,   150,    -1,   108,   215,   150,    -1,
     194,   108,   215,   150,    -1,    75,   214,   150,    -1,   195,
      75,   214,   150,    -1,    76,   215,   150,    -1,   195,    76,
     215,   150,    -1,   104,   214,   150,    -1,   196,   104,   214,
     150,    -1,   103,   215,   150,    -1,   196,   103,   215,   150,
      -1,    58,   198,   150,    -1,   197,    59,   199,   150,    -1,
     197,    60,   215,   150,    -1,   197,    61,   215,   150,    -1,
     197,    63,   150,    -1,   197,    62,   214,   214,   215,   215,
     215,   215,   158,   211,   159,   150,    -1,   197,    62,   214,
     214,   215,   215,   215,   215,   158,   211,   159,   214,   158,
     210,   159,   150,    -1,    15,    -1,    83,    -1,    35,    -1,
      16,    -1,    51,    -1,    21,    -1,   134,   150,   201,    -1,
     200,   201,    -1,   214,   214,   135,   150,    -1,   214,   214,
     214,   135,   150,    -1,   203,   158,   204,   159,   150,    -1,
     203,   158,   150,   204,   159,   150,    -1,   203,   158,   150,
     204,   150,   159,   150,    -1,   203,   150,   158,   204,   159,
     150,    -1,   203,   150,   158,   150,   204,   159,   150,    -1,
     203,   150,   158,   150,   204,   150,   159,   150,    -1,    92,
      -1,    89,    -1,    90,    -1,    91,    -1,    88,    -1,   205,
     206,    -1,   205,   206,   150,    -1,   204,   205,   206,    -1,
     204,   205,   206,   150,    -1,   206,    -1,   206,   150,    -1,
     204,   206,    -1,   204,   206,   150,    -1,   208,    -1,   204,
     208,    -1,   214,   160,    -1,   215,   156,   207,   157,   215,
     154,   215,    -1,   215,   156,   207,   157,   215,   155,   215,
      -1,   215,   156,   207,   157,   215,    -1,   215,   156,   207,
     154,   215,    -1,   215,   156,   207,   155,   215,    -1,   215,
     156,   207,    -1,   207,   157,   215,   154,   215,    -1,   207,
     157,   215,   155,   215,    -1,   207,   157,   215,    -1,   207,
     154,   215,    -1,   207,   155,   215,    -1,   207,    -1,   215,
     156,   207,   157,   215,   154,   215,   209,    -1,   215,   156,
     207,   157,   215,   155,   215,   209,    -1,   215,   156,   207,
     157,   215,   209,    -1,   215,   156,   207,   154,   215,   209,
      -1,   215,   156,   207,   155,   215,   209,    -1,   215,   156,
     207,   209,    -1,   207,   157,   215,   154,   215,   209,    -1,
     207,   157,   215,   155,   215,   209,    -1,   207,   157,   215,
     209,    -1,   207,   154,   215,   209,    -1,   207,   155,   215,
     209,    -1,   207,   209,    -1,    53,   161,   214,   162,    -1,
      55,   161,   214,   162,    -1,    54,   161,   214,   162,    -1,
      31,   161,   214,   162,   163,   202,    -1,    56,   164,   214,
     165,   214,   165,   214,   166,    -1,    56,   164,   214,   165,
     214,   166,    -1,   214,    -1,   150,   214,    -1,   210,   214,
      -1,   210,   150,    -1,   214,    -1,   150,   214,    -1,   211,
     214,    -1,   211,   150,    -1,   214,   175,   215,   215,    -1,
     150,   214,   175,   215,   215,    -1,   212,   214,   175,   215,
     215,    -1,   212,   150,    -1,   214,   174,   215,   215,    -1,
     150,   214,   174,   215,   215,    -1,   213,   214,   174,   215,
     215,    -1,   213,   150,    -1,   151,    -1,   151,    -1,   152,
      -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   108,   108,   116,   117,   121,   122,   123,   124,   125,
     126,   127,   131,   133,   138,   141,   144,   147,   150,   153,
     156,   159,   165,   167,   169,   171,   173,   175,   177,   179,
     181,   183,   185,   187,   189,   191,   193,   195,   197,   199,
     202,   205,   208,   211,   214,   217,   219,   221,   223,   225,
     227,   229,   231,   233,   235,   237,   239,   241,   243,   248,
     249,   250,   251,   252,   253,   254,   255,   256,   257,   258,
     259,   260,   261,   262,   263,   264,   265,   269,   270,   271,
     272,   273,   274,   275,   279,   281,   286,   288,   290,   292,
     297,   299,   304,   307,   310,   313,   319,   322,   325,   328,
     331,   334,   340,   341,   342,   343,   344,   345,   346,   347,
     348,   349,   350,   351,   352,   353,   354,   355,   356,   360,
     361,   362,   363,   364,   365,   366,   367,   368,   372,   373,
     374,   375,   379,   383,   385,   390,   392,   394,   396,   403,
     404,   405,   409,   410,   414,   416,   421,   423,   425,   427,
     429,   431,   433,   438,   441,   444,   447,   453,   456,   458,
     461,   466,   469,   471,   474,   479,   482,   484,   487,   492,
     495,   497,   500,   505,   507,   509,   511,   516,   518,   520,
     522,   524,   526,   532,   541,   542,   546,   547,   548,   549,
     553,   555,   560,   562,   567,   569,   571,   573,   575,   577,
     582,   583,   584,   585,   586,   590,   601,   612,   628,   644,
     656,   668,   685,   701,   709,   720,   725,   728,   731,   734,
     737,   740,   743,   746,   749,   752,   755,   758,   761,   764,
     767,   770,   773,   776,   779,   782,   785,   788,   791,   794,
     800,   802,   804,   809,   814,   816,   821,   823,   825,   827,
     832,   834,   836,   838,   843,   845,   847,   849,   854,   856,
     858,   860,   865,   870,   872
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "FORMALG", "FAILURE", "RNDVAR",
  "RNDTYPE", "MONTEALG", "MONTEINT", "MONTEFLT", "PMAFLAG", "STCRND",
  "ACCELERATION", "AERO", "ANALYSIS", "ANALYTIC", "ADJOINT", "AREA",
  "AOATTACK", "AUTOMATIC", "AVEFLAG", "CENTRAL", "CONSTRAINT", "CONTINUE",
  "CONVECT", "CTRLCOST", "CRITERIA", "CTE", "CTE1", "CTE2", "DAMPING",
  "DEFOPR", "DENSITY", "DISPLEVEL", "DISPLACMENT", "DIRECT", "DSGVAR",
  "DX", "DY", "DZ", "ELCPOW", "ELCRHS", "ELCVOLT", "EMOD1", "EMOD2", "END",
  "ENERGY", "EQUALITY", "ESTAT", "ESTATCRIT", "FAILPROB", "FORWARD",
  "FREQUENCY", "FUNCRIT", "FUNVAR", "FUNOPR", "GENERATE", "GMOD", "GRDTYP",
  "GRDMTH", "GRDRELINC", "GRDABSINC", "GRDFILTER", "GRDSMI", "INEQUALITY",
  "INIT", "INERTIA", "INTFORCE", "KINETIC", "LAYER", "LCASE", "MASS",
  "MIXED", "MMAALG", "GCMALG", "MMAINT", "MMAFLT", "MOMAER", "NLPALG",
  "NLPINT", "NLPFLT", "NODELE", "NODPOS", "NUMERIC", "OBJECTIVE", "OCMALG",
  "OCMFLT", "OCMINT", "OPCOSINUS", "OPKSFUNC", "OPMULTI", "OPSINUS",
  "OPSUM", "PERMTVY", "PHI", "PMACRITERIA", "PULLIN", "RX", "RY", "RZ",
  "RELINDEX", "RGDVELOC", "SALALG", "SALFLT", "SALINT", "SBOOM", "SDTEMP",
  "SLPALG", "SLPFLT", "SLPINT", "SOLVER", "STCVAR", "STRESS", "STRLEVEL",
  "STRVOL", "S1L", "S1M", "S1U", "S2L", "S2M", "S2U", "SXL", "SXM", "SXU",
  "SYL", "SYM", "SYU", "SZL", "SZM", "SZU", "STCANA", "STCDYN", "STCDAM",
  "STCTRA", "STCELEVAR", "STCVARTYP", "TEMP", "THICK_rel", "THPOWER",
  "TIME", "THERMCOND", "TRADIATE", "VELOCITY", "VML", "VMM", "VMU",
  "VOLTAGE", "XMACH", "YAWANGLE", "YMODUS", "NewLine", "IntConstant",
  "DblConstant", "String", "'+'", "'-'", "'*'", "'^'", "'{'", "'}'", "':'",
  "'['", "']'", "'='", "'('", "','", "')'", "$accept", "FinalizedData",
  "All", "Component", "CriteriaSet", "Criteria", "CritTyp", "StrTyp",
  "DispTyp", "TimeCase", "FailureCriteriaSet", "RndVariableSet",
  "RndVariable", "StcVariableSet", "AccTyp", "CompTyp", "FldAttrTyp",
  "ElecTyp", "StcAnalysisSet", "StcAnalysisData", "StcDynTyp",
  "StcDampTyp", "SolverSet", "SolverTyp", "MonteData", "NlpData",
  "OcmData", "SlpData", "MmaData", "SalData", "GradData", "GradTyp",
  "GradMethod", "StcEleVarSet", "StcEleVarInfo", "Function", "FuncTyp",
  "StatementSet", "StatementNumber", "Statement", "OprTyp", "FunctionDef",
  "Generate", "IntegerList", "DirectIntegerList", "DispLevelList",
  "StrLevelList", "Integer", "Float", 0
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
     405,   406,   407,   408,    43,    45,    42,    94,   123,   125,
      58,    91,    93,    61,    40,    44,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,   167,   168,   169,   169,   170,   170,   170,   170,   170,
     170,   170,   171,   171,   172,   172,   172,   172,   172,   172,
     172,   172,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   173,
     173,   173,   173,   173,   173,   173,   173,   173,   173,   174,
     174,   174,   174,   174,   174,   174,   174,   174,   174,   174,
     174,   174,   174,   174,   174,   174,   174,   175,   175,   175,
     175,   175,   175,   175,   176,   176,   177,   177,   177,   177,
     178,   178,   179,   179,   179,   179,   180,   180,   180,   180,
     180,   180,   181,   181,   181,   181,   181,   181,   181,   181,
     181,   181,   181,   181,   181,   181,   181,   181,   181,   182,
     182,   182,   182,   182,   182,   182,   182,   182,   183,   183,
     183,   183,   184,   185,   185,   186,   186,   186,   186,   187,
     187,   187,   188,   188,   189,   189,   190,   190,   190,   190,
     190,   190,   190,   191,   191,   191,   191,   192,   192,   192,
     192,   193,   193,   193,   193,   194,   194,   194,   194,   195,
     195,   195,   195,   196,   196,   196,   196,   197,   197,   197,
     197,   197,   197,   197,   198,   198,   199,   199,   199,   199,
     200,   200,   201,   201,   202,   202,   202,   202,   202,   202,
     203,   203,   203,   203,   203,   204,   204,   204,   204,   204,
     204,   204,   204,   204,   204,   205,   206,   206,   206,   206,
     206,   206,   206,   206,   206,   206,   206,   206,   206,   206,
     206,   206,   206,   206,   206,   206,   206,   206,   206,   206,
     207,   207,   207,   208,   209,   209,   210,   210,   210,   210,
     211,   211,   211,   211,   212,   212,   212,   212,   213,   213,
     213,   213,   214,   215,   215
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     2,     1,     2,     1,     1,     1,     1,     1,
       1,     1,     3,     2,     3,     5,     4,     6,     4,     5,
       6,     7,     1,     4,     1,     4,     2,     3,     3,     4,
       3,     3,     1,     1,     2,     2,     1,     3,     1,     6,
       9,     7,     7,     8,     8,     2,     2,     2,     2,     2,
       2,     2,     1,     3,     2,     2,     1,     1,     2,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     2,     2,     4,     6,     3,     5,
       3,     2,     7,     5,     8,     6,     7,     6,     6,     5,
       7,     6,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     1,     3,     2,     3,     3,     4,     7,     1,
       1,     1,     1,     1,     4,     3,     3,     6,     6,     6,
       6,     6,     6,     3,     3,     4,     4,     3,     4,     3,
       4,     3,     4,     3,     4,     3,     4,     3,     4,     3,
       4,     3,     4,     3,     4,     3,     4,     3,     4,     4,
       4,     3,    12,    16,     1,     1,     1,     1,     1,     1,
       3,     2,     4,     5,     5,     6,     7,     6,     7,     8,
       1,     1,     1,     1,     1,     2,     3,     3,     4,     1,
       2,     2,     3,     1,     2,     2,     7,     7,     5,     5,
       5,     3,     5,     5,     3,     3,     3,     1,     8,     8,
       6,     6,     6,     4,     6,     6,     4,     4,     4,     2,
       4,     4,     4,     6,     8,     6,     1,     2,     2,     2,
       1,     2,     2,     2,     4,     5,     5,     2,     4,     5,
       5,     2,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned short yydefact[] =
{
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       3,     5,     6,     7,     8,     9,    11,    10,     0,     0,
       0,     0,     0,     0,     0,     1,     2,     4,   262,    13,
       0,     0,    91,     0,     0,     0,     0,     0,   134,     0,
     191,     0,     0,    90,     0,    12,     0,   133,   190,     0,
       0,    38,    33,     0,     0,    56,     0,     0,    22,     0,
       0,     0,     0,     0,    32,    24,     0,     0,     0,    52,
       0,    36,     0,     0,     0,     0,     0,    57,     0,     0,
       0,   204,   201,   202,   203,   200,    88,     0,     0,     0,
     141,   140,     0,   139,   143,     0,   142,     0,     0,     0,
     145,     0,     0,    86,     0,   144,     0,    34,     0,     0,
       0,    55,    54,     0,    46,    47,    45,    26,    78,    79,
      80,    81,    82,    83,    50,    77,     0,     0,    35,     0,
      49,    48,    58,     0,     0,     0,    62,    63,    64,    65,
      66,    67,    68,    69,    70,    71,    72,    73,    74,    75,
      76,    59,    60,    61,     0,    51,     0,     0,     0,     0,
       0,    14,     0,     0,   263,   264,     0,     0,     0,     0,
     129,     0,   128,   130,     0,     0,   135,   136,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    31,     0,     0,
      28,     0,     0,   246,    37,     0,    53,    27,     0,     0,
       0,    30,     0,     0,    84,    85,     0,    16,     0,    18,
      89,     0,     0,     0,     0,     0,     0,   262,     0,     0,
     209,   227,   213,     0,     0,     0,   131,    99,   103,   106,
     107,   105,   112,   113,   114,     0,   108,   115,   116,   117,
     111,   109,   110,   132,   104,     0,     0,   102,   137,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   146,   192,
       0,    87,    98,     0,     0,     0,     0,    29,   247,   249,
      23,   248,    25,     0,     0,     0,    15,     0,     0,     0,
      19,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     211,   214,   205,   210,     0,     0,     0,   239,   215,     0,
      93,     0,     0,     0,    97,   101,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   193,    96,   100,
       0,     0,     0,     0,     0,    20,     0,    17,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   194,   207,   212,
     206,   225,   226,   224,   221,    95,     0,   126,   127,   123,
     120,   121,   122,   125,   124,   118,   119,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   153,   154,     0,     0,     0,     0,
       0,     0,     0,     0,    39,     0,    21,     0,     0,   197,
       0,   240,   242,   241,     0,   195,   208,   237,   238,     0,
       0,   236,     0,     0,     0,   233,    92,     0,   138,     0,
       0,     0,     0,     0,   150,   152,     0,     0,     0,     0,
     147,     0,     0,     0,     0,   148,     0,     0,     0,     0,
     151,     0,     0,     0,     0,   149,   155,   156,     0,   257,
      42,     0,     0,    41,     0,     0,     0,     0,     0,     0,
     245,     0,   198,     0,   196,   222,   223,   219,   220,   218,
      94,   169,   171,   184,   185,     0,     0,     0,     0,     0,
       0,     0,     0,   157,   159,     0,     0,   163,   161,     0,
       0,   175,   173,     0,     0,   167,   165,     0,     0,     0,
       0,     0,     0,   261,    44,     0,     0,    43,     0,     0,
     199,   243,   234,   235,   231,   232,     0,     0,   230,   177,
     170,   172,   187,   189,   186,   188,     0,     0,     0,     0,
     181,   158,   160,   164,   162,   176,   174,   168,   166,     0,
       0,   254,     0,     0,     0,    40,   244,   216,   217,   178,
     179,   180,     0,   255,   256,     0,     0,   258,   228,   229,
       0,   259,   260,     0,     0,     0,     0,     0,     0,   250,
     251,   253,     0,   252,   182,     0,     0,     0,     0,   183
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,     8,     9,    10,    11,    29,    79,   154,   124,   162,
      12,    13,    32,    14,   245,   355,   174,   246,    15,    38,
      92,    95,    16,   100,   258,   364,   367,   373,   360,   370,
     414,   465,   516,    17,    40,    86,    87,   218,   219,   220,
     221,   222,   163,   192,   558,   379,   445,   125,   224
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -181
static const short yypact[] =
{
     100,  -118,  -106,  -103,   -99,   -84,   -71,   -61,    94,   116,
    -181,   -64,   -64,   -64,   -64,   313,   -64,   -64,   -64,   -64,
     -64,   -64,   -64,   313,   -64,  -181,  -181,  -181,  -181,  -181,
     605,   249,  -181,    97,   -29,     5,    -7,   -64,  -181,    34,
    -181,   -64,   363,  -181,   -26,  -181,    34,  -181,  -181,   -64,
     -64,  -181,  -181,    -6,   -64,  -181,   -64,   -64,   -36,   -34,
     -64,   -64,   215,   -64,  -181,    -3,   -64,   -64,   -64,  -181,
     -64,  -181,   -64,   -64,     7,   405,   -64,  -181,   -64,    15,
     -76,  -181,  -181,  -181,  -181,  -181,  -181,   -81,   -76,    50,
    -181,  -181,   -17,  -181,  -181,    52,  -181,   -64,   -13,     3,
    -181,   -63,   -76,  -181,    50,  -181,   215,  -181,   -64,   -64,
     215,  -181,  -181,    72,  -181,  -181,  -181,  -181,  -181,  -181,
    -181,  -181,  -181,  -181,  -181,  -181,   215,    72,  -181,   215,
    -181,  -181,  -181,   405,   -64,   -64,  -181,  -181,  -181,  -181,
    -181,  -181,  -181,  -181,  -181,  -181,  -181,  -181,  -181,  -181,
    -181,  -181,  -181,  -181,   -76,  -181,   215,   -64,    -1,   -64,
     -76,  -181,     2,    16,  -181,  -181,   401,    12,   269,   -76,
    -181,   215,  -181,  -181,   401,   150,  -181,  -181,    36,   -76,
     357,   257,    56,    73,   401,   401,   150,  -181,   -76,   -76,
     -64,   -64,    40,  -181,  -181,   231,  -181,  -181,   -76,   -76,
     -76,  -181,   -37,   -64,  -181,  -181,   -64,  -181,    66,  -181,
    -181,   289,    74,    76,    99,   108,   331,   125,   224,   -33,
     136,    39,  -181,   129,   141,   106,  -181,  -181,  -181,  -181,
    -181,  -181,  -181,  -181,  -181,   -64,  -181,  -181,  -181,  -181,
    -181,  -181,  -181,  -181,  -181,   401,   401,  -181,  -181,   -76,
     158,   161,   169,   171,   198,   202,   -64,   -76,   285,  -181,
     213,  -181,  -181,   401,   401,   145,   201,  -181,  -181,  -181,
    -181,  -181,  -181,   293,   306,   -64,  -181,   256,   226,   -22,
    -181,   331,   252,   -64,   -64,   -64,   -64,    77,   260,   -33,
     262,  -181,   264,  -181,   -76,   -76,   -76,  -181,  -181,   450,
    -181,   272,   -76,   321,  -181,  -181,   -64,   241,   241,   253,
     259,   232,   248,   274,   279,   -64,   -76,  -181,  -181,  -181,
     266,   266,   292,   303,   351,  -181,   -64,  -181,   323,    81,
     334,   301,   324,   352,   353,   322,   367,  -181,   368,  -181,
    -181,   457,   457,   -16,    87,  -181,   -21,  -181,  -181,  -181,
    -181,  -181,  -181,  -181,  -181,  -181,  -181,   369,   -64,   -76,
     209,   209,   -64,   -76,    89,   -76,   -64,    -5,   -76,   -64,
      -4,   -76,   -64,   153,  -181,  -181,   385,   386,   -64,   243,
     215,   246,   276,   276,   379,   341,  -181,   380,   388,  -181,
     377,  -181,  -181,  -181,   391,  -181,  -181,  -181,  -181,   -76,
     -76,  -181,   -76,   -76,   -76,  -181,  -181,   392,  -181,   393,
     394,     8,   -64,   -76,   438,   438,   402,   403,   -64,   -76,
     438,   404,   407,   -76,   -64,   438,   408,   410,   -76,   -64,
     438,   411,   412,   -76,   -64,   438,  -181,  -181,   215,  -181,
    -181,   215,   -76,  -181,   -64,   297,   405,   316,    72,   -64,
    -181,   414,  -181,   401,  -181,   457,   457,   457,   457,    18,
    -181,  -181,  -181,  -181,  -181,   416,   417,   419,   247,   -76,
     -76,   -64,   420,  -181,  -181,   422,   423,  -181,  -181,   424,
     425,  -181,  -181,   426,   427,  -181,  -181,   428,   429,   -76,
     -76,   -76,   405,  -181,  -181,   405,   -76,  -181,   319,   415,
    -181,  -181,  -181,  -181,  -181,  -181,   -76,   -76,  -181,  -181,
    -181,  -181,  -181,  -181,  -181,  -181,   430,   432,   433,   -64,
    -181,  -181,  -181,  -181,  -181,  -181,  -181,  -181,  -181,   -76,
     -76,  -181,   -76,   -76,   -76,  -181,  -181,   457,   457,  -181,
    -181,  -181,   -76,  -181,  -181,   -76,   -76,  -181,  -181,  -181,
     -76,  -181,  -181,   -76,   -76,   431,   359,   -64,   326,  -181,
    -181,  -181,   361,  -181,  -181,   434,    72,   329,   435,  -181
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
    -181,  -181,  -181,   575,  -181,   565,  -181,  -115,   -46,  -181,
    -181,  -181,   568,  -181,   409,  -181,   484,   421,  -181,   567,
    -181,  -181,  -181,   545,  -181,  -181,  -181,  -181,   286,  -181,
     101,  -181,  -181,  -181,   569,   -25,  -181,  -180,  -111,  -104,
     298,   -78,  -129,  -123,  -181,   275,   218,   -11,   200
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -264
static const short yytable[] =
{
      30,    31,    33,    34,   195,    39,    41,    42,    33,    44,
      30,    46,    94,    41,   108,   114,   206,   103,   197,   158,
     213,   214,   215,   463,    93,    96,    97,   134,    90,   157,
     101,   282,    18,   208,   158,   158,   287,    98,   106,   107,
     158,    99,   109,   110,    19,   111,   112,    20,   115,   116,
     117,    21,   126,   411,   411,   128,   129,   130,   158,   131,
     187,   132,   133,   135,   190,   155,    22,   156,   170,   167,
      91,   158,   182,   277,   158,   164,   165,   168,   175,    23,
     194,   423,   424,   196,   178,   159,   179,    28,    28,    24,
     183,   464,   297,   186,    25,   158,   301,   188,   189,   428,
     429,   329,   193,    88,     1,     2,    89,   289,   212,   104,
     201,     3,   212,   276,   290,   292,   193,    28,   164,   165,
       1,     2,   113,   198,   199,   226,     4,     3,   327,   406,
     213,   214,   215,   176,   213,   214,   215,   180,   399,   400,
     291,   210,     4,   158,    28,    28,   202,   411,   204,   227,
     328,   171,   207,   181,   160,   127,    28,   223,    28,   261,
     262,    26,   158,   203,   247,   161,   209,   228,   418,   419,
     211,   289,   506,   507,   229,   247,   289,   230,   290,   267,
     268,   271,   231,   290,   271,   338,   248,   232,   233,   234,
     269,    28,   278,   294,   295,   279,   296,   172,   173,   270,
     223,    28,   177,    28,   291,   223,   259,   223,   260,   291,
       5,   411,   397,   398,   401,   405,   280,   407,   289,   235,
     304,   305,   191,    28,   303,   290,     5,   335,   217,   165,
       6,   387,   217,   165,     7,   283,   336,   284,   318,   319,
     388,   402,   403,   236,   404,   313,     6,   237,   238,   239,
       7,   291,   118,   119,   120,   212,   300,   164,   165,    80,
     285,   433,   434,   512,   324,   256,   257,   411,   513,   286,
     223,   223,   331,   332,   333,   334,   223,   213,   214,   215,
     166,  -263,   514,   212,   412,   413,   293,   240,   169,   298,
     241,   242,   356,   315,   316,   357,   243,   299,   515,   244,
     212,    28,   184,   320,   376,   213,   214,   215,   307,   380,
     380,   308,   121,   122,   123,   385,   358,   359,   223,   309,
     212,   310,   213,   214,   215,   498,   502,   503,   504,   505,
     508,   496,   362,   363,   442,   368,   369,    81,    82,    83,
      84,    85,   213,   214,   215,   365,   366,   409,   311,   347,
     348,   416,   312,   349,   200,   422,   371,   372,   427,   321,
     205,   432,   212,   317,   350,   351,    28,   438,   441,   225,
     441,   446,   446,   102,   322,   217,   165,   532,   352,   249,
     533,   269,    28,   288,   213,   214,   215,   323,   265,   266,
     272,   326,   489,   439,    28,   490,   439,    28,   273,   274,
     275,   466,   440,   217,   165,   443,   325,   475,   548,   549,
     337,   330,   339,   480,   340,   353,   378,    28,   484,   216,
     217,   165,   345,   488,   374,   302,   444,    28,   501,   375,
     250,   251,   384,   492,   495,   252,   495,   193,   499,   281,
     217,   165,   253,   567,    35,    36,    37,   493,    28,   306,
     382,    81,    82,    83,    84,    85,   494,   314,   354,   254,
     519,   383,   415,   390,   255,   420,   493,    28,   425,   269,
      28,   430,    28,   386,   435,   497,   561,    28,   535,   269,
      28,   394,   217,   165,   389,   562,   391,   271,   568,    81,
      82,    83,    84,    85,   341,   342,   343,   468,   469,   470,
     471,   472,   346,   213,   214,   215,   449,   450,   542,   557,
      28,   564,    28,   158,   392,   393,   377,   395,   396,   408,
     136,   137,   138,   139,   140,   141,   142,   143,   144,   145,
     146,   147,   148,   149,   150,   436,   437,   448,   452,   451,
     453,   454,   460,   461,   462,   559,   560,   563,   151,   152,
     153,   565,   473,   474,   477,   193,   271,   478,   481,   410,
     482,   485,   486,   417,   500,   421,   509,   510,   426,   511,
     520,   431,   521,   522,   523,   524,   525,   526,   527,   528,
     539,   536,   540,   541,    27,   569,    45,    43,   185,   556,
      47,   105,   566,    48,   361,   263,   381,   344,     0,   455,
     456,   447,   457,   458,   459,     0,     0,   264,     0,     0,
       0,     0,     0,   467,     0,     0,     0,    49,    50,   476,
       0,     0,     0,   479,     0,     0,     0,     0,   483,     0,
      51,     0,     0,   487,     0,    52,     0,     0,    53,    54,
       0,     0,   491,     0,     0,    55,    56,    57,     0,     0,
       0,    58,     0,    59,     0,    60,     0,    61,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   517,
     518,    62,    63,    64,     0,     0,    65,     0,     0,     0,
       0,     0,    66,     0,     0,     0,     0,    67,     0,   529,
     530,   531,     0,     0,     0,     0,   534,     0,     0,     0,
      68,    69,     0,     0,     0,    70,   537,   538,     0,     0,
      71,    72,     0,     0,     0,     0,     0,    73,    74,    75,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   543,
     544,     0,   545,   546,   547,     0,     0,     0,     0,     0,
       0,    76,   550,    77,     0,   551,   552,    78,     0,     0,
     553,     0,     0,   554,   555
};

static const short yycheck[] =
{
      11,    12,    13,    14,   127,    16,    17,    18,    19,    20,
      21,    22,    19,    24,    20,    49,    14,    42,   133,    56,
      53,    54,    55,    15,    35,    36,    37,    20,    23,    14,
      41,   211,   150,   162,    56,    56,   216,     3,    49,    50,
      56,     7,    53,    54,   150,    56,    57,   150,    59,    60,
      61,   150,    63,    58,    58,    66,    67,    68,    56,    70,
     106,    72,    73,    74,   110,    76,   150,    78,    18,   150,
      65,    56,   135,   202,    56,   151,   152,   158,    89,   150,
     126,    86,    87,   129,    95,    70,    97,   151,   151,   150,
     101,    83,   221,   104,     0,    56,   225,   108,   109,   103,
     104,   281,   113,     6,     4,     5,   135,   218,    31,   135,
     156,    11,    31,   150,   218,   219,   127,   151,   151,   152,
       4,     5,   158,   134,   135,   171,    26,    11,   150,   150,
      53,    54,    55,   150,    53,    54,    55,   150,   154,   155,
     218,   166,    26,    56,   151,   151,   157,    58,   159,   174,
     279,   101,   150,   150,   139,   158,   151,   168,   151,   184,
     185,    45,    56,   164,   175,   150,   150,    17,    79,    80,
     158,   282,   154,   155,    24,   186,   287,    27,   282,   190,
     191,   192,    32,   287,   195,   289,   150,    37,    38,    39,
     150,   151,   203,   154,   155,   206,   157,   147,   148,   159,
     211,   151,   150,   151,   282,   216,   150,   218,   135,   287,
     110,    58,   341,   342,   343,   344,   150,   346,   329,    69,
     245,   246,   150,   151,   235,   329,   110,   150,   151,   152,
     130,   150,   151,   152,   134,   161,   159,   161,   263,   264,
     159,   154,   155,    93,   157,   256,   130,    97,    98,    99,
     134,   329,    37,    38,    39,    31,   150,   151,   152,    10,
     161,   108,   109,    16,   275,     8,     9,    58,    21,   161,
     281,   282,   283,   284,   285,   286,   287,    53,    54,    55,
      80,   156,    35,    31,    75,    76,   150,   137,    88,   160,
     140,   141,   303,     8,     9,   306,   146,   156,    51,   149,
      31,   151,   102,   158,   315,    53,    54,    55,   150,   320,
     321,   150,    97,    98,    99,   326,    75,    76,   329,   150,
      31,   150,    53,    54,    55,   448,   455,   456,   457,   458,
     459,   446,    79,    80,   380,   103,   104,    88,    89,    90,
      91,    92,    53,    54,    55,    86,    87,   358,   150,    28,
      29,   362,   150,    32,   154,   366,   108,   109,   369,   158,
     160,   372,    31,   150,    43,    44,   151,   378,   379,   169,
     381,   382,   383,    10,    81,   151,   152,   492,    57,   179,
     495,   150,   151,   159,    53,    54,    55,    81,   188,   189,
     159,   165,   438,   150,   151,   441,   150,   151,   198,   199,
     200,   412,   159,   151,   152,   159,   150,   418,   537,   538,
     150,   159,   150,   424,   150,    94,   150,   151,   429,   150,
     151,   152,   150,   434,   150,   225,   150,   151,   453,   150,
      73,    74,    81,   444,   445,    78,   447,   448,   449,   150,
     151,   152,    85,   566,   131,   132,   133,   150,   151,   249,
     158,    88,    89,    90,    91,    92,   159,   257,   137,   102,
     471,   158,   361,   162,   107,   364,   150,   151,   367,   150,
     151,   370,   151,   150,   373,   159,   150,   151,   159,   150,
     151,   159,   151,   152,   150,   159,   162,   498,   159,    88,
      89,    90,    91,    92,   294,   295,   296,    59,    60,    61,
      62,    63,   302,    53,    54,    55,   165,   166,   519,   150,
     151,   150,   151,    56,   162,   162,   316,   150,   150,   150,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   150,   150,   158,   150,   159,
     163,   150,   150,   150,   150,   556,   557,   558,   143,   144,
     145,   562,   150,   150,   150,   566,   567,   150,   150,   359,
     150,   150,   150,   363,   150,   365,   150,   150,   368,   150,
     150,   371,   150,   150,   150,   150,   150,   150,   150,   150,
     150,   166,   150,   150,     9,   150,    21,    19,   104,   158,
      23,    46,   158,    24,   308,   186,   321,   299,    -1,   399,
     400,   383,   402,   403,   404,    -1,    -1,   186,    -1,    -1,
      -1,    -1,    -1,   413,    -1,    -1,    -1,    12,    13,   419,
      -1,    -1,    -1,   423,    -1,    -1,    -1,    -1,   428,    -1,
      25,    -1,    -1,   433,    -1,    30,    -1,    -1,    33,    34,
      -1,    -1,   442,    -1,    -1,    40,    41,    42,    -1,    -1,
      -1,    46,    -1,    48,    -1,    50,    -1,    52,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   469,
     470,    66,    67,    68,    -1,    -1,    71,    -1,    -1,    -1,
      -1,    -1,    77,    -1,    -1,    -1,    -1,    82,    -1,   489,
     490,   491,    -1,    -1,    -1,    -1,   496,    -1,    -1,    -1,
      95,    96,    -1,    -1,    -1,   100,   506,   507,    -1,    -1,
     105,   106,    -1,    -1,    -1,    -1,    -1,   112,   113,   114,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   529,
     530,    -1,   532,   533,   534,    -1,    -1,    -1,    -1,    -1,
      -1,   136,   542,   138,    -1,   545,   546,   142,    -1,    -1,
     550,    -1,    -1,   553,   554
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     4,     5,    11,    26,   110,   130,   134,   168,   169,
     170,   171,   177,   178,   180,   185,   189,   200,   150,   150,
     150,   150,   150,   150,   150,     0,    45,   170,   151,   172,
     214,   214,   179,   214,   214,   131,   132,   133,   186,   214,
     201,   214,   214,   179,   214,   172,   214,   186,   201,    12,
      13,    25,    30,    33,    34,    40,    41,    42,    46,    48,
      50,    52,    66,    67,    68,    71,    77,    82,    95,    96,
     100,   105,   106,   112,   113,   114,   136,   138,   142,   173,
      10,    88,    89,    90,    91,    92,   202,   203,     6,   135,
      23,    65,   187,   214,    19,   188,   214,   214,     3,     7,
     190,   214,    10,   202,   135,   190,   214,   214,    20,   214,
     214,   214,   214,   158,    49,   214,   214,   214,    37,    38,
      39,    97,    98,    99,   175,   214,   214,   158,   214,   214,
     214,   214,   214,   214,    20,   214,   115,   116,   117,   118,
     119,   120,   121,   122,   123,   124,   125,   126,   127,   128,
     129,   143,   144,   145,   174,   214,   214,    14,    56,    70,
     139,   150,   176,   209,   151,   152,   215,   150,   158,   215,
      18,   101,   147,   148,   183,   214,   150,   150,   214,   214,
     150,   150,   135,   214,   215,   183,   214,   175,   214,   214,
     175,   150,   210,   214,   175,   210,   175,   174,   214,   214,
     215,   175,   214,   164,   214,   215,    14,   150,   209,   150,
     202,   158,    31,    53,    54,    55,   150,   151,   204,   205,
     206,   207,   208,   214,   215,   215,   175,   202,    17,    24,
      27,    32,    37,    38,    39,    69,    93,    97,    98,    99,
     137,   140,   141,   146,   149,   181,   184,   214,   150,   215,
      73,    74,    78,    85,   102,   107,     8,     9,   191,   150,
     135,   202,   202,   181,   184,   215,   215,   214,   214,   150,
     159,   214,   159,   215,   215,   215,   150,   209,   214,   214,
     150,   150,   204,   161,   161,   161,   161,   204,   159,   205,
     206,   208,   206,   150,   154,   155,   157,   209,   160,   156,
     150,   209,   215,   214,   202,   202,   215,   150,   150,   150,
     150,   150,   150,   214,   215,     8,     9,   150,   202,   202,
     158,   158,    81,    81,   214,   150,   165,   150,   209,   204,
     159,   214,   214,   214,   214,   150,   159,   150,   206,   150,
     150,   215,   215,   215,   207,   150,   215,    28,    29,    32,
      43,    44,    57,    94,   137,   182,   214,   214,    75,    76,
     195,   195,    79,    80,   192,    86,    87,   193,   103,   104,
     196,   108,   109,   194,   150,   150,   214,   215,   150,   212,
     214,   212,   158,   158,    81,   214,   150,   150,   159,   150,
     162,   162,   162,   162,   159,   150,   150,   209,   209,   154,
     155,   209,   154,   155,   157,   209,   150,   209,   150,   214,
     215,    58,    75,    76,   197,   197,   214,   215,    79,    80,
     197,   215,   214,    86,    87,   197,   215,   214,   103,   104,
     197,   215,   214,   108,   109,   197,   150,   150,   214,   150,
     159,   214,   175,   159,   150,   213,   214,   213,   158,   165,
     166,   159,   150,   163,   150,   215,   215,   215,   215,   215,
     150,   150,   150,    15,    83,   198,   214,   215,    59,    60,
      61,    62,    63,   150,   150,   214,   215,   150,   150,   215,
     214,   150,   150,   215,   214,   150,   150,   215,   214,   175,
     175,   215,   214,   150,   159,   214,   174,   159,   210,   214,
     150,   202,   209,   209,   209,   209,   154,   155,   209,   150,
     150,   150,    16,    21,    35,    51,   199,   215,   215,   214,
     150,   150,   150,   150,   150,   150,   150,   150,   150,   215,
     215,   215,   174,   174,   215,   159,   166,   215,   215,   150,
     150,   150,   214,   215,   215,   215,   215,   215,   209,   209,
     215,   215,   215,   215,   215,   215,   158,   150,   211,   214,
     214,   150,   159,   214,   150,   214,   158,   210,   159,   150
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
#line 109 "Relparser.y"
    { 
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->inputcheck();
	  return 0; 
        ;}
    break;

  case 12:
#line 132 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addCriteria(yyvsp[0].crit); ;}
    break;

  case 13:
#line 134 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addCriteria(yyvsp[0].crit); ;}
    break;

  case 14:
#line 139 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-2].ival; yyval.crit.data=yyvsp[-1].cdata; yyval.crit.t=-1; yyval.crit.anaId=0;
          yyval.crit.igen=0; ;}
    break;

  case 15:
#line 142 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-4].ival; yyval.crit.data=yyvsp[-3].cdata; yyval.crit.t=-1; yyval.crit.anaId=yyvsp[-1].ival-1;
          yyval.crit.igen=0; ;}
    break;

  case 16:
#line 145 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-3].ival; yyval.crit.data=yyvsp[-2].cdata; yyval.crit.t=yyvsp[-1].fval; yyval.crit.anaId=0;
          yyval.crit.igen=0; ;}
    break;

  case 17:
#line 148 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-5].ival; yyval.crit.data=yyvsp[-4].cdata; yyval.crit.t=yyvsp[-3].fval; yyval.crit.anaId=yyvsp[-1].ival-1;
          yyval.crit.igen=0; ;}
    break;

  case 18:
#line 151 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-3].ival; yyval.crit.data=yyvsp[-2].cdata; yyval.crit.t=-1; yyval.crit.anaId=0; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 19:
#line 154 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-4].ival; yyval.crit.data=yyvsp[-3].cdata;   yyval.crit.t=yyvsp[-2].fval; yyval.crit.anaId=0; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 20:
#line 157 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-5].ival; yyval.crit.data=yyvsp[-4].cdata;   yyval.crit.t=-1; yyval.crit.anaId=yyvsp[-2].ival-1; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 21:
#line 160 "Relparser.y"
    { yyval.crit.initialize(); yyval.crit.num=yyvsp[-6].ival; yyval.crit.data=yyvsp[-5].cdata;   yyval.crit.t=yyvsp[-4].fval; yyval.crit.anaId=yyvsp[-2].ival-1; 
	  yyval.crit.igen=1; yyval.crit.gena=yyvsp[-1].gen.a; yyval.crit.genb=yyvsp[-1].gen.b; yyval.crit.gens=yyvsp[-1].gen.s; ;}
    break;

  case 22:
#line 166 "Relparser.y"
    { yyval.cdata.typ =  0; ;}
    break;

  case 23:
#line 168 "Relparser.y"
    { yyval.cdata.typ =  0; yyval.cdata.ip = yyvsp[-1].ilist.ip; yyval.cdata.ipsize = yyvsp[-1].ilist.ipsize; ;}
    break;

  case 24:
#line 170 "Relparser.y"
    { yyval.cdata.typ =  1; ;}
    break;

  case 25:
#line 172 "Relparser.y"
    { yyval.cdata.typ =  1; yyval.cdata.ip = yyvsp[-1].ilist.ip; yyval.cdata.ipsize = yyvsp[-1].ilist.ipsize; ;}
    break;

  case 26:
#line 174 "Relparser.y"
    { yyval.cdata.typ =  2; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 27:
#line 176 "Relparser.y"
    { yyval.cdata.typ =  3; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 28:
#line 178 "Relparser.y"
    { yyval.cdata.typ =  4; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; yyval.cdata.i[2]=0; ;}
    break;

  case 29:
#line 180 "Relparser.y"
    { yyval.cdata.typ =  4; yyval.cdata.i[0]=yyvsp[-2].ival-1; yyval.cdata.i[1]=yyvsp[-1].ival; yyval.cdata.i[2]=yyvsp[0].ival-1;;}
    break;

  case 30:
#line 182 "Relparser.y"
    { yyval.cdata.typ =  5; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 31:
#line 184 "Relparser.y"
    { yyval.cdata.typ =  6; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 32:
#line 186 "Relparser.y"
    { yyval.cdata.typ =  7; ;}
    break;

  case 33:
#line 188 "Relparser.y"
    { yyval.cdata.typ =  8; ;}
    break;

  case 34:
#line 190 "Relparser.y"
    { yyval.cdata.typ =  9; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 35:
#line 192 "Relparser.y"
    { yyval.cdata.typ =  9; yyval.cdata.i[0]=yyvsp[0].ival+2; ;}
    break;

  case 36:
#line 194 "Relparser.y"
    { yyval.cdata.typ =  9; yyval.cdata.i[0]=6; ;}
    break;

  case 37:
#line 196 "Relparser.y"
    { yyval.cdata.typ = 10; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 38:
#line 198 "Relparser.y"
    { yyval.cdata.typ = 11; ;}
    break;

  case 39:
#line 200 "Relparser.y"
    { yyval.cdata.typ = 12; yyval.cdata.i[0]=yyvsp[-4].ival; yyval.cdata.d[0]=yyvsp[-3].fval; yyval.cdata.d[1]=yyvsp[-2].fval; yyval.cdata.i[1]=yyvsp[-1].ival; yyval.cdata.i[2]=yyvsp[0].ival; 
          yyval.cdata.ipsize=0; ;}
    break;

  case 40:
#line 203 "Relparser.y"
    { yyval.cdata.typ = 12; yyval.cdata.i[0]=yyvsp[-7].ival; yyval.cdata.d[0]=yyvsp[-6].fval; yyval.cdata.d[1]=yyvsp[-5].fval; yyval.cdata.i[1]=yyvsp[-4].ival; yyval.cdata.i[2]=yyvsp[-3].ival;
          yyval.cdata.ip  =yyvsp[-1].ilist.ip; yyval.cdata.ipsize= yyvsp[-1].ilist.ipsize;;}
    break;

  case 41:
#line 206 "Relparser.y"
    { yyval.cdata.typ =13; yyval.cdata.i[0]=yyvsp[-5].ival; yyval.cdata.i[1]=yyvsp[-4].ival; yyval.cdata.d[0]=yyvsp[-3].fval; 
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 42:
#line 209 "Relparser.y"
    { yyval.cdata.typ =13; yyval.cdata.i[0]=1; yyval.cdata.i[1]=yyvsp[-4].ival; yyval.cdata.d[0]=yyvsp[-3].fval;
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 43:
#line 212 "Relparser.y"
    { yyval.cdata.typ =20; yyval.cdata.i[0]=yyvsp[-6].ival; yyval.cdata.i[1]=yyvsp[-5].ival; yyval.cdata.d[0]=yyvsp[-4].fval; yyval.cdata.i[2]=yyvsp[-3].ival;
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 44:
#line 215 "Relparser.y"
    { yyval.cdata.typ =20; yyval.cdata.i[0]=1; yyval.cdata.i[1]=yyvsp[-5].ival; yyval.cdata.d[0]=yyvsp[-4].fval; yyval.cdata.i[2]=yyvsp[-3].ival;
          yyval.cdata.dl  =yyvsp[-1].dpllist.dl; yyval.cdata.dlsize=yyvsp[-1].dpllist.dlsize; ;}
    break;

  case 45:
#line 218 "Relparser.y"
    { yyval.cdata.typ =  15; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 46:
#line 220 "Relparser.y"
    { yyval.cdata.typ =  16; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 47:
#line 222 "Relparser.y"
    { yyval.cdata.typ =  16; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 48:
#line 224 "Relparser.y"
    { yyval.cdata.typ =  17; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 49:
#line 226 "Relparser.y"
    { yyval.cdata.typ =  18; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 50:
#line 228 "Relparser.y"
    { yyval.cdata.typ =  19; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 51:
#line 230 "Relparser.y"
    { yyval.cdata.typ =  21; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 52:
#line 232 "Relparser.y"
    { yyval.cdata.typ = 22; ;}
    break;

  case 53:
#line 234 "Relparser.y"
    { yyval.cdata.typ = 23; yyval.cdata.i[0]=yyvsp[-1].ival-1; yyval.cdata.i[1]=yyvsp[0].ival; ;}
    break;

  case 54:
#line 236 "Relparser.y"
    { yyval.cdata.typ =  24; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 55:
#line 238 "Relparser.y"
    { yyval.cdata.typ =  25; yyval.cdata.i[0]=yyvsp[0].ival-1; ;}
    break;

  case 56:
#line 240 "Relparser.y"
    { yyval.cdata.typ =  26; ;}
    break;

  case 57:
#line 242 "Relparser.y"
    { yyval.cdata.typ = 27;  ;}
    break;

  case 58:
#line 244 "Relparser.y"
    { yyval.cdata.typ = 28; yyval.cdata.i[0] = yyvsp[0].ival-1; ;}
    break;

  case 59:
#line 248 "Relparser.y"
    { yyval.ival = 0;  ;}
    break;

  case 60:
#line 249 "Relparser.y"
    { yyval.ival = 1;  ;}
    break;

  case 61:
#line 250 "Relparser.y"
    { yyval.ival = 2;  ;}
    break;

  case 62:
#line 251 "Relparser.y"
    { yyval.ival = 10; ;}
    break;

  case 63:
#line 252 "Relparser.y"
    { yyval.ival = 11; ;}
    break;

  case 64:
#line 253 "Relparser.y"
    { yyval.ival = 12; ;}
    break;

  case 65:
#line 254 "Relparser.y"
    { yyval.ival = 20; ;}
    break;

  case 66:
#line 255 "Relparser.y"
    { yyval.ival = 21; ;}
    break;

  case 67:
#line 256 "Relparser.y"
    { yyval.ival = 22; ;}
    break;

  case 68:
#line 257 "Relparser.y"
    { yyval.ival = 30; ;}
    break;

  case 69:
#line 258 "Relparser.y"
    { yyval.ival = 31; ;}
    break;

  case 70:
#line 259 "Relparser.y"
    { yyval.ival = 32; ;}
    break;

  case 71:
#line 260 "Relparser.y"
    { yyval.ival = 40; ;}
    break;

  case 72:
#line 261 "Relparser.y"
    { yyval.ival = 41; ;}
    break;

  case 73:
#line 262 "Relparser.y"
    { yyval.ival = 42; ;}
    break;

  case 74:
#line 263 "Relparser.y"
    { yyval.ival = 50; ;}
    break;

  case 75:
#line 264 "Relparser.y"
    { yyval.ival = 51; ;}
    break;

  case 76:
#line 265 "Relparser.y"
    { yyval.ival = 52; ;}
    break;

  case 77:
#line 269 "Relparser.y"
    { yyval.ival = yyvsp[0].ival-1; ;}
    break;

  case 78:
#line 270 "Relparser.y"
    { yyval.ival = 0;   ;}
    break;

  case 79:
#line 271 "Relparser.y"
    { yyval.ival = 1;   ;}
    break;

  case 80:
#line 272 "Relparser.y"
    { yyval.ival = 2;   ;}
    break;

  case 81:
#line 273 "Relparser.y"
    { yyval.ival = 3;   ;}
    break;

  case 82:
#line 274 "Relparser.y"
    { yyval.ival = 4;   ;}
    break;

  case 83:
#line 275 "Relparser.y"
    { yyval.ival = 5;   ;}
    break;

  case 84:
#line 280 "Relparser.y"
    { yyval.fval = yyvsp[0].ival-1; ;}
    break;

  case 85:
#line 282 "Relparser.y"
    { yyval.fval = yyvsp[0].fval-1; ;}
    break;

  case 86:
#line 287 "Relparser.y"
    { int num=yyvsp[-1].ival; /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addFailureCriteria(num,yyvsp[0].fall,0,0); ;}
    break;

  case 87:
#line 289 "Relparser.y"
    { int num=yyvsp[-3].ival; /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addFailureCriteria(num,yyvsp[0].fall,1,yyvsp[-1].fval); ;}
    break;

  case 88:
#line 291 "Relparser.y"
    { int num=yyvsp[-1].ival; /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addFailureCriteria(num,yyvsp[0].fall,0,0); ;}
    break;

  case 89:
#line 293 "Relparser.y"
    { int num=yyvsp[-3].ival; /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addFailureCriteria(num,yyvsp[0].fall,1,yyvsp[-1].fval); ;}
    break;

  case 90:
#line 298 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addRndvar(yyvsp[0].absvar); ;}
    break;

  case 91:
#line 300 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addRndvar(yyvsp[0].absvar); ;}
    break;

  case 92:
#line 305 "Relparser.y"
    { yyval.absvar.num=yyvsp[-6].ival; yyval.absvar.dist=yyvsp[-5].ival; yyval.absvar.mean=yyvsp[-4].fval; yyval.absvar.sdev=yyvsp[-3].fval; yyval.absvar.low=yyvsp[-2].fval; yyval.absvar.upp=yyvsp[-1].fval; 
          yyval.absvar.igen=0; ;}
    break;

  case 93:
#line 308 "Relparser.y"
    { yyval.absvar.num=yyvsp[-4].ival; yyval.absvar.dist=yyvsp[-3].ival; yyval.absvar.mean=yyvsp[-2].fval; yyval.absvar.sdev=yyvsp[-1].fval; yyval.absvar.low=-1.0e10; yyval.absvar.upp=1.0e10; 
	  yyval.absvar.igen=0; ;}
    break;

  case 94:
#line 311 "Relparser.y"
    { yyval.absvar.num=yyvsp[-7].ival; yyval.absvar.dist=yyvsp[-6].ival; yyval.absvar.mean=yyvsp[-5].fval; yyval.absvar.sdev=yyvsp[-4].fval; yyval.absvar.low=yyvsp[-3].fval; yyval.absvar.upp=yyvsp[-2].fval; 
	  yyval.absvar.igen=1; yyval.absvar.gena=yyvsp[-1].gen.a; yyval.absvar.genb=yyvsp[-1].gen.b ; yyval.absvar.gens=yyvsp[-1].gen.s; ;}
    break;

  case 95:
#line 314 "Relparser.y"
    { yyval.absvar.num=yyvsp[-5].ival; yyval.absvar.dist=yyvsp[-4].ival; yyval.absvar.mean=yyvsp[-3].fval ; yyval.absvar.sdev=yyvsp[-2].fval; yyval.absvar.low=-1.0e10; yyval.absvar.upp=1.0e10; 
	  yyval.absvar.igen=1; yyval.absvar.gena=yyvsp[-1].gen.a; yyval.absvar.genb=yyvsp[-1].gen.b ; yyval.absvar.gens=yyvsp[-1].gen.s; ;}
    break;

  case 96:
#line 320 "Relparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall,1); ;}
    break;

  case 97:
#line 323 "Relparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall,1); ;}
    break;

  case 98:
#line 326 "Relparser.y"
    { int num=yyvsp[-3].ival; int loc1=yyvsp[-1].ival-1; int loc2=0;
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addStcvar(num,yyvsp[-2].ival,loc1,loc2,yyvsp[0].fall,1); ;}
    break;

  case 99:
#line 329 "Relparser.y"
    { int num=yyvsp[-3].ival; int loc1=yyvsp[-1].ival-1; int loc2=0;
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addStcvar(num,yyvsp[-2].ival,loc1,loc2,yyvsp[0].fall,1); ;}
    break;

  case 100:
#line 332 "Relparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall,1); ;}
    break;

  case 101:
#line 335 "Relparser.y"
    { int num=yyvsp[-4].ival; int loc1=yyvsp[-2].ival-1; int loc2=yyvsp[-1].ival;
	  /*dynamic_cast<Domain_opt*>(domain)->*/relpro->addStcvar(num,yyvsp[-3].ival,loc1,loc2,yyvsp[0].fall,1); ;}
    break;

  case 102:
#line 340 "Relparser.y"
    { yyval.ival = yyvsp[0].ival-1; ;}
    break;

  case 103:
#line 341 "Relparser.y"
    { yyval.ival = 0; ;}
    break;

  case 104:
#line 342 "Relparser.y"
    { yyval.ival = 1; ;}
    break;

  case 105:
#line 343 "Relparser.y"
    { yyval.ival = 3; ;}
    break;

  case 106:
#line 344 "Relparser.y"
    { yyval.ival = 4; ;}
    break;

  case 107:
#line 345 "Relparser.y"
    { yyval.ival =10; ;}
    break;

  case 108:
#line 346 "Relparser.y"
    { yyval.ival = 5; ;}
    break;

  case 109:
#line 347 "Relparser.y"
    { yyval.ival = 5; ;}
    break;

  case 110:
#line 348 "Relparser.y"
    { yyval.ival = 2; ;}
    break;

  case 111:
#line 349 "Relparser.y"
    { yyval.ival = 6; ;}
    break;

  case 112:
#line 350 "Relparser.y"
    { yyval.ival = 0; ;}
    break;

  case 113:
#line 351 "Relparser.y"
    { yyval.ival = 1; ;}
    break;

  case 114:
#line 352 "Relparser.y"
    { yyval.ival = 2; ;}
    break;

  case 115:
#line 353 "Relparser.y"
    { yyval.ival = 3; ;}
    break;

  case 116:
#line 354 "Relparser.y"
    { yyval.ival = 4; ;}
    break;

  case 117:
#line 355 "Relparser.y"
    { yyval.ival = 5; ;}
    break;

  case 118:
#line 356 "Relparser.y"
    { yyval.ival = 100000 * yyvsp[-1].ival + yyvsp[0].ival; ;}
    break;

  case 119:
#line 360 "Relparser.y"
    { yyval.ival = yyvsp[0].ival-1; ;}
    break;

  case 120:
#line 361 "Relparser.y"
    { yyval.ival = 1; ;}
    break;

  case 121:
#line 362 "Relparser.y"
    { yyval.ival = 2; ;}
    break;

  case 122:
#line 363 "Relparser.y"
    { yyval.ival = 4; ;}
    break;

  case 123:
#line 364 "Relparser.y"
    { yyval.ival = 7; ;}
    break;

  case 124:
#line 365 "Relparser.y"
    { yyval.ival = 8; ;}
    break;

  case 125:
#line 366 "Relparser.y"
    { yyval.ival = 9; ;}
    break;

  case 126:
#line 367 "Relparser.y"
    { yyval.ival = 10; ;}
    break;

  case 127:
#line 368 "Relparser.y"
    { yyval.ival = 11; ;}
    break;

  case 128:
#line 372 "Relparser.y"
    { yyval.ival = 1; ;}
    break;

  case 129:
#line 373 "Relparser.y"
    { yyval.ival = 2; ;}
    break;

  case 130:
#line 374 "Relparser.y"
    { yyval.ival = 3; ;}
    break;

  case 131:
#line 375 "Relparser.y"
    { yyval.ival = 4 + yyvsp[0].ival; ;}
    break;

  case 132:
#line 379 "Relparser.y"
    { yyval.ival = 100; ;}
    break;

  case 133:
#line 384 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structrel->addAnalysisData(yyvsp[0].anadat);  ;}
    break;

  case 134:
#line 386 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structrel->addAnalysisData(yyvsp[0].anadat); ;}
    break;

  case 135:
#line 391 "Relparser.y"
    { yyval.anadat.typ=0; yyval.anadat.ival[0]=yyvsp[-1].ival; ;}
    break;

  case 136:
#line 393 "Relparser.y"
    { yyval.anadat.typ=1; yyval.anadat.ival[0]=yyvsp[-1].ival; ;}
    break;

  case 137:
#line 395 "Relparser.y"
    { yyval.anadat.typ=1; yyval.anadat.ival[0]=yyvsp[-2].ival; yyval.anadat.ival[1]=yyvsp[-1].ival; ;}
    break;

  case 138:
#line 397 "Relparser.y"
    { yyval.anadat.typ=2; yyval.anadat.ival[0]=yyvsp[-5].ival; yyval.anadat.ival[1]=yyvsp[-4].ival; 
	            yyval.anadat.rval[0]=yyvsp[-3].fval; yyval.anadat.rval[1]=yyvsp[-2].fval;
                    yyval.anadat.ival[2]=yyvsp[-1].ival; ;}
    break;

  case 139:
#line 403 "Relparser.y"
    { yyval.ival=yyvsp[0].ival; ;}
    break;

  case 140:
#line 404 "Relparser.y"
    { yyval.ival=0;  ;}
    break;

  case 141:
#line 405 "Relparser.y"
    { yyval.ival=1;  ;}
    break;

  case 142:
#line 409 "Relparser.y"
    { yyval.ival=yyvsp[0].ival; ;}
    break;

  case 143:
#line 410 "Relparser.y"
    { yyval.ival=-1; ;}
    break;

  case 144:
#line 415 "Relparser.y"
    { int num=yyvsp[-1].ival; /*dynamic_cast<Domain_opt*>(domain)->*/relpro->relsol->addSolver(num,yyvsp[0].solv.typ,*(yyvsp[0].solv.param),*(yyvsp[0].solv.pgrad)); ;}
    break;

  case 145:
#line 417 "Relparser.y"
    { int num=yyvsp[-1].ival; /*dynamic_cast<Domain_opt*>(domain)->*/relpro->relsol->addSolver(num,yyvsp[0].solv.typ,*(yyvsp[0].solv.param),*(yyvsp[0].solv.pgrad)); ;}
    break;

  case 146:
#line 422 "Relparser.y"
    { yyval.solv.typ=0;   yyval.solv.param=&yyvsp[0].nlp; yyval.solv.pgrad=0; ;}
    break;

  case 147:
#line 424 "Relparser.y"
    { yyval.solv.typ=100; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 148:
#line 426 "Relparser.y"
    { yyval.solv.typ=101; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 149:
#line 428 "Relparser.y"
    { yyval.solv.typ=102; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 150:
#line 430 "Relparser.y"
    { yyval.solv.typ=103; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 151:
#line 432 "Relparser.y"
    { yyval.solv.typ=104; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 152:
#line 434 "Relparser.y"
    { yyval.solv.typ=105; yyval.solv.param=&yyvsp[-1].nlp; yyval.solv.pgrad=&yyvsp[0].grad; ;}
    break;

  case 153:
#line 439 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 154:
#line 442 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 155:
#line 445 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 156:
#line 448 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 157:
#line 454 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 158:
#line 457 "Relparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 159:
#line 459 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 160:
#line 462 "Relparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 161:
#line 467 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 162:
#line 470 "Relparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 163:
#line 472 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 164:
#line 475 "Relparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 165:
#line 480 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 166:
#line 483 "Relparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 167:
#line 485 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 168:
#line 488 "Relparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 169:
#line 493 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 170:
#line 496 "Relparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 171:
#line 498 "Relparser.y"
    { yyval.nlp.initialize(); 
          yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 172:
#line 501 "Relparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 173:
#line 506 "Relparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 174:
#line 508 "Relparser.y"
    { yyval.nlp.ival[yyvsp[-2].ival] = yyvsp[-1].ival; yyval.nlp.iflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 175:
#line 510 "Relparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 176:
#line 512 "Relparser.y"
    { yyval.nlp.rval[yyvsp[-2].ival] = yyvsp[-1].fval; yyval.nlp.rflag[yyvsp[-2].ival] = 1; ;}
    break;

  case 177:
#line 517 "Relparser.y"
    { yyval.grad.typ=yyvsp[-1].ival; yyval.grad.filter=0; ;}
    break;

  case 178:
#line 519 "Relparser.y"
    { yyval.grad.mth=yyvsp[-1].ival; ;}
    break;

  case 179:
#line 521 "Relparser.y"
    { yyval.grad.epstyp=0 ; yyval.grad.epsval=yyvsp[-1].fval; ;}
    break;

  case 180:
#line 523 "Relparser.y"
    { yyval.grad.epstyp=1 ; yyval.grad.epsval=yyvsp[-1].fval; ;}
    break;

  case 181:
#line 525 "Relparser.y"
    { dynamic_cast<Domain_opt*>(domain)->setSemiSAflag(1); ;}
    break;

  case 182:
#line 527 "Relparser.y"
    { yyval.grad.filter=1;  yyval.grad.filterTyp=yyvsp[-9].ival; yyval.grad.filterScale=yyvsp[-8].ival;  
          yyval.grad.radius=yyvsp[-7].fval; yyval.grad.maxCount=yyvsp[-6].fval;  yyval.grad.minExp=yyvsp[-5].fval;  yyval.grad.maxExp=yyvsp[-4].fval;
          yyval.grad.numFilCrit=yyvsp[-2].ilist.ipsize;      yyval.grad.filcritList=yyvsp[-2].ilist.ip; 
	  yyval.grad.numGroups=0;              
	  yyval.grad.numFilGrps=0;               yyval.grad.filGrpsList=0; ;}
    break;

  case 183:
#line 533 "Relparser.y"
    { yyval.grad.filter=1;  yyval.grad.filterTyp=yyvsp[-13].ival; yyval.grad.filterScale=yyvsp[-12].ival;  
          yyval.grad.radius=yyvsp[-11].fval; yyval.grad.maxCount=yyvsp[-10].fval;  yyval.grad.minExp=yyvsp[-9].fval;  yyval.grad.maxExp=yyvsp[-8].fval;
          yyval.grad.numFilCrit=yyvsp[-6].ilist.ipsize;      yyval.grad.filcritList=yyvsp[-6].ilist.ip; 
	  yyval.grad.numGroups=yyvsp[-4].ival;              
	  yyval.grad.numFilGrps=yyvsp[-2].ilist.ipsize;      yyval.grad.filGrpsList=yyvsp[-2].ilist.ip; ;}
    break;

  case 184:
#line 541 "Relparser.y"
    { yyval.ival=0; ;}
    break;

  case 185:
#line 542 "Relparser.y"
    { yyval.ival=1; ;}
    break;

  case 186:
#line 546 "Relparser.y"
    { yyval.ival=0; ;}
    break;

  case 187:
#line 547 "Relparser.y"
    { yyval.ival=1; ;}
    break;

  case 188:
#line 548 "Relparser.y"
    { yyval.ival=2; ;}
    break;

  case 189:
#line 549 "Relparser.y"
    { yyval.ival=3; ;}
    break;

  case 190:
#line 554 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structrel->addOptInf(yyvsp[0].sevdata.var,yyvsp[0].sevdata.ele1,yyvsp[0].sevdata.ele2,yyvsp[0].sevdata.typ); ;}
    break;

  case 191:
#line 556 "Relparser.y"
    { /*dynamic_cast<Domain_opt*>(domain)->*/structrel->addOptInf(yyvsp[0].sevdata.var,yyvsp[0].sevdata.ele1,yyvsp[0].sevdata.ele2,yyvsp[0].sevdata.typ); ;}
    break;

  case 192:
#line 561 "Relparser.y"
    { yyval.sevdata.var=yyvsp[-3].ival-1; yyval.sevdata.ele1=yyvsp[-2].ival-1; yyval.sevdata.ele2=-1; yyval.sevdata.typ=yyvsp[-1].ival; ;}
    break;

  case 193:
#line 563 "Relparser.y"
    { yyval.sevdata.var=yyvsp[-4].ival-1; yyval.sevdata.ele1=yyvsp[-3].ival-1; yyval.sevdata.ele2=yyvsp[-2].ival-1; yyval.sevdata.typ=yyvsp[-1].ival; ;}
    break;

  case 194:
#line 568 "Relparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-4].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 195:
#line 570 "Relparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-5].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 196:
#line 572 "Relparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-6].ival; yyval.fall.fdata=yyvsp[-3].func; ;}
    break;

  case 197:
#line 574 "Relparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-5].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 198:
#line 576 "Relparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-6].ival; yyval.fall.fdata=yyvsp[-2].func; ;}
    break;

  case 199:
#line 578 "Relparser.y"
    { initFuncall(yyval.fall); yyval.fall.typ=yyvsp[-7].ival; yyval.fall.fdata=yyvsp[-3].func; ;}
    break;

  case 200:
#line 582 "Relparser.y"
    { yyval.ival=0; ;}
    break;

  case 201:
#line 583 "Relparser.y"
    { yyval.ival=1; ;}
    break;

  case 202:
#line 584 "Relparser.y"
    { yyval.ival=2; ;}
    break;

  case 203:
#line 585 "Relparser.y"
    { yyval.ival=3; ;}
    break;

  case 204:
#line 586 "Relparser.y"
    { yyval.ival=4; ;}
    break;

  case 205:
#line 591 "Relparser.y"
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

  case 206:
#line 602 "Relparser.y"
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

  case 207:
#line 613 "Relparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
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

  case 208:
#line 629 "Relparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
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

  case 209:
#line 645 "Relparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=1; yyval.func->numgen=0; 
	  yyval.func->numstat[0]=-1;  
	  yyval.func->oprtyp[0]=yyvsp[0].sum.oprtyp; yyval.func->oprnum[0]=yyvsp[0].sum.oprnum;
	  yyval.func->a[0]=yyvsp[0].sum.a; yyval.func->p[0]=yyvsp[0].sum.p; yyval.func->b[0]=yyvsp[0].sum.b; 
	  if ( yyvsp[0].sum.igen ) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"no number for summands assigned\n");
	    fprintf (stderr," *** This is required for generating summands\n");
	    return -1;
	  }
	;}
    break;

  case 210:
#line 657 "Relparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=1; yyval.func->numgen=0; 
	  yyval.func->numstat[0]=-1;  
	  yyval.func->oprtyp[0]=yyvsp[-1].sum.oprtyp; yyval.func->oprnum[0]=yyvsp[-1].sum.oprnum;
	  yyval.func->a[0]=yyvsp[-1].sum.a; yyval.func->p[0]=yyvsp[-1].sum.p; yyval.func->b[0]=yyvsp[-1].sum.b; 
	  if ( yyvsp[-1].sum.igen ) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"no number for summands assigned\n");
	    fprintf (stderr," *** This is required for generating summands!\n");
	    return -1;
	  }
	;}
    break;

  case 211:
#line 669 "Relparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	  yyval.func->numstat[num]=-1;  
	  yyval.func->oprtyp[num]=yyvsp[0].sum.oprtyp; yyval.func->oprnum[num]=yyvsp[0].sum.oprnum; 
	  yyval.func->a[num]=yyvsp[0].sum.a; yyval.func->p[num]=yyvsp[0].sum.p; yyval.func->b[num]=yyvsp[0].sum.b;
	  if ( yyvsp[0].sum.igen ) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"no number for summands assigned\n");
	    fprintf (stderr," *** This is required for generating summands!\n");
	    return -1;
	  }
	;}
    break;

  case 212:
#line 686 "Relparser.y"
    { int num=yyval.func->numopr; yyval.func->numopr++; 
	  if (yyval.func->numopr > MAXOPR) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	  yyval.func->numstat[num]=-1;  
	  yyval.func->oprtyp[num]=yyvsp[-1].sum.oprtyp; yyval.func->oprnum[num]=yyvsp[-1].sum.oprnum; 
	  yyval.func->a[num]=yyvsp[-1].sum.a; yyval.func->p[num]=yyvsp[-1].sum.p; yyval.func->b[num]=yyvsp[-1].sum.b;
	  if ( yyvsp[-1].sum.igen ) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"number of defined summands > MAXOPR");
	    return -1;
	  }
	;}
    break;

  case 213:
#line 702 "Relparser.y"
    { yyval.func=buildFuncdata(); yyval.func->numopr=0; yyval.func->numgen=0; int num=yyvsp[0].fdef.num;
	  if (num > MAXOPR) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"number of subfunctions > MAXOPR");
	    return -1;
	  }
	  yyval.func->subfunc[num]=yyvsp[0].fdef.falldata;;}
    break;

  case 214:
#line 710 "Relparser.y"
    { int num=yyvsp[0].fdef.num;
	  if (num > MAXOPR) {
	    fprintf (stderr," *** Error reading reliability analysis input file: ");
	    fprintf (stderr,"number of subfunctions > MAXOPR");
	    return -1;
	  }
	  yyval.func->subfunc[num]=yyvsp[0].fdef.falldata; ;}
    break;

  case 215:
#line 721 "Relparser.y"
    { yyval.ival = yyvsp[-1].ival; ;}
    break;

  case 216:
#line 726 "Relparser.y"
    { yyval.sum.a=yyvsp[-6].fval; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 217:
#line 729 "Relparser.y"
    { yyval.sum.a=yyvsp[-6].fval; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 218:
#line 732 "Relparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=yyvsp[0].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 219:
#line 735 "Relparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 220:
#line 738 "Relparser.y"
    { yyval.sum.a=yyvsp[-4].fval; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 221:
#line 741 "Relparser.y"
    { yyval.sum.a=yyvsp[-2].fval; yyval.sum.oprtyp=yyvsp[0].opr.oprtyp; yyval.sum.oprnum=yyvsp[0].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 222:
#line 744 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 223:
#line 747 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-4].opr.oprtyp; yyval.sum.oprnum=yyvsp[-4].opr.oprnum; yyval.sum.p=yyvsp[-2].fval; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 224:
#line 750 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=yyvsp[0].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 225:
#line 753 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 226:
#line 756 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-2].opr.oprtyp; yyval.sum.oprnum=yyvsp[-2].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[0].fval;  
	  yyval.sum.igen=0; ;}
    break;

  case 227:
#line 759 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[0].opr.oprtyp; yyval.sum.oprnum=yyvsp[0].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=0; ;}
    break;

  case 228:
#line 762 "Relparser.y"
    { yyval.sum.a=yyvsp[-7].fval; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;  ;}
    break;

  case 229:
#line 765 "Relparser.y"
    { yyval.sum.a=yyvsp[-7].fval; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 230:
#line 768 "Relparser.y"
    { yyval.sum.a=yyvsp[-5].fval; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=yyvsp[-1].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 231:
#line 771 "Relparser.y"
    { yyval.sum.a=yyvsp[-5].fval; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 232:
#line 774 "Relparser.y"
    { yyval.sum.a=yyvsp[-5].fval; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 233:
#line 777 "Relparser.y"
    { yyval.sum.a=yyvsp[-3].fval; yyval.sum.oprtyp=yyvsp[-1].opr.oprtyp; yyval.sum.oprnum=yyvsp[-1].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 234:
#line 780 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;    ;}
    break;

  case 235:
#line 783 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-5].opr.oprtyp; yyval.sum.oprnum=yyvsp[-5].opr.oprnum; yyval.sum.p=yyvsp[-3].fval; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 236:
#line 786 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=yyvsp[-1].fval; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 237:
#line 789 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 238:
#line 792 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-3].opr.oprtyp; yyval.sum.oprnum=yyvsp[-3].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=-yyvsp[-1].fval;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;   ;}
    break;

  case 239:
#line 795 "Relparser.y"
    { yyval.sum.a=1.0; yyval.sum.oprtyp=yyvsp[-1].opr.oprtyp; yyval.sum.oprnum=yyvsp[-1].opr.oprnum; yyval.sum.p=1.0; yyval.sum.b=0.0;  
	  yyval.sum.igen=1; yyval.sum.gena=yyvsp[0].gen.a; yyval.sum.genb=yyvsp[0].gen.b; yyval.sum.gens=yyvsp[0].gen.s;  ;}
    break;

  case 240:
#line 801 "Relparser.y"
    { yyval.opr.oprtyp=0; yyval.opr.oprnum=yyvsp[-1].ival-1; ;}
    break;

  case 241:
#line 803 "Relparser.y"
    { yyval.opr.oprtyp=1; yyval.opr.oprnum=yyvsp[-1].ival-1; ;}
    break;

  case 242:
#line 805 "Relparser.y"
    { yyval.opr.oprtyp=2; yyval.opr.oprnum=yyvsp[-1].ival-1; ;}
    break;

  case 243:
#line 810 "Relparser.y"
    { yyval.fdef.num=yyvsp[-3].ival-1; yyval.fdef.falldata=yyvsp[0].fall; ;}
    break;

  case 244:
#line 815 "Relparser.y"
    { yyval.gen.a=yyvsp[-5].ival-1 ; yyval.gen.b = yyvsp[-3].ival-1 ; yyval.gen.s = yyvsp[-1].ival ; ;}
    break;

  case 245:
#line 817 "Relparser.y"
    { yyval.gen.a=yyvsp[-3].ival-1 ; yyval.gen.b = yyvsp[-1].ival-1 ; yyval.gen.s = 1 ; ;}
    break;

  case 246:
#line 822 "Relparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival-1); ;}
    break;

  case 247:
#line 824 "Relparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival-1); ;}
    break;

  case 248:
#line 826 "Relparser.y"
    { yyval.ilist.add(yyvsp[0].ival-1); ;}
    break;

  case 249:
#line 828 "Relparser.y"
    { ;}
    break;

  case 250:
#line 833 "Relparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival); ;}
    break;

  case 251:
#line 835 "Relparser.y"
    { yyval.ilist.ipsize=0; yyval.ilist.mxsize=0; yyval.ilist.add(yyvsp[0].ival); ;}
    break;

  case 252:
#line 837 "Relparser.y"
    { yyval.ilist.add(yyvsp[0].ival); ;}
    break;

  case 253:
#line 839 "Relparser.y"
    { ;}
    break;

  case 254:
#line 844 "Relparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 255:
#line 846 "Relparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 256:
#line 848 "Relparser.y"
    { yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 257:
#line 850 "Relparser.y"
    { ;}
    break;

  case 258:
#line 855 "Relparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 259:
#line 857 "Relparser.y"
    { yyval.dpllist.dlsize=0; yyval.dpllist.mxsize=0; yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 260:
#line 859 "Relparser.y"
    { yyval.dpllist.add(yyvsp[-3].ival-1,yyvsp[-2].ival,yyvsp[-1].fval,yyvsp[0].fval); ;}
    break;

  case 261:
#line 861 "Relparser.y"
    { ;}
    break;

  case 262:
#line 866 "Relparser.y"
    { yyval.ival = yyvsp[0].ival; ;}
    break;

  case 263:
#line 871 "Relparser.y"
    { yyval.fval = yyvsp[0].ival; ;}
    break;

  case 264:
#line 873 "Relparser.y"
    { yyval.fval = yyvsp[0].fval; ;}
    break;


    }

/* Line 1000 of yacc.c.  */
#line 3358 "Relparser.C"

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


#line 875 "Relparser.y"


#endif

