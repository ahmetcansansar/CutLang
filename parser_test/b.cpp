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
     DEF = 258,
     CMD = 259,
     ELE = 260,
     MUO = 261,
     LEP = 262,
     PHO = 263,
     JET = 264,
     BJET = 265,
     QGJET = 266,
     NUMET = 267,
     METLV = 268,
     PHI = 269,
     ETA = 270,
     ABSETA = 271,
     PT = 272,
     PZ = 273,
     NBF = 274,
     DR = 275,
     DPHI = 276,
     NELE = 277,
     NMUO = 278,
     NLEP = 279,
     NPHO = 280,
     NJET = 281,
     NBJET = 282,
     NQGJET = 283,
     HT = 284,
     METMWT = 285,
     MWT = 286,
     MET = 287,
     ALL = 288,
     LEPSF = 289,
     FILLHISTOS = 290,
     NB = 291,
     ID = 292,
     SIN = 293,
     COS = 294,
     TAN = 295,
     INT = 296,
     OR = 297,
     AND = 298,
     LT = 299,
     GT = 300,
     LE = 301,
     GE = 302,
     EQ = 303,
     NE = 304,
     IRG = 305,
     ERG = 306,
     Unary = 307
   };
#endif
/* Tokens.  */
#define DEF 258
#define CMD 259
#define ELE 260
#define MUO 261
#define LEP 262
#define PHO 263
#define JET 264
#define BJET 265
#define QGJET 266
#define NUMET 267
#define METLV 268
#define PHI 269
#define ETA 270
#define ABSETA 271
#define PT 272
#define PZ 273
#define NBF 274
#define DR 275
#define DPHI 276
#define NELE 277
#define NMUO 278
#define NLEP 279
#define NPHO 280
#define NJET 281
#define NBJET 282
#define NQGJET 283
#define HT 284
#define METMWT 285
#define MWT 286
#define MET 287
#define ALL 288
#define LEPSF 289
#define FILLHISTOS 290
#define NB 291
#define ID 292
#define SIN 293
#define COS 294
#define TAN 295
#define INT 296
#define OR 297
#define AND 298
#define LT 299
#define GT 300
#define LE 301
#define GE 302
#define EQ 303
#define NE 304
#define IRG 305
#define ERG 306
#define Unary 307




/* Copy the first part of user declarations.  */
#line 1 "parse.y"
 
#include <math.h>
#include "stdlib.h"
#include <iostream>
#include <string>
#include <map>
#include <iterator>
extern int yylex();
extern int yyparse();
void yyerror(const char *s) { std::cout << s << std::endl; } 
int cutcount;
using namespace std;
string tmp;
int pnum;
map<string,string> vars;
map<string,string> parts;
map<int,string> cuts;


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
#line 19 "parse.y"
{
	double real;
	char* s;
}
/* Line 193 of yacc.c.  */
#line 224 "b.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 237 "b.cpp"

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
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   331

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  70
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  16
/* YYNRULES -- Number of rules.  */
#define YYNRULES  68
/* YYNRULES -- Number of states.  */
#define YYNSTATES  137

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   307

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      68,    69,    54,    52,    65,    53,     2,    55,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    58,     2,
       2,     2,     2,    67,     2,     2,     2,     2,     2,    64,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      63,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    57,    66,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    61,
       2,     2,     2,    62,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    59,     2,    60,     2,     2,     2,     2,
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
      45,    46,    47,    48,    49,    50,    51,    56
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     6,     9,    10,    15,    20,    25,    30,
      35,    40,    45,    50,    55,    60,    65,    70,    73,    76,
      77,    84,    87,    89,    93,    97,   101,   105,   109,   113,
     117,   121,   125,   127,   130,   132,   135,   136,   139,   142,
     145,   151,   153,   155,   157,   161,   165,   169,   173,   177,
     181,   186,   191,   195,   199,   203,   207,   211,   215,   219,
     223,   226,   231,   236,   241,   245,   247,   249,   251
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      71,     0,    -1,    72,    80,    -1,    72,    73,    -1,    -1,
       3,    37,    58,    77,    -1,     3,    37,    58,    85,    -1,
      59,    77,    60,    61,    -1,    59,    77,    60,    62,    -1,
      59,    77,    60,    63,    -1,    59,    77,    60,    64,    -1,
      59,    77,    60,    14,    -1,    59,    77,    60,    15,    -1,
      59,    77,    60,    16,    -1,    59,    77,    60,    17,    -1,
      59,    77,    60,    18,    -1,    59,    77,    60,    19,    -1,
      75,    20,    -1,    75,    21,    -1,    -1,    59,    77,    76,
      65,    77,    60,    -1,    77,    78,    -1,    78,    -1,     5,
      66,    79,    -1,     6,    66,    79,    -1,     7,    66,    79,
      -1,     8,    66,    79,    -1,     9,    66,    79,    -1,    10,
      66,    79,    -1,    11,    66,    79,    -1,    12,    66,    79,
      -1,    13,    66,    79,    -1,    37,    -1,    53,    41,    -1,
      41,    -1,    80,    81,    -1,    -1,     4,    84,    -1,     4,
      33,    -1,     4,    82,    -1,    84,    67,    83,    58,    83,
      -1,    84,    -1,    33,    -1,    82,    -1,    85,    44,    85,
      -1,    85,    45,    85,    -1,    85,    46,    85,    -1,    85,
      47,    85,    -1,    85,    48,    85,    -1,    85,    49,    85,
      -1,    85,    50,    85,    85,    -1,    85,    51,    85,    85,
      -1,    84,    43,    84,    -1,    84,    42,    84,    -1,    68,
      84,    69,    -1,    85,    52,    85,    -1,    85,    53,    85,
      -1,    85,    54,    85,    -1,    85,    55,    85,    -1,    85,
      57,    85,    -1,    53,    85,    -1,    39,    68,    85,    69,
      -1,    38,    68,    85,    69,    -1,    40,    68,    85,    69,
      -1,    68,    85,    69,    -1,    36,    -1,    41,    -1,    74,
      -1,    37,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    44,    44,    46,    47,    49,    67,    86,    92,    98,
     104,   110,   116,   122,   128,   134,   140,   146,   152,   159,
     159,   167,   175,   187,   194,   198,   202,   206,   210,   214,
     218,   222,   226,   243,   244,   246,   247,   249,   250,   251,
     253,   258,   259,   261,   263,   267,   271,   275,   279,   283,
     287,   291,   295,   299,   303,   308,   312,   316,   320,   324,
     328,   332,   336,   340,   344,   348,   349,   350,   352
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "DEF", "CMD", "ELE", "MUO", "LEP", "PHO",
  "JET", "BJET", "QGJET", "NUMET", "METLV", "PHI", "ETA", "ABSETA", "PT",
  "PZ", "NBF", "DR", "DPHI", "NELE", "NMUO", "NLEP", "NPHO", "NJET",
  "NBJET", "NQGJET", "HT", "METMWT", "MWT", "MET", "ALL", "LEPSF",
  "FILLHISTOS", "NB", "ID", "SIN", "COS", "TAN", "INT", "OR", "AND", "LT",
  "GT", "LE", "GE", "EQ", "NE", "IRG", "ERG", "'+'", "'-'", "'*'", "'/'",
  "Unary", "'^'", "':'", "'{'", "'}'", "'m'", "'q'", "'P'", "'E'", "','",
  "'_'", "'?'", "'('", "')'", "$accept", "input", "definitions",
  "definition", "function", "list", "@1", "particules", "particule",
  "index", "commands", "command", "ifstatement", "action", "condition",
  "e", 0
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
     305,   306,    43,    45,    42,    47,   307,    94,    58,   123,
     125,   109,   113,    80,    69,    44,    95,    63,    40,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    70,    71,    72,    72,    73,    73,    74,    74,    74,
      74,    74,    74,    74,    74,    74,    74,    74,    74,    76,
      75,    77,    77,    78,    78,    78,    78,    78,    78,    78,
      78,    78,    78,    79,    79,    80,    80,    81,    81,    81,
      82,    83,    83,    83,    84,    84,    84,    84,    84,    84,
      84,    84,    84,    84,    84,    85,    85,    85,    85,    85,
      85,    85,    85,    85,    85,    85,    85,    85,    85
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     2,     2,     0,     4,     4,     4,     4,     4,
       4,     4,     4,     4,     4,     4,     4,     2,     2,     0,
       6,     2,     1,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     2,     1,     2,     0,     2,     2,     2,
       5,     1,     1,     1,     3,     3,     3,     3,     3,     3,
       4,     4,     3,     3,     3,     3,     3,     3,     3,     3,
       2,     4,     4,     4,     3,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       4,     0,    36,     1,     0,     3,     2,     0,     0,    35,
       0,    38,    65,    68,     0,     0,     0,    66,     0,     0,
       0,    67,     0,    39,    37,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    32,     0,     5,    22,     6,
       0,     0,     0,    60,    32,    19,     0,     0,    17,    18,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    21,     0,     0,     0,
       0,     0,    54,    64,    53,    52,    42,    43,     0,    41,
      44,    45,    46,    47,    48,    49,     0,     0,    55,    56,
      57,    58,    59,    34,     0,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    62,    61,    63,    11,    12,    13,
      14,    15,    16,     7,     8,     9,    10,     0,     0,     0,
      50,    51,    33,     0,    40,    56,    20
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,     2,     5,    21,    22,    81,    37,    38,   105,
       6,     9,    87,    88,    89,    25
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -42
static const yytype_int16 yypact[] =
{
     -42,     2,     6,   -42,   -25,   -42,    13,   -35,   119,   -42,
      53,   -42,   -42,   -42,   -41,   -40,   -27,   -42,   178,   217,
     202,   -42,    -5,   -42,   -38,   266,   -32,   -31,     1,     4,
      12,    16,    18,    31,    42,   -15,   178,   217,   -42,    22,
     178,   178,   178,    23,   -42,   117,   -36,   227,   -42,   -42,
     202,   202,   129,   178,   178,   178,   178,   178,   178,   178,
     178,   178,   178,   178,   178,   178,   -28,   -28,   -28,   -28,
     -28,   -28,   -28,   -28,   -28,   -33,   -42,   195,   233,   240,
      86,    45,   -42,   -42,    64,   -42,   -42,   -42,    30,   -38,
      22,    22,    22,    22,    22,    22,   153,   153,    41,    41,
      23,    23,    23,   -42,    70,   -42,   -42,   -42,   -42,   -42,
     -42,   -42,   -42,   -42,   -42,   -42,   -42,   -42,   -42,   -42,
     -42,   -42,   -42,   -42,   -42,   -42,   -42,   217,   129,   178,
      22,    22,   -42,   126,   -42,    14,   -42
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -42,   -42,   -42,   -42,   -42,   -42,   -42,   -18,   -34,   257,
     -42,   -42,   105,   -14,    65,   -10
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -69
static const yytype_int16 yytable[] =
{
      39,    45,     3,    76,    50,    51,    50,    51,    43,     4,
      47,    76,     7,   103,   -60,    48,    49,     8,   -60,    61,
      62,    63,    64,    10,    65,   104,    75,    40,    41,    52,
      77,    78,    79,    82,    66,    67,    83,   -68,   -68,   -68,
     -68,    42,   -68,    90,    91,    92,    93,    94,    95,    96,
      97,    98,    99,   100,   101,   102,   -60,   -60,    26,    27,
      28,    29,    30,    31,    32,    33,    34,    68,   -60,   -60,
      69,    65,   -60,    24,    61,    62,    63,    64,    70,    65,
      65,   -60,    71,   -60,    72,    46,   130,   131,   128,    12,
      35,    14,    15,    16,    17,    63,    64,    73,    65,    76,
     117,   118,   119,   120,   121,   122,    18,    51,    74,   133,
     127,   132,    19,    23,   134,    84,    85,     0,     0,   135,
       0,    36,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    26,    27,    28,    29,    30,    31,    32,    33,    34,
       0,     0,     0,     0,     0,     0,     0,   123,   124,   125,
     126,     0,    11,     0,    44,    12,    13,    14,    15,    16,
      17,     0,    86,    44,     0,    12,    13,    14,    15,    16,
      17,     0,    18,     0,     0,     0,     0,    80,    19,     0,
       0,     0,    18,     0,     0,     0,   136,    20,    19,    12,
      13,    14,    15,    16,    17,     0,     0,    20,     0,     0,
       0,     0,     0,     0,     0,    61,   129,    63,    64,     0,
      65,     0,    19,     0,    12,    13,    14,    15,    16,    17,
       0,    36,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    18,     0,     0,     0,     0,     0,    19,    12,    13,
      14,    15,    16,    17,     0,     0,    36,    61,    62,    63,
      64,     0,    65,     0,    44,    18,     0,     0,     0,     0,
       0,    19,     0,     0,   114,     0,     0,     0,     0,     0,
      20,    53,    54,    55,    56,    57,    58,    59,    60,    61,
      62,    63,    64,     0,    65,    61,    62,    63,    64,     0,
      65,     0,    61,    62,    63,    64,    83,    65,     0,     0,
       0,     0,   115,     0,     0,     0,     0,     0,     0,   116,
      53,    54,    55,    56,    57,    58,    59,    60,    61,    62,
      63,    64,     0,    65,   106,   107,   108,   109,   110,   111,
     112,   113
};

static const yytype_int16 yycheck[] =
{
      10,    19,     0,    37,    42,    43,    42,    43,    18,     3,
      20,    45,    37,    41,     0,    20,    21,     4,     4,    52,
      53,    54,    55,    58,    57,    53,    36,    68,    68,    67,
      40,    41,    42,    69,    66,    66,    69,    52,    53,    54,
      55,    68,    57,    53,    54,    55,    56,    57,    58,    59,
      60,    61,    62,    63,    64,    65,    42,    43,     5,     6,
       7,     8,     9,    10,    11,    12,    13,    66,    54,    55,
      66,    57,    58,     8,    52,    53,    54,    55,    66,    57,
      57,    67,    66,    69,    66,    20,    96,    97,    58,    36,
      37,    38,    39,    40,    41,    54,    55,    66,    57,   133,
      14,    15,    16,    17,    18,    19,    53,    43,    66,   127,
      65,    41,    59,     8,   128,    50,    51,    -1,    -1,   129,
      -1,    68,     5,     6,     7,     8,     9,    10,    11,    12,
      13,     5,     6,     7,     8,     9,    10,    11,    12,    13,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    61,    62,    63,
      64,    -1,    33,    -1,    37,    36,    37,    38,    39,    40,
      41,    -1,    33,    37,    -1,    36,    37,    38,    39,    40,
      41,    -1,    53,    -1,    -1,    -1,    -1,    60,    59,    -1,
      -1,    -1,    53,    -1,    -1,    -1,    60,    68,    59,    36,
      37,    38,    39,    40,    41,    -1,    -1,    68,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    52,    53,    54,    55,    -1,
      57,    -1,    59,    -1,    36,    37,    38,    39,    40,    41,
      -1,    68,     5,     6,     7,     8,     9,    10,    11,    12,
      13,    53,    -1,    -1,    -1,    -1,    -1,    59,    36,    37,
      38,    39,    40,    41,    -1,    -1,    68,    52,    53,    54,
      55,    -1,    57,    -1,    37,    53,    -1,    -1,    -1,    -1,
      -1,    59,    -1,    -1,    69,    -1,    -1,    -1,    -1,    -1,
      68,    44,    45,    46,    47,    48,    49,    50,    51,    52,
      53,    54,    55,    -1,    57,    52,    53,    54,    55,    -1,
      57,    -1,    52,    53,    54,    55,    69,    57,    -1,    -1,
      -1,    -1,    69,    -1,    -1,    -1,    -1,    -1,    -1,    69,
      44,    45,    46,    47,    48,    49,    50,    51,    52,    53,
      54,    55,    -1,    57,    67,    68,    69,    70,    71,    72,
      73,    74
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    71,    72,     0,     3,    73,    80,    37,     4,    81,
      58,    33,    36,    37,    38,    39,    40,    41,    53,    59,
      68,    74,    75,    82,    84,    85,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    37,    68,    77,    78,    85,
      68,    68,    68,    85,    37,    77,    84,    85,    20,    21,
      42,    43,    67,    44,    45,    46,    47,    48,    49,    50,
      51,    52,    53,    54,    55,    57,    66,    66,    66,    66,
      66,    66,    66,    66,    66,    85,    78,    85,    85,    85,
      60,    76,    69,    69,    84,    84,    33,    82,    83,    84,
      85,    85,    85,    85,    85,    85,    85,    85,    85,    85,
      85,    85,    85,    41,    53,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    69,    69,    69,    14,    15,    16,
      17,    18,    19,    61,    62,    63,    64,    65,    58,    53,
      85,    85,    41,    77,    83,    85,    60
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
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
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
        case 5:
#line 49 "parse.y"
    {
                                        pnum=0;
                                        map<string, string>::iterator it ;
                                        string name = (yyvsp[(2) - (4)].s);
                                        it = parts.find(name);
                        
                                        if(it != parts.end()) {
                                                cout <<name<<" : " ;
                                                yyerror("Particule already defined");
                                                YYERROR;//stops parsing if variable already defined
                                                
                                        }
                                         
                                         string phrase= (yyvsp[(4) - (4)].s);
                                         parts.insert(make_pair(name,phrase));
                                         //cout<<"\ndef "<<$2<<":"<<$4<<endl;
                                         
				;}
    break;

  case 6:
#line 67 "parse.y"
    {
                                        pnum=0;
                                        map<string, string>::iterator it ;
                                        string name = (yyvsp[(2) - (4)].s);
                                        it = vars.find(name);
                        
                                        if(it != vars.end()) {
                                                cout <<name<<" : " ;
                                                yyerror("Variable already defined");
                                                YYERROR;//stops parsing if variable already defined
                                                
                                        }
                                         
                                         string phrase= (yyvsp[(4) - (4)].s);
                                         vars.insert(make_pair(name,phrase));
                                         //cout<<"\ndef "<<$2<<":"<<$4<<endl;
                                         
				;}
    break;

  case 7:
#line 86 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }m";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 8:
#line 92 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }q";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 9:
#line 98 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }P";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 10:
#line 104 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }E";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 11:
#line 110 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }Phi";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 12:
#line 116 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }Eta";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 13:
#line 122 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }AbsEta";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 14:
#line 128 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }Pt";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 15:
#line 134 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }Pz";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 16:
#line 140 "parse.y"
    {     
                                        string s=(yyvsp[(2) - (4)].s);
                                        tmp="{ "+s+" }ndf";                        
                                        (yyval.s)=strdup(tmp.c_str());

                                ;}
    break;

  case 17:
#line 146 "parse.y"
    {     
                                        
                                        string s=(yyvsp[(1) - (2)].s);                                       
                                        s=s+"dR";                        
                                        (yyval.s)=strdup(s.c_str());
                                ;}
    break;

  case 18:
#line 152 "parse.y"
    {    
                                        
                                        string s=(yyvsp[(1) - (2)].s);                                       
                                        s=s+"dPhi";                        
                                        (yyval.s)=strdup(s.c_str());
                                ;}
    break;

  case 19:
#line 159 "parse.y"
    { pnum=0; ;}
    break;

  case 20:
#line 159 "parse.y"
    { 
                                                        string s=(yyvsp[(2) - (6)].s);
                                                        string s2=(yyvsp[(5) - (6)].s);
                                                        s="{ "+s+" , "+s2+" }";                        
                                                        (yyval.s)=strdup(s.c_str());

                                                        ;}
    break;

  case 21:
#line 167 "parse.y"
    {                                                 
                                                char s [512];
                                                strcpy(s,(yyval.s)); 
                                                strcat(s," ");
                                                strcat(s,(yyvsp[(2) - (2)].s));
                                                strcpy((yyval.s),s);                                       

                                        ;}
    break;

  case 22:
#line 175 "parse.y"
    {if (pnum==0){
                                                (yyval.s)=strdup((yyvsp[(1) - (1)].s));
                                                pnum++;                                                
                                        }
                                        else{                                                
                                                char s [512];
                                                strcpy(s,(yyval.s)); 
                                                strcat(s," ");
                                                strcat(s,(yyvsp[(1) - (1)].s));
                                                strcpy((yyval.s),s);
                                        };}
    break;

  case 23:
#line 187 "parse.y"
    {
                            //do something with name and index? Maybe put them in lists
                                                        
                            tmp="ele_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                            (yyval.s)=strdup(tmp.c_str());
                                                        
                            ;}
    break;

  case 24:
#line 194 "parse.y"
    {       tmp="muo_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 25:
#line 198 "parse.y"
    {       tmp="lep_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 26:
#line 202 "parse.y"
    {       tmp="pho_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 27:
#line 206 "parse.y"
    {       tmp="jet_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 28:
#line 210 "parse.y"
    {       tmp="bjet_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 29:
#line 214 "parse.y"
    {       tmp="qgjet_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 30:
#line 218 "parse.y"
    {       tmp="numet_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 31:
#line 222 "parse.y"
    {       tmp="metlv_"+to_string((int)(yyvsp[(3) - (3)].real));                        
                                (yyval.s)=strdup(tmp.c_str());
                                
                        ;}
    break;

  case 32:
#line 226 "parse.y"
    { //we want the original defintions as well
                map<string, string>::iterator it ;
                it = parts.find((yyvsp[(1) - (1)].s));
     
                if(it == parts.end()) {
                        cout <<(yyvsp[(1) - (1)].s)<<" : " ;
                        yyerror("Particule not defined");
                        YYERROR;//stops parsing if variable not found
                        
                }
                else {
                        tmp= it->second ;
                        (yyval.s)=strdup(tmp.c_str());
                }

               ;}
    break;

  case 33:
#line 243 "parse.y"
    {(yyval.real)=-(yyvsp[(2) - (2)].real);;}
    break;

  case 34:
#line 244 "parse.y"
    {(yyval.real)= (yyvsp[(1) - (1)].real);;}
    break;

  case 37:
#line 249 "parse.y"
    {;;}
    break;

  case 38:
#line 250 "parse.y"
    {;;}
    break;

  case 39:
#line 251 "parse.y"
    {;;}
    break;

  case 40:
#line 253 "parse.y"
    { string s1=(yyvsp[(1) - (5)].s); string s3=(yyvsp[(3) - (5)].s);string s4=(yyvsp[(5) - (5)].s);
                        tmp=s1+" [] "+s3+" "+s4;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 41:
#line 258 "parse.y"
    {(yyval.s)=(yyvsp[(1) - (1)].s);;}
    break;

  case 42:
#line 259 "parse.y"
    {tmp= " all " ;
                        (yyval.s)=strdup(tmp.c_str());;;}
    break;

  case 43:
#line 261 "parse.y"
    {(yyval.s)=(yyvsp[(1) - (1)].s);;}
    break;

  case 44:
#line 263 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" < "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 45:
#line 267 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" > "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 46:
#line 271 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" <= "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 47:
#line 275 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" >= "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 48:
#line 279 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" == "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 49:
#line 283 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" != "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 50:
#line 287 "parse.y"
    { string s1=(yyvsp[(1) - (4)].s); string s3=(yyvsp[(3) - (4)].s);string s4=(yyvsp[(4) - (4)].s);
                        tmp=s1+" [] "+s3+" "+s4;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 51:
#line 291 "parse.y"
    { string s1=(yyvsp[(1) - (4)].s); string s3=(yyvsp[(3) - (4)].s);string s4=(yyvsp[(4) - (4)].s);
                        tmp=s1+" ][ "+s3+" "+s4; 
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 52:
#line 295 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" and "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 53:
#line 299 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
                        tmp=s1+" or "+s3;   
                        (yyval.s)=strdup(tmp.c_str()); 
                        ;}
    break;

  case 54:
#line 303 "parse.y"
    { string s3=(yyvsp[(2) - (3)].s);
                                tmp=" ( "+s3+" ) ";   
                                (yyval.s)=strdup(tmp.c_str()); 
                                ;}
    break;

  case 55:
#line 308 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
               tmp=s1+" + "+s3;   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 56:
#line 312 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
               tmp=s1+" - "+s3;   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 57:
#line 316 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
               tmp=s1+" * "+s3;   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 58:
#line 320 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
               tmp=s1+" / "+s3;   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 59:
#line 324 "parse.y"
    { string s1=(yyvsp[(1) - (3)].s); string s3=(yyvsp[(3) - (3)].s);
               tmp=s1+" ^ "+s3;   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 60:
#line 328 "parse.y"
    { string s1=(yyvsp[(2) - (2)].s);
               tmp=" -"+s1;   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 61:
#line 332 "parse.y"
    { string s3=(yyvsp[(3) - (4)].s);
               tmp=" cos( "+s3+" ) ";   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 62:
#line 336 "parse.y"
    { string s3=(yyvsp[(3) - (4)].s);
               tmp=" sin( "+s3+" ) ";   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 63:
#line 340 "parse.y"
    { string s3=(yyvsp[(3) - (4)].s);
               tmp=" tan( "+s3+" ) ";   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 64:
#line 344 "parse.y"
    { string s3=(yyvsp[(2) - (3)].s);
               tmp=" ( "+s3+" ) ";   
               (yyval.s)=strdup(tmp.c_str()); 
               ;}
    break;

  case 65:
#line 348 "parse.y"
    { tmp=to_string((yyvsp[(1) - (1)].real)); (yyval.s)=strdup(tmp.c_str()); ;}
    break;

  case 66:
#line 349 "parse.y"
    { tmp=to_string((int)(yyvsp[(1) - (1)].real)); (yyval.s)=strdup(tmp.c_str()); ;}
    break;

  case 67:
#line 350 "parse.y"
    {(yyval.s)=(yyvsp[(1) - (1)].s); pnum=0;;}
    break;

  case 68:
#line 352 "parse.y"
    { //we want the original defintions as well
                map<string, string>::iterator it ;
                it = vars.find((yyvsp[(1) - (1)].s));
     
                if(it == vars.end()) {
                        cout <<(yyvsp[(1) - (1)].s)<<" : " ;
                        yyerror("Variable not defined");
                        YYERROR;//stops parsing if variable not found
                        
                }
                else {
                        tmp= it->second ;
                        (yyval.s)=strdup(tmp.c_str());
                }
               ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2145 "b.cpp"
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


#line 368 "parse.y"

int main(void) {
        yyparse(); 
        cout<<"\n Variables results: \n";
	std::map<std::string, string>::iterator it = vars.begin();
        while(it != vars.end())
        {
                std::cout<<it->first<<" :: "<<it->second<<std::endl;
                it++;
        }
        cout<<"\n Particles results: \n";
	 it = parts.begin();
        while(it != parts.end())
        {
                std::cout<<it->first<<" :: "<<it->second<<std::endl;
                it++;
        }
	cout<<"\n CUTS : \n";
	std::map<int, string>::iterator iter = cuts.begin();
        while(iter != cuts.end())
        {
                cout<<iter->first<<" :: "<<iter->second<<endl;
                iter++;
        }		
                }


                // simple calculator
// e : e '+' e  { $$ = $1 + $3 ;string s=$2;
//                                         tmp="{ "+s+" }m";                        
//                                         $$=strdup(tmp.c_str()); }
//    | e '-' e { $$ = $1 - $3 ; }
//    | e '*' e { $$ = $1 * $3 ; }
//    | e '/' e { $$ = $1 / $3*1.0 ; }
//    | e '^' e { $$ = pow($1,$3) ;  } 	
//    |'-' e %prec Unary { $$ = - $2 ; }
//    | COS '(' e ')' { $$ = cos($3) ;}
//    | SIN '(' e ')' { $$ = sin($3) ;}
//    | TAN '(' e ')' { $$ = tan($3) ;}
//    |'(' e ')' { $$ = $2 ;}
//    | NB { $$ = $1 ;} 
//    | INT { $$ = $1 ;}	
//    | function 
//    ;
