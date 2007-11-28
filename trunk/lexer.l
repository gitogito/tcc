%{
#include "sim.h"
#include "y.tab.h"
#include "mem.h"

#define yywrap()	1

extern int yydebug;

long lineno = 1;
%}

%option nounput

ws	[ \t]+
word	[a-zA-Z_][a-zA-Z_0-9]*
digit	[0-9]
num1	{digit}+\.?([eE][-+]?{digit}+)?
num2	{digit}*\.{digit}+([eE][-+]?{digit}+)?
number	{num1}|{num2}

%%

^{ws}?(\#.*)?\n {
    /* ignore blank line (including comment) */
    ++lineno;
}

{ws} {
    /* ignore white spaces */
}

\#.* {
    /* ignore comment */
}

---+		return TK_LINE;
active		return TK_ACTIVE;
box		return TK_BOX;
circle		return TK_CIRCLE;
edge		return TK_EDGE;
ellipse		return TK_ELLIPSE;
fix		return TK_FIX;
heat		return TK_HEAT;
lambda		return TK_LAMBDA;
noactive	return TK_NOACTIVE;
polygon		return TK_POLYGON;
rect		return TK_RECT;
triangle	return TK_TRIANGLE;
sweep		return TK_SWEEP;
world		return TK_WORLD;

{number} {
    yylval.val = atof(yytext);
    return TK_NUMBER;
}

{word} {
    yylval.str = estrdup(yytext);
    return TK_WORD;
}

:{word} {
    yylval.str = estrdup(yytext + 1);
    return TK_SYMBOL;
}

\*\* {
    return TK_POW;
}

. {
    return yytext[0];
}

\n {
    /* ignore */
    ++lineno;
}

%%
