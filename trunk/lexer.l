%{
#include "sim.h"
#include "y.tab.h"
#include "mem.h"

extern int yydebug;

long lineno = 1;
%}

%option nounput
%option noyywrap

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

---+		return TK_DASH;
active		return TK_ACTIVE;
box		return TK_BOX;
circleperi	return TK_CIRCLEPERI;
circle		return TK_CIRCLE;
ellipseperi	return TK_ELLIPSEPERI;
ellipse		return TK_ELLIPSE;
fix		return TK_FIX;
heat		return TK_HEAT;
lambda		return TK_LAMBDA;
line		return TK_LINE;
inactive	return TK_INACTIVE;
polygon		return TK_POLYGON;
rect		return TK_RECT;
triangle	return TK_TRIANGLE;
sphere		return TK_SPHERE;
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
