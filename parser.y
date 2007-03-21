%{

#include <string.h>
#include "tc.h"
#include "sim.h"
#include "parser.h"

#define YYDEBUG 1

int yyerror(const char *s);
int yylex();

extern long lineno;

enum {
    ST_START,
    ST_WORLD,
    ST_ACTIVE,
    ST_FIX,
    ST_HEATFLOW,
    ST_LAMBDA,
};

static int state = ST_START;
static double value;
static int dir = -1;

%}

%union {
    double val;
    char *str;
    struct point_t point;
    struct vector2d_t vector2d;
    Obj *obj;
}

%token <val> TK_NUMBER
%token <str> TK_SYMBOL

%token	TK_LINES TK_ACTIVE TK_BOX TK_FIX TK_HEATFLOW TK_LAMBDA TK_RECT TK_WORLD

%type <point>		point
%type <point>		vector
%type <vector2d>	vector2d
%type <val>		expr
%type <obj>		obj

%left		'+' '-'
%left		'*' '/'
%right		NEG
%right		TK_POW

%%

input: 
    /* empty */

  | input command

command:
    TK_LINES TK_WORLD point ',' vector ',' expr ',' expr ',' expr
    {
	if (config_parser->world != NULL)
	    yyerror("world is already defined");
        config_parser->world = world_new($3.x, $3.y, $3.z, $5.x, $5.y, $5.z, $7, $9, $11);
	state = ST_WORLD;
    }

  | TK_LINES TK_ACTIVE
    {
	state = ST_ACTIVE;
    }

  | TK_LINES TK_FIX expr
    {
	state = ST_FIX;
	value = $3;
    }

  | TK_LINES TK_HEATFLOW TK_SYMBOL ',' expr
    {
	state = ST_HEATFLOW;

	if (strcmp($3, ":LEFT") == 0)
	    dir = DIR_LEFT;
	else if (strcmp($3, ":RIGHT") == 0)
	    dir = DIR_RIGHT;
	else if (strcmp($3, ":FRONT") == 0)
	    dir = DIR_FRONT;
	else if (strcmp($3, ":BACK") == 0)
	    dir = DIR_BACK;
	else if (strcmp($3, ":BELOW") == 0)
	    dir = DIR_BELOW;
	else if (strcmp($3, ":ABOVE") == 0)
	    dir = DIR_ABOVE;
	else
	    yyerror("unknown direction");

	value = $5;
    }

  | TK_LINES TK_LAMBDA expr
    {
	state = ST_LAMBDA;
	value = $3;
    }

  | obj
    {
	Obj *obj;

	obj = $1;
	switch (state) {
	case ST_ACTIVE:
	    obj->uval.i = 1;
	    aryobj_push(config_parser->active_obj_ary, obj);
	    break;
	case ST_FIX:
	    obj->uval.d = value;
	    aryobj_push(config_parser->fix_obj_ary, obj);
	    break;
	case ST_HEATFLOW:
	    obj->uval.h = heatflow_new(dir, value);
	    aryobj_push(config_parser->heatflow_obj_ary, obj);
	    break;
	case ST_LAMBDA:
	    obj->uval.d = value;
	    aryobj_push(config_parser->lambda_obj_ary, obj);
	    break;
	default:
	    yyerror("object must not be here");
	}
    }

point:
    '(' expr ',' expr ',' expr ')'
    {
	$$.x = $2;
	$$.y = $4;
	$$.z = $6;
    }

vector:
    '<' expr ',' expr ',' expr '>'
    {
	$$.x = $2;
	$$.y = $4;
	$$.z = $6;
    }

vector2d:
    '<' expr ',' expr '>'
    {
	$$.v1 = $2;
	$$.v2 = $4;
    }

expr:
    expr '+' expr { $$ = $1 + $3; }

  | expr '-' expr { $$ = $1 - $3; }

  | expr '*' expr { $$ = $1 * $3; }

  | expr '/' expr { $$ = $1 / $3; }

  | '(' expr ')' { $$ = $2; }

  | '-' TK_NUMBER  %prec NEG { $$ = - $2; }

  | TK_NUMBER { $$ = $1; }

obj:
    TK_BOX point ',' vector
    {
	$$ = obj_new(OBJ_BOX);
	$$->uobj.box = box_new(config_parser->world, $2.x, $2.y, $2.z, $4.x, $4.y, $4.z);
    }

  | TK_RECT point ',' TK_SYMBOL ',' vector2d
    {
	int axis;

	$$ = obj_new(OBJ_RECT);
	if (strcmp($4, ":X") == 0)
	    axis = AXIS_X;
	else if (strcmp($4, ":Y") == 0)
	    axis = AXIS_Y;
	else if (strcmp($4, ":Z") == 0)
	    axis = AXIS_Z;
	else
	    yyerror("unknown axis");
	$$->uobj.rect = rect_new(config_parser->world, $2.x, $2.y, $2.z, axis, $6.v1, $6.v2);
    }

%%

int yyerror(const char *s)
{
    warn_exit("%s at line %ld", s, lineno);
}
