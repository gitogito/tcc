%{

#include <string.h>
#include <math.h>
#include "tc.h"
#include "sim.h"
#include "mem.h"

#define YYDEBUG 1

int yyerror(const char *s);
int yylex();

extern long lineno;

typedef struct var_t {
    char *name;
    double val;
    struct var_t *next;
} var_t;

static var_t *varlist = NULL;

static var_t *get_var(char *varname, int create_p)
{
    var_t *v;

    for (v = varlist; v != NULL; v = v->next) {
	if (strcmp(v->name, varname) == 0) {
	    return v;
	}
    }
    if (!create_p) {
	return NULL;
    }
    v = varlist;
    varlist = EALLOC(var_t);
    varlist->name = varname;
    varlist->next = v;
    return varlist;
}

static void var_assign(char *varname, double value)
{
    var_t *v;

    v = get_var(varname, 1);
    v->val = value;
}

static int symbol_to_axis(char *symbol)
{
    if (strcmp(symbol, "X") == 0)
	return AXIS_X;
    else if (strcmp(symbol, "Y") == 0)
	return AXIS_Y;
    else if (strcmp(symbol, "Z") == 0)
	return AXIS_Z;
    else
	yyerror("unknown axis");
    /* NOTREACHED */
    bug("NOTREACHED");
    return -1;
}

%}

%union {
    double		val;
    char		*str;
    Point		point;
    Vector2d		point2d;
    Vector2d_ary	*point2d_ary;
    Vector2d		vector2d;
    Vector2d_ary	*vector2d_ary;
    Obj			*obj;
    AryObj		*aryobj;
}

%token <val> TK_NUMBER
%token <str> TK_WORD TK_SYMBOL

%token	TK_ACTIVE TK_BOX TK_CIRCLE TK_DASH TK_ELLIPSE TK_EDGE TK_FIX TK_HEAT
        TK_LAMBDA TK_LINE TK_NOACTIVE TK_POLYGON TK_RECT TK_SWEEP TK_TRIANGLE
	TK_WORLD

%type <point>		point
%type <point>		vector
%type <point2d>		point2d
%type <point2d_ary>	point2d_ary
%type <vector2d>	vector2d
%type <vector2d_ary>	vector2d_ary
%type <val>		expr var_assign
%type <obj>		obj
%type <aryobj>		objs

%left		'+' '-'
%left		'*' '/'
%right		NEG
%right		TK_POW

%%

input:
    world commands

  | var_assigns world commands

var_assigns:
    var_assign
    {}

  | var_assigns var_assign

var_assign:
    TK_WORD '=' expr
    {
	var_assign($1, $3);
	$$ = $3;
    }

world:
    TK_DASH TK_WORLD point ',' vector ',' expr ',' expr ',' expr
    {
	if (world != NULL)
	    warn_exit("world is already defined at line %ld", lineno);
	world = world_new($3.x, $3.y, $3.z, $5.x, $5.y, $5.z, $7, $9, $11);
    }

  | TK_DASH TK_WORLD point ',' point ',' expr ',' expr ',' expr
    {
	if (world != NULL)
	    warn_exit("world is already defined at line %ld", lineno);
	world = world_new($3.x, $3.y, $3.z,
	    $5.x - $3.x, $5.y - $3.y, $5.z - $3.z, $7, $9, $11);
    }

commands:
    command

  | commands command

command:
    TK_DASH TK_ACTIVE obj
    {
	$3->uval.i = 1;
	aryobj_push(config->active_obj_ary, $3);
    }

  | TK_DASH TK_NOACTIVE obj
    {
	$3->uval.i = 0;
	aryobj_push(config->active_obj_ary, $3);
    }

  | TK_DASH TK_FIX expr obj
    {
	$4->uval.d = $3;
	aryobj_push(config->fix_obj_ary, $4);
    }

  | TK_DASH TK_FIX TK_HEAT expr obj
    {
	$5->uval.d = $4;
	aryobj_push(config->fixheat_obj_ary, $5);
    }

  | TK_DASH TK_HEAT TK_FIX expr obj
    {
	$5->uval.d = $4;
	aryobj_push(config->fixheat_obj_ary, $5);
    }

  | TK_DASH TK_HEAT expr obj
    {
	$4->uval.d = $3;
	aryobj_push(config->heat_obj_ary, $4);
    }

  | TK_DASH TK_LAMBDA expr obj
    {
	$4->uval.d = $3;
	obj_offset($4);
	aryobj_push(config->lambda_obj_ary, $4);
    }

point:
    '(' expr ',' expr ',' expr ')'
    {
	$$.x = $2;
	$$.y = $4;
	$$.z = $6;
    }

point2d:
    '(' expr ',' expr ')'
    {
	$$.x = $2;
	$$.y = $4;
    }

point2d_ary:
    point2d
    {
	$$ = vector2d_ary_new();
	vector2d_ary_push($$, $1);
    }

  | point2d_ary ',' point2d
    {
	$$ = $1;
	vector2d_ary_push($$, $3);
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
	$$.x = $2;
	$$.y = $4;
    }

vector2d_ary:
    vector2d
    {
	$$ = vector2d_ary_new();
	vector2d_ary_push($$, $1);
    }

  | vector2d_ary ',' vector2d
    {
	$$ = $1;
	vector2d_ary_push($$, $3);
    }

expr:
    expr '+' expr { $$ = $1 + $3; }

  | expr '-' expr { $$ = $1 - $3; }

  | expr '*' expr { $$ = $1 * $3; }

  | expr '/' expr { $$ = $1 / $3; }

  | expr TK_POW expr { $$ = pow($1, $3); }

  | '(' expr ')' { $$ = $2; }

  | '-' expr  %prec NEG { $$ = - $2; }

  | TK_NUMBER { $$ = $1; }

  | TK_WORD
    {
	var_t *v;

	v = get_var($1, 0);
	if (v == NULL)
	    warn_exit("can't find variable: '%s' at line %ld", $1, lineno);
	$$ = v->val;
    }

  | var_assign { $$ = $1; }

obj:
    TK_BOX point ',' vector
    {
	$$ = obj_new(OBJ_BOX);
	$$->uobj.box = box_new($2.x, $2.y, $2.z, $4.x, $4.y, $4.z);
    }

  | TK_BOX point ',' point
    {
	$$ = obj_new(OBJ_BOX);
	$$->uobj.box = box_new($2.x, $2.y, $2.z, $4.x - $2.x, $4.y - $2.y, $4.z - $2.z);
    }

  | TK_RECT TK_SYMBOL ',' point ',' vector2d
    {
	int axis;

	$$ = obj_new(OBJ_RECT);
	axis = symbol_to_axis($2);
	$$->uobj.rect = rect_new($4.x, $4.y, $4.z, axis, $6.x, $6.y);
    }

  | TK_RECT TK_SYMBOL ',' point ',' point2d
    {
	int axis;

	$$ = obj_new(OBJ_RECT);
	axis = symbol_to_axis($2);
	$$->uobj.rect = rect_new($4.x, $4.y, $4.z, axis, $6.x - $4.x, $6.y - $4.y);
    }

  | TK_TRIANGLE TK_SYMBOL ',' point ',' vector2d ',' vector2d
    {
	int axis;

	$$ = obj_new(OBJ_TRIANGLE);
	axis = symbol_to_axis($2);
	$$->uobj.triangle = triangle_new($4.x, $4.y, $4.z,
	    axis, $6.x, $6.y, $8.x, $8.y);
    }

  | TK_TRIANGLE TK_SYMBOL ',' point ',' point2d ',' point2d
    {
	int axis;

	$$ = obj_new(OBJ_TRIANGLE);
	axis = symbol_to_axis($2);
	$$->uobj.triangle = triangle_new($4.x, $4.y, $4.z,
	    axis, $6.x - $4.x, $6.y - $4.y, $8.x - $4.x, $8.y - $4.y);
    }

  | TK_CIRCLE TK_SYMBOL ',' point ',' expr
    {
	int axis;

	$$ = obj_new(OBJ_CIRCLE);
	axis = symbol_to_axis($2);
	$$->uobj.circle = circle_new($4.x, $4.y, $4.z, axis, $6);
    }

  | TK_ELLIPSE TK_SYMBOL ',' point ',' expr ',' expr
    {
	int axis;

	$$ = obj_new(OBJ_ELLIPSE);
	axis = symbol_to_axis($2);
	$$->uobj.ellipse = ellipse_new($4.x, $4.y, $4.z, axis, $6, $8);
    }

  | TK_POLYGON TK_SYMBOL ',' point ',' vector2d_ary
    {
	int axis;

	$$ = obj_new(OBJ_POLYGON);
	axis = symbol_to_axis($2);
	$$->uobj.polygon = polygon_new($4.x, $4.y, $4.z, axis, $6);
    }

  | TK_POLYGON TK_SYMBOL ',' point ',' point2d_ary
    {
	int axis;

	$$ = obj_new(OBJ_POLYGON);
	axis = symbol_to_axis($2);
	$$->uobj.polygon = polygon_new2($4.x, $4.y, $4.z, axis, $6);
    }

  | TK_LINE TK_SYMBOL ',' point ',' vector2d_ary
    {
	int axis;

	$$ = obj_new(OBJ_LINE);
	axis = symbol_to_axis($2);
	$$->uobj.line = line_new($4.x, $4.y, $4.z, axis, $6);
    }

  | TK_LINE TK_SYMBOL ',' point ',' point2d_ary
    {
	int axis;

	$$ = obj_new(OBJ_LINE);
	axis = symbol_to_axis($2);
	$$->uobj.line = line_new2($4.x, $4.y, $4.z, axis, $6);
    }

  | TK_SWEEP TK_SYMBOL ',' expr ',' obj
    {
	int axis;

	$$ = obj_new(OBJ_SWEEP);
	axis = symbol_to_axis($2);
	$$->uobj.sweep = sweep_new(axis, $4, $6);
    }

  | TK_EDGE obj
    {
	$$ = obj_new(OBJ_EDGE);
	$$->uobj.edge = edge_new($2);
    }
  | '[' objs ']'
    {
	$$ = obj_new(OBJ_OBJARY);
	$$->uobj.objary = objary_new($2);
    }

objs:
    obj
    {
	$$ = aryobj_new();
	aryobj_push($$, $1);
    }
  | objs obj
    {
	$$ = $1;
	aryobj_push($$, $2);
    }

%%

int yyerror(const char *s)
{
    warn_exit("%s at line %ld", s, lineno);
    return 0;	/* NOTREACHED */
}
