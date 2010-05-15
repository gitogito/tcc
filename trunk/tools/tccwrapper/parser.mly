%{
open Syntax

let pi = 4.0 *. atan(1.0)
%}

%token EOF
%token ACTIVE BOX CIRCLE CIRCLEPERI ELLIPSE ELLIPSEPERI FIX
       HEAT INACTIVE LAMBDA LINE PI POLYGON RECT SWEEP
       TRIANGLE WORLD
%token DASH LPAR RPAR LPARA RPARA LPARL RPARL COMMA
%token POW PLUS MINUS MULTI DIV EQ
%token <float> NUMBER
%token <string> WORD
%token <string> SYMBOL

%left  PLUS MINUS
%left  MULTI DIV
%right NEG
%right POW

%start input

%type <Syntax.input> input

%%

input:
    world commands {
      `Input (`Nil, $1, $2) }
  | var_assigns world commands {
      `Input ($1, $2, $3) }
;

var_assigns:
    var_assign
    {
      match $1 with
        `Var_assign (s, e) -> `Var_assigns (`Var_assign (s, e), `Nil)
      | _ -> failwith "not reached"
    }
  | var_assigns var_assign
    {
      match $2 with
        `Var_assign (s, e) -> `Var_assigns (`Var_assign (s, e), $1)
      | _ -> failwith "not reached"
    }
;

var_assign:
    WORD EQ expr { `Var_assign ($1, $3) }
;

world:
    DASH WORLD point COMMA vector COMMA expr COMMA expr COMMA expr
    {
      `World ($3, to_point $3 $5, $7, $9, $11)
    }
  | DASH WORLD point COMMA point COMMA expr COMMA expr COMMA expr
    {
      `World ($3, $5, $7, $9, $11)
    }
;

commands:
    command { `Commands ($1, `Nil) }
  | commands command { `Commands ($2, $1) }
;

command:
    DASH ACTIVE obj { `Active $3 }
  | DASH INACTIVE obj { `Inactive $3 }
  | DASH FIX expr obj { `Fix ($3, $4) }
  | DASH FIX HEAT expr obj { `Fixheat ($4, $5) }
  | DASH HEAT FIX expr obj { `Fixheat ($4, $5) }
  | DASH HEAT expr obj { `Heat ($3, $4) }
  | DASH LAMBDA expr obj { `Lambda ($3, $4) }
;

point:
    LPAR expr COMMA expr COMMA expr RPAR
    {
      `Point ($2, $4, $6)
    }
;

point2d:
    LPAR expr COMMA expr RPAR
    {
      `Point2d ($2, $4)
    }
;

point2d_list:
    point2d { `Point2d_list ($1, `Nil) }
  | point2d_list COMMA point2d { `Point2d_list ($3, $1) }
;

vector:
    LPARA expr COMMA expr COMMA expr RPARA
    {
      `Vector ($2, $4, $6)
    }
;

vector2d:
    LPARA expr COMMA expr RPARA
    {
      `Vector2d ($2, $4)
    }
;

vector2d_list:
    vector2d { `Vector2d_list ($1, `Nil) } 
  | vector2d_list COMMA vector2d { `Vector2d_list ($3, $1) }
;

expr:
    expr PLUS expr { `Plus ($1, $3) }
  | expr MINUS expr { `Minus ($1, $3) }
  | expr MULTI expr { `Multi ($1, $3) }
  | expr DIV expr { `Div ($1, $3) }
  | expr POW expr { `Pow ($1, $3) }
  | LPAR expr RPAR { `Paren $2 }
  | MINUS expr %prec NEG { `Neg $2 }
  | NUMBER { `Number $1 }
  | PI { `Number pi }
  | WORD { `Var $1 }
  | var_assign { $1 }
;

obj:
    BOX point COMMA vector
    {
      `Box ($2, to_point $2 $4)
    }
  | BOX point COMMA point
    {
      `Box ($2, $4)
    }
  | RECT SYMBOL COMMA point COMMA vector2d
    {
      `Rect ($2, $4, to_point2d $2 $4 $6)
    }
  | RECT SYMBOL COMMA point COMMA point2d
    {
      `Rect ($2, $4, $6)
    }
  | TRIANGLE SYMBOL COMMA point COMMA vector2d COMMA vector2d
    {
      `Triangle ($2, $4, to_point2d $2 $4 $6, to_point2d $2 $4 $8)
    }
  | TRIANGLE SYMBOL COMMA point COMMA point2d COMMA point2d
    {
      `Triangle ($2, $4, $6, $8)
    }
  | CIRCLE SYMBOL COMMA point COMMA expr
    {
      `Circle ($2, $4, $6)
    }
  | CIRCLEPERI SYMBOL COMMA point COMMA expr COMMA expr COMMA expr
    {
      `Circleperi ($2, $4, $6, $8, $10)
    }
  | ELLIPSE SYMBOL COMMA point COMMA expr COMMA expr
    {
      `Ellipse ($2, $4, $6, $8)
    }
  | ELLIPSEPERI SYMBOL COMMA point COMMA expr COMMA expr COMMA expr COMMA expr
    {
      `Ellipseperi ($2, $4, $6, $8, $10, $12)
    }
  | POLYGON SYMBOL COMMA point COMMA vector2d_list
    {
      `Polygon ($2, $4, to_point2d_list $2 $4 $6)
    }
  | POLYGON SYMBOL COMMA point COMMA point2d_list
    {
      `Polygon ($2, $4, $6)
    }
  | LINE SYMBOL COMMA point COMMA vector2d_list
    {
      `Line ($2, $4, to_point2d_list $2 $4 $6)
    }
  | LINE SYMBOL COMMA point COMMA point2d_list
    {
      `Line ($2, $4, $6)
    }
  | SWEEP SYMBOL COMMA expr COMMA obj
    {
      `Sweep ($2, $4, $6)
    }
  | LPARL objs RPARL
    {
      `Objs $2
    }
;

objs:
    obj { `Objs ($1, `Nil) }
  | objs obj { `Objs ($2, $1) }
;
