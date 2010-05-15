{
open Parser

let lineno = ref 1
}

let digit = ['0'-'9']
let num1 = digit+ '.'? (['e' 'E'] ['-' '+']? digit+)?
let num2 = digit* '.' digit+ (['e' 'E'] ['-' '+']? digit+)?
let number = num1 | num2

let alpha = ['a'-'z' 'A'-'Z']
let word = (alpha | '_') (alpha | '_' | digit)*

rule token = parse
  [' ' '\t']            { token lexbuf }
| '#' [^'\n']* '\n'     { lineno := !lineno + 1; token lexbuf }
| '\n'                  { lineno := !lineno + 1; token lexbuf }
| "active"              { ACTIVE }
| "box"                 { BOX }
| "circle"              { CIRCLE }
| "circleperi"          { CIRCLEPERI }
| "ellipse"             { ELLIPSE }
| "ellipseperi"         { ELLIPSEPERI }
| "fix"                 { FIX }
| "heat"                { HEAT }
| "lambda"              { LAMBDA }
| "line"                { LINE }
| "inactive"            { INACTIVE }
| "pi"                  { PI }
| "polygon"             { POLYGON }
| "rect"                { RECT }
| "sweep"               { SWEEP }
| "triangle"            { TRIANGLE }
| "world"               { WORLD }
| number as s           { NUMBER (float_of_string s) }
| word as s             { WORD s }
| ':' word as s         { SYMBOL s }
| "---" '-'*            { DASH }
| '('                   { LPAR }
| ')'                   { RPAR }
| '<'                   { LPARA }
| '>'                   { RPARA }
| '['                   { LPARL }
| ']'                   { RPARL }
| ','                   { COMMA }
| "**"                  { POW }
| '+'                   { PLUS }
| '-'                   { MINUS }
| '*'                   { MULTI }
| '/'                   { DIV }
| '='                   { EQ }
| eof                   { EOF }
| _
    {
      failwith
        (Printf.sprintf
           "unknown token '%s' at line %d" (Lexing.lexeme lexbuf) !lineno)
    }
