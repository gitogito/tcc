open Syntax

let rec print_var_assigns = function
    `Nil -> ()
  | `Var_assigns (va, vas) ->
      (* var_assigns has vars in reversed order, so va is printed after vas *)
      print_var_assigns vas;
      match va with `Var_assign (s, e) ->
        Printf.printf "%s = %s\n" s (string_of_expr e)

let print_world = function
    `World (pnt1, pnt2, dx, dy, dz) ->
      print_endline "--- world";
      print_string (string_of_point pnt1);
      print_string ", ";
      print_string (string_of_point pnt2);
      Printf.printf ", %s, %s, %s\n"
        (string_of_expr dx)
        (string_of_expr dy)
        (string_of_expr dz)

let rec print_objs = function
    `Nil -> ()
  | `Objs (obj, objs) ->
      print_objs objs;
      print_obj obj

and print_obj = function
    `Box (p1, p2) ->
      Printf.printf "box %s, %s\n" (string_of_point p1) (string_of_point p2)
  | `Rect (s, p, p2d) ->
      Printf.printf "rect %s, %s, %s\n" s (string_of_point p)
        (string_of_point2d p2d)
  | `Triangle (s, p, p2d1, p2d2) ->
      Printf.printf "triangle %s, %s, %s, %s\n" s (string_of_point p)
        (string_of_point2d p2d1) (string_of_point2d p2d2)
  | `Circle (s, p, e) ->
      Printf.printf "circle %s, %s, %s\n" s (string_of_point p)
        (string_of_expr e)
  | `Circleperi (s, p, e1, e2, e3) ->
      Printf.printf "circleperi %s, %s, %s, %s, %s\n" s (string_of_point p)
        (string_of_expr e1) (string_of_expr e2) (string_of_expr e3)
  | `Ellipse (s, p, e1, e2) ->
      Printf.printf "ellipse %s, %s, %s, %s\n" s (string_of_point p)
        (string_of_expr e1) (string_of_expr e2)
  | `Ellipseperi (s, p, e1, e2, e3, e4) ->
      Printf.printf "ellipseperi %s, %s, %s, %s, %s, %s\n" s (string_of_point p)
        (string_of_expr e1) (string_of_expr e2)
        (string_of_expr e3) (string_of_expr e4)
  | `Polygon (s, p, p2d_list) ->
      Printf.printf "polygon %s, %s, %s\n" s (string_of_point p)
        (string_of_point2d_list p2d_list)
  | `Line (s, p, p2d_list) ->
      Printf.printf "line %s, %s, %s\n" s (string_of_point p)
        (string_of_point2d_list p2d_list)
  | `Sweep (s, e, obj) ->
      Printf.printf "sweep %s, %s, " s (string_of_expr e);
      print_obj obj
  | `Objs objs ->
      print_endline "["
      print_objs objs
      print_endline "]"

let print_command = function
    `Active obj ->
      print_endline "--- active";
      print_obj obj
  | `Inactive obj ->
      print_endline "--- inactive";
      print_obj obj
  | `Fix (expr, obj) ->
      Printf.printf "--- fix %s\n" (string_of_expr expr);
      print_obj obj
  | `Fixheat (expr, obj) ->
      Printf.printf "--- fix heat %s\n" (string_of_expr expr);
      print_obj obj
  | `Heat (expr, obj) ->
      Printf.printf "--- heat %s\n" (string_of_expr expr);
      print_obj obj
  | `Lambda (expr, obj) ->
      Printf.printf "--- lambda %s\n" (string_of_expr expr);
      print_obj obj

let rec print_commands = function
    `Nil -> ()
  | `Commands (a, l) ->
      print_command a;
      print_newline ();
      print_commands l

let print_input = function
    `Input (var_assigns, world, commands) ->
      print_var_assigns var_assigns;
      print_newline ();

      print_world world;
      print_newline ();

      print_commands commands

let () =
  let lexbuf = Lexing.from_channel stdin in
  let input = Parser.input Lexer.token lexbuf in
    print_input input
