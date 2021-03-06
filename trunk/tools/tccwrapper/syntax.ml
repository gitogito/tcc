(* All lists are LIFO (last in, first out). Take care of order *)

type 'a mylist = Nil | Cons of 'a * 'a mylist

type var_assign = Var_assign of string * expr
and expr =
  Plus of expr * expr
| Minus of expr * expr
| Multi of expr * expr
| Div of expr * expr
| Pow of expr * expr
| Paren of expr
| Neg of expr
| Number of float
| Var of string
| Var_assign_expr of var_assign

type point = Point of expr * expr * expr
type vector = Vector of expr * expr * expr
type point2d = Point2d of expr * expr
type vector2d = Vector2d of expr * expr

type obj =
    Box of point * point
  | Rect of string * point * point2d
  | Triangle of string * point * point2d * point2d
  | Circle of string * point * expr
  | Circleperi of string * point * expr * expr * expr
  | Ellipse of string * point * expr * expr
  | Ellipseperi of string * point * expr * expr * expr * expr
  | Polygon of string * point * point2d mylist
  | Line of string * point * point2d mylist
  | Sweep of string * expr * obj
  | Objs of obj mylist

type world = World of point * point * expr * expr * expr
type command =
    Active of obj
  | Inactive of obj
  | Fix of expr * obj
  | Fixheat of expr * obj
  | Heat of expr * obj
  | Lambda of expr * obj

type input = Input of var_assign mylist * world * command mylist

(* return point from vector relative to the point *)
let to_point point vector =
  match point, vector with
    Point (px, py, pz), Vector (vx, vy, vz) ->
      Point (Plus (px, vx), Plus (py, vy), Plus (pz, vz))

(* return point2d from vector2d relative to the point *)
let to_point2d axis point vector2d =
  match point, vector2d with
    Point (px, py, pz), Vector2d (v1, v2) ->
      match axis with
        ":X" -> Point2d (Plus (py, v1), Plus (pz, v2))
      | ":Y" -> Point2d (Plus (pz, v1), Plus (px, v2))
      | ":Z" -> Point2d (Plus (px, v1), Plus (py, v2))
      | _ -> failwith "unknown axis"

(* return point2d list from vector2d_list relative to the point *)
let rec to_point2d_list axis point vector2d_list =
  match vector2d_list with
    Nil -> Nil
  | Cons (v, vl) ->
      Cons (to_point2d axis point v, to_point2d_list axis point vl)

let rec string_of_expr = function
    Plus (e1, e2) -> Printf.sprintf "%s + %s" (string_of_expr e1)
                        (string_of_expr e2)
  | Minus (e1, e2) -> Printf.sprintf "%s - %s" (string_of_expr e1)
                         (string_of_expr e2)
  | Multi (e1, e2) -> Printf.sprintf "%s * %s" (string_of_expr e1)
                         (string_of_expr e2)
  | Div (e1, e2) -> Printf.sprintf "%s / %s" (string_of_expr e1)
                       (string_of_expr e2)
  | Pow (e1, e2) -> Printf.sprintf "%s ** %s" (string_of_expr e1)
                       (string_of_expr e2)
  | Paren e -> Printf.sprintf "(%s)" (string_of_expr e)
  | Neg e -> Printf.sprintf "-%s" (string_of_expr e)
  | Number e -> Printf.sprintf "%g" e
  | Var e -> e
  | Var_assign_expr Var_assign (s, e) ->
      Printf.sprintf "%s = %s\n" s (string_of_expr e)

let string_of_point = function
    Point (x, y, z) ->
      Printf.sprintf "(%s, %s, %s)"
        (string_of_expr x)
        (string_of_expr y)
        (string_of_expr z)

let string_of_point2d = function
    Point2d (p1, p2) ->
      Printf.sprintf "(%s, %s)"
        (string_of_expr p1)
        (string_of_expr p2)

(* return "," in the beginning of string *)
let rec string_of_point2d_list = function
    Nil -> ""
  | Cons (p2d, p2dl) ->
      string_of_point2d_list p2dl ^ ", " ^ string_of_point2d p2d
