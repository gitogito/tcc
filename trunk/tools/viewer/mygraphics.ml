open Printf

type t = {
  width : int; height : int;
  xa : float; ya : float;
  xb : float; yb : float;
}

let get_int_xy gr x y =
  (int_of_float ((x -. gr.xa) /. (gr.xb -. gr.xa) *. float gr.width),
   int_of_float ((y -. gr.ya) /. (gr.yb -. gr.ya) *. float gr.height))

let open_graph width height xa ya xb yb =
  Graphics.open_graph (sprintf " %dx%d" width height);
  { width; height; xa; ya; xb; yb }

let moveto gr x y =
  let nx, ny = get_int_xy gr x y in
  Graphics.moveto nx ny

let draw_poly_line gr ary =
  Graphics.draw_poly_line
    (Array.map
       (fun (x, y) ->
          (int_of_float ((x -. gr.xa) /. (gr.xb -. gr.xa) *. float gr.width),
           int_of_float ((y -. gr.ya) /. (gr.yb -. gr.ya) *. float gr.height)))
       ary)

let fill_rect gr xa ya xb yb =
  let x = int_of_float ((xa -. gr.xa) /. (gr.xb -. gr.xa) *. float gr.width) in
  let y = int_of_float ((ya -. gr.ya) /. (gr.yb -. gr.ya) *. float gr.height) in
  let w = int_of_float ((xb -. xa) /. (gr.xb -. gr.xa) *. float gr.width) in
  let h = int_of_float ((yb -. ya) /. (gr.yb -. gr.ya) *. float gr.height) in
  Graphics.fill_rect x y w h
