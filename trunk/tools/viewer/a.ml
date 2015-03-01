open Printf
open Scanf
open ExtLib

type world = {
  nx : int; xa : float; xb : float;
  ny : int; ya : float; yb : float;
  nz : int; za : float; zb : float;
}

type data = { v : float; act : bool }

type axis = X | Y | Z

let get_world ic =
  let aux () =
    let s = input_line ic in
    sscanf s "# %d %g %g" (fun n a b -> (n, a, b))
  in
  let nx, xa, xb = aux () in
  let ny, ya, yb = aux () in
  let nz, za, zb = aux () in
  { nx; xa; xb;
    ny; ya; yb;
    nz; za; zb }

let get_line_data ic =
  let s = input_line ic in
  sscanf s "%g %g %g %g %d" (fun x y z v act -> (x, y, z, v, act))

let get_data ic =
  let world = get_world ic in
  let ary =
    Array.init world.nx (fun _ ->
        Array.init world.ny (fun _ ->
            Array.init world.nz (fun _ -> { v = 0.0; act = false })))
  in
  for iz = 0 to world.nz - 1 do
    for iy = 0 to world.ny - 1 do
      for ix = 0 to world.nx - 1 do
        let _, _, _, v, act = get_line_data ic in
        let act = (act = 1) in
        ary.(ix).(iy).(iz) <- { v; act }
      done;
      let s = input_line ic in
      if s <> "" then
        failwith "invalid input"
    done
  done;
  (world, ary)

let get_min_max world ary =
  let mi = ref ary.(0).(0).(0).v in
  let ma = ref ary.(0).(0).(0).v in
  for iz = 0 to world.nz - 1 do
    for iy = 0 to world.ny - 1 do
      for ix = 0 to world.nx - 1 do
        if ary.(ix).(iy).(iz).v < !mi then mi := ary.(ix).(iy).(iz).v;
        if ary.(ix).(iy).(iz).v > !ma then ma := ary.(ix).(iy).(iz).v
      done
    done
  done;
  (!mi, !ma)

let get_xyz i n a b =
  a +. (float i) /. float (n - 1) *. (b -. a)

let get_2d_array world ary axis n =
  let n1, a1, b1, n2, a2, b2 =
    match axis with
    | X -> (world.ny, world.ya, world.yb, world.nz, world.za, world.zb)
    | Y -> (world.nx, world.xa, world.xb, world.nz, world.za, world.zb)
    | Z -> (world.nx, world.xa, world.xb, world.ny, world.ya, world.yb)
  in
  let ary_2d =
    Array.init n1 (fun _ ->
        Array.init n2 (fun _ -> 0.0))
  in
  for i2 = 0 to n2 - 1 do
    for i1 = 0 to n1 - 1 do
      let v =
        match axis with
        | X -> ary.(n).(i1).(i2).v
        | Y -> ary.(i1).(n).(i2).v
        | Z -> ary.(i1).(i2).(n).v
      in
      ary_2d.(i1).(i2) <- v
    done
  done;
  (n1, a1, b1, n2, a2, b2, ary_2d)

let draw gr world v_min v_max ary axis n =
  let n1, _, _, n2, _, _, ary_2d = get_2d_array world ary axis n in
  Graphics.clear_graph ();

  for i2 = 0 to n2 - 1 do
    for i1 = 0 to n1 - 1 do
      let xa = float i1 /. float n1 in
      let xb = float (i1 + 1) /. float n1 in
      let ya = float i2 /. float n2 in
      let yb = float (i2 + 1) /. float n2 in
      let gray = int_of_float (256.0 *. (ary_2d.(i1).(i2) -. v_min) /. (v_max -. v_min)) in
      let gray = if gray <= 255 then gray else 255 in
      Graphics.set_color (Graphics.rgb gray gray gray);
      Mygraphics.fill_rect gr xa ya xb yb
    done
  done;

  Graphics.set_color (Graphics.rgb 128 128 255);
  for i1 = 1 to n1 do
    Mygraphics.draw_poly_line gr [| (float i1 /. float n1, 0.0); (float i1 /. float n1, 1.0) |]
  done;
  for i2 = 0 to n2 do
    Mygraphics.draw_poly_line gr [| (0.0, float i2 /. float n2); (1.0, float i2 /. float n2) |]
  done;

  let axis_s = match axis with X -> "X" | Y -> "Y" | Z -> "Z" in
  let n_xyz, a, b =
    match axis with
    | X -> (world.nx, world.xa, world.xb)
    | Y -> (world.ny, world.ya, world.yb)
    | Z -> (world.nz, world.za, world.zb)
  in
  let xyz = get_xyz n n_xyz a b in
  Graphics.set_color (Graphics.rgb 0 0 0);
  Mygraphics.moveto gr 0.005 1.005;
  Graphics.draw_string (sprintf "%s %g" axis_s xyz);

  Graphics.synchronize ()

let print_2d_array world ary axis n =
  let n1, a1, b1, n2, a2, b2, ary_2d = get_2d_array world ary axis n in
  for i2 = 0 to n2 - 1 do
    for i1 = 0 to n1 - 1 do
      let xyz1, xyz2 = (get_xyz i1 n1 a1 b1, get_xyz i2 n2 a2 b2) in
      printf "%g\t%g\t%g\n" xyz1 xyz2 ary_2d.(i1).(i2)
    done;
    print_newline ()
  done

let gui world ary axis n =
  let v_min, v_max = get_min_max world ary in
  let axis = ref axis in
  let n = ref n in
  let n_max = ref world.nx in
  let open Graphics in
  let gr = Mygraphics.open_graph 500 550 0.0 0.0 1.0 (550.0 /. 500.0) in
  Graphics.auto_synchronize false;
  draw gr world v_min v_max ary !axis !n;
  loop_at_exit
    [Button_down; Key_pressed]
    (fun status ->
       if status.keypressed then begin
         begin
           match status.key with
           | 'x' ->
               axis := X;
               n := 0;
               n_max := world.nx
           | 'y' ->
               axis := Y;
               n := 0;
               n_max := world.ny
           | 'z' ->
               axis := Z;
               n := 0;
               n_max := world.nz
           | 'j' ->
               if !n > 0 then decr n
           | 'k' ->
               if !n < !n_max - 1 then incr n
           | 'p' ->
               print_2d_array world ary !axis !n;
               fprintf stderr "print done\n";
               flush stderr
           | 'q' ->
               raise Exit
           | _ ->
               ()
         end;
         draw gr world v_min v_max ary !axis !n
       end)

let cui world ary axis n =
  print_2d_array world ary axis n

let () =
  let world, ary = get_data stdin in
  match Array.length Sys.argv with
  | 1 ->
      gui world ary X 0
  | 3 ->
      let axis =
        match String.lowercase Sys.argv.(1) with
        | "x" -> X
        | "y" -> Y
        | "z" -> Z
        | _ -> failwith "invalid axis"
      in
      let n = int_of_string Sys.argv.(2) in
      cui world ary axis n
  | _ ->
      failwith "invalid argc"
