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

let get_xyz i n a b =
  a +. (float i) /. float (n - 1) *. (b -. a)

let print_2d axis n =
  let world, ary = get_data stdin in
  let flag_newline = ref false in
  for iz = 0 to world.nz - 1 do
    for iy = 0 to world.ny - 1 do
      flag_newline := false;
      for ix = 0 to world.nx - 1 do
        let i =
          match axis with
          | X -> ix
          | Y -> iy
          | Z -> iz
        in
        if n = i then begin
          flag_newline := true;
          let xyz1, xyz2 =
            match axis with
            | X -> (get_xyz iy world.ny world.ya world.yb,
                    get_xyz iz world.nz world.za world.zb)
            | Y -> (get_xyz ix world.nx world.xa world.xb,
                    get_xyz iz world.nz world.za world.zb)
            | Z -> (get_xyz ix world.nx world.xa world.xb,
                    get_xyz iy world.ny world.ya world.yb)
          in
          printf "%g\t%g\t%g\n" xyz1 xyz2 ary.(ix).(iy).(iz).v
        end
      done;
      if !flag_newline && (axis = Y || axis = Z) then print_newline ()
    done;
    if axis = X then print_newline ()
  done

let () =
  let axis =
    match String.uppercase Sys.argv.(1) with
    | "X" -> X
    | "Y" -> Y
    | "Z" -> Z
    | _ -> failwith "invalid axis"
  in
  let n = int_of_string Sys.argv.(2) in
  print_2d axis n
