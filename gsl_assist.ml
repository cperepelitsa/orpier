(*  Orpie -- a fullscreen RPN calculator for the console
 *  Copyright (C) 2003-2004, 2005, 2006-2007 Paul Pelzl
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License, Version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  Please send bug reports, patches, etc. to Paul Pelzl at 
 *  <pelzlpj@eecs.umich.edu>.
 *)

(* gsl_assist.ml
 *
 * Contains a bunch of helper functions that make ocamlgsl a bit easier
 * to work with. *)


open Gsl


let complex_zero = Gsl_complex.complex ~re:0. ~im:0.
let complex_one = Gsl_complex.complex ~re:1. ~im:0.


(* simplifying functions for use with ocamlgsl bindings *)
let cmpx_of_int i   = {Gsl_complex.re=Big_int.float_of_big_int i; Gsl_complex.im=0.0}
let cmpx_of_float f = {Gsl_complex.re=f; Gsl_complex.im=0.0}
let cmat_of_fmat fm =
   let rows, cols = Matrix.dims fm and
   f_array = Matrix.to_array fm in
   let c_array = Array.map cmpx_of_float f_array in
   Matrix_complex.of_array c_array rows cols



(* 1-norm of matrix *)
let one_norm mat =
   let n, m = Matrix.dims mat in
   let maxval = ref (-1.0) in
   let sum = ref 0.0 in
   for j = 0 to pred m do
      sum := 0.0;
      for i = 0 to pred n do
         sum := !sum +. (abs_float mat.{i, j})
      done;
      if !sum > !maxval then
         maxval := !sum
      else
         ()
   done;
   !maxval



(* solve a complex linear system using LU decomposition *)
let solve_complex_LU ?(protect=true) mat b =
  let mA = Vectmat.cmat_convert ~protect mat in
  let vB = (`CV (Vector_complex.copy b)) in
  let (len, _) = Vectmat.dims mA in
  let p = Permut.create len in
  let _ = Linalg.complex_LU_decomp mA p in
  let x = Vector_complex_flat.create len in
  Linalg.complex_LU_solve mA p vB (`CVF x);
  Vector_complex_flat.to_complex_array x



(* arch-tag: DO_NOT_CHANGE_a19e0df2-6d6b-4925-87eb-be2a2926ffbb *)
