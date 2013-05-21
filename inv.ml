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

open Rpc_stack
open Gsl
open Gsl_assist


let inv (stack : rpc_stack) (evaln : int -> unit) =
   evaln 1;
   let gen_el = stack#pop () in
   match gen_el with
   |RpcFloatUnit (el, uu) ->
      stack#push (RpcFloatUnit (1.0 /. el, Units.div Units.empty_unit uu))
   |RpcComplexUnit (el, uu) ->
      stack#push (RpcComplexUnit 
         (Gsl_complex.inverse el, Units.div Units.empty_unit uu))
   |RpcFloatMatrixUnit (el, uu) ->
      let new_unit = Units.pow uu (~-. 1.0) in
      let n, m = (Matrix.dims el) in
      if n = m then
         let copy_el1 = Matrix.copy el
         and copy_el2 = Matrix.copy el
         and perm     = Permut.create m
         and inv      = Matrix.create m m
         and vv       = Matrix.create m m
         and ss       = Vector.create m
         and work     = Vector.create m in
         begin
            (* test for singular matrix *)
            (* first factor out matrix norm, since the GSL SVD algorithm has
             * issues with large magnitude matrices *)
            let norm = Gsl_assist.one_norm copy_el1 in
            Matrix.scale copy_el1 (1.0 /. norm);
            (* now compute condition number as largest singular value
             * divided by smallest singular value *)
            Linalg._SV_decomp (`M copy_el1) (`M vv) (`V ss) (`V work);
            let condition_number = 
               (Vector.get ss 0) /. (Vector.get ss (pred m)) 
            in
            (* if the condition number is too large for machine precision,
             * then abort with error *)
            if condition_number > 1e14 then
               (stack#push gen_el;
               raise (Invalid_argument "cannot invert ill-conditioned matrix"))
            else
               let _ = Linalg._LU_decomp (`M copy_el2) perm in
               (Linalg._LU_invert (`M copy_el2) perm (`M inv);
               stack#push (RpcFloatMatrixUnit (inv, new_unit)))
         end
      else
         (stack#push gen_el;
         raise (Invalid_argument "cannot invert non-square matrix"))
   |RpcComplexMatrixUnit (el, uu) ->
      let new_unit = Units.pow uu (~-. 1.0) in
      let n, m = (Matrix_complex.dims el) in
      if n = m then
         let copy_el = Vectmat.cmat_convert ~protect:true (`CM el) and
         perm = Permut.create m and
         inv = Matrix_complex.create m m in
         try
            let _ = Linalg.complex_LU_decomp copy_el perm in
            Linalg.complex_LU_invert copy_el perm (`CM inv);
            stack#push (RpcComplexMatrixUnit (inv, new_unit))
         with Error.Gsl_exn _ ->
            (stack#push gen_el;
            raise (Invalid_argument "cannot invert singular matrix"))
      else
         (stack#push gen_el;
         raise (Invalid_argument "cannot invert non-square matrix"))
   |_ ->
      (stack#push gen_el;
      raise (Invalid_argument "inversion is undefined for this data
      type"))



(* arch-tag: DO_NOT_CHANGE_d8ce074c-3d77-4448-b3c6-9e239b853aad *)
