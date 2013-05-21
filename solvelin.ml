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

(* solve a linear system Ax = b, with input nxn matrix A and output nx1
 * matrix b *)
let solve_linear (stack : rpc_stack) (evaln : int -> unit) =
   evaln 2;
   let gen_el2 = stack#pop () in
   let gen_el1 = stack#pop () in
   match gen_el1 with
   |RpcFloatMatrixUnit (el1, u1) ->
      begin match gen_el2 with
      |RpcFloatMatrixUnit (el2, u2) ->
         let n1, m1 = Matrix.dims el1 in
         if n1 <> m1 then
            (stack#push gen_el2;
            stack#push gen_el1;
            raise (Invalid_argument "multiplier matrix must be square"))
         else
            let n2, m2 = Matrix.dims el2 in
            if m2 <> 1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument "resultant matrix must be a column")
            end else if n2 <> m1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument ("dimensions of multiplier and " ^ 
               "resultant matrices do not match"))
            end else begin
               let b = Matrix.to_array el2 in
               let x = Linalg.solve_LU (`M el1) (`A b) in
               let x_mat = Matrix.of_array x m1 1 in
               stack#push (RpcFloatMatrixUnit (x_mat, Units.div u2 u1))
            end
      |RpcComplexMatrixUnit (el2, u2) ->
         let n1, m1 = Matrix.dims el1 in
         if n1 <> m1 then
            (stack#push gen_el2;
            stack#push gen_el1;
            raise (Invalid_argument "multiplier matrix must be square"))
         else
            let n2, m2 = Matrix_complex.dims el2 in
            if m2 <> 1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument "resultant matrix must be a column")
            end else if n2 <> m1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument ("dimensions of multiplier and" ^ 
               "resultant matrices do not match"))
            end else begin
               let a_cpx = Gsl_assist.cmat_of_fmat el1 in
               let b_arr = Matrix_complex.to_array el2 in
               let b_vec = Vector_complex.of_array b_arr in
               let x = Gsl_assist.solve_complex_LU (`CM a_cpx) b_vec in
               let x_mat = Matrix_complex.of_complex_array x m1 1 in
               stack#push (RpcComplexMatrixUnit (x_mat, Units.div u2 u1))
            end
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "both arguments of solve_linear must be matrices"))
      end
   |RpcComplexMatrixUnit (el1, u1) ->
      begin match gen_el2 with
      |RpcFloatMatrixUnit (el2, u2) ->
         let n1, m1 = Matrix_complex.dims el1 in
         if n1 <> m1 then begin
            stack#push gen_el2;
            stack#push gen_el1;
            raise (Invalid_argument "multiplier matrix must be square")
         end else
            let n2, m2 = Matrix.dims el2 in
            if m2 <> 1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument "resultant matrix must be a column")
            end else if n2 <> m1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument ("dimensions of multiplier and" ^ 
               "resultant matrices do not match"))
            end else begin
               let b_cpx = Gsl_assist.cmat_of_fmat el2 in
               let b_arr = Matrix_complex.to_array b_cpx in
               let b_vec = Vector_complex.of_array b_arr in
               let x = Gsl_assist.solve_complex_LU (`CM el1) b_vec in
               let x_mat = Matrix_complex.of_complex_array x m1 1 in
               stack#push (RpcComplexMatrixUnit (x_mat, Units.div u2 u1))
            end
      |RpcComplexMatrixUnit (el2, u2) ->
         let n1, m1 = Matrix_complex.dims el1 in
         if n1 <> m1 then begin
            stack#push gen_el2;
            stack#push gen_el1;
            raise (Invalid_argument "multiplier matrix must be square")
         end else
            let n2, m2 = Matrix_complex.dims el2 in
            if m2 <> 1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument "resultant matrix must be a column")
            end else if n2 <> m1 then begin
               stack#push gen_el1;
               stack#push gen_el2;
               raise (Invalid_argument ("dimensions of multiplier and" ^ 
               "resultant matrices do not match"))
            end else begin
               let b_arr = Matrix_complex.to_array el2 in
               let b_vec = Vector_complex.of_array b_arr in
               let x = Gsl_assist.solve_complex_LU (`CM el1) b_vec in
               let x_mat = Matrix_complex.of_complex_array x m1 1 in
               stack#push (RpcComplexMatrixUnit (x_mat, Units.div u2 u1))
            end
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "both arguments of solve_linear must be matrices"))
      end
   |_ ->
      (stack#push gen_el1;
      stack#push gen_el2;
      raise (Invalid_argument "both arguments of solve_linear must be matrices"))



(* arch-tag: DO_NOT_CHANGE_c11268a8-d37a-4573-98db-e985cc338e38 *)
