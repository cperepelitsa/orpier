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
open Big_int

let mult (stack : rpc_stack) (evaln : int -> unit) =
   evaln 2;
   let gen_el2 = stack#pop () in
   let gen_el1 = stack#pop () in
   match gen_el1 with
   |RpcInt el1 -> (
      match gen_el2 with
      |RpcInt el2 ->
         stack#push (RpcInt (mult_big_int el1 el2))
      |RpcFloatUnit (el2, uu2) ->
         let f_el1 = float_of_big_int el1 in
         stack#push (RpcFloatUnit (f_el1 *. el2, uu2))
      |RpcComplexUnit (el2, uu2) ->
         let c_el1 = {
            Gsl_complex.re = float_of_big_int el1;
            Gsl_complex.im = 0.0
         } in
         stack#push (RpcComplexUnit (Gsl_complex.mul c_el1 el2, uu2))
      |RpcFloatMatrixUnit (el2, uu) ->
         let result = Matrix.copy el2 in
         Matrix.scale result (float_of_big_int el1);
         stack#push (RpcFloatMatrixUnit (result, uu))
      |RpcComplexMatrixUnit (el2, uu) ->
         let c_el1 = cmpx_of_int el1 in
         let result = Matrix_complex.copy el2 in
         Matrix_complex.scale result c_el1;
         stack#push (RpcComplexMatrixUnit (result, uu))
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "incompatible types for multiplication"))
      )
   |RpcFloatUnit (el1, uu1) -> (
      match gen_el2 with
      |RpcInt el2 ->
         let f_el2 = float_of_big_int el2 in
         stack#push (RpcFloatUnit (el1 *. f_el2, uu1))
      |RpcFloatUnit (el2, uu2) ->
         stack#push (RpcFloatUnit (el1 *. el2, Units.mult uu1 uu2))
      |RpcComplexUnit (el2, uu2) ->
         let c_el1 = c_of_f el1 in
         stack#push (RpcComplexUnit (Gsl_complex.mul c_el1 el2, Units.mult uu1 uu2))
      |RpcFloatMatrixUnit (el2, uu2) ->
         let uprod = Units.mult uu1 uu2 in
         let result = Matrix.copy el2 in
         Matrix.scale result el1;
         stack#push (RpcFloatMatrixUnit (result, uprod))
      |RpcComplexMatrixUnit (el2, uu2) ->
         let uprod = Units.mult uu1 uu2 in
         let result = Matrix_complex.copy el2 in
         Matrix_complex.scale result (c_of_f el1);
         stack#push (RpcComplexMatrixUnit (result, uprod))
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "incompatible types for multiplication"))
      )
   |RpcComplexUnit (el1, uu1) -> (
      match gen_el2 with
      |RpcInt el2 ->
         let c_el2 = cmpx_of_int el2 in
         stack#push (RpcComplexUnit (Gsl_complex.mul el1 c_el2, uu1))
      |RpcFloatUnit (el2, uu2) ->
         let c_el2 = c_of_f el2 in
         stack#push (RpcComplexUnit (Gsl_complex.mul el1 c_el2, Units.mult uu1 uu2))
      |RpcComplexUnit (el2, uu2) ->
         stack#push (RpcComplexUnit (Gsl_complex.mul el1 el2, Units.mult uu1 uu2))
      |RpcFloatMatrixUnit (el2, uu2) ->
         let c_el2 = cmat_of_fmat el2 in
         Matrix_complex.scale c_el2 el1;
         stack#push (RpcComplexMatrixUnit (c_el2, Units.mult uu1 uu2))
      |RpcComplexMatrixUnit (el2, uu2) ->
         let result = Matrix_complex.copy el2 in
         Matrix_complex.scale result el1;
         stack#push (RpcComplexMatrixUnit (result, Units.mult uu1 uu2))
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "incompatible types for multiplication"))
      )
   |RpcFloatMatrixUnit (el1, uu1) -> (
      match gen_el2 with
      |RpcInt el2 ->
         let result = Matrix.copy el1 in
         Matrix.scale result (float_of_big_int el2);
         stack#push (RpcFloatMatrixUnit (result, uu1))
      |RpcFloatUnit (el2, uu2) ->
         let result = Matrix.copy el1 in
         Matrix.scale result el2;
         stack#push (RpcFloatMatrixUnit (result, Units.mult uu1 uu2))
      |RpcComplexUnit (el2, uu2) ->
         let c_el1 = cmat_of_fmat el1 in
         Matrix_complex.scale c_el1 el2;
         stack#push (RpcComplexMatrixUnit (c_el1, Units.mult uu1 uu2))
      |RpcFloatMatrixUnit (el2, uu2) ->
         let n1, m1 = (Matrix.dims el1)
         and n2, m2 = (Matrix.dims el2) in
         if m1 = n2 then
            let result = Matrix.create n1 m2 in
            Blas.gemm Blas.NoTrans Blas.NoTrans
            1.0 el1 el2 0.0 result;
            stack#push (RpcFloatMatrixUnit (result, Units.mult uu1 uu2))
         else
            (stack#push gen_el1;
            stack#push gen_el2;
            raise (Invalid_argument "incompatible matrix dimensions for multiplication"))
      |RpcComplexMatrixUnit (el2, uu2) ->
         let n1, m1 = (Matrix.dims el1)
         and n2, m2 = (Matrix_complex.dims el2) in
         if m1 = n2 then
            let c_el1 = cmat_of_fmat el1
            and result = Matrix_complex.create n1 m2 in
            Blas.Complex.gemm Blas.NoTrans Blas.NoTrans
            complex_one c_el1 el2 complex_zero result;
            stack#push (RpcComplexMatrixUnit (result, Units.mult uu1 uu2))
         else
            (stack#push gen_el1; 
            stack#push gen_el2;
            raise (Invalid_argument "incompatible matrix dimensions for multiplication"))
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "incompatible types for multiplication"))
      )
   |RpcComplexMatrixUnit (el1, uu1) -> (
      match gen_el2 with
      |RpcInt el2 ->
         let c_el2 = cmpx_of_int el2 in
         let result = Matrix_complex.copy el1 in
         Matrix_complex.scale result c_el2;
         stack#push (RpcComplexMatrixUnit (result, uu1))
      |RpcFloatUnit (el2, uu2) ->
         let result = Matrix_complex.copy el1 in
         Matrix_complex.scale result complex_one;
         stack#push (RpcComplexMatrixUnit (result, Units.mult uu1 uu2))
      |RpcComplexUnit (el2, uu2) ->
         let result = Matrix_complex.copy el1 in
         Matrix_complex.scale result complex_one;
         stack#push (RpcComplexMatrixUnit (result, Units.mult uu1 uu2))
      |RpcFloatMatrixUnit (el2, uu2) ->
         let n1, m1 = (Matrix_complex.dims el1)
         and n2, m2 = (Matrix.dims el2) in
         if m1 = n2 then
            let c_el2 = cmat_of_fmat el2
            and result = Matrix_complex.create n1 m2 in
            Blas.Complex.gemm Blas.NoTrans Blas.NoTrans
            complex_one el1 c_el2 complex_zero result;
            stack#push (RpcComplexMatrixUnit (result, Units.mult uu1 uu2))
         else
            (stack#push gen_el1;
            stack#push gen_el2;
            raise (Invalid_argument "incompatible matrix dimensions for multiplication"))
      |RpcComplexMatrixUnit (el2, uu2) ->
         let n1, m1 = (Matrix_complex.dims el1)
         and n2, m2 = (Matrix_complex.dims el2) in
         if m1 = n2 then
            let result = Matrix_complex.create n1 m2 in
            Blas.Complex.gemm Blas.NoTrans Blas.NoTrans
            complex_one el1 el2 complex_zero result;
            stack#push (RpcComplexMatrixUnit (result, Units.mult uu1 uu2))
         else
            (stack#push gen_el1;
            stack#push gen_el2;
            raise (Invalid_argument "incompatible matrix dimensions for multiplication"))
      |_ ->
         (stack#push gen_el1;
         stack#push gen_el2;
         raise (Invalid_argument "incompatible types for multiplication"))
      )
   |_ ->
      (stack#push gen_el1;
      stack#push gen_el2;
      raise (Invalid_argument "incompatible types for multiplication"))


(* arch-tag: DO_NOT_CHANGE_5fc03e41-d1d3-40da-8b68-9a85d96148d0 *)
