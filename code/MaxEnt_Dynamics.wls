#!/usr/bin/env wolframscript
(* ::Package:: *)

(*This script lets me play with the dynamics using the
maximum entropy state*)
Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];
Needs["ThesisTools`"]


(*
===========================================
=CONTRUCTION OF THE MAXENT STATE AND ERROR=
===========================================
*)

(*There are two basic parameters for the coarse graining problem.
The z coordinate of the state and the swap probability*)
swapP = 0.8;
zcoord = 0.9;
targetstate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;

(*We then have to construct the maximum entropy state. Since we
have acces to the rz coordinate, we obtain the lagrange mult by
numerically approximating the inverse value*)
lambdastep=0.05;
lowerlambda=-4;
upperlambda=0;

rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]];

ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];

Print["The maximum entropy state is: ", ZMaxEnt//MatrixForm]
Print["Fidelity between the coarse maximum entropy state and the target state: ", fidelity[CGKraus[ZMaxEnt,swapP],targetstate]]


(*
===========================================
===============SWAP DYNAMICS===============
===========================================
*)
swappedMaxEnt=SWAP[1] . ZMaxEnt . SWAP[1];
finalstate=CGKraus[swappedMaxEnt,swapP];
Print["Sphere contracts by a factor of: ", realfactor=With[{bv=densityMatrixToPoint[{targetstate,finalstate},gellMannBasis[1]]},bv[[2,3]]/bv[[1,3]]]]
Print["Theoretical contraction factor is: ", thefactor=SWAPContractionFactor[1,swapP,lagrangemult]]
Print["Relative error: ", 100*Abs[realfactor-thefactor]/thefactor, "%"]


(* ::Input:: *)
(*\.b4*)
