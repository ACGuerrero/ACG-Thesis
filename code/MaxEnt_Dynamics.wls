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
swapP = 0.99;
zcoord = 0.8;
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

(*
swappedMaxEnt=SWAP[1] . ZMaxEnt . SWAP[1];
finalstate=CGKraus[swappedMaxEnt,swapP];
Print["Sphere contracts by a factor of: ", realfactor=With[{bv=densityMatrixToPoint[{targetstate,finalstate},gellMannBasis[1]]},bv[[2,3]]/bv[[1,3]]]]
Print["Theoretical contraction factor is: ", thefactor=SWAPContractionFactor[1,swapP,lagrangemult]]
Print["Relative error: ", 100*Abs[realfactor-thefactor]/thefactor, "%"]

*)


(* ::Input:: *)
(*\.b4*)


(*
===========================================
================  general  ================
===========================================
*)
(*We create mixed states with a radius of zcoord, then we
obtain their corresponding MaxEnt states*)
MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[zcoord,1500];


H=Pi/4*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]];
unitary=MatrixExp[I*H];

NewMaxEnt=unitary . ZMaxEnt . Dagger[unitary];
NewMaxEntsNotInZ=unitary . # . Dagger[unitary]&/@MaxEntsNotInZ;
CoarseMaxEntsNotInZ=Chop[CGKraus[#,swapP]]&/@MaxEntsNotInZ;
CoarseNewMaxEntsNotInZ=Chop[CGKraus[#,swapP]]&/@NewMaxEntsNotInZ;
(*
differences=ToSphericalCoordinates[densityMatrixToPoint[CoarseMaxEntsNotInZ,gellMannBasis[1]]]-ToSphericalCoordinates[densityMatrixToPoint[CoarseNewMaxEntsNotInZ,gellMannBasis[1]]];
UnitaryEvolutionAss=ApplyUnitaryButSlowly[unitary,steps,assignements];
CoarseEvolutionAss=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolutionAss[[i]]],{i,1,Length[UnitaryEvolutionAss]}];
Histogram[#,{-2*Pi,2*Pi,0.05}]&/@Transpose[differences]*)
PlotTwoCoarseSets[CoarseMaxEntsNotInZ,CoarseNewMaxEntsNotInZ,{"Initial states","Final states"},"Evolution of the Assignement Map"]


(*Export["../figures/"<>"MaxEnt_SWAP_t0vst1_"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",
	CompareTwoCoarseSets[
		coarseassign,
		evolcoarsemaxents,
		{ToExpression["\\rho", TeXForm, HoldForm],"\[ScriptCapitalC](S"<>ToString[ToExpression["\\varrho", TeXForm, HoldForm]]<>"S)"},"Effective evolution using MaxEnt"
	]
]

Export["../figures/"<>"MapAss_SWAP_t0vst1_"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",
	CompareTwoCoarseSets[
		coarsemaxents,
		evolcoarseassign,
		{ToExpression["\\rho", TeXForm, HoldForm],"\[ScriptCapitalC](S"<>ToString[ToExpression["\\varrho", TeXForm, HoldForm]]<>"S)"},"Effective evolution using AssMap"
	]
]
*)
