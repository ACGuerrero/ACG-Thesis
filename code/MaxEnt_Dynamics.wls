#!/usr/bin/env wolframscript
(* ::Package:: *)

(*This script lets me play with the dynamics using the
maximum entropy state*)
Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
Needs["ThesisTools`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];


(*
===========================================
=CONTRUCTION OF THE MAXENT STATE AND ERROR=
===========================================
*)

(*There are two basic parameters for the coarse graining problem.
The z coordinate of the state and the swap probability*)
swapP = 0.01;
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

(*
swappedMaxEnt=SWAP[1] . ZMaxEnt . SWAP[1];
finalstate=CGKraus[swappedMaxEnt,swapP];
Print["Sphere contracts by a factor of: ", realfactor=With[{bv=densityMatrixToPoint[{targetstate,finalstate},gellMannBasis[1]]},bv[[2,3]]/bv[[1,3]]]]
Print["Theoretical contraction factor is: ", thefactor=SWAPContractionFactor[1,swapP,lagrangemult]]
Print["Relative error: ", 100*Abs[realfactor-thefactor]/thefactor, "%"]

*)


(*
===========================================
================  general  ================
===========================================
*)
(*We create mixed states with a radius of zcoord, then we
obtain their corresponding MaxEnt states*)
MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[zcoord,1500];
(*H=Pi/4*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]];
unitary=MatrixExp[I*H];*)
unitary=CNOT
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


lambdastep=0.05;
lowerlambda=-4;
upperlambda=0;
unitary=swapGate;
Do[Do[
	rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
	lagrangemult=LagrangeMultFromZCoord[rzlambdatable,z][[1]];
	ZMaxEnt=CGMaxEntStateLM[lagrangemult,p];
	MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[z,1500];
	NewMaxEntsNotInZ=unitary . # . Dagger[unitary]&/@MaxEntsNotInZ;
	CoarseMaxEntsNotInZ=Chop[CGKraus[#,p]]&/@MaxEntsNotInZ;
	CoarseNewMaxEntsNotInZ=Chop[CGKraus[#,p]]&/@NewMaxEntsNotInZ;
	Export["../figures/SWAP/"<>"MaxEnt_SWAP_t0vst1"<>"_p="<>ToString[NumberForm[p, {3, 2}]]<>"_z="<>ToString[NumberForm[z, {3, 2}]]<>".png",
	PlotTwoCoarseSets[
		CoarseMaxEntsNotInZ,CoarseNewMaxEntsNotInZ,
		{"Initial states","Final states"},"Effective evolution SWAP, z="<>ToString[z]<>", p="<>ToString[p]]
	];
	,{z,0.05,0.95,0.05}],{p,0.05,0.95,0.05}]


NumberForm[0.1, {3, 2}]


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
