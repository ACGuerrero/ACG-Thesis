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
swapP = 0.4;
zcoord = 0.9;
ztargetstate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;

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


testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{1.,1.,3.}] . Rest[gellMannBasis[1]])/2;


H=Pi/4*PauliMatrix[3];
unitary=KroneckerProduct[MatrixExp[I*H],IdentityMatrix[2]];
testMaxEnt=MaxEntForStateNotInZ[testTargetState,ZMaxEnt];
newtestMaxEnt=unitary . testMaxEnt . Dagger[unitary];
neweffective=CGKraus[newtestMaxEnt,swapP];
Norm/@densityMatrixToPoint[{testTargetState,neweffective},gellMannBasis[1]]
With[
{vecs=densityMatrixToPoint[{testTargetState,swapP*partialTraceB[testMaxEnt],(1-swapP)*partialTraceA[testMaxEnt],neweffective,swapP*partialTraceB[newtestMaxEnt],(1-swapP)*partialTraceA[newtestMaxEnt]},gellMannBasis[1]]},
	ShowWithBlochSphere[
		{
		ListPointPlot3D[
			{#}&/@vecs,
			PlotStyle->{Red,Green,Blue,Darker[Red],Darker[Green],Darker[Blue]},
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},
			AxesOrigin->{0,0,0},
			Ticks->None,
			AxesLabel->{"x","y","z"},
			PlotLegends->{HoldForm[Subscript[\[Rho],i]],HoldForm[Subscript[\[Rho],iA]],HoldForm[Subscript[\[Rho],iB]],HoldForm[Subscript[\[Rho],f]],HoldForm[Subscript[\[Rho],fA]],HoldForm[Subscript[\[Rho],fB]]}
		],
		Graphics3D[
		{Green,Arrowheads[Small],Arrow[{vecs[[3]],vecs[[2]]+vecs[[3]]}]}
		],
		Graphics3D[
		{Darker[Green],Arrowheads[Small],Arrow[{vecs[[3]],vecs[[3]]+vecs[[5]]}]}
		],
		Graphics3D[
		{Blue,Arrowheads[Small],Arrow[{{0,0,0},vecs[[3]]}]}
		],
		Graphics3D[
		{Darker[Green],Arrowheads[Small],Arrow[{{0,0,0},vecs[[5]]}]}
		]
		}
	]
]


Ham[t_]:=t*PauliMatrix[3];
U[t_]:=KroneckerProduct[MatrixExp[I*Ham[t]],IdentityMatrix[2]];
n=10;
points=densityMatrixToPoint[Table[CGKraus[U[t] . testMaxEnt . Dagger[U[t]],swapP],{t,0,Pi,N[Pi/n]}],gellMannBasis[1]];
points2=swapP*densityMatrixToPoint[(partialTraceB[#]&/@Table[U[t] . testMaxEnt . Dagger[U[t]],{t,0,Pi,N[Pi/n]}]),gellMannBasis[1]];
figure=With[
	{vecs=densityMatrixToPoint[{testTargetState,swapP*partialTraceB[testMaxEnt],(1-swapP)*partialTraceA[testMaxEnt]},gellMannBasis[1]]},
	ShowWithBlochSphere[{
		ListPointPlot3D[
			{#}&/@vecs,
			PlotStyle->{Red,Green,Blue},
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},
			AxesOrigin->{0,0,0},
			Ticks->None,
			AxesLabel->{"x","y","z"},
			PlotLegends->{\[Rho],HoldForm[Subscript[\[Rho],A]],HoldForm[Subscript[\[Rho],B]]}
		],
		Graphics3D[
		{Blue,Line[{{0,0,0},vecs[[3]]}]}
		],
		Graphics3D[
		{Darker[Green],Sequence@@(Line[{vecs[[3]],#}]&/@points)}
		],
		Graphics3D[
		{Darker[Green],Sequence@@(Line[{{0,0,0},#}]&/@points2)}
		],
		ListPointPlot3D[{points,points2},PlotStyle->{Red,Darker[Green]}]
	}]
]
Export["../figures/"<>"Ux1_H=sigmaz_sequence_"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>".png",figure]


MatrixExp[I*t*PauliMatrix[3]]


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
