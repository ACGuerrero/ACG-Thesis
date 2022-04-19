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
swapP = 0.8;
zcoord = 0.6;
ztargetstate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;

(*We then have to construct the maximum entropy state. Since we
have acces to the rz coordinate, we obtain the lagrange mult by
numerically approximating the inverse value*)
lambdastep=0.05;
lowerlambda=0.1;
upperlambda=4.;

rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]]

ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];

Print["The maximum entropy state is: ", ZMaxEnt//MatrixForm]
Print["Fidelity between the coarse maximum entropy state and the target state: ", fidelity[CGKraus[ZMaxEnt,swapP],ztargetstate]]


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
Clear[unitary]
(*We create mixed states with a radius of zcoord, then we
obtain their corresponding MaxEnt states*)
MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[zcoord,1500];
(*H=Pi/4*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]];
unitary=MatrixExp[I*H];*)
unitary[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];
NewMaxEnt=unitary . ZMaxEnt . Dagger[unitary[1]];
NewMaxEntsNotInZ=unitary[1] . # . Dagger[unitary[1]]&/@MaxEntsNotInZ;
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
r=zcoord;
swapP = 0.6;
rltable=RzLambdaTable[r,lowerlambda,upperlambda,lambdastep];
l=LagrangeMultFromZCoord[rltable,r][[1]];
unit1[t_]:=MatrixExp[-I*t*Pi*(PauliMatrix[1]+PauliMatrix[2])/Sqrt[2]];
unit2[t_]:=IdentityMatrix[2];
rot1[t_]:=SU2ToSO3[unit1[t]];
rot2[t_]:=SU2ToSO3[unit2[t]];
transformation[t_]:=(swapP*Tanh[-l*swapP]*rot1[t]+(1-swapP)*Tanh[-l*(1-swapP)]*rot2[t])
t=0.5;
figure=TransformedSphereMesh[1,transformation[t]]
Export["../figures/"<>"U1xU2_H1=Pi(sx+sy)_H2=Id"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"t="<>ToString[t]<>".png",figure]


swapP = 0.7;
zcoord = 0.8;
lambdastep=0.05;
lowerlambda=0.1;
upperlambda=4.;
rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]]
ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];
MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[zcoord,1500];
steps=100;
(*unitary[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];*)
unitary[t_]:=MatrixPower[CNOT,t]
end=N[Pi];
pngset=Table[
	ShowWithBlochSphere[{
		ListPointPlot3D[
			densityMatrixToPoint[Chop[CGKraus[#,swapP]&/@(unitary[t] . # . Dagger[unitary[t]]&/@MaxEntsNotInZ)],gellMannBasis[1]],
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},
			AxesOrigin->{0,0,0},
			Ticks->None,
			AxesLabel->{"x","y","z"},ViewPoint->Front]
	}],
	{t,0,end,end/steps}];
PNGSetToGif[set_]:=Export["../figures/"<>"CNOT_effevol_FRONT"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>".gif",
Flatten[{set, Table[set[[i]], {i, Length[set]}]}]]
PNGSetToGif[pngset]


H=Pi*(IdentityMatrix[4]-
	KroneckerProduct[PauliMatrix[3],IdentityMatrix[2]]-
	KroneckerProduct[IdentityMatrix[2],PauliMatrix[1]]+
	KroneckerProduct[PauliMatrix[3],PauliMatrix[1]])/4;
MatrixExp[-I*H]//MatrixForm


(*
CoarseEvolutionPNGSet[zcoord_,swapP_,]:= With[
{coarseev=Table[Map[coarseGraining2[#,swapP]&,ApplyUnitaryButSlowly[unitary,steps,assignements][[i]]],{i,1,Length[ApplyUnitaryButSlowly[unitary,steps,assignements]]}];
},
Table[
Labeled[
Show[
ListPointPlot3D[densityMatrixToPoint[coarseev[[i]],gellMannBasis[1]],BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]
],
{"t="<>ToString[i],"Coarse evolution for p="<>ToString[swapP]<>", z="<>ToString[zcoord]},
{Top,Bottom}], 
{i,Length[coarseev]}]];

PNGSetToGif[set_]:=Export["../figures/"<>"coarse_swap_evol_"<>ToString[steps]<>"steps"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".gif",
Flatten[{gif, Table[gif[[i]], {i, Length[gif] }]}]]
*)


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


testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{-2.,-5,2}] . Rest[gellMannBasis[1]])/2;
H=Pi/4*(PauliMatrix[1]+PauliMatrix[2]);
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
];


testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{-0.3,-0.3,1.0}] . Rest[gellMannBasis[1]])/2;
testMaxEnt=MaxEntForStateNotInZ[testTargetState,ZMaxEnt];
(*Ham1=1*PauliMatrix[3];
Ham2=10*PauliMatrix[3];
U1[t_]:=SWAP[t];
(*U1[t_]:=KroneckerProduct[MatrixExp[t*I*Ham1],MatrixExp[t*I*Ham2]];*)
U1[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];*)
(*U1[t_]:=MatrixPower[CNOT,t]*)
U1[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],IdentityMatrix[2]]];
n=20;
swapP = 0.9;
end=1;
points=densityMatrixToPoint[Table[CGKraus[U1[t] . testMaxEnt . Dagger[U1[t]],swapP],{t,0,end,N[end/n]}],gellMannBasis[1]];
points2=swapP*densityMatrixToPoint[(partialTraceB[#]&/@Table[U1[t] . testMaxEnt . Dagger[U1[t]],{t,0,end,N[end/n]}]),gellMannBasis[1]];
points3=(1-swapP)*densityMatrixToPoint[(partialTraceA[#]&/@Table[U1[t] . testMaxEnt . Dagger[U1[t]],{t,0,end,N[end/n]}]),gellMannBasis[1]];
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
		Graphics3D[{Darker[Blue],Thickness[0.001],Opacity[0.2],Line/@Transpose[{points2,points}]}],
		Graphics3D[{Darker[Blue],Thickness[0.001],Opacity[0.5],Line/@Transpose[{Table[{0,0,0},Length[points3]],points3}]}],
		Graphics3D[{Black,Line[{{0,0,0},vecs[[1]]}]}],
		ListPointPlot3D[{points,points2,points3},PlotStyle->{{Darker[Red],Opacity[0.5]},{Darker[Green],Opacity[0.2]},{Darker[Blue],Opacity[0.2]}}]
	}]
]


Export["../figures/"<>"U1xU2_H1=100(sx-sy)_H2=sz"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_sequence"<>".png",figure];

Show[Graphics3D[{Darker[Green],Opacity[1],Thickness[0.001],Line/@Transpose[{points3,points}]}],
ListPointPlot3D[points3]],
Graphics3D[{Darker[Red],Opacity[1],Thickness[0.001],Line/@Transpose[{Table[{0,0,0},Length[points]],points}]}];


Export["../figures/"<>"Ux1_H=1x+2y_sequence"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>".png",figure];


Chop[partialTraceB[testother]/2]//MatrixForm


testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{0.3,-0.6,0.80}] . Rest[gellMannBasis[1]])/2;
testMaxEnt=MaxEntForStateNotInZ[testTargetState,ZMaxEnt];
testother=KroneckerProduct[testTargetState,IdentityMatrix[2]];
(*
U1[t_]:=MatrixExp[I*t*Pi*KroneckerProduct[PauliMatrix[3],IdentityMatrix[2]]];
U1[t_]:=SWAP[t];
(*U1[t_]:=KroneckerProduct[MatrixExp[t*I*Ham1],MatrixExp[t*I*Ham2]];*)
U1[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];*)
(*U1[t_]:=MatrixPower[CNOT,t]*)
Ham1=(-PauliMatrix[3]+PauliMatrix[2]);
Ham2=10*(-PauliMatrix[2]+2PauliMatrix[1]);
U1[t_]:=KroneckerProduct[MatrixExp[t*I*Pi*Ham1],MatrixExp[t*I*Pi*Ham2]];
n=5000;
swapP = 0.9;
end=1.;
points=densityMatrixToPoint[Table[CGKraus[U1[t] . testMaxEnt . Dagger[U1[t]],swapP],{t,0,end,N[end/n]}],gellMannBasis[1]];
points2=densityMatrixToPoint[((partialTraceB[#]/2)&/@Table[U1[t] . testother . Dagger[U1[t]],{t,0,end,N[end/n]}]),gellMannBasis[1]];
figure=With[
	{vecs=densityMatrixToPoint[{testTargetState,Chop[partialTraceB[testother]]/2},gellMannBasis[1]]},
	ShowWithBlochSphere[{
		ListPointPlot3D[
			{#}&/@vecs,
			PlotStyle->{Red,Black},
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},
			AxesOrigin->{0,0,0},
			Ticks->None,
			AxesLabel->{"x","y","z"},
			PlotLegends->{\[Rho],HoldForm[Subscript[\[Rho],A]],HoldForm[Subscript[\[Rho],B]]}
		],
		ListPointPlot3D[{points,points2},PlotStyle->{{Darker[Red],Opacity[1],PointSize[0.005]},{Black,PointSize[0.005]}}]
	}]
]


Export["../figures/"<>"U1xU2_H1=(sy-sz)_H2=10(2sx-sy)"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>".png",%];


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





l[r_,p_,lo_,up_,st_]:=First@LagrangeMultFromZCoord[
	RzLambdaTable[p,lo,up,st],
	r];
r=0.8;
sp=0.95;
runi=Normalize[{-1,1,1}];
rvec=r*runi;
lag=l[r,sp,low=0.,up=6.,ste=0.002];
ra=Tanh[lag*sp]*runi;
rb=Tanh[lag*(1-sp)]*runi;
transformation2[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	{rb*(sp*ra*runi[[1]]^2+(1-sp)*runi[[1]]),
	rb*ra*(sp*runi[[1]]*runi[[2]]+(1-sp)*runi[[2]]*runi[[3]]),
	ra*(sp*runi[[3]]+(1-sp)*rb*runi[[3]]^2)}
	];
Export["../figures/CNOT/"<>"sphere_CNOT_t=0.0"<>"_z="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh[r,(IdentityMatrix[3] . #)&,0.5]]
Export["../figures/CNOT/"<>"sphere_CNOT_t=1.0"<>"_z="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh[r,transformation2,0.5]]


Clear[unitary]
swapP=0.8;
zcoord=0.9;
lambdastep=0.05;
lowerlambda=0.1;
upperlambda=4.;
rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]]
ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];
MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[zcoord,1500];
unitary[t_]:=CNOT
NewMaxEnt=unitary . ZMaxEnt . Dagger[unitary[1]];
NewMaxEntsNotInZ=unitary[1] . # . Dagger[unitary[1]]&/@MaxEntsNotInZ;
CoarseMaxEntsNotInZ=Chop[CGKraus[#,swapP]]&/@MaxEntsNotInZ;
CoarseNewMaxEntsNotInZ=Chop[CGKraus[#,swapP]]&/@NewMaxEntsNotInZ;
PlotTwoCoarseSets[CoarseMaxEntsNotInZ,CoarseNewMaxEntsNotInZ,{"Initial states","Final states"},"Evolution of the MaxEnt"]



