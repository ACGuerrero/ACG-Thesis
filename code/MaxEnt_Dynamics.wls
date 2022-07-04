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
zcoord = 0.8;
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
==  Numerical evolution of MaxEnt states  ==
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


(*
===========================================
==  Numerical evolution of MaxEnt states GIF!  ==
===========================================
*)
swapP = 0.5;
zcoord = 1.;
lambdastep=0.05;
lowerlambda=0.1;
upperlambda=4.;
rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]]
ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];
MaxEntsNotInZ=MaxEntForStateNotInZ[#,ZMaxEnt]&/@UniformMixedStates[zcoord,1500];
steps=100;
(*unitary[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];*)
unitary[t_]:=MatrixExp[-I*t*{{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}}]
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


(*
===========================================
 Effective state as sum of partial traces
===========================================
*)


(*There are two basic parameters for the coarse graining problem.
The z coordinate of the state and the swap probability*)
swapP = 0.8;
zcoord = 0.8;
ztargetstate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;

(*We then have to construct the maximum entropy state. Since we
have acces to the rz coordinate, we obtain the lagrange mult by
numerically approximating the inverse value*)
lambdastep=0.05;
lowerlambda=0.1;
upperlambda=4.;

rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]];

ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];


testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{4.,4,7}] . Rest[gellMannBasis[1]])/2;
testMaxEnt=MaxEntForStateNotInZ[testTargetState,ZMaxEnt];

EffAsSum=With[
{vecs=densityMatrixToPoint[{testTargetState,swapP*partialTraceB[testMaxEnt],(1-swapP)*partialTraceA[testMaxEnt]},gellMannBasis[1]]},
	ShowWithBlochSphere[
		{
		ListPointPlot3D[
			{#}&/@vecs,
			PlotStyle->{{Red,PointSize[0.02]},{Darker[Green],PointSize[0.001]},{Darker[Blue],PointSize[0.001]}},
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{0,1.},{0,1.},{0,1.}},
			AxesOrigin->{0,0,0},
			Ticks->Automatic,
			AxesLabel->{"x","y","z"},
			LabelStyle->Black,
			PlotLegends->Placed[{HoldForm[\[Rho]],HoldForm[Subscript[\[Rho],A]],HoldForm[Subscript[\[Rho],B]]},{0.7,0.7}],
 ViewPoint -> {0.4, -1.7, 0.5}
		],
		Graphics3D[
		{Darker[Green],Arrowheads[0.03], Thickness[0.003], Arrow[{{0,0,0},vecs[[2]]}]}
		],
		Graphics3D[
		{Darker[Blue],Arrowheads[0.03],Thickness[0.003],Arrow[{vecs[[2]],vecs[[2]]+vecs[[3]]}]}
		]
		}
	]
];


Export["../figures/"<>"effective_state_as_sum"<>".png",EffAsSum];


(*
===========================================
==== Trajectories of rho, rhoA and rhoB ===
===========================================
*)
testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{-0.5,-0.5,.90}] . Rest[gellMannBasis[1]])/2;
testMaxEnt=MaxEntForStateNotInZ[testTargetState,ZMaxEnt];
(*Ham1=1*PauliMatrix[3];
Ham2=10*PauliMatrix[3];
U1[t_]:=SWAP[t];
(*U1[t_]:=KroneckerProduct[MatrixExp[t*I*Ham1],MatrixExp[t*I*Ham2]];*)
U1[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];*)
(*U1[t_]:=MatrixPower[CNOT,t]*)
U1[t_]:=MatrixExp[I*t*Pi*KroneckerProduct[PauliMatrix[3],IdentityMatrix[2]]];
n=400;
swapP = 0.7;
end=1;
points=densityMatrixToPoint[Table[CGKraus[U1[t] . testMaxEnt . Dagger[U1[t]],swapP],{t,0,end,N[end/n]}],gellMannBasis[1]];
points2=swapP*densityMatrixToPoint[(partialTraceB[#]&/@Table[U1[t] . testMaxEnt . Dagger[U1[t]],{t,0,end,N[end/n]}]),gellMannBasis[1]];
points3=(1-swapP)*densityMatrixToPoint[(partialTraceA[#]&/@Table[U1[t] . testMaxEnt . Dagger[U1[t]],{t,0,end,N[end/n]}]),gellMannBasis[1]];
figure=With[
	{vecs=densityMatrixToPoint[{testTargetState,swapP*partialTraceB[testMaxEnt],(1-swapP)*partialTraceA[testMaxEnt]},gellMannBasis[1]]},
	ShowWithBlochSphere[{
		ListPointPlot3D[
			{#}&/@vecs,
			PlotStyle->{{Darker[Red],PointSize[0.0002]},{Darker[Green]},{Blue,PointSize[0.02]}},
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},
			AxesOrigin->{0,0,0},
			AxesStyle->Black,
			Ticks->None,
			AxesLabel->{"x","y","z"},
			LabelStyle->Black,
			PlotLegends->Placed[{\[Rho][t],HoldForm[Subscript[\[Rho],A]][t],HoldForm[Subscript[\[Rho],B]][t]},{0.9,0.7}],
			ViewPoint->{2, -2, 1}
		],
		Graphics3D[
		{Darker[Green], Thickness[0.003], Line[{{0,0,0},vecs[[2]]}]}
		],
		Graphics3D[{Blue,Thickness[0.002],Opacity[0.2],Line/@Transpose[{points2,points}]}],
		Graphics3D[{Blue,Thickness[0.002],Opacity[0.2],Line/@Transpose[{Table[{0,0,0},Length[points3]],points3}]}],
		ListPointPlot3D[{points,points2,points3},PlotStyle->{{Darker[Red],Opacity[0.5]},{Darker[Green],Opacity[0.2]},{Darker[Blue],Opacity[0.3]}}]
	}]
]
Export["../figures/"<>"U1xU2_H1=sz_H2=Id"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_sequence"<>".png",figure];


Export["../figures/"<>"U1xU2_H1=100(sx-sy)_H2=sz"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_sequence"<>".png",figure];

Show[Graphics3D[{Darker[Green],Opacity[1],Thickness[0.001],Line/@Transpose[{points3,points}]}],
ListPointPlot3D[points3]],
Graphics3D[{Darker[Red],Opacity[1],Thickness[0.001],Line/@Transpose[{Table[{0,0,0},Length[points]],points}]}];


Export["../figures/"<>"Ux1_H=1x+2y_sequence"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>".png",figure];


(*
===========================================
Small oscilations around no error trayectory
===========================================
*)

zcoord=0.8;
testTargetState=(IdentityMatrix[2]+zcoord*Normalize[{1,1,2}] . Rest[gellMannBasis[1]])/2;
testMaxEnt=MaxEntForStateNotInZ[testTargetState,ZMaxEnt];
testother=KroneckerProduct[testTargetState,IdentityMatrix[2]];
(*
U1[t_]:=MatrixExp[I*t*Pi*KroneckerProduct[PauliMatrix[3],IdentityMatrix[2]]];
U1[t_]:=SWAP[t];
(*U1[t_]:=KroneckerProduct[MatrixExp[t*I*Ham1],MatrixExp[t*I*Ham2]];*)
U1[t_]:=MatrixExp[I*t*KroneckerProduct[PauliMatrix[3],PauliMatrix[3]]];*)
(*U1[t_]:=MatrixPower[CNOT,t]*)
Ham1=(-PauliMatrix[3]);
Ham2=15*(-PauliMatrix[2]+2PauliMatrix[1])/Sqrt[2];
U1[t_]:=KroneckerProduct[MatrixExp[t*I*Pi*Ham1],MatrixExp[t*I*Pi*Ham2]];
n=6000;
swapP = 0.9;
end=1.;
points=densityMatrixToPoint[Table[CGKraus[U1[t] . testMaxEnt . Dagger[U1[t]],swapP],{t,0,end,N[end/n]}],gellMannBasis[1]];
points2=densityMatrixToPoint[((partialTraceB[#]/2)&/@Table[U1[t] . testother . Dagger[U1[t]],{t,0,end,N[end/n]}]),gellMannBasis[1]];
figure=With[
	{vecs=densityMatrixToPoint[{testTargetState,Chop[partialTraceB[testother]]/2},gellMannBasis[1]]},
	ShowWithBlochSphere[{
		ListPointPlot3D[
			{#}&/@vecs,
			PlotStyle->{{Darker[Red],PointSize[0.002]},{Black,PointSize[0.002]}},
			BoxRatios->{1, 1, 1},
			Boxed->False,
			PlotRange->{{-1.1,1.1},{-1.1,1.1},{-1.1,1.1}},
			AxesOrigin->{0,0,0},
			Ticks->None,
			AxesLabel->{"x","y","z"},
			LabelStyle->Black,
			PlotLegends->Placed[{\[Rho][t],HoldForm[Subscript[\[Rho],A][t]],HoldForm[Subscript[\[Rho],B]]},{0.9,0.7}],
			ViewPoint -> {-0.4, -0.4, 0.5}
		],
		ListPointPlot3D[{points,points2},PlotStyle->{{Darker[Red],Opacity[1],PointSize[0.005]},{Black,PointSize[0.005]}}]
	}]
]

Export["../figures/"<>"U1xU2_H1=(sz)_H2=15(sx-sy)"<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_far.png",%];


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
sp=0.05;
Do[
sp=i;
runi=Normalize[{-1,1,1}];
rvec=r*runi;
lag=l[r,sp,low=0.,up=6.,ste=0.002];
ra=Tanh[lag*sp]*runi;
rb=Tanh[lag*(1-sp)]*runi;
transformationCNOTC[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	{rb*(sp*ra*runi[[1]]^2+(1-sp)*runi[[1]]),
	rb*ra*(sp*runi[[1]]*runi[[2]]+(1-sp)*runi[[2]]*runi[[3]]),
	ra*(sp*runi[[3]]+(1-sp)*rb*runi[[3]]^2)}
	];
	Export["../figures/CNOT/"<>"sphere_CNOT_t=1.0"<>"_z="<>ToString[r]<>"_p="<>ToString[i*100]<>".png",FunctionSphereMesh2[r,transformationCNOTC,0.5]],
{i,0.0,1.0,0.01}
]


transformationCNOT[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	r*(Cos[freq*tim]^4+Sin[freq*tim]^4)
	+2*sp*Sin[freq*tim]*Cos[freq*tim]*ra*(Sin[freq*tim]*Cos[freq*tim]*(rb*runi[[1]]*runi+(1-rb*runi[[1]])*{-runi[[1]],-runi[[2]],runi[[3]]})-(Cos[freq*tim]^2-Sin[freq*tim]^2)*(1+rb*runi[[1]])*{-runi[[2]],runi[[1]],0})
	+2*(1-sp)*Sin[freq*tim]*Cos[freq*tim]*rb*(Sin[freq*tim]*Cos[freq*tim]*(ra*runi[[3]]*runi+(1-ra*runi[[3]])*{runi[[1]],-runi[[2]],-runi[[3]]})-(Cos[freq*tim]^2-Sin[freq*tim]^2)*(1+ra*runi[[3]])*{0,-runi[[1]],runi[[2]]})
	];


FunctionSphereMesh2[r_,function_,opa_:1.]:=Show[
	ParametricPlot3D[function[{r*Sin[t] Cos[p], r*Sin[t] Sin[p], r*Cos[t]}], {t, 0, \[Pi]}, {p, 0, 2 \[Pi]}, Lighting -> {"Neutral", White},
	PlotStyle -> {GrayLevel[.5],Opacity[opa]}, PlotTheme -> None,
			BoxRatios->{1, 1, 1},
			PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},
			AxesLabel->{Style["x", 16],Style["y", 16],Style["z", 16]},
			LabelStyle->Black,
			PlotLegends->Placed[Style["p="<>ToString[sp],16],{0.65,0.9}]],
	ParametricPlot3D[function[{r*Cos[t],r*Sin[t],0}],{t,0,2*Pi},PlotStyle->Red],
	ListPointPlot3D[(Map[function[#]&,{{0,0,r},{0,0,-r}}]),PlotStyle->{Red,PointSize[0.03]}]]


set=Table[Import["../figures/CNOT/"<>"sphere_CNOT_t=1.0"<>"_z="<>ToString[r]<>"_p="<>ToString[i*100]<>".png"],{i,0.0,1.0,0.02}];


PNGSetToGif[set_]:=Export["../figures/"<>"CNOT_GIF_p_evol"<>"_r="<>ToString[r]<>".gif",
Flatten[{set, Table[set[[i]], {i, Length[set]}]}]]
PNGSetToGif[set]


(* ::Input:: *)
(*\.b4*)


ToString[0.02*100]


r=0.9;
sp=0.1;
runi=Normalize[{1,0.01,0.01}];
rvec=r*runi;
lag=l[r,sp,low=0.,up=6.,ste=0.002];
ra=Tanh[lag*sp]*runi;
rb=Tanh[lag*(1-sp)]*runi;
freq=Pi/4;
tim=0.0;
transformationBitFli[r_]:=With[{
	runi=Normalize[r]},
	{runi[[1]],
	(1-2*sp)*runi[[2]],
	(1-2*sp)*runi[[3]]}
	];
transformationCNOT[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	r*(Cos[freq*tim]^4+Sin[freq*tim]^4)
	+2*sp*Sin[freq*tim]*Cos[freq*tim]*ra*(Sin[freq*tim]*Cos[freq*tim]*(rb*runi[[1]]*runi+(1-rb*runi[[1]])*{-runi[[1]],-runi[[2]],runi[[3]]})-(Cos[freq*tim]^2-Sin[freq*tim]^2)*(1+rb*runi[[1]])*{-runi[[2]],runi[[1]],0})
	+2*(1-sp)*Sin[freq*tim]*Cos[freq*tim]*rb*(Sin[freq*tim]*Cos[freq*tim]*(ra*runi[[3]]*runi+(1-ra*runi[[3]])*{runi[[1]],-runi[[2]],-runi[[3]]})-(Cos[freq*tim]^2-Sin[freq*tim]^2)*(1+ra*runi[[3]])*{0,-runi[[1]],runi[[2]]})
	];
Export["../../figures/CNOT/"<>"sphere_CNOT_t="<>ToString[tim]<>"_r="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh[r,transformationCNOT,0.5]]


l[r_,p_,lo_,up_,st_]:=First@LagrangeMultFromZCoord[
	RzLambdaTable[p,lo,up,st],
	r];
r=0.9;
sp=0.5;
runi=Normalize[{-1,1,1}];
rvec=r*runi;
lag=l[r,sp,low=0.,up=6.,ste=0.002];
ra=Tanh[lag*sp]*runi;
rb=Tanh[lag*(1-sp)]*runi;
transformationCNOTC[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	{rb*(sp*ra*runi[[1]]^2+(1-sp)*runi[[1]]),
	rb*ra*(sp*runi[[1]]*runi[[2]]+(1-sp)*runi[[2]]*runi[[3]]),
	ra*(sp*runi[[3]]+(1-sp)*rb*runi[[3]]^2)}
	];
Export["../../figures/CNOT/"<>"sphere_CNOT_t=1.0"<>"_z="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh2[r,transformationCNOTC,0.5]];


l[r_,p_,lo_,up_,st_]:=First@LagrangeMultFromZCoord[
	RzLambdaTable[p,lo,up,st],
	r];
r=0.9;
sp=0.5;
runi=Normalize[{-1,1,1}];
rvec=r*runi;
lag=l[r,sp,low=0.,up=6.,ste=0.002];
ra=Tanh[lag*sp]*runi;
rb=Tanh[lag*(1-sp)]*runi;
transformationCNOTC[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	{rb*(sp*ra*runi[[1]]^2+(1-sp)*runi[[1]]),
	rb*ra*(sp*runi[[1]]*runi[[2]]+(1-sp)*runi[[2]]*runi[[3]]),
	ra*(sp*runi[[3]]+(1-sp)*rb*runi[[3]]^2)}
	];
Export["../../figures/CNOT/"<>"sphere_CNOT_t=1.0"<>"_z="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh2[r,transformationCNOTC,0.5]];


Norm[transformationCNOT[runi]]


r=0.9;
sp=0.5;
omega=Pi/4;
tim=1.;
transformationIsingp2[r_]:=With[{
	runi=Normalize[r]},
	{r[[1]]*Cos[2*omega*tim]-r[[3]]*r[[2]]*Sin[2*omega*tim],
	r[[2]]*Cos[2*omega*tim]+r[[3]]*r[[1]]*Sin[2*omega*tim],
	r[[3]]}
	];
Export["../../figures/Ising/"<>"sphere_Ising_t="<>ToString[tim]<>"_z="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh2[r,transformationIsingp2,1.]];


{r[[1]]*Cos[2*omega*t]-r[[3]]*r[[2]]*Sin[2*omega*t],
	r[[2]]*Cos[2*omega*t]+r[[3]]*r[[1]]*Sin[2*omega*t],
	r[[3]]}


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



