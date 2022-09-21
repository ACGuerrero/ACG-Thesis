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
r=0.9;
sp=0.9;
lag=l[r,sp,low=0.,up=6.,ste=0.002];
Do[
tim=i;
omega=Pi/4;
transformationIsing[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	{r[[1]]*Cos[omega*tim]^2-r[[1]]*Sin[omega*tim]^2-Sin[omega*tim]*Cos[omega*tim]*(sp*rb*runi[[3]]*(-ra*runi[[2]])+(1-sp)*ra*runi[[3]]*(-rb*runi[[2]])),
	r[[2]]*Cos[omega*tim]^2-r[[2]]*Sin[omega*tim]^2-Sin[omega*tim]*Cos[omega*tim]*(sp*rb*runi[[3]]*(ra*runi[[1]])+(1-sp)*ra*runi[[3]]*(rb*runi[[1]])),
	r[[3]]}
	];
	Export["../figures/Ising/"<>"sphere_Ising_t="<>ToString[i*100]<>"_r="<>ToString[r]<>"_p="<>ToString[sp]<>".png",FunctionSphereMesh2[r,transformationIsing,0.5]],
{i,0.0,1.0,0.1}
]


set=Table[Import["../figures/Ising/"<>"sphere_Ising_t="<>ToString[i*100]<>"_r="<>ToString[r]<>"_p="<>ToString[sp]<>".png"],{i,0.0,1.0,0.1}];
PNGSetToGif[set_]:=Export["../figures/"<>"Ising_GIF_t_evol"<>"_r="<>ToString[r]<>"_p="<>ToString[sp]<>".gif",
Flatten[{set, Table[set[[i]], {i, Length[set]}]}]]
PNGSetToGif[set]


transformationCNOT[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	r*(Cos[freq*tim]^4+Sin[freq*tim]^4)
	+2*sp*Sin[freq*tim]*Cos[freq*tim]*ra*(Sin[freq*tim]*Cos[freq*tim]*(rb*runi[[1]]*runi+(1-rb*runi[[1]])*{-runi[[1]],-runi[[2]],runi[[3]]})-(Cos[freq*tim]^2-Sin[freq*tim]^2)*(1+rb*runi[[1]])*{-runi[[2]],runi[[1]],0})
	+2*(1-sp)*Sin[freq*tim]*Cos[freq*tim]*rb*(Sin[freq*tim]*Cos[freq*tim]*(ra*runi[[3]]*runi+(1-ra*runi[[3]])*{runi[[1]],-runi[[2]],-runi[[3]]})-(Cos[freq*tim]^2-Sin[freq*tim]^2)*(1+ra*runi[[3]])*{0,-runi[[1]],runi[[2]]})
	];
transformationCNOTC[r_]:=With[{
	runi=Normalize[r],
	ra=Tanh[lag*sp],
	rb=Tanh[lag*(1-sp)]},
	{rb*(sp*ra*runi[[1]]^2+(1-sp)*runi[[1]]),
	rb*ra*(sp*runi[[1]]*runi[[2]]+(1-sp)*runi[[2]]*runi[[3]]),
	ra*(sp*runi[[3]]+(1-sp)*rb*runi[[3]]^2)}
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
