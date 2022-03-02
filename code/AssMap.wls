#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];
Needs["ThesisTools`"]


(*There are two basic parameters for the coarse graining problem.
The z coordinate of the state and the swap probability*)
swapP = 0.2;
zcoord = 0.5;
targetstate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@ ASSIGNMENT MAP CONSTRUCION @@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)

(*First we generate or import the data. The GenerateMHData function
won't create anything if the file exists already*)
n = 5000;
beta = 150;
delta = 0.6;
steps=20;
zcoarsestate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;
data=GenerateMHData[n, beta, delta, swapP, zcoord];
(*We create mixed states with a radius of zcoord, then we
obtain their assignement map using the loaded data*)
AssMap=Chop[Total[data]/Length[data]];
Print["The assigned state is: ", AssMap//MatrixForm]
Print["Fidelity between the assigned state and the target state: ", fidelity[CGKraus[AssMap,swapP],targetstate]]
AssMapsNotInZ=AssignementMapForStateNotInZ[#,data]&/@UniformMixedStates[zcoord,1000];


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@ COMPARISON BETWEEN MAXENT AND ASSMAP @@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
lambdastep=0.05;
lowerlambda=-4;
upperlambda=0;
rzlambdatable=RzLambdaTable[swapP,lowerlambda,upperlambda,lambdastep];
lagrangemult=LagrangeMultFromZCoord[rzlambdatable,zcoord][[1]];
ZMaxEnt=CGMaxEntStateLM[lagrangemult,swapP];
Print["The maximum entropy state is: ", ZMaxEnt//MatrixForm]
Print["Fidelity between the coarse maximum entropy state and the target state: ", fidelity[CGKraus[ZMaxEnt,swapP],targetstate]]
(*As a function of n*)
fidn=Table[
	Table[
		{n,With[
			{data=GenerateMHData[n, beta, delta, p, zcoord,False]},
		fidelity[Chop[Total[data]/Length[data]],ZMaxEnt]]},
	{n,100,5000,100}],
{p,0.1,0.5,0.1}];
(*As a function of p*)
fidp=Table[
		{p,With[
			{data=GenerateMHData[1000, beta, delta, p, zcoord,False]},
		fidelity[Chop[Total[data]/Length[data]],ZMaxEnt]]},
	{p,0.1,0.5,0.05}];


Export["../figures/fidelityMaxEntAssMap_vs_n.png",Labeled[ListLinePlot[fid,PlotLegends->{"p=0.1","p=0.2","p=0.3","p=0.4","p=0.5"},Frame -> True, FrameLabel -> {Style["n",FontSize->12], Style["Fidelity",FontSize->12]}],"Fidelity between MaxEnt and AssMap",Top]]


Export["../figures/fidelityMaxEntAssMap_vs_p.png",Labeled[ListLinePlot[fidp,Frame -> True, FrameLabel -> {Style["p",FontSize->12], Style["Fidelity",FontSize->12]}],"Fidelity between MaxEnt and AssMap",Top]]


CoarseNewAssMapsNotInZ[[1]]


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@ DYNAMICS @@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
unitary=swapGate;
NewAssMap=unitary . AssMap . Dagger[unitary];
NewAssMapsNotInZ=unitary . # . Dagger[unitary]&/@AssMapsNotInZ;
CoarseAssMapsNotInZ=Chop[CGKraus[#,swapP]]&/@AssMapsNotInZ;
CoarseNewAssMapsNotInZ=Chop[CGKraus[#,swapP]]&/@NewAssMapsNotInZ;
differences=ToSphericalCoordinates[densityMatrixToPoint[CoarseAssMapsNotInZ,gellMannBasis[1]]]-ToSphericalCoordinates[densityMatrixToPoint[CoarseNewAssMapsNotInZ,gellMannBasis[1]]];
(*UnitaryEvolutionAss=ApplyUnitaryButSlowly[unitary,steps,assignements];
CoarseEvolutionAss=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolutionAss[[i]]],{i,1,Length[UnitaryEvolutionAss]}];*)
Histogram[#,{-2*Pi,2*Pi,0.05}]&/@Transpose[differences]
PlotTwoCoarseSets[CoarseAssMapsNotInZ,CoarseNewAssMapsNotInZ,{"Initial states","Final states"},"Evolution of the Assignement Map"]


Export["../figures/"<>"MaxEnt_SWAP_t0vst1_"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",
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


(*This is by default all commented. These are commands that you can 
run in a new line to visualize special parts of the dynamics*)

(*
--------> CREATE A GIF OF THE COARSE EVOLUTION

gif = Table[
Labeled[
Show[
ListPointPlot3D[densityMatrixToPoint[CoarseEvolution[[i]],gellMannBasis[1]],BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]
],
{"t="<>ToString[i],"Coarse SWAP for p="<>ToString[swapP]<>", z="<>ToString[zcoord]},
{Top,Bottom}], 
{i,Length[CoarseEvolution]}];
Export["../figures/"<>"coarse_swap_evol_"<>ToString[steps]<>"steps"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".gif",
Flatten[{gif, Table[gif[[i]], {i, Length[gif] }]}]]

--------> CREATE A GIF OF THE FINE EVOLUTION
gif = Table[
Labeled[
visualizeBipartiteSystem[UnitaryEvolution[[i]]],
{"t="<>ToString[i],"Fine SWAP for p="<>ToString[swapP]<>", z="<>ToString[zcoord]},
{Top,Bottom}],
{i,1,Length[UnitaryEvolution]}];
Export["../figures/"<>"swap_evol_"<>ToString[steps]<>"steps"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".gif",
Flatten[{gif, Table[gif[[i]], {i, Length[gif] }]}]]

--------> CREATE PNG FOR EACH STEP FINE
Table[
Export["../figures/"<>"swap_evol"<>"_step"<>"_t="<>StringTake["0"<>ToString[i],-2]<>"_n"<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".gif",
Labeled[
visualizeBipartiteSystem[SWAPEvolution[[i]]],
"t="<>ToString[i],
Top]
]
],
{i,1,Length[SWAPEvolution]}];


--------> CREATE PNG FOR EACH STEP COARSE
Table[
Export["../figures/"<>"coarse_swap_evol"<>"_step"<>"_t="<>StringTake["0"<>ToString[i],-2]<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",
Labeled[
Show[
ListPointPlot3D[densityMatrixToPoint[CGSWAPEvolution[[i]],gellMannBasis[1]],BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]
],
"t="<>ToString[i],
Top]
],
{i,Length[CGSWAPEvolution]}
]


--------> Compare MaxEnt and Assignement

*)
