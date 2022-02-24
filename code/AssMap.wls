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
fid=Table[
	Table[
		{n,With[
			{data=GenerateMHData[n, beta, delta, p, zcoord,False]},
		fidelity[Chop[Total[data]/Length[data]],ZMaxEnt]]},
	{n,100,5000,100}],
{p,0.1,0.5,0.1}];


Export["../figures/fidelityMaxEntAssMap_vs_n.png",Labeled[ListLinePlot[fid,PlotLegends->{"p=0.1","p=0.2","p=0.3","p=0.4","p=0.5"},Frame -> True, FrameLabel -> {Style["n",FontSize->12], Style["Fidelity",FontSize->12]}],"Fidelity between MaxEnt and AssMap",Top]]


Labeled[ListLinePlot[fid,Frame -> True, FrameLabel -> {Style["n",FontSize->12], Style["Fidelity",FontSize->12]}],"p="<>ToString[swapP],Top]


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@ THIRD PART - MAXENT DYNAMICS @@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
unitary=swapGate;

newmaxents=unitary . # . Dagger[unitary]&/@maxents;
coarsemaxents=coarseGraining2[#,swapP]&/@maxents;
evolcoarsemaxents=coarseGraining2[#,swapP]&/@newmaxents;
(*UnitaryEvolutionMaxEnt=ApplyUnitaryButSlowly[unitary,steps,assignements];*)
CoarseEvolutionMaxEnt=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolutionAss[[i]]],{i,1,Length[UnitaryEvolutionMaxEnt]}];

swapP = 0.1;
rztable=RzLambdaTable[swapP,-4,4,lambdastep=0.0001]
Show[ListPlot[Table[
With[
{bv=densityMatrixToPoint[{(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2,
CGKraus[swapGate . CGMaxEntStateLM[LagrangeMultFromZCoord[rztable,zcoord][[1]],swapP] . swapGate,swapP]},gellMannBasis[1]]},
{zcoord,bv[[2,3]]/bv[[1,3]]}],
{zcoord,0.1,1,0.01}]],Plot[SWAPContractionFactor[1,swapP,LagrangeMultFromZCoord[rztable,zcoord][[1]]],{zcoord,0.1,1}]]


Plot[Evaluate[Table[SWAPContractionFactor[1,p,LagrangeMultFromZCoord[rztable,zcoord][[1]]]/.{p->i},{i,0,0.5,0.1}]],{zcoord,0.1,1},AxesLabel->{Style[z,FontSize->18],Style[a,FontSize->18]},PlotLegends->{Table["p="<>ToString[i],{i,0,0.5,0.1}]}]


Show[Join[Table[Plot[SWAPContractionFactor[1,p,LagrangeMultFromZCoord[RzLambdaTable[swapP,-4,0,lambdastep=0.01],zcoord][[1]]],{zcoord,0.1,1}],{p,0,0.5,0.1}],Table[Plot[SWAPContractionFactor[1,p,LagrangeMultFromZCoord[RzLambdaTable[swapP,0,4,lambdastep=0.01],zcoord][[1]]],{zcoord,-1,-0.1}],{p,0,0.5,0.1}]]]


(* ::InheritFromParent:: *)
(*Graphics[{{RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[{Opacity[0.126], Thickness[Small]}], {}, {RGBColor[0.987148, 0.8073604000000001, 0.49470040000000004`], EdgeForm[{Opacity[0.126], Thickness[Small]}], Rectangle[{0.065, 0}, {0.07, 1000.}, RoundingRadius -> 0]}, {}, {}}, {}}, AspectRatio -> GoldenRatio^(-1), Axes -> {True, True}, AxesLabel -> {None, None}, AxesOrigin -> {-0.008, 0}, FrameLabel -> {{None, None}, {None, None}}, FrameTicks -> {{Automatic, Automatic}, {Automatic, Automatic}}, GridLines -> {None, None}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], PlotRange -> {{0., 0.4}, {All, All}}, PlotRangePadding -> {{Scaled[0.02], Scaled[0.02]}, {Scaled[0.02], Scaled[0.05]}}, Ticks -> {Automatic, Automatic}]*)


(coarsemaxents-evolcoarsemaxents)[[2]]//Norm


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@ THIRD PART - ASSIGNEMENT DYNAMICS @@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)

newassign=unitary . # . Dagger[unitary]&/@assignements;
coarseassign=coarseGraining2[#,swapP]&/@assignements;
evolcoarseassign=coarseGraining2[#,swapP]&/@newassign;
(*UnitaryEvolutionAss=ApplyUnitaryButSlowly[unitary,steps,assignements];*)
CoarseEvolutionAss=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolutionAss[[i]]],{i,1,Length[UnitaryEvolutionAss]}];


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
