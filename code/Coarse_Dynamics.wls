(* ::Package:: *)

Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];
Needs["ThesisTools`"]


(*
This script lets me play with the dynamics of a system.
The first part generates the data. The second part is all
about the visualization.
*)

(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@ FIRST PART - MAXENT @@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)



A=MatrixExp[-l3*GObsMaxEnt[p,3]]/Tr[MatrixExp[-l3*GObsMaxEnt[p,3]]];
A//MatrixForm
(partialTraceB[A]) //MatrixForm
(partialTraceA[A]) //MatrixForm


(*There are two basic parameters for the coarse graining problem.
The z coordinate of the state and the swap probability*)
swapP = 0.5;
zcoord = 0.5;

(*Given those two, we can obtain a maximum entropy state*)
lambdastep=0.05;
rzlambda=RzLambdaTable[lambdastep,swapP];
ZMaxEnt=CGMaxEntStateLM[LagrangeMultFromZCoord[rzlambda,zcoord][[1]],swapP];
ZMaxEnt//MatrixForm


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@ SECOND PART - ASSIGNEMENT @@@@@@@@@@@@@@@@@@*)
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
mixedstates=UniformMixedStates[zcoord,1000];
assignements=AssignementMapForStateNotInZ[#,data]&/@mixedstates;
maxents=MaxEntForStateNotInZ[#,ZMaxEnt]&/@mixedstates;



(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@ THIRD PART - MAXENT DYNAMICS @@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
unitary=CNOT;

newmaxents=unitary . # . Dagger[unitary]&/@maxents;
coarsemaxents=coarseGraining2[#,swapP]&/@maxents;
evolcoarsemaxents=coarseGraining2[#,swapP]&/@newmaxents;
(*UnitaryEvolutionMaxEnt=ApplyUnitaryButSlowly[unitary,steps,assignements];*)
CoarseEvolutionMaxEnt=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolutionAss[[i]]],{i,1,Length[UnitaryEvolutionMaxEnt]}];


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@ THIRD PART - ASSIGNEMENT DYNAMICS @@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)

newassign=unitary . # . Dagger[unitary]&/@assignements;
coarseassign=coarseGraining2[#,swapP]&/@assignements;
evolcoarseassign=coarseGraining2[#,swapP]&/@newassign;
(*UnitaryEvolutionAss=ApplyUnitaryButSlowly[unitary,steps,assignements];*)
CoarseEvolutionAss=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolutionAss[[i]]],{i,1,Length[UnitaryEvolutionAss]}];


CompareTwoCoarseSets[set1_,set2_,legend_]:=Show[
ListPointPlot3D[
{densityMatrixToPoint[set1,gellMannBasis[1]],densityMatrixToPoint[set2,gellMannBasis[1]]},
BoxRatios->{1,1,1},
PlotRange->{{-1.,1.},{-1.`,1.`},{-1.,1.}},
PlotLegends->legend,
PlotLabel->Style["Avg dist="<>ToString[Mean[Table[Norm[set1[[i]]-set2[[i]],"Frobenius"],{i,Length[set1]}]]],15]
],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]
]


Export["../figures/"<>"AssVsMaxEnt_CNOT_init_"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",CompareTwoCoarseSets[coarseassign,coarsemaxents,{"AssMap","MaxEnt"}]]
Export["../figures/"<>"AssVsMaxEnt_CNOT_evol_"<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",CompareTwoCoarseSets[evolcoarseassign,evolcoarsemaxents,{"AssMap","MaxEnt"}]]


CompareTwoCoarseSets[evolcoarseassign,evolcoarsemaxents,{"AssMap","MaxEnt"}]


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
