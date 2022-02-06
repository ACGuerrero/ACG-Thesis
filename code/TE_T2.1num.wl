(* ::Package:: *)

Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
Needs["ThesisTools`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];


(*The following are actually pretty good parameters*)
n=1000;
beta=100;
delta=0.6;
swapP=0.3;
zcoord=0.8;
steps=20

zcoarsestate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;

(*Let's import the previously created data*)
filename="MHstates_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta];
data=Get["MH_data/"<>filename<>"_data"<>".m"];

(*Mixed states*)
mixedstates=UniformMixedStates[zcoord,1000];

(*Assignement maps and swap*)
assignements=AssignementMapForStateNotInZ[#,data]&/@mixedstates;
swappedassign=swapGate . # . Dagger[swapGate]&/@assignements;
(*
(*New coarse*)
swappedmixed=coarseGraining2[#,swapP]&/@swappedassign;

Export["../figures/coolimage.png",Show[ListPointPlot3D[densityMatrixToPoint[#,gellMannBasis[1]]&/@{mixedstates,swappedmixed},BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]]]
*)
SWAPEvolution=ApplyUnitaryButSlowly[swapGate,steps,assignements];
CGSWAPEvolution=Table[Map[coarseGraining2[#,swapP]&,SWAPEvolution[[i]]],{i,1,Length[SWAPEvolution]}];
(*
Table[Export["../figures/"<>"swap_evol"<>"_step"<>"_t="<>StringTake["0"<>ToString[i],-2]<>"_n"<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",visualizeBipartiteSystem[SWAPEvolution[[i]]]],{i,1,Length[SWAPEvolution]}];
*)
Table[
Export["../figures/"<>"coarse_swap_evol"<>"_step"<>"_t="<>StringTake["0"<>ToString[i],-2]<>"_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta]<>".png",
Show[
ListPointPlot3D[densityMatrixToPoint[CGSWAPEvolution[[i]],gellMannBasis[1]],BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]]],
{i,Length[CGSWAPEvolution]}
]





CGSWAPEvolution//Length


CGSWAPEvolution[[1]]


Export["../figures/TEST.png",Show[ListPointPlot3D[densityMatrixToPoint[CGSWAPEvolution[[1]],gellMannBasis[1]],BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]]]

