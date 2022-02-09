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
(*@@@@@@@@@@@@@@@@@@@@@@@@@ FIRST PART @@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)


(*There are two basic parameters for the coarse graining problem.
The z coordinate of the state and the swap probability*)
swapP = 0.5;
zcoord = 0.8;

(*Given those two, we can obtain a maximum entropy state*)
lambdastep=0.05
zobs[l_,p_]:=-(1/2) Sech[l*p] Sech[l-p*l] (Sinh[l]+(1-2 p) Sinh[l-2 p*l])
rzlambda=Table[{zobs[l,swapP],l},{l,-4,0,lambdastep}]
nearpos[haystack_,value_]:= With[{ f = Nearest[haystack -> Range@Length@haystack]},f[value, 1]];
lagrangemult=With[{data=Transpose[rzlambda]},data[[2,nearpos[data[[1]],zcoord]]]];
G[i_,p_]:=p*KroneckerProduct[PauliMatrix[i],IdentityMatrix[2]]+(1-p)*KroneckerProduct[IdentityMatrix[2],PauliMatrix[i]]
A=-\[Lambda]3*G[3,p];
ExpMat=MatrixExp[A];
Z=Tr[ExpMat];
MaxEnt=ExpMat/Z;

(*First we generate or import the data. The GenerateMHData function
won't create anything if the file exists already*)
n = 1000;
beta = 100;
delta = 0.6;
steps=20;
zcoarsestate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;
GenerateMHData[n, beta, delta, swapP, zcoord]

(*This imports the data*)
filename="MHstates_n="<>ToString[n]<>"_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta];
data=Get["MH_data/"<>filename<>"_data"<>".m"];

(*We create mixed states with a radius of zcoord, then we
obtain their assignement map using the loaded data*)
mixedstates=UniformMixedStates[zcoord,1000];
assignements=AssignementMapForStateNotInZ[#,data]&/@mixedstates;

(*Now the dynamics. We apply an unitary to all the assignements*)
unitary=swapGate;
newedassign=unitary . # . Dagger[unitary]&/@assignements;
newcoarse=coarseGraining2[#,swapP]&/@newedassign;

UnitaryEvolution=ApplyUnitaryButSlowly[unitary,steps,assignements];
CoarseEvolution=Table[Map[coarseGraining2[#,swapP]&,UnitaryEvolution[[i]]],{i,1,Length[UnitaryEvolution]}];

(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@ SECOND PART @@@@@@@@@@@@@@@@@@@@@@@@@*)
(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)


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

*)






