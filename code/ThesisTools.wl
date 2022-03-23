(* ::Package:: *)

BeginPackage["ThesisTools`"]
Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
Pauli2Basis::usage="something"
CNOT::usage="something"
MatrixToLatex::usage="something"
ShowOperatorsWithSphere::usage="something"
KetsToBlochVector::usage="something"
UnitaryStep::usage="something"
ApplyUnitaryButSlowly::usage="something"
UniformMixedStates::usage="something"
MixedToPure::usage="something"
UB2Q::usage="something"
UnitaryToRotateFineStates::usage="something"
AssignementMapForStateNotInZ::usage="something"
RotatedFineStates::usage="something"
testDistribution::usage="something"
metropolisHastingsSample::usage="something"
cgfidelitylist::usage="something"
cgfrobeniuslist::usage="something"
GenerateMHData::usage="GenerateMHData[n,beta,delta,swapP,zcoord,print_:True] generates a set of pure states using the Metropolis Hastings algorithm. If the data exists, it imports it."
GObsMaxEnt::usage="something"
CGMaxEntStateLM::usage="something"
NearestPosition::usage="something"
GObsMaxEnt::usage="something"
CGMaxEntStateLM::usage="something"
ZCoordFromLagrangeMult::usage="ZCoordFromLagrangeMult[\[Lambda],p] obtains the z coordinate of the effective state given the lagrange multiplier"
RzLambdaTable::usage="RzLambdaTable[p,low,up,step] a table of {z,\[Lambda]} with \[Lambda] ranging from low to up with step step."
LagrangeMultFromZCoord::usage="Finds the lagrange multiplier associated with a p and z values. To be used along RzLambdaTable like LagrangeMultFromZCoord[RzLambdaTable[p,low,up,step],z]"
MaxEntForStateNotInZ::usage="something"
CGKrausOp::usage="something"
CGKraus::usage="CGKraus[rho,p] applies the coarse graining map using its Kraus Operators"
SWAPContractionFactor::usage="SWAPContractionFactor[t,p,\[Lambda]] gives the contraction factor of the coarse system resulting from applying the swap gate to a MaxEnt state characterized by p,\[Lambda]."
SWAP::usage="SWAP[t] applies the operator at a time t. t=1 is the full swap gate, while t=0 is the identity operator."
PlotTwoCoarseSets::usage="PlotTwoCoarseSets[set1,set2,legend,title] takes two sets of two level density operators and plots their bloch vectors."
PlotTwoCoarseSetsWLine::usage="PlotTwoCoarseSets[set1,set2,legend,title] takes two sets of two level density operators and plots their bloch vectors, with a line joining corresponding points."
ShowWithBlochSphere::usage="Acts like Show function, but appends a Graphics3D showing the bloch sphere. Also, argument is a list."
Begin["`Private`"]


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@GENERAL TOOLS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)


Pauli2Basis=Flatten[Table[Pauli[{i,j}],{i,0,3},{j,0,3}],1];
CNOT={{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}};
NearestPosition[haystack_,value_]:= With[{ f = Nearest[haystack -> Range@Length@haystack]},f[value, 1]];

MatrixToLatex[matrix_]:=
"\\begin{pmatrix}"<>Fold[
(#1<>"\\"<>"\\"<>#2)&,
Map[
Fold[(ToString[#1]<>"&"<>ToString[#2])&,#]&,matrix
]
]<>"\\end{pmatrix}"

ShowOperatorsWithSphere[operators_]:=
Show[
ListPointPlot3D[
densityMatrixToPoint[operators,gellMannBasis[1]],BoxRatios->{1, 1, 1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}}],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]]

KetsToBlochVector[kets_,basis_]:=
densityMatrixToPoint[ketsToDensity[kets],basis]

UnitaryStep[unitary_,t_]:=MatrixPower[unitary,t];

ApplyUnitaryButSlowly[unitary_,steps_,rhos_]:=
With[
{stepsize=1/steps},
Table[
MapMonitored[
UnitaryStep[unitary,i*stepsize] . # . Dagger[UnitaryStep[unitary,i*stepsize]]&,
rhos
],
{i,0,steps}
]
]

UniformMixedStates[r_,n_]:=
With[
{basis=gellMannBasis[1]},
Map[
(IdentityMatrix[2]+# . Rest[basis])/2&,
r*KetsToBlochVector[Table[RandomState[2],n],basis]
]
]

MixedToPure[mixedstate_]:=
{Cos[#[[1]]/2],Exp[I*#[[2]]]*Sin[#[[1]]/2]}&[
Rest[ToSphericalCoordinates[densityMatrixToPoint[{mixedstate},gellMannBasis[1]][[1]]]]
]

UB2Q[qubit2_,qubit1_]:=
With[
{u1={{qubit1[[1]],-Conjugate[qubit1[[2]]]},{qubit1[[2]],Conjugate[qubit1[[1]]]}}
,u2={{Conjugate[qubit2[[1]]],Conjugate[qubit2[[2]]]},{-qubit2[[2]],qubit2[[1]]}}},u1 . u2]

UnitaryToRotateFineStates[cstate_]:=
KroneckerProduct[#,#]&[
UB2Q[{1,0},MixedToPure[cstate]]
]

AssignementMapForStateNotInZ[cstate_,fdata_]:=
With[
{preassignement=Total[fdata/Length[fdata]],
Ubig=UnitaryToRotateFineStates[cstate]
},
Ubig . preassignement . Dagger[Ubig]
]

RotatedFineStates[cstate_,fdata_]:=
With[
{preassignement=Total[fdata/Length[fdata]],
Ubig=UnitaryToRotateFineStates[cstate]
},
Ubig . #&/@fdata
]

(*All about the MaxEnt state*)
GObsMaxEnt[p_,i_]:=p*KroneckerProduct[PauliMatrix[i],IdentityMatrix[2]]+(1-p)*KroneckerProduct[IdentityMatrix[2],PauliMatrix[i]];
CGMaxEntStateLM[lambda_,p_]:=With[{ExpMat=MatrixExp[-lambda*GObsMaxEnt[p,3]]},ExpMat/Tr[ExpMat]]
ZCoordFromLagrangeMult[l_,p_]:=-(p*Tanh[l*p]+(1-p)*Tanh[l*(1-p)]);
RzLambdaTable[p_,low_,up_,step_]:=Transpose[Table[{ZCoordFromLagrangeMult[l,p],l},{l,low,up,step}]];
LagrangeMultFromZCoord[data_,zcoord_]:=data[[2,NearestPosition[data[[1]],zcoord]]]
MaxEntForStateNotInZ[cstate_,ZMaxEnt_]:=
With[
{Ubig=UnitaryToRotateFineStates[cstate]
},
Ubig . ZMaxEnt . Dagger[Ubig]
]

CGKrausOp[p_]:={
Sqrt[p]KroneckerProduct[IdentityMatrix[2],{1,0}],
Sqrt[p]KroneckerProduct[IdentityMatrix[2],{0,1}],
Sqrt[1-p]KroneckerProduct[IdentityMatrix[2],{1,0}] . swapGate,
Sqrt[1-p]KroneckerProduct[IdentityMatrix[2],{0,1}] . swapGate
};

CGKraus[rho_,p_]:=With[{kr=CGKrausOp[p]},Total[(# . rho . ConjugateTranspose[#])&/@kr]];

SWAPContractionFactor[t_,p_,\[Lambda]_]:=((E^(-I \[Pi] t - \[Lambda] - 
  2 p \[Lambda]) (2 E^(
    I \[Pi] t + 2 p \[Lambda]) (-1 + E^(2 \[Lambda])) - 
   E^(2 \[Lambda]) (1 + E^(2 I \[Pi] t)) (-1 + 2 p) + 
   E^(4 p \[Lambda]) (1 + E^(2 I \[Pi] t)) (-1 + 
      2 p)))/(4 (Sinh[\[Lambda]] + (1 - 2 p) Sinh[\[Lambda] - 
      2 p \[Lambda]])))//Chop;

SWAP[t_]:={{1,0,0,0},{0,(1+Exp[I*Pi*t])/2,(1-Exp[I*Pi*t])/2,0},{0,(1-Exp[I*Pi*t])/2,(1+Exp[I*Pi*t])/2,0},{0,0,0,1}}

PlotTwoCoarseSets[set1_,set2_,legend_,title_]:=Show[
ListPointPlot3D[
{densityMatrixToPoint[set1,gellMannBasis[1]],densityMatrixToPoint[set2,gellMannBasis[1]]},
BoxRatios->{1,1,1},
PlotRange->{{-1.,1.},{-1.`,1.`},{-1.,1.}},
PlotLegends->legend,
PlotLabel->title,
AxesLabel->{"x","y","z"}
],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]
]

PlotTwoCoarseSetsWLine[set1_,set2_,legend_,title_]:=
With[{points={densityMatrixToPoint[set1,gellMannBasis[1]],densityMatrixToPoint[set2,gellMannBasis[1]]}},
Show[
ListPointPlot3D[points,BoxRatios->{1,1,1},PlotRange->{{-1.,1.},{-1.,1.},{-1.,1.}},PlotLegends->legend,PlotLabel->title],
Graphics3D[{Opacity[0.2],ColorData[1][1],Thickness[0.005],Line/@Transpose[points]}],
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True]
]
]


ShowWithBlochSphere[args_]:=Show[
Sequence@@args,
Graphics3D[{Opacity[0.2],GrayLevel[0.9],Sphere[]},BoxRatios->1,Axes->True,
AxesLabel->{"x","y","z"}]
]


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


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@TOOLS FROM ADANERICK@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)

testDistribution[beta_,targetstate_,state_,swapP_]:= Exp[-beta*Norm[targetstate - coarseGraining2[state,swapP],"Frobenius"]]

metropolisHastingsSample[size_,\[Beta]_,\[Delta]_,swapP_,initialstate_,targetstate_]:= Module[{n = 0, X = initialstate, Y, U, \[Alpha], statelist = {}},
	While[n < size,
		(*Y = ketsToDensity[randomKets[4,1]][[1]];*)
		U = randomSmallEvolution[4,\[Delta]];
		Y = U . X . ConjugateTranspose[U];
		\[Alpha] = Min[testDistribution[\[Beta],targetstate,Y,swapP]/testDistribution[\[Beta],targetstate,X,swapP],1];
		X = RandomChoice[{\[Alpha], 1 - \[Alpha]}->{Y,X}];
		If[X == Y, AppendTo[statelist,X];n++]
	];
	Return[statelist]
	]

cgfidelitylist[bigstateslist_,targetstate_,swapP_]:=fidelity[coarseGraining2[#,swapP],targetstate]&/@bigstateslist
cgfrobeniuslist[bigstateslist_,targetstate_,swapP_]:=Norm[coarseGraining2[#,swapP]-targetstate,"Frobenius"]&/@bigstateslist

CreateData[n_,beta_,delta_,swapP_,zcoord_]:=Module[{
targetstate,
data,
filename,
thertestdist,
time
},
targetstate=(IdentityMatrix[2]+zcoord*PauliMatrix[3])/2;
{time,data}=Timing[metropolisHastingsSample[n+400,beta,delta,swapP,ketsToDensity[randomKets[4,1]][[1]],targetstate][[401;;]]];
thertestdist=cgfrobeniuslist[data,targetstate,swapP];
Print["Data created."];
Print["Process took "<>ToString[time/60]<>" minutes"];
Print["Mean distance between images and target: ", Mean[thertestdist]];
Print["Number of states in the file: ",Length[data]];
Return[data]
]

GenerateMHData[n_,beta_,delta_,swapP_,zcoord_,print_:True]:=
Module[{
targetstate,
data,
filename,
thertestdist,
time,
size
},
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];
filename="MHstates_z="<>ToString[zcoord]<>"_p="<>ToString[swapP]<>"_beta="<>ToString[beta]<>"_delta="<>ToString[delta];
If[FileExistsQ["/home/acastillo/Documents/tesis-adan/code/MH_data/"<>filename<>"_data"<>".m"],
	If[print,Print["Set of data exists already."]];
	data=Get["MH_data/"<>filename<>"_data"<>".m"];
	size=Length[data];
	If[print,Print["Data imported from file. File contains "<>ToString[size]<>" elements."]];
	If[size>=n,
		If[print,Print["File has enough data"]];
		Return[data[[;;n]]],
		If[print,Print["Set of data doesn't have enough data. Creating more data."]];
		data=Join[CreateData[n-Length[data],beta,delta,swapP,zcoord],data];
		Export["MH_data/"<>filename<>"_data"<>".m",data];
		Return[data]
	],
	If[print,Print["Set of data doesn't exist. Creating data."]];
	data=CreateData[n,beta,delta,swapP,zcoord];
	Export["MH_data/"<>filename<>"_data"<>".m",data];
	Return[data]
]
];
End[]
EndPackage[]
