(* ::Package:: *)

BeginPackage["ThesisTools`"]
Needs["CoolTools`"]
Needs["Carlos`"]
Needs["Quantum`"]
Pauli2Basis::usage="something"
MatrixToLatex::usage="something"
ShowOperatorsWithSphere::usage="something"
KetsToBlochVector::usage="something"
Begin["`Private`"]


(*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@GENERAL TOOLS@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*)


Pauli2Basis=Flatten[Table[Pauli[{i,j}],{i,0,3},{j,0,3}],1];

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

UnitaryStep[unitary_,t_]:=
With[{unitaryH=-I*MatrixLog[unitary]},MatrixExp[I*t*unitaryH]]

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

End[]
EndPackage[]
