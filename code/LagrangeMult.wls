#!/usr/bin/env wolframscript
(* ::Package:: *)

$Assumptions=Element[{\[Lambda]3,p,rz},Reals] && 0<=p<=1 && rz\[Element]Reals && 0<=rz<=1;
Z=MatrixExp[-\[Lambda]3*G[3,p]];
MaxEnt=Z/Tr[Z];
G[i_,p_]:=p*KroneckerProduct[PauliMatrix[i],IdentityMatrix[2]]+(1-p)*KroneckerProduct[IdentityMatrix[2],PauliMatrix[i]]
Z=MatrixExp[-\[Lambda]3*G[3,p]];
MaxEnt=Z/Tr[Z];
obsz=Tr[G[3,p] . Z]/Tr[Z]//FullSimplify
Solve[obsz==rz,\[Lambda]3,Reals, Assumptions->\[Lambda]3>0]



