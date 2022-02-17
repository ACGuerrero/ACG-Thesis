#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["ThesisTools`"]
Needs["CoolTools`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];


(*Reconstrucci\[OAcute]n del MaxEnt*)
$Assumptions=Element[{p,rz,\[Lambda],Z},Reals] && rz\[Element]Reals && 0<=rz<=1&&0<=p<=1;
A=-\[Lambda]*GObsMaxEnt[p,3];
ExpMat=MatrixExp[A];
Z=Tr[ExpMat];
Clear[Z]
MaxEnt=ExpMat/Z;
Print["In terms of the Lagrange multiplier, the initial effective state is: ", (Refine[CGKraus[MaxEnt,p],{Element[p, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand)//MatrixForm]
Print["In terms of the Lagrange multiplier, the evolved effective state is: ", (Refine[CGKraus[swapGate . MaxEnt . swapGate,p],{Element[p, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand)//MatrixForm]



