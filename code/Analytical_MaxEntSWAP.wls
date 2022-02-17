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
MaxEnt=ExpMat/Z;
rhoz0=(Refine[CGKraus[MaxEnt,p],{Element[p, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand);
rhozt=(Refine[CGKraus[swapGate . MaxEnt . swapGate,p],{Element[p, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand);
Print[Subscript[\[Psi],"max"],"=", (MaxEnt//FullSimplify)//MatrixForm]
Print["Notice that:"]
Print[Subscript[\[Psi],"max"],"=", (partialTraceB[MaxEnt]//FullSimplify)//MatrixForm,"xox",(partialTraceA[MaxEnt]//FullSimplify)//MatrixForm]

Print[Subscript[\[Rho],z],"(0)=",(rhoz0//FullSimplify)//MatrixForm]
Print[Subscript[\[Rho],z],"(t)=", (rhozt//FullSimplify)//MatrixForm]
Print["For the case p=0.5, the following results apply"]
Print[Subscript[\[Psi],"max"],"=", (MaxEnt/.p->1/2)//MatrixForm]
Print["Notice that:"]
Print[Subscript[\[Psi],"max"],"=", partialTraceB[(MaxEnt/.p->1/2)]//MatrixForm,"xox",partialTraceA[(MaxEnt/.p->1/2)]//MatrixForm]
Print["So, not only the state is separable but it is of the form rho xox rho."]
Print[Subscript[\[Rho],z],"(0)=",(rhoz0/.p->1/2)//MatrixForm]
Print[Subscript[\[Rho],z],"(t)=", (rhozt/.p->1/2)//MatrixForm]
Print["This might be rewritten using the expresion of the lagrange multiplier"]
lambda05=Log[(1-Subscript[r,z])/(1+Subscript[r,z])];
Print[Subscript[\[Psi],"max"],"=", (MaxEnt//.{\[Lambda]->lambda05,p->1/2}//FullSimplify)//MatrixForm]
Print[Subscript[\[Rho],z],"(0)=",(rhoz0//.{\[Lambda]->lambda05,p->1/2}//FullSimplify)//MatrixForm]
Print[Subscript[\[Rho],z],"(t)=", (rhozt//.{\[Lambda]->lambda05,p->1/2}//FullSimplify)//MatrixForm]



