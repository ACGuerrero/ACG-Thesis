#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["ThesisTools`"]
Needs["CoolTools`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];


(*Reconstrucci\[OAcute]n del MaxEnt*)
$Assumptions=Element[{p,rz,\[Lambda],Z,t},Reals] && rz\[Element]Reals && 0<=rz<=1&&0<=p<=1&&0<=t<=1;
A=-\[Lambda]*GObsMaxEnt[p,3];
ExpMat=MatrixExp[A];
Z=Tr[ExpMat];
1/Z//FullSimplify
MaxEnt=ExpMat/Z;

(*Obtenci\[OAcute]n del factor de contracci\[OAcute]n*)
rhoz0=(Refine[CGKraus[MaxEnt,p],{Element[p, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand);
Print["The bloch vector of the last operator is: "]
rhoz0bv=densityMatrixToPoint[{rhoz0},gellMannBasis[1]][[1]]//FullSimplify
rhozt=(Refine[CGKraus[swapGate . MaxEnt . swapGate,p],{Element[p, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand);
Print["The bloch vector of the last operator is: ", rhoztbv]
rhoztbv=densityMatrixToPoint[{rhozt},gellMannBasis[1]][[1]]//FullSimplify
Print["The factor separating this states is: " (contractionfactor=rhoztbv[[3]]/rhoz0bv[[3]])//FullSimplify];
plot1=Plot3D[contractionfactor,{p,0,1},{\[Lambda],-4,4},AxesLabel->{Style[p,FontSize->18],Style[\[Lambda],FontSize->18],Style[a,FontSize->18]},ColorFunction->GrayLevel,ViewPoint -> {3, 2.2, 2}]
Export["../figures/ContractionFactorSWAP_3D_lambda-4to4.png",plot1]
plot3=Plot[Evaluate[Table[contractionfactor/.{p->i},{i,0,0.5,0.1}]],{\[Lambda],-8,0},AxesLabel->{Style[\[Lambda],FontSize->18],Style[a,FontSize->18]},PlotLegends->{Table["p="<>ToString[i],{i,0,0.5,0.1}]}]
Export["../figures/ContractionFactorSWAP_2D_lambda-8to0.png",plot3]


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


Print["To develop the same results for a power of the unitary we need a analytical expression."]
Print["The eigensystem of the swap unitary:"]
Eigensystem[swapGate]





Print["To develop the same results for a power of the unitary we need a analytical expression."]
Print["The eigensystem of the swap unitary:"]
P={{1,0,0,0},{0,1,1,0},{0,-1,1,0},{0,0,0,1}}
P . swapGate . P


SWAP[t_]:={{1,0,0,0},{0,(1+Exp[I*Pi*t])/2,(1-Exp[I*Pi*t])/2,0},{0,(1-Exp[I*Pi*t])/2,(1+Exp[I*Pi*t])/2,0},{0,0,0,1}}


SwapedMaxEntT=Refine[SWAP[t] . MaxEnt . ConjugateTranspose[SWAP[t]],
{Element[t, Reals]}];
CoarseEvT=Refine[CGKraus[SwapedMaxEntT,p],
{Element[p, Reals]}]


rhoztt=(Refine[CGKraus[SWAP[t] . MaxEnt . ConjugateTranspose[SWAP[t]],p],{Element[p, Reals],Element[t, Reals], Element[\[Lambda], Reals], 0<=p<=1}]//Expand);
Print["The bloch vector of the last operator is: ", 
rhozttbv=densityMatrixToPoint[{rhoztt},gellMannBasis[1]][[1]]//FullSimplify]

Print["The factor separating this states is: ", (contractionfactor=rhozttbv[[3]]/rhoz0bv[[3]])//FullSimplify];


(* ::Input:: *)
(*\.b4*)


SWAPContractionFactor[1,0.3,1]
