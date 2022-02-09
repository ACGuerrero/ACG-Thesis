#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["ThesisTools`"]
Clear[p]
$Assumptions=Element[{\[Lambda]3,p,rz},Reals] && 0<=p<=1 && rz\[Element]Reals && 0<=rz<=1;
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];
(*Observable fino*)
G[i_,p_]:=p*KroneckerProduct[PauliMatrix[i],IdentityMatrix[2]]+(1-p)*KroneckerProduct[IdentityMatrix[2],PauliMatrix[i]]

(*Reconstrucci\[OAcute]n del MaxEnt*)
A=-\[Lambda]3*G[3,p];
ExpMat=MatrixExp[A];
Z=Tr[ExpMat];
MaxEnt=ExpMat/Z;
(*
obsz=Tr[G[3,p] . Z]/Tr[Z]//FullSimplify
Solve[obsz==rz,\[Lambda]3,Reals, Assumptions->\[Lambda]3>0]
*)
(expr=Tr[G[3,p] . MaxEnt]//FullSimplify);
plot1=Plot3D[expr,{p,0,1},{\[Lambda]3,-8,8},AxesLabel->{Style[p,FontSize->18],Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},ColorFunction->GrayLevel,ViewPoint -> {3, 2.2, 2}];
Export["../figures/LagrangeMult_lambda-8to8.png",plot1]
plot2=Plot3D[expr,{p,0,1},{\[Lambda]3,-4,0},AxesLabel->{Style[p,FontSize->18],Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},ColorFunction->GrayLevel,ViewPoint -> {3, 2.2, 2}];
Export["../figures/LagrangeMult_lambda-4to0.png",plot2]
plot3=Plot[Evaluate[Table[expr/.{p->i},{i,0,0.5,0.1}]],{\[Lambda]3,-4,0},AxesLabel->{Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},PlotLegends->{Table["p="<>ToString[i],{i,0,0.5,0.1}]}];
Export["../figures/rz_has_inverse_lambda-4to0.png",plot3]
plot3=Plot[Evaluate[Table[expr/.{p->i},{i,0,0.5,0.1}]],{\[Lambda]3,-4,4},AxesLabel->{Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},PlotLegends->{Table["p="<>ToString[i],{i,0,0.5,0.1}]}];
Export["../figures/rz_has_inverse_lambda-4to4.png",plot3];  


ContourPlot[expr, {p,0,1},{\[Lambda]3,0,-1}, Contours -> Table[i,{i,0.0,0.4,0.02}], 
   ContourStyle -> Black, ContourShading -> None,ContourLabels->True, 
   FrameLabel -> {"p", "\[Lambda]3"}]


ContourPlot[expr, {p,0,1},{\[Lambda]3,0,1}, Contours -> -Table[i,{i,0.0,0.4,0.02}], 
   ContourStyle -> Black, ContourShading -> None,ContourLabels->True, 
   FrameLabel -> {"p", "\[Lambda]3"}]
