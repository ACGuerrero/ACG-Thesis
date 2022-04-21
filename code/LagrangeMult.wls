#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["ThesisTools`"]
Needs["CoolTools`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];


(*Reconstrucci\[OAcute]n del MaxEnt*)
$Assumptions=Element[{p,rz,\[Lambda],Z,t},Reals] && rz\[Element]Reals && 0<=rz<=1&&0<=p<=1&&0<=t<=1;
A=Subscript[\[Lambda],3]*GObsMaxEnt[p,3];
ExpMat=MatrixExp[A];
Z=Tr[ExpMat];
1/Z//FullSimplify;
MaxEnt=ExpMat/Z;
rz=ZCoordFromLagrangeMult[\[Lambda],p];
plot3=Plot[Evaluate[Table[rz/.{p->i},{i,0,0.5,0.1}]],{\[Lambda],-4,0},AxesLabel->{Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},PlotLegends->{Table["p="<>ToString[i],{i,0,0.5,0.1}]}];
trb=(IdentityMatrix[2]+Tanh[p*\[Lambda]]*PauliMatrix[3])/2;
tra=(IdentityMatrix[2]+Tanh[\[Lambda]*(1-p)]*PauliMatrix[3])/2;
eff=p*trb+(1-p)*tra;


plot4=Plot[
	Evaluate[
		Table[ZCoordFromLagrangeMult[\[Lambda],p]/.{p->i},{i,0.0,0.5,0.25}]
		],
		{\[Lambda],0,8},
	PlotStyle->{
		{Thickness[0.003],Black},
		{Thickness[0.003],Dashing[Large],Red},
		{Thickness[0.003],DotDashed,Blue}},
	AxesStyle->Black,PlotRange->{{0,8},{0,1}},AxesLabel->{Style[\[Lambda],FontSize->15],Style[Subscript[r,\[Rho]],FontSize->15]},
	PlotLegends->Placed[{Table["p="<>ToString[i],{i,0.0,0.5,0.25}]},{0.8,0.5}]]


	Export["../figures/r(lambda).png",plot4]


Tr[SubMaxEnt1=MatrixExp[-\[Lambda]*p*PauliMatrix[3]]]//FullSimplify


$Assumptions=Element[{p,rx,ry,rz,\[Lambda]1,\[Lambda]2,\[Lambda]3,Z,t},Reals] && rz\[Element]Reals && 0<=rz<=1&&0<=p<=1&&0<=t<=1;
A=-(\[Lambda]1*GObsMaxEnt[p,1]+\[Lambda]2*GObsMaxEnt[p,2]+\[Lambda]3*GObsMaxEnt[p,3]);
ExpMat=MatrixExp[A];
Z=Tr[ExpMat];
MaxEnt=(ExpMat/Z)/.{Sqrt[\[Lambda]1^2+\[Lambda]2^2+\[Lambda]3^2]->\[Lambda]}//FullSimplify;


rx=Tr[GObsMaxEnt[p,1] . MaxEnt]/.{Sqrt[\[Lambda]1^2+\[Lambda]2^2+\[Lambda]3^2]->\[Lambda]}//FullSimplify;
ry=Tr[GObsMaxEnt[p,2] . MaxEnt]/.{Sqrt[\[Lambda]1^2+\[Lambda]2^2+\[Lambda]3^2]->\[Lambda]}//FullSimplify;
rz=Tr[GObsMaxEnt[p,3] . MaxEnt]/.{Sqrt[\[Lambda]1^2+\[Lambda]2^2+\[Lambda]3^2]->\[Lambda]}//FullSimplify;


(* ::Input:: *)
(*\[DoubleDot]*)


(*
ContourPlot[expr, {p,0,1},{\[Lambda]3,0,-1}, Contours -> Table[i,{i,0.0,0.4,0.02}], 
   ContourStyle -> Black, ContourShading -> None,ContourLabels->True, 
   FrameLabel -> {"p", "\[Lambda]3"}]
plot1=Plot3D[expr,{p,0,1},{\[Lambda]3,-8,8},AxesLabel->{Style[p,FontSize->18],Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},ColorFunction->GrayLevel,ViewPoint -> {3, 2.2, 2}];
Export["../figures/LagrangeMult_lambda-8to8.png",plot1]
plot2=Plot3D[expr,{p,0,1},{\[Lambda]3,-4,0},AxesLabel->{Style[p,FontSize->18],Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},ColorFunction->GrayLevel,ViewPoint -> {3, 2.2, 2}];
Export["../figures/LagrangeMult_lambda-4to0.png",plot2]
plot3=Plot[Evaluate[Table[expr/.{p->i},{i,0,0.5,0.1}]],{\[Lambda]3,-4,0},AxesLabel->{Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},PlotLegends->{Table["p="<>ToString[i],{i,0,0.5,0.1}]}];
Export["../figures/rz_has_inverse_lambda-4to0.png",plot3]

plot3=Plot[Evaluate[Table[expr/.{p->i},{i,0.5,1,0.1}]],{\[Lambda]3,-4,4},AxesLabel->{Style[Subscript[\[Lambda],3],FontSize->18],Style[Subscript[r,z],FontSize->18]},PlotLegends->{Table["q="<>ToString[i],{i,0,0.5,0.1}]}];
Export["../figures/rz_has_inverse_lambda-4to4.png",plot3];  
*)
data=Table[{rz[l,0.8],l},{l,-4,0,0.05}];
data//MatrixForm
ListLinePlot[data]
nearpos[haystack_,value_]:= With[{ f = Nearest[haystack -> Range@Length@haystack]},f[value, 1]];
nearpos[Transpose[data][[1]],0.8]
Transpose[data][[1,nearpos[Transpose[data][[1]],0.8]]]
Transpose[data][[2,nearpos[Transpose[data][[1]],0.8]]]


ContourPlot[expr, {p,0,1},{\[Lambda]3,0,-1}, Contours -> Table[i,{i,0.0,0.4,0.02}], 
   ContourStyle -> Black, ContourShading -> None,ContourLabels->True, 
   FrameLabel -> {"p", "\[Lambda]3"}]
