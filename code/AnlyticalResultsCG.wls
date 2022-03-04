#!/usr/bin/env wolframscript
(* ::Package:: *)

Needs["ThesisTools`"]
Needs["CoolTools`"]
Needs["Quantum`"]
SetDirectory["/home/acastillo/Documents/tesis-adan/code"];


??CGKrausOp


krausop=CGKrausOp[p];
unitary=KroneckerProduct[krausop[[1]],{{1},{0},{0},{0}}]+KroneckerProduct[krausop[[2]],{{0},{1},{0},{0}}]+KroneckerProduct[krausop[[3]],{{0},{0},{1},{0}}]+KroneckerProduct[krausop[[4]],{{0},{0},{0},{1}}]


unitary


Clear[CGUnitary]
CGUnitary[p_]:={
{Sqrt[p],0,0,0,-Sqrt[1-p],0,0,0},
{0,Sqrt[p],0,0,0,-Sqrt[1-p],0,0},
{Sqrt[1-p],0,0,0,Sqrt[p],0,0,0},
{0,0,Sqrt[1-p],0,0,0,Sqrt[p],0},
{0,0,Sqrt[p],0,0,0,-Sqrt[1-p],0},
{0,0,0,Sqrt[p],0,0,0,-Sqrt[1-p]},
{0,Sqrt[1-p],0,0,0,Sqrt[p],0,0},
{0,0,0,Sqrt[1-p],0,0,0,Sqrt[p]}}


CGUnitary[p]//MatrixForm


CGUnitary[0.5] . ConjugateTranspose[CGUnitary[0.5]]//MatrixForm


p=0.5;
u=CGUnitary[p];
state=ketsToDensity[{RandomState[4]}][[1]];
bigstate=KroneckerProduct[state,{{1,0},{0,0}}];
fidelity[CGKraus[state,p],PartialTrace[u . bigstate . Dagger[u],1]]
