(* ::Package:: *)

(* ::Input::Initialization:: *)
dim = 2; 
fieldType = "True" (* grav *); 
nc = 1000; 
charges = RandomReal[{-1,1},nc];


centers = {{-5,0},{5,0}};
positions = Join[RandomVariate[MultinormalDistribution[centers[[1]],IdentityMatrix[2]],Floor[nc/2]],RandomVariate[MultinormalDistribution[centers[[2]],IdentityMatrix[2]],Ceiling[nc/2]]];
sources = MapThread[Join[#1,{#2}]&,{positions,charges}];

Needs["CCompilerDriver`"];
lib=CreateLibrary[{"./mathematica.cpp"},"FMM_Forces",
    "CompileOptions"->{"-Wall","-Wno-unknown-pragmas","-pedantic","-fno-builtin",
        "-std=c++17","-fopenmp"},
    "SystemCompileOptions"->{"-O3", "-fPIC"},
    "Defines" -> {"DIM" -> ToString[dim], "FIELD_TYPE" -> fieldType},
    "Libraries" -> {"stdc++fs", "gsl", "gslcblas"},

    "ShellOutputFunction"->Print,"ShellCommandFunction"->Print
];
FMMForces = LibraryFunctionLoad[lib,"FMM_Forces",
    {{Real,1}(* positions*),
    {Real,1}(* charges *),
    Integer (* items per leaf *), 
    Real (* eps*)},{Real,1}
];

{t, forces} = AbsoluteTiming[FMMForces[Flatten[positions],charges,5,0.01]];
Print[t]



