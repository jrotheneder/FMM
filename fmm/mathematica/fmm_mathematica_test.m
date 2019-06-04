(* ::Package:: *)

(* ::Input::Initialization:: *)
dim = 2; 
nc = 20000; 
charges = RandomReal[{-1,1},nc];

centers = {{-5,0},{5,0}};
positions = Join[RandomVariate[MultinormalDistribution[centers[[1]],IdentityMatrix[2]],Floor[nc/2]],RandomVariate[MultinormalDistribution[centers[[2]],IdentityMatrix[2]],Ceiling[nc/2]]];
sources = MapThread[Join[#1,{#2}]&,{positions,charges}];

Needs["CCompilerDriver`"];
lib=CreateLibrary[{"./mathematica.cpp"},"FMM2D_Forces","CompileOptions"->{"-Wall","-Wno-unknown-pragmas","-pedantic","-fno-builtin","-std=c++17","-fopenmp","-lstdc++fs"},"SystemCompileOptions"->{"-O3", "-lgsl", "-lcblas","-fPIC"},"ShellOutputFunction"->Print,"ShellCommandFunction"->Print];
FMM2D\[LetterSpace]Forces=LibraryFunctionLoad[lib,"FMM2D_Forces",{{Real,1}(* positions*),{Real,1}(* charges *),Integer (* items per leaf *), Real (* eps*)},{Real,1}];

Print[AbsoluteTiming[FMM2D\[LetterSpace]Forces[Flatten[positions],charges,5,0.01]][[1]]]



