(* ::Package:: *)

BeginPackage["MathematicaFMM`"];


Needs["CCompilerDriver`"];

(*Begin["`Private`"];*)
	
	sourceFile = "MathematicaFMM.cpp";
	compileOpts = {"-Wall","-Wno-unknown-pragmas","-pedantic","-fno-builtin","-std=c++17","-fopenmp"};
	sysCompileOpts = {"-O3","-fPIC"}; 
	libs = {"stdc++fs","gsl","gslcblas"}; 
	
	lib2DGrav=CreateLibrary[{sourceFile},"FMM2DGrav",
	"CompileOptions"-> compileOpts,
	"SystemCompileOptions"->sysCompileOpts,
	"Defines"->{"DIM"->"2","FIELD_TYPE"->"True"}, (* True for gravitational, false for Coulomb *)
	"Libraries"->libs,"ShellOutputFunction"->Print,"ShellCommandFunction"->Print];
	
	lib3DGrav=CreateLibrary[{sourceFile},"FMM3DGrav",
	"CompileOptions"-> compileOpts,
	"SystemCompileOptions"->sysCompileOpts,
	"Defines"->{"DIM"->"3","FIELD_TYPE"->"True"}, (* True for gravitational, false for Coulomb *)
	"Libraries"->libs,"ShellOutputFunction"->Print,"ShellCommandFunction"->Print];
	
	(*
	lib2DCoul=CreateLibrary[{sourceFile},"FMM2DCoul",
	"CompileOptions"\[Rule] compileOpts,
	"SystemCompileOptions"\[Rule]sysCompileOpts,
	"Defines"\[Rule]{"DIM"\[Rule]"2","FIELD_TYPE"\[Rule]"False"}, (* True for gravitational, false for Coulomb *)
	"Libraries"\[Rule]libs,"ShellOutputFunction"\[Rule]Print,"ShellCommandFunction"\[Rule]Print];
	
	lib3DCoul=CreateLibrary[{sourceFile},"FMM3DCoul",
	"CompileOptions"\[Rule] compileOpts,
	"SystemCompileOptions"\[Rule]sysCompileOpts,
	"Defines"\[Rule]{"DIM"\[Rule]"3","FIELD_TYPE"\[Rule]"False"}, (* True for gravitational, false for Coulomb *)
	"Libraries"\[Rule]libs,"ShellOutputFunction"\[Rule]Print,"ShellCommandFunction"\[Rule]Print];
	*)
	
	FMM2DGravParticleForces=LibraryFunctionLoad[lib2DGrav,"FMMParticleForces",{{Real,1}(*positions*),{Real,1}(*charges*),Integer (*items per leaf*),Real (*eps*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	FMM2DGravParticlePotentialEnergies=LibraryFunctionLoad[lib2DGrav,"FMMParticlePotentialEnergies",{{Real,1}(*positions*),{Real,1}(*charges*),Integer (*items per leaf*),Real (*eps*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	FMM2DGravEvaluatePotentials=LibraryFunctionLoad[lib2DGrav,"FMMEvaluatePotentials",{{Real,1}(*positions*),{Real,1}(*charges*),Integer (*items per leaf*),Real (*eps*),Real (*force_eps*),{Real,1}(*eval points*),"Boolean" (*tree to file?*)},{Real,1} (*output array of potentials*)];

	Direct2DGravParticleForces=LibraryFunctionLoad[lib2DGrav,"DirectParticleForces",{{Real,1}(*positions*),{Real,1}(*charges*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	Direct2DGravParticlePotentialEnergies=LibraryFunctionLoad[lib2DGrav,"DirectParticlePotentialEnergies",{{Real,1}(*positions*),{Real,1}(*charges*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	Direct2DGravEvaluatePotentials=LibraryFunctionLoad[lib2DGrav,"DirectEvaluatePotentials",{{Real,1}(*positions*),{Real,1}(*charges*),Real (*force_eps*),{Real,1}(*eval points*)},{Real,1} (*output array of potentials*)];

	FMM3DGravParticleForces=LibraryFunctionLoad[lib3DGrav,"FMMParticleForces",{{Real,1}(*positions*),{Real,1}(*charges*),Integer (*items per leaf*),Real (*eps*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	FMM3DGravParticlePotentialEnergies=LibraryFunctionLoad[lib3DGrav,"FMMParticlePotentialEnergies",{{Real,1}(*positions*),{Real,1}(*charges*),Integer (*items per leaf*),Real (*eps*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	FMM3DGravEvaluatePotentials=LibraryFunctionLoad[lib3DGrav,"FMMEvaluatePotentials",{{Real,1}(*positions*),{Real,1}(*charges*),Integer (*items per leaf*),Real (*eps*),Real (*force_eps*),{Real,1}(*eval points*),"Boolean" (*tree to file?*)},{Real,1} (*output array of potentials*)];

	Direct3DGravParticleForces=LibraryFunctionLoad[lib3DGrav,"DirectParticleForces",{{Real,1}(*positions*),{Real,1}(*charges*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	Direct3DGravParticlePotentialEnergies=LibraryFunctionLoad[lib3DGrav,"DirectParticlePotentialEnergies",{{Real,1}(*positions*),{Real,1}(*charges*),Real (*force_eps*)},{Real,1} (*output array of sources*)];
	Direct3DGravEvaluatePotentials=LibraryFunctionLoad[lib3DGrav,"DirectEvaluatePotentials",{{Real,1}(*positions*),{Real,1}(*charges*),Real (*force_eps*),{Real,1}(*eval points*)},{Real,1} (*output array of potentials*)];
	
(*End[];*)


FMM2DGravParticleForces::usage = "FMM2DGravParticleForces[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, items_per_leaf, acc_eps, force_smoothing_eps] computes gravitational interaction forces for a set of particles in 2D in parallel via the FMM.";
FMM2DGravParticlePotentialEnergies::usage = "FMM2DGravParticlePotentialEnergies[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, items_per_leaf, acc_eps, force_smoothing_eps] computes gravitational potential energies for a set of particles in 2D in parallel via the FMM";
FMM2DGravEvaluatePotentials::usage = "FMM2DGravEvaluatePotentials[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, items_per_leaf, acc_eps, force_smoothing_eps, eval_points = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, tree_to_file = True|False] computes gravitational potentials at a set of points by a set of charges in 2D in parallel via the FMM"; 

Direct2DGravParticleForces::usage = "Direct2DGravParticleForces[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, force_smoothing_eps] computes gravitational forces for a set of particles in 2D in parallel by the direct method";
Direct2DGravParticlePotentialEnergies::usage = "Direct2DGravParticlePotentialEnergies[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, force_smoothing_eps] computes gravitational potential energies for a set of particles in 2D in parallel by the direct method";
Direct2DGravEvaluatePotentials::usage = "Direct2DGravEvaluatePotentials[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, force_smoothing_eps, eval_points = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}] computes gravitational potentials at a set of points by a set of charges in 2D in parallel by the direct method";

FMM3DGravParticleForces::usage = "FMM3DGravParticleForces[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(z\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\),\!\(\*SubscriptBox[\(z\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, items_per_leaf, acc_eps, force_smoothing_eps] computes gravitational interaction forces for a set of particles in 3D in parallel via the FMM.";
FMM3DGravParticlePotentialEnergies::usage = "FMM3DGravParticlePotentialEnergies[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(z\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\),\!\(\*SubscriptBox[\(z\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, items_per_leaf, acc_eps, force_smoothing_eps] computes gravitational potential energies for a set of particles in 3D in parallel via the FMM";
FMM3DGravEvaluatePotentials::usage = "FMM3DGravEvaluatePotentials[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(z\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\),\!\(\*SubscriptBox[\(z\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, items_per_leaf, acc_eps, force_smoothing_eps, eval_points = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}, tree_to_file = True|False] computes gravitational potentials at a set of points by a set of charges in 3D in parallel via the FMM"; 

Direct3DGravParticleForces::usage = "Direct3DGravParticleForces[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(z\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\),\!\(\*SubscriptBox[\(z\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, force_smoothing_eps] computes gravitational forces for a set of particles in 3D in parallel by the direct method";
Direct3DGravParticlePotentialEnergies::usage = "Direct3DGravParticlePotentialEnergies[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(z\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\),\!\(\*SubscriptBox[\(z\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, force_smoothing_eps] computes gravitational potential energies for a set of particles in 3D in parallel by the direct method";
Direct3DGravEvaluatePotentials::usage = "Direct3DGravEvaluatePotentials[positions = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(z\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\),\!\(\*SubscriptBox[\(z\), \(2\)]\)...}, charges = {\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...}, force_smoothing_eps, eval_points = {\!\(\*SubscriptBox[\(x\), \(1\)]\),\!\(\*SubscriptBox[\(y\), \(1\)]\),\!\(\*SubscriptBox[\(x\), \(2\)]\),\!\(\*SubscriptBox[\(y\), \(2\)]\)...}] computes gravitational potentials at a set of points by a set of charges in 3D in parallel by the direct method";


EndPackage[]
