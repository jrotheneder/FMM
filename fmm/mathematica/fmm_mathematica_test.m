dim = 2; 
fieldType = "True" (* grav *); 
nc = 200;
charges = RandomReal[{-1,1},nc];

centers = {{-5,0},{5,0}};
positions = Join[RandomVariate[MultinormalDistribution[centers[[1]],IdentityMatrix[2]],Floor[nc/2]],RandomVariate[MultinormalDistribution[centers[[2]],IdentityMatrix[2]],Ceiling[nc/2]]];
sources = MapThread[Join[#1,{#2}]&,{positions,charges}];

Needs["MathematicaFMM`"];

FMMForces =  FMM2DGravParticleForces;
FMMPotentials =  FMM2DGravParticlePotentialEnergies;
FMMPointEvalPotentials =  FMM2DGravEvaluatePotentials;

DirectPotentials =  Direct2DGravParticlePotentialEnergies; 
DirectForces =  Direct2DGravParticleForces;
DirectPointEvalPotentials = Direct2DGravEvaluatePotentials;

{tf, forces} = AbsoluteTiming[
    FMMParticleForces[Flatten[positions],charges,20,0.01,0]];
{tp, potentials} = AbsoluteTiming[
    FMMParticlePotentials[Flatten[positions],charges,20,0.01]];
{tpd, refpotentials} = AbsoluteTiming[DirectPotentials[Flatten[positions],charges]];

Print[tp]
Print[tpd]
Print[tf]
Print[potentials]
Print[tpx, ", max rel dev is ", Max[Abs[potentials-refpotentials]]]



