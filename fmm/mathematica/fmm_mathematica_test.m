dim = 2; 
fieldType = "True" (* grav *); 
nc = 1000; 
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
    FMMParticlePotentials[Flatten[positions],charges,20,0.01,0]];
{tpd, refpotentials} = AbsoluteTiming[DirectPotentials[Flatten[positions],charges,0]];
{tpx, evalpots} = AbsoluteTiming[
    FMMPointEvalPotentials[Flatten[positions],charges, 50, 0.001, 0,
    Flatten[positions], False]
];

Print[tp]
Print[tpd]
Print[tf]
Print[tpx, ", mean dev is ", Mean[refpotentials - charges * evalpots]]



