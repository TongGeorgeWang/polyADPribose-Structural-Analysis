seq = 'AUCUUUCUUUUCACCUCCUCCUUUCCCUUA'
fname = 'D_Seq_1.mat'
Make_ssRNA_mol_prms(seq,fname)
%%
seq = 'AAGAAUAAAAGAGAAGCCACCCCACCCAGA'
fname = 'D_Seq_2.mat'
Make_ssRNA_mol_prms(seq,fname)
%%
seq = 'AGUCCCAUAAGUCAUCCGUCCAGUUCCAGU'
fname = 'D_Seq_3.mat'
Make_ssRNA_mol_prms(seq,fname)
%%
load D_Seq_1.mat
mol_prms = iterative_refinement(mol_prms,'SAM23Mg',1);
load D_Seq_2.mat
mol_prms = iterative_refinement(mol_prms,'SAM23Mg',1);
%%
load D_Seq_3.mat
mol_prms = iterative_refinement(mol_prms,'SAM23Mg',1);