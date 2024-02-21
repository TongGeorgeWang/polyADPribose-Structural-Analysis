function [mc_output, gajoe_output, torsion_analysis] = ...
    batch_run(mol_prms,prms_in)

dirname =           prms_in.dirname;
num_structures =    prms_in.num_structures;
num_mc_steps =      prms_in.num_mc_steps;
num_burnin =        prms_in.num_burnin;
postfix =           prms_in.postfix;
datadir =           prms_in.datadir;
curvenames =        prms_in.curvenames;
generations =       prms_in.generations;

if isdir(dirname)
    error('ABORT: directory already exists.');
else
    mkdir(dirname);
end

matfilename = [dirname '\' postfix '.mat'];

save(matfilename,'mol_prms','prms_in');

mc_output = mc_pool(mol_prms,num_structures,num_mc_steps,num_burnin);

save(matfilename,'mc_output','-append');

process_pool_inline(mol_prms,mc_output,[],postfix,dirname);

gajoe_output = gajoeit_inline(dirname,postfix,datadir,curvenames,generations);

save(matfilename,'gajoe_output','-append');

torsion_analysis = reweight_torsions(mol_prms,mc_output,gajoe_output);

save(matfilename,'torsion_analysis','-append');

end