function mol_prms = iterative_refinement(mol_prms,prefix,niter)
    run_prms = struct(...
        'dirname','',...
        'num_structures',500,... % did 500 for SAM-II % should be 1000
        'num_mc_steps',   50,... % should be 50
        'num_burnin',     10,... % should be 10
        'postfix','X',...
        'datadir','saxs_data',...
        'curvenames',[prefix '.dat'],...
        'generations',50);
    t0 = tic;
    fprintf(1,'\nREFINEMENT, TIME IS: %s (%.2f hours elapsed)\n',datestr(now),toc(t0)/3600);
    
    i = 1;
    while isdir(sprintf('%s_%02d',prefix,i-1))
        i = i+1;
    end
    imin = i;
    
    for i=imin:(imin+niter-1)
        run_prms.dirname = sprintf('%s_%02d',prefix,i-1);
        [~, ~, torsion_analysis] = ...
            batch_run(mol_prms,run_prms);
        mol_prms.suite_weights = torsion_analysis.ws_new;
        fprintf(1,'\nREFINEMENT, TIME IS: %s (%.2f hours elapsed)\n',datestr(now),toc(t0)/3600);
    end
end