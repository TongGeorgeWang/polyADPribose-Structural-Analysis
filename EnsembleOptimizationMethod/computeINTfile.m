clear; close all 
%goal: running crysol and gajoe from Kush's PAR15 models

%% to do before running this script
% copy the pdb and rename the frames
%rename_frames.m
%%
%move pdbs from PDB folder to folder where processing is to be done

load('X.mat') 
prms_in.dirname = 'PARna'; 
prms_in.num_structures = 2066; 
prms_in.curvenames = 'PARna.dat'; 


%% from batch_run.m
dirname =           prms_in.dirname;
num_structures =    prms_in.num_structures;
num_mc_steps =      prms_in.num_mc_steps;
num_burnin =        prms_in.num_burnin;
postfix =           prms_in.postfix;
datadir =           prms_in.datadir;
curvenames =        prms_in.curvenames;
generations =       prms_in.generations;
%%
crysol_options = '/lm 15 /fb 18 /sm .3 /ns 61';
pdbdir = prms_in.dirname;
cmd_cleanup = {};

%%  
for i=0:prms_in.num_structures
    fnames = [num2str(i) '_X.pdb']; 
    [fprefix] = [num2str(i) '_X'];
    logfilename = [fprefix '00.log'];
    intfilename = [fprefix '00.int'];
    almfilename = [fprefix '00.alm'];

    cmdlist = {'echo off'};
    cmdlist = [cmdlist, cmd_cleanup];
    cmdlist = [cmdlist,'erase ' logfilename ' ' intfilename ' ' almfilename];
    cmdlist = [cmdlist,['crysol ' fnames ' ',crysol_options]];
    remote_dos(pdbdir,cmdlist,3);
    [rg,dmax] = read_crysol_logfile(pdbdir,logfilename);
    write_size_list(pdbdir,i,rg,dmax);
    [q,iq] = read_crysol_intfile(pdbdir,intfilename);
    write_intensity_file(pdbdir,postfix,i,q,iq);
    cmd_cleanup = {['erase ' logfilename ' ' intfilename ' ' almfilename ' ']};
 %   cmd_cleanup = {['erase ' logfilename ' ' intfilename ' ' almfilename ' ' fnames{1}]};

end
remote_dos(pdbdir,cmd_cleanup,0);


%%
%gajoe_output = gajoeit_inline(dirname,postfix,datadir,curvenames,generations);
%%