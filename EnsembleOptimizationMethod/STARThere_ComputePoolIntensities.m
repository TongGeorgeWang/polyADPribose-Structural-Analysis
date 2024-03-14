clear; close all 


%% Generate CRYSOL-derived theoretical scattering intensities for a pdb structural ensemble
%   Input: pdb ensemble in a subfolder
%       The PDBs should be named in ascending order from 0 (ie '0.pdb',
%       '1.pdb', '2.pdb', ...)
%   Output: .int file containing CRYSOL-computed I(q) arrays of every structure of the ensemble; 
%       this is a downstream input for GAJOE
% 
%   CRYSOL is a program from ATSAS. The required executable and its
%   dependencies should be included in the pdb subfolder. This program has
%   been properly tested with older versions of ATSAS (up to ATSAS2.8.4).
%
%   originally written by either Derrick, Alex P, and/or Suzette; modified
%   and commented by George W heavily because it was very confusing
%

%% Specify these parameters

dirname = 'PDBsKanja_12'; % Name of the subfolder with the pdb ensembles
fileTypes = '.pdb'; % ie '.pdb', '.cif'
num_structures = 406; % Number of pdbs in the pdb ensemble subfolder 
crysol_options = ' /lm 15 /fb 18 /sm .3 /ns 61'; % see ATSAS documentation


%% Execute CRYSOL and save theoretical I(q) for each structure in the pool 

pdbdir = dirname;
cmd_cleanup = {};

for i=0:num_structures-1
    fnames = [num2str(i) fileTypes]; 
    logfilename = [num2str(i) '00.log'];
    intfilename = [num2str(i) '00.int'];
    almfilename = [num2str(i) '00.alm'];
    
    % CRYSOL will be executed in a separate DOS (command line) batch; set up the batch here
    cmdlist = {'echo off'};
    cmdlist = [cmdlist, cmd_cleanup];
    cmdlist = [cmdlist,['crysol ' fnames crysol_options]];
    
    % Run the batch commands that were set up above (ie, run CRYSOL)
    remote_dos(pdbdir,cmdlist,3)
    
    % Write the size_list file with Rg and Dmax of each structure and a
    % .int file with the I(q) of each structure
    [rg,dmax] = read_crysol_logfile(pdbdir,logfilename);
    write_size_list(pdbdir,i,rg,dmax);
    [q,iq] = read_crysol_intfile(pdbdir,intfilename);
    write_intensity_file(pdbdir,i,q,iq);

    % Extra command to delete the crysol outputs, to keep the directory tidier
    cmd_cleanup = {['erase ' logfilename ' ' intfilename ' ' almfilename ' ']};

end

remote_dos(pdbdir,cmd_cleanup,0);


