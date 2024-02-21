function [fnames,fullfnames] = export_pdb(mol_prms,mc_output,idx,postfix,dirname)

if nargin==4 || isempty(dirname)
    fpath = '';
else
    fpath = [dirname '\'];
    if ~isdir(dirname)
        mkdir(dirname);
    end
end

pdbfilefmtstr = ['%05d_' postfix '.pdb'];
fnames = cell(size(idx));
fullfnames = cell(size(idx));

num_mc_struct = size(mc_output.c_all,1);
num_extras = 0;
if (isfield(mol_prms, 'pool_extras'))
    num_extras = numel(mol_prms.pool_extras);
end

for k = 1:length(idx)
    fnames{k} = sprintf(pdbfilefmtstr,idx(k));
    fullfnames{k} = [fpath fnames{k}];
    export_1pdb(fullfnames{k},mol_prms,mc_output,idx(k),num_mc_struct,num_extras);
end

end

function export_1pdb(fname,mol_prms,mc_output,idx,num_mc_struct,num_extras)
% Number the structures in the pool as follows:
%   The structure corresponding to the nth row of mc_output shall be the
%   nth structure.
%   If there are additional structures in pool_extras, those shall be
%   numbered starting at (num_mc_struct + 1).

    if (idx <= num_mc_struct)
        if isfield(mol_prms, 'seq')
            if isfield(mol_prms, 'startIdx')
                atoms = clist2atoms(mc_output.c_all(idx,:),...
                    mol_prms.suites,mol_prms.template_atoms, mol_prms.seq,...
                    mol_prms.startIdx);
            else
                atoms = clist2atoms(mc_output.c_all(idx,:),...
                    mol_prms.suites,mol_prms.template_atoms, mol_prms.seq);
            end
        else
            atoms = clist2atoms(mc_output.c_all(idx,:),...
                mol_prms.suites,mol_prms.template_atoms);
        end
    elseif isfield(mol_prms, 'pool_extras') && idx-num_mc_struct <= num_extras
        atoms = mol_prms.pool_extras{idx-num_mc_struct};
    else
        error('No structure with given index %d found.', idx)
    end
    
    quick_pdb_writer(fname, atoms);
end