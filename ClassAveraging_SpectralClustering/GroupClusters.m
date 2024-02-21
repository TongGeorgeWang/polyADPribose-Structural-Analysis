
%% Copy all PDBs into subfolders in the working directory based on cluster number 
%   Only execute this part when you are satisfied with the clustering situation 
%
%   Requires you to input the 'clusters' variable [Nstructuresx2 matrix]
%   from 'RNA_classAverage.m' script and load it into the workspace before
%   running. 
%
%   GW - 2023 March  

%Nstructures = 75; %PAR22na
%Nstructures = 78; %PAR22naMg
%Nstructures = 297; %PAR15na
Nstructures = 582; %PAR15naMg

folderName = 'PAR15naMg';
basename = 'PAR15naMgPDB ('; % preceeding the number
postname = ').pdb'; % following the number


%% Do the file copying 
nClusters = numel(unique(clusters(:,2)));
clusterPDBsAll = clusters(:,1);
if ~exist([folderName,'/PDBs_SpectralClustered'],'dir')
    mkdir([folderName,'/PDBs_SpectralClustered'])
end

if ~exist([folderName,'/PDBs_SpectralClustered/Cluster1'],'dir')
    for i = 0:(nClusters-1)
        mkdir([folderName,'/PDBs_SpectralClustered/Cluster',num2str(i)])
        
        clusterPDBindices = find(clusters(:,2)==i);
        clusterPDBnumbers = clusterPDBsAll(clusterPDBindices);
        
        for j = 1:numel(clusterPDBnumbers)
            copyfile([folderName,'/',basename,num2str(clusterPDBnumbers(j)),postname],...
                     [folderName,'/PDBs_SpectralClustered/Cluster',num2str(i)])
        end
        
    end
end











