
%% Cluster EOM-selected structures using the hierarchical method 
%   GW August 2022
%   Follows the steps described in: https://www.mathworks.com/help/stats/hierarchical-clustering.html
%
%   Feb 2023 - made plots publication friendly and also fixed that
%   bug that caused dendrogram colouring to be incorrect
%

clear; close all

%% User entered parameters for where the GAJOE output is
dir = 'ExampleInputAndOutputData_PAR15na';
gajoe_output = 'best_curve001.txt';
suffix = '_X.pdb';
%xlim_Dmax = [35 130]; ylim_Rg = [10 40]; xlim_R = [0 120]; % for PAR15
xlim_Dmax = [50 150]; ylim_Rg = [15 50]; xlim_R = [0 140]; % for PAR22
fontSize = 20;
axisFontSize = 22;

% Play around with the clustering parameters; the number of clustercolors
% should agree with the number of clusters resulting from the distance cutoff
clusterDistanceCutoff = 300;
% Dendrogram and scatterplot may plot colours in different orders; may need
% to declare different orders for the clusters to match

% adjust colours - will need to play around with them
clusterColors_scatter = [1 .5 0; 0 1 0; 0 0 1; 1 0 0; 0 1 1]; % PAR15na
%clusterColors_scatter = [1 0 0; 0 0 1; 0 1 0; 1 .5 0; 0 1 1]; % PAR15naMg
%clusterColors_scatter = [1 0 0; 0 0 1; 0 1 1; 1 .5 0; 0 1 0]; % PAR22na
%clusterColors_scatter = [0 0 1; 0 1 1; 1 0 0; 1 .5 0; 0 1 0]; % PAR22naMg
clusterColors_dend = [0 1 1; 0 0 1; 1 .5 0; 1 0 0; 0 1 0];


% For reference, [1 0 0]=red; [0 1 0]=green; [0 0 1]=blue; [0 1 1]=cyan; [1 0 1]=magenta;
% -- -- -- -- -- -- -- -- -- -- -- --


%% Import variables 
% R (Ree), Rg, Dmax, compute E (fret eff) from R
load([dir,'/Ree_values.mat'])
%E = fret_eff(R,60)';
Sizelist = importdata([dir,'/Size_list.txt', ' ']);
Sizelist(:,3) = R;
pool.rg = Sizelist(:,2);

% GAJOE results
best_curve = importdata([dir,'/',gajoe_output], ' ',2);
bestcurve = best_curve.data;  
gaj = length(bestcurve(:,1));


%% Delineate GAJOE structures
Dmax_gaj = []; Rg_gaj = []; R_gaj = []; %E_gaj = [];
for n = 1:gaj
    D_Rep = repmat(Sizelist(bestcurve(n,1),3),bestcurve(n,2),1); 
    Rg_Rep = repmat(Sizelist(bestcurve(n,1),2),bestcurve(n,2),1);
    R_Rep = repmat(R(bestcurve(n,1)),bestcurve(n,2),1);
    %E_Rep = repmat(E(bestcurve(n,1)),bestcurve(n,2),1);
    Dmax_gaj = [Dmax_gaj; D_Rep]; Rg_gaj = [Rg_gaj; Rg_Rep];
    R_gaj = [R_gaj; R_Rep];% E_gaj = [E_gaj; E_Rep];
    clear D_Rep Rg_Rep R_rep %E_Rep
end
gajoe.dmax = Dmax_gaj; gajoe.rg = Rg_gaj; % This is the structural pool that we wish to cluster 
gajoe.r = R_gaj;


%% Plot all EOM-selected points in 2D Euclidean space 
figure('Name',dir); hold on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',axisFontSize)
grid on; box on
xlabel('$R_{EE}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex') 
ylabel('$R_{g}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex')
title('All EOM structures','FontSize',fontSize)
xlim(xlim_R); ylim(ylim_Rg)

scatter(gajoe.r, gajoe.rg, 20,'k','filled')
hold off


%% Perform hierarchical clustering using Ward's method

X = [gajoe.r, gajoe.rg];

Y = pdist(X); % Compute paired distances between all points
Ysquare = squareform(Y); 

Z = linkage(Y,'ward');
I = inconsistent(Z); % Inconsistency matrix (quality control metric)

T = cluster(Z,'cutoff',clusterDistanceCutoff,'Criterion','distance');


%% Compute centroid of each cluster and take structure closest to centroid as the representative structure of that cluster 

clusteredData = sortrows([X,T],3);
nClusters = numel(unique(clusteredData(:,3))); 
for i = 1:nClusters
    clusterIdx = find(clusteredData(:,3) == i);
    dataInCurrentCluster = [clusteredData(clusterIdx,1), clusteredData(clusterIdx,2)];     
    clusterCentroid(i,:) = [mean(dataInCurrentCluster(:,1)) , mean(dataInCurrentCluster(:,2))];
    
    nStructuresInCluster = numel(dataInCurrentCluster(:,1));
    
    meanR = mean(dataInCurrentCluster(:,1));
    stDevR = std(dataInCurrentCluster(:,1));
    
    meanRg = mean(dataInCurrentCluster(:,2));
    stDevRg = std(dataInCurrentCluster(:,2));
    
    disp(['Cluster ',num2str(i),' mean R = ',num2str(meanR),'; stDev = ',num2str(stDevR),' Angstroms; # of structures in cluster: ',num2str(nStructuresInCluster)])
    disp(['Cluster ',num2str(i),' mean Rg = ',num2str(meanRg),'; stDev = ',num2str(stDevRg),' Angstroms; # of structures in cluster: ',num2str(nStructuresInCluster)])
    
    representativeStructureIndex(i,:) = dsearchn(dataInCurrentCluster, clusterCentroid(i,:));
    representativeStructure(i,:) = dataInCurrentCluster(representativeStructureIndex(i,:),:);
    
    dataInAllClusters{i} = unique(dataInCurrentCluster,'rows'); 
end


%% Make dendrogram and scatter plots

figure('Name',dir)
set(gcf,'color','w')
box on
set(gca,'LineWidth',2)
set(gca,'FontSize',axisFontSize)
H = dendrogram(Z,0,'ColorThreshold',clusterDistanceCutoff);
xlabel('Node','FontSize',fontSize) 
ylabel('Distance','FontSize',fontSize)
set(gca,'FontName', 'Arial')
title('Dendrogram of linkages')

%// Controlling the colours of the dendrogram
lineColours = cell2mat(get(H,'Color'));
colourList = unique(lineColours, 'rows');
newLineColours = lineColours;

%// Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines             
for colour = 2:size(colourList,1)
    %// Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    %// Replace the colour for those lines
    newLineColours(idx, :) = repmat(clusterColors_dend(colour-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for line = 1:size(H,1)
    set(H(line), 'Color', newLineColours(line,:));
end


figure('Name',dir); hold on
set(gcf, 'color', 'w')
grid on; box on
set(gca,'LineWidth',2)
set(gca,'FontSize',axisFontSize)
xlabel('$R_{EE}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex','FontName', 'Arial') 
ylabel('$R_{g}$ ($\AA$)','FontSize',fontSize,'Interpreter','latex','FontName', 'Arial')
%title('All EOM structures, clustered','FontSize',fontSize)
xlim(xlim_R); ylim(ylim_Rg)

colororder()
scatter(clusteredData(:,1), clusteredData(:,2), 20, clusteredData(:,3),'filled')
scatter(representativeStructure(:,1),representativeStructure(:,2),200,'k*') % highlight representative structures
colormap(clusterColors_scatter)
hold off

disp(' ')


%% Get and print representative structure PDB numbers

representativeStructure = sort(representativeStructure,1); % sort so PDB #s will be in order
for i = 1:numel(representativeStructure(:,1))
    representativeStructureNumbers(i) = find(...
        Sizelist(:,3)==representativeStructure(i,1) & ...
        Sizelist(:,2)==representativeStructure(i,2) );
end

disp(['Representative structures: ',num2str(representativeStructureNumbers)])
disp(' ')


%% Print list of PDB numbers of all structures in each cluster 
for j = 1:nClusters
    clusterPDBindices = ismember(Sizelist,dataInAllClusters{j});
    clusterPDBnumbers{j} = find(clusterPDBindices(:,2)==1 & clusterPDBindices(:,3)==1);
    
    suffixArray = repmat(suffix,numel(clusterPDBnumbers{j}),1);
    writecell( strtrim( cellstr( num2str(clusterPDBnumbers{j}, ['%d',suffix] ) ) ) ,...
        [dir,'/cluster',num2str(j),'_PDBlist.txt'])
    % A bit of song and dance is required to remove whitespaces in the names
end
