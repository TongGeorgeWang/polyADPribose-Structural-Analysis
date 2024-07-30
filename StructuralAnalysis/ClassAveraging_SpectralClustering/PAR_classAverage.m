clear; close all

%% Assign structures within an ensemble to unique classes based on nearest-neighbor RMSD minimization 
%   This workflow requires you to first generate RMSD values of all possible pairwise alignments
%       using PYMOL (can use 'alignAll.pml' script)
%
%   Designed for ensembles of disordered nucleic acid conformers
%   
%   NOTE: due to the variable nature of K-means clustering (the final step
%   of this algorithm), a slightly different clustering answer may be
%   yielded upon each run. We ran the K-means clustering until the graph
%   was clustered in a sensible way. 
%
%   GW - March 2023

Nstructures = 297; %PAR15na


folderName = 'ExampleInputData_PAR15na';
rmsdFile = 'rmsd.txt';

% Threshold RMSD value to be considered in the same class
rmsdThreshold = 11; 

nKmeans = 4; % Number of clusters to use for final K-means on the node/edge graph

%% Load RMSDs and set up matrix 
RMSD_load = readmatrix([folderName,'/',rmsdFile]);
RMSD = RMSD_load(:,2);
RMSD = reshape(RMSD,[Nstructures, Nstructures]);

% Plot heatmap
figure('Renderer', 'painters', 'Position', [10 10 580 540])
h = heatmap(RMSD);
h.GridVisible = 'off';
colormap('jet')
caxis([0,25])  
set(gcf,'color','w')
set(gca,'FontSize',20)
Labels = 1:Nstructures; CustomLabels = string(Labels);
CustomLabels(mod(Labels,5) ~= 0) = " "; % Replace all but the fifth elements by spaces 
h.XDisplayLabels = CustomLabels; h.YDisplayLabels = CustomLabels;
xlabel('structure'); ylabel('structure')
%annotation('textarrow',[1,1],[0.5,0.5],'string','RMSD','Interpreter','latex','Fontsize',20, ...
 %     'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');

  
%% Define binary RMSD pairings (connectivity matrix)
RMSD_binary = zeros(Nstructures, Nstructures); 
for i = 1:Nstructures
    for j = 1:Nstructures
        if RMSD(i,j) <= rmsdThreshold
            RMSD_binary(i,j) = 1;
        else 
            RMSD_binary(i,j) = 0;
        end
    end
end

% Plot connectivity matrix
figure('Renderer', 'painters', 'Position', [10 10 580 540])
h = heatmap(RMSD_binary);
h.GridVisible = 'off';
colormap('jet')
caxis([0,1])  
set(gcf,'color','w')
set(gca,'FontSize',20)
h.XDisplayLabels = CustomLabels; h.YDisplayLabels = CustomLabels;
xlabel('structure'); ylabel('structure')


%% Perform spectral clustering on connectivity matrix; plot graph

figure; hold all
G = graph(RMSD_binary,'omitselfloops','upper');
h = plot(G,'-db','LineWidth',1,'MarkerSize',5);
set(gcf,'color','w')
set(gca,'FontSize',20)
set(gca,'LineWidth',2)
box on 
grid off
h.EdgeColor = [0.5 0.5 0.5];
h.Marker = 'O';

C = SpecClust(RMSD_binary, nKmeans);

% Color graph nodes by cluster
clusters = [1:Nstructures; C']';
clusters = sortrows(clusters,2); 
nClusters = numel(unique(clusters(:,2)));
colors = colormap(jet);
colorSpacing = floor(numel(colors(:,1)) / (nKmeans+1) );
for j = 1:nClusters+1
    clusterIndices{j} = clusters(clusters(:,2)==j-1);
    highlight(h,clusterIndices{j},'NodeColor',colors(j*colorSpacing,:))
end

 


