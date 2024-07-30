
%% Compute number of base stacking events of PAR sequence
%   GW - January 2023 v3: 
%           added weights
%           made pairwise distance maps publication quality 
%
%   GW - August 2022
%   
%   Method: index atoms in the adenine bases
%           fit a plane to each base
%           take the normal vector of this plane
%           Criteria for base stacking:
%               i) normal vectors are reasonably colinear (<45 degrees) 
%               ii) normal vectors are in close proximity (~3.4A is the formal limit, but probably need to go a bit higher)
%
%   Identities of base atoms:
%       N9, C8, N7, C5, C4, N3, C2, N1, C6
%

clear; close all

%% Enter these variables
folderName = 'ExampleInputData_PAR15na';

nBaseAtoms = 9; 
nBases = 15;

baseStackDistanceThreshold = 5; % in Angstroms
baseStackAngleThreshold = 45; % in degrees 

weightsFilename = 'weights.txt'; %comment this line out if no weights


%% Import data
%PDBnumbers = strsplit( fileread([folderName,'/PDBlist.txt']) ); 
PDBnumbers = load([folderName,'/PDBlist.txt']);
nStructures = numel(PDBnumbers);

if exist('weightsFilename','var')
    w = load([folderName,'/',weightsFilename]);
    w = w(:,2)./sum(w(:,2));
end


%% Sample atoms in adenine bases

for name = 1:nStructures
    
    pdb = pdbread([folderName,'/',num2str(PDBnumbers(name)),'_X.pdb']);
    pdb = pdb.Model.Atom;
    
    %% Index base atoms along the chain
    for i = 1:pdb(1,end).AtomSerNo
        atomNames{i} = pdb(1,i).AtomName;
    end
    
    wherePoints = [find(strcmpi(atomNames,'N9')), find(strcmpi(atomNames,'C8')), find(strcmpi(atomNames,'N7')),...
        find(strcmpi(atomNames,'C5')), find(strcmpi(atomNames,'C4')), find(strcmpi(atomNames,'N3')),...
        find(strcmpi(atomNames,'C2')), find(strcmpi(atomNames,'N1')), find(strcmpi(atomNames,'C6'))];
    wherePoints = sort(wherePoints);
    
    % Get coordinates of base atoms
    for i = 1:numel(wherePoints)
        Ps_x(i,1) = pdb(1,wherePoints(i)).X ;
        Ps_y(i,1) = pdb(1,wherePoints(i)).Y ;
        Ps_z(i,1) = pdb(1,wherePoints(i)).Z ;
    end
    S(:,1) = Ps_x;
    S(:,2) = Ps_y;
    S(:,3) = Ps_z;
    
    
    %% Individuate each base, fit a plane to each base
    for i = 1:nBases
        S_indv = S(nBaseAtoms*(i-1)+1:nBaseAtoms*i,:);
        tails(i,:) = mean(S_indv,1); % Base points of normal vectors to each plane - xx
        X = [ones(size(S_indv(:,1))), S_indv(:,1), S_indv(:,2)];
        y = S_indv(:,3);
        [b{name},bint{name},r{name},rint{name},stats{name}] = regress(y,X); % Fit plane using bi-linear regression
        coeff = b{name};
        
        x1fit = min(S_indv(:,1)):0.01:max(S_indv(:,1));
        x2fit = min(S_indv(:,2)):0.01:max(S_indv(:,2));
        [X1FIT{i},X2FIT{i}] = meshgrid(x1fit,x2fit);
        
        YFIT{i} = coeff(1) + coeff(2)*X1FIT{i} + coeff(3)*X2FIT{i};
        [nx,ny,nz] = surfnorm(X1FIT{i},X2FIT{i},YFIT{i}); % Determine normal vector for each plane
        normVectors{i} = [nx(1,1),ny(1,1),nz(1,1)];
        heads(i,:) = tails(i,:) + 5.*normVectors{i};
        
    end
    
    
    %% Determine pairwise distances between each average point in each base & pairwise angles between normals
    
    %baseDistances = pdist(tails);
    for i=1:numel(normVectors)
        for j=1:numel(normVectors)
            
            % Distances
            connection = tails(i,:) - tails(j,:);
            baseDistances(i,j) = abs(norm(connection));
            
            % Angles
            v1 = normVectors{i};
            v2 = normVectors{j};
            baseAngles(i,j) = acos(dot(v1,v2,2)./...
                (sqrt(v1(:,1).^2 + v1(:,2).^2 +v1(:,3).^2) .* sqrt(v2(:,1).^2 + v2(:,2).^2 +v2(:,3).^2) )); % Based on definition of dot product
            
            if i==j % Set zeros to NaN to avoid cointing them as base stacking events
                baseDistances(i,j)=NaN;
                baseAngles(i,j)=NaN;
            end
                    
        end
    end
    
    % Save baseDistances for downstream averaging
    baseDistancesAll{name} = baseDistances; 
    
    % Delineate if base-stacked if bases are <=3.4A apart, and normals are reasonably colinear
    baseAngles = real(baseAngles);
    baseAnglesDegrees = baseAngles .* (180./pi);
    obtuse = find(baseAnglesDegrees>90);
    baseAnglesDegrees(obtuse) = baseAnglesDegrees(obtuse)-180;
    
    [whereBaseStacked_row{name}, whereBaseStacked_column{name}] = find(baseDistances<=baseStackDistanceThreshold & abs(baseAnglesDegrees)<=baseStackAngleThreshold);
    %nBaseStacked(name) = numel(whereBaseStacked_row{name});
    nBaseStacked(name) = numel(unique(whereBaseStacked_row{name})); % only include unique stacking events 

    disp(['# of base stacking events in structure ',num2str(name),' = ',num2str(nBaseStacked(name))])
    
end

if exist('weightsFilename','var')
    mean_nBaseStacked = sum(w.*nBaseStacked')./sum(w); % weighted mean
else
    mean_nBaseStacked = mean(nBaseStacked',1); % arithmetic mean
end

meanFracBaseStacked = mean_nBaseStacked ./ nBases .* 100;
stErrFracBaseStacked = ( std(nBaseStacked)./sqrt(numel(nBaseStacked)) ) ./ nBases .* 100; 
disp(['On average, ',num2str(meanFracBaseStacked),' +/- ',num2str(stErrFracBaseStacked),' % of bases are stacked.'])


%% Plot positions of base atoms and planes for the last iteration for visualization purposes
figure; hold all
set(gcf,'color','w')
title('Base positions of PAR molecule')
axis equal

for j = 1:nBases
    H = surf(X1FIT{j},X2FIT{j},YFIT{j});
    H.EdgeColor = 'none';
    % Visualize normal vectors from tail points
    xNorm = [tails(j,1),heads(j,1)]; yNorm = [tails(j,2),heads(j,2)]; zNorm = [tails(j,3),heads(j,3)];
    
    % Colour bases that are stacked red; unstacked black 
    if ismember(j,whereBaseStacked_row{name})
        H.FaceColor = [1 0.2 0.2];
        plot3(xNorm,yNorm,zNorm,'r','LineWidth',3) 
    else
        H.FaceColor = [0.8 0.8 0.8];
        plot3(xNorm,yNorm,zNorm,'k','LineWidth',3) 
    end
    
end

plot3(S(:,1),S(:,2),S(:,3),'K*') % Plot base positions
plot3(tails(:,1),tails(:,2),tails(:,3),'K--') % Plot connectivity for visualization purposes 
plot3(tails(1,1),tails(1,2),tails(1,3),'rO','MarkerSize',10) % Mark where 5' and 3' - most bases are 
plot3(tails(end,1),tails(end,2),tails(end,3),'r^','MarkerSize',10)


%% Plot heatmap of pairwise base stacking instances, averaged across all iterations

if exist('weightsFilename','var')
    w3D = reshape(w,[1,1,numel(w)]);
    meanBaseDistances = sum(w3D.*cat(3, baseDistancesAll{:}), 3)./sum(w); % weighted mean
else
    meanBaseDistances = mean(cat(3, baseDistancesAll{:}), 3); % arithmetic mean
end

meanBaseDistances(isnan(meanBaseDistances))=0; % set NaNs to 0

figure('Renderer', 'painters', 'Position', [10 10 650 650])
h = heatmap(meanBaseDistances);
colormap('jet')
caxis([0,90])  
set(gcf,'color','w')
set(gca,'FontSize',20)
xlabel('base'); ylabel('base')
annotation('textarrow',[1,1],[0.5,0.5],'string','$\AA$','Interpreter','latex','Fontsize',20, ...
      'HeadStyle','none','LineStyle','none','HorizontalAlignment','center');
title(['N = ', num2str(nStructures)])


%% Plot heatmap of pairwise base stacking instances of last iteration
% %baseStackingIndex = baseDistances .* abs(baseAnglesDegrees) ./ (baseStackDistanceThreshold.*baseStackAngleThreshold); 
% baseStackingIndex = baseDistances;
% figure
% heatmap(baseStackingIndex)
% colormap('jet')
% %caxis([0,1])
% set(gcf,'color','w')
% xlabel('base'); ylabel('base')
% title('Pairwise distances between bases in PAR molecule')


%% Plot where base stacks are occurring (between which pairs of bases) 
whereBaseStacked_array(:,1) = vertcat(whereBaseStacked_column{:});
whereBaseStacked_array(:,2) = vertcat(whereBaseStacked_row{:});
pairwiseBaseStackingMatrix = zeros(nBases);

idx_BS = sub2ind(size(pairwiseBaseStackingMatrix),...
            whereBaseStacked_array(:,1),whereBaseStacked_array(:,2) );

%for n = 1:numel(whereBaseStacked_array(:,1))
%    pairwiseBaseStackingMatrix(idx_BS(n)) = pairwiseBaseStackingMatrix(idx_BS(n))+1; 
%end
pairwiseBaseStackingMatrix(idx_BS) = pairwiseBaseStackingMatrix(idx_BS)+1; 

pairwiseBaseStackingMatrix(1:nBases+1:end) = NaN; % set diagonals to NaN

figure('Renderer', 'painters', 'Position', [10 10 420 350])
heatmap(pairwiseBaseStackingMatrix)
colormap('jet')
%caxis()
set(gcf,'color','w')
set(gca,'FontSize',20)
xlabel('base'); ylabel('base')



