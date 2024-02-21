
%% Compute tortuosity of PAR sequence 
%
%   GW - updated January 2023 to v3
%   Now computes weighted tortuosities 
%   To have ~equilength chain segments, changed atom sampling:  
%       Sample mean of each P atom pair
%       Sample O atom in between ribose sugars
%   
%
%   GW - August 2022
%   Method 1: sample coordinates of P atoms (PA & PB) and O atoms in between the two
%   bases (O1D) as 'points'
%       Take pairwise lengths of adjacent points
%       C = total arc length = total length of pairwise vectors
%       L = end-to-end distance = length of 1st & last coordinates
%
%       Tortuosity (T) = C/L
%
%   Method 2: sample coordinates of P atoms (PA & PB) and O atoms in between the two
%   bases (O1D)
%       Take pairwise lengths of adjacent points
%       SA = sum of all angles between pairwise vectors
%       L = end-to-end distance 
%       
%       Tortuosity (T) = SA/L
%

clear; close all

%% Enter these variables
folderName = 'ExampleInputData_PAR15na';


weightsFilename = 'weights.txt'; %comment this line out if no weights


%% Import data and weights, if defined
PDBnumbers = strsplit( fileread([folderName,'/PDBlist.txt']) ); 
nStructures = numel(PDBnumbers);

if exist('weightsFilename','var')
    w = load([folderName,'/',weightsFilename]);
    w = w(:,2)./sum(w(:,2));
end

%% Sample phosphorus atom coordinates and draw pairwise vectors between adjacent phosphorus atoms 

for name = 1:nStructures
    pdb = pdbread([folderName,'/',PDBnumbers{name},'_X.pdb']);
    pdb = pdb.Model.Atom;
    
    %% Index phosphorus atoms along the chain
    for i = 1:pdb(1,end).AtomSerNo
        atomNames{i} = pdb(1,i).AtomName;
    end
    
    %% Old way (Aug 2022)
%     wherePoints = [find(strcmpi(atomNames,'PA')), find(strcmpi(atomNames,'PB')), find(strcmpi(atomNames,'O1D'))];
%     wherePoints = sort(wherePoints);
%     
%     % Get coordinates of phosphorus atoms
%     for i = 1:numel(wherePoints)
%         Ps_x(i,1) = pdb(1,wherePoints(i)).X ;
%         Ps_y(i,1) = pdb(1,wherePoints(i)).Y ;
%         Ps_z(i,1) = pdb(1,wherePoints(i)).Z ;
%     end
%     S(:,1) = Ps_x;
%     S(:,2) = Ps_y;
%     S(:,3) = Ps_z;
    
    %% New way (Jan 2023)
    wherePhos = [find(strcmpi(atomNames,'PA')), find(strcmpi(atomNames,'PB'))]; wherePhos = sort(wherePhos);
    for i = 1:2:numel(wherePhos)
        Ps_x(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).X, pdb(1,wherePhos(i+1)).X] ) ;
        Ps_y(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).Y, pdb(1,wherePhos(i+1)).Y] ) ;
        Ps_z(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).Z, pdb(1,wherePhos(i+1)).Z] ) ;
    end
    
    whereOs = find(strcmpi(atomNames,'O1D')); whereOs = sort(whereOs);
    for j = 1:numel(whereOs)
        Os_x(1,j) = pdb(1,whereOs(j)).X;
        Os_y(1,j) = pdb(1,whereOs(j)).Y;
        Os_z(1,j) = pdb(1,whereOs(j)).Z;
    end
    
    % % Alternate P and O
    S(:,1) = reshape([Ps_x(1+rem(0:max(numel(Ps_x),numel(Os_x))-1, numel(Ps_x))) ; Os_x(1+rem(0:max(numel(Ps_x),numel(Os_x))-1, numel(Os_x)))],1,[]);
    S(:,2) = reshape([Ps_y(1+rem(0:max(numel(Ps_y),numel(Os_y))-1, numel(Ps_y))) ; Os_y(1+rem(0:max(numel(Ps_y),numel(Os_y))-1, numel(Os_y)))],1,[]);
    S(:,3) = reshape([Ps_z(1+rem(0:max(numel(Ps_z),numel(Os_z))-1, numel(Ps_z))) ; Os_z(1+rem(0:max(numel(Ps_z),numel(Os_z))-1, numel(Os_z)))],1,[]);
    
    
    S_1 = [S; 0 0 0];
    S_2 = [0 0 0; S];
    
    % Set up vectors
    v = S_2 - S_1;
    v = v(2:end-1,:); % truncate first and last points
    
    % Save coordinates for downstream averaging
    S_avgX(:,name) = S(:,1);
    S_avgY(:,name) = S(:,2);
    S_avgZ(:,name) = S(:,3);
    
    
    %% Compute tortuosity index (T) via method 1
    
    magnitudes = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    C = sum(magnitudes);
    vEE = S(end,:)-S(1,:); % end-to-end vector
    L = sqrt(vEE(1).^2 + vEE(2).^2 +vEE(3).^2);
    
    T_method1(name) = C./L ;
    
    
    %% Compute tortuosity index (T) via method 2
    
    v1 = v(1:end-1,:);
    v2 = v(2:end,:);
    angles = acos(dot(v1,v2,2)./...
        (sqrt(v1(:,1).^2 + v1(:,2).^2 +v1(:,3).^2) .* sqrt(v2(:,1).^2 + v2(:,2).^2 +v2(:,3).^2) )); % Based on definition of dot product
    %angles = atan2(norm(cross(v1,v2)),dot(v1,v2)); % More robust calculation at small angles
    anglesDegrees = angles .* (180./pi);
    SA = sum(anglesDegrees);
    
    vEE = S(end,:)-S(1,:); % end-to-end vector
    
    T_method2(name) = SA./L ;
    
end

% Get weighted mean tortuosities across all structures in pool 
if exist('weightsFilename','var')
    T_method1_ensembleMean = sum(w.*T_method1')./sum(w); % weighted mean
    T_method2_ensembleMean = sum(w.*T_method2')./sum(w); 
else
    T_method1_ensembleMean = mean(T_method1,1); % arithmetic mean
    T_method2_ensembleMean = mean(T_method2,1); 
end
T_method1_ensembleStd = std(T_method1',0,1);
T_method2_ensembleStd = std(T_method2',0,1);

disp(['Mean tortuosity via method 1 = ',num2str(T_method1_ensembleMean)])
disp(['StD tortuosity via method 1 = ',num2str(T_method1_ensembleStd)])
disp('')

disp(['Mean tortuosity via method 2 = ',num2str(T_method2_ensembleMean)])
disp(['StD tortuosity via method 2 = ',num2str(T_method2_ensembleStd)])
disp('')


%% Plot backbone positions of all samples (averaged) as a wire frame model, to visualize
figure; hold all
set(gcf,'color','w')
title('Wireframe of PAR backbone (mean positions)')

if exist('weightsFilename','var')
    meanS = [ sum(w.*S_avgX')./sum(w);... 
              sum(w.*S_avgY')./sum(w);... 
              sum(w.*S_avgZ')./sum(w) ]'; % weighted mean
else
    meanS = [mean(S_avgX,2), mean(S_avgY,2), mean(S_avgZ,2)]; % arithmetic mean
end

% Colour point positions by spatial variance across the set of structures (kind of like a crystallographic B-factor map)
varS = [var(S_avgX,0,2), var(S_avgY,0,2), var(S_avgZ,0,2)];
varMetricS = sqrt(varS(:,1).^2 + varS(:,2).^2 + varS(:,3).^2);
scatter3(meanS(:,1),meanS(:,2),meanS(:,3),[],varMetricS,'filled')
a = colorbar; a.Label.String = 'Spatial variance (Angstroms)'; a.Label.FontSize = 15;
plot3(meanS(:,1),meanS(:,2),meanS(:,3),'K-')
plot3(meanS(1,1),meanS(1,2),meanS(1,3),'rs','MarkerFaceColor','red','MarkerSize',10) % show 5' end
plot3(meanS(end,1),meanS(end,2),meanS(end,3),'r^','MarkerFaceColor','red','MarkerSize',10) % show 3' end

%% Plot backbone positions of last sample as a wire frame model, to visualize 
% figure; hold all
% set(gcf,'color','w')
% title('Wireframe of PAR backbone')
% plot3(S(:,1),S(:,2),S(:,3),'KO-') 

%% Plot boxplot of tortuosity values from both methods
figure; hold all
set(gcf,'color','w')
xlabel('Tortuosity index')
boxplot(T_method1',...
    'Labels',{'Arc length method'},...
    'Whisker',1)

figure; hold all
set(gcf,'color','w')
xlabel('Tortuosity index')
boxplot(T_method2',...
    'Labels',{'Sum of angles method'},...
    'Whisker',1)






