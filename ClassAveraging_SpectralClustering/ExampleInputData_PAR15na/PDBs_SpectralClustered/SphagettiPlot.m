clear
%close all

%% Plot all structures within a spectral cluster to assess similarity and gross features
%   Requires you to organize PDBs into subfolders based on spectral
%   clustering determination
%
%   Uses phosphorus point sampling method from tortuosity scripts  
%
%   GW - 2023 March


%% Enter these variables relating to the PDB inputs
folderName = 'Cluster3/Aligned';
basename = 'PAR15naPDB_';
%PDBnumbers = [1,10,11,12,125,13,131,132,133,134,135,137,138,139,14,140,141,142,143,144,145,146,147,148,149,15,150,151,152,153,154,155,156,157,158,159,16,160,161,17,18,190,237,261,262,3,30,4,5,6,7,8,9]; %C0
%PDBnumbers = [136,162,188,19,191,192,2,216,217,218,220,228,23,238,239,24,245,246,247,248,249,25,250,26,260,263,264,267,268,269,27,270,271,272,277,278,279,28,280,281,282,283,284,285,29,293,296,297,31,32]; %C1
%PDBnumbers = [100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,126,127,128,129,130,163,265,266,289,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99]; %C2
PDBnumbers = [164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,189,193,194,195,196,197,198,199,20,200,201,202,203,204,205,206,207,208,209,21,210,211,212,213,214,215,219,22,221,222,223,224,225,226,227,229,230,231,232,233,234,235,236,240,241,242,243,244,251,252,253,254,255,256,257,258,259,273,274,275,276,286,287,288,290,291,292,294,295]; %C3

nAtoms = 813;

% For some reason the PAR subunit indices get reordered at this point in
% the workflow; need to manually look in the pdb and reorder.
trueOrder = [1;2;8;9;10;11;12;13;14;15;3;4;5;6;7];

%%
nStructures = numel(PDBnumbers);
figure; hold all
set(gcf,'color','w')
set(gcf,'Position',[10 10 1150 1200])


%% Sample phosphorus atom coordinates and draw pairwise vectors between adjacent phosphorus atoms 
for name = 1:nStructures
    pdb = pdbread([folderName,'/',basename,num2str(PDBnumbers(name)),'.pdb']);
    pdb = pdb.Model.Atom;
    
    pdbMatrix = zeros([nAtoms,4]);
    % Need to resort PDB because bases are out of order
    pdbMatrix(:,1) = [pdb(1:nAtoms).resSeq]';
    pdbMatrix(:,2) = [pdb(1:nAtoms).X]';
    pdbMatrix(:,3) = [pdb(1:nAtoms).Y]';
    pdbMatrix(:,4) = [pdb(1:nAtoms).Z]';
    
    pdbMatrixSorted = sort(pdbMatrix,1);
    
    order = unique([pdb(1:nAtoms).resSeq]','stable');
    
    %% Index phosphorus atoms along the chain
    for i = 1:pdb(1,end).AtomSerNo
        atomNames{i} = pdb(1,i).AtomName;
    end
    
    %% Sample average of each pair of P atoms
    wherePhos = [find(strcmpi(atomNames,'PA')), find(strcmpi(atomNames,'PB'))]; wherePhos = sort(wherePhos);
    for i = 1:2:numel(wherePhos)
        Ps_x(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).X, pdb(1,wherePhos(i+1)).X] ) ;
        Ps_y(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).Y, pdb(1,wherePhos(i+1)).Y] ) ;
        Ps_z(1,(i+1)/2) = mean( [pdb(1,wherePhos(i)).Z, pdb(1,wherePhos(i+1)).Z] ) ;
    end
    
    %Reorder
    Ps_x_sorted = Ps_x(trueOrder); Ps_y_sorted = Ps_y(trueOrder); Ps_z_sorted = Ps_z(trueOrder);

% If you want to sample middle O atoms too, uncomment this section (will require some unfucking)
%     whereOs = find(strcmpi(atomNames,'O1D')); %whereOs = sort(whereOs);
%     for j = 1:numel(whereOs)
%         Os_x(1,j) = pdb(1,whereOs(j)).X;
%         Os_y(1,j) = pdb(1,whereOs(j)).Y;
%         Os_z(1,j) = pdb(1,whereOs(j)).Z;
%     end
%     
%     %Reorder
%     Os_x_sorted = Os_x(trueOrder); Os_y_sorted = Os_y(trueOrder); Os_z_sorted = Os_z(trueOrder);
%     
%     % % Alternate P and O
%     S(:,1) = reshape([Ps_x_sorted(1+rem(0:max(numel(Ps_x_sorted),numel(Os_x_sorted))-1, numel(Ps_x_sorted))) ; Os_x(1+rem(0:max(numel(Ps_x_sorted),numel(Os_x_sorted))-1, numel(Os_x_sorted)))],1,[]);
%     S(:,2) = reshape([Ps_y_sorted(1+rem(0:max(numel(Ps_y_sorted),numel(Os_y_sorted))-1, numel(Ps_y_sorted))) ; Os_y(1+rem(0:max(numel(Ps_y_sorted),numel(Os_y_sorted))-1, numel(Os_y_sorted)))],1,[]);
%     S(:,3) = reshape([Ps_z_sorted(1+rem(0:max(numel(Ps_z_sorted),numel(Os_z_sorted))-1, numel(Ps_z_sorted))) ; Os_z(1+rem(0:max(numel(Ps_z_sorted),numel(Os_z_sorted))-1, numel(Os_z_sorted)))],1,[]);
% 
%     % Save coordinates for downstream averaging
%     S_allX(:,name) = S(:,1);
%     S_allY(:,name) = S(:,2);
%     S_allZ(:,name) = S(:,3);
    
    S_allX(:,name) = Ps_x_sorted;
    S_allY(:,name) = Ps_y_sorted;
    S_allZ(:,name) = Ps_z_sorted;
    %plot3(S_allX(:,name),S_allY(:,name),S_allZ(:,name),'.-','Color',[0.85 0.85 0.85])  
    
end

%% Plot backbone positions of all samples (averaged) 

%figure; hold all
%set(gcf,'color','w')
%title('Wireframe of mean backbone position')

meanS = [mean(S_allX,2), mean(S_allY,2), mean(S_allZ,2)];

% if exist('weightsFilename','var')
%     meanS = [ sum(w.*S_avgX')./sum(w);... 
%               sum(w.*S_avgY')./sum(w);... 
%               sum(w.*S_avgZ')./sum(w) ]'; % weighted mean
% else
%     meanS = [mean(S_avgX,2), mean(S_avgY,2), mean(S_avgZ,2)]; % arithmetic mean
% end

% Colour point positions by spatial variance across the set of structures (kind of like a crystallographic B-factor map)
varS = [var(S_allX,0,2), var(S_allY,0,2), var(S_allZ,0,2)];
varMetricS = sqrt(varS(:,1).^2 + varS(:,2).^2 + varS(:,3).^2);
scatter3(meanS(:,1),meanS(:,2),meanS(:,3),300,varMetricS,'filled','MarkerEdgeColor','K','LineWidth',1)
format long g
a = colorbar; a.Label.String = 'Spatial variance (Angstroms)'; a.Label.FontSize = 35;
colormap(pink)
set(gca,'ColorScale','log')
a.FontSize = 30;
caxis([0, 300])
plot3(meanS(:,1),meanS(:,2),meanS(:,3),'K-','LineWidth',1.5)
plot3(meanS(1,1),meanS(1,2),meanS(1,3),'rs','MarkerFaceColor','red','MarkerSize',20) % show 5' end
plot3(meanS(end,1),meanS(end,2),meanS(end,3),'r^','MarkerFaceColor','red','MarkerSize',10) % show 3' end


% Save coordinates 
wireframe = [meanS(:,1),meanS(:,2),meanS(:,3)];

