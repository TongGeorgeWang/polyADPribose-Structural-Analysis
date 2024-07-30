
%% Compute OCF using Alex P's method (adapted from 'JJ_OCF.m') 
%
%   GW - updated January 2023 to v3
%   Now computes weighted OCF  
%   Also added computation of correlation length 
%   To have ~equilength chain segments, changed atom sampling:  
%       Sample mean of each P atom pair
%       Sample O atom in between ribose sugars
%
%   GW - August 2022 

clear; close all 

%% Fill out these variables 

folderName = 'ExampleInputData_PAR15na';

nAtoms = 813;
weightsFilename = 'weights.txt'; %comment this line out if no weights


%% Import data and weights, if defined
%PDBnumbers = strsplit( fileread([folderName,'/PDBlist.txt']) );
PDBnumbers = load([folderName,'/PDBlist.txt']);
nStructures = numel(PDBnumbers);

if exist('weightsFilename','var')
    w = load([folderName,'/',weightsFilename]);
    w = w(:,2)./sum(w(:,2));
end

for name = 1:nStructures 
%pdb = pdbread([folderName,'/',PDBnumbers{name},'_X.pdb']);
pdb = pdbread([folderName,'/',num2str(PDBnumbers(name)),'_X.pdb']);
pdb = pdb.Model.Atom;

%% Index phosphorus atoms along the chain
for i = 1:pdb(1,end).AtomSerNo
    atomNames{i} = pdb(1,i).AtomName; 
end
wherePhos = [find(strcmpi(atomNames,'PA')), find(strcmpi(atomNames,'PB'))];
%wherePhos = find(strcmpi(atomNames,'PA'));
%wherePhos = find(strcmpi(atomNames,'PB'));
wherePhos = sort(wherePhos);

%% Get vectors at phosphorus atoms - index each P atom separately
% for i = 1:numel(wherePhos)
%     Ps_x(i,1) = pdb(1,wherePhos(i)).X ;
%     Ps_y(i,1) = pdb(1,wherePhos(i)).Y ;
%     Ps_z(i,1) = pdb(1,wherePhos(i)).Z ;
% end 
% 
% S(:,1) = Ps_x;
% S(:,2) = Ps_y;
% S(:,3) = Ps_z;

%% Get vectors at phosphorus atoms - take average of each pair of P atoms
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


%% Compute vectors between sampled atoms
S_1 = [S; 0 0 0];
S_2 = [0 0 0; S];

% Set up vectors
v = S_2 - S_1; 
v = v(2:end-1,:); 
% Normalize 
vl = sqrt(sum(v.^2,2));
vn = v./repmat(vl,[1,size(v,2),1]);
nv = zeros(size(vn,1)-1,size(vn,3));

% Compute P-O displacements (b value in correlation length computation)
bVec = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
bMean(name) = mean(bVec);

%% Compute average of cos(theta) vs separation, for each chain
for j = 1:(size(vn,1)-1)
    nv(j,:) = squeeze(mean((vn(1:(end-j),1,:).*vn((j+1):end,1,:) + ...
        vn(1:(end-j),2,:).*vn((j+1):end,2,:) + ...
        vn(1:(end-j),3,:).*vn((j+1):end,3,:)),1)); 
end

OCF(name,:) = nv(:);
end

% Get weighted mean OCF across all structures in pool 
if exist('weightsFilename','var')
    OCF_ensembleMean = sum(w.*OCF)./sum(w); % weighted mean
else
    OCF_ensembleMean = mean(OCF,1); % arithmetic mean
end
OCF_ensembleVar = var(OCF,0,1);

% OCF_ensembleMean = mean(OCF,1);
% OCF_ensembleVar = var(OCF,0,1);
% OCF_ensembleStErr = std(OCF,0,1) ./ sqrt(numel(OCF)) ; 


%% Compute correlation length
if exist('weightsFilename','var')
    b = sum(w.*bMean')./sum(w); % weighted mean
else
    b = mean(bMean); % arithmetic mean
end
b_err = var(bMean);
l_OCF = b.*sum(OCF_ensembleMean)
l_OCF_err = l_OCF.*sqrt( (b_err./b).^2 + (mean(OCF_ensembleVar)./l_OCF).^2 ) % propagate error multiplicatively 


%% Visualize |i-j| vs OCF_ij, AKA <cos(th)ij>
x = 1:1:numel(OCF_ensembleMean);
figure; hold all 
errorbar(x,OCF_ensembleMean,OCF_ensembleVar,'k-o')
plot(x,0.*x,'k--','LineWidth',3)
set(gcf,'Color','w')
grid on; box on
xlabel('|i-j|','FontSize',22); ylabel('$<cos(\theta)_{ij}>$', 'Interpreter','latex','FontSize',22)
title(['Weighted OCF for ',folderName,'; N = ',num2str(nStructures)],'Interpreter','none')

   
   