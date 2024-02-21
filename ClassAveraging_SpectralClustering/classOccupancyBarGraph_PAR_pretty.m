clear; close all

%% Plot data for PAR manuscript figures as bar graphs with errorbars 
%   GW - May 2023
%
%   Values updated to reflect Round 2 of data processing 
%
% Color guide:
%   Blue: #0800ff
%   Dark turquoise: #12abde
%   Light turquoise: #12dbd5
%   Light green: #83f25e
%   Yellow green: #89c910
%   Yellow: #cee310
%   Light orange: #f7ad23
%   Dark orange: #d95a0b
%   Mahogany: #91220f
%
%   Requires either the hex2rgb.m script (from FileExchange), or for you to
%   manually specify the colors in rgb form. 

fontSize = 30;
barWidth = 0.5;
sideways = 0; % To flip or not to flip 


%% PAR15na

figure('Name','PAR15na')
names = {'0','1','2','3'};

x = [1:4]; 
occ = [16.8,34.0,17.8,31.3]; 

colors = {hex2rgb('#12dbd5'),hex2rgb('#cee310'),hex2rgb('#0800ff'),hex2rgb('#d95a0b')};
b = bar(occ,barWidth);
b.FaceColor = 'flat';
for k = 1:numel(x)
    b.CData(k,:) = colors{k};
end
b.LineWidth = 2;

%box on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',fontSize)
ylim([14,35])
set(gca,'xticklabel',names)
ylabel('Prevalence (%)','FontSize',fontSize)

ylim([0,40])

if sideways == 1
    camroll(-90)
end


%% PAR15naMg

figure('Name','PAR15naMg')
names = {'0','1','2','3','4'};

x = [1:5]; 
occ = [10.0,29.4,22.2,33.8,4.6]; 

colors = {hex2rgb('#d95a0b'),hex2rgb('#0800ff'),hex2rgb('#12dbd5'),hex2rgb('#83f25e'),hex2rgb('#cee310')};
b = bar(occ,barWidth);
b.FaceColor = 'flat';
for k = 1:numel(x)
    b.CData(k,:) = colors{k};
end
b.LineWidth = 2;

box on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',fontSize)
ylim([14,35])
set(gca,'xticklabel',names)
ylabel('Prevalence (%)','FontSize',fontSize)

ylim([0,40])

if sideways == 1
    camroll(-90)
end



%% PAR22na

figure('Name','PAR22na')
names = {'0','1','2','3','4'};

x = [1:5]; 
occ = [21.3,25.3,17.3,20.0,16.0]; 

colors = {hex2rgb('#83f25e'),hex2rgb('#0800ff'),hex2rgb('#cee310'),hex2rgb('#d95a0b'),hex2rgb('#12dbd5')};
b = bar(occ,barWidth);
b.FaceColor = 'flat';
for k = 1:numel(x)
    b.CData(k,:) = colors{k};
end
b.LineWidth = 2;

box on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',fontSize)
ylim([14,35])
set(gca,'xticklabel',names)
ylabel('Prevalence (%)','FontSize',fontSize)

ylim([0,40])

if sideways == 1
    camroll(90)
end


%% PAR22naMg

figure('Name','PAR22naMg')
names = {'0','1','2','3','4'};

x = [1:5]; 
occ = [21.8,24.4,9.0,23.1,12.8]; 

colors = {hex2rgb('#89c910'),hex2rgb('#12abde'),hex2rgb('#12dbd5'),hex2rgb('#91220f'),hex2rgb('#f7ad23')};
b = bar(occ,barWidth);
b.FaceColor = 'flat';
for k = 1:numel(x)
    b.CData(k,:) = colors{k};
end
b.LineWidth = 2;

box on
set(gcf,'color','w')
set(gca,'LineWidth',2)
set(gca,'FontSize',fontSize)
ylim([14,35])
set(gca,'xticklabel',names)
ylabel('Prevalence (%)','FontSize',fontSize)

ylim([0,40])

if sideways == 1
    camroll(90)
end


