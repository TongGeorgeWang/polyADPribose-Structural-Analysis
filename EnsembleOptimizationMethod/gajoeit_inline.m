function gajoe_output = gajoeit_inline(pdbdir,postfix,datadir,curvenames,generations)

postfix = ['_' postfix]; % changed conventions on this...

print_darwin();
t0 = tic;
print_etime(t0);

if nargin < 1
    error('GAJOEIT: one argument at least is required');
end
if nargin < 5 || isempty(generations)
    generations = 200;
end
if (nargin < 4) || isempty(curvenames)
    curvenames = '*.dat';
end
if nargin < 3 || isempty(datadir)
    datadir = pdbdir;
end
if nargin < 2 || isempty(postfix)
    postfix = ['_' pdbdir];
end

if iscellstr(curvenames)
    fn = curvenames;
else
    fn = cellstr(ls([datadir '\' curvenames]));
end
clear curvenames

if isempty(fn{1})
    error('GAJOEIT: data files not found in data directory');
end

fprintf(1,'\nGAJOEIT input options\n\n');
fprintf(1,'%25s : %-s\n',...
    'pdb directory',pdbdir,...
    'pdb file postfix',postfix,...
    'directory with SAXS data',datadir);
fprintf(1,'%25s : %s\n','data files',fn{1});
if length(fn)>1
    fprintf(1,[sprintf('%25s','') '   %s\n'],fn{2:end});
end
fprintf('%25s : %d\n','number of generations',generations);

for i=1:length(fn) % copy files over
    copyfile([datadir '\' fn{i}],...
            [pdbdir '\' fn{i}]);
    [~,curvenames{i}] = fileparts(fn{i});
end

cmdlist = {'echo off'};

for i=1:length(curvenames)
    curvename = curvenames{i};

    % GAJOE 1.3
    cmdlist = [cmdlist,['copy Size_list.txt Size_list' curvename '.txt']]; % GAJOE requires several RANCH format output files as input
    cmdlist = [cmdlist,'echo 0 > mycommands']; % program mode? 0 = Ensemble Optimization
    cmdlist = [cmdlist,['echo ' curvename '.dat >> mycommands']]; % experimental data file name
    cmdlist = [cmdlist,'echo 1 >> mycommands'];                   % angular units. 1 = inverse angstrom, 0 = inverse nm
    cmdlist = [cmdlist,['echo jun' postfix '00.int >> mycommands']]; % name of the intensities file
    cmdlist = [cmdlist,'echo y >> mycommands']; % "Do you want to change the settings?"
    cmdlist = [cmdlist,'echo 1000 >> mycommands']; % number of generations
    cmdlist = [cmdlist,'echo 50 >> mycommands'];   % number of ensembles
    cmdlist = [cmdlist,'echo 20 >> mycommands'];   % number of curves per ensemble
    cmdlist = [cmdlist,'echo 10 >> mycommands'];   % max number of mutations per ensemble
    cmdlist = [cmdlist,'echo 20 >> mycommands'];   % number of crossings per generation
    cmdlist = [cmdlist,['echo ' num2str(generations) ' >> mycommands']]; % number of times you want the process repeated
    cmdlist = [cmdlist,'echo y >> mycommands'];    % Allow repetitions?
    cmdlist = [cmdlist,'echo n >> mycommands'];    % Allow constant subtraction?
    cmdlist = [cmdlist,'echo n >> mycommands'];    % Create optional files?
    cmdlist = [cmdlist,'echo n >> mycommands'];    % Create analysis files?
    cmdlist = [cmdlist,'..\gajoe13 < mycommands'];
%     cmdlist = [cmdlist,'erase mycommands'];
end

fprintf(1,'launch dos window, evaluate program\n');

remote_dos(pdbdir,cmdlist,2);

fprintf(1,'\nDONE!\n');
print_etime(t0);

gajoe_output = gajoe_readall(pdbdir); %comment out when running crysol_code_sp.m (or similar code)

print_etime(t0);
end

function print_darwin()
darwin_ascii = [...
'               _....._   ';...
'       ''.  .-''`       `. ';...
'         ><  G A J O E  )';...
'       .''  `-.,_     _,'' ';...
'                `L`L`    '];
disp(darwin_ascii);
fprintf(1,'\n');
end