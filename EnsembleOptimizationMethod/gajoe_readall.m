function out_data = gajoe_readall(gajoedir)

% intensity file names
intensity_file_name = ls([gajoedir '\jun*.int']);

% load intensities
fprintf(1,'\nreading intensity file %s ...\n',intensity_file_name);
[q,iq] = read_intensity_file([gajoedir '\' intensity_file_name]);
fprintf(1,'\b done!\n');

out_data.pool.q = q;
out_data.pool.iq = iq;

% log file  
logfile_names = cellstr(ls([gajoedir '\*.log']));
fprintf(1,'reading result: %2d of %-2d\n',0,length(logfile_names))  
%figure;clf;  
for i=1:length(logfile_names)
    fprintf('\b\b\b\b\b\b\b\b\b%2d of %-2d\n',i,length(logfile_names));
    % read the log file
    gajoe_log = read_gajoe_log([gajoedir '\' logfile_names{i}]);
    
    % generate file namescc
    runidx = str2double(logfile_names{i}(3:5));
    GAdir = [gajoedir '\' sprintf('GA%03d',runidx)];
    
    % SAXS data file
    data_filename = [gajoedir '\' gajoe_log.Experimental_data_file_name];
    
    % size file
    [~,data_prefix] = fileparts(data_filename);
    size_file_name = [gajoedir '\Size_list' data_prefix '.txt'];
    
    % fit file
    fit_file_name = [GAdir '\' gajoe_log.Fit_to_the_experimental_data];
    
    % ensembles file
    ensembles_file_name = [GAdir '\selected_ensem' sprintf('%03d',runidx),'.txt'];
    
    % read saxs data, size list, best fit, and ensembles
    saxs_data = read_saxs_data(data_filename);
    size_list = read_size_list(size_file_name);
    [chi_total,best_fit] = read_fit_file(fit_file_name);
    ensembles = read_ensembles(ensembles_file_name);
    
    % generate rg histogram
%     edgevec_rg = linspace(min([size_list.rg]),max([size_list.rg]),30);
%     edgevec_dmax = linspace(min([size_list.dmax]),max([size_list.dmax]),30);
%     rg_hist = histc([size_list(ensembles(:)).rg],edgevec_rg)/...
%         numel(ensembles);
%     dmax_hist = hist([size_list(ensembles(:)).dmax],edgevec_dmax)/...
%         numel(ensembles);
%     
%     subplot(1,2,1);stairs(edgevec_rg,rg_hist);hold all;
%     subplot(1,2,2);stairs(edgevec_dmax,dmax_hist);hold all;

    out_data.GArun(runidx) = struct(...
        'gajoe_log',gajoe_log,...
        'saxs_data',saxs_data,...
        'chi_total',chi_total,...
        'best_fit',best_fit,...
        'ensembles',ensembles);
    
    out_data.pool.rg = [size_list.rg];
    out_data.pool.dmax = [size_list.dmax];
    
end
fprintf('done!\n\n');
end

function ensembles_out = read_ensembles(ensembles_file_name)

nskip0 = 4;
nskip1 = 2;
m = importdata(ensembles_file_name,' ',nskip0);
nrows = size(m.data,1);
ncols = size(m.data,2);

fid = fopen(ensembles_file_name);

ensembles = zeros(size(m.data));

j=0;
while ~feof(fid)
j = j+1;
for i=1:nskip0
    fgetl(fid);
end
    z = fscanf(fid,'%d',[ncols,nrows])';
    if all(size(z)==[nrows,ncols])
        ensembles(:,:,j) = z;
    else
        break;
    end
for i=1:nskip1
    fgetl(fid);
end
end
fclose(fid);

ensembles_out = squeeze(ensembles(:,1,:));

end

function [chi_total,best_fit] = read_fit_file(fit_file_name)
%fit_file_name = 'gajoe_dA30_c2_00\GA001\profilesda30_mg00001.fit';
nlines = 0;
fid = fopen(fit_file_name);
fgetl(fid); nlines = nlines + 1; % skip one line
while 1
    newl = fgetl(fid); nlines = nlines + 1;
    a = regexpi(newl,'CYCLE:\s+(\d+) Chi-total:\s+([\d\.]+)','tokens');
    if isempty(a),break;end
    a = a{1};
    if isempty(a{1})||isempty(a{2}),break; end
    chi_total(str2double(a{1})) = str2double(a{2});
end
fclose(fid);
m = importdata(fit_file_name,' ',nlines-1);
best_fit = struct('q',m.data(:,1),'int',m.data(:,3));
end

function [prms,msgs] = read_gajoe_log(logfilename)

fid = fopen(logfilename);

prms = struct();
msgs = {};

while ~feof(fid)
    newl = fgetl(fid);
    s2 = strtrim(regexp(newl, '(\.*) : ', 'split'));
    if length(s2)==1
        if ~isempty(regexp(s2{1},'^Files in the pdb folder:','start'))
            break;
        end
        if isempty(regexp(s2{1},'(^Ensemble optimization)|(^Files created:)'))
            msgs = [msgs;s2];
        end
    elseif length(s2)==2
        sfn = regexprep(s2{1},' ','_');
        %fprintf('%s = %s\n',s2{:});
        num_match = regexp(s2{2},'^(\d+)$','match');
        if isempty(num_match)
            prms.(sfn) = s2{2};
        else
            prms.(sfn) = str2double(s2{2});
        end
    end
    %fprintf(1,'%s\n',newl);
end
% s2 = strtrim(regexp(newl, '(\.+) : ', 'split'));
% fprintf('%s = %s\n',s2{:});
fclose(fid);

end


function [q,iq] = read_intensity_file(int_filename)

fid = fopen(int_filename);

a = regexpi(fgetl(fid),'S values\s+(\d+)','tokens');
npts = str2double(a{1});

q = fscanf(fid, '%f', npts);

fgetl(fid);
fgetl(fid);

i = 0;
while ~feof(fid)
    i = i+1;
    iq(:,i) = fscanf(fid, '%f', npts);
    
    fgetl(fid);
    fgetl(fid);
end

fclose(fid);
end



function data_out = read_saxs_data(saxs_datfile_name)
m = dlmread(saxs_datfile_name);
data_out = struct('q',m(:,1),'int',m(:,2),'err',m(:,3));
end


function size_list = read_size_list(size_list_filename)
m = dlmread(size_list_filename,'');
size_list = struct('rg',num2cell(m(:,2)),'dmax',num2cell(m(:,3)));
end
