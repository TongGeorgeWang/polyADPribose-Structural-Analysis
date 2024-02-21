function [rg,dmax] = read_crysol_logfile(pdbdir,logfilename)

rg_regexp = 'rg from the slope of net intensity \.+\s*:\s*(?<rg>[\d\.]+)';
dmax_regexp = 'Envelope\s+diameter\s*:\s*(?<dmax>[\d\.]+)';

dmax = NaN; rg = NaN;
if ~isempty(pdbdir)
    fid = fopen([pdbdir '\' logfilename]);
else
    fid = fopen(logfilename);
end
while ~feof(fid);
    s = fgetl(fid);
    match_dmax = regexpi(s,dmax_regexp,'names');
    if ~isempty(match_dmax)
        dmax = str2double(match_dmax.dmax);
        break;
    end
end
while ~feof(fid);
    s = fgetl(fid);
    match_rg = regexpi(s,rg_regexp,'names');
    if ~isempty(match_rg)
        rg = str2double(match_rg.rg);
        break;
    end
end
fclose(fid);
end