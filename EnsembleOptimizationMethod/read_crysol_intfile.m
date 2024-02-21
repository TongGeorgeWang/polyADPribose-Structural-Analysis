function [q,iq] = read_crysol_intfile(pdbdir,intfilename)
if isempty(pdbdir)
    fn = intfilename;
else
    fn = [pdbdir '\' intfilename];
end
a = importdata(fn,' ',1);
q = a.data(:,1);
iq = a.data(:,2);
end