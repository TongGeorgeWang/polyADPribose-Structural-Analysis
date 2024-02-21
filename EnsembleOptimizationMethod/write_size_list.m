function write_size_list(pdbdir,idx,rg,dmax)
if isempty(pdbdir)
    fn = 'Size_list.txt';
else
    fn = [pdbdir '\Size_list.txt'];
end

if isempty(ls(fn)) || idx==1
    fid = fopen(fn,'w');
else
    fid = fopen(fn,'a');
end
for i=1:length(idx)
    fprintf(fid,'%6d %7.2f %7.2f\n',idx(i),rg(i),dmax(i));
end
fclose(fid);

end