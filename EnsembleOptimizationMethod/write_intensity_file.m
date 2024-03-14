function write_intensity_file(pdbdir,idx,q,iq)
if isempty(pdbdir)
    fn = 'crysolIntensitiesOfPool.int';
else
    fn = [pdbdir '\crysolIntensitiesOfPool.int'];
end

if isempty(ls(fn)) || idx == 1
    fid = fopen(fn,'w');
    fprintf(fid,'    S values %5d\n',length(q));
    fprintf(fid,'%14.6E\n',q);
else
    fid = fopen(fn,'a');
end

fprintf(fid,' Curve no. %5d\n',idx);
fprintf(fid,'%14.6E\n',iq);

fclose(fid);

end