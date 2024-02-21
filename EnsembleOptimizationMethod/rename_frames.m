%rename frames to correct format for gajoe script 

startNum = 0;
endNum = 2066;
basename = 'f';

for n = startNum:endNum
    oldname = [basename num2str(n) '.pdb'];
    newname = [num2str(n) '_X.pdb'];
    copyfile(oldname,newname);
end

%%
delete f*
%%