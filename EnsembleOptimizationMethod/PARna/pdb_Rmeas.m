%get end to end distance
%'C:\Users\Suzette\Documents\2021\Single-strand code\MATLAB code\Dseq2_10'
%
nAtoms = 813; % User enters the total number of atoms in the PDB files (accessed by looking in PDB file with notepad)

%%
for n = 1:2066
    namestr = [num2str(n) '_X.pdb']; n
    r_X =  pdbread(namestr);
    x1 = r_X.Model.Atom(1,1).X; %/00001_X///A`1/O5': 1  O5'   A     1      55.990  58.820  37.040
    y1 = r_X.Model.Atom(1,1).Y;
    z1 = r_X.Model.Atom(1,1).Z;
    x2 = r_X.Model.Atom(1,nAtoms).X; %/00001_X///A`30/O3': 972  O3'   A    30      57.530  55.500 127.220
    y2 = r_X.Model.Atom(1,nAtoms).Y;
    z2 = r_X.Model.Atom(1,nAtoms).Z;
    R(n) = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 );
end

%%
R = R';
save Ree_values R
%%

