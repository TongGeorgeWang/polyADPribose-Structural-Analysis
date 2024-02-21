function suite_weights = calc_suite_weights(suites,torsion_names,torsion_weights,Kc2)

%torsion_names = unique({suites.torsion_type});

nTij = zeros(length(suites),length(torsion_names));
for i=1:length(suites)
    for j=1:length(torsion_names)
        nTij(i,j) = strcmp(suites(i).torsion_type,torsion_names{j});
    end
end

nc2i = strcmp({suites.sugar_5p_pucker},'c2');

wTj = torsion_weights;
if ~iscolumn(wTj), wTj = wTj'; end

%Kc2 = 1.9;
wc3 = Kc2*(1-nc2i)*nTij*wTj;
wc2 = (nc2i)*nTij*wTj;
if wc3~=0 && wc2~=0
    suite_weights = nTij*wTj.*(wc3/wc2).^nc2i';
else
    suite_weights = nTij*wTj;
end
suite_weights = suite_weights/sum(suite_weights);
end