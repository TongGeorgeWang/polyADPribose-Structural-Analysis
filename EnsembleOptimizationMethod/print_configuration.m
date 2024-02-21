function print_configuration(c)
maxline = 60;                                        
ndigits = floor(log10(max(c)))+1;
nperline = floor((maxline-4)/(ndigits+1));
fmstr = repmat([' %',num2str(ndigits),'d'],1,nperline);
fprintf(1,'configuration:');

if length(c) <= nperline
    fprintf(1,['\n [' fmstr],c);
    %fprintf(1,' ]\n');
else
    fprintf(1,['\n [' fmstr],c(1:nperline));
    fprintf(1,[' ...\n  ' fmstr],c((nperline+1):end));
end
fprintf(1,']\n');
end