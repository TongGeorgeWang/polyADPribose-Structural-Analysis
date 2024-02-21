function print_etime(t0)
tt = toc(t0);
if tt < 2*60
    fprintf(1,'ETIME = %.2f seconds\n',tt);
elseif tt < 2*3600
    fprintf(1,'ETIME = %.2f minutes\n',tt/60);
else
    fprintf(1,'ETIME = %.2f hours\n',tt/3600);
end
end