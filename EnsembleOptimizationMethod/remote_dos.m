function [a,b] = remote_dos(p,cmd,waitflag)
if nargin==2
    waitflag = 0;
end

if ~iscellstr(cmd)
    cmd = {cmd};
end

if waitflag==0
    fid = fopen(sprintf('%s/remote_dos_tmp.bat',p),'w');
    fprintf(fid,'cd %s',p);
    fprintf(fid,'\n%s',cmd{:});
    fclose(fid);
    [a,b] = dos([p '\remote_dos_tmp.bat']);
elseif waitflag==1
    fid = fopen(sprintf('%s/remote_dos_tmp.bat',p),'w');
    fprintf(fid,'cd %s',p);
    fprintf(fid,'\n%s',cmd{:});
    fclose(fid);
    dos([p '\remote_dos_tmp.bat &']);
    a = 0; b = '';
elseif waitflag==2
    if dos_is_busy(p)
        error('dos is already executing in that directory?');
    end
    %cmd = ['echo busy > remote_dos_status',cmd];
    cmd = [cmd 'erase remote_dos_status'];
    cmd = [cmd 'exit'];
    fid = fopen(sprintf('%s/remote_dos_tmp.bat',p),'w');
    fprintf(fid,'cd %s',p);
    fprintf(fid,'\n%s',cmd{:});
    fclose(fid);
    dlmwrite([p '/remote_dos_status'],1);
    dos([p '\remote_dos_tmp.bat &']);
    % 12 hour timeout, check once every second
    wait_for_dos(1,43200,p);
elseif waitflag==3
    if dos_is_busy(p)
        error('dos is already executing in that directory?');
    end
    %cmd = ['echo busy > remote_dos_status',cmd];
    cmd = [cmd 'erase remote_dos_status'];
    cmd = [cmd 'exit'];
    fid = fopen(sprintf('%s/remote_dos_tmp.bat',p),'w');
    fprintf(fid,'cd %s',p);
    fprintf(fid,'\n%s',cmd{:});
    fclose(fid);
    dlmwrite([p '/remote_dos_status'],1);
    dos([p '\remote_dos_tmp.bat &']);
    % 10 minute timeout, check 10 times per second
    wait_for_dos(.1,600,p);
end

%delete(sprintf('%s/remote_dos_tmp.bat',p)); <- creates problems sometimes

end

function wait_for_dos(period,timeout,p)
t = timer('Period',period,...
    'TimerFcn',@(x,y) mycbfun(x,p),...
    'ExecutionMode','fixedDelay',...
    'TasksToExecute',ceil(timeout/period));
start(t);
wait(t);

end

function mycbfun(t,p)
    if ~dos_is_busy(p), stop(t); end
end

function y = dos_is_busy(p)
    y = ~isempty(ls([p '/remote_dos_status']));
end