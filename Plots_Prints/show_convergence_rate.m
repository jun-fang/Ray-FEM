function r = show_convergence_rate(N,err,xlab,k,opt,str)
%% show convergence rate of an error sequence
%
%  r = show_convergence_rate(N,err) finds the number r such that err = N^r and plots the
%  err vs N in loglog scale.
% 
%  show_convergence_rate(N,err,'h'): show it wrt mesh size h;
%  show_convergence_rate(N,err,'omega'): show it wrt frequency omega;
% 
%  r = showrate(N,err,[],k) finds the number r such that err(k:end)=N(k:end)^r.
%
%  The function accepts standard plotting properting. For example, r =
%  showrate(N,err,[],[],'r') will plot the error curve in red. 
%
%  based on Long Chen's showrate.m

if ~exist('k','var') || isempty(k)
    k = 1; 
end
if ~exist('str','var')
    str = 'Error';
end
err(err == 0) = 1e-16; % Prevent the case err = 0, log(err) = -Inf.
p = polyfit(log(N(k:end)),log(err(k:end)),1);
r = single(p(1));
s = 0.75*err(1)/N(1)^r;
if exist('opt','var')
    h = loglog(N,err,opt);
    set(h,'linewidth',2);
    hold on
    h = loglog(N,s*N.^r,opt);
    set(h,'linewidth',1,'linestyle','--','marker','none');    
else
    loglog(N,err,'-*','linewidth',2);
    hold on
    loglog(N,s*N.^r,'k-.','linewidth',1)
end


if ~exist('xlab','var')
    axis tight;
    xlabel('Number of unknowns','FontSize', 20); ylabel('Error','FontSize', 20);
    title(['Rate of convergence: O(N^{' num2str(r) '})'],'FontSize', 20);
    h_legend = legend(str,['CN^{' num2str(r,2) '}'],'LOCATION','northwest');
    set(h_legend,'FontSize', 14);
elseif xlab == 'omega'
    axis tight;
    xlabel('frequency \omega','FontSize', 20); ylabel('Error','FontSize', 20);
    title(['Error \sim O(\omega^{' num2str(r) '})'],'FontSize', 20);
    h_legend = legend(str,['C \omega^{' num2str(r,2) '}'],'LOCATION','northwest');
    set(h_legend,'FontSize', 14);
elseif xlab == 'h'
    axis tight;
    xlabel('mesh size h','FontSize', 20); ylabel('Error','FontSize', 20);
    title(['Error \sim O(h^{' num2str(r) '})'],'FontSize', 20);
    h_legend = legend(str,['Ch^{' num2str(r,2) '}'],'LOCATION','northwest');
    set(h_legend,'FontSize', 14);
end

