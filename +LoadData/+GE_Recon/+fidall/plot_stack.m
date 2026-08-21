function plot_stack(spec,par)
%PLOT_STACK  Display spectra as stacked plot
%  plot_stack(spec,par)
%
%  7/2021 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('par','var'), par = []; end
if isempty(par),        par = struct; end
if ~isreal(spec)
    warning('~isreal(spec): taking abs'); 
    spec = abs(spec);
end


%% frequency axis
[n1,n2,n3]=size(spec);
if isfield(par,'sample_frequency')
    hz = (-n2/2:n2/2-1)/n2*par.sample_frequency;
    ylstr = 'spec [Hz]';    
else
    hz = (-n2/2:n2/2-1);
    ylstr = 'spec [pts]';
end


%% plotting
for l3=1:n3
    if n3>1
        figure(100+l3);
    else
        clf
    end
    x = repmat(hz,[n1 1]).';
    z = repmat((1:n1).',[1 n2]).';
    plot3(z,x,spec(:,:,l3).');
    axis tight
    view(-87,21);
    ylabel(ylstr);

    t_str = sprintf('n1=%g  n2=%g',n1,n2);
    if n3>1, t_str = sprintf('%s   l3=%g',t_str,l3); end
    title(t_str);
    
    if n3>1, set(gcf,'name',sprintf('Slice%d',l3)); end
    grid on
end

end      % plot_stack.m
