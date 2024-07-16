function fa = vfa(n,flip_final,fname,verb)
% VFA Design variable flip angle list neglecting relaxation
%  fa = vfa(n,flip_final)
%                                                          [unit] (default)
%          n  Number of excitation pulses
% flip_final  Flip angle of last RF pulse                  [deg]  (90)
%      fname  File name for storing as vap_flip_XXX.fdl           ('')
%       verb  Verbose mode                                        (2)
%             0=off; 1=verbose; 2=verbose+plotting
%
%         fa  Array of flip angles                         [deg]
%
% Literature: MRI book, Haacke; p464; Eq18.37.
% 2/2021 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('flip_final','var'), flip_final = []; end
if isempty(flip_final),        flip_final = 90; end
if ~exist('fname','var'),      fname = ''; end
if ~exist('verb','var'),       verb = []; end
if isempty(verb),              verb = 2; end
if flip_final>90, warning('flip_final(=%g)>90',flip_final); end
flip_final = flip_final*pi/180;        % convert to [rad]


%% calculate variable flip angle list
fa = zeros(1,n);
fa(end) = flip_final;
for l=(n-1):-1:1
    fa(l) = atan(sin(fa(l+1)));
end


%% check for constant signal
if verb>0
    signal = zeros(1,n);
    mz=1;
    for l=1:n
        signal(l) = sin(fa(l))*mz;
        mz = mz*cos(fa(l));
    end
    fprintf('SNR = %g [%%]\n',sum(signal)/sqrt(n)*100);
    
    if verb>1
        subplot(2,1,1);
        plot(fa*180/pi,'b-x');
        ylabel('flip [deg]');
        axis([0 n 0 flip_final*180/pi]);
        grid on
        subplot(2,1,2);
        plot(signal,'b-x');
        xlabel('excitation'); ylabel('signal (Mz(1)=1)');
        axis([0 n 0 max(signal)*1.5]);
        grid on
    end
end


%% convert to [deg]
fa = fa*180/pi;


%% store as vap file
if ~isempty(fname)
    if islogical(fname)
        if fname
            fname = sprintf('vap_vfa%d_%d.fdl',n,round(flip_final*180/pi));
        else
            return
        end
    end
    if ischar(fname)
        if verb>0, fprintf('Saving vfa to ''%s''\n',fname); end
        write_fdl(fname,fa,'flip');
    else
        fprintf('~ischar(fname); skipping saving');
    end
end


end      % vfa.m
