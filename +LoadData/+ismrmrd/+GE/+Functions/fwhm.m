function [lw,out] = fwhm(spec,bw,freqs,verb,plt,frac_max)
%FWHM Calculate spectral linewidth (full width at half maximum)
%[lw,out] = fwhm(spec,bw,freqs,verb,plt,frac_max)
%        spec   Spectra (2D; spectral dimension along 2)
%               if complex -> fwhm of magnitude - corrected by 1/sqrt(3)
%          bw   bandwidth [Hz]
%       freqs   (approximate) positions of maxima
%               if empty->take global maxima                    ([])
%        verb   verbose output                                  (true)
%         plt   plotting                                        (true)
%    frac_max   Calculate width at max fraction (half max=0.5)  (0.5)
%
%          lw   linewidth [Hz] (scaled to bw)
%         out   output structure
%                    .fp  pos. half freq
%                    .fm  neg. half max freq
%               .ind_max  index to maximum
%                 .freqs  frequency of maximum
%
% 11/2020 Rolf Schulte
if nargin<2, help(mfilename); return; end

%% input
if ~exist('bw','var'),    bw = []; end
if isempty(bw),           bw = 1; end
if ~exist('freqs','var'), freqs = []; end
if ~exist('verb','var'),  verb = []; end
if isempty(verb),         verb = true; end
if ~exist('plt','var'),   plt = []; end
if isempty(plt),          plt = true; end
si = size(spec);
if length(si)~=2, error('spec must be 2D'); end
if si(2)==1,      error('size(spec,2)==1'); end
if length(bw)~=1, error('bw not scalar value'); end
if ~exist('frac_max','var'), frac_max = []; end
if isempty(frac_max),        frac_max = 0.5; end
if frac_max>0.99, error('frac_max(=%g)>0.99',frac_max); end
if frac_max<0.01, error('frac_max(=%g)<0.01',frac_max); end


x_ax = -(-si(2)/2:(si(2)/2-1))/si(2)*bw;
mag_corr = false;
if ~isreal(spec)
    if verb
        fprintf('Attention: spec complex: taking magnitude; correcting lw by 1/sqrt(3)\n');
    end
    spec = abs(spec);
    mag_corr = true;
end

%%
if all(abs(spec)<1d-10)
   warning('spec==0: exiting fwhm.m');
   lw = 0;
   if nargout>1
       out.fp = 0; out.fm = 0; out.ind_max = 1; out.freqs = 0;
   end
   return
end


%% locating + refining indices of maxima
if isempty(freqs)       % determine global maxima of spectra
    [~,ind_max] = max(spec,[],2);
    freqs = zeros(size(ind_max));
    for l1=1:si(1)
        freqs(l1,1) = x_ax(ind_max(l1,1));
        if verb
            fprintf('l1=%g: ind_max=%g; freq=%g\n',l1,ind_max(l1,1),freqs(l1,1)); 
        end
    end
else                    % if freqs given, look for nearest local maxima
    ind_max = zeros(size(freqs));
    for l1=1:size(freqs,1)
        for l2=1:size(freqs,2)
            [~,im] = min(abs(x_ax-freqs(l1,l2)));
            ll = 0;
            while true
                if (im==1)
                    fprintf('Warning1: im==1; l1=%g,l2=%g; break\n',l1,l2); 
                    break; 
                end
                if (im==si(2))
                    fprintf('Warning1: im==%g; l1=%g,l2=%g; break\n',si(2),l1,l2); 
                    break; 
                end
                if (ll>si(2))
                    fprintf('Warning1: excessive iterations (%g); l1=%g,l2=%g; break\n',ll,l1,l2); 
                    break; 
                end
                [~,ii] = max(spec(l1,(im+(-1:1))));
                if (ii==2), break; end
                im = im+ii-2;
                ll = ll+1;
            end
            if verb
                fprintf('l1=%g; l2=%g: ind_max=%g; freq=%g\n',...
                    l1,l2,im,x_ax(im));
            end
            ind_max(l1,l2) = im;
            freqs(l1,l2) = x_ax(im);
        end
    end
end

%% actual linewidth determination
lw = zeros(size(ind_max));   % linewidths
fp = lw;                     % freqs at half maximum in plus direction 
fm = lw;                     % freqs at half maximum in minus direction 
for l1=1:size(ind_max,1)
    for l2=1:size(ind_max,2)
        im = ind_max(l1,l2);
        val = spec(l1,im);   % value of full maximum
        ind_hp = im;         % index to >half maximum
        ll = 0;
        while true           % walk up
            if (ind_hp==si(2))
                fprintf('Warning2: ind_hp==%g; l1=%g,l2=%g; break\n',si(2),l1,l2); 
                break; 
            end
            if (ll>si(2))
                fprintf('Warning2: excessive iterations (%g); l1=%g,l2=%g; break\n',ll,l1,l2); 
                break;
            end
            if (spec(l1,ind_hp)>(val*frac_max))
                ind_hp = ind_hp+1;
            else
                break; 
            end
            ll = ll+1;
        end
        ind_hm = im;         % index to <half maximum
        ll = 0;
        while true           % walk down
            if (ind_hm==1)
                fprintf('Warning2: ind_hm==1; l1=%g,l2=%g; break\n',l1,l2); 
                break; 
            end
            if (ll>si(2)) 
                fprintf('Warning2: excessive iterations (%g); l1=%g,l2=%g; break\n',ll,l1,l2); 
                break;
            end
            if (spec(l1,ind_hm)>(val*frac_max))
                ind_hm = ind_hm-1;
            else
                break; 
            end
            ll = ll+1;
        end
        
        % interpolate for intermediate positions
        % ilp = (im:(ind_hp+1));
        % fp(l1,l2) = interp1(spec(l1,ilp),x_ax(ilp),val/2,'spline');
        if ind_hp==1
            warning('ind_hp==1: adding 1');
            ind_hp = 2;
        end
        if ind_hp==si(2)
            warning('ind_hp==si(2)(=%d): subtracting 1',si(2));
            ind_hp = si(2)-1;
        end
        ilp = ind_hp + (-1:1);
        
        fp(l1,l2) = interp1(spec(l1,ilp),x_ax(ilp),val*frac_max,'linear');
        % ilm = ((ind_hm-1):im);
        % fm(l1,l2) = interp1(spec(l1,ilm),x_ax(ilm),val/2,'spline');
        if ind_hm==1
            warning('ind_hm==1: adding 1');
            ind_hm = 2;
        end
        if ind_hm==si(2)
            warning('ind_hm==si(2)(=%d): subtracting 1',si(2));
            ind_hm = si(2)-1;
        end
        ilm = ind_hm + (-1:1);
        
        
        fm(l1,l2) = interp1(spec(l1,ilm),x_ax(ilm),val*frac_max,'linear');
        lw(l1,l2) = abs(fp(l1,l2)-fm(l1,l2));
        if mag_corr
            % correction for magnitude mode
            % factor of 1/sqrt(3) derived for Lorentzian peaks
            % fits Gaussians approximately as well
            lw(l1,l2) = lw(l1,l2)/sqrt(3);  
        end
        if verb
            fprintf('l1=%g; l2=%g: lw=%g; val_max=%g; ind_max=%g; ind_hp=%g; ind_hm=%g\n',...
                l1,l2,lw(l1,l2),val,im,ind_hp,ind_hm);
        end
    end
    if plt
        if size(ind_max,1)>1, figure(10+l1); clf; end
        plot(x_ax,spec(l1,:),'b-x',...
            freqs(l1,:),spec(l1,ind_max(l1,:)),'g*',...
            fp(l1,:),spec(l1,ind_max(l1,:))*frac_max,'r*',...
            fm(l1,:),spec(l1,ind_max(l1,:))*frac_max,'m*');
        grid on; axis tight;
    end
end

%% populate output structure
if nargout>1
    out.fp = fp; out.fm = fm; out.ind_max = ind_max; out.freqs = freqs;
end

end      % main function fwhm.m
