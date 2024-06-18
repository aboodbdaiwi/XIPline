function fidall_b1shimming(b1file,cenX,cenY,radius,fac_thresh,method)
% FIDALL_B1SHIMMING Wrapper function to evaluate optimimal dual drive 
%   settings using B1maps acquired with fidall
%fidall_b1shimming(b1file,cenX,cenY,radius,fac_thresh,method)
%    b1file  Filename of fidall-BS matlab file
%      cenX  Centre of mask in x                      [pixel]
%      cenY  Centre of mask in y                      [pixel]
%    radius  Radius of circular mask                  [pixel]
%            skip masking if 0
%fac_thresh  Factor of threshold for segmenting mask  (default=1)
%    method  b1optimiser method (see help there)      (default=2)

%     b1fit  Resulting b1map
%   txAtten  Amplitude difference between Q and I channel Tx  [dB/10]
%            txAtten = -200*log10(x(1));
%            Prescan UI: Amplitude Attenuation
%   txPhase  Phase difference between Q and I channel Tx      [deg]
%            txPhase = x(2)*180/pi; 
%            Prescan UI: Phase Delay
%
% See also B1OPTIMISER, BLOSE_B1MAP and SEGMENT_SPHERE by Rolf Schulte
% 9/2016 Guido Kudielka, modified RFS

if (nargin<1), help(mfilename); return; end;

if ~exist('cenX','var'),       cenX = []; end
if ~exist('cenY','var'),       cenY = []; end
if ~exist('radius','var'),     radius = []; end
if ~exist('fac_thresh','var'), fac_thresh = []; end
if isempty(fac_thresh),        fac_thresh = 1; end
if ~exist('method','var'),     method = []; end
if isempty(method),            method = 2; end

if ~exist('b1file','var'), error('filename for b1 data is missing'); end
if ~exist(b1file,'file'),  warning('file (%s) not found',b1file); end
verb = 2; 

%% load data and generate b1 map
b1data = load(b1file);
b1map = blosi_b1map(b1data.bb,b1data.h,1,0,verb);
bbabs = sqrt(mean(mean(conj(b1data.bb).*b1data.bb,6),3));
figure(20); clf; imagesc(bbabs); drawnow;
if isempty(cenX) || isempty(cenY),
    fprintf('Enter centre of circular mask\n');
    cenX = input('x=');
    cenY = input('y=');
end
if isempty(radius),
    fprintf('Enter radius of circular mask\n');
    radius = input('radius=');
end



%% generate mask
if radius>0
    mask = segment_sphere(b1map(:,:,1),[cenY cenX],radius,fac_thresh);
else
    mask = [];
end


%% evaluate optimised b1
figure(21); clf;
b1optimiser(b1map,mask,method,[],verb);

end  % end of main function (fidall_b1shimming.m)
