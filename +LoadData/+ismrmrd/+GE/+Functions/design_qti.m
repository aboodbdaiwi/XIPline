function design_qti(nrep,fname,fa,tr,ti,tmax,trec)
% DESIGN_QTI Design QTI flip angle scheme for MRP and add to ak_grad mat
% design_qti(nrep,fname,fa,tr,ti,tmax,trec)
%  nrep   #exc per QTI cycle
% fname   Waveform file name
%    fa   Parameters for QTI flip angle modulation
%         
%    tr
%    ti
%  tmax
%  trec
%
%  2/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
if ~exist('fa','var'),   fa = []; end
if isempty(fa),          fa = 1.5; end
if length(fa)<2,         fa = [fa 70]; end
if length(fa)<3,         fa = [fa 0.8]; end
if ~exist('tr','var'),   tr = []; end
if isempty(tr),          tr = 10; end
if ~exist('ti','var'),   ti = []; end
if isempty(ti),          ti = 18; end
if ~exist('tmax','var'), tmax = []; end
if isempty(tmax),        tmax = 4.2d3; end
if ~exist('trec','var'), trec = []; end
if isempty(trec),        trec = 7d3; end


%% QTI flip angle scheme
tall = nrep*tr;
tt = (1:nrep)*tr;
if tall<8d3,  warning('tall(=%g[ms])<8000[ms]',tall); end
if tmax>tall, warning('tmax(=%g[ms])>tall(=%g[ms])',tmax,tall); end
if trec>tall, warning('trec(=%g[ms])>tall(=%g[ms])',trec,tall); end
if tmax>trec, warning('tmax(=%g[ms])>trec(=%g[ms])',tmax,trec); end
[~,imax] = min(abs(tt-tmax));
[~,irec] = min(abs(tt-trec));
flip_arr = [linspace(fa(1),fa(2),imax),linspace(fa(2),fa(3),irec-imax),...
    fa(3)*ones(1,nrep-irec)];
if size(flip_arr,2)~=nrep
    warning('size(flip_arr,2)(=%d)~=nrep(=%d)',size(flip_arr,2),nrep); 
end


%% plotting
clf;
plot(tt,flip_arr); 
title('QTI Flip Angle Scheme');
xlabel('time [ms]'); ylabel('Flip angle [deg]');
axis([0 tall 0 ceil(fa(2))]); grid on



%% loading waveform file
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
    load([fname '.mat']);
    if exist('method','var')
        warning('method existing; overwriting');
        clear method
    end
end


%% method structure
method.ExcPulse = {2};
method.InvPulse = {10,'adiabatic'};
method.FISP_inversion_enable = 'Yes';
method.InversionTime = ti;
method.VariableTR = zeros(1,nrep);
method.VariableFlip = flip_arr;
method.VariableTE = zeros(1,nrep);
method.VariablePhase = zeros(1,nrep);
method.Method = '';
method.threeD = 1;


%% actual saving of method file
if ~isempty(fname)
    write_fdl([fname '_flip.fdl'],flip_arr,'flip');
    clear tall nrep tt tr imax irec fa flip_arr te tmax trec ti
    save([fname '.mat']);    % add method structure to waveform file
end


end      % design_qti.m
