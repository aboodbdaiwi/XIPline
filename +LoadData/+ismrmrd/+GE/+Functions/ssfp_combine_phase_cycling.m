function bc = ssfp_combine_phase_cycling(b0,b1,h0,h1,te,t1,t2,fname,f0add)
%SSFP_COMBINE_PHASE_CYCLING  Combine phase-cycled bSSFP images using 
%   simulated frequeny response
%   bc = ssfp_combine_phase_cycling(b0,b1,h0,h1,te,t1,t2,fname)
%   b0   Images acquired with different TEs (x,y,z,te)
%   b1   bSSFP images acquired with phase-cycling (x,y,z,npc)
%        b0/b1 can be mat filename with reconstructed data
%   h0   Header structure
%   h1   Header structure
%   te   Echo times                                   [s]  ([0 2 5]*1d-3)
%   t1   T1 relaxation time                           [s]  (490d-3)
%   t2   T2 relaxation time                           [s]  (160d-3)
%fname   Save results to fname (png+mat+dicom)
%f0add   Add f0 to df0 map                            [Hz] (0)
%
%   bc   Combined bSSFP images
%
% 12/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('te','var'), te = []; end
if isempty(te),        te = [0 2 5]*1d-3; end
if ~exist('t1','var'), t1 = []; end
if isempty(t1),        t1 = 490d-3; end
if ~exist('t2','var'), t2 = []; end
if isempty(t2),        t2 = 160d-3; end
if ~exist('fname','var'), fname = []; end
if ~exist('f0add','var'), f0add = []; end
if isempty(f0add),        f0add = 0; end


%% load if mat filename given
if ischar(b0)
    if ~exist('b0','file')
        x0 = load(b0);
        b0 = x0.bb;
        h0 = x0.h;
    else
        warning('file(=%s) not found',b0);
    end
end
if ischar(b1)
    if ~exist('b1','file')
        x1 = load(b1);
        b1 = x1.bb;
        h1 = x1.h;
    else
        warning('file(=%s) not found',b1);
    end
end


%% checks
sub_check_header(h0,h1);     % ensure same geometry
te = te(:).';
if max(te)>0.1, warning('max(te)(=%g)>0.1',max(te)); end
if size(b0,4)~=length(te)
    warning('size(b0,4)(=%d)~=length(te)(=%d)',size(b0,4),length(te))
end
n0 = size(b0);
n1 = size(b1);
if length(n0)~=4, warning('b0 not 4D; dim=%d',length(n0)); end
if length(n1)~=4, warning('b1 not 4D; dim=%d',length(n1)); end


%% calculate B0 map
fov0 = h0.rdb_hdr.fov;
fov1 = h1.rdb_hdr.fov;
if abs(fov0-fov1)>0.1
    fprintf('Different FOV (%g;%g[mm]); interpolating b0\n',fov0,fov1);
    nn = [round(n0(1:3)*fov1/fov0) n0(4)];
    b0 = truma(b0,false,nn);
end
if any(n0(1:3)~=n1(1:3))
    fprintf('Different matrix sizes; interpolating b0\n');
    b0 = truma(b0,true,n1(1:3));
end
df0 = calc_b0map(b0,te);
if h0.image.specnuc~=h1.image.specnuc
    fprintf('Different nuclei (%d;%d)\n',h0.image.specnuc,h1.image.specnuc);
    grat = gyrogamma(h0.image.specnuc)/gyrogamma(h1.image.specnuc);
    fprintf('\tscaling df0 map by gamma ratio = %g\n',grat);
    df0 = df0/grat;
end
if abs(f0add)>1d-10
    fprintf('Adding %g [Hz] to df0\n',f0add);
    df0 = df0 + f0add;
end


%% calculate bSSFP frequency response
npc = size(b1,4);
if npc~=4, warning('npc(=%d)~=4',npc); end
flip = h1.image.mr_flip; 
tr = h1.image.tr*1d-6;
ww = ssfp_freq_resp(flip,tr,t1,t2,df0,npc);
bc = sum(b1.*ww,4);


%% plotting
figure(30); clf
imagesc_ind3d(abs(bc));


%% saving
if ~isempty(fname)
    print([fname '.png'],'-dpng','-r600','-painters');
    save([fname '.mat'],'bc','df0','ww',...
        'b0','b1','h0','h1','te','t1','t2','npc','flip','tr');
    inp.SeriesNumber = 100*h1.series.se_no;
    write_dicom(fname,bc,h1,inp);
end


end      % ssfp_combine_phase_cycling.m


%% sub-functions
% make sure geometrical prescription is identical
function sub_check_header(h0,h1)
fld = {...
    {'series','start_loc'},...
    {'series','end_loc'},...
    {'rdb_hdr','fov'},...
    {'image','slthick'},...
    {'image','scanspacing'},...
    {'rdb_hdr','user26'},...
    {'rdb_hdr','user27'},...
    {'rdb_hdr','user28'},...
    {'data_acq_tab','rotate'},...
    {'data_acq_tab','transpose'}};

for l=1:length(fld)
    if isfield(h0,fld{l}{1}) && isfield(h1,fld{l}{1})
        if isfield(h0.(fld{l}{1}),fld{l}{2}) && isfield(h1.(fld{l}{1}),fld{l}{2})
            v0 = h0.(fld{l}{1}).(fld{l}{2});
            v1 = h1.(fld{l}{1}).(fld{l}{2});
            % fprintf('%g %g\n',v0,v1);
            if abs(v0-v1)>1d-10
                warning('h0.%s.%s(=%g)~=h1.%s.%s(=%g)',...
                    (fld{l}{1}),(fld{l}{2}),v0,...
                    (fld{l}{1}),(fld{l}{2}),v1);
                fprintf('Function requires same geometry to work!!!\n');
            end
        else
            warning('field(=%s) not existing',fld{l}{2});
        end
    else
        warning('field(=%s) not existing',fld{l}{1});
    end
end
end      % sub_check_header
