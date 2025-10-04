function traj = philipsradialcoords(delay,trajtype,fullfilename)

%% Function to calculate trajectories from information passed in a Philips sin file
% Largely copied from the PhilipsRadialCoords node in GPI
% This only calculates trajectories for 3D kooshball (i.e. not good for
% spiral, stack of stars, etc... it may be useful to add to this and have
% that functionality in the long run, but keep it simple for now.
%
% Also has functionality for calculating Golden Means trajectories as
% implemented in Peter Niedbalski's addition to the CPIR complete patch

%delay    = trajectory delay required. If not supplied, the value given in
%           the .sin file is used
%trajtype = 0 (default) - calculates Philip's native Interleaved Linear
%                         trajectory order
%         = 1           - calculates Golden Means Trajectories
%         = 2           - calculated Haltoned Spiral Trajectories
%filename = name of the .sin file that contains values of importance
%           If not specified, user is directed to select the file
%=========================================================================

if nargin == 0
    [filename,path]=uigetfile('*.sin','Select sin file');
    trajtype = 0;
    delay = -200; %Set delay to some absurd value so that I know when it has been supplied or not
    fullfilename = fullfile(path,filename);
elseif nargin == 1
    [filename,path]=uigetfile('*.sin','Select sin file');
    fullfilename = fullfile(path,filename);
    trajtype = 0;
elseif nargin == 2
    [filename,path]=uigetfile('*.sin','Select sin file');
    fullfilename = fullfile(path,filename);
end

%Open and read in the text from the sin file
sinFile=fopen(char(fullfilename));
sinRead=textscan(sinFile,'%s','delimiter','\n');
sinRead=sinRead{1};

%% Get lines I care about in this file:
% non_cart_max_encoding_nrs
% non_cart_min_encoding_nrs 
% nr_echoes
% non_cart_fid_delay
% non_cart_fid_slope
% non_cart_fid_smoother
% non_cart_fid_trailing_sample

for index=1:size(sinRead,1)
    testStr = char(sinRead{index});
    if contains(testStr,'non_cart_min_encoding_nrs')
        readline = sinRead{index};
        min_encodes = readline(43:end);       
    end

    if contains(testStr,'non_cart_max_encoding_nrs')
        readline = sinRead{index};
        max_encodes = readline(43:end);
    end

    if contains(testStr,'nr_echoes')
        readline = sinRead{index};
        nr_echoes = readline(43:end);
    end
    if contains(testStr,'non_cart_fid_delay')
        readline = sinRead{index};
        fid_delay = readline(43:end);
    end
    if contains(testStr,'non_cart_fid_slope')
        readline = sinRead{index};
        fid_slope = readline(43:end);
    end
    if contains(testStr,'non_cart_fid_smoother')
        readline = sinRead{index};
        fid_smoother = readline(43:end);
    end
    if contains(testStr,'non_cart_fid_trailing_sample')
        readline = sinRead{index};
        fid_trailing = readline(43:end);
    end
end

max_encodes_num = str2num(max_encodes);
min_encodes_num = str2num(min_encodes);
echoes = str2num(nr_echoes);
if delay == -200
    delay = str2num(fid_delay);
end
slope = str2num(fid_slope);
smoother = str2num(fid_smoother);
ntrail = ceil(abs(str2num(fid_trailing)));

nsamp = max_encodes_num(1) - min_encodes_num(1) + 1;
nprof = max_encodes_num(2)+1;
interleaves = max_encodes_num(3)+1;

%% Calculate the gradient shape:
FILTER_TAPS = 5; % Don't know what this is... set this way in GPI
int_ratio_1 = round(FILTER_TAPS * smoother - delay);
int_ratio_2 = ceil(delay - smoother);
int_ratio_3 = ceil(delay - smoother) + ceil(slope);

start_slope = max([0, (int_ratio_2 + int_ratio_1)]+1);
end_slope = max([0, (int_ratio_3 + int_ratio_1)])+1;

grad = zeros(1,nsamp+ntrail+int_ratio_1);
grad(1:start_slope) = 0;
grad(start_slope:end_slope) = (int_ratio_2:int_ratio_3) - (delay-smoother);

if slope > 0
    grad(end_slope:end) = slope;
else
    grad(end_slope:end) = 1;
end

if smoother > 0
    filt = exp(-(0:(ceil(FILTER_TAPS*smoother)))/smoother);
    if length(filt) > 0
        grad_filt = conv(grad,filt);
    else
        grad_filt = grad;
    end
else
    grad_filt = grad;
end

%% Sanity Check - see if I got a reasonable looking gradient
% figure('Name','Gradient Shape Check')
% plot(1:length(grad_filt),grad_filt)

%% Get the kr location from the gradient shape
kr = cumtrapz(grad_filt(int_ratio_1:(int_ratio_1+nsamp+ntrail)));
norm_fact = 0.5 * ((nsamp+ntrail-1)/((nsamp+ntrail)*kr(nsamp+ntrail)));
kr = kr(1:nsamp)*norm_fact;

%% Let's do another sanity check:
% figure('Name', 'Trajectory Shape Check')
% plot(1:length(kr),kr,'*')
    
%I'll probably have to add a bunch of stuff to deal with multiple echoes...
%for now, just skip it, since we haven't yet done any scans in that fashion

%% Calculate rotation value
phi = 2*pi*((1:interleaves)-1) / interleaves;
phi = phi';
phi = repmat(phi,[1,nprof]);

if trajtype == 0
    z = (2* ((1:nprof) -1) +1 -nprof)/nprof;
    z = repmat(z,[interleaves,1]);
    z(2:2:end,:) = -z(2:2:end,:);
    phi = (phi + sqrt(nprof*pi/interleaves).*asin(z));
elseif trajtype == 1
    phi1 = 0.46557123;
    phi2 = 0.6823278;
    totproj = nprof*interleaves;
    z_a = ((0:(totproj-1))*phi1 - floor((0:(totproj-1))*phi1))*2-1;
    phi_a = ((0:(totproj-1))*phi2 - floor((0:(totproj-1))*phi2))*2*pi;
    
   % z_a1 = reshape(z_a,nprof,interleaves)';
   % phi_a1 = reshape(phi_a,nprof,interleaves)';
    
    z_a1 = reshape(z_a,interleaves,nprof);
    phi_a1 = reshape(phi_a,interleaves,nprof);
    
    golden_angles = mod((0:(nprof-1))*111.246,360);
    [~,sort_index] = sort(golden_angles);
    sorting = 1:nprof;
    New_Indices = sorting(sort_index);
    
    z = z_a1(:,New_Indices);
    phi = phi_a1(:,New_Indices);
%Haltoned Spiral
elseif trajtype == 2
    %calc spiral order
    M_PI = 3.14159265358979323846;
    PrevAngle = 0;
    for proj = 0:nprof-1
        currZ = -1.0 + 2.0 * proj/nprof;
        arch_z(proj+1) = currZ;
        if (proj == 0)
            arch_azi(proj+1) = 0;
        else
            arch_azi(proj+1) = mod(PrevAngle + 3.6/(sqrt(nprof*(1-currZ*currZ))), 2.0*M_PI);
        end
        PrevAngle = arch_azi(proj+1);
    end
    
    %calc halton order
    p1 = 2;
    %p2 = 3;
    for proj = 0:nprof-1
        z = haltonnumber(proj+1,p1) * 2 - 1;
        %phi = 2 * M_PI * haltonnumber(proj+1,p2);
        halt_polar(proj+1) = acos(z);
        %halt_azi(proj+1) = phi; %not used here
    end
    
    %sort spiral via halton
    [~,sort_index] = sort(halt_polar);
    z = arch_z(:,sort_index);
    phi = arch_azi(:,sort_index);    
else
    error('trajtype must be either 0, 1, or 2');
end

cp = cos(phi);
sp = sin(phi);
st = sqrt(1-z.*z);

coords = zeros(interleaves,nprof,nsamp,3);
outdims = [interleaves, nprof, nsamp];
for i = 1:interleaves
    for j = 1:nprof
        coords(i,j,:,1) = cp(i,j)*st(i,j)*kr;
        coords(i,j,:,2) = sp(i,j)*st(i,j)*kr;
        coords(i,j,:,3) = z(i,j)*kr;
    end
end

traj = coords;

end

function result = haltonnumber(index, base)
result = 0;
f = 1.0;
i = index;
while i>0
    f = f/base;
    result = result + f*mod(i,base);
    i = fix(i/base);
end

end
