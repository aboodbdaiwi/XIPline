function spec = csi_orientation(spec,h)
%CSI_ORIENTATION  Orient images (xy plane) from xyz to AP-RL-SI coordinates
%
% spec = csi_orientation(spec,h)
%
%     spec  Spectrum        [#s,#x,#y,#z,#t,#c]
%        h  p-file header 
%           h.data_acq_tab.rotate
%           h.data_acq_tab.transpose
%
% 11/2015 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% misc parameters + checks
rot90fac = h.data_acq_tab.rotate;
trnsps =   h.data_acq_tab.transpose;
if any(diff(rot90fac)), warning('rot90fac'); disp(rot90fac); end
if any(diff(trnsps)),   warning('trnsps');   disp(tsnsps); end
rot90fac = rot90fac(1);
trnsps = trnsps(1);
nn = size(spec); nn = [nn ones(1,6-length(nn))];


%% determine which operation to follow
% values found empirically by comparing to rot90/fliplr in recon_grid3d
ori = 0;                % default=do nothing
switch trnsps
    case 0,
        switch rot90fac
            case 0, ori = 1; %2;
            case 1, ori = 7; %4;
            case 2, ori = 2; %1;
            case 3, ori = 4; %7;
            otherwise, warning('h.data_acq_tab.rotate(=%g) not 0,1,2 or 3',rot90fac)
        end
    case 3,
        switch rot90fac
            case 0, ori = 5; %6;
            case 1, ori = 0; %3;
            case 2, ori = 6; %5;
            case 3, ori = 3; %0;
            otherwise, warning('h.data_acq_tab.rotate(=%g) not 0,1,2 or 3',rot90fac)
        end
    otherwise,
        warning('h.data_acq_tab.transpose(=%g) not 0 or 3: unknown image orientation',...
            trnsps(1));
end

%% perform actual orientation by looping through x and y indices
if ori>0,
    if ori>3,
        tmp = zeros(nn(1),nn(3),nn(2),nn(4),nn(5),nn(6));
    else
        tmp = zeros(nn);
    end

    for l2=1:nn(2),
        for l3=1:nn(3),
            switch ori
                case 0,%ll2 = l2;         ll3 = l3;
                case 1, ll2 = nn(2)-l2+1; ll3 = l3;
                case 2, ll2 = l2;         ll3 = nn(3)-l3+1;
                case 3, ll2 = nn(2)-l2+1; ll3 = nn(3)-l3+1;
                case 4, ll3 = l2;         ll2 = l3;
                case 5, ll3 = nn(2)-l2+1; ll2 = l3;
                case 6, ll3 = l2;         ll2 = nn(3)-l3+1;
                case 7, ll3 = nn(2)-l2+1; ll2 = nn(3)-l3+1;
            end
            tmp(:,ll2,ll3,:,:,:) = spec(:,l2,l3,:,:,:);
        end
    end
    spec = tmp;
end

end

