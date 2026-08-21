function gcs = get_gradcoil(gcn)
%GET_GRADCOIL Translate gradient coil number into string
% gcs = get_gradcoil(gcn)
%  gcn  Gradient Coil number
%       or header with mrconfig field (from ScanArchive)
%  gcs  Gradient Coil string
%
%  8/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% extract gradcoil number from h.mrconfig.GradCoil
gcs = '';
if isstruct(gcn)
    h = gcn;
    if isfield(h,'mrconfig')
        if isfield(h.mrconfig,'GradCoil')
            gcn = h.mrconfig.GradCoil;
        else
            warning('No field h.mrconfig.GradCoil');
            return;
        end
    else
        warning('No field h.mrconfig');
        return;
    end
end


%% lookup list from MrConfiguration.py
switch gcn
    case 2, gcs = 'BRM';
    case 3 ,gcs = 'TRM';
    case 5, gcs = 'miniCRM';
    case 6, gcs = 'CRM';
    case 7, gcs = 'ORM';
    case 8, gcs = 'XRM';
    case 9, gcs = 'MFO';
    case 10, gcs = 'MFOCPG';
    case 11, gcs = 'XRMW';
    case 12, gcs = 'VRMW';
    case 13, gcs = 'HRMW';
    case 15, gcs = 'ESPRM';
    case 16, gcs = 'HRMB';
    case 17, gcs = 'IRMW';
    case 19, gcs = 'JXG';
    otherwise, warning('gradcoil num=%d unknown',gcn);
end

end      % get_gradcoil.m
