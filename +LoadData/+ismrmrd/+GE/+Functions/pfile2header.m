function h = pfile2header(hdr)
%PFILE2HEADER Convert pfile header structure from Orchestra to original
%   h = pfile2header(hdr)
% hdr   Header structure as returned by GERecon
%       hdr=GERecon('Pfile.Header',pfile);
%
%   h   Header structure as returned by p-file matlab reading routines
%
% See also GERecon, read_MR_headers.
% 1/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% list for renaming field names
rnm = {{'RawHeader','rdb_hdr'},...
    {'DataAcqTable','data_acq_tab'},...
    {'PrescanHeader','psc'},...
    {'ExamData','exam'},...
    {'SeriesData','series'},...
    {'ImageData','image'},...
    {'GradData','grad_data'},...
    {'CttEntry','cttentry'}};


%% actual assignment to new field names
for l=1:length(rnm)
    if isfield(hdr,rnm{l}{1})
        h.(rnm{l}{2}) = hdr.(rnm{l}{1});
    else
        warning('hdr.%s field not existing',rnm{l}{1});
    end
end

%% remove whitespace from psd name
if isfield(h.image,'psdname')
    h.image.psdname = deblank(h.image.psdname);
end
if isfield(h.image,'psd_iname')
    h.image.psd_iname = deblank(h.image.psd_iname);
end



end  % main function pfile2header
