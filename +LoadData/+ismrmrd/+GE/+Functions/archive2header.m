function h = archive2header(arc,verb,rdcsz)
%ARCHIVE2HEADER Convert archive to pfile header structure
%    h = archive2header(arc)
%  arc   Archive structure as returned by GERecon
%        archive=GERecon('Archive.Load','XXX.h5');
% verb   Verbose mode                                           (false)
%rdcsz   Reduce header size by removing not uncommon fields     (true)
%
%    h   Header structure as returned by p-file matlab reading routines
%
% See also GERecon, read_MR_headers.
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default input
if ~exist('verb','var'),  verb = []; end
if isempty(verb),         verb = false; end
if ~exist('rdcsz','var'), rdcsz = []; end
if isempty(rdcsz),        rdcsz = true; end


%% conversion
h = arc.DownloadData;
h = subfun_remove_trailing_str(h,'rdb_hdr_',verb);
h.rdb_hdr = h.rec; 
h = rmfield(h,'rec');
h.rdb_hdr.dab = [h.rdb_hdr.dab(1).start_rcv h.rdb_hdr.dab(1).stop_rcv];


%% remove not used fields
if rdcsz
    flds = {'per_pass','unlock_raw','nex_tab','nex_abort_tab','tool'};
    
    for ll=1:length(flds)
        if isfield(h,flds{ll})
            h = rmfield(h,flds{ll});
        end
    end
    h.data_acq_tab = h.data_acq_tab(1);
end


end      % archive2header.m


%% sub-functions
function h = subfun_remove_trailing_str(h,str,verb)

fn = fieldnames(h);
if verb, fprintf('length(fn)=%g\n',length(fn)); end
for l=1:length(fn)
    xx = h.(fn{l});
    if size(xx,1)>1
        if verb, fprintf('Attention: size(xx,1)(=%g)>1\n',size(xx,1)); end
    end
    if size(xx,2)>1
        if verb, fprintf('Attention: size(xx,2)(=%g)>1\n',size(xx,2)); end
    end
    if isstruct(xx) && (size(xx,2)==1)
        h.(fn{l}) = subfun_remove_trailing_str(xx,str,verb);
    end
    if regexpi(fn{l},str)==1
        nf = fn{l}((length(str)+1):end);
        if ~isempty(nf)
            if ~isempty(regexp(nf(1),'\d','once'))
                if verb
                    fprintf('Attention: field (''%s'') starting with number: adding X\n',nf);
                end
                nf = ['X' nf];
            end
            h.(nf) = h.(fn{l});
            h = rmfield(h,fn{l});
        else
            if verb
                fprintf('Attention: new field name empty: ''%s'' -> ''%s''\n',...
                    fn{l},nf);
            end
        end
    else
        if verb
            fprintf('Attention: cannot remove ''%s'' from ''%s''\n',str,fn{l});
        end
    end
end

end      % subfun_remove_trailing_str
