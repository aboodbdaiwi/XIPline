function [freq,str] = mns_freq_converter(f0,nucleus,molecule,verb)
%MNS_FREQ_CONVERTER  Determine MNS frequencies from 1H centre frequency f0
%    freq = mns_freq_converter(f0,nucleus,molecule,verb)
%      f0   Proton centre frequency                             [Hz]
% nucleus   Desired x-nucleus: 13C, 23Na, 129Xe, 31P
%molecule   Desired molecule (if []->all)                            ([])
%    verb   Verbose mode                                             (true)
%    freq   Determined MNS frequencies                          [Hz]
%
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end
if ~exist('nucleus','var'),  nucleus = []; end
if isempty(nucleus),         nucleus = 'all'; end
if ~exist('molecule','var'), molecule = []; end
if isempty(molecule),        molecule = ''; end
str = sprintf('B0 = %g[T]\n',f0/gyrogamma(1)*2*pi); 
str = sprintf('%s    1H H2O      = %9d [Hz]\n',str,round(f0)); 


%% MNS/1H ratios
% 13C lactate, pyruvate, bicarbonate courtesy to Esben Hansen
% 129Xe courtesy to Bastian Driehuys multi-centre study guid
% nucleus molecule ratio validated
mns_tab = {...
    {'13C','lactate',1/3.976224812,1},...
    {'13C','pyruvate',1/3.976271443,1},...
    {'13C','bicarb',1/3.976311538,1},...
    {'13C','urea',0.251419,0},...
    {'13C','oil',0.25145,0},...
    {'13C','peg',0.251465,0},...
    {'23Na','',0.264506,0},...
    {'129Xe','air',  0.2766032749,0},...
    {'129Xe','lung', 0.2766029834,0},...        % 1ppm below air
    {'129Xe','blood',0.2766632951,0},...        % 218ppm above lung
    {'129Xe','paren',0.2766577631,0},...        % +198ppm; lung parenchyma
    {'129Xe','mid',0.2766605291,0},...          % 208ppm: middle dissolved
    {'31P','',0.404792,0}};
nmt = length(mns_tab);


%% actual frequency conversion
val_str{1} = 'not validated';
val_str{2} = '';
if ~isempty(regexpi(nucleus,'all')) || isempty(nucleus)
    freq = zeros(1,nmt);
    for lmt=1:nmt
        f0mns = round(f0*mns_tab{lmt}{3});
        freq(1,lmt) = f0mns;
        str = sprintf('%s%6s %-8s = %9d [Hz] \t %s\n',str,...
            mns_tab{lmt}{1},mns_tab{lmt}{2},f0mns,val_str{mns_tab{lmt}{4}+1});
    end
else
    freq = [];
    for lmt=1:nmt
        if ~isempty(regexpi(mns_tab{lmt}{1},nucleus))
            if ~isempty(regexpi(mns_tab{lmt}{2},molecule)) || isempty(molecule)
                f0mns = round(f0*mns_tab{lmt}{3});
                freq = [freq f0mns];
                str = sprintf('%s%6s %-8s = %9d [Hz] \t %s\n',str,...
                    mns_tab{lmt}{1},mns_tab{lmt}{2},f0mns,val_str{mns_tab{lmt}{4}+1});
            end
        end
    end
end


%% print output
if verb
    fprintf(str);
end

end      % main function mns_freq_converter.m
