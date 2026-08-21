function [freq,str] = mns_freq_converter(f0,nucleus,molecule,verb)
%MNS_FREQ_CONVERTER  Determine MNS frequencies from 1H centre frequency f0
%    freq = mns_freq_converter(f0,nucleus,molecule,verb)
%      f0   Proton centre frequency                             [Hz]
% nucleus   Desired x-nucleus: 13C, 23Na, 129Xe, 31P
%molecule   Desired molecule (if []->all)                            ([])
%    verb   Verbose mode                                             (true)
%    freq   Determined MNS frequencies                          [Hz]
%
%  8/2022 Rolf Schulte
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
    {'2H','D2O',      0.153506071852142, 653.62},...
    {'13C','lactate', 1/3.976224812,    1070.63},...
    {'13C','pyruvate',1/3.976271443,    1070.63},...
    {'13C','bicarb',  1/3.976311538,    1070.63},...
    {'13C','urea',    1/3.976302499380294,1070.63},...
    {'13C','oil',     0.25145,          1070.63},...
    {'13C','peg',     0.251465,         1070.63},...
    {'23Na','',       0.264517727834842,1126.09},...
    {'31P','Ph buffer', 0.404792,       1723.51},...
    {'129Xe','gas bag', 0.2766032749,  -1177.60},...
    {'129Xe','gas lung',0.2766029834,  -1177.60},...  % 1ppm below air
    {'129Xe','blood', 0.2766632951,    -1177.60},...  % 218ppm above lung
    {'129Xe','tissue',0.2766577631,    -1177.60},...  % +198ppm; lung parenchyma
    {'129Xe','mid',  0.2766605291,     -1177.60}};    % 208ppm: middle dissolved
nmt = length(mns_tab);
GAM = 4257.59;                                  % GE 1H gyromagnetic ratio 


%% actual frequency conversion
if ~isempty(regexpi(nucleus,'all')) || isempty(nucleus)
    freq = zeros(1,nmt);
    for lmt=1:nmt
        f0mns = round(f0*mns_tab{lmt}{3});
        freq(1,lmt) = f0mns;
        nu = (1 - mns_tab{lmt}{3}*GAM/abs(mns_tab{lmt}{4}))*1d6;
        str = sprintf('%s%6s %-8s = %9d [Hz] \t %8g [ppm]\n',str,...
            mns_tab{lmt}{1},mns_tab{lmt}{2},f0mns,nu);
    end
else
    freq = [];
    for lmt=1:nmt
        if ~isempty(regexpi(mns_tab{lmt}{1},nucleus))
            if ~isempty(regexpi(mns_tab{lmt}{2},molecule)) || isempty(molecule)
                f0mns = round(f0*mns_tab{lmt}{3});
                freq = [freq f0mns];
                nu = (1 - mns_tab{lmt}{3}*GAM/abs(mns_tab{lmt}{4}))*1d6;
                str = sprintf('%s%6s %-8s = %9d [Hz] \t %8g [ppm]\n',str,...
                    mns_tab{lmt}{1},mns_tab{lmt}{2},f0mns,nu);
            end
        end
    end
end


%% print output
if verb
    fprintf(str);
end

end      % main function mns_freq_converter.m
