function [T1,T2] = nist_phantom(model,plane,B0,verb)
%NIST_PHANTOM T1 and T2 values of NIST-ISMRM QalibreMD phantoms
%[T1,T2] = nist_phantom(model,plane,B0,verb) or tt = nist...
%  model   1=model 130, 2=essential Rev A, 3=TO5                        (1)
%  plane   1=T1 plane, 2=T2 plane, no values for PD plane               (2)
%     B0   Field strength (currently only 3T)                           (3)
%   verb   Verbose mode                                              (true)
%     T1   T1 relaxation times [ms]
%     T2   T2 relaxation times [ms]
%     tt   Relaxation times [2,#vials] [ms]
% Literature: 
% model130 SN>=0042: nist_ismrm_SystemPhantomManual_RevL.pdf p24
% essential: Essential-System-Phantom-Manual_Rev-A-compressed.pdf p21+p
% 11/2022 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% default input
if ~exist('model','var'), model = []; end
if isempty(model),        model = 1; end
if ~exist('plane','var'), plane = []; end
if isempty(plane),        plane = 2; end
if ~exist('B0','var'),    B0 = []; end
if isempty(B0),           B0 = 3; end
if B0~=3, error('only 3T implemeted'); end
if ~exist('verb','var'),  verb = []; end
if isempty(verb),         verb = true; end


%% table with values
switch model
    case 1    % model130 SN>=0042: nist_ismrm_SystemPhantomManual_RevL.pdf p24
        switch plane
            case 1
                T1 = [1838 1398 998.3 725.8 509.1 367 258.7 184.7 130.8 90.9 ...
                    64.2 46.28 32.65 22.95];
                T2 = [1354 1035 728.3 524.4 368.6 266.7 189.3 134.1 93.8 65.7 ...
                    46.8 33.11 23.69 16.73];
            case 2
                T1 = [2756 2281 1961 1552 1341 1017 782.1 589.7 443.8 299.8 ...
                    237.8 170.5 121.8 86.9];
                T2 = [645.8 423.6 286 184.8 134.1 94.4 62.51 44.98 30.95 20.1...
                    15.4 10.85 7.59 5.35];
            otherwise, error('plane not found');
        end
    case 2    % essential: Essential-System-Phantom-Manual_Rev-A-compressed.pdf p21+p
        % values for 3T 20°C
        switch plane
            case 1  % p21
                T1 = [1883.97 1330.16 987.27 690.08 484.97 ...
                    341.58    240.86  174.95 121.08 85.75 ...
                    60.21     42.89   30.40  21.44];
                T2 = [1489.41 1026.78 736.97 521.27 354.38 ...
                    248.03    174.99  126.33 86.91  60.86 ...
                    42.44     30.66   21.76  15.30];
            case 2 % p25
                T1 = [2478.19 2185.54 1901.28 1549.98 1197.57 ...
                    1026.43   805.1   599.96  431.22  292.87 ...
                    226.53    158.17  116.7   82.46];
                T2 = [552.73  379.48  267.29  175.05  112.66 ...
                    88.89     63.42   44.24   29.88   19.40 ...
                    14.74     10.52   7.27    5.10];
            otherwise, error('plane not found');
        end
    case 3    % Eurospin TO5; leedstestobjects.com
        % temperature=292K
        T1a = [200 299 296 463 450 444 604 596 754 745 ...
            903 1448 966 1034 1160 1276 1262 1415];
        T2a = [52 73 113 53 94 154 95 136 116 157 ...
            137 390 224 167 214 204 184 174];
        % temperature=296K
        T1b = [223 334 331 516 502 496 674 666 841 831 ...
            1007 1615 1078 1153 1293 1422 1407 1576];
        T2b = [50 70 111 49 89 151 89 129 109 149 ...
            128 373 212 156 201 190 171 161];
        T1 = round((T1a+T1b)/2);  % mean is T=21C
        T2 = round((T2a+T2b)/2);
    otherwise, error('model (=%s) not found');
end


%% print info
if verb
    switch model
        case 1, fprintf('model 130: ');
        case 2, fprintf('essential Rev A: ');
        case 3, fprintf('TO5\n');
        otherwise, warning('model unknown');
    end
    if model~=3, fprintf('T%d plane\n',plane); end
    nn = 1:length(T1);
    fprintf('vial  ');
    for ll=nn, fprintf('%4d ',ll); end
    fprintf('\nT1    ');
    for ll=nn, fprintf('%4d ',round(T1(ll))); end
    fprintf('\nT2    ');
    for ll=nn, fprintf('%4d ',round(T2(ll))); end
    fprintf('\nT1/T2 ');
    for ll=nn, fprintf('%.3g ',T1(ll)./T2(ll)); end
    fprintf('\n');
end


%% adapt output arguments
switch nargout
    case 0, clear T1 T2
    case 1, T1 = [T1;T2];
    case 2
    otherwise
        error('nargout(=%d) ~= 0,1,2',nargout);
end


end      % nist_phantom.m
