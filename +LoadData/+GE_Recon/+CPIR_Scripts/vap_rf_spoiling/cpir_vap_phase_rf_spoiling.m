% RF phase spoiling incrementation, CPIR https://mriquestions.com/spoiling---what-and-how.html
MAX_PARS = 16382*16; %fid all limits

phase = 0.5 * 117 * ((1:MAX_PARS).^2 + (1:MAX_PARS));
phase = mod(phase,360);

filename = sprintf('vap_phase_rf_spoiling_%d.fdl',MAX_PARS);
write_fdl(filename,phase,'phase');

% function write_fdl(fname,x,what)
%WRITE_FDL  Write vector x into binary file readable for fidall
%  write_fdl(fname,x,what)
%  fname  output filename
%         use meaningful name and link it to 'vap_XXX<cv20>.fdl')
%      x  binary vector
%   what  'te','tr'                                    [us]
%         'phase','flip','phi','theta'                 [deg]
%         'freq'                                       [Hz]
%         'dudri'        [0=I+Q on, 1=Ion+Qoff, 2=Ioff+Qon,3=I+Q off]
%         'agrad'
%         scaling factor
%
%  4/2021 Rolf Schulte