% Calibration, XeCTC and CPIR
% In theory, can just use default spectroscopy fidall2 mode and use custom 
% rf freqs

% Set params
n_diss=500;
n_gas=20;
n_exc=n_diss+n_gas;

% XeCTC
diss_ppm=218;
diss_hz=7696;
freqs = [diss_hz*ones(n_diss,1); zeros(n_gas,1)];

filename = sprintf('cal_xectc_nexc%d_ndiss%d_dfreq%d_freq.fdl',n_exc,n_diss,diss_ppm);
write_fdl(filename,freqs,'freq');

% CPIR
diss_ppm=202;
diss_hz=7143;
freqs = [diss_hz*ones(n_diss,1); zeros(n_gas,1)];

filename = sprintf('cal_cpir_nexc%d_ndiss%d_dfreq%d_freq.fdl',n_exc,n_diss,diss_ppm);
write_fdl(filename,freqs,'freq');