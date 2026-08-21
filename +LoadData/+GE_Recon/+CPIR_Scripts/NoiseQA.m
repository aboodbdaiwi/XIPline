%% Noise QA
%Script to check for noise

%% Select ScanArchive
[file,location] = uigetfile('*','Select ScanArchive file');
fname = fullfile(location,file);

%% Load ScanArchive
%[d,h,archive] = read_archive(fname,verb,do_single,no_add,do_slice_order)
%     fname   filename                         (default=last h5 one in dir)
%      verb   verbose mode: 0=off;1=simple,2=more,3=all                 (1)
verb = 1;
% do_single   d in single precision                                  (true)
do_single = false;
%    no_add   force control.operation=0                             (false)
no_add = false;
%do_slice_order  perform slice ordering      (true; false for 3drad/burzte)
do_slice_order = true;

[data,header,archive]=read_archive(fname,verb,do_single,no_add,do_slice_order);
% data = [nframes,frame_size,nphases,nechoes,nslices,nreceivers]

clear verb do_single no_add do_slice_order file location

%% Plot Data
noise = squeeze(data);
data_size = size(noise);
noise = shiftdim(noise,1); % make readout first dim
noise = reshape(noise, size(noise,1),[]);
plot_noise = noise;

if size(noise, 2) > 1000 % for faster plotting
    projs = round(linspace(1,size(noise,2),1000));
    plot_noise = plot_noise(:,projs);
end

figure('Name','Noise Data','WindowState','maximized')
t=tiledlayout(3,3);
title(t,['Data Size = ', num2str(data_size)])

% full readout
nexttile
plot(real(plot_noise),'r','LineStyle','none','Marker','.','MarkerSize',5)
title('Real')
nexttile
plot(imag(plot_noise),'g','LineStyle','none','Marker','.','MarkerSize',5)
title('Imaginary')
nexttile
plot(abs(plot_noise),'b','LineStyle','none','Marker','.','MarkerSize',5)
title('Magnitude')

% zoom readout
nexttile
plot(real(plot_noise(1:64,:)),'r','LineStyle','none','Marker','.','MarkerSize',5)
nexttile
plot(imag(plot_noise(1:64,:)),'g','LineStyle','none','Marker','.','MarkerSize',5)
nexttile
plot(abs(plot_noise(1:64,:)),'b','LineStyle','none','Marker','.','MarkerSize',5)

% histogram
nexttile
histogram(real(noise(:)))
nexttile
histogram(imag(noise(:)))
nexttile
histogram(abs(noise(:)))

%% QQ plot
figure('Name','QQ Plots','WindowState','maximized')
t=tiledlayout(1,3);

% full readout
nexttile
qqplot(real(noise(:)))
title('Real')
nexttile
qqplot(imag(noise(:)))
title('Imaginary')
nexttile
pd = makedist("Rician");
qqplot(abs(noise(:)),pd)
title('Magnitude')

