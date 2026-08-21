%% Info - update
% traj_mat_file='vent_cpir_floret\floret_h1_fov375_mtx97_arms12p6_hub3_0p4_intlv2304_kdt2_gmax21_smax148_dur1p6_M8.mat';
% traj_wav_file='vent_cpir_floret\floret_h1_fov375_mtx97_arms12p6_hub3_0p4_intlv2304_kdt2_gmax21_smax148_dur1p6_M8.wav';
traj_mat_file='interleaved_radial3D_129xe_fov400_mtx64_intlv2002_kdt10_gmax31_smax106_dur1p2_coca_rew2_1p25.mat';
traj_wav_file='interleaved_radial3D_129xe_fov400_mtx64_intlv2002_kdt10_gmax31_smax106_dur1p2_coca_rew2_1p25.wav';

%% Read/Analyze
wf = load_waveform(traj_wav_file);
[grad,bw,fov,desc,N,params,n_kpts] = read_ak_wav(traj_wav_file);
[snr,lw,sz,dur,mdcf,psf_vol] = analyse_trajectory(traj_mat_file);
[PThresh,pt,PTmax,gmax,smax,t,f] = pns(grad,'HRMW');

%% Plot Grads
figure();
fig=tiledlayout(3,1);
t_grad = (0:N.gpts-1)*params(6)*1d-6*1d3;

nexttile
hold on
yline(0,'k--');
for p=1:size(grad,2)
    plot(t_grad,grad(:,p,1)*1d3)
end
ylabel('X grad')
hold off

nexttile
hold on
yline(0,'k--');
for p=1:size(grad,2)
    plot(t_grad,grad(:,p,2)*1d3)
end
ylabel('Y grad')
hold off

nexttile
hold on
yline(0,'k--');
if size(grad,3) == 3    
    for p=1:size(grad,2)
        plot(t_grad,grad(:,p,3)*1d3)
    end
else
    plot(t_grad,zeros(N.gpts,1))
end
ylabel('Z grad')
xlabel('Time (ms)','FontSize',12)
hold off

title(fig, traj_mat_file(1:end-4),'Interpreter','none','FontSize',10);
ylabel(fig, 'Grad (mT/m)','FontSize',12);

colororder(turbo(max(N.intl)));

%% Plot Coords
num_trajs = size(wf.ind,1);
wf_read = size(wf.k,2)/num_trajs;
ktraj = reshape(wf.k,wf_read,num_trajs,3);
num_plot = 25;
if size(ktraj,2) < num_plot
    num_plot = size(ktraj,2);
end
figure();
hold on
for p=1:num_plot
    plot3(ktraj(:,p,1),ktraj(:,p,2),ktraj(:,p,3),'LineWidth',1)
end
hold off
axis equal
xlabel('kx')
ylabel('ky')
zlabel('kz')
colororder(turbo(max(num_plot)));