function Dis_Traj = get_allinone_distraj(Dis_Fid)

%% Bunch of hardcoding. Not pretty, but easy.
gamma = 11.777953;
RUT = 100*1e3; %in ns
MaxGrad = 9.475937;
ADC_Dur = 716.8*1e3;%in ns

Resolution = [0.00625 0.00625 0.00625];

%Dwell = 22.4/2;

Dur = ADC_Dur-RUT;
Grad = linspace(0,1,RUT);
Grad((RUT):ADC_Dur) = 1;
Grad = Grad * MaxGrad;

Pts = size(Dis_Fid,1);
ADC_Dur = ADC_Dur*1e-9; %Convert from ns to s
RUT = RUT*1e-9; %Convert from ns to s
Dw = ADC_Dur/Pts;

Arm_untimed = cumtrapz(Grad);
Grad_Time = 0:1e-9:(ADC_Dur-1e-9); %in s

Time = 0:Dw:(ADC_Dur-Dw);
Arm = interp1(Grad_Time,Arm_untimed,Time); %Now we are in mT*s/m

Arm = Arm*gamma/1000;
k_loc = Arm;
NPts = size(Dis_Fid,1);
NPro = size(Dis_Fid,2);

ind = zeros(1,NPro);
for i = 1:NPro
    ind(i) = Halton_rand(i-1,2);
end
[~,newind] = sort(ind);

l_kz = ((2*(newind-1))+1-NPro)/NPro;
l_alph = sqrt(NPro*pi)*asin(l_kz);

traj = zeros(3,NPts,NPro);
for i = 1:NPro
%     traj(1,:,i) = k_loc*sin(Angs(i,2))*cos(Angs(i,1));
%     traj(2,:,i) = k_loc*sin(Angs(i,2))*sin(Angs(i,1));
%     traj(3,:,i) = k_loc*cos(Angs(i,2));
    traj(1,:,i) = k_loc*sqrt(1-(l_kz(i))^2)*cos(l_alph(i));
    traj(2,:,i) = k_loc*sqrt(1-(l_kz(i))^2)*sin(l_alph(i));
    traj(3,:,i) = k_loc*l_kz(i);
end

kFOV_desired = 1./(Resolution);
kMax_desired = kFOV_desired/2;
max_k = max(kMax_desired); %Here, we are in 1/m
Dis_Traj = traj/max_k/2;
Dis_Traj = Dis_Traj/1000;

rad = sqrt(Dis_Traj(1,:,:).^2 + Dis_Traj(2,:,:).^2 + Dis_Traj(3,:,:).^2);
%Dis_Traj = traj/max(rad(:))/2;
hold_traj = Dis_Traj;
hold_traj(1,:,:) = Dis_Traj(2,:,:);
hold_traj(2,:,:) = Dis_Traj(3,:,:);
hold_traj(3,:,:) = Dis_Traj(1,:,:);

Dis_Traj = hold_traj;
