% change
k_plot = k2;
ro_length = size(ind2,2);
points = 1000000;

%run
k_plot = permute(k_plot,[3 2 1]);

figure()
x = k_plot(1,1:points);
y = k_plot(2,1:points);
z = k_plot(3,1:points);
c = 1:length(z);
surface([x;x], [y;y], [z;z], [c;c],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
axis equal
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])
colormap(jet)

figure()
x = k_plot(1,1:ro_length:end);
y = k_plot(2,1:ro_length:end);
z = k_plot(3,1:ro_length:end);
c = 1:length(z);
surface([x;x], [y;y], [z;z], [c;c],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 2);
axis equal
xlim([-0.5 0.5])
ylim([-0.5 0.5])
zlim([-0.5 0.5])
colormap(jet)