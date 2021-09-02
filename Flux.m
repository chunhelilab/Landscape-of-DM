clear;
clc;
load firing_rate

time_step = 0.005;
t = 0 : length(firing_rate) * time_step;
a = 1;
[xd, yd, d, pr, flux_x, flux_y, F_data, flux_x1, flux_y1] = flux_2D(firing_rate, t, a);

%% visualize
load U U
U(find(U > 15)) = 15;

x_range = [0, 30];
y_range = [0, 30];
g = 51;
x_step = (x_range(2) - x_range(1)) / (g - 1);
y_step = (y_range(2) - y_range(1)) / (g - 1);

figure(1)
surf(x_range(1) : x_step : x_range(2) , y_range(1) : y_step : y_range(2) , U-20)
axis equal
shading interp
xlabel('r_1 (Sps/s)');
ylabel('r_2 (Sps/s)');
zlabel('U');
set(gca,'FontSize',12)
view([0, 90]);
xlim([x_range(1), x_range(2)])
ylim([y_range(1), y_range(2)])

hold on
q = quiver(xd(1, :)+a/2, yd(:, 1)+a/2, flux_x, flux_y);
q.Color = 'w';
q.LineWidth = 1;
q.AutoScaleFactor = 0.9;