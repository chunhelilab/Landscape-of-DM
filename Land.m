%%%w+ = 1.61, mu_A = mu_B = 58, bistable
%% read firing rate
clear;
clc;
T = 500;  %%total simulation time
stable_time = 0; %%%% stable state (time after stable_time) data is used for landscape
time_step = 0.005;  %%%firing rate is caculated for time step 
time_window = 0.05;  %%%firing rate in the time window is averaged
filenumber = 150;
rate_number = (T - time_window) / time_step;
firing_rate = zeros(rate_number * filenumber, 2);

for idx = 1 : filenumber
    filename = ['firing_rate', num2str(idx-1), '.txt'];
    [S1_rate, S2_rate] = textread(filename, '%n%n','delimiter', ' ','headerlines', 1);  %%%skip the first row
    firing_rate(rate_number * (idx-1) + 1: rate_number * idx, 1) = S1_rate;
    firing_rate(rate_number * (idx-1) + 1: rate_number * idx, 2) = S2_rate;
end
save firing_rate firing_rate

%% the trajectory of firing rate
% figure(1)
% plot_time =  [stable_time : time_step : T - time_window - time_step];
% plot(plot_time, S1_rate)
% hold on
% plot(plot_time, S2_rate)
% legend('S1_{rate}','S2_{rate}')
% xlabel('t(s)')
% ylabel('Spikes/s')


%% calculate landscape 
x_range = [0, 30];
y_range = [0, 30];
g = 51;
x_step = (x_range(2) - x_range(1)) / (g - 1);
y_step = (y_range(2) - y_range(1)) / (g - 1);
for m = 1 : g
    for n = 1:g
        pps(m, n) = sum(firing_rate(:, 1) >= x_range(1) + (m -1) * x_step & firing_rate(:, 1) < x_range(1) + m * x_step & firing_rate(:, 2) >= y_range(1) + (n -1) * y_step & firing_rate(:, 2) < y_range(1) + n * y_step);
    end
end


pps = pps' / sum(sum(pps));   
U = -log(pps);
save U U

%% the landscape of firing rate
load U U
U(find(U > 15)) = 15;

figure(2)
surf(x_range(1) : x_step : x_range(2) , y_range(1) : y_step : y_range(2) , U)
shading interp
xlabel('r_1 (Sps/s)');
ylabel('r_2 (Sps/s)');
zlabel('U');
set(gca, 'FontSize', 12)
view([-29, 59]);
xlim([x_range(1), x_range(2)])
ylim([y_range(1), y_range(2)])

