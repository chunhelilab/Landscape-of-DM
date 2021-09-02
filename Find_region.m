%%%% select the region of steady state
clear;
clc;
load U
U(find(U > 15)) = 15;

x_range = [0, 30];
y_range = [0, 30];
g = 51;
x_step = (x_range(2) - x_range(1)) / (g - 1);
y_step = (y_range(2) - y_range(1)) / (g - 1);

U1 = U(1 : (g-1)/2, :);   %%B wins
U2 = U((g-1)/2 + 1 : end, :);   %%A wins

%% center of ellipitic region
mu = zeros(2);
[mu(1, 1), mu(2, 1)] = find(U1 == min(min(U1)));
[mu(1, 2), mu(2, 2)] = find(U2 == min(min(U2)), 1, 'first');
mu(1, 1) = mu(1, 1).* x_step;
mu(2, 1) = mu(2, 1).* y_step;
mu(1, 2) = (mu(1, 2) + (g-1)/2 - 1).* x_step;
mu(2, 2) = mu(2, 2).* y_step;

%% major and minor axies of ellipitic region
sig=[5 5; 2.5 2.5] * 0.3;

figure(3)
surf(x_range(1) : x_step : x_range(2) , y_range(1) : y_step : y_range(2) , U-20)
axis equal
shading interp
xlabel('r_1 (Sps/s)');
ylabel('r_2 (Sps/s)');
zlabel('U');
set(gca,'FontSize',12)
view([0, 90]);
xlim([2, 23])
ylim([2, 23])

%% region1 for DS2
hold on
ecc1 = axes2ecc(sig(1, 1), sig(2, 1)); % eccentricity ratio
[elat1,elon1] = ellipse1(mu(1, 1), mu(2, 1),[sig(1, 1) ecc1], 90);  %90: direction angle
plot(elat1,elon1, '--m', 'LineWidth', 1.5)

%% region2 for DS1
hold on
ecc2 = axes2ecc(sig(1, 2), sig(2, 2)); 
[elat2,elon2] = ellipse1(mu(1, 2), mu(2, 2),[sig(1, 2) ecc2], 0);
plot(elat2,elon2, '--g', 'LineWidth', 1.5)

save region mu sig
