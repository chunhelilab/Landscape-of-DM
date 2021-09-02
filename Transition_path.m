clear;
load firing_rate
load region
%% find the transition way
%Use the sigma principle to specify areas of two states
%Two elliptical regions
num = zeros(1, size(firing_rate, 1));

num((firing_rate(:, 1) - mu(1, 1)).^2 / (sig(1, 1).^2) + ...
    (firing_rate(:, 2) - mu(2, 1)).^2 / (sig(2, 1).^2) <= 1) = 1;  %1: region1, high r2, low r1
num((firing_rate(:, 1) - mu(1, 2)).^2 / (sig(1, 2).^2) + ...
    (firing_rate(:, 2) - mu(2, 2)).^2 / (sig(2, 2).^2) <= 1) = 3;   %2: region2, high r1, low r2

index(1, :) = find((num==1) | (num==3));     %the index of region1 and region2
index(2, :) = num(find((num==1) | (num==3)));      %the value of region1 and region2

transition1 = find((index(2, 1:end-1) - index(2, 2:end) == -2));  %find the transition point: 1¡ú2
transition2 = find((index(2, 1:end-1) - index(2, 2:end) == 2));  %find the transition point: 2¡ú1

way_1to2 = zeros(length(transition1), 2);
way_2to1 = zeros(length(transition2), 2);

way_1to2(:, 1) = index(1, [transition1]);  %every row is a path from 1 to 2; 1st column is the start point in region1, 2st column is the end point in region2
way_1to2(:, 2) = index(1, [transition1 + 1]);
way_2to1(:, 1) = index(1, [transition2]);  %every row is a path from 2 to 1; 1st column is the start point in region2, 2st column is the end point in region1
way_2to1(:, 2) = index(1, [transition2 + 1]);

%% plot single way
figure(1)
for i = 1 : size(way_1to2, 1)
    plot(firing_rate(way_1to2(i, 1):way_1to2(i, 2), 1), ...
        firing_rate(way_1to2(i, 1):way_1to2(i, 2), 2), 'm-');
    hold on
end
xlabel('r_1')
ylabel('r_2')

figure(2)
for i = 1 : (size(way_2to1,1))
    plot(firing_rate(way_2to1(i, 1):way_2to1(i, 2), 1), ...
        firing_rate(way_2to1(i, 1):way_2to1(i, 2), 2), 'g-');
    hold on
end
xlabel('r_1')
ylabel('r_2')


%%
%Divide each path into almost equidistantly into (n+1) points

n = 10000;
WAY_1to2 = zeros(n + 1, 2, size(way_1to2, 1));
WAY_2to1 = zeros(n + 1, 2, size(way_2to1, 1));

%Calculate WAY_1to2
for i = 1:size(way_1to2, 1)
    tmp_num = way_1to2(i, 1):way_1to2(i, 2);
    len = 1:(way_1to2(i, 2) - way_1to2(i, 1));
    
    %%%total geometrical distance 
    for j = 1:size(len, 2)
        len(j) = sqrt((firing_rate(tmp_num(j), 1) - firing_rate(tmp_num(j + 1), 1))^2 + ...
                      (firing_rate(tmp_num(j), 2) - firing_rate(tmp_num(j + 1), 2))^2);
    end
    
    len_per = floor(n * len / sum(len));
    rem = n - sum(len_per);
    len_per(1, 1:rem) = len_per(1, 1:rem) + 1;
    
    len_num = zeros(1, size(len_per, 2) + 1);
    len_num(1, 2:end) = cumsum(len_per, 2);
    
    for j = 1:size(len_per, 2)
        for k = len_num(j):(len_num(j+1) - 1)
            lambda = (k - len_num(j)) / (len_num(j + 1) - len_num(j));
            WAY_1to2(k + 1, 1, i) = (1 - lambda) * firing_rate(tmp_num(j), 1) ...
                                    + lambda * firing_rate(tmp_num(j + 1), 1);
            WAY_1to2(k + 1, 2, i) = (1 - lambda) * firing_rate(tmp_num(j), 2) ...
                                     + lambda * firing_rate(tmp_num(j + 1), 2);
        end
    end
    WAY_1to2(n + 1, 1, i) = firing_rate(tmp_num(size(tmp_num, 2)), 1);
    WAY_1to2(n + 1, 2, i) = firing_rate(tmp_num(size(tmp_num, 2)), 2);
end

%Calculate WAY_2to1
for i = 1:size(way_2to1, 1)
    tmp_num = way_2to1(i, 1):way_2to1(i, 2);
    len = 1:(way_2to1(i, 2) - way_2to1(i, 1));
    
    %%%total geometrical distance 
    for j = 1:size(len, 2)
        len(j) = sqrt((firing_rate(tmp_num(j), 1) - firing_rate(tmp_num(j + 1), 1))^2 + ...
                      (firing_rate(tmp_num(j), 2) - firing_rate(tmp_num(j + 1), 2))^2);
    end
    
    len_per = floor(n * len / sum(len));
    rem = n - sum(len_per);
    len_per(1, 1:rem) = len_per(1, 1:rem) + 1;
    
    len_num = zeros(1, size(len_per, 2) + 1);
    len_num(1, 2:end) = cumsum(len_per, 2);
    
    for j = 1:size(len_per, 2)
        for k = len_num(j):(len_num(j+1) - 1)
            lambda = (k - len_num(j)) / (len_num(j + 1) - len_num(j));
            WAY_2to1(k + 1, 1, i) = (1 - lambda) * firing_rate(tmp_num(j), 1) ...
                + lambda * firing_rate(tmp_num(j + 1), 1);
            WAY_2to1(k + 1, 2, i) = (1 - lambda) * firing_rate(tmp_num(j), 2) ...
                + lambda * firing_rate(tmp_num(j + 1), 2);
        end
    end
    WAY_2to1(n + 1, 1, i) = firing_rate(tmp_num(size(tmp_num, 2)), 1);
    WAY_2to1(n + 1, 2, i) = firing_rate(tmp_num(size(tmp_num, 2)), 2);
    
end
%%
%Calculate the mean path

mean_1to2 = mean(WAY_1to2, 3);
mean_2to1 = mean(WAY_2to1, 3);

%%
%plot average way
figure(3)
hold on
f1 = plot(mean_1to2(:, 1), mean_1to2(:, 2), 'm.-');
hold on;
f2 = plot(mean_2to1(:, 1), mean_2to1(:, 2), 'g.-');


    

        

