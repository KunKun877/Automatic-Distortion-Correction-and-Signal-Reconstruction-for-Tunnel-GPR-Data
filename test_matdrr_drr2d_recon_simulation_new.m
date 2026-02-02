clc;clear;
close all;
%% 加载数据
d = load('data/field_final_resample.mat').field;
% rate = 0.45;
% mask=drr_genmask(d,1-rate,'c',201415);
mask = load("data/mask_final.mat").mask;
missing_columns = sum(mask(1, :));  % 统计缺失列数（假设整列相同）
total_columns = size(mask, 2);      % 总列数
rate = 1-(missing_columns / total_columns);
fprintf('缺失列占比: %.2f%%\n', rate*100);
d1 = d.*mask;

%% 可选：手动添加大量缺失
% s1 = 55;
% s2 = 65;
% d(:,s1:s2) = 0;
% mask(:,s1:s2) = 0;


%% 进行降采样，增加程序运行速度
k = 10;
sampled_rows = 1:k:size(d, 1);
d = d(sampled_rows,:);
mask = mask(sampled_rows,:);

%% 画图
figure();
subplot(211);
imagesc(d);
colormap('gray');
colorbar;
title("丢失数据图像");
subplot(212);
imagesc(mask);
colormap('gray');
colorbar;
title("丢失数据位置");


%% 
N=3;Niter = 80; mode=0; verb=1; win_t = 40; win_x = 40;
dr = drr_spatial_recon(d, mask, 20, 10, Niter, eps, verb, mode, 0, win_t, win_x);
dd = d-dr;
timewindow = 15;
dt = timewindow*1e-9/size(dr,1);
traces = 1:size(dr,2);
time = dt:dt:dt*(size(dr,1)-1)+dt;
figure;
imagesc([], time*1e9,dr);
xlabel('Trace number','fontsize',20);
ylabel('Time(ns)','fontsize',20);
set(gca,'linewidth',1,'fontsize',20);
colormap('gray');

figure();
subplot(311);
imagesc(d);
colormap('gray');
colorbar;
subplot(312);
imagesc(dr);
colormap('gray');
colorbar;
subplot(313);
imagesc(dd);
colormap('gray');
colorbar;     

%%
figure();
subplot(211);
imagesc(d);
colormap('gray');
subplot(212);
imagesc(dr);
colormap('gray');

% field1 = load('../data/result/f_1_500.mat').field;
% field2 = load('../data/result/f_501_1000.mat').field;
% field3 = load('../data/result/f_1001_1400.mat').field;
% 
% field = [field1,field2,field3];
% 
% figure();
% imagesc(field);
% colormap('gray');
% colorbar();


%% 平滑处理  一维测试
% field = load('../data/result/f_1_500.mat').field;
% figure();
% imagesc(field);
% colormap('gray');
% 
% flow = field(85,:);
% flow2 = field(260,:);
% 
% filter_width = 5;
% s_flow = movmean(flow, filter_width);
% s_flow2 = movmean(flow2, filter_width);
% 
% figure();
% subplot(211);
% plot(flow, 'LineWidth', 2)
% subplot(212);
% plot(s_flow, 'LineWidth', 2)
% 
% figure();
% subplot(211);
% plot(flow2, 'LineWidth', 2)
% subplot(212);
% plot(s_flow2, 'LineWidth', 2)

%% 平滑处理 二维测试
field = dr;
% field = load('../data/result/f.mat').field;

for i = 1:10
    s_field = field;
    filter_width = 3;
    for i = 1:size(field, 1) % 遍历图像的每一行
        s_field(i, :) = movmean(field(i, :), filter_width);
    end
    field_t = field.*(mask) + s_field.*(1-mask);
end

figure();
subplot(311);
imagesc(field);
colormap('gray');
subplot(312);
imagesc(field_t);
colormap('gray');
subplot(313);
imagesc(d);
colormap('gray');
field = field_t;

figure();
subplot(211);
imagesc(d);
colormap('gray');
hold on;
missing_values = find(all(d == 0, 1));
y_limits = ylim;
x_limits = xlim;
short_y_length = 0.05 * (y_limits(2) - y_limits(1));
short_y_limits = [y_limits(1), y_limits(1) + short_y_length];
for col = missing_values
    plot([col, col], short_y_limits, 'Color', 'blue', 'LineWidth', 2);
end
subplot(212);
imagesc(field);
colormap('gray');



% ss = field_t;
% [U,S,V]=svd(ss,'econ');
% s = diag(S);
% figure();
% plot(s);








