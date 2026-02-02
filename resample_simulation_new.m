clc;
clear;
close all;

%% 读取矩阵
field = load("data/similation_file/No100_2.mat").field;
% [field, time, traces, dt] = read_B_scan_out("data/similation_file/Concrete_rebars_L_merged.out",'Ey');
% field = field-repmat(mean(field,2),1,size(field, 2));
% field = field(:,1:274);
% [field, time, traces, dt] = read_B_scan_out("data\similation_file\No5_1_bar_merged.out",'Ex');
% field = field-repmat(mean(field,2),1,size(field, 2));
% [field, time, traces, dt] = read_B_scan_out("data/similation_file/Concrete_rebars_merged.out",'Ey');
% [field_empty, time_empty, traces_empty, dt_empty] = read_B_scan_out("data/similation_file/Concrete_rebars_void_merged.out",'Ey');
% field = field - field_empty;

timewindow = 15;
dt = timewindow*1e-9/size(field,1);
traces = 1:size(field,2);
time = dt:dt:dt*(size(field,1)-1)+dt;
figure;
imagesc([], time*1e9,field);
xlabel('Trace number','fontsize',20);
ylabel('Time(ns)','fontsize',20);
set(gca,'linewidth',1,'fontsize',20);
% title('RAW','fontname','Times New Roman','linewidth',4,'fontsize',35);
colormap('gray');

% figure();
% imagesc(field);
% colormap('gray');


%% 速度曲线设定及位置标号
N = size(field,2);
v_non_uniform = ones(1,N);
x_known = [1,   19, 24, 37,  43, 88, 98, 167]; % 已知点的位置（索引）
v_known = [0.7, 1,  1,  1.3, 1,  1,   1.4, 1.4];  % 已知点的值
% x_known = [1,   30, 40, 60,  70, 120, 135, 274]; % 已知点的位置（索引）
% v_known = [0.7, 1,  1,  1.3, 1,  1,   1.7, 1.7];  % 已知点的值
% x_known = [1,   30, 45, 65,  80, 100, 130, 150, 155, 180, 200, 274]; % 已知点的位置（索引）
% v_known = [0.7, 1,  1,  1.3, 1,  1,   0.8, 1,   1,   1.2, 1,   1];  % 已知点的值
% x_known = [1, 10, 20, 30, 35, 45, 60, 70, 168]; % 已知点的位置（索引）
% v_known = [0.7, 1, 1, 0.8, 1, 1, 1.4, 1, 1];  % 已知点的值
% x_known = [1, 15, 25, 50, 80, 100, 130, 155, 168]; % 已知点的位置（索引）
% v_known = [0.7, 1, 1, 0.6, 1, 1, 1.8, 1, 1];  % 已知点的值
% x_known = [1, 7, 24 36, 54, 75, 109, 168]; % 已知点的位置（索引）
% v_known = [1.4, 1.4, 1, 1.4, 1.4, 2.6, 2.6, 2.6];  % 已知点的值
x_all = 1:x_known(end);
v = interp1(x_known, v_known, x_all, 'linear', 'extrap');
v_non_uniform = v(1:end);
% v_non_uniform = v_non_uniform.*1.2;

% 速度曲线添加随机扰动
v_smooth = v;
envelope_ratio = 0.02; % 15%的波动范围
upper_env = v_smooth * (1 + envelope_ratio);
lower_env = v_smooth * (1 - envelope_ratio);
v_enveloped = zeros(size(v_smooth));
% --- 关键修改：设置随机种子 ---
seed_value = 117;  % 设定一个固定的整数（你可以随便改，比如 1, 42, 2024）
rng(seed_value);   % 锁定随机数生成器的状态
% ---------------------------
for i = 1:length(v_smooth)
    v_enveloped(i) = lower_env(i) + (upper_env(i)-lower_env(i)) * rand();
end
v_enveloped = smoothdata(v_enveloped, 'movmean', 3);
v_non_uniform = v_enveloped;
% v_non_uniform = v_non_uniform.*1.2;

%% 可选：带插值的非均匀B扫构造
range_float = zeros(1, 10000); % 预分配大一点，稍后截断
range_float(1) = 1;
current_pos = 1;
idx = 1;

while true
    % 累加浮点数步长
    current_pos = current_pos + v_non_uniform(idx);
    
    % 边界检查：如果超出原始图像范围，停止
    if current_pos > size(field, 2)
        break;
    end
    
    idx = idx + 1;
    range_float(idx) = current_pos;
end
range_float = range_float(1:idx); 

% 原始坐标轴 (1, 2, 3, ..., N)
original_grid = 1:size(field, 2);

% 核心代码：一行搞定线性插值
% 'linear' 实现了你说的线性采样
% 'pchip' 或 'spline' 可以实现更平滑的三次样条插值（保留更多高频细节，但可能有振铃效应，推荐 linear）
field_non_uniform = interp1(original_grid, field', range_float, 'linear')';

% 3. 绘图验证
timewindow = 15;
dt = timewindow*1e-9/size(field_non_uniform,1);
time = dt:dt:dt*(size(field_non_uniform,1)-1)+dt;

figure;
imagesc([], time*1e9, field_non_uniform);
xlabel('Trace number', 'fontsize', 20);
ylabel('Time(ns)', 'fontsize', 20);
set(gca, 'linewidth', 1, 'fontsize', 20);
colormap('gray');

sss = fun_trace_detect(field_non_uniform);
figure();
plot(sss);
ylim([0.93 1]);
%% 可选：给不均匀B扫图加入一段数据缺失
start_col = 57;    
end_col = 66; 
field_non_uniform(:,start_col:end_col) = 0;





%% 速度画图


% --- 1. 定义符合真实场景的物理参数 ---
gpr_acq_freq = 20; % GPR采集频率 (Hz), 即每秒采集20道
avg_speed_mps = 0.2; % 设定的平均扫描速度 (m/s)

% --- 2. 将无量纲曲线转换为带物理单位的曲线 ---
% 生成真实的速度曲线 (m/s)
v_real_mps = v_non_uniform * avg_speed_mps;
% 生成真实的时间轴 (s)
time_axis_s = (x_all - 1) / gpr_acq_freq;

% --- 3. 绘图 ---
figure();
% figure('Name', 'Synthesized Velocity Curve with Physical Units', 'NumberTitle', 'off');
plot(time_axis_s, v_real_mps, '-b', 'LineWidth', 2);
% title('Synthesized Velocity Curve for Simulation');
set(gca, 'FontSize', 30);
xlabel('Time [s]', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Velocity[m/s]','FontSize', 20, 'FontWeight', 'bold');

grid on;
box on;
set(gca, 'FontSize', 12);



% figure();
% imagesc(field_non_uniform);
% colormap('gray');
% title("非均匀图像");


%% 速度曲线求直方图
data = v_real_mps;
% data = v_non_uniform;
nbins = 75; % 你可以根据需要调整这个数值
[counts, edges] = histcounts(data, nbins);
binWidth = edges(2) - edges(1);
pdfEstimate = counts / (sum(counts) * binWidth);
centers = (edges(1:end-1) + edges(2:end)) / 2;
[maxPdfValue, maxIndex] = max(pdfEstimate);
mostFrequentSpeed = centers(maxIndex) + binWidth/2;
maxBinEdges = [edges(maxIndex), edges(maxIndex + 1)];
disp("最大占比速度为：");
disp(mostFrequentSpeed);

figure();
subplot(311);
plot(data, '-b', 'LineWidth', 1); % 使用蓝色实线表示插值结果
title('v');
subplot(312);
bar(centers, pdfEstimate, 'histc'); % 'histc' 参数使得柱状图的宽度与bins匹配
title('pdf of v');
xlabel('Velocity (normalized)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
hold on;
plot(mostFrequentSpeed, maxPdfValue, 'ro', 'MarkerFaceColor', 'r');
text(mostFrequentSpeed, maxPdfValue, sprintf('%.2f (%.2f-%.2f)', mostFrequentSpeed, maxBinEdges(1), maxBinEdges(2)), ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');
v_non_uniform_new = data/mostFrequentSpeed;
% v_non_uniform_new = v_non_uniform_new/2;
subplot(313);
plot(v_non_uniform_new, '-b', 'LineWidth', 1); % 使用蓝色实线表示插值结果
title('normalized v');
%% 可选：带插值的数据删减、重排操作
% 1. 获取图像的实际尺寸作为基准
% field_non_uniform 是 [Time x Traces]，我们需要 Traces 数量
num_traces = size(field_non_uniform, 2); 

% 2. 检查并修正速度向量 v_non_uniform_new 的长度
% 必须保证速度向量长度足以覆盖所有道，或者如果不一致，进行截断/补齐
v_len = length(v_non_uniform_new);

if v_len < num_traces
    % 如果速度向量比图像短（少见，防止报错），用最后一个速度值补齐
    warning('速度向量长度小于图像宽度，正在补齐...');
    v_pad = repmat(v_non_uniform_new(end), 1, num_traces - v_len);
    v_correction = [v_non_uniform_new, v_pad];
elseif v_len > num_traces
    % 如果速度向量比图像长，直接截断
    v_correction = v_non_uniform_new(1:num_traces);
else
    % 长度一致
    v_correction = v_non_uniform_new;
end

% 3. 计算每一道非均匀数据的“真实浮点坐标”
% 强制初始化为 num_traces 长度，确保与 field_non_uniform' 的行数严格一致
range_float_pos = zeros(1, num_traces);
range_float_pos(1) = 1;
current_pos = 1;

% 循环计算位置 (注意循环次数是 num_traces - 1)
for i = 1:num_traces-1
    current_pos = current_pos + v_correction(i);
    range_float_pos(i+1) = current_pos;
end

% 4. 确定“要保留的整数道” (Target Indices)
target_indices_raw = round(range_float_pos);

% 边界处理：防止 round 后索引小于1
target_indices_raw(target_indices_raw < 1) = 1;

% 去重：得到最终需要计算数值的整数位置 (如 [1, 2, 5])
target_unique = unique(target_indices_raw);

% 5. 初始化稀疏矩阵
% 最大的道数由我们算出的 target_unique 决定
max_target_idx = max(target_unique);
field_resample = zeros(size(field_non_uniform, 1), max_target_idx); 

% 6. 核心步骤：插值重构数值
% 这里的 field_non_uniform' 行数是 num_traces
% 这里的 range_float_pos 长度也是 num_traces
% 维度完美匹配，不会再报错
interpolated_traces = interp1(range_float_pos, field_non_uniform', target_unique, 'linear')';

% 处理可能的 NaN (如果插值点在范围外)
interpolated_traces(isnan(interpolated_traces)) = 0;

% 7. 填入数据
field_resample(:, target_unique) = interpolated_traces;

% 8. 更新 missing_values (供后续 mask 和绘图使用)
missing_values = find(all(field_resample == 0, 1));
% full_sequence = 1:size(field_resample, 2);
% missing_values = setdiff(full_sequence, target_unique);





%%  旧的重构方法
% v_non_uniform_new = v_non_uniform_new(1:size(field_non_uniform,2));
% % 求更新的速度对应的标号
% range = 1;
% range_non_uniform_new(1) = 1;
% for i= 1:size(v_non_uniform_new,2)-1
%     range = range + v_non_uniform_new(i)*1;
%     range_non_uniform_new(i+1) = range;
% end
% range_non_uniform_new = round(range_non_uniform_new);
% 
% 
% % 重新排列
% field_resample = zeros(size(field_non_uniform,1),range_non_uniform_new(end));
% for i = 1:size(field_non_uniform,2)
%     field_resample(:,range_non_uniform_new(i)) = field_non_uniform(:,i);
% end
% missing_values = find(all(field_resample == 0, 1));
% % unique_data = unique(range_non_uniform_new);  % 去重后的数据
% % full_sequence = 1:range_non_uniform_new(end);  % 完整的整数序列
% % missing_values = setdiff(full_sequence, unique_data);


%%
% 把数据丢失也加进去
% missing_values = union(missing_values, [60, 170]);
% field_resample(:,[60,170]) = 0;

figure();
subplot(211);
imagesc(field_non_uniform);
colormap('gray');
title("非均匀图像");
subplot(212);
imagesc(field_resample);
colormap('gray');
title("重构图像");
hold on;
for col = missing_values
    % 绘制从顶部到底部的垂直线
    plot([col, col], [1,100], 'Color', 'blue', 'LineWidth', 2);
end


timewindow = 15;
dt = timewindow*1e-9/size(field_resample,1);
traces = 1:size(field_resample,2);
time = dt:dt:dt*(size(field_resample,1)-1)+dt;
figure;
imagesc([], time*1e9,field_resample);
xlabel('Trace number','fontsize',20);
ylabel('Time(ns)','fontsize',20);
set(gca,'linewidth',1,'fontsize',20);
colormap('gray');
hold on;
y_limits = ylim;
x_limits = xlim;
% for col = missing_values
%     plot([col, col], y_limits, 'Color', 'blue', 'LineWidth', 2);
% end
short_y_length = 0.05 * (y_limits(2) - y_limits(1));
short_y_limits = [y_limits(1), y_limits(1) + short_y_length];
for col = missing_values
    plot([col, col], short_y_limits, 'Color', 'blue', 'LineWidth', 2);
end
%% 存储数据
mask = ones(size(field_resample));
mask(:,missing_values) = 0;
figure();
subplot(211);
imagesc(field_resample);
colormap('gray');
hold on;
y_limits = ylim;
x_limits = xlim;
short_y_length = 0.05 * (y_limits(2) - y_limits(1));
short_y_limits = [y_limits(1), y_limits(1) + short_y_length];
for col = missing_values
    plot([col, col], short_y_limits, 'Color', 'blue', 'LineWidth', 2);
end
subplot(212);
imagesc(mask);
colormap('gray');

% field = field_non_uniform;
% save("data/field_non_uniform_new.mat","field");
% field = field_resample;
% save("data/field_final_resample.mat","field");
% save("data/mask_final_new.mat","mask");





