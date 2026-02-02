clc;clear;
close all;
%% 
d = load(  'data/field_final_resample.mat').field;
mask = load("data/mask_final.mat").mask;
% d = load('data/field_resample.mat').field;
% mask = load("data/mask.mat").mask;
figure();
subplot(211);
imagesc(d);
colormap('gray');
subplot(212);
imagesc(mask);
colormap('gray');
missing_columns = sum(mask(1, :));  % 统计缺失列数（假设整列相同）
total_columns = size(mask, 2);      % 总列数
rate = 1-(missing_columns / total_columns);
fprintf('缺失列占比: %.2f%%\n', rate*100);
k = 10;
sampled_rows = 1:k:size(d, 1);
d = d(sampled_rows,:);
mask = mask(sampled_rows,:);
lx = 40;ly = 40;

%% 可选：手动添加大量缺失
% s1 = 55;
% s2 = 65;
% d(:,s1:s2) = 0;
% mask(:,s1:s2) = 0;

%% 
d = d.*mask;
% figure();
% imagesc(d);
% colormap('gray');
d = P_H(d,lx,ly);
mask = P_H(mask,lx,ly);
% disp('正在对图像矩阵进行SVD分解...');
% [U, S, V] = svd(d, 'econ'); % 'econ' 可以在矩阵非方时提高效率
% disp('SVD分解完成。');
[L,appr_rank] = create_low_rank_direct_3(d, mask);
% [L,appr_rank] = create_low_rank_direct(d, 0.33);
disp(['找到的合适的秩序号为：',num2str(appr_rank)]);
% figure();
% plot(L);
figure();
plot(L, 'b-', 'LineWidth', 2.5);
xlim([0 200]); 
set(gca, 'FontSize', 25, 'LineWidth', 1.5);
xlabel('Rank', 'FontSize', 25);
ylabel('Frobenius norm', 'FontSize', 25);
function [dout]=P_H(din,lx,ly)
% forming trajectory matrix (also known as Hankel matrix in SSA context)
% This is a standard, robust implementation.
[nx,ny]=size(din);
lxx = nx - lx + 1; % Number of patches vertically
lyy = ny - ly + 1; % Number of patches horizontally
dout = zeros(lx*ly, lxx*lyy, 'like', din);
patch_idx = 1;
for j = 1:lyy
    for i = 1:lxx
        patch = din(i:(i+lx-1), j:(j+ly-1));
        dout(:, patch_idx) = patch(:);
        patch_idx = patch_idx + 1;
    end
end
end
function [L , appr_rank] = create_low_rank_direct_3(img, mask)
    % 1. 执行SVD
    disp('正在对图像矩阵进行SVD分解...');
    [U, S, V] = svd(img, 'econ');
    disp('SVD分解完成。');
    appr_rank = size(S, 1);
    norm_img = norm(img.*(1-mask), 'fro');
    disp(['原始图像缺失部分范数：',num2str(norm_img)]);
    L = zeros(1,(size(S, 1) -1));
    for i = 1:200
        img_appr_temp = U(:, 1:i) * S(1:i, 1:i) * V(:, 1:i)';
        L(i) = norm(img_appr_temp.*(1-mask), 'fro');
    end
end

