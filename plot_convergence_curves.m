%% 收敛梯度对比图
data_dir = 'E:\科研\湿润问题变式1\compare';
momentum_values = {'0_0.015','0.1_0.015','0.2_0.015','0.3_0.015','0.4_0.015'};
beta_labels     = {'0', '0.1', '0.2', '3', '0.4'};
% momentum_values = {'0.4_0.015'};
% beta_labels     = { '0.4'};
colors     = [0.00 0.45 0.70;
              0.85 0.33 0.10;
              0.47 0.67 0.19;
              0.49 0.18 0.56;
              0.93 0.69 0.13];
markers = {'o', 's', '^', 'd', 'v'};
linestyles = {'-', '--', ':', '-.', '-'};

figure('Color','w','Position',[100,100,900,500]);
hold on;

for k = 1:length(momentum_values)
    fname = fullfile(data_dir, sprintf('grad1_1_%s.mat', momentum_values{k}));
    if ~isfile(fname)
        fprintf('文件不存在: %s，跳过\n', fname);
        continue;
    end
    data = load(fname);
    vars = fieldnames(data);
    grad = data.(vars{1});
    grad = grad(:);
    steps = (1:length(grad)) * 200;

    % 稀疏标记索引
    marker_spacing = max(1, floor(length(steps)/15));
    idx = 1:marker_spacing:length(steps);

    % 一次画线+标记，图例自动显示两者
    semilogy(steps, grad, ...
        [linestyles{k}, markers{k}], ...
        'Color',           colors(k,:), ...
        'LineWidth',       2, ...
        'MarkerSize',      14, ...
        'MarkerFaceColor', 'w', ...
        'MarkerIndices',   idx, ...          % ← 只在稀疏点显示标记
        'DisplayName',     sprintf('\\gamma = %s', beta_labels{k}));
end

yline(1e-6, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

xlabel('Iteration steps',  'FontSize', 20);
ylabel('$\|\nabla_{\mathbf{P}}\mathcal{E}_h\|$', 'FontSize', 10, 'Interpreter', 'latex');
% title('不同动量项下的收敛曲线对比', 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 24);

set(gca, 'YScale',         'log', ...
         'YMinorTick',     'off', ...
         'FontSize',       24, ...       % ← 可以随意调大
         'GridAlpha',      0.8, ...
         'MinorGridAlpha', 0.15);

% ↓ 关键：手动固定刻度位置，防止字体变大后自动减少
yticks([1e-6, 1e-4, 1e-2, 1e0, 1e2]);

box on;
ylim([1e-7, 1e3]);
hold off;