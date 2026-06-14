%% ── 加载四个已有结果 ──
base_path = 'E:\科研\湿润问题变式1\结果1initial\';

dt_list = [0.0002, 0.0001, 0.00005, 2.5e-5,1.25e-5,5e-6];   % 对应 dt1~dt4

phi_res = cell(6, 1);
for i = 1:6
    fname = fullfile(base_path, sprintf('phi_result1_dt%d.mat', i));
    tmp = load(fname, 'phi_p');
    phi_res{i} = tmp.phi_p(:);   % 统一成列向量
    fprintf('dt%d (dt=%.5f) 加载完成，长度 %d\n', i, dt_list(i), numel(phi_res{i}));
end

%% ── 以最细 dt（dt4）为参考解 ──
phi_ref = phi_res{6};

%% ── 计算误差和收敛阶 ──
n_compare = 5;   % dt1, dt2, dt3 对比 dt4
err_L2   = zeros(n_compare, 1);
err_Linf = zeros(n_compare, 1);

for lv = 1:n_compare
    diff         = phi_res{lv} - phi_ref;
    err_L2(lv)   = sqrt(mean(diff.^2));
    err_Linf(lv) = max(abs(diff));
end

%% ── 打印结果表 ──
fprintf('\n%-12s  %-14s  %-8s  %-14s  %-8s\n', ...
        'dt', 'L2误差', 'L2阶', 'Linf误差', 'Linf阶');

for lv = 1:n_compare
    if lv == 1
        fprintf('%-12.5f  %-14.4e  %-8s  %-14.4e  %-8s\n', ...
            dt_list(lv), err_L2(lv), '---', err_Linf(lv), '---');
    else
        p_L2   = log(err_L2(lv-1)   / err_L2(lv))   / log(dt_list(lv-1) / dt_list(lv));
        p_Linf = log(err_Linf(lv-1) / err_Linf(lv)) / log(dt_list(lv-1) / dt_list(lv));
        fprintf('%-12.5f  %-14.4e  %-8.3f  %-14.4e  %-8.3f\n', ...
            dt_list(lv), err_L2(lv), p_L2, err_Linf(lv), p_Linf);
    end
end

%% ── 收敛曲线 ──
figure;
loglog(dt_list(1:5), err_L2,   'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'L2误差');
hold on;
loglog(dt_list(1:5), err_Linf, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'L∞误差');

% 参考斜率线（过第一个点）
dt_ref_line = dt_list(1:4);
loglog(dt_ref_line, err_L2(1)*(dt_ref_line/dt_ref_line(1)).^1, ...
    'k--', 'LineWidth', 1.2, 'DisplayName', '1阶参考');
loglog(dt_ref_line, err_L2(1)*(dt_ref_line/dt_ref_line(1)).^2, ...
    'k:',  'LineWidth', 1.2, 'DisplayName', '2阶参考');

xlabel('\Deltat', 'FontSize', 13);
ylabel('误差', 'FontSize', 13);
title('时间收敛阶测试', 'FontSize', 14);
legend('Location', 'northwest');
grid on;

% 在图上标注各点的收敛阶
for lv = 2:n_compare
    p_L2 = log(err_L2(lv-1)/err_L2(lv)) / log(dt_list(lv-1)/dt_list(lv));
    text(dt_list(lv), err_L2(lv)*1.3, sprintf('p=%.2f', p_L2), ...
        'FontSize', 10, 'Color', 'b', 'HorizontalAlignment', 'center');
end
