theta = [102, 103, 104, 105, 106, 107];
Ec = [0.0123, 0.0153, 0.0183, 0.0214, 0.0244, 0.0276];

w = [0.05, 0.06, 0.07, 0.08, 0.09];
Ed = [0.0158, 0.0123, 0.0100, 0.0068, 0.0042];

h = [0.10, 0.11, 0.12, 0.13, 0.14, 0.15];
Ee = [0.0046, 0.0064, 0.0079, 0.0092, 0.0109, 0.0123];

phi = [60, 90, 120, 150, 180];
Ef = [0.0141, 0.0135, 0.0128, 0.0127, 0.0123];

cyan = [0.30, 0.75, 0.93];

% 1. 画布尺寸稍大，容纳更大字号
figure('Position', [100, 100, 800, 650], 'Color', 'white');

% 2. 子图位置（保持接近正方形）
ax_positions = [
    0.10,  0.60,  0.34, 0.34;   % 子图1
    0.56,  0.60,  0.34, 0.34;   % 子图2
    0.10,  0.12,  0.34, 0.34;   % 子图3
    0.56,  0.12,  0.34, 0.34;   % 子图4
];

%% 子图1：Young angle (theta)
ax1 = axes('Position', ax_positions(1,:));
plot(theta, Ec, '-o', 'Color', cyan, 'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', cyan, 'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Young angle (\theta)', 'FontSize', 14);
ylabel('Energy barrier', 'FontSize', 14);
xlim([102, 107]);
xticks(102:1:107);
ylim([0.010, 0.030]);
yticks(0.010:0.005:0.030);
set(gca, 'Box', 'on', 'FontSize', 12);
text(0.05, 0.95, '(1)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

%% 子图2：Bottom Gap (w)
ax2 = axes('Position', ax_positions(2,:));
plot(w, Ed, '-o', 'Color', cyan, 'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', cyan, 'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Bottom Gap (w)', 'FontSize', 14);
ylabel('Energy barrier', 'FontSize', 14);
xlim([0.05, 0.09]);
xticks(0.05:0.01:0.09);
ylim([0.003, 0.017]);
yticks(0.004:0.003:0.016);
set(gca, 'Box', 'on', 'FontSize', 12);
text(0.05, 0.95, '(2)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

%% 子图3：Bottom height (h)
ax3 = axes('Position', ax_positions(3,:));
plot(h, Ee, '-o', 'Color', cyan, 'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', cyan, 'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Bottom depth (d)', 'FontSize', 14);
ylabel('Energy barrier', 'FontSize', 14);
xlim([0.10, 0.15]);
xticks(0.10:0.01:0.15);
ylim([0.004, 0.014]);
yticks(0.004:0.002:0.014);
ax3.YAxis.Exponent = 0;
ytickformat('%.3f');
set(gca, 'Box', 'on', 'FontSize', 12);
text(0.05, 0.95, '(3)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');

%% 子图4：Central angle (phi)
ax4 = axes('Position', ax_positions(4,:));
plot(phi, Ef, '-o', 'Color', cyan, 'MarkerFaceColor', 'white', ...
    'MarkerEdgeColor', cyan, 'MarkerSize', 7, 'LineWidth', 1.5);
xlabel('Central angle (\phi)', 'FontSize', 14);
ylabel('Energy barrier', 'FontSize', 14);
xlim([60, 180]);
xticks(60:30:180);
ylim([0.012, 0.015]);
yticks(0.012:0.001:0.015);
set(gca, 'Box', 'on', 'FontSize', 12);
text(0.05, 0.95, '(4)', 'Units', 'normalized', 'FontSize', 18, 'FontWeight', 'bold');