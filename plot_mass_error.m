clear; clc; close all;

% 定义三个文件的绝对路径
files = {
    'E:\科研\湿润问题变式1\结果1initial\mass_saddle_1.mat',
    'E:\科研\湿润问题变式1\结果1initial\mass_saddle_2.mat',
    'E:\科研\湿润问题变式1\结果1initial\mass_saddle_3.mat'
};

% 循环处理每个文件
for i = 1:length(files)
    % 1. 加载数据
    if ~isfile(files{i})
        warning(['找不到文件: ', files{i}]);
        continue;
    end
    data = load(files{i});
    
    % 自动提取 mat 文件中的变量（假设每个文件里只有一个误差数组）
    fields = fieldnames(data);
    err_array = data.(fields{1});
    
    % 确保数据是列向量
    err_array = err_array(:);
    
    % 2. 构造 X 轴（每隔 200 个点计数）
    % 第1个数据对应步数200，第2个对应400，以此类推...
    steps = (1:length(err_array)) * 200;
    max_step = steps(end); % 获取当前数据的最大步数
    
    % 3. 创建图形并绘制
    fig = figure('Color', 'w', 'Position', [100+i*50, 100+i*50, 600, 450]);
    plot(steps, err_array, 'b-', 'LineWidth', 1.5);
    
    % 4. 设置坐标轴格式
    set(gca, 'FontSize', 20, 'LineWidth', 1);
    xlabel('Iteration steps', 'FontSize', 20);
    ylabel('Error of the Mass', 'FontSize', 20);
    
    % 限制 X 轴的范围，严格从 0 到当前最大步数（不留多余空白）
    xlim([0, max_step]);
    
    % （可选）强制 X 轴使用科学计数法，就像您截图里右下角的 x 10^4
    ax = gca;
    ax.XAxis.Exponent = 4; 
    
%     % 5. 自动保存为高质量的 PNG 图片，准备插入 LaTeX
%     output_filename = sprintf('mass_error_%d.png', i);
%     exportgraphics(fig, output_filename, 'Resolution', 300);
%     
%     fprintf('成功生成并保存图片: %s (最大步数: %d)\n', output_filename, max_step);
end