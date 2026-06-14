function [u, X, Y, psi, abs_grad_psi, dx, dy, params] = droplet_rectangular_domain_pillar_pencil()
%% ====================== 1. 定义矩形计算域 ======================
x_left = -0.1;    x_right = 1.1;
y_down = -0.1;    y_up = 0.9;
Nx = 241; Ny = 201;
% Nx = 121; Ny = 101;
dx = (x_right - x_left) / (Nx - 1);
dy = (y_up - y_down) / (Ny - 1);
x = linspace(x_left, x_right, Nx);
y = linspace(y_down, y_up, Ny);
[X, Y] = meshgrid(x, y);

%% ====================== 2. 定义柱子几何参数 ======================
W_total   = 1.0;
H_total   = 0.8;
h         = 0.15;
w         = 0.06;
w_p       = 0.06;
alpha_deg = 90;
dx_pillar = h * cotd(alpha_deg);

pillar_lefts  = [];
pillar_rights = [];
current_x = 0.05;

while (current_x + w_p) <= W_total + 1e-6
    pillar_lefts(end+1)  = current_x + dx_pillar;
    pillar_rights(end+1) = current_x + dx_pillar + w_p;
    current_x = current_x + w_p + w;
end
n_pillars = length(pillar_lefts);

fprintf('柱数: %d,  跨度/间隙严格保持: %.4f\n', n_pillars, w);

%% ====================== 3 & 4. 极严密 2D 几何构造“铅笔形/尖顶” ψ 相场 ======================
% 铅笔形几何特征：下半部为直立矩形，上半部为三角形尖顶
h_tri = 0.06;               % 三角形尖顶的高度 (设定为柱宽大小，看起来比例最协调)
h_straight = h - h_tri;     % 下方直立矩形底座的高度
tol = 1e-6;

% 初始化全域为流体区域 (psi = 1)
psi = ones(Ny, Nx);

% 切除计算域外部及地板 (psi = 0)
x_phys_left  = 0;
x_phys_right = 1;
y_phys_up    = H_total;
psi( Y < 0 | Y > y_phys_up | X < x_phys_left - tol | X > x_phys_right + tol ) = 0;

for i = 1:n_pillars
    x_L = pillar_lefts(i);
    x_R = pillar_rights(i);
    xc = (x_L + x_R) / 2;
    
    % 严格确保柱间距精确的横向边界约束
    valid_X = (X > x_L + tol) & (X < x_R - tol);
    
    % 1. 判定矩形底座区域
    in_rect = valid_X & (Y >= 0) & (Y <= h_straight);
    
    % 2. 判定尖顶三角形区域：当 Y 从 h_straight 升至 h 时，允许的宽度从 w_p 线性收缩到 0
    in_tri = valid_X & (Y > h_straight) & (Y <= h) & ...
             (abs(X - xc) < (w_p / 2) * (1 - (Y - h_straight) / h_tri) + tol);
    
    % 满足任意一部分则判定为固体柱体 (psi = 0)
    psi(in_rect | in_tri) = 0;
end

%% ====================== 5. 构造 delta 函数 ======================
abs_grad_psi = zeros(Ny, Nx);

psi_below    = [zeros(1,Nx);    psi(1:end-1,:)];
psi_above    = [psi(2:end,:);   zeros(1,Nx)   ];
psi_left_nb  = [zeros(Ny,1),    psi(:,1:end-1)];
psi_right_nb = [psi(:,2:end),   zeros(Ny,1)   ];

in_phys = (X >= x_phys_left - tol) & (X <= x_phys_right + tol);

bottom_interface = (psi > 0.5) & (psi_below    < 0.5) & in_phys;
top_interface    = (psi > 0.5) & (psi_above    < 0.5) & in_phys;
left_interface   = (psi > 0.5) & (psi_left_nb  < 0.5) & in_phys & ~bottom_interface;
right_interface  = (psi > 0.5) & (psi_right_nb < 0.5) & in_phys & ~bottom_interface;

abs_grad_psi(bottom_interface) = abs_grad_psi(bottom_interface) + 1/dy;
abs_grad_psi(top_interface)    = abs_grad_psi(top_interface)    + 1/dy;
abs_grad_psi(left_interface)   = abs_grad_psi(left_interface)   + 1/dx;
abs_grad_psi(right_interface)  = abs_grad_psi(right_interface)  + 1/dx;

fprintf('delta函数积分 = %.4f\n', sum(abs_grad_psi(:))*dx*dy);

%% ====================== 6. 初始化液滴 ======================
u = -ones(Ny, Nx);
cx = 0.5; r = 0.25; cy = 0.1 + r -0.02;
inside_circle = (X - cx).^2 + (Y - cy).^2 <= r^2;
u(inside_circle) = 1;

% 确保固体内部（尖顶柱子内部）严格标记为 -1
u(psi == 0) = -1;

%% ── 精确质量补偿（贪心逐点） ──────────────────────────────────────────
target_mass  = -0.350425;
current_mass = sum(u(:) .* psi(:)) * dx * dy;
M_deficit    = target_mass - current_mass;   

if abs(M_deficit) < 1e-14
    fprintf('初始质量完美匹配，无需补偿。\n');
else
    cell_mass = 2 * dx * dy;   

    if M_deficit > 0
        candidates = find(u == -1 & psi == 1);
        dist2 = (X(candidates) - cx).^2 + (Y(candidates) - cy).^2;
        [~, ord] = sort(dist2, 'ascend');   
        candidates = candidates(ord);
        flip_sign  = +1;   
    else
        candidates = find(u == 1 & psi == 1);
        dist2 = (X(candidates) - cx).^2 + (Y(candidates) - cy).^2;
        [~, ord] = sort(dist2, 'descend');  
        candidates = candidates(ord);
        flip_sign  = -1;   
    end

    remaining = abs(M_deficit);
    n_full    = 0;

    for k = 1:length(candidates)
        idx = candidates(k);

        if remaining >= cell_mass - 1e-15
            u(idx)  = flip_sign;
            remaining = remaining - cell_mass;
            n_full    = n_full + 1;
        else
            u(idx) = u(idx) + sign(M_deficit) * remaining / (dx * dy);
            u(idx) = max(-1, min(1, u(idx)));   
            remaining = 0;
            break;
        end

        if remaining < 1e-15, break; end
    end

    final_mass = sum(u(:) .* psi(:)) * dx * dy;
    fprintf('整格翻转: %d 个  |  末格精确补偿: 1 个\n', n_full);
end
%% ====================== 7. 参数打包输出 ======================
params = struct();
params.h             = h;
params.w             = w;
params.w_p           = w_p;
params.n_pillars     = n_pillars;
params.pillar_lefts  = pillar_lefts;
params.pillar_rights = pillar_rights;
params.dx_pillar     = dx_pillar;
params.x_phys_left   = x_phys_left;
params.x_phys_right  = x_phys_right;
params.H_total       = H_total;

%% ====================== 8. 可视化 (无缝抗锯齿铅笔形绘制) ======================
figure('Color','w','Position',[100,100,700,400]);
pcolor(X, Y, u);
shading interp;
colormap(parula);
hold on;

for i = 1:n_pillars
    x_L = pillar_lefts(i);
    x_R = pillar_rights(i);
    xc = (x_L + x_R) / 2;
    
    % 构建铅笔形多边形的 5 个顶点 (首尾闭合，兼容倾斜 dx_pillar)
    % 顺序：左底角 -> 右底角 -> 右直道顶 -> 尖顶 -> 左直道顶 -> 闭合左底角
    x_poly = [x_L - dx_pillar, x_R - dx_pillar, x_R, xc, x_L, x_L - dx_pillar];
    y_poly = [0,               0,               h_straight, h, h_straight, 0  ];
    
    % 填充并应用抗锯齿平滑渲染，保证线条圆润没有像素狗牙
    fill(x_poly, y_poly, [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5);
end

axis equal;
axis([x_phys_left, x_phys_right, 0, H_total]);
axis off;
set(gca, 'Position', [0, 0, 1, 1]);
colorbar;
title(sprintf('矩形计算域 | %d根尖顶铅笔柱', n_pillars));
hold off;
end