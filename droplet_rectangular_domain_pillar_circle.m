function [u, X, Y, psi, abs_grad_psi, dx, dy, params] = droplet_rectangular_domain_pillar_circle()
%% ====================== 1. 定义矩形计算域 ======================
x_left = -0.1;    x_right = 1.1;
y_down = -0.1;    y_up = 0.9;
% Nx = 241; Ny = 201;
% Nx = 121; Ny = 101;
Nx = 401; Ny = 334;
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

%% ====================== 3 & 4. 支持【弧度调节】的极严密 2D 相场 ======================
tol = 1e-6;

% ▼▼▼ 核心可调参数：圆心角 ▼▼▼
arc_angle_deg = 180;  % 【在这里修改弧度！】180 为标准半圆；往小调(如 120, 90)弧顶就会变平
% ▲▲▲ 核心可调参数：圆心角 ▲▲▲

% 安全限制，防止出现数学除 0 错误
arc_angle_deg = max(1, min(180, arc_angle_deg)); 
alpha_rad = arc_angle_deg * pi / 180;

% 通过三角函数，动态计算满足条件的 圆半径、拱高 和 真实的圆心 Y 坐标
R = (w_p / 2) / sin(alpha_rad / 2);        % 动态外接圆半径 (弧越平，半径越大)
h_dome = R * (1 - cos(alpha_rad / 2));     % 弧顶纯凸起的高度
h_straight = h - h_dome;                   % 柱子直立（矩形）部分的高度
y_c = h - R;                               % 真实的圆心 Y 坐标 (随弧度变化)

% 初始化全域为流体区域 (psi = 1)
psi = ones(Ny, Nx);

% 切除计算域外部及地板，使其成为固体 (psi = 0)
x_phys_left  = 0;
x_phys_right = 1;
y_phys_up    = H_total;
psi( Y < 0 | Y > y_phys_up | X < x_phys_left - tol | X > x_phys_right + tol ) = 0;

for i = 1:n_pillars
    x_L = pillar_lefts(i);
    x_R = pillar_rights(i);
    xc = (x_L + x_R) / 2;   
    
    valid_X = (X > x_L + tol) & (X < x_R - tol);
    
    % 1. 直立矩形柱体判定
    in_rect = valid_X & (Y >= 0) & (Y <= h_straight);
    
    % 2. 严密的二维圆缺方程 (圆心变成了 y_c)
    in_semi = valid_X & ((X - xc).^2 + (Y - y_c).^2 <= R^2 + tol) & (Y > h_straight);
    
    psi(in_rect | in_semi) = 0;
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

% 你现在可以随意修改这里的中心了，不会再变形了
cx = 0.5; 
r = 0.25; 
cy = 0.1 + r -0.05;

inside_circle = (X - cx).^2 + (Y - cy).^2 <= r^2;
u(inside_circle) = 1;

% 覆盖固体部分
u(psi == 0) = -1;

%% ── 精确质量补偿（已解开锁定） ──────────────────────────────────────────
current_mass = sum(u(:) .* psi(:)) * dx * dy;

% 【重要】直接让目标质量等于你画的圆，停止强行扭曲液滴形状！
target_mass  =  -0.350425; 
M_deficit    = target_mass - current_mass;   

if abs(M_deficit) < 1e-14
    fprintf('初始质量完美匹配，无需补偿。\n');
else
    % 如果将来有强制设定的需求，这里的补偿算法依然可用
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
    for k = 1:length(candidates)
        idx = candidates(k);
        if remaining >= cell_mass - 1e-15
            u(idx)  = flip_sign;
            remaining = remaining - cell_mass;
        else
            u(idx) = u(idx) + sign(M_deficit) * remaining / (dx * dy);
            u(idx) = max(-1, min(1, u(idx)));   
            break;
        end
    end
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
params.arc_angle_deg =arc_angle_deg;
%% ====================== 8. 可视化 (匹配动态圆弧的抗锯齿画图) ======================
figure('Color','w','Position',[100,100,700,400]);
pcolor(X, Y, u);
shading interp;
colormap(parula);
hold on;

for i = 1:n_pillars
    x_L = pillar_lefts(i);
    x_R = pillar_rights(i);
    xc = (x_L + x_R) / 2;
    
    % 动态计算画图所需的弧度范围：保证画图始终从柱子的右顶点画到左顶点
    theta = linspace(pi/2 - alpha_rad/2, pi/2 + alpha_rad/2, 200); 
    
    x_semi = xc + R * cos(theta);
    y_semi = y_c + R * sin(theta); % 动态圆心 y_c
    
    % 无缝拼接成一个连续的一体化多边形
    x_poly = [x_L - dx_pillar, x_R - dx_pillar, x_semi, x_L - dx_pillar];
    y_poly = [0,               0,               y_semi, 0              ];
    
    % 填充多边形 (带有抗锯齿保护的边框)
    fill(x_poly, y_poly, [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5);
end

axis equal;
axis([x_phys_left, x_phys_right, 0, H_total]);
axis off;
set(gca, 'Position', [0, 0, 1, 1]);
colorbar;
title(sprintf('矩形计算域 | %d根柱子 (可调圆弧 %d°)', n_pillars, arc_angle_deg));
hold off;
end