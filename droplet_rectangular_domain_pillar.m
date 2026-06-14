function [u, X, Y, psi, abs_grad_psi, dx, dy, params] = droplet_rectangular_domain_pillar()
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

%% ====================== 3. 构造固体表面高度函数 y_surface(X) ======================
y_surface = zeros(Ny, Nx);
tol = 1e-6;

for i = 1:n_pillars
    pillar_col = (X > pillar_lefts(i) + tol) & (X < pillar_rights(i) - tol);
    y_surface(pillar_col) = h-tol;
end

%% ====================== 4. 构造 ψ ======================
x_phys_left  = 0;
x_phys_right = 1;
y_phys_up    = H_total;

psi = zeros(Ny, Nx);
psi( Y >= y_surface & Y <= y_phys_up & ...
     X >= x_phys_left - tol & X <= x_phys_right + tol ) = 1;

%% ====================== 5. 构造 delta 函数 ======================
abs_grad_psi = zeros(Ny, Nx);

psi_below    = [zeros(1,Nx);    psi(1:end-1,:)];
psi_above    = [psi(2:end,:);   zeros(1,Nx)   ];
psi_left_nb  = [zeros(Ny,1),    psi(:,1:end-1)];
psi_right_nb = [psi(:,2:end),   zeros(Ny,1)   ];

% 两侧都加容差，保证 x=0 和 x=1 的格点都被包含
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
cx = 0.5; r = 0.25; cy = 0.1 + r -0.08;
inside_circle = (X - cx).^2 + (Y - cy).^2 <= r^2;
u(inside_circle) = 1;
u(Y < y_surface) = -1;

%% ── 精确质量补偿（贪心逐点） ──────────────────────────────────────────
target_mass  = -0.350425;
% target_mass  = -0.5;
current_mass = sum(u(:) .* psi(:)) * dx * dy;
M_deficit    = target_mass - current_mass;   % 正 = 需要增加质量

if abs(M_deficit) < 1e-14
    fprintf('初始质量完美匹配，无需补偿。\n');
else
    cell_mass = 2 * dx * dy;   % 单个格点从 ±1 翻转贡献的质量变化

    if M_deficit > 0
        % 质量不足：在液滴外边缘（u=-1）由近到远翻转为 +1
        candidates = find(u == -1 & psi == 1);
        dist2 = (X(candidates) - cx).^2 + (Y(candidates) - cy).^2;
        [~, ord] = sort(dist2, 'ascend');   % 从边界向外
        candidates = candidates(ord);
        flip_sign  = +1;   % -1 → +1
    else
        % 质量过多：在液滴内边缘（u=+1）由远到近翻转为 -1
        candidates = find(u == 1 & psi == 1);
        dist2 = (X(candidates) - cx).^2 + (Y(candidates) - cy).^2;
        [~, ord] = sort(dist2, 'descend');  % 从边界向内
        candidates = candidates(ord);
        flip_sign  = -1;   % +1 → -1
    end

    remaining = abs(M_deficit);
    n_full    = 0;

    for k = 1:length(candidates)
        idx = candidates(k);

        if remaining >= cell_mass - 1e-15
            % 整格翻转
            u(idx)  = flip_sign;
            remaining = remaining - cell_mass;
            n_full    = n_full + 1;
        else
            % 最后一个格点：精确补齐，允许中间值
            u(idx) = u(idx) + sign(M_deficit) * remaining / (dx * dy);
            u(idx) = max(-1, min(1, u(idx)));   % 防浮点溢出
            remaining = 0;
            break;
        end

        if remaining < 1e-15, break; end
    end

    final_mass = sum(u(:) .* psi(:)) * dx * dy;
    fprintf('整格翻转: %d 个  |  末格精确补偿: 1 个\n', n_full);
    fprintf('目标质量: %.10f\n实际质量: %.10f\n残差:     %.2e\n', ...
            target_mass, final_mass, abs(final_mass - target_mass));
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

%% ====================== 8. 可视化 ======================
figure('Color','w','Position',[100,100,700,400]);
pcolor(X, Y, u);
shading interp;
colormap(parula);
hold on;

for i = 1:n_pillars
    fill([pillar_lefts(i), pillar_rights(i), pillar_rights(i), pillar_lefts(i)], ...
         [0, 0, h, h], [0.7 0.7 0.7], 'EdgeColor', 'none');
end

axis equal;
axis([x_phys_left, x_phys_right, 0, H_total]);
axis off;
set(gca, 'Position', [0, 0, 1, 1]);
colorbar;
title(sprintf('矩形计算域 | %d根柱子', n_pillars));
hold off;
end