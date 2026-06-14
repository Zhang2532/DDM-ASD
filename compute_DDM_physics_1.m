function [E, G] = compute_DDM_physics_1(u_1D, psi_binary, abs_grad_psi, dx, dy, epsilon, theta)
[Ny, Nx] = size(psi_binary);
u   = reshape(u_1D, Ny, Nx);
psi = psi_binary;

%% === 1. 周期前向差分梯度（Ny×Nx，与节点同维）===
% 循环边界：最后一列/行的右/上邻居绕回第一列/行
grad_u_x = (u(:, [2:end, 1]) - u) / dx;   % Ny×Nx
grad_u_y = (u([2:end, 1], :) - u) / dy;   % Ny×Nx

%% === 2. 通量（直接用 psi 乘，无维度问题）===
flux_x = psi .* grad_u_x;   % Ny×Nx
flux_y = psi .* grad_u_y;   % Ny×Nx

%% === 3. 周期后向差分散度（与前向差分对偶）===
div_psi_grad_u = (flux_x - flux_x(:, [end, 1:end-1])) / dx ...
               + (flux_y - flux_y([end, 1:end-1], :)) / dy;

%% === 4. 势能及其导数 ===
u2      = u .* u;
F_u     = 0.25 * (u2 - 1).^2;
dF_du   = u .* (u2 - 1);
d2F_du2 = 3*u2 - 1;
c       = (sqrt(2)/2) * cos(theta);
W_u     = -c * (u - u2.*u/3);
dW_du   = -c * (1 - u2);
d2W_du2 =  2*c * u;

%% === 5. 能量 ===
E_diffuse = 0.5 * epsilon * dx * dy * ...
    (psi(:)' * (grad_u_x(:).^2 + grad_u_y(:).^2));
E_bulk    = (dx * dy / epsilon) * (psi(:)' * F_u(:));
E_surf    =  dx * dy            * (abs_grad_psi(:)' * W_u(:));
E = E_diffuse + E_bulk + E_surf;

%% === 6. 梯度 G0 ===
G0 = -epsilon * div_psi_grad_u ...
     + (1/epsilon) * psi .* dF_du ...
     + abs_grad_psi .* dW_du;

%% === 7. 质量守恒投影 ===
psi_flat  = psi(:);
denom     = dx * dy * (psi_flat' * psi_flat);
numerator = dx * dy * (G0(:)' * psi_flat);
lambda    = -numerator / denom;
G_2D = G0 + lambda * psi;

G    = -G_2D(:);
end

%% =========================================================
