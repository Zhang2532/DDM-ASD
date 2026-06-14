function [E, G, H_op] = compute_DDM_physics(u_1D, psi_binary, abs_grad_psi, dx, dy, epsilon, theta, c0_mass, dt)
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

%% === 8. Hessian 算子 ===
H_op = @(w_1D) apply_hessian(w_1D, psi, abs_grad_psi, ...
    d2F_du2, d2W_du2, dx, dy, epsilon, Ny, Nx, denom);
end

%% =========================================================
function Hw = apply_hessian(w_1D, psi, abs_grad_psi, ...
                             d2F_du2, d2W_du2, dx, dy, epsilon, Ny, Nx, denom)
% w = reshape(w_1D, Ny, Nx);
% 
% % 周期前向差分
% grad_w_x = (w(:, [2:end, 1]) - w) / dx;
% grad_w_y = (w([2:end, 1], :) - w) / dy;
% 
% % 通量（直接乘 psi）
% flux_x = psi .* grad_w_x;
% flux_y = psi .* grad_w_y;
% 
% % 周期后向差分散度
% div_psi_grad_w = (flux_x - flux_x(:, [end, 1:end-1])) / dx ...
%                + (flux_y - flux_y([end, 1:end-1], :)) / dy;
% 
% Hw0 = -epsilon * div_psi_grad_w ...
%       + (1/epsilon) * psi .* d2F_du2 .* w ...
%       + abs_grad_psi .* d2W_du2 .* w;
% 
% lambda_H = -(dx * dy * (Hw0(:)' * psi(:))) / denom;
% Hw = Hw0(:) + lambda_H * psi(:);
w = reshape(w_1D, Ny, Nx);
psi_flat = psi(:);

% === 第一步：投影输入 Pw（对应 P_0）===
alpha = dx * dy * (psi_flat' * w(:)) / denom;
w_proj = w - alpha * psi;   % Pw，满足 <ψ, Pw>=0

% === 第二步：H 作用在 Pw 上 ===
grad_w_x = (w_proj(:, [2:end, 1]) - w_proj) / dx;
grad_w_y = (w_proj([2:end, 1], :) - w_proj) / dy;
flux_x = psi .* grad_w_x;
flux_y = psi .* grad_w_y;
div_psi_grad_w = (flux_x - flux_x(:, [end, 1:end-1])) / dx ...
               + (flux_y - flux_y([end, 1:end-1], :)) / dy;

Hw0 = -epsilon * div_psi_grad_w ...
      + (1/epsilon) * psi .* d2F_du2 .* w_proj ...
      + abs_grad_psi .* d2W_du2 .* w_proj;

% === 第三步：投影输出 P(Hw0)（对应 P_V^T）===
lambda_H = -(dx * dy * (Hw0(:)' * psi_flat)) / denom;
Hw = Hw0(:) + lambda_H * psi_flat;

end