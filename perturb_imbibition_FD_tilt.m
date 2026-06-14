function phi_out_1D = perturb_imbibition_FD_tilt(phi_in_1D, X, Y, psi, dx, dy, gap_index, params)
    [Ny, Nx] = size(X);
    phi_out = reshape(phi_in_1D, Ny, Nx);

    %% 1. 读取几何参数
    h             = params.h;
    w             = params.w;
    w_p           = params.w_p;
    n_pillars     = params.n_pillars;
    pillar_lefts  = params.pillar_lefts;
    pillar_rights = params.pillar_rights;
    dx_pillar     = params.dx_pillar;

    if gap_index < 1 || gap_index >= n_pillars
        error('gap_index 超出范围，有效范围为 1 到 %d', n_pillars - 1);
    end

    %% 2. 间隙边界（随高度线性变化，与斜面一致）
    % pillar_lefts/rights 存的是顶部坐标，底部需减去 dx_pillar
    frac = min(max(Y / h, 0), 1);   % Ny×Nx，clamp [0,1]
    x_gap_left_2D  = pillar_rights(gap_index)    - dx_pillar .* (1 - frac);
    x_gap_right_2D = pillar_lefts(gap_index + 1) - dx_pillar .* (1 - frac);

    %% 3. 扰动区域（间隙底部，高度方向 Y <= h/6）
    in_gap      = (X >= x_gap_left_2D) & (X <= x_gap_right_2D);
    in_lower    = (Y <= h/0.5 );
    target_mask = in_gap & in_lower & (psi > 0.5);
%     target_mask = in_gap & in_lower ;
    if ~any(target_mask(:))
        warning('未找到指定的扰动区域！请检查网格或间隙位置。');
        phi_out_1D = phi_out(:);
        return;
    end

    %% 4. 扰动前质量
    C0_before = sum(phi_out(:) .* psi(:)) * dx * dy;

    %% 5. 施加扰动
    phi_temp = phi_out;
    phi_temp(target_mask) = -1;

    %% 6. 质量损失
    C0_temp   = sum(phi_temp(:) .* psi(:)) * dx * dy;
    mass_lost = C0_before - C0_temp;

    %% 7. 补偿区域
%     x_min = 0.5 - 0.11;   x_max = 0.5 + 0.19;
%      x_min = 0.5 - 0.2;   x_max = 0.5 + 0.1;
%  x_min = 0.5 - 0.22;   x_max = 0.5 + 0.08;
% x_min = 0.5 - 0.15;   x_max = 0.5 + 0.15;
%     y_min = 3.8 * 0.15;   y_max = 3.8 * 0.15 + 0.05;
%     y_min = 3.5 * 0.15;   y_max = 3.5 * 0.15 + 0.05;
%  y_min = 4 * 0.15;   y_max = 4 * 0.15 + 0.05;
%  x_min = 0.27  ;   x_max = 0.27 + 0.05;
%  y_min = 0.33 - 0.15;   y_max = 0.33 + 0.15;
% x_min = 0.2-0.06;x_max = x_min + 0.05;
% x_min = 0.68+0.07;x_max = x_min + 0.05;
x_min = 0.75-0.05;x_max = x_min + 0.05;
y_min = 0.2 ;y_max = y_min + 0.2;
    comp_mask = (X >= x_min) & (X <= x_max) & ...
                (Y >= y_min) & (Y <= y_max) & (psi > 0.5);

    if ~any(comp_mask(:))
        warning('没有找到用于补偿质量的区域！质量将不再守恒。');
    else
        vol_comp = sum(psi(comp_mask)) * dx * dy;
        phi_temp(comp_mask) = phi_temp(comp_mask) + (mass_lost / vol_comp);
        phi_temp(comp_mask) = min(1, phi_temp(comp_mask));
    end

    %% 8. 检查最终质量
    phi_out  = phi_temp;
    C0_after = sum(phi_out(:) .* psi(:)) * dx * dy;
    fprintf('扰动节点数 = %d\n', sum(target_mask(:)));
    fprintf('补偿节点数 = %d\n', sum(comp_mask(:)));
    fprintf('质量变化: %.6f → %.6f (差值: %.2e)\n', C0_before, C0_after, C0_after - C0_before);

    %% 9. 可视化（平行四边形柱子）
    figure('Name', sprintf('Gap %d 扰动结果', gap_index), 'Color', 'w');
    pcolor(X, Y, phi_out);
    shading interp;
    colormap(parula);
    caxis([-1, 1]);
    colorbar;
    hold on;

    for i = 1:n_pillars
        px = [pillar_lefts(i)-dx_pillar, pillar_rights(i)-dx_pillar, ...
              pillar_rights(i),           pillar_lefts(i)];
        py = [0, 0, h, h];
        fill(px, py, [0.65 0.65 0.65], 'EdgeColor', 'k', 'LineWidth', 0.8);
    end

    axis equal;
    axis([0, 1.1, 0, 0.8]);
    axis off;
    title(sprintf('Gap %d Imbibition 扰动及上方质量补偿', gap_index));
    hold off;

    %% 10. 转回 1D
    phi_out_1D = phi_out(:);
end