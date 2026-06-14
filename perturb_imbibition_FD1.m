function phi_out_1D = perturb_imbibition_FD1(phi_in_1D, X, Y, psi, dx, dy, gap_index, params)

    [Ny, Nx] = size(X);
    phi_out = reshape(phi_in_1D, Ny, Nx);

    %% 1. 从 params 读取几何参数
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

    %% 2. 查间隙位置
    x_gap_left  = pillar_rights(gap_index)  - dx_pillar;
    x_gap_right = pillar_lefts(gap_index+1) - dx_pillar;

    fprintf('=== 执行 Imbibition 扰动 ===\n');
    fprintf('目标 Gap %d : X = [%.4f, %.4f]\n', gap_index, x_gap_left, x_gap_right);

    %% 3. 确定扰动区域（间隙上半部分）
    in_gap           = (X >= x_gap_left  - 1e-6) & (X <= x_gap_right + 1e-6);
    y_perturb_bottom = h / 7;
    in_height        = (Y >= y_perturb_bottom) & (Y <= h);
    target_mask      = in_gap & in_height ;

    if ~any(target_mask(:))
        warning('未找到指定的扰动区域！');
        phi_out_1D = phi_out(:);
        return;
    end

    %% 4. 扰动前质量
    C0_before = sum(phi_out(:) .* psi(:)) * dx * dy;

    %% 5. 施加扰动：间隙上半段设为液体
    phi_temp = phi_out;
    phi_temp(target_mask) = 1;

    %% 6. 多出的质量
    C0_temp     = sum(phi_temp(:) .* psi(:)) * dx * dy;
    mass_gained = C0_temp - C0_before;

    %% 7. 内部抽水补偿
%     x_min = 0.5 - 0.15;   x_max = 0.5 + 0.15;
%     x_min = 0.5 - 0.1;   x_max = 0.5 + 0.2;%w=0.08
    x_min = 0.5 - 0.1;   x_max = 0.5 + 0.07;%w=0.08
    y_min = 2.0 * h;       y_max = 3.0 * h-0.05;
    comp_mask = (X >= x_min) & (X <= x_max) & ...
                (Y >= y_min) & (Y <= y_max) & (psi > 0.5);

    if ~any(comp_mask(:))
        warning('没有找到抽水区域！质量将不再守恒。');
    else
        vol_comp = sum(psi(comp_mask)) * dx * dy;
        phi_temp(comp_mask) = phi_temp(comp_mask) - (mass_gained / vol_comp);
        phi_temp(comp_mask) = max(-1, phi_temp(comp_mask));
    end

    %% 8. 检查质量
    phi_out  = phi_temp;
    C0_after = sum(phi_out(:) .* psi(:)) * dx * dy;
    fprintf('注入节点数 = %d\n', sum(target_mask(:)));
    fprintf('抽水节点数 = %d\n', sum(comp_mask(:)));
    fprintf('质量变化: %.6f → %.6f (差值: %.2e)\n', C0_before, C0_after, C0_after - C0_before);

    %% 9. 可视化
    figure('Name', sprintf('Gap %d 部分侵入扰动', gap_index), 'Color', 'w');
    pcolor(X, Y, phi_out);
    shading interp;
    colormap(parula);
    caxis([-1, 1]);
    colorbar;
    hold on;
    for i = 1:n_pillars
        fill([pillar_lefts(i), pillar_rights(i), pillar_rights(i), pillar_lefts(i)], ...
             [0, 0, h, h], [0.7 0.7 0.7], 'EdgeColor', 'none');
    end
    plot([x_min, x_max, x_max, x_min, x_min], [y_min, y_min, y_max, y_max, y_min], ...
         'r-', 'LineWidth', 1.5);
    axis equal;
    axis([0, 1, 0, 0.8]);
    axis off;
    title(sprintf('间隙 Gap %d 上半段注入液体，红框为抽水区', gap_index));
    hold off;

    %% 10. 转回 1D
    phi_out_1D = phi_out(:);
end