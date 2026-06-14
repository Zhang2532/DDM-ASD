%% MAIN SCRIPT FOR DDM-ASD SADDLE POINT SEARCH
% This script serves as the main driver for executing the Diffuse-Domain
% Accelerated Saddle Dynamics (DDM-ASD) algorithm to find an index-k
% saddle point for the wetting transition problem.
%
% Workflow:
% 1. Define geometry and initial droplet state.
% 2. Set simulation parameters (saddle index, dt, momentum, etc.).
% 3. Apply a perturbation to create an initial guess.
% 4. Run the main DDM-ASD iteration loop.
% 5. Monitor convergence and physical properties in real-time.
% 6. Compute Morse index of the final state.

clc; clear;

%% 1. CHOOSE GEOMETRY AND INITIALIZE DROPLET
% Uncomment the desired geometry setup function.
[u, X, Y, psi, abs_grad_psi, dx, dy, params] = droplet_rectangular_domain_pillar();
% [u, X, Y, psi, abs_grad_psi, dx, dy, params] = droplet_rectangular_domain_pillar_circle();
% [u, X, Y, psi, abs_grad_psi, dx, dy, params] = droplet_rectangular_domain_pillar_tilt();

% --- Option to load a previously saved state ---
use_saved = 0;
if use_saved && isfile('path/to/your/phi_result.mat')
    load('path/to/your/phi_result.mat', 'phi_p');
    u = reshape(phi_p, size(u,1), size(u,2));
end

%% 2. SET SIMULATION PARAMETERS
k         = 1;         % Index of the saddle point to search
l         = 0.5;       % Initial step size for finite difference Hessian
phi_initial = u(:);
dt        = 0.0001;    % Time step size (beta in the paper)
gm        = 0.1;       % Momentum parameter (gamma in the paper)
epsilonf  = 1e-6;      % Gradient convergence tolerance
epsilon_l = 1e-6;      % Lower bound for 'l'
maxstep   = 10^6;      % Maximum number of iterations

% --- Physical Parameters ---
epsilon = 0.015;
theta   = deg2rad(102);

% --- Initialization ---
n = 0;
phi_p = phi_initial;
v = [];
c0_mass = sum(u(:) .* psi(:)) * dx * dy;
psi_flat = psi(:);
psi_mask = (psi_flat == 1);

% --- History Logging ---
grad_history   = [];
energy_history = [];
mass_history   = [];
mass_diff_history = [];
step_history   = [];
mass_prev = c0_mass;

%% 3. CREATE INITIAL GUESS VIA PERTURBATION
% This step "kicks" the system away from a local minimum towards a saddle.
phi_p = perturb_imbibition_FD(phi_initial, X, Y, psi, dx, dy, 5, params);

%% 4. COMPUTE INITIAL EIGENVECTORS (UNSTABLE DIRECTIONS)
N_total = numel(phi_p);
[~, ~, H_op] = compute_DDM_physics(phi_p, psi, abs_grad_psi, dx, dy, epsilon, theta, c0_mass, dt);
opts.issym = true;
opts.maxit = 1000;
opts.tol   = 1e-8;
if k > 0
    [V_init, ~] = eigs(H_op, N_total, k, 'smallestreal', opts);
    v = V_init;
end

% --- Get initial gradient ---
[~, f] = compute_DDM_physics_1(phi_p, psi, abs_grad_psi, dx, dy, epsilon, theta);

%% 5. MAIN DDM-ASD ITERATION LOOP
fprintf('Starting DDM-ASD search for an index-%d saddle point...\n', k);
tic;
while norm(f) > epsilonf
    phi_old = phi_p;

    % --- Step 1: Compute search direction g (Eq. 2.23) ---
    g = f;
    if k > 0
        g = f - 2 * v * (v' * f);
    end

    % --- Step 2: Update phase-field phi_p (Eq. 2.19) ---
    if n == 0
        phi_p = phi_p + dt * g;
    else
        phi_p = phi_p + dt * g + gm * (phi_p - phi_old1);
    end
    phi_old1 = phi_old;
    
    v_old = v;

    % --- Step 3: Update eigenvectors v (Eq. 2.19 & 2.24) ---
    for i = 1:k
        % 3a. Approximate Hessian-vector product H*v_i using finite differences
        [~, F_plus]  = compute_DDM_physics_1(phi_p + l * v_old(:,i), psi, abs_grad_psi, dx, dy, epsilon, theta);
        [~, F_minus] = compute_DDM_physics_1(phi_p - l * v_old(:,i), psi, abs_grad_psi, dx, dy, epsilon, theta);
        u_i = -(F_plus - F_minus) / (2 * l); % This is H(phi, v_i)

        % 3b. Compute eigenvector update direction d_i
        d_i = -u_i + (v(:, i)' * u_i) * v(:, i);
        if i > 1
             d_i = d_i + 2 * v(:, 1:i-1) * (v(:, 1:i-1)' * u_i);
        end

        % 3c. Update v_i
        v(:, i) = v_old(:, i) + dt * d_i;

        % 3d. Modified Gram-Schmidt Orthonormalization
        if i > 1
            coeffs = v(:, 1:i-1)' * v(:, i);
            v(:, i) = v(:, i) - v(:, 1:i-1) * coeffs;
        end
        v(:, i) = v(:, i) / norm(v(:,i));
    end

    % --- Step 4: Update finite-difference step 'l' and compute new gradient 'f' ---
    l = max(l / (1 + dt), epsilon_l);
    [~, f] = compute_DDM_physics_1(phi_p, psi, abs_grad_psi, dx, dy, epsilon, theta);
    
    n = n + 1;

    % --- Monitoring & Visualization (every 200 steps) ---
    if mod(n, 200) == 0 || n == 1
        [E_cur, ~] = compute_DDM_physics_1(phi_p, psi, abs_grad_psi, dx, dy, epsilon, theta);
        mass_cur = sum(phi_p(:) .* psi(:)) * dx * dy;
        mass_diff = mass_cur - mass_prev;
        mass_prev = mass_cur;
        
        grad_history(end+1)   = norm(f);
        energy_history(end+1) = E_cur;
        mass_history(end+1)   = mass_cur;
        mass_diff_history(end+1) = mass_diff;
        step_history(end+1)   = n;

        fprintf('Step: %d, |g|: %.3e, E: %.4f, Mass Err: %.2e\n', n, norm(f), E_cur, mass_diff);
        % Add plotting routines here if desired...
    end

    if n >= maxstep
        fprintf('Maximum number of iterations reached.\n');
        break;
    end
end
toc;

%% 6. FINAL ANALYSIS: COMPUTE MORSE INDEX
fprintf('Convergence reached. Final gradient norm: %.3e\n', norm(f));
fprintf('Computing Morse index of the final state...\n');
[~, ~, H_op] = compute_DDM_physics(phi_p, psi, abs_grad_psi, dx, dy, epsilon, theta, c0_mass, dt);
n_eig = min(10, N_total - 1);
[~, D] = eigs(H_op, N_total, n_eig, 'smallestreal', opts);
eig_vals = sort(diag(D));
morse_index = sum(eig_vals < -1e-10);

fprintf('Morse Index = %d\n', morse_index);
fprintf('Smallest eigenvalues: '); fprintf('%.4e  ', eig_vals(1:min(5,end))'); fprintf('\n');

%% 7. FINAL VISUALIZATION
plot_phase_field(X, Y, phi_p, 'pillar', true, 'params', params);
title(sprintf('Final State (Index-%d Saddle), Step %d', k, n));