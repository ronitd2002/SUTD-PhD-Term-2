function Q1_steak_PDE()
%% ============================================================
%% Question 1: PDE Problem — Cooking the Perfect Steak
%% 99.500 Applied Mathematics for Engineering, April 2026
%%
%% HOW TO RUN:
%%   Simply type  Q1_steak_PDE  in the MATLAB command window,
%%   or press the Run button. No arguments needed.
%%
%% TESTED ON: MATLAB R2025b
%% ============================================================
%
% PDE:  dT/dt = kappa * d2T/dz2,   z in [0, beta],  t > 0
%
% Boundary conditions:
%   (BC1)  T(0,t) = Tp                              [Dirichlet, pan side]
%   (BC2) -alpha * dT/dz|_{z=beta}
%              = gamma*(T(beta,t) - delta*Tp)       [Robin,     top side]
%
% Initial condition:
%   T(z,0) = T0  (uniform)
%
% Numerical Method: Method of Lines (MOL)
%   - Space : linspace(0,2,50) as specified in the question
%   - Time  : ode15s  (stiff, variable-step, high-order solver)
%   - Robin BC algebraically eliminated inside each RHS evaluation
%
% ============================================================

clc; close all;

%% ---- Physical parameters ----
kappa = 0.1;    % thermal diffusivity
alpha = 0.05;   % material conductivity
beta  = 2;      % steak thickness
gam   = 3;      % Robin BC coefficient  ('gam' avoids clash with built-in gamma())
delta = 0.9;    % proportionality coefficient
T0    = 5;      % uniform initial temperature
Ti    = 10;     % target internal (centre) temperature

%% ---- Spatial grid (as specified in question) ----
z  = linspace(0, beta, 50);   % 50 points in [0, 2]
N  = length(z);                % = 50
dz = z(2) - z(1);             % uniform spacing = 2/49

%% ---- Pan-side temperatures to investigate ----
Tp_list  = [15, 25, 35, 45];
clr_list = {'b', 'r', [0.1 0.7 0.1], 'm'};

cook_times = NaN(1, numel(Tp_list));
T_cook     = cell(1, numel(Tp_list));   % temperature profile at cook time

%% ============================================================
%% Solve PDE for each Tp (Parts a and b)
%% ============================================================
for k = 1:numel(Tp_list)

    Tp = Tp_list(k);

    % Initial condition
    T_init    = T0 * ones(N, 1);
    T_init(1) = Tp;     % apply Dirichlet BC at t = 0

    % ODE solver options
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.2);

    % Integrate.  500 time units is more than sufficient for all Tp.
    [t_sol, T_sol] = ode15s( ...
        @(t, T) steak_rhs(t, T, kappa, dz, N, Tp, alpha, gam, delta), ...
        [0, 500], T_init, opts);

    % Re-apply BCs to stored solution (ode15s does not store algebraic values)
    T_sol(:, 1) = Tp;
    for idx = 1:numel(t_sol)
        T_sol(idx, N) = robin_TN(T_sol(idx, N-1), Tp, dz, alpha, gam, delta);
    end

    % Minimum temperature across the spatial domain at each time step
    Tmin_t = min(T_sol, [], 2);

    % First time index where the coldest point reaches Ti
    idx_cross = find(Tmin_t >= Ti, 1, 'first');

    if isempty(idx_cross)
        fprintf('Tp = %2d : Target Ti = %g not reached within t = 500\n', Tp, Ti);
        continue
    end

    % Linear interpolation between the two bracketing time steps
    t1 = t_sol(idx_cross-1);   Tm1 = Tmin_t(idx_cross-1);
    t2 = t_sol(idx_cross);     Tm2 = Tmin_t(idx_cross);
    t_star = t1 + (Ti - Tm1)/(Tm2 - Tm1) * (t2 - t1);

    cook_times(k) = t_star;
    T_cook{k}     = T_sol(idx_cross, :);   % profile at crossing index
end

%% ============================================================
%% Part (a): Print cooking time table
%% ============================================================
fprintf('\n');
fprintf('=====================================================\n');
fprintf('  Part (a): Cooking Time Table  (Ti = %g)\n', Ti);
fprintf('=====================================================\n');
fprintf('  %-8s  %-20s\n', 'Tp', 'Cook time  t*');
fprintf('  %s\n', repmat('-', 1, 32));
for k = 1:numel(Tp_list)
    if isnan(cook_times(k))
        fprintf('  %-8d  %-20s\n', Tp_list(k), 'Not reached');
    else
        fprintf('  %-8d  %-20.4f\n', Tp_list(k), cook_times(k));
    end
end
fprintf('=====================================================\n\n');

%% ============================================================
%% Part (b): All temperature profiles on one graph
%% ============================================================
figure('Name', 'Part (b): Temperature Profiles at Cook Time', ...
       'NumberTitle', 'off', 'Position', [50 50 950 620]);
hold on;

for k = 1:numel(Tp_list)
    if ~isnan(cook_times(k))
        plot(z, T_cook{k}, 'Color', clr_list{k}, 'LineWidth', 2.5, ...
             'DisplayName', sprintf('T_p = %d   (t* = %.2f)', ...
             Tp_list(k), cook_times(k)));
    end
end

yline(Ti, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Target  T_i = 10');

xlabel('z  (distance through steak;  z = 0: pan side)', 'FontSize', 13);
ylabel('Temperature  T(z, t*)', 'FontSize', 13);
title({'Part (b) — Temperature Profile at the Moment T_{min} = T_i = 10', ...
       sprintf('(kappa=%.1f, alpha=%.2f, beta=%g, gamma=%g, delta=%.1f, T0=%g)', ...
               kappa, alpha, beta, gam, delta, T0)}, 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 12);
grid on; box on;

saveas(gcf, 'Q1b_temperature_profiles.png');
fprintf('Part (b) figure saved: Q1b_temperature_profiles.png\n\n');

%% ============================================================
%% Part (c): Top 3 Key Findings
%% ============================================================
fprintf('=====================================================\n');
fprintf('  Part (c): Top 3 Key Findings\n');
fprintf('=====================================================\n\n');

fprintf('  FINDING 1 — Higher Tp drastically shortens cook time\n');
fprintf('  -------------------------------------------------------\n');
fprintf('  Cook time decreases strongly and nonlinearly with Tp.\n');
fprintf('  Pan temperature is the dominant control variable.\n\n');

fprintf('  FINDING 2 — Steep temperature gradient always at pan side\n');
fprintf('  -----------------------------------------------------------\n');
fprintf('  At cook time, z=0 is locked at Tp while the centre\n');
fprintf('  is just reaching Ti=10. Higher Tp creates a steeper\n');
fprintf('  gradient, confirming the Dirichlet BC drives heating.\n\n');

fprintf('  FINDING 3 — Top surface is warmer than the centre\n');
fprintf('  ---------------------------------------------------\n');
fprintf('  With delta=0.9, the Robin BC drives the top surface\n');
fprintf('  toward ~0.9*Tp (always above Ti by cook time).\n');
fprintf('  The steak has an asymmetric doneness profile:\n');
fprintf('  hottest at pan side, intermediate at top, coldest inside.\n\n');
fprintf('=====================================================\n\n');

%% ============================================================
%% Part (d): What happens if heat is not turned off?
%% ============================================================
fprintf('=====================================================\n');
fprintf('  Part (d): Continued heating beyond T_i\n');
fprintf('=====================================================\n\n');
fprintf('  All temperatures continue rising toward steady state:\n');
fprintf('      kappa * d2T/dz2 = 0  =>  T_ss(z) linear in z\n');
fprintf('  with T_ss(0)=Tp and Robin condition fixing T_ss(beta).\n');
fprintf('  In practice: the pan-side exterior burns/chars while\n');
fprintf('  the interior becomes dry and overcooked. The longer\n');
fprintf('  past t* the steak is left, the worse the gradient:\n');
fprintf('  exterior far above safe temperature, interior still\n');
fprintf('  rising toward steady state.\n\n');
fprintf('=====================================================\n\n');

%% ============================================================
%% Part (e): Complete standalone code for Tp = 45
%% ============================================================
fprintf('=====================================================\n');
fprintf('  Part (e): Standalone solution for Tp = 45\n');
fprintf('=====================================================\n');

% Parameters
kappa_e = 0.1;  alpha_e = 0.05;  beta_e  = 2;
gam_e   = 3;    delta_e = 0.9;   T0_e    = 5;   Ti_e = 10;
Tp_e    = 45;

% Spatial grid
z_e  = linspace(0, beta_e, 50);
N_e  = length(z_e);
dz_e = z_e(2) - z_e(1);

% Initial condition
T_init_e    = T0_e * ones(N_e, 1);
T_init_e(1) = Tp_e;

% Solve
opts_e = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.2);
[t_e, T_e] = ode15s( ...
    @(t, T) steak_rhs(t, T, kappa_e, dz_e, N_e, Tp_e, alpha_e, gam_e, delta_e), ...
    [0, 100], T_init_e, opts_e);

% Post-process BCs
T_e(:, 1) = Tp_e;
for ii = 1:numel(t_e)
    T_e(ii, N_e) = robin_TN(T_e(ii, N_e-1), Tp_e, dz_e, alpha_e, gam_e, delta_e);
end

% Find cook time
Tmin_e  = min(T_e, [], 2);
idx_e   = find(Tmin_e >= Ti_e, 1, 'first');
t1_e    = t_e(idx_e-1);   Tm1_e = Tmin_e(idx_e-1);
t2_e    = t_e(idx_e);     Tm2_e = Tmin_e(idx_e);
tstar_e = t1_e + (Ti_e - Tm1_e)/(Tm2_e - Tm1_e) * (t2_e - t1_e);

fprintf('\n  Tp = 45 : Cook time  t* = %.4f\n\n', tstar_e);

% Figure for Part (e)
figure('Name', 'Part (e): Tp = 45 Full Solution', ...
       'NumberTitle', 'off', 'Position', [80 80 1100 480]);

subplot(1, 2, 1);
contourf(z_e, t_e, T_e, 40);
colorbar;
hold on;
yline(tstar_e, 'w--', 'LineWidth', 2.0);
text(0.05, tstar_e * 1.05, sprintf(' t* = %.2f', tstar_e), ...
     'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('z  (position)', 'FontSize', 12);
ylabel('Time  t',       'FontSize', 12);
title('T(z,t) Space-Time Contour  (T_p = 45)', 'FontSize', 12);

subplot(1, 2, 2);
plot(z_e, T_e(idx_e, :), 'r', 'LineWidth', 2.5, ...
     'DisplayName', sprintf('T(z, t* = %.2f)', tstar_e));
yline(Ti_e, 'k--', 'LineWidth', 1.8, 'DisplayName', 'T_i = 10');
xlabel('z  (position)', 'FontSize', 12);
ylabel('Temperature  T', 'FontSize', 12);
title(sprintf('Profile at Cook Time  t* = %.2f   (T_p = 45)', tstar_e), ...
      'FontSize', 12);
legend('Location', 'southeast', 'FontSize', 11);
grid on; box on;

saveas(gcf, 'Q1e_Tp45.png');
fprintf('Part (e) figures saved: Q1e_Tp45.png\n');

end   % <<<  end of main function Q1_steak_PDE


%% ============================================================
%% LOCAL HELPER FUNCTIONS  (must be after the main function)
%% ============================================================

function dTdt = steak_rhs(~, T, kappa, dz, N, Tp, alpha, gam, delta)
% Method-of-Lines RHS for the 1D steak heat conduction PDE.
%
% Discretisation:
%   Node 1  (z = 0)    : Dirichlet BC  ->  T(1) = Tp  (no ODE update)
%   Node N  (z = beta) : Robin BC solved algebraically via robin_TN()
%   Nodes 2..N-1       : Central 2nd-order FD for d2T/dz2

    dTdt = zeros(N, 1);

    % Enforce Dirichlet at pan side
    T(1) = Tp;

    % Enforce Robin at top surface (algebraic expression)
    T(N) = robin_TN(T(N-1), Tp, dz, alpha, gam, delta);

    % Interior nodes: standard central finite difference
    for j = 2:N-1
        dTdt(j) = kappa * (T(j+1) - 2*T(j) + T(j-1)) / dz^2;
    end

    % Boundary nodes are algebraically constrained — zero time derivatives
    dTdt(1) = 0;
    dTdt(N) = 0;
end


function TN = robin_TN(TN1, Tp, dz, alpha, gam, delta)
% Solves for T at the top-surface node (z = beta) from the Robin BC:
%
%   -alpha * (T_N - T_{N-1}) / dz  =  gamma * (T_N - delta * Tp)
%
% Rearranging for T_N:
%   -alpha*T_N/dz + alpha*T_{N-1}/dz  =  gamma*T_N - gamma*delta*Tp
%   T_N * (alpha/dz + gamma)           =  gamma*delta*Tp + alpha*T_{N-1}/dz
%
%   => T_N = (gamma*delta*Tp + (alpha/dz)*T_{N-1}) / (alpha/dz + gamma)

    TN = (gam * delta * Tp + (alpha / dz) * TN1) / (alpha / dz + gam);
end