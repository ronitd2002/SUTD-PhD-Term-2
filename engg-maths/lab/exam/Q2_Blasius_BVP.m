function Q2_Blasius_BVP()

clc; close all;

%% ---- Problem setup ----
L_values = [5, 20, 35, 50];
nL       = numel(L_values);

% Physical reference
FPP0_REF = 0.332057;

fprintf('\n=== Blasius BVP Results ===\n\n');
fprintf('  %-4s  %-8s  %-12s  %-10s  %-10s\n', ...
        'L', 'Guess', 'Status', 'maxres', 'f''''(0)');
fprintf('  %s\n', repmat('-', 1, 55));

figure('Name','Blasius Solutions','Position',[50 50 1200 1000]);

%% ============================================================
%% Loop over all cases
%% ============================================================

for k = 1:nL
    L      = L_values(k);
    N_mesh = 100*L + 1;
    eta0   = linspace(0, L, N_mesh);

    for g = 1:2

        %% ---- Initial guess ----
        if g == 1
            solinit = bvpinit(eta0, [0;0;0]);
            glabel  = '0';
        else
            solinit = bvpinit(eta0, @eta_guess);
            glabel  = 'eta';
        end

        %% ---- Solver options ----
        opts = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',10000);

        %% ---- Solve safely ----
        warnState = warning('off','MATLAB:bvp4c:RelTolNotMet');
        try
            sol = bvp4c(@blasius_ode, @blasius_bc, solinit, opts);
            success = true;
        catch
            sol = [];
            success = false;
        end
        warning(warnState);

        %% ---- Evaluate solution ----
        if ~success
            fprintf('  %-4d  %-8s  %-12s\n', L, glabel, 'FAILED');
            continue;
        end

        eta_plot = linspace(0, L, max(1000, N_mesh));
        Y        = deval(sol, eta_plot);

        f   = Y(1,:);
        fp  = Y(2,:);
        fpp = Y(3,:);

        %% ---- Residual (approximate diagnostic) ----
        deta  = eta_plot(2)-eta_plot(1);
        fpp3  = gradient(fpp, deta);
        res   = abs(2*fpp3 + f.*fpp);
        maxres = max(res);

        fpp0 = fpp(1);

        %% ---- Physical validation ----
        is_physical = (maxres < 0.05) && ...
                      (abs(fpp0 - FPP0_REF) < 0.05) && ...
                      (fp(end) > 0.9);

        if is_physical
            status = 'OK';
            color  = [0 0.45 0.74]; % blue
        else
            status = 'NON-PHYS';
            color  = [0.85 0.15 0.1]; % red
        end

        fprintf('  %-4d  %-8s  %-12s  %-10.3e  %-10.4f\n', ...
                L, glabel, status, maxres, fpp0);

        %% ---- Plot ----
        subplot(nL,2,(k-1)*2 + g);
        plot(eta_plot, f, 'Color', color, 'LineWidth', 2);
        xlabel('\eta'); ylabel('f(\eta)');
        title(sprintf('L=%d, f0=%s (%s)', L, glabel, status));
        grid on;

        % Show full domain (your concern)
        xlim([0 L]);

    end
end

sgtitle('Blasius BVP Solutions (Blue = Physical, Red = Failed)','FontWeight','bold');

%% ============================================================
%% PART (c): Clean solution (L=50, good guess)
%% ============================================================

fprintf('\n=== Part (c): L=50, f0=eta ===\n');

L = 50;
eta = linspace(0,L,100*L+1);

solinit = bvpinit(eta, @eta_guess);

opts = bvpset('RelTol',1e-6,'AbsTol',1e-8,'NMax',10000);

warnState = warning('off','MATLAB:bvp4c:RelTolNotMet');
sol = bvp4c(@blasius_ode, @blasius_bc, solinit, opts);
warning(warnState);

eta_fine = linspace(0,L,5000);
Y = deval(sol, eta_fine);

f   = Y(1,:);
fp  = Y(2,:);
fpp = Y(3,:);

deta = eta_fine(2)-eta_fine(1);
fpp3 = gradient(fpp,deta);
maxres = max(abs(2*fpp3 + f.*fpp));

fprintf('maxres = %.3e\n', maxres);
fprintf('f''''(0) = %.6f\n', fpp(1));
fprintf('f''(L) = %.6f\n\n', fp(end));

figure;
subplot(1,3,1);
plot(eta_fine,f,'b','LineWidth',2);
title('f(\eta)'); grid on;

subplot(1,3,2);
plot(eta_fine,fp,'r','LineWidth',2); hold on;
yline(1,'k--');
title("f'(\eta)"); grid on;

subplot(1,3,3);
plot(eta_fine,fpp,'m','LineWidth',2); hold on;
yline(FPP0_REF,'k--');
title("f''(\eta)"); grid on;

sgtitle('Blasius Solution (L=50, good guess)');

end

%% ============================================================
%% ODE + BC
%% ============================================================

function dydeta = blasius_ode(~,y)
dydeta = [y(2);
          y(3);
         -0.5*y(1).*y(3)];
end

function res = blasius_bc(ya,yb)
res = [ya(1);
       ya(2);
       yb(2)-1];
end

function y0 = eta_guess(eta)
y0 = [eta; 1; 0];
end