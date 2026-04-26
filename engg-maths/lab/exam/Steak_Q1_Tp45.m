function Steak_Q1_Tp45

% Parameters
kappa = 0.1;
alpha = 0.05;
beta  = 2;
gamma = 3;
delta = 0.9;
T0    = 5;
Tp    = 45;
Ti    = 10;

% Mesh
z = linspace(0,2,50);
t = linspace(0,5,1000);   % suitable time mesh

% Symmetry parameter (Cartesian 1D)
m = 0;

% Solve PDE
sol = pdepe(m, ...
    @(z,t,u,DuDz) pdefun(z,t,u,DuDz,kappa), ...
    @(z) icfun(z,T0), ...
    @(zl,ul,zr,ur,t) bcfun(zl,ul,zr,ur,t,Tp,alpha,kappa,gamma,delta), ...
    z, t);

u = sol;   % temperature matrix, rows=time, cols=space

% Find first time when minimum temperature reaches Ti
minT = min(u,[],2);                  % minimum over space for each time
idx  = find(minT >= Ti, 1, 'first'); % first time index

if isempty(idx)
    disp('Target temperature not reached in chosen time interval.');
else
    t_target = t(idx);
    fprintf('For Tp = %g, target Ti = %g is reached at t = %.4f\n', Tp, Ti, t_target);

    % Plot profile at target time
    figure
    plot(z, u(idx,:), 'LineWidth', 2)
    xlabel('z')
    ylabel('Temperature T(z,t)')
    title(['Temperature profile at t = ', num2str(t_target), ' for T_p = ', num2str(Tp)])
    grid on
end

% Surface plot
figure
surf(z,t,u)
xlabel('z')
ylabel('t')
zlabel('T')
title(['Temperature evolution for T_p = ', num2str(Tp)])

end

% PDE function
function [c,f,s] = pdefun(~,~,~,DuDz,kappa)
c = 1;
f = kappa*DuDz;
s = 0;
end

% Initial condition
function u0 = icfun(~,T0)
u0 = T0;
end

% Boundary conditions
function [pl,ql,pr,qr] = bcfun(~,ul,~,ur,~,Tp,alpha,kappa,gamma,delta)

% Left boundary: T(0,t) = Tp
pl = ul - Tp;
ql = 0;

% Right boundary:
% -alpha*dT/dz = gamma*(T(beta,t) - delta*Tp)
% Since f = kappa*dT/dz:
% gamma*(ur - delta*Tp) + (alpha/kappa)*f = 0
pr = gamma*(ur - delta*Tp);
qr = alpha/kappa;

end