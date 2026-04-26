function SpherePDE_allDeltas

% -----------------------------
% Parameters
% -----------------------------
alpha = 1.0;     % diffusion coefficient in PDE
beta  = 0.0;     % source constant
gamma = 0.0;     % nonlinear coefficient multiplying T^2
k     = 1.0;     % conductivity coefficient in boundary condition

r0 = 1.0;        % outer radius
T0 = 0.0;        % initial temperature at r = 0
T1 = 1.0;        % reference / ambient temperature at r = r0

deltas = [1 5 10];   % values of delta to test

% -----------------------------
% Mesh
% -----------------------------
r = linspace(0,r0,50);
t = linspace(0,1,50);

% -----------------------------
% Symmetry parameter for sphere
% -----------------------------
m = 2;

% -----------------------------
% Loop over deltas
% -----------------------------
for j = 1:length(deltas)

    delta = deltas(j);

    sol = pdepe(m, ...
        @(r,t,u,DuDr) pdefun(r,t,u,DuDr,alpha,beta,gamma), ...
        @(r) icfun(r,T0,T1,r0), ...
        @(rl,ul,rr,ur,t) bcfun(rl,ul,rr,ur,t,k,delta,T1), ...
        r, t);

    u = sol;

    % 3D surface plot
    figure
    surf(r,t,u)
    xlabel('r')
    ylabel('t')
    zlabel('T')
    title(['Sphere PDE solution for \delta = ', num2str(delta)])

    % Final-time profile
    figure
    plot(r,u(end,:),'LineWidth',2)
    xlabel('r')
    ylabel('T(r,1)')
    title(['Temperature profile at final time for \delta = ', num2str(delta)])
    grid on

end

end

% ----------------------------------
% PDE function
% c*u_t = r^{-m} d/dr ( r^m f ) + s
% ----------------------------------
function [c,f,s] = pdefun(~,~,u,DuDr,alpha,beta,gamma)
c = 1;
f = alpha*DuDr;
s = beta + gamma*u.^2;
end

% ----------------------------------
% Initial condition
% T(r,0) = T0 + (r/r0)(T1-T0)
% ----------------------------------
function u0 = icfun(r,T0,T1,r0)
u0 = T0 + (r/r0)*(T1 - T0);
end

% ----------------------------------
% Boundary conditions
% At r=0: dT/dr = 0
% At r=r0: -k dT/dr = delta*(T-T1)
%
% p + q*f = 0
% ----------------------------------
function [pl,ql,pr,qr] = bcfun(~,~,~,ur,~,k,delta,T1)

% Left boundary: r = 0
% dT/dr = 0  --> f = alpha*dT/dr = 0
pl = 0;
ql = 1;

% Right boundary: r = r0
% -k*dT/dr = delta*(ur - T1)
%
% Since f = alpha*dT/dr,
% dT/dr = f/alpha
% so one exact form would be:
% -(k/alpha)f - delta*(ur-T1) = 0
%
% If alpha = 1, this simplifies directly.
pr = delta*(ur - T1);
qr = k;

end