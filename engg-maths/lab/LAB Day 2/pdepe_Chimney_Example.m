% Parameters
z = linspace(0, 1, 1e3); % discretize space
t = linspace(0, 1e3, 1e3); % discretize time
m = 0;

% Solve PDE
sol = pdepe(m, @pdefun, @icfun, @bcfun, z, t);

% Extract solution
u = sol;

% Plot solution
figure(1)
surf(z,t,u)
xlabel('Chimney position/ m'); % Label for x-axis
ylabel('Time/ s'); % Label for y-axis
zlabel('CO2 concentration/ mol/m3'); % Label for z-axis
shading interp
colorbar

% Find the closest point to z = 0.5
desired_z = 0.5; % Desired spatial position
[~, z_idx] = min(abs(z - desired_z)); % Index of closest point

% Extract the solution at the closest z
z_sol = sol(:, z_idx); % Extract solution for all time at closest z

% Plot the solution over time
figure(2)
plot(t, z_sol, 'LineWidth', 2);
xlabel('Time/ s');
ylabel('CO2 xoncentration/ mol/m3');
title(['CO2 concentration over time at z = ', num2str(z(z_idx))]);
grid on;

function [c,f,s] = pdefun(z,t,u,DuDz)
    D = 1e-3; % Define diffusion coefficient 
    c = 1;
    f = D*DuDz;
    s = 0;
end

function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
    pl = ul - 0.4;
    ql = 0;
    pr = ur - 0.016;
    qr = 0;
end

function [u0] = icfun(x);
    u0 = 0;
end