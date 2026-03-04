% Parameters
r = 0.5;  % Growth rate
y0 = 1;   % Initial condition`
tspan = [0 10];  % Time span for the solution
tspan_analyt = linspace(0, 10, 1000);

% Function
growthODE = @(t, y) r * y;

% Solve
[t, y] = ode45(growthODE, tspan, y0);

% Analytical solution
y_analyt = @(t) exp(r * t);

% Solve the analytical solution using tspan_analyt
y_sol = y_analyt(tspan_analyt);

% Plot
figure;
plot(t, y, 'o', 'LineWidth', 1); % Numerical solution as discrete points
hold on;
plot(tspan_analyt, y_sol, 'b-', 'LineWidth', 1); % Analytical solution as a continuous line
hold off;
xlabel('Time (t)');
ylabel('y(t)');
title('Exponential Growth Solution using ode45');
legend('numerical', 'analytical', 'Location', 'southeast');
grid on;