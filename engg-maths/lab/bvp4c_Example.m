% Define the ODE as a function handle
ode = @(x, y) [y(2); -y(1)]; % System of first-order equations

% Define the boundary conditions as a function handle
bc = @(ya, yb) [ya(1) - 1; yb(1) - 1]; % Now y(0.5*pi) = 1

% Initial guess for the solution
solinit = bvpinit(linspace(0, 0.5 * pi, 100), [1, 0]); 

% Solve the BVP using bvp4c
sol = bvp4c(ode, bc, solinit);

% Evaluate the solution at desired points
x = linspace(0, 0.5 * pi, 100); % Updated interval for evaluation
x_analyt = linspace(0, 0.5 * pi, 10);
y = deval(sol, x);

% Analytical solution
y_analyt = @(x) cos(x) + sin(x);

% Solve the analytical solution using tspan_analyt:
y_sol = y_analyt(x_analyt);

% Plot the solution
figure;
plot(x, y(1, :));
hold on;
plot(x_analyt, y_sol, 'o');
hold off;
title('Solution to the BVP with y(0.5*pi) = 1');
xlabel('x');
ylabel('y(x)');
title('Exponential Growth Solution using ode45');
legend('numerical', 'analytical', 'Location', 'southeast');
grid on;