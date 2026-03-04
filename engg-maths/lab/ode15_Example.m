% Parameters
k = 1e4;      
y0 = 0;         
tspan = [0 0.02]; 
tspan_analyt = linspace(0, 0.02, 1000);

% Function
stiffODE = @(t, y) -k * y + 3000;

% Use ode15s to solve the ODE
[t_ode15s, y_ode15s] = ode15s(stiffODE, tspan, y0);

% Analytical solution
y_analyt = @(t) 0.3*(1-exp(-10000 * t));
y_sol = y_analyt(t_ode15s);

% Plot the results
figure;
plot(t_ode15s, y_ode15s, 'o', 'LineWidth', 1);
hold on;
plot(t_ode15s, y_sol, '-r', 'LineWidth', 1);
xlabel('Time (t)');
ylabel('y(t)');
legend('numerical', 'analytical', 'Location', 'southeast');
grid on;
