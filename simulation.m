clear
clc

% constants
c1 = 26.525;
c2 = 0.0619;
c3 = 3.343e-3;
b  = 11.58e-5;
g  = 9.807;

% time range
tspan = [0 3500];

% initial conditions
y0 = [0;0];   % [z , dz/dt]

[t,y] = ode45(@(t,y) ode_system(t,y,c1,c2,c3,b,g), tspan, y0);

z = y(:,1);
v = y(:,2);

% Plot 1: z vs t
figure
plot(t,z,'LineWidth',2)
xlabel('t')
ylabel('z')
title('z vs t')
grid on
saveas(gcf,'z_vs_t.png')

% Plot 2: dz/dt vs t
figure
plot(t,v,'LineWidth',2)
xlabel('t')
ylabel('dz/dt')
title('dz/dt vs t')
grid on
saveas(gcf,'velocity_vs_t.png')

% Plot 3: dz/dt vs z
figure
plot(z,v,'LineWidth',2)
xlabel('z')
ylabel('dz/dt')
title('dz/dt vs z')
grid on
saveas(gcf,'velocity_vs_z.png')