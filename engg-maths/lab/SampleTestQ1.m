clear all
clc
t0=0;
tb=3500; 
z0=[0,0]; 
options=odeset('RelTol',1e-5,'AbsTol',[1e-7 1e-7]);
[t,z1]=ode45(@sample1,[t0,tb],z0,options);
figure(1)
plot(z1(:,1),z1(:,2)); 
ylabel('dz/dt');
xlabel('t');
title('dz/dt vs z');
grid on

figure(2)
plot(t,z1(:,1)); 

ylabel('z');
xlabel('t');
title('z vs t');
grid on

figure(3)
plot(t,z1(:,2)); 
ylabel('dz/dt');
xlabel('t');
title('dz/dt vs t');
grid on

function dydx=sample1(t,z)
c1=26.525;
c2=6.19e-2;
c3=3.343e-3;
b=11.58e-5;
g=9.807;
dydx=zeros(2,1);
dydx(1)=z(2);
dydx(2)=c1*exp(-b*z(1))-sign(z(2))*c2*(c3*sqrt(abs(z(2)))+0.5407*exp((-b*z(1))/2)*abs(z(2))^2)-g;
end