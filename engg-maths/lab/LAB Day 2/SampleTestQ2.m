function SampleTestQ2

m = 1;

% set the lower, upper limits and the mesh number
x = linspace(0,1,20); 

% set the lower, upper limits and the mesh number
t = linspace(0,1,20); 

sol = pdepe(m,@pdefun,@icfun,@bcfun,x,t);

% Extract solution
u = sol;

% Use a surface plot to show the complete solution
surf(x,t,u)
xlabel('Coordinate x')
ylabel('Coordinate t')
% A profile at a certain time
figure
plot(x,u(end,:)) % this is to show at the final time
xlabel('Coordinate x')
ylabel('Temperature T') % can change accordingly

% define the PDE
function [c,f,s] = pdefun(x,t,u,DuDx) 
c = 1.0;
f = DuDx;
s = 0;

% define IC
function u0 = icfun(x) 
u0 = 1.0;

% define BCs
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t) 
delta = 10;
pl = 0;
ql = 1;
qr = 1;
pr=delta*ur;