function SampleTestQ2_allDeltas

m = 1;

x = linspace(0,1,10000);
t = linspace(0,1,10000);

deltas = [1 5 10];

for k = 1:length(deltas)
    
    delta = deltas(k);
    
    sol = pdepe(m,@pdefun,@icfun,@(xl,ul,xr,ur,t) bcfun(xl,ul,xr,ur,t,delta),x,t);
    u = sol;
    
    % Surface plot
    figure
    surf(x,t,u)
    title(['Surface plot for \delta = ', num2str(delta)])
    xlabel('x')
    ylabel('t')
    zlabel('T')
    
    % Profile at final time
    figure
    plot(x,u(end,:),'LineWidth',2)
    title(['T vs x at t=1 for \delta = ', num2str(delta)])
    xlabel('x')
    ylabel('T')
    
end

% PDE
function [c,f,s] = pdefun(~,~,~,DuDx)
c = 1;
f = DuDx;
s = 0;

% IC
function u0 = icfun(~)
u0 = 1;

% BC
function [pl,ql,pr,qr] = bcfun(~,~,~,ur,~,delta)

% x = 0
pl = 0;
ql = 1;

% x = 1
pr = delta*ur;
qr = 1;