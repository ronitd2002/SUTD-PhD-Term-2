function dydt = ode_system(t,y,c1,c2,c3,b,g)

z = y(1);
v = y(2);

acc = c1*exp(-b*z) ...
      - sign(v)*c2*(c3*abs(v)^(1/2) + 0.5407*exp(-b*z/2)*abs(v)^2) ...
      - g;

dydt = [v; acc];

end