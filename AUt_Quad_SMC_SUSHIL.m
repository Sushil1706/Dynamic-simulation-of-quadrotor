clc; clear all; close all;

% Simulation times, in seconds.
start_time = 0;
end_time = 25;
dt = 0.005;
times = start_time:dt:end_time;

m = 1;
I = [7.5*10^-3 0 0
    0 7.5*10^-3 0
    0 0 1.3*10^-2];
L = 0.232;
kdr = 0;

g = 9.81;
k1 = 4.905;
k2 = 0.08164;
ksmc = [20 20 1 4 1.25 1.25];
lambda = [20 20 1 4 4.25 4.25];
j = 0;


x = zeros(12,1);
xd = zeros(12,1);
xddot = zeros(12,1);
e = zeros(4001,6);
f = zeros(4001,6);

xdp = zeros(1,6);
xddotp = zeros(1,6);
omegab = zeros(3,1);
vb = zeros(3,1);

w = pi/10;

for t = times
t
j = j + 1;
xd(5) = pi/6;
xd(7) = 3;
xd(9) = cos(w*t);
xd(11) = sin(w*t);

xddot(1) = (xd(1) - xdp(1))/dt;
xddot(3) = (xd(3) - xdp(2))/dt;
xddot(5) = (xd(5) - xdp(3))/dt;
xddot(7) = (xd(7) - xdp(4))/dt;
xddot(9) = (xd(9) - xdp(5))/dt;
xddot(11) = (xd(11) - xdp(6))/dt;

xddtot(2) = (xddot(1) - xddotp(1))/dt;
xddtot(4) = (xddot(3) - xddotp(2))/dt;
xddtot(6) = (xddot(5) - xddotp(3))/dt;
xddtot(8) = (xddot(7) - xddotp(4))/dt;
xddtot(10) = (xddot(9) - xddotp(5))/dt;
xddtot(12) = (xddot(11) - xddotp(6))/dt;

data = smc(m, I, g, x, xd, xddot, lambda, ksmc);

% Ft = [(data(1)+m*g)/(cos(x(1))*cos(x(3)));
%       data(2);
%       data(3);
%       data(4)];
%   
%  M = [k1 k1 k1 k1;
%       -L*k1 L*k1 -L*k1 L*k1;
%       L*k1 -L*k1 -L*k1 L*k1;
%       -k2 -k2 k2 k2];
%  
% delta = M\Ft;
% 
% T = thrust(delta, k1);
% tau = torques(delta, L, k2, k1);


Fb = [0; 0; data(1)];
tau = [data(2); data(3); data(4)];
xd(1) = data(5);
xd(3) = data(6);
T = (Fb+[0; 0; m*g])/(cos(x(1))*cos(x(3)));


omegadotb = angular_acceleration(tau, omegab, I);
omegab = omegab + omegadotb*dt;

thetadotv = omegab2thetadotv(omegab, x);
x(2) = thetadotv(1);
x(4) = thetadotv(2);
x(6) = thetadotv(3);
x(1) = x(1) + x(2)*dt;
x(3) = x(3) + x(4)*dt;
x(5) = x(5) + x(6)*dt;

ab = acceleration(T, m, g, kdr, x);
vb = vb + ab*dt;
vi = vb2xdot(vb, x);
x(10) = vi(1);
x(12) = vi(2);
x(8) = vi(3);
x(9) = x(9) + vi(1)*dt;
x(11) = x(11) + vi(2)*dt;
x(7) = x(7) + vi(3)*dt;


xdp = [xd(1) xd(3) xd(5) xd(7) xd(9) xd(11)];
xddotp = [xddot(1) xddot(3) xddot(5) xddot(7) xddot(9) xddot(11)];



T01=[cos(x(5)), -sin(x(5)), 0, x(9);
     sin(x(5))*cos(x(1)), cos(x(5))*cos(x(1)), -sin(x(1)), -x(7)*sin(x(1));
     sin(x(5))*sin(x(1)), cos(x(5))*cos(x(1)), cos(x(1)), x(7)*cos(x(1));
     0, 0, 0, 1];
 
% T01=[cos(x(1)), -sin(x(1)), 0, x(9);
%      sin(x(1))*cos(x(5)), cos(x(1))*cos(x(5)), -sin(x(5)), -x(7)*sin(x(5));
%      sin(x(1))*sin(x(5)), cos(x(1))*cos(x(5)), cos(x(5)), x(7)*cos(x(5));
%      0, 0, 0, 1];
 T12 =[-sin(x(3)), -cos(x(3)), 0, 0;
     0, 0, 1, x(11);
     -cos(x(3)), 0, 0, 0;
     0, 0, 0, 1];
 T02=T01*T12;
 
 T23 =[1 0 0 0;
     0 1 0 0;
     0 1 1 L;
     0 0 0 1];
 
 T24 = [1 0 0 0;
     0 0 1 L;
     0 0 -1 -L;
     0 0 0 1];
 
  T25 = [1 0 0 0;
     0 -1 0 0;
     0 -1 -1 -L;
     0 0 0 1];
 
  T26 = [1 0 0 0;
     0 0 -1 -L;
     0 0 1 L;
     0 0 0 1];
 
 Ta=T02*T23;
 px1(j) = Ta(1,4);
 py1(j) = Ta(2,4);
 pz1(j) = Ta(3,4);
 
 Tb=T02*T24;
 px2(j) = Tb(1,4);
 py2(j) = Tb(2,4);
 pz2(j) = Tb(3,4);
 
 Tc=T02*T25;
 px3(j) = Tc(1,4);
 py3(j) = Tc(2,4);
 pz3(j) = Tc(3,4);
 
 Td=T02*T26;
 px4(j) = Td(1,4);
 py4(j) = Td(2,4);
 pz4(j) = Td(3,4);
 
% px1(j) = x(9) - L*sin(x(5)) - x(11)*sin(x(5));
% py1(j) = L*cos(x(1))*cos(x(5)) - x(7)*sin(x(1)) + x(11)*cos(x(1))*cos(x(5));
% pz1(j) = x(7)*cos(x(1)) + L*cos(x(1))*cos(x(5)) + x(11)*cos(x(1))*cos(x(5));
% 
% px2(j) =  x(9) + L*sin(x(5)) - x(11)*sin(x(5)) - L*cos(x(5))*cos(x(3));
% py2(j) = x(11)*cos(x(1))*cos(x(5)) - L*cos(x(1))*cos(x(5)) - x(7)*sin(x(1)) - L*cos(x(1))*cos(x(3))*sin(x(5));
% pz2(j) = x(7)*cos(x(1)) - L*cos(x(1))*cos(x(5)) + x(11)*cos(x(1))*cos(x(5)) - L*cos(x(3))*sin(x(1))*sin(x(5));
% 
% px3(j) =  x(9) + L*sin(x(5)) - x(11)*sin(x(5));
% py3(j) = x(11)*cos(x(1))*cos(x(5)) - L*cos(x(1))*cos(x(5)) - x(7)*sin(x(1));
% pz3(j) = x(7)*cos(x(1)) - L*cos(x(1))*cos(x(5)) + x(11)*cos(x(1))*cos(x(5));
% 
% px4(j) =  x(9) - L*sin(x(5)) - x(11)*sin(x(5)) + L*cos(x(5))*cos(x(3));
% py4(j) = x(11)*cos(x(1))*cos(x(5)) + L*cos(x(1))*cos(x(3))*sin(x(5));
% pz4(j) = x(11)*cos(x(1))*cos(x(5)) + L*cos(x(3))*sin(x(1))*sin(x(5));

px(j) = x(9);
pxd(j)= xd(9);
py(j) = x(11);
pyd(j) = xd(11);
pz(j) = x(7);
pzd(j) = xd(7);

e(j,1) = x(1);
e(j,2) = x(3);
e(j,3) = x(5);
e(j,4) = x(7);
e(j,5) = x(9);
e(j,6) = x(11);

f(j,1) = xd(1);
f(j,2) = xd(3);
f(j,3) = xd(5);
f(j,4) = xd(7);
f(j,5) = xd(9);
f(j,6) = xd(11);

% pause(dt);

end
for j=1:length(px)
    plot3([px1(j) px(j) px3(j)],[py1(j) py(j) py3(j)],[pz1(j) pz(j) pz3(j)], 'r')
   
    hold on
    plot3([px2(j) px(j) px4(j)],[py2(j) py(j) py4(j)],[pz2(j) pz(j) pz4(j)],'b')
    plot3(px(1:j),py(1:j),pz(1:j))
    plot3(pxd(1:j),pyd(1:j),pzd(1:j),'r')
    
    hold off
    axis([-2 2 -2 2 -1 4])
    grid on
    pause(dt)
end

figure(1)
plot(e(:,3))
hold on;
plot(f(:,3))

figure(2)
plot(e(:,4))
hold on;
plot(f(:,4))

figure(3)
plot(e(:,5))
hold on;
plot(f(:,5))

figure(4)
plot(e(:,6))
hold on;
plot(f(:,6))

% function T = thrust(inputs, k1)
% T = [0; 0; k1 * sum(inputs)];
% end
% 
% function tau = torques(inputs, L, k2, k1)
% % Inputs are values for ?i
% 
% M = [-L*k1 L*k1 -L*k1 L*k1;
%      L*k1 -L*k1 -L*k1 L*k1;
%      -k2 -k2 k2 k2];
% 
% tau = M * inputs;
% end

function ab = acceleration(T, m, g, kdr, x) 
theta = [x(1); x(3); x(5)];
gravity = [0; 0; -g];
% R = rotation(angles);
R = [cos(theta(3))*cos(theta(2)),                                              cos(theta(2))*sin(theta(3)),                                              -sin(theta(2));
     cos(theta(3))*sin(theta(1))*sin(theta(2))-cos(theta(1))*sin(theta(3)),    cos(theta(1))*cos(theta(3))+sin(theta(1))*sin(theta(3))*sin(theta(2)),    cos(theta(2))*sin(theta(1));
     sin(theta(1))*sin(theta(3))+cos(theta(1))*cos(theta(3))*sin(theta(2)),    cos(theta(1))*sin(theta(3))*sin(theta(2))-cos(theta(3))*sin(theta(1)),    cos(theta(1))*cos(theta(2))];

xdot = [x(8); x(10); x(12)];

gb = R * gravity;
% T = (Fb+[0; 0; m*g])/(cos(x(1))*cos(x(3)));
Fd = -kdr * xdot;

ab = gb + T/m + kdr;
end

function omegadotb = angular_acceleration(tau, omegab, I)
omegadotb = I\(tau - cross(omegab, I * omegab));

end

function thetadotv = omegab2thetadotv(omegab, x)
theta = [x(1); x(3); x(5)];
r = [1 sin(theta(1))*tan(theta(2)) cos(theta(1))*tan(theta(2));
    0 cos(theta(1)) -sin(theta(1));
    0 sin(theta(1))*sec(theta(2)) cos(theta(1))*sec(theta(2))];

thetadotv = r * omegab;
end

% function omegab = thetadotv2omegab(thetadotv, x)
% theta = [x(1); x(3); x(5)];
% r = [1 sin(theta(1))*tan(theta(2)) cos(theta(1))*tan(theta(2));
%     0 cos(theta(1)) -sin(theta(1));
%     0 sin(theta(1))*sec(theta(2)) cos(theta(1))*sec(theta(2))];
% omegab = r\thetadotv;
% end

function vi = vb2xdot(vb, x)
theta = [x(1); x(3); x(5)];
R = [cos(theta(3))*cos(theta(2)),                                              cos(theta(2))*sin(theta(3)),                                              -sin(theta(2));
     cos(theta(3))*sin(theta(1))*sin(theta(2))-cos(theta(1))*sin(theta(3)),    cos(theta(1))*cos(theta(3))+sin(theta(1))*sin(theta(3))*sin(theta(2)),    cos(theta(2))*sin(theta(1));
     sin(theta(1))*sin(theta(3))+cos(theta(1))*cos(theta(3))*sin(theta(2)),    cos(theta(1))*sin(theta(3))*sin(theta(2))-cos(theta(3))*sin(theta(1)),    cos(theta(1))*cos(theta(2))];
 
vi = R\vb;
end


function data = smc(m, I, g, x, xd, xddot, lambda, ksmc)
diff = x - xd;

e = [diff(1) diff(3) diff(5) diff(7) diff(9) diff(11)];
edot = [diff(2) diff(4) diff(6) diff(8) diff(10) diff(12)];
s = edot + lambda.*e;
u1 = (m/(cos(x(1))*cos(x(3))))*(xddot(8)+g-lambda(4)*(x(8)-xd(8))-ksmc(4)*sign(s(4)));

ux = (m/u1)*(xddot(10)-lambda(5)*(x(10)-xd(10))-ksmc(5)*sign(s(5)));
uy = (m/u1)*(xddot(12)-lambda(6)*(x(12)-xd(12))-ksmc(6)*sign(s(6)));

xd(10);

p = ux*sin(xd(5))-uy*cos(xd(5));
if p > 1/2
    p = 1/2;
elseif p < -1/2
    p = -1/2;
end

q = (uy*sin(xd(5))+ux*cos(xd(5)))/((1 - (p^2))^(1/2));

if q > 1/2
    q = 1/2;
elseif q < -1/2
    q = -1/2;
end


xd(1) = asin(p);
xd(3) = asin(q);

diff = x - xd;
e = [diff(1) diff(3) diff(5) diff(7) diff(9) diff(11)];
edot = [diff(2) diff(4) diff(6) diff(8) diff(10) diff(12)];
s = edot + lambda.*e;

u2 = I(1,1)*(xddot(2)-((I(2,2)-I(3,3))/I(1,1))*x(4)*x(6)-lambda(1)*(x(2)-xd(2))-ksmc(1)*sign(s(1)));
u3 = I(2,2)*(xddot(4)-((I(3,3)-I(1,1))/I(2,2))*x(2)*x(6)-lambda(2)*(x(4)-xd(4))-ksmc(2)*sign(s(2)));
u4 = I(3,3)*(xddot(6)-((I(1,1)-I(2,2))/I(3,3))*x(2)*x(4)-lambda(3)*(x(6)-xd(6))-ksmc(3)*sign(s(3)));
u = [u1;
     u2;
     u3;
     u4];

fx = [x(2); ((I(2,2)-I(3,3))/I(1,1))*x(4)*x(6);
      x(4); ((I(3,3)-I(1,1))/I(2,2))*x(2)*x(6);
      x(6); ((I(1,1)-I(2,2))/I(3,3))*x(2)*x(4);
      x(8); -g;
      x(10); 0;
      x(12); 0];
  
 gx = [0,                   0,           0,          0;
       0,                   1/I(1,1),    0,          0;
       0,                   0,           0,          0;
       0,                   0,           1/I(2,2),   0;
       0,                   0,           0,          0;
       0,                   0,           0,          1/I(3,3);
       0,                   0,           0,          0;
       cos(x(1))*cos(x(3))/m,   0,           0,          0;
       0,                   0,           0,          0;
       ux/m,                0,           0,          0;
       0,                   0,           0,          0;
       uy/m,                0,           0,          0];

xdot = fx + gx*u;

F = xdot(8)*m;
tau = I*[xdot(2); xdot(4); xdot(6)];

data = [F; tau; xd(1); xd(3)];
end
