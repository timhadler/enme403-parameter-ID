% Task 2 code for ENME403 Paramter ID Assignment
% Integral method
% 27/05/2020
% Tim Hadler

clc, clear

% Load project data
load('StudentID_number44663394.mat')

Ag = Project_Data(:, 1);        % External acceleration
V = Project_Data(:, 2);
Vd = Project_Data(:, 3);
Vdd = Project_Data(:, 4);
Z = Project_Data(:, 5);

m = Mass;
c = Damping;
xdata = 1:2001;

dt = 0.005;

syms x th
sign = piecewise(x < 0, -1, x > 0, 1, x == 0, 0);

P = zeros(2001, 1);

for i = 1:2001
    P(i) = 0.5*Vd(i)*(subs(sign, x, Vd(i)*Z(i))+1)*(abs(Z(i)))^2;
end

a1 = 1;
b1 = 2001;
a2 = 1;
b2 = 1000;
a3 = 1001;
b3 = 2001;

Z1 = Z(b1) - Z(a1);
Z2 = Z(b2) - Z(a2);
Z3 = Z(b3) - Z(a3);

x1 = -dt/2*(P(a1) + P(b1) + 2*((sum(P(a1+1:b1-1)))));
x2 = -dt/2*(P(a2) + P(b2) + 2*((sum(P(a2+1:b2-1)))));
x3 = -dt/2*(P(a3) + P(b3) + 2*((sum(P(a3+1:b3-1)))));

y1 = Z1 - dt/2*((Vd(a1) + Vd(b1)) + 2*sum(V(a1+1:b1-1)));
y2 = Z2 - dt/2*((Vd(a2) + Vd(b2)) + 2*sum(V(a2+1:b2-1)));
y3 = Z3 - dt/2*((Vd(a3) + Vd(b3)) + 2*sum(V(a3+1:b3-1)));

Y = [y1;y2;y3];
X = [x1;x2;x2];

th = (X'*X)\X'*Y;



% Checking
ft = @(ag, vdd, vd) -m*ag - m*vdd-c*vd; 


a = 0.2637;
k0 = 1.1067;

dy = 1/sqrt(abs(th));
dzdv = @(i) 1 - 0.5*(subs(sign, x, Vd(i)*Z(i)) + 1)*(abs(Z(i)/dy))^2;
f = @(dzdv) a*k0 + (1 - a)*k0*dzdv;

i = 2000;
myf = f(dzdv(i))
F = ft(Ag(i), Vdd(i), Vd(i))


