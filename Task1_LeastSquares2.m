 % Task 1 code for ENME403 Paramter ID Assignment
% Least Squares
% 23/05/2020
% Tim Hadler

clc, clear, close all

% Load project data
load('StudentID_number44663394.mat')

Ag = Project_Data(:, 1);        % External acceleration
V = Project_Data(:, 2);
Vd = Project_Data(:, 3);
Vdd = Project_Data(:, 4);
Z = Project_Data(:, 5);

m = Mass;
c = Damping;


% Equations
f = @(ag, vdd, vd) -m*ag - m*vdd-c*vd;         % EOM

% Matricies
A = [V Z];
F = zeros(2001, 1);
myAg = zeros(2001, 1);

for i = 1:2001
    F(i) = f(Ag(i), Vdd(i), Vd(i));
end

% LS solve for theta
% th1 = alpha*k0, th2 = (1-alpha)*k0
th = inv(A'*A)*A'*F;

k0 = sum(th);
alpha = th(1) / k0;

for i = 1:2001
    myAg(i) = -Vdd(i) - c/m*Vd(i) - alpha*k0/m*V(i) - (1-alpha)*k0/m*Z(i);
end