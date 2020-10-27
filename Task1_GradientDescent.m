% Task 1 code for ENME403 Paramter ID Assignment
% Gradient Descent
% 17/05/2020
% Tim Hadler

clc, clear;

% Load project data
load('StudentID_number44663394.mat')

Ag = Project_Data(:, 1);
V = Project_Data(:, 2);
Vdot = Project_Data(:, 3);
Vddot = Project_Data(:, 4);
Z = Project_Data(:, 5);

% Symbolic variables, v derivatives
% input acceleration ag, hysteretic element z
% unknown parameters theta1 = k0, initial sprink coeff
% theta2 = alpha, post-yielding ratio
syms v vdot vddot ag z theta1 theta2 

% Known parameters
m = Mass;
c = Damping;

% Algorithm parameters
n = length(Project_Data);   % number of data points
lr = 0.01;                     % Learning rate
tol = 1e-5;                    % Error tolerance
theta0 = [1;1];             % Initial guess vector


theta = [];           % Storage for theta vectors
grad_J = [];          % Storage for cost function gradient per iteration


%EOM: ag = -vddot - c/m*vdot - alpha*k0/m*v - (1-alpha)*k0/m*z + e
%EOM rearranged for measurment error e. 
e = ag + vddot + c/m*vdot + theta2*theta1/m*v + (1-theta2)*theta1/m*z;

% Cost function J
% J = sum of e^2
J = 0;
j=e^2;
for i = 1:n
    J = J + subs(j, {ag, v, vdot, vddot, z}, {Ag(i), V(i), Vdot(i), Vddot(i), Z(i)});
end

% Gradient of cost function J
% Partial derivatives of J, w respect to unknown parameters

pd_J_theta1 = diff(J, theta1);
pd_J_theta2 = diff(J, theta2);
            
theta(:, 1) = theta0 - lr * [subs(pd_J_theta1, {theta1, theta2}, {theta0(1), theta0(2)});
                             subs(pd_J_theta2, {theta1, theta2}, {theta0(1), theta0(2)})];
for k = 1:500
    grad_J(:, k) = [subs(pd_J_theta1, {theta1, theta2}, {theta(1, k), theta(2, k)});
                    subs(pd_J_theta2, {theta1, theta2}, {theta(1, k), theta(2, k)})
                    ];
                
    
    
    theta(:, k+1) = theta(:, k) - lr * grad_J(:, k);
    
    if (k~=1)
        if ((abs(theta(k) - theta(k-1))) < tol) 
            break
        end
    end
end
k
final_theta = theta(:, k+1)
