% Task 1 code for ENME403 Paramter ID Assignment
% Least Squares
% 23/05/2020
% Tim Hadler

clc, clear

% Load project data
load('StudentID_number44663394.mat')

Ag = Project_Data(:, 1);
V = Project_Data(:, 2);
Vdot = Project_Data(:, 3);
Vddot = Project_Data(:, 4);
Z = Project_Data(:, 5);

% Breakpoints
bp = zeros(1, 13);       
bp(1) =1;               
bp(13) = 2001;

% Data
xdata = linspace(1, 2001, 2001).';
ydata = Ag;

% Segment function
fun = @(th, xdata) th(1) + th(2) * xdata;

% Possible range for each bp, length(data) // #segments
rangebp = 154;

for i=2:12
    bp(i) = (i-1)*rangebp;
end

Ri=0;
R=0;
theta=zeros(12,2);
for n=1:12
    theta(n, :) = lsqcurvefit(fun, [1, 1], xdata(bp(n):bp(n+1)), ydata(bp(n):bp(n+1)));
    plc = 0;
    if n~= 1
        plc = (theta(n-1,1)-theta(n,1))/(theta(n,2)-theta(n-1,2));
        if plc > bp(n) && plc < bp(n+1)
            bp(n) = plc;
        else
            pp=666
            n
            break
        end
    end
   
    
    Ri=0;
    for k = xdata(round(bp(n))):(xdata(round(bp(n+1)))-1)
        Ri = Ri + (ydata(k) - theta(n,1)-theta(n,2) * xdata(k))^2;
    end
    
    plot(xdata(bp(n):bp(n+1)), theta(n,1) + theta(n,2)*xdata(bp(n):bp(n+1)));
    hold on
end



plot(xdata, ydata);
 