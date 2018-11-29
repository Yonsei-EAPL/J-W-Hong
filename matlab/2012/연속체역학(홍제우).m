%% history
% 2013-05 Je-Woo Hong, Term-Project

%% statement
% This program is made for the thermal transport problem,
%   to solve the diffusion equation and advection equation.
% Mother equation is 
%   " dT/dt = -u dT/dx + k d2T/dx2 ".
% FTCS (forward-Time Central-space) Method and Leapfrog Method are used.
% FTCS method is using,
%   " (T_i^n+1 - T_i^n)/dt = k (T_i+1^n - 2T_i^n + T_i-1^n)/(dx^2)  "
% Leapfrog method is using, 
%   " (T_i^n+1 - T_i^n-1)/2dt = k (T_i+1^n - 2T_i^n + T_i-1^n+1)/(dx^2) "


%% assinging constant values
L = 120;    % length of advection
alpha = 3.5*10^(-7);    % advection constant
T0 = 274.15;    % initial temperature
dx = 1;    % lattice of length
t = 100;    % total time
dt = 0.5;   % time interval
v= 1;   % velocity of temperature
w = 8;  % width of temperature
rga = 1;
rgb = L/dx;    % number of lattice
C = v*dt/(2*dx);    % constant of advection calculation
D = alpha*dt/dx^2;    % constant of diffusion calculation
m = t/dt;    % the number of times of calculation
temp = zeros(rgb,2);    % temporary bin for FTCS
temp2 = zeros(rgb,3);   % temporary bin for Leapfrog
result_FTCS = zeros(m,rgb,2);   % final result for FTCS(time,distance,temperature)
result_Leapfrog = zeros(m,rgb,2);   % final result for Leapfrog


%% initial condition
for i = rga:rgb
    for j = 1:m
        result_FTCS(j,i,1) = i*dx;
        result_Leapfrog(j,i,1) = i*dx;
    end
end
for i = rga+1:rga+1+w
    temp(i,1) = T0;
    temp2(1,1) = T0;
end



%% calculation FTCS method
for i = 1:m
    for j = rga+1:rgb-1
        temp(j,2) = D*temp(j+1,1)+(1-C-2*D)*temp(j,1)+(C+D)*temp(j-1,1);
    end
    % updating
    for j = rga:rgb
        temp(j,1) = temp(j,2);
        result_FTCS(i,j,2) = temp(j,2);
    end
    if mod(i,10)==0
        figure(1)
        plot(temp(:,2))
        hold on
    end
end


%% calculation Leapfrog method
% first step using upwind scheme
for i = rga+1:rgb-1
    temp2(i,2) = (1-C)*temp2(i,1)+C*temp2(i-1,1);
end
for i = 1:m
    for j = rga+1:rgb-1
        temp2(j,3) = temp2(j,1) - C*(temp2(j+1,2)+temp2(j-1,2)) + 2*D*(temp2(j+1,2)-2*temp2(j,2)+temp2(j-1,2));
    end
    % updating
    for j = rga:rgb
        temp2(j,1) = temp2(j,2);
        temp2(j,2) = temp2(j,3);
        result_Leapfrog(i,j,2) = temp2(j,3);
    end
    if mod(i,10)==0
        figure(2)
        plot(temp2(:,3))
        hold on
    end
end

%% clear other variable
clear L alpha T0 dx t dt v w rga rgb C D m temp temp2
clear i j 
