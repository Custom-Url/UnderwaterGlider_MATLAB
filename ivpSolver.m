function [t,zRK4] = ivpSolver(t0,z0,dt,tend,deltaV,offset)
% ivpSolver    Solve an initial value problem (IVP) and plot the result
% 
%     [T,Z] = ivpSolver(T0,Z0,DT,TE) computes the IVP solution using a step 
%     size DT, beginning at time T0 and initial state Z0 and ending at time 
%     TEND. The solution is output as a time vector T and a matrix of state 
%     vectors Z.

% Set initial conditions
t(1) = t0;
zRK4(:,1) = z0;

% Continue stepping until the end time is exceeded
n=1;
while t(n) <= tend
    % Increment the time vector by one time step
    t(n+1) = t(n) + dt;

    % Apply RK4s method for one time step
    zRK4(:,n+1) = stepRungeKutta(t(n), zRK4(:,n), dt,deltaV,offset);

    n = n+1;
end
