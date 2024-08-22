function znext = stepRungeKutta(t,z,dt,deltaV,offset)
% stepEuler    Compute one step using the Euler method
% 
%     ZNEXT = stepEuler(T,Z,DT) computes the state vector ZNEXT at the next
%     time step T+DT

% Calculate the state derivative from the current state
dz = stateDeriv(t, z, deltaV,offset);

% Calculate the next state vector from the previous one using RungeKutta
% method
% update equation
A = dt * stateDeriv(t, z, deltaV,offset);
B= dt * stateDeriv(t+dt/2, z+A/2, deltaV,offset);
C= dt * stateDeriv(t+dt/2, z+B/2, deltaV,offset);
D= dt * stateDeriv(t+dt, z+C, deltaV,offset);
znext = z + (A +2*B + 2*C + D)/6;