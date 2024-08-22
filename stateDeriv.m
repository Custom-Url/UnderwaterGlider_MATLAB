function dz = stateDeriv(t,z,deltaV,offset)
%   Calculate the state derivative for the underwater glider system
% 
%     DZ = stateDeriv(T,Z) computes the derivative DZ = [V; A] of the 
%     state vector Z = [X; V], where X is displacement, V is velocity,
%     and A is acceleration.

%   Sets constants for stateDeriv
Gravity = 9.81;
Mass = 3.9;
CDrag = 0.78;
Area = pi/400;
Density = 997;
TimeInterval = 600;

k = ceil(t/TimeInterval);
%   Rounds (t/TimeInterval) up to a constant 'k' which represents an interval of either sinking or climbing

direction = mod(k,2);
%   Calculates if 'k' is odd or even to decide what type of interval it represents
%   Direction = 1 during an odd interval / glider is sinking

if direction == 0
    CLift = -2.76; 
    FS = 1;
else
    CLift = 2.76;
    FS = -1;
end 

%   Buoyancy equation using outputs of shooting.m
Buoyancy = Density * Gravity * (0.0039 +  offset + (deltaV * FS));

%   ODE
dz1 = z(3);
dz2 = z(4);
dz3 = (Density * Area * ((z(3)^2)+(z(4)^2))^0.5) * ((-z(4)*CLift) - (z(3)*CDrag)) /(2*Mass);
dz4 = (Buoyancy/Mass) - Gravity +  (Density * Area * ((z(3)^2)+(z(4)^2))^0.5) * ((z(3)*CLift) + (-z(4)*CDrag)) /(2*Mass);

dz = [dz1; dz2; dz3; dz4];

end