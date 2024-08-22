function [offset, deltaV] = Shooting(t0,z0,dt,tend,x,y,floor)
%   Shooting        Computes offset and detlaV for a target point
%
%       [OFFSET,DELTAV] = Shooting(T0,Z0,DT,TEND,X,Y,FLOOR) computes the BVP solution
%       using a step size DT, beginning at time T0 and initial state Z0 and ending at time
%       TEND. The target point is X,Y and the solution avoids going below FLOOR.

%   Set initial conditons for calculating offset
offset = 0;
off1 = -6.0000e-05;
off2 = 6.0000e-05;
Error1 = 1;
Error2 = 1;

%   Calculates a displacement offset so glider has neutral bouyancy (accurate to within 1mm)
while  abs(Error1) && abs(Error2) > 0.001
    %   Runs ivpSolver for an initial offset guess and calculates Error1
    [t,zRK4] = ivpSolver(t0,z0,dt,tend,0,off1);
    Error1 = zRK4(2,100) - z0(2);
    %   Runs ivpSolver for a secondary offset guess and calculates Error2
    [t,zRK4] = ivpSolver(t0,z0,dt,tend,0,off2);
    Error2 = zRK4(2,100) - z0(2);
    %   Calculates next guess, offset, using previous offsets and errors
    offset = off2 - (Error2*(off2 - off1)/(Error2-Error1));
    %   Replaces the furthest previous guess with the new offset
    if abs(Error1) < abs(Error2)
        off2 = offset;
    else
        off1 = offset;
    end
end

%   Checks conditions are suitable for the glider to be neutrally bouyant
if abs(offset) > 6e-05
    Result = 'UNSUITABLE BUOYANCY CONDITIONS:';
    fprintf('%s Required offset displacement is %d ml.\n', Result,round(offset * 1e06));
else
    Result = 'SUITABLE BUOYANCY CONDITIONS:';
    fprintf('%s Required offset displacement is %d ml.\n', Result,round(offset * 1e06));

    %   Create new variables and arrays for second shooting method to calculate deltaV
    %   Arrays hold values for each iteration of initial guesses
    Iteration = 1;
    deltaV_array = [0,0,0,0,0,0,];
    Accuracy_array = [0,0,0,0,0,0,];
    Time_array = [0,0,0,0,0,0,];

    %   Iteration variable allows 6 initial guesses (0ml-10ml, 10ml-20ml ect)
    %   to ensure all possible solutions are found
    while Iteration < 7
        %   Set initial conditons for calculating displacement
        dv1 = (Iteration - 1) * 1.0000e-05;
        dv2 = Iteration * 1.0000e-05;
        Error1 = 1;
        Error2 = 1;
        Counter = 0;

        %   Calculates displacement required to go through TARGET (aims for accuracy within 5mm)
        while  abs(Error1) && abs(Error2) > 0.005
            %   Runs ivpSolver for an initial displacement guess and calculates Error1
            [t,zRK4] = ivpSolver(t0,z0,dt,tend,dv1,offset);
            Xvalues = zRK4(1,:);
            Xdiff = abs(x-Xvalues);
            [m,i] = min(Xdiff);
            Error1 = zRK4(2,i) - y;
            %   Calculates value for Accuracy - Distance between target point and closest glider position
            Accuracy = ((Error1^2) + (m^2))^0.5;
            %   Runs ivpSolver for a secondary displacement guess and calculates Error2
            [t,zRK4] = ivpSolver(t0,z0,dt,tend,dv2,offset);
            Xvalues = zRK4(1,:);
            Xdiff = abs(x-Xvalues);
            [m,i] = min(Xdiff) ;
            Error2 = zRK4(2,i) - y;
            %   Calculates next guess, displacement, using previous displacements and errors
            deltaV = dv2 - (Error2*(dv2 - dv1)/(Error2-Error1));

            %   Replaces the furthest previous guess with the new offset
            if abs(Error1) < abs(Error2)
                dv2 = deltaV;
            else
                dv1 = deltaV;
                Accuracy = ((Error1^2) + (m^2))^0.5;
            end
            %   Counter ends the while loop, if Error has not been suitably reduced, after 20 iterations
            Counter = Counter + 1;
            if Counter == 20
                %   Sets Error below threshold to force end the while loop
                Error1 = 0.001;
                Error2 = 0.001;
            end
        end
        %   deltaV < 0 is impossible results therefore deltaV = Inf to stop value being considered later
        if deltaV < 0
            deltaV = Inf;
        end
        %   Adds values to previously created arrays for later selection
        deltaV_array(Iteration) = deltaV;
        Accuracy_array(Iteration) = Accuracy;
        Time_array(Iteration) = ceil(i/(600/dt));
        Iteration = Iteration + 1;
    end

    %   Selects most battery efficient deltaV
    Efficiency = abs(deltaV_array) .* Time_array;
    [a,b] = min(Efficiency);
    deltaV = deltaV_array(b);
    [t,zRK4] = ivpSolver(t0,z0,dt,tend,deltaV,offset);

    %   Checks to see if route exceeds glider's displacement capacity and attempts to find a superior route
    if abs(offset) + deltaV > 6e-05
        [deltaV,b] = min(deltaV_array(deltaV_array>0));
        [t,zRK4] = ivpSolver(t0,z0,dt,tend,deltaV,offset);
        if abs(offset) + deltaV < 6e-05
            disp 'ALTERNATIVE DELTAV WITH LESS THAN MAXIMUM EFFICIECNCY FOUND TO REDUCE DISPLACEMENT'
        end
    end

    %   Checks to see if route collides with floor and attempts to find a superior route
    if min(zRK4(2,:)) < floor                 %Example Case: Shooting(0,[0;0;0;0],0.1,6000,400,-20,-23.5);
        [deltaV,b] = min(deltaV_array(deltaV_array>0));
        [t,zRK4] = ivpSolver(t0,z0,dt,tend,deltaV,offset);
        if min(zRK4(2,:)) > floor
            disp 'ALTERNATIVE ROUTE WITH LESS THAN MAXIMUM EFFICIENCY FOUND TO AVOID FLOOR COLLISION'
        end
    end

    %   Calls Accuracy again to ensure value matches current deltaV selection
    Accuracy = Accuracy_array(b);


    %   Below section relates to output of a figure and text output in the command window to improve usability
    %   Text ouptut in the command window relating to success/failure point of the simulation
    if min(zRK4(2,:)) < floor
        disp 'ROUTE CANNOT BE FOUND: PATHING COLLIDES WITH FLOOR'
    else
        if abs(offset) + deltaV > 6e-05
            Result = 'EXCEEDS DISPLACEMENT CAPACITY:';
        else
            Result = 'MOST ENERGY EFFICIENT ROUTE FOUND:';
        end

        fprintf('%s Required total displacement is %d ml.\n',Result,round((deltaV + abs(offset)) * 1e06));

        if abs(offset) + deltaV < 6e-05
            fprintf('Accuracy is %g m.\n\n',round(Accuracy, 3, "significant"));
            %   Checks route Accuracy and attempts to find a superior route if Accuracy > 1m
            if Accuracy > 1
                disp 'Target May Be Outside of Nominal Pathing or Other Routes Being Less Efficient'

                %   Produces reduced Accuracy_array with outliers set to Inf to exclude consideration as improved route
                Average_Accuracy = filloutliers(Accuracy_array,Inf, 'percentiles', [0,70]);
                %   New Efficiency metric where Accuracy is now considered and new most efficient route is selected
                Efficiency = abs(deltaV_array) .* Time_array .* Average_Accuracy;
                [a,b] = min(Efficiency);
                deltaV = deltaV_array(b);
                Accuracy = Accuracy_array(b);
                %   Checks if new deltaV has Accuracy < 1m and outputs the suitability of this route
                if  Accuracy < 1
                    if abs(offset) + deltaV < 6e-05
                        fprintf('%d ml Displacement May Be Less Energy Efficient But Has a Accuracy of %d mm.\n',round((deltaV + offset) * 1e06),round(Accuracy*1e3));
                    else
                        disp 'More Accurate Route Exceeds Displacement Capacity '
                    end
                end
            else

                %     Plot the result
                figure(1)
                patch([-1 1 1 -1 -1]*1e12,[floor-1e12 floor-1e12 0 0 -1],[0 0 0 0 0],'y');  %   Adds 'water' to plot
                patch([-1 1 1 -1 -1]*1e12,[0 0 floor floor -1],[0 0 0 0 0],'c');    %   Adds 'floor'to plot

                hold on
                plot(zRK4(1,:),zRK4(2,:),'r','LineWidth',1.5)   %   Plot figure1

                text(x,y,'\leftarrow TARGET')           %   Labels TARGET on figure1
                xlim([0 max(zRK4(1,:))])
                ylim([(floor-10) (max(zRK4(2,:))+10)])  %   Sets autoscale parameters
                xlabel('Horizontal Displacement')
                ylabel('Vertical Displacement')
                legend('Floor','Water','Glider Path')
                hold off
            end
        end
    end
end