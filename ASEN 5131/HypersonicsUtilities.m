classdef HypersonicsUtilities
    properties(Constant)
        inverseScaleHeight = 1.387e-4; % /m, aTilde
        airDensity_seaLevel = 1.225; % kg/m^3, rho_SL
        earthSurfaceGravityAcceleration = 9.81; % m/s^2, g
        earthRadius = 6378e3; % m
        velocityRatio_decelMax = exp(-.5); % v/v_0
        velocityRatio_heatingMax = exp(-1/6); % v/v_0
    end
    
    methods(Static)
        %% trajectories
        function gammaDot = flightPathAngleRate(velocity, lift, weight, radialDistance, flightPathAngle)
            g = HypersonicsUtilities.earthSurfaceGravityAcceleration;
            gammaDot = (g/velocity)*(lift/weight - (1 - velocity^2/g/radialDistance)*cos(flightPathAngle));
        end
        
        function vDot = acceleration(thrust, drag, m, flightPathAngle)
            g = HypersonicsUtilities.earthSurfaceGravityAcceleration;
            vDot = (thrust/m - drag/m - g*sin(flightPathAngle));
        end
        
        function vDot = accelerationW(thrust, drag, weight, flightPathAngle)
            m = weight/HypersonicsUtilities.earthSurfaceGravityAcceleration;
            vDot = HypersonicsUtilities.acceleration(thrust, drag, m, flightPathAngle);
        end
        
        function rDot = radialDistanceRate(velocity, flightPathAngle)
            rDot = velocity*sin(flightPathAngle);
        end
        
        function thetaDot = earthAngleRate(velocity, flightPathAngle, radialDistance)
            thetaDot = velocity*cos(flightPathAngle)/radialDistance;
        end
        
        function h = altitude(radialDistance)
            h = radialDistance - HypersonicsUtilities.earthRadius;
        end
        
        function rho = airDensity(altitude)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            rho_SL = HypersonicsUtilities.airDensity_seaLevel;
            rho = rho_SL*exp(-aTilde*altitude);
        end
        
        function h = altitudeFromAirDensity(airDensity)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            rho_SL = HypersonicsUtilities.airDensity_seaLevel;
            h = log(airDensity/rho_SL)/(-aTilde);
        end
        
        function d = drag(dragCoefficient, area, airDensity, velocity)
            d = .5.*dragCoefficient.*area.*airDensity.*velocity.^2;
        end
        
        function l = lift(liftCoefficient, area, airDensity, velocity)
            l = .5.*liftCoefficient.*area.*airDensity.*velocity.^2;
        end
        
        function b = ballisticCoefficient(weight, dragCoefficient, area)
            b = weight/dragCoefficient/area;
        end
        
        function lb = ballisticCoefficientLift(weight, liftCoefficient, area)
            lb = weight/liftCoefficient/area;
        end
        
        function qDot = heatTransfer(k, airDensity, velocity)
            qDot = k.*sqrt(airDensity).*velocity.^3;
        end
        
        %% ballistic trajectory
        function vDot = decelerationBallisticSimple(m, D)
            vDot = -D/m;
        end
        
        function vr = velocityRatio(airDensity, beta)
            vr = exp(-airDensity/beta);
        end
        
        function beta = betaFromWeight(weight, flightPathAngle, dragCoefficient, area)
            mass = weight/HypersonicsUtilities.earthSurfaceGravityAcceleration;
            beta = HypersonicsUtilities.beta(mass, flightPathAngle, dragCoefficient, area);
        end
        
        function beta = beta(mass, flightPathAngle, dragCoefficient, area)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            beta = -2*mass*aTilde*sin(flightPathAngle)/dragCoefficient/area;
        end
        
        function rho = airDensity_maxHeating(beta)
            rho = beta/6;
        end
        
        function gs = maxG(velocity_0, flightPathAngle)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            g = HypersonicsUtilities.earthSurfaceGravityAcceleration;
            gs = -velocity_0^2*aTilde*sin(flightPathAngle)/2/g*exp(-1);
        end
        
        function qDot_max = maxHeating(k, beta, velocity_0)
            qDot_max = k*sqrt(beta/6)*exp(-.5)*velocity_0^3;
        end

        %% aerothermaldynamics
        function KEh0 = kineticEnergyRatio(M, gamma)
            KEh0 = M.^2./(2/(gamma - 1) + M.^2);
        end
    end
end