classdef HypersonicsUtilities
    properties(Constant)
        inverseScaleHeight = 1.1387e-4 % /m, aTilde
        airDensity_seaLevel = 1.225 % kg/m^3, rho_SL
        earthSurfaceGravityAcceleration = 9.81 % m/s^2, g
        velocityRatio_decelMax = exp(-.5); % v/v_0
        velocityRatio_heatingMax = exp(-1/6); % v/v_0
    end

    methods(Static)
        %% trajectories
        function rho = airDensity(altitude)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            rho_SL = HypersonicsUtilities.airDensity_seaLevel;
            rho = rho_SL*exp(-aTilde*altitude);
        end

        function h = altitude(airDensity)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            rho_SL = HypersonicsUtilities.airDensity_seaLevel;
            h = log(airDensity/rho_SL)/(-aTilde);
        end

        function d = drag(dragCoefficient, area, airDensity, velocity)
            d = .5*dragCoefficient*area*airDensity.velocity;
        end

        function l = lift(liftCoefficient, area, airDensity, velocity)
            l = .5*liftCoefficient*area*airDensity*velocity;
        end

        function b = ballisticCoefficient(weight, dragCoefficient, area)
            b = weight/dragCoefficient/area;
        end

        function lb = ballisticCoefficientLift(weight, liftCoefficient, area)
            lb = weight/liftCoefficient/area;
        end

        function qDot = heatTransfer(k, airDensity, velocity)
            qDot = k*sqrt(airDensity)*velocity^3;
        end

        %% ballistic trajectory
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

        function gs = maxG(velocity_0, flightPathAngle)
            aTilde = HypersonicsUtilities.inverseScaleHeight;
            g = HypersonicsUtilities.earthSurfaceGravityAcceleration;
            gs = -velocity_0^2*aTilde*sin(flightPathAngle)/2/g*exp(-1);
        end

        function qDot_max = maxHeating(k, beta, velocity_0)
            qDot_max = k*sqrt(beta/6)*exp(-.5)*velocity_0^3;
        end
    end
end