classdef HypersonicsUtilities
    properties(Constant)
        airDensity_seaLevel = 1.225; % kg/m^3, rho_SL
        airPressure_seaLevel = 1.01325e5; % N/m^2
        inverseScaleHeight = 1.387e-4; % /m, aTilde
        earthSurfaceGravityAcceleration = 9.81; % m/s^2, g
        earthRadius = 6378e3; % m
        gasConstant_air = 287; % J/kgK
        ratioSpecificHeats_air = 1.4;
        velocityRatio_decelMax = exp(-.5); % v/v_0
        velocityRatio_heatingMax = exp(-1/6); % v/v_0

        pascalsPerAtmosphere = 101325;
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
        function KEh0 = kineticEnergyRatio(M_1, gamma)
            KEh0 = M_1.^2./(2/(gamma - 1) + M_1.^2);
        end

        function Re = reynoldsNumber(rho, u, L, mu)
            Re = rho*u*L/mu;
        end

        function mu = sutherlandsLawViscosity(T)
            mu = 1.458e-6*T^1.5/(T + 110.4);
        end

        function mu = sutherlandsLawThermalConductivity(T)
            mu = 1.993e-3*T^1.5/(T + 112);
        end

        function a = speedOfSound(gamma, R, T)
            a = sqrt(gamma*R*T);
        end

        function c_p = specificHeatAtConstantPressure(gamma, R)
            c_p = gamma/(gamma - 1)*R;
        end

        function h_0 = enthalpyTotal(h, u)
            h_0 = h + .5*u^2;
        end

        function rho = densityRealGasEffect(p, z, R, T)
            rho = p/(z*R*T);
        end

        %% normal shock relations
        function ratio = densityVelocity12Ratio(gamma, M_1)
            ratio = (gamma + 1)*M_1^2/((gamma - 1)*M_1^2 + 2);
        end

        function ratio = pressureRatio(gamma, M_1)
            ratio = (2*gamma*M_1^2 - (gamma - 1))/(gamma + 1);
        end

        function M_2 = machAfterShock(gamma, M_1)
            M_2 = sqrt((2 + (gamma - 1)*M_1^2)/(2*gamma*M_1^2 - (gamma - 1)));
        end
    end
end