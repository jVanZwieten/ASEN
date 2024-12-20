classdef HypersonicsUtilities
    properties(Constant)
        airDensity_seaLevel = 1.225; % kg/m^3, rho_SL
        airPressure_seaLevel = 1.01325e5; % N/m^2
        inverseScaleHeight = 1.387e-4; % /m, aTilde
        earthSurfaceGravityAcceleration = 9.81; % m/s^2, g
        earthRadius = 6378e3; % m
        gasConstant_air = 287; % J/kgK
        gasConstant_universal = 8314.5; % J/kgK
        ratioSpecificHeats_air = 1.4;
        specificHeatAtConstantPressure_air = HypersonicsUtilities.specificHeatAtConstantPressure(HypersonicsUtilities.ratioSpecificHeats_air, HypersonicsUtilities.gasConstant_air);
        velocityRatio_decelMax = exp(-.5); % v/v_0
        velocityRatio_heatingMax = exp(-1/6); % v/v_0

        pascalsPerAtmosphere = 101325;

        stefanBoltzmannConstant = 5.67e-8; % W/m^2K^4
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
        
        function qDot_max = maxHeatingBallistic(k, beta, velocity_0)
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

        function a = speedOfSoundAir(T)
            gamma = HypersonicsUtilities.ratioSpecificHeats_air;
            R = HypersonicsUtilities.gasConstant_air;
            a = HypersonicsUtilities.speedOfSound(gamma, R, T);
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

        function alphaOverAlpha1 = massFractionFunction(rho, rho_d, theta_d, T)
            alphaOverAlpha1 = rho_d/rho*exp(-theta_d/T);
        end

        function R_s = speciesGasConstant(molarMass_s)
            R = HypersonicsUtilities.gasConstant_universal;
            R_s = R/molarMass_s;
        end

        function c_vVib = specificHeatVibrational(R, theta_v, T)
            a = theta_v/T;
            c_vVib = R*a^2*exp(a)/(exp(a) - 1)^2;
        end

        function c_p = specificHeatConstantPressureDiatomic(R, theta_v, T)
            c_p = 7/2*R + HypersonicsUtilities.specificHeatVibrational(R, theta_v, T);
        end

        function c_p = specificHeatConstantPressureMonoatomic(R)
            c_p = 5/2*R;
        end

        function p = idealGasPressure(rho, T)
            R = HypersonicsUtilities.gasConstant_air;
            p = rho*R*T;
        end

        function rho = idealGasDensity(p, T)
            R = HypersonicsUtilities.gasConstant_air;
            rho = p/R/T;
        end

        function T_0OverT = stagnationTemperatureRatio(M, gamma)
            T_0OverT = 1 + (gamma - 1)/2*M^2;
        end

        function T_0OverT = stagnationTemperatureRatioAir(M)
            gamma = HypersonicsUtilities.ratioSpecificHeats_air;
            T_0OverT = HypersonicsUtilities.stagnationTemperatureRatio(M, gamma);
        end

        %% liftingReentry
        function qDot_max = maxHeatingLifting(k, B, LDRatio, v_0)
            qDot_max = k*sqrt(8*B/27/LDRatio)*v_0^2;
        end
        
        function v = velocityAtMaxHeatingLifting(v_0)
            v = v_0*sqrt(2/3);
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

        %% boundary layer
        function qDot_w = detraEtAllHeatFlux(rho_inf, u_inf, R_N)
            qDot_w = 5.156e-5*sqrt(rho_inf)*u_inf^3.15/sqrt(R_N);
        end

        function dudx_e = velocityGradient(R_N, p_e, p_inf, rho_e)
            dudx_e = 1/R_N*sqrt(2*(p_e - p_inf)/rho_e);
        end

        function qDot_w = fayRiddellHeatFlux(k, Pr_w, rho_e, mu_e, rho_w, mu_w, dudx_e, h_0e, h_w, F)
            qDot_w = k/Pr_w^.6*(rho_e*mu_e)^.4*(rho_w*mu_w)^.1*sqrt(dudx_e)*(h_0e - h_w)*F;
        end
        
        function F = fayRiddelCase1(Le, h_dOverH_0e)
            F = 1 + (Le^0.52 - 1)*h_dOverH_0e;
        end
                
        function F = fayRiddelCase2(Le, h_dOverH_0e)
            F = 1 + (Le^0.63 - 1)*h_dOverH_0e;
        end
                
        function F = fayRiddelCase3(h_dOverH_0e)
            F = 1 - h_dOverH_0e;
        end

        function TStar = EckertReferenceTemperature(T_e, T_0e, T_w, r)
            TStar = .5*(T_e + T_w) + .22*r*(T_0e - T_e);
        end
        
        function h_aw = hawFromR(r, h_0e, h_e)
            h_aw = r*(h_0e - h_e) + h_e;
        end

        function C_h = reynoldsAnalogy(C_fx, Pr)
            C_h = .5*C_fx*Pr^(-2/3);
        end

        function qDot_w = heatFluxFromHeatingCoefficient(C_h, rhoStar, u_e, h_aw, h_w)
            qDot_w = C_h*rhoStar*u_e*(h_aw - h_w);
        end

        function qDot_rad = radiationHeatFlux(T_w, emissivity)
            sigma = HypersonicsUtilities.stefanBoltzmannConstant;
            qDot_rad = emissivity*sigma*(T_w^4);
        end
    end
end