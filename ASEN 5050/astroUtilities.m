classdef astroUtilities
    methods(Static)
        function H = AngularMomentum(R, V)
            H = cross(R, V);
        end

        function h = angularMomentumFromae(a, e, mu)
            h = sqrt(a*(1 - e^2)*mu);
        end
        
        function omega = argumentOfPeriapsisFromNE(N, E)
            omega = astroUtilities.argumentOfPeriapsisFromne(Utilities.UnitVector(N), Utilities.UnitVector(E));
        end
        
        function omega = argumentOfPeriapsisFromne(n, e)
            omega = acos(dot(n, e));
        end

        function r = conicEquation(a, e, nu)
            r = a*(1 - e^2)/(1 + e*cos(nu));
        end
        
        function C = DirectionCosineMatrix(RAAN, AOP, trueAnomaly, inclination)
            theta = AOP + trueAnomaly;
            C = [cos(RAAN)*cos(theta) - sin(RAAN)*cos(inclination)*sin(theta), -cos(RAAN)*sin(theta) - sin(RAAN)*cos(inclination)*cos(theta), sin(RAAN)*sin(inclination);
                sin(RAAN)*cos(theta) + cos(RAAN)*cos(inclination)*sin(theta), -sin(RAAN)*sin(theta) + cos(RAAN)*cos(inclination)*cos(theta), -cos(RAAN)*sin(inclination);
                sin(inclination)*sin(theta), sin(inclination)*cos(theta), cos(inclination)];
        end
        
        function E = EccentricityFromRV(R, V, mu)
            E = (1/mu)*((norm(V)^2 - mu/norm(R))*R - dot(R, V)*V);
        end
        
        function E = EccentricityFromRVH(R, V, H, mu)
            E = cross(V, H)/mu - R/norm(R);
        end
        
        function e = eccentricityFromSemiLatusRectum(p, a)
            e = sqrt(1 - p/a);
        end
        
        function e = eccentricityFromNuInfinity(nu_inf)
            e = -1/cos(nu_inf);
        end
        
        function phi = flightPathAngleFromENu(e, nu)
            phi = atan(e*sin(nu)/(1 + e*cos(nu)));
        end

        function phi = flightPathAngleFromhrv(h, r, v)
            phi = acos(h/r/v);
        end
        
        function i = inclinationFromH(H)
            i = astroUtilities.inclinationFromHhat(Utilities.UnitVector(H));
        end
        
        function i = inclinationFromHhat(Hhat)
            i = acos(Hhat(3));
        end
        
        function N = LineOfNodesFromH(H)
            N = [-H(2); H(1); 0];
        end
        
        function r_p = periapsisFromae(a, e)
            r_p = a*(1 - e);
        end
        
        function R_p = PeriaspisFromaE(a, E)
            r_p = astroUtilities.periapsisFromae(a, norm(E));
            R_p = PeriapsisFromrE(r_p, E);
        end
        
        function R_p = PeriapsisFromrE(r_p, E)
            Ehat = Utilities.UnitVector(E);
            R_p = PeriapsisFromrEhat(r_p, Ehat);
        end

        function R_p = PeriapsisFromrEhat(r_p, Ehat)
            R_p = r_p*Ehat;
        end

        function P = period(a, mu)
            P = 2*pi/sqrt(mu)*sqrt(a^3); % s
        end
        
        function Omega = RAANFromN(N)
            Omega = astroUtilities.RAANFromNhat(N/norm(N));
        end
        
        function Omega = RAANFromNhat(Nhat)
            Omega = acos(Nhat(1));
        end
        
        function p = semiLatusRectumFromH(h, mu)
            p = h^2/mu;
        end
        
        function a = semiMajorAxisFromspecificNRG(nrg_mech, mu)
            a = -mu/(2*nrg_mech);
        end
        
        function specNRG = specificNRGFromA(a,mu)
            specNRG = -mu/(2*a);
        end
        
        function specMechNRG = specificNRGFromRV(R, V, mu)
            specMechNRG = astroUtilities.specificNRGFromrv(norm(R), norm(V), mu);
        end
        
        function specMechNRG = specificNRGFromrv(r, v, mu)
            specMechNRG = v^2/2 - mu/r;
        end
        
        function nu = trueAnomalyFromRE(R, E)
            nu = acos(dot(E, R)/(norm(E)*norm(R)));
        end

        function nu = trueAnomalyFromhre(h, r, e, mu)
            nu = acos((h^2/r/mu - 1)/e);
        end
        
        function del = turningAngleFromEccentricity(e)
            del = 2*asin(1/e);
        end
        
        function v = velocityAtApsisFromHR(h, r)
            v = h/r;
        end
        
        function v = velocityAtRFromSpecificNRG(specNRG, r, mu)
            v = sqrt(2*(specNRG + mu/r));
        end
    end
end