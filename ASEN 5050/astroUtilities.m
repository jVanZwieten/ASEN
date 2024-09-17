classdef astroUtilities
    methods(Static)
        function specMechNRG = specificMechanicalNRGFromA(a,mu)
            specMechNRG = -mu/(2*a);
        end

        function specMechNRG = specificMechanicalNRGFromRV(R, V, mu)
            specMechNRG = astroUtilities.specificMechanicalNRGFromrv(norm(R), norm(V), mu);
        end

        function specMechNRG = specificMechanicalNRGFromrv(r, v, mu)
            specMechNRG = v^2/2 - mu/r;
        end

        function v = vFromSpeceificMechanicalNRG(nrg_mech, mu, r)
            v = sqrt(2*(nrg_mech + mu/r));
        end

        function p = semiLatusRectumFromH(h, mu)
            p = h^2/mu;
        end

        function e = eccentricityFromSemiLatusRectum(p, a)
            e = sqrt(1 - p/a);
        end

        function E = EccentricityFromRVH(R, V, H, mu)
            E = cross(V, H)/mu - R/norm(R);
        end

        function E = EccentricityFromRV(R, V, mu)
            E = (1/mu)*((norm(V)^2 - mu/norm(R))*R - dot(R, V)*V);
        end

        function a = semiMajorAxisFromSpecificMechanicalNRG(nrg_mech, mu)
            a = -mu/(2*nrg_mech);
        end

        function e = eccentricityFromTrueAnomalyInfinity(theta)
            e = -1/cos(theta);
        end

        function del = turningAngleFromEccentricity(e)
            del = 2*asin(1/e);
        end

        function r_p = periapsisFromAE(a, e)
            r_p = a*(1 - e);
        end

        function v = velocityAtRFromMechNRG(nrg_mech, r, mu)
            v = sqrt(2*(nrg_mech + mu/r));
        end

        function H = angularMomentumFromRV(R, V)
            H = cross(R, V);
        end

        function i = inclinationFromH(H)
            i = astroUtilities.inclinationFromHhat(H/norm(H));
        end

        function i = inclinationFromHhat(Hhat)
            i = acos(Hhat(3));
        end

        function N = LineOfNodesFromH(H)
            N = [-H(2); H(1); 0];
        end

        function Omega = RAANFromN(N)
            Omega = astroUtilities.RAANFromNhat(N/norm(N));
        end

        function Omega = RAANFromNhat(Nhat)
            Omega = acos(Nhat(1));
        end

        function omega = argumentOfPeriapsisFromNE(N, E)
            omega = acos(dot(N, E)/(norm(N)*norm(E)));
        end

        function nu = trueAnomalyFromRE(R, E)
            nu = acos(dot(E, R)/(norm(E)*norm(R)));
        end

        function C = DirectionCosineMatrix(RAAN, AOP, trueAnomaly, inclination)
            theta = AOP + trueAnomaly;
            C = [cos(RAAN)*cos(theta) - sin(RAAN)*cos(inclination)*sin(theta), -cos(RAAN)*sin(theta) - sin(RAAN)*cos(inclination)*cos(theta), sin(RAAN)*sin(inclination);
            sin(RAAN)*cos(theta) + cos(RAAN)*cos(inclination)*sin(theta), -sin(RAAN)*sin(theta) + cos(RAAN)*cos(inclination)*cos(theta), -cos(RAAN)*sin(inclination);
            sin(inclination)*sin(theta), sin(inclination)*cos(theta), cos(inclination)];
        end

        function phi = flightPathAngleFromENu(e, nu)
            phi = atan(e*sin(nu)/(1 + e*cos(nu)));
        end

        function R_p = RpFromAE(a, E)
            r_p = astroUtilities.periapsisFromAE(a, norm(E));
            R_p = PeriapsisFromRE(r_p, E);
        end

        function R_p = PeriapsisFromRE(r_p, E)
            R_p = r_p*E/norm(E);
        end

        function v = velocityAtApsisFromHR(h, r)
            v = h/r;
        end
    end
end