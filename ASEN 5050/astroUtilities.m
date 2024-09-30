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

        function r = conicEquationp(p, e, nu)
            r = p/(1 + e*cos(nu));
        end
        
        function C = DirectionCosineMatrix(RAAN, AOP, trueAnomaly, inclination)
            theta = AOP + trueAnomaly;
            C = [cos(RAAN)*cos(theta) - sin(RAAN)*cos(inclination)*sin(theta), -cos(RAAN)*sin(theta) - sin(RAAN)*cos(inclination)*cos(theta), sin(RAAN)*sin(inclination);
                sin(RAAN)*cos(theta) + cos(RAAN)*cos(inclination)*sin(theta), -sin(RAAN)*sin(theta) + cos(RAAN)*cos(inclination)*cos(theta), -cos(RAAN)*sin(inclination);
                sin(inclination)*sin(theta), sin(inclination)*cos(theta), cos(inclination)];
        end

        function E = eccentricAnomalyFromaer(a, e, r)
            E = acos((a - r)/(a*e));
        end

        function E = eccentricAnomalyFromnu(nu, e)
            E = 2*atan(sqrt((1 - e)/(1 + e))*tan(nu/2));
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

        function f = fFunction(r, p, deltaNu)
            f = 1 - r/p*(1 - cos(deltaNu));
        end
        
        function fDot = fDotFunction(mu, p, deltaNu, r, r_0)
            fDot = sqrt(mu/p)*tan(deltaNu/2)*((1 - cos(deltaNu))/p - 1/r - 1/r_0);
        end

        function g = gFunction(r, r_0, mu, p, deltaNu)
            g = r*r_0/sqrt(mu*p)*sin(deltaNu);
        end

        function gDot = gDotFunction(r_0, p, deltaNu)
            gDot = 1 - r_0/p*(1 - cos(deltaNu));
        end

        function i = inclinationFromH(H)
            i = astroUtilities.inclinationFromHhat(Utilities.UnitVector(H));
        end
        
        function i = inclinationFromHhat(Hhat)
            i = acos(Hhat(3));
        end

        function E = KeplerNewtonSolvert(t, a, e, mu)

            n = astroUtilities.meanMotionFroma(a, mu);
            M = n*t;

            E = KeplerNewtonSolver(M, a, e, mu)
        end

        function E = KeplerNewtonSolver(M, a, e, mu)
            tolerance = 1e-6;
            E = M; % initial guess

            delta = inf;
            while(abs(delta) > tolerance)
                delta = astroUtilities.g(E, e, M)/astroUtilities.dgdE(E, e);
                E = E - delta;
            end
        end

        function gE = g(E, e, M)
            gE = astroUtilities.meanAnomalyFromE(E, e) - M;
        end

        function dgdE = dgdE(E, e)
            dgdE = 1 - e*cos(E);
        end
        
        function N = LineOfNodesFromH(H)
            N = [-H(2); H(1); 0];
        end

        function M = meanAnomalyFromE(E, e)
            M = E + e*sin(E);
        end

        function n = meanMotionFroma(a, mu)
            n = sqrt(mu/a^3);
        end

        function R = positionfg(R_0, V_0, f, g)
            R = f*R_0 + g*V_0;
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
        
        function p = semiLatusRectumFromh(h, mu)
            p = h^2/mu;
        end

        function a = semiMajorAxisFrompe(p, e)
            a = p/(1 - e^2);
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

        function dt = timeBetweenEccentricAnomalies(E_1, E_2, e, n)
            dt = 1/n*((E_2 - e*sin(E_2)) - (E_1 - e*sin(E_1)));
        end

        function dt = timeBetweenEccentricAnomaliesa(E_1, E_2, a, e, mu)
            n = astroUtilities.meanMotionFroma(mu, a);
            dt = astroUtilities.timeBetweenEccentricAnomalies(E_1, E_2, n, e);
        end

        function dt = timeBetweenTrueAnomalies(nu_1, nu_2, e, n)
            E_1 = astroUtilities.eccentricAnomalyFromnu(nu_1, e);
            E_2 = astroUtilities.eccentricAnomalyFromnu(nu_2, e);
            dt = astroUtilities.timeBetweenEccentricAnomalies(E_1, E_2, e, n);
        end

        function dt = timeBetweenTrueAnomaliesa(nu_1, nu_2, a, e, mu)
            n = astroUtilities.meanMotionFroma(mu, a);
            dt = astroUtilities.timeBetweenTrueAnomalies(nu_1, nu_2, e, n);
        end

        function t = timePeriapsisToEccentricAnomaly(E, e, n)
            t = (E - e*sin(E))/n;
        end

        function nu = trueAnomalyFromRE(R, E)
            nu = acos(dot(E, R)/(norm(E)*norm(R)));
        end

        function nu = trueAnomalyFromhre(h, r, e, mu)
            nu = acos((h^2/r/mu - 1)/e);
        end

        function nu = trueAnomalyFrompre(p, r, e)
            nu = acos((p/r - 1)/e);
        end

        function nu = trueFromEccentricAnomaly(E, e)
            nu = 2*atan(sqrt((1 + e)/(1 - e))*tan(E/2));
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

        function V = velocityfDotgDot(R_0, V_0, fDot, gDot)
            V = fDot*R_0 + gDot*V_0;
        end
    end
end