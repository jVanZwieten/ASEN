classdef lambertUtilities
    methods(Static)
        function [c, s] = lambertGeometricQuantities(r_1, r_2, deltaNu)
            c = sqrt(r_1^2 + r_2^2 - 2*r_1*r_2*cos(deltaNu));
            s = (r_1 + r_2 + c)/2;
        end

        function [TOF, alpha, beta, n] = lambertsEquation(c, s, a, mu, longTof, longWay)
            n = sqrt(mu/a^3);

            alpha = 2*asin(sqrt(s/(2*a)));
            if(longTof) alpha = 2*pi - alpha; end

            beta = 2*asin(sqrt((s - c)/(2*a)));
            if(longWay) beta = -beta; end
            
            TOF = lambertUtilities.lambertsEquationGeometric(alpha, beta, n);
        end

        function TOF = lambertsEquationGeometric(alpha, beta, n)
            TOF = ((alpha - beta) - (sin(alpha) - sin(beta)))/n;
        end

        function TOF_min = lambertsEquationMin(c, s, longWay, mu)
            a_min = s/2
            n_min = sqrt(mu/a_min^3)
            alpha_min = pi
            beta_min = 2*asin(sqrt((s - c)/s))
            if(longWay) beta_min = -beta_min; end
            
            TOF_min = lambertUtilities.lambertsEquationGeometric(alpha_min, beta_min, n_min);
        end
        
        function TOF_parabolic = lambertsEquationParabolic(c, s, mu, longWay)
            longWayFactor = -1 + 2*longWay;
            TOF_parabolic = (sqrt(2/mu)*(s^(3/2) + longWayFactor*(s - c)^(3/2)))/3;
        end

        function TOF = lambertsEquationHyperbolic(c, s, a, mu, longWay)
            longWayFactor = -1 + 2*longWay;

            n = sqrt(mu/a^3);
            alpha = 2*asinh(sqrt(s/(2*abs(a))));
            beta = 2*asinh(sqrt((s - c)/(2*abs(a))));
            
            TOF = (sinh(alpha) - alpha + longWayFactor*(sinh(beta) - beta))/n;
        end

        function p = semiLatusRectum(a, c, s, r_1, r_2, alpha, beta)
            p = 4*a*(s - r_1)*(s - r_2)/c^2*sin((alpha + beta)/2)^2;
        end
    end
end