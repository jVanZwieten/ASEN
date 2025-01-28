classdef CR3BPUtilities
    properties(Constant)
    end
    methods(Static)
        function mStar = characteristicMass(m1, m2)
            mStar = m1 + m2;
        end

        function mu = massRatio(m1, m2)
            mu = m2/(m1 + m2);
        end

        function tStar = characteristicTime(lStar, mStar)
            Gtil = CelestialParameters.universalGravitationalConstant;
            tStar = sqrt(lStar^3/(Gtil*mStar));
        end

        function dXdt = cr3bpEom(t, X, mu)
            x = X(1);
            y = X(2);
            z = X(3);
            xDot = X(4);
            yDot = X(5);
            zDot = X(6);

            R = [x; y; z];
            r_1 = norm(R + [mu; 0; 0]);
            r_2 = norm(R + [mu - 1; 0; 0]);

            xDotDot = 2*yDot + x - (1 - mu)*(x + mu)/r_1^3 - mu*(x - 1 + mu)/r_2^3;
            yDotDot = -2*xDot + y - (1 - mu)*y/r_1^3 - mu*y/r_2^3;
            zDotDot = -(1 - mu)*z/r_1^3 - mu*z/r_2^3;

            dXdt = [xDot; yDot; zDot; xDotDot; yDotDot; zDotDot];
        end
    end
end