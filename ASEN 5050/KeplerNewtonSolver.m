classdef KeplerNewtonSolver
    methods(Static)
        function E = KeplerNewtonSolvert(t, a, e, mu)

            n = astroUtilities.meanMotionFroma(a, mu);
            M = n*t;

            E = KeplerNewtonSolver(M, a, e, mu);
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
    end
end