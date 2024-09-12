classdef astroUtilities
    methods(Static)

        function nrg_mech = mechanicalNRGFromA(a,mu)
            nrg_mech = -mu/(2*a);
        end

        function nrg_mech = mechanicalNRGFromRV(R, V, mu)
            nrg_mech = astroUtilities.mechanicalNRGFromrv(norm(R), norm(V), mu);
        end

        function nrg_mech = mechanicalNRGFromrv(r, v, mu)
            nrg_mech = v^2/2 - mu/r;
        end

        function v = vFromMechanicalNRG(nrg_mech, mu, r)
            v = sqrt(2*(nrg_mech + mu/r));
        end

        function p = semiLatusRectumFromH(h, mu)
            p = h^2/mu;
        end

        function e = eccentricityFromSemiLatusRectum(p, a)
            e = sqrt(1 - p/a);
        end

        function e = eccentricityFromRVH(R, V, H, mu)
            e = cross(V, H)/mu - R/norm(R);
        end

        function a = semiMajorAxisFromMechanicalNRG(nrg_mech, mu)
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
    end
end