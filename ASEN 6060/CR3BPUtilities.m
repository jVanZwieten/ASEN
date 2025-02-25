classdef CR3BPUtilities
    properties(Constant)
        Gtil = CelestialParameters.universalGravitationalConstant;

        %% Earth-Luna
        GMtil_earth = CelestialParameters.gravityParameter_earth;
        Mtil_earth = CR3BPUtilities.GMtil_earth/CR3BPUtilities.Gtil;

        GMtil_luna = CelestialParameters.gravityParameter_luna;
        Mtil_luna = CR3BPUtilities.GMtil_luna/CR3BPUtilities.Gtil;

        mStar_earthLuna = CR3BPUtilities.characteristicMass(CR3BPUtilities.Mtil_earth, CR3BPUtilities.Mtil_luna);
        mu_earthLuna = CR3BPUtilities.massRatio(CR3BPUtilities.Mtil_earth, CR3BPUtilities.Mtil_luna);

        lStar_earthLuna = CelestialParameters.semiMajorAxis_luna;
        tStar_earthLuna = CR3BPUtilities.characteristicTime(CR3BPUtilities.lStar_earthLuna, CR3BPUtilities.mStar_earthLuna);

        lagrangePoints_earthLuna = CR3BPUtilities.lagrangePoints(CR3BPUtilities.mu_earthLuna);

        %% Sol-Earth
        GMtil_sol = CelestialParameters.gravityParameter_sol;
        Mtil_sol = CR3BPUtilities.GMtil_sol/CR3BPUtilities.Gtil;

        mStar_solEarth = CR3BPUtilities.characteristicMass(CR3BPUtilities.Mtil_sol, CR3BPUtilities.Mtil_earth);
        mu_solEarth = CR3BPUtilities.massRatio(CR3BPUtilities.Mtil_sol, CR3BPUtilities.Mtil_earth);

        lStar_solEarth = CelestialParameters.semiMajorAxis_earth;
        tStar_solEarth = CR3BPUtilities.characteristicTime(CR3BPUtilities.lStar_solEarth, CR3BPUtilities.mStar_solEarth);
    end
    methods(Static)
        function mStar = characteristicMass(m1, m2)
            mStar = m1 + m2;
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

            rVec = [x; y; z];
            [r_1, r_2] = CR3BPUtilities.r_primaries(rVec, mu);

            xDotDot = 2*yDot + x - (1 - mu)*(x + mu)/r_1^3 - mu*(x - 1 + mu)/r_2^3;
            yDotDot = -2*xDot + y - (1 - mu)*y/r_1^3 - mu*y/r_2^3;
            zDotDot = -(1 - mu)*z/r_1^3 - mu*z/r_2^3;

            dXdt = [xDot; yDot; zDot; xDotDot; yDotDot; zDotDot];
        end

        function dXdt = cr3bpEomEquilibriumLinearized(t, X, UstarMat)
            xi = X(1);
            eta = X(2);
            delta = X(3);
            xiDot = X(4);
            etaDot = X(5);
            deltaDot = X(6);


            UstarPlanarMat = UstarMat(1:2, 1:2);
            OmegaMat = [0 2; -2 0];
            planarAcceleration = [UstarPlanarMat, OmegaMat]*[xi; eta; xiDot; etaDot];

            Ustar_zz = UstarMat(3, 3);
            deltaDotDot = Ustar_zz*delta;

            dXdt = [xiDot; etaDot; deltaDot; planarAcceleration; deltaDotDot];
        end

        function [t, X] = integrateCr3bp(xVec_0, mu, tSpan, tolerancePower)
            if(length(tSpan) == 1)
                tSpan = [0, tSpan];
            end

            solverOptions = odeset('stats', 'off', 'RelTol', 10^tolerancePower, 'AbsTol', 10^tolerancePower);
            [t, X] = ode45(@(t, X) CR3BPUtilities.cr3bpEom(t, X, mu), tSpan, xVec_0, solverOptions);
        end

        function [t, X, t_f, xVec_f, index_f] = integrateCr3bpUntil(xVec_0, mu, tSpan, tolerancePower, terminationEvent)
            solverOptions = odeset('Events', terminationEvent, 'stats', 'off', 'RelTol', 10^tolerancePower, 'AbsTol', 10^tolerancePower);
            [t, X, t_f, xVec_f, index_f] = ode45(@(t, X) CR3BPUtilities.cr3bpEom(t, X, mu), tSpan, xVec_0, solverOptions);
        end

        function [t, X, t_f, xVec_f, index_f] = integrateLinearizedCr3bp(rVec_equilibrium, deltaXVec_0, mu, tSpan, tolerancePower)
            if(length(rVec_equilibrium) == 2)
                rVec_equilibrium = [rVec_equilibrium; 0];
            end
            assert(length(rVec_equilibrium) == 3)
            if(length(deltaXVec_0) == 4)
                deltaXVec_0 = [deltaXVec_0(1:2); 0; deltaXVec_0(3:4); 0];
            end
            assert(length(deltaXVec_0) == 6)

            UstarMat = CR3BPUtilities.potentialDoubleDerivativeAtMat(rVec_equilibrium, mu);
            solverOptions = odeset('stats', 'off', 'RelTol', 10^tolerancePower, 'AbsTol', 10^tolerancePower);
            [t, X] = ode45(@(t, X) CR3BPUtilities.cr3bpEomEquilibriumLinearized(t, X, UstarMat), tSpan, deltaXVec_0, solverOptions);
        end

        function C = jacobiConstant(XVec, mu)
            rVec = XVec(1:3);

            if length(XVec) == 6
                vVec = XVec(4:6);
            else
                vVec = zeros(3, 1);
            end

            x = rVec(1);
            y = rVec(2);
            [r_1, r_2] = CR3BPUtilities.r_primaries(rVec, mu);

            C = x^2 + y^2 + 2*(1 - mu)/r_1 + 2*mu/r_2 - dot(vVec, vVec);
        end

        function LVecs = lagrangePoints(mu)
            equilibriumEquation = @(x) x - (1 - mu)*(x + mu)/abs(x + mu)^3 - mu*(x - 1 + mu)/abs(x - 1 + mu)^3;

            guessDeltaFactor = 0.1;
            initialGuess = [1 - mu*(1 + guessDeltaFactor), 1 - mu*(1 - guessDeltaFactor), -mu*(1 + guessDeltaFactor)];

            options = optimoptions('fsolve', 'Display', 'off');
            roots = NaN(1, length(initialGuess));
            for i = 1:length(initialGuess)
                roots(i) = fsolve(equilibriumEquation, initialGuess(i), options);
            end

            triangularPoints = {[.5 - mu; sqrt(3)/2; 0], [.5 - mu; -sqrt(3)/2; 0]};

            LVecs = cellfun(@(x) [x; 0; 0], num2cell(roots), 'UniformOutput', false);
            LVecs = [LVecs, triangularPoints];
        end

        function lagrangePlot(lagrangePoints)
            figure
            hold on
            for i = 1:length(lagrangePoints)
                plot(lagrangePoints{i}(1), lagrangePoints{i}(2), 'o')
            end
            xlabel('x')
            ylabel('y')
            legend('L1', 'L2', 'L3', 'L4', 'L5')
            hold off
        end

        function [lagrangePointStabilitiesT] = lagrangeStability(mu)
            lagrangePoints = CR3BPUtilities.lagrangePoints(mu);

            lagrangePointStabilitiesT = table('Size', [5, 6], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'InPlaneEigenvalue1', 'InPlaneEigenvalue2', 'InPlaneEigenvalue3', 'InPlaneEigenvalue4', 'OutPlaneEigenvalue1', 'OutPlaneEigenvalue2'}, 'RowNames', {'L1', 'L2', 'L3', 'L4', 'L5'});
            for i = 1:length(lagrangePoints)
                lagrangePoint = lagrangePoints{i};
                [Ustar_xx, Ustar_xy, Ustar_yy, Ustar_zz] = CR3BPUtilities.potentialDoubleDerivativeAt(lagrangePoint, mu);

                a = -4 + Ustar_xx + Ustar_yy;
                b = sqrt((4 - Ustar_xx - Ustar_yy)^2 - 4*(Ustar_xx*Ustar_yy - Ustar_xy^2));
                Lambdas = [(a + b)/2, (a - b)/2];

                lambdas_inPlane = [];
                for Lambda = Lambdas
                    lambdas_inPlane = [lambdas_inPlane, sqrt(Lambda), -sqrt(Lambda)];
                end

                lambdas_outPlane = [sqrt(Ustar_zz), -sqrt(Ustar_zz)];

                lagrangePointStabilitiesT(i, :) = {lambdas_inPlane(1), lambdas_inPlane(2), lambdas_inPlane(3), lambdas_inPlane(4), lambdas_outPlane(1), lambdas_outPlane(2)};
            end
        end

        function xiDotVec = librationPlanarOrbitVelocity(rVec_equilibrium, xiVec, lambda_3, mu)
            xi_0 = xiVec(1);
            eta_0 = xiVec(2);

            Ustar_xx = CR3BPUtilities.potentialDoubleDerivativeAt(rVec_equilibrium, mu);

            alpha_3 = (lambda_3^2 - Ustar_xx)/2/lambda_3;

            xiDot = lambda_3/alpha_3*eta_0;
            etaDot = alpha_3*lambda_3*xi_0;

            xiDotVec = [xiDot; etaDot];
        end

        function P = librationPeriod(eigenvalue)
            P = 2*pi/abs(imag(eigenvalue));
        end

        function A = linearizedA(rVec_equilibrium, mu)
            UstarMat = CR3BPUtilities.potentialDoubleDerivativeAtMat(rVec_equilibrium, mu);
            Omega = [0, 2, 0; -2, 0, 0; 0, 0, 0];
            A = [zeros(3), eye(3); UstarMat, Omega];
        end

        function A_2D = linearizedA2d(rVec_equilibrium, mu)
            UstarMat = CR3BPUtilities.potentialDoubleDerivativeAtMat(rVec_equilibrium, mu);
            Omega = [0, 2; -2, 0];
            A_2D = [zeros(2), eye(2); UstarMat(1:2, 1:2), Omega];
        end

        function mu = massRatio(m1, m2)
            mu = m2/(m1 + m2);
        end

        function plotLinearizedLibrationOrbit(X)
            figure
            plot(X(:, 1), X(:, 2), 'DisplayName', 'Trajectory')
            hold on
            plot(0, 0, 'rd', 'MarkerFaceColor', 'r', 'DisplayName', 'L1')
            legend
            xlabel('\xi')
            ylabel('\eta')
            axis equal
            hold off
        end

        function plotCr3bpTrajectory(X, mu, p1Name, p2Name)
            figure
            hold on
            colorMap = flipud(jet(size(X, 1)));
            scatter3(X(:,1), X(:,2), X(:,3), 2, colorMap, 'filled', 'DisplayName', 'Trajectory');
            if ~isempty(p1Name)
                plot3(-mu, 0, 0, 'bo', 'MarkerSize', 10, 'DisplayName', p1Name);
            end
            if ~isempty(p2Name)
                plot3(1 - mu, 0, 0, 'ro', 'MarkerSize', 5, 'DisplayName', p2Name);
            end
            view(3)
            legend
            xlabel('x');
            ylabel('y');
            zlabel('z');
        end

        function plotZeroVelocityCurve(C, mu, p1Name, p2Name)
            xRange = linspace(-1.5, 1.5, 400);
            yRange = linspace(-1.5, 1.5, 400);
            [X, Y] = meshgrid(xRange, yRange);

            Z = zeros(size(X));
            for i = 1:length(xRange)
                for j = 1:length(yRange)
                    rVec = [X(i, j); Y(i, j); 0];
                    Z(i, j) = CR3BPUtilities.jacobiConstant(rVec, mu);
                end
            end

            figure;
            contourf(X, Y, Z, [-Inf C, max(Z(:))], 'LineColor', 'none', 'DisplayName', 'Forbidden Region');
            colormap([0.5 0.5 0.5; 1 1 1]); % Set colormap to grey for forbidden region and white for allowable region
            caxis([-1 1]);
            xlabel('x');
            ylabel('y');
            axis equal;
            grid on;
            hold on;

            if ~isempty(p1Name)
                plot(-mu, 0, 'bo', 'MarkerSize', 10, 'DisplayName', p1Name);
            end
            if ~isempty(p2Name)
                plot(1 - mu, 0, 'ro', 'MarkerSize', 6, 'DisplayName', p2Name);
            end

            legend
        end

        function [r_1, r_2] = r_primaries(rVec, mu)
            [rVec_1, rVec_2] = CR3BPUtilities.rVec_primaries(rVec, mu);
            r_1 = norm(rVec_1);
            r_2 = norm(rVec_2);
        end

        function [rVec_1, rVec_2] = rVec_primaries(rVec, mu)
            if(length(rVec) == 2)
                rVec = [rVec; 0];
            end
            assert(length(rVec) == 3)

            rVec_1 = rVec - [-mu; 0; 0];
            rVec_2 = rVec - [1 - mu; 0; 0];
        end

        function [Ustar_xx, Ustar_xy, Ustar_yy, Ustar_zz] = potentialDoubleDerivativeAt(rVec, mu)
            x = rVec(1);
            y = rVec(2);
            [r_1, r_2] = CR3BPUtilities.r_primaries(rVec, mu);

            Ustar_xx = 1 - (1 - mu)/r_1^3 - mu/r_2^3 + 3*(1 - mu)*(x + mu)^2/r_1^5 + 3*mu*(x - 1 + mu)^2/r_2^5;
            Ustar_xy = 3*(1 - mu)*(x + mu)*y/r_1^5 + 3*mu*(x - 1 + mu)*y/r_2^5;
            Ustar_yy = 1 - (1 - mu)/r_1^3 - mu/r_2^3 + 3*(1 - mu)*y^2/r_1^5 + 3*mu*y^2/r_2^5;
            Ustar_zz = -(1 - mu)/r_1^3 - mu/r_2^3;
        end

        function UstarMat = potentialDoubleDerivativeAtMat(rVec_equilibrium, mu)
            [Ustar_xx, Ustar_xy, Ustar_yy, Ustar_zz] = CR3BPUtilities.potentialDoubleDerivativeAt(rVec_equilibrium, mu);
            UstarMat = [Ustar_xx, Ustar_xy, 0; Ustar_xy, Ustar_yy, 0; 0, 0, Ustar_zz];
        end
    end
end