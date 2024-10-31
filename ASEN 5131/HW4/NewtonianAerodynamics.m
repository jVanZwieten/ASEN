function [body, aerodynamicCoefficients] = NewtonianAerodynamics(shapeFileName, alphaA, referenceArea)
    % NEWTONIANAERODYNAMICS Computes the pressure, lift, and drag coefficients of the shape defined by the input file for all input angles of attack alpha.
    %
    %   shapeFileName: should be in .dat format
    %   alphaA: array of angles of attack, alphas, to calculate aerodynamics on (radians)
    %
    %   XYZ body-centered cartesian frame:
    %   X - roll axis
    %   Y - yaw axis
    %   Z - pitch axis
    %
    %   airflow velocity set to the +roll (x) axis, rotated -pitch (z) axis

    verticies = readmatrix(shapeFileName);
    n_triangles = size(verticies, 1)/3;

    triangles = Triangle3d.empty(n_triangles, 0);
    triangleConnectionMatrix = zeros(n_triangles,3);
    for i = 1:n_triangles
        j = (i - 1)*3 + 1;
        triangles(i) = Triangle3d(verticies(j, :)', verticies(j + 1, :)', verticies(j + 2, :)');
        triangleConnectionMatrix(i, :) = [ j, j + 1, j + 2 ];
    end
    
    body = table(triangles');
    body.nHats = [triangles.nHat]';
    body.areas = [triangles.A]';
    
    aerodynamicCoefficients = table([], [], [], 'VariableNames', {'alpha', 'C_L', 'C_D'});

    for alpha = alphaA
        uHat = airflowVector(alpha);

        alphaDegrees = rad2deg(alpha);
        alphaName = num2str(abs(alphaDegrees));
        C_pColumnName = strcat('C_p_', alphaName);

        body.(C_pColumnName) = unmodifiedNewtonianCoeffeicentOfPressure(uHat, body.nHats')';
        [C_L, C_D] = aeroforceCoeffiecients(uHat, body.(C_pColumnName), body.nHats, body.areas, referenceArea);
        
        aerodynamicCoefficients = [aerodynamicCoefficients; {alphaDegrees, C_L, C_D}];

        figure
        trisurf(triangleConnectionMatrix, verticies(:,1), verticies(:,3), verticies(:,2), body.(C_pColumnName));
        title('Cp');
        xlabel('x (m)'); ylabel('z (m)'); zlabel('y (m)');
        colorbar; caxis([ 0 2 ]);
    end
    
    function uHat = airflowVector(alpha)
        % AIRFLOWVECTOR Returns the unit vector, uHat, representing the direction of airflow with respect to the body for a given angle of attack.
        %   alpha: angle of attack of the body, the difference in orientation between the body and airflow (radians)
    
        uHat = [cos(alpha); sin(alpha); 0];
        assert(norm(uHat) == 1)
    end
    
    function C_p = unmodifiedNewtonianCoeffeicentOfPressure(uHat, nHat)
        % UNMODIFIEDNEWTONIANCOEFFEICENTOFPRESSURE Returns the coefficient of pressure, C_p, of an element of analysis based on unmodified Newtonian aerodynamic analysis.
        %   Assumes that element is flat plate with orientation nHat.
        %   If angle(uHat, nHat) > pi/2, the element is not facing the flow, and is assumed to be shadowed by another element, therefore C_p = 0.
        %
        %   uHat: unit vector representing direction of airflow, 3x1
        %   nHat: unit vector representing the outward orientation of the element, pointing outward from the element face, 3xn
    
        cosphi = uHat' * nHat;
        cosphi = min(cosphi, 0); % defends against shadowing
        C_p = 2*cosphi.^2;
        assert(all(C_p <= 2 & C_p >= 0))
    end
    
    function [C_L, C_D] = aeroforceCoeffiecients(uHat, C_ps, nHats, As, A_r)
        FV_pNet = sum(-nHats .* C_ps .* As, 1)';
    
        LAxis = liftAxis(uHat);
    
        F_L = dot(FV_pNet, LAxis);
        C_L = F_L/A_r;
        
        F_D = dot(FV_pNet, uHat);
        C_D = F_D/A_r;
    end
    
    function LHat = liftAxis(uHat)
        yAxis = [0; 1; 0];
        LHat = cross(uHat, cross(yAxis, uHat));
        LHat = unitVector(LHat);
    end

    function Ahat = unitVector(A)
        % UNITVECTOR returns the unit vector(s) of input
        % If n vectors, input as m x n
        Ahat = A./vecnorm(A);
    end
end