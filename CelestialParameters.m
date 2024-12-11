classdef CelestialParameters
    properties(Constant)
        kmPerAu = 149597870.7; % km/AU
        universalGravitationalConstant = 6.6743e-20; % km^3/(kg*s^2)
        
        %% Sun
        gravityParameter_sun = 1.32712428e11; % km^3/s^2
        
        %% Earth
        accelerationGravity_earthSurface = 9.81; % m/s^2
        eccentricity_earth = sqrt(2*CelestialParameters.flattening_earth - CelestialParameters.flattening_earth^2); % WGS-84 fundamental parameter
        flattening_earth = 1/298.257223563; % WGS-84 fundamental parameter
        gravityParameter_earth = 3.986012e5 % km^3/s^2
        j2_earth = 1082.64e-6;
        radius_earth = 6378137; % m, WGS-84 fundamental parameter
        semiMajorAxis_earthAu = 1.0000010178; % AU
        rotationRate_earth = 7.2921151467e-5; % rad/s

        %% Moon
        radius_Moon = 1738; % km
        gravityParameter_moon = 4902.799; % km^3/s^2
       
        %% Mars
        gravityParameter_mars = 4.305e4; % km^3/s^2
        j2_mars = 1960e-6; % Wahrheit 2009
        radius_Mars = 3397.2; % km
        rotationPeriod_marsDays = 1.025957; % days
        rotationPeriod_mars = CelestialParameters.rotationPeriod_marsDays*24*60*60; % s
        rotationRate_mars = 2*pi/CelestialParameters.rotationPeriod_mars; % rad/s
        semiMajorAxis_mars = 1.523679342*CelestialParameters.kmPerAu; % km

        %% Jupiter
        gravityParameter_jupiter = 1.268e8; % km^3/s^2
        semiMajorAxis_jupiter = 5.202603191*CelestialParameters.kmPerAu; % km
        
        %% Saturn
        gravityParameter_saturn = 3.794e7; % km^3/s^2
        semiMajorAxis_saturnAu = 9.554909595; % AU
    end

    methods(Static)
        function km = au2km(au)
            km = au*CelestialParameters.kmPerAu;
        end

        function au = km2au(km)
            au = km/CelestialParameters.kmPerAu;
        end

        function mu = gravitationalParameter(m)
            mu = CelestialParameters.universalGravitationalConstant*m;
        end
    end
end