classdef GnssUtilities
    methods(Static)
        function Rn = autoCorrelateFunction(CA)
            Rn = [1:length(CA); ones(1, length(CA))];
            Rn(2, :) = arrayfun(@(x) GnssUtilities.autoCorrelateN(CA, x), Rn(1, :));
        end

        function R = autoCorrelateN(CA, n)
            x = GnssUtilities.CaToX(CA);
            xShiftN = Utilities.leftShift(x, n);
            R = GnssUtilities.correlate(x, xShiftN);
        end

        function AzElRg = AzElRangeToObserver(REcef_observer, RGeod_observer, REcef_satellites)
            REcef_observerSatellites = REcef_satellites - REcef_observer;

            C_Ecef2Enu = GnssUtilities.TranformationMatrix_Ecef2Enu(RGeod_observer);
            REnu_observerSatellites = C_Ecef2Enu*REcef_observerSatellites;

            [e, n, u] = Utilities.decompose(REnu_observerSatellites);
            azimuth = atan2(e, n);
            range = vecnorm(REnu_observerSatellites);
            elevation = asin(u./range);

            AzElRg = [azimuth
                elevation
                range];
        end

        function AzElRg = AzElRangeToObserverEcef(REcef_observer, REcef_satellites)
            RGeod_observer = GnssUtilities.Ecef2Geodetic(REcef_observer);
            AzElRg = GnssUtilities.AzElRangeToObserver(REcef_observer, RGeod_observer, REcef_satellites);
        end

        function AzElRg =  AzElRangeToObserverGeodetic(RGeod_observer, REcef_satellites)
            REcef_observer = GnssUtilities.Geodetic2EcefSurface(RGeod_observer);
            AzElRg = GnssUtilities.AzElRangeToObserver(REcef_observer, RGeod_observer, REcef_satellites);
        end

        function R = correlate(x1, x2)
            R = sum(x1.*x2)/length(x1);
        end

        function Rn = crossCorrelateFunction(CA1, CA2)
            x1 = GnssUtilities.CaToX(CA1);
            x2 = GnssUtilities.CaToX(CA2);
            Rn =  GnssUtilities.crossCorrelateFunctionX(x1, x2);
        end

        function Rn = crossCorrelateFunctionX(x1, x2)
            assert(length(x1)==length(x2))
            Rn = [1:length(x1); ones(1, length(x1))];
            Rn(2, :) = arrayfun(@(x) GnssUtilities.crossCorrelateNX(x1, x2, x), Rn(1, :));
        end

        function R = crossCorrelateN(CA1, CA2, n)
            x1 =  GnssUtilities.CaToX(CA1);
            x2 =  GnssUtilities.CaToX(CA2);
            R = GnssUtilities.crossCorrelateNX(x1, x2, n);
        end

        function R = crossCorrelateNX(x1, x2, n)
            x2 = Utilities.leftShift(x2, n);
            R = GnssUtilities.correlate(x1, x2);
        end

        function x = CaToX(CA)
            x = (CA==0)*1 + (CA==1)*-1;
        end

        function RGeod = Ecef2Geodetic(REcef)
            e = Utilities.eccentricity_earth;

            [x, y, z] = Utilities.decompose(REcef);

            p = norm([x y]);
            r = norm(REcef);

            longitute = atan2(y, x);

            latitude = asin(z/r); % initial guess
            dLatitude = 1;
            while(dLatitude > 1e-8)
                latitude_previous = latitude;
                radiusOfCurveInMeridian = GnssUtilities.radiusOfCurveInMeridian(latitude_previous);
                latitude = atan((z + radiusOfCurveInMeridian*e^2*sin(latitude_previous))/p);

                dLatitude = abs(latitude - latitude_previous);
            end

            RGeod = [latitude; longitute];
        end

        function azElRanges = FilterAboveHorizon(azElRanges)
            azElRanges = azElRanges(:, azElRanges(2, :) >= 0);
        end

        function [G, G1, G2] = generateCA(phaseSelector, bits)
            G1 = ones(1, 10);
            G2 = ones(1, 10);
            [G, G1, G2] = GnssUtilities.generateCAFromRegisters(phaseSelector, G1, G2, bits);
        end

        function [G, G1, G2] = generateCAFromRegisters(phaseSelector, G1, G2, bits)
            G = NaN(1, bits);

            for i = 1:bits
                G1i = G1(10);

                S = [G2(phaseSelector(1)) G2(phaseSelector(2))];
                G2i = mod(sum(S), 2);

                XGi = mod(sum([G1i G2i]), 2);
                G(i) = XGi;

                G1Feedback = [G1(3) G1(10)];
                G1 = [mod(sum(G1Feedback), 2) G1(1:9)];

                G2Feedback = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
                G2 = [mod(sum(G2Feedback), 2) G2(1:9)];
            end
        end

        function REcef = Geodetic2EcefSurface(RGeod)
            [latitudeGeodetic, longitude] = Utilities.decompose(RGeod);
            e_earth = Utilities.eccentricity_earth;
            radiusOfCurveInMeridian = GnssUtilities.radiusOfCurveInMeridian(latitudeGeodetic);

            REcef = [
                radiusOfCurveInMeridian*cos(latitudeGeodetic)*cos(longitude)
                radiusOfCurveInMeridian*cos(latitudeGeodetic)*sin(longitude)
                radiusOfCurveInMeridian*(1 - e_earth^2)*sin(latitudeGeodetic)];
        end

        function A_i = ionosphericDelayAmplitude(alphas, latGeomag_IPP)
            latGeomag_IPPSC = Utilities.rad2semicircle(latGeomag_IPP);
            A_i = alphas(1) + alphas(2)*latGeomag_IPPSC + alphas(3)*latGeomag_IPPSC^2 + alphas(4)*latGeomag_IPPSC^3; % s
            A_i = max(A_i, 0);
        end

        function P_i = ionosphericDelayPeriod(betas, latGeomag_IPP)
            latGeomag_IPPSC = Utilities.rad2semicircle(latGeomag_IPP);
            P_i = betas(1) + betas(2)*latGeomag_IPPSC + betas(3)*latGeomag_IPPSC^2 + betas(4)*latGeomag_IPPSC^3; % s
            P_i = min(P_i, 72000);
        end

        function delay_ionosphericL1Gps = Klobuchar(RGeod_observer, azEl_sat, t_GPS, alphas, betas)
            assert(length(alphas) == 4)
            assert(length(betas) == 4)
            A1 = 5e-9;

            RGeomag_Ipp = GnssUtilities.RGeomagnetic(RGeod_observer, azEl_sat);
            amplitude_ionosphericDelay = ionosphericDelayAmplitude(alphas, RGeomag_Ipp(1));
            period_ionosphericDelay = ionosphericDelayPeriod(betas, RGeomag_Ipp(1));
            
            t_Ipp = GnssUtilities.localTime(t_GPS, long_observer);
            phase_ionosphericDelay = 2*pi*(t_Ipp - 504000)/period_ionosphericDelay;

            slantFactor = 1 + 16*(0.53 - Utilities.rad2semicircle(el_sat));

            if(phase_ionosphericDelay >= pi/2)
                delay_ionosphericL1Gps = slantFactor*A1;
            else
                delay_ionosphericL1Gps = slantFactor*(A1 + amplitude_ionosphericDelay*cos(phase_ionosphericDelay));
            end
        end

        function [latitude, longitude] = latGeodeticLongDeg(R_Ecef)
            rads = GnssUtilities.latGeodeticLong(R_Ecef);
            [latitude, longitude] = rad2deg(rads(:));
        end

        function RGeomag = RGeomagnetic(RGeod_observer, azEl_sat)
            el_earthCenteredSC = 0.0137/(Utilities.rad2semicircle(azEl_sat(2)) + 0.11); % semicircles

            lat_IppSC = Utilities.rad2semicircle(latGeodet_observer) + el_earthCenteredSC*cos(azEl_sat(1)); % semicircles
            lat_IppSC = Utilities.clamp(lat_IppSC, 0.416, -0.416);

            long_IppSC = Utilities.rad2semicircle(long_observer) + el_earthCenteredSC*sin(azEl_sat(1))/cos(semicircle2rad(lat_IppSC)); % semicircles
            long_Ipp = Utilities.semicircle2rad(long_IppSC);

            latGeomagSC = lat_IppSC + 0.064*cos(long_Ipp - 1.617); % semicircles
            latGeomag = Utilities.semicircle2rad(latGeomagSC);

            RGeomag = [latGeomag; long_Ipp];
        end

        function t = localTime(t_GPS, longitude)
            sPerRad = Utilities.days2s(1)/(2*pi);
            t = sPerRad*longitude + t_GPS;

            t = GnssUtilities.normalizeSOfDay(t);
        end

        function t = normalizeSOfDay(t)
            secondsOfADay = Utilities.days2s(1);
            while(t >= secondsOfADay)
                t = t - secondsOfADay;
            end
            while(t < 0)
                t = t + secondsOfADay;
            end
        end

        function C_earth = radiusOfCurveInMeridian(latitudeGeodetic)
            r_earth = Utilities.radius_earth;
            e_earth = Utilities.eccentricity_earth;
            C_earth = r_earth/sqrt(1 - e_earth^2*sin(latitudeGeodetic)^2);
        end

        function SatellitesInView = SatellitesInView(AzElRg_satellites, minElevation)
            elevationRow = 4;
            SatellitesInView = AzElRg_satellites(:, AzElRg_satellites(elevationRow, :) > minElevation);
        end

        function SatellitesInView = SatellitesInViewAtEcef(satellitePositions, REcef_observer, minElevation)
            AzElRg_satellites = [satellitePositions(1:2, :); GnssUtilities.AzElRangeToObserverEcef(REcef_observer, satellitePositions(3:5, :))];
            SatellitesInView = GnssUtilities.SatellitesInView(AzElRg_satellites, minElevation);
        end

        function SatellitesInView = SatellitesInViewAtGeodetic(satellitePositions, RGeod_observer, minElevation)
            AzElRg_satellites = [satellitePositions(1:2, :); GnssUtilities.AzElRangeToObserverLatLong(RGeod_observer, satellitePositions)];
            SatellitesInView = GnssUtilities.SatellitesInView(AzElRg_satellites, minElevation);
        end

        function C_Ecef2Enu = TranformationMatrix_Ecef2Enu(RGeod)
            [lat, long] = Utilities.decompose(RGeod);

            C_Ecef2Enu = [-sin(long) cos(long) 0
                -sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat)
                cos(lat)*cos(long) cos(lat)*sin(long) sin(lat)];
        end
    end
end