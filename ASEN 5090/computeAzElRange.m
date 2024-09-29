function [azimuth, elevation, range] = computeAzElRange(userECEF, satECEF)
    REcef = satECEF - userECEF;
    [lat, long] = GnssUtilities.latGeodeticLong(REcef);

    C_ECEF2ENU = ECEF2ENU(lat, long);
    [e, n, u] = Utilities.decompose(C_ECEF2ENU*REcef);

    azimuth = atan2(e, n);
    range = norm([e n u]);
    elevation = asin(u/range);
end