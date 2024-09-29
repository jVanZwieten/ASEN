function C_ECEF2ENU = ECEF2ENU(lat_deg,lon_deg)
    lat = deg2rad(lat_deg);
    long = deg2rad(lon_deg);
    C_ECEF2ENU = [-sin(long) cos(long) 0
        -sin(lat)*cos(long) -sin(lat)*sin(long) cos(lat)
        cos(lat)*cos(long) cos(lat)*sin(long) sin(lat)];
end