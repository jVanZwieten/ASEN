function [health,satPos,satClkCorr] = alm2pos(almanac,t_input,prn)

%==========================================================================
%==========================================================================
% [health,x] = alm2pos(ephem_all,t_input,prn)
%
% Calculates the position from an ephemeris 
%  matrix (see read_GPSbroadcast.m).  The input ephem_all can 
%  be generated by the read_GPSyuma.m or read_GPSbroadcast.m function.
%
% Modified by Y. Wang 9/20/2023 for init. satClkCorr & cleaning up
% Modified by P. Axelrad 9/10/2018 to remove extra functionality
% Author: Ben K. Bradley
% Date: 07/19/2009
%
%
% INPUT:               Description                                  Units
%
%  ephem_all    - matrix of gps satellite orbit parameters           (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: blank (zero), delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: blank (zero), i_rate, rate of change of inclination angle, rad/s
%                 col11: blank (zero), Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: blank (zero), Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: blank (zero), Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: blank (zero), Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: blank (zero), Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: blank (zero), Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris/almanac (seconds into GPS week)
%                 col18: blank (zero), IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe)
%                 col20: blank (zero), Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: blank (zero), Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: blank (zero), Timing Group Delay (TGD), seconds
%                 col25: health, satellite health (0=good and usable)
%
%
%  t_input      - GPS times to calculate values at                 [WN TOW] (nx2)
%  prn          - PRN to compute values for (one satellite only)                       
%
%
% OUTPUT:       
%    
%  health       - health of satellite (0=good)                              (nx1)
%  x            - position of satellite (ECEF)                  [x y z]   m (nx3)
%                                     
%

% Coupling:
%
%   mean2eccentric.m
%
% References:
% 
%   [1] Interface Control Document: IS-GPS-200D
%         < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%   [2] Zhang, J., et.all. "GPS Satellite Velocity and Acceleration
%         Determination using the Broadcast Ephemeris". The Journal of
%         Navigation. (2006), 59, 293-305.
%            < http://journals.cambridge.org/action/displayAbstract;jsess ...
%                ionid=C6B8C16A69DD7C910989C661BAB15E07.tomcat1?fromPage=online&aid=425362 >
%
%   [3] skyplot.cpp by the National Geodetic Survey
%          < http://www.ngs.noaa.gov/gps-toolbox/skyplot/skyplot.cpp >
%
%
%
%  2015/01/22  B.K. Bradley - the capability to look for updated ephem
%                              entries that occur at odd times within each
%                              2hr window has been commented out in this 
%                              function and added to read_GPSbroadcast.m
%                              instead. This moves the computational
%                              overhead to the reading which only occurs
%                              once.
%
%  2021/09/10 P. Axelrad - removed variables not included in the almanac
%  model.
%
%==========================================================================
%==========================================================================

% Load GPS Accepted WGS-84 Constants 
muE = 3.986005e14;     % WGS-84 value, m^3/s^2
wE  = 7.2921151467e-5; % WGS-84 value, rad/s 
c   = 2.99792458e8;    % SI speed of light, m/s

% Initialize Output Variables for Speed 
sz         = size(t_input,1);
satPos     = ones(sz,3) * NaN;
health     = ones(sz,1) * NaN; 
satClkCorr = ones(sz,1) * NaN; 

% Pull out ephemerides for the selected PRN
kk  = find(almanac(:,1) == prn);  
sat_alm = almanac(kk,:);        

% If no matching PRN found, returning data will be NaNs
if isempty(kk),return,end 

% Start Main Calculation Loop 
for tt = 1:sz % loop through all input times
    % Pull out variables from the almanac matrix
    %======================================================================   
    Toe = sat_alm(17);
    gps_wk = sat_alm(19);
    dt  = (t_input(tt,1) - gps_wk)*604800 + (t_input(tt,2) - Toe); % seconds difference from epoch
    a   = sat_alm(5)^2;           % semimajor axis, sqrt(a) = gps_ephem_all(:,5) (meters)
    ecc = sat_alm(4);             % eccentricity
    n0  = sqrt(muE/a^3);          % nominal mean motion (rad/s)
    n   = n0 + sat_alm(3);        % corrected mean motion, delta_n = gps_ephem_all(:,3)
    M   = sat_alm(2) + n*dt;      % mean anomaly, M0 = gps_ephem_all(:,2)
    inc = sat_alm(7) ;            % inclination
    perigee  = sat_alm(8);     % argument of perigee

    % Compute true and eccentric anomaly...
    %======================================================================        
    % Compute Eccentric Anomaly, rad
    E    = mean2eccentric(M,ecc);
    cosE = cos(E);  
    sinE = sin(E);
    % Compute true anomaly, rad
    nu    = atan2( sqrt(1 - ecc*ecc).*sinE,  cosE-ecc ); 
    % Compute the argument of latitude, rad 
    u = nu + perigee;  % true anomaly + argument of perigee

    % Compute radius and inclination
    %======================================================================   
    r   = a * (1 - ecc*cosE) ;                        % corrected radius  

    cosu = cos(u);    
    sinu = sin(u);  

    % Compute satellite position in orbital plane
    %======================================================================
    xo = r * cosu;    % satellite x-position in orbital plane
    yo = r * sinu;    % satellite y-position in orbital plane

    % Corrected longitude of ascending node for node rate and Earth rotation
    %======================================================================
    node = sat_alm(6) + (sat_alm(9) - wE)*dt -  (wE * Toe); 

    % Calculate GPS Satellite Position in ECEF (m)
    %======================================================================
    cosi = cos(inc);    sini = sin(inc);
    coso = cos(node);   sino = sin(node);

    % Satellite position in ECEF (m)
    satPos(tt,1) = xo*coso - yo*cosi*sino;  %x-position  
    satPos(tt,2) = xo*sino + yo*cosi*coso;  %y-position 
    satPos(tt,3) = yo*sini;                 %z-position

    % Satellite clocl bias (m)
    satClkCorr(tt,1) = c*((sat_alm(23)*dt + sat_alm(22))*dt + sat_alm(21)); % meters

    % Keep track of health of each satellite
    %======================================================================      
    health(tt,1) = sat_alm(25); % satellite health (0.00 is useable)


end % END of t_input loop =================================================
%==========================================================================    











    
    