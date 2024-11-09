function C = Cecef2enu(lat_deg, lon_deg)

lat = lat_deg*pi/180; lon = lon_deg*pi/180;

% From Misra & Enge p.137
C = [-sin(lon)           cos(lon)            0;
     -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
     cos(lat)*cos(lon)  cos(lat)*sin(lon)  sin(lat)];

%Calt = C1(pi/2-lat)*C3(lon+pi/2)

end 
 
% function C=C1(th)
% sth=sin(th);cth=cos(th);
% C=[1 0 0 ; 0 cth sth ; 0 -sth cth];
% end
% 
% 
% function C=C3(th)
% sth=sin(th);cth=cos(th);
% C=[cth sth 0; -sth cth 0; 0 0 1];
% end