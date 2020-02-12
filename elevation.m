% elavation.m	
%
% this GPS utility calculates the elevation from a
% reference position specified in ECEF coordinates (e.g. antenna
% location) to another position specified in ECEF coordinates
% (e.g. satellite location)
%
% input: 'satLoc' matrix which rows contain an SV id number, a GPS
%		time (seconds), and the ECEF coordinates (meters) of the 
%		satellite location
%						[ svID GPStime ECEFx ECEFy ECEFz ;
%						  svID GPStime ECEFx ECEFy ECEFz ;
%											...
%						  svID GPStime ECEFx ECEFy ECEFz ]
%	 		*'obsLoc' vector which contains the ECEF coordinates
%		(meters) of a reference position which elevation and azimuth
%		will be calculated from
%						[ ECEFx ECEFy ECEFz ]
%           *not used since defined in constant.m
%
% output: 'elevation' matrix which rows contain the an SV id number, a
%		GPS time (seconds)and the elevation look angles
%		(degrees) to the satellite
%						[ svID GPStime elevation ;
%						  svID GPStime elevation ;
%											...
%		  				  svID GPStime elevation ]
%
function elevation = elevation(satLoc,obsLoc);
% define constants
% 	constant;
% define satellite locations in ECEF coordinates
	satX = satLoc(:,3);	% meters
	satY = satLoc(:,4);	% meters
	satZ = satLoc(:,5);	% meters
% define observation location in ECEF coordinates
	obsX = obsLoc(1);		% meters
	obsY = obsLoc(2);		% meters
	obsZ = obsLoc(3);		% meters
% compute unit vector from observation station to satellite position
	r = sqrt((satX-obsX).^2+(satY-obsY).^2+(satZ-obsZ).^2); 		
	dx = (satX-obsX)./r;
	dy = (satY-obsY)./r;
	dz = (satZ-obsZ)./r;
% compute the observation latitude and longitude
% wgs84 = wgs84Ellipsoid('meters')
AAA = ecef2geodetic(obsX,obsY,obsZ);

% compute look-angles from observation station to satellite position 
   north = -cos(AAA(2))*sin(AAA(1)).*dx-sin(AAA(2))*sin(AAA(1)).*dy+cos(AAA(1)).*dz; 
   east = -sin(AAA(2)).*dx+cos(AAA(2)).*dy; 
   vertical = cos(AAA(2))*cos(AAA(2)).*dx+sin(AAA(2))*cos(AAA(1)).*dy+sin(AAA(1)).*dz; 
% compute elevation
	elevation = rad2deg((pi/2-acos(vertical)));     % degrees
	return;