function [ion_delay] = klubochar_ion(alpha,beta,reciver_loc,toe,ct)
wgs84 = referenceEllipsoid('wgs84');
[lat,long,h_geo] = ecef2geodetic(wgs84,reciver_loc(1),reciver_loc(2),reciver_loc(3)) ;
recever_geodetic(1,1:3) = [lat,long,h_geo] ;
e_reciver = [-sind(recever_geodetic(2)),cosd(recever_geodetic(2)),0];
n_reciver = [-cosd(recever_geodetic(2))*sind(recever_geodetic(1)) ,-sind(recever_geodetic(2))*sind(recever_geodetic(1)),cosd(recever_geodetic(1))];
u_reciver = [cosd(recever_geodetic(2))*cosd(recever_geodetic(1)) , sind(recever_geodetic(2))*cosd(recever_geodetic(1)) , sind(recever_geodetic(1)) ];

phi_u = recever_geodetic(1)*pi/180 ;
lambda_u = recever_geodetic(2)*pi/180 ;
h_gps = 350e3;
phi_p = 78.3 *pi/180 ;
lambda_p = 291 * pi/180 ;
satr_eph = length(toe) ;
for i = 1: satr_eph
    ro(i,1:3) = (ct(i,:) - reciver_loc)/norm(reciver_loc - ct(i,:)) ; 
    elevation(i,1) = asin(sum(ro(i,:).*u_reciver));
    azimuth(i,1) = atan2(sum(ro(i,:).*e_reciver),sum(ro(i,:).*n_reciver)) ;
end
si = 0.5*pi - elevation - asin( (wgs84.SemimajorAxis/(wgs84.SemimajorAxis+h_gps))*cos(elevation)) ;
phi_ipp = asin(sin(phi_u)*cos(si)+cos(phi_u)*sin(si).*cos(azimuth)) ;
lambda_ipp = lambda_u + (si.*sin(azimuth))./(cos(phi_ipp)) ;
phi_megnatic = asin(sin(phi_ipp).*sin(phi_p)+cos(phi_ipp).*cos(phi_p).*cos(lambda_ipp-lambda_p));
t_ipp = 43200*lambda_ipp/pi + toe;
for i = 1:length(t_ipp)
    while t_ipp(i) > 86400
        t_ipp(i) = t_ipp(i) - 86400 ;
    end
    while t_ipp(i) < 0
        t_ipp(i) = t_ipp(i) + 86400 ;
    end
end
for i = 1 : length(phi_megnatic)
    amplitude(i,1) = alpha(1)*(phi_megnatic(i)/pi)^0 +...
                     alpha(2)*(phi_megnatic(i)/pi)^1+...
                     alpha(3)*(phi_megnatic(i)/pi)^2+...
                     alpha(4)*(phi_megnatic(i)/pi)^3 ;
    if amplitude(i,1)<0
        amplitude(i,1) = 0 ;
    end
    period(i,1) = beta(1)*(phi_megnatic(i)/pi)^0 +...
                  beta(2)*(phi_megnatic(i)/pi)^1+...
                  beta(3)*(phi_megnatic(i)/pi)^2+...
                  beta(4)*(phi_megnatic(i)/pi)^3 ;
    if period(i,1)<72000
        period(i,1) = 72000 ;
    end
end
phase_ion = 2*pi*(t_ipp-50400)./period;
slant_factor = (1-(wgs84.SemimajorAxis.*cos(elevation)./(wgs84.SemimajorAxis+h_gps)  ).^2).^-0.5 ;
for i =1:length(phase_ion)
    if abs(phase_ion(i))<0.5*pi
        ion_delay(i,1) = (5e-9 +amplitude(i)*cos(phase_ion(i)))*slant_factor(i) ; 
    elseif 0.5*pi<=abs(phase_ion(i))
        ion_delay(i,1) = 5e-9 * slant_factor(i) ;
    end
end
end

