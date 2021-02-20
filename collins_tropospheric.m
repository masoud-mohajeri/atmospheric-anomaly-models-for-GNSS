function [trop_collins] = collins_tropospheric(reciver_loc,ct,day,Hemisphere)
%% 'n' for Northern hemisphere
%% 's' for Southern Hemisphere
wgs84 = referenceEllipsoid('wgs84');
[lat,long,h_geo] = ecef2geodetic(wgs84,reciver_loc(1),reciver_loc(2),reciver_loc(3)) ;
recever_geodetic(1,1:3) = [lat,long,h_geo] ;
e_reciver = [-sind(recever_geodetic(2)),cosd(recever_geodetic(2)),0];
n_reciver = [-cosd(recever_geodetic(2))*sind(recever_geodetic(1)) ,-sind(recever_geodetic(2))*sind(recever_geodetic(1)),cosd(recever_geodetic(1))];
u_reciver = [cosd(recever_geodetic(2))*cosd(recever_geodetic(1)) , sind(recever_geodetic(2))*cosd(recever_geodetic(1)) , sind(recever_geodetic(1)) ];
phi = recever_geodetic(1) ;
satr_eph = length(ct) ;
for i = 1: satr_eph
    ro(i,1:3) = (ct(i,:) - reciver_loc)./sqrt((reciver_loc(1) - ct(i,1))^2 +...
                                              (reciver_loc(2) - ct(i,2))^2 +...
                                              (reciver_loc(3) - ct(i,3))^2 ) ;
    elevation(i,1) = asin(sum(ro(i,:).*u_reciver));
    azimuth(i,1) = atan2(sum(ro(i,:).*e_reciver),sum(ro(i,:).*n_reciver)) ;
end
[p0,delta_p,t0,delta_t,e0,delta_e,beta0,delta_beta,lambda0,delta_lambda] = collins_table(phi) ;
if Hemisphere == 'n'
    daymin = 28 ;
elseif Hemisphere == 's'
    daymin = 211 ;
end
p = p0 - delta_p*cos(2*pi*(day - daymin )/(365.25)) ;
t = t0 - delta_t*cos(2*pi*(day - daymin )/(365.25)) ;
e = e0 - delta_e*cos(2*pi*(day - daymin )/(365.25)) ;
beta = beta0 - delta_beta*cos(2*pi*(day - daymin )/(365.25)) ;
lambda = lambda0 - delta_lambda*cos(2*pi*(day - daymin )/(365.25)) ;
k1_collins =  77.604 ;
k2_collins =  382000 ;
rd_collins = 287.054 ;
gm_collins = 9.784 ;
g_collins = 9.80665 ;
trz0_dry = (k1_collins * rd_collins * p * 1e-6)/gm_collins ;
trz0_wet = (k2_collins * rd_collins*e*1e-6)/(((lambda+1)*gm_collins -beta*rd_collins)*t) ;
trz_dry = ((1- (beta*h_geo/t))^(g_collins/(rd_collins*beta)) )*trz0_dry;
trz_wet = (((1- (beta*h_geo/t))^((((lambda+1)*g_collins)/(rd_collins*beta))-1)) )*trz0_wet ;
for i = 1:satr_eph
    trop_collins(i,1) = (trz_dry + trz_wet)*(1.001/sqrt(0.002001 +(sin(elevation(i)))^2));
end
end

function [p0,delta_p,t0,delta_t,e0,delta_e,beta0,delta_beta,lambda0,delta_lambda] = collins_table(phi)

z = [ 1013.25  299.65   26.31   6.3e-3   2.77 ;
      1017.25  294.15   21.79   6.05e-3  3.15 ; 
      1015.75  283.15   11.66   5.58e-3  2.57 ;
      1011.75  272.15   6.78    5.39e-3  1.81 ;
      1013.00  263.65   4.11    4.53e-3  1.55 ] ;
del=[ 0.00     0.00     0.00    0.00e-3  0.00 ;
      -3.75    7.00     8.85    0.25e-3  0.33 ;
      -2.25    11.00    7.24    0.32e-3  0.46 ;
      -1.75    15.00    5.36    0.81e-3  0.74 ;
      -0.50    14.50    3.39    0.62e-3  0.30 ] ;


if abs(phi) <= 15
    p0 =           z(1,1) ;
    t0 =           z(1,2) ;
    e0 =           z(1,3);
    beta0 =        z(1,4) ;
    lambda0 =      z(1,5) ;
    delta_p =      del(1,1) ;
    delta_t =      del(1,2) ;
    delta_e =      del(1,3) ;
    delta_beta =   del(1,4) ;
    delta_lambda = del(1,5) ;
elseif 15 < abs(phi) && abs(phi) <= 30
    p0 =           interpolate(15,z(1,1),30,z(2,1),phi) ;
    t0 =           interpolate(15,z(1,2),30,z(2,2),phi) ;
    e0 =           interpolate(15,z(1,3),30,z(2,3),phi) ;
    beta0 =        interpolate(15,z(1,4),30,z(2,4),phi) ;
    lambda0 =      interpolate(15,z(1,5),30,z(2,5),phi) ;
    delta_p =      interpolate(15,del(1,1),30,del(2,1),phi) ;
    delta_t =      interpolate(15,del(1,2),30,del(2,2),phi) ;
    delta_e =      interpolate(15,del(1,3),30,del(2,3),phi) ;
    delta_beta =   interpolate(15,del(1,4),30,del(2,4),phi) ;
    delta_lambda = interpolate(15,del(1,5),30,del(2,5),phi) ;
elseif 30 < abs(phi) && abs(phi) <= 45
    p0 =           interpolate(30,z(2,1),45,z(3,1),phi) ;
    t0 =           interpolate(30,z(2,2),45,z(3,2),phi) ;
    e0 =           interpolate(30,z(2,3),45,z(3,3),phi) ;
    beta0 =        interpolate(30,z(2,4),45,z(3,4),phi) ;
    lambda0 =      interpolate(30,z(2,5),45,z(3,5),phi) ;
    delta_p =      interpolate(30,del(2,1),45,del(3,1),phi) ;
    delta_t =      interpolate(30,del(2,2),45,del(3,2),phi) ;
    delta_e =      interpolate(30,del(2,3),45,del(3,3),phi) ;
    delta_beta =   interpolate(30,del(2,4),45,del(3,4),phi) ;
    delta_lambda = interpolate(30,del(2,5),45,del(3,5),phi) ;
elseif 45 < abs(phi) && abs(phi) <= 60
    p0 =           interpolate(45,z(3,1),60,z(4,1),phi) ;
    t0 =           interpolate(45,z(3,2),60,z(4,2),phi) ;
    e0 =           interpolate(45,z(3,3),60,z(4,3),phi) ;
    beta0 =        interpolate(45,z(3,4),60,z(4,4),phi) ;
    lambda0 =      interpolate(45,z(3,5),60,z(4,5),phi) ;
    delta_p =      interpolate(45,del(3,1),60,del(4,1),phi) ;
    delta_t =      interpolate(45,del(3,2),60,del(4,2),phi) ;
    delta_e =      interpolate(45,del(3,3),60,del(4,3),phi) ;
    delta_beta =   interpolate(45,del(3,4),60,del(4,4),phi) ;
    delta_lambda = interpolate(45,del(3,5),60,del(4,5),phi) ;
elseif 60 < abs(phi) && abs(phi) <= 75
    p0 =           interpolate(60,z(4,1),75,z(5,1),phi) ;
    t0 =           interpolate(60,z(4,2),75,z(5,2),phi) ;
    e0 =           interpolate(60,z(4,3),75,z(5,3),phi) ;
    beta0 =        interpolate(60,z(4,4),75,z(5,4),phi) ;
    lambda0 =      interpolate(60,z(4,5),75,z(5,5),phi) ;
    delta_p =      interpolate(60,del(4,1),75,del(5,1),phi) ;
    delta_t =      interpolate(60,del(4,2),75,del(5,2),phi) ;
    delta_e =      interpolate(60,del(4,3),75,del(5,3),phi) ;
    delta_beta =   interpolate(60,del(4,4),75,del(5,4),phi) ;
    delta_lambda = interpolate(60,del(4,5),75,del(5,5),phi) ;
elseif 75 < abs(phi)
    p0 =           z(5,1) ;
    t0 =           z(5,2) ;
    e0 =           z(5,3);
    beta0 =        z(5,4) ;
    lambda0 =      z(5,5) ;
    delta_p =      del(5,1) ;
    delta_t =      del(5,2) ;
    delta_e =      del(5,3) ;
    delta_beta =   del(5,4) ;
    delta_lambda = del(5,5) ;
end

end

function [yp] = interpolate(x1,y1,x2,y2,xp)
yp = ((y2-y1)/(x2-x1))*(xp-x1)+y1 ;
end

