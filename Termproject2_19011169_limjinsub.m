load("nav.mat");
figure;

colors = {'b', 'g', 'r'}; % Colors for GPS, QZSS, and BDS respectively

for n = 1:1:3
    if n == 1
        satellite = nav.GPS;
        sat_name = 'GPS';
    elseif n == 2
        satellite = nav.QZSS;
        sat_name = 'QZSS';
    else 
        satellite = nav.BDS;
        sat_name = 'BDS';
    end

    my_lat = 37;
    my_lon = 127;
    a = satellite.a/10^3;
    e = satellite.e;
    i = satellite.i;
    omega = satellite.omega;
    M0 = satellite.M0;
    OMEGA = satellite.OMEGA;

    mu = 3.986004418*10^5;

    period = 2*pi*sqrt(mu/a^3);

    if M0 < 0
        M0 = M0 + 2*pi;
    end

    t = 0:60:86390; % 0부터 86399까지 1초 간격으로 시간 범위 설정

    r_perifocal = zeros(3, length(t));
    v_perifocal = zeros(3, length(t));
    r_eci = zeros(3, length(t));
    v_eci = zeros(3, length(t));
    r_ecef = zeros(3, length(t));
    v_ecef = zeros(3, length(t));
    r_enu = zeros(3, length(t));
    r_lon = zeros(1, length(t));
    r_lat = zeros(1, length(t));
    r_az = zeros(1, length(t));
    r_el = zeros(1, length(t));

    for p = 1:1:length(t)
        M = M0 + mean_motion(a) * t(p);
        M = rad2deg(M);
        M = mod(M, 360);

        second = rem(t(p), 60);
        minute = floor(t(p)/60);
        minute = mod(minute, 60);
        hour = floor(t(p)/3600);
        toc = satellite.toc;
        toc_G = satellite.toc;
        toc(1,4) = hour;
        toc(1,5) = minute;
        toc(1,6) = second;

        r_perifocal(:,p) = solveRangeInPerifocalFrame(a, e, mean2true(M,e));
        v_perifocal(:,p) = solveVelocityInPerifocalFrame(a, e, mean2true(M,e));
    
        r_eci(:,p) = PQW2ECI(omega, i, OMEGA)*r_perifocal(:,p);
        v_eci(:,p) = PQW2ECI(omega, i, OMEGA)*v_perifocal(:,p);

        r_ecef(:,p) = ECI2ECEF_DCM(toc_G)*r_eci(:,p);
        v_ecef(:,p) = ECI2ECEF_DCM(toc_G)*v_eci(:,p);

        r_my = [6400*cosd(my_lon)*sind(my_lat); 6400*cosd(my_lat)*sind(my_lon); 6400*sind(my_lat)];
        r_enu(:,p) = eci2enu(my_lon, my_lat, toc)*(r_eci(:,p)-r_my);

        r_lon(1,p) = rad2deg(lon(r_ecef(:,p)));
        r_lat(1,p) = rad2deg(lat(r_ecef(:,p)));
        r_az(1,p) = azimuth(r_enu(:,p));
        r_el(1,p) = elevation(r_enu(:,p), my_lat);  %북위 37도
        
    end
   
    subplot(4,2,1), plot3(r_perifocal(1,:), r_perifocal(2,:), r_perifocal(3,:), colors{n});
    hold on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Perifocal Frame');
    grid on;

    subplot(4,2,2), plot3(r_eci(1,:), r_eci(2,:), r_eci(3,:), colors{n});
    hold on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('ECI Frame');
    grid on;

    subplot(4,2,3), plot3(r_ecef(1,:), r_ecef(2,:), r_ecef(3,:), colors{n});
    hold on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('ECEF Frame');
    grid on;

    subplot(4,2,4), plot3(r_enu(1,:), r_enu(2,:), r_enu(3,:), colors{n});
    hold on;
    xlabel('East');
    ylabel('North');
    zlabel('Up');
    title('ENU Frame');
    grid on;

    subplot(4,2,5), geoplot(r_lat, r_lon, colors{n});
    hold on;
    title('Geographic Plot');
    grid on;

    subplot(4,2,6), polarplot(deg2rad(r_az), r_el, colors{1});
    hold on;
    title([sat_name 'GPS']);
    grid on;
    subplot(4,2,7), polarplot(deg2rad(r_az), r_el, colors{2});
    hold on;
    title([sat_name 'QZSS']);
    grid on;
    subplot(4,2,8), polarplot(deg2rad(r_az), r_el, colors{3});
    hold on;
    title([sat_name 'BDS']);
    grid on;
end


%% PQW to ECI
function [rotation_matrix] = PQW2ECI(arg_prg, inc_angle, RAAN) 

arg_prg = deg2rad(arg_prg);
inc_angle = deg2rad(inc_angle);
RAAN = deg2rad(RAAN);
cos_AP = cos(arg_prg);
sin_AP = sin(arg_prg);
cos_INC = cos(inc_angle);
sin_INC = sin(inc_angle);
cos_RAAN = cos(RAAN);
sin_RAAN = sin(RAAN);

% Rotation matrix components
R11 = cos_RAAN * cos_AP - sin_RAAN * sin_AP * cos_INC;
R12 = -cos_RAAN * sin_AP - sin_RAAN * cos_AP * cos_INC;
R13 = sin_RAAN * sin_INC;
R21 = sin_RAAN * cos_AP + cos_RAAN * sin_AP * cos_INC;
R22 = -sin_RAAN * sin_AP + cos_RAAN * cos_AP * cos_INC;
R23 = -cos_RAAN * sin_INC;
R31 = sin_AP * sin_INC;
R32 = cos_AP * sin_INC;
R33 = cos_INC;

% Construct the rotation matrix
rotation_matrix = [R11, R12, R13;
                   R21, R22, R23;
                   R31, R32, R33];
end

%% r_PQW
function [rangeInPQW] = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)  % 각도 unit : deg, 거리 unit : km로 통일

% Convert input angles from degrees to radians
true_anomaly_rad = deg2rad(true_anomaly);

% Calculate the range (position vector) in the perifocal frame
r_pqw = semimajor_axis * (1 - eccentricity^2) ./ (1 + eccentricity * cos(true_anomaly_rad));

% Construct the position vector in the perifocal frame
rangeInPQW = [r_pqw .* cos(true_anomaly_rad);
              r_pqw .* sin(true_anomaly_rad);
              zeros(size(true_anomaly_rad))];
end

%% v_PQW
function [velocityInPQW] = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly) % 각도 unit : deg, 거리 unit : km로 통일

true_anomaly_rad = deg2rad(true_anomaly);

% Gravitational constant (μ) in km^3/s^2
mu = 3.986004418*10^5;

% Calculate the range (position vector) in the perifocal frame
p = semimajor_axis * (1 - eccentricity^2);

% Calculate the magnitude of the velocity in the perifocal frame
v_pqw = sqrt(mu ./ p);

% Construct the velocity vector in the perifocal frame
velocityInPQW = [-v_pqw .* sin(true_anomaly_rad);
                 v_pqw .* (eccentricity + cos(true_anomaly_rad));
                 zeros(size(true_anomaly_rad))];
end

%% eci to ecef
function [DCM] = ECI2ECEF_DCM(time)

standard_time = [2000,1,1,12,0,0];
standard_time = datetime(standard_time);
standard_jd = juliandate(standard_time);

time = datetime(time);
jd = juliandate(time);
GMST = mod(280.4606 + 360.9856473*(jd - standard_jd), 360);
GMST_rad = deg2rad(GMST);

DCM = [cos(GMST_rad), sin(GMST_rad), 0;
       -sin(GMST_rad), cos(GMST_rad), 0;
       0, 0, 1];
end

%% eci to enu
function [rotation_matrix] = eci2enu(my_longitude, my_latitude, toc_G)
    standard_time = [2000,1,1,12,0,0];
    standard_time = datetime(standard_time);
    standard_jd = juliandate(standard_time);

    toc_G = datetime(toc_G);
    jd = juliandate(toc_G);
    GMST = mod(280.4606 + 360.9856473*(jd - standard_jd), 360);
    GMST_rad = deg2rad(GMST);

    rotation_matrix = [-sin(GMST_rad-my_longitude), cos(GMST_rad-my_longitude), 0;
                       -sin(my_latitude)*cos(GMST_rad-my_longitude), -sin(my_latitude)*sin(GMST_rad-my_longitude), cos(my_latitude);
                       cos(my_latitude)*cos(GMST_rad-my_longitude), cos(my_latitude)*sin(GMST_rad-my_longitude), sin(my_latitude)];
end


%% latitude
function [latitude] = lat(ecef)
% Calculate the magnitude of the xy components
magnitude_xy = sqrt(ecef(1,:).^2 + ecef(2,:).^2);

% Check if magnitude_xy is close to zero and set latitude accordingly
latitude = zeros(size(ecef(1,:)));
zero_idx = abs(magnitude_xy) < 1e-8;
latitude(zero_idx) = sign(ecef(3, zero_idx)) * pi/2;
nonzero_idx = ~zero_idx;
latitude(nonzero_idx) = 50*atan2(ecef(3, nonzero_idx), magnitude_xy(nonzero_idx));
end


%% longitude
function [longitude] = lon(ecef)
longitude = atan2(ecef(2,:), ecef(1,:));
end


%% az
function az = azimuth(ENU)
% Calculate the azimuth angle using the atan2 function
az = rad2deg(atan2(ENU(2,:), ENU(1,:)));

% Adjust the azimuth angle to the range [0, 360)
az = mod(az, 360);
end

%% el
function el = elevation(r_enu, el_mask)
el = rad2deg(atan2(r_enu(3,:), sqrt(r_enu(1,:).^2 + r_enu(2,:).^2)));
el(el < el_mask) = NaN;
end

%% mean anomaly to true anomaly
function [true_anomaly] = mean2true(mean_anomaly, eccentricity)
mean_anomaly = deg2rad(mean_anomaly);
eccentricity_anomaly = mean_anomaly;
delta = 1;

while abs(delta) > 10^(-8)
    delta = (eccentricity_anomaly - eccentricity * sin(eccentricity_anomaly) - mean_anomaly) / (1 - eccentricity * cos(eccentricity_anomaly));
    eccentricity_anomaly = eccentricity_anomaly - delta;
end

true_anomaly = atan2(sqrt(1 - eccentricity^2) * sin(eccentricity_anomaly), cos(eccentricity_anomaly) - eccentricity);
true_anomaly = rad2deg(true_anomaly);
end

%% mean motion
function [mean_motion] = mean_motion(semimajor_axis)
% Gravitational constant (μ) in km^3/s^2
mu = 3.986004418*10^5;

% Calculate the mean motion
mean_motion = sqrt(mu / semimajor_axis^3);
end
