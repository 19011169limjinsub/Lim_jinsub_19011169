YYYY = input('YYYY: ');
MM = input('MM: ');
DD = input('DD: ');
hh = input('hh: ');
mm = input('mm: ');
ss = input('ss: ');
time = datetime(YYYY, MM, DD, hh, mm, ss);

n = input('Number of satellites: ');

ENU = zeros(n, 3);
for i = 1:n
    fprintf('Satellite %d:\n', i);
    ENU(i, 1) = input('E (km): ');
    ENU(i, 2) = input('N (km): ');
    ENU(i, 3) = input('U (km): ');
end

el_mask = input('el_mask (deg): ');

result = ECI2ECEF_DCM(time);
disp('ECI2ECEF_DCM =');
disp(result);

az = zeros(1, n);
el = zeros(1, n);

for i = 1:n
    result = azimuth(ENU(i, :));
    az(i) = result;
end

disp('azimuth angle =');
disp(az);

for i = 1:n
    result = elevation(ENU(i, :), el_mask);
    el(i) = result;
end

disp('elevation angle =');
disp(el);

function DCM = ECI2ECEF_DCM(time)
    % Convert the given time to Julian date
    JD = juliandate(datetime(time));

    % Compute the Greenwich Mean Sidereal Time (GMST) in degrees
    GMST = mod(280.4606 + 360.9856473*(JD - 2451545), 360);

    % Convert GMST to radians
    GMST_rad = deg2rad(GMST);

    % Calculate the Direction Cosine Matrix (DCM) for ECI to ECEF transformation
    DCM = [cos(GMST_rad), sin(GMST_rad), 0;
           -sin(GMST_rad), cos(GMST_rad), 0;
           0, 0, 1];
end

function az = azimuth(ENU)
    % Compute the azimuth angle in degrees from the given ENU coordinates

    % Calculate the azimuth angle using the atan2 function
    az = rad2deg(atan2(ENU(2), ENU(1)));

    % Adjust the azimuth angle to the range [0, 360)
    az = mod(az, 360);
end

function el = elevation(ENU, el_mask)
    % Compute the elevation angle in degrees from the given ENU coordinates

    % Calculate the elevation angle using the atan2 function
    el = rad2deg(atan2(ENU(3), sqrt(ENU(1).^2 + ENU(2).^2)));

    % Set elevation angles below the mask to NaN
    el(el < el_mask) = NaN;
end
