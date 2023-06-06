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
    az = rad2deg(atan2(ENU(:, 2), ENU(:, 1)));

    % Adjust the azimuth angle to the range [0, 360)
    az = mod(az, 360);
end

function el = elevation(ENU, el_mask)
    % Compute the elevation angle in degrees from the given ENU coordinates

    % Calculate the elevation angle using the atan2 function
    el = rad2deg(atan2(ENU(:, 3), sqrt(ENU(:, 1).^2 + ENU(:, 2).^2)));

    % Set elevation angles below the mask to NaN
    el(el < el_mask) = NaN;
end




