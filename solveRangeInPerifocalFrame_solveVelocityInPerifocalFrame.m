semimajor_axis = input('semimajor_axis = ');
eccentricity = input('eccentricity = ');
true_anomaly = input('true_anomaly = ');

result = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly);
disp('rangeInPQW =');
disp(result);
result = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly);
disp('velocityInPQW =');
disp(result);

function rangeInPQW = solveRangeInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)
    % Convert input angles from degrees to radians
    true_anomaly_rad = deg2rad(true_anomaly);
    
    % Calculate the range (position vector) in the perifocal frame
    r_pqw = semimajor_axis * (1 - eccentricity^2) / (1 + eccentricity * cos(true_anomaly_rad));
    
    % Construct the position vector in the perifocal frame
    rangeInPQW = [r_pqw * cos(true_anomaly_rad);
                  r_pqw * sin(true_anomaly_rad);
                  0];
end


function velocityInPQW = solveVelocityInPerifocalFrame(semimajor_axis, eccentricity, true_anomaly)
    % Convert input angles from degrees to radians
    true_anomaly_rad = deg2rad(true_anomaly);
    
    % Gravitational constant (μ) in km^3/s^2
    mu = 3.986004418*10^5;
    
    % Calculate the range (position vector) in the perifocal frame
    p = semimajor_axis * (1 - eccentricity^2);
    
    % Calculate the magnitude of the velocity in the perifocal frame
    v_pqw = sqrt(mu / p);
    
    % Construct the velocity vector in the perifocal frame
    velocityInPQW = [-v_pqw * sin(true_anomaly_rad);
                     v_pqw * (eccentricity + cos(true_anomaly_rad));
                     0];
end