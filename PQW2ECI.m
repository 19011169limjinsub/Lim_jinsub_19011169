function R = perifocal_to_eci(arg_prg, inc_angle, RAAN)
    % Perifocal frame to ECI frame transformation
    % Inputs:
    %   arg_prg: Argument of periapsis (rad)
    %   inc_angle: Inclination angle (rad)
    %   RAAN: Right Ascension of the Ascending Node (rad)
    % Outputs:
    %   R: Rotation matrix (3x3)
    
    arg_prg = input('Argument of periapsis (rad): ');
    inc_angle = input('Inclination angle (rad): ');
    RAAN = input('Right Ascension of the Ascending Node (rad): ');

    % Compute necessary trigonometric values
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
    R = [R11, R12, R13;
         R21, R22, R23;
         R31, R32, R33];
    disp('Rotation matrix:');
    disp(R);
end
