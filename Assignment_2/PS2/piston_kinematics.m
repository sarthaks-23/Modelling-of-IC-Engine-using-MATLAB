function V = piston_kinematics(bore, stroke, rod, CR, theta)

    V_swept = (3.14/4) * bore^2 * stroke; % Swept volume
    V_clearance = V_swept / (CR - 1); % Clearance volume
    V = zeros(size(theta));
    for i = 1:length(theta)
        r = stroke / 2; 
        R = rod / r; 
        V(i) = V_clearance * (1 + 0.5 * (CR - 1) * (R + 1 - cosd(theta(i)) - sqrt(R^2 - sind(theta(i))^2))); 
    end
end
