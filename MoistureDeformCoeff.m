% Function Definition for the Coefficients of Moisture Deformation
function [betaMatrix] = MoistureDeformCoeff(properties, theta)
    m = cosd(theta);
    n = sind(theta);

    beta1 = properties.b1;
    beta2 = properties.b2;

    betax = beta1.*m.^2 + beta2.*n.^2;
    betay = beta1.*n.^2 + beta2.*m.^2;
    betaxy = 2.*(beta1 - beta2).*m.*n;

    betaMatrix = [betax betay betaxy];
end
