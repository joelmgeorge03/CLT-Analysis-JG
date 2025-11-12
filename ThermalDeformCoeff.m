% Function Definition for the Coefficients of Thermal Deformation
function [alphaMatrix] = ThermalDeformCoeff(properties, theta)
    m = cosd(theta);
    n = sind(theta);
    
    alpha1 = properties.a1;
    alpha2 = properties.a2;

    alphax = alpha1.*m.^2 + alpha2.*n.^2;
    alphay = alpha1.*n.^2 + alpha2.*m.^2;
    alphaxy = 2.*(alpha1 - alpha2).*m.*n;

    alphaMatrix = [alphax alphay alphaxy];
end