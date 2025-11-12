%% Function Definition for the Transformed Reduced Stiffness Matrix, [Qbar]: 
function [QbarMatrix] = TransformedRSM(properties, theta)
    E1 = properties.E1 ; 
    E2 = properties.E2 ; 
    v12 = properties.v12 ; 
    G12 = properties.G12 ; 
    v21 = v12*(E2/E1) ;

    Q11 = E1/(1 - v12*v21);
    Q12 = (v12*E2)/(1 - v12*v21);
    Q22 = E2/(1 - v12*v21);
    Q66 = G12;

    m = cosd(theta); 
    n = sind(theta);

    Qbar11 = Q11.*m.^4 + 2.*(Q12 + 2.*Q66).*n.^2.*m.^2 + Q22.*n.^4; 
    Qbar12 = (Q11 + Q22 - 4.*Q66).*n.^2.*m.^2 + Q12.*(n.^4 + m.^4); 
    Qbar16 = (Q11 - Q12 - 2.*Q66).*n.*m.^3 + (Q12 - Q22 + 2.*Q66).*n.^3.*m; 
    Qbar22 = Q11.*n.^4 + 2.*(Q12 + 2.*Q66).*n.^2.*m.^2 + Q22.*m.^4; 
    Qbar26 = (Q11 - Q12 - 2.*Q66).*n.^3.*m + (Q12 - Q22 + 2.*Q66).*n.*m.^3; 
    Qbar66 = (Q11 + Q22 - 2.*Q12 - 2.*Q66).*n.^2.*m.^2 + Q66.*(n.^4 + m.^4);

    QbarMatrix = [Qbar11 Qbar12 Qbar16;
                  Qbar12 Qbar22 Qbar26;
                  Qbar16 Qbar26 Qbar66];
end