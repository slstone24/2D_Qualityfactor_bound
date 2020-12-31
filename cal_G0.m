% G0^EE_zz
function G0 = cal_G0(xy, omega)
    % calculate Green's function
    G0_r = @(r) omega^2 * 1j/4 .* besselh(0, omega * r); 
    
    Nxy = length(xy);
    G0 = zeros(Nxy,Nxy);
    for i = 1:Nxy
        dr = xy(i,:) - xy;
        dr = sqrt(dr(:,1).^2 + dr(:,2).^2);
        G0(i,:) = G0_r(dr);
    end
    G0(1:(1+size(G0,1)):end) = omega^2 * 1j/4; 
end
