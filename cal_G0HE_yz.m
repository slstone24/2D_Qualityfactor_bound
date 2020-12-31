%G0^HE_yz
function G0HE_yz = cal_G0HE_yz(xy, omega)
    ImG0HE_yz_x1x2 = @(x1,x2) 1/4 .* besselj(1, norm(x1-x2)) ...
        .* (x1(1) - x2(1)) ./ norm(x1-x2); % unit: k^2
    
    Nxy = length(xy);
    G0HE_yz = zeros(Nxy,Nxy);
    for i = 1:Nxy
        x1 = xy(i,:);
        for j = 1:Nxy
           x2 = xy(j,:);
           G0HE_yz(i,j) = omega^2 * ImG0HE_yz_x1x2(omega * x1,omega * x2);
        end
    end
    % set diagonal to zero (as they should be)
    G0HE_yz(1:1+size(G0HE_yz,1):end) = 0; 
end
