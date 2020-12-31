function [A,beta,c,pftor] = obj_Q(geo,opt,G0,G0Hx,G0Hy,xi_Mat,e_rad,dS)
    switch opt.obj
        case 'Qext'
            A = zeros(size(G0));
            beta = e_rad;
        case 'Qsca'
            A = imag(G0);
            beta = zeros(size(e_rad)); 
        case 'Qabs'
            A = imag(xi_Mat);
            beta = zeros(size(e_rad));
        case 'Qfactor'
            A = (geo.epsr)*(G0)'*G0 + (G0Hx+G0Hy)'*(G0Hx+G0Hy);
            beta = e_rad;
    end
    c = 0;
%     [~, dS] = get_2D_grid(geo);
    switch geo.geo % different ext/.. coefficient for sphere and thin film
        case 'sphere'
            pftor = dS / (2*geo.R); 
        case 'film'
            pftor = dS / geo.Dy;
    end
end
