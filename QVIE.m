function Q_VIE = QVIE(epsr,rigeig,xy,h,omega_mode,flag)
   
    chi = epsr - 1;
    p_current = (epsr - 1) .* rigeig; %polarization current from radiated field
    p = p_current(flag(:) == 1); %extract just those points within cylinder
    
    % % Greens matrices
    dS = h^2;
    G0EE = dS * cal_G0(xy,omega_mode);
    G0Hx = dS * cal_G0HE_xz(xy, omega_mode);
    G0Hy = dS * cal_G0HE_yz(xy, omega_mode);

    % % incident field
    e_rad = G0EE * p;
    e_tot = p / chi;
    hx_rad = G0Hx * p;
    hy_rad = G0Hy * p;
    
    % QVIE
    U_stored = 1/4 * (epsr * norm(e_rad)^2 + 1 * (norm(hx_rad)^2 + norm(hy_rad)^2)) * dS;
    P_rad = 1/2 * imag(omega_mode'*p'*G0EE*p) * dS;
    Q_VIE = real(omega_mode) * (U_stored / P_rad);

    % display results
    error = norm(e_rad - e_tot) / norm(e_rad) * 100; %calculate percent error
    fprintf('\nError in VIE = %4.2f%%\n', error); 
    disp(['Prad = ', num2str(P_rad)]);
    disp(['Ustore = ', num2str(U_stored)])
    fprintf('Q_VIE = %4.2f\n', Q_VIE)
end

