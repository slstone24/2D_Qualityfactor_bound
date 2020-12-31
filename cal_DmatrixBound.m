function cs = cal_DmatrixBound(epsr,ND,xy,w,h)
    % input: epsr, material permittivity
    %        geo, geometric information
    %        opt, other information (polarization, objective, ...)
    %
    % output: fbd, upper bound of the objective, it stores three different
    %         upper bounds, fbd.ImG0 is the upper bound calculated by real
    %         energy conservation, fbd.ReG0 is the upper bound calculated
    %         by both real and reactive energy conservation, fbd.D is the
    %         upper bound calculated by both real and reactive energy
    %         conservation and additional local D-matrix energy conservation
    
    % ------------- Init ----------------
    
    % user input
    dS = h^2;
    xi = get_xi(epsr);
    G0 = dS * cal_G0(xy, w);
    Ngrid = length(G0);
    xi_Mat = xi*eye(Ngrid);
    S = G0 + xi_Mat;
    
    % D matrices
    D = cell(1,ND);
    D{1} = eye(size(S)) * 1i;
    D{2} = eye(size(S)) * 1;
    status = cell(1,ND);
    popt = cell(1,ND);

    %% ImG0 bound

    [status{1},popt{1}] = bound_Ds(S,{D{1}}); %ImG0 bound
    [status{2},popt{2}] = bound_Ds(S,D); %ReG0 bound
  
    %% iterate D matrix bounds (optmize over Ddiff) 

    for i = 3:ND
        tic
        einc = zeros(size(G0,1),1);
        D{i} = get_Dopt(popt{i-1},S,einc,w);
        [status{i}, popt{i}] = bound_Ds(S,D);
        t = toc;
        fprintf(' D matrix: %d / %d: %4.2f seconds, +%5.3f%% \n', i, ND,t)
        fprintf(' status: %s \n',status{i})
        if strcmp('Infeasible',status{i}) == 1 || strcmp('Failed',status{i}) == 1
            break
        else
            continue
        end
    end
    
    cs = status;
end



