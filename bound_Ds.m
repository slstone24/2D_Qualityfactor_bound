function [status, popt] = bound_Ds(S, D)
    % init
    D = D(~cellfun('isempty',D)); % remove empty D entry
    ND = length(D);
    
    % construct matrices for SDP
    C = cell(1,ND);
    for i = 1:ND
        C{i} = Mat_real(D{i}*(S));
    end
    
    % solve SDP by cvx
    n = length(S);
    cvx_begin quiet
        variable X(n,n) hermitian
        minimize 0;
        subject to
            trace(X) == 1;
            for i = 1:ND
                C{i}(:)'*X(:) == 0;
            end
            X == hermitian_semidefinite(n);
    cvx_end
    
    % show feasability
    status = cvx_status;
    
    % extract optimal polarization current
    popt = extract_p_opt(X);
end
