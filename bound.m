function q = bound(n,cs,wi,wr,omega_vals,ND)
    
    for j = 1:length(wi)
        for k = 1:ND
            if ~any(strcmp(cs{n,j},'Failed')) || ~any(strcmp(cs{n,j},'Infeasible'))
                fprintf('found bd')
                q = -0.5 * (wr(n) / imag(wi(j)));
            return
            end
            fprintf('not feasible, qbd = qana')
            q = - 0.5 * (wr(n) / imag(omega_vals(n)));
        end
    end

end