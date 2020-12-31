% In this function, we try to find the mode profiles with certain resonant 
% frequency, omega. Mathematically, we try to solve the eigenprobelm:
% cur*cur*E = omega^2*eps*E

function [rigeig,rigeigval] = ModeSolverTE(h,dim,BC,epsilon,num_mode,omega)
    %% setup and get operator 
    M = round(dim(2)/h)+1;
    N = round(dim(1)/h)+1;
    
    Matx = ones(M, N); %epsx, staggered
    Maty = ones(M, N); %epsy, staggered
    Matz = ones(M, N); %muz, not staggered
    [pmlx, pmly, pmlz] = PML(dim, h, BC);
    [Maxwell] = Eigen_Maxwell_Operator_Construct(dim, h, BC, (Matx .* pmlx).', (Maty .* pmly).', (Matz .* pmlz).');
    B = spdiags(reshape(epsilon',N*M,1),0,N*M,N*M);
    A = sparse(Maxwell);

    Operator = B\A;
    display('finding eigenvectors/eigenvalues')
    [rigeig,rigeigval] = eigs(Operator,num_mode,(omega)^2);
    display('save eigenvectors/eigenvalues')
            
end


