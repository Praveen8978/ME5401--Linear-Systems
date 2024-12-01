function [F,K] = decoupling(A,B,C)
    % the relative degree calculated
    sigma = [1 1];

    Bstar=[C(1,:)*A^(sigma(1)-1)*B;
       C(2,:)*A^(sigma(2)-1)*B];

    if det(Bstar) == 0
        error('Bstar matrix is singular')
    end

    F = inv(Bstar);
    
    % calculating C**
    % assuming H(s) = diag(1/(s+1))
    phi = A + eye(3);

    C_star = [C(1,:)*phi;
              C(2,:)*phi];
    K= Bstar\C_star;
end