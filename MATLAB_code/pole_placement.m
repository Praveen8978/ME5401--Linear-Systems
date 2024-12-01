function K = pole_placement(A,B,P)
    % this function is used to place poles by pole placement 
    % use Full rank method .i.e. controllable canonical form to calculate
    % poles
    
    % required characteristic polynomial
    lambda1 = P(1);
    lambda2 = P(2);
    lambda3 = P(3);

    syms s
    polynomial=(s-lambda1)*(s-lambda2)*(s-lambda3);
    Ad_cof=double(coeffs(polynomial));
    
    % controlability check:
    ctrb = [B A*B A*A*B];
    if rank(ctrb) < 3
        error('it is uncontrollable!');
    end
    
    for i=1:6
        WR(i) = rank(ctrb(:,1:i));
    end
    
    rho = [2 1]; % The controllability Index is determined by inspection of
                 % WR
    WA = [ctrb(:,1)'; ctrb(:,3)'; ctrb(:,2)']';

    M = inv(WA);
    M1 = M(rho(2),:);
    M2 = M(rho(1)+rho(2),:);

    % to find T - transformation matrix
    T = [M1; M1*A; M2];
    T_inv = inv(T);
    
    A_bar = T*A*T_inv;
    B_bar = T*B;
    
    
    % rounding off the low values to zero
    A_bar(abs(A_bar)<10^(-10))=0;
    B_bar(abs(B_bar)<10^(-10))=0;

    A_bar = round(A_bar, 2);
    B_bar = round(B_bar, 2);
    
    % solving for K
    K_bar = sym('k', [2 3]); 
    closed_loop_matrix= A_bar - B_bar*K_bar

    % let the desired closed loop matrix be:
    Ad = [0 1 0;
        0 0 1;
        -Ad_cof(1:3)];

    % solving the equation
    equation = Ad == closed_loop_matrix;

    K_num=solve(equation)
    K_ans=struct2array(K_num);
    K_ans=double(K_ans);
    Kbar=[K_ans(1:3);
         K_ans(4:6)];

    % tranforming K back to original
    K = Kbar*T;
    
end
