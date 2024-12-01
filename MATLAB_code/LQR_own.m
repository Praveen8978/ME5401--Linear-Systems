function [K, V, U, P, M] = LQR_own(A,B,Q,R)
M = [ A   -B/R*B';
     -Q   -A'];

[evec, eval] = eig(M);
eval = sum(eval);
evec_stable = evec(:,find(real(eval)<0));
% P = evec_stable(,:)
V = evec_stable(1:3,:);
U = evec_stable(4:6,:);
P = U/V;

K = real(R \ B' * P);

end