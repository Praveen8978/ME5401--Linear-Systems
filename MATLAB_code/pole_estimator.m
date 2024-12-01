function [pole1, pole2] = pole_estimator(zeta, omega)
pole1 = -zeta * omega + omega * sqrt(zeta^2 - 1);
pole2 = -zeta * omega - omega * sqrt(zeta^2 - 1);
end