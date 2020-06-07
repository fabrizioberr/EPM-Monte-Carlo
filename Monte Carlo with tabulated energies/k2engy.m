function energy = k2engy(kx,ky,kz,iv)
% energy dispersion relation, referred to GAMMA conduction band edge
%k input in m^-1
M0 = 9.1095D-31; % electron mass, kg
Q = 1.6e-19; %
HBAR = 1.054D-34; % reduced Planck constant, J*s
if iv == 1
    meff = 0.067*M0;
    alpha = 0.64/Q;
else
     meff = (1.2*0.2^2)^(1/3)*M0;
     alpha = 0.46/Q;
end
gamma =HBAR^2*(kx^2+ky^2+kz^2)/(2*meff);
energy = 2*gamma./(1+sqrt(1+4*alpha*gamma));
end

