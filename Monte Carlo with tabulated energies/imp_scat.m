function [W_imp] = imp_scat(Ef,T,iv,nI)
%for impurity elastic scattering in montecarlo in NON-PARABOLIC approx.
% iv identifies the zone: 1 -- > Gamma, 2 --> L
%nI impurity concentration in m^-3
HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
M0 = 9.1095D-31; % electron mass, kg
eps0 = 8.854D-12; % vacuum permittivity constant, F/m
eps_s = 12.9*eps0; % static dielectric constant, F/m

q0=sqrt((Q^2*nI)/(kB*T*eps_s));
Z = 1;
if iv==1
    meff = 0.067*M0;
    alpha = 0.64/Q; % 
    N=@(E) ((2*meff)^1.5)*sqrt(E.*(1+alpha*E))/(4*pi^2*HBAR^3).*(1+2*alpha*E);

else
    meff=(1.2*0.2^2)^(1/3)*M0;  %
    alpha = 0.46/Q;
    N=@(E) ((2*meff)^(1.5))*sqrt(E.*(1+alpha*E))/(4*pi^(2)*HBAR^(3)).*(1+2*alpha*E);
end
k=@(E) sqrt(2*meff*E.*(1+alpha*E))/HBAR; 

W_imp=(2*pi*nI*Z^2*Q^4*N(Ef))/((HBAR*eps_s^2)*(q0^2*(4*k(Ef)^2+q0^2)));
end

