function [W_aco] = aco_scat(Ef,T,iv)
%for elastic acoustisc scattering in montecarlo in NON-PARABOLIC approx.
% iv identifies the zone: 1 -- > Gamma, 2 --> L
%costants from appendix tomizawa

HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
M0 = 9.1095D-31; % electron mass, kg
rho = 5320; % mass density, kg/mˆ3
v_l = 5240; % longitudinal sound velocity, m/s
Daco = 7*Q; % acoustic deformation potential, J


     if iv==1
        meff = 0.067*M0; % effective mass, kg
        alpha = 0.64/Q; % nonparabolicity factor, 1/J
        N=@(E) ((2*meff)^(1.5))*sqrt(E.*(1+alpha*E))/(4*pi^(2)*HBAR^(3)).*(1+2*alpha*E);

     else 
        meff = (1.2*0.2^2)^(1/3)*M0; % effective mass, kg
        alpha = 0.46/Q; % nonparabolicity factor, 1/J
        N=@(E) ((2*meff)^(1.5))*sqrt(E.*(1+alpha*E))/(4*pi^(2)*HBAR^(3)).*(1+2*alpha*E);

     end
    W_aco=2*pi*Daco^2*kB*T/(HBAR*v_l^2*rho)*N(Ef);

end

