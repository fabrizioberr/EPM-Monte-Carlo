function k = engy2k(E,iv)
% inverse energy dispersion relation, referred to GAMMA conduction band edge
%k input in m^-1. input energy in J
M0 = 9.1095D-31; % electron mass, kg
HBAR = 1.054D-34; % reduced Planck constant, J*s
meff_G = 0.067*M0;
meff_L = (1.2*0.2^2)^(1/3)*M0;
k =(iv==1) * sqrt(2*meff_G*E/HBAR^2) + ...
(iv==2) * sqrt(2*meff_L*E/HBAR^2);
end

