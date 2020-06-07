function [W_em,W_ab] = iv_scat(Ef,T,iv)
%for intervalley inelastic optical scattering in montecarlo in NON-PARABOLIC approx.
% iv identifies the zone: 1 -- > Gamma, 2 --> L
%costants from and Tomizawa Hess
%EDIT:
%-LINES 37,40 : isequal(iv,1) since from L to Gamma the reference must
%not be shifted (Ef is already with respect to gamma reference) 

HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
M0 = 9.1095D-31; % electron mass, kg
rho = 5320; % mass density, kg/mˆ3

 % delta energy gamma-L valleys, J
Dij=1e+11*Q; %coupling constant gamma-L [J/m]
omega_ij=0.0278*Q/HBAR; % da libro HESS pag.108
Nij=1/(exp((HBAR*omega_ij)/(kB*T))-1);
W_em = 0;
W_ab = 0;
deltaE_ij=0.29*Q;

if iv==1
    meff_ij=(1.2*0.2^2)^(1/3)*M0;  %final effective mass in L
    alpha = 0.46/Q; % final nonparabolicity factor, 1/J
    Zj = 4;
else
    meff_ij = 0.067*M0;
    alpha = 0.64/Q; % final nonparabolicity factor, 1/J
    Zj = 1;

end
coeff=(pi*Dij^2*Zj)/(rho*omega_ij);

N=@(E) ((2*meff_ij)^1.5)*sqrt(E.*(1+alpha*E))/(4*pi^2*HBAR^3).*(1+2*alpha*E);

    if (Ef-HBAR*omega_ij-deltaE_ij>0 && iv == 1)  %it is different from pop_scat since in pop_scat electron must remains in L valley. Here, it can change valley
        W_em=coeff*N(Ef-HBAR*omega_ij-deltaE_ij)*(Nij+1);
    end
    if (Ef>0 && iv == 2) 
        W_em = coeff*N(Ef-HBAR*omega_ij+deltaE_ij)*(Nij+1);
    end
    if (Ef+HBAR*omega_ij-deltaE_ij>0 && iv == 1)
        W_ab=coeff*N(Ef+HBAR*omega_ij-deltaE_ij)*Nij;
    end
    if (Ef>0 && iv == 2) 
        W_ab = coeff*N(Ef+HBAR*omega_ij+deltaE_ij)*Nij;
    end
end

