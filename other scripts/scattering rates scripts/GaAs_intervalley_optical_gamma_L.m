clear all
close all

% constants
c_light = 2.99792458e+8; % light velocity, m/s
H = 6.626070040e-34; % Planck constant, J*s
HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
T = 300; % temperature, K
eps0 = 8.854D-12; % vacuum permittivity constant, F/m
M0 = 9.1095D-31; % electron mass, kg
rho = 5320; % mass density, kg/mˆ3
eps_s = 12.9*eps0; % static dielectric constant, F/m
eps_infty = 10.9*eps0; % high-frequency dielectric constant, F/m
v_l = 5240; % longitudinal sound velocity, m/s
meff = 0.067*M0; % effective mass, kg
alpha = 0.64/Q; % nonparabolicity factor, 1/J
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J
Dij=1e+11*Q; %coupling constant gamma-L [J/m]
meff_ij=(1.2*0.2^2)^(1/3)*M0;  
deltaE_ij=0.29*Q; % delta energy gamma-L valleys, J
omega_ij=0.0278*Q/HBAR; % da libro HESS pag.108
Nij=1/(exp((HBAR*omega_ij)/(kB*T))-1);
%%
coeff=(pi*Dij^(2)*4)/(rho*omega_ij);
N=@(E) ((2*meff_ij)^(1.5))*sqrt(E)/(4*pi^(2)*HBAR^(3));

nE = 300; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
W_gamma_L_emi=zeros(1,200);
W_gamma_L_abs=zeros(1,200);
for ie = 1:nE, 
    E = vE(ie);
    if (E-HBAR*omega_ij-deltaE_ij>=0)
    W_gamma_L_emi(ie)=coeff*N(E-HBAR*omega_ij-deltaE_ij)*(Nij+1);
    end
    if (E+HBAR*omega_ij-deltaE_ij>=0)
    W_gamma_L_abs(ie)=coeff*N(E+HBAR*omega_ij-deltaE_ij)*(Nij);
    end
end


plot(vE/Q, W_gamma_L_emi,'k','linewidth',2)
hold on
plot(vE/Q,W_gamma_L_abs,'r','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Intervalley (Gamma-L) scattering rate, 1/s'), xlabel('Energy, eV')
legend('emission','absorption')
