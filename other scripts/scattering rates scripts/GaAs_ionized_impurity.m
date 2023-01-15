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
wpop=hwpop/HBAR;
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J
eps_p=1/(1/eps_infty-1/eps_s);
nI=[1e16 1e18 ]*1e6; %density of impurity
Z=1;
q0=sqrt((Q^2*nI)/(kB*T*eps_s));
%%
N=@(E) ((2*meff)^(1.5))*sqrt(E)/(4*pi^(2)*HBAR^(3));
k=@(E) sqrt(2*meff*E)/HBAR; 

nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
for ie = 1:nE 
    E = vE(ie);
    W_impurity_16(ie)=(2*pi*nI(1)*Z^2*Q^4*N(E))/((HBAR*eps_s^(2))*(q0(1)^(2)*(4*k(E)^2+q0(1)^(2))));
    W_impurity_18(ie)=(2*pi*nI(2)*Z^2*Q^4*N(E))/((HBAR*eps_s^(2))*(q0(2)^(2)*(4*k(E)^2+q0(2)^(2))));
   
end

figure(1)
plot(vE/Q, W_impurity_16,'k','linewidth',2)
hold on
plot(vE/Q,W_impurity_18,'r','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('scattering rate, 1/s'), xlabel('Energy, eV')
legend('1x10^{16}','1x10^{18}')