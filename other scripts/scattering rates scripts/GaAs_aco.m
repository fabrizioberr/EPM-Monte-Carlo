clear all
close all

% constants
c_light = 2.99792458e+8; % light velocity, m/s
H = 6.626070040e-34; % Planck constant, J*s
HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
%T = 300; % temperature, K
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
%%
nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
Waco_abs_300 = zeros(1,nE);
Waco_emi_300 = zeros(1,nE);
Waco_abs_77 = zeros(1,nE);
Waco_emi_77 = zeros(1,nE);
T=[300 77];
for ie = 1:nE, 
    E = vE(ie);
    [Waco_emi_300(ie), Waco_abs_300(ie)] = aco_scat_inelastic(E,T(1));
    [Waco_emi_77(ie), Waco_abs_77(ie)] = aco_scat_inelastic(E,T(2));
end

figure(1), hold on
semilogy(vE/Q,Waco_emi_300,'r.-','linewidth',2)
semilogy(vE/Q,Waco_abs_300,'b.-','linewidth',2)
semilogy(vE/Q,Waco_emi_300+Waco_abs_300,'k.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Acoustic scattering rate (300 K), 1/s'), xlabel('Energy, eV')
legend('emission','absorption','total')

figure(2), hold on
semilogy(vE/Q,Waco_emi_77,'r.-','linewidth',2)
semilogy(vE/Q,Waco_abs_77,'b.-','linewidth',2)
semilogy(vE/Q,Waco_emi_77+Waco_abs_77,'k.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Acoustic scattering rate, 1/s (77 K)'), xlabel('Energy, eV')
legend('emission','absorption','total')

figure(3), hold on
semilogy(vE/Q,Waco_emi_300+Waco_abs_300,'r.-','linewidth',2)
semilogy(vE/Q,Waco_emi_77+Waco_abs_77,'b.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Acoustic scattering rate, 1/s'), xlabel('Energy, eV')
legend('300 K','77 K')

%% elastic approximation (emission rate = absorption rate)
N=@(E) ((2*meff)^(1.5))*sqrt(E)/(4*pi^(2)*HBAR^(3));

nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
for ie = 1:nE, 
    E = vE(ie);
    Waco_elastic_300(ie)=(2*pi*Daco^(2)*kB*T(1))/(HBAR*v_l^(2)*rho)*N(E);
    Waco_elastic_77(ie)=(2*pi*Daco^(2)*kB*T(2))/(HBAR*v_l^(2)*rho)*N(E);
end

figure(4), hold on
semilogy(vE/Q,Waco_elastic_300,'k.-','linewidth',2)
semilogy(vE/Q,Waco_emi_300+Waco_abs_300,'r.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Acoustic scattering rate, 1/s'), xlabel('Energy, eV')
legend('elastic','inelastic')
title('300 K')

figure(5), hold on
semilogy(vE/Q,Waco_elastic_77,'k.-','linewidth',2)
semilogy(vE/Q,Waco_emi_77+Waco_abs_77,'r.-','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Acoustic elastic scattering rate, 1/s'), xlabel('Energy, eV')
legend('elastic','inelastic')
title('77 K')


