%% 
clear
syms x k q0
f = (1-x)/(2*k^2*(1-x)+q0^2)^2; %integrand of eq.2.2
res=int(f,x);
simplify(res) %result of eq.2.2

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
eps_s = 13.5*eps0; % static dielectric constant, F/m
eps_infty = 10.9*eps0; % high-frequency dielectric constant, F/m
v_l = 5240; % longitudinal sound velocity, m/s
meff = 0.067*M0; % effective mass, kg
alpha = 0.64/Q; % nonparabolicity factor, 1/J
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
Daco = 7*Q; % acoustic deformation potential, J
egap = 1.424*Q; % energy gap, J
%% mu_aco
N=@(E) ((2*meff)^(1.5))*sqrt(E)/(4*pi^(2)*HBAR^(3));
f_m=@(E,T) exp(-E/(kB*T));

nE = 1000; % number of energy points
vE = linspace(1e-4,1,nE)*Q; % energy axis, J
vE = vE';
Tmin = 8;
Tmax = 1000;
T = linspace(Tmin,Tmax,100);

Waco_elastic = 2*pi*Daco^2*kB/(HBAR*v_l^(2)*rho)*N(vE)*T;
tau_aco(:,:) = 1./Waco_elastic(:,:);
mu_aco = zeros(1,length(T));
coeff_mobility = @(T) 2*Q/(meff*3*kB*T);
for i = 1:length(T)
    mu_aco(i) = 2*Q./(meff*3*kB*T(i)).*trapz(vE,f_m(vE,T(i)).*vE.*tau_aco(:,i).*N(vE))...
        ./trapz(vE,f_m(vE,T(i)).*N(vE));
end
figure,
loglog(T,mu_aco,'r--','DisplayName','acoustic elast. d.p.')
hold on
%% mu impurity
nI=1e14*1e6; %density of impurity
Z=1;
k=sqrt(2*meff*vE)/HBAR; 
q0=sqrt((Q^2*nI)./(kB*T*eps_s));

x=linspace(-1,1,100);
integral_symbolic = zeros(length(vE),length(T));
for i = 1:length(vE)
    k_sym=k(i);
    for j = 1:length(T)
        q0_sym = q0(j);
         integral_symbolic(i,j) = -1/(4*k_sym^4)+1/(4*k_sym^4)*...
             log(1+4*k_sym^2/q0_sym^2); %see integral at the beginning of the code
    end
end
coeff_impurity=pi*nI*Z^2*Q^4/(HBAR*eps_s^2)*N(vE);
tau_imp = (coeff_impurity.*integral_symbolic(:,1:end)).^(-1);
mu_imp = zeros(1,length(T));
for i = 1:length(T)
    mu_imp(i) = 2*Q./(meff*3*kB*T(i)).*trapz(vE,f_m(vE,T(i)).*vE.*tau_imp(:,i).*N(vE))...
        ./trapz(vE,f_m(vE,T(i)).*N(vE));
end
loglog(T,mu_imp,'b--','DisplayName','ionized impurity')
hold on
%% Matthiessen’s rule
mu_0 = (mu_imp.^(-1)+mu_aco.^(-1)).^-1;
loglog(T,mu_0,'k','DisplayName','total (except Jacoboni one)');
legend
xlabel('Temperature [K]');
ylabel('mobility [m^2/(s \cdot V)]')
title(['nI = ' num2str(nI) ' m^{-3}'])
%% Jacoboni approximated formula
energy0 = HBAR.^2*q0.^2/(2*meff);
F_bh = log(1+12*kB*T./energy0)-12*kB*T./energy0./(1+12*kB*T./energy0);
mu_jacoboni = 64*eps_s^2/sqrt(meff)*sqrt(pi)./(nI*Q^3*F_bh).*(2*kB*T).^(3/2);
loglog(T,mu_jacoboni,'g--','DisplayName','ionized impurity [Jacoboni, 11.49]')