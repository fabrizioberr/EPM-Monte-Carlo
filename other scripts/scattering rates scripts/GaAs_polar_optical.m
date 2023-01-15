clear all
close all
% constants
H = 6.626070040e-34; % Planck constant, J*s
HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
T = 300; % temperature, K
eps0 = 8.854D-12; % vacuum permittivity constant, F/m
M0 = 9.1095D-31; % electron mass, kg
eps_s = 12.9*eps0; % static dielectric constant, F/m
eps_infty = 10.9*eps0; % high-frequency dielectric constant, F/m
meff = 0.067*M0; % effective mass, kg
alpha = 0.64/Q; % nonparabolicity factor, 1/J
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
wpop=hwpop/HBAR;
egap = 1.424*Q; % energy gap, J
eps_p=1/(1/eps_infty-1/eps_s);
%%
gamma=@(E) E*(1+alpha*E); %for non parabolic
d_gamma=@(E) 1+2*alpha*E; %derivative of gamma(E)

A_e=@(E) (2*(1+alpha*E)*(1+alpha*(E-hwpop))+alpha*(gamma(E)+gamma(E-hwpop)))^2;
A_a=@(E) (2*(1+alpha*E)*(1+alpha*(E+hwpop))+alpha*(gamma(E)+gamma(E+hwpop)))^2;
B_e=@(E) -2*alpha*gamma(E)^(0.5)*gamma(E-hwpop)^(0.5)*(4*(1+alpha*E)*(1+alpha*(E-hwpop))+alpha*(gamma(E)+gamma(E-hwpop)));
B_a=@(E) -2*alpha*gamma(E)^(0.5)*gamma(E+hwpop)^(0.5)*(4*(1+alpha*E)*(1+alpha*(E+hwpop))+alpha*(gamma(E)+gamma(E+hwpop)));
C_e=@(E) 4*(1+2*alpha*E)*(1+2*alpha*(E-hwpop))*(1+alpha*E)*(1+alpha*(E-hwpop));
C_a=@(E) 4*(1+2*alpha*E)*(1+2*alpha*(E+hwpop))*(1+alpha*E)*(1+alpha*(E+hwpop));

k=@(E) sqrt(2*meff*E)/HBAR; %FOR PARABOLIC
qmin_abs=@(E) k(E)*abs(1-sqrt(1+hwpop/E)); 
qmin_emi=@(E) k(E)*abs(1-sqrt(1-hwpop/E));
qmax_abs=@(E) k(E)*abs(1+sqrt(1+hwpop/E));
qmax_emi=@(E) k(E)*abs(1+sqrt(1-hwpop/E));

coeff1=(Q^2*sqrt(meff)*wpop)/(4*pi*eps_p*sqrt(2)*HBAR);
coeff2=(Q^2*wpop)/(8*pi*eps_p);
B=bernoulli(sym(1:10));
z=sym('z','real');
x=sym('x','real');
nphon=1/z-1/2;
for m=1:5
    nphon = nphon + B(2*m)*z^(2*m-1)/factorial(2*m);
end
z=(hwpop)/(kB*T);
Nq=double(subs(nphon)); %well approximates the exact form (line below)
%Nq=1/(exp((hwpop)/(kB*T))-1); 

nE = 200; % number of energy points
vE = linspace(0,1,nE)*Q; % energy axis, J
Wpop_inelastic_abs=zeros(1,nE);
Wpop_inelastic_emi=zeros(1,nE);
Wpop_parabolic_abs=zeros(1,nE);
Wpop_parabolic_emi=zeros(1,nE);
for ie = 1:nE
    E = vE(ie);
    Wpop_inelastic_abs(ie)=(coeff1/sqrt(gamma(E)))*d_gamma(E+hwpop)*(Nq)*(1/C_a(E))*(A_a(E)*log(abs((gamma(E)^(0.5)+gamma(E+hwpop)^(0.5))/(gamma(E)^(0.5)-gamma(E+hwpop)^(0.5))))+B_a(E));
    Wpop_parabolic_abs(ie)=(coeff2*(k(E)/E))*(Nq)*log(qmax_abs(E)/qmin_abs(E));
    if (E>=hwpop)
    Wpop_inelastic_emi(ie)=(coeff1/sqrt(gamma(E)))*d_gamma(E-hwpop)*(Nq+1)*(1/C_e(E))*(A_e(E)*log(abs((gamma(E)^(0.5)+gamma(E-hwpop)^(0.5))/(gamma(E)^(0.5)-gamma(E-hwpop)^(0.5))))+B_e(E));
    Wpop_parabolic_emi(ie)=(coeff2*(k(E)/E))*(Nq+1)*log(qmax_emi(E)/qmin_emi(E));
    end
end

figure(1)
plot(vE/Q, Wpop_inelastic_emi,'b','linewidth',2)
hold on
plot(vE/Q,Wpop_inelastic_abs,'r','linewidth',2)
plot(vE/Q,Wpop_inelastic_abs+Wpop_inelastic_emi,'k','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Polar scattering rate (non parabolic), 1/s'), xlabel('Energy, eV')
legend('emission','absorption','total')

figure(2)
plot(vE/Q, Wpop_parabolic_emi,'b','linewidth',2)
hold on
plot(vE/Q,Wpop_parabolic_abs,'r','linewidth',2)
plot(vE/Q,Wpop_parabolic_abs+ Wpop_parabolic_emi,'k','linewidth',2)
set(gca,'FontSize',14,'FontName','Arial','box','on')
ylabel('Polar scattering rate (parabolic), 1/s'), xlabel('Energy, eV')
legend('emission','absorption','total')