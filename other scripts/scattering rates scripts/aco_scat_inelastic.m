function [Waco_emi,Waco_abs] = aco_scat(E,T)

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
%% acoustic deformation potential scattering - GaAs (inelastic 1.120)
B=bernoulli(sym(1:10));
z=sym('z','real');
x=sym('x','real');
nphon=1/z-1/2;
for m=1:5
    nphon = nphon + B(2*m)*z^(2*m-1)/factorial(2*m);
end
gamma=@(E) E*(1+alpha*E);
coeff=@(E) (meff^(0.5)*(kB*T)^3*Daco^2)/(2^(2.5)*pi*HBAR^(4)*v_l^(4)*rho*sqrt(gamma(E)));%coefficient definition
E_s=meff*v_l^(2)/2;
C=4*sqrt(E_s)/(kB*T*(1-4*alpha*E_s));
if gamma(E)<(E_s/(1-4*alpha*E_s))
    x1a=C*(sqrt(E_s)*(1+2*alpha*E)-sqrt(gamma(E)));
    x2a=C*(sqrt(E_s)*(1+2*alpha*E)+sqrt(gamma(E)));
    x1e=0;
    x2e=0;
else 
    x1a=0;
    x2a=C*(sqrt(E_s)*(1+2*alpha*E)+sqrt(gamma(E)));
    x1e=0;
    x2e=-C*(sqrt(E_s)*(1+2*alpha*E)-sqrt(gamma(E)));
end

F1=nphon*z^2;
F2=nphon*z^3;
G1=(nphon+1)*z^2;
G2=(nphon+1)*z^3;

cost1=1+2*alpha*E;
cost2=2*alpha*kB*T;
Waco_abs=coeff(E)*(cost1*(int(F1,0,x2a)-int(F1,0,x1a))+cost2*(int(F2,0,x2a)-int(F2,0,x1a)));
Waco_emi=coeff(E)*(cost1*(int(G1,0,x2e)-int(G1,0,x1e))-cost2*(int(G2,0,x2e)-int(G2,0,x1e)));
Waco_abs=double(Waco_abs);
Waco_emi=double(Waco_emi);
end

