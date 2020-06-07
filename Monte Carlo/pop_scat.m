function [W_em,W_ab] = pop_scat(Ef,T,iv,Nq)
%for polar optical scattering in montecarlo in NON-PARABOLIC approx.
%%%%%ACTHUNG: COPY IN MAIN THE FOLLOWING
% B=bernoulli(sym(1:10));
% z=sym('z','real');
% x=sym('x','real');
% nphon=1/z-1/2;
% for m=1:5
%     nphon = nphon + B(2*m)*z^(2*m-1)/factorial(2*m);
% end
% z=(hwpop)/(kB*T);
% Nq=double(subs(nphon)); %well approximates the exact form (line below)

HBAR = 1.054D-34; % reduced Planck constant, J*s
Q = 1.6021766208e-19; % elementary charge, C

M0 = 9.1095D-31; % electron mass, kg
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
wpop=hwpop/HBAR;
eps0 = 8.854D-12; % vacuum permittivity constant, F/m
eps_s = 12.9*eps0; % static dielectric constant, F/m
eps_infty = 10.9*eps0; % high-frequency dielectric constant, F/m
eps_p=1/(1/eps_infty-1/eps_s);
if iv ==1
    alpha = 0.64/Q;
    meff = 0.067*M0;
    gamma=@(E) E*(1+alpha*E);
else
    alpha = 0.46/Q;
    meff=(1.2*0.2^2)^(1/3)*M0;
    
    gamma=@(E) (E)*(1+alpha*E);

end
d_gamma=@(E) 1+2*alpha*E;

A_e=@(E) (2*(1+alpha*E)*(1+alpha*(E-hwpop))+alpha*(gamma(E)+gamma(E-hwpop)))^2;
A_a=@(E) (2*(1+alpha*E)*(1+alpha*(E+hwpop))+alpha*(gamma(E)+gamma(E+hwpop)))^2;
B_e=@(E) -2*alpha*gamma(E)^(0.5)*gamma(E-hwpop)^(0.5)*(4*(1+alpha*E)*(1+alpha*(E-hwpop))+alpha*(gamma(E)+gamma(E-hwpop)));
B_a=@(E) -2*alpha*gamma(E)^(0.5)*gamma(E+hwpop)^(0.5)*(4*(1+alpha*E)*(1+alpha*(E+hwpop))+alpha*(gamma(E)+gamma(E+hwpop)));
C_e=@(E) 4*(1+2*alpha*E)*(1+2*alpha*(E-hwpop))*(1+alpha*E)*(1+alpha*(E-hwpop));
C_a=@(E) 4*(1+2*alpha*E)*(1+2*alpha*(E+hwpop))*(1+alpha*E)*(1+alpha*(E+hwpop));

coeff1=(Q^2*sqrt(meff)*wpop)/(4*pi*eps_p*sqrt(2)*HBAR);
if Ef>hwpop
    W_em=(coeff1/sqrt(gamma(Ef)))*d_gamma(Ef-hwpop)*(Nq+1)*(1/C_e(Ef))*(A_e(Ef)...
        *log(abs((gamma(Ef)^(0.5)+gamma(Ef-hwpop)^(0.5))/(gamma(Ef)^(0.5)-gamma(Ef-hwpop)^(0.5))))+B_e(Ef));
else
    W_em = 0;
end
if Ef>0
    W_ab=(coeff1/sqrt(gamma(Ef)))*d_gamma(Ef+hwpop)*(Nq)*(1/C_a(Ef))*(A_a(Ef)...
        *log(abs((gamma(Ef)^(0.5)+gamma(Ef+hwpop)^(0.5))/(gamma(Ef)^(0.5)-gamma(Ef+hwpop)^(0.5))))+B_a(Ef));
else
    W_ab = 0;
end
end

