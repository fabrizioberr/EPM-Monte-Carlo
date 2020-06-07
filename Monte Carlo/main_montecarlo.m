%main montecarlo
clear
close all
tic
T = 300; % temperature, K
HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
eps0 = 8.854D-12; % vacuum permittivity constant, F/m
eps_s = 12.9*eps0; % static dielectric constant, F/m
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
dEgammaL = 0.29*Q;

tau_sim = 1e-10; % simulation time, s
tau = 1e-16; % flight time, s
kx = 0; ky = 0; kz = 0; % initial state (step 1)
F = 0.1e5; % electric field along the x-axis, V/m
iv = 1; % valley index (1 -> G, 2 ->L)

nI = 1e20; %impurity concentration in m^-3
time = 0:tau:tau_sim; % time axis, s
nt = length(time);
DeltaE = zeros(1,nt);
mean_v_tau = zeros(1,nt);
E_tau = zeros(1,nt);
q0=sqrt((Q^2*nI)/(kB*T*eps_s));
%for polar optical scattering
B=bernoulli(sym(1:10));
z=sym('z','real');
x=sym('x','real');
nphon=1/z-1/2;
for m=1:5
    nphon = nphon + B(2*m)*z^(2*m-1)/factorial(2*m);
end
z=(hwpop)/(kB*T);
Nq=double(subs(nphon)); %well approximates the exact form

R = @(kx,ky,kz) [ky/sqrt(kx^2+ky^2), ...
kx*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)), ....
kx/sqrt(kx^2+ky^2+kz^2); ...
-kx/sqrt(kx^2+ky^2), ...
ky*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)), ...
ky/sqrt(kx^2+ky^2+kz^2); ...
0, ...
-sqrt(kx^2+ky^2)/sqrt(kx^2+ky^2+kz^2), ...
kz/sqrt(kx^2+ky^2+kz^2)];

for ii = 1:nt
t = time(ii);
Ei = k2engy(kx,ky,kz,iv);
% move particle in k-space (step 3)
kx = kx -Q*F*tau/HBAR;
% collect statistics (step 4)
Ef = k2engy(kx,ky,kz,iv);
DeltaE(ii) = Ef - Ei; % energy variation during drift
E_tau(ii) = (Ei + Ef)/2; % average energy during drift
mean_v_tau(ii) = -DeltaE(ii)/(Q*F*tau);
% choose scattering mechanim and select final state (step 5,6)
[W(1),W(2)] = pop_scat(Ef,T,iv,Nq); %emission and absorption, respectively NAN FOR ABSORPTION
[W(3)] = aco_scat(Ef,T,iv);
[W(4)] = imp_scat(Ef,T,iv,nI);
[W(5),W(6)] = iv_scat(Ef,T,iv);
L = cumsum(W);
r = rand; % pick a random number
norm_ki = norm([kx,ky,kz]); %initial k vector of electron AFTER drift, BEFORE scattering
phi = 2*pi*rand;
r_theta = rand;

flag_scat = 1;
if(r < L(1)*tau)  % select polar optical emission
    f = 2*sqrt(Ef*(Ef-0.0354*Q))/(sqrt(Ef)-sqrt(Ef-0.0354*Q)).^2; %MODIFIFIED
    cos_theta = ((1+f)-(1+2*f)^rand)/f;%MODIFIED 
    norm_kf = engy2k(Ef-0.0354*Q,iv); %energy optical phonon subtracted
        
elseif(r < L(2)*tau)  % select polar optical absorption
    f = 2*sqrt(Ef*(Ef+0.0354*Q))/(sqrt(Ef)-sqrt(Ef+0.0354*Q)).^2;%MODIFIED
    cos_theta = ((1+f)-(1+2*f)^rand)/f; %MODIFIED
    norm_kf = engy2k(Ef+0.0354*Q,iv); %energy optical phonon
elseif(r < L(3)*tau)  % select acoustic scattering
     cos_theta = 1-2*rand;
     norm_kf = norm_ki;
elseif (r < L(4)*tau) %impurity scattering
    r_theta=rand;
    cos_theta = 1-2*r_theta/(1+(1-r_theta)*(2*norm_ki/q0)^2);
    norm_kf = norm_ki;
elseif (r < L(5)*tau) % intervalley emission optical
    cos_theta = 1-2*rand;
    iv = 1.*(iv==2) + 2.*(iv ==1);
            norm_kf = engy2k(Ef-0.0278*Q+dEgammaL*isequal(iv,1)-dEgammaL*isequal(iv,2),iv); %energy optical phonon
elseif (r < L(6)*tau) % intervalley absorption optical
      cos_theta = 1-2*rand;
      iv = 1.*(iv==2) + 2.*(iv ==1);
     norm_kf = engy2k(Ef+0.0278*Q+dEgammaL*isequal(iv,1)-dEgammaL*isequal(iv,2),iv); %energy optical phonon
else % no scattering has occured
    flag_scat = 0;
end 
if flag_scat
    theta = acos(cos_theta);
    k_out = R(kx,ky,kz)*[norm_kf*sin(theta)*cos(phi); norm_kf*sin(theta)*sin(phi); norm_kf*cos(theta)];
    kx = k_out(1);
    ky = k_out(2);
    kz = k_out(3);
end

end % end of for loop
toc
%% post-processing: convergence, diffusivity, ..
mean_E_t = cumsum(E_tau)*tau./time;
figure(1),
plot(time*1e12,mean_E_t/Q*1000,'DisplayName','Average cumulative energy');
ylabel('Electron energy [meV]');
xlabel('time [ps]');
legend

mean_v_t = -cumsum(DeltaE)./(Q*F*time);
figure(2),
plot(time*1e12,mean_v_t*100,'DisplayName','Average cumulative velocity');
ylabel('velocity [cm/s]');
xlabel('time [ps]');
legend
%% autocorrelation function
figure(3),
delta_v= mean_v_tau-mean_v_t;
delta_v=delta_v(2:end); %since first element is NaN by dividing 0
MAXLAG = 20000;
[Cv,lags] = xcorr(delta_v,delta_v,MAXLAG);
plot(lags(MAXLAG+1:end)*tau*1e12,Cv(MAXLAG+1:end)/length(delta_v),'ro-');
xlabel('time, ps')
ylabel('Autocorrelation function, m^2/s^2')

%diffusivity
vector_D = trapz(lags(MAXLAG+1:end)*tau,Cv(MAXLAG+1:end)/length(delta_v))