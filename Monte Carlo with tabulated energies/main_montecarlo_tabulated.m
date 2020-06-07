%main montecarlo
clear
close all
T = 300; % temperature, K
HBAR = 1.054D-34; % reduced Planck constant, J*s
kB = 1.3806488e-23; % Boltzmann constant, J/K
Q = 1.6021766208e-19; % elementary charge, C
eps0 = 8.854D-12; % vacuum permittivity constant, F/m
eps_s = 12.9*eps0; % static dielectric constant, F/m
hwpop = 0.0354*Q; % longitudinal optical phonon energy, J
homega_ij = 0.0278*Q; % longitudinal optical scattering for intervalley
dEgammaL = 0.29*Q;

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
%rotation matrix
R = @(kx,ky,kz) [ky/sqrt(kx^2+ky^2), ...
kx*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)), ....
kx/sqrt(kx^2+ky^2+kz^2); ...
-kx/sqrt(kx^2+ky^2), ...
ky*kz/(sqrt(kx^2+ky^2+kz^2)*sqrt(kx^2+ky^2)), ...
ky/sqrt(kx^2+ky^2+kz^2); ...
0, ...
-sqrt(kx^2+ky^2)/sqrt(kx^2+ky^2+kz^2), ...
kz/sqrt(kx^2+ky^2+kz^2)];

% preparing scattering rates
nI = 1e20; %impurity concentration in m^-3
q0=sqrt((Q^2*nI)/(kB*T*eps_s));
N = 1e3; 
E_max = 1.75; % in eV
vE = Q*linspace(0,E_max,N);
Wmatrix = zeros(7,length(vE),2); %1st dim. different scatt.rates, 2nd dim. energy value and 3rd dim the
%for impact ionization coefficients
Eth = 1.7*Q; %energy J, 
S = 2.5e15; % s^-1
Egap = 1.441*Q; %for GaAs, @300 K
idx_gap = find((vE-Egap).*(vE-Egap>0),1); %finds the idx in vE corresponding to energy greater than Egap
idx_Eth = find((vE-Eth).*(vE-Eth>0),1)-1; %finds the idx in vE corresponding to largest energy SMALLER than Eth
for i = idx_gap:idx_Eth
    Wmatrix(7,i,1) = S*((vE(i)-Eth)/Eth)^2; 
    Wmatrix(7,i,2) = Wmatrix(7,i,1); %equal scattering rates because it depends on the kinetic energy.
end

% GAMMA valley
iv = 1;
for i =1:length(vE)
    [Wmatrix(1,i,iv),Wmatrix(2,i,iv)] = pop_scat(vE(i),T,iv,Nq); %emission and absorption, respectively
    [Wmatrix(3,i,iv)] = aco_scat(vE(i),T,iv);
    [Wmatrix(4,i,iv)] = imp_scat(vE(i),T,iv,nI);
    [Wmatrix(5,i,iv),Wmatrix(6,i,iv)] = iv_scat(vE(i),T,iv);
end
%L valley
iv = 2;
for i =1:length(vE)
    [Wmatrix(1,i,iv),Wmatrix(2,i,iv)] = pop_scat(vE(i),T,iv,Nq); %emission and absorption, respectively
    [Wmatrix(3,i,iv)] = aco_scat(vE(i),T,iv);
    [Wmatrix(4,i,iv)] = imp_scat(vE(i),T,iv,nI);
    [Wmatrix(5,i,iv),Wmatrix(6,i,iv)] = iv_scat(vE(i),T,iv);
end
Lmatrix = cumsum(Wmatrix,1);
%% MONTECARLO
close all
tic
nF = 1;
vF = linspace(25,66,nF)*1e6;

tau_sim = 1e-8; % simulation time, s
tau = 1e-16; % flight time, s
time = 0:tau:tau_sim; % time axis, s
nt = length(time);

mean_E_final = zeros(1,nF);
mean_v_final = zeros(1,nF);
occup_gamma = zeros(1,nF); %counts the fraction of time the electron is gamma for each EL.FIELD cycle
vector_std_E = zeros(1,nF);
vector_std_v = zeros(1,nF);
vector_D=zeros(1,nF); %diffusivity

v_mean_alpha = zeros(1,nF); %impact ioniz.
v_std_alpha = zeros(1,nF); %impact ioniz.
for jj = 1:nF
    F = vF(jj); % electric field along the x-axis, V/m
    DeltaE = zeros(1,nt);
    mean_v_tau = zeros(1,nt);
    E_tau = zeros(1,nt);

    kx = 0; ky = 0; kz = 0; % initial state (step 1)
    iv = 1;
    flag_gamma = 0; %counts how many times the electron is gamma for each TIME cycle
    tot_path = 0;
    count_impact = 0;
    mean_free_path = zeros(1,nt);
    for ii = 1:nt
        flag_gamma = flag_gamma + isequal(iv,1);
        t = time(ii);
        Ei = k2engy(kx,ky,kz,iv);
        % move particle in k-space (step 3)
        kx = kx -Q*F*tau/HBAR;
        % collect statistics (step 4)
        Ef = k2engy(kx,ky,kz,iv);
        DeltaE(ii) = Ef - Ei; % energy variation during drift
        E_tau(ii) = (Ei + Ef)/2; % average energy during drift
        mean_v_tau(ii) = - DeltaE(ii)/(Q*F*tau); % diffusion: we use this also for istantaneous velocity
        tot_path = tot_path + abs(mean_v_tau(ii)*tau);
        % choose scattering mechanim and select final state (step 5,6)
        [value, idx] = min(abs(vE-Ef));
        Ef = vE(idx);
        r = rand; % pick a random number
        norm_ki = norm([kx,ky,kz]); %initial k vector of electron AFTER drift, BEFORE scattering
        phi = 2*pi*rand;
        r_theta = rand;

        flag_scat = 1;
        if(r < Lmatrix(1,idx(1),iv)*tau)  % select polar optical emission
            f = 2*sqrt(Ef*(Ef-0.0354*Q))/(sqrt(Ef)-sqrt(Ef-0.0354*Q)).^2; 
            cos_theta = ((1+f)-(1+2*f)^r_theta)/f; 
            norm_kf = engy2k(Ef-0.0354*Q,iv); %energy optical phonon subtracted

        elseif(r < Lmatrix(2,idx(1),iv)*tau)  % select polar optical absorption
            f = 2*sqrt(Ef*(Ef+0.0354*Q))/(sqrt(Ef)-sqrt(Ef+0.0354*Q)).^2;
            cos_theta = ((1+f)-(1+2*f)^r_theta)/f; 
            norm_kf = engy2k(Ef+0.0354*Q,iv); %energy optical phonon

        elseif(r < Lmatrix(3,idx(1),iv)*tau)  % select acoustic scattering
             cos_theta = 1-2*r_theta;
             norm_kf = norm_ki;

        elseif (r < Lmatrix(4,idx(1),iv)*tau) %impurity scattering
            cos_theta = 1-2*r_theta/(1+(1-r_theta)*(2*norm_ki/q0)^2);
            norm_kf = norm_ki;

        elseif (r < Lmatrix(5,idx(1),iv)*tau) % intervalley emission optical
            cos_theta = 1-2*r_theta;
            iv = 1.*(iv==2) + 2.*(iv ==1);
            norm_kf = engy2k(Ef-0.0278*Q+dEgammaL*isequal(iv,1)-dEgammaL*isequal(iv,2),iv); %energy optical phonon
        elseif (r < Lmatrix(6,idx(1),iv)*tau) % intervalley absorption optical
              cos_theta = 1-2*r_theta;
              iv = 1.*(iv==2) + 2.*(iv ==1);
              norm_kf = engy2k(Ef+0.0278*Q+dEgammaL*isequal(iv,1)-dEgammaL*isequal(iv,2),iv); %energy optical phonon
        elseif (r < Lmatrix(7,idx(1),iv)*tau) %impact ionization
              cos_theta = 1-2*r_theta; %being isotropic, [30,pag.5]
              iv = randi([1 2],1); %see [30,pag.5];
              norm_kf = engy2k(Ef-Egap,iv);
              if count_impact == 0
                  mean_free_path(1) = tot_path;
                  count_impact = 1;
              else
                  count_impact = count_impact+1;
                  mean_free_path(count_impact) = tot_path-sum(mean_free_path(1:count_impact-1));
              end
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
    %average energy
    mean_E_t = cumsum(E_tau)*tau./time;
    mean_E_final(jj) = mean_E_t(end);
    vector_std_E(jj)= std(mean_E_t(floor(nt/10):end)); %we pick last 90% of simulation to evaluate std
   
    %average velocity
    mean_v_t = -cumsum(DeltaE)./(Q*F*time); % diffusion: we use it for avg time up to simulation time step
    mean_v_final(jj) = mean_v_t(end); 
    vector_std_v(jj) = std(mean_v_t(floor(nt/10):end)); %we pick last 90% of simulation to evaluate std
    occup_gamma(jj) = flag_gamma/nt;
%     %impact ionization
%     mean_free_path = mean_free_path((mean_free_path>0)); % to remove the 0 tail
%     v_alpha = 1./mean_free_path; %only last 90% elements
%     v_mean_alpha(jj) = mean(v_alpha);
%    v_std_alpha(jj) = std(v_alpha); 
    mean_free_path = mean_free_path((mean_free_path>0)); % to remove the 0 tail
    v_mean_alpha(jj) = 1/mean(mean_free_path);
    v_std_alpha(jj) = std(1./mean_free_path);
%     
    disp(['Progress ' num2str(jj/nF*100) '%'])
end
toc
%% plotting mean energy vs electric field, velocity vs electric field, occupancies in valleys
figure(1),
%plot(vF/1e5,mean_E_final/Q*1000,'DisplayName','Average energy');
errorbar(vF/1e5,mean_E_final/Q*1000,vector_std_E/Q*1000,'DisplayName','Average energy');
ylabel('Electron energy [meV]');
xlabel('Electric field [kV/cm]');
legend

figure(2),
%plot(vF/1e5,abs(mean_v_final*100),'DisplayName','Average velocity');
errorbar(vF/1e5,abs(mean_v_final*100),abs(vector_std_v*100),'DisplayName','Average velocity');
ylabel('velocity [cm/s]');
xlabel('Electric field [kV/cm]');
legend
%% impact ionization coefficient
figure(6),
%plot(1./vF*1e8,v_alpha,'DisplayName','Electrons');
errorbar(1./vF*1e8,v_mean_alpha,v_std_alpha,'DisplayName','Electrons');
xlabel('E^{-1} [10^{-8} m/V]');
ylabel('Ionization rate (m^{-1})');
legend