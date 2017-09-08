function [ I,B,Fc ] = solar_flux(P, T, rho,qt,co2_ppm,ch4_ppm,Ap ,emiss,bowen,THETA)
% for radiation model see:
% http://www.bartonlevenson.com/NewPlanetTemps.html

Sflux=1365; 
rhov=rho.*qt./4; % kg/m^3 - note this is scaled down by a factor of 4 
                % (because water vapour is distributed in height)
R_v=461;    %Specific gas constant for water vapour [J/kg K]
sigma=5.67*10^(-8); % Boltzman constant [w/m^2/k^4]



% solar declination angle++++++++++++++++++++++++++++++++++++++++++++++++++
delta=0;
%--------------------------------------------------------------------------

% Insolation
costheta_s=cos(THETA).*cos(delta).*1./pi; % integrated over 1 day
% costheta_s=costheta_s.*2./(pi./2);
I=Sflux.*costheta_s; 

A=0.306;
F=I.*(1-A);


% Partial pressure
e_h2o=T.*rhov.*R_v;
e_co2=P.*co2_ppm.*1e-6;
e_ch4=P.*ch4_ppm.*1e-6;

% Optical depth 
tau_h2o=0.087*(e_h2o).^(1/2);
tau_co2=0.029*(e_co2).^(1/2);
tau_ch4=0.725*(e_ch4).^(1/2);

tau=tau_h2o+tau_co2+tau_ch4;

% Emissivity 
epsilon=1./(1+tau.*0.75);

% epsilon=epsilon.*max((T-265)./30.,0.0001);

% raw greenhouse temperature:
T0 = (F./(epsilon.*sigma)).^(1/4);


% visible optical thickness parametrisation:
tvis=0.46.*(tau-0.723).^0.411; % changed to 0.46
tvis(find(tau<0.724))=0;

% atmospheric cooling loss:
Labs=F.*(1-exp(-tvis));


% sensible and latent heat:
F0=emiss.*sigma.*T0.^4;
Fsi=F-Labs;
Fabs=(1-Ap).*Fsi+emiss.*(F0-F);
Fc=0.369.*Fabs.*tau./(-0.6+2.*tau);

% surface flux:
Fs=F0-Labs-Fc;

I=Fs;
B=emiss.*sigma.*T.^4;

% if(length(find(imag(I))) | length(find(imag(B))))
%     I;
% end

% Surface temperature
% T_s = (I.*(ones(size(Ap))-Ap)./(epsilon.*4.*sigma)).^(1/4);
% T_s=max(T_s, 260);

