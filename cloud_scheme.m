function [qc,qv,Tc]=cloud_scheme(theta,P,qt)
% [qc,qv,Tc]=cloud_scheme(theta,P,qt)
% Given an array of potential temperature, pressure and total water content
% calculate the cloud water, water vapour and associated cloud temperature
% after the latent heating due to condensation has taken place

global L_v;
global R_v;
L_v=2.5e6;
R_a=287;
c_p=1000;
R_v=461;

T_A=theta*((P./100000).^(R_a/c_p));

fun=@(Tcl)T_A+(L_v/c_p)*(qt-(R_a/R_v)*svp(Tcl)./P)-Tcl;

Tc=fzero(fun,200);

if Tc < T_A
    qv=qt;
    qc=0.;
    Tc=T_A;
else
    qv=(R_a/R_v)*svp(Tc)./P;
    qc=qt-qv;
end


% Clausius-Clapeyron equation
function es=svp(T)
global L_v;
global R_v;

es=(610.7).*exp((L_v/R_v)*(1/273.15-1/T));
