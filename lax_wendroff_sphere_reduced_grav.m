function [unew, vnew, h_new, PT_new, qt_new] = ...
   lax_wendroff_sphere_reduced_grav(dphi, dtheta, dt, g, u, v, h, H,Re, THETA, F, PT, qt)

% This function performs one timestep of the Lax-Wendroff scheme
% applied to the shallow water equations
% Ripa system: http://www4.ncsu.edu/~acherto/papers/Ripa.pdf

% note, need to put accelerations into flux term:
% code this up properly, with proper BCs, etc

% calculate spatial differences
dx=dphi.*(Re).*cos(THETA);
dy=dtheta.*(Re);
dy1=dtheta.*(Re).*cos(THETA);
dy2=dy.*ones(size(THETA));

v1=v.*cos(THETA); 

% PT(:,[1:4])=repmat(PT(:,5),[1 4]);
% PT(:,end-3:end)=repmat(PT(:,end-4),[1 4]);


PT_ref=100.;

rg=(PT-PT_ref)./PT_ref;

% First work out mid-point values in time and space
uh = u.*h;
vh = v.*h;
vh1 = v1.*h;


% continuity equation (calculate mid-point values at 0.5*dt):
h_mid_xt = 0.5.*(h(2:end,:)+h(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(uh(2:end,:)-uh(1:end-1,:));
h_mid_yt = 0.5.*(h(:,2:end)+h(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(vh1(:,2:end)-vh1(:,1:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%POTENTIAL TEMPERATURE=PT
hPT = h.*rg;
uhPT = u.*hPT;
vhPT = v.*hPT;
v1hPT = v1.*hPT;
vPT1 =v1.*rg;
uPT=u.*rg;

hPT_mid_xt = 0.5.*(hPT(2:end,:)+hPT(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(uhPT(2:end,:)-uhPT(1:end-1,:))... 
  +0.125.*dt.*(F(2:end,:)+F(1:end-1,:)).*(vh(2:end,:)+vh(1:end-1,:)); 
hPT_mid_yt = 0.5.*(hPT(:,2:end)+hPT(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(v1hPT(:,2:end)-v1hPT(:,1:end-1)) ...
  +0.125.*dt.*(F(:,2:end)+F(:,1:end-1)).*(vh(:,2:end)+vh(:,1:end-1)); 


PT_mid_xt = 0.5.*(rg(2:end,:)+rg(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(uPT(2:end,:)-uPT(1:end-1,:));
PT_mid_yt = 0.5.*(rg(:,2:end)+rg(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(vPT1(:,2:end)-vPT1(:,1:end-1));

% rg_mid_xt=(PT_mid_xt-PT_ref)./PT_ref;
% rg_mid_yt=(PT_mid_yt-PT_ref)./PT_ref;
rg_mid_xt=PT_mid_xt;
rg_mid_yt=PT_mid_yt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cloud water Content=qt
hqt = h.*qt;
uhqt = u.*hqt;
vhqt = v.*hqt;
v1hqt = v1.*hqt;
vqt1 =v1.*qt;
uqt=u.*qt;

hqt_mid_xt = 0.5.*(hqt(2:end,:)+hqt(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(uhqt(2:end,:)-uhqt(1:end-1,:))... 
 ;% +0.125.*dt.*(F(2:end,:)+F(1:end-1,:)).*(vh(2:end,:)+vh(1:end-1,:)); 
hqt_mid_yt = 0.5.*(hqt(:,2:end)+hqt(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(v1hqt(:,2:end)-v1hqt(:,1:end-1)) ...
  +0.125.*dt.*(F(:,2:end)+F(:,1:end-1)).*(vh(:,2:end)+vh(:,1:end-1)); 


qt_mid_xt = 0.5.*(qt(2:end,:)+qt(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(uqt(2:end,:)-uqt(1:end-1,:));
qt_mid_yt = 0.5.*(qt(:,2:end)+qt(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(vqt1(:,2:end)-vqt1(:,1:end-1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v-phi, or u momentum equation (calculate mid-point values at 0.5*dt):
Ux = uh.*u+0.5.*rg.*g.*h.^2;
Uy = uh.*v1;
uh_mid_xt = 0.5.*(uh(2:end,:)+uh(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(Ux(2:end,:)-Ux(1:end-1,:))... 
  +0.125.*dt.*(F(2:end,:)+F(1:end-1,:)).*(vh(2:end,:)+vh(1:end-1,:)); 
uh_mid_yt = 0.5.*(uh(:,2:end)+uh(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(Uy(:,2:end)-Uy(:,1:end-1)) ...
  +0.125.*dt.*(F(:,2:end)+F(:,1:end-1)).*(vh(:,2:end)+vh(:,1:end-1)); 



% v-theta, or v momentum equation (calculate mid-point values at 0.5*dt):
Vx = uh.*v;
% Vy = vh1.*v+0.5.*g.*h.^2.*cos(THETA);
Vy = vh1.*v;
Vy2 = 0.5.*rg.*g.*h.^2;
vh_mid_xt = 0.5.*(vh(2:end,:)+vh(1:end-1,:)) ...
  -(0.5.*dt./(0.5.*(dx(2:end,:)+dx(1:end-1,:)))).*(Vx(2:end,:)-Vx(1:end-1,:))...
  -0.125.*dt.*(F(2:end,:)+F(1:end-1,:)).*(uh(2:end,:)+uh(1:end-1,:)); 

vh_mid_yt = 0.5.*(vh(:,2:end)+vh(:,1:end-1)) ...
  -(0.5.*dt./(0.5.*(dy1(:,2:end)+dy1(:,1:end-1)))).*(Vy(:,2:end)-Vy(:,1:end-1)) ...
  -(0.5.*dt./(dy)).*(Vy2(:,2:end)-Vy2(:,1:end-1)) ...
    -0.125.*dt.*(F(:,2:end)+F(:,1:end-1)).*(uh(:,2:end)+uh(:,1:end-1)); 

% calculate mid-point value of cos (theta)
c_mid_yt=cos(0.5.*(THETA(:,2:end)+THETA(:,1:end-1)));


% Now use the mid-point values to predict the values at the next timestep
% continuity:
h_new = h(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(uh_mid_xt(2:end,2:end-1)-uh_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).* ...
  (vh_mid_yt(2:end-1,2:end).*c_mid_yt(2:end-1,2:end)-vh_mid_yt(2:end-1,1:end-1).*c_mid_yt(2:end-1,1:end-1));
% h_new=h_new+dt.*h_tend.*0.5.*(h(2:end-1,2:end-1)+h_new);

% h_new=max(h_new,10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uhPT_mid_xt = (uh_mid_xt.*hPT_mid_xt)./h_mid_xt;
vhPT_mid_yt = (vh_mid_yt.*hPT_mid_yt)./h_mid_yt;

% uhPT_mid_xt = uh_mid_xt.*PT_mid_xt;
% vhPT_mid_yt = vh_mid_yt.*PT_mid_yt;



%hPT
hPT_new = hPT(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(uhPT_mid_xt(2:end,2:end-1)-uhPT_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).* ...
  (vhPT_mid_yt(2:end-1,2:end).*c_mid_yt(2:end-1,2:end)-vhPT_mid_yt(2:end-1,1:end-1).*c_mid_yt(2:end-1,1:end-1));
% h_new=h_new+dt.*h_tend.*0.5.*(h(2:end-1,2:end-1)+h_new);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uhqt_mid_xt = (uh_mid_xt.*hqt_mid_xt)./h_mid_xt;
vhqt_mid_yt = (vh_mid_yt.*hqt_mid_yt)./h_mid_yt;

% uhPT_mid_xt = uh_mid_xt.*PT_mid_xt;
% vhPT_mid_yt = vh_mid_yt.*PT_mid_yt;



%hPT
hqt_new = hqt(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(uhqt_mid_xt(2:end,2:end-1)-uhqt_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).* ...
  (vhqt_mid_yt(2:end-1,2:end).*c_mid_yt(2:end-1,2:end)-vhqt_mid_yt(2:end-1,1:end-1).*c_mid_yt(2:end-1,1:end-1));
% h_new=h_new+dt.*h_tend.*0.5.*(h(2:end-1,2:end-1)+h_new);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u-momentum equation:
Ux_mid_xt = uh_mid_xt.*uh_mid_xt./h_mid_xt + 0.5.*rg_mid_xt.*g.*h_mid_xt.^2;
Uy_mid_yt = uh_mid_yt.*vh_mid_yt./h_mid_yt.*c_mid_yt;
uh_new = uh(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(Ux_mid_xt(2:end,2:end-1)-Ux_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).*(Uy_mid_yt(2:end-1,2:end)-Uy_mid_yt(2:end-1,1:end-1));
%  + dt.*u_tendency.*0.5.*(h(2:end-1,2:end-1)+h_new); 


% v-momentum equation:
Vx_mid_xt = uh_mid_xt.*vh_mid_xt./h_mid_xt;
% Vy_mid_yt = vh_mid_yt.*vh_mid_yt./h_mid_yt.*c_mid_yt + 0.5.*g.*h_mid_yt.^2.*c_mid_yt;
Vy_mid_yt = vh_mid_yt.*vh_mid_yt./h_mid_yt.*c_mid_yt;
Vy_mid_yt2 = 0.5.*rg_mid_yt.*g.*h_mid_yt.^2;
vh_new = vh(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(Vx_mid_xt(2:end,2:end-1)-Vx_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).*(Vy_mid_yt(2:end-1,2:end)-Vy_mid_yt(2:end-1,1:end-1)) ...
  - (dt./(dy)).* ... ...
  (Vy_mid_yt2(2:end-1,2:end)-Vy_mid_yt2(2:end-1,1:end-1));

% Calculate new PT
PT_new = (hPT_new./h_new).*PT_ref+PT_ref;
% New reduced gravity
rg_new=(PT_new-PT_ref)./PT_ref;


% add on Coriolis and contribution of orography to pressure gradient:
% g=0;
uh_new = uh_new  +dt.*.5.*(F(2:end-1,2:end-1).*v(2:end-1,2:end-1) - ...
    0.5.*g.*(rg(2:end-1,2:end-1)+rg_new).*(H(3:end,2:end-1)-H(1:end-2,2:end-1))./ ...
    (dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)) ).* ...
    (h(2:end-1,2:end-1)+h_new);   

vh_new = vh_new  -dt.*.5.*(F(2:end-1,2:end-1).*u(2:end-1,2:end-1) + ...
    0.5.*g.*(rg(2:end-1,2:end-1)+rg_new).*(H(2:end-1,3:end)-H(2:end-1,1:end-2))./ ...
    (dy2(2:end-1,2:end-1)+dy2(2:end-1,1:end-2))   ).* ...
    (h(2:end-1,2:end-1)+h_new);


%  + dt.*v_tendency.*0.5.*(h(2:end-1,2:end-1)+h_new); 
% re-calculate u and v. and PT
unew = uh_new./h_new;
vnew = vh_new./h_new;

qt_new = hqt_new./h_new;
