function [unew, vnew, h_new] = ...                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   function [u_new, v_new, h_new, PT_new] = ...
   lax_wendroff_sphere_barotropic(dphi, dtheta, dt, g, u, v, h, H,Re, THETA, F)

% This function performs one timestep of the Lax-Wendroff scheme
% applied to the shallow water equations

% note, need to put accelerations into flux term:
% code this up properly, with proper BCs, etc

% calculate spatial differences
dx=dphi.*(Re).*cos(THETA);
dy=dtheta.*(Re);
dy1=dtheta.*(Re).*cos(THETA);
dy2=dy.*ones(size(THETA));

v1=v.*cos(THETA); 



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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v-phi, or u momentum equation (calculate mid-point values at 0.5*dt):
Ux = uh.*u+0.5.*g.*h.^2;
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
Vy2 = 0.5.*g.*h.^2;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u-momentum equation:
Ux_mid_xt = uh_mid_xt.*uh_mid_xt./h_mid_xt + 0.5.*g.*h_mid_xt.^2;
Uy_mid_yt = uh_mid_yt.*vh_mid_yt./h_mid_yt.*c_mid_yt;
uh_new = uh(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(Ux_mid_xt(2:end,2:end-1)-Ux_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).*(Uy_mid_yt(2:end-1,2:end)-Uy_mid_yt(2:end-1,1:end-1));
%  + dt.*u_tendency.*0.5.*(h(2:end-1,2:end-1)+h_new); 


% v-momentum equation:
Vx_mid_xt = uh_mid_xt.*vh_mid_xt./h_mid_xt;
% Vy_mid_yt = vh_mid_yt.*vh_mid_yt./h_mid_yt.*c_mid_yt + 0.5.*g.*h_mid_yt.^2.*c_mid_yt;
Vy_mid_yt = vh_mid_yt.*vh_mid_yt./h_mid_yt.*c_mid_yt;
Vy_mid_yt2 = 0.5.*g.*h_mid_yt.^2;
vh_new = vh(2:end-1,2:end-1) ...
  - (dt./(0.5.*(dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)))).*(Vx_mid_xt(2:end,2:end-1)-Vx_mid_xt(1:end-1,2:end-1)) ...
  - (dt./(0.5.*(dy1(2:end-1,2:end-1)+dy1(2:end-1,1:end-2)))).*(Vy_mid_yt(2:end-1,2:end)-Vy_mid_yt(2:end-1,1:end-1)) ...
  - (dt./(dy)).* ... ...
  (Vy_mid_yt2(2:end-1,2:end)-Vy_mid_yt2(2:end-1,1:end-1));

% add on Coriolis and contribution of orography to pressure gradient:
% g=0;
uh_new = uh_new  +dt.*.5.*(F(2:end-1,2:end-1).*v(2:end-1,2:end-1) - ...
    g.*(H(3:end,2:end-1)-H(1:end-2,2:end-1))./ ...
    (dx(2:end-1,2:end-1)+dx(1:end-2,2:end-1)) ).* ...
    (h(2:end-1,2:end-1)+h_new);   

vh_new = vh_new  -dt.*.5.*(F(2:end-1,2:end-1).*u(2:end-1,2:end-1) + ...
    g.*(H(2:end-1,3:end)-H(2:end-1,1:end-2))./ ...
    (dy2(2:end-1,2:end-1)+dy2(2:end-1,1:end-2))   ).* ...
    (h(2:end-1,2:end-1)+h_new);


%  + dt.*v_tendency.*0.5.*(h(2:end-1,2:end-1)+h_new); 
% re-calculate u and v. and PT
unew = uh_new./h_new;
vnew = vh_new./h_new;
