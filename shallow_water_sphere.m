% Implementation of the shallow water equations in a spherical coordinate
% system. Paul Connolly, University of Manchester, September 2017

% MODEL SETTINGS+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IDEALISED_JET=2;
EARTH_WINDS=3;
NO_WINDS=4;

FLAT_TOPO=1;
%--------------------------------------------------------------------------

% USER SETTINGS++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


add_random_height_noise=true;
initially_geostrophic=true;
viscous_dissipation=false;
viscous_dissipation_days=20*365; %0.25;
dissipate_h=true;
initial_winds=IDEALISED_JET; % set to IDEALISED_JET for a jet
nudge=true; % nudge mean model height field to that observed
restart=false;

topo=FLAT_TOPO;

vis=5e6; % turbulent dissipation coefficient

dt_mins              = 0.5;   % Timestep (minutes)
output_interval_mins = 60;    % Time between outputs (minutes)
forecast_length_days = 10;     % Total simulation length (days)

% constants for Earth
g=10.44;       % gravity
rho=.19;       % density [ in[kg/m^3] ]
Re=5.4155760e7;     % radius of planet [m]
f=2.*pi./((10.44).*3600); % radians per second - planet's rotation rate


scale_height=60e3; % scale height of earth's atmosphere (density and pressure fall by 1/e)

% set-up of model grid

dtheta=0.2.*pi./180;  % the north-south (latitude coordinate)
dphi=1.*pi./180;    % the east-west (longitude coordinate)

phi=[0:dphi:2.*pi-dphi]; % the longitude mesh 
theta=dtheta./2-85.*pi./180:dtheta:85.*pi./180-dtheta./2; % the latitude mesh
if initial_winds==IDEALISED_JET
    %theta=dtheta./2+15.*pi./180:dtheta:85.*pi./180-dtheta./2; % the latitude mesh
    theta=dtheta./2+65.*pi./180:dtheta:86.5.*pi./180-dtheta./2; % the latitude mesh
end
%--------------------------------------------------------------------------


% MODEL SET UP FROM INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++
[THETA,PHI]=meshgrid(theta,phi);    %grid of [lat,long]

nx=length(phi); % number of points in east-west direction
ny=length(theta); % number of points in north-south direction

if topo==FLAT_TOPO
    H = zeros(nx, ny); % terrain height
end


F=2.*f.*sin(THETA); % Coriolis parameter
F=max(abs(F),1e-5).*sign(F);        

dt = dt_mins*60.0; % Timestep (s)
output_interval = output_interval_mins*60.0; % Time between outputs (s)
forecast_length = forecast_length_days*24.0*3600.0; % Forecast length (s)
nt = fix(forecast_length/dt)+1; % Number of timesteps          
timesteps_between_outputs = fix(output_interval/dt);
noutput = ceil(nt/timesteps_between_outputs); % Number of output frames       


if restart == false
    % Initialize the wind to rest
    u=zeros(nx, ny);
    v=zeros(nx, ny);

    switch initial_winds
        case IDEALISED_JET
            u(:,:)=normpdf(THETA, 77.*pi./180,1.*pi./180);
            u(:,:)=100.*u(:,:)./max(u(1,:));
            v(:,:)=0.;
            max_wind = 2000;
            %%%%%%%
        case NO_WINDS
            v(:,:)=0;
            u(:,:)=0;
            max_wind = 300;
            %%%%%%%%%%%
    end
    u_nudge=u(1,:);     
    % SET UP HEIGHT FIELD FROM OBSERVED WINDS (GEOSTROPHIC BALANCE) +++++++
    height=zeros(nx, ny);

    height(:,end)=scale_height;         
    for i=ny-1:-1:1            
       height(:,i)=height(:,i+1)+ ...
           0.5.*(F(:,i+1)+F(:,i)).* ...
           0.5.*(u(:,i+1)+u(:,i)).* ...
           (Re)./g.*dtheta;
    end
    
    
    % ARRAYS TO HOLD DEL^2 TERMS+++++++++++++++++++++++++++++++++++++++++++
    del2a=zeros(nx+2,ny+2);
    del2b=zeros(nx+2,ny+2);
    del2c=zeros(nx+2,ny+2);

    % DEFINE X, Y, DX, DY
    x=PHI.*cos(THETA).*Re;
    y=THETA.*Re;
    dx=dphi.*cos(THETA).*Re;
    dy=dtheta.*Re.*ones(size(THETA));
    h_tend=0;
    u_accel=0.;
    v_accel=0.;
    %----------------------------------------------------------------------



    % ADD SOME RANDOM NOISE TO HEIGHT FIELD++++++++++++++++++++++++++++++++
    if add_random_height_noise
        ind=find(180./pi.*THETA>74 & 180./pi.*THETA<80 );
        height(ind) = height(ind) + 1000.*0.6e5.*rand(size(height(ind)))./height(ind);
    end
    %----------------------------------------------------------------------



    % SET WINDS TO BE IN GEOSTROPHIC BALANCE (CENTRED DIFFERENCE)++++++++++
    if initially_geostrophic
        % calculate spatial differences
        dx=dphi.*(Re).*cos(THETA);
        dya=dtheta.*(Re);
        dy1=dtheta.*(Re).*cos(THETA);
        dy2=dy.*ones(size(THETA));
       % Centred spatial differences to compute geostrophic wind
       u(:,2:end-1)=-0.5.*g.*(height(:,3:end)-height(:,1:end-2)) ./ ...
           (dya.*(F(:,2:end-1)));

       v(2:end-1,:)=g.*(height(3:end,:)-height(1:end-2,:)) ./ ...
           ((dx(2:end-1,:)+dx(1:end-2,:)).*(F(2:end-1,:)));

       % periodic:
       v(1,:)=g.*(height(2,:)-height(end,:)) ./ ...
           ((dx(1,:)+dx(end,:)).*(F(1,:)));;
       v(end,:)=g.*(height(1,:)-height(end-1,:)) ./ ...
           ((dx(end,:)+dx(1,:)).*(F(end,:)));;

       u(:,1)=u(:,2);
       u(:,end)=u(:,end-1);


       % Don't allow the initial wind speed to exceed max_wind anywhere
       u(find(u>max_wind)) = max_wind;
       u(find(u<-max_wind)) = -max_wind;
       v(find(v>max_wind)) = max_wind;
       v(find(v<-max_wind)) = -max_wind;
           
       
       h_nudge=height-H;
       
    end
    %----------------------------------------------------------------------
   

    % Define h as the depth of the fluid (whereas "height" is the height of
    % the upper surface)
    h = height - H;

else
    u=u_save(:,:,end);
    v=v_save(:,:,end);
    h=h_save(:,:,end);
end


%--------------------------------------------------------------------------




% Initialize the 3D arrays where the output data will be stored
u_save = zeros(nx, ny, noutput);
v_save = zeros(nx, ny, noutput);
h_save = zeros(nx, ny, noutput);
t_save = zeros(1, noutput);


% Index to stored data
i_save = 1;




% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MAIN LOOP
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for n = 1:nt
  % Every fixed number of timesteps we store the fields
  if mod(n-1,timesteps_between_outputs) == 0           
    max_u = sqrt(max(u(:).*u(:)+v(:).*v(:)));
    disp(['Time = ' num2str((n-1)*dt/3600) ...
	  ' hours (max ' num2str(forecast_length_days*24) ...
		   '); max(|u|) = ' num2str(max_u)]);
    u_save(:,:,i_save) = u;
    v_save(:,:,i_save) = v;
    h_save(:,:,i_save) = h;
    t_save(i_save) = (n-1).*dt;
    i_save = i_save+1;

    
  end
   
  
   
  %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
  % solve for the dynamics (barotropic version):
  [unew, vnew, h_new] = ...
      lax_wendroff_sphere_barotropic(dphi, dtheta, dt, g, ...
      [u(end-1:end,:); u; u(1:2,:)], ...
      [v(end-1:end,:);v;v(1:2,:)], ...
      [h(end-1:end,:);h;h(1:2,:)],...
      [H(end-1:end,:);H;H(1:2,:)],Re,...
      [THETA(end-1:end,:);THETA;THETA(1:2,:)],...
      [F(end-1:end,:);F;F(1:2,:)]); 
      
      
      
  %------------------------------------------------------------------------
  unew=unew([2:end-1],:);
  vnew=vnew([2:end-1],:);
  h_new=h_new([2:end-1],:);
  
  
  % apply viscosity:
  if viscous_dissipation==true & (n-1).*dt<viscous_dissipation_days.*86400
      % gradient of wind:
      u1=[unew(end,:);unew;unew(1,:)]; u1=[u1(:,1),u1,u1(:,end)];
      [uy,ux]=gradient(u1);     
      ux=ux./[dx(end,:); dx(1:end,1:end); dx(1,:)];uy=uy./[dy(end,:);dy(1:end,1:end);dy(1,:)];
      
      v1=[vnew(end,:);vnew;vnew(1,:)]; v1=[v1(:,1),v1,v1(:,end)];
      [vy,vx]=gradient(v1);
      vx=vx./[dx(end,:);dx(1:end,1:end);dx(1,:)];vy=vy./[dy(end,:);dy(1:end,1:end);dy(1,:)];
      

      % divergence of grad:
      del2a=divergence([y(end,:);y(1:end,1:end);y(1,:)],[x(end,:); x(1:end,1:end); x(1,:)],uy,ux);
      del2b=divergence([y(end,:);y(1:end,1:end);y(1,:)],[x(end,:); x(1:end,1:end); x(1,:)],vy,vx);

      % add viscous term to wind fields:
      unew=unew+dt.*vis.*del2a(2:end-1,2:end-1);
      vnew=vnew+dt.*vis.*del2b(2:end-1,2:end-1);
      
      if dissipate_h == true
          % gradient of h:
          h1=[h_new(end,:);h_new;h_new(1,:)]; h1=[h1(:,1),h1,h1(:,end)];
          [hy,hx]=gradient(h1);
          hx=hx./[dx(end,:); dx(1:end,1:end); dx(1,:)];hy=hy./[dy(end,:);dy(1:end,1:end);dy(1,:)];
          % divergence of grad:
          del2c=divergence([y(end,:);y(1:end,1:end);y(1,:)],[x(end,:); x(1:end,1:end); x(1,:)],hy,hx);
          % add viscous term to wind fields:
          h_new=h_new+dt.*vis.*del2c(2:end-1,2:end-1);

      end
  end
  
  
  
  
  
  
  % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  % set the boundary conditions
  % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if nudge==true & (n-1).*dt<viscous_dissipation_days.*86400
      unew(:,:)=unew(:,:)+(u_nudge(:,2:end-1)-unew(:,:) )./6e5.*dt;
  end
  u(:,2:end-1) = unew(:,[1:end]);
  v(:,2:end-1) = vnew(:,[1:end]);
  
   
  % add nudge terms
%   if nudge==true & (n-1).*dt<viscous_dissipation_days.*86400
%       h_new(:,:)=h_new(:,:)+(h_nudge(:,2:end-1)-h_new(:,:) ).*1e-4.*dt;
%   end
  
 
  h(:,2:end-1) = h_new(:,:);
%   h(:,[1 end])=h(:,[2 end-1]);
  
  %%%%%%%
  
  
  %%%%%%%%%%%%%%%%%%%%%%% 
     


 


  
  
  % Don't allow the initial wind speed to exceed max_wind anywhere
  u(find(u>max_wind)) = max_wind;
  u(find(u<-max_wind)) = -max_wind;
  v(find(v>max_wind)) = max_wind;
  v(find(v<-max_wind)) = -max_wind;
  % -----------------------------------------------------------------------
  
end
% -------------------------------------------------------------------------
% END OF MAIN LOOP
% -------------------------------------------------------------------------
