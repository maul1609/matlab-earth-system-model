% Implementation of the shallow water equations in a spherical coordinate
% system. Paul Connolly, University of Manchester, September 2017

% MODEL SETTINGS+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
IDEALISED_JET=2;
EARTH_WINDS=3;
NO_WINDS=4;

FLAT_TOPO=1;
EARTH_TOPO=2;
%--------------------------------------------------------------------------

% USER SETTINGS++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
interactive_theta=true;
radiation=false;
cloud_diabatic=false; % include the latent heating due to 
                      % condensation in dynamics


add_random_height_noise=false;
initially_geostrophic=true;
viscous_dissipation=true;
smagorinksky=true;
cvis=0.2;
viscous_dissipation_days=0.25;
dissipate_h=false;
initial_winds=EARTH_WINDS; % set to IDEALISED_JET for a jet
nudge=true; % nudge mean model height field to that observed
restart=false;

topo=EARTH_TOPO;

vis=5e6; % turbulent dissipation coefficient

dt_mins              = 0.5;   % Timestep (minutes)
output_interval_mins = 60;    % Time between outputs (minutes)
forecast_length_days = 2;     % Total simulation length (days)

% constants for Earth
g=9.81;       % gravity
rho=.5;       % density [ in[kg/m^3] ]
Re=6.4e6;     % radius of planet [m]
f=2.*pi./((24).*3600); % radians per second - planet's rotation rate

R_a=287;    % Specific gas constant for air [J/kg K]
R_v=461;    % Specific gas constant for water vapour [J/kg K]
c_p=1000;   % specific heat [J/kg K]
L_v=2.5e6;  % latent heat of vapourisation [J/kg]
co2_ppm=800;
ch4_ppm=1.75;
tau_evap=86400;
tau_precip=1000;
tau_rad=86400;

scale_height=1.0.*10e3; % scale height of earth's atmosphere (density and pressure fall by 1/e)

% set-up of model grid

dtheta=1.*pi./180;  % the north-south (latitude coordinate)
dphi=2.*pi./180;    % the east-west (longitude coordinate)

phi=[0:dphi:2.*pi-dphi]; % the longitude mesh 
theta=dtheta./2-85.*pi./180:dtheta:85.*pi./180-dtheta./2; % the latitude mesh
if initial_winds==IDEALISED_JET
    theta=dtheta./2+15.*pi./180:dtheta:85.*pi./180-dtheta./2; % the latitude mesh
end
%--------------------------------------------------------------------------


% MODEL SET UP FROM INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++
[THETA,PHI]=meshgrid(theta,phi);    %grid of [lat,long]

nx=length(phi); % number of points in east-west direction
ny=length(theta); % number of points in north-south direction

if topo==FLAT_TOPO
    H = zeros(nx, ny); % terrain height
    Apland=0.05.*ones(size(H));  
    emiss=0.95.*ones(size(Apland));
    bowen=ones(size(Apland));;
elseif topo==EARTH_TOPO
    load earth_topo_1degx1deg.mat
    [r,c]=size(Hfield);
    Hfield=Hfield(r-nx+1:r,c-ny+1:c);
    Apland=0.3.*ones(size(Hfield));  
    bowen=ones(size(Apland));;
    bowen(find(Hfield(:)<20))=0.1;
    bowen(find(Hfield(:)>=20))=1;
    emiss=0.95.*ones(size(Apland));
    emiss(find(Hfield(:)<20))=1;
    Apland(find(Hfield(:)<20))=0.1;
    latar=[1:nx]; % for albedo at poles
    latar=repmat(latar,[ny 1]);
    Apland((latar(:)<20 | latar(:)>160))=0.1;
    
    
    H=(max(Hfield./g,0));           
end
qt=0.02.*ones(nx,ny);
qc_arr=ones(nx,ny);
load cloud_scheme_data THETA_A P_A QT_A QC_A QV_A TC_A;


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
            u(:,:)=normpdf(THETA,55.*pi./180,1.*pi./180);
            u(:,:)=100.*u(:,:)./max(u(1,:));
            v(:,:)=0.;
            max_wind = 2000;
            %%%%%%%%%%%
%             T=300.*ones(nx,ny);
%             T(find(THETA.*180./pi>45))=290;
            T=50.*cos(THETA)+250;
            %%%%%%%
        case EARTH_WINDS
            v(:,:)=0;
            u(:,:)=0;
            max_wind = 300;
            %%%%%%%%%%%
            T=300.*ones(nx,ny);
            %%%%%%%
        case NO_WINDS
            v(:,:)=0;
            u(:,:)=0;
            max_wind = 300;
            %%%%%%%%%%%
            T=300.*ones(nx,ny);
            T(find(THETA.*180./pi>55))=270;
            %%%%%%%
    end
    PT=T;
    T_s=T(:, :);
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
    
    if initial_winds==EARTH_WINDS
        load ecmwf_msl_pressure_1degx1deg msl_p;
        height=msl_p./g./rho;       %msl_p = mean sea level pressure
        h_nudge=height-H;
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
        ind=find(180./pi.*THETA>50 & 180./pi.*THETA<60 );
        height(ind) = height(ind) + 50.*rand(size(height(ind))).*(abs(F(ind))./3e-4);
        T(ind)=T(ind)+1.*rand(size(T(ind)));
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

%%%%%%%%%%%
% Find original pressure field
P=h.*rho*g;
%%%%%%

%--------------------------------------------------------------------------




% Initialize the 3D arrays where the output data will be stored
u_save = zeros(nx, ny, noutput);
v_save = zeros(nx, ny, noutput);
h_save = zeros(nx, ny, noutput);
th_save = zeros(nx, ny, noutput);
qt_save=zeros(nx,ny,noutput); % Create qt
qc_save=zeros(nx,ny,noutput); % Create qc
T_s_save=zeros(nx,ny,noutput); 
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
    th_save(:,:,i_save) = PT;
    qt_save(:,:,i_save)=qt;
    qc_save(:,:,i_save)=qc_arr;
    T_s_save(:,:,i_save)=T_s;
    t_save(i_save) = (n-1).*dt;
    i_save = i_save+1;

    
  end
   
  
   
  %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  % Call the Lax-Wendroff scheme to move forward one timestep. Note, that
  % this is also coded to do the accelerations, but not the viscous term
  if interactive_theta
      
      if cloud_diabatic % calculate the new potential temperature for use 
                        % dynamics
        % Calculate the new PT due to latent heat++++++++++++++++++++++++++
        % Find new pressure field
        P=h.*rho*g;
        %Find new T
        T=PT.*(P./100000).^(R_a/c_p);
        %%%%%%%
        % cloud look-up table
        qc_arr=interp3(THETA_A,P_A,QT_A,   permute(QC_A,[2 1 3 ])  ,...
        max(min(PT,399),201),max(min(P,1.09e5),20001),...
        max(min(qt,0.039),1.1e-4),'linear');
        qt=max(qt,qc_arr); % ensure that total water > cloud water
        T=T+qc_arr.*L_v./c_p;
        PT=T./(P./100000).^(R_a/c_p);
        %------------------------------------------------------------------
      end
      
      
      % Solve for the dynamics (with reduced gravity term):
      [unew, vnew, h_new, PT_new, qt_new] = ...
          lax_wendroff_sphere_reduced_grav(dphi, dtheta, dt, g, ...
          [u(end-1:end,:); u; u(1:2,:)], [v(end-1:end,:);v;v(1:2,:)], ...
          [h(end-1:end,:);h;h(1:2,:)],...
          [H(end-1:end,:);H;H(1:2,:)],...
          Re,[THETA(end-1:end,:);THETA;THETA(1:2,:)],...
          [F(end-1:end,:);F;F(1:2,:)], ...
          ([PT(end-1:end,:);PT;PT(1:2,:)]), ...
          [qt(end-1:end,:);qt;qt(1:2,:)]); 
    
    
    
    
    
      if cloud_diabatic % calculate the new potential temperature (after 
                        % condensation - take away)
        % Readjust the PT - evaporate the cloud++++++++++++++++++++++++++++
        %Find new pressure field
        P_new=h_new.*rho*g;
        %Find new T
        T_new=PT_new.*(P_new./100000).^(R_a/c_p);
        %%%%%%%
        % cloud look-up table
        qc_arr=interp3(THETA_A,P_A,QT_A,   permute(QC_A,[2 1 3 ])  ,...
        max(min(PT,399),201),max(min(P,1.09e5),20001),...
        max(min(qt,0.039),1.1e-4),'linear');
        qt=max(qt,qc_arr); % ensure that total water > cloud water
        T_new=T_new-qc_new.*L_v./c_p;
%         T_new(2:end-1,:)=T_new(2:end-1,:)-qc_arr(:,2:end-1).*L_v./c_p;
        PT_new=T_new./(P_new./100000).^(R_a/c_p);
        %------------------------------------------------------------------
      end
  else
      
      
      % solve for the dynamics (barotropic version):
      [unew, vnew, h_new, PT_new, qt_new] = ...
          lax_wendroff_sphere_barotropic(dphi, dtheta, dt, g, ...
          [u(end-1:end,:); u; u(1:2,:)], ...
          [v(end-1:end,:);v;v(1:2,:)], ...
          [h(end-1:end,:);h;h(1:2,:)],...
          [H(end-1:end,:);H;H(1:2,:)],Re,...
          [THETA(end-1:end,:);THETA;THETA(1:2,:)],...
          [F(end-1:end,:);F;F(1:2,:)], ...
          ([PT(end-1:end,:);PT;PT(1:2,:)]), ...
          [qt(end-1:end,:);qt;qt(1:2,:)]); 
      
      
  end
      
  %------------------------------------------------------------------------
  unew=unew([2:end-1],:);
  vnew=vnew([2:end-1],:);
  h_new=h_new([2:end-1],:);
  PT_new=PT_new([2:end-1],:);
  qt_new=qt_new([2:end-1],:);
  
  qt_new(find(qt_new(:)<1e-6))=1e-6;
  
  % add nudge terms
  if nudge==true & (n-1).*dt<viscous_dissipation_days.*86400
      h_new(:,:)=h_new(:,:)+(h_nudge(:,2:end-1)-h_new(:,:) ).*1e-4.*dt;
  end
  
  % apply viscosity:
  if viscous_dissipation==true & (n-1).*dt<viscous_dissipation_days.*86400
      % gradient of wind:
      u1=[unew(end,:);unew;unew(1,:)]; 
      u1(:,:)=(u1(:,:)+[u(end,2:end-1);u(:,2:end-1);u(1,2:end-1)])./2;
      u1=[u1(:,1),u1,u1(:,end)];
%       [uy,ux]=gradient(u1);     
%       ux=ux./[dx(end,:); dx(1:end,1:end); dx(1,:)];uy=uy./[dy(end,:);dy(1:end,1:end);dy(1,:)];
      
      v1=[vnew(end,:);vnew;vnew(1,:)]; 
      v1(:,:)=(v1(:,:)+[v(end,2:end-1);v(:,2:end-1);v(1,2:end-1)])./2;
      v1=[v1(:,1),v1,v1(:,end)];
%       [vy,vx]=gradient(v1);
%       vx=vx./[dx(end,:);dx(1:end,1:end);dx(1,:)];vy=vy./[dy(end,:);dy(1:end,1:end);dy(1,:)];
      

%       % divergence of grad:
%       % del^2 of u
%       del2a=divergence([y(end,:);y(1:end,1:end);y(1,:)],[x(end,:); x(1:end,1:end); x(1,:)],uy,ux);
%       % del^2 of v
%       del2b=divergence([y(end,:);y(1:end,1:end);y(1,:)],[x(end,:); x(1:end,1:end); x(1,:)],vy,vx);

      del2a = 1./(Re.^2.*cos(THETA(:,2:end-1))).* ...
          (cos(THETA(:,2:end-1)).*(u1(2:end-1,3:end)-u1(2:end-1,2:end-1))./dtheta - ... 
           (cos(THETA(:,1:end-2)).*(u1(2:end-1,2:end-1)-u1(2:end-1,1:end-2))./dtheta )./dtheta) + ...
           1./(Re.^2.*cos(THETA(:,2:end-1)).^2).* ...
          ((u1(3:end,2:end-1)-u1(2:end-1,2:end-1))./dphi - ... 
           ((u1(2:end-1,2:end-1)-u1(1:end-2,2:end-1))./dphi )./dphi);
       
      del2b = 1./(Re.^2.*cos(THETA(:,2:end-1))).* ...
          (cos(THETA(:,2:end-1)).*(v1(2:end-1,3:end)-v1(2:end-1,2:end-1))./dtheta - ... 
           (cos(THETA(:,1:end-2)).*(v1(2:end-1,2:end-1)-v1(2:end-1,1:end-2))./dtheta )./dtheta) + ...
           1./(Re.^2.*cos(THETA(:,2:end-1)).^2).* ...
          ((v1(3:end,2:end-1)-v1(2:end-1,2:end-1))./dphi - ... 
           ((v1(2:end-1,2:end-1)-v1(1:end-2,2:end-1))./dphi )./dphi);

      % add viscous term to wind fields:
      if(smagorinksky==true)
          vis2=cvis.^2.*Re.*(Re).*cos(THETA(:,2:end-1)).*dphi.*dtheta.* ...
              sqrt(  ((u1(3:end,2:end-1)-u1(1:end-2,2:end-1)) ./(Re.*cos(THETA(:,2:end-1)).*2.*dphi)).^2 + ...
                    ((v1(2:end-1,3:end)-v1(2:end-1,1:end-2)) ./(Re.*2.*dtheta)).^2 + ...
                    0.5.*(...
                    ((u1(2:end-1,3:end)-u1(2:end-1,1:end-2)) ./(Re.*2.*dtheta)) + ...
                    ((v1(3:end,2:end-1)-v1(1:end-2,2:end-1)) ./(Re.*cos(THETA(:,2:end-1)).*2.*dphi)) ...
                    ).^2 ...
                );
%           unew=unew+dt.*vis2.*del2a(2:end-1,2:end-1);
%           vnew=vnew+dt.*vis2.*del2b(2:end-1,2:end-1);
          unew=unew+dt.*vis2.*del2a;
          vnew=vnew+dt.*vis2.*del2b;
      else
          unew=unew+dt.*vis.*del2a;
          vnew=vnew+dt.*vis.*del2b;
      end
      
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
  u(:,2:end-1) = unew(:,[1:end]);
  v(:,2:end-1) = vnew(:,[1:end]);
  
   
  
  PT(:,2:end-1) = PT_new(:,[1:end]);
  qt(:,2:end-1) = qt_new(:,[1:end]);
 
  h(:,2:end-1) = h_new(:,:);
%   h(:,[1 end])=h(:,[2 end-1]);
  
  %%%%%%%
  
  %Find new pressure field
  P=h.*rho*g;
  
  %Find new T - needed for q cloud
  T=PT.*(P./100000).^(R_a/c_p);
  %%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%% 
     
  if radiation
      % cloud look-up table
      qc_arr=interp3(THETA_A,P_A,QT_A,   permute(QC_A,[2 1 3 ])  ,...
        max(min(PT,399),201),max(min(P,1.09e5),20001),...
        max(min(qt,0.039),1.1e-4),'linear');
      qt=max(qt,qc_arr); % ensure that total water > cloud water

      % local albedo:
      Ap=Apland.*(qc_arr<1e-4)+(qc_arr>=1e-4).*0.35;

      [ I,B ,Fc] = solar_flux(P, T, rho,qt-qc_arr, co2_ppm,ch4_ppm,...
          Ap,emiss,bowen,THETA) ;
        
      I=real(I);B=real(B);
      %  shortwave flux, minus the IR emission:
      T=T+(I-B)*dt./tau_rad;
      
      % source of water
      % https://en.wikipedia.org/wiki/Bowen_ratio
      % say Qh+Qe=Fc and adjust Fc to equal Qh, so Qh=Fc-Qe
      % define bowen ratio B=Qh/Qe,: therefore Qh=Fc-1/B*Qh, so Qh+Qh/B=Fc,
      % Qh(1+1./B)=Fc, Qh=Fc/(1+1/B)
      % Fc=Fc./(1+1./bowen);
      lhf=Fc./(1+bowen);
      qt=qt+lhf./L_v./1000.*dt;

      % sink of water (e.g. precipitation):
      qt=max(qt-qt.*min(qc_arr./3e-3.*dt./tau_precip,0.9),1e-6);
      %%%%%%%%%%%%%%%%%%%%%%%   
  end

  % equation for potential temperature:
  PT = T./((P./100000.).^(R_a/c_p));
  
 


  
  
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
