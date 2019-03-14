% an example plotting program with m_map
plot1='vort';
print1='no';

figure('renderer','painters');%maxfigsize
if viscous_dissipation == false
    vis=0.;
end
warning off;
[r,c,p]=size(u_save);
lon=length(phi);
phi2=[phi(lon/2+1:lon)-360.*pi./180 phi(1:lon/2)];
PHI2=[PHI(lon/2+1:lon,:)-360.*pi./180; PHI(1:lon/2,:)];
for i=480:p
    vorticity = zeros(size(u_save(:,:,1)));
    vorticity(2:end-1,2:end-1) = (u_save(2:end-1,1:end-2,i)-u_save(2:end-1,3:end,i)) ...
        + (v_save(3:end,2:end-1,i)-v_save(1:end-2,2:end-1,i));
%     m_proj('robinson','long',[-180 180],'lat',[-90 90]);
    m_proj('ortho','long',180,'lat',80);
    hold off;
    switch plot1
        case 'h'
            F1=h_save([lon/2+1:lon 1:lon/2],:,i); 
        case 'vort'
            F1=vorticity([lon/2+1:lon 1:lon/2],:); 
        otherwise
            disp(['error ']);
            return;
    end
    m_pcolor(phi2.*180./pi,theta.*180./pi,F1');shading flat
    m_grid('fontsize',6,'xticklabels',[],'xtick',[-180:30:360],... 
        'ytick',[-80:20:60 78 80],'yticklabels',[],'linest','-','color',[0 0 0],'linewidth',0.1);
    hold on;
%     m_gshhs_c('color','k');
%     m_quiver(PHI2(1:3:end,1:3:end)'.*180./pi,THETA(1:3:end,1:3:end)'.*180./pi,...
%         u_save([91:3:180 1:3:90],1:3:end,i)',v_save([91:3:180 1:3:90],1:3:end,i)')
    title({['Time (hrs): ',num2str(t_save(i)./3600)],
        [num2str(i),' of ',num2str(p),'; viscosity: ',num2str(vis)]});
    colorbar('horiz')

    pause; %(0.1);
    
    if strcmp(print1,'yes')
        eval(['print -dpng -r25 aframe',num2str(i,'%03d'),'.png']);
    end
end
