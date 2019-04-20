function fourier_wave_number(fileName)
% uses findpeaks and fourier analysis to track the current wave number and
% the rotation speed (i.e. using phase information from the fft).

n_files=length(fileName); 

load(fileName{1},'t_save');
np=length(t_save);

wave_number=zeros(n_files,np,1);
rotation_rate=zeros(n_files,np-1,1);
phase=zeros(n_files,np,1);
u_jet2=zeros(n_files,1);
% figure;
for j=1:n_files
    load(fileName{j},'t_save');

    dt_sec=(t_save(2)-t_save(1));
    
    load(fileName{j},'u_save');
    load(fileName{j},'v_save');
    load(fileName{j},'u_jet');
    
    u_jet2(j)=u_jet;
    
    
    
    [r,c,p]=size(v_save);
    
    Fs=r;
    T=1/Fs;
    L=r;
    t = (0:L-1)*T;

    f = Fs*(0:(L/2))/L;

    for i=1:np
        % Calculate vorticity
        vorticity = zeros(size(u_save(:,:,1)));
        vorticity(2:end-1,2:end-1) = (u_save(2:end-1,1:end-2,i)-u_save(2:end-1,3:end,i)) ...
            + (v_save(3:end,2:end-1,i)-v_save(1:end-2,2:end-1,i));
        %--
        X=mean(-v_save(:,30:35,i),2)';
        Y = fft(X);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);



        [pk,locs]=findpeaks([X X(1:3)]);
        if (numel(locs)>0)
            if (locs(1) == locs(end)-c)
                locs(end)=[];
            end
        end

        ind=numel(locs);
        
        phs=unwrap(angle(Y)).*180./pi;
%         phs=(angle(Y)).*180./pi;

        disp(['wave number is : ',num2str(f(ind+1)), '; phase: ',num2str(phs(ind+1))]);
        phase(j,i)=phs(ind+1);
        wave_number(j,i)=f(ind+1);

    end

    phaseold(j,:)=phase(j,:);
%     phase(j,:)=unwrap(phaseold(j,:).*pi./180).*180./pi./wave_number(j,:);
    phase(j,:)=unwrap(phaseold(j,:).*pi./180).*180./pi./(wave_number(j,:));
    rotation_rate(j,:)=-diff(phase(j,:),1,2) ./ ...
        dt_sec.*86400.*365.25./360; % full rotations per year
end


cmap_lev=64;
map=jet(cmap_lev);
u1=linspace(u_jet2(1),u_jet2(end),cmap_lev);

mean_rot=zeros(n_files,9);
std_rot=zeros(n_files,9);
for j=1:n_files
    for i=1:9
        ind=find(wave_number(j,:)==i);
        mean_rot(j,i)=mean(rotation_rate(j,ind(1:end-2)));
        std_rot(j,i)=std(rotation_rate(j,ind(1:end-2)));

        h=plot(i,mean_rot(j,i),'ko','markersize',10);hold on;
        if(n_files > 1)
            row=interp1(u1,1:cmap_lev,u_jet(j),'nearest');
        else
            row=1;
        end
        set(h,'markerfacecolor',map(row,:));

    end
end
xlabel('wave number (per full rotation)');
ylabel('mean rotation rate (rotations per year)');
h=colorbar
ylabel(h,'speed of jet (m/s)')
caxis([u_jet2(1) u_jet2(end)+0.01]);



