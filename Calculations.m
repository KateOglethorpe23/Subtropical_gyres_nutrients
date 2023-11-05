% Kate Oglethorpe
% Sep 2022

% This M-file calculates the nitrate supply by mesoscale eddy stirring and 
% microscale turbulence to the winter mixed layer and thermocline of the 
% subtropical gyres. 

close all
clear

%% Load variable data

%variables
N_file = 'woa18_all_n00_01.nc';              %ntirate (micromol kg-1)
N_z = ncread(N_file, 'n_an');
T_file = 'woa18_decav_t00_01.nc';            %temperature (degC)
T_z = ncread(T_file, 't_an');             
S_file = 'woa18_decav_s00_01.nc';            %salinity (unitless)
S_z = ncread(S_file, 's_an');
w_mld_file = 'woa18_A5B7_M0200_01.nc';       %winter mixed layer depth (m)
w_mld_z = ncread(w_mld_file, 'M_an');  
load WOA18_K.mat                             %eddy diffusivity (m2 s-1)
K_iso_z = WOA18_K.K_WOA;
load global_kappa_epsilon_for_Katherine      %diapycnal diffusivity (m2 s-1)    
K_dia_z = A4.K;
load chl_final.mat %from 'Chl_Feb_Aqua_Modis.txt' -- chlorophyll conc (X)

%dimensions
lon_woa = ncread(T_file, 'lon');             %N, T, S, w_mld 
lat_woa = ncread(T_file, 'lat');                     
z_woa = ncread(T_file, 'depth');             
lon_Groeskamp = WOA18_K.x;                   %K_iso 
lat_Groeskamp = WOA18_K.y;                   
z_Groeskamp = WOA18_K.z;                     
lon_Whalen = A4.lon;                         %K_dia 
lat_Whalen = A4.lat;                         
z_Whalen = A4.z;                                    

%% Edit data

%reorder all dimensions excluding K_iso
T_z = permute(T_z, [2,1,3]);                	 
S_z = permute(S_z, [2,1,3]);                     
N_z = permute(N_z, [2,1,3]);
w_mld_z = permute(w_mld_z, [2, 1]);
K_iso_z = permute(K_iso_z, [2,1,3]);

% %infill missing K_dia at 100m
% K_dia_new = nan(length(lat_Whalen),...       
%     length(lon_Whalen),length(z_Whalen)+1);                     
% K_dia_new(:,:,1) = K_dia_z(:,:,1); 
% K_dia_new(:,:,2:10) = K_dia_z;
% z_Whalen_new = nan(length(z_Whalen)+1,1); 
% z_Whalen_new(1) = 0;
% z_Whalen_new(2:10) = z_Whalen;

%shift K_iso lonxlat grid
K_iso_2nd = K_iso_z(:,[1:length(lon_Groeskamp)/2],:);
K_iso_1st = K_iso_z(:,[(length(lon_Groeskamp)/2)+1:length(lon_Groeskamp)],:);
K_iso_new = [K_iso_1st K_iso_2nd];

% %shift Kdia lonxlat grid
% lon_2nd = lon_Whalen([1:length(lon_Whalen)/2]);
% lon_1st = -flip(lon_2nd);
% lon_new = [lon_1st lon_2nd];   
% K_dia_2nd = K_dia_new(:,[1:length(lon_Whalen)/2],:);
% K_dia_1st = K_dia_new(:,[(length(lon_Whalen)/2):length(lon_Whalen)-1],:);
% K_dia_new = [K_dia_1st K_dia_2nd];

lon_2nd = lon_Whalen([1:length(lon_Whalen)/2]);
lon_1st = -flip(lon_2nd);
lon_new = [lon_1st lon_2nd];  
K_dia_2nd = K_dia_z(:,[1:length(lon_Whalen)/2],:);
K_dia_1st = K_dia_z(:,[(length(lon_Whalen)/2):length(lon_Whalen)-1],:);
K_dia_new = [K_dia_1st K_dia_2nd];

%map K_dia on 1 deg res
K_dia_1deg = interp3(lon_new, lat_Whalen,z_Whalen, K_dia_new, ...
    double(lon_woa'), double(lat_woa),z_Whalen);

%% Extract data for top 1000m of world's subtropical gyres

%define dims
glob_lon = lon_woa                            
glob_lat = lat_woa
z_1000 = -z_woa(1:47); 
% z_Whalen = -z_Whalen_new(1:6);

z_Whalen = -z_Whalen(1:4);

%top 1000m
S_z = S_z(:,:,1:47);       
T_z = T_z(:,:, 1:47);
N_z = N_z(:,:, 1:47);
K_iso_z = K_iso_new(:,:, 1:47);
K_dia_z = K_dia_1deg(:,:, 1:4);


%% Load mask for subtropical gyres

%masks (based on mask of chl < 0.1 mg m^-^3)
load NA_mask.mat
load SA_mask.mat
load SP1_mask.mat
load SP2_mask.mat
load NP1_mask.mat
load NP2_mask.mat
load I_mask.mat

%outlines
load NA_lat_chl.mat
load NA_lon_chl.mat
load SA_lat_chl.mat
load SA_lon_chl.mat
load I_lat_chl.mat
load I_lon_chl.mat
load NP1_lat_chl.mat
load NP1_lon_chl.mat
load NP2_lat_chl.mat
load NP2_lon_chl.mat
load SP1_lat_chl.mat
load SP1_lon_chl.mat
load SP2_lat_chl.mat
load SP2_lon_chl.mat

%plot
figure
m_proj('miller');
m_pcolor(glob_lon, glob_lat, chl_final)
hold on
m_line(NA_lon_chl, NA_lat_chl, 'color','k','linewi',0.3)
m_hatch(NA_lon_chl, NA_lat_chl,'single',30,2,'color','k'); 
m_line(SA_lon_chl, SA_lat_chl, 'color','k','linewi',0.3)
m_hatch(SA_lon_chl, SA_lat_chl,'single',30,2,'color','k'); 
m_line(I_lon_chl, I_lat_chl, 'color','k','linewi',0.3)
m_hatch(I_lon_chl, I_lat_chl,'single',30,2,'color','k'); 
m_line(NP1_lon_chl, NP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP1_lon_chl, NP1_lat_chl,'single',30,2,'color','k'); 
m_line(SP1_lon_chl, SP1_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP1_lon_chl, SP1_lat_chl,'single',30,2,'color','k'); 
m_line(NP2_lon_chl, NP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(NP2_lon_chl, NP2_lat_chl,'single',30,2,'color','k'); 
m_line(SP2_lon_chl, SP2_lat_chl, 'color','k','linewi',0.3)
m_hatch(SP2_lon_chl, SP2_lat_chl,'single',30,2,'color','k'); 
m_grid('fontsize',9);
m_coast('color','k')
set(gcf,'color','w')
c = colorbar
cmocean('algae')
caxis([0 6e4])
c.Label.String = 'Chl-a concentration (mg m^-^3)'
set(gca,'Fontsize',9)

%% Plot missing K_dia data 0-500m
% 
% %plot
% figure
% m_proj('miller','lat',[-60 60],'lon',[-180 180]);
% m_pcolor(glob_lon, glob_lat, mean(K_dia_z(:,:,[1:3]),3))
% hold on
% m_grid('fontsize',9);
% m_coast('color','k')
% set(gcf,'color','w')
% c = colorbar
% %cmocean('thermal')
% caxis([0 6e-5])
% c.Label.String = 'm^2 s^-^1'
% set(gca,'Fontsize',9)
% title('Average diapycnal diffusivity 0-500m')

%% Convert N conc (mg m-3) to N (mol m-3)

pressure (dbar) of upper ocean of gyre
p_z = nan([length(glob_lat) length(glob_lon) 47]);                        
for lat =1:length(glob_lat);
   for lon = 1:length(glob_lon);               
       p_z(lat,lon,:) = gsw_p_from_z(z_1000,glob_lat(lat));
   end
end

density (kg m^-^3) of upper ocean gyre
[gamma_n_z] = eos80_legacy_gamma_n(S_z,T_z,p_z,glob_lon,glob_lat); % density (kg m^-3)
%save('gamma_n_z.mat','gamma_n_z');
%load gamma_n_z.mat

density = gamma_n_z + 1000;          
nitrate_z= N_z.* density; 
nitrate_z= nitrate_z./ 1e6; 

%% Map variable data on density surfaces

%map nitrate on [20.98:0.02:28.82] density surfaces (mol m-3)
N = nan([length(glob_lat) length(glob_lon) 393]);                       
for lat =1:length(glob_lat)                         
    for lon = 1:length(glob_lon)               
        p_orig = squeeze(gamma_n_z(lat, lon, :)); 
        N_orig = squeeze(nitrate_z(lat, lon, :));
        p_new = [20.98:0.02:28.82]';                       
        if sum(isnan(p_orig)) < 30 
            a = find(~isnan(p_orig));   %remove nans
            N_no_NaNs = N_orig(a);
            p_no_nans = p_orig(a);
            [uniq_N, idx_N] =  unique(N_no_NaNs,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_N, idx_p);
            uniq_N = N_no_NaNs(idx);
            uniq_p = p_no_nans(idx);
            N_new = interp1(uniq_p,uniq_N,p_new); 
            N(lat, lon, :) = permute(N_new,[3 2 1]); 
        else
            N(lat, lon, :) = nan(1,1,393);
        end
    end
end
%save('N2.mat','N')

%map nitrate on [20.99:0.02:28.81] density surfaces (mol m-3)
N_offset = nan([length(glob_lat) length(glob_lon) 392]);                
for lat =1:length(glob_lat)                         
    for lon = 1:length(glob_lon)               
        p_orig = squeeze(gamma_n_z(lat, lon, :)); 
        N_orig = squeeze(nitrate_z(lat, lon, :));
        p_new = [20.99:0.02:28.81]'; 
        if sum(isnan(p_orig)) < 30 
            a = find(~isnan(p_orig));
            N_no_NaNs = N_orig(a);
            p_no_nans = p_orig(a);
            [uniq_N, idx_N] =  unique(N_no_NaNs,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_N, idx_p);
            uniq_N = N_no_NaNs(idx);
            uniq_p = p_no_nans(idx);
            N_new = interp1(uniq_p,uniq_N,p_new);  
            N_offset(lat, lon, :) = permute(N_new,[3 2 1]); 
        else
            N_offset(lat, lon, :) = nan(1,1,392);
        end
    end
end

%3D array; lon x lat x depth
A = repmat(z_1000,1,length(glob_lon));                    
glob_z = repmat(A,[1 1 length(glob_lat)]);
glob_z = permute(glob_z, [3 2 1]);

%map depth on density surfaces [20.98:0.02:28.82] (m)
Z = nan([length(glob_lat) length(glob_lon) 393]);                       
for lat =1:length(glob_lat)                        
    for lon = 1:length(glob_lon)               
        p_orig = squeeze(gamma_n_z(lat, lon,:)); 
        z_orig = squeeze(glob_z(lat, lon,:));
        p_new = [20.98:0.02:28.82]';
        if sum(isnan(p_orig)) < 30 
            a = find(~isnan(p_orig));
            z_no_NaNs = z_orig(a);
            p_no_nans = p_orig(a);
            [uniq_z, idx_z] =  unique(z_no_NaNs,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_z, idx_p);
            uniq_z = z_no_NaNs(idx);
            uniq_p = p_no_nans(idx);
            z_new = interp1(uniq_p,uniq_z,p_new); 
            Z(lat, lon,:) = permute(z_new,[3 2 1]);
        else
            Z(lat, lon,:) = nan(1,1,393);
        end
    end
end
%save('Z.mat','Z')

%map depth on density surfaces [20.99:0.02:28.81] (m)
Z_offset = nan([length(glob_lat) length(glob_lon) 392]);                
for lat =1:length(glob_lat)                         
    for lon = 1:length(glob_lon)               
        p_orig = squeeze(gamma_n_z(lat, lon,:)); 
        z_orig = squeeze(glob_z(lat, lon,:));
        p_new = [20.99:0.02:28.81]';
        if sum(isnan(p_orig)) < 30 
            a = find(~isnan(p_orig));   %remove nans
            z_no_NaNs = z_orig(a);
            p_no_nans = p_orig(a);
            [uniq_z, idx_z] =  unique(z_no_NaNs,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_z, idx_p);
            uniq_z = z_no_NaNs(idx);
            uniq_p = p_no_nans(idx);
            z_new = interp1(uniq_p,uniq_z,p_new); 
            Z_offset(lat, lon,:) = permute(z_new,[3 2 1]);
        else
            Z_offset(lat, lon,:) = nan(1,1,392);
        end
    end
end

%map depth on density surfaces [20.97:0.02:28.83] (m)
Z_offset_2 = nan([length(glob_lat) length(glob_lon) 394]);              
for lat =1:length(glob_lat)                         
    for lon = 1:length(glob_lon)               
        p_orig = squeeze(gamma_n_z(lat, lon,:)); 
        z_orig = squeeze(glob_z(lat, lon,:));
        p_new = [20.97:0.02:28.83]';
        if sum(isnan(p_orig)) < 30 
            a = find(~isnan(p_orig));
            z_no_NaNs = z_orig(a);
            p_no_nans = p_orig(a);
            [uniq_z, idx_z] =  unique(z_no_NaNs,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_z, idx_p);
            uniq_z = z_no_NaNs(idx);
            uniq_p = p_no_nans(idx);
            z_new = interp1(uniq_p,uniq_z,p_new); 
            Z_offset_2(lat, lon,:) = permute(z_new,[3 2 1]);
        else
            Z_offset_2(lat, lon,:) = nan(1,1,394);
        end
    end
end

%get density of wml depth (kg m-3)
w_mld_dens = nan(length(glob_lat),length(glob_lon));
for lat =1:length(glob_lat);
    for lon =1:length(glob_lon);
        if sum(isnan(w_mld_z(lat,lon))) > 0;
            w_mld_dens(lat,lon) = NaN;
        else
            depth = w_mld_z(lat,lon);
            p = gsw_p_from_z(-depth,glob_lat(lat));
            [val idx] = min(abs(z_1000--depth));
            w_mld_dens(lat,lon) = eos80_legacy_gamma_n(S_z(lat,lon,idx),T_z(lat,lon,idx),p,glob_lon(lon),glob_lat(lat));
%             w_mld_dens(lat,lon) = round(dens/0.02)*0.02;
        end
    end
end
%save('w_mld_dens.mat','w_mld_dens')

%map Kiso on density surfaces (m2 s-1)
K_iso = nan([length(glob_lat) length(glob_lon) 393]);
for lat =1:length(glob_lat);
    for lon = 1:length(glob_lon);       
        p_orig = reshape((gamma_n_z(lat ,lon, :)),[47,1]); 
        K_iso_orig = reshape((K_iso_z(lat ,lon,:)),[47,1]);
        p_new = [20.98:0.02:28.82]'; 
        if sum(isnan(p_orig)) < 30 ;
            a = find(~isnan(p_orig));
            K_iso_no_nans = K_iso_orig(a);
            p_no_nans = p_orig(a);
            [uniq_K_iso, idx_K_iso] =  unique(K_iso_no_nans,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_K_iso, idx_p);
            uniq_K_iso = K_iso_no_nans(idx);
            uniq_p = p_no_nans(idx);
            K_iso_new = interp1(uniq_p,uniq_K_iso,p_new); 
            K_iso(lat ,lon,:) = permute(K_iso_new,[3 2 1]); 
        else
            K_iso(lat ,lon,:) = nan(1,1,393);
        end
    end
end
%save('K_iso.mat','K_iso')

%map Kdia on density surfaces (m2 s-1)
K_dia = nan([length(glob_lat) length(glob_lon) 393]);
for lat =1:length(glob_lat)
    for lon = 1:length(glob_lon)       
        A = find(z_1000 == -300 | z_1000 == -500 ...
            | z_1000 == -700 | z_1000 == -900); 
        p_orig = squeeze(gamma_n_z(lat ,lon,A)); 
        K_dia_orig = squeeze(K_dia_z(lat ,lon,:));
        p_new =  [20.98:0.02:28.82]'; 
        if sum(isnan(p_orig)) < 3 
            a = find(~isnan(p_orig));
            K_dia_no_nans = K_dia_orig(a);
            p_no_nans = p_orig(a);
            [uniq_K_dia, idx_K_dia] =  unique(K_dia_no_nans,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_K_dia, idx_p);
            uniq_K_dia = K_dia_no_nans(idx);
            uniq_p = p_no_nans(idx);
            K_dia_new = interp1(uniq_p,uniq_K_dia,p_new);
            shallowest_gamma_n = round(gamma_n_z(lat, lon, 1) / 0.02) * 0.02;
            shallowest_gamma_n_idx = find(round(p_new,2)==round(shallowest_gamma_n,2));
            m_300_gamma_n = find(~isnan(K_dia_new), 1);
            K_dia_newer =  nan(1, length(p_new))';
            K_dia_newer(shallowest_gamma_n_idx:(m_300_gamma_n-1)) = K_dia_new(find(~isnan(K_dia_new), 1));
            K_dia_newer(m_300_gamma_n:end) = K_dia_new(m_300_gamma_n:end);
            K_dia(lat ,lon,:) = permute(K_dia_newer,[3 2 1]);    
        else
            K_dia(lat ,lon,:) = nan(1,1,393);
        end
    end
end
%save('K_dia.mat','K_dia')

%get offset diffusivity (m2 s-1)
K_dia_offset = nan([length(glob_lat) length(glob_lon) 392]);
for lat =1:length(glob_lat)
    for lon =1:length(glob_lon)       
        A = find(z_1000 == -300 | z_1000 == -500 ...
            | z_1000 == -700 | z_1000 == -900); 
        p_orig = squeeze(gamma_n_z(lat ,lon,A)); 
        K_dia_orig = squeeze(K_dia_z(lat ,lon,:));
        p_new = [20.99:0.02:28.81]'; 
        if sum(isnan(p_orig)) < 3 
            a = find(~isnan(p_orig));
            K_dia_no_nans = K_dia_orig(a);
            p_no_nans = p_orig(a);
            [uniq_K_dia, idx_K_dia] =  unique(K_dia_no_nans,'stable'); %remove duplicates
            [uniq_p, idx_p] =  unique(p_no_nans,'stable');
            idx = intersect(idx_K_dia, idx_p);
            uniq_K_dia = K_dia_no_nans(idx);
            uniq_p = p_no_nans(idx);
            K_dia_new = interp1(uniq_p,uniq_K_dia,p_new);
            %infill shallow diffusities (<300m)
            shallowest_gamma_n = round(gamma_n_z(lat, lon, 1) / 0.02) * 0.02;
            if mod(shallowest_gamma_n / 0.01, 2) == 0
                shallowest_gamma_n = shallowest_gamma_n + 0.01;
            end
            shallowest_gamma_n_idx = find(round(p_new,2)==round(shallowest_gamma_n,2));
            m_300_gamma_n = find(~isnan(K_dia_new), 1);
            K_dia_newer =  nan(1, length(p_new))';
            K_dia_newer(shallowest_gamma_n_idx:(m_300_gamma_n-1)) = K_dia_new(find(~isnan(K_dia_new), 1));
            K_dia_newer(m_300_gamma_n:end) = K_dia_new(m_300_gamma_n:end);
            K_dia_offset(lat ,lon,:) = permute(K_dia_newer,[3 2 1]); 
        else
            K_dia_offset(lat ,lon,:) = nan(1,1,392);
        end
    end
end
%save('K_dia_offset.mat','K_dia_offset')

%% Calculate eddy diffusive N flux 

%calculate isopycnal N gradient on [20.98:0.02:28.82] (mol m-4)
N_grad_iso_lon = nan([180,360,393]);
N_grad_iso_lat = nan([180,360,393]);
for p = 1:393
    for lat = 2:179
        for lon = 2:359  
            dN_lon = N(lat,(lon+1),p) - N(lat,(lon-1),p);
            dN_lat = N((lat+1),lon,p) - N((lat-1),lon,p);
            lon_dist = (sw_dist([glob_lat(lat) glob_lat(lat)],[glob_lon(lon+1) glob_lon(lon-1)],'km'))*1000; %horizontal dist
            lat_dist = (sw_dist([glob_lat(lat+1) glob_lat(lat-1)],[glob_lon(lon) glob_lon(lon)], 'km'))*1000;
            z_change_lon = Z(lat,(lon+1),p) - Z(lat,(lon-1),p); %dist along isopycnal
            z_change_lat = Z((lat+1),lon,p) - Z((lat-1),lon,p);
            iso_dist_lon = sqrt((lon_dist).^2 + (z_change_lon).^2);
            iso_dist_lat = sqrt((lat_dist).^2 + (z_change_lat).^2);
            N_grad_iso_lon(lat,lon,p) = dN_lon / iso_dist_lon; %gradients
            N_grad_iso_lat(lat,lon,p) = dN_lat / iso_dist_lat;            
        end
    end
end
%save('N_grad_iso_lon.mat','N_grad_iso_lon')
%save('N_grad_iso_lat.mat','N_grad_iso_lat')

%convert units
K_iso_yr = K_iso*60*60*24*365; %(m2 s-1 -> m2 yr-1)

%calculate isopycnal N flux, F_iso, on [20.98:0.02:28.82] (mol m2 yr-1)
F_iso_lat = -K_iso_yr.*N_grad_iso_lat;     
F_iso_lon = -K_iso_yr.*N_grad_iso_lon;     
F_iso = sqrt((F_iso_lon).^2 + (F_iso_lat).^2); 
%save('F_iso_lat.mat','F_iso_lat')
%save('F_iso_lon.mat','F_iso_lon')

%% Calculate diapycnal diffusive N flux 

%calculate diapycnal N gradient, dN/dz*, on [20.99:0.02:28.81] (mol m-4)
N_grad_dia = nan([length(glob_lat),length(glob_lon),392]);
for lat =1:length(glob_lat);
    for lon = 1:length(glob_lon);
        for p = 1:392
            dN = N(lat,lon,p)- N(lat,lon,p+1);
            h = Z(lat, lon, p) - Z(lat, lon, p+1);
            N_grad_dia(lat, lon, p) = dN / h;
        end
    end
end
%save('N_grad_dia.mat','N_grad_dia')

%convert units 
K_dia_yr = K_dia_offset*60*60*24*365; %(m2s-1 -> m2yr-1)

%calculate diapycnal diffusive N flux, F_dia_diff on [20.99:0.02:28.81] %(mol m2 yr-1)
F_dia_diff = -K_dia_yr.*N_grad_dia;
%save('F_dia_diff.mat','F_dia_diff')

%% Calculate diapycnal advective N flux 

%calculate diapycnal density gradient, dp/dz, on [20.98:0.02:28.82]
%which needs z on [20.97:0.02:28.83] (kg m-4)
dens_grad_dia = nan(length(glob_lat),length(glob_lon),393);
for lat =1:length(glob_lat);
    for lon =1:length(glob_lon);
        for p = 1:392;
            if sum(isnan(Z_offset_2(lat, lon, :)))>392;
                dens_grad_dia(lat, lon, p) = NaN;
            else
                dp = -0.02;
                dist = Z_offset_2(lat, lon, p) - Z_offset_2(lat, lon, p+1); 
                dens_grad_dia(lat, lon, p) = dp/dist;
            end
        end
    end
end

%calculate diapycnal diffusive density flux on [20.98:0.02:28.82] (kg m-2 s-1)
F_dia_dens = - K_dia.*dens_grad_dia;
%save('F_dia_dens.mat','F_dia_dens')
%save('dens_grad_dia.mat','dens_grad_dia')

%calculate diapycnal velocity, w*,on [20.99:0.02:28.81] (m s-1)
%from contrasts of the diapycnal diffusive density flux 
w = nan(length(glob_lat),length(glob_lon),392);
for lat =1:length(glob_lat);
    for lon = 1:length(glob_lon);
        for p = 1:392;
            if sum(isnan(K_dia_offset(lat, lon, :)))>392;
                w(lat, lon, p) = NaN;
            else
                F_dia_dens_grad = F_dia_dens(lat, lon, p) - F_dia_dens(lat, lon, p+1);
                dp = -0.02;
                w(lat, lon, p) = -F_dia_dens_grad./dp; %negative because upwards is positive (not upwards to high density)
            end
        end
    end
end
%save('w.mat','w')

%calculate diapcynal diffusive N flux, F_dia_adv, on [20.99:0.02:28.81]
%(mol m-2 yr-1)
F_dia_adv = w*60*60*24*365.*N_offset;
%save('F_dia_adv.mat','F_dia_adv')

%% Calculate convergence of the isopycnal diffusive N flux

%convergence of Fiso on [20.98:0.02:28.82]
conv_F_iso = nan(length(glob_lat), length(glob_lon), 393);
for lat =2:length(glob_lat)-1;
    for lon = 2:length(glob_lon)-4;
        for p = 1:393;
            dF_lon = F_iso_lon(lat,(lon+1),p) - F_iso_lon(lat,(lon-1),p);    
            dF_lat = F_iso_lat((lat+1),lon,p) - F_iso_lat((lat-1),lon,p); 
            lon_dist = (sw_dist([glob_lat(lat) glob_lat(lat)],[glob_lon(lon+1) glob_lon(lon-1)],'km'))*1000; %horizontal dist
            lat_dist = (sw_dist([glob_lat(lat+1) glob_lat(lat-1)],[glob_lon(lon) glob_lon(lon)], 'km'))*1000;
            z_change_lon = Z(lat,(lon+1),p) - Z(lat,(lon-1),p);
            z_change_lat = Z((lat+1),lon,p) - Z((lat-1),lon,p);
            iso_dist_lon = sqrt((lon_dist)^2 + (z_change_lon)^2);          
            iso_dist_lat = sqrt((lat_dist)^2 + (z_change_lat)^2);         
            conv_F_iso_lon = dF_lon / iso_dist_lon;    
            conv_F_iso_lat = dF_lat / iso_dist_lat;    
            conv_F_iso(lat, lon, p) = -(conv_F_iso_lon + conv_F_iso_lat);
        end
    end
end
%save('conv_F_iso.mat','conv_F_iso')

%% Calculate convergence of the diapycnal diffusive N flux (don't need?)

%convergence of F_dia_diff excludlatng top density layer on [25.00:0.02:28.80]
conv_F_dia_diff_main = nan(length(glob_lat), length(glob_lon), 392);
for lat =1:length(glob_lat);
    for lon = 1:length(glob_lon);
        for p = 1:391
            dF = F_dia_diff(lat,lon,p+1) - F_dia_diff(lat,lon,p);   
            h = -(Z(lat, lon, p+1) -  Z(lat, lon, p)); 
            conv_F_dia_diff_main(lat,lon,p) = dF/h; 
        end
    end
end

%convergence of F_dia_diff for top density layer on [25.00:0.02:28.80]
conv_F_dia_diff_top = nan(length(glob_lat),length(glob_lon));
for lat =1:length(glob_lat);
    for lon = 1:length(glob_lon);
        if sum(isnan(F_dia_diff(lat, lon, :)))>0; 
            dF_top = F_dia_diff(lat, lon, 1);
            h_top = -Z(lat, lon, 1);
            conv_F_dia_diff_top(lat,lon) = dF_top / h_top;
        else
            conv_F_dia_diff_top(lat,lon) = NaN;
        end
    end
end

%total convergence of F_dia_diff on [25.00:0.02:28.80]
conv_F_dia_diff(:,:,1) = conv_F_dia_diff_top;
conv_F_dia_diff(:,:,[2:393]) = conv_F_dia_diff_main;
%save('conv_F_dia_diff.mat','conv_F_dia_diff')

% Calculate N supply by mesoscale eddies (isopycnal diffusion)

% WML
wml_F_iso = nan(length(glob_lat), length(glob_lon)); 
for lat = [find(glob_lat == -40.5):find(glob_lat == 40.5)]
    for lon = 1:length(glob_lon)
        if sum(isnan(Z(lat, lon, :))) > 392
            wml_F_iso(lat, lon) = NaN;
        else
            iso_in_layer = find(-Z(lat, lon, :) <= w_mld_z(lat, lon));
            
            if iso_in_layer > 0
                iso_mid_except_top = squeeze((Z(lat, lon, iso_in_layer(2:end) - 1) + Z(lat, lon, iso_in_layer(2:end))) / 2);
                iso_mid_top = (Z(lat, lon, iso_in_layer(1)) / 2);    
                iso_mid = [iso_mid_top; iso_mid_except_top];               
                iso_mid_offset = [0; iso_mid(1:end-1)];
                iso_thickness = iso_mid_offset - iso_mid;
                layer_supply = squeeze(conv_F_iso(lat, lon, iso_in_layer)) .* iso_thickness;
                wml_F_iso(lat, lon) = sum(layer_supply, 'omitnan');
            else
                dens = find(~isnan(Z(lat, lon, :)));
                wml_F_iso(lat, lon) = conv_F_iso(lat, lon, dens(1)) * -w_mld_z(lat, lon);
            end
        end
    end
end

% Thermocline
therm_F_iso = nan(length(glob_lat), length(glob_lon)); 
for lat = [find(glob_lat == -40.5):find(glob_lat == 40.5)]
    for lon = 1:length(glob_lon)
        if sum(isnan(Z(lat, lon, :))) > 392
            therm_F_iso(lat, lon) = NaN;
        else
            iso_in_layer = find(-Z(lat, lon, [1:find(dens_grid == 27)]) > w_mld_z(lat, lon));
            
            if iso_in_layer > 0
                iso_mid_except_top = squeeze((Z(lat, lon, iso_in_layer(2:end) - 1) + Z(lat, lon, iso_in_layer(2:end))) / 2);
                iso_mid_top = (Z(lat, lon, iso_in_layer(1)) / 2);    
                iso_mid = [iso_mid_top; iso_mid_except_top];               
                iso_mid_offset = [0; iso_mid(1:end-1)];
                iso_thickness = iso_mid_offset - iso_mid;
                layer_supply = squeeze(conv_F_iso(lat, lon, iso_in_layer)) .* iso_thickness;
                therm_F_iso(lat, lon) = sum(layer_supply, 'omitnan');
            else
                therm_F_iso(lat, lon) = NaN;
            end
        end
    end
end

% Calculate N supply by microscale turbulence (diapycnal diffusion)

% WML
wml_F_dia_diff = nan(length(glob_lat), length(glob_lon)); 
for lat = [find(glob_lat == -40.5):find(glob_lat == 40.5)]
    for lon = 1:length(glob_lon)
        if sum(isnan(Z(lat, lon, :))) > 392
            wml_F_dia_diff(lat, lon) = NaN;
        elseif sum(isnan(K_dia(lat, lon, :))) > 392
            wml_F_dia_diff(lat, lon) = NaN;
        else
            densities = find(-Z(lat, lon, :) <= w_mld_z(lat, lon));
            if densities > 0
                wml_F_dia_diff(lat, lon) = F_dia_diff(lat, lon, densities(end));
            else
                densities_no_nans = find(~isnan(Z(lat, lon, :)));
                wml_F_dia_diff(lat, lon) = F_dia_diff(lat, lon, densities_no_nans(1));
            end
        end
    end
end

% Thermocline
therm_F_dia_diff = nan(length(glob_lat), length(glob_lon)); 
for lat = [find(glob_lat == -40.5):find(glob_lat == 40.5)]
    for lon = 1:length(glob_lon)
        if sum(isnan(Z(lat, lon, :))) > 392
            wml_F_dia_diff(lat, lon) = NaN;
        elseif sum(isnan(K_dia(lat, lon, :))) > 392
            therm_F_dia_diff(lat, lon) = NaN;
        else
            densities = find(-Z(lat, lon, [1:find(dens_grid == 27)]) <= w_mld_z(lat, lon));
            fluxes_no_nans = find(~isnan(F_dia_diff(lat, lon, :)));
            A = ismember(densities, fluxes_no_nans);
            densities = densities(A);
            if densities > 0
                therm_F_dia_diff(lat, lon) = -F_dia_diff(lat, lon, densities(1)) + F_dia_diff(lat, lon, densities(end));
            else
                therm_F_dia_diff(lat, lon) = NaN;
            end
        end
    end
end

% Calculate N supply by microscale turbulence (diapcynal advection)

% WML
wml_F_dia_adv = nan(length(glob_lat), length(glob_lon)); 
for lat = [find(glob_lat == -40.5):find(glob_lat == 40.5)]
    for lon = 1:length(glob_lon)
        if sum(isnan(Z(lat, lon, :))) > 392
            wml_F_dia_adv(lat, lon) = NaN;
        elseif sum(isnan(K_dia(lat, lon, :))) > 392
            wml_F_dia_adv(lat, lon) = NaN;
        else
            densities = find(-Z(lat, lon, :) <= w_mld_z(lat, lon));
            if densities > 0
                wml_F_dia_adv(lat, lon) = F_dia_adv(lat, lon, densities(end));
            else
                densities_no_nans = find(~isnan(Z(lat, lon, :)));
                wml_F_dia_adv(lat, lon) = F_dia_adv(lat, lon, densities_no_nans(1));
            end
        end
    end
end

% Thermocline
therm_F_dia_adv = nan(length(glob_lat), length(glob_lon)); 
for lat = [find(glob_lat == -40.5):find(glob_lat == 40.5)]
    for lon = 1:length(glob_lon)
        if sum(isnan(Z(lat, lon, :))) > 392
            therm_F_dia_adv(lat, lon) = NaN;
        elseif sum(isnan(K_dia(lat, lon, :))) > 392
            therm_F_dia_adv(lat, lon) = NaN;
        else
            densities = find(-Z(lat, lon, [1:find(dens_grid == 27)]) <= w_mld_z(lat, lon));
            fluxes_no_nans = find(~isnan(F_dia_adv(lat, lon, :)));
            A = ismember(densities, fluxes_no_nans);
            densities = densities(A);
            if densities > 0
                therm_F_dia_adv(lat, lon) = -F_dia_adv(lat, lon, densities(1)) + F_dia_adv(lat, lon, densities(end));
            else
                therm_F_dia_adv(lat, lon) = NaN;
            end
        end
    end
end

%% Save N supply estimates

% save('wml_F_iso.mat','wml_F_iso')
% save('therm_F_iso.mat','therm_F_iso')
% save('wml_F_dia_diff.mat','wml_F_dia_diff')
% save('therm_F_dia_diff.mat','therm_F_dia_diff')
% save('wml_F_dia_adv.mat','wml_F_dia_adv')
% save('therm_F_dia_adv.mat','therm_F_dia_adv')

