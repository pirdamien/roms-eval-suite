
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc'
rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/mean_2000_2005/'
rep_out = './Fig/'
file_avg = [rep_in 'pacmed_avg.nc'];
file_bgc = [rep_in 'pacmed_bgc_avg.nc'];

theta_s = 6.0; theta_b = 6.0;
hc = 250; sc_type = 'new2012';
NZ = 100 ; epsilon = 1e-3 ;

load('../colormap_IsleOfDogs.dat');
colormap_mine=colormap_IsleOfDogs(:,2:4) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon = ncread(grid_file,'lon_rho') ;
lat = ncread(grid_file,'lat_rho') ;
h   = ncread(grid_file,'h') ;
mask= ncread(grid_file,'mask_rho') ;
pm  = ncread(grid_file,'pm') ;
pn  = ncread(grid_file,'pn') ;
lon(lon<0) = lon(lon<0)+360 ;
[NX,NY]=size(h);
for z=1:NZ
mask3d(:,:,z) = mask ;
end

lon_rshp = [140:1:270] ; 
lat_rshp = [-30:1:30] ;
[lat_rshp lon_rshp] = meshgrid(lat_rshp,lon_rshp) ;
dep_rshp = [[0:1:10] [12:2:100] [105:5:300] [310:10:500]] ; 
%dep_rshp = [0:25:500] ;

    zeta = squeeze(mean(ncread(file_avg,'zeta'),3)) ;
    temp = squeeze(mean(ncread(file_avg,'temp'),4)) ;
    salt = squeeze(mean(ncread(file_avg,'salt'),4)) ;
    nitr = squeeze(mean(ncread(file_bgc, 'NO3'),4)) ;
    phos = squeeze(mean(ncread(file_bgc, 'PO4'),4)) ;
    oxyg = squeeze(mean(ncread(file_bgc,  'O2'),4)) ;
    silc = squeeze(mean(ncread(file_bgc,'SiO3'),4)) ;
    ntro = squeeze(mean(ncread(file_bgc, 'N2O'),4)) ;
    iron = squeeze(mean(ncread(file_bgc,  'Fe'),4)) ;

    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);
    temp(mask3d==0) = NaN ; salt(mask3d==0) = NaN ;
    nitr(mask3d==0) = NaN ; phos(mask3d==0) = NaN ;
    oxyg(mask3d==0) = NaN ; silc(mask3d==0) = NaN ;
    ntro(mask3d==0) = NaN ; iron(mask3d==0) = NaN ;
    for z=1:length(dep_rshp)
        disp(['---> ' num2str(z) '/' num2str(length(dep_rshp))])
        var=permute(vinterp(permute(temp,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        temp_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;             
        var=permute(vinterp(permute(salt,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        salt_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(nitr,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        nitr_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(phos,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        phos_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(oxyg,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        oxyg_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(silc,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        silc_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(ntro,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        ntro_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(iron,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        iron_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
    end

temp_180160 = squeeze(nanmean(temp_rshp(40: 61,:,:),1)) ; temp_160120 = squeeze(nanmean(temp_rshp(61:101,:,:),1)) ;
salt_180160 = squeeze(nanmean(salt_rshp(40: 61,:,:),1)) ; salt_160120 = squeeze(nanmean(salt_rshp(61:101,:,:),1)) ;
nitr_180160 = squeeze(nanmean(nitr_rshp(40: 61,:,:),1)) ; nitr_160120 = squeeze(nanmean(nitr_rshp(61:101,:,:),1)) ;
phos_180160 = squeeze(nanmean(phos_rshp(40: 61,:,:),1)) ; phos_160120 = squeeze(nanmean(phos_rshp(61:101,:,:),1)) ;
oxyg_180160 = squeeze(nanmean(oxyg_rshp(40: 61,:,:),1)) ; oxyg_160120 = squeeze(nanmean(oxyg_rshp(61:101,:,:),1)) ;
silc_180160 = squeeze(nanmean(silc_rshp(40: 61,:,:),1)) ; silc_160120 = squeeze(nanmean(silc_rshp(61:101,:,:),1)) ;
ntro_180160 = squeeze(nanmean(ntro_rshp(40: 61,:,:),1)) ; ntro_160120 = squeeze(nanmean(ntro_rshp(61:101,:,:),1)) ;
iron_180160 = squeeze(nanmean(iron_rshp(40: 61,:,:),1)) ; iron_160120 = squeeze(nanmean(iron_rshp(61:101,:,:),1)) ;
for i=1:size(lat_rshp,2)
    for j=1:length(dep_rshp)
        dep2d(i,j) = -dep_rshp(j) ;
        lat2d(i,j) = lat_rshp(1,i);
    end
end

temp_eq = squeeze(temp_rshp(:,31,:)) ; 
salt_eq = squeeze(salt_rshp(:,31,:)) ; 
nitr_eq = squeeze(nitr_rshp(:,31,:)) ; 
phos_eq = squeeze(phos_rshp(:,31,:)) ; 
oxyg_eq = squeeze(oxyg_rshp(:,31,:)) ; 
silc_eq = squeeze(silc_rshp(:,31,:)) ; 
ntro_eq = squeeze(ntro_rshp(:,31,:)) ; 
iron_eq = squeeze(iron_rshp(:,31,:)) ; 
for i=1:size(lon_rshp,1)
    for j=1:length(dep_rshp)
        depeq(i,j) = -dep_rshp(j) ;
        loneq(i,j) = lon_rshp(i,1);
    end
end


%%%%%%%%%%%%%% TEMPERATURE & SALINITY %%%%%%%%%%%%%%%

file_t = '/data/project1/data/WOA18/temperature/woa18_decav_t00_04.nc' ;
file_s = '/data/project1/data/WOA18/salinity/woa18_decav_s00_04.nc' ;
lon_obs = ncread(file_t,'lon');
lat_obs = ncread(file_t,'lat');
dep_obs = ncread(file_t,'depth',1,37);
temp_obs = squeeze(ncread(file_t,'t_an',[1 1 1 1],[inf inf 37 1]));
salt_obs = squeeze(ncread(file_s,'s_an',[1 1 1 1],[inf inf 37 1]));

temp_180160_obs = squeeze(nanmean(temp_obs( 1: 80,:,:),1)) ; 
temp_160120_obs = squeeze(nanmean(temp_obs(81:240,:,:),1)) ;
salt_180160_obs = squeeze(nanmean(salt_obs( 1: 80,:,:),1)) ;
salt_160120_obs = squeeze(nanmean(salt_obs(81:240,:,:),1)) ;
for i=1:length(lat_obs)
    for j=1:length(dep_obs)
        dep2d_obs(i,j) = -dep_obs(j) ;
        lat2d_obs(i,j) = lat_obs(i);
    end
end
temp_eq_obs = squeeze(nanmean(temp_obs(:,360:361,:),2)) ;
salt_eq_obs = squeeze(nanmean(salt_obs(:,360:361,:),2)) ;
lon_obs = [lon_obs(721:end)' lon_obs(1:720)'+360]';
temp_eq_obs = [temp_eq_obs(721:end,:)' temp_eq_obs(1:720,:)']'; 
salt_eq_obs = [salt_eq_obs(721:end,:)' salt_eq_obs(1:720,:)']';
for i=1:length(lon_obs)
    for j=1:length(dep_obs)
        depeq_obs(i,j) = -dep_obs(j) ;
        loneq_obs(i,j) = lon_obs(i);
    end
end


figure
%pcolor(lat2d,dep2d,temp_180160) ; shading flat ; 
contourf(lat2d,dep2d,temp_180160,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([6 30]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_180160,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Temperature on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Temp_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,temp_160120) ; shading flat ; 
contourf(lat2d,dep2d,temp_160120,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([6 30]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_160120,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Temperature on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Temp_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,temp_eq) ; shading flat ; 
contourf(loneq,depeq,temp_eq,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([6 30]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,temp_eq,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Temperature along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Temp_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_obs,dep2d_obs,temp_180160_obs) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,temp_180160_obs,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ; 
caxis([6 30]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_180160_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ; 
title('Temperature on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Temp_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_obs,dep2d_obs,temp_160120_obs) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,temp_160120_obs,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([6 30]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_160120_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;                    
title('Temperature on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Temp_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_obs,depeq_obs,temp_eq_obs) ; shading flat ; 
contourf(loneq_obs,depeq_obs,temp_eq_obs,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([6 30]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,temp_eq_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Temperature along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Temp_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close



figure
%pcolor(lat2d,dep2d,salt_180160) ; shading flat ; 
contourf(lat2d,dep2d,salt_180160,[-0.025:0.05:40.025],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([33.8 36.3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_180160,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Salinity on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Salt_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,salt_160120) ; shading flat ; 
contourf(lat2d,dep2d,salt_160120,[-0.025:0.05:40.025],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([33.8 36.3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_160120,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Salinity on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Salt_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,salt_eqs) ; shading flat ; 
contourf(loneq,depeq,salt_eq,[-0.025:0.05:40.025],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([34.4 35.4]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,salt_eq,[33:0.2:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Salinity along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Salt_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close

figure
%pcolor(lat2d_obs,dep2d_obs,salt_180160_obs) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,salt_180160_obs,[-0.025:0.05:40.025],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([33.8 36.3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_180160_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Salinity on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Salt_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_obs,dep2d_obs,salt_160120_obs) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,salt_160120_obs,[-0.025:0.05:40.025],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([33.8 36.3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_160120_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Salinity on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Salt_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_obs,depeq_obs,salt_eq_obs) ; shading flat ; 
contourf(loneq_obs,depeq_obs,salt_eq_obs,[-0.025:0.05:40.025],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([34.4 35.4]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,salt_eq_obs,[33:0.2:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Salinity along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Salt_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close

%%%%%%%%%%%%%%%%%%%%%%%%% NO3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dep2d_vobs lat2d_vobs depeq_vobs loneq_vobs
file_t = '/data/project7/pdamien/DATA/BGC/WOA18/nitr/woa18_all_n00_01.nc' ;
lon_vobs = ncread(file_t,'lon');
lat_vobs = ncread(file_t,'lat');
dep_vobs = ncread(file_t,'depth',1,37);
var_obs = squeeze(ncread(file_t,'n_an',[1 1 1 1],[inf inf 37 1]));

var_180160_obs = squeeze(nanmean(var_obs( 1:20,:,:),1)) ;
var_160120_obs = squeeze(nanmean(var_obs(21:60,:,:),1)) ;
for i=1:length(lat_vobs)
    for j=1:length(dep_vobs)
        dep2d_vobs(i,j) = -dep_vobs(j) ;
        lat2d_vobs(i,j) = lat_vobs(i);
    end
end
var_eq_obs = squeeze(nanmean(var_obs(:,90:91,:),2)) ;
lon_vobs = [lon_vobs(181:end)' lon_vobs(1:180)'+360]';
var_eq_obs = [var_eq_obs(181:end,:)' var_eq_obs(1:180,:)']';
for i=1:length(lon_vobs)
    for j=1:length(dep_vobs)
        depeq_vobs(i,j) = -dep_vobs(j) ;
        loneq_vobs(i,j) = lon_vobs(i);
    end
end

figure
%pcolor(lat2d,dep2d,nitr_180160) ; shading flat ; 
contourf(lat2d,dep2d,nitr_180160,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 38]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_180160,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Nitrate on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Nitr_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,nitr_160120) ; shading flat ; 
contourf(lat2d,dep2d,nitr_160120,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 38]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_160120,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Nitrate on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Nitr_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,nitr_eq) ; shading flat ; 
contourf(loneq,depeq,nitr_eq,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 38]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,temp_eq,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Nitrate along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Nitr_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close

figure
%pcolor(lat2d_vobs,dep2d_vobs,var_180160_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_180160_obs,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 38]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_180160_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Nitrate on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Nitr_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_vobs,dep2d_vobs,var_160120_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_160120_obs,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 38]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_160120_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Nitrate on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Nitr_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_vobs,depeq_vobs,var_eq_obs) ; shading flat ; 
contourf(loneq_vobs,depeq_vobs,var_eq_obs,[-0.25:0.5:40.25],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 38]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,temp_eq_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Nitrate along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Nitr_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close

%%%%%%%%%%%%%%%%%%%%%%%%% PO4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dep2d_vobs lat2d_vobs depeq_vobs loneq_vobs
file_t = '/data/project7/pdamien/DATA/BGC/WOA18/phos/woa18_all_p00_01.nc' ;
lon_vobs = ncread(file_t,'lon');
lat_vobs = ncread(file_t,'lat');
dep_vobs = ncread(file_t,'depth',1,37);
var_obs = squeeze(ncread(file_t,'p_an',[1 1 1 1],[inf inf 37 1]));

var_180160_obs = squeeze(nanmean(var_obs( 1:20,:,:),1)) ;
var_160120_obs = squeeze(nanmean(var_obs(21:60,:,:),1)) ;
for i=1:length(lat_vobs)
    for j=1:length(dep_vobs)
        dep2d_vobs(i,j) = -dep_vobs(j) ;
        lat2d_vobs(i,j) = lat_vobs(i);
    end
end
var_eq_obs = squeeze(nanmean(var_obs(:,90:91,:),2)) ;
lon_vobs = [lon_vobs(181:end)' lon_vobs(1:180)'+360]';
var_eq_obs = [var_eq_obs(181:end,:)' var_eq_obs(1:180,:)']';
for i=1:length(lon_vobs)
    for j=1:length(dep_vobs)
        depeq_vobs(i,j) = -dep_vobs(j) ;
        loneq_vobs(i,j) = lon_vobs(i);
    end
end


figure
%pcolor(lat2d,dep2d,phos_180160) ; shading flat ; 
contourf(lat2d,dep2d,phos_180160,[-0.05:0.1:5.05],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_180160,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Phosphate on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Phos_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,phos_160120) ; shading flat ; 
contourf(lat2d,dep2d,phos_160120,[-0.05:0.1:5.05],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_160120,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Phosphate on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Phos_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,phos_eq) ; shading flat ; 
contourf(loneq,depeq,phos_eq,[-0.05:0.1:5.05],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 3]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,temp_eq,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Phosphate along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Phos_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close

figure
%pcolor(lat2d_vobs,dep2d_vobs,var_180160_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_180160_obs,[-0.05:0.1:5.05],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_180160_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Phosphate on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Phos_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_vobs,dep2d_vobs,var_160120_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_160120_obs,[-0.05:0.1:5.05],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 3]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_160120_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Phosphate on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Phos_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_vobs,depeq_vobs,var_eq_obs) ; shading flat ; 
contourf(loneq_vobs,depeq_vobs,var_eq_obs,[-0.05:0.1:5.05],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 3]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,temp_eq_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Phosphate along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Phos_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close



%%%%%%%%%%%%%%%%%%%%%%%%% O2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dep2d_vobs lat2d_vobs depeq_vobs loneq_vobs
file_t = '/data/project7/pdamien/DATA/BGC/WOA18/oxyg/woa18_all_o00_01.nc' ;
lon_vobs = ncread(file_t,'lon');
lat_vobs = ncread(file_t,'lat');
dep_vobs = ncread(file_t,'depth',1,37);
var_obs = squeeze(ncread(file_t,'o_an',[1 1 1 1],[inf inf 37 1]));

var_180160_obs = squeeze(nanmean(var_obs( 1:20,:,:),1)) ;
var_160120_obs = squeeze(nanmean(var_obs(21:60,:,:),1)) ;
for i=1:length(lat_vobs)
    for j=1:length(dep_vobs)
        dep2d_vobs(i,j) = -dep_vobs(j) ;
        lat2d_vobs(i,j) = lat_vobs(i);
    end
end
var_eq_obs = squeeze(nanmean(var_obs(:,90:91,:),2)) ;
lon_vobs = [lon_vobs(181:end)' lon_vobs(1:180)'+360]';
var_eq_obs = [var_eq_obs(181:end,:)' var_eq_obs(1:180,:)']';
for i=1:length(lon_vobs)
    for j=1:length(dep_vobs)
        depeq_vobs(i,j) = -dep_vobs(j) ;
        loneq_vobs(i,j) = lon_vobs(i);
    end
end

figure
%pcolor(lat2d,dep2d,oxyg_180160) ; shading flat ; 
contourf(lat2d,dep2d,oxyg_180160,[-2.5:5:302.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 250]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_180160,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Oxygen on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Oxyg_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,oxyg_160120) ; shading flat ; 
contourf(lat2d,dep2d,oxyg_160120,[-2.5:5:302.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 250]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_160120,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Oxygen on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Oxyg_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,oxyg_eq) ; shading flat ; 
contourf(loneq,depeq,oxyg_eq,[-2.5:5:302.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 250]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,salt_eq,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Oxygen along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Oxyg_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_vobs,dep2d_vobs,var_180160_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_180160_obs,[-2.5:5:302.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 250]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_180160_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Oxygen on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Oxyg_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_vobs,dep2d_vobs,var_160120_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_160120_obs,[-2.5:5:302.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 250]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_160120_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Oxygen on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Oxyg_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_vobs,depeq_vobs,var_eq_obs) ; shading flat ; 
contourf(loneq_vobs,depeq_vobs,var_eq_obs,[-2.5:5:302.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 250]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,salt_eq_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Oxygen along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Oxyg_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


%%%%%%%%%%%%%%%%%%%%%%%%% SILICATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dep2d_vobs lat2d_vobs depeq_vobs loneq_vobs
file_t = '/data/project7/pdamien/DATA/BGC/WOA18/silc/woa18_all_i00_01.nc' ;
lon_vobs = ncread(file_t,'lon');
lat_vobs = ncread(file_t,'lat');
dep_vobs = ncread(file_t,'depth',1,37);
var_obs = squeeze(ncread(file_t,'i_an',[1 1 1 1],[inf inf 37 1]));

var_180160_obs = squeeze(nanmean(var_obs( 1:20,:,:),1)) ;
var_160120_obs = squeeze(nanmean(var_obs(21:60,:,:),1)) ;
for i=1:length(lat_vobs)
    for j=1:length(dep_vobs)
        dep2d_vobs(i,j) = -dep_vobs(j) ;
        lat2d_vobs(i,j) = lat_vobs(i);
    end
end
var_eq_obs = squeeze(nanmean(var_obs(:,90:91,:),2)) ;
lon_vobs = [lon_vobs(181:end)' lon_vobs(1:180)'+360]';
var_eq_obs = [var_eq_obs(181:end,:)' var_eq_obs(1:180,:)']';
for i=1:length(lon_vobs)
    for j=1:length(dep_vobs)
        depeq_vobs(i,j) = -dep_vobs(j) ;
        loneq_vobs(i,j) = lon_vobs(i);
    end
end

figure
%pcolor(lat2d,dep2d,silc_180160) ; shading flat ; 
contourf(lat2d,dep2d,silc_180160,[-0.5:1:70.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 50]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_180160,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Silicate on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Silc_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,silc_160120) ; shading flat ; 
contourf(lat2d,dep2d,silc_160120,[-0.5:1:70.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 50]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_160120,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Silicate on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Silc_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,silc_eq) ; shading flat ; 
contourf(loneq,depeq,silc_eq,[-0.5:1:70.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 50]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,salt_eq,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Silicate along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Silc_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_vobs,dep2d_vobs,var_180160_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_180160_obs,[-0.5:1:70.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 50]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_180160_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Silicate on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Silc_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_vobs,dep2d_vobs,var_160120_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_160120_obs,[-0.5:1:70.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 50]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_160120_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Silicate on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Silc_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_vobs,depeq_vobs,var_eq_obs) ; shading flat ; 
contourf(loneq_vobs,depeq_vobs,var_eq_obs,[-0.5:1:70.5],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0 50]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,salt_eq_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Silicate along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Silc_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


%%%%%%%%%%%%%%%%%%%%%%%%% Nitrous Oxide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dep2d_vobs lat2d_vobs depeq_vobs loneq_vobs
file_t = '/data/project7/pdamien/DATA/BGC/Ncycle/n2ofromnn.nc' ;
lon_vobs = ncread(file_t,'lon');
lat_vobs = ncread(file_t,'lat');
dep_vobs = ncread(file_t,'DEPTH',1,37);
var_obs = squeeze(ncread(file_t,'N2O',[1 1 1 1],[inf inf 37 inf]));
var_obs = squeeze(mean(var_obs,4)) ; 

var_180160_obs = squeeze(nanmean(var_obs( 1:20,:,:),1)) ;
var_160120_obs = squeeze(nanmean(var_obs(21:60,:,:),1)) ;
for i=1:length(lat_vobs)
    for j=1:length(dep_vobs)
        dep2d_vobs(i,j) = -dep_vobs(j) ;
        lat2d_vobs(i,j) = lat_vobs(i);
    end
end
var_eq_obs = squeeze(nanmean(var_obs(:,90:91,:),2)) ;
lon_vobs = [lon_vobs(181:end)' lon_vobs(1:180)'+360]';
var_eq_obs = [var_eq_obs(181:end,:)' var_eq_obs(1:180,:)']';
for i=1:length(lon_vobs)
    for j=1:length(dep_vobs)
        depeq_vobs(i,j) = -dep_vobs(j) ;
        loneq_vobs(i,j) = lon_vobs(i);
    end
end

figure
%pcolor(lat2d,dep2d,ntro_180160) ; shading flat ; 
contourf(lat2d,dep2d,ntro_180160,[-0.0005:0.001:0.0805],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.005 0.045]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_180160,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('N2O on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'N2O_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,ntro_160120) ; shading flat ; 
contourf(lat2d,dep2d,ntro_160120,[-0.0005:0.001:0.0805],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.005 0.045]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,temp_160120,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('N2O on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'N2O_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,ntro_eq) ; shading flat ; 
contourf(loneq,depeq,ntro_eq,[-0.0005:0.001:0.0805],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.005 0.045]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,temp_eq,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('N2O along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'N2O_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_vobs,dep2d_vobs,var_180160_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_180160_obs,[-0.0005:0.001:0.0805],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.005 0.045]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_180160_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('N2O on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'N2O_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_vobs,dep2d_vobs,var_160120_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_160120_obs,[-0.0005:0.001:0.0805],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.005 0.045]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,temp_160120_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('N2O on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'N2O_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_vobs,depeq_vobs,var_eq_obs) ; shading flat ; 
contourf(loneq_vobs,depeq_vobs,var_eq_obs,[-0.0005:0.001:0.0805],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.005 0.045]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,temp_eq_obs,[0:5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('N2O along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'N2O_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


%%%%%%%%%%%%%%% IRON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear dep2d_vobs lat2d_vobs depeq_vobs loneq_vobs
file_t = '/data/project7/pdamien/DATA/BGC/MappedIron/Monthly_dFe.nc' ;
lon_vobs = ncread(file_t,'Longitude');
lat_vobs = ncread(file_t,'Latitude');
dep_vobs = ncread(file_t,'Depth',1,20);
var_obs = squeeze(ncread(file_t,'dFe_RF',[1 1 1 13],[inf inf 20 1]))*0.001;
for z=1:length(dep_vobs)
    var_obs(:,:,z) = inpaint_nans(squeeze(var_obs(:,:,z)),2) ;
end

var_180160_obs = squeeze(nanmean(var_obs( 1:19,:,:),1)) ;
var_160120_obs = squeeze(nanmean(var_obs(20:59,:,:),1)) ;
for i=1:length(lat_vobs)
    for j=1:length(dep_vobs)
        dep2d_vobs(i,j) = -dep_vobs(j) ;
        lat2d_vobs(i,j) = lat_vobs(i);
    end
end
var_eq_obs = squeeze(var_obs(:,73,:)) ;
lon_vobs = [lon_vobs(179:end)' lon_vobs(1:187)'+360]';
var_eq_obs = [var_eq_obs(179:end,:)' var_eq_obs(1:187,:)']';
for i=1:length(lon_vobs)
    for j=1:length(dep_vobs)
        depeq_vobs(i,j) = -dep_vobs(j) ;
        loneq_vobs(i,j) = lon_vobs(i);
    end
end

figure
%pcolor(lat2d,dep2d,iron_180160) ; shading flat ; 
contourf(lat2d,dep2d,iron_180160,[-0.05e-4:0.1e-4:10.05e-4],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.5e-4 5.5e-4]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_180160,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Fe on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'IRON_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,iron_160120) ; shading flat ; 
contourf(lat2d,dep2d,iron_160120,[-0.05e-4:0.1e-4:10.05e-4],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.5e-4 5.5e-4]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,salt_160120,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Fe on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'IRON_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,iron_eq) ; shading flat ; 
contourf(loneq,depeq,iron_eq,[-0.05e-4:0.1e-4:10.05e-4],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.5e-4 5.5e-4]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,salt_eq,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Fe along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'IRON_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_vobs,dep2d_vobs,var_180160_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_180160_obs,[-0.05e-4:0.1e-4:10.05e-4],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.5e-4 5.5e-4]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_180160_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Fe on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'IRON_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_vobs,dep2d_vobs,var_160120_obs) ; shading flat ; 
contourf(lat2d_vobs,dep2d_vobs,var_160120_obs,[-0.05e-4:0.1e-4:10.05e-4],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.5e-4 5.5e-4]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,salt_160120_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Fe on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'IRON_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_vobs,depeq_vobs,var_eq_obs) ; shading flat ; 
contourf(loneq_vobs,depeq_vobs,var_eq_obs,[-0.05e-4:0.1e-4:10.05e-4],'edgecolor','none')
colorbar ; colormap(colormap_mine) ;
caxis([0.5e-4 5.5e-4]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,salt_eq_obs,[33:0.5:40],'Color','k','LineWidth',1,'ShowText','on') ;
title('Fe along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'IRON_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
















