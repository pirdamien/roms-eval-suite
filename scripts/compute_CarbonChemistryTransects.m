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
load('./colormap_calvacanti.dat');
colormap_min2=colormap_calvacanti(:,2:4) ;

par1type = 1 ; % first input parameter - Alk
par2type = 2 ; % second input parameter - 2 for DIC, 3 for pH
pHscale  = 2 ; % 1 = total pH, 2 = sea water scale
k1k2c    = 14; % Millero et al, 2010 sea water scale
kso4c    = 1 ; % bisulfate ion dissociation Dickson (1990) J. Chem. Thermodyn.
%kfconstant = 2 ;
%kbors    = 1 ; % boron:salt relationship Uppstrom 1979

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon = ncread(grid_file,'lon_rho') ;
lat = ncread(grid_file,'lat_rho') ;
h   = ncread(grid_file,'h') ;
mask= ncread(grid_file,'mask_rho') ;
pm  = ncread(grid_file,'pm') ;
pn  = ncread(grid_file,'pn') ;
lon(lon<0) = lon(lon<0)+360 ;
[NX,NY]=size(mask) ; 

lon_rshp = [140:1:270] ;
lat_rshp = [-30:1:30] ;
[lat_rshp lon_rshp] = meshgrid(lat_rshp,lon_rshp) ;
dep_rshp = [[0:1:10] [12:2:100] [105:5:300] [310:10:500]] ;
%dep_rshp = [5 25 50 100 250 500] ;
[NX,NY]=size(lat_rshp) ;

for t=1:12

    disp(num2str(t))

    zeta = ncread(file_avg,'zeta',[1 1 t],[inf inf 1]) ;
    temp = ncread(file_avg,'temp',[1 1 1 t],[inf inf inf 1]) ;
    salt = ncread(file_avg,'salt',[1 1 1 t],[inf inf inf 1]) ;
    po4 = ncread(file_bgc ,'PO4' ,[1 1 1 t],[inf inf inf 1]) ;
    sil = ncread(file_bgc ,'SiO3',[1 1 1 t],[inf inf inf 1]) ;
    dic = ncread(file_bgc ,'DIC' ,[1 1 1 t],[inf inf inf 1]) ;
    alk = ncread(file_bgc ,'Alk' ,[1 1 1 t],[inf inf inf 1]) ;
    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);

    for z=1:length(dep_rshp)

        disp(['---> ' num2str(z) '/' num2str(length(dep_rshp))])

        var=permute(vinterp(permute(temp,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        tem_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(salt,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        sal_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(po4,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        po4_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(sil,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        sil_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(dic,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        dic_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
        var=permute(vinterp(permute(alk,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        alk_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;
     
        pres_1 = abs(sw_pres(ones(NX,NY)*dep_rshp(z),lat_rshp)) ;
        rho_1 = real(sw_dens(sal_rshp(:,:,z),tem_rshp(:,:,z),pres_1)) ;
        alk_1 = real(alk_rshp(:,:,z) ./(rho_1*0.001)) ;
        dic_1 = real(dic_rshp(:,:,z) ./(rho_1*0.001)) ;
        sil_1 = real(sil_rshp(:,:,z) ./(rho_1*0.001)) ;
        po4_1 = real(po4_rshp(:,:,z) ./(rho_1*0.001)) ;
        pres_vec = reshape(pres_1,1,NX*NY);
        alk_vec  = reshape(alk_1 ,1,NX*NY);
        dic_vec  = reshape(dic_1 ,1,NX*NY);
        sil_vec  = reshape(sil_1 ,1,NX*NY);
        po4_vec  = reshape(po4_1 ,1,NX*NY);
        tem_vec  = reshape(tem_rshp(:,:,z) ,1,NX*NY);
        sal_vec  = reshape(sal_rshp(:,:,z) ,1,NX*NY);
    [DATA,HEADERS,NICEHEADERS]= ...
  CO2SYS(alk_vec,dic_vec,par1type,par2type, ...
         sal_vec,tem_vec,tem_vec,pres_vec,pres_vec, ...
         sil_vec,po4_vec,pHscale,k1k2c,kso4c);
        OmgAr_vec = DATA(:,16) ;
        pH_vec    = DATA(:,33) ;
        OmgAr_rshp(:,:,z) = real(reshape(OmgAr_vec,NX,NY)) ;
        pH_rshp(:,:,z)    = real(reshape(pH_vec   ,NX,NY)) ;

     end


dic_180160(:,:,t) = squeeze(nanmean(dic_rshp(40: 61,:,:),1)) ; dic_160120(:,:,t) = squeeze(nanmean(dic_rshp(61:101,:,:),1)) ;
alk_180160(:,:,t) = squeeze(nanmean(alk_rshp(40: 61,:,:),1)) ; alk_160120(:,:,t) = squeeze(nanmean(alk_rshp(61:101,:,:),1)) ;
omg_180160(:,:,t) = squeeze(nanmean(OmgAr_rshp(40: 61,:,:),1)) ; omg_160120(:,:,t) = squeeze(nanmean(OmgAr_rshp(61:101,:,:),1)) ;
ph_180160(:,:,t)  = squeeze(nanmean(pH_rshp(40: 61,:,:),1)) ; ph_160120(:,:,t) = squeeze(nanmean(pH_rshp(61:101,:,:),1)) ;

dic_eq(:,:,t) = squeeze(dic_rshp(:,31,:)) ;
alk_eq(:,:,t) = squeeze(alk_rshp(:,31,:)) ;
omg_eq(:,:,t) = squeeze(OmgAr_rshp(:,31,:)) ;
ph_eq(:,:,t)  = squeeze(pH_rshp(:,31,:)) ;

end

for i=1:size(lat_rshp,2)
    for j=1:length(dep_rshp)
        dep2d(i,j) = -dep_rshp(j) ;
        lat2d(i,j) = lat_rshp(1,i);
    end
end

for i=1:size(lon_rshp,1)
    for j=1:length(dep_rshp)
        depeq(i,j) = -dep_rshp(j) ;
        loneq(i,j) = lon_rshp(i,1);
    end
end

%%% Obs

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'Depth');
dic_obs = squeeze(ncread(file,'TCO2'));
dic_obs_eq = squeeze(nanmean(dic_obs(:,90:91,:),2)) ; 
dic_obs_180160 = squeeze(nanmean(dic_obs(161:180,:,:),1)) ;
dic_obs_160120 = squeeze(nanmean(dic_obs(181:220,:,:),1)) ;

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TAlk.nc' ;
alk_obs = squeeze(ncread(file,'TAlk'));
alk_obs_eq = squeeze(nanmean(alk_obs(:,90:91,:),2)) ;
alk_obs_180160 = squeeze(nanmean(alk_obs(161:180,:,:),1)) ;
alk_obs_160120 = squeeze(nanmean(alk_obs(181:220,:,:),1)) ;

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.OmegaA.nc' ;
omg_obs = squeeze(ncread(file,'OmegaA'));
omg_obs_eq = squeeze(nanmean(omg_obs(:,90:91,:),2)) ;
omg_obs_180160 = squeeze(nanmean(omg_obs(161:180,:,:),1)) ;
omg_obs_160120 = squeeze(nanmean(omg_obs(181:220,:,:),1)) ; 

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.pHtsinsitutp.nc' ;
ph_obs = squeeze(ncread(file,'pHtsinsitutp'));
ph_obs_eq = squeeze(nanmean(ph_obs(:,90:91,:),2)) ;
ph_obs_180160 = squeeze(nanmean(ph_obs(161:180,:,:),1)) ;
ph_obs_160120 = squeeze(nanmean(ph_obs(181:220,:,:),1)) ; 

for i=1:length(lat_obs)
    for j=1:length(dep_obs)
        dep2d_obs(i,j) = -dep_obs(j) ;
        lat2d_obs(i,j) = lat_obs(i);
    end
end
for i=1:length(lon_obs)
    for j=1:length(dep_obs)
        depeq_obs(i,j) = -dep_obs(j) ;
        loneq_obs(i,j) = lon_obs(i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DIC %%%%%%%%%%%%%%%%%%%

dic_180160 = squeeze(mean(dic_180160,3)) ; 
dic_160120 = squeeze(mean(dic_160120,3)) ;
dic_eq     = squeeze(mean(dic_eq    ,3)) ;

figure
%pcolor(lat2d,dep2d,dic_180160) ; shading flat ; 
contourf(lat2d,dep2d,dic_180160,[1700:5:2500],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([1900 2300]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,dic_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('DIC on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'DIC_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,dic_160120) ; shading flat ; 
contourf(lat2d,dep2d,dic_160120,[1700:5:2500],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([1900 2300]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,dic_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('DIC on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'DIC_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,dic_eq) ; shading flat ; 
contourf(loneq,depeq,dic_eq,[1700:5:2500],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([1900 2300]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,dic_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('DIC along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'DIC_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_obs,dep2d_obs,dic_obs_180160) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,dic_obs_180160,[1700:5:2500],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([1900 2300]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,dic_obs_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('DIC on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'DIC_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_obs,dep2d_obs,dic_obs_160120) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,dic_obs_160120,[1700:5:2500],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([1900 2300]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,dic_obs_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('DIC on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'DIC_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_obs,depeq_obs,dic_obs_eq) ; shading flat ; 
contourf(loneq_obs,depeq_obs,dic_obs_eq,[1700:5:2500],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([1900 2300]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,dic_obs_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('DIC along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'DIC_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Alk %%%%%%%%%%%%%%%%%%%

alk_180160 = squeeze(mean(alk_180160,3)) ;
alk_160120 = squeeze(mean(alk_160120,3)) ;
alk_eq     = squeeze(mean(alk_eq    ,3)) ;

figure
%pcolor(lat2d,dep2d,alk_180160) ; shading flat ; 
contourf(lat2d,dep2d,alk_180160,[1700:5:2600],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([2200 2500]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,alk_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('Alk on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Alk_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,alk_160120) ; shading flat ; 
contourf(lat2d,dep2d,alk_160120,[1700:5:2600],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([2200 2500]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,alk_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('Alk on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Alk_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,alk_eq) ; shading flat ; 
contourf(loneq,depeq,alk_eq,[1700:5:2600],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([2200 2500]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,alk_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('Alk along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Alk_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_obs,dep2d_obs,alk_obs_180160) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,alk_obs_180160,[1700:5:2600],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([2200 2500]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,alk_obs_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('Alk on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Alk_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_obs,dep2d_obs,alk_obs_160120) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,alk_obs_160120,[1700:5:2600],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([2200 2500]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,alk_obs_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('Alk on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Alk_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_obs,depeq_obs,alk_obs_eq) ; shading flat ; 
contourf(loneq_obs,depeq_obs,alk_obs_eq,[1700:5:2600],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([2200 2500]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,alk_obs_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('Alk along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Alk_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Omg %%%%%%%%%%%%%%%%%%%

omg_180160 = squeeze(mean(omg_180160,3)) ;
omg_160120 = squeeze(mean(omg_160120,3)) ;
omg_eq     = squeeze(mean(omg_eq    ,3)) ;

figure
%pcolor(lat2d,dep2d,omg_180160) ; shading flat ; 
contourf(lat2d,dep2d,omg_180160,[0:0.05:5],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([0 4.5]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,dic_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('\Omega_{Ar} on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Omg_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,omg_160120) ; shading flat ; 
contourf(lat2d,dep2d,omg_160120,[0:0.05:5],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([0 4.5]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,dic_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('\Omega_{Ar} on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Omg_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,omg_eq) ; shading flat ; 
contourf(loneq,depeq,omg_eq,[0:0.05:5],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([0 4.5]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,dic_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('\Omega_{Ar} along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Omg_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_obs,dep2d_obs,omg_obs_180160) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,omg_obs_180160,[0:0.05:5],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([0 4.5]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,dic_obs_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('\Omega_{Ar} on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Omg_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_obs,dep2d_obs,omg_obs_160120) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,omg_obs_160120,[0:0.05:5],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([0 4.5]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,dic_obs_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('\Omega_{Ar} on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Omg_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_obs,depeq_obs,omg_obs_eq) ; shading flat ; 
contourf(loneq_obs,depeq_obs,omg_obs_eq,[0:0.05:5],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([0 4.5]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,dic_obs_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('\Omega_{Ar} along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'Omg_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pH %%%%%%%%%%%%%%%%%%%

ph_180160 = squeeze(mean(ph_180160,3)) ;
ph_160120 = squeeze(mean(ph_160120,3)) ;
ph_eq     = squeeze(mean(ph_eq    ,3)) ;

figure
%pcolor(lat2d,dep2d,ph_180160) ; shading flat ; 
contourf(lat2d,dep2d,ph_180160,[7:0.005:9],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([7.56 8.18]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,alk_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('pH on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'pH_180160_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d,dep2d,ph_160120) ; shading flat ; 
contourf(lat2d,dep2d,ph_160120,[7:0.005:9],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([7.56 8.18]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d,dep2d,alk_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('pH on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'pH_160120_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq,depeq,ph_eq) ; shading flat ; 
contourf(loneq,depeq,ph_eq,[7:0.005:9],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([7.56 8.18]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq,depeq,alk_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('pH along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'pH_eq_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close


figure
%pcolor(lat2d_obs,dep2d_obs,ph_obs_180160) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,ph_obs_180160,[7:0.005:9],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([7.56 8.18]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,alk_obs_180160,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('pH on transect 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'pH_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(lat2d_obs,dep2d_obs,ph_obs_160120) ; shading flat ; 
contourf(lat2d_obs,dep2d_obs,ph_obs_160120,[7:0.005:9],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([7.56 8.18]) ; ylim([-500 0]) ; xlim([-30 30])
hold on ; contour(lat2d_obs,dep2d_obs,alk_obs_160120,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('pH on transect 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'pH_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close
figure
%pcolor(loneq_obs,depeq_obs,ph_obs_eq) ; shading flat ; 
contourf(loneq_obs,depeq_obs,ph_obs_eq,[7:0.005:9],'edgecolor','none')
colorbar ; colormap(colormap_min2) ;
caxis([7.56 8.18]) ; ylim([-500 0]) ; xlim([140 270])
hold on ; contour(loneq_obs,depeq_obs,alk_obs_eq,[1700:50:2500],'Color','k','LineWidth',1,'ShowText','on') ;
title('pH along the equator','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 3] ;
name_file = 'pH_eq_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg') ; print(gcf,'-painters',[rep_out name_file],'-djpeg') ;
print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
close





