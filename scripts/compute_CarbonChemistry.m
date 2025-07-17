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
    tem_surf(:,:,t) = squeeze(temp(:,:,end));
    sal_surf(:,:,t) = squeeze(salt(:,:,end));
    po4_surf(:,:,t) = squeeze(po4(:,:,end));
    sil_surf(:,:,t) = squeeze(sil(:,:,end));
    dic_surf(:,:,t) = squeeze(dic(:,:,end));
    alk_surf(:,:,t) = squeeze(alk(:,:,end));
    test=permute(vinterp(permute(temp,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     tem_200m(:,:,t) = test ;
    test=permute(vinterp(permute(salt,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     sal_200m(:,:,t) = test ;
    test=permute(vinterp(permute(po4,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     po4_200m(:,:,t) = test ;
    test=permute(vinterp(permute(sil,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     sil_200m(:,:,t) = test ;
    test=permute(vinterp(permute(dic,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     dic_200m(:,:,t) = test ;
    test=permute(vinterp(permute(alk,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     alk_200m(:,:,t) = test ;

    %%%% surf %%%%
    pres = abs(sw_pres(squeeze(z3d_v1(end,:,:))',lat)) ;
    rho = real(sw_dens(sal_surf(:,:,t),tem_surf(:,:,t),pres)) ;
    alk = real(alk_surf(:,:,t) ./(rho*0.001)) ;
    dic = real(dic_surf(:,:,t) ./(rho*0.001)) ;
    sil = real(sil_surf(:,:,t) ./(rho*0.001)) ;
    po4 = real(po4_surf(:,:,t) ./(rho*0.001)) ;
        pres_vec = reshape(pres,1,NX*NY);
        alk_vec  = reshape(alk ,1,NX*NY);
        dic_vec  = reshape(dic ,1,NX*NY);
        sil_vec  = reshape(sil ,1,NX*NY);
        po4_vec  = reshape(po4 ,1,NX*NY);
        tem_vec  = reshape(tem_surf(:,:,t) ,1,NX*NY);
        sal_vec  = reshape(sal_surf(:,:,t) ,1,NX*NY);
    [DATA,HEADERS,NICEHEADERS]= ...
  CO2SYS(alk_vec,dic_vec,par1type,par2type, ...
         sal_vec,tem_vec,tem_vec,pres_vec,pres_vec, ...
         sil_vec,po4_vec,pHscale,k1k2c,kso4c);
        OmgAr_vec = DATA(:,16) ;
        pH_vec    = DATA(:,33) ;
    OmgAr_surf(:,:,t) = real(reshape(OmgAr_vec,NX,NY)) ;
    pH_surf(:,:,t)    = real(reshape(pH_vec   ,NX,NY)) ;

    %%%% 200 %%%%
    pres = abs(sw_pres(ones(NX,NY)*200,lat)) ;
    rho = real(sw_dens(sal_200m(:,:,t),tem_200m(:,:,t),pres)) ;
    alk = real(alk_200m(:,:,t) ./(rho*0.001)) ;
    dic = real(dic_200m(:,:,t) ./(rho*0.001)) ;
    sil = real(sil_200m(:,:,t) ./(rho*0.001)) ;
    po4 = real(po4_200m(:,:,t) ./(rho*0.001)) ;
        pres_vec = reshape(pres,1,NX*NY);
        alk_vec  = reshape(alk ,1,NX*NY);
        dic_vec  = reshape(dic ,1,NX*NY);
        sil_vec  = reshape(sil ,1,NX*NY);
        po4_vec  = reshape(po4 ,1,NX*NY);
        tem_vec  = reshape(tem_200m(:,:,t) ,1,NX*NY);
        sal_vec  = reshape(sal_200m(:,:,t) ,1,NX*NY);
    [DATA,HEADERS,NICEHEADERS]= ...
  CO2SYS(alk_vec,dic_vec,par1type,par2type, ...
         sal_vec,tem_vec,tem_vec,pres_vec,pres_vec, ...
         sil_vec,po4_vec,pHscale,k1k2c,kso4c);
        OmgAr_vec = DATA(:,16) ;
        pH_vec    = DATA(:,33) ;
    OmgAr_200m(:,:,t) = real(reshape(OmgAr_vec,NX,NY)) ;
    pH_200m(:,:,t)    = real(reshape(pH_vec   ,NX,NY)) ;

end
clear zeta temp salt po4 sil dic alk

%%%%% DIC %%%%%

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TCO2.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
obs_surf = squeeze(ncread(file,'TCO2',[1 1 1 ],[inf inf 1]));
obs_200m = squeeze(ncread(file,'TCO2',[1 1 10],[inf inf 1]));
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;

figure
var = squeeze(nanmean(dic_surf,3)) ; var(mask==0) = NaN ; %var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([1900 2300])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface DIC','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanDICsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([1900 2300])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface DIC','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanDICsrf_GLODAP' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(dic_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([1900 2300])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean DIC at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanDIC200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([1900 2300])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean DIC at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanDIC200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%% Alk %%%%%

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.TAlk.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
obs_surf = squeeze(ncread(file,'TAlk',[1 1 1 ],[inf inf 1]));
obs_200m = squeeze(ncread(file,'TAlk',[1 1 10],[inf inf 1]));
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;

figure
var = squeeze(nanmean(alk_surf,3)) ; var(mask==0) = NaN ; %var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([2200 2500])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface Alk','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanALKsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([2200 2500])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface Alk','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanALKsrf_GLODAP' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(alk_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([2200 2500])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean Alk at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanALK200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[1700:5:2500],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([2200 2500])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean Alk at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanALK200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%% Omg Arg %%%%%

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.OmegaA.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
obs_surf = squeeze(ncread(file,'OmegaA',[1 1 1 ],[inf inf 1]));
obs_200m = squeeze(ncread(file,'OmegaA',[1 1 10],[inf inf 1]));
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;


figure
var = squeeze(nanmean(OmgAr_surf,3)) ; var(mask==0) = NaN ; %var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.05:5],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([0 4.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface \Omega_{Ar}','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOMGsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.05:5],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([0 4.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface \Omega_{Ar}','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOMGsrf_GLODAP' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(nanmean(OmgAr_200m,3)) ; var(mask==0) = NaN ; %var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.05:5],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([0 4.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean \Omega_{Ar} at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOMG200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.05:5],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([0 4.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean \Omega_{Ar} at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOMG200_GLODAP' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%% pH %%%%%

file = '/data/project7/pdamien/DATA/BGC/GLODAP/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.pHtsinsitutp.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
obs_surf = squeeze(ncread(file,'pHtsinsitutp',[1 1 1 ],[inf inf 1]));
obs_200m = squeeze(ncread(file,'pHtsinsitutp',[1 1 10],[inf inf 1]));
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;


figure
var = squeeze(nanmean(pH_surf,3)) ; var(mask==0) = NaN ; %var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[7:0.005:9],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([7.56 8.18])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface pH','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPHsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[7:0.005:9],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([7.56 8.18])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface pH','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPHsrf_GLODAP' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(nanmean(pH_200m,3)) ; var(mask==0) = NaN ; %var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[7:0.005:9],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([7.56 8.18])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean pH at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPH200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[7:0.005:9],'edgecolor','none') ;
colorbar ; colormap(colormap_min2) ; caxis([7.56 8.18])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean pH at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPH200_GLODAP' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close





