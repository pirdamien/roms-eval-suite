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


for t=1:12
    zeta = ncread(file_avg,'zeta',[1 1 t],[inf inf 1]) ;
    no3 = ncread(file_bgc ,'NO3' ,[1 1 1 t],[inf inf inf 1]) ;
    po4 = ncread(file_bgc ,'PO4' ,[1 1 1 t],[inf inf inf 1]) ;
    oxy = ncread(file_bgc , 'O2' ,[1 1 1 t],[inf inf inf 1]) ;
    sil = ncread(file_bgc ,'SiO3',[1 1 1 t],[inf inf inf 1]) ;
    n2o = ncread(file_bgc ,'N2O' ,[1 1 1 t],[inf inf inf 1]) ;
    irn = ncread(file_bgc , 'Fe' ,[1 1 1 t],[inf inf inf 1]) ;
    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);
    no3_surf(:,:,t) = squeeze(no3(:,:,end));
    po4_surf(:,:,t) = squeeze(po4(:,:,end));
    oxy_surf(:,:,t) = squeeze(oxy(:,:,end));
    sil_surf(:,:,t) = squeeze(sil(:,:,end));
    n2o_surf(:,:,t) = squeeze(n2o(:,:,end));
    irn_surf(:,:,t) = squeeze(irn(:,:,end));
    test=permute(vinterp(permute(no3,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     no3_200m(:,:,t) = test ;
    test=permute(vinterp(permute(po4,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     po4_200m(:,:,t) = test ;
    test=permute(vinterp(permute(oxy,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     oxy_200m(:,:,t) = test ;
    test=permute(vinterp(permute(sil,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     sil_200m(:,:,t) = test ;
    test=permute(vinterp(permute(n2o,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     n2o_200m(:,:,t) = test ;
    test=permute(vinterp(permute(irn,[3 2 1]),-abs(z3d_v1),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     irn_200m(:,:,t) = test ;
end
clear zeta no3 po4 oxy sil n2o irn

%%%%% NO3 %%%%%

file = '/data/project7/pdamien/DATA/BGC/WOA18/nitr/woa18_all_n00_01.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'depth');
files = dir('/data/project7/pdamien/DATA/BGC/WOA18/nitr/woa18_all_n*.nc') ;
files = files(2:13) ;
for t=1:12
    obs_surf(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'n_an',[1 1 1 1],[inf inf 1 1]));
    obs_200m(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'n_an',[1 1 25 1],[inf inf 1 1]));
end
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(181:end,:)' lon_obs(1:180,:)'+360]' ;
lat_obs = [lat_obs(181:end,:)' lat_obs(1:180,:)']' ;
for t=1:12
obs_surf(:,:,t) = [squeeze(obs_surf(181:end,:,t))' squeeze(obs_surf(1:180,:,t))']' ;
obs_200m(:,:,t) = [squeeze(obs_200m(181:end,:,t))' squeeze(obs_200m(1:180,:,t))']' ;
end

figure
var = squeeze(nanmean(no3_surf,3)) ; var(mask==0) = NaN ; var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.2:30],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 10])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface nitrate','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNITsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.2:30],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 10])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface nitrate','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNITsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(no3_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 40])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean nitrate at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNIT200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 40])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean nitrate at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNIT200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(no3_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.2:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly nitrate variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNITsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_surf,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.2:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly nitrate variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNITsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(no3_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.2:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly nitrate variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNIT200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_200m,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.2:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly nitrate variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNIT200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

clear obs_surf obs_200m

%%%%% PO4 %%%%%

file = '/data/project7/pdamien/DATA/BGC/WOA18/phos/woa18_all_p00_01.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'depth');
files = dir('/data/project7/pdamien/DATA/BGC/WOA18/phos/woa18_all_p*.nc') ;
files = files(2:13) ;
for t=1:12
    obs_surf(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'p_an',[1 1 1 1],[inf inf 1 1]));
    obs_200m(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'p_an',[1 1 25 1],[inf inf 1 1]));
end
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(181:end,:)' lon_obs(1:180,:)'+360]' ;
lat_obs = [lat_obs(181:end,:)' lat_obs(1:180,:)']' ;
for t=1:12
obs_surf(:,:,t) = [squeeze(obs_surf(181:end,:,t))' squeeze(obs_surf(1:180,:,t))']' ;
obs_200m(:,:,t) = [squeeze(obs_200m(181:end,:,t))' squeeze(obs_200m(1:180,:,t))']' ;
end

figure
var = squeeze(nanmean(po4_surf,3)) ; var(mask==0) = NaN ; var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.05:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface phosphate','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPO4srf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.05:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface phosphate','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPO4srf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(po4_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1:5],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean phosphate at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPO4200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.1:5],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean phosphate at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPO4200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(po4_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.01:0.6],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly phosphate variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdPO4srf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_surf,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.01:0.6],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly phosphate variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdPO4srf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(po4_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.01:0.6],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly phosphate variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdPO4200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_200m,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.01:0.6],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly phosphate variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdPO4200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

clear obs_surf obs_200m

%%%%% OXY %%%%%

file = '/data/project7/pdamien/DATA/BGC/WOA18/oxyg/woa18_all_o00_01.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'depth');
files = dir('/data/project7/pdamien/DATA/BGC/WOA18/oxyg/woa18_all_o*.nc') ;
files = files(2:13) ;
for t=1:12
    obs_surf(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'o_an',[1 1 1 1],[inf inf 1 1]));
    obs_200m(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'o_an',[1 1 25 1],[inf inf 1 1]));
end
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(181:end,:)' lon_obs(1:180,:)'+360]' ;
lat_obs = [lat_obs(181:end,:)' lat_obs(1:180,:)']' ;
for t=1:12
obs_surf(:,:,t) = [squeeze(obs_surf(181:end,:,t))' squeeze(obs_surf(1:180,:,t))']' ;
obs_200m(:,:,t) = [squeeze(obs_200m(181:end,:,t))' squeeze(obs_200m(1:180,:,t))']' ;
end

figure
var = squeeze(nanmean(oxy_surf,3)) ; var(mask==0) = NaN ; var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:5:400],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([180 360])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface oxygen','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOXYsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:5:400],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([180 360])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface oxygen','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOXYsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(oxy_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:10:350],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 280])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean oxygen at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOXY200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:10:350],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 280])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean oxygen at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanOXY200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(oxy_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2:50],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 35])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly oxygen variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdOXYsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_surf,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2:50],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 35])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly oxygen variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdOXYsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(oxy_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 20])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly oxygen variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdOXY200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_200m,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 20])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly oxygen variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdOXY200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

clear obs_surf obs_200m


%%%%% SIL %%%%%

file = '/data/project7/pdamien/DATA/BGC/WOA18/silc/woa18_all_i00_01.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'depth');
files = dir('/data/project7/pdamien/DATA/BGC/WOA18/silc/woa18_all_i*.nc') ;
files = files(2:13) ;
for t=1:12
    obs_surf(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'i_an',[1 1 1 1],[inf inf 1 1]));
    obs_200m(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'i_an',[1 1 25 1],[inf inf 1 1]));
end
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(181:end,:)' lon_obs(1:180,:)'+360]' ;
lat_obs = [lat_obs(181:end,:)' lat_obs(1:180,:)']' ;
for t=1:12
obs_surf(:,:,t) = [squeeze(obs_surf(181:end,:,t))' squeeze(obs_surf(1:180,:,t))']' ;
obs_200m(:,:,t) = [squeeze(obs_200m(181:end,:,t))' squeeze(obs_200m(1:180,:,t))']' ;
end

figure
var = squeeze(nanmean(sil_surf,3)) ; var(mask==0) = NaN ; var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 8])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface silicate','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSILsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.1:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 8])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface silicate','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSILsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(sil_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:1:80],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 60])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean silicate at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSIL200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:1:80],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 60])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean silicate at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSIL200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(sil_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);%
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1:10],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly silicate variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSILsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_surf,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.1:10],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly silicate variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSILsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = permute(squeeze(std(permute(sil_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 8])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly silicate variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSIL200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_200m,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.1:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 8])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly silicate variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSIL200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

clear obs_surf obs_200m


%%%%% N2O %%%%%

file = '/data/project7/pdamien/DATA/BGC/Ncycle/n2ofromnn.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'DEPTH');
obs_surf = squeeze(ncread(file,'N2O',[1 1 1 1],[inf inf 1 inf]));
obs_200m = squeeze(ncread(file,'N2O',[1 1 25 1],[inf inf 1 inf]));
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(181:end,:)' lon_obs(1:180,:)'+360]' ;
lat_obs = [lat_obs(181:end,:)' lat_obs(1:180,:)']' ;
for t=1:12
obs_surf(:,:,t) = [squeeze(obs_surf(181:end,:,t))' squeeze(obs_surf(1:180,:,t))']' ;
obs_200m(:,:,t) = [squeeze(obs_200m(181:end,:,t))' squeeze(obs_200m(1:180,:,t))']' ;
end

figure
var = squeeze(nanmean(n2o_surf,3)) ; var(mask==0) = NaN ; var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.0005:0.03],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0.005 0.018])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface n2o','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanN2Osrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.0005:0.03],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0.005 0.018])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface n2o','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanN2Osrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(n2o_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.002:0.1],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0.0 0.06])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean n2o at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanN2O200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.002:0.1],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.06])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean n2o at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanN2O200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(n2o_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);%
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1e-3:10e-3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 3e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly n2o variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdN2Osrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_surf,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.1e-3:10e-3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 3e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly n2o variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdN2Osrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(n2o_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.2e-3:20e-3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 6e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly n2o variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdN2O200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_200m,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.2e-3:20e-3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 6e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly n2o variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdN2O200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

clear obs_surf obs_200m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IRON %%%%%%%%%%%%%%%%%%%%%%%%

file = '/data/project7/pdamien/DATA/BGC/MappedIron/Monthly_dFe.nc' ;
lon_obs = ncread(file,'Longitude');
lat_obs = ncread(file,'Latitude');
dep_obs = ncread(file,'Depth');
obs_surf = squeeze(ncread(file,'dFe_RF',[1 1 1 1],[inf inf 1 12]));
obs_200m = 0.25*squeeze(ncread(file,'dFe_RF',[1 1 16 1],[inf inf 1 12])) + ...
           0.75*squeeze(ncread(file,'dFe_RF',[1 1 17 1],[inf inf 1 12])) ;
 
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(179:end,:)' lon_obs(1:178,:)'+360]' ;
lat_obs = [lat_obs(179:end,:)' lat_obs(1:178,:)']' ;
for t=1:12
obs_surf(:,:,t) = inpaint_nans([squeeze(obs_surf(179:end,:,t))' squeeze(obs_surf(1:178,:,t))']',2) ;
obs_200m(:,:,t) = inpaint_nans([squeeze(obs_200m(179:end,:,t))' squeeze(obs_200m(1:178,:,t))']',2) ;
end
obs_surf = obs_surf*0.001;
obs_200m = obs_200m*0.001;


figure
var = squeeze(nanmean(irn_surf,3)) ; var(mask==0) = NaN ; var(var<0)=epsilon ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface Fe','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanIRNsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(obs_surf,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean surface Fe','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanIRNsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(irn_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean Fe at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanIRN200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(nanmean(obs_200m,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean Fe at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanIRN200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(irn_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);%
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.5e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly Fe variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdIRNsrf_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_surf,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.5e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly Fe variance at surface','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdIRNsrf_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = permute(squeeze(std(permute(irn_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.1e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly Fe variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdIRN200_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(obs_200m,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[[0:0.005e-3:0.1e-3] [0.2e-3:0.1e-3:10e-3]],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.1e-3])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Monthly Fe variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdIRN200_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

clear obs_surf obs_200m






