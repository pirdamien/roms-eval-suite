
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     MLD Model                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc'
rep_out = './Fig/'
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

load('MLD_model_JAMSTEC.mat')


figure
var = squeeze(nanmean(MLD_mod.MLD(:,:,:),3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 80])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanMLD_romsJAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (MLD_mod.MLD(:,:,12)+MLD_mod.MLD(:,:,1)+MLD_mod.MLD(:,:,2))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('winter mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'winterMLD_romsJAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (MLD_mod.MLD(:,:,6)+MLD_mod.MLD(:,:,7)+MLD_mod.MLD(:,:,8))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('summer mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'summerMLD_romsJAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(MLD_mod.MLD,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:100],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly MLD variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdMLD_romsJAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('MLD_model_GOSML.mat')


figure
var = squeeze(nanmean(MLD_mod.MLD(:,:,:),3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 80])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanMLD_romsGOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (MLD_mod.MLD(:,:,12)+MLD_mod.MLD(:,:,1)+MLD_mod.MLD(:,:,2))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('winter mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'winterMLD_romsGOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (MLD_mod.MLD(:,:,6)+MLD_mod.MLD(:,:,7)+MLD_mod.MLD(:,:,8))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('summer mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'summerMLD_romsGOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(MLD_mod.MLD,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:2.5:100],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly MLD variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdMLD_romsGOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% dataset_MLD/GOSML/

file = 'dataset_MLD/GOSML/mixed_layer_properties_mean.nc' ; 

lon_obs = ncread(file,'longitude') ; 
lat_obs = ncread(file,'latitude') ; 
mld_obs = ncread(file,'depth_mean') ; 
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ; 

figure
var = squeeze(nanmean(mld_obs,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 80])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanMLD_GOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (mld_obs(:,:,12)+mld_obs(:,:,1)+mld_obs(:,:,2))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('winter mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'winterMLD_GOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (mld_obs(:,:,6)+mld_obs(:,:,7)+mld_obs(:,:,8))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('summer mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'summerMLD_GOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(mld_obs,[3 2 1]))),[2 1]) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:100],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly MLD variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdMLD_GOSML' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get the metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mod1 = squeeze(nanmean(MLD_mod.MLD(:,:,:),3)) ;
mod2 = squeeze( 1/3 * (MLD_mod.MLD(:,:,12)+MLD_mod.MLD(:,:,1)+MLD_mod.MLD(:,:,2))) ;
mod3 = squeeze( 1/3 * (MLD_mod.MLD(:,:,6)+MLD_mod.MLD(:,:,7)+MLD_mod.MLD(:,:,8))) ;
mod4 = permute(squeeze(std(permute(MLD_mod.MLD,[3 2 1]))),[2 1]) ;

obs1 = squeeze(nanmean(mld_obs,3)) ;
obs2 = squeeze( 1/3 * (mld_obs(:,:,12)+mld_obs(:,:,1)+mld_obs(:,:,2))) ;
obs3 = squeeze( 1/3 * (mld_obs(:,:,6)+mld_obs(:,:,7)+mld_obs(:,:,8))) ;
obs4 = permute(squeeze(std(permute(mld_obs,[3 2 1]))),[2 1]) ;

obs1_interp = griddata(lon_obs,lat_obs,obs1,lon,lat) ;  
obs2_interp = griddata(lon_obs,lat_obs,obs2,lon,lat) ;
obs3_interp = griddata(lon_obs,lat_obs,obs3,lon,lat) ;
obs4_interp = griddata(lon_obs,lat_obs,obs4,lon,lat) ;

mod1(isnan(obs1_interp)) = NaN ; obs1_interp(isnan(mod1)) = NaN ;
mod2(isnan(obs2_interp)) = NaN ; obs2_interp(isnan(mod2)) = NaN ;
mod3(isnan(obs3_interp)) = NaN ; obs3_interp(isnan(mod3)) = NaN ;
mod4(isnan(obs4_interp)) = NaN ; obs4_interp(isnan(mod4)) = NaN ;


vmod = mod2 ; 
vobs = obs2_interp ; 
% Compute Standard Deviations
STD_mod = std(vmod(~isnan(vmod))); STD_obs = std(vobs(~isnan(vobs)));;
% Compute Relative Standard Deviation
rel_STD = STD_mod / STD_obs;
% Compute RMSE (Root Mean Square Error)
RMSE = sqrt(mean((vmod(~isnan(vmod)) - vobs(~isnan(vobs))).^2));
% Compute Relative RMSE
rel_RMSE = RMSE / STD_obs;
% Compute Pearson Correlation Coefficient
COR = corr(vmod(~isnan(vmod)), vobs(~isnan(vobs)));
disp(['Relative STD: ' num2str(rel_STD)]);
disp(['Relative RMSE: ' num2str(rel_RMSE)]);
disp(['Pearson COR: ' num2str(COR)]);
disp(['Obs STD: ' num2str(STD_obs)]);
%disp(['Mod STD: ' num2str(STD_mod)]);
disp(['RMSE: ' num2str(RMSE)]);
% Check Taylor Diagram Equation
check_value = rel_RMSE - sqrt(rel_STD^2 + 1 - 2 * rel_STD * COR);


load('../redblue_update.mat')


figure
var = mod1 - obs1_interp ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-100.5:1:100.5],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-50 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('\Delta annual mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanMLD_DELTA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = mod2 - obs2_interp  ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-100.5:1:100.5],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-50 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('\Delta winter mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'winterMLD_DELTA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = mod3 - obs3_interp  ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-100.5:1:100.5],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-50 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('\Delta summer mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'summerMLD_DELTA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = mod4 - obs4_interp ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-100.5:1:100.5],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-30 30])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('\Delta monthly MLD variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdMLD_DELTA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close











%%%%%%%%%% JAMSTEC/

rep = 'dataset_MLD/JAMSTEC/' ;
list = dir([rep '*.nc']) ;
file = [list(1).folder '/' list(1).name]
lon_obs = ncread(file,'LONGITUDE') ;
lon_obs(lon_obs<0) = lon_obs(lon_obs<0)+360 ; 
%lon_obs = [lon_obs(181:end)' lon_obs(1:180)']' ;  
lat_obs = ncread(file,'LATITUDE') ;
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;

for t=1:12
    file = [rep 'ml_MLD_CLIM_' sprintf('%02d',t) '.nc'] ;
    var = ncread(file,'MLD') ; 
%    var = [var(181:end,:)' var(1:180,:)']' ;
    mld_obs(:,:,t) = inpaint_nans(var,2); ;
%    mld_obs(:,:,t) =  var ; 
end

figure
var = squeeze(nanmean(mld_obs,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 80])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanMLD_JAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze( 1/3 * (mld_obs(:,:,12)+mld_obs(:,:,1)+mld_obs(:,:,2))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('winter mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'winterMLD_JAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze( 1/3 * (mld_obs(:,:,6)+mld_obs(:,:,7)+mld_obs(:,:,8))) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:300],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 120])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('summer mean MLD','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'summerMLD_JAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(mld_obs,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:2.5:100],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 50])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly MLD variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdMLD_JAMSTEC' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

