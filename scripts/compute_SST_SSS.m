
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc' 
rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/mean_2000_2005/' 
rep_out = './Fig/'
file = [rep_in 'pacmed_avg.nc']; 

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
    zeta = ncread(file,'zeta',[1 1 t],[inf inf 1]) ;  
    temp = ncread(file,'temp',[1 1 1 t],[inf inf inf 1]) ;
    salt = ncread(file,'salt',[1 1 1 t],[inf inf inf 1]) ;
    [z3d_v1,Cw] = zlevs4(h,squeeze(mean(zeta,3)), theta_s, theta_b, hc, NZ , 'r',sc_type); 
    temp_surf(:,:,t) = squeeze(temp(:,:,end));
    salt_surf(:,:,t) = squeeze(salt(:,:,end));
    test=permute(vinterp(permute(temp,[3 2 1]),-abs(permute(z3d_v1,[1 3 2])),-abs(200)),[2 1]) ;
    test(find(test<0)) = epsilon;
    temp_200m(:,:,t) = test ;
    test=permute(vinterp(permute(salt,[3 2 1]),-abs(permute(z3d_v1,[1 3 2])),-abs(200)),[2 1]) ;
    test(find(test<0)) = epsilon;
    salt_200m(:,:,t) = test ;
end
clear temp zeta salt

file = '/data/project1/data/WOA18/temperature/woa18_decav_t00_04.nc' ;
lon_aviso = ncread(file,'lon');
lat_aviso = ncread(file,'lat');
dep_aviso = ncread(file,'depth');

files = dir('/data/project1/data/WOA18/temperature/woa18_decav_t*.nc') ; 
files = files(2:13) ; 
for t=1:12
    temp_surf_aviso(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'t_an',[1 1 1 1],[inf inf 1 1]));
    temp_200m_aviso(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'t_an',[1 1 25 1],[inf inf 1 1]));
end

files = dir('/data/project1/data/WOA18/salinity/woa18_decav_s*.nc') ;
files = files(2:13) ;
for t=1:12
    salt_surf_aviso(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'s_an',[1 1 1 1],[inf inf 1 1]));
    salt_200m_aviso(:,:,t) = squeeze(ncread([files(t).folder '/' files(t).name],'s_an',[1 1 25 1],[inf inf 1 1]));
end

[lat2_aviso lon2_aviso] = meshgrid(lat_aviso,lon_aviso) ; 
lon2_aviso = [lon2_aviso(721:end,:)' lon2_aviso(1:720,:)'+360]' ; 
lat2_aviso = [lat2_aviso(721:end,:)' lat2_aviso(1:720,:)']' ;
for t=1:12
temp_surf_aviso(:,:,t) = [squeeze(temp_surf_aviso(721:end,:,t))' squeeze(temp_surf_aviso(1:720,:,t))']' ;
temp_200m_aviso(:,:,t) = [squeeze(temp_200m_aviso(721:end,:,t))' squeeze(temp_200m_aviso(1:720,:,t))']' ;
salt_surf_aviso(:,:,t) = [squeeze(salt_surf_aviso(721:end,:,t))' squeeze(salt_surf_aviso(1:720,:,t))']' ;
salt_200m_aviso(:,:,t) = [squeeze(salt_200m_aviso(721:end,:,t))' squeeze(salt_200m_aviso(1:720,:,t))']' ;
end

%%%%

figure
var = squeeze(nanmean(temp_surf,3)) ; var(mask==0) = NaN ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 31])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean SST','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSST_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(temp_surf_aviso,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 31])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean SST','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSST_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(temp_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.25:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 6])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly SST variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSST_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(temp_surf_aviso,[3 2 1]))),[2 1]) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:0.25:15],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 6])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly SST variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSST_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%

figure
var = squeeze(nanmean(temp_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 22])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('T at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanT200m_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(temp_200m_aviso,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 22])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('T at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanT200m_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(temp_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly T variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdT200m_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(temp_200m_aviso,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:0.1:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly T variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdT200m_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get the metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = permute(squeeze(std(permute(temp_surf,[3 2 1]))),[2 1]) ; mod1(mskm==0) = NaN ;
mod2 = permute(squeeze(std(permute(temp_200m,[3 2 1]))),[2 1]) ; mod2(mskm==0) = NaN ;
obs1 = permute(squeeze(std(permute(temp_surf_aviso,[3 2 1]))),[2 1]) ;
obs2 = permute(squeeze(std(permute(temp_200m_aviso,[3 2 1]))),[2 1]) ;
obs1_interp = griddata(double(lon2_aviso),double(lat2_aviso),obs1,lon,lat) ;
obs2_interp = griddata(double(lon2_aviso),double(lat2_aviso),obs2,lon,lat) ;
mod1(isnan(obs1_interp)) = NaN ; obs1_interp(isnan(mod1)) = NaN ;
mod2(isnan(obs2_interp)) = NaN ; obs2_interp(isnan(mod2)) = NaN ;

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



%%%

figure
var = squeeze(nanmean(salt_surf,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.25:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([31.5 37.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean SSS','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSSS_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(salt_surf_aviso,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:0.25:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([31.5 37.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean SSS','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanSSS_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(salt_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.05:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly SSS variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSSS_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(salt_surf_aviso,[3 2 1]))),[2 1]) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:0.05:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly SSS variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdSSS_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%


figure
var = squeeze(nanmean(salt_200m,3)) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([33.5 36.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean S at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanS200m_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(nanmean(salt_200m_aviso,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
%m_contourf(lon2_aviso,lat2_aviso,var,[0:0.1:40],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([33.5 36.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean S at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanS200m_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(salt_200m,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.005:0.3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.15])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly S variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdS200m_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(salt_200m_aviso,[3 2 1]))),[2 1]) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2_aviso,lat2_aviso,var) ; shading flat ;
m_contourf(lon2_aviso,lat2_aviso,var,[0:0.005:0.3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.15])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly S variance at 200m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdS200m_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get the metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = permute(squeeze(std(permute(salt_surf,[3 2 1]))),[2 1]) ; mod1(mskm==0) = NaN ;
mod2 = permute(squeeze(std(permute(salt_200m,[3 2 1]))),[2 1]) ; mod2(mskm==0) = NaN ;
obs1 = permute(squeeze(std(permute(salt_surf_aviso,[3 2 1]))),[2 1]) ;
obs2 = permute(squeeze(std(permute(salt_200m_aviso,[3 2 1]))),[2 1]) ;
obs1_interp = griddata(double(lon2_aviso),double(lat2_aviso),obs1,lon,lat) ;
obs2_interp = griddata(double(lon2_aviso),double(lat2_aviso),obs2,lon,lat) ;
mod1(isnan(obs1_interp)) = NaN ; obs1_interp(isnan(mod1)) = NaN ;
mod2(isnan(obs2_interp)) = NaN ; obs2_interp(isnan(mod2)) = NaN ;

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


%%%%%%%%%

%sal_surf_aviso_interp = griddata(double(lon2d1_aviso_reshape), ...
%                                 double(lat2d1_aviso_reshape), ... 
%                                 double(sali_surf_aviso_reshape), ...
%                                 double(lon), ...
%                                 double(lat),'linear') ;








