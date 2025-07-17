
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc' 
rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/mean_2000_2005/' 
rep_out = './Fig/'
file_avg = [rep_in 'pacmed_avg.nc'];
file_bgc = [rep_in 'pacmed_bgc_avg.nc']; 
file_dia = [rep_in 'pacmed_bgc_dia_avg.nc'];

theta_s = 6.0; theta_b = 6.0;
hc = 250; sc_type = 'new2012';
NZ = 100 ; epsilon = 1e-3 ;

co2_file = '/data/project7/pdamien/DATA/BGC/pCO2/OCADS/MPI_SOM_FFN_2022_NCEI_OCADS.nc' ;
n2o_file = '/data/project7/pdamien/DATA/BGC/Ncycle/n2oFlux-Yang2020.nc' ; 

load('../colormap_IsleOfDogs.dat');
colormap_mine=colormap_IsleOfDogs(:,2:4) ;
load('../redblue_update.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon = ncread(grid_file,'lon_rho') ;
lat = ncread(grid_file,'lat_rho') ;
h   = ncread(grid_file,'h') ;
mask= ncread(grid_file,'mask_rho') ;
pm  = ncread(grid_file,'pm') ;
pn  = ncread(grid_file,'pn') ;
lon(lon<0) = lon(lon<0)+360 ;

for t=1:12
    disp(num2str(t)) 
    CO2flx(:,:,t) = ncread(file_dia,'FG_CO2',[1 1 t],[inf inf 1]) ; % mmol/m2/s 
    N2Oflx(:,:,t) = ncread(file_dia,'FG_N2O',[1 1 t],[inf inf 1]) ; % mmol/m2/s
end

CO2flx_obs = ncread(co2_file,'fgco2_smoothed',[1 1 241],[inf inf 12]) ; % mol/m2/yr
CO2flx_lon = ncread(co2_file,'lon') ;
CO2flx_lat = ncread(co2_file,'lat') ;
% convert to mmol/m2/s
CO2flx_obs = CO2flx_obs / (365.25*24*60*60) * 1000 ; 
[CO2flx_lat,CO2flx_lon] = meshgrid(CO2flx_lat,CO2flx_lon) ;
CO2flx_lon = [CO2flx_lon(181:end,:)' CO2flx_lon(1:180,:)'+360]' ; 
CO2flx_lat = [CO2flx_lat(181:end,:)' CO2flx_lat(1:180,:)']' ;
CO2flx_obs = permute([permute(CO2flx_obs(181:end,:,:),[2 1 3]) ...
                      permute(CO2flx_obs(  1:180,:,:),[2 1 3])],[2 1 3]) ;


N2Oflx_obs1= ncread(n2o_file,'n2oFlux_EnsMean_g-pm2-pyr') ; % g/m2/yr
N2Oflx_obs2= ncread(n2o_file,'n2oFluxSeas_g-pm2-pyr')     ; % g/m2/yr
N2Oflx_lon = ncread(n2o_file,'longitude_2d') ;
N2Oflx_lat = ncread(n2o_file,'latitude_2d') ;
% convert to mmol/m2/s
N2Oflx_obs1 = N2Oflx_obs1 * 1000 / 44.013 / (365.25*24*60*60) ; 
N2Oflx_obs2 = N2Oflx_obs2 * 1000 / 44.013 / (365.25*24*60*60) ;
N2Oflx_lon = [N2Oflx_lon(721:end,:)' N2Oflx_lon(1:720,:)'+360]' ;
N2Oflx_lat = [N2Oflx_lat(721:end,:)' N2Oflx_lat(1:720,:)']' ;
N2Oflx_obs1= permute([permute(N2Oflx_obs1(721:end,:,:),[2 1 3]) ...
                      permute(N2Oflx_obs1(  1:720,:,:),[2 1 3])],[2 1 3]) ;
N2Oflx_obs2= permute([permute(N2Oflx_obs2(721:end,:,:),[2 1 3]) ...
                      permute(N2Oflx_obs2(  1:720,:,:),[2 1 3])],[2 1 3]) ;



%%%%%%%%%%%%%%%%%%%%%%% CO2 %%%%%%%%%%%%%%%%%%%%%

figure
var = squeeze(mean(CO2flx,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-10.025e-4:0.05e-4:10.025e-4],'edgecolor','none') ;
cb1 = colorbar ; colormap(C) ;  caxis([-2e-4 2e-4]) ; 
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean CO_2 flux','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanCO2flx_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(mean(CO2flx_obs,3)) ; var(abs(var)>1e10)=NaN ; var=-var ;                                              
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(CO2flx_lon,CO2flx_lat,var) ; shading flat ;
m_contourf(CO2flx_lon,CO2flx_lat,var,[-10.025e-4:0.05e-4:10.025e-4],'edgecolor','none') ;
cb1 = colorbar ; colormap(C) ;  caxis([-2e-4 2e-4]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean CO_2 flux','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanCO2flx_Obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = permute(squeeze(std(permute(CO2flx,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.05e-4:5e-4],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1.2e-4]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly CO_2 flux variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdCO2flx_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(CO2flx_obs,[3 2 1]))),[2 1]) ; 
var(abs(var)>1e4)=NaN ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(CO2flx_lon,CO2flx_lat,var) ; shading flat ;
m_contourf(CO2flx_lon,CO2flx_lat,var,[0:0.05e-4:5e-4],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1.2e-4]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly CO_2 flux variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdCO2flx_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%% N2O %%%%%%%%%%%%%%%

figure
var = squeeze(mean(N2Oflx,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-10.01e-7:0.02e-7:10.01e-7],'edgecolor','none') ;
cb1 = colorbar ; colormap(C) ;  caxis([-1e-7 1e-7]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean N_2O flux','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanN2Oflx_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(mean(N2Oflx_obs1,3)) ; var(abs(var)>1e10)=NaN ; var=-var ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(N2Oflx_lon,N2Oflx_lat,var) ; shading flat ;
m_contourf(N2Oflx_lon,N2Oflx_lat,var,[-10.01e-7:0.02e-7:10.01e-7],'edgecolor','none') ;
cb1 = colorbar ; colormap(C) ;  caxis([-1e-7 1e-7]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean N_2O flux','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanN2Oflx_Obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close



figure
var = permute(squeeze(std(permute(N2Oflx,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.05e-7:5e-7],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1.0e-7]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly N_2O flux variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdN2Oflx_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close



figure
var = permute(squeeze(std(permute(N2Oflx_obs1,[3 2 1]))),[2 1]) ;
var(abs(var)>1e4)=NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(N2Oflx_lon,N2Oflx_lat,var) ; shading flat ;
m_contourf(N2Oflx_lon,N2Oflx_lat,var,[0:0.05e-7:5e-7],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1.0e-7]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly N_2O flux variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdN2Oflx_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get the metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = permute(squeeze(std(permute(CO2flx,[3 2 1]))),[2 1]) ; ; mod1(mskm==0) = NaN ;
mod2 = permute(squeeze(std(permute(N2Oflx,[3 2 1]))),[2 1]) ; ; mod2(mskm==0) = NaN ;
obs1 = permute(squeeze(std(permute(CO2flx_obs,[3 2 1]))),[2 1]) ; obs1(abs(obs1)>1e4)=NaN ;
obs2 = permute(squeeze(std(permute(N2Oflx_obs1,[3 2 1]))),[2 1]) ; obs2(abs(obs2)>1e4)=NaN ;
obs1_interp = griddata(double(CO2flx_lon),double(CO2flx_lat),double(obs1),lon,lat) ;
obs2_interp = griddata(double(N2Oflx_lon),double(N2Oflx_lat),double(obs2),lon,lat) ;
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






