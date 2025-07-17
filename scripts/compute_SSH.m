
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc'
rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/Y*/'
rep_out = './Fig/'
list = [rep_in 'pacmed_avg.*.nc'];

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

list_file = dir(list) ;
list_file = list_file(12:83) ; 
for t=1:length(list_file)
    file = [list_file(t).folder '/' list_file(t).name] 
    zeta_month(:,:,t) = ncread(file,'zeta') ;
end

load('./monthly_AVISO_2000_2005.mat')
for i=1:size(monthly_AVISO.sla_monthly,1)
    for j=1:size(monthly_AVISO.sla_monthly,2)
        lon2d_aviso(i,j) = monthly_AVISO.lon_aviso(i) ;
        lat2d_aviso(i,j) = monthly_AVISO.lat_aviso(j) ; 
    end
end

figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,squeeze(nanmean(zeta_month,3))) ; shading flat ;
m_contourf(lon,lat,squeeze(nanmean(zeta_month,3)),[-1:0.05:2],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([-0.4 1.1])
hold on 
m_contour(lon,lat,squeeze(nanmean(zeta_month,3)),[0 0.5],'--','Color','k') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean SSH [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = [rep_out 'meanSSH_roms'] ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon2d_aviso,lat2d_aviso,squeeze(nanmean(monthly_AVISO.adt_monthly,3))-0.35) ; shading flat ;
m_contourf(lon2d_aviso,lat2d_aviso, ...
           squeeze(nanmean(monthly_AVISO.adt_monthly,3))-0.35,[-1:0.05:2],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([-0.4 1.1])
hold on 
m_contour(lon2d_aviso,lat2d_aviso, ...
          squeeze(nanmean(monthly_AVISO.adt_monthly,3))-0.35,[0 0.5],'--','Color','k') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('mean SSH [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = [rep_out 'meanSSH_aviso'] ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,permute(squeeze(std(permute(zeta_month,[3 2 1]))),[2 1])) ; shading flat ;
m_contourf(lon,lat,permute(squeeze(std(permute(zeta_month,[3 2 1]))),[2 1]),[0:0.02:0.5],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.2])
hold on
m_contour(lon,lat,permute(squeeze(std(permute(zeta_month,[3 2 1]))),[2 1]),[0 0.1],'--','Color','k') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly SSH variance [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = [rep_out 'stdSSH_roms'] ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,permute(squeeze(std(permute(zeta_month,[3 2 1]))),[2 1])) ; shading flat ;
m_contourf(lon2d_aviso,lat2d_aviso, ...
       permute(squeeze(std(permute(monthly_AVISO.adt_monthly,[3 2 1]))),[2 1]),[0:0.02:0.5],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.2])
hold on
m_contour(lon2d_aviso,lat2d_aviso, ...
       permute(squeeze(std(permute(monthly_AVISO.adt_monthly,[3 2 1]))),[2 1]),[0 0.1],'--','Color','k') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly SSH variance [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = [rep_out 'stdSSH_aviso'] ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get the metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = squeeze(nanmean(zeta_month,3)) ; mod1(mskm==0) = NaN ;
mod2 = permute(squeeze(std(permute(zeta_month,[3 2 1]))),[2 1]) ; mod2(mskm==0) = NaN ;
obs1 = squeeze(nanmean(monthly_AVISO.adt_monthly,3))-0.35 ;
obs2 = permute(squeeze(std(permute(monthly_AVISO.adt_monthly,[3 2 1]))),[2 1]) ;
obs1_interp = griddata(double(lon2d_aviso),double(lat2d_aviso),obs1,lon,lat) ;
obs2_interp = griddata(double(lon2d_aviso),double(lat2d_aviso),obs2,lon,lat) ;
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


load('../redblue_update.mat')

figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,squeeze(nanmean(zeta_month,3))) ; shading flat ;
m_contourf(lon,lat,mod1-obs1_interp+0.05,[-2.02:0.01:2.02],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-0.5 0.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('\Delta mean SSH [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = ['meanSSH_DELTA'] ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,permute(squeeze(std(permute(zeta_month,[3 2 1]))),[2 1])) ; shading flat ;
m_contourf(lon,lat,mod2-obs2_interp,[[-1.01:0.005:1.01]],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-0.1 0.1])
hold on
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('\Delta monthly SSH variance [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = ['stdSSH_DELTA'] ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close




