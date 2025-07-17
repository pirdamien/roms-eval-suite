
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

chl_file = '/data/project7/pdamien/DATA/BGC/Aqua-MODIS/*.L3m.MC.CHL.chlor_a.4km.nc' ;

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
    disp(num2str(t)) 
    chl = ncread(file_dia,'TOT_CHL',[1 1 1 t],[inf inf inf 1]) ;
    chl_surf(:,:,t) = squeeze(chl(:,:,end));
end

rep = dir(chl_file);
rep = [rep(8:13)' rep(1:6)']';
for t=1:length(rep(:,1))
    file = [rep(t,1).folder '/' rep(t,1).name]
    var = ncread(file,'chlor_a') ;
    chl_obs(:,:,t) = inpaint_nans(var,2) ;
end
lon_obs = ncread(file,'lon') ;
lat_obs = ncread(file,'lat') ;
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;

for t=1:12
chl_obs(:,:,t) = [squeeze(chl_obs(4321:end,:,t))' squeeze(chl_obs(1:4320,:,t))']' ; 
end
lon_obs = [lon_obs(4321:end,:)' lon_obs(1:4320,:)'+360]' ;
lat_obs = [lat_obs(4321:end,:)' lat_obs(1:4320,:)']' ;



figure
var = squeeze(mean(chl_surf,3)) ; var(mask==0) = NaN ; var(var<0)=0.00001 ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon,lat,real(log10(var))) ; shading flat ;
%m_contourf(lon,lat,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; 
caxis([log10(0.02) log10(2)]) ; 
set(cb1,'XTick',[log10(0.02) log10(0.05) log10(0.1) log10(0.5) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.02' '0.05' '0.1' '0.5' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean surface Chl','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanCHL_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(mean(chl_obs,3)) ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,real(log10(var))) ; shading flat ;
%m_contourf(lon_obs,lat_obs,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ;
caxis([log10(0.02) log10(2)]) ;
set(cb1,'XTick',[log10(0.02) log10(0.05) log10(0.1) log10(0.5) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.02' '0.05' '0.1' '0.5' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean surface Chl','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanCHL_MODIS' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close



figure
var = permute(squeeze(std(permute(chl_surf,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
var(var<1e-6)=1e-6 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon,lat,real(log10(var))) ; shading flat ;
%m_contourf(lon,lat,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; 
caxis([log10(0.001) log10(2)]) ;
set(cb1,'XTick',[log10(0.001) log10(0.01) log10(0.1) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.001' '0.01' '0.1' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly Chl variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdCHL_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = permute(squeeze(std(permute(chl_obs,[3 2 1]))),[2 1]) ;
var(var<1e-6)=1e-6 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,real(log10(var))) ; shading flat ;
%m_contourf(lon_obs,lat_obs,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ;
caxis([log10(0.001) log10(2)]) ;
set(cb1,'XTick',[log10(0.001) log10(0.01) log10(0.1) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.001' '0.01' '0.1' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly Chl variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdCHL_MODIS' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close




clear chl_obs
chl_file = '/data/project7/pdamien/DATA/BGC/SeaWifs/Monthly_clim/*L3m.MC.CHL.chlor_a.9km.nc' ;

rep = dir(chl_file);
rep = [rep(5:12)' rep(1:4)']';
for t=1:length(rep(:,1))
    file = [rep(t,1).folder '/' rep(t,1).name]
    var = ncread(file,'chlor_a') ;
    chl_obs(:,:,t) = inpaint_nans(var,2) ;
%     chl_obs(:,:,t) = var ;
end
lon_obs = ncread(file,'lon') ;
lat_obs = ncread(file,'lat') ;
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;

for t=1:12
chl_obs(:,:,t) = [squeeze(chl_obs(2161:end,:,t))' squeeze(chl_obs(1:2160,:,t))']' ;
end
lon_obs = [lon_obs(2161:end,:)' lon_obs(1:2160,:)'+360]' ;
lat_obs = [lat_obs(2161:end,:)' lat_obs(1:2160,:)']' ;

figure
var = squeeze(mean(chl_obs,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,real(log10(var))) ; shading flat ;
%m_contourf(lon_obs,lat_obs,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ;
caxis([log10(0.02) log10(2)]) ;
set(cb1,'XTick',[log10(0.02) log10(0.05) log10(0.1) log10(0.5) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.02' '0.05' '0.1' '0.5' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean surface Chl','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanCHL_SeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(chl_obs,[3 2 1]))),[2 1]) ;
var(var<1e-6)=1e-6 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,real(log10(var))) ; shading flat ;
%m_contourf(lon_obs,lat_obs,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ;
caxis([log10(0.001) log10(2)]) ;
set(cb1,'XTick',[log10(0.001) log10(0.01) log10(0.1) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.001' '0.01' '0.1' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly Chl variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdCHL_SeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%%%%%%%%%%%%

clear chl_obs
chl_file = '/data/project7/pdamien/DATA/BGC/SeaWifs/Monthly/SEASTAR_SEAWIFS_GAC.' ;


for t=1:12
    rep = dir([chl_file '*' num2str(t,'%2.2d') '01*.nc']);
    for i=1:length(rep(:,1))
        file = [rep(i,1).folder '/' rep(i,1).name]
        if i==1        
        var = ncread(file,'chlor_a')/length(rep(:,1)) ;
        else 
        var = var + ncread(file,'chlor_a')/length(rep(:,1)) ;
        end
        chl_obs(:,:,t) = inpaint_nans(var,2) ;
    end    
end
lon_obs = ncread(file,'lon') ;
lat_obs = ncread(file,'lat') ;
[lat_obs,lon_obs] = meshgrid(lat_obs,lon_obs) ;

for t=1:12
chl_obs(:,:,t) = [squeeze(chl_obs(2161:end,:,t))' squeeze(chl_obs(1:2160,:,t))']' ;
end
lon_obs = [lon_obs(2161:end,:)' lon_obs(1:2160,:)'+360]' ;
lat_obs = [lat_obs(2161:end,:)' lat_obs(1:2160,:)']' ;

figure
var = squeeze(mean(chl_obs,3)) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,real(log10(var))) ; shading flat ;
%m_contourf(lon_obs,lat_obs,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ;
caxis([log10(0.02) log10(2)]) ;
set(cb1,'XTick',[log10(0.02) log10(0.05) log10(0.1) log10(0.5) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.02' '0.05' '0.1' '0.5' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual mean surface Chl','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanCHL_SeaWifs20002005' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(chl_obs,[3 2 1]))),[2 1]) ;
var(var<1e-6)=1e-6 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,real(log10(var))) ; shading flat ;
%m_contourf(lon_obs,lat_obs,real(log10(var)),[-3:0.05:1],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ;
caxis([log10(0.001) log10(2)]) ;
set(cb1,'XTick',[log10(0.001) log10(0.01) log10(0.1) log10(1) log10(2)])
set(cb1,'XTickLabel',{'0.001' '0.01' '0.1' '1' '2'})
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly Chl variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdCHL_SeaWifs20002005' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get the metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = permute(squeeze(std(permute(chl_surf,[3 2 1]))),[2 1]) ; mod1(mskm==0) = NaN ;
obs1 = permute(squeeze(std(permute(chl_obs,[3 2 1]))),[2 1])  ; 
mod1(mod1<1e-6)=1e-6 ; obs1(obs1<1e-6)=1e-6 ;
obs1_interp = griddata(double(lon_obs),double(lat_obs),obs1,lon,lat) ;
mod1(isnan(obs1_interp)) = NaN ; obs1_interp(isnan(mod1)) = NaN ;

vmod = mod1 ;
vobs = obs1_interp ;
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















