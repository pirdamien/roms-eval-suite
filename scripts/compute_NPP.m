
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

npp1_file = '/data/project7/pdamien/DATA/BGC/NPP/Standard_VGPM_SeaWifs/' ;
npp2_file = '/data/project7/pdamien/DATA/BGC/NPP/CBPM_SeaWifs/' ;

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
    prd = squeeze(ncread(file_dia,'TOT_PROD',[1 1 1 t],[inf inf inf 1])) ;
    ssh = squeeze(ncread(file_avg,'zeta',[1 1 t],[inf inf 1])) ;
    [z_w,Cw] = zlevs4(h',ssh', theta_s, theta_b, hc, NZ , 'w',sc_type);
    dz = diff(z_w); dz=permute(dz,[3 2 1]) ; 
    prd_int = sum(prd.*dz,3) * 12 * (24*60*60) ; % z-int + convert from mmolC to mgC + from s-1 to day-1 ; 
    prd_roms(:,:,t) = prd_int;
end

%%%%

for t=1:12
    rep = dir([npp1_file 'npp.*' num2str(t,'%2.2d') '.hdf']);
    for i=1:length(rep(:,1))
        file = [rep(i,1).folder '/' rep(i,1).name]
        var = hdfread(file, 'npp') ; %mgC m-2 day-1% 
        var(var<0) = NaN ; var=fliplr(var') ;
        if i==1
        varm = var/length(rep(:,1)) ;
        else
        varm = varm + var/length(rep(:,1)) ;
        end
     end
     npp_obs1(:,:,t) = inpaint_nans(double(varm),2) ;
end
lon_obs1 = [-180+0.5*360/size(varm,1) : 360/size(varm,1) : 180-0.5*360/size(varm,1)] ;
lat_obs1 = [- 90+0.5*180/size(varm,2) : 180/size(varm,2) :  90-0.5*180/size(varm,2)] ;
[lat_obs1,lon_obs1] = meshgrid(lat_obs1,lon_obs1) ;

for t=1:12
npp_obs1(:,:,t) = [squeeze(npp_obs1(1081:end,:,t))' squeeze(npp_obs1(1:1080,:,t))']' ;
end
lon_obs1 = [lon_obs1(1081:end,:)' lon_obs1(1:1080,:)'+360]' ;
lat_obs1 = [lat_obs1(1081:end,:)' lat_obs1(1:1080,:)']' ;

%%%%

for t=1:12
    rep = dir([npp2_file 'cbpm.*' num2str(t,'%2.2d') '.hdf']);
    for i=1:length(rep(:,1))
        file = [rep(i,1).folder '/' rep(i,1).name]
        var = hdfread(file, 'npp') ; %mgC m-2 day-1% 
        var(var<0) = NaN ; var=fliplr(var') ;
        if i==1
        varm = var/length(rep(:,1)) ;
        else
        varm = varm + var/length(rep(:,1)) ;
        end
     end
     npp_obs2(:,:,t) = inpaint_nans(double(varm),2) ;
end
lon_obs2 = [-180+0.5*360/size(varm,1) : 360/size(varm,1) : 180-0.5*360/size(varm,1)] ;
lat_obs2 = [- 90+0.5*180/size(varm,2) : 180/size(varm,2) :  90-0.5*180/size(varm,2)] ;
[lat_obs2,lon_obs2] = meshgrid(lat_obs2,lon_obs2) ;

for t=1:12
npp_obs2(:,:,t) = [squeeze(npp_obs2(1081:end,:,t))' squeeze(npp_obs2(1:1080,:,t))']' ;
end
lon_obs2 = [lon_obs2(1081:end,:)' lon_obs2(1:1080,:)'+360]' ;
lat_obs2 = [lat_obs2(1081:end,:)' lat_obs2(1:1080,:)']' ;

%%%%%

for t=1:12
    rep = dir([npp3_file 'cafe.*' num2str(t,'%2.2d') '.hdf']);
    for i=1:length(rep(:,1))
        file = [rep(i,1).folder '/' rep(i,1).name]
        var = hdfread(file, 'npp') ; %mgC m-2 day-1% 
        var(var<0) = NaN ; var=fliplr(var') ;
        if i==1
        varm = var/length(rep(:,1)) ;
        else
        varm = varm + var/length(rep(:,1)) ;
        end
     end
     npp_obs3(:,:,t) = inpaint_nans(double(varm),2) ;
end
lon_obs3 = [-180+0.5*360/size(varm,1) : 360/size(varm,1) : 180-0.5*360/size(varm,1)] ;
lat_obs3 = [- 90+0.5*180/size(varm,2) : 180/size(varm,2) :  90-0.5*180/size(varm,2)] ;
[lat_obs3,lon_obs3] = meshgrid(lat_obs3,lon_obs3) ;

for t=1:12
npp_obs3(:,:,t) = [squeeze(npp_obs3(1081:end,:,t))' squeeze(npp_obs3(1:1080,:,t))']' ;
end
lon_obs3 = [lon_obs3(1081:end,:)' lon_obs3(1:1080,:)'+360]' ;
lat_obs3 = [lat_obs3(1081:end,:)' lat_obs3(1:1080,:)']' ;


%%%%

figure
var = squeeze(mean(prd_roms,3)) ; var(mask==0) = NaN ; var(var<0)=0.00001 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1000]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual NPP','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNPP_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = squeeze(mean(npp_obs1,3)) ; var(var<0)=0.00001 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs1,lat_obs1,var) ; shading flat ;
m_contourf(lon_obs1,lat_obs1,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1000]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual NPP','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNPP_VGPMSeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(mean(npp_obs2,3)) ; var(var<0)=0.00001 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs2,lat_obs2,var) ; shading flat ;
m_contourf(lon_obs2,lat_obs2,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1000]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual NPP','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNPP_CBPMSeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = squeeze(mean(npp_obs3,3)) ; var(var<0)=0.00001 ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs2,lat_obs2,var) ; shading flat ;
m_contourf(lon_obs3,lat_obs3,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 1000]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('annual NPP','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanNPP_CAFESeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var = permute(squeeze(std(permute(prd_roms,[3 2 1]))),[2 1]) ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 400]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly NPP variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNPP_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(npp_obs1,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs1,lat_obs1,var) ; shading flat ;
m_contourf(lon_obs1,lat_obs1,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 400]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly NPP variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNPP_VGPMSeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(npp_obs2,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs2,lat_obs2,var) ; shading flat ;
m_contourf(lon_obs2,lat_obs2,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 400]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly NPP variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNPP_CBPMeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = permute(squeeze(std(permute(npp_obs3,[3 2 1]))),[2 1]) ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs2,lat_obs2,var) ; shading flat ;
m_contourf(lon_obs3,lat_obs3,var,[0:20:4000],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 400]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('monthly NPP variance','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'stdNPP_CAFESeaWifs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%%%%%%%%%%%%%%%%%%%%%

mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = permute(squeeze(std(permute(prd_roms,[3 2 1]))),[2 1]) ; mod1(mskm==0) = NaN ;
mod2 = permute(squeeze(std(permute(prd_roms,[3 2 1]))),[2 1]) ; mod2(mskm==0) = NaN ;
mod3 = permute(squeeze(std(permute(prd_roms,[3 2 1]))),[2 1]) ; mod3(mskm==0) = NaN ;
obs1 = permute(squeeze(std(permute(npp_obs1,[3 2 1]))),[2 1]) ;
obs2 = permute(squeeze(std(permute(npp_obs2,[3 2 1]))),[2 1]) ;
obs3 = permute(squeeze(std(permute(npp_obs3,[3 2 1]))),[2 1]) ;
obs1_interp = griddata(double(lon_obs1),double(lat_obs1),obs1,lon,lat) ;
obs2_interp = griddata(double(lon_obs2),double(lat_obs2),obs2,lon,lat) ;
obs3_interp = griddata(double(lon_obs3),double(lat_obs3),obs3,lon,lat) ;
mod1(isnan(obs1_interp)) = NaN ; obs1_interp(isnan(mod1)) = NaN ;
mod2(isnan(obs2_interp)) = NaN ; obs2_interp(isnan(mod2)) = NaN ;
mod3(isnan(obs3_interp)) = NaN ; obs3_interp(isnan(mod3)) = NaN ;

vmod = mod3 ;
vobs = obs3_interp ;
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













