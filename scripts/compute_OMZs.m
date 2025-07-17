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
[NX,NY]=size(lon) ; 

if 0

for t=1:12

    t

    zeta = ncread(file_avg,'zeta',[1 1 t],[inf inf 1]) ;
    oxy = ncread(file_bgc , 'O2' ,[1 1 1 t],[inf inf inf 1]) ;
    [z3d_r,Cr] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);
    [z3d_w,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'w',sc_type);
    dz = diff(z3d_w); dz=permute(dz,[3 2 1]) ;
    z3d_r = permute(z3d_r,[3 2 1]) ;
    z3d_w = permute(z3d_w,[3 2 1]) ; 

    omz_mask = oxy < 5 ;
    THICKomz_005(:,:,t) = squeeze(sum(omz_mask.*dz,3)) ;     
    omz_mask = oxy < 10 ;                            
    THICKomz_010(:,:,t) = squeeze(sum(omz_mask.*dz,3)) ;
    omz_mask = oxy < 25 ;                            
    THICKomz_025(:,:,t) = squeeze(sum(omz_mask.*dz,3)) ;
    omz_mask = oxy < 50 ;
    THICKomz_050(:,:,t) = squeeze(sum(omz_mask.*dz,3)) ;
    omz_mask = oxy < 80 ;
    THICKomz_080(:,:,t) = squeeze(sum(omz_mask.*dz,3)) ;
    omz_mask = oxy < 100 ;
    THICKomz_100(:,:,t) = squeeze(sum(omz_mask.*dz,3)) ;

    test=permute(vinterp(permute(oxy,[3 2 1]),-abs(permute(z3d_r,[3 2 1])),-abs(200)),[2 1]) ;
    test(test<0) = epsilon ;     oxy_200m(:,:,t) = test ;
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ; 
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:)).*pm(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1)).*pn(:,2:end-1) ;
    Hgrad_200(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ; 
  
    test=permute(vinterp(permute(oxy,[3 2 1]),-abs(permute(z3d_r,[3 2 1])),-abs(300)),[2 1]) ;
    test(test<0) = epsilon ;     oxy_300m(:,:,t) = test ;
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:)).*pm(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1)).*pn(:,2:end-1) ;
    Hgrad_300(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    test=permute(vinterp(permute(oxy,[3 2 1]),-abs(permute(z3d_r,[3 2 1])),-abs(400)),[2 1]) ;
    test(test<0) = epsilon ;     oxy_400m(:,:,t) = test ;
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:)).*pm(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1)).*pn(:,2:end-1) ;
    Hgrad_400(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    test=permute(vinterp(permute(oxy,[3 2 1]),-abs(permute(z3d_r,[3 2 1])),-abs(500)),[2 1]) ;
    test(test<0) = epsilon ;     oxy_500m(:,:,t) = test ;
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:)).*pm(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1)).*pn(:,2:end-1) ;
    Hgrad_500(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    minOXY(:,:,t) = zeros(NX,NY)*NaN ;
    depMIN(:,:,t) = zeros(NX,NY)*NaN ;
    for i=1:NX
    for j=1:NY
        vec = squeeze(oxy(i,j,:)) ; 
        [O2_min, idx_min] = min(vec) ; 
        if ~isempty(idx_min)
           minOXY(i,j,t) = O2_min ; 
           depMIN(i,j,t) = z3d_r(i,j,idx_min) ;
        end
    end
    end 

end

    OMZs.THICKomz_005 = squeeze(mean(THICKomz_005,3)) ; 
    OMZs.THICKomz_010 = squeeze(mean(THICKomz_010,3)) ;
    OMZs.THICKomz_025 = squeeze(mean(THICKomz_025,3)) ;
    OMZs.THICKomz_050 = squeeze(mean(THICKomz_050,3)) ;
    OMZs.THICKomz_080 = squeeze(mean(THICKomz_080,3)) ;
    OMZs.THICKomz_100 = squeeze(mean(THICKomz_100,3)) ;
    OMZs.oxy_200m = squeeze(mean(oxy_200m,3)) ;
    OMZs.oxy_300m = squeeze(mean(oxy_300m,3)) ;
    OMZs.oxy_400m = squeeze(mean(oxy_400m,3)) ;
    OMZs.oxy_500m = squeeze(mean(oxy_500m,3)) ;
    OMZs.Hgrad_200 = squeeze(mean(Hgrad_200,3)) ;
    OMZs.Hgrad_300 = squeeze(mean(Hgrad_300,3)) ;
    OMZs.Hgrad_400 = squeeze(mean(Hgrad_400,3)) ;
    OMZs.Hgrad_500 = squeeze(mean(Hgrad_500,3)) ;
    OMZs.minOXY = squeeze(mean(minOXY,3)) ;
    OMZs.depMIN = squeeze(mean(depMIN,3)) ;

    save('OMZs.mat','OMZs','-v7.3')

    return
end

load('OMZs.mat')

file = '/data/project7/pdamien/DATA/BGC/WOA18/oxyg/woa18_all_o00_01.nc' ;
lon_obs = ncread(file,'lon');
lat_obs = ncread(file,'lat');
dep_obs = ncread(file,'depth');
files = dir('/data/project7/pdamien/DATA/BGC/WOA18/oxyg/woa18_all_o*.nc') ;
files = files(2:13) ;

depw = [0 0.5*(dep_obs(1:end-1)+dep_obs(2:end))' dep_obs(end)]' ; 
dz = diff(depw) ; 
for i=1:length(lon_obs)
for j=1:length(lat_obs)
    dz3d(i,j,:) = dz ;
end
end

[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;
lon_obs = [lon_obs(181:end,:)' lon_obs(1:180,:)'+360]' ;
lat_obs = [lat_obs(181:end,:)' lat_obs(1:180,:)']' ;
[NX,NY]=size(lon_obs) ;

dx = ones(NX,NY) * NaN ; 
dy = ones(NX,NY) * NaN ;
for i=2:NX-1
for j=2:NY-1
    dx(i,j) = sw_dist([0.5*(lat_obs(i-1,j)+lat_obs(i,j)) 0.5*(lat_obs(i+1,j)+lat_obs(i,j))], ...
                      [0.5*(lon_obs(i-1,j)+lon_obs(i,j)) 0.5*(lon_obs(i+1,j)+lon_obs(i,j))],'km')*1000 ; 
    dy(i,j) = sw_dist([0.5*(lat_obs(i,j-1)+lat_obs(i,j)) 0.5*(lat_obs(i,j+1)+lat_obs(i,j))], ...
                      [0.5*(lon_obs(i,j-1)+lon_obs(i,j)) 0.5*(lon_obs(i,j+1)+lon_obs(i,j))],'km')*1000 ;
end
end

for t=1:12

    oxy_obs = squeeze(ncread(file,'o_an')) ; 
    oxy_obs(:,:,1:57) = squeeze(ncread([files(t).folder '/' files(t).name],'o_an')) ; 
    oxy_obs = permute([permute(oxy_obs(181:end,:,:),[2 1 3]) permute(oxy_obs(1:180,:,:),[2 1 3])],[2 1 3]) ;

    omz_mask = oxy_obs < 5 ;
    THICKomz_005_obs(:,:,t) = squeeze(sum(omz_mask.*dz3d,3)) ;
    omz_mask = oxy_obs < 10 ;
    THICKomz_010_obs(:,:,t) = squeeze(sum(omz_mask.*dz3d,3)) ;
    omz_mask = oxy_obs < 25 ;
    THICKomz_025_obs(:,:,t) = squeeze(sum(omz_mask.*dz3d,3)) ;
    omz_mask = oxy_obs < 50 ;
    THICKomz_050_obs(:,:,t) = squeeze(sum(omz_mask.*dz3d,3)) ;
    omz_mask = oxy_obs < 80 ;
    THICKomz_080_obs(:,:,t) = squeeze(sum(omz_mask.*dz3d,3)) ;
    omz_mask = oxy_obs < 100 ;
    THICKomz_100_obs(:,:,t) = squeeze(sum(omz_mask.*dz3d,3)) ;

    test = squeeze(oxy_obs(:,:,25)) ;
    obs_200m(:,:,t) = test ; 
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:))./dx(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1))./dy(:,2:end-1) ;
    Hgrad_200_obs(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    test = squeeze(oxy_obs(:,:,29)) ;
    obs_300m(:,:,t) = test ;                    
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:))./dx(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1))./dy(:,2:end-1) ;
    Hgrad_300_obs(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    test = squeeze(oxy_obs(:,:,33)) ;
    obs_400m(:,:,t) = test ;                    
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:))./dx(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1))./dy(:,2:end-1) ;
    Hgrad_400_obs(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    test = squeeze(oxy_obs(:,:,37)) ;
    obs_500m(:,:,t) = test ;                    
    testi = 0.5*(test(1:end-1,:)+test(2:end,:)) ;
    dodx = zeros(NX,NY)*NaN ; dodx(2:end-1,:) = (testi(2:end,:)-testi(1:end-1,:))./dx(2:end-1,:) ;
    testj = 0.5*(test(:,1:end-1)+test(:,2:end)) ;
    dody = zeros(NX,NY)*NaN ; dody(:,2:end-1) = (testj(:,2:end)-testj(:,1:end-1))./dy(:,2:end-1) ;
    Hgrad_500_obs(:,:,t) = sqrt( dodx.*dodx + dody.*dody ) ;

    minOXY_obs(:,:,t) = zeros(NX,NY)*NaN ;
    depMIN_obs(:,:,t) = zeros(NX,NY)*NaN ;
    for i=1:NX
    for j=1:NY
        vec = squeeze(oxy_obs(i,j,:)) ;
        [O2_min, idx_min] = min(vec) ;
        if ~isempty(idx_min)
           minOXY_obs(i,j,t) = O2_min ;
           depMIN_obs(i,j,t) = dep_obs(idx_min) ;
        end
    end
    end

end

msk_obs = oxy_obs(:,:,1) ; 
msk_obs = msk_obs*0+1 ; 


%%%%%%%% Thickness

figure
var1 = 0.4*(OMZs.oxy_200m+OMZs.oxy_300m+OMZs.oxy_400m+OMZs.oxy_500m) ; 
var = OMZs.THICKomz_050 ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:100:5000],'edgecolor','none') ;
colorbar ; colormap(parula) ; caxis([0 1000])
m_contour(lon,lat,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('OMZ thickness','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'OMZthickness_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close



figure
var1 = 0.4*squeeze(mean(obs_200m+obs_300m+obs_400m+obs_500m,3)) ;
var = squeeze(mean(THICKomz_050_obs,3)) ; var = var.*msk_obs ; 
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:100:5000],'edgecolor','none') ;
colorbar ; colormap(parula) ; caxis([0 1000])
m_contour(lon_obs,lat_obs,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('OMZ thickness','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'OMZthickness_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%%%%%%%%% Min Oxygen

figure
var1 = 0.4*(OMZs.oxy_200m+OMZs.oxy_300m+OMZs.oxy_400m+OMZs.oxy_500m) ;
var = OMZs.minOXY ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:5:500],'edgecolor','none') ;
colorbar ; colormap(parula) ; caxis([0 150])
%m_contour(lon,lat,var,[10 50 100 150],'Color','k','Showtext','on') ;
m_contour(lon,lat,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Minimum O_2','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'MinO2_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close



figure
var1 = 0.4*squeeze(mean(obs_200m+obs_300m+obs_400m+obs_500m,3)) ;
var = squeeze(mean(minOXY_obs,3)) ; var = var.*msk_obs ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:5:500],'edgecolor','none') ;
colorbar ; colormap(parula) ; caxis([0 150])
%m_contour(lon_obs,lat_obs,var,[10 50 100 150],'Color','k','Showtext','on') ;
m_contour(lon_obs,lat_obs,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Minimum O_2','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'MinO2_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


mskm= ncread('/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc','mask_rho') ;
mod1 = OMZs.minOXY ; mod1(mskm==0) = NaN ;
obs1 = squeeze(mean(minOXY_obs,3)) ; obs1 = obs1.*msk_obs ;
obs1_interp = griddata(double(lon_obs),double(lat_obs),obs1,lon,lat) ;
mod1(isnan(obs1_interp)) = NaN ; obs1_interp(isnan(mod1)) = NaN ;




%%%%%%%%%%%% Depth Min Oxygen


figure
var1 = 0.4*(OMZs.oxy_200m+OMZs.oxy_300m+OMZs.oxy_400m+OMZs.oxy_500m) ;
var = OMZs.depMIN ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[-5000:20:0],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([-1200 0])
%m_contour(lon,lat,var,[10 50 100 150],'Color','k','Showtext','on') ;
m_contour(lon,lat,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Depth of Minimum O_2','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'MinO2Depth_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


figure
var1 = 0.4*squeeze(mean(obs_200m+obs_300m+obs_400m+obs_500m,3)) ;
var = -squeeze(mean(depMIN_obs,3)) ; var = var.*msk_obs ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[-5000:20:0],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([-1200 0])
%m_contour(lon_obs,lat_obs,var,[10 50 100 150],'Color','k','Showtext','on') ;
m_contour(lon_obs,lat_obs,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Depth of Minimum O_2','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'MinO2Depth_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%%%% H GRAD %%%%%%%%%%%%%%%%%%%%%%%%%

figure
var1 = 0.4*(OMZs.oxy_200m+OMZs.oxy_300m+OMZs.oxy_400m+OMZs.oxy_500m) ;
%var = 0.4*(OMZs.Hgrad_200+OMZs.Hgrad_300+OMZs.Hgrad_400+OMZs.Hgrad_500) ;
var = OMZs.Hgrad_300 ;
var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon,lat,var) ; shading flat ;
%m_contourf(lon,lat,var,[0:100:5000],'edgecolor','none') ;
colorbar ; colormap(parula) ; caxis([0 15e-5])
m_contour(lon,lat,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Horizontal O_2 gradient at 300m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'HgradO2_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close



figure
var1 = 0.4*squeeze(mean(obs_200m+obs_300m+obs_400m+obs_500m,3)) ;
var = squeeze(mean(Hgrad_300_obs,3)) ; var = var.*msk_obs ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
%m_contourf(lon_obs,lat_obs,var,[0:100:5000],'edgecolor','none') ;
colorbar ; colormap(parula) ; caxis([0 15e-5])
m_contour(lon_obs,lat_obs,var1,[10 50 100 150],'Color','k','Showtext','on') ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('Horizontal O_2 gradient at 300m','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'HgradO2_WOA' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%% HISTO 

var = OMZs.THICKomz_005.*mask./(pm.*pn) ;
vol005_roms = sum(var(~isnan(var)))/1e9 ; 
var = OMZs.THICKomz_010.*mask./(pm.*pn) ;
vol010_roms = sum(var(~isnan(var)))/1e9 ;
var = OMZs.THICKomz_025.*mask./(pm.*pn) ;
vol025_roms = sum(var(~isnan(var)))/1e9 ;
var = OMZs.THICKomz_050.*mask./(pm.*pn) ;
vol050_roms = sum(var(~isnan(var)))/1e9 ;
var = OMZs.THICKomz_080.*mask./(pm.*pn) ;
vol080_roms = sum(var(~isnan(var)))/1e9 ;
var = OMZs.THICKomz_100.*mask./(pm.*pn) ;
vol100_roms = sum(var(~isnan(var)))/1e9 ;

mask_obs = griddata(lon,lat,mask,double(lon_obs),double(lat_obs),'nearest') ;
mask_nan = griddata(lon,lat,mask,double(lon_obs),double(lat_obs)) ; 
mask_obs(isnan(mask_nan)) = 0 ; 

var = squeeze(mean(THICKomz_005_obs,3)).*mask_obs.*(dx.*dy) ;
vol005_obs = sum(var(~isnan(var)))/1e9 ;
var = squeeze(mean(THICKomz_010_obs,3)).*mask_obs.*(dx.*dy) ;
vol010_obs = sum(var(~isnan(var)))/1e9 ;
var = squeeze(mean(THICKomz_025_obs,3)).*mask_obs.*(dx.*dy) ;
vol025_obs = sum(var(~isnan(var)))/1e9 ;
var = squeeze(mean(THICKomz_050_obs,3)).*mask_obs.*(dx.*dy) ;
vol050_obs = sum(var(~isnan(var)))/1e9 ;
var = squeeze(mean(THICKomz_080_obs,3)).*mask_obs.*(dx.*dy) ;
vol080_obs = sum(var(~isnan(var)))/1e9 ;
var = squeeze(mean(THICKomz_100_obs,3)).*mask_obs.*(dx.*dy) ;
vol100_obs = sum(var(~isnan(var)))/1e9 ;


threshold = [5 10 25 50 80 100] ; 
vol_roms = [vol005_roms vol010_roms vol025_roms vol050_roms vol080_roms vol100_roms] ; 
vol_obs = [vol005_obs vol010_obs vol025_obs vol050_obs vol080_obs vol100_obs] ;

x_labels = {'5', '10', '25', '50', '80', '100'};
figure;
b = bar([vol_obs', vol_roms'], 'grouped'); % Ensure column vectors for correct grouping
b(1).FaceColor = [0.3022 0.4473 0.5198];  % Light blue for vol1
b(2).FaceColor = [0.9075 0.5288 0.3265];  % Light red for vol2
xticklabels(x_labels);
xticks(1:length(x_labels)); % Ensure labels align with bars
xlabel('O_2 threshold'); ylabel('km^3');
legend({'WOA', 'roms'}, 'Location', 'northwest'); % Adjust legend
set(gca, 'FontSize', 12); % Adjust font size for readability
box on; % Add a box around the plot
grid on; % Optional, adds grid lines
title('low O_2 volume','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 2 3] ;
name_file = 'OMZhisto' ;
print(gcf,'-painters',[rep_out name_file],'-depsc')
print(gcf,'-painters',[rep_out name_file],'-depsc')
close














