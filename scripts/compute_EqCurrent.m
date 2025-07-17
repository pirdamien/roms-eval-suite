
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_outputs/PACHUG/pachug_grd.nc'
rep_in = '/data/project3/pdamien/ROMS_outputs/PACHUG/'
rep_out = './Fig/'
file = [rep_in 'pachug_avg.2017-2021.nc'];

theta_s = 6.0; theta_b = 6.0;
hc = 250; sc_type = 'new2012';
NZ = 100 ; epsilon = 1e-3 ;

%load('../colormap_IsleOfDogs.dat');
%colormap_mine=colormap_IsleOfDogs(:,2:4) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon = ncread(grid_file,'lon_rho') ;
lat = ncread(grid_file,'lat_rho') ;
h   = ncread(grid_file,'h') ;
mask= ncread(grid_file,'mask_rho') ;
pm  = ncread(grid_file,'pm') ;
pn  = ncread(grid_file,'pn') ;
angle  = ncread(grid_file,'angle') ;
lon(lon<0) = lon(lon<0)+360 ;
[NX,NY]=size(h);
for z=1:NZ
angle3d(:,:,z) = angle ;
mask3d(:,:,z) = mask ;
end

lon_rshp = [140:1:270] ; 
lat_rshp = [-20:1:20] ;
[lat_rshp lon_rshp] = meshgrid(lat_rshp,lon_rshp) ;
dep_rshp = [[0:1:10] [12:2:100] [105:5:300] [310:10:500]] ; 

%%%% Compute transects %%%

figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
m_plot([140 270],[0 0],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
%179E-160W
m_plot([179 200],[20 20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
m_plot([179 200],[-20 -20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
m_plot([179 179],[-20 20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
m_plot([200 200],[-20 20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
%160W_120W
m_plot([200 240],[20 20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
m_plot([200 240],[-20 -20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
m_plot([200 200],[-20 20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
m_plot([240 240],[-20 20],'LineWidth',2,'Color',[0, 0.4470, 0.7410]) ;
title('Transect locations','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'Transects' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%for t=1:12
%    u = squeeze(ncread(file,'u',[1 1 1 t],[inf inf inf 1])) ;
%    v = squeeze(ncread(file,'v',[1 1 1 t],[inf inf inf 1])) ;
%    zeta = squeeze(ncread(file,'zeta',[1 1 t],[inf inf 1])) ;
%    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);    
%    uct(1:NX,1:NY,1:NZ)=NaN ; vct(1:NX,1:NY,1:NZ)=NaN ;
%    uct(2:end-1,:,:) = (u(2:end,:,:)+u(1:end-1,:,:))/2 ;
%    vct(:,2:end-1,:) = (v(:,2:end,:)+v(:,1:end-1,:))/2 ;
%    uE = cos(angle) .*  uct - sin(angle) .*  vct ;
%    vN = sin(angle) .*  uct + cos(angle) .*  vct ;
%    uE(mask3d==0) = NaN ;
%    vN(mask3d==0) = NaN ;
%    for z=1:length(dep_rshp)
%        var=permute(vinterp(permute(uE,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
%        uE_rshp(:,:,z,t) = griddata(lon,lat,var,lon_rshp,lat_rshp) ; 
%    end
%end

lon = lon(1400:2800,200:1200) ;
lat = lat(1400:2800,200:1200) ;
h   =   h(1400:2800,200:1200) ;
mask= mask(1400:2800,200:1200) ;
pm  = pm(1400:2800,200:1200) ;
pn  = pn(1400:2800,200:1200) ;
angle  = angle(1400:2800,200:1200) ;
angle3d = angle3d(1400:2800,200:1200,:) ;
mask3d = mask3d(1400:2800,200:1200,:) ;



    u = squeeze(mean(ncread(file,'u'),4)) ;
    v = squeeze(mean(ncread(file,'v'),4)) ;
    zeta = squeeze(mean(ncread(file,'zeta'),3)) ;
    uct(1:NX,1:NY,1:NZ)=NaN ; vct(1:NX,1:NY,1:NZ)=NaN ;
    uct(2:end-1,:,:) = (u(2:end,:,:)+u(1:end-1,:,:))/2 ;
    vct(:,2:end-1,:) = (v(:,2:end,:)+v(:,1:end-1,:))/2 ;

    zeta = zeta(1400:2800,200:1200) ;
    uct = uct(1400:2800,200:1200,:) ;
    vct = vct(1400:2800,200:1200,:) ;

    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);
    uE = cos(angle) .*  uct - sin(angle) .*  vct ;
    vN = sin(angle) .*  uct + cos(angle) .*  vct ;
    uE(mask3d==0) = NaN ;
    vN(mask3d==0) = NaN ;
    for z=1:length(dep_rshp)
	disp([num2str(z) '/' num2str(length(dep_rshp))])
        var=permute(vinterp(permute(uE,[3 2 1]),-abs(z3d_v1),-abs(dep_rshp(z))),[2 1]) ;
        uE_rshp(:,:,z) = griddata(lon,lat,var,lon_rshp,lat_rshp) ;             
    end

uE_180160 = squeeze(nanmean( uE_rshp(40: 61,:,:),1)) ;
uE_160120 = squeeze(nanmean( uE_rshp(61:101,:,:),1)) ;
for i=1:size(lat_rshp,2)
    for j=1:length(dep_rshp)
        dep2d(i,j) = -dep_rshp(j) ;
        lat2d(i,j) = lat_rshp(1,i);
    end
end

%data1 = '/data/project1/data/Tropical_Currents_Cravatte_2017/Mean_zonal_currents_179E-160W.cdf';
%uEobs_180160 = ncread(data1,'U_SADCPD');
%dep_obs = ncread(data1,'DEPTH');
%lat_obs  = ncread(data1,'LATI');
%for i=1:length(lat_obs)
%    for j=1:length(dep_obs)
%        lat2d_obs1(i,j) = lat_obs(i) ;
%        dep2d_obs1(i,j) = dep_obs(j) ;
%    end
%end

%data1 = '/data/project1/data/Tropical_Currents_Cravatte_2017/Mean_zonal_currents_160W_120W.cdf';
%uEobs_160120 = ncread(data1,'U_SADCPD');
%dep_obs = ncread(data1,'DEPTH');
%lat_obs  = ncread(data1,'LATI');
%for i=1:length(lat_obs)
%    for j=1:length(dep_obs)
%        lat2d_obs2(i,j) = lat_obs(i) ;
%        dep2d_obs2(i,j) = dep_obs(j) ;
%    end
%end

load('../redblue_update.mat')

figure
%pcolor(lat2d,dep2d,uE_180160) ; shading flat ; 
contourf(lat2d,dep2d,uE_180160,[-1.01:0.02:1.01],'edgecolor','none')
colorbar ; colormap(C) ;
caxis([-0.5 0.50]) ; ylim([-500 0]) ; xlim([-20 20])
hold on ; contour(lat2d,dep2d,uE_180160,[-99 0],'Color','k','LineWidth',1) ; 
hold on ; contour(lat2d,dep2d,uE_180160,[-0.5 0.5],'--','Color','k','LineWidth',1) ;
title('Zonal u 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'uE_180160_romsPACHUG' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-depsc')
print(gcf,'-painters',[rep_out name_file],'-depsc')
close

figure
%pcolor(lat2d_obs1,-dep2d_obs1,uEobs_180160/100) ; shading flat ; 
contourf(lat2d_obs1,-dep2d_obs1,uEobs_180160/100,[-1.01:0.02:1.01],'edgecolor','none')
colorbar ; colormap(C) ;
caxis([-0.5 0.50]) ; ylim([-500 0]) ; xlim([-20 20])
hold on ; contour(lat2d_obs1,-dep2d_obs1,uEobs_180160/100,[-99 0],'Color','k','LineWidth',1) ; 
hold on ; contour(lat2d_obs1,-dep2d_obs1,uEobs_180160/100,[-0.5 0.5],'--','Color','k','LineWidth',1) ;
title('Zonal u 180-160','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'uE_180160_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-depsc')
print(gcf,'-painters',[rep_out name_file],'-depsc')
close

figure
%pcolor(lat2d,dep2d,uE_160120) ; shading flat ; 
contourf(lat2d,dep2d,uE_160120,[-1.01:0.02:1.01],'edgecolor','none')
colorbar ; colormap(C) ;
caxis([-0.5 0.50]) ; ylim([-500 0]) ; xlim([-20 20])
hold on ; contour(lat2d,dep2d,uE_160120,[-99 0],'Color','k','LineWidth',1) ;
hold on ; contour(lat2d,dep2d,uE_160120,[-0.5 0.5],'--','Color','k','LineWidth',1) ;
title('Zonal u 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'uE_160120_romsPACHUG' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-depsc')
print(gcf,'-painters',[rep_out name_file],'-depsc')
close

figure
%pcolor(lat2d_obs2,-dep2d_obs2,uEobs_160120/100) ; shading flat ; 
contourf(lat2d_obs2,-dep2d_obs2,uEobs_160120/100,[-1.01:0.02:1.01],'edgecolor','none')
colorbar ; colormap(C) ;
caxis([-0.5 0.50]) ; ylim([-500 0]) ; xlim([-20 20])
hold on ; contour(lat2d_obs2,-dep2d_obs2,uEobs_160120/100,[-99 0],'Color','k','LineWidth',1) ;
hold on ; contour(lat2d_obs2,-dep2d_obs2,uEobs_160120/100,[-0.5 0.5],'--','Color','k','LineWidth',1) ;
title('Zonal u 160-120','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'uE_160120_obs' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-depsc')
print(gcf,'-painters',[rep_out name_file],'-depsc')
close
















