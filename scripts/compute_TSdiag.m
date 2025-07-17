
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

init_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/pacmed_iniY1999M01D01.nc' ;

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


% SWP : South West Pacific
zone1 = h*0 ; zone1(800:808,200:208)=1 ;
for z=1:NZ
mask1(:,:,z) =  mask(800:808,200:208) ;
end
% NWP : North West Pacific
zone2 = h*0 ; zone2(560:568,540:548)=1 ;
for z=1:NZ
mask2(:,:,z) =  mask(560:568,540:548) ;
end
% CCS : California Current System
zone3 = h*0 ; zone3(1180:1188,720:728)=1 ;
for z=1:NZ
mask3(:,:,z) =  mask(1180:1188,720:728) ;
end
% SPG : South Pacific Gyre
zone4 = h*0 ; zone4(1350:1358,210:218)=1 ;
for z=1:NZ
mask4(:,:,z) =  mask(1350:1358,210:218) ;
end
% SG : Subpolar  Gyre
zone5 = h*0 ; zone5(810:818,760:768)=1 ;
for z=1:NZ
mask5(:,:,z) =  mask(810:818,760:768) ;
end
% PC : Peru Current
zone6 = h*0 ; zone6(1630:1638,380:388)=1 ;
for z=1:NZ
mask6(:,:,z) =  mask(1630:1638,380:388) ;
end




temp_vecbound = [0:0.2:30] ; 
temp_veccentr = 0.5*(temp_vecbound(1:end-1)+temp_vecbound(2:end)) ;  

temp1_vecbound = [0:1:30] ;
temp1_veccentr = 0.5*(temp1_vecbound(1:end-1)+temp1_vecbound(2:end)) ;

%%%% Init file

temp_init_mod = ncread(init_file,'temp') ;
salt_init_mod = ncread(init_file,'salt') ;
temp_init_mod = squeeze(temp_init_mod) ;
salt_init_mod = squeeze(salt_init_mod) ;

temp1_init_mod = temp_init_mod(800:808,200:208,:) ;
salt1_init_mod = salt_init_mod(800:808,200:208,:) ;
temp1_init_mod(mask1==0)=NaN;
salt1_init_mod(mask1==0)=NaN;
vart= temp1_init_mod(:) ; vars= salt1_init_mod(:) ;
for z=1:length(temp_veccentr) 
    [A,b] = find(vart>=temp_vecbound(z) & vart<temp_vecbound(z+1)) ;
    salt1_init_mod_p05(z) = prctile(vars(A),05) ;
    salt1_init_mod_p95(z) = prctile(vars(A),95) ; 
    salt1_init_mod_mn(z) = mean(vars(A)) ;
end
temp1_init_mean = squeeze(nanmean(nanmean(temp1_init_mod,1),2)); 
salt1_init_mean = squeeze(nanmean(nanmean(salt1_init_mod,1),2));

temp2_init_mod = temp_init_mod(560:568,540:548,:) ;
salt2_init_mod = salt_init_mod(560:568,540:548,:) ;
temp2_init_mod(mask2==0)=NaN;
salt2_init_mod(mask2==0)=NaN;
vart= temp2_init_mod(:) ; vars= salt2_init_mod(:) ;
for z=1:length(temp_veccentr)
    [A,b] = find(vart>=temp_vecbound(z) & vart<temp_vecbound(z+1)) ;
    salt2_init_mod_p05(z) = prctile(vars(A),05) ;
    salt2_init_mod_p95(z) = prctile(vars(A),95) ;   
    salt2_init_mod_mn(z) = mean(vars(A)) ;
end
temp2_init_mean = squeeze(nanmean(nanmean(temp2_init_mod,1),2));
salt2_init_mean = squeeze(nanmean(nanmean(salt2_init_mod,1),2));

temp3_init_mod = temp_init_mod(1180:1188,720:728,:) ;
salt3_init_mod = salt_init_mod(1180:1188,720:728,:) ;
temp3_init_mod(mask3==0)=NaN;
salt3_init_mod(mask3==0)=NaN;
vart= temp3_init_mod(:) ; vars= salt3_init_mod(:) ;
for z=1:length(temp_veccentr)
    [A,b] = find(vart>=temp_vecbound(z) & vart<temp_vecbound(z+1)) ;
    salt3_init_mod_p05(z) = prctile(vars(A),05) ;
    salt3_init_mod_p95(z) = prctile(vars(A),95) ;   
    salt3_init_mod_mn(z) = mean(vars(A)) ;
end
temp3_init_mean = squeeze(nanmean(nanmean(temp3_init_mod,1),2));
salt3_init_mean = squeeze(nanmean(nanmean(salt3_init_mod,1),2));

temp4_init_mod = temp_init_mod(1350:1358,210:218,:) ;
salt4_init_mod = salt_init_mod(1350:1358,210:218,:) ;
temp4_init_mod(mask4==0)=NaN;
salt4_init_mod(mask4==0)=NaN;
vart= temp4_init_mod(:) ; vars= salt4_init_mod(:) ;
for z=1:length(temp_veccentr)
    [A,b] = find(vart>=temp_vecbound(z) & vart<temp_vecbound(z+1)) ;
    salt4_init_mod_p05(z) = prctile(vars(A),05) ;
    salt4_init_mod_p95(z) = prctile(vars(A),95) ;
    salt4_init_mod_mn(z) = mean(vars(A)) ;
end
temp4_init_mean = squeeze(nanmean(nanmean(temp4_init_mod,1),2));
salt4_init_mean = squeeze(nanmean(nanmean(salt4_init_mod,1),2));

temp5_init_mod = temp_init_mod(810:818,760:768,:) ;
salt5_init_mod = salt_init_mod(810:818,760:768,:) ;
temp5_init_mod(mask5==0)=NaN;
salt5_init_mod(mask5==0)=NaN;
vart= temp5_init_mod(:) ; vars= salt5_init_mod(:) ;
for z=1:length(temp_veccentr)
    [A,b] = find(vart>=temp_vecbound(z) & vart<temp_vecbound(z+1)) ;
    salt5_init_mod_p05(z) = prctile(vars(A),05) ;
    salt5_init_mod_p95(z) = prctile(vars(A),95) ;
    salt5_init_mod_mn(z) = mean(vars(A)) ;
end
temp5_init_mean = squeeze(nanmean(nanmean(temp5_init_mod,1),2));
salt5_init_mean = squeeze(nanmean(nanmean(salt5_init_mod,1),2));

temp6_init_mod = temp_init_mod(1630:1638,380:388,:) ;
salt6_init_mod = salt_init_mod(1630:1638,380:388,:) ;
temp6_init_mod(mask6==0)=NaN;
salt6_init_mod(mask6==0)=NaN;
vart= temp6_init_mod(:) ; vars= salt6_init_mod(:) ;
for z=1:length(temp_veccentr)
    [A,b] = find(vart>=temp_vecbound(z) & vart<temp_vecbound(z+1)) ;
    salt6_init_mod_p05(z) = prctile(vars(A),05) ;
    salt6_init_mod_p95(z) = prctile(vars(A),95) ;
    salt6_init_mod_mn(z) = mean(vars(A)) ;
end
temp6_init_mean = squeeze(nanmean(nanmean(temp6_init_mod,1),2));
salt6_init_mean = squeeze(nanmean(nanmean(salt6_init_mod,1),2));

%%%%% Mod profile

temp = ncread(file,'temp') ; 
salt = ncread(file,'salt') ;

    vart = squeeze(temp(800:808,200:208,:,:));
    vars = squeeze(salt(800:808,200:208,:,:));
    for t=1:size(vart,4)
         test = squeeze(vart(:,:,:,t)) ; test(mask1==0)=NaN;
         temp1_mod(:,:,:,t) = test ;
         test = squeeze(vars(:,:,:,t)) ; test(mask1==0)=NaN;
         salt1_mod(:,:,:,t) = test ;
    end

    vart = squeeze(temp(560:568,540:548,:,:));
    vars = squeeze(salt(560:568,540:548,:,:));
    for t=1:size(vart,4)
         test = squeeze(vart(:,:,:,t)) ; test(mask2==0)=NaN;
         temp2_mod(:,:,:,t) = test ;
         test = squeeze(vars(:,:,:,t)) ; test(mask2==0)=NaN;
         salt2_mod(:,:,:,t) = test ;
    end

    vart = squeeze(temp(1180:1188,720:728,:,:));
    vars = squeeze(salt(1180:1188,720:728,:,:));
    for t=1:size(vart,4)
         test = squeeze(vart(:,:,:,t)) ; test(mask3==0)=NaN;
         temp3_mod(:,:,:,t) = test ;
         test = squeeze(vars(:,:,:,t)) ; test(mask3==0)=NaN;
         salt3_mod(:,:,:,t) = test ;
    end

    vart = squeeze(temp(1350:1358,210:218,:,:));
    vars = squeeze(salt(1350:1358,210:218,:,:));
    for t=1:size(vart,4)
         test = squeeze(vart(:,:,:,t)) ; test(mask4==0)=NaN;
         temp4_mod(:,:,:,t) = test ;
         test = squeeze(vars(:,:,:,t)) ; test(mask4==0)=NaN;
         salt4_mod(:,:,:,t) = test ;
    end

    vart = squeeze(temp(810:818,760:768,:,:));
    vars = squeeze(salt(810:818,760:768,:,:));
    for t=1:size(vart,4)
         test = squeeze(vart(:,:,:,t)) ; test(mask5==0)=NaN;
         temp5_mod(:,:,:,t) = test ;
         test = squeeze(vars(:,:,:,t)) ; test(mask5==0)=NaN;
         salt5_mod(:,:,:,t) = test ;
    end

    vart = squeeze(temp(1630:1638,380:388,:,:));
    vars = squeeze(salt(1630:1638,380:388,:,:));
    for t=1:size(vart,4)
         test = squeeze(vart(:,:,:,t)) ; test(mask6==0)=NaN;
         temp6_mod(:,:,:,t) = test ;
         test = squeeze(vars(:,:,:,t)) ; test(mask6==0)=NaN;
         salt6_mod(:,:,:,t) = test ;
    end

%%%%% obs %%%%%%

file = '/data/project1/data/WOA18/temperature/woa18_decav_t00_04.nc' ;
lon_aviso = ncread(file,'lon') ;
lat_aviso = ncread(file,'lat') ;
tem_aviso = ncread(file,'t_an') ;
file = '/data/project1/data/WOA18/salinity/woa18_decav_s00_04.nc' ;
sal_aviso = ncread(file,'s_an') ;

tem_aviso(tem_aviso>1e4) = NaN;
sal_aviso(sal_aviso>1e4) = NaN;

for i=1:length(lat_aviso)
lon2d1_aviso(:,i) = lon_aviso ;
end
for i=1:length(lon_aviso)
lat2d1_aviso(i,:) = lat_aviso ;
end

tem_aviso_reshape = tem_aviso(721:end,:,:) ;
tem_aviso_reshape(721:1440,:,:) = tem_aviso(1:720,:,:) ;
lat2d_aviso_reshape = lat2d1_aviso(721:end,:) ;
lat2d_aviso_reshape(721:1440,:) = lat2d1_aviso(1:720,:) ;
lon2d_aviso_reshape = lon2d1_aviso(721:end,:) ;
lon2d_aviso_reshape(721:1440,:) = lon2d1_aviso(1:720,:)+360 ;
sal_aviso_reshape = sal_aviso(721:end,:,:) ;
sal_aviso_reshape(721:1440,:,:) = sal_aviso(1:720,:,:) ;

for z=1:size(sal_aviso_reshape,3)
    z
    sal1_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(sal_aviso_reshape(480:1200,:,z))), ...
                                double(lon(800:808,200:208)),double(lat(800:808,200:208)),'linear') ;
    tem1_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(tem_aviso_reshape(480:1200,:,z))), ...
                                double(lon(800:808,200:208)),double(lat(800:808,200:208)),'linear') ;
    sal2_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(sal_aviso_reshape(480:1200,:,z))), ...
                                double(lon(560:568,540:548)),double(lat(560:568,540:548)),'linear') ;
    tem2_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(tem_aviso_reshape(480:1200,:,z))), ...
                                double(lon(560:568,540:548)),double(lat(560:568,540:548)),'linear') ;
    sal3_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(sal_aviso_reshape(480:1200,:,z))), ...
                                double(lon(1180:1188,720:728)),double(lat(1180:1188,720:728)),'linear') ;
    tem3_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(tem_aviso_reshape(480:1200,:,z))), ...
                                double(lon(1180:1188,720:728)),double(lat(1180:1188,720:728)),'linear') ;
    sal4_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(sal_aviso_reshape(480:1200,:,z))), ...
                                double(lon(1350:1358,210:218)),double(lat(1350:1358,210:218)),'linear') ;
    tem4_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(tem_aviso_reshape(480:1200,:,z))), ...
                                double(lon(1350:1358,210:218)),double(lat(1350:1358,210:218)),'linear') ;
    sal5_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(sal_aviso_reshape(480:1200,:,z))), ...
                                double(lon(810:818,760:768)),double(lat(810:818,760:768)),'linear') ;
    tem5_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(tem_aviso_reshape(480:1200,:,z))), ...
                                double(lon(810:818,760:768)),double(lat(810:818,760:768)),'linear') ;
    sal6_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(sal_aviso_reshape(480:1200,:,z))), ...
                                double(lon(1630:1638,380:388)),double(lat(1630:1638,380:388)),'linear') ;
    tem6_obs(:,:,z) =  griddata(double(lon2d_aviso_reshape(480:1200,:)),double(lat2d_aviso_reshape(480:1200,:)), ... 
                                double(squeeze(tem_aviso_reshape(480:1200,:,z))), ...
                                double(lon(1630:1638,380:388)),double(lat(1630:1638,380:388)),'linear') ;
end

vart= tem1_obs(:) ; vars= sal1_obs(:) ;
for z=1:length(temp1_veccentr)
    [A,b] = find(vart>=temp1_vecbound(z) & vart<temp1_vecbound(z+1)) ;
    salt1_obs_p05(z) = prctile(vars(A),05) ;
    salt1_obs_p95(z) = prctile(vars(A),95) ;
    salt1_obs_mn(z) = mean(vars(A)) ;
end
temp1_obs_mean = squeeze(nanmean(nanmean(tem1_obs,1),2));
salt1_obs_mean = squeeze(nanmean(nanmean(sal1_obs,1),2));


vart= tem2_obs(:) ; vars= sal2_obs(:) ;
for z=1:length(temp1_veccentr)
    [A,b] = find(vart>=temp1_vecbound(z) & vart<temp1_vecbound(z+1)) ;
    salt2_obs_p05(z) = prctile(vars(A),05) ;
    salt2_obs_p95(z) = prctile(vars(A),95) ;
    salt2_obs_mn(z) = mean(vars(A)) ;
end
temp2_obs_mean = squeeze(nanmean(nanmean(tem2_obs,1),2));
salt2_obs_mean = squeeze(nanmean(nanmean(sal2_obs,1),2));


vart= tem3_obs(:) ; vars= sal3_obs(:) ;
for z=1:length(temp1_veccentr)
    [A,b] = find(vart>=temp1_vecbound(z) & vart<temp1_vecbound(z+1)) ;
    salt3_obs_p05(z) = prctile(vars(A),05) ;
    salt3_obs_p95(z) = prctile(vars(A),95) ;
    salt3_obs_mn(z) = mean(vars(A)) ;
end
temp3_obs_mean = squeeze(nanmean(nanmean(tem3_obs,1),2));
salt3_obs_mean = squeeze(nanmean(nanmean(sal3_obs,1),2));

vart= tem4_obs(:) ; vars= sal4_obs(:) ;
for z=1:length(temp1_veccentr)
    [A,b] = find(vart>=temp1_vecbound(z) & vart<temp1_vecbound(z+1)) ;
    salt4_obs_p05(z) = prctile(vars(A),05) ;
    salt4_obs_p95(z) = prctile(vars(A),95) ;
    salt4_obs_mn(z) = mean(vars(A)) ;
end
temp4_obs_mean = squeeze(nanmean(nanmean(tem4_obs,1),2));
salt4_obs_mean = squeeze(nanmean(nanmean(sal4_obs,1),2));

vart= tem5_obs(:) ; vars= sal5_obs(:) ;
for z=1:length(temp1_veccentr)
    [A,b] = find(vart>=temp1_vecbound(z) & vart<temp1_vecbound(z+1)) ;
    salt5_obs_p05(z) = prctile(vars(A),05) ;
    salt5_obs_p95(z) = prctile(vars(A),95) ;
    salt5_obs_mn(z) = mean(vars(A)) ;
end
temp5_obs_mean = squeeze(nanmean(nanmean(tem5_obs,1),2));
salt5_obs_mean = squeeze(nanmean(nanmean(sal5_obs,1),2));

vart= tem6_obs(:) ; vars= sal6_obs(:) ;
for z=1:length(temp1_veccentr)
    [A,b] = find(vart>=temp1_vecbound(z) & vart<temp1_vecbound(z+1)) ;
    salt6_obs_p05(z) = prctile(vars(A),05) ;
    salt6_obs_p95(z) = prctile(vars(A),95) ;
    salt6_obs_mn(z) = mean(vars(A)) ;
end
temp6_obs_mean = squeeze(nanmean(nanmean(tem6_obs,1),2));
salt6_obs_mean = squeeze(nanmean(nanmean(sal6_obs,1),2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure ; hold on ; 
%X_plot = [temp1_veccentr(~isnan(salt1_obs_p05)), fliplr(temp1_veccentr(~isnan(salt1_obs_p95)))] ;
%Y_plot = [salt1_obs_p05(~isnan(salt1_obs_p05)) , ...
%   fliplr(salt1_obs_p95(~isnan(salt1_obs_p95)))] ;
%fill(Y_plot, X_plot , 1,'facecolor',[0.8500, 0.3250, 0.0980],'edgecolor','none','facealpha', 0.1);
%X_plot = [temp_veccentr(~isnan(salt1_init_mod_p05)), fliplr(temp_veccentr(~isnan(salt1_init_mod_p95)))] ;
%Y_plot = [salt1_init_mod_p05(~isnan(salt1_init_mod_p05)) , ...
%   fliplr(salt1_init_mod_p95(~isnan(salt1_init_mod_p95)))] ;
%fill(Y_plot, X_plot , 1,'facecolor',[0, 0.4470, 0.7410],'edgecolor','none','facealpha', 0.1);
p1 = scatter(salt1_mod(:),temp1_mod(:),10,'k','.') ; 
p2 = plot(salt1_obs_mean,temp1_obs_mean,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
p3 = plot(salt1_init_mean,temp1_init_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
%p2 = plot(salt1_obs_mn,temp1_veccentr,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
%plot(salt1_obs_p05,temp1_veccentr,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
%plot(salt1_obs_p95,temp1_veccentr,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
%p3 = plot(salt1_init_mod_mn,temp_veccentr,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
%plot(salt1_init_mod_p05,temp_veccentr,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1)
%plot(salt1_init_mod_p95,temp_veccentr,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1) 
grid on ; legend([p2 p3 p1],'WOA','Glorys','ROMS','Location','SouthEast') ; box on ;
title('South West Pacific Water Masses') ; xlabel('Salinity') ; ylabel('Temperature') ;
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = [rep_out 'SWP_WaterMasses'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
print(gcf,'-painters',name_file,'-djpeg') ; print(gcf,'-painters',name_file,'-djpeg') ;
close


figure ; hold on ;
%X_plot = [temp1_veccentr(~isnan(salt2_obs_p05)), fliplr(temp1_veccentr(~isnan(salt2_obs_p95)))] ;
%Y_plot = [salt2_obs_p05(~isnan(salt2_obs_p05)) , ...
%   fliplr(salt2_obs_p95(~isnan(salt2_obs_p95)))] ;
%fill(Y_plot, X_plot , 1,'facecolor',[0.8500, 0.3250, 0.0980],'edgecolor','none','facealpha', 0.1);
%X_plot = [temp_veccentr(~isnan(salt2_init_mod_p05)), fliplr(temp_veccentr(~isnan(salt2_init_mod_p95)))] ;
%Y_plot = [salt2_init_mod_p05(~isnan(salt2_init_mod_p05)) , ...
%   fliplr(salt2_init_mod_p95(~isnan(salt2_init_mod_p95)))] ;
%fill(Y_plot, X_plot , 1,'facecolor',[0, 0.4470, 0.7410],'edgecolor','none','facealpha', 0.1);
p1 = scatter(salt2_mod(:),temp2_mod(:),10,'k','.') ;
p2 = plot(salt2_obs_mean,temp2_obs_mean,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
p3 = plot(salt2_init_mean,temp2_init_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
%p2 = plot(salt2_obs_mn,temp1_veccentr,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
%plot(salt2_obs_p05,temp1_veccentr,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
%plot(salt2_obs_p95,temp1_veccentr,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
%p3 = plot(salt2_init_mod_mn,temp_veccentr,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
%plot(salt2_init_mod_p05,temp_veccentr,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1)
%plot(salt2_init_mod_p95,temp_veccentr,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1)
grid on ; legend([p2 p3 p1],'WOA','Glorys','ROMS','Location','SouthEast') ; box on ;
title('North West Pacific Water Masses') ; xlabel('Salinity') ; ylabel('Temperature') ;
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = [rep_out 'NWP_WaterMasses'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
print(gcf,'-painters',name_file,'-djpeg') ; print(gcf,'-painters',name_file,'-djpeg') ;
close



figure ; hold on ;
%X_plot = [temp1_veccentr(~isnan(salt3_obs_p05)), fliplr(temp1_veccentr(~isnan(salt3_obs_p95)))] ;
%Y_plot = [salt3_obs_p05(~isnan(salt3_obs_p05)) , ...
%   fliplr(salt3_obs_p95(~isnan(salt3_obs_p95)))] ;
%fill(Y_plot, X_plot , 1,'facecolor',[0.8500, 0.3250, 0.0980],'edgecolor','none','facealpha', 0.1);
%X_plot = [temp_veccentr(~isnan(salt3_init_mod_p05)), fliplr(temp_veccentr(~isnan(salt3_init_mod_p95)))] ;
%Y_plot = [salt3_init_mod_p05(~isnan(salt3_init_mod_p05)) , ...
%   fliplr(salt3_init_mod_p95(~isnan(salt3_init_mod_p95)))] ;
%fill(Y_plot, X_plot , 1,'facecolor',[0, 0.4470, 0.7410],'edgecolor','none','facealpha', 0.1);
p1 = scatter(salt3_mod(:),temp3_mod(:),10,'k','.') ;
p2 = plot(salt3_obs_mean,temp3_obs_mean,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
p3 = plot(salt3_init_mean,temp3_init_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
%p2 = plot(salt3_obs_mn,temp1_veccentr,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
%plot(salt3_obs_p05,temp1_veccentr,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
%plot(salt3_obs_p95,temp1_veccentr,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1)
%p3 = plot(salt3_init_mod_mn,temp_veccentr,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
%plot(salt3_init_mod_p05,temp_veccentr,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1)
%plot(salt3_init_mod_p95,temp_veccentr,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1)
grid on ; legend([p2 p3 p1],'WOA','Glorys','ROMS','Location','SouthWest') ; box on ;
title('California Current System Water Masses') ; xlabel('Salinity') ; ylabel('Temperature') ;
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = [rep_out 'CCS_WaterMasses'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
print(gcf,'-painters',name_file,'-djpeg') ; print(gcf,'-painters',name_file,'-djpeg') ;
close

figure ; hold on ;
p1 = scatter(salt4_mod(:),temp4_mod(:),10,'k','.') ;
p2 = plot(salt4_obs_mean,temp4_obs_mean,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
p3 = plot(salt4_init_mean,temp4_init_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
grid on ; legend([p2 p3 p1],'WOA','Glorys','ROMS','Location','SouthEast') ; box on ;
title('South Pacific Gyre Water Masses') ; xlabel('Salinity') ; ylabel('Temperature') ;
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = [rep_out 'SPG_WaterMasses'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
print(gcf,'-painters',name_file,'-djpeg') ; print(gcf,'-painters',name_file,'-djpeg') ;
close

figure ; hold on ;
p1 = scatter(salt5_mod(:),temp5_mod(:),10,'k','.') ;
p2 = plot(salt5_obs_mean,temp5_obs_mean,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
p3 = plot(salt5_init_mean,temp5_init_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
grid on ; legend([p2 p3 p1],'WOA','Glorys','ROMS','Location','SouthWest') ; box on ;
title('Subpolar Gyre Water Masses') ; xlabel('Salinity') ; ylabel('Temperature') ;
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = [rep_out 'SG_WaterMasses'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
print(gcf,'-painters',name_file,'-djpeg') ; print(gcf,'-painters',name_file,'-djpeg') ;
close

figure ; hold on ;
p1 = scatter(salt6_mod(:),temp6_mod(:),10,'k','.') ;
p2 = plot(salt6_obs_mean,temp6_obs_mean,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2) ;
p3 = plot(salt6_init_mean,temp6_init_mean,'Color',[0, 0.4470, 0.7410],'LineWidth',2) ;
grid on ; legend([p2 p3 p1],'WOA','Glorys','ROMS','Location','SouthEast') ; box on ;
title('Peru Current Water Masses') ; xlabel('Salinity') ; ylabel('Temperature') ;
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = [rep_out 'PC_WaterMasses'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
print(gcf,'-painters',name_file,'-djpeg') ; print(gcf,'-painters',name_file,'-djpeg') ;
close


