
%% plot diagnostic P output

close all
clear all

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc'
rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/mean_2000_2005/'
rep_out = './Fig/'
file_avg = [rep_in 'pacmed_avg.nc'];
file_bgc = [rep_in 'pacmed_bgc_avg.nc'];
file_dia = [rep_in 'pacmed_bgc_dia_avg.nc'];
file_flx = [rep_in 'pacmed_flx3.nc'];

theta_s = 6.0; theta_b = 6.0;
hc = 250; sc_type = 'new2012';
NZ = 100 ; epsilon = 1e-3 ;

load('../colormap_IsleOfDogs.dat');
colormap_mine=colormap_IsleOfDogs(:,2:4) ;
load('./colormap_calvacanti.dat');
colormap_min2=colormap_calvacanti(:,2:4) ;

lon = ncread(grid_file,'lon_rho') ;
lat = ncread(grid_file,'lat_rho') ;
h = ncread(grid_file,'h') ;
mask = ncread(grid_file,'mask_rho') ;
pm = ncread(grid_file,'pm') ;
pn = ncread(grid_file,'pn') ;
lon(lon<0) = lon(lon<0) + 360 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1

%% load diagnostic output

load diagP_output.mat

% unpack
grid1 = output.grid;
AREA = repmat(grid1.Areat,[1 1 length(grid1.zt)]);
M3d = output.M3d;
iocn = find(M3d==1);
[ny,nx,nz] = size(M3d);
MSKS = output.MSKS;
R3d = output.R3d;
regnames = output.regnames;
nreg = length(regnames);
Pflux = output.Pflux;

%% make basin profiles

% trim the continental shelves/slopes
Zbot = repmat(sum(M3d.*grid1.DZT3d,3),[1 1 nz]);
mask2k = double(Zbot>2000);
mask2k(mask2k==0) = NaN;
R3d = mask2k.*R3d;

% loop over regions
regs = [1:12];
for k = 1:length(regs)
    mask = double(R3d==regs(k));
    Fprof(:,k) = squeeze(nansum(nansum(Pflux.*mask.*AREA,1),2)./nansum(nansum(mask.*AREA,1),2));
end

lon_data = [1:2:359] ; 
lat_data = [-90:2:90] ; 

for i=1:length(lat_data)
    for j=1:length(lon_data)
        lon2d_data(i,j) = lon_data(j) ; 
        lat2d_data(i,j) = lat_data(i) ; 
    end
end

% figure
% pcolor(lon2d_data,lat2d_data,squeeze(R3d(:,:,1))) ; shading flat ;
% colorbar

depth_interp =[10 20 30 40 50 60 70 80 90 100 120 140 160 180 200 ...
    220 240 260 280 300 350 400 450 500 600 700 800 900 1000 1500 2000 ...
    2500 3000 4000 5000];

region_interp = griddata(lon2d_data,lat2d_data,squeeze(R3d(:,:,1)),double(lon),double(lat)) ;
% for z=1:length(depth_interp)
%     region3d_interp(:,:,z) = region_interp ;
% end


figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_contourf(lon2d_data,lat2d_data,squeeze(R3d(:,:,1)),[0:1:20],'edgecolor',[0.6 0.6 0.6],'LineWidth',2) ;
colormap(parula) ; caxis([0 12])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('POC regions','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'POC_regions' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

for t=1:12
   
    zeta = squeeze(ncread(file_avg,'zeta',[1 1 t],[inf inf 1])) ;   
    POC  = squeeze(ncread(file_dia,'POC_FLUX_IN',[1 1 1 t],[inf inf inf 1])) ;
    PARinc = squeeze(ncread(file_flx,'PARinc',[1 1 t],[inf inf 1])) ;   
    PAR  = squeeze(ncread(file_flx,'PAR',[1 1 1 t],[inf inf inf 1])) ; 

    PARw(:,:,1:NZ-1) = 0.5.*(PAR(:,:,1:end-1)+PAR(:,:,2:end)) ; 
    PARw(:,:,NZ) = PARinc ; 
   
    % compute euphotic layer : 
    eupho_mask = PAR < 0.01*squeeze(PARw(:,:,end)) ;
    [~, eupho_index] = min(eupho_mask, [], 3); % size (lon, lat)
    % Handle cases where no depth satisfies the condition (all false in mask)
    eupho_index(~any(eupho_mask,3)) = NaN;
 
    [z3d_v1,Cw] = zlevs4(h,zeta, theta_s, theta_b, hc, NZ, 'w',sc_type);
     z3d_v1 = permute(z3d_v1(2:NZ+1,:,:),[2 3 1]) ;
     
     for z=1:length(depth_interp) 
     disp(['time = ' num2str(t) ' , z = ' num2str(z) '/' num2str(length(depth_interp))])    
     POC_interp = vinterp(permute(POC,[3 2 1]),-abs(permute(z3d_v1,[3 2 1])), ...
              -abs(depth_interp(z)) ) ; 
     POC_interp = permute(POC_interp,[2 1]) ;     
     POCmean_roms(z,t,1) = nanmean(POC_interp(region_interp==4)) ;
     POCmean_roms(z,t,2) = nanmean(POC_interp(region_interp==7)) ;
     POCmean_roms(z,t,3) = nanmean(POC_interp(region_interp==10)) ;
     POCmean_roms(z,t,4) = nanmean(POC_interp(region_interp==11)) ;
     POCmean_roms(z,t,5) = nanmean(POC_interp(region_interp==2)) ;
     POCstd_roms(z,t,1) = nanstd(POC_interp(region_interp==4)) ;
     POCstd_roms(z,t,2) = nanstd(POC_interp(region_interp==7)) ;
     POCstd_roms(z,t,3) = nanstd(POC_interp(region_interp==10)) ;
     POCstd_roms(z,t,4) = nanstd(POC_interp(region_interp==11)) ;
     POCstd_roms(z,t,5) = nanstd(POC_interp(region_interp==2)) ;
     end

     POC_interp = vinterp(permute(POC,[3 2 1]),-abs(permute(z3d_v1,[3 2 1])),-abs(100)) ;
     POC100(:,:,t) = POC_interp' ;
     POC_interp = vinterp(permute(POC,[3 2 1]),-abs(permute(z3d_v1,[3 2 1])),-abs(75)) ; 
     POC75(:,:,t)  = POC_interp' ; 
     clear POCeuphotic
     for i=1:size(lon,1)
     for j=1:size(lon,2)
         if isnan(eupho_index(i,j))
            POCeuphotic(i,j,t) = NaN ; 
         else
            POCeuphotic(i,j,t) = POC(i,j,eupho_index(i,j)) ;
         end
     end
     end

end
     
romsPOC.POCmean_roms =  POCmean_roms ;
romsPOC.POCstd_roms =  POCstd_roms ;
romsPOC.depth = depth_interp ;
romsPOC.name = {'ST.Pac';'Trop.Pac';'N.Pac';'Pac.OMZ';'SubAnt'};
romsPOC.POC100 = POC100 ; 
romsPOC.POC75 = POC75 ;
romsPOC.POCeuphotic = POCeuphotic ;

save('romsPOC.mat','romsPOC')

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('romsPOC.mat')
rep_out = './Fig/'

load diagP_output.mat

% unpack
grid1 = output.grid;
AREA = repmat(grid1.Areat,[1 1 length(grid1.zt)]);
M3d = output.M3d;
iocn = find(M3d==1);
[ny,nx,nz] = size(M3d);
MSKS = output.MSKS;
R3d = output.R3d;
regnames = output.regnames;
nreg = length(regnames);
Pflux = output.Pflux;

%% make basin profiles

% trim the continental shelves/slopes
Zbot = repmat(sum(M3d.*grid1.DZT3d,3),[1 1 nz]);
mask2k = double(Zbot>2000);
mask2k(mask2k==0) = NaN;
R3d = mask2k.*R3d;

% loop over regions
regs = [1:12];
for k = 1:length(regs)
    mask = double(R3d==regs(k));
    Fprof(:,k) = squeeze(nansum(nansum(Pflux.*mask.*AREA,1),2)./nansum(nansum(mask.*AREA,1),2));
end

%% plot profiles

% set up figure
%f1 = figure;
%set(f1,'color','w','position',[100 100 1800 600]);

I = find(grid1.zw>100);

% loop over regions

    indexes = [4 7 10 11 2] ; 

    for t=1:length(indexes) 
    i = indexes(t)
    figure 
    hold on 
    X_plot = [romsPOC.depth , fliplr(romsPOC.depth)] ;
    Y_plot = [(squeeze(nanmean(romsPOC.POCmean_roms(:,:,t),2)) + ... 
         squeeze(nanmean(romsPOC.POCstd_roms(:,:,t),2)))'*24*60*60 , ...
        fliplr(squeeze(nanmean(romsPOC.POCmean_roms(:,:,t),2))'- ...
         squeeze(nanmean(romsPOC.POCstd_roms(:,:,t),2))')*24*60*60] ;        
    fill(Y_plot, X_plot , 1,'facecolor',[0, 0.4470, 0.7410],'edgecolor','none','facealpha', 0.1);
    plot(squeeze(nanmean(romsPOC.POCmean_roms(:,:,t),2))*24*60*60, ...
        romsPOC.depth,'Color',[0, 0.4470, 0.7410],'linewidth',2);
    plot((squeeze(nanmean(romsPOC.POCmean_roms(:,:,t),2))+ ...
         squeeze(nanmean(romsPOC.POCstd_roms(:,:,t),2)))*24*60*60, ...
        romsPOC.depth,'--','Color',[0, 0.4470, 0.7410],'linewidth',1);
    plot((squeeze(nanmean(romsPOC.POCmean_roms(:,:,t),2))- ...
          squeeze(nanmean(romsPOC.POCstd_roms(:,:,t),2)))*24*60*60, ...
        romsPOC.depth,'--','Color',[0, 0.4470, 0.7410],'linewidth',1);
    plot(Fprof(I,i)/0.00943/365.24,grid1.zw(I),'Color',[0.8500, 0.3250, 0.0980],'linewidth',2);
    set(gca,'Ydir','reverse') ; ylim([0 2000])
    xlim([0 12]) ; grid on ; xlabel('POC flux [mmol m^{-2} d^{-1}]')
    title(regnames{i},'fontsize',16) ; set(gca,'fontsize',14) ; box on ; 
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 2] ;
    name_file = [rep_out 'meanPOCprofile_' num2str(i) '.svg'] ;
    plot2svg(name_file)
    close
     end


figure
var = squeeze(mean(romsPOC.POC100,3))*24*60*60 ; % mmol/m2/d  
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.5:100],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 20]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('POC export [mmolC m^{-2} d^{-1}]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPOC_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%% Clements 2023 %%%%% 

file_obs = './dataset_POC/clements2023/Euphotic_Export.nc' ; 
poc_obs = ncread(file_obs,'Flux') ; 
lon_obs = ncread(file_obs,'lon') ;
lat_obs = ncread(file_obs,'lat') ;
[lat_obs lon_obs] = meshgrid(lat_obs,lon_obs) ;


figure
var = squeeze(nanmean(poc_obs,3))/ 12.011 ; % mmol/m2/d  
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.5:100],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 20]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('POC export [mmolC m^{-2} d^{-1}]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPOC_obs_clements' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close


%%%%% Nowicki 2022 %%%%% 

file_obs = './dataset_POC/Nowicki2022/biopump_model_output.nc' ;
poc_obs = ncread(file_obs,'POCflux') ;
lon_obs = ncread(file_obs,'LON') ;
lat_obs = ncread(file_obs,'LAT') ;
dep_obs = ncread(file_obs,'DEPTH') ;

poc_obs = 0.85 * poc_obs(:,:,3) + 0.15 * poc_obs(:,:,4) ; 
lon_obs = 0.85 * lon_obs(:,:,3) + 0.15 * lon_obs(:,:,4) ;
lat_obs = 0.85 * lat_obs(:,:,3) + 0.15 * lat_obs(:,:,4) ;
dep_obs = 0.85 * dep_obs(:,:,3) + 0.15 * dep_obs(:,:,4) ;

figure
var = squeeze(nanmean(poc_obs,3)) / 365.25 ; % mmol/m2/d  
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon_obs,lat_obs,var) ; shading flat ;
m_contourf(lon_obs,lat_obs,var,[0:0.5:100],'edgecolor','none') ;
cb1 = colorbar ; colormap(colormap_mine) ; caxis([0 20]) ;
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('POC export [mmolC m^{-2} d^{-1}]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'meanPOC_obs_Nowicki' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

