
clear all


grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc'
rep_out = './Fig/'
load('../colormap_IsleOfDogs.dat');
colormap_mine=colormap_IsleOfDogs(:,2:4) ;

lon = ncread(grid_file,'lon_rho') ;
lat = ncread(grid_file,'lat_rho') ;
h   = ncread(grid_file,'h') ;
mask= ncread(grid_file,'mask_rho') ;
pm  = ncread(grid_file,'pm') ;
pn  = ncread(grid_file,'pn') ;
angle = ncread(grid_file,'angle') ;
cosa = cos(angle);sina = sin(angle);
lon(lon<0) = lon(lon<0)+360 ;


% first step define index for M2 and S1

lon(900,450)
lat(900,450)
rep = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/HIStidal/' ;
files_his = dir([rep '*his*.nc']) ; 
for t=1:length(files_his)
    file = [files_his(t,1).folder '/' files_his(t,1).name] ;
    if t==1
       zeta = squeeze(ncread(file,'zeta',[900 450 1],[1 1 inf])) ;
    else
       zeta = [zeta' squeeze(ncread(file,'zeta',[900 450 1],[1 1 inf]))']' ;
    end
end
timeseries_SF1 = squeeze(zeta) - mean(squeeze(zeta)) ;

if 0

figure ; subplot(2,1,1) ; plot(timeseries_SF1) ;  grid on ;
xlabel('hour') ; ylabel('ssh') ; title('TimeSeries') 
Y1 = fft(timeseries_SF1(1:end-1));
P21 = abs(Y1/(length(timeseries_SF1)-1));
P11 = P21(1:(length(timeseries_SF1)-1)/2+1);
P11(2:end-1) = 2*P11(2:end-1);
f = (0:((length(timeseries_SF1)-1)/2))/(length(timeseries_SF1)-1);
subplot(2,1,2) ; plot(f,P11) ; hold on ;
title('Spectrum of ssh - SF')
xlabel('f (1/h)')
ylabel('|P1(f)|')
set(gca, 'YScale', 'log')
grid on ; hold on
plot([1/24 1/24],get(gca,'Ylim'),'r--')
plot([1/12 1/12],get(gca,'Ylim'),'k--')

figure
plot(P11) ; grid on ; grid minor

return
end

% Compute Pamp (probably there is a better script to build here)

istrM2 = 115 ; 
iendM2 = 130 ;
istrK1 = 55 ;
iendK1 = 70 ;

clear avg

M2amp(1:size(lat,1),1:size(lat,2))=NaN ;
K1amp(1:size(lat,1),1:size(lat,2))=NaN ;

ndomx = 10 ; ndomy = 5 ;
[Lp,Mp] = size(mask);
Tp = length(zeta);

szx = floor(Lp/ndomx);
szy = floor(Mp/ndomy);

icmin = [0:ndomx-1]*szx;
jcmin = [0:ndomy-1]*szy;
icmax = [1:ndomx]*szx;
jcmax = [1:ndomy]*szy;
icmin(1) = 1;
jcmin(1) = 1;
icmax(end) = Lp;
jcmax(end) = Mp;

for domx = 1:ndomx
      for domy = 1:ndomy
        [ domx domy]
        icb = icmin(domx);
        ice = icmax(domx);
        jcb = jcmin(domy);
        jce = jcmax(domy);
        % read zeta block
        for t=1:length(files_his)
        file = [files_his(t,1).folder '/' files_his(t,1).name] ;
        if t==1
            zeta_piece = squeeze(ncread(file,'zeta',[icb jcb 1],[ice-icb+1 jce-jcb+1 inf])) ;
        else
            zeta_piece = permute([permute(zeta_piece,[1 3 2]) ... 
                         permute(squeeze(ncread(file,'zeta',[icb jcb 1],[ice-icb+1 jce-jcb+1 inf])), ...
                         [1 3 2])],[1 3 2]) ;
        end
        end
        % fft
        zk=fft(zeta_piece,[],3);
        % reshape spectrum
        zk=abs(zk/(Tp));
        zk = zk(:,:,1:(Tp-1)/2);
        zk(:,:,2:end-1) = 2*zk(:,:,2:end-1);
        % maximum in a certain frequency segment 
        M2amp(icb:ice,jcb:jce) = max(zk(:,:,istrM2:iendM2),[],3) ;
        K1amp(icb:ice,jcb:jce) = max(zk(:,:,istrK1:iendK1),[],3) ;
      end
end

save('Zamp_M2filter_12km.mat','M2amp')
save('Zamp_K1filter_12km.mat','K1amp')

return


figure
var = M2amp ; var(mask==0) = NaN ;  
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.025:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('M2 amplitude [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'TidesM2_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = K1amp ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.025:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('K1 amplitude [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'TidesK1_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%% 

file_tide = '/data/project8/pdamien/ROMS_outputs/PACMED12KMnew/FORCINGS/pacmed12_tidestpxo9_19950101.nc' ; 
M2re_tpxo = squeeze(ncread(file_tide,'ssh_Re',[1 1 1],[inf inf 1])) ;
M2im_tpxo = squeeze(ncread(file_tide,'ssh_Im',[1 1 1],[inf inf 1])) ;
K1re_tpxo = squeeze(ncread(file_tide,'ssh_Re',[1 1 5],[inf inf 1])) ;
K1im_tpxo = squeeze(ncread(file_tide,'ssh_Im',[1 1 5],[inf inf 1])) ;

M2amp_tpxo = abs(M2re_tpxo+ 1i*M2im_tpxo);
K1amp_tpxo = abs(K1re_tpxo+ 1i*K1im_tpxo);

figure
var = M2amp_tpxo ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.025:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 1])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('M2 amplitude [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'TidesM2_tpxo' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
var = K1amp_tpxo ; var(mask==0) = NaN ;
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
%m_pcolor(lon,lat,var) ; shading flat ;
m_contourf(lon,lat,var,[0:0.025:3],'edgecolor','none') ;
colorbar ; colormap(colormap_mine) ; caxis([0 0.5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('K1 amplitude [m]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'TidesK1_tpxo' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Internal tides %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rep = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/Y*/' ;
files_dia = dir([rep 'pacmed_dia.*.nc']) ;
files_dia = files_dia(12:end) ; 

[NX,NY] = size(lon) ; 
up = zeros(NX-1,NY) ; vp = zeros(NX,NY-1) ; 
for t=1:length(files_dia)
    file = [files_dia(t,1).folder '/' files_dia(t,1).name] ;
    up = up + ncread(file,'up')/length(files_dia) ; 
    vp = vp + ncread(file,'vp')/length(files_dia) ;
end
upc = lon.*NaN ; vpc = lon.*NaN ;
upc(2:end-1,:) = 0.5.*(up(1:end-1,:)+up(2:end,:)) ; 
vpc(:,2:end-1) = 0.5.*(vp(:,1:end-1)+vp(:,2:end)) ;
UPe = upc.*cosa - vpc.*sina;
VPn = vpc.*cosa + upc.*sina;

load('../redblue_update.mat')

figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon,lat,UPe) ; shading flat ;
%m_contourf(lon,lat,UPe,[-20.1:0.2:20.1],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-5 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('depth-integrated baroclinic zonal pressure flux [W.m^{-2}]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'baroclinicUP_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

figure
hold on
m_proj('robinson','lon',[102 293],'lat',[-56 67]);
m_pcolor(lon,lat,VPn) ; shading flat ;
%m_contourf(lon,lat,VPn,[-20.1:0.2:20.1],'edgecolor','none') ;
colorbar ; colormap(C) ; caxis([-5 5])
m_gshhs_l('patch',[.8 .8 .8]);
m_grid('tickdir','in','yaxislocation','left','xaxislocation','bottom','ticklen',.01,'fontsize',12);
m_plot(lon(1,:),lat(1,:),'LineWidth',1,'Color','k') ;
m_plot(lon(end,:),lat(end,:),'LineWidth',1,'Color','k') ;
m_plot(lon(:,1),lat(:,1),'LineWidth',1,'Color','k') ;
m_plot(lon(:,end),lat(:,end),'LineWidth',1,'Color','k') ;
title('depth-integrated baroclinic meridional pressure flux [W.m^{-2}]','fontsize',16) ; set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 6 5] ;
name_file = 'baroclinicVP_roms' ;
print(gcf,'-painters',[rep_out name_file],'-djpeg')
print(gcf,'-painters',[rep_out name_file],'-djpeg')
close

















