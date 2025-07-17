
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rep25 = '/data/project8/pdamien/ROMS_outputs/PACMED25KMnew' ; 
grid25 = [rep25 '/FORCINGS/pacmed25_grd.nc'] ;
file25_avg = dir([rep25 '/RUN/LOOP*/Y*/pacmed_avg.*.nc' ]) ;
file25_bgc = dir([rep25 '/RUN/LOOP*/Y*/pacmed_bgc_avg.*.nc' ]) ; 
file25_avg = [file25_avg(2:52)' file25_avg(79:109)']' ; 
file25_bgc = [file25_bgc(1:51)' file25_bgc(78:108)']' ; 

rep12 = '/data/project8/pdamien/ROMS_outputs/PACMED12KMnew' ;
grid12 = [rep12 '/FORCINGS/pacmed12_grd.nc'] ;
file12_avg = dir([rep12 '/RUNS/Y*/pacmed_avg.*.nc' ]) ;
file12_bgc = dir([rep12 '/RUNS/Y*/pacmed_bgc_avg.*.nc' ]) ;

rep_out = './Fig/'

theta_s = 6.0; theta_b = 6.0;
hc = 250; sc_type = 'new2012';
NZ = 100 ; epsilon = 1e-3 ;

load('../colormap_IsleOfDogs.dat');
colormap_mine=colormap_IsleOfDogs(:,2:4) ;

if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i25_1 = 407 ; j25_1 = 103 ; % SWP : South West Pacific 
i25_2 = 286 ; j25_2 = 273 ; % NWP : North West Pacific
i25_3 = 599 ; j25_3 = 363 ; % CCS : California Current System
i25_4 = 685 ; j25_4 = 108 ; % SPG : South Pacific Gyre
i25_5 = 412 ; j25_5 = 383 ; % SG : Subpolar  Gyre
i25_6 = 827 ; j25_6 = 193 ; % PC : Peru Current

i12_1 = 804 ; j12_1 = 204 ; % SWP : South West Pacific 
i12_2 = 564 ; j12_2 = 544 ; % NWP : North West Pacific
i12_3 = 1184; j12_3 = 724 ; % CCS : California Current System
i12_4 = 1354; j12_4 = 214 ; % SPG : South Pacific Gyre
i12_5 = 814 ; j12_5 = 764 ; % SG : Subpolar  Gyre
i12_6 = 1634; j12_6 = 384 ; % PC : Peru Current

%%
%i = i12_6 ; j = j12_6 ;
%target_lon = lon12(i,j) ;  
%target_lat = lat12(i,j) ;
%% Assuming lon and lat are 2D matrices of size (Nx, Ny)
%dist = (lon25 - target_lon).^2 + (lat25 - target_lat).^2;  % Compute squared distance
%[~, idx] = min(dist(:));  % Find the index of the minimum distance
%% Convert linear index to row and column indices
%[i25, j25] = ind2sub(size(lon25), idx);
%% Output the indices
%disp(['Closest index: (', num2str(i25), ', ', num2str(j25), ')']);
%disp([num2str(lon25(i25,j25)) ',' num2str(lat25(i25,j25))])
%disp([num2str(lon12(i,j)) ',' num2str(lat12(i,j))])
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lon25 = ncread(grid25,'lon_rho') ; lon25(lon25<0) = lon25(lon25<0)+360 ; 
lat25 = ncread(grid25,'lat_rho') ;
h25   = ncread(grid25,'h') ;

lon12 = ncread(grid12,'lon_rho') ; lon12(lon12<0) = lon12(lon12<0)+360 ;
lat12 = ncread(grid12,'lat_rho') ;
h12   = ncread(grid12,'h') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SpinUpTS(1).h25 = h25(i25_1,j25_1) ; 
SpinUpTS(2).h25 = h25(i25_2,j25_2) ;
SpinUpTS(3).h25 = h25(i25_3,j25_3) ;
SpinUpTS(4).h25 = h25(i25_4,j25_4) ;
SpinUpTS(5).h25 = h25(i25_5,j25_5) ;
SpinUpTS(6).h25 = h25(i25_6,j25_6) ;

SpinUpTS(1).h12 = h12(i12_1,j12_1) ;
SpinUpTS(2).h12 = h12(i12_2,j12_2) ;
SpinUpTS(3).h12 = h12(i12_3,j12_3) ;
SpinUpTS(4).h12 = h12(i12_4,j12_4) ;
SpinUpTS(5).h12 = h12(i12_5,j12_5) ;
SpinUpTS(6).h12 = h12(i12_6,j12_6) ;

ind = 0 ; 
for t=1:length(file25_avg)

    ind = ind + 1 ;

    file = [file25_avg(t).folder '/' file25_avg(t).name] 

    SpinUpTS(1).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ; 
    SpinUpTS(2).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(3).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(4).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(5).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(6).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;

    var = ncread(file,'zeta') ; 
    SpinUpTS(1).zeta(ind) = var(i25_1,j25_1) ; SpinUpTS(2).zeta(ind) = var(i25_2,j25_2) ;
    SpinUpTS(3).zeta(ind) = var(i25_3,j25_3) ; SpinUpTS(4).zeta(ind) = var(i25_4,j25_4) ;
    SpinUpTS(5).zeta(ind) = var(i25_5,j25_5) ; SpinUpTS(6).zeta(ind) = var(i25_6,j25_6) ;
    var = ncread(file,'temp') ;
    SpinUpTS(1).temp(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).temp(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).temp(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).temp(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).temp(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).temp(ind,:) = var(i25_6,j25_6,:) ;
    var = ncread(file,'salt') ;
    SpinUpTS(1).salt(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).salt(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).salt(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).salt(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).salt(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).salt(ind,:) = var(i25_6,j25_6,:) ;

    file = [file25_bgc(t).folder '/' file25_bgc(t).name]

    var = ncread(file,'NO3') ;
    SpinUpTS(1).nitr(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).nitr(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).nitr(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).nitr(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).nitr(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).nitr(ind,:) = var(i25_6,j25_6,:) ;
    var = ncread(file,'Fe') ;
    SpinUpTS(1).iron(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).iron(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).iron(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).iron(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).iron(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).iron(ind,:) = var(i25_6,j25_6,:) ;
    var = ncread(file,'O2') ;
    SpinUpTS(1).oxyg(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).oxyg(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).oxyg(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).oxyg(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).oxyg(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).oxyg(ind,:) = var(i25_6,j25_6,:) ;
    var = ncread(file,'DIC') ;
    SpinUpTS(1).dic(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).dic(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).dic(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).dic(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).dic(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).dic(ind,:) = var(i25_6,j25_6,:) ;
    var = ncread(file,'Alk') ;
    SpinUpTS(1).alk(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).alk(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).alk(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).alk(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).alk(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).alk(ind,:) = var(i25_6,j25_6,:) ;
    var = ncread(file,'SiO3') ;
    SpinUpTS(1).silc(ind,:) = var(i25_1,j25_1,:) ; SpinUpTS(2).silc(ind,:) = var(i25_2,j25_2,:) ;
    SpinUpTS(3).silc(ind,:) = var(i25_3,j25_3,:) ; SpinUpTS(4).silc(ind,:) = var(i25_4,j25_4,:) ;
    SpinUpTS(5).silc(ind,:) = var(i25_5,j25_5,:) ; SpinUpTS(6).silc(ind,:) = var(i25_6,j25_6,:) ;

end

for t=1:length(file12_avg)

    ind = ind + 1 ;

    file = [file12_avg(t).folder '/' file12_avg(t).name]

    SpinUpTS(1).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(2).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(3).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(4).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(5).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;
    SpinUpTS(6).time(ind) = ncread(file,'ocean_time')/(24*60*60) + datenum(1995,1,1) ;

    var = ncread(file,'zeta') ;
    SpinUpTS(1).zeta(ind) = var(i12_1,j12_1) ; SpinUpTS(2).zeta(ind) = var(i12_2,j12_2) ;
    SpinUpTS(3).zeta(ind) = var(i12_3,j12_3) ; SpinUpTS(4).zeta(ind) = var(i12_4,j12_4) ;
    SpinUpTS(5).zeta(ind) = var(i12_5,j12_5) ; SpinUpTS(6).zeta(ind) = var(i12_6,j12_6) ;
    var = ncread(file,'temp') ;
    SpinUpTS(1).temp(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).temp(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).temp(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).temp(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).temp(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).temp(ind,:) = var(i12_6,j12_6,:) ;
    var = ncread(file,'salt') ;
    SpinUpTS(1).salt(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).salt(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).salt(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).salt(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).salt(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).salt(ind,:) = var(i12_6,j12_6,:) ;

    file = [file12_bgc(t).folder '/' file12_bgc(t).name]

    var = ncread(file,'NO3') ;
    SpinUpTS(1).nitr(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).nitr(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).nitr(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).nitr(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).nitr(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).nitr(ind,:) = var(i12_6,j12_6,:) ;
    var = ncread(file,'Fe') ;
    SpinUpTS(1).iron(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).iron(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).iron(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).iron(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).iron(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).iron(ind,:) = var(i12_6,j12_6,:) ;
    var = ncread(file,'O2') ;
    SpinUpTS(1).oxyg(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).oxyg(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).oxyg(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).oxyg(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).oxyg(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).oxyg(ind,:) = var(i12_6,j12_6,:) ;
    var = ncread(file,'DIC') ;
    SpinUpTS(1).dic(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).dic(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).dic(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).dic(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).dic(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).dic(ind,:) = var(i12_6,j12_6,:) ;
    var = ncread(file,'Alk') ;
    SpinUpTS(1).alk(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).alk(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).alk(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).alk(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).alk(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).alk(ind,:) = var(i12_6,j12_6,:) ;
    var = ncread(file,'SiO3') ;
    SpinUpTS(1).silc(ind,:) = var(i12_1,j12_1,:) ; SpinUpTS(2).silc(ind,:) = var(i12_2,j12_2,:) ;
    SpinUpTS(3).silc(ind,:) = var(i12_3,j12_3,:) ; SpinUpTS(4).silc(ind,:) = var(i12_4,j12_4,:) ;
    SpinUpTS(5).silc(ind,:) = var(i12_5,j12_5,:) ; SpinUpTS(6).silc(ind,:) = var(i12_6,j12_6,:) ;

end

save('SpinUpTS.mat','SpinUpTS')

return 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PLOTS %%%%%%%%%%%%%%%%%%

load('SpinUpTS.mat')

for i=1:6

    zeta = SpinUpTS(i).zeta ; 
    h    = [ones(1,82)*SpinUpTS(i).h25 ones(1,20)*SpinUpTS(i).h12] ; 
    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);
    z3d_v1 = permute(z3d_v1,[2 1]) ; 

    time = SpinUpTS(i).time ; 
    year = datestr(time,'yyyy') ; 
    tm4plt = time ; 
    tm4plt(26:end) = tm4plt(26:end)-datenum(str2num(year(26,:)),1,1)+time(25) ; 
    tm4plt(52:end) = tm4plt(52:end)-datenum(str2num(year(52,:)),1,1)+time(51) ;    
    tm4plt(78:end) = tm4plt(78:end)-datenum(str2num(year(78,:)),1,1)+time(77) ;  
    for z=1:NZ
        tm2d(:,z) = tm4plt ; 
    end

    deprange = [-2400 0] ; 

    figure ; 
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).temp) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).temp,[0:1:40],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([0 30]) ;
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).temp,[2 5 10 20 30],'k','ShowText','on') 
    ylim(deprange) ; ylabel('depth [m]') 
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Temperature time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Temp_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close

    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).salt) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).salt,[30:0.05:40],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([34.0 36.5]) ;
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).salt,[34 34.5 35 35.5 36 36.5],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Salinity time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Salt_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close

    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).nitr) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).nitr,[-1:1:45],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([0 40]) ;
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).nitr,[20 30 35 40],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Nitrate time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Nitr_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close


    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).oxyg) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).oxyg,[-1:2:400],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([0 300]) 
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).oxyg,[0 20 50 100 200],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Oxygen time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Oxyg_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close

    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).dic) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).dic,[1500:5:2500],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([1920 2370])
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).dic,[1900 2000 2100 2200 2300],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('DIC time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['DIC_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close

    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).alk) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).alk,[1800:5:2800],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([2270 2430])
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).alk,[2100 2200 2250 2300 2350 2400 2500 2600],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Alkalinity time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Alk_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close

    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).iron) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).iron,[-0.2e-4:0.2e-4:20e-4],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ; caxis([0 10e-4])
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).iron,[1e-4 5e-4 10e-4],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Iron time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Iron_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close

    figure ;
%    pcolor(tm2d,z3d_v1,SpinUpTS(i).silc) ; shading flat ; 
    contourf(tm2d,z3d_v1,SpinUpTS(i).silc,[0:5:300],'edgecolor','none')
    colorbar ; colormap(colormap_mine) ;  caxis([0 200])
    hold on ; contour(tm2d,z3d_v1,SpinUpTS(i).silc,[20 50 100 150],'k','ShowText','on')
    ylim(deprange) ; ylabel('depth [m]')
    xticks([tm4plt(13) tm4plt(26) tm4plt(39) tm4plt(52) tm4plt(65) tm4plt(78) tm4plt(91)]) ;
    xticklabels({year(13,:) year(26,:) year(39,:) year(52,:) year(65,:) year(78,:) year(91,:)});
    title('Silicate time series','fontsize',14) ; set(gca,'fontsize',14) ;
    fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 10 2] ;
    name_file = ['Silc_SpinUp_' num2str(i)] ;
    print(gcf,'-painters',[rep_out name_file],'-depsc') ; print(gcf,'-painters',[rep_out name_file],'-depsc') ;
    close


end
