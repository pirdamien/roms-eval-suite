
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     MLD Model                         %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid_file = '/data/project3/pdamien/ROMS_pdamien/config/pacmed12km/grid/pacmed_12km_grd.nc'
%rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/Y*/'
%list = [rep_in 'pacmed_avg.*.nc'];

rep_in = '/data/project3/pdamien/ROMS_outputs/PACMED12KM/mean_2000_2005/'
file = [rep_in 'pacmed_avg.nc'];

rep_out = './Fig/'

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


%list_file = dir(list) ;
%list_file = list_file(12:83) ;

%for t=1:length(list_file)
for t=1:12

%    file = [list_file(t).folder '/' list_file(t).name]
    zeta   = squeeze(ncread(file,'zeta',[1 1 t],[inf inf 1])) ;
    temper = squeeze(ncread(file,'temp',[1 1 1 t],[inf inf inf 1])) ;
    salini = squeeze(ncread(file,'salt',[1 1 1 t],[inf inf inf 1])) ;
    [z3d_v1,Cw] = zlevs4(h',zeta', theta_s, theta_b, hc, NZ , 'r',sc_type);
    z3d_v1 = permute(z3d_v1,[3 2 1]) ; 

    MLD1=zeta.*NaN ; MLD2=zeta.*NaN ;
    for i=1:size(zeta,1)
%    tic
    disp([num2str(t) ' - ' num2str(i)])
    for j=1:size(zeta,2)
    %for i=1:5
    %for j=1:5
        if ( h(i,j)>20 )
        if ( mask(i,j)==1 )
             temp = squeeze(temper(i,j,:)) ; temp=flipud(temp) ; 
             sal  = squeeze(salini(i,j,:)) ; sal =flipud(sal ) ;
             depth= squeeze(z3d_v1(i,j,:)) ; depth=-flipud(depth) ; 
%             pres = sw_pres(depth,lat(i,j)) ; 
             pres = depth ;
             MLD1(i,j) = findmldGOSML(temp,sal,pres)   ; %%% Holte and Talley 2009 routine ,density algo to compare with GOSML 
%             MLD2(i,j) = findmldJAMSTEC(temp,sal,pres) ; %%% method MLD 1 (min(dens,temp)) criteron suited for JAMSTEC
        end
        end
    end 
%    toc
    end

MLD_mod.MLD(:,:,t)   = MLD1 ;
%MLD_mod.MLD(:,:,t) = MLD2 ;
%MLD_mod.year(t)  = str2num(list_file(t).folder(end-6:end-3)) ;
%MLD_mod.month(t) = str2num(list_file(t).folder(end-1:end)) ;

save('MLD_model_GOSML.mat','MLD_mod','-v7.3')

end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



