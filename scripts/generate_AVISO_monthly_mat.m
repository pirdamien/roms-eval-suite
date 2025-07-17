
clear all

file = '/data/project1/data/AVISO/data/2000/dt_global_allsat_phy_l4_20000101_20190101.nc' ;   
lon_aviso = ncread(file,'longitude') ; 
lat_aviso = ncread(file,'latitude');

ind = 0 ;

for an=2000:2005
    
    an

for m=1:12   
    
    if m<10
        month = ['0' num2str(m)];
    else
        month = num2str(m);
    end
    
    rep = dir(['/data/project1/data/AVISO/data/' num2str(an) '/*' num2str(an) num2str(month) '*.nc']) ;
    
    for d=1:length(rep(:,1)) 
        
    file = ['/data/project1/data/AVISO/data/' num2str(an) '/' rep(d,1).name];  
    ncid = netcdf.open(file) ;
    varid_sla = netcdf.inqVarID(ncid,'sla') ;
    sla(:,:,d) = double(netcdf.getVar(ncid,varid_sla)) ;
    varid_adt = netcdf.inqVarID(ncid,'adt') ;
    adt(:,:,d) = double(netcdf.getVar(ncid,varid_adt)) ;
    netcdf.close(ncid)
    
    end
  
    ind=ind+1;
    test = nanmean(sla,3) ; test((test)<-2e5)=NaN ;
    sla_monthly(:,:,ind) = test*0.0001 ;
    test = nanmean(adt,3) ; test((test)<-2e5)=NaN ;
    adt_monthly(:,:,ind) = test*0.0001 ;
    date_monthly(ind) = datenum(an,m,1) ;
    
    clear sla
    clear adt
    
end

end

monthly_AVISO.lon_aviso = lon_aviso ;
monthly_AVISO.lat_aviso = lat_aviso ;
monthly_AVISO.sla_monthly = sla_monthly ;
monthly_AVISO.adt_monthly = adt_monthly ;
monthly_AVISO.date_monthly = date_monthly ;

save('monthly_AVISO_2000_2005.mat','monthly_AVISO','-v7.3')


return

figure
pcolor(squeeze(nanmean(monthly_AVISO.sla_monthly,3))') ; shading flat ; colorbar
caxis([-0.1 0.1])

test = permute(squeeze(std(permute(monthly_AVISO.sla_monthly,[3 2 1]))),[2 1]);
figure
pcolor(test') ; shading flat ; colorbar
caxis([0 0.4])


figure
pcolor(squeeze(nanmean(monthly_AVISO.adt_monthly,3))') ; shading flat ; colorbar
caxis([-1.5 1.5])


