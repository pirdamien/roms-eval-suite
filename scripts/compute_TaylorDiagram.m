
clear all 

% For each variable, here is a vector with 
% [relSTD relRMS COR obsSTD RMSE] 

% temp 
TEM_srfmean = [0.989 0.077 0.998  6.852  0.528] ;
TEM_200mean = [0.987 0.193 0.983  4.908  0.948] ;
TEM_srfmvar = [1.028 0.197 0.983  1.166  0.230] ;
TEM_200mvar = [0.939 0.772 0.718  0.248  0.191] ;
% salt 
SAL_srfmean = [0.997 0.280 0.961  1.026  0.288] ;
SAL_200mean = [1.005 0.215 0.978  0.635  0.136] ;
SAL_srfmvar = [1.082 0.806 0.705  0.161  0.130] ;
SAL_200mvar = [0.799 0.901 0.615  0.026  0.023] ;
% MLD 
MLD_mean     = [1.083 0.769 0.731 11.796  9.069] ;
MLD_wint     = [1.018 0.593 0.830 23.819 14.134] ;
MLD_summ     = [1.016 0.450 0.900 28.383 12.783] ;
MLD_mvar     = [1.171 0.907 0.757 10.927  9.917] ;
MLD_wintp100 = [1.018 4.170 0.830 23.819 99.310] ; %MLD+100 to test ; 
% SSH 
SSH_mean    = [0.976 0.291 0.984  0.284  0.082] ;
SSH_mvar    = [0.976 0.565 0.844  0.035  0.020] ;

% NO3 
NITR_srfmean = [0.768 0.818 0.668  3.844  3.147] ;
NITR_200mean = [0.865 0.389 0.930 10.217  3.972] ;
% PO4 surf and 200m , mean and monthly std
PHOS_srfmean = [0.692 0.788 0.812  0.340  0.268] ;
PHOS_200mean = [0.756 0.404 0.940  0.766  0.309] ;
%  O2 surf and 200m , mean and monthly std
OXYG_srfmean = [0.992 0.198 0.993 32.909  6.510] ;
OXYG_200mean = [0.833 0.596 0.926 66.728 39.757] ;
% SIO2 surf and 200m , mean and monthly std
SILC_srfmean = [0.480 0.768 0.697  6.837  5.253] ;
SILC_200mean = [1.120 0.640 0.852 18.129 11.600] ;
%  Fe surf and 200m , mean and monthly std
IRON_srfmean = [3.889 4.185 0.128  4.53e-5  0.0001896] ;
IRON_200mean = [1.411 1.076 0.676  0.0001738  0.00018711] ;
% N2O surf and 200m , mean and monthly std
NTOX_srfmean = [0.823 0.558 0.940  0.002  0.001] ;
NTOX_200mean = [0.478 0.864 0.835  0.012  0.010] ;

% DIC surf and 200m , mean
DIC_srfmean = [1.157 1.164 0.507 54.644 63.595] ;
DIC_200mean = [1.117 0.663 0.873 71.483 47.358] ;
% Alk surf and 200m , mean
ALK_srfmean = [0.367 1.191 -0.213 53.789 64.059] ;
ALK_200mean = [0.557 1.286 -0.305 30.612 39.383] ;
% Omg surf and 200m , mean
OMG_srfmean = [0.984 0.444 0.930  0.617  0.274] ;
OMG_200mean = [0.832 0.371 0.960  0.756  0.280] ;
% pH surf and 200m , mean
PHC_srfmean = [0.686 0.744 0.670  0.036  0.027] ;
PHC_200mean = [0.638 0.652 0.924  0.154  0.100] ;

% [relSTD relRMS COR obsSTD RMSE]
% CHL surf , mean and monthly std
CHL_srfmean = [0.830 0.864 0.569  0.403  0.348] ;
CHL_srfmvar = [1.066 1.005 0.537  0.238  0.239] ;
% FGCo2 , mean and std
CO2flx_srfmean = [1.135 0.764 0.751  3.907e-5  2.984e-5] ;
CO2flx_srfmvar = [1.395 1.138 0.720  1.905e-5  2.168e-5] ;
% FGN2O , mean and std
N2Oflx_srfmean = [0.342 0.895 0.483  2.292e-8  2.051e-8] ;
N2Oflx_srfmvar = [0.766 0.892 0.517  1.150e-8  1.025e-8] ;
% NPP , mean and std
NPP_VGPMmean = [0.645 0.914 0.455 261.82 239.34] ;
NPP_CBPMmean = [0.595 0.873 0.528 284.01 247.74] ;
NPP_CAFEmean = [1.237 1.230 0.591 136.58 165.25] ;

NPP_VGPMmvar = [0.694 0.726 0.688 136.91  99.37] ;
NPP_CBPMmvar = [0.742 0.688 0.747 128.13  88.13] ;
NPP_CAFEmvar = [0.381 1.423 0.296  60.09  85.52] ;

% POC , mean and std
POC_mean1 = [0.547 0.834 0.602 3.695 3.082] ; % Clements 2023
POC_mean2 = [1.361 1.257 0.698 1.570 1.973] ; % Nowicki 2022

return

%%% Physics Diagram
%VAR2plot = [TEM_srfmean' TEM_200mean' SAL_srfmean' SAL_200mean' ...
%            TEM_srfmvar' TEM_200mvar' SAL_srfmvar' SAL_200mvar' ...
%            SSH_mean' SSH_mvar' ...
%            MLD_mean' MLD_wint' MLD_summ' MLD_mvar']' ; 
%VAR2plot = [TEM_srfmean' TEM_200mean' SAL_srfmean' SAL_200mean' ...
%            TEM_srfmvar' TEM_200mvar' SAL_srfmvar' SAL_200mvar' ...
%            SSH_mean' SSH_mvar' ...
%            MLD_wint' MLD_summ']' ;
VAR2plot = [TEM_srfmean' TEM_200mean' SAL_srfmean' SAL_200mean' ...
            SSH_mean' SSH_mvar' ...
            MLD_wint' MLD_summ']' ;
[pp tt axl] = MYtaylordiag(squeeze(VAR2plot(:,1)),squeeze(VAR2plot(:,2)),squeeze(VAR2plot(:,3)));
ii=1;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ; 
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=2;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=3;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=4;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=5;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=6;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=7;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = '>' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4660, 0.6740, 0.1880] ;
ii=8;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = '<' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4660, 0.6740, 0.1880] ;
legend('','','mean Temp srf','mean Temp 200m','mean Salt srf','mean Salt 200m' , ...
            'mean SSH','m.var SSH','mean MLD winter','mean MLD summer')
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = ['FigArpae/TaylorDiag_Physics'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
close


%%% Nutrient & Oxygen Diagram
VAR2plot = [NITR_srfmean' NITR_200mean' ...
            PHOS_srfmean' PHOS_200mean' ...
            OXYG_srfmean' OXYG_200mean' ...
            SILC_srfmean' SILC_200mean' ...
            NTOX_srfmean' NTOX_200mean']' ;
[pp tt axl] = MYtaylordiag(squeeze(VAR2plot(:,1)),squeeze(VAR2plot(:,2)),squeeze(VAR2plot(:,3)));
ii=1;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=2;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=3;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=4;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=5;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=6;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=7;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4660, 0.6740, 0.1880] ;
ii=8;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4660, 0.6740, 0.1880] ;
ii=9;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4940, 0.1840, 0.5560] ;
ii=10;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4940, 0.1840, 0.5560] ;
legend('','','mean NO_3 srf','mean NO_3 200m','mean PO_4 srf','mean PO_4 200m' , ...
             'mean O_2 srf','mean O_2 200m','mean SiO_3 srf','mean SiO_3 200m' , ...
             'mean N_2O srf','mean N_2O 200m')
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = ['FigArpae/TaylorDiag_Nutrients'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
close



%%% Carbon Chemistry Diagram
VAR2plot = [DIC_srfmean' DIC_200mean' ...
            OMG_srfmean' OMG_200mean' ...
            PHC_srfmean' PHC_200mean']' ;
[pp tt axl] = MYtaylordiag(squeeze(VAR2plot(:,1)),squeeze(VAR2plot(:,2)),squeeze(VAR2plot(:,3)));
ii=1;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=2;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=3;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=4;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=5;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=6;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'diamond' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
legend('','','mean DIC srf','mean DIC 200m','mean \Omega_{Ar} srf','mean \Omega_{Ar} 200m' , ...
             'mean pH srf','mean pH 200m')
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = ['FigArpae/TaylorDiag_CarbonSys'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
close


%%% BGC Cycling Diagram
VAR2plot = [CHL_srfmean' CHL_srfmvar' ...
            NPP_VGPMmean' NPP_VGPMmvar' ...
            NPP_CBPMmean' NPP_CBPMmvar' ...
            CO2flx_srfmean' CO2flx_srfmvar' ...
            N2Oflx_srfmean' N2Oflx_srfmvar' ...
            POC_mean1' POC_mean2' ...
            NPP_CAFEmean' NPP_CAFEmvar']' ;
[pp tt axl] = MYtaylordiag(squeeze(VAR2plot(:,1)),squeeze(VAR2plot(:,2)),squeeze(VAR2plot(:,3)));
ii=1;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=2;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.9290, 0.6940, 0.1250] ;
ii=3;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=4;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.8500, 0.3250, 0.0980] ;
ii=5;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.6350, 0.0780, 0.1840] ;
ii=6;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.6350, 0.0780, 0.1840] ;
ii=7;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4660, 0.6740, 0.1880] ;
ii=8;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4660, 0.6740, 0.1880] ;
ii=9;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4940, 0.1840, 0.5560] ;
ii=10;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.4940, 0.1840, 0.5560] ;
ii=11;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.3010, 0.7450, 0.9330] ;
ii=12;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0, 0.4470, 0.7410] ;
ii=13;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'square' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.93, 0.16, 0.31] ;
ii=14;
set(tt(ii),'fontsize',9,'fontweight','bold','String','')
set(pp(ii),'markersize',8) ; pp(ii).Marker = 'hexagram' ;
pp(ii).MarkerEdgeColor = 'k' ; pp(ii).MarkerFaceColor = [0.93, 0.16, 0.31] ;
legend('','','mean Chl','m.var Chl','mean NPP VGPM','m.var NPP VGPM' , ...
            'mean NPP CBPM','m.var NPP CBPM','mean CO_2 flx','m.var CO_2 flx' , ...
            'mean N_2O flx','m.var N_2O flx', 'mean POC vs Clements 2023' , 'mean POC vs Nowicki 2022' , ...
            'mean NPP CAFE','m.var NPP CAFE')
set(gca,'fontsize',14) ;
fig=gcf; fig.PaperUnits = 'inches' ; fig.PaperPosition = [0 0 5 5] ;
name_file = ['FigArpae/TaylorDiag_BGCflx'] ;
print(gcf,'-painters',name_file,'-depsc') ; print(gcf,'-painters',name_file,'-depsc') ;
close








% 
%VAR2plot = [MLD_wint' MLD_wintp100']' ;
%[pp tt axl] = MYtaylordiag(squeeze(VAR2plot(:,1)),squeeze(VAR2plot(:,2)),squeeze(VAR2plot(:,3)));



