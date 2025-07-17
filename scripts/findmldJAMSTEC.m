function [mixeddp] = findmldJAMSTEC(temp,sal,pres)

% Calculate the index of the reference value
m = length(sal); 
starti = min(find((pres-10).^2==min((pres-10).^2)));  
pres = pres(starti:m);
sal = sal(starti:m);
temp = temp(starti:m);
starti = 1;
m = length(sal); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the potential density anomaly, with a reference pressure of 0 
pden = sw_pden(sal,temp,pres,0)-1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Search for the first level that exceeds the potential density threshold
mldepthdens = m;
for j = starti:m
    if abs(pden(starti)-pden(j))>.03 
        mldepthdens = j;
        break;
    end
end
% Interpolate to exactly match the potential density threshold
clear pdenseg presseg presinterp pdenthreshold
presseg = [pres(mldepthdens-1) pres(mldepthdens)];
pdenseg = [pden(starti)-pden(mldepthdens-1) pden(starti) - pden(mldepthdens)];
P = polyfit(presseg,pdenseg,1);
presinterp = presseg(1):.5:presseg(2);
pdenthreshold = polyval(P,presinterp);
% The potential density threshold MLD value:
mldepthdens = presinterp(max(find(abs(pdenthreshold)<.03)));

%% Search for the first level that exceeds the temperature threshold
%mldepthptmp = m;
%for j = starti:m
%    if abs(temp(starti)-temp(j))>.2
%        mldepthptmp = j;
%        break;
%    end
%end
%% Interpolate to exactly match the temperature threshold
%clear tempseg presseg presinterp tempthreshold
%presseg = [pres(mldepthptmp-1) pres(mldepthptmp)];
%tempseg = [temp(starti)-temp(mldepthptmp-1) temp(starti) - temp(mldepthptmp)];
%P = polyfit(presseg,tempseg,1);
%presinterp = presseg(1):.5:presseg(2);
%tempthreshold = polyval(P,presinterp);
%% The temperature threshold MLD value:
%mldepthptmp = presinterp(max(find(abs(tempthreshold)<.2)));

mixeddp = mldepthdens ;
%mixeddp = min([mldepthdens mldepthptmp]) ; 

return
