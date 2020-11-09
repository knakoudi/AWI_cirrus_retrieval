% AWI_cirrus_retrieval: config_constrained_Klett
% *AWI (Alfred Wegener Institute)

% This routine contains the constrained Klett configuration algorithm 
% The configuration parameters of convergence range, reference profile and
% reference value are determined (see Nakoudi et al.
% (2020): An extended lidar-based cirrus cloud retrieval scheme: first
% application over an Arctic site, submitted to Optics Express // Flowchart Fig. 4 and section 2.3.2)

% Last update: 09-11-2020
% Author: Konstantina Nakoudi (konstantina.nakoudi@awi.de)

%% Prerequisites: 
% 1. The cirrus geometrical properties should already be determined 
% 2. Air density and Temperature profiles are needed for estimating the molecular extinction and backscatter coefficients
% 3. If needed, the Pr2 signals can be averaged temporally. Physically
% meaningful averaging periods can be defined by the Lanzante method (Lanzante, 1996)
%https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/(SICI)1097-0088(199611)16:11%3C1197::AID-JOC89%3E3.0.CO;2-L
% 4. A routine for the classical Klett-Fernald retrievals

% If needed select specific time bins for analysis
t1= 1; %default
t2= %your value

% Convergence Range determination
% This is a 500m zone, where the variance of the median Pr2 signal is minimum
% The Klett solution convergence will be performed within the convergence range

% function: conv_range
%% Input parameters:
% filename: yymmdd
% H: Height vector
% TimeXXX: Time vector
% PXXX: Lidar Pr2 signal
% Wvl: operating wavelength
% t1 and t2: selected time bins

%% Output parameters:
% listvalidXXX: valid cirrus profiles
% listvalidXXX_all: valid cirrus and cloud-free profiles
% listvalidXXX_final: well-correlated cirrus and cloud-free profiles
% listvalidXXX_ci: well-correlated cirrus profiles 
% convrXXX: convergence range bins 

% function [listvalidXXX,listvalidXXX_all,listvalidXXX_final,listvalidXXX_ci,convrXXX]=  conv_range(filename,H,TimeXXX,PXXX,WvlXXX,t1,t2);
% 
% % Select the high quality cirrus and cloud-free profiles
% listvalid=find(cloud_flag_final'==1);
% cloudfree = find(cloud_free_flag'==1);

%  select specific time bins
% for i= 1:length(listvalid); if listvalid(i)<t1 || listvalid(i)>t2;  listvalid(i)=nan;  end; end
% 
% % cirrus & cloud-free profs
% listvalid_all1 = union(listvalid,cloudfree);
% listvalid_all = listvalid_all1(~isnan(listvalid_all1)); 
% 
% % Determination of convergence zone 
% hstep = 66; % (66 range bins = 500 m for 7.5 m vertical resoltion \\ adapt for different vertical resolution) 
% hstepnr = fix(length(H)/hstep);
% H500 = zeros(hstepnr,1);
% P_hmedian= zeros(hstepnr,size(PXXX,2)); % median Pr2 signal within 500-m zones
% for i=1:hstepnr
%     P_hmedian(i,:) = median(PXXX(((i*hstep)-hstep+1):(i*hstep),:),1,'omitnan'); 
%     H500(i) = H((i*hstep));
% end
%
% Constrain convergence range: above incomplete overlap range and 1km below lowest Cbase (hbelcl)
% Hmin1 = 81; % [height bins] --> full overlap 
% hbelcl = min(Cbase_dyn(t1:t2))-133;
% diff1 = nan(1,length(H500));   diff2 = nan(1,length(H500)); 
% for i = 1:length(H500)
%     diff1(1,i) = abs(H(Hmin1)-H500(i));
%     [mini,wo1] =  min(diff1);   Hwo1 = H500(wo1+1);
%     diff2(1,i) = abs(H(hbelcl)-H500(i));
%     [mini,wo2] =  min(diff2);    Hwo2 = H500(wo2);
% end
% 
% % Candidate convergence range bins
% candid_r = Hmin1:hbelcl; 
% 
% % Exclude range bins with low- or mid-level clouds 
% if any(~isnan(lowcloud_flag)) || any(~isnan(midcloud_flag))
% lowclbins = min(Cbase_dyn_low) : max(Ctop_dyn_low) ;  % affected range bins 
% for j=1:length(candid_r); if any(candid_r(j) == lowclbins);  candid_r(j)=nan ; end; end
% candid_r = candid_r(~isnan(candid_r)); 
% 
% for i = 1:length(H500)
%     diff3(1,i) = abs(H(lowclbins(1))-H500(i)); 
%     [mini,wo3] =  min(diff3);   Hwo3 = H500(wo3); H500(wo3)=nan; % min low-cloud range
%     
%     diff4(1,i) = abs(H(lowclbins(end))-H500(i)); 
%     [mini,wo4] =  min(diff4);   Hwo4 = H500(wo4); H500(wo4)=nan; % max low-cloud range
%     
%     % exclude the H500 range bins between min-max
%     H500(wo3:wo4) = nan; 
% end
% end 
% 
% % Temporal correlation of signals below the cirrus cloud
% P_hmedian_tempmed = median(P_hmedian,2,'omitnan'); % temporal median of P_hmedian 
% P_tempmed = median(PXXX,2,'omitnan'); % temporal median of PXXX
% r3 = nan(size(PXXX,2),1);
% lowsigcorr = nan(size(PXXX,2),1);
% for i = listvalid_all
% ri = corrcoef(PXXX(candid_r,i),P_tempmed(candid_r ));  r3(i) = ri(1,2);
% if r3(i) < 0.98 && ~isnan(r3(i))
%    disp('Signal  profile not well-correlated: reject prof'); 
%    lowsigcorr(i) = 1;  
%    %BSR532ref(:,i)=nan; BetaAer532ref(:,i)=nan;
% else   lowsigcorr(i) = 0;
% end
% end
% 
% % Well-correlated cirrus and cloud-free profs 
% test = find(lowsigcorr==0);
% listvalid_final1 = intersect(listvalid_all',test);
% listvalid_final = listvalid_final1(~isnan(listvalid_final1)); 
% 
% % Flag profiles: 
%%  Well-correlated cirrus profiles 
% listvalid_ci1 = intersect(listvalid',test);
% for i= 1:length(listvalid_ci1); if listvalid_ci1(i) < t1 || listvalid_ci1(i) > t2;  listvalid_ci1(i)=nan; end; end
% listvalid_ci = listvalid_ci1(~isnan(listvalid_ci1)); 
% 
% % Signal temporal variance in 500m-zones
% varP_temp = nan(length(H500(wo1:wo2)),1);
% for j= (wo1+1) : wo2
%     varP_temp(j,1) = var(P_hmedian(j,listvalid_final),0,2,'omitnan');
%
%     % if a candidate range contains low- or mid-cloud, exclude it 
%     if any(~isnan(lowcloud_flag)) || any(~isnan(midcloud_flag));
%     varP_temp(wo3:wo4,1)=nan;  
%     end  
%     [mini,idmin] = (min(varP_temp));
% end
% 
% % Finally, check whether the variance is equally low at a higher range
% % If yes, select the higher range as the error due to backward Klett integration is lower at higher ranges 
% if idmin < length(varP_temp)
%     idx = length(varP_temp) - idmin;
%     for i = idx:1
%         if abs (real(log10(varP_temp(i+idx))) - real(log10(mini))) < 1
%            disp('equally low variance found above')
%            mini = i+idx
%         else
%             disp('covergence zone remains as is')
%         end
%     end 
% end
% 
% % Convergence range 
% convrXXX = ((idmin*hstep)-hstep+1):(idmin*hstep); 
% disp(['reference zone for ',num2str(Wvl),' defined'])
% 
% return
% end % function

% Reference profile and Reference value determination

%  Initial guess  Klett retrieval 

%  initial conditions
% BSR calibration value  
BSRAtFitXXX = 1.01; 
dimen=size(PXXX);  BSRAtFitXXXarr=zeros(dimen(2),1)+BSRAtFitXXX;  eps =0.05; BSRAtFitXXXerr= eps./2; 

% Calibration range 
dH = dh;      Hmin1 = 81; %full overlap     
hbelcl = min(Cbase_dyn)-65; % 1km below min(Cbase)
FitRangeC=[13000 15000];    FitSel=find(H>FitRangeC(1) & H<FitRangeC(2)); 

% check SNR at Calibration range
SNRfitXXX = mean(SNRXXX(FitSel,:),1,'omitnan'); 
for i=1:n_final; if SNRfitXXX(i)<3; disp('WARNING! low SNR at calibration range'); end; end; 

% check SNR at Ctop
SNRCtopXXX= nan(n_final,1); for i=1:n_final ; SNRCtopXXX(i) = mean(SNRXXX(Ctop_dyn(i)-3:Ctop_dyn(i),i),1);  end 
for i=1:n_final; if SNRCtopXXX(i)<3; disp('WARNING! low SNR at Ctop'); end; end; 

% Check if Ctop < FitRange
for t=1:n_final;if H(Ctop_dyn(t,1))>FitRangeC(1);disp('Ci inside Calibr. Window'); end;Ctop_max=max(Ctop_dyn); Ctop_max=H(Ctop_max); end

% Initial guess LR array: LRci and LRaer zones 
LRXXXarr=ones(size(PXXX)).*val_LRXXX;  %Predefine LR array 
LRXXXarrerr = ones(size(LRXXXarr))+15;  % LR error
BeRaXXX_avg=Density(:,1:n_final).*raybckwq (WvlXXX,'p','p', Temp(:,1:n_final), Density(:,1:n_final)); %molecular bsc

% Use LRci only in WCT-derived cirrus bins 
WCTcirrusarr= zeros(size(PXXX_avg));
for i= 1: n_final
 if ~isnan(Cbase_dyn(i)) &&  ~isnan(Ctop_dyn(i))  
 WCTcirrus = Cbase_dyn(i): Ctop_dyn(i); 
 WCTcirrusarr(WCTcirrus,i) = WCTcirrus; LRXXXarr(WCTcirrus,i) = val_LRXXXWo+0; 
 end
end

% classical Klett-Fernald retrieval
% BetaAerXXXKlett_classic: particulate backscatter coefficient at XXX nm
% BSRXXXKlett_classic: backscatter ratio at XXX nm
BetaAerXXXKlett_classic=zeros(size(PXXX));  BSRXXXKlett_classic =zeros(size(PXXX)); 
Perr= PXXX./SNRXXX; % signal noise
for i= 1 : n_final 
FitSel = find(H>FitRangeC(1)& H<FitRangeC(2)); %FitSel=Selected points for applying the Bound.Cond. 
% // replace with your routine for classical Klett-Fernald retrievals
%[Beta] = klettinv_ableit4( BSRAtFitXXXarr(i), FitRangeC,H,PXXX(:,i),Perr(:,i),LRXXXarr(:,i),AlRayXXX_avg(:,i),BeRaXXX_avg(:,i)); //
Betaaer(:)=Beta-BeRaXXX_avg(:,i);   Btemp(:)=Beta./BeRaXXX_avg(:,i);
BSRXXXKlett_classic(:,i)=Btemp; 
end

% Initial guess COD time-series
CODXXX_first_avg = nan(size(PXXX,2),1);
for i=1:n_final; 
if ~isnan(Cbase_dyn(i)) &&  ~isnan(Ctop_dyn(i)) 
CODXXX_first_avg(i)=trapz(H(Cbase_dyn(i):Ctop_dyn(i)),BetaAerXXXKlett_classic(Cbase_dyn(i):Ctop_dyn(i),i)...
.*LRXXXarr(Cbase_dyn(i):Ctop_dyn(i),i),1); 
end
end

% Initial guess integrated backscatter (b_int) time-series
BetaintXXX_first_avg = nan(size(PXXX,2),1);
for i= 1: n_final; 
if ~isnan(Cbase_dyn(i)) &&  ~isnan(Ctop_dyn(i))  
BetaintXXX_first_avg(i) = trapz(H(Cbase_dyn(i):Ctop_dyn(i)),BetaAerXXXKlett_classic(Cbase_dyn(i):Ctop_dyn(i),i)); 
end
end
figure(); subplot(2,1,1); plot(CODXXX_first_avg,'.'); subplot(2,1,2); plot(BetaintXXX_first_avg,'.')

% Reference profile: profile with minium b_int within the cirrus range
for i=1:length(BetaintXXX_first_avg); if BetaintXXX_first_avg(i)<=0; BetaintXXX_first_avg(i)=nan; end; end; igood = find(negXXX_avg==0);
[BintXXXavg_min,minidx] = min(min(BetaintXXX_first_avg(igood),1,'omitnan')); minidxXXX = igood(minidx)

%  BSRref (median and standard deviation): BSR within the convergence range of reference profile
BSRXXXref = median(median(BSRXXXKlettold(convrXXX,minidxXXX),'omitnan'))
BSRXXXref_std = std(std(BSRXXXKlettold(convrXXX,minidxXXX),'omitnan'))
figure(); plot(BSRXXXKlettold(:,minidxXXX),H)

% Now we are ready to perform the constrained Klett retrievals
