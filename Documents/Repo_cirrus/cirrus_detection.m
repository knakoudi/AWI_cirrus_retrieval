% AWI_cirrus_retrieval: cirrus_detection
% *AWI (Alfred Wegener Institute)
% This routine contains the cirrus detection algorithm via the dynamic WCT method 
% The |WCT/std| and SNR related criteria are applied (from Nakoudi et al.
% (2020): An extended lidar-based cirrus cloud retrieval scheme: first
% application over an Arctic site, submitted to Optics Express)

% //This algorithm has been developed and tested on lidar data obtained over
% the AWIPEV research base, Ny-Alesund, Svalbard
% Some parameters may need tuning if the algorithm is applied on data from
% mid-latitudes or the tropics//

% Last update: 09-11-2020
% Author: Konstantina Nakoudi (konstantina.nakoudi@awi.de)

% Prerequisites:           
%% The WCT profile needs to be calculated by Equation (6): Kokkalis et al. (2020) Boundary-Layer Meteorology 
% https://doi.org/10.1007/s10546-020-00514-z

%% It is recommended to screen low- and mid-level clouds (e.g. via the WCT method) and check the lidar signal quality aloft

%% //The lidar range-corrected signal (Pr2) should be provided in dimensions [range bins, time bins]. 
%% Otherwise, the routine should be adapted accordingly //

%% Cirrus Detection %%

% Configuration and thresholds

% WCT dilation parameter 
n = 12;  % 90 m for a system with 7.5 m vertical resolution as the Koldewey Aerosol Raman Lidar // adapt accordingly
step = 7.5 % [m] // lidar system vertical resolution  // adapt accordingly
alpha = n*step; % dilation in [m]
ihalf = (alpha/(2*step));   % should be ~=1 and integer

% SNR ratio thresholds 
% This threshold compares the cirrus marginal layers' SNR to the adjacent areas' SNR 
SNRratio_thresXXX = 1.1;  % for Cbase of 532S channel // modify for different channels/lidar systems
SNRratio_rev_thresXXX = 1.2; % for Ctop of 532S channel //modify for different channels/lidar systems

% function: WCT_dyn_cbase
% This function determines the Cbase using dynamic WCT thresholds
%% Input parameters:
% Pr2_XXX_norm - normalized Pr2 signal for comparability with previous studies
% H - range vector
% SNRXXX - Signal-to-noise ratio of XXX channel (perpendicular polarization detection capability)
% WCTXXX_d90 - WCT profile of XXX signal, calculated with dilation alpha = 90 m 
% ihalf - vertical moving window for WCT, std and SNR calculations
% cloud_flag - low- and mid-level cloud screening 
% SNRratio_thresXXX - SNR ratio threshold for the Cbase at the XXX channel

%% Output parameters:
% Pr2_std - standard deviation of Pr2 signal within alpha/2 below each range bin
% WCT_STDXXX_RATIO - ratio of WCT to std within alpha/2 below each range bin
% Cbase_dyn - Cbase time-series
% cbase_flag - //1 for sucessful Cbase detection // 0 for no Cbase detection
% WCT_cbase_dyn - WCT value at Cbase
% WCT_STD_RATIO_cbase_dyn - ratio of WCT to std at Cbase
% SNRmed_in - median SNR at marginal cirrus layer above Cbase
% SNRmed_out - median SNR at adjacent area below Cbase
% SNRratioXXX - SNR ratio at Cbase

% function [Pr2_std, WCT_STD_RATIO, Cbase_dyn,cbase_flag,WCT_cbase_dyn,WCT_STD_RATIO_cbase_dyn,...SNRmed_in,SNRmed_out,SNRratio]...
%     = WCT_dyn_cbase(Pr2_norm, H,SNR, WCT,ihalf,cloud_flag,SNRratio_thres)
% 
% %% Calculate Signal std in overlapping fragments of alpha/2 
% %% Same bins as WCT 
% %% For Cbase we use the ihalf bins below each bin (we look for upward signal changes)
% disp('Calculating Cbase dynamically ')
% %% Search for Cbase above 5km
% cuth = find(H< (5000));   Pr2_norm(cuth,:) = nan; % avoid low- and mid-level clouds
% cuth2 = find(H > 12500); % cut signal above 12 km (avoid noise) % adjust accordingly for mid-latitudes and tropics
% Pr2_norm(cuth2,:) = nan;
% 
% % Standard deviation of Pr2 signal within a/2 below each range bin  
% Pr2_std = nan(size(Pr2_norm));
% start=ihalf+1;
% ende= size(H,1)-ihalf-1;
%  for t = 1: size(Pr2_norm,2) %time steps
%  for z =start:ende %height interval
%  Pr2_std(z,t) = std(Pr2_norm(((z-ihalf):z-1),t),'omitnan');
%  end
%  end
%  
%  % SNR  within a/2 below each range bin  (candidate outer zone)
%  SNRmed_out =  nan(size(Pr2_norm));
%  for t = 1: size(Pr2_norm,2) %tsteps
%  for z =start:ende; %height interval
%  SNRmed_out(z,t) = median(SNR(((z-ihalf):z-1),t),'omitnan');
%  end
%  end
%  
%  % SNR within a/2 above each range bin (candidate marginal zone)
%  SNRmed_in =  nan(size(Pr2_norm));
%  for t = 1: size(Pr2_norm,2) %tsteps
%  for z =start:ende; %height interval
%  SNRmed_in(z,t) = median(SNR((z+1:(z+ihalf)),t),'omitnan');
%  end
%  end
%  
%  % SNR ratio 
%  SNRratio = SNRmed_in ./ SNRmed_out ; 
%  
%  % At each range bin calculate WCT/std (within a/2 below)
%  %% At Cbase (signal increase) the ratio will be negative
% WCT_STD_RATIO = nan(size(Pr2_norm)); 
% for z = 2 : size(Pr2_norm,1) %vertical bins
% WCT_STD_RATIO(z,:) = WCT(z,:) ./ Pr2_std(z,:) ; 
% end
% 
% % Cbase selection
% % 1st Condition: Search upwards for WCT/std < -1  (localmin)
% %% These are local minimums that indicate signal increase higher than the std of signal
% localmin_ratio = nan(size(Pr2_norm));
% for t = 1: size(Pr2_norm,2) %tsteps   
% q = find( WCT_STD_RATIO(:,t) < -1) ;  %index
% localmin_ratio(q,t) = WCT_STD_RATIO(q,t);    
% end
% 
% %% Find WCT @ candidate Cbase range bins
% localmin_idx = nan(size(Pr2_norm));
% localmin_idx (isnan(localmin_ratio)==1) = 0; % no possible Cbase
% localmin_idx (isnan(localmin_ratio)==0) = 1; % possible Cbase
% 
% WCT_localmin = nan(size(Pr2_norm));
% for t = 1: size(Pr2_norm,2) %tsteps  
% q2 =  find(localmin_idx(:,t)==1)  ;  %index   
% WCT_localmin(q2,t) =  WCT(q2,t); 
% end
% 
% WCT_cbase_all = nan(size(Pr2_norm));
% Cbase_possible = zeros(size(Pr2_norm));
% for t = 1: size(Pr2_norm,2) %tsteps      
% q3 =  find(isnan(WCT_localmin(:,t))== 0) ; 
% WCT_cbase_all(q3,t) =  WCT_localmin(q3,t); 
% Cbase_possible(q3-1,t) = 1;
% end
% 
% %% 2nd Condition: SNR > 2 at five consecutive bins above possible Cbase
% nbins = 5;
% for t = 1: size(Pr2_norm,2) %tsteps    
% q44 = find(Cbase_possible(:,t) == 1); %index of possible Cbases
% if isempty(q44)==0 %&& length(q4)>1
% for j = 1: length(q44)    
% y1(j,:) = q44(j): q44(j) + nbins ;
% snr(:,j)= SNR(y1(j,:),t); %  for each q4(possible Cbase)
% 
% if any(snr(:,j) < 2)==1; Cbase_possible(q44(j),t) = 0;
% else
%     Cbase_possible(q44(j),t) = 1;
% end  % snr
% end % q4 size
% end % isempty(q4)
% end % tsteps
% 
% %% 3rd condition: SNR ratio threshold
% for t = 1: size(Pr2_norm,2) %tsteps    
% q4 = find(Cbase_possible(:,t) == 1); %index of possible Cbases
% if isempty(q4)==0 %&& length(q4)>1
% for j = 1: length(q4)    
%     
% snr_all(:,j)= SNRratio(q4(j,:),t); %  for each q4(possible Cbase)
% if snr_all(:,j) > SNRratio_thres ; Cbase_possible(q4(j),t) = 1;
% else
%     Cbase_possible(q4(j),t) = 0;
% end  % snr
% end % q4 size
% end % isempty(q4)
% end % tsteps
% 
% %% Select the nearest range bin that fullfils the three conditions
% WCT_cbase_dyn = nan(size(Pr2_norm,2),1); % WCT @ base
% Cbase_dyn = nan(size(Pr2_norm,2),1);  % Cbase bin
% 
% for t = 1: size(Pr2_norm,2) %tsteps      
% q5 =  find(Cbase_possible(:,t)==1,1,'first') ;  %index  
% if isempty(q5)==1; continue; end % no cloud detected
% WCT_cbase_dyn(t,:) =  WCT_cbase_all(q5+1,t); % WCT@CCbase
% Cbase_dyn(t,:) = q5; %Cbase  
% end
% 
% % Ancillary parameters
% %% WCT at Cbase
% WCT_cbase_dyn = nan(size(Cbase_dyn));
% for t = 1: size(Pr2_norm,2) %tsteps
%     if isnan(Cbase_dyn(t))==1; continue; end 
% WCT_cbase_dyn(t) =  WCT(Cbase_dyn(t)+1,t);   
% end
% 
% %% WCT/std at Cbase
% WCT_STD_RATIO_cbase_dyn =  nan(size(Cbase_dyn));
% for t = 1: size(Pr2_norm,2) %tsteps
%     if isnan(Cbase_dyn(t))==1; continue; end 
% WCT_STD_RATIO_cbase_dyn (t) =  WCT_STD_RATIO (Cbase_dyn(t)+1,t);   
% end
% 
% %% Cloud base flag
% cbase_flag = nan(size(Pr2_norm,2),1); 
% cbase_flag(isnan(Cbase_dyn)==0) = 1; % Cbase was sucessfully detected
% cbase_flag(isnan(Cbase_dyn)==1) = 0; % no Cbase detected
% end

% function: WCT_dyn_ctop
% This function determines the Ctop (calculated from the reversed signal: top to bottom)
%% Input parameters:
% Pr2__norm_rvs - reversed normalized Pr2 signal for comparability with previous studies
% H - height vetor
% SNRrev - reversed SNR of XXX channel
% WCT_dnw - WCT profile of reversed XXX signal, calculated with dilation alpha = 90 m
% ihalf - vertical moving window for WCT, std and SNR calculations
% cloud_flag - low- and mid-level cloud screening 
% SNRratio_rev_thres -  SNR ratio threshold for the Ctop at the XXX channel

%% Output parameters:
% Pr2_std_for_ctop - standard deviation of Pr2 signal within alpha/2 below each range bin
% Ctop_dyn -  Ctop time-series
% ctop_flag - //1 for sucessful Ctop detection // 0 for no Ctop detection
% WCT_defliped - defliped WCT profile
% WCT_STD_defliped - defliped WCT/std profile
% WCT_ctop_defliped - WCT value at Ctop
% WCT_STD_ctop_defliped - ratio of WCT to std at Ctop
% SNRratio_defliped - SNR ratio at Ctop

% function [Pr2_std_for_ctop, Ctop_dyn,ctop_flag,...
%       WCT_defliped,WCT_STD_defliped,WCT_ctop_defliped,WCT_STD_ctop_defliped,~, SNRratio_defliped]...
%     = WCT_dyn_ctop(Pr2__norm_rvs, H, SNRrev , WCT_dnw,ihalf,cloud_flag,SNRratio_rev_thres)
% 
% % ! ! ! Pr2__norm_rvs and WCT_dnw are reversed ! ! ! 
% %       (1st element --> end of signal)
% disp('Calculating Ctop dynamically above 5km')
% %% Search for Ctop above 5km
% cuth = find(H< (5000));
% Pr2__norm_rvs(cuth,:) = nan;
% cuth2 = find(H > 12000); % cut signal above 12 km (avoid noise)
% Pr2__norm_rvs(cuth2,:) = nan;
% 
% % Standard deviation of Pr2 signal within a/2 above below range bin (reversed signal) 
% Pr2_std_for_ctop = nan(size(Pr2__norm_rvs));
% start=ihalf+1;
% ende= size(H,1)-ihalf-1;
%  for t = 1: size(Pr2__norm_rvs,2) %tsteps
%  for z =start:ende; %height interval
%  Pr2_std_for_ctop(z,t) = std(Pr2__norm_rvs(((z-ihalf):z-1),t),'omitnan');
%  end
%  end
%  
%  % SNR  within a/2 below each range bin  (candidate outer zone)
%  SNRmed_out =  nan(size(Pr2__norm_rvs));
%  for t = 1: size(Pr2__norm_rvs,2) %tsteps
%  for z =start:ende; %height interval
%  SNRmed_out(z,t) = median(SNRrev(((z-ihalf):z-1),t),'omitnan');
%  end
%  end
%  
%  % SNR  within a/2 above each range bin  (candidate marginal zone)
%  SNRmed_in =  nan(size(Pr2__norm_rvs));
%  for t = 1: size(Pr2__norm_rvs,2) %tsteps
%  for z =start:ende; %height interval
%  SNRmed_in(z,t) = median(SNRrev((z+1:(z+ihalf)),t),'omitnan');
%  end
%  end
%  
%  % SNR ratio
%  SNRratio_rev = SNRmed_in ./ SNRmed_out ; 
%  
%  % At each range bin calculate WCT/std (within a/2 below)
%  %% At Ctop (signal increase) the ratio will be negative
% WCT_STD_RATIO_for_ctop = nan(size(Pr2__norm_rvs)); 
% for z = 2 : size(Pr2__norm_rvs,1) %vertical bins
% WCT_STD_RATIO_for_ctop(z,:) = WCT_dnw(z,:) ./ Pr2_std_for_ctop(z,:) ; 
% end
% 
% % Ctop selection
% % 1st Condition: Search upwards for WCT/std < -1  (localmin)
% %% These are local minimums that indicate signal increase higher than the std of signal
% localmin_ratio = nan(size(Pr2__norm_rvs));
% for t = 1: size(Pr2__norm_rvs,2) %tsteps 
%   clear q
% q = find( WCT_STD_RATIO_for_ctop(:,t) < - 1 ) ;  %index
% localmin_ratio(q,t) = WCT_STD_RATIO_for_ctop(q,t);    
% end
% 
% %% Find WCT @ all localmin positions (@ localmin_idx)
% localmin_idx = nan(size(Pr2__norm_rvs));
% localmin_idx (isnan(localmin_ratio)==1) = 0; % no possible Ctop
% localmin_idx (isnan(localmin_ratio)==0) = 1; % possible Ctop
% 
% WCT_localmin = nan(size(Pr2__norm_rvs));
% for t = 1: size(Pr2__norm_rvs,2) %tsteps
%   clear q2
% q2 =  find(localmin_idx(:,t)==1)  ;  %index   
% WCT_localmin(q2,t) =  WCT_dnw(q2,t); 
% end
% 
% % 1st Condition: WCT/std
% % Find all positions that fullfil the wct threshold
% WCT_ctop_possible = nan(size(Pr2__norm_rvs)); % WCT @ all possible Ctop bins
% Ctop_possible = zeros(size(Pr2__norm_rvs));
% for t = 1: size(Pr2__norm_rvs,2) %tsteps      
% clear q3
% q3 =  find(isnan(WCT_localmin(:,t))== 0) ; 
% if isempty(q3)==0
% WCT_ctop_possible(q3,t) =  WCT_localmin(q3,t); end
% Ctop_possible(q3-1,t) = 1; % 1 bin before signal increase
% end
% 
% % 2nd Condition: SNR > 2 for five consecutive bins after possible Ctop
% nbins = 5;
% Ctop_possible22 = zeros(size(Pr2__norm_rvs));
% for t = 1: size(Pr2__norm_rvs,2) %tsteps
%  clear q44
% q44 = find(Ctop_possible(:,t) == 1); %index of possible Ctops
% if isempty(q44)==0
% for j = 1: length(q44)
%     clear y; clear snr
% y(j,:) = q44(j): q44(j) + nbins;
% snr(:,j)= SNRrev(y(j,:),t); %  for each q4(possible Ctop)
% 
% if any(snr(:,j) < 2)==1; Ctop_possible22(q44(j),t) = 0;
% else
%     Ctop_possible22(q44(j),t) = 1;
% end  % snr
% end % q4 size
% end % isempty(q4)
% end % tsteps
% 
% %% 3rd condition: SNR ratio threshold
% Ctop_possible2 = zeros(size(Pr2__norm_rvs));
% for t = 1: size(Pr2__norm_rvs,2) %tsteps 
%    clear q4
%  q4 = find(Ctop_possible22(:,t) == 1); %index of possible Ctops
% if isempty(q4)==0 %&& length(q4)>1
% clear snr_all
% for j = 1: length(q4)     
% snr_all(:,j)= SNRratio_rev(q4(j,:),t); %  for each q4(possible Ctop)
% if snr_all(:,j) > SNRratio_rev_thres ; Ctop_possible2(q4(j),t) = 1; % && snr_all(:,j) < 10;
% else
%     Ctop_possible2(q4(j),t) = 0;
% end  % snr
% end % q4 size
% end % isempty(q4)
% end % tsteps
% 
% %% 4th Condition: SNR ratio should increase for 3bins below Ctop 
% Ctop_possible3 = zeros(size(Pr2__norm_rvs));
% for t = 1: size(Pr2__norm_rvs,2) %tsteps   
%  clear q5
%  q5 = find(Ctop_possible2(:,t) == 1); %index of possible Ctops
% if isempty(q5)==0 %&& length(q5)>1
% clear snr_below; clear below;
% for j = 1: length(q5)    
% below(j,:) = q5(j) : q5(j)+3;
% snr_below(:,j)= diff(SNRratio_rev(below(j,:),t)); %  for each q5(possible Ctop)
% if all(snr_below(:,j)>0) ==1 ; Ctop_possible3(q5(j),t) = 1; 
% else
%     Ctop_possible3(q5(j),t) = 0;
% end  % snr
% end % q5 size
% end % isempty(q5)
% end % tsteps
% 
% %% Select only the first element(downwards) that fullfils the 4 conditions
% WCT_ctop_dyn = nan(size(Pr2__norm_rvs,2),1); % WCT @ top
% WCT_STD_ctop_dyn = nan(size(Pr2__norm_rvs,2),1); % WCT/STD @ top
% SNR_ctop_dyn = nan(size(Pr2__norm_rvs,2),1); % WCT/STD @ top
% Ctop_dyn_dnw = nan(size(Pr2__norm_rvs,2),1);  % Ctop bin
% 
% for t = 1: size(Pr2__norm_rvs,2) %tsteps      
% q6 =  find(Ctop_possible3(:,t)==1,1,'first') ;  %index 
% if isempty(q6)==1; continue; end % no cloud detected
% WCT_ctop_dyn(t,:) =  WCT_dnw(q6+1,t); % WCT@CCtop
% WCT_STD_ctop_dyn(t,:) =  WCT_STD_RATIO_for_ctop(q6+1,t);
% SNR_ctop_dyn(t,:) =  SNRrev(q6,t);
% Ctop_dyn_dnw(t,:) = q6; %Ctop  
% end
% 
% %% Cloud top flag
% ctop_flag = nan(size(Pr2__norm_rvs,2),1);
% ctop_flag(isnan(Ctop_dyn_dnw)==0) = 1;
% ctop_flag(isnan(Ctop_dyn_dnw)==1) = 0;
% 
% % // Before plotting and saving bring all the affected quantities to original height vector!
% % so that they are compatible with the original Signal, WCT,SNR, Cbase etc. //
% 
% Pr2_defliped = flip(Pr2__norm_rvs,1);
% WCT_defliped = flip(WCT_dnw); 
% WCT_STD_defliped = flip(WCT_STD_RATIO_for_ctop);
% SNRratio_defliped = flip(SNRratio_rev,1);
% Ctop_dyn = length(H)-(Ctop_dyn_dnw-1);
% WCT_ctop_defliped =  WCT_ctop_dyn;
% WCT_STD_ctop_defliped = WCT_STD_ctop_dyn;
% 
% end%function

% Now we can configure the constrained Klett retrieval (see config_constrained_Klett.m)
