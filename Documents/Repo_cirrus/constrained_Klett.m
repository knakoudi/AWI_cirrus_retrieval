% AWI_cirrus_retrieval: constrained_Klett
% *AWI (Alfred Wegener Institute)

% This routine contains the constrained Klett algorithm for cirrus cloud optical retrievals
% The algorithm is described in Nakoudi et al. (2020): An extended lidar-based cirrus cloud retrieval scheme: first
% application over an Arctic site, submitted to Optics Express // Flowchart Fig. 4 and section 2.3.2

% Last update: 09-11-2020
% Authors:  Christoph Ritter (christoph.ritter@awi.de)
%                Konstantina Nakoudi (konstantina.nakoudi@awi.de)

%% Prerequisites: 
%  Defined cirrus geometrical boundaries (see cirrus_detection.m)
% configuration of constrained Klett parameters (see config_constrained_Klett.m)
%  A routine for the classical Klett-Fernald retrievals

%function: constrained_Klett_func
%% Input parameters:
% PXXX: Lidar signal at XXX channel
% SNRXX: Signal-to-Noise ratio of XXX channel
% listvalidXXX: valid cirrus profiles
% BSRXXXref: BSR reference value at XXX wavelelngth
% convrXXX: convergence range bins for XXX wanelength
% LRci_initial: initial guess for cirrus LR
% LR_part: initial guess for particulate LR
% AlRay: molecular extinction coefficient at XXX wavelength
% BeRa: molecular backscatter coefficient at XXX wavelength
% Density: air density 
% Cbase_dyn: cloud base from dynamic WCT
% Ctop_dyn: cloud top from dynamic WCT
% BSRcal_XXX: BSR calibration value at far-range
% BSRcal_XXX_err: BSR calibration value error at far-range
% CalRange: Calibation range for backward Klett-Fernald retrieval
% H: Height vector
% dH: vertical resolution

%% Output parameters:
%  BSRXXX_Klettconstr: BSR as derived from constrained  Klett at XXX wavelength 
%  BSRXXXerr_Klettconstr: BSR error as derived from constrained  Klett at XXX wavelength 
% BetaXXX_Klettconstr: particulate backscatter coefficient as derived from constrained  Klett at XXX wavelength 
% BetaXXXerr_Klettconstr: particulate backscatter coefficient error as derived from constrained  Klett at XXX wavelength 
% dBetaXXXdP: particulate backscatter coefficient error due to signal noise
% dBetaXXXdR: particulate backscatter coefficient error due to BSR uncertainty
% dBetaXXXdLR: particulate backscatter coefficient error due to LR uncertainty
% LRXXX_iter: LR array of indicidual iterations
% iter: total number of iterations
% LRXXX_final: final LR array
% LRXXX_ci: final cirrus LR
% CODXXXci: cirrus COD
% goodprofXXX: profiles with reasonable cirrus LR as determined by  Ansmann
% et al., 1992; Chen et al., 2002; Borovoi et al., 2014; Okamoto et al., 2019; Okamoto et al., 2020

% function [BSRXXX_Klettconstr, BSRXXXerr_Klettconstr, BetaXXX_Klettconstr, BetaXXXerr_Klettconstr, dBetaXXXdP, dBetaXXXdR,...
%     dBetaXXXdLR, LRXXX_iter, iter, LRXXX_final, LRXXX_ci, CODXXXci, goodprofXXX] ...
%          = constrained_Klett_func(PXXX, SNRXXX, listvalidXXX, BSRXXXref,convrXXX, LRci_initial, LR_part, AlRay, BeRa, Density,...
%                                           Cbase_dyn, Ctop_dyn, BSRcal_XXX, BSRcal_XXX_err, CalRange, H, dH)
 
% Physically reasonable LR limits as derived from:
% A. Ansmann, U. Wandinger, M. Riebesell, C. Weitkamp, and W. Michaelis, ?Independent measurement of extinction
%and backscatter profiles in cirrus clouds by using a combined raman elastic-backscatter lidar,? Appl. Opt. 31, 7113?7131 (1992).
%
% W.-N. Chen, C.-W. Chiang, and J.-B. Nee, ?Lidar ratio and depolarization ratio for cirrus clouds,? Appl. Opt. 41,
% 6470?6476 (2002).
%
% A. Borovoi, A. Konoshonkin, and N. Kustova, ?Backscatter ratios for arbitrary oriented hexagonal ice crystals of
% cirrus clouds,? Opt. letters 39, 5788?5791 (2014).
%
% H. Okamoto, K. Sato, A. Borovoi, H. Ishimoto, K. Masuda, A. Konoshonkin, and N. Kustova, ?Interpretation of
% lidar ratio and depolarization ratio of ice clouds using spaceborne high-spectral-resolution polarization lidar,? Opt.
% Express 27, 36587?36600 (2019).
%
% H. Okamoto, K. Sato, A. Borovoi, H. Ishimoto, K. Masuda, A. Konoshonkin, and N. Kustova, ?Wavelength
% dependence of ice cloud backscatter properties for space-borne polarization lidar applications,? Opt. Express 28,
% 29178?29191 (2020).

LRlowerlim = 5;  % maximmum cirrus LR limit            
LRupperlim = 90;  % maximmum cirrus LR limit

Hcalcrange = 32000; % maximum range for retrieval 
LRarrerr = ones(size(PXXX))+15;  % LR uncertainty 
Hmin1 = 81; % full-overlap range

% pre-allocation
BetaXXX_Klettconstr=zeros(size(PXXX)); BetaXXXerr_Klettconstr=zeros(size(PXXX)); 
BSRXXX_Klettconstr =zeros(size(PXXX)); BSRXXXerr_Klettconstr=zeros(size(PXXX)); %predefine new solution-->iter
dBetaXXXdLR=zeros(size(PXXX));dBetaXXXdR=zeros(size(PXXX)); 
dBetaXXXdP=zeros(size(PXXX));

LRXXX_final = LR_part.*ones(size(PXXX));  

%% Initiation of iteration
conv_perc = 0.3; % convergence percentange 
itmax= 30; % maximum number of iterations // 2-4 iterations are usually enough
LRXXX_iter=zeros(itmax, size(PXXX,2));  

% flag_iter = 1; --> too high LR more than 3 consecutive times
% flag_iter = -1; --> too low LR more than 3 consecutive times
flag_iter = nan(size(PXXX,2),1); 
flag_ok = nan(size(PXXX,2),1);

% Use LRci only in WCT-derived cirrus bins 
WCTcirrus = [];    WCTcirrusarr= zeros(size(PXXX));  iter=zeros(size(PXXX,2),1);
Pschwelle=1e-8;   
Perr= PXXX./SNRXXX; % signal noise
for i= listvalidXXX'  
if ~isnan(Cbase_dyn(i)) &&  ~isnan(Ctop_dyn(i))  
WCTcirrus = Cbase_dyn(i): Ctop_dyn(i); 
WCTcirrusarr(WCTcirrus,i) = WCTcirrus; 
LRXXX_final(WCTcirrus,i) = LRci_initial; 
end
if ~isempty(WCTcirrus)

Sel=find(H >= 0 & H <= Hcalcrange & PXXX(:,i)>Pschwelle);%Sel==Selected points for calculation  %PC(:,i)>Pschwelle excludes NaN 
problem2=any(Density(Sel,i)) ; if problem2==0, Sel=[]; end

mincalclength=2000; mincalcsteps=round(mincalclength ./ dH);
if length(Sel) >mincalcsteps,      
iter=0; condi=1;  hz=0; tz=0;  % hz / tz Zaehler f?r hohe - tiefe Werte
refer=zeros(itmax,1);  refer2=zeros(itmax,1);

%  While condi is True --> iterate through all timesteps
while condi
iter = iter+1;    
% initial guess Klett retrieval with LR1 
% //replace with your routine for classical Klett-Fernald retrievals
% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRcal_XXX(i), CalRange, H(Sel), PXXX(Sel,i), Perr(Sel,i), ...
%     LRXXX_final(Sel,i), AlRay(Sel,i), BeRa(Sel,i));//
Betaaer(Sel)=Beta-BeRa(Sel,i);  %Baer=Btot-Bray
Betaaererr(Sel,i)=abs(dBdR.*BSRcal_XXX_err)+abs(dBdLR.*LRarrerr(Sel,i))+abs(dBdP.*Perr(Sel,i)); %Beta error 
Btemp(Sel)=Beta./BeRa(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa(Sel,i); %Temporary BSR and BSR error

% initial guess Klett retrieval with LR2 = LR1 + 1 sr
% // replace with your routine for classical Klett-Fernald retrievals 
% [Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4( BSRcal_XXX(i), CalRange, H(Sel), PXXX(Sel,i), Perr(Sel,i), ...
%     LRXXX_final(Sel,i)+1, AlRay(Sel,i), BeRa(Sel,i));  //
%Btemp2(Sel)=Beta2./BeRa(Sel,i); 
% Btemperr(Sel)=Betaaererr(Sel,i)./BeRa(Sel,i); %Temporary BSR and BSR error

BSR1(iter) = median(Btemp(convrXXX),'omitnan');  % median BSR1 within convergence range (using the BETA1)
BSR2(iter) = median(Btemp2(convrXXX),'omitnan'); % median BSR2 within convergence range (using the BETA2)
dBSRdLR = BSR2(iter)-BSR1(iter); % difference of the two BSR median values --> dependence on LR
dBSR = BSRXXXref - BSR1(iter); % difference between the BSRref  and BSR1 
dLR=dBSR./dBSRdLR; %  BSRref-BSR1/BSR2-BSR1 // see Eq. 3 from Nakoudi et al. (2020)

% Is the convergence to BSRref is satisfactory?
idcir = find(WCTcirrusarr(:,i)~=0);
if BSR1(iter) ./ BSRXXXref > 1+conv_perc./100,     %  BSR1/BSRref > 100.3% 
    LRXXX_final(idcir,i) = LRXXX_final(idcir,i) +dLR; % tune LR by dLR>0
    %if median LR exceeds upper limit, then make LR equal to upper limit
    if median(LRXXX_final(idcir,i)) > LRupperlim 
        LRXXX_final(idcir,i) = LRupperlim;              
        hz=hz+1;
    end
elseif BSR1(iter) ./ BSRXXXref < 1-conv_perc./100,      %  BSR1/BSRref < 99.7% 
    LRXXX_final(idcir,i) = LRXXX_final(idcir,i) +dLR; % tune LR by dLR<0
    %if median LR exceeds lower limit, then make LR equal to lower limit
    if median(LRXXX_final(idcir,i)) < LRlowerlim  
        LRXXX_final(idcir,i) = LRlowerlim;       
        tz=tz+1;
    end
else
    condi=0; disp('Finally, a reasonable BSR was found (very close to BSRXXXref)'); flag_ok(i)=1; 
end

% flag problematic profiles
if iter >=itmax, condi=0;  disp('no convergence'); flag_noconv(i)=1;  end  %  max iterations were reached without satisfactory convergence
if hz>3, condi = 0;  disp('high LR for 3 times'); flag_iter(i)=1; end    % too high LR was found more than 3 times
if tz>3, condi = 0; disp('low LR for 3 times'); flag_iter(i)=-1; end      % too low LR was found more than 3 times
LRXXX_iter(iter,i) = LRXXX_final(idcir(1),i);   % for diagnostic reasons --> how LR changes in each iteration
end % while
   
iter(i) = iter; % number of iterations  
if LRXXX_final(idcir,i) > LRupperlim; disp('very high LR'); flag_iter(i)=1; end  
if LRXXX_final(idcir,i) < LRlowerlim; disp('very low LR'); flag_iter(i)=-1; end

%  constrained Klett solution and uncertainty characterization
BSRXXX_Klettconstr(Sel,i)=Btemp(Sel); BSRXXXerr_Klettconstr(Sel,i)=Btemperr(Sel); 
BetaXXX_Klettconstr(Sel,i)=Betaaer(Sel); BetaXXXerr_Klettconstr(Sel,i)=Betaaererr(Sel,i); 
dBetaXXXdP(Sel,i)=dBdP; 
dBetaXXXdR(Sel,i)=dBdR; 
dBetaXXXdLR(Sel,i)=dBdLR;

if mean(BSRXXX_Klettconstr(Hmin1:Cbase_dyn(i),i))<1; disp('WARNING:BSR<1'); flag_negb(i)=1; end
end 
end
end 

% cirrus LR
LRXXX_ci = nan(size(PXXX,2),1); 
for i = listvalidXXX'
if ~isnan(Cbase_dyn(i)) &&  ~isnan(Ctop_dyn(i)) 
    LRXXX_ci(i) = mean(LRXXX_final(Cbase_dyn(i):Ctop_dyn(i),i)); 
end
end

% Estimate particulate backscatter coefficient at channel with
% perdendicular polarization (only for Lidar system with polarization detection capabilitites)
BetaAerXXXtot = BetaXXX_Klettconstr + BetaXXX_S; 

 % cirrus COD // see Eq. 4 from Nakoudi et al. (2020)
CODXXXci = nan(n_final,1);
for i= listvalidXXX'
CODXXXci(i) = trapz(H(Cbase_dyn(i):Ctop_dyn(i)), LRXXXci(i)'.*BetaAerXXXtot(Cbase_dyn(i):Ctop_dyn(i),i)); 
end

% Flag profiles with reasonable LR
 goodprofXX = nan(size(PXXX,2),1); 
for i = listvalidXXX'; if flag_ok(i)==1;  goodprofXX(i) = 1;  end; end

