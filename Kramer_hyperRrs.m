%Sasha Kramer
%sjkramer@bu.edu
%UCSB IGPMS - MBARI - Boston University

%%%Script to run hyperspectral reflectance model from Kramer et al. (2022)
%Runs with: betasw_ZHH2009.m, gsm_invert.m, gsm_cost.m
%Inputs: measured Rrs, temp, sal, asw, A&B coeffs for aph

%Rrs is a function of IOPs, specifically absorption and backscattering. 
%Here, you will build the absorption and backscattering terms from their 
%component parts, and then match the modeled spectra to the measured 
%spectra to optimize the model fit.

%Make sure all functions are in your path / directory!

%Map the directory where you will load your data:
cd /Users/skramer/Documents/UCSB/Research/Data/HPLC_Aph_Rrs

%Load your data
load Kramer_rrs_testdata.mat
%includes Rrs, chl, T, S, latlon 

%The model uses reflectance = f(IOPs) from Gordon et al. (1988) which uses
%below the surface reflectance (rrs = Lu(0-)/Ed(0-)). If you are using 
%above surface reflectance (Rrs = Lu(0+)/Ed(0+)), you will need to convert
%your reflectances before running this model, using the equation below from
%Lee et al. (2002):
rrs = Rrs./(0.52 + 1.7*Rrs);
rrs = rrs'; %putting matrix in right configuration for gsm_invert.m
Rrs = Rrs';

%Define wavelengths
wave = 400:1:700;

%%First, define total absorption as a sum of seawater absorption (asw), phytoplankton absorption (aph) and CDOM plus other detrital matter (acdm)
%a_tot = asw + aphstar + acdm
load asw_all.mat
asw = asw_all(51:end,2); %subset for your wavelengths
clear asw_all

%aph = import A & B coefficients from aph_A_B_Coeffs_Sasha_RSE_paper.txt
%aph = A.*chl.^B; %inversion model will solve for chl as an output
load AB_coeffs.mat 
A = A(51:end); %subset for your wavelengths
B = B(51:end); %subset for your wavelengths

%acdm slope is a function of Rrs (just above surface):
%acdm_s = -(0.01447 + 0.00033*Rrs490/Rrs555); %You will need to define Rrs490 and Rrs555 based on your Rrs data, using test values below
acdm_s = -(0.01447 + 0.00033*Rrs(91,:)./Rrs(156,:));
acdm = exp(acdm_s'.*(wave-443))';
clear acdm_s

%%Then, define backscattering as a sum of seawater backscattering (bbsw) and backscattering by particles (bbp)
%bb_tot = bbsw + bbp

%bsw comes from Zhang et al. (2009): run for each spectrum (example shown with test data matrix)
%[~,~,bsw] = betasw_ZHH2009(wave,T,[],S);
for i = 1:17
    [~,~,bsw(i,:)] = betasw_ZHH2009(wave,T(i),[],S(i));
end
bbsw = 0.5*bsw';
clear bsw i S T

%bbp slope is a function of rrs (just below surface):
%bbp_s = 2.0*(1.-1.2*exp(-0.9*rrs440./rrs555)); %You will need to define rrs440 and rrs555 based on your rrs data, using test values below
bbp_s = 2.0*(1.-1.2*exp(-0.9*rrs(41,:)./rrs(156,:)));
bbp = (443./wave').^bbp_s;
clear bbp_s

%Put IOPs together: run for each spectrum (example shown with test data matrix)
%[IOPs,output] = gsm_invert(rrs,asw,bbsw,bbp,A,B,acdm);
for i = 1:17
    [IOPs(i,:),~] = gsm_invert(rrs(:,i)',asw',bbsw(:,i)',bbp(:,i)',A',B',acdm(:,i)');
end
clear i
%outputs = chl, acdm443, bbp443

%Reconstruct Rrs: run for each spectrum (example shown with test data matrix)
%a = asw + A.*IOPs(:,1).^B + IOPs(:,2)*acdm;
%bb = bbsw + IOPs(:,3)*bbp;
%rrsP = bb./(a+bb);
for i = 1:17
    a(:,i) = asw + (A.*IOPs(i,1).^B) + (IOPs(i,2)*acdm(:,i));
    bb(:,i) = bbsw(:,i) + IOPs(i,3)*bbp(:,i);
end
clear i
rrsP = bb./(a+bb);

g = [0.0949 0.0794]; %coefficients from Gordon et al. (1988)
modrrs = (g(1) + g(2)*rrsP).*rrsP;

%Convert rrs back to Rrs
modRrs = (0.52*modrrs)./(1-(1.7*modrrs));

%Residual between measured and modeled (to use for Kramer_hyperRrs.m)
rrsD = rrs-modrrs;
RrsD = Rrs-modRrs;
