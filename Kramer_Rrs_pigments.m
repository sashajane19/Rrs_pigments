%Sasha Kramer
%sjkramer@bu.edu
%UCSB IGPMS - MBARI - Boston University

%%%Script to model pigments from hyperspectral reflectance residual
%Modified from aph version written by Dylan Catlett: https://github.com/dcat4/bioOptix_and_PFTs

%Map the directory where you will load your data:
cd /Users/skramer/Documents/UCSB/Research/Data/HPLC_Aph_Rrs

%Load your samples (formatted here as a .mat file):
load HGSM145_20200918.mat

hgsm_wl = HGSM145.wl; %wavelengths
hgsm_diff = HGSM145.Rrsdiff'; %Rrs residual

%Set location to save modeled pigment results
cd /Users/skramer/Documents/UCSB/Research/Data/HPLC_Aph_Rrs/modelrun/20210810

%Possible model inputs for different interations
%Paper uses just 2nd deriv of Rrs residual - other results in supplement
diffD = diff(hgsm_diff,1,2); %1st deriv residual
diffD2 = diff(hgsm_diff,2,2); %2nd deriv residual

%Set parameters that are constant for all modeling runs:
n_permutations = 100; %the number of independent model validations to do (each formulates and validates a model)
max_pcs = 30; %max number of spectral pc's to incorporate into the model - 30
mdl_pick_metric = 'MAE'; %pick the metric by which to evaluate pigment model fit - R2, RMSE, avg, med, ens, bias, MAE
k = 5;
pft_index = 'pigment';
ofn_suffix = '_rrsD2_1nm_MAE_20210810.mat'; % the suffix of the name of the file (change for your output)

%Loop through and train models for each of the pigments listed below.
%Make sure they each have a match to a name in the "vars" variable
pigs2mdl = {'Tchla','Zea','DVchla','ButFuco','HexFuco','Allo','MVchlb','Neo','Viola','Fuco','Chlc12','Chlc3','Perid'};
vars = {'Tchla','Tchlb','Tchlc','ABcaro','ButFuco','HexFuco','Allo','Diadino','Diato','Fuco','Perid','Zea','MVchla','DVchla','Chllide','MVchlb','DVchlb','Chlc12','Chlc3','Lut','Neo','Viola','Phytin','Phide','Pras'};

%Set spectra to use in the model:
daph = diffD2; %"daph" is the variable name from Dylan's code for consistency, but here we are using d2rrs

%% Start modeling...
%this script keeps track of how long it takes ("tic" and "toc")
%when the script is done running, handel's messiah plays - to turn off, comment out lines 55-56
tic
for i = 1:length(pigs2mdl)
    this1 = pigs2mdl{i};
    jj = strcmp(vars,this1);
    coli = find(jj == 1);
    pft = Global_RHPLC(:,coli); 
    output_file_name = [this1,ofn_suffix];
    [coefficients, intercepts, summary_gofs, all_gofs] = rrsModelTrain(daph, pft, pft_index, n_permutations, max_pcs, k, mdl_pick_metric, output_file_name);
end
toc

load handel
sound(y,Fs)
