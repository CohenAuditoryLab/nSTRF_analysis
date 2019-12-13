% Defining STRF params
Params.T1=0;
Params.T2=0.15;
Params.Fss=24414.0625;
Params.SPL=80;
Params.MdB=30;
Params.ModType='dB';
Params.Sound='MR';
Params.NBlocks=400;
Params.UF=10;
Params.sprtype='float';
Params.p=0.001;
Params.SModType='dB';

% full path to .mat containing spike times
spikeClusters = '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat';

% full path to .spr sprfile
sprfile = '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/DNR_Cortex_96k5min_4_50.spr';

% full path to .mat Trig
Trig = '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/AudiResp_16_24-190326-154559_triggers.mat';

% Cassius_190324 saved as "spikeTimeRip"
% Cassius_190326 data saved as "st_clu" 
version = 'st_clu';

% 'y' to analyze all clusters
% 'n' to analyze specified numClusters amt of clusters
analyzeAll = 'n';
numClusters = 2;

% collect indices of individual clusters
getIndividual(Params, spikeClusters, sprfile, Trig, version, analyzeAll, numClusters);

% collect indices of 2STRF clusters
getPaired(Params, spikeClusters, sprfile, Trig, version, analyzeAll, numClusters);