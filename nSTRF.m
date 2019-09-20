%function nSTRF(spikes, cluster1, cluster2, time_window)
%
%   FILE NAME   : nSTRF.m
%   DESCRIPTION : 
%
%INPUT PARAMS
%   spikes      : full path to cluster.mat variable file containing all spikes
%   cluster1    : index of first cluster to get spike series from
%   cluster2    : index of second cluster to get spike series from
%   time_window : fixed time window to chunk spikes with
%   sprfile     : full path to spectral profile file
%
%RETURNED VARIABLES
%
%   
%
% (C) Shannon Lin, Edited Sept 2019


% TASK
% Calculate STRF for each cluster, digitize them to 0 or 1 (0 being not significant, 1 being significant from Monty?s function)
% 
% Calculate OR of 2 STRfs
% Calculate AND of 2 STRFs
% 
% Calculate nSTRF
% Find xCorr btw nSTRF/OR STRF
% Find xCorr btw nSTRF/AND STRF

% Tested by running function as follows: 
% nSTRF('/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/spike_times_ripple_clust.mat', 1, 2, 200, '/Users/shannon1/Documents/F19/neuroResearch/Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50.spr','/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/AudiResp_16_24-190326-154559_triggers.mat')

% spike_clusters='/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/spike_times_ripple_clust.mat'
% sprfile='/Users/shannon1/Documents/F19/neuroResearch/Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50_param.mat'
% Trig_190326='/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/AudiResp_16_24-190326-154559_triggers.mat'
% Trig_190324='/Users/shannon1/Documents/F19/neuroResearch/Cassius-190324/AudiResp_24_24-190324-175634_triggers.mat'
function nSTRF(spike_clusters, cluster1, cluster2, time_window, sprfile, Trig)
    % load spike trains from index cluster1 and index cluster2 from
    % spike series in specified input spikeClusters
    spikeTimeRipClusStruct = load(spike_clusters);
    spikeTimeRipClus = spikeTimeRipClusStruct.spikeTimeRipClus;
    spike1 = spikeTimeRipClus{cluster1, 2};
    spike2 = spikeTimeRipClus{cluster2, 2};
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
    % Compute STRF of each neuron
    s_bin=0.15;
    T1=0;
    T2=s_bin;
    Fss=24414.0625;
    % convert spike times from ms to sec
    spike1 = double(spike1)/1000.0;
    spet1=spike1 * Fss;
    SPL=80;
    MdB=30;
    ModType='dB';
    Sound='MR';
    NBlocks=100;
    UF=10;
    sprtype='float';
    try
        [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
    catch me
        disp(me);
        STRF1A=[];STRF2A=[];No1A=[];Wo1A=[];No1A=[];Wo1A=[];
    end
    if isempty(STRF1A)
        disp('Error when computing STRF of each neuron, exiting function')
        return;
    end
    STRF1 = STRF1A;
    STRF2 = STRF2A;
    No1   =  No1A;
    Wo1   =  Wo1A;
    No2   =  No2A;
    Wo2   =  Wo2A;
%     % Get STRF by averaging STRF params retrieved from TrigA and TrigB
%     STRF1 = (STRF1A+STRF1B)/2;
%     STRF2 = (STRF2A+STRF2B)/2;
%     No1   =  No1A+No1B;
%     Wo1   =  (Wo1A+Wo1B)/2;
%     No2   =  No2A+No2B;
%     Wo2   =  (Wo2A+Wo2B)/2;
%     
%     
%     % Compute significant STRF of each neuron (mark 1 as sig, 0 as not)
%     p=0.001;
%     SModeType='dB';
%     [STRF1s,Tresh1]=wstrfstat(STRF1,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
%     [STRF1s,Tresh2]=wstrfstat(STRF2,p,No2,Wo2,PP,MdB,ModType,Sound,SModType);
% 
%     % Compute nSTRF using different bin sizes
%     %   What does this mean?
%     %   From findOverlap function, I have overlap spike times
%     overlapSpike = findOverlap(spikes, cluster1, cluster2, time_window);
% 
%     % Compute from OR/AND STRF from STRF1, STRF2
%     
%     % Compute xCorr btw OR/nSTRF and AND/nSTRF
%     
%     % Compute the cross correlation between the 2 spike trains and calc
%     % significance of corr coeff (calculate the CC between 0-10 ms)
%  
end