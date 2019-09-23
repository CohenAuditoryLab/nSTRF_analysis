%function [nSTRF,orSTRF,andSTRF,nCorrOr,nCorrAnd] = pairedSTRFanalysis(spike_clusters,cluster1,cluster2,time_window,sprfile,Trig)
%
%   FILE NAME   : pairedSTRFanalysis.m
%   DESCRIPTION : This file computes the significant STRF of spike times
%   from 2 clusters, calculates the OR and AND STRF btw the sig STRF of
%   the 2 clusters, computes the nSTRF btw the 2 clusters, and then finds
%   the cross correlation btw the OR/nSTRF and AND/nSTRF
%   
%
%INPUT PARAMS
%   spike_clusters  : full path to cluster.mat variable file containing all spikes
%   cluster1        : index of first cluster to get spike series from
%   cluster2        : index of second cluster to get spike series from
%   time_window     : fixed time window to chunk spikes with
%   sprfile         : full path to spectral profile file
%   Trig            : full path to trigger data
%
%RETURNED VARIABLES
%   nSTRF           : nSTRF between 2 specified spike clusters
%   orSTRF          : orSTRF btw the 2 sig STRF of specified clusters
%   andSTRF         : andSTRF btw the 2 sig STRF of specified clusters  
%   nCorrOr         : cross correlation between nSTRF/orSTRF
%   nCorrAnd        : cross correlation between nSTRF/andSTRF
%
% (C) Shannon Lin, Edited Sept 2019

% Tested by running function as follows: 
% pairedSTRFanalysis('/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/spike_times_ripple_clust.mat', 1, 2, 200, '/Users/shannon1/Documents/F19/neuroResearch/Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50.spr','/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/AudiResp_16_24-190326-154559_triggers.mat')

% spike_clusters='/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/spike_times_ripple_clust.mat'
% sprfile='/Users/shannon1/Documents/F19/neuroResearch/Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50_param.mat'
% Trig_190326='/Users/shannon1/Documents/F19/neuroResearch/Cassius-190326/AudiResp_16_24-190326-154559_triggers.mat'
% Trig_190324='/Users/shannon1/Documents/F19/neuroResearch/Cassius-190324/AudiResp_24_24-190324-175634_triggers.mat'

function [nSTRF,orSTRF,andSTRF,nCorrOr,nCorrAnd] = pairedSTRFanalysis(spike_clusters,cluster1,cluster2,time_window,sprfile,Trig)

    % Load spike times at index cluster1 and cluster2 in spike_clusters
    spikeTimeRipClusStruct = load(spike_clusters);
    spikeTimeRipClus = spikeTimeRipClusStruct.spikeTimeRipClus;
    spike1 = spikeTimeRipClus{cluster1, 2};
    spike2 = spikeTimeRipClus{cluster2, 2};
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
    
    % Compute STRF of each cluster
    s_bin=0.15;
    T1=0;
    T2=s_bin;
    Fss=24414.0625;
    spike1 = double(spike1)/1000.0; % convert spike times from ms to sec
    spike2 = double(spike2)/1000.0; % convert spike times from ms to sec
    spet1=spike1 * Fss;
    spet2=spike2 * Fss;
    SPL=80;
    MdB=30;
    ModType='dB';
    Sound='MR';
    NBlocks=100;
    UF=10;
    sprtype='float';
    try
        [taxis,faxis,clusOneSTRF1A,clusOneSTRF2A,clusOnePP,clusOneWo1A,clusOneWo2A,clusOneNo1A,clusOneNo2A,clusOneSPLN]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        [taxis,faxis,clusOneSTRF1B,clusOneSTRF2B,clusOnePP,clusOneWo1B,clusOneWo2B,clusOneNo1B,clusOneNo2B,clusOneSPLN]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
        [taxis,faxis,clusTwoSTRF1A,clusTwoSTRF2A,clusTwoPP,clusTwoWo1A,clusTwoWo2A,clusTwoNo1A,clusTwoNo2A,clusTwoSPLN]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        [taxis,faxis,clusTwoSTRF1B,clusTwoSTRF2B,clusTwoPP,clusTwoWo1B,clusTwoWo2B,clusTwoNo1B,clusTwoNo2B,clusTwoSPLN]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
    catch me
        disp('Error when computing STRF1, exiting function')
        disp(me);
        return;
    end
    % Average STRF1 from TrigA and TrigB for both clusters
    %   Ignoring STRF2 for now
    clusOneSTRF = (clusOneSTRF1A+clusOneSTRF1B)/2;
    clusOneNo1 = clusOneNo1A + clusOneNo1B;
    clusOneWo1 = (clusOneWo1A + clusOneWo1B)/2;
    clusTwoSTRF = (clusTwoSTRF1A+clusTwoSTRF1B)/2;
    clusTwoNo1 = clusTwoNo1A + clusTwoNo1B;
    clusTwoWo1 = (clusTwoWo1A + clusTwoWo1B)/2;
   
    % Compute significant STRF of each cluster (mark 1 as sig, 0 as not)
    p=0.001;
    SModType='dB';
    [clusOneSTRF1s,clusOneTresh1]=wstrfstat(clusOneSTRF,p,clusOneNo1,clusOneWo1,clusOnePP,MdB,ModType,Sound,SModType);
    [clusTwoSTRF1s,clusTwoTresh1]=wstrfstat(clusTwoSTRF,p,clusTwoNo1,clusTwoWo1,clusTwoPP,MdB,ModType,Sound,SModType);

    % Compute nSTRF using different bin sizes
    %   Currently only calculating nSTRF for specified time_window
    [nSTRF]=findOverlap(spike_clusters, cluster1, cluster2, time_window);
    assignin('base', 'nSTRF', nSTRF);

    % Compute from OR/AND STRF from STRF1, STRF2
    orSTRF = clusOneSTRF1s | clusTwoSTRF1s;
    andSTRF = clusOneSTRF1s & clusTwoSTRF1s;
    orSTRF = double(orSTRF);
    andSTRF = double(andSTRF);
    assignin('base', 'orSTRF', orSTRF);
    assignin('base', 'andSTRF', andSTRF);
%     O = any(orSTRF(:) > 0);
%     A = any(andSTRF(:) > 0);
%     assignin('base', 'O', O);
%     assignin('base', 'A', A);
    
    % Compute xCorr btw OR/nSTRF and AND/nSTRF
    nCorrOr = xcorr2(nSTRF, orSTRF);
    nCorrAnd = xcorr2(nSTRF, andSTRF);
    assignin('base', 'nCorrOR', nCorrOr);
    assignin('base', 'nCorrAND', nCorrAnd); 
%     oCorr = any(nCorrOr(:) > 0);
%     aCorr = any(nCorrAnd(:) > 0);
%     assignin('base', 'oCorr', oCorr);
%     assignin('base', 'aCorr', aCorr);
end