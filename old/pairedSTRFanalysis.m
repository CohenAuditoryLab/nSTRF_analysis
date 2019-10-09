% function [clusOneSTRF,clusTwoSTRF,nSTRF,orSTRF,andSTRF,nCorrOr,nCorrAnd] = pairedSTRFanalysis(spike_clusters,cluster1,cluster2,time_window,sprfile,Trig)
%
%   FILE NAME   : pairedSTRFanalysis.m
%   DESCRIPTION : This file computes the significant STRF of spike times
%   from 2 clusters, calculates the OR and AND STRF btw the sig STRF of
%   the 2 clusters, computes the nSTRF btw the 2 clusters, and then finds
%   the cross correlation btw the OR/nSTRF and AND/nSTRF
%   
%
% INPUT PARAMS
%   spike_clusters  : full path to cluster.mat variable file containing all spikes
%   cluster1        : first cluster to get spike series from
%   cluster2        : second cluster to get spike series from
%   time_window     : fixed time window to chunk spikes with
%   sprfile         : full path to spectral profile file
%   Trig            : full path to trigger data
%   version         : 'st_clu' or 'spikeTimeRip', refers to what struct exists in
%                       spike_clusters
%
% RETURNED VARIABLES
%   clusOneSTRF     : STRF of cluster 1 spike time
%   clusTwoSTRF     : STRF of cluster 2 spike time
%   nSTRF           : nSTRF between 2 specified spike clusters
%   orSTRF          : orSTRF btw the 2 sig STRF of specified clusters
%   andSTRF         : andSTRF btw the 2 sig STRF of specified clusters  
%   nCorrOr         : cross correlation between nSTRF/orSTRF
%   nCorrAnd        : cross correlation between nSTRF/andSTRF
%
% (C) Shannon Lin, Edited Sept 2019

% Tested by running function as follows: 
% pairedSTRFanalysis('/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat', 6, 9, 500, '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/DNR_Cortex_96k5min_4_50.spr','/Users/shannon1/Documents/F19/neuroResearch/nSTRF/AudiResp_16_24-190326-154559_triggers.mat', 'st_clu')

function [clusOneSTRF,clusTwoSTRF,nSTRF,orSTRF,andSTRF,nCorrOr,nCorrAnd] = pairedSTRFanalysis(spike_clusters,cluster1,cluster2,time_window,sprfile,Trig, version)
    % Load spike times at index cluster1 and cluster2 in spike_clusters
    % along with overlap spike
    [spikeTimeRipClus, spike1, spike2, overlap] = findOverlap(spike_clusters, cluster1, cluster2, time_window, version);
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;

    % Compute STRF of each cluster
    s_bin=0.15;
    T1=0;
    T2=s_bin;
    Fss=24414.0625;
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
        [oneTaxis,oneFaxis,clusOneSTRF1A,clusOneSTRF2A,clusOnePP,clusOneWo1A,clusOneWo2A,clusOneNo1A,clusOneNo2A,clusOneSPLN]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        [oneTaxis,oneFaxis,clusOneSTRF1B,clusOneSTRF2B,clusOnePP,clusOneWo1B,clusOneWo2B,clusOneNo1B,clusOneNo2B,clusOneSPLN]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
        [twoTaxis,twoFaxis,clusTwoSTRF1A,clusTwoSTRF2A,clusTwoPP,clusTwoWo1A,clusTwoWo2A,clusTwoNo1A,clusTwoNo2A,clusTwoSPLN]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        [twoTaxis,twoFaxis,clusTwoSTRF1B,clusTwoSTRF2B,clusTwoPP,clusTwoWo1B,clusTwoWo2B,clusTwoNo1B,clusTwoNo2B,clusTwoSPLN]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
    catch me
        disp('Error when computing STRF1/STRF2 of individual clusters, exiting function')
        disp(me);
        return;
    end
    
    % Average STRF1 from TrigA and TrigB for both clusters
    clusOneSTRF = (clusOneSTRF1A+clusOneSTRF1B)/2;
    assignin('base', 'clusOneSTRF', clusOneSTRF);
    clusOneNo1 = clusOneNo1A + clusOneNo1B;
    clusOneWo1 = (clusOneWo1A + clusOneWo1B)/2;
    clusTwoSTRF = (clusTwoSTRF1A+clusTwoSTRF1B)/2;
    assignin('base', 'clusTwoSTRF', clusTwoSTRF);
    clusTwoNo1 = clusTwoNo1A + clusTwoNo1B;
    clusTwoWo1 = (clusTwoWo1A + clusTwoWo1B)/2;
     
    % Compute significant STRF of each cluster (mark 1 as sig, 0 as not)
    p=0.001;
    SModType='dB';
    [clusOneSTRF1s,clusOneTresh1]=wstrfstat(clusOneSTRF,p,clusOneNo1,clusOneWo1,clusOnePP,MdB,ModType,Sound,SModType);
    assignin('base', 'clusOneSTRF1s', clusOneSTRF1s);
    [clusTwoSTRF1s,clusTwoTresh1]=wstrfstat(clusTwoSTRF,p,clusOneNo1,clusOneWo1,clusOnePP,MdB,ModType,Sound,SModType);
    assignin('base', 'clusTwoSTRF1s', clusTwoSTRF1s);
    
    % Compute nSTRF of overlap spike time
    if (isempty(overlap))
        disp('No overlap spike time between spike1 and spike2 so nSTRF cannot be calculated, exiting function')
        return;
    end
    spetOverlap=overlap * Fss;
    try
        [nTaxis,nFaxis,nSTRF1A,nSTRF2A,nPP,nWo1A,nWo2A,nNo1A,nNo2A,nSPLN]=rtwstrfdbint(sprfile,T1,T2,spetOverlap',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        assignin('base', 'nSTRF1A', nSTRF1A);
        [nTaxis,nFaxis,nSTRF1B,nSTRF2B,nPP,nWo1B,nWo2B,nNo1B,nNo2B,nSPLN]=rtwstrfdbint(sprfile,T1,T2,spetOverlap',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
        assignin('base', 'nSTRF1B', nSTRF1B);
    catch me
        disp('Error when computing nSTRF1A/B, exiting function')
        disp(me);
        return;
    end
    
    % Average STRF1 from TrigA and TrigB for overlap spike times
    nSTRF = (nSTRF1A+nSTRF1B)/2;
    assignin('base', 'nSTRF', nSTRF);
    nNo1 = nNo1A + nNo1B;
    nWo1 = (nWo1A + nWo1B)/2;
    assignin('base', 'nSTRF', nSTRF);
    
    % Compute from OR/AND STRF from STRF1, STRF2
    orSTRF = clusOneSTRF1s | clusTwoSTRF1s;
    andSTRF = clusOneSTRF1s & clusTwoSTRF1s;
    orSTRF = double(orSTRF);
    andSTRF = double(andSTRF);
    assignin('base', 'orSTRF', orSTRF);
    assignin('base', 'andSTRF', andSTRF);
     
    % Compute significant nSTRF
    [nSTRF1s,nTresh1]=wstrfstat(nSTRF,p,nNo1,nWo1,nPP,MdB,ModType,Sound,SModType);
    assignin('base', 'nSTRF1s', nSTRF1s);
    
    % === TO DO ==== 
    % cross correlation btw n/OR, n/AND (check out xcorr, xcorr2)
    % Compute correlation btw OR/nSTRF and AND/nSTRF
%     nCorrOr = nSTRF1s | orSTRF;
%     nCorrAnd = nSTRF1s & orSTRF;
%     nCorrOr = double(nCorrOr);
%     nCorrAnd = double(nCorrAnd);
    assignin('base', 'nCorrOr', nCorrOr);
    assignin('base', 'nCorrAnd', nCorrAnd);
    
    % plot STRF1, STRF2, nSTRF, orSTRF, andSTRF
    subplot(2,5,1)
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),clusOneSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['STRF no ' int2str(cluster1)]);
    
    subplot(2,5,2)
    pcolor(twoTaxis,log2(twoFaxis/twoFaxis(1)),clusTwoSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['STRF no ' int2str(cluster2)]);
    
    subplot(2,5,3)
    pcolor(nTaxis,log2(nFaxis/nFaxis(1)),nSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['nSTRF ' int2str(cluster1) ':' int2str(cluster2) '-' num2str(double(time_window))]);
    
    subplot(2,5,4)
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),orSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title('orSTRF');
    
    subplot(2,5,5)
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),andSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title('andSTRF');
    
    subplot(2,5,6)
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),nCorrOr);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title('nCorrOr');
    
    subplot(2,5,7)
    pcolor(twoTaxis,log2(twoFaxis/twoFaxis(1)),nCorrAnd);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title('nCorrAnd');
    
    subplot(2,5,8)
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),clusOneSTRF1s);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['STRF1s no ' int2str(cluster1)]);
    
    subplot(2,5,9)
    pcolor(twoTaxis,log2(twoFaxis/twoFaxis(1)),clusTwoSTRF1s);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['STRF2s no ' int2str(cluster2)]);
    
    % close all opened files
    fclose all;
end

