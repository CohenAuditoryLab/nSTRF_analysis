% function pairedSTRFanalysis(spike_clusters,cluster1,cluster2,time_window,sprfile,Trig)
%
%   FILE NAME   : pairedSTRFanalysis.m
%   DESCRIPTION : This file computes the significant STRFs of spike times
%       from 2 clusters, the OR/AND STRF btw the two STRFs,
%       coincident spike times at calculated optimal bin size, and 
%       cross correlation btw coincidents/OR-AND STRF 
%
% INPUT PARAMS
%   spikeClusters   : full path to cluster.mat file containing all spikes
%   cluster1        : index of first cluster to get spike series from
%   cluster2        : index of second cluster to get spike series from
%   sprfile         : full path to spectral profile file
%   Trig            : full path to trigger data
%   version         : 'st_clu' or 'spikeTimeRip', which struct is in 
%                       spikeClusters
%
% RETURNED VARIABLES
%   N/A, data saved to workspace
%
% (C) Shannon Lin, Edited Nov 2019

% Note to self (run function as follows):
% pairedSTRFanalysis('/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat', 9, 9, '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/DNR_Cortex_96k5min_4_50.spr','/Users/shannon1/Documents/F19/neuroResearch/nSTRF/AudiResp_16_24-190326-154559_triggers.mat', 'st_clu')

function pairedSTRFanalysis(spikeClusters,cluster1,cluster2,sprfile,Trig,version)
    % Find and load spike times corresponding to cluster1 and cluster2 in spike_clusters
    spikeTimeRipClusStruct = load(spikeClusters);
    if (strcmp(version, 'spikeTimeRip'))
        allSpikeTimes = spikeTimeRipClusStruct.spikeTimeRipClus;
    elseif (strcmp(version, 'st_clu'))
        allSpikeTimes = spikeTimeRipClusStruct.st_clu;
    else
        disp('Invalid version name, please input "st_clu" or "spikeTimeRip"')
        return;
    end
    assignin('base', 'spikeTimeRipClus', allSpikeTimes);
    index1 = find(cell2mat(allSpikeTimes(:,1))==cluster1);
    index2 = find(cell2mat(allSpikeTimes(:,1))==cluster2);
    if (isempty(index1) || isempty(index2)) 
        disp('Cluster number specified does not exist in data, see spikeTimeRipClus in workspace')
        return
    end
    spike1 = allSpikeTimes{index1, 2};
    spike2 = allSpikeTimes{index2, 2};
    assignin('base', 'spike1', spike1);
    assignin('base', 'spike2', spike2);
    
    % Define STRF parameters
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
    s_bin=0.15;
    T1=0;
    T2=s_bin;
    Fss=24414.0625;
    spet1=spike1 * Fss;
    spet2=spike2 * Fss;
    assignin('base', 'spet1', spet1);
    assignin('base', 'spet2', spet2);
    SPL=80;
    MdB=30;
    ModType='dB';
    Sound='MR';
    NBlocks=100;
    UF=10;
    sprtype='float';
    
    % Find optimal bin size
    optimalBinSize = findBin(Fss, spet1, spet2);
    assignin('base', 'time_window', optimalBinSize);
    
    % Find coincident spike times for spike1 and spike2
    [coin1, coin2] = findOverlap2(spike1, spike2, optimalBinSize);

    % Compute STRF of each cluster
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
    
    % Average STRF1 from TrigA and TrigB for both spike trains
    clusOneSTRF = (clusOneSTRF1A+clusOneSTRF1B)/2;
    assignin('base', 'clusOneSTRF', clusOneSTRF);
    clusOneNo1 = clusOneNo1A + clusOneNo1B;
    clusOneWo1 = (clusOneWo1A + clusOneWo1B)/2;
    clusTwoSTRF = (clusTwoSTRF1A+clusTwoSTRF1B)/2;
    assignin('base', 'clusTwoSTRF', clusTwoSTRF);
    clusTwoNo1 = clusTwoNo1A + clusTwoNo1B;
    clusTwoWo1 = (clusTwoWo1A + clusTwoWo1B)/2;
     
    % Compute significant STRF of each cluster
    p=0.05;
    SModType='dB';
    [clusOneSTRF1s,clusOneTresh1]=wstrfstat(clusOneSTRF,p,clusOneNo1,clusOneWo1,clusOnePP,MdB,ModType,Sound,SModType);
    clusOneSTRF1s = binarizeSTRFs(clusOneSTRF1s);
    assignin('base', 'clusOneSTRF1s', clusOneSTRF1s);
    [clusTwoSTRF1s,clusTwoTresh1]=wstrfstat(clusTwoSTRF,p,clusOneNo1,clusOneWo1,clusOnePP,MdB,ModType,Sound,SModType);
    clusTwoSTRF1s = binarizeSTRFs(clusTwoSTRF1s);
    assignin('base', 'clusTwoSTRF1s', clusTwoSTRF1s);
    
    % Compute STRF coincident spike times
    if (~(isempty(coin1) && (isempty(coin2))))
        spet1=coin1 * Fss;
        spet2=coin2 * Fss;
        try
            [n1Taxis,n1Faxis,n1STRF1A,n1STRF2A,n1PP,n1Wo1A,n1Wo2A,n1No1A,n1No2A,n1SPLNA]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            [n2Taxis,n2Faxis,n2STRF1A,n2STRF2A,n2PP,n2Wo1A,n2Wo2A,n2No1A,n2No2A,n2SPLNA]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            assignin('base', 'n1STRF1A', n1STRF1A);
            assignin('base', 'n2STRF1A', n2STRF1A);
            [n1Taxis,n1Faxis,n1STRF1B,n1STRF2B,n1PP,n1Wo1B,n1Wo2B,n1No1B,n1No2B,n1SPLNB]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            [n2Taxis,n2Faxis,n2STRF1B,n2STRF2B,n2PP,n2Wo1B,n2Wo2B,n2No1B,n2No2B,n2SPLNB]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            assignin('base', 'n1STRF1B', n1STRF1B);
            assignin('base', 'n2STRF1B', n2STRF1B);
        catch me
            disp('Error when computing nSTRF1A/B, exiting function')
            disp(me);
            return;
        end
        % Average STRF1 from TrigA and TrigB for coincident spike times
        coin1STRF = (n1STRF1A+n1STRF1B)/2;
        assignin('base', 'coin1STRF', coin1STRF);
        n1No1 = n1No1A + n1No1B;
        n1Wo1 = (n1Wo1A + n1Wo1B)/2;
        
        coin2STRF = (n2STRF1A+n2STRF1B)/2;
        assignin('base', 'coin2STRF', coin2STRF);
        n2No1 = n2No1A + n2No1B;
        n2Wo1 = (n2Wo1A + n2Wo1B)/2;
    end
    
    % Compute significant coinSTRFs
    [coin1STRF1s,n1Tresh1]=wstrfstat(coin1STRF,p,n1No1,n1Wo1,n1PP,MdB,ModType,Sound,SModType);
    coin1STRF1s = binarizeSTRFs(coin1STRF1s);
    assignin('base', 'coin1STRF1s', coin1STRF1s);
    [coin2STRF1s,n2Tresh1]=wstrfstat(coin2STRF,p,n2No1,n2Wo1,n2PP,MdB,ModType,Sound,SModType);
    coin2STRF1s = binarizeSTRFs(coin2STRF1s);
    assignin('base', 'coin2STRF1s', coin2STRF1s);

    % Compute OR/AND STRF from STRFs
    orSTRF = clusOneSTRF1s | clusTwoSTRF1s;
    andSTRF = clusOneSTRF1s & clusTwoSTRF1s;
    orSTRF = double(orSTRF);
    andSTRF = double(andSTRF);
    assignin('base', 'orSTRF', orSTRF);
    assignin('base', 'andSTRF', andSTRF);
    
    % Compute cross correlation between coin1/2STRFs, OR/AND STRF
    coin1CorrOr = corr2(coin1STRF1s, orSTRF);
    coin1CorrAnd = corr2(coin1STRF1s, andSTRF);
    assignin('base', 'coin1CorrOr', coin1CorrOr);
    assignin('base', 'coin1Corr', coin1CorrAnd);
    % Index of zero lag is (0+rows, 0+cols) 
    zeroLagRowIndex = size(coin1CorrOr, 1);
    zeroLagColIndex = size(coin1CorrOr, 2);
    coin1CorrOrZeroLag = coin1CorrOr(zeroLagRowIndex, zeroLagColIndex);
    coin1CorrAndZeroLag = coin1CorrAnd(zeroLagRowIndex, zeroLagColIndex);
    assignin('base', 'coin1CorrOrZeroLag', coin1CorrOrZeroLag);
    assignin('base', 'coin1CorrAndZeroLag', coin1CorrAndZeroLag);
    coin2CorrOr = corr2(coin2STRF1s, orSTRF);
    coin2CorrAnd = corr2(coin2STRF1s, andSTRF);
    assignin('base', 'coin2CorrOr', coin2CorrOr);
    assignin('base', 'coin2Corr', coin2CorrAnd);
    coin2CorrOrZeroLag = coin2CorrOr(zeroLagRowIndex, zeroLagColIndex);
    coin2CorrAndZeroLag = coin2CorrAnd(zeroLagRowIndex, zeroLagColIndex);
    assignin('base', 'coin2CorrOrZeroLag', coin2CorrOrZeroLag);
    assignin('base', 'coin2CorrAndZeroLag', coin2CorrAndZeroLag);
    
    % plot STRF, STRFs, coinSTRF, coinSTRFs
    rows = 6;
    cols = 2;
    num = 1;
    figure();
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),clusOneSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster1) ' STRF']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(twoTaxis,log2(twoFaxis/twoFaxis(1)),clusTwoSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster2) ' STRF']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(n1Taxis,log2(n1Faxis/n1Faxis(1)),coin1STRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster1) ' coinSTRF']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(n2Taxis,log2(n2Faxis/n2Faxis(1)),coin2STRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster2) ' coinSTRF']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),clusOneSTRF1s);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster1) ' STRFs']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(twoTaxis,log2(twoFaxis/twoFaxis(1)),clusTwoSTRF1s);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster2) ' STRFs']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(n1Taxis,log2(n1Faxis/n1Faxis(1)),coin1STRF1s);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster1) ' coinSTRFs']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(n2Taxis,log2(n2Faxis/n2Faxis(1)),coin2STRF1s);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title(['#' int2str(cluster2) ' coinSTRFs']);
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),orSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title('orSTRF');
    
    subplot(rows,cols,num)
    num = num + 1;
    pcolor(oneTaxis,log2(oneFaxis/oneFaxis(1)),andSTRF);
    colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    title('andSTRF');
    
    xAxis = 0.1;
    yAxis = 1;
    offset = -0.3;
    corrOr1 = subplot(rows, cols, num);
    text(xAxis, yAxis,['coin' int2str(cluster1) 'CorrOrZeroLag: ' num2str(coin1CorrOrZeroLag)]);
    yAxis = yAxis + offset;
    set (corrOr1, 'visible', 'off')
    
    corrAnd1 = subplot(rows, cols, num);
    text(xAxis, yAxis,['coin' int2str(cluster1) 'CorrAndZeroLag: ' num2str(coin1CorrAndZeroLag)]);
    yAxis = yAxis + offset;
    set (corrAnd1, 'visible', 'off')
    
    corrOr2 = subplot(rows, cols, num);
    text(xAxis, yAxis,['coin' int2str(cluster2) 'CorrOrZeroLag: ' num2str(coin2CorrOrZeroLag)]);
    yAxis = yAxis + offset;
    set (corrOr2, 'visible', 'off')
    
    corrAnd2 = subplot(rows, cols, num);
    text(xAxis, yAxis,['coin' int2str(cluster2) 'CorrAndZeroLag: ' num2str(coin2CorrAndZeroLag)]);
    set (corrAnd2, 'visible', 'off')
  
    % close all opened files
    fclose all;
end

