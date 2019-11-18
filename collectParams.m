% function [clusData] = collectParams(spikeClusters,sprfile,Trig,version)
%
%   FILE NAME   : collectParams.m
%   DESCRIPTION : This file save various cluster parameters in nxn matrix
%       where n = number of clusters in inputted spikeClusters .mat file;
%       (i,i) in outputted matrix has fields from strfparam, (i, j)
%       contains pairedSTRFanalysis data with cluster i and j
%       
%
% INPUT PARAMS
%   spikeClusters   : full path to cluster.mat file containing all spikes
%   sprfile         : full path to spectral profile file
%   Trig            : full path to trigger data
%   version         : 'st_clu' or 'spikeTimeRip', which struct is in 
%                       spikeClusters
%
% RETURNED VARIABLES
%   N/A

% VARIABLES SAVED TO WORKSPACE
% clusData          : contains individual params
% pairData          : contains cross correlation parameters
%
% (C) Shannon Lin, Edited Nov 2019

% Tested function as follows:
% collectParams('/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat', '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/DNR_Cortex_96k5min_4_50.spr','/Users/shannon1/Documents/F19/neuroResearch/nSTRF/AudiResp_16_24-190326-154559_triggers.mat', 'st_clu')

function collectParams(spikeClusters,sprfile,Trig,version)

    % Load in spikeClusters and compute num of clusters are in there
    spikeTimeRipClusStruct = load(spikeClusters);
    if (strcmp(version, 'spikeTimeRip'))
        spikeTimeRipClus = spikeTimeRipClusStruct.spikeTimeRipClus;
    elseif (strcmp(version, 'st_clu'))
        spikeTimeRipClus = spikeTimeRipClusStruct.st_clu;
    else
        disp('Invalid version name, please input "st_clu" or "spikeTimeRip"')
        return;
    end
    assignin('base', 'spikeTimeRipClus', spikeTimeRipClus);
    field = sprintf('spikeTimeRipClusStruct.%s', version);
    numClusters = size(eval(field),1);
    
    % Define STRF parameters
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
    s_bin=0.15;
    T1=0;
    T2=s_bin;
    Fss=24414.0625;
    SPL=80;
    MdB=30;
    ModType='dB';
    Sound='MR';
    NBlocks=100;
    UF=10;
    sprtype='float';
    % Define STRFs params
    p=0.05;
    SModType='dB';
    
    % Store STRF, STRFs data of each cluster to a struct
    for k=1:3
        % Compute spike event times
        spike = spikeTimeRipClus{k, 2};
        spet = spike * Fss;
        try
            [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
        catch me
            disp('Error when computing STRF1/STRF2 of cluster %d, exiting function', spikeTimeRipClus{k, 1})
            disp(me);
            return;
        end
        
        % Average STRF1 from TrigA and TrigB
        STRF = (STRF1A+STRF1B)/2;
        No1 = No1A + No1B;
        Wo1 = (Wo1A + Wo1B)/2;
        
        % Compute significant STRF
        [STRF1s,Tresh1]=wstrfstat(STRF,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
        % Convert to 0/1s
        STRF1s = binarizeSTRFs(STRF1s);
        
        % Assign fields in struct
        clusParam.num = spikeTimeRipClus{k, 1};
        clusParam.spike = spike;
        clusParam.spet = spet;
        clusParam.STRF1s = STRF1s;
        [RFParam]=strfparam(taxis,faxis,STRF,No1,PP,Sound);
        for fn = fieldnames(RFParam)'
           clusParam.(fn{1}) = RFParam.(fn{1});
        end
        
        % Save in clusData
        clusData(k, 1) = clusParam;
    end
    assignin('base', 'clusData', clusData);
    
    for i=1:2
        for j=1:2
            if (i == j)
                pairedClusData(i, j) = clusData(i, 1);
            else
                pairData.clusOne = spikeTimeRipClus{i, 1};
                pairData.clusOne = spikeTimeRipClus{j, 1};
                struct1 = clusData(i, 1);
                struct2 = clusData(j, 1);
                % Find optimal bin size
                optimalBinSize = findBin(Fss, struct1.spet, struct2.spet);
                pairData.optimalBinSize = optimalBinSize;
                % Find coincident spike times for spike1 and spike2
                [coin1, coin2] = findOverlap2(struct1.spike, struct2.spike, optimalBinSize);
                pairData.coin1 = coin1;
                pairData.coin2 = coin2;

                % Compute STRF coincident spike times
                if (~(isempty(coin1) && (isempty(coin2))))
                    spet1=coin1 * Fss;
                    spet2=coin2 * Fss;
                    try
                        [n1Taxis,n1Faxis,n1STRF1A,n1STRF2A,n1PP,n1Wo1A,n1Wo2A,n1No1A,n1No2A,n1SPLNA]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                        [n1Taxis,n1Faxis,n1STRF1B,n1STRF2B,n1PP,n1Wo1B,n1Wo2B,n1No1B,n1No2B,n1SPLNB]=rtwstrfdbint(sprfile,T1,T2,spet1',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                        [n2Taxis,n2Faxis,n2STRF1A,n2STRF2A,n2PP,n2Wo1A,n2Wo2A,n2No1A,n2No2A,n2SPLNA]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                        [n2Taxis,n2Faxis,n2STRF1B,n2STRF2B,n2PP,n2Wo1B,n2Wo2B,n2No1B,n2No2B,n2SPLNB]=rtwstrfdbint(sprfile,T1,T2,spet2',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                    catch me
                        disp('Error when computing nSTRF1A/B, exiting function')
                        disp(me);
                        return;
                    end
                    % Average STRF1 from TrigA and TrigB for coincident spike times
                    coin1STRF = (n1STRF1A+n1STRF1B)/2;
                    n1No1 = n1No1A + n1No1B;
                    n1Wo1 = (n1Wo1A + n1Wo1B)/2;
                    pairData.coin1STRF = coin1STRF;
                    
                    coin2STRF = (n2STRF1A+n2STRF1B)/2;
                    n2No1 = n2No1A + n2No1B;
                    n2Wo1 = (n2Wo1A + n2Wo1B)/2;
                    pairData.coin2STRF = coin2STRF;
                end

                % Compute significant coinSTRFs
                [coin1STRFs,n1Tresh1]=wstrfstat(coin1STRF,p,n1No1,n1Wo1,n1PP,MdB,ModType,Sound,SModType);
                coin1STRFs = binarizeSTRFs(coin1STRFs);
                pairData.coin1STRFs = coin1STRFs;
                [coin2STRFs,n2Tresh1]=wstrfstat(coin2STRF,p,n2No1,n2Wo1,n2PP,MdB,ModType,Sound,SModType);
                coin2STRFs = binarizeSTRFs(coin2STRFs);
                pairData.coin2STRFs = coin2STRFs;

                [~, index1] = find([clusData.num] == i);
                [~, index2] = find([clusData.num] == j);
                clusOneSTRF1s = clusData(index1).STRF1s;
                clusTwoSTRF1s = clusData(index2).STRF1s;
                % Compute OR/AND STRF from STRFs
                orSTRF = clusOneSTRF1s | clusTwoSTRF1s;
                andSTRF = clusOneSTRF1s & clusTwoSTRF1s;
                orSTRF = double(orSTRF);
                andSTRF = double(andSTRF);
                pairData.orSTRF = orSTRF;
                pairData.andSTRF = andSTRF;

                % Compute cross correlations (STRF1/OR,STRF2/OR,STRF1/AND,STRF2/AND)
                coin1CorrOr = corr2(coin1STRFs, orSTRF);
                coin1CorrAnd = corr2(coin1STRFs, andSTRF);
                
                % Index of zero lag is (0+rows, 0+cols) 
                zeroLagRowIndex = size(coin1CorrOr, 1);
                zeroLagColIndex = size(coin1CorrOr, 2);
                coin1CorrOrZeroLag = coin1CorrOr(zeroLagRowIndex, zeroLagColIndex);
                coin1CorrAndZeroLag = coin1CorrAnd(zeroLagRowIndex, zeroLagColIndex);
                pairData.coin1CorrOrZeroLag = coin1CorrOrZeroLag;
                pairData.coin1CorrAndZeroLag = coin1CorrAndZeroLag;
                
                coin2CorrOr = corr2(coin2STRFs, orSTRF);
                coin2CorrAnd = corr2(coin2STRFs, andSTRF);
                coin2CorrOrZeroLag = coin2CorrOr(zeroLagRowIndex, zeroLagColIndex);
                coin2CorrAndZeroLag = coin2CorrAnd(zeroLagRowIndex, zeroLagColIndex);
                pairData.coin2CorrOrZeroLag = coin2CorrOrZeroLag;
                pairData.coin2CorrAndZeroLag = coin2CorrAndZeroLag;
                
                coin1CorrCoin2 = corr2(coin1STRFs, coin2STRFs);
                coin1CorrCoin2ZeroLag = coin1CorrCoin2(zeroLagRowIndex, zeroLagColIndex);
                pairData.coin1CorrCoin2ZeroLag = coin1CorrCoin2ZeroLag;

                % Using Monty's code to calculate xcorr
                RSTRF = strfcorrcorrected(clusOneSTRF1A,clusOneSTRF1B,clusOneSTRF1s,clusTwoSTRF1A,clusTwoSTRF1B,clusTwoSTRF1s,n2Taxis,n2Faxis,n2PP);
                montyR = RSTRF.R;
                montyCoin1CorrCoin2ZeroLag = montyR(zeroLagRowIndex, zeroLagColIndex);
                pairData.montyCoin1CorrCoin2ZeroLag = montyCoin1CorrCoin2ZeroLag;
                pairData.MontyRSTRF = RSTRF;
                
                pairedClusData(i, j) = pairData;
            end
        end
    end
   
end

