% function [pairedData] = getPaired(Params, spikeTimeRipClus,sprfile,Trig,numClusters)
%
%   FILE NAME   : getPaired.m
%   DESCRIPTION : This file saves pairwise analysis data to a 
%       pairedData struct database
%       
%
% INPUT PARAMS
%   Params           : specified params to compute STRF
%   spikeTimeRipClus : struct containing spike times
%   sprfile          : full path to spectral profile file
%   Trig             : full path to trigger data
%   numClusters      : number of clusters to analyze
%   clusData         : contains individual cluster indices data
%
% RETURNED VARIABLES
%   pairedData          : contains nSTRF indices data
%       
%
% (C) Shannon Lin, Edited Dec 2019

function [pairedData] = getPaired(Params, spikeTimeRipClus,sprfile,Trig,numClusters,clusData)
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
    
    % loading params from input
    T1=Params.T1;
    T2=Params.T2;
    Fss=Params.Fss;
    SPL=Params.SPL;
    MdB=Params.MdB;
    ModType=Params.ModType;
    Sound=Params.Sound;
    NBlocks=Params.NBlocks;
    UF=Params.UF;
    sprtype=Params.sprtype;
    p=Params.p;
    SModType=Params.SModType;
    
    entered = 0;
    numIter = 1;
    for i=1:numClusters
        for j=1:numClusters
            if (i == j)
                continue;
            else
                clusNo1 = spikeTimeRipClus{i, 1};
                clusNo2 = spikeTimeRipClus{j, 1};
                if (entered == 1)
                    % skip if pair already analyzed
                    allClusOnes = [pairedData.clusOne];
                    allClusTwos = [pairedData.clusTwo];
                    index1 = find(allClusOnes == clusNo2);
                    index2 = find(allClusTwos == clusNo1);
                    if (index1 > 0)
                        if (index2 > 0)
                            continue;
                        end
                    end
                end
                entered = 1;
                pairData.clusOne = clusNo1;
                pairData.clusTwo = clusNo2;
                % Retrieve STRF params from clusData struct
                struct1 = clusData(i, 1);
                struct2 = clusData(j, 1);
                
                % Find optimal bin size to calculate coincident spike times
                optimalBinSize = findBin(Fss, struct1.spet, struct2.spet);
                pairData.optimalBinSize = optimalBinSize;
                
                % Find coincident spike times for spike1 and spike2
                [coin1, coin2] = findOverlap2(struct1.spike, struct2.spike, optimalBinSize);
                
                % Compute STRF coincident spike times
                if (~(isempty(coin1) && (isempty(coin2))))
                    spet1=(coin1 * Fss)';
                    spet2=(coin2 * Fss)';
                    nSTRFOne.spet = spet1;
                    nSTRFTwo.spet = spet2;
                    try
                        [~,~,n1STRF1A,~,~,n1Wo1A,~,n1No1A,~,SPLN1]=rtwstrfdbint(sprfile,T1,T2,spet1,TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                        [n1Taxis,n1Faxis,n1STRF1B,~,n1PP,n1Wo1B,~,n1No1B,~,~]=rtwstrfdbint(sprfile,T1,T2,spet1,TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                        [~,~,n2STRF1A,~,~,n2Wo1A,~,n2No1A,~,SPLN2]=rtwstrfdbint(sprfile,T1,T2,spet2,TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                        [n2Taxis,n2Faxis,n2STRF1B,~,n2PP,n2Wo1B,~,n2No1B,~,~]=rtwstrfdbint(sprfile,T1,T2,spet2,TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
                    catch me
                        X = ['Error when computing nSTRF1A/B btw clus no ', int2str(i), ' and ', int2str(j)];
                        disp(X)
                        disp(me);
                        continue;
                    end
                    % Average STRF1 from TrigA and TrigB for coincident spike times
                    coin1STRF = (n1STRF1A+n1STRF1B)/2;
                    n1Wo1 = (n1Wo1A + n1Wo1B)/2;
                    n1No1 = n1No1A + n1No1B;
                    % save nSTRF params for both coincident spikes
                    nSTRFOne.STRF = coin1STRF;
                    nSTRFOne.taxis = n1Taxis;
                    nSTRFOne.faxis = n1Faxis;
                    nSTRFOne.PP = n1PP;
                    nSTRFOne.Wo1 = n1Wo1;
                    nSTRFOne.Wo1 = n1Wo1;
                    nSTRFOne.No1 = n1No1;
                    nSTRFOne.SPLN = SPLN1;
                    coin2STRF = (n2STRF1A+n2STRF1B)/2;
                    n2No1 = n2No1A + n2No1B;
                    n2Wo1 = (n2Wo1A + n2Wo1B)/2;
                    nSTRFTwo.STRF = coin2STRF;
                    nSTRFTwo.taxis = n2Taxis;
                    nSTRFTwo.faxis = n2Faxis;
                    nSTRFTwo.PP = n2PP;
                    nSTRFTwo.Wo1 = n2Wo1;
                    nSTRFTwo.Wo1 = n2Wo1;
                    nSTRFTwo.No1 = n2No1;
                    nSTRFTwo.SPLN = SPLN2;
                end 

                % Save indices for coin1STRF and coin2STRF
                [RFParam1]=strfparam(n1Taxis,n1Faxis,coin1STRF,n1No1,n1PP,Sound);
                for fn = fieldnames(RFParam1)'
                   nSTRFOne.(fn{1}) = RFParam1.(fn{1});
                end
                [RFParam2]=strfparam(n2Taxis,n2Faxis,coin2STRF,n2No1,n2PP,Sound);
                for fn = fieldnames(RFParam2)'
                   nSTRFTwo.(fn{1}) = RFParam2.(fn{1});
                end
                
                % Compute significant coinSTRFs
                [coin1STRFs,Tresh1]=wstrfstat(coin1STRF,p,n1No1,n1Wo1,n1PP,MdB,ModType,Sound,SModType);
                coin1STRFsBinary = binarizeSTRFs(coin1STRFs);
                nSTRFOne.STRFs = coin1STRFs;
                nSTRFOne.Tresh1 = Tresh1; 
                nSTRFOne.STRFsBinary = coin1STRFsBinary;
                [coin2STRFs,Tresh2]=wstrfstat(coin2STRF,p,n2No1,n2Wo1,n2PP,MdB,ModType,Sound,SModType);
                coin2STRFsBinary = binarizeSTRFs(coin2STRFs);
                nSTRFTwo.STRFs = coin2STRFs;
                nSTRFTwo.Tresh2 = Tresh2; 
                nSTRFTwo.STRFsBinary = coin2STRFsBinary;
                
                % Save struct of cluster one and two data
                pairData.nSTRFOne = nSTRFOne;
                pairData.nSTRFTwo = nSTRFTwo;
                
                % retrieve STRFs from clusData struct
                [~, index1] = find([clusData.clusterNo] == clusNo1);
                [~, index2] = find([clusData.clusterNo] == clusNo2);

                clusOneSTRF1s = clusData(index1).STRF1sBinary;
                clusTwoSTRF1s = clusData(index2).STRF1sBinary;
                
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
                % Compute cross corrs (coin1/OR, coin2/OR, coin1/coin2)
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

                % Using Monty's code to calculate coin1/coin2 xcorr
                % (TO DO: why is it a diff # than what I get?)
                clusOneSTRF1A = clusData(index1).STRF1A;
                clusOneSTRF1B = clusData(index1).STRF1B;
                clusTwoSTRF1A = clusData(index2).STRF1A;
                clusTwoSTRF1B = clusData(index2).STRF1B;
                RSTRF = strfcorrcorrected(clusOneSTRF1A,clusOneSTRF1B,clusOneSTRF1s,clusTwoSTRF1A,clusTwoSTRF1B,clusTwoSTRF1s,n2Taxis,n2Faxis,n2PP);
                montyR = RSTRF.R;
                montyCoin1CorrCoin2ZeroLag = montyR(zeroLagRowIndex, zeroLagColIndex);
                pairData.montyCoin1CorrCoin2ZeroLag = montyCoin1CorrCoin2ZeroLag;
                pairData.MontyRSTRF = RSTRF;

                % Integrating Monty's ShuffleXCorr.m code
                STRF1s = struct1.STRF1s;
                STRF2s = struct2.STRF1s;
                taxis = pairData.nSTRFTwo.taxis;
                faxis = pairData.nSTRFTwo.faxis;
                PP = pairData.nSTRFTwo.PP;
                % Predict Correlation From the STRFs - Chen 2012
                [T,R,Rcc,RR]=strf2xcorr(taxis,faxis,STRF1s,STRF2s,PP,'y');
                pairData.T = T;
                pairData.R = R;
                pairData.Rcc = Rcc;
                pairData.RR = RR;
                pairedData(numIter, 1) = pairData;
                numIter = numIter + 1;
                assignin('base', 'pairedData', pairedData);
            end
        end
    end
end