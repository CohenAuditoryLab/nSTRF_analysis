% function [clusData] = getIndividual(Params, spikeTimeRipClus,sprfile,Trig,numClusters)
%
%   FILE NAME   : getIndividual.m
%   DESCRIPTION : This file save individual cluster parameters in
%       struct clusData
%       
%
% INPUT PARAMS
%   Params           : specified params to compute STRF
%   spikeTimeRipClus : struct containing spike times
%   sprfile          : full path to spectral profile file
%   Trig             : full path to trigger data
%   numClusters      : number of clusters to analyze
%
% RETURNED VARIABLES
%   clusData          : contains individual cluster indices data
%
% (C) Shannon Lin, Edited Dec 2019

function [clusData] = getIndividual(Params, spikeTimeRipClus,sprfile,Trig,numClusters)
    % Define STRF parameters
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
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
    
    % Store STRF, STRFs data of each cluster to a struct
    for k=1:numClusters
        % Compute spike event times
        spike = spikeTimeRipClus{k, 2};
        spet = spike * Fss;
        try
            [~,~,STRF1A,STRF2A,~,Wo1A,Wo2A,No1A,No2A,~]=rtwstrfdbint(sprfile,T1,T2,spet',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
        catch me
            disp('Error when computing STRF1/STRF2 of cluster %d, exiting function', spikeTimeRipClus{k, 1})
            disp(me);
            return;
        end
        
        % Average STRF from TrigA and TrigB
        STRF1 = (STRF1A+STRF1B)/2;
        No1 = No1A + No1B;
        Wo1 = (Wo1A + Wo1B)/2;
        STRF2 = (STRF2A+STRF2B)/2;
        No2 = No2A + No2B;
        Wo2 = (Wo2A + Wo2B)/2;
        
        % Compute significant STRF
        [STRF1s,Tresh1]=wstrfstat(STRF1,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
        [STRF2s,Tresh2]=wstrfstat(STRF2,p,No2,Wo2,PP,MdB,ModType,Sound,SModType);
        % Convert to 0/1s
        STRF1sBinary = binarizeSTRFs(STRF1s);
        STRF2sBinary = binarizeSTRFs(STRF2s);
        
        % Assign fields in struct
        clusParam.clusterNo = spikeTimeRipClus{k, 1};
        clusParam.spike = spike;
        clusParam.spet = spet;
        clusParam.taxis = taxis;
        clusParam.faxis = faxis;
        clusParam.PP = PP;
        clusParam.SPLN = SPLN;
        clusParam.No1 = No1;
        clusParam.Wo1 = Wo1;
        clusParam.STRF1A = STRF1A;
        clusParam.STRF1B = STRF1B;
        clusParam.STRF1 = STRF1;
        clusParam.STRF1s = STRF1s;
        clusParam.STRF1sBinary = STRF1sBinary;
        clusParam.Tresh1 = Tresh1;
        clusParam.No2 = No2;
        clusParam.Wo2 = Wo2;
        clusParam.STRF2A = STRF2A;
        clusParam.STRF2B = STRF2B;
        clusParam.STRF2 = STRF2;
        clusParam.STRF2s = STRF2s;
        clusParam.STRF2sBinary = STRF2sBinary;
        clusParam.Tresh2 = Tresh2;
        
        [RFParam]=strfparam(taxis,faxis,STRF1,Wo1,PP,Sound);
        for fn = fieldnames(RFParam)'
           clusParam.(fn{1}) = RFParam.(fn{1});
        end
        
        % Save in clusData
        clusData(k, 1) = clusParam;
    end
    assignin('base', 'clusData', clusData);
end

