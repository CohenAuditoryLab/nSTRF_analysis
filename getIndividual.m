% function [clusData] = getIndividual(Params, spikeTimeRipClus,sprfile,Trig,numClusters)
%
%   FILE NAME   : getIndividual.m
%   DESCRIPTION : This file save individual cluster parameters in
%       struct clusData
%       
%
% INPUT PARAMS
%   Params           : struct of STRF params, must include
%                           T1,T2,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,
%                           sprtype,p,SModType
%   spikeTimeRipClus : struct containing spike times
%   sprfile          : full path to spectral profile file
%   Trig             : full path to trigger data
%   numClusters      : number of clusters to analyze
%
% RETURNED VARIABLES
%   clusData          : contains individual cluster indices data
%
% (C) Shannon Lin, Edited Dec 2019

function [clusData] = getIndividual(Params,spikeTimeRipClus,sprfile,Trig,numClusters)
    % Retrieve trig
    trigStruct = load(Trig);
    TrigA = trigStruct.TrigA;
    TrigB = trigStruct.TrigB;
    % Define STRF parameters
    if (~isfield(Params,'T1') || ~isfield(Params,'T2') || ~isfield(Params,'Fss') || ~isfield(Params,'SPL') ...
        || ~isfield(Params,'MdB') || ~isfield(Params,'ModType') || ~isfield(Params,'Sound') || ...
        ~isfield(Params,'NBlocks') || ~isfield(Params,'UF') || ~isfield(Params,'sprtype') || ...
        ~isfield(Params,'p') || ~isfield(Params,'SModType'))
        disp('Check specifications above for necessary Params needed to calculate STRF, exiting func')
        return; 
    end
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
    for i=1:numClusters
        % Compute spike event times
        spike = spikeTimeRipClus{i, 2};
        spet = spike * Fss;
        try
            [~,~,STRF1A,STRF2A,~,Wo1A,Wo2A,No1A,No2A,~]=rtwstrfdbint(sprfile,T1,T2,spet',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
            [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
        catch me
            disp('Error when computing STRF1/STRF2 of cluster %d, exiting function', spikeTimeRipClus{i, 1})
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
        clusParam.clusterNo = spikeTimeRipClus{i, 1};
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

        % Save params from strfparam
        [RFParam]=strfparam(taxis,faxis,STRF1,Wo1,PP,Sound);
        for fn = fieldnames(RFParam)'
           clusParam.(fn{1}) = RFParam.(fn{1});
        end
        
        % Collecting indices from gstrfmodel.m
        % Errs parameter squiggled out, gstrfmodel.m has error in code
        % gstrfmodel.m for STRF1
        [STRF1m,STRF1am,STRF1bm,STRF1cm,x0,w,sf0,spectrop,t0,c,tf0,q,k,belta,SIs,SIt,SI,~,alpha_d,~]=gstrfmodel(STRF1,taxis,faxis,PP,Tresh1,'nsvd','nsvd','n');
        clusParam.STRF1m = STRF1m;
        clusParam.STRF1am = STRF1am;
        clusParam.STRF1bm = STRF1bm;
        clusParam.STRF1cm = STRF1cm;
        clusParam.x0_1 = x0;
        clusParam.w_1 = w;
        clusParam.sf0_1 = sf0;
        clusParam.spectrop_1 = spectrop;
        clusParam.t0_1 = t0;
        clusParam.c_1 = c;
        clusParam.tf0_1 = tf0;
        clusParam.q_1 = q;
        clusParam.k_1 = k;
        clusParam.belta_1 = belta;
        clusParam.SIs_1 = SIs;
        clusParam.SIt_1 = SIt;
        clusParam.SI_1 = SI;
        clusParam.alpha_d_1 = alpha_d;
        
        % gstrmodel.m for STRF2
        [STRF2m,STRF2am,STRF2bm,STRF2cm,x0,w,sf0,spectrop,t0,c,tf0,q,k,belta,SIs,SIt,SI,~,alpha_d,N]=gstrfmodel(STRF2,taxis,faxis,PP,Tresh2,'nsvd','nsvd','n');
        clusParam.STRF2m = STRF2m;
        clusParam.STRF2am = STRF2am;
        clusParam.STRF2bm = STRF2bm;
        clusParam.STRF2cm = STRF2cm;
        clusParam.x0_2 = x0;
        clusParam.w_2 = w;
        clusParam.sf0_2 = sf0;
        clusParam.spectrop_2 = spectrop;
        clusParam.t0_2 = t0;
        clusParam.c_2 = c;
        clusParam.tf0_2 = tf0;
        clusParam.q_2 = q;
        clusParam.k_2 = k;
        clusParam.belta_2 = belta;
        clusParam.SIs_2 = SIs;
        clusParam.SIt_2 = SIt;
        clusParam.SI_2 = SI;
        clusParam.alpha_d_2 = alpha_d;
        
        % Save in clusData
        % Can't figure out how to preallocate for speed
        clusData(i, 1) = clusParam;
    end
    assignin('base', 'clusData', clusData);
end

