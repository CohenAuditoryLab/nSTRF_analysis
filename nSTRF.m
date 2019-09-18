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
%
%RETURNED VARIABLES
%
%   
%
% (C) Shannon Lin, Edited Sept 2019

function nSTRF(spike_clusters, cluster1, cluster2, time_window)
    % load spike trains from index cluster1 and index cluster2 from
    % spike series in specified input spikeClusters
    spikeTimeRipClusStruct = load(spike_clusters);
    spikeTimeRipClus = spikeTimeRipClusStruct.spikeTimeRipClus;
    spike1 = spikeTimeRipClus{cluster1, 2};
    spike2 = spikeTimeRipClus{cluster2, 2};
    % Compute STRF of each neuron
    %   NEED INPUT PARAMS to rtwstrfdbint function
    %   Should I load them in as an input argument to nSTRF function?
    %   I think spike1 and spike2 should be the respective "spet" input vars?
    try
        [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint(SpecFile,T1,T2,spet,Trig,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
    catch me
        STRF1A=[];STRF2A=[];No1A=[];Wo1A=[];No1A=[];Wo1A=[];
    end

    try
        [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(SpecFile,T1,T2,spet,Trig,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
    catch me
        STRF1B=[];STRF2B=[];No1B=[];Wo1B=[];No1B=[];Wo1B=[];
    end
    if isempty(STRF1A) && isempty(STRF1B)
        disp('Error when computing STRF of each neuron, exiting function')
        return;
    end
    STRF1 = (STRF1A+STRF1B)/2;
    STRF2 = (STRF2A+STRF2B)/2;
    No1   =  No1A+No1B;
    Wo1   =  (Wo1A+Wo1B)/2;
    No2   =  No2A+No2B;
    Wo2   =  (Wo2A+Wo2B)/2;
    
    % Compute significant STRF of each neuron + track if pixels are
    % excitatory or inhibitatory
    %   NEED MORE INPUT PARAMS to calculate STRFs, probably same as the ones used when
    %   calculating STRF?
    %   Excitatory or inhibitatory, what's that mean?
    [STRF1s,Tresh1]=wstrfstat(STRF1,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
    [STRF1s,Tresh2]=wstrfstat(STRF2,p,No2,Wo2,PP,MdB,ModType,Sound,SModType);

    % Compute nSTRF using different bin sizes
    %   What does this mean?
    %   From findOverlap function, I have overlap spike times
    overlapSpike = findOverlap(spikes, cluster1, cluster2, time_window);

    % Compute from individual STRFs the OR and AND STRF
    %   Is it the logical OR btw pixels of STRF1s and STRF2s?
    %   Similarly, are we calculating logical AND btw pixels of STRF1s and STRF2s?
    
    % Compare OR and AND STRF with nSTRF
    %   Compare as in cross correlation?
    
    % Compute the cross correlation between the 2 spike trains and calc
    % significance of corr coeff (calculate the CC between 0-10 ms)
    %   What are the 2 spike trains? Is it CC btw STRF1s/STRF2s?
    %   Or is it CC btw nSTRF/OR STRF and btw nSTRF/AND STRF?
 
end