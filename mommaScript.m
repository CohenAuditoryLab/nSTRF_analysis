% function [clusData, pairedData] = mommaScript(spikeClusters,sprfile,Trig,version,numClusters)
%
%   FILE NAME   : mommaScript.m
%   DESCRIPTION : This file collects individual and pairwise STRF indices
%       
%
% INPUT PARAMS
%   spikeClusters : full path to .mat containing spike times
%   sprfile       : full path to spectral profile file
%   Trig          : full path to trigger data
%   version       : 'st_clu' or 'spikeTimeRip', which struct is in 
%                       spikeClusters (Cassius_190326 is "st_clu",
%                       Cassius_190324 is "spikeTimeRip")
%   numClusters   : optional argument, can specify number of clusters 
%                       to analyze, default is analyze all clusters in .mat
%
% RETURNED VARIABLES
%   clusData          : contains individual cluster indices data
%   pairedData        : contains nSTRF indices data
%
% (C) Shannon Lin, Edited Dec 2019

% Note to self, run as follows:
% mommaScript('/Users/shannon1/Documents/Research/Cohen/nSTRF/spike_times_ripple_clust_new.mat', '/Users/shannon1/Documents/Research/Cohen/nSTRF/DNR_Cortex_96k5min_4_50.spr','/Users/shannon1/Documents/Research/Cohen/nSTRF/AudiResp_16_24-190326-154559_triggers.mat', 'st_clu', 2)
function [clusData, pairedData] = mommaScript(spikeClusters,sprfile,Trig,version,numClusters)
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
    assignin('base', 'Params', Params);

    % Load in spikeClusters and compute total num of clusters in there
    spikeTimeRipClusStruct = load(spikeClusters);
    % Cassius_190324 saved as "spikeTimeRip"
    % Cassius_190326 data saved as "st_clu" 
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

    % number of clusters to analyze is total number of clusters in
    % spikeClusters
    if nargin<5
         numClusters = size(eval(field),1);
    end
    numClusters = int8(numClusters);
    
    % collect indices of individual clusters
    clusData = getIndividual(Params, spikeTimeRipClus, sprfile, Trig, numClusters);
    % collect indices of nSTRF pairwise clusters
    pairedData = getPaired(Params, spikeTimeRipClus, sprfile, Trig, numClusters, clusData,'n');
end