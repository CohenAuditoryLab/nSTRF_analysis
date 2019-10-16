% function [allSpikeTimes, spike1, spike2, overlap]=findOverlap(spikeClusters, cluster1, cluster2, time_window, version)
%
%   FILE NAME   : findOverlap.m
%   DESCRIPTION : overlap spike times from 2 clusters of spike trains
%
% INPUT PARAMS
%   spike1              : first spike time series
%   spike2              : second spike time series
%
% RETURNED VARIABLES
%   optimalBinSize      : optimal bin size for finding coincident spike
%                           times
%
% (C) Shannon Lin, Edited Oct 2019

function [optimalBinSize] = findBin(spike1, spike2)
    [corr, lags] = xcorr(spike1, spike2);
%     assignin('base', 'corr', corr);
%     assignin('base', 'lags', lags);
%     plot(lags, corr);
    [~, peakIndex] = max(corr);
    optimalBinSize = abs(lags(1, peakIndex) / 1e3);
end