% function [optimalBinSize] = findBin(spike1, spike2)
%
%   FILE NAME   : findBin.m
%   DESCRIPTION : finds optimal time bin size between 2 spike trains by
%       finding at what lag the peak cross correlation occurs at
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