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

function [optimalBinSize] = findBin(Fss, spet1, spet2)
    % Original sampling rate - ideally use the sampling rate for your data
    % Fs=12000; 
    Fs = Fss;
    % The desired sampling rate for correlation analysis - the spike train 
    % is resampled at Fsd prior to computing the correlation
    % Can change to 250 for 4 ms resolution
    Fsd=1000;    

    MaxTau = 250;  %Correlation lag - up to 250 ms
    T = 600;       %Experiment duration - Im assuming this is the correct value

    Zero='n';   % Do not remove correlation value at zero lag - only useful for autocorrelograms
    Mean='n';   %Do not remove mean correaltion value - DC
    Disp='y';
    
    % R=correlation matrix, optimalBinSize time of maximum peak
    [R, optimalBinSize] = xcorrspike(spet1,spet2,Fs,Fsd,MaxTau,T,Zero,Mean,Disp);
    assignin('base', 'R', R);
%     [corr, lags] = xcorr(spike1, spike2);
%     assignin('base', 'corr', corr);
%     assignin('base', 'lags', lags);
%     [~, peakIndex] = max(corr);
%     optimalBinSize = abs(lags(1, peakIndex) / 1e3);
%     
%     plot(lags, corr);
%     title('xcorr between spike trains 6 and 9 from Cassius190326')
%     xlabel('lags (in ms)')
%     ylabel('corr')
%     text(-1500, 1.6e8,['optimalBinSize = ' num2str(optimalBinSize) 'sec']);
end