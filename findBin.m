% function [optimalBinSize] = findBin(Fss, spet1, spet2)
%
%   FILE NAME   : findBin.m
%   DESCRIPTION : finds optimal time bin size between 2 spike event times
%       finding at what lag the peak cross correlation occurs at
%
% INPUT PARAMS
%   Fss                : sampling rate
%   spet1              : first spike event time series
%   spet2              : second spike event time series
%
% RETURNED VARIABLES
%   R                   : correlation matrix
%   optimalBinSize      : optimal bin size for finding coincident spike
%                           times

function [R, optimalBinSize] = findBin(Fss, spet1, spet2)
    % Monty's code:
    % The desired sampling rate for correlation analysis - the spike train 
    % is resampled at Fsd prior to computing the correlation
    % Can change to 250 for 4 ms resolution
    Fsd=1000;    

    % Correlation lag - up to 250 ms
    MaxTau = 250;  
    % Experiment duration
    T = 600;       

    % Don't remove corr value at zero lag - only useful for autocorrelograms
    Zero='n'; 
    % Don't remove mean correlation value - DC
    Mean='n';  
    Disp='n';
    
    % R=correlation matrix, optimalBinSize time of maximum peak
    [R, optimalBinSize] = xcorrspike(spet1,spet2,Fss,Fsd,MaxTau,T,Zero,Mean,Disp);

    % Shannon's optimalBinSize calculation, using Monty's function instead^
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