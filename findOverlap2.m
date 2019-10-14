% function [allSpikeTimes, spike1, spike2, coin1, coin2]=findOverlap2(spikeClusters, cluster1, cluster2, time_window, version)
%
%   FILE NAME   : findOverlap2.m
%   DESCRIPTION : overlap spike times from 2 clusters of spike trains
%
% INPUT PARAMS
%   spike1          : first spike time series
%   spike2          : second spike time series
%   time_window     : optimal bin size to find coincident spike times
%
% RETURNED VARIABLES
%   coin1               : spike1 coincident spike times with spike2
%   coin2               : spike2 coincident spike times with spike1
%
% (C) Shannon Lin, Edited Oct 2019

function [coin1, coin2]=findOverlap2(spike1, spike2, time_window)
    % find total number of time_window chunks to loop through
    max1 = max(spike1);
    max2 = max(spike2);
    if (max1 >= max2) 
        largestTime = max1;
    else
        largestTime = max2;
    end
    numIterations = uint64(largestTime/time_window + 1);
    [dim1, ~] = size(spike1);
    [dim2, ~] = size(spike2);
    coin1 = zeros(dim1, 1);
    coin2 = zeros(dim2, 1);
    
    % mark respective spike times if both fired within same time window
    for i = 1:numIterations
        lb = time_window*(double(i)-1.0);
        ub = time_window*double(i);
        spike1Contains = find(spike1 <= ub & spike1 >= lb);
        spike2Contains = find(spike2 <= ub & spike2 >= lb);
        if (~isempty(spike1Contains) && ~isempty(spike2Contains))
            coin1 = cat(1, coin1, findAllCoins(spike1, lb, ub));
            coin2 = cat(1, coin2, findAllCoins(spike2, lb, ub));
        end
    end
    
    coin1 = coin1(coin1 ~= 0);
    coin2 = coin2(coin2 ~= 0);
    assignin('base', 'coin1', coin1);
    assignin('base', 'coin2', coin2);
    
    function coincidentSpikes = findAllCoins(spike, lb, ub)
        [dim, ~] = size(spike);
        coincidentSpikes = zeros(dim, 1);
        for j = 1:dim
            if (spike(j) >= lb && spike(j) <= ub)
                coincidentSpikes(j) = spike(j);
            end
        end
        coincidentSpikes = coincidentSpikes(coincidentSpikes ~= 0);
    end
end

