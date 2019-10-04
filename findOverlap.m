% function [overlap]=findOverlap(spikes, cluster1, cluster2, time_window)
%
%   FILE NAME   : findOverlap.m
%   DESCRIPTION : overlap spike times from 2 clusters of spike trains
%
% INPUT PARAMS
%   spikeClusters   : full path to cluster.mat variable file containing all spikes
%   cluster1        : first cluster to get spike series from
%   cluster2        : second cluster to get spike series from
%   time_window     : fixed time window to chunk spikes with
%   version         : 'st_clu' or 'spikeTimeRip', refers to what struct exists in
%                       spike_clusters
%
% RETURNED VARIABLES
%
%   spikeTimeRipClus    : mx2 cell that holds cluster spike times 
%   spike1              : spike time at specified cluster 1
%   spike2              : spike time at specified cluster 2
%   overlap             : mxn double array that contains overlap times of spikes
%
% (C) Shannon Lin, Edited Sept 2019
%
% EXAMPLE
% If I ran the script as such:
%   findOverlap('/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat', 6, 9, 250, 'st_clu')
% I would be returned with an array of overlap spike times (individual
%   elements in array represent the middle of the time window with overlap)
%   eg: 255375 in the array means both clusters spiked between the
%   window of 255250 and 255500)

function [allSpikeTimes, spike1, spike2, overlap]=findOverlap(spikeClusters, cluster1, cluster2, time_window, version)
    spikeTimeRipClusStruct = load(spikeClusters);
    if (strcmp(version, 'spikeTimeRip'))
        allSpikeTimes = spikeTimeRipClusStruct.spikeTimeRipClus;
    elseif (strcmp(version, 'st_clu'))
        allSpikeTimes = spikeTimeRipClusStruct.st_clu;
    else
        disp('Invalid version name, please input "st_clu" or "spikeTimeRip"')
        return;
    end
    assignin('base', 'spikeTimeRipClus', allSpikeTimes);
    index1 = find(cell2mat(allSpikeTimes(:,1))==cluster1);
    index2 = find(cell2mat(allSpikeTimes(:,1))==cluster2);
    if (isempty(index1) || isempty(index2)) 
        disp('Cluster number specified does not exist in data, see below for valid cluster numbers')
        return
    end
    spike1 = allSpikeTimes{index1, 2};
    spike2 = allSpikeTimes{index2, 2};
    assignin('base', 'spike1', spike1);
    assignin('base', 'spike2', spike2);
    
    % find total number of time_window chunks to loop through
    max1 = max(spike1);
    max2 = max(spike2);
    if (max1 > max2) 
        largestTime = max1;
    else
        largestTime = max2;
    end
    numIterations = uint64(largestTime/time_window + 1);
    overlap = zeros(numIterations, 1);
    
    % fill in overlap array with center of time_window chunk if both
    % neurons fired within a particular time_window
    for i = 1:numIterations
        lb = time_window*(i-1);
        ub = time_window*i;
        spike1Contains = find(spike1 <= ub & spike1 >= lb);
        spike2Contains = find(spike2 <= ub & spike2 >= lb);
        if (~isempty(spike1Contains) && ~isempty(spike2Contains))
            overlap(i) = (lb+ub)/2;
        end
    end
    
    % shrink array to contain only nonzero numbers, save to workspace
    overlap = overlap(overlap ~= 0);
    assignin('base', 'overlap', overlap);
end