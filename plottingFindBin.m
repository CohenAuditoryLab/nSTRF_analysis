spikeClusters = '/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat';
version = 'st_clu';
cluster1 = 6;
cluster2 = 9;

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
    disp('Cluster number specified does not exist in data, see spikeTimeRipClus in workspace')
    return
end
spike1 = allSpikeTimes{index1, 2};
spike2 = allSpikeTimes{index2, 2};
assignin('base', 'spike1', spike1);
assignin('base', 'spike2', spike2);

optimalBinSize = findBin(spike1, spike2)