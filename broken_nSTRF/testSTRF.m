% Change these paths to run script

sprfile = '/Users/shannon1/Documents/F19/neuroResearch/broken_nSTRF/DNR_Cortex_96k5min_4_50.spr';
spike_clusters = '/Users/shannon1/Documents/F19/neuroResearch/broken_nSTRF/spike_times_ripple_clust_old.mat';
% spike_clusters = '/Users/shannon1/Documents/F19/neuroResearch/broken_nSTRF/spike_times_ripple_clust_new.mat';
Trig = '/Users/shannon1/Documents/F19/neuroResearch/broken_nSTRF/AudiResp_16_24-190326-154559_triggers.mat';
trigStruct = load(Trig);
TrigA = trigStruct.TrigA;
TrigB = trigStruct.TrigB;

cluster1 = 1;
cluster2 = 2;
time_window=500;

[spikeTimeRipClus, spike1, spike2, overlap] = findOverlap(spike_clusters, cluster1, cluster2, time_window, 'old');    
% [spikeTimeRipClus, spike1, spike2, overlap] = findOverlap(spike_clusters, cluster1, cluster2, time_window, 'new'); 
if (isempty(overlap))
    disp('No overlap spike time between spike1 and spike2 so nSTRF cannot be calculated, exiting function')
    return;
end
s_bin=0.15;
T1=0;
T2=s_bin;
Fss=24414.0625;
overlap = double(overlap)/1000.0; % convert spike times ms to sec
spetOverlap=overlap * Fss;
SPL=80;
MdB=30;
ModType='dB';
Sound='MR';
NBlocks=100;
UF=10;
sprtype='float';
try
%     [nTaxis,nFaxis,nSTRF1A,nSTRF2A,nPP,nWo1A,nWo2A,nNo1A,nNo2A,nSPLN]=rtwstrfdbint(sprfile,T1,T2,spetOverlap',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
%     assignin('base', 'nSTRF1A', nSTRF1A);
    [nTaxis,nFaxis,nSTRF1B,nSTRF2B,nPP,nWo1B,nWo2B,nNo1B,n2B,nSPLN]=rtwstrfdbint(sprfile,T1,T2,spetOverlap',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype); 
    assignin('base', 'nSTRF1B', nSTRF1B);
catch me
    disp('Error when computing nSTRF, exiting function')
    disp(me);
    return;
end