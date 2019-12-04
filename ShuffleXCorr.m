

%Loading SPET and Trigger Files
load('/Users/shannon1/Documents/F19/neuroResearch/nSTRF/AudiResp_16_24-190326-154559_triggers.mat')
load('/Users/shannon1/Documents/F19/neuroResearch/nSTRF/spike_times_ripple_clust_new.mat')

%Batch to compute xcorr between recording locations
%Need to add shuffling procedure 
Fs= 24414.0625  %Original sampling rate - ideally use the sampling rate for your data - but the way its done here it wont matter
Fsd=250         % The desired sampling rate for correlation analysis - the spike train is resampled at Fsd prior to computing the correlation
                % you can change to 250 for 4 ms resolution
MaxTau=.50      %Correlation lag - up to 250 ms
T=300           %Experiment duration - Im assuming this is the correct value
Zero='n'        % Do not remove correlation value at zero lag - only useful for autocorrelograms
Mean='n'        %Do not remove mean correaltion value - DC
Disp='y'

%Finding spike times for different trials
K=1 %Select first cluster
L=2 %Select second cluster
spet1=round(cell2mat(st_clu(K,2))*Fs)';             %Spike even time vector for first unit (in sample number not real time)
spet2=round(cell2mat(st_clu(L,2))*Fs)';             %Spike event time vector for second unit
[spet1A,spet1B]=spet2spetab(spet1,TrigA,TrigB,Fs);
[spet2A,spet2B]=spet2spetab(spet2,TrigA,TrigB,Fs);
[spet1As]=shufflespet(spet1A);
[spet1Bs]=shufflespet(spet1B);
[spet2As]=shufflespet(spet2A);
[spet2Bs]=shufflespet(spet2B);


%Trial Shuffled spike train cross covariance
Disp='n';
[R1A2B]=xcovspike(spet1A,spet2B,Fs,Fsd,MaxTau,Disp);
[R1B2A]=xcovspike(spet1A,spet2B,Fs,Fsd,MaxTau,Disp);
R12shuf=(R1A2B+R1B2A)/2;

%Spike train cross covariance
[R1A2A]=xcovspike(spet1A,spet2A,Fs,Fsd,MaxTau,Disp);
[R1B2B]=xcovspike(spet1B,spet2B,Fs,Fsd,MaxTau,Disp);
R12=(R1A2A+R1B2B)/2;

%Random Control - Trial shuffled Spike train cross covariance
[R1A2Bs]=xcovspike(spet1As,spet2Bs,Fs,Fsd,MaxTau,Disp);
[R1B2As]=xcovspike(spet1As,spet2Bs,Fs,Fsd,MaxTau,Disp);
R12shuf_s=(R1A2B+R1B2A)/2;

%Random Control - Spike train cross covariance
[R1A2As]=xcovspike(spet1As,spet2As,Fs,Fsd,MaxTau,Disp);
[R1B2Bs]=xcovspike(spet1Bs,spet2Bs,Fs,Fsd,MaxTau,Disp);
R12s=(R1A2As+R1B2Bs)/2;


%%%%%%%%%%%%%%%%%%%%%%%% COMPUTING STRFS %%%%%%%%%%%%%%%%%%%%%%

%Generating Receptive Fields
SpecFile='DNR_Cortex_96k5min_4_50.spr'
T1=0;
T2=.25; % changing to 0.15 gives me trouble, need it for getPaired.m tho
Fss=Fs
SPL=70
MdB=30
ModType='dB'
SModType='dB'
Sound='MR';
NBlocks=400
sprtype='float';
p=0.001

% T1=0;
% T2=0.15;
% Fss=24414.0625;
% SPL=80;
% MdB=30;
% ModType='dB';
% SModType='dB';
% Sound='MR';
% NBlocks=100;
% sprtype='float';
% p=0.05;

%Computing STRF for first cluster
[taxis,faxis,STRF1A,SSS,PP,Wo1A,Wo2,No1A,No2,SPLN]=rtwstrfdb(SpecFile,T1,T2,spet1,TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,sprtype);
[taxis,faxis,STRF1B,SSS,PP,Wo1B,Wo2,No1B,No2,SPLN]=rtwstrfdb(SpecFile,T1,T2,spet1,TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,sprtype);
STRF1=(STRF1A+STRF1B)/2;
[STRF1s,Tresh]=wstrfstat(STRF1,p,No1A+No1B,(Wo1A+Wo1B)/2,PP,MdB,ModType,Sound,SModType);

%Computing STRF for second cluster
[taxis,faxis,STRF2A,SSS,PP,Wo2A,Wo2,No2A,No2,SPLN]=rtwstrfdb(SpecFile,T1,T2,spet2,TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,sprtype);
[taxis,faxis,STRF2B,SSS,PP,Wo2B,Wo2,No2B,No2,SPLN]=rtwstrfdb(SpecFile,T1,T2,spet2,TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,sprtype);
STRF1=(STRF2A+STRF2B)/2;
[STRF2s,Tresh]=wstrfstat(STRF1,p,No2A+No2B,(Wo2A+Wo2B)/2,PP,MdB,ModType,Sound,SModType);

% STRF1s = pairedData.nSTRFOne.STRFs;
% STRF2s = pairedData.nSTRFTwo.STRFs;
%Predict Correlation From the STRFs - Chen 2012
[T,R,Rcc,RR]=strf2xcorr(taxis,faxis,STRF1s,STRF2s,PP,'y');


%%%%%%%%%%%%%%% PLOTTING RESULTS %%%%%%%%%%%%%%%%

%Receptive field for cluster 1
figure
subplot(331)
Max=max(max(abs(STRF1s)));
imagesc(taxis*1000,log2(faxis/faxis(1)),STRF1s),caxis([-Max Max])
set(gca,'YDir','normal')
xlabel('Delay (ms)')
ylabel('Freq. (Oct)')
title(['STRF clust ' num2str(cell2mat(st_clu(K,1)))])

%Receptive field for cluster 2
subplot(334)
Max=max(max(abs(STRF2s)));
imagesc(taxis*1000,log2(faxis/faxis(1)),STRF2s),caxis([-Max Max])
set(gca,'YDir','normal')
xlabel('Delay (ms)')
ylabel('Freq. (Oct)')
title(['STRF clust ' num2str(cell2mat(st_clu(L,1)))])

%Receptive field covariance
[RSTRF] = strfcorr(STRF1s,STRF2s,taxis,faxis,PP,100,4);
subplot(335)
Max=max(max(abs(RR)));
imagesc(T*1000,log2(faxis/faxis(1)),RR),caxis([-Max Max])
set(gca,'YDir','normal')
Pos=get(gca,'Position');
set(gca,'Position',[Pos(1) Pos(2)+.15 Pos(3) Pos(4)])
xlabel('Lag (ms)')
ylabel('Freq. (Oct)')
title('STRF cross-corr')

%Predicted spike train covariance
subplot(336)
Max=max(max(abs(RR)));
plot(T*1000,Rcc/30,'r') %Rescaled amplitude - note that as per Chen et al 2012, the predicted value is an upper bound. The shape prediction is actually very good. 
xlim([-MaxTau MaxTau])
Pos=get(gca,'Position');
set(gca,'Position',[Pos(1) Pos(2)+.15 Pos(3) Pos(4)])
title('Predicted corr (red) / orig. corr (black)')

%Plotting Spike Train Covariance
hold on
N=(length(R12)-1)/2;
Tau=(-N:N)/Fsd*1000;
plot(Tau,R12,'k')
xlim([-250 250])





