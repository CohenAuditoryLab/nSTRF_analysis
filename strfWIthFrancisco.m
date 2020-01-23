% Revisions from v5:
%   Packaged into compact function that has data, paramfile, sprfile routes 
%   as inputs and outputs STRFs along with plots of parameters
% Note that inputted data must include TrigA, TrigB, and rez_merged


%SnipTimeStamp (goes from 1-64)+SnipTimeStamp1(1-32)=>ML
%SnipTimeStamp1 (goes from 33-128)=>AL

% load data_for_Shannon.mat

function [STRF1, STRF2]= strf_analysis(data, paramfile, sprfile)

load(data);

% Data must include TrigA, TrigB, and rez_merged which should be variables
% in workspace after loading data in
if (exist('TrigA', 'var') ~= 1)
    error('Data must include TrigA')
elseif (exist('TrigB', 'var') ~= 1)
    error('Data must include TrigB')
elseif (exist('rez_merged', 'var')~= 1)
    error('Data must include rez_merged')
end

fs=24414.0625;
block=2;
% paramfile= '/Users/shannon1/Documents/S19/Cohen_research/STRF/Moving_Ripple_generation/DNR_Cortex_96k5min_param.mat';

result = load(paramfile);

warning('off','MATLAB:interp1:UsePCHIP')
% sprfile='/Users/shannon1/Documents/S19/Cohen_research/STRF/Moving_Ripple_generation/DNR_Cortex_96k5min.spr';
s_bin=0.15;
cs=1;

% myKsDir='C:\work\ToSort_Sam\SAM-190503_b2_xpz5AL';
% spikes = loadKSdir(myKsDir);
spikes=rez_merged.st3;
index_clusters=unique(spikes(:,5));
conta=1;

% initialize array of params
% hard-coded 379 to match # with images saved 
% MATLAB indexing starts at 1 and not 0 so params array is off by 1 (index
% 1 corresponds with picture 0)
bestFm_params = zeros(379, 1);
bTMF_params = zeros(379, 1);
cTMF_params = zeros(379, 1);
bestRD_params = zeros(379, 1);
bSMF_params = zeros(379, 1);
cSMF_params = zeros(379, 1);
PLI_params = zeros(379, 1);
DSI_params = zeros(379, 1);
BW_params = zeros(379, 1);
Delay_params = zeros(379, 1);
Duration_params = zeros(379, 1);
FmUpperCutoff_params = zeros(379, 1);
RDUpperCutoff_params = zeros(379, 1);

% initialize array of PeakBF values
PeakBF_params = zeros(379, 1);
save_strf = '';
for s=1:length(index_clusters)
    index_temp=find(spikes(:,5)==index_clusters(s));
    spet=spikes(index_temp,1);
    cluster_number=index_clusters(s);
        try
            [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint(sprfile,0,s_bin,spet',TrigA,fs,80,30,'dB','MR',1700,5,'float');
        catch me
        %disp(me);
            STRF1A=[];STRF2A=[];No1A=[];Wo1A=[];No1A=[];Wo1A=[];
        end
    
        try
            [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(sprfile,0,s_bin,spet',TrigB,fs,80,30,'dB','MR',1700,5,'float');
        catch me
            %disp(me);
            STRF1B=[];STRF2B=[];No1B=[];Wo1B=[];No1B=[];Wo1B=[];
        end
    
        if ~isempty(STRF1A) && ~isempty(STRF1B)
        STRF1 = (STRF1A+STRF1B)/2;
        STRF2 = (STRF2A+STRF2B)/2;
        No1   =  No1A+No1B;
        Wo1   =  (Wo1A+Wo1B)/2;
        No2   =  No2A+No2B;
        Wo2   =  (Wo2A+Wo2B)/2;
        
        [STRF1s,Tresh1]=wstrfstat(STRF1,0.001,No1,Wo1,PP,30,'dB','MR','dB');
        [STRF2s,Tresh2]=wstrfstat(STRF2,0.001,No2,Wo2,PP,30,'dB','MR','dB');
        %Getting STRF param
        [RF1P]=strfparam(taxis,faxis,STRF1s,Wo1,PP,'MR',500,4,0.5,0.05,0.1,'n');
        
        % Update collection of .PeakBF values
        PeakBF_params(cluster_number+1, 1) = RF1P.PeakBF;
        
        % Update collection of .BestFm, bTMF, cTMF values
        % BestFm returns 1 or 2 values, fill in matrix accordingly
        [rowFm, colFm] = size(RF1P.BestFm);
        if (colFm == 2) 
            bestFm_params(cluster_number+1, 1) = (abs(RF1P.BestFm(1,1)) + abs(RF1P.BestFm(1,2)))/2;
        else
            bestFm_params(cluster_number+1, 1) = RF1P.BestFm;
        end     
        bTMF_params(cluster_number+1, 1) = RF1P.bTMF;
        cTMF_params(cluster_number+1, 1) = RF1P.cTMF;
        
        % Update collection of various params
        % BestRD returns 1 or 2 values, fill in matrix accordingly
        [rowRD, colRD] = size(RF1P.BestRD);
        if (colRD == 2)
            bestRD_params(cluster_number+1, 1) = (abs(RF1P.BestRD(1,1)) + abs(RF1P.BestRD(1,2)))/2;
        else
            bestRD_params(cluster_number+1, 1) = RF1P.BestRD;
        end     
        bSMF_params(cluster_number+1, 1) = RF1P.bSMF;
        cSMF_params(cluster_number+1, 1) = RF1P.cSMF;
        
        % Update collection of .PLI and .DSI values
        PLI_params(cluster_number+1, 1) = RF1P.PLI;
        DSI_params(cluster_number+1, 1) = RF1P.DSI;
        
        BW_params(cluster_number+1, 1) = RF1P.BW;
        Delay_params(cluster_number+1, 1) = RF1P.Delay;
        Duration_params(cluster_number+1, 1) = RF1P.Duration;
        FmUpperCutoff_params(cluster_number+1, 1) = RF1P.FmUpperCutoff;
        RDUpperCutoff_params(cluster_number+1, 1) = RF1P.RDUpperCutoff;

        [RDHist1,FMHist1,RDHist2,FMHist2,Time1,Time2] =rtfhist(paramfile,spet',TrigA,fs);
        [Xa,Ya,Na]=hist2(FMHist1,RDHist2,20,20,'n');
        [RDHist1b,FMHist1b,RDHist2b,FMHist2b,Time1b,Time2b]=rtfhist(paramfile,spet',TrigB,fs);
        [Xb,Yb,Nb]=hist2(FMHist1,RDHist2,20,20,'n');
    
        elseif isempty(STRF1A) && ~isempty(STRF1B)
        STRF1 = STRF1B;
        STRF2 = STRF2B;
        No1   =  No1B;
        Wo1   =  Wo1B;
        No2   =  No2B;
        Wo2   =  Wo2B;
        elseif isempty(STRF1B) && ~isempty(STRF1A)
        STRF1 = STRF1A;
        STRF2 = STRF2A;
        No1   =  No1A;
        Wo1   =  Wo1A;
        No2   =  No2A;
        Wo2   =  Wo2A;
        elseif isempty(STRF1B) && isempty(STRF1A)
        STRF1 = [];
        STRF2 = [];
        No1   =  [];
        Wo1   =  [];
        No2   =  [];
        Wo2   =  [];
        taxis =  [];
        faxis =  [];
        PP=[];
        SPLN=[];
        end
    
    STRFDataxpz2_AL(conta) = struct('cluster',cluster_number,'No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN,'Xa',Xa,'Ya',Ya,'Na',Na,'Xb',Xb,'Yb',Yb,'Nb',Nb);
    figure;
    taxis = (taxis)*1e3;
    faxis = (faxis)*1e3;

    if ~isempty(STRF1)
        subplot(1,4,1)
        pcolor(taxis,log2(faxis/faxis(1)),STRF1);
        
        colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
        subplot(1,4,2)
        pcolor(RF1P.Fm,RF1P.RD,RF1P.RTF);
        axis([-50 50 0 3.9]);colormap jet;
        subplot(1,4,3)
        hist2(FMHist1,RDHist1,40,40,'y'); 
        [x,y,n] = hist2(FMHist1,RDHist1,40,40,'y');
        subplot(1,4,4) 
        contour(x, y, n);

        % create new folder named strf_analysis, cd into it, and save path
        % to save STRF contours into the folder
        if ~exist('strf_analysis','dir')
            mkdir('strf_analysis')
            addpath('strf_analysis','-end') 
            cd strf_analysis
            save_strf = pwd;
        end
        print(['STRF_contour' num2str(cluster_number)],'-djpeg');
        close all
        conta=conta+1;
     end
end
cd(save_strf)
save(save_strf,'STRFDataxpz2_AL');

subplot(1,4,1)
[PLI_2, PLI_4, PLI_6, PLI_8] = average_param(PeakBF_params, PLI_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [PLI_2 PLI_4 PLI_6 PLI_8];
bar(x,y)
title('PLI vs PeakBF');

subplot(1,4,2)
[DSI_2, DSI_4, DSI_6, DSI_8] = average_param(PeakBF_params, DSI_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [DSI_2 DSI_4 DSI_6 DSI_8];
bar(x,y)
title('DSI vs PeakBF');

% plot various parameters vs. PeakBF
subplot(1,4,3)
[bTMF_2, bTMF_4, bTMF_6, bTMF_8] = average_param(PeakBF_params, bTMF_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [bTMF_2 bTMF_4 bTMF_6 bTMF_8];
bar(x,y)
title('bTMF vs PeakBF');

subplot(1,4,4)
[bSMF_2, bSMF_4, bSMF_6, bSMF_8] = average_param(PeakBF_params, bSMF_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [bSMF_2 bSMF_4 bSMF_6 bSMF_8];
bar(x,y)
title('bSMF vs PeakBF');

figure();

subplot(1,4,1)
[BW_2, BW_4, BW_6, BW_8] = average_param(PeakBF_params, BW_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [BW_2 BW_4 BW_6 BW_8];
bar(x,y)
title('BW vs PeakBF');

subplot(1,4,2)
[Delay_2, Delay_4, Delay_6, Delay_8] = average_param(PeakBF_params, Delay_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [Delay_2 Delay_4 Delay_6 Delay_8];
bar(x,y)
title('Delay vs PeakBF');

subplot(1,4,3)
[Duration_2, Duration_4, Duration_6, Duration_8] = average_param(PeakBF_params, Duration_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [Duration_2 Duration_4 Duration_6 Duration_8];
bar(x,y)
title('Duration vs PeakBF');

subplot(1,4,4)
[FmUpperCutoff_2, FmUpperCutoff_4, FmUpperCutoff_6, FmUpperCutoff_8] = average_param(PeakBF_params, FmUpperCutoff_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [FmUpperCutoff_2 FmUpperCutoff_4 FmUpperCutoff_6 FmUpperCutoff_8];
bar(x,y)
title('FmUpperCutoff vs PeakBF');

figure();

subplot(1,4,1)
[RDUpperCutoff_2, RDUpperCutoff_4, RDUpperCutoff_6, RDUpperCutoff_8] = average_param(PeakBF_params, RDUpperCutoff_params);
x = categorical({'0-2','2-4','4-6', '6-8'});
y = [RDUpperCutoff_2 RDUpperCutoff_4 RDUpperCutoff_6 RDUpperCutoff_8];
bar(x,y)
title('RDUpperCutoff vs PeakBF');

subplot(1,4,2)
scatter(PLI_params, DSI_params);
title('DSI vs PLI');

subplot(1,4,3)
scatter(DSI_params, PLI_params);
title('PLI vs DSI');
