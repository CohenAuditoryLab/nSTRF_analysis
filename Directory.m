% Shannon's summary of function descriptions from Monty's Keck toolbox
% (inputs in parenthesi are optional input args)
%
% gstrfmodel: fitting STRF by Gabor functions
%   inputs: STRF,taxis,faxis,PP,Tresh,TreshType,Method,display
%   outputs: STRFm,STRFam,STRFbm,STRFcm,x0,w,sf0,spectrop,t0,c,tf0,q,k,belta,SIs,SIt,SI,Errs,alpha_d,N
%   
% strfcorr: computes correlation between two STRFs to determine, the level of similarity. 
%   inputs: STRF1,STRF2,taxis,faxis,PP (MaxDelay,MaxFreqShift)
%   outputs: R, tau, dX, optimal delay, frequency shift, SI, SI00, SIt, SIf
%
% strfcorrcorrected: like strfcorr except takes into account internal noise
%   inputs: STRF1A,STRF1B,STRF1s,STRF2A,STRF2B,STRF2s,taxis,faxis,PP
%       (NoiseFlag,MaxDelay,MaxFreqShift)
%   outputs: R, tau, dX, optimal delay, frequency shift, SI, SI00, SIt, SIf
% 
% batchdbsplind: Computes SI for the contrast-intensity response curves
%   inputs: BatchFile,Header,T,p
%   outputs: SI,SInoLin,FileList,SpikeType
