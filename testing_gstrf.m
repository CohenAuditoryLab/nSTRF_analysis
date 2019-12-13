struct = clusData(1, 1);
[STRFm,STRFam,STRFbm,STRFcm,x0,w,sf0,spectrop,t0,c,tf0,q,k,belta,SIs,SIt,SI,Errs,alpha_d,N]=gstrfmodel(struct.STRF1,struct.taxis,struct.faxis,struct.PP,struct.Tresh1, 'nsvd');

% [STRFm,STRFam,STRFbm,STRFcm,x0,w,sf0,spectrop,t0,c,tf0,q,k,belta,SIs,SIt,SI,Errs,alpha_d,N]=gstrfmodel(struct.STRF1,struct.taxis,struct.faxis,struct.PP,struct.Tresh,'nsvd','nsvd','n');