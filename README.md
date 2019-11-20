**Overall goal of project:**
Analyze pairwise correlations between neurons in various layers. Used Cassius_190324 and Cassius_190326 data in Shannon folder 
uploaded to PennBox by Jaejin.

**pairedSTRFanalysis.m:**
This file computes the significant STRFs of spike times from 2 clusters, the OR/AND STRF btw the two STRFs,
coincident spike times at calculated optimal bin size, and 
cross correlation between coincident spikes and OR/AND STRF. Coincident spike times between spike1 and spike2 is computed by
recording the original spike1/2 times in periods of overlap within some bin size. The optimal bin size was computed by
finding the time lag at which cross correlation was at its peak.

