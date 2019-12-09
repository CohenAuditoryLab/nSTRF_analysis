***Overall goal of project:***

Analyze pairwise correlations between neurons in various layers. Used Cassius_190324 and Cassius_190326 data in Shannon folder 
uploaded to PennBox by Jaejin (also contains sprfile, trig files, etc).

**Main 3 final files:**
getIndividual.m: computes and saves various pairwise parameters
getPaired.m: computes and saves various individual cluster parameters
mommaScript.m: script to run the above two scripts

**Contents.m**

Monty's summary of files in his code base, helpful in locating particular functions.

**getIndividual.m**

This file records individual cluster parameters in struct clusData. Can be easily modified to hold more parameters if of interest.

clusData parameters include: 
cluster number, spike train, spike event time, STRF1A, STRF1B, STRF1, STRF1s, STRF1sBinary, and all parameters from Monty's strfparam including delay, duration, best frequency (in octaves and Hz), spectral bandwidth (octaves and Hz), delay at STRF peak, best frequency at peak, delay measurement at peak of temporal envelope, best frequency at peak of spectral envelope, various envelope duration measurements, best modulation rate, best ripple density, best temporal modulation frequency, best spectral modulation freq, temporal modulation freq centroid, spectral modulation freq centroid, spectral MTF bandwidth, temporal MTF bandwidth, temporal modulation freq upper/lower cutoff, spectral modulation freq upper/lower cutoff, direction selectivity inde, peak response from ripple density plot, temporal envelope, spectral envelope, phase locking index.

**getPaired.m**
This file records pairwise analysis indices in pairedData struct.

pairedData parameters include: 
clusOne (cluster one number), ClusTwo, optimalBinSize, nSTRFOne (STRF of first cluster coincident spike time), nSTRFTwo, orSTRF, andSTRF, coin1CorrOrZeroLag, coin1CorrAndZeroLag, coin2CorrOrZeroLag, coin2CorrAndZeroLag, coin1CorrCoin2ZeroLag, montyCoin1CorrCoin2ZeroLag (xcorr between coincident first spike and coincident second spike using Monty's calculations), MontyRSTRF (struct containing various other indices), T, R, Rcc, RR (last four from Monty's ShuffleXCorr.m)

**pairedSTRFanalysis.m:**

This file computes the significant STRFs of spike times from 2 clusters, the OR/AND STRF btw the two STRFs,
coincident spike times at calculated optimal bin size, and cross correlation between coincident spikes and OR/AND STRF.

**findOverlap2.m:**

This file computes the coincident spike times. Coincident spike times between spike1 and spike2 is computed by
recording the original spike1/2 times in periods of overlap within some bin size. The optimal bin size is computed using findBin.m

**findBin.m:**

This file computes optimal bin size to use when recording coincident spike times. Optimal bin size is the time lag when cross correlation spike between inputted spet1 and spet2 is at its peak, calculated using Monty's xcorrspike.m

**binarizeSTRFs.m:**

This file computes the siginificant STRF using Monty's wstrfstat.m and then converts the outputted STRFs to 0/1 where 1 is 
significant, 0 else. Significance is determined as any pixel greater than threshold, which I defined as 0.05/total number of pixels. The result of this file was used to calculate OR/AND STRF between two STRFs using logical or/and between respective pixels.

**Monty's code**

rtwstrfdbint.m, strfcorrcorrected.m	, strfparam.m, wstrfstat.m, xcorrspike.m

