***Overall goal of project:***

Analyze pairwise correlations between neurons in various layers. Used Cassius_190324 and Cassius_190326 data in Shannon folder 
uploaded to PennBox by Jaejin (also contains sprfile, trig files, etc).

**Contents.m**

Monty's summary of files in his code base, helpful in locating particular functions.

**collectParams.m**

Main file of the project. This file records individual cluster parameters in struct clusData including cluster number, spike train, spike event time, STRF1A, STRF1B, STRF1s, and all parameters from strfparam including . This file also records pairwise 
analysis indices in pairedData struct including


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

