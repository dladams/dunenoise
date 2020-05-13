## Iceberg in dunenoise

David Adams<br>
May 2020

### Introduction

Dunenoise is one in series of DUNE analysis packages

### Noise distributions

Top-level fcl is provided to evaluate sample and (50 sample) integrated noise for each channel.
To evaluate the noise before and after CNR (coherent noise removal):
<pre>
duneproc ibtaiNoiseVsChan iceberg005044
duneproc ibcnrNoiseVsChan iceberg005044
</pre>
The ouput of the first job is in the run directory ibtaiNoiseVsChan/iceberg005044 and includes plots of sample noise at
<pre>
chmet_nsgrms_all_run005044.png
</pre>
and integrated noise at
<pre>
chmet_nsgrms50_all_run005044.png
</pre>
The underlying histograms may be found in the corresponding tpad files (replace png with tpad in the file names).

Plots showing the distributions of sample and integrated noise before and after CNR can be created from the those histograms
using the [drawNoise](../Script/drawNoise) script:
<pre>
drawNoise 5044
</pre>
An example plot is [here](noise_tai-tai-50-cnr-cnr-50_zcGood-uvGood_run005044.png).
