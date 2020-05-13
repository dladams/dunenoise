## Iceberg in dunenoise

David Adams<br>
May 2020

### Introduction

Dunenoise is one in series of DUNE analysis packages.
It depends on [duneproc](https://github.com/dladams/duneproc) which must also be installed.
The examples below make use of the duneproc command from that package.
See the Iceberg tutorial in that package for some guidance on package installation, use of the duneproc command
and construction of dataset descriptions used with that command.

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

### DFT power plots

Top-level fcl is also provided to create DFT (discrete Fourier transform) distributions of noise power vs. frequency, e.g.
<pre>
duneproc ibtaiNoiseDft iceberg005044
duneproc ibcnrNoiseDft iceberg005044
</pre>
to repectively create the plots before and after CNR.
The run directory ibtaiNoiseDft/iceberg005044 will hold the plot file [dftpowt_run005044.png](dftpowt_run005044.png) with a separate DFT power distribution
for each of the four wire planes (z1, z2, u and v).
The corresponding tpad file provides acccess to the underlying histograms.

A signal finder is applied and inverted and then regions of 1000 contiguous samples are transformed to obtain the
final power distributions which are averaged over all such regions in each wire plane.
The distribution is binned---100 bins for 500 frequencies---and the limits chosen so the zero frequency term is recorded as underflow
and appears as a solid bar in the plots.

The script [ibMergeDftPlots](../Script/ibMergeDftPlots) overlays the before and after CNR plots, e.g.
<pre>
ibMergeDftPlots 5044
</pre>
This produces the plot file [dft5044.png](dft5044.png).

The zero frequncy contributions are not shown on these plots butit is possible to make unbinned plots of the first 100 frequencies (i.e. below 200 kHz) that includes this term:
<pre>
duneproc ibtaiNoiseDft/ibsetDftFmax200 iceberg005044
duneproc ibcnrNoiseDft/ibsetDftFmax200 iceberg005044
ibMergeDftPlots 5044 200
</pre>
with output now in [dft5044_fmax200.png](dft5044_fmax200.png).
