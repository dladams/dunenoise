## Protodune in dunenoise

David Adams<br>
July 2020

### Introduction

Dunenoise is one in series of DUNE analysis packages.
It depends on [duneproc](https://github.com/dladams/duneproc) which must also be installed.
The examples below make use of the duneproc command from that package.
See the protoDUNE tutorials e.g. [tutorial 2](https://github.com/dladams/duneproc/blob/master/doc/tutorial02.md) in that package for some guidance on package installation, use of the duneproc command
and construction of dataset descriptions used with that command.

### Noise distributions

Code for this has yet to be added. Examples may be found in the coce used to generate
the noise plots for the protoDUNE-SP performance paper. See ~dladams/dev/dudev03/np04-prod/jul19.

### DFT power plots

The code for the performance note is partly transferred. For more, see ~dladams/dev/dudev02/np04-proc/sep2019/dft1000.

The current code (2020-07-20) is that used to produce the DFT spectra for the July 2020 charge injection
runs. For results, see my [web area](https://internal.dunescience.org/people/dladams/protodune/studies/2020-07-20_wire_plane_injection).

Top-level fcl is also provided to create DFT (discrete Fourier transform) distributions of noise power vs. frequency, e.g.
<pre>
duneproc pdspTaiNoiseDft/pdspSetDftFmax100 np04_run011780_evts000000-001000
</pre>
produces the raw DFT plots and then the Root script pdmergeDftPlots may be used to pretty up the result:
<pre>
root> pdMergeDftPlots("./pdspTaiNoiseDft/pdspSetDftFmax100/np04_run011780_evts000000-00100_proc000005/dftpowtlog_run011780.tpad", "5 kHz on Z:", "", 2e-4)
</pre>

To make these plots, a signal finder is applied and inverted and then regions of 2000 contiguous samples are transformed
to obtain the final power distributions which are averaged over all such regions in each wire plane.
Bad and noisy channels are excluded.
