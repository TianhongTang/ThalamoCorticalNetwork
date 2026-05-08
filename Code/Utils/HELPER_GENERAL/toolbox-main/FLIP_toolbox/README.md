# (variable) Frequency-based Layer Identification Procedure (FLIP and vFLIP) for Spectrolaminar analysis

<https://doi.org/10.5061/dryad.9w0vt4bnp>

## **FLIP Algorithm File Descriptions**

### Getting Started

There are three files available for download: FLIPAnalysis.m, main_FLIP_script, and data.mat
The only file necessary to perform the FLIP algorithm on an existing non-normalized powermatrix is FLIPAnalysis.m. The provided tutorial details how to utilize FLIP. Two supplementary files are provided:
data.mat contains the data used in the paper, and main_FLIP_script.m is a script that interact with the FLIPAnalysis function.
main_FLIP_script is helpful in demonstrating how FLIP can be used in data analysis.

### File and function Descriptions

#### Published\_main\_FLIP\_script.pdf

This PDF contains the author's run of main_FLIP_script.m and associated figures produced by that script, so that users can verify they are getting the expected figure outputs when they run main_FLIP_script.m themselves.

#### data.mat

This data file contains a struct that encompasses the data recorded in experimentation.
There are five fields contained in the struct: CSD, relpow, meta, example1_vlPFC_lfp, and example2_7A_lfp

*   CSD: Current Source Density data contained in a three-dimensional matrix, the first dimension represents each probe, the second dimension represents each channel, and the third dimension represents each time point (from -0.1s pre-stimulus to +0.5s post-stimulus) included in the analysis. The units of CSD are z-score change from baseline (using the standard error across trials to perform the z-score at each point in time). Positive values denote sources and negative values denote sinks.
*   relpow: Relative Power data is contained in a three-dimensional matrix, like CSD, the first dimension represents each probe, the second dimension represents each channel, and the third dimension represents each frequency (from 1-150 Hz) included in the analysis. The units of the relpow data are normalized units between 0 and 1 (1=the channel with the highest power at that frequency).
*   example1_vIPFC_lfp: a sample LFP dataset, a three-dimensional matrix, the first dimension represents probe channels, the second dimension represents trials, and the third dimension represents time. Units of these local field potentials (LFP) are in microvolts by time (sampling rate = 1,000 Hz)
*   example2_7A_lfp: also a sample LFP dataset, a three-dimensional matrix, the first dimension represents probe channels, the second dimension represents trials, and the third dimension represents time. Units of these local field potentials (LFP) are in microvolts by time (sampling rate = 1,000 Hz)
*   meta: Metadata about each probe used in the analysis. The rows of meta correspond to the rows of CSD and relpow. Filename: The filename of the associated session, including the date of recording. probenumber: The probenumber on that session. study_number refers to study 1/2 reported in Mendoza-Halliday et al., 2023. monkey_number identifies the research subject. brain_area_num identifies the brain area numerically. interchannel_distance identifies the distance between contacts on a probe in units of microns. brain_area denotes each area's name. crossover_default is the cross identified by the FLIP algorithm with default settings. crossover_VFLIP is cross-identified by the vFLIP algorithm. CSD_sink_channel is the channel corresponding to the "early" current sink identified with a manual inspection method.

#### FLIPAnalysis.m

This file contains the FLIP and vFLIP algorithms.

#### main\_FLIP\_script.m

This file contains a few sample interactions with the FLIPAnalysis function.
To use this file, ensure the data.mat and FLIPAnalysis.m files are downloaded. Additionally, ensure it is in the same folder as FLIPAnalysis.m or add FLIPAnalysis.m to the path.

#### FLIP Algorithm Introduction.pdf

A step-by-step tutorial that reviews how to use the FLIP and vFLIP algorithms for spectrolaminar analysis.
