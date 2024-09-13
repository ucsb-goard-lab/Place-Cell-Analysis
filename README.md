# Place Cell Analysis Pipeline
Pipeline to analyze Place Cells: including data extraction, preprocessing, and analysis of place cells.

## The pipeline consists of three main steps:
#### Data Extraction:
The pipeline starts by adding all necessary functions to the path. It then extracts data from the suite2p output and converts it into a format that is compatible with the pipeline. The data is then deconcatenated based on the number of environments.
#### Preprocessing:
The pipeline performs local neuropil subtraction to remove the background fluorescence from the raw fluorescence signal. It also detects times when the mouse is moving and removes times when the mouse is not moving from the data.
#### Place Cell Analysis:
The pipeline analyzes the place cells in the data. It identifies active cells, computes spatial information, and bins and smooths the activity data. It also identifies valid place cells and speed cells.

## Usage
To run the pipeline, open your data folder with the output from suite2p and the neurotar, run the master script RemappingAnalysisMaster.m, and the program will output a bunch of data, including all the DFF information, all the place cell dynamics, and remapping info. 

For instance, the input data from one recording will have 3 neurotar files from the 3 different environments, and one big suite2p file with all the data from each of the three recordings concatenated together.

## Contributors
Nora Wolcott: `HPC_Analysis_Pipeline_Method3`, `DeConcatenateEnvironments_v2`, `PC_reliability_checker_WTR_v2`, `Moving_v3`, `Active_Cells`, `getLaps`, `suite2p2data`, `Save_Data`, `plotOneD`, `plotDonutHeatmap`,`plotAllLapByLap`,`getExtractorInpt`, 
       
Will Redman: `Spike_Max`, `Speed_Cells`, `SpatialInfoComputer`, `PC_reliability_checker_WTR_v2`, `OneD_Track_Anaysis_v2`, `Normalizer`, `Moving_v3`, `Active_Cells`, `DFF_transients`, `suite2p2data`,

Santi Acosta: `NewNeurotarExtractor`,

Kevin Sit: `colormapMaker`,

James Roney: `subroutine_test_r_HPC`, `subroutine_find_corr_HPC`, 

Pengcheng Zhou and ported from the Python implementation from Johannes Friedrich: `foopsi_oasisAR1`, `estimate_time_constant`, `deconvolveCa`, `GetSn`, `oasisAR1`,

Eftychios Pnevmatikakis: `foopsi_oasisAR1`, `deconvolveCa`, `oasisAR1`

## References
- Hainmueller T and Bartos M (2018). Parallel emergence of stable and dynamic memory engrams in the hippocampus. Nature 558, 292–296.
- Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging.
- Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data.
