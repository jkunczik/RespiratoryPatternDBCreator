# RespiratoryPatternDBCreator

This repository contains the code for the following paper:

> Kunczik, J.; Hubbermann, K.; Mösch, L.; Follmann, A.; Czaplik, M.; Barbosa Pereira, C. Breathing Pattern Monitoring by Using Remote Sensors. Sensors 2022, 22, 8854. https://doi.org/10.3390/s22228854

The code can be used to repdroduce the dataset used in the publication. If the code was useful for your research, please consider citing the original paper. 

## Code description

| module               | description                                                                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| dataset_augmentation | augments the original dataset by randomly combining various original signals with each other and then scrambling the resulting signal's frequency |
| feature_extraction   | code to extract features, useful for breathing pattern classification                                                                             |
| labeling             | small GUI to label the signals                                                                                                                    |
| signal_alignment     | aligns all recorded signals with each other                                                                                                       |
| wfdb_export          | convert data to wfdb for publication on PhysioNet                                                                                                                                                  |

The `main`script combines all  presented processing steps and algorithms. 

## Installation
The code was written and tested for MATLAB R2022a. To run the code, download or clone this repository and add its root folder, as well as all subfolders to the search path. Run `main.m` to run the complete processing pipeline.
