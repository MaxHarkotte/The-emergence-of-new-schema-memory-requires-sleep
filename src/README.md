# The-emergence-of-new-schema-memory-requires-sleep
All video files analyzed with the following scripts can be found at: osf.io/u9h7d/
Tracking and behavioral scoring were done in AnyMaze. To view Videos, please download AnyMaze from: any-maze.com

Processed data can be found under ~/dat

- '00-EncodingClean.csv' contains all tracking and behavioral scoring data from the encoding episodes
- '00-Sleep.csv' contains the sleep scoring off all animals that were sleeping during the post-encoding interval
- '00-TestClean.csv' contains all tracking and behavioral scoring from the retrieval phase 
- '00-TestClean_Rearing.csv' contains behavioral scoring of rearing behavior from the retrieval phase of the main experiment
- 'cfos.csv' contains positive c-Fos cell counts from each subregion

All analyses were conducted with R version 4.4.1 (2024-06-14 ucrt). 

- '01_Behavioral_Analyses.R' was used to run statistics and plot on all behavioral tracking data 
- '02_cFos_Analysis.R' was used to run statistics and correlational network analyses on cFos-positive cells 

