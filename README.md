# dream-beat

Model Presentation : [click here](https://www.synapse.org/#!Synapse:syn21827452/wiki/602762)
Dream BEAT-AML challenge.

The goal of the Beat AML DREAM Challenge is to define patient subpopulations tailored to individual treatments by discovering (genomic and transcriptomic) biomarkers of drug sensitivity, as evaluated on an unpublished cohort of patients from the BeatAML project.

We developed a transfer learning and ensemble learning approach where we learnt characteristics from published cohorts and trained the new data and new features available from the unpublished cohort for prognostic modelling . 

Our results were in the top 5 but most importantly the framework developed could be applied for any dataset. 


 [Methods and Workflow](Methods_and_Workflow_SubChallenges.pdf) for explainations.


Challenge 1: Ex-vivo drug sensitivity prediction with kernel multi-task regression.

Challenge 2: Outcome prediction from knowledge transfer and federated learning of large scale AML genomic study.



Repo structure:


- **data/**: raw data and processed dataframes (features and responses) for predictions.


- **src/**: generic source functions as *predictors* and *evaluation* for subchallenges 1 and 2 separately.


- **prediction/**: scripts to run predictions. Please create subfolders if necessary with the dates YYMMDD.


- **analysis/**: folder for exploratory analysis. Please create subfolders if necessary with the dates YYMMDD.


- **docker/**: ....
