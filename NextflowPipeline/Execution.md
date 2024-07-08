## Execution of the pipeline in IFB cluster

To execute the pipeline in IFB cluster, three commands should be executed in the following order:

First, move to the desired folder where you want to obtain the results of the execution of the pipeline using this command:

`cd PATH_TO_YOUR_OUTPUT_FOLDER`


Then, create a screen to visualize the execution state of the pipeline using this command:

`screen -S SCREEN_NAME`


Finally, run the following command to execute the pipeline:

`sbatch PATH_TO_JOB_EXECUTION PATH_TO_PIPELINE_MAIN PATH_TO_CONFIGURATION_FILE`

With:

PATH_TO_JOB_EXECUTION: The path containing the .job pipeline execution file.

PATH_TO_PIPELINE_MAIN: The path containing the .nf pipeline execution file.

PATH_TO_CONFIGURATION_FILE: The path containing the .config pipeline configuration file.



