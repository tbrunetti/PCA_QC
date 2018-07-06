# PCA_QC
dynamic correlation analysis of PCs with metadata for QC metrics

A way to look at whether the first 10 or n PCs in your dataset correlate with any measured metadata metrics

Can be used with almost any continuous data matrix that has corresponding metadata.  Primarily built for use with RNA-seq data.  Run can be performed on a local standard computer; no HPC required.  If running on HPC, please use in interactive mode.

## Required R libraries  
* factoextra
* ggfortify
* ggplot2
* data.table
* gplots

## Output Files  
3 output files are generated:  
1. output_name_pca_matrix_output.csv  --provided so you can easily test your regression by looking at residuals
2. output_name_metadata_sorted_output.csv  --provided so you can easily test your regression by looking at residuals
3. output_name.pdf 


## Required Inputs
* counts in raw counts,  FPKM, or TPM units for all genes/isoform used in experiment for every sample (as few or as many as one would like to consider in analysis)
* metadata file where metrics are listed for every sample in experiment (as few or as many as one would like to consider in analysis) 

The counts data should be formatted as follows in **CSV format**:  
* __1st Line__: sample name or ID  
* __1st Column__: name of gene or isoform  
* all cells should be populated with a non-negative integer or float     

![alt text](https://github.com/tbrunetti/PCA_QC/blob/develop/counts_file_format.png "counts file format")  

The metadata file should be formatted as follows in **CSV format**:  
* __1st Line__: C or N, indicating if the column variable should be considerd C for categorical or N for continuous/numeric  
* __2nd Line__: the metadata variable name (can have an few or as many as one would like to consider in analysis)  
* __1st Column__: the sample name or ID, matching those in the counts data.  If they are sorted in the same order, analysis will be faster
* **if a cell is left empty**, it will be populated as NA and exlcuded from analysis  

![alt text](https://github.com/tbrunetti/PCA_QC/blob/develop/metadata_file_format.png "metadata file format")  


## Running PCA-dynamic.R with test data

To initiate the script, run the following on the command line:
```
Rscript PCA-dynamic.R
```

This will prompt the user for file information and output directory and prefix information.  Below is an example for the provided test data:

```
Full path to counts file (no whitespaces): test_set/input_files/test_counts_file.csv
You entered the following:  test_set/input_files/test_counts_file.csv
Is this correct? (Y/N/Q) Y
Full path to meta data file (no whitespaces): test_set/input_files/test_metadata_file_with_types.csv
Is this correct? (Y/N/Q) Y
Name of output file (no whitespaces): test_set/results/test_data_output
Is this correct? (Y/N/Q) Y
```

This will print out some output (for logging and error purposes) to the screen/standard out.  All of the actual output files will be located in the output file location the user specified.  In this case is is test_set/results/ with all output files containing the prefix **test_data_ouput**.



