#PCA_QC
dynamic correlation analysis of PCs with metadata for QC metrics

A way to look at whether the first 10 or n PCs in your dataset correlate with any measured metadata metrics

Can be used with almost any continuous data matrix that has corresponding metadata.  Primarily built for use with RNA-seq data.  Run can be performed on a local standard copmuter; no HPC required.

##Required R libraries  
* factoextra
* ggfortify
* ggplot2
* data.table
* gplots

##Output Files  
3 output files are generated:  
1. output_name_pca_matrix_output.csv  
2. output_name__metadata_sorted_output.csv  
3. output_name.pdf  


##Required Inputs
* counts in raw counts,  FPKM, or TPM units for all genes/isoform used in experiment for every sample
* metadata file where metrics are listed for every sample in experiment

##Running PCA-dynamic.R with test data

To initiate the script, run the following on the command line:
```
Rscript PCA-dynamic.R
```

This will prompt the user for file information and output directory and prefix information.  Below is an example for the provided test data:

```
Full path to counts file: test_set/input_files/test_counts_file.csv
You entered the following:  test_set/input_files/test_counts_file.csv
Is this correct? (Y/N/Q) Y
Full path to meta data file: test_set/input_files/test_metadata_file_with_types.csv
Is this correct? (Y/N/Q) Y
Name of output file (no whitespaces): test_set/results/test_data_output
Is this correct? (Y/N/Q) Y
```

This will print a some output (for logging and error purposes) to the screen.  All of the actualy output files will be in the output file location the user specified.  In this case is is test_set/results/ with all output files containing the prefix **test_data_ouput**.



