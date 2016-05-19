# TSScall
Developed by Christopher Lavender, National Institute of Environmental Health Sciences. Based on work by Adam Burkholder, Brian Bennett, and David Fargo.
## Overview
TSScall identifies transcription start sites (TSSs) from Start-seq data (Nechaev et al. *Science*, 2010). Operating both with and without a reference annotation, TSScall allows for rapid annotation of TSSs across an entire genome.
## TSS calling procedure
TSScall first divides a genome into a series of windows. If a reference annotation is used, windows are first generated at TSSs in the reference.  Windows are then made in other areas of the genome where Start-seq coverage is present but TSSs are annotated in the reference.

TSScall divides each window into 200-nt bins. A TSS is called within the bin with the highest total read counts and at the individual nt position with the highest number of reads.  In unannotated windows, this process is performed iteratively: after calling a TSS, all reads within a distance threshold are removed from the window, and calling repeats until no reads are left.
## Requirements
TSScall was developed in Python 3.4.4.

TSScall requires input Start-seq data to first be processed. Two bedGraph files are required, one for the forward strand and one for the reverse strand. TSScall expects each bedGraph file to contain the 5' end counts of reads.

TSScall also requires a chromosome sizes file.  The genome size is taken from this file, with the genome size used in read threshold determination.
## Parameters
### --call_method
Selects call method:

*bin_winner* (Default):

To call a TSS, each search window is divided into 200-nt bins.  In the bin with the most reads, the TSS is called at the position with the most reads.

*global*:

The position within a given search window with the most reads is called as the TSS.
### --annotation, -a ANNOTATION
Sets annotation fie.  Annotation file must be in GTF format.
### --detail_file DETAIL_FILE
If a detail_file is set, additional information not found in the output BED file will be added to a tab-delimited TXT file.
### --fdr FDR
If set, determines read threshold for TSS calling by FDR.  Threshold will be the minimum read count with an expected FDR less than the set FDR value.

This is the default method.  Sets FDR at 0.001
### --false_positives FALSE_POSITIVES
If set, determines read threshold for TSS calling by allowed false positives.  Threshold will be the minimum read count with the expected false positives less than the set value.
### --set_read_threshold SET_READ_THRESHOLD
Sets the read threshold at user-defined count.
### --bidirectional_threshold BIDIRECTIONAL_THRESHOLD
TSScall associates TSSs into bidirectional pairs.  These associations are noted in the detail file.  This parameter sets the allowed distance threshold for bidirectional TSS association.

Default: 1000
### --cluster_threshold CLUSTER_THRESHOLD
TSScall associates TSSs into clusters for super enhancer determination.  This parameter sets the allowed distance threshold for cluster association.

Default: 1000
### --annotation_search_window ANNOTATION_SEARCH_WINDOW
Search windows are made about each TSS in an annotation file.  This parameter determines the size of the annotation search window by setting the maximum distance away from an annotated TSS that a TSS may be called.

Default: 1000
### --annotation_join_distance ANNOTATION_JOIN_DISTANCE
In some annotations, TSSs may be very close in space.  If within the annotation join distance, annotation search windows will be merged, resulting in a single call for the merged windows.

Default: 200
### --nutss_search_window NUTSS_SEARCH_WINDOW
Sets the distance between which unannotated TSSs may be called.  No two nuTSSs may be closer than the nuTSS search window apart.

Default: 250
### --nutss_filter_size NUTSS_FILTER_SIZE
Prior to calling of unannotated TSSs, Start-seq reads are filtered based on annotated TSSs and TSSs called from annotation.  If a Start-seq read is within the nuTSS filter size of a TSS in the annotation or called from the annotation, it is filtered prior to unannotated TSS calling.

Default: 750
### --cluster_bed CLUSTER_BED
Writes clusters to unstranded BED file CLUSTER_BED
