# TSScall

Developed by Christopher Lavender, National Institute of Environmental Health Sciences. Based on work by Adam Burkholder, Brian Bennett, and David Fargo.

## Overview

TSScall identifies transcription start sites (TSSs) from Start-seq data (Nechaev et al. *Science*, 2010). TSScall allows for rapid annotation of TSSs across an entire genome.

## Obtaining TSScall

TSScall may be downloaded directly from this GitHub repository. TSScall was designed to be run as a standalone Python script that does not require installation or outside dependencies. The only requirement is a standard installation of Python (2.7 or 3).

## Using TSScall

TSScall is run from the commandline.

``python TSScall.py [OPTIONS] FORWARD_BEDGRAPH REVERSE_BEDGRAPH CHROMOSOME_SIZES OUTPUT_FILE``

`FORWARD_BEDGRAPH` and `REVERSE_BEDGRAPH` are bedGraph files describing the coverage of 5' ends from Start-seq reads.  Stranded coverage is required. `FORWARD_BEDGRAPH` and `REVERSE_BEDGRAPH` describe the coverage for the forward (+) strand and the reverse (-) strand, respectively.

`CHROMOSOME_SIZES` is a text file describing the size of each chromosome/contig in the reference genome. A sample file in this format can be found [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes). For well-established genomes, we recommend generating this file using the USCS `fetchChromosomeSizes` script, found [here](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/).

Usage information may be found by running `python TSScall.py --help`.

### Optional arguments

`--annotation, -a ANNOTATION_FILE`

Sets annotation fie.  The annotation file must be in GTF format.

Setting the annotation file affects how TSScall identifies TSSs. If an annotation file is used, TSScall will search for TSSs within regions centered on annotation-defined TSSs.

`--detail_file DETAIL_FILE`

If a detail_file is set, additional information not found in the output BED file will be added to a tab-delimited TXT file.

`--cluster_bed CLUSTER_BED`

If set, TSScall will associate TSSs into clusters and write these clusters to an unstranded BED file, CLUSTER_BED.

`--calling_method [bin_winner|global]`

Sets the calling method for identifying TSSs. Either `bin_winner` or `global` may be specified.

`bin_winner` is the default method used. Here, bins of read coverage are used to select TSSs. This method favors selecting TSSs where the coverage may be locally difuse.

If `global` is set, only individual positions are considered. This method is appropriate when the 5' end coverage is consistently focused on single nucleotide positions.

`--fdr FDR`

`--false_positives FALSE_POSITIVES`

`--set_read_threshold SET_READ_THRESHOLD`

`--fdr`, `--false_positives`, and `--set_read_threshold` set the method used for the determining the read threshold for a valid TSS. Only one of these three options may be set at a time.

If `--fdr` is set, the read threshold is determined by using the associated FDR.  The read threshold for a valid TSS will be the minimum read count with an expected FDR less than the set FDR value. This is the default method, with a default FDR of 0.001.

If `--false_positives` is set, the read threshold is determined by considering FALSE_POSITIVES to be the maximum number of allowed false positives.  The read threshold will be the minimum read count where the number of expected false positives will be less than the specified value.

If `--set_read_threshold` is set, the read threshold is set the user-defined count. This is useful for setting the threshold by eye or for using outside methods.

`--bidirectional_threshold BIDIRECTIONAL_THRESHOLD`

TSScall associates TSSs into bidirectional pairs. These associations are noted in the detail file. This parameter sets the allowed distance threshold for bidirectional TSS association. The default value is 1000.

`--cluster_threshold CLUSTER_THRESHOLD`

TSScall associates TSSs into clusters by proximity. These clusters may be helpful for super enhancer determination. This parameter sets the allowed distance threshold for cluster association. The default value is 1000.

`--annotation_search_window ANNOTATION_SEARCH_WINDOW`

Search windows are made about each TSS in an annotation file.  This parameter determines the size of the annotation search window by setting the maximum distance away from an annotated TSS that a TSS may be called. The default value is 1000.

`--annotation_join_distance ANNOTATION_JOIN_DISTANCE` 

In some annotations, TSSs may be very close in space.  If within the annotation join distance, annotation search windows will be merged, resulting in a single call for the merged windows. The default value is 200.

`--nutss_search_window NUTSS_SEARCH_WINDOW`

Sets the distance between which unannotated TSSs may be called.  No two nuTSSs may be closer than the nuTSS search window apart. The default value is 250.

`--nutss_filter_size NUTSS_FILTER_SIZE`

Prior to calling of unannotated TSSs, Start-seq reads are filtered based on annotated TSSs and TSSs called from annotation.  If a Start-seq read is within the nuTSS filter size of a TSS in the annotation or called from the annotation, it is filtered prior to unannotated TSS calling. The default value is 750.

## Utilities for annotation of results

Within the /utils directory are scripts for annotation of TSScall results. Additional scripts may be present in the /utils directory; these scripts may be under development or designed for niche functions.

We recommend follow-up analysis with TSSclassify to associate called TSSs with gene models.

### TSSclassify

```Usage```

### TSScompare

```Usage```
