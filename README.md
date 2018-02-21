# TSScall

Developed by Christopher Lavender, National Institute of Environmental Health Sciences. Based on work by Adam Burkholder, Brian Bennett, and David Fargo.

## Overview

TSScall identifies transcription start sites (TSSs) from Start-seq data (Nechaev et al. *Science*, 2010). TSScall allows for rapid annotation of TSSs across an entire genome.

## Obtaining TSScall

TSScall may be downloaded directly from this GitHub repository. TSScall was designed to be run as a standalone Python script that does not require installation or outside dependencies. The only requirement is a standard installation of Python (2.7 or 3).

## Using TSScall

TSScall is run from the commandline.

``python TSScall.py [OPTIONS] FORWARD_BEDGRAPH REVERSE_BEDGRAPH CHROMOSOME_SIZES OUTPUT_FILE``

`FORWARD_BEDGRAPH` and `REVERSE_BEDGRAPH` are bedGraph files describing the coverage of 5' ends from Start-seq reads at single-nucleotide resolution.  Stranded information is required. `FORWARD_BEDGRAPH` and `REVERSE_BEDGRAPH` describe the coverage for the forward (+) strand and the reverse (-) strand, respectively.

5' end coverage may be calculated using the [bedtools utility package](http://bedtools.readthedocs.io/en/latest/). Coverage files may be generated from an alignment in BAM format.

```
bedtools genomecov -ibam BAM_FILE -bg -5 -strand + > FORWARD_BEDGRAPH
bedtools genomecov -ibam BAM_FILE -bg -5 -strand - > REVERSE_BEDGRAPH
```

`CHROMOSOME_SIZES` is a text file describing the size of each chromosome/contig in the reference genome. A sample file in this format can be found [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes). For well-established genomes, we recommend generating this file using the USCS `fetchChromosomeSizes` script, found [here](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/).

Usage information may be found by running `python TSScall.py --help`.

### Optional arguments

`--annotation, -a ANNOTATION_FILE`

Sets annotation fie.  The annotation file must be in GTF format.

Setting the annotation file affects how TSScall identifies TSSs. If an annotation file is used, TSScall will first search for TSSs within regions centered on annotation-defined TSSs before searching for unannotated TSSs.

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

`--utss_search_window NUTSS_SEARCH_WINDOW`

Sets the distance between which unannotated TSSs may be called.  No two uTSSs may be closer than the uTSS search window apart. The default value is 250.

`--utss_filter_size NUTSS_FILTER_SIZE`

If an annotation is specified, TSScall will first search for TSSs near annotated start sites of gene models before calling additional unannotated TSSs. Prior to calling unannotated TSSs, reads are filtered and not considered in unannotated TSS calling if they are near annotated TSSs. Reads are filtered if they are within the uTSS filter size of either an annotated TSS or TSS called from annotation. The default uTSS filter size is 750.

### Output files

#### BED output file

The standard output file for TSScall is a list of TSSs in [BED format](http://genome.ucsc.edu/FAQ/FAQformat#format1). Entries are sorted by chromosome and position. TSSs are arbitrarily named by order of strand and position. By default, a TSS name will have the 'uTSS' prefix. If a TSS corresponds to a transcript in an input annotation, its name will have the 'obsTSS' prefix. The associated transcript can be found in the optional detail file.  

No information is stored in the score column (column 5). A placeholder value of '0' is used to maintain the BED file format.

#### Detail file

The optional detail file contains additional information about each TSS in a tab-delimited TXT file. Column headers are as follows:

* **TSS ID** Name assigned by TSS call.
* **Type** Type of TSS called. Available types are ...
* **Transcripts** Semi-colon delimited list of transcript IDs if the TSS was called from annotation.
* **Gene ID** Semi-colon delimited list of gene IDs if the TSS was called from annotation.
* **Strand**
* **Chromosome**
* **Position**
* **Reads** Coverage of 5' ends at the called TSS.
* **Divergent?** Boolean field describing whether or not a divergent TSS was called.
* **Divergent partner** If divergent, gives the assigned ID of the divergent TSS.
* **Divergent distance** The distance from the TSS to its divergent partner.
* **Convergent?** Boolean field describing whether or not a convergent TSS was called.
* **Convergent partner** If convergent, gives the assigned ID of the convergent TSS.
* **Convergent distance** The distance from the TSS to its convergent partner.
* **TSS cluster** TSScall will associate TSSs into clusters based on proximity. Gives the assigned ID of the cluster.
* **TSSs in associated cluster** Gives the number of TSSs in the associated cluster. May be helpful in identifying super enhancers.

The detail file also includes any additional information found in the attribute field of the [GTF file](https://useast.ensembl.org/info/website/upload/gff.html). Additional information will have headers consistent with the GTF annotation.

#### Cluster BED file

TSScall will associate TSSs into clusters based on proximity. This file gives the location of TSS clusters in [BED format](http://genome.ucsc.edu/FAQ/FAQformat#format1).

## Utilities for annotation of results

Within the /utils directory are scripts for annotation of TSScall results. Additional scripts may be present in the /utils directory; these scripts may be under development or designed for niche functions.

We recommend follow-up analysis with TSSclassify to associate called TSSs with gene models.

### TSSclassify

```perl TSSclassify.pl DETAIL_FILE GTF_ANNOTATION_FILE```

TSSclassify is a perl script that annotates a TSScall detail file with information from a GTF annotation.

### TSScompare

```python TSScompare.py BED_FILE_1 BED_FILE_2```

TSScompare considers two BED files and finds TSSs shared between them. An optional `--distance_threshold INTEGER` or `-d INTEGER` flag may be used to set a distance threshold; any TSSs within this distance threshold will be considered shared. Additional usage information can be found using `python TSScompare.py --help`.
