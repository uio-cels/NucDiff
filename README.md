# NucDiff manual
## 1 Introduction
NucDiff locates and categorizes differences between two closely related nucleotide sequences. It is able to deal with very fragmented genomes, structural rearrangements and various local differences. These features make NucDiff to be perfectly suitable to compare assemblies with each other or with available reference genomes.   

NucDiff provides information about the types of differences and their locations. It is possible to upload the results into genome browser for visualization and further inspection. It was written in Python and uses the NUCmer package from MUMmer[1] for sequence comparison. 




## 2 Download and installation
### 2.1 Prerequisites
NucDiff can be run on Linux and Mac OS. It uses MUMmer v3.23 which should be installed and be in the PATH before running NucDiff.  

The MUMmer tarball can be downloaded at http://sourceforge.net/projects/mummer/ .


### 2.2 Installation 
To install NucDiff, download the source code tarball and unpack it like this:

```
$wget  <reference  in the internet>
$tar -xzf <name of tarball>.tar.gz
```


## 3 Running NucDiff
### 3.1 Command line syntax and input arguments
To run NucDiff, run `nucdiff.py` script with valid input arguments:

```
$ python  nucdiff.py [-h] [--reloc_dist [Int]]
                          [--nucmer_opt [NUCMER_OPT]]
                          [--filter_opt [FILTER_OPT]]
                          [--delta_file [Delta_file.delta]] [--proc [Int]]
                          [--ref_name_full [{yes,no}]]
                          [--query_name_full [{yes,no}]] [--version]
                          Reference.fasta Query.fasta Output_dir Prefix
```

Positional arguments:
*  **Reference.fasta** - fasta file with the reference sequences
*  **Query.fasta** - fasta file with the query sequences
*  **Output_dir** - path to the directory where all intermediate and final results will be stored
*  **Prefix** - name that will be added to all generated files including the ones created by NUCmer

Optional arguments:
*  **-h, --help** - show this help message and exit
*  **--reloc_dist** - minimum distance between two relocated fragments [10000]
*  **--nucmer_opt** - nucmer run options. By default, NUCmer will be run with its default parameters values, except the --maxmatch parameter. --maxmatch is hard coded and cannot be changed. To change any other parameter values, type parameter names and new values inside single or double quotation marks.
*  **--filter_opt** - delta-filter run options. By default, it will be run with -q parameter only.
*  **--delta_file** - path to the already existing delta file (NUCmer output file)
*  **--proc** - number of processes to be used [1]
*  **--ref_name_full** - print full reference names in output files ('yes' value). In case of 'no', everything after the first space will be ignored. ['no']
*  **--query_name_full** - print full query names in output files ('yes' value). In case of 'no', everything after the first space will be ignored. ['no']
*  **--version** - show program's version number and exit



### 3.2 Running examples
Running example with NucDiff and NUCmer predefined parameters values, except NUCmer --maxmatch parameter.  --maxmatch is hard coded and cannot be changed neither to --mum nor to --mumreference:

```
$python nucdiff.py my_reference.fasta my_query.fasta my_output_dir my_prefix
```

Running example when user needs to change NUCmer and NucDiff default parameter values:

```
$python nucdiff.py --proc 5 --ref_name_full 'yes' --query_name_full 'yes' --nucmer_opt '-c 200 -l 250' my_reference.fasta my_query.fasta my_output_dir my_prefix
```


## 4 Method overview

### 4.1 NucDiff steps

The NucDiff workflow is shown in Figure 1. The detailed description of all steps can be found in [2].

![](Figures_manual/workflow.png)

Figure 1: The NucDiff workflow

### 4.2 Types of differences

All types of differences are classified into 3 groups: Global, Local and Structural (Figure 2). 


![](Figures_manual/types_of_differences.png)

Figure 2: Classification of the types of differences with group names found in coloured boxes with capitalised names and the specific types found in white boxes with lowercase names.

#### 4.2.1 Global differences

Global differences affect the whole query sequence. This group consists of only one type, called unaligned sequence.
* **unaligned sequence** is a query sequence that has no matches of length equal to or longer than a given number of bases (65 by default) with the reference genome.  


#### 4.2.2 Local differences

Local differences involve various types of insertions, deletions and substitutions. NucDiff distinguishes between six types of insertions (the insertion subgroup in Figure 2):
* **simple insertion** - an insertion of bases in the query sequence that were not present anywhere on the reference genome
* **duplication** - an insertion in the query sequence of an extra copy of some reference sequence not adjacent to this region, creating an interspersed repeat, or increasing the copy number of an interspersed repeat  
* **tandem duplication** - an insertion of an extra copy of some reference sequence region adjacent to this region in the query sequence  
* **inserted gap** - an insertion of unknown bases (Ns) in the query sequence in a region which is continuous (without a gap) in the reference, or which results in an elongation of the region of unknown bases in the reference
* **unaligned beginning** - unaligned bases in the beginning of a query sequence
* **unaligned end** - unaligned bases at the end of query sequence

There are several types of deletions (the Deletion subgroup in Figure 2):
* **simple deletion** - a deletion of some bases, present in the reference sequence, from a query sequence
* **collapsed repeat** - a deletion of one copy of an interspersed repeat from the reference sequence in a query sequence
* **collapsed tandem repeat** - a deletion of one or more tandem repeat units from the reference sequence in a query sequence

And, last, there are two types of substitutions (the substitution subgroup in Figure 2):
* **substitution** - a substitution of some reference sequence region with another sequence of the same exact length not present anywhere in the reference genome  (note that this sequence is not categorised as unaligned sequence because it is within a fragment that overlaps between query and reference). SNPs can be considered as a subcategory of simple substitutions.
* **gap** - a substitution of some reference sequence region with unknown sequence (Ns) of the same length. If the query has an enlarged gap, then this will be classified as gap+inserted gap differences, while a shortened gap is classified as gap+simple deletion differences.

#### 4.2.3 Structural differences

NucDiff detects several structural differences. These can be grouped into intra- and inter-chromosomal differences, and  some of these contain groups of types:
* **Translocation** -  an inter-chromosomal structural rearrangement, when two regions located on the different reference sequences, are placed adjacent to each other (**simple translocation**), overlap (**translocation with overlap**), or contain an unmapped region between them (**translocation with inserted gap** in case of inserted Ns, and **translocation with insertion** otherwise) in the same query sequence
* **Relocation** - an  intra-chromosomal structural rearrangement, when two regions located in different parts of the same reference sequence are placed adjacent to each other (**simple relocation**), overlap (**relocation with overlap**),  or or contain an unmapped region between them (**relocation with inserted gap** in case of inserted Ns, and **relocation with insertion** otherwise) in the same query sequence
* **reshuffling** - an intra-chromosomal structural rearrangement, which occurs when several neighbouring reference sequence regions are placed in a different order in a query sequence
* **inversion** -  an intra-chromosomal structural rearrangement, which occurs when a query sequence region is the reverse complement of a reference sequence region

The Translocation is a part of the Inter-chromosomal subgroup, while Relocation, reshuffling and inversion belong to the Intra-chromosomal subgroup in Figure 2.

## 5. NucDiff output
NucDiff puts its output in the directory `<output_dir>/results`. The output consists of 4 files: `<prefix>_query_coord.gff`, `<prefix>_ref_coord.gff`, `<prefix>_mapped_blocks.gff`,     `<prefix>_stat.out`. 


### 5.1 < prefix>_query_coord.gff
The ```<prefix>_query_coord.gff``` is a file in gff3 format containing information about the differences between query and reference sequences in query-based coordinates. An example of such a file is given below: 

```
##gff-version 3
##sequence-region	scf7180000001043	1	155799
scf7180000001043	.	Differences		90		173		.	+	.	Name=relocation-overlap;length=84
scf7180000001043	.	Differences		90		90		.	+	.	Name=relocation-overlap_st;length=84;ref_sequence=NC_010079;ref_coord=2181532-2181532
scf7180000001043	.	Differences		148		148		.	+	.	Name=substitution;length=1;ref_sequence=NC_010079;ref_coord=1961731-1961731
scf7180000001043	.	Differences		169		169		.	+	.	Name=substitution;length=1;ref_sequence=NC_010079;ref_coord=1961752-1961752
scf7180000001043	.	Differences		173		173		.	+	.	Name=relocation-overlap_end;length=84;ref_sequence=NC_010079;ref_coord=1961756-1961756
scf7180000001043	.	Differences		21217	21764	.	+	.	Name=gap;length=548;ref_sequence=NC_010079;ref_coord=2202659-2203206
scf7180000001043	.	Differences		21764	21764	.	+	.	Name=deletion;length=31;ref_sequence=NC_010079;ref_coord=2203207-2203237
scf7180000001043	.	Differences		21768	21768	.	+	.	Name=substitution;length=1;ref_sequence=NC_010079;ref_coord=2203241-2203241
scf7180000001043	.	Differences		40375	40375	.	+	.	Name=insertion;length=1;ref_sequence=NC_010079;ref_coord=2221847-2221847
scf7180000001043	.	Differences		40376	40395	.	+	.	Name=inserted_gap;length=20;ref_sequence=NC_010079;ref_coord=2221847-2221847
scf7180000001043	.	Differences		40396	40405	.	+	.	Name=insertion;length=10;ref_sequence=NC_010079;ref_coord=2221847-2221847
scf7180000001043	.	Differences		66534	66534	.	+	.	Name=deletion;length=1;ref_sequence=NC_010079;ref_coord=2247977-2247977
...
```

Information about each type of differences is represented in the following way:
* **Unaligned sequence**:
   * the name of the query sequence (1 column)
   * the query sequence start being equal to 1 (4 column)
   * the query sequence end being equal to the length of the sequence (5 column)
   * the difference type being equal to "unaligned_sequence" (9 column, Name)
   * the query sequence length (9 column, length)  
* **All types of Insertions**:
   * the name of the query sequence where the insertion was detected (1 column) 
   * the query starting and ending positions of the inserted region (4, 5 columns) 
   * the insertion type. The "simple insertion" will be represented as "insertion". (9 column, Name)
   * the length of the inserted region (9 column, length)
   * the name of the reference sequence where the insertion was detected (9 column, ref_sequence)
   * the reference coordinate of the query base preceding the inserted bases (9 column, ref_coord)
* **All types of Deletions**:
   * the name of query sequence where the deletion was detected (1 column)
   * the query coordinate of the reference base preceding the deleted bases (4, 5 columns)
   * the deletion type (9 column, Name)
   * the length of the deleted region (9 column, length)
   * the name of the reference sequence where the deletion was detected (9 column, ref_sequence)
   * the reference starting and ending coordinates of the deleted region (9 column, ref_coord)
* **All types of Substitutions**:
   * the name of query sequence where the substitution was detected (1 column)
   * the query starting and ending positions of the substituted region (4, 5 columns)
   * the substitution type (9 column, Name)
   * the length of the substituted region (9 column,length)
   * the name of the reference sequence where the substitution was detected (9 column, ref_sequence)
   * the reference starting and ending coordinates of the substituted region (9 column, ref_coord)
* **Simple translocations/relocations**:
<br />   For each pair of translocated/relocated fragments:
   * first line: 
      * the name of the query sequence where the difference was detected (1 column) 
	  * the query ending position of the first fragment (4 column) 
	  * the query starting position of the second fragment (5 column) 
	  * the difference type being equal to "translocation" or "relocation" (9 column, Name)
	  * the length of the region between two translocated/relocated fragments being equal to 0 (9 column, length)
   * second line:
      * the name of the query sequence where the difference was detected (1 column)  
	  * the query ending position of the first fragment (4, 5 columns)
      * the entry type being equal to "translocation_end" or "relocation_end" (9 column, Name)
      * the length of the region between two fragments being equal to 0 (9 column, length) 
      * the name of the reference sequence where the end of the first fragment corresponds to (9 column, ref_sequence)
      * the reference coordinate of the end of the first fragment (9 column, ref_coord)	  
   * third line:
      * the name of the query sequence where the difference was detected (1 column)  
	  * the query start position of the second fragment (4, 5 columns)
      * the entry type being equal to "translocation_st" or "relocation_st" (9 column, Name)
      * the length of the region between two fragments being equal to 0 (9 column, length) 
      * the name of the reference sequence where the start of the second fragment corresponds to (9 column, ref_sequence)
      * the reference coordinate of the start of the second fragment (9 column, ref_coord)	  
* **Translocations/Relocations with insertion**
<br />   For each pair of translocated/relocated fragments with insertion:
   * first line: 
      * the name of the query sequence where the difference was detected (1 column) 
	  * the query starting and ending positions of the inserted region (4, 5 columns) 
	  * the difference type being equal to "translocation-insertion" or "relocation-insertion" (9 column, Name)
	  * the length of the inserted region (9 column, length)
   * second line:
      * the name of the query sequence, where the difference was detected (1 column)  
	  * the query starting position of the inserted region (4, 5 columns)
      * the entry type being equal to "translocation-insertion_st" or "relocation-insertion_st" (9 column, Name)
      * the length of the inserted region (9 column, length) 
      * the name of the reference sequence where the query base preceding the start of the inserted region corresponds to (9 column, ref_sequence)
      * the reference coordinate of the query base preceding the start of the inserted region (9 column, ref_coord)	  
   * third line:
      * the name of the query sequence where the difference was detected (1 column) 
	  * the query ending position of the inserted region (4, 5 columns)
      * the entry type being equal to "translocation-insertion_end" or "relocation-insertion_end" (9 column, Name)
      * the length of overlapped region (9 column, length) 
      * the name of the reference sequence where the query base next to the end of the inserted region corresponds to (9 column, ref_sequence)
      * the reference coordinate of the query base next to the end of the inserted region  (9 column, ref_coord)	  
* **Translocations/Relocations with overlap**
<br />   For each pair of translocated/relocated fragments with overlap:
   * first line: 
      * the name of the query sequence where the difference was detected (1 column) 
	  * the query starting and ending positions of the overlapped region (4, 5 columns) 
	  * the difference type being equal to "translocation-overlap" or "relocation-overlap" (9 column, Name)
	  * the length of the overlapped region (9 column, length)
   * second line:
      * the name of the query sequence where the difference was detected (1 column)  
	  * the query starting position of the overlapped region (4, 5 columns)
      * the entry type being equal to "translocation-overlap_st" or "relocation-overlap_st" (9 column, Name)
      * the length of the overlapped region (9 column, length) 
      * the name of the reference sequence where the start of the overlapped region corresponds to (9 column, ref_sequence)
      * the reference coordinate of the overlapped region start (9 column, ref_coord)	  
   * third line:
      * the name of the query sequence where the difference was detected (1 column) 
	  * the query ending position of the overlapped region (4, 5 columns)
      * the entry type being equal to "translocation-overlap_end" or "relocation-overlap_end" (9 column, Name)
      * the length of the overlapped region (9 column, length) 
      * the name of the reference sequence where the end of the overlapped region corresponds to (9 column, ref_sequence)
      * the reference coordinate of the overlapped region end (9 column, ref_coord)	  
* **Reshufflings**:
   * the name of the query sequence where the reshuffled region was detected (1 column) 
   * the query starting and ending positions of the reshuffled region (4, 5 columns) 
   * the difference type being equal to "reshuffling-part_x_gr_y" where part corresponds to the order number in the reference sequence , and gr corresponds to the reference reshuffled region number (9 column, Name)
   * the length of the reshuffled region (9 column, length)
   * the name of the reference sequence where the reshuffled region corresponds to (9 column, ref_sequence)
   * the reference coordinates of the reshuffled region (9 column, ref_coord)   
* **Inversions**:
   * the name of the query sequence where the inversion was detected (1 column)  
   * the query starting and ending positions of the inverted region (4, 5 columns)
   * the difference type being equal to "inversion" (9 column, Name)
   * the length of the inverted region (9 column, length) 
   * the name of the reference sequence where the inverted region corresponds to (9 column, ref_sequence)
   * the reference starting and ending coordinates of the inverted region (9 column, ref_coord)	  
	
There can be 3 types of  "translocation/relocation with insertion" difference, depending on the bases inside the inserted region:
* "translocation-insertion" or "relocation-insertion" - the inserted region consists of ATGC's only
* "translocation-inserted_gap" or "relocation-inserted_gap" - the inserted region consists of N's only  
* "translocation-insertion_ATGCN" or "relocation-insertion_ATGCN" - the inserted region consists of ATGC's and N's

All these cases are handled in the same way as it was described for "translocation/relocation with insertion" case with a proper name for each difference type.


### 5.2 < prefix>_ref_coord.gff
The ```<prefix>_ref_coord.gff``` is a file in gff3 format, containing information about the differences between query and reference sequences in the reference-based coordinates.
Here is an example of such a file:

```
##gff-version 3
##sequence-region	NC_010079	1	2872915
NC_010079	.	Differences		113168	113168	.	+	.	Name=substitution;length=1;query_seq=scf7180000001044;query_coord=690816-690816
NC_010079	.	Differences		155856	155856	.	+	.	Name=inserted_gap;length=8;query_seq=scf7180000001044;query_coord=733505-733512
NC_010079	.	Differences		155857	155869	.	+	.	Name=gap;length=13;query_seq=scf7180000001044;query_coord=733513-733525
NC_010079	.	Differences		169962	169962	.	+	.	Name=insertion;length=35;query_seq=scf7180000001044;query_coord=747691-747725
NC_010079	.	Differences		170013	170013	.	+	.	Name=relocation_st;length=1;query_seq=scf7180000001044;query_coord=432669-432669
NC_010079	.	Differences		170013	170013	.	+	.	Name=relocation-overlap_st;length=51;query_seq=scf7180000001044;query_coord=432669-432669
NC_010079	.	Differences		170026	170026	.	+	.	Name=substitution;length=1;query_seq=scf7180000001044;query_coord=432656-432656
NC_010079	.	Differences		170026	170026	.	+	.	Name=substitution;length=1;query_seq=scf7180000001044;query_coord=747789-747789
NC_010079	.	Differences		170064	170064	.	+	.	Name=inserted_gap;length=245;query_seq=scf7180000001044;query_coord=747828-748072
NC_010079	.	Differences		170065	170439	.	+	.	Name=deletion;length=375;query_seq=scf7180000001044;query_coord=747827-747827
NC_010079	.	Differences		170229	170229	.	+	.	Name=relocation-inserted_gap;length=20;query_seq=scf7180000001044;query_coord=432433-432452
NC_010079	.	Differences		170229	170229	.	+	.	Name=relocation_end;length=1;query_seq=scf7180000001044;query_coord=432453-432453
NC_010079	.	Differences		170579	170580	.	+	.	Name=substitution;length=2;query_seq=scf7180000001044;query_coord=748212-748213
NC_010079	.	Differences		170585	170585	.	+	.	Name=substitution;length=1;query_seq=scf7180000001044;query_coord=748218-748218
NC_010079	.	Differences		222814	222814	.	+	.	Name=query_st;length=1;query_seq=scf7180000001029;query_coord=1279-1279
...
```

Information about each type of differences is represented in the following way:
* **All types of Insertions**:
   * the name of the reference sequence where the insertion was detected (1 column)
   * the reference coordinate of the query base preceding the inserted bases (4, 5 columns)   
   * the insertion type. The "simple insertion" will be represented as "insertion". (9 column, Name)
   * the length of the inserted region (9 column, length)
   * the name of the query sequence where the insertion was detected (9 column, query)
   * the query starting and ending positions of the inserted region (9 column, query_coord)
* **All types of Deletions**:
   * the name of reference sequence where the deletion was detected (1 column)
   * the reference starting and ending coordinates of the deleted region (4, 5 columns)
   * the deletion type. The "simple deletion" will be represented as "deletion". (9 column, Name)
   * the length of the deleted region (9 column, length)
   * the name of the query sequence where the deletion was detected (9 column, query)
   * the query coordinate of the reference base preceding the deleted bases (9 column, query_coord)
* **All types of Substitutions**:
   * the name of reference sequence where the substitution was detected (1 column)
   * the reference starting and ending positions of the substituted region (4, 5 columns)
   * the substitution type (9 column, Name)
   * the length of the substituted region (9 column,length)
   * the name of the query sequence where the substitution was detected (9 column, query)
   * the query starting and ending coordinates of the substituted region (9 column, query_coord)
* **Simple translocations/relocations**:
<br />    For each translocated/relocated fragment:
    * first line:
	  * the name of the reference sequence where the fragment is located (1 column)
      * the reference coordinate of the query starting position of the fragment (4, 5 columns)
	  * the entry type being equal to "translocation_st" ("relocation_st"), if the query fragment starting position corresponds to the reference fragment starting position , or "translocation_end" ("relocation_end"), if the query fragment starting position corresponds to the reference fragment ending position. Basically, it denotes the start or the end of the corresponding reference fragment depending on the query fragment orientation comparing to the reference fragment. (9 column, Name)
	  * the length of entry being equal to 1 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query fragment starting position of the fragment (9 column, query_coord) 
 	* second line:
	  * the name of the reference sequence where the fragment is located (1 column)
      * the reference coordinate of the query ending position of the fragment (4, 5 columns)
	  * the entry type being equal to "translocation_st" ("relocation_st"), if the query fragment starting position corresponds to the reference fragment ending position , or "translocation_end" ("relocation_end"), if the query fragment starting position corresponds to the reference fragment starting position (9 column, Name)
	  * the length of entry being equal to 1 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query fragment ending position of the fragment (9 column, query_coord) 
	* third line: 
      * the name of the reference sequence where the fragment is located (1 column)
	  * the reference coordinate of the corresponding query misjoin point (either start or end of reference fragment) (4, 5 columns)
	  * the difference type being equal to "translocation" or "relocation" (9 column, Name)
	  * the length of the region between two translocated/relocated fragments, being equal to 0 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query coordinate of the misjoin point corresponding to the given reference fragment (either start or end of the query fragment) (9 column, query_coord) 
* **Translocations/Relocations with insertion**:
<br />    For each translocated/relocated fragment:
    * first line:
	  * the name of the reference sequence where the fragment is located (1 column)
      * the reference coordinate of the query starting position of the fragment (4, 5 columns)
	  * the entry type being equal to "translocation_st" ("relocation_st"), if the query fragment starting position corresponds to the reference fragment starting position , or "translocation_end" ("relocation_end"), if the query fragment starting position corresponds to the reference fragment ending position. Basically, it denotes the start or the end of the corresponding reference fragment depending on the query fragment orientation comparing to the reference fragment. (9 column, Name)
	  * the length of entry being equal to 1 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query fragment starting position of the fragment (9 column, query_coord) 
 	* second line:
	  * the name of the reference sequence where the fragment is located (1 column)
      * the reference coordinate of the query ending position of the fragment (4, 5 columns)
	  * the entry type being equal to "translocation_st" ("relocation_st"), if the query fragment starting position corresponds to the reference fragment ending position , or "translocation_end" ("relocation_end"), if the query fragment starting position corresponds to the reference fragment starting position (9 column, Name)
	  * the length of entry being equal to 1 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query fragment ending position of the fragment (9 column, query_coord) 
	* third line: 
      * the name of the reference sequence where the fragment is located (1 column)
	  * the reference coordinate of the base preceding the insertion or the base next after insertion (depending on the location of insertion region relatively to the corresponding reference fragment) (4, 5 column)
	  * the difference type being equal to "translocation_insertion" or "translocation_inserted_gap" or "translocation_ATGCN" based on the nature of inserted region ("relocation_insertion" or "relocation_inserted_gap" or "relocation_ATGCN" in case of relocation) (9 column, Name)
	  * the length of the inserted region (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query coordinate of the inserted region (9 column, query_coord) 
* **Translocations/Relocations with overlap**:
<br />    For each translocated/relocated fragment:
    * first line:
	  * the name of the reference sequence where the fragment is located (1 column)
      * the reference coordinate of the query starting position of the fragment (4, 5 columns)
	  * the entry type being equal to "translocation_st" ("relocation_st"), if the query fragment starting position corresponds to the reference fragment starting position , or "translocation_end" ("relocation_end"), if the query fragment starting position corresponds to the reference fragment ending position. Basically, it denotes the start or the end of the corresponding reference fragment depending on the query fragment orientation comparing to the reference fragment. (9 column, Name)
	  * the length of entry being equal to 1 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query fragment starting position of the fragment (9 column, query_coord) 
 	* second line:
	  * the name of the reference sequence where the fragment is located (1 column)
      * the reference coordinate of the query ending position of the fragment (4, 5 columns)
	  * the entry type being equal to "translocation_st" ("relocation_st"), if the query fragment starting position corresponds to the reference fragment ending position , or "translocation_end" ("relocation_end"), if the query fragment starting position corresponds to the reference fragment starting position (9 column, Name)
	  * the length of entry being equal to 1 (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query fragment ending position of the fragment (9 column, query_coord) 
	* third line: 
      * the name of the reference sequence where the fragment is located (1 column)
	  * the reference coordinates of the overlapped region (4, 5 columns)
	  * the difference type being equal to "translocation_overlap" or "relocation_overlap" (9 column, Name)
	  * the length of the overlapped region (9 column, length)
	  * the name of the query sequence where the fragment is located (9 column, query) 
	  * the query coordinate of the overlapped region (9 column, query_coord) 
* **Reshufflings**:
<br />    For each reshuffled region:
   * the name of the reference sequence where the reshuffled region was detected (1 column) 
   * the reference starting and ending positions of the reshuffled region (4, 5 columns) 
   * the difference type being equal to "reshuffling-part_x_gr_y", where part correspond to the order number in the reference sequence , and gr corresponds to the reference reshuffled region number (9 column, Name)
   * the length of the reshuffled region (9 column, length)
   * the name of the query sequence where the reshuffled region corresponds to (9 column, query)
   * the query coordinates of the reshuffled region (9 column, query_coord)   
* **Inversions**:
   * the name of the reference sequence where the inversion was detected (1 column)  
   * the reference starting and ending positions of the inverted region (4, 5 columns)
   * the difference type being equal to "inversion" (9 column, Name)
   * the length of the inverted region (9 column, length) 
   * the name of the query sequence where the inverted region corresponds to (9 column, ref_sequence)
   * the query starting and ending coordinates of the inverted region (9 column, ref_coord)	  

Some additional information is output as well:
* Reference coordinates of query sequence start and end position ("query_st", "query_end" entry type) with the corresponding query coordinates
* Unmapped reference regions ("uncovered_region" entry type)    
   

### 5.3 < prefix>_mapped_blocks.gff
The ```<prefix>_mapped_blocks.gff```  is a gff file that contains information about the coverage of reference sequences by query sequences or query sequences blocks, formed by splitting the query sequences at the points of translocations, relocations, and reshufflings.
The example of such a file is given below:

```
##gff-version 3
##sequence-region	gi|169786809|ref|NC_010394.1|	1	23319
gi|169786809|ref|NC_010394.1|	.	MappedBlock		1		30463	.	+	.	Name=scf7180000000886;length=30463;query_length=62487;query_coord=1-30431
gi|169786809|ref|NC_010394.1|	.	MappedBlock		29271	36827	.	+	.	Name=scf7180000000990;length=7557;query_length=7557;query_coord=1-7557
gi|169786809|ref|NC_010394.1|	.	MappedBlock		29619	30706	.	+	.	Name=scf7180000000885;length=1088;query_length=2047;query_coord=960-2047
gi|169786809|ref|NC_010394.1|	.	MappedBlock		35280	44675	.	+	.	Name=scf7180000000991;length=9396;query_length=9396;query_coord=1-9396
gi|169786809|ref|NC_010394.1|	.	MappedBlock		44416	45430	.	+	.	Name=scf7180000001074;length=1015;query_length=1017;query_coord=1-1015
gi|169786809|ref|NC_010394.1|	.	MappedBlock		45355	46063	.	+	.	Name=scf7180000000860;length=709;query_length=709;query_coord=1-709
gi|169786809|ref|NC_010394.1|	.	MappedBlock		46149	65070	.	+	.	Name=scf7180000000909;length=18922;query_length=18922;query_coord=1-18922
gi|169786809|ref|NC_010394.1|	.	MappedBlock		52610	53442	.	+	.	Name=scf7180000000797;length=833;query_length=1364;query_coord=532-1364
gi|169786809|ref|NC_010394.1|	.	MappedBlock		72450	92655	.	+	.	Name=scf7180000000952;length=20206;query_length=20207;query_coord=1-20207
```
Each line contains following information:
* the name of the reference sequence (1 column)
* the reference starting and ending coordinates of the mapped block (4,5  columns)
* the name of the query sequence (9 column, Name)
* the length of the mapped block (9 column, length)
* the length of the whole query sequence (9 column, query_length)
* the query coordinates of the mapped block (9 column, query_coord). 

### 5.4 < prefix>_stat.out
The ```<prefix>_stat.out``` file contains information about the number of the found differences between the query and reference sequences. The following information is provided in the file:
* The total number of all reshufflings, inversions, wrong sequences and all types of insertions, deletions, substitutions, translocations, and relocations (Total number).
* The number of all insertions, deletions, substitutions, translocations, and relocations separately for each group. 
* The total  number of reshufflings, reshuffled fragments (reshuffled blocks), inversions, and wrong sequences.
* The number of reference regions that were not covered by any assembly sequence (Uncovered ref regions num) and the total length of these regions (Uncovered ref regions len)
* A detailed information about the number of insertions, deletions, substitutions, relocations and translocations separately for each type.

As an additional information, the number of sequences in query and reference genomes is provided. 

The example of the file is given below:
```
Total number	800
Insertions	143
Deletions	91
Substitutions	525
Translocations	0
Relocations	26
Reshufflings	4
Reshuffled blocks	14
Inversions	10
Unaligned sequences	0

Uncovered ref regions num	13
Uncovered ref regions len	45120

DETAILED INFORMATION:	
substitution	483
gap	42

insertion	73
duplication	5
tandem_duplication	0
unaligned_beginning	0
unaligned_end	2
inserted_gap	63

deletion	77
collapsed_repeat	9
tandem_collapsed_repeat	5

translocation	0
translocation-insertion	0
translocation-insertion_ATGCN	0
translocation-inserted_gap	0
translocation-overlap	0

circular_genome_start	1
relocation	1
relocation-insertion	1
relocation-insertion_ATGCN	1
relocation-inserted_gap	9
relocation-overlap	14


ADDITIONAL INFORMATION:
query sequences 10
reference sequences 2

``` 

## References
[1] Kurtz S et al. Versatile and open software for comparing large genomes. Genome Biol. 2004;5(2):R12. doi 10.1186/gb-2004-5-2-r12.

[2] Khelik K et al. NucDiff: in-depth characterization and annotation of differences between two sets of DNA sequences. 2016 (unpublished article)
