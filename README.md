# Nexus
Pan-genome compacted de Bruijn graphs with support for approximate pattern matching using search schemes

<!--- TODO: reference paper 
Columba was introduced in our [paper](https://doi.org/10.1016/j.isci.2021.102687). If you find this code useful in your research, please cite: 
```
@article{RENDERS2021102687,
title = {Dynamic partitioning of search patterns for approximate pattern matching using search schemes},
journal = {iScience},
volume = {24},
number = {7},
pages = {102687},
year = {2021},
issn = {2589-0042},
doi = {https://doi.org/10.1016/j.isci.2021.102687},
url = {https://www.sciencedirect.com/science/article/pii/S2589004221006556},
author = {Luca Renders and Kathleen Marchal and Jan Fostier},
keywords = {Algorithms, Bioinformatics, Computer science, High-performance computing in bioinformatics},
```--->

These instructions will get you a copy of the project up and running on your local machine.

#  Prerequisites

This package requires a number of packages to be installed on your system. Required: CMake (3.0 or higher); Google's Sparsehash; gcc (GCC 6.4 or a more recent version) 

How to install these packages:

As a root, execute the following commands:

on Redhat / Fedora distributions
```bash
yum install cmake
yum install sparsehash-devel
``` 

on Ubuntu / Debian distributions
```bash
apt-get install cmake
apt-get install libsparsehash-dev
```  
# Installing Nexus:

The installation is now simple. First, clone Nexus from the GitHub address

    git clone "https://github.com/biointec/nexus.git"

From this directory, run the following commands:
```bash
mkdir build
cd build
cmake ..
make 
```
# Usage:
Nexus aligns reads to a compressed pan-genome de Bruijn graph. To do this you need to build the implicit representation of this graph, along with the underlying bidirectional FM-index, based on the input data. Currently we only support input data with an alphabet of length 6 (for DNA: A, C, G, T + separation characters: $, %). Both separation character must be present to guarantee correct functionality. In other words, at least two genome strains should be included.

Nexus integrates code from three external repositories:
* [Columba](https://github.com/biointec/columba) is base of the search schemes implementation
* [radixSA64](https://github.com/mariusmni/radixSA64) for building suffix arrays
* [sux](https://github.com/vigna/sux) for bit vector rank and select support

## Building the Representation
To build the implicit representation of the compressed de Bruijn graph, along with the underlying bidirectional FM-index, only the search text (the pan-genome) is required. We make use of [radixSA64](https://github.com/mariusmni/radixSA64) for building suffix arrays. This is integrated in our building code.

To build the implicit representation of the compressed de Bruijn graph, run the following command in the `build` folder. 
```bash
./nexusBuild [basefile] [k]
```

The basefile parameter demonstrates where the search text can be found. The k parameter determines the minimum node length of the compressed de Bruijn graph. 

Options:

```
[options]
  -s  --sa-sparseness	Suffix array sparseness factors to be used. This option can be repeated multiple times for multiple versions of the suffix array. This option takes values in {1, 2, 4, 8, 16, 32, 64, 128, 256}. Use "all" to use all options. [default = 1]
  -c  --cp-sparseness	Sparseness factor that indicates how many checkpoints must be stored to identify nodes. This option can be repeated multiple times for multiple versions of the checkpoint sparseness. Use "none" to use no checkpoints. [default = 128]
  -p  --progress	Report extra progress updates
```
### Example 1
After installing Nexus, the nexus directory should look like this:

    .
    ├── build
    ├── cmake
    ├── radixSA64
    ├── search_schemes
    ├── src
    └── sux
    
In this example we will build the implicit representation of the compressed de Bruijn graph for a pan-genome of four E. coli strains. To do this, we will create an `example` folder.
To create this folder, navigate to the nexus folder. Here, enter the following command
```bash
mkdir example
```
To this new directory, copy the example file found [here](https://github.com/biointec/nexus/releases/download/v1.0.0/EscherichiaColi4Strains.txt). These are four E. coli strains (accession numbers FM180568, CP000247, CP001671, and CP000468) where all non-ACGT characters were replaced accordingly. They were appended to each other with the '%' separation character in between them. The sentinel character '$' was also appended to the end of this file.


To build the implicit representation of the compressed de Bruijn graph for a minimum node length of `k = 20`, navigate to the `build` folder and run the following command:
```bash
./nexusBuild ../example/EscherichiaColi4Strains 20
```
The index files are then written to the same folder. Your directory structure will now look like:
 ```
    .
    ├── build
    ├── cmake
    ├── example 
    |   ├── EscherichiaColi4Strains.B.left 
    |   ├── EscherichiaColi4Strains.brt 
    |   ├── EscherichiaColi4Strains.DBG 
    |   ├── EscherichiaColi4Strains.rev.bwt 
    |   ├── EscherichiaColi4Strains.sa.bv.1 
    |   ├── EscherichiaColi4Strains.B.right.128 
    |   ├── EscherichiaColi4Strains.bwt 
    |   ├── EscherichiaColi4Strains.left.map 
    |   ├── EscherichiaColi4Strains.right.map.128 
    |   ├── EscherichiaColi4Strains.txt 
    |   ├── EscherichiaColi4Strains.B.right.full.128 
    |   ├── EscherichiaColi4Strains.cct 
    |   ├── EscherichiaColi4Strains.rev.brt 
    |   └── EscherichiaColi4Strains.sa.1
    ├── radixSA64
    ├── search_schemes
    ├── src
    └── sux
```

Congratulations! You have used Nexus to build the implicit representation of the compressed de Bruijn graph of a pan-genome of four E. coli strains!

---


## Using the Representation
Nexus can align reads in a fasta (`.FASTA`, `.fasta`, `.fa`) or fastq (`.fq`, `.fastq`) format. Part of this code (specifically, the search schemes implementation) is based on part of the code of Nexus is based on [Columba 1.0](https://github.com/biointec/columba/releases/tag/v1.0).
To align your reads, use the following format: 

```bash
./nexus [options] basefilename readfile.[ext]
```

options:

```
 [options]
 -sfr --strain-free	strain-free matching
 -e   --max-ed		maximum edit distance [default = 0]
 -s   --sa-sparseness	suffix array sparseness factor [default = 1]
 -c   --cp-sparseness	sparseness factor that indicates how many checkpoints must be stored to identify nodes. Use "none" to use no checkpoints. [default = 128]
 -f   --filter		filtering type that should be used to filter the occurrences. This option is only valid in case of strain-free matching. Options:
	linear		linear filtering is efficient but does not filter out all redundant occurrences. Additionally, in some exceptional cases, a non-optimal replacement occurrence can be chosen. This is the default option.
	complete	complete filtering leads to a set of occurrences with no redundancy. This option is very slow however and thus not recommended.
 -p   --partitioning		Add flag to do uniform/static/dynamic partitioning. Dynamic partitioning cannot be used with strain-free matching. [default = static]
 -m   --metric		Add flag to set distance metric (editnaive/editopt/hamming) [default = editopt]
 -ss  --search-scheme	Choose the search scheme
  options: 
	kuch1	Kucherov k + 1
	kuch2	Kucherov k + 2
	kianfar	Optimal Kianfar scheme
	manbest	Manual best improvement for Kianfar scheme (only for ed = 4)
	pigeon	Pigeonhole scheme
	01*0	01*0 search scheme
	naive	naive backtracking
	custom	custom search scheme, the next parameter should be a path to the folder containing this search scheme

[ext]
	one of the following: fq, fastq, FASTA, fasta, fa
Following input files are required:
	<base filename>.txt: input text T
	<base filename>.cct: character counts table
	<base filename>.sa.[saSF]: sparse suffix array, with suffix array sparseness factor [saSF] elements
	<base filename>.sa.bv.[saSF]: bitvector indicating which elements of the suffix array are stored.
	<base filename>.bwt: BWT of T
	<base filename>.rev.bwt: BWT of the reverse of T
	<base filename>.brt: Prefix occurrence table of T
	<base filename>.rev.brt: Prefix occurrence table of the reverse of T
	<base filename>.DBG: variable k and the compressed de Bruijn graph.
	<base filename>.B.left: bitvector B_left for the compressed de Bruijn graph.
	<base filename>.B.right.[cpSF]: bitvector B_right for the compressed de Bruijn graph, with checkpoint sparseness factor [cpSF].
	<base filename>.B.right.full.[cpSF]: bitvector B_right_full for the compressed de Bruijn graph, with checkpoint sparseness factor [cpSF].
	<base filename>.left.map: node identifier mapping corresponding to B_left.
	<base filename>.right.map.[cpSF]: node identifier mapping corresponding to B_right, with checkpoint sparseness factor [cpSF].



    
```

### Strain-Fixed versus Strain-Free Matching

The `-sfr` or `--strain-free` option indicates a crucial difference in functionality. If this option is not added, Nexus defaults to strain-fixed matching. This means that only substrings of the original strains in the pan-genome can be found as occurrences. If the option is included, Nexus switches to strain-free matching. In this case, occurrences can be non-contiguous with respect to the strains in the pan-genome. In other words, an occurrence does not necessarily entirely occur in one of the original strains of the pan-genome. The choice of which type of matching is used can influence the other options as is mentioned in their description above. 

The matching duration, the number of index nodes visited, the number of graph nodes visited, and the number of reported/unique matches (strain-fixed) or node paths (strain-free) will be printed to stdout along with some additional metrics. For strain-fixed matching specifically, the number of reported node paths is included as an additional metric since multiple matches can correspond to the same node path. 

The matches will be written to a custom output file in the folder where your readfile was. This output file for strain-fixed matching will be a tab-separated file with the fields: `Identifier` (identifies the read), `SubgraphID` (identifies the node paths that were found for this read), `Path` (the node path corresponding to this occurrence), `Strain` (the number of the strain in which this occurrence lies), `Position` (the position of this occurrence in the pan-genome), `Length` (the length of this occurrence), `ED` (the edit distance of this occurrence) and `reverseComplement` (1 if this occurrence was found on the reverse complement of the reference, 0 otherwise). Similarly for strain-free matching, a tab-separated file is generated with the fields: `Identifier`, `SubgraphID`, `Path`, `DistanceFromLeftEnd`, `Length`, `ED` and `reverseComplement`. The `Strain` and `Position` fields are not present since an occurrence does not necessarily belong to one strain/correspond to a certain position in the original pan-genome anymore. The `DistanceFromLeftEnd` field is new and indicates the distance from the start of the match to the start of the first node of the node path, partly substituting the previous `Position` field. For each optimal alignment under the maximum given edit distance a line will be present. This output file will be called `readfile_output.txt`.






### Example 2
Consider the final directory structure from [example 1](##Example-1). 
Copy this [file](https://github.com/biointec/nexus/releases/download/v1.0.0/EscherichiaColi4Strains.reads.fasta) to this directory. 
This file contains 100 000 reads of length 100 all sampled from the pan-genome of four E. coli strains (i.e., each read is sampled from one of the four strains). Thus, each read will have at least one exact occurrence.
If you want to align these reads in a strain-fixed way using the Pigeonhole scheme up to an edit distance of 3 to our reference pan-genome, run the following command in the `build` folder:
```bash
./nexus -e 3 -ss pigeon ../example/EscherichiaColi4Strains ../example/EscherichiaColi4Strains.reads.fasta
```

After this operation your directory structure will look like this:
``` 
    .
    ├── build
    ├── cmake
    ├── example 
    |   ├── EscherichiaColi4Strains.B.left 
    |   ├── EscherichiaColi4Strains.brt 
    |   ├── EscherichiaColi4Strains.DBG 
    |   ├── EscherichiaColi4Strains.rev.bwt 
    |   ├── EscherichiaColi4Strains.sa.bv.1 
    |   ├── EscherichiaColi4Strains.B.right.128 
    |   ├── EscherichiaColi4Strains.bwt 
    |   ├── EscherichiaColi4Strains.left.map 
    |   ├── EscherichiaColi4Strains.reads.fasta 
    |   ├── EscherichiaColi4Strains.reads.fasta_output.txt 
    |   ├── EscherichiaColi4Strains.right.map.128 
    |   ├── EscherichiaColi4Strains.txt 
    |   ├── EscherichiaColi4Strains.B.right.full.128 
    |   ├── EscherichiaColi4Strains.cct 
    |   ├── EscherichiaColi4Strains.rev.brt 
    |   └── EscherichiaColi4Strains.sa.1
    ├── radixSA64
    ├── search_schemes
    ├── src
    └── sux
```

The results can be found in `EscherichiaColi4Strains.reads.fasta_output.txt`.

To align the reads in a strain-free way, run
```bash
./nexus -sfr -e 3 -ss pigeon ../example/EscherichiaColi4Strains ../example/EscherichiaColi4Strains.reads.fasta
```

Note that strain-free matching requires more runtime and reports a higher number of node paths. 

---
**NOTE**

This second alignment of reads will overwrite the `EscherichiaColi4Strains.reads.fasta_output.txt`. Before running a second time, make sure to back up the original file if you would like to keep it stored.

---

Congratulations! You are now able to use Nexus to align reads to a pan-genome of four E. coli strains!

## Custom Search Schemes

The search scheme can either be one of the hardcoded search schemes present in Nexus or you can provide a custom search scheme. In the `search_schemes` folder a number of search schemes is already present. 

To make your own search scheme you need to create a folder containing at least a file called `name.txt`, which contains the name of your scheme on the first line. 
For every maximum edit/hamming distance a subfolder should be present, which contains at least the file `searches.txt`. In this file the searches of your scheme are written line per line. Each line consists of three space-separated arrays: pi, L and U. Each array is written between curly braces {} and the values are comma-separated.

---
**NOTE**

The pi array should be zero-based! The connectivity property should always be satisfied. The L and U array cannot decrease.

---
### Static Partitioning
If you want to provide optimal static partitioning, you can create a file named `static_partitioning.txt` in the folder of the maximum edit/hamming distance this partitioning is for. This file should contain one line with percentages (values between 0 and 1) separated by spaces. The ith percentage corresponds to the starting position (relative to the size of the pattern) of the (i + 1)th part (again this is zero based). The starting position of the first part is always zero and should **not** be provided.

### Dynamic Partitioning
Similarly, to provide values for dynamic partitioning you can create a file called `dynamic_partitioning.txt`. This file should contain two lines. The first line contains percentages (again between 0 and 1) that correspond to the seeding positions, relative to the size of the pattern, of all parts, except the first and last part. 
The second line should contain space-separated integers corresponding to the weights of each part.

---
**NOTE**

Dynamic partitioning is only possible in the case of strain-fixed matching.

---

### Folder Structure Example
Consider a search scheme which supports maximal edit/hamming distances 1, 2 and 4. For distance 1 no static or dynamic partitioning values are known. For distance 2 only static partitioning values are known and for distance 4 both static and dynamic partitioning values are known. The folder structure of this search scheme should look like this:

``` .
    ├── 1
    |   ├── searches.txt
    ├── 2 
    |   ├── searches.txt
    |   ├── static_partitioning.txt
    ├── 4
    |   ├── dynamic_partitioning.txt
    |   ├── searches.txt
    |   ├── static_partitioning.txt
    └── name.txt
```

### Example 3
Consider the pigeonhole search scheme for maximum edit distance 4. The `searches.txt` file should look like this:

```
{0,1,2,3,4} {0,0,0,0,0} {0,4,4,4,4}
{1,2,3,4,0} {0,0,0,0,0} {0,4,4,4,4}
{2,3,4,1,0} {0,0,0,0,0} {0,4,4,4,4}
{3,4,2,1,0} {0,0,0,0,0} {0,4,4,4,4}
{4,3,2,1,0} {0,0,0,0,0} {0,4,4,4,4}
```

### Example 4
The adapted search schemes based on those by Kucherov can be found in the directories:
`search_schemes/kuch_k+1_adapted` and `search_schemes/kuch_k+2_adapted`.

### Other Examples
In the `search_schemes` folder the hardcoded search schemes of Nexus are available as custom search schemes. 

## Pattern Matching to the Underlying Linear Reference.

As part of the code of Nexus is based on [Columba 1.0](https://github.com/biointec/columba/releases/tag/v1.0), Nexus also implements lossless approximate pattern matching directly to the underlying linear reference. In this case, no node path in the pan-genome graph is found. We refer to this functionality as Nexus's built-in version of Columba.

Just as all of the other functionality, it is assumed here that the linear reference is a pan-genome containing DNA characters `A`, `C`, `G` and `T`; and separation characters `$` and `%`. Both separation characters must be present to garantee correct performance.

Matching a read file to the linear reference text can be done by running the following command from the `build` folder.

```bash
./columba [options] basefilename readfile.[ext]
```

options:

```
 [options]
  -e  --max-ed		maximum edit distance [default = 0]
  -s  --sa-sparseness	suffix array sparseness factor [default = 1]
  -p  --partitioning 	Add flag to do uniform/static/dynamic partitioning [default = dynamic]
  -m   --metric	Add flag to set distance metric (editnaive/editopt/hamming) [default = editopt]
  -ss --search-scheme	Choose the search scheme
  options:
	kuch1	Kucherov k + 1
	kuch2	Kucherov k + 2
	kianfar	 Optimal Kianfar scheme
	manbest	 Manual best improvement for Kianfar scheme (only for ed = 4)
	pigeon	 Pigeonhole scheme
	01*0	01*0 search scheme
	custom	custom search scheme, the next parameter should be a path to the folder containing this search scheme

[ext]
	one of the following: fq, fastq, FASTA, fasta, fa
Following input files are required:
	<base filename>.txt: input text T
	<base filename>.cct: character counts table
	<base filename>.sa.[saSF]: suffix array sample every [saSF] elements
	<base filename>.bwt: BWT of T
	<base filename>.brt: Prefix occurrence table of T
	<base filename>.rev.brt: Prefix occurrence table of the reverse of T


    
```

As a result, a custom output file is written in the folder where your readfile was. This `.txt` file contains for each occurrence the following fields: `identifier`, `Position`, `Length`, `ED` and `reverseComplement`. `identifier` identifies the read, the other fields have been discussed above.

# Visualizing Subgraphs
We also foresee functionality to visualize node paths along with their surrounding neighborhood.

The visualization functionality is implemented such that every connection between two nodes is shown separately. This means that it is possible that multiple edges have the same source and target nodes if this particular combination occurs multiple times in the pan-genome reference. Note that in the `Lossless Approximate Pattern Matching on Pan-genome de Bruijn Graphs` paper, edges are defined differently: each edge must have a different combination of source and target nodes. In other words, all parallel connections are considered together as one edge. We choose not to simplify the visualized subgraphs as such, to provide a more complete overview.

## Prerequisites
### Cytoscape
The visualization files that are generated by this implementation can be opened in Cytoscape. You can download the latest version of Cytoscape [here](https://cytoscape.org/download.html).

### yFiles Layout Algorithms
Additionally, for smooth visualization we recommend installing `yFiles Layout Algorithms` in cytoscape. To do so:
* Select `Apps` from the top menu bar
* Select `App Manager...`
* In the `Install Apps` pane, make sure that `Cytoscape App Store` is selected as `Download Site`
* In the search bar, enter `yFiles Layout Algorithms`
* Click `Install`

### Styles
To import a clear graph style, you can create a styles file using the following command (run from the `build` folder): 

```bash
./createStyles [numberOfStrains]
```

As a parameter, enter the number of strains that your pan-genome graph contains. This way, the edge colors are spread out nicely.

After running this program, a file `PanGenomeSubgraph.xml` is created in the `build` folder. In Cytoscape, do the following:
* Select `File` from the top menu bar
* Hover `Import`
* Select `Styles from File...`
* Select the `PanGenomeSubgraph.xml` file you just generated
* Select `Open`

## Visualizing a Read

It is possible to input a DNA read, which is then matched to the pan-genome graph. The resulting occurrences are then visualized. To do so, execute the following in the `build` directory: 

```
./visualizeRead [options] basefilename read
```

There are pattern matching options, which are identical to the ones described above for matching read files. Additionally, two new visualization options are included.

```
 [Pattern matching options]
 -sfr --strain-free	strain-free matching
 -e   --max-ed		maximum edit distance [default = 0]
 -s   --sa-sparseness	suffix array sparseness factor [default = 1]
 -c   --cp-sparseness	sparseness factor that indicates how many checkpoints must be stored to identify nodes. Use "none" to use no checkpoints. Choose a value that was also used during the building process. [default = 128]
 -f   --filter		filtering type that should be used to filter the occurrences. This option is only valid in case of strain-free matching. Options:
	linear		linear filtering is efficient but does not filter out all redundant occurrences. Additionally, in some exceptional cases, a non-optimal replacement occurrence can be chosen. This is the default option.
	complete	complete filtering leads to a set of occurrences with no redundancy. This option is very slow however and thus not recommended.
 -p   --partitioning		Add flag to do uniform/static/dynamic partitioning. Dynamic partitioning cannot be used with strain-free matching. [default = static]
 -m   --metric		Add flag to set distance metric (editnaive/editopt/hamming) [default = editopt]
 -ss  --search-scheme	Choose the search scheme
  options:
	kuch1	Kucherov k + 1
	kuch2	Kucherov k + 2
	kianfar	Optimal Kianfar scheme
	manbest	Manual best improvement for kianfar scheme (only for ed = 4)
	pigeon	Pigeonhole scheme
	01*0	01*0 search scheme
	naive	naive backtracking
	custom	custom search scheme, the next parameter should be a path to the folder containing this search scheme

 [Visualization options]
 -d   --visualization-depth		Depth of the visualized neighborhood around the paths of interest [default = 3]
 -o   --output-files		Prefix of the output files that will be created during the visualization process [default = basefilename]

Following input files are required:
	<base filename>.txt: input text T
	<base filename>.cct: character counts table
	<base filename>.sa.[saSF]: sparse suffix array, with suffix array sparseness factor [saSF] elements
	<base filename>.sa.bv.[saSF]: bitvector indicating which elements of the suffix array are stored.
	<base filename>.bwt: BWT of T
	<base filename>.rev.bwt: BWT of the reverse of T
	<base filename>.brt: Prefix occurrence table of T
	<base filename>.rev.brt: Prefix occurrence table of the reverse of T
	<base filename>.DBG: variable k and the compressed de Bruijn graph.
	<base filename>.B.left: bitvector B_left for the compressed de Bruijn graph.
	<base filename>.B.right.[cpSF]: bitvector B_right for the compressed de Bruijn graph, with checkpoint sparseness factor [cpSF].
	<base filename>.B.right.full.[cpSF]: bitvector B_right_full for the compressed de Bruijn graph, with checkpoint sparseness factor [cpSF].
	<base filename>.left.map: node identifier mapping corresponding to B_left.
	<base filename>.right.map.[cpSF]: node identifier mapping corresponding to B_right, with checkpoint sparseness factor [cpSF].

```

This procedure outputs two files: `outputfilename_SubgraphOverview.tsv` and `outputfilename_SubgraphEdges.tsv`. 

Note however that if your read does not occur in the graph respecting your search parameters, these two files will be empty except for the header.

`outputfilename_SubgraphOverview.tsv` contains an overview of the occurrences corresponding to the input read. This overview is similar to what is reported when mapping a readfile to the graph. For strain-fixed matching, `SubgraphID`, `Path`, `Strain`, `Position`, `Length` and `ED` are reported. For strain-free matching, `SubgraphID`, `Path`, `DistanceFromLeftEnd`, `Length` and `ED` are reported. Note that only the read itself is matched to the graph, not its reverse complement.

`outputfilename_SubgraphEdges.tsv` contains all edges for all subgraphs and can be used to visualize the subgraphs in Cytoscape. This can be done as follows:
* either drag the `outputfilename_SubgraphEdges.tsv` into the `Network` pane or import it as a `Network from File` via the `File` dropdown in the top menu bar. Check that the columns are interpreted by Cytoscape as follows:
  * EdgeKey = Edge Attribute
  * Source = Source Node
  * OmegaShort = Source Node Attribute
  * OmegaFull = Source Node Attribute
  * PartOfPath = Source Node Attribute
  * Target = Target Node
  * OmegaShort = Target Node Attribute
  * OmegaFull = Target Node Attribute
  * PartOfPath = Target Node Attribute
  * Color = Edge Attribute
* Select the network that is shown, then navigate to `Style` in the left menu bar
* In the styles dropdown menu, choose `PanGenomeSubgraph` which you imported earlier
* Select `Layout` in the top menu
* Select `yFiles Hierarchic Layout`

You should now see a clear set of subgraphs corresponding to your input read. 

### Example 5

Use again the pan-genome from [example 1](##Example-1).

If you want to align and visualize a read (e.g., `CGGCATCCAGGTCGTTAATGATGATAGTTGGTCTGGACATTTTTACTCCATGTCGTCGGTACTGCGAGTGTCGCAGATAAACATACCCAAAAGAAAACCC`) in a strain-fixed way using the Pigeonhole scheme up to an edit distance of 3 to our reference pan-genome, run the following command in the `build` folder:
```bash
./visualizeRead -e 3 -ss pigeon ../example/EscherichiaColi4Strains CGGCATCCAGGTCGTTAATGATGATAGTTGGTCTGGACATTTTTACTCCATGTCGTCGGTACTGCGAGTGTCGCAGATAAACATACCCAAAAGAAAACCC
```

The default neighborhood depth of 3 is used. The results can now be found in `EscherichiaColi4Strains_SubgraphOverview.tsv` and `EscherichiaColi4Strains_SubgraphEdges.tsv`. You can now visualize the subgraph in Cytoscape as described above. This should give the following result:

![Example of subgraph](https://github.com/biointec/nexus/blob/media/images/exampleSubgraph.png?raw=true)

The darker nodes are part of the node path. Lighter nodes are part of the neighborhood.

To align the read in a strain-free way, run
```bash
./visualizeRead -sfr -e 3 -ss pigeon ../example/EscherichiaColi4Strains CGGCATCCAGGTCGTTAATGATGATAGTTGGTCTGGACATTTTTACTCCATGTCGTCGGTACTGCGAGTGTCGCAGATAAACATACCCAAAAGAAAACCC
```

No additional node paths are found for this specific read.

## Visualizing a Node Path

You can also directly visualize a node path of interest. To do so, execute the following in the `build` directory: 

```
./visualizePath [options] basefilename path
```

The path should be a comma-separated list of nodes (no whitespaces).

In this case, there are similar but less options:

```
 [Visualization options]
 -d   --visualization-depth		Depth of the visualized neighborhood around the paths of interest [default = 3]
 -o   --output-files		Prefix of the output files that will be created during the visualization process [default = basefilename]
 -s   --sa-sparseness	suffix array sparseness factor [default = 1]
 -c   --cp-sparseness	sparseness factor that indicates how many checkpoints must be stored to identify nodes. Use "none" to use no checkpoints. Choose a value that was also used during the building process.[default = 128]

Following input files are required:
	<base filename>.txt: input text T
	<base filename>.cct: character counts table
	<base filename>.sa.[saSF]: sparse suffix array, with suffix array sparseness factor [saSF] elements
	<base filename>.sa.bv.[saSF]: bitvector indicating which elements of the suffix array are stored.
	<base filename>.bwt: BWT of T
	<base filename>.rev.bwt: BWT of the reverse of T
	<base filename>.brt: Prefix occurrence table of T
	<base filename>.rev.brt: Prefix occurrence table of the reverse of T
	<base filename>.DBG: variable k and the compressed de Bruijn graph.
	<base filename>.B.left: bitvector B_left for the compressed de Bruijn graph.
	<base filename>.B.right.[cpSF]: bitvector B_right for the compressed de Bruijn graph, with checkpoint sparseness factor [cpSF].
	<base filename>.B.right.full.[cpSF]: bitvector B_right_full for the compressed de Bruijn graph, with checkpoint sparseness factor [cpSF].
	<base filename>.left.map: node identifier mapping corresponding to B_left.
	<base filename>.right.map.[cpSF]: node identifier mapping corresponding to B_right, with checkpoint sparseness factor [cpSF].

```

This time, only one output file is generated: `outputfilename_SubgraphEdges.tsv`. The contents are analogous to what was described above. This file can also be visualized in the same way as described above.

### Example 6

Use again the pan-genome from [example 1](##Example-1).

We now want to visualize the subgraph around node path `{551,73827}` (which was one of the occurrences from Example 5) with the default neighborhood depth of 3. To do so, run the following command in the `build` folder:
```bash
./visualizePath ../example/EscherichiaColi4Strains 551,73827
```

The result can now be found in `EscherichiaColi4Strains_SubgraphEdges.tsv`. You can now visualize the subgraph in Cytoscape as described above. This should show one of the subgraphs from the image above.
# Reproducing Results
The results from our paper 'Lossless Approximate Pattern Matching on Pan-genome de Bruijn Graphs'<!-- [Lossless Approximate Pattern Matching on Pan-genome de Bruijn Graphs](TODO)--> can be reproduced by using the following instructions.

First download and compile the source code of [Nexus 1.0.0](https://github.com/biointec/nexus/releases/tag/v1.0.0).

## Dataset

### Reference
We used two human genome strains (GRCh37 and GRCh38) to build a reference pan-genome. From each strain we removed all N's and substituted them with a random nucleotide.
Then we removed the first line, concatenated the chromosomes, removed newline characters, and added a separation character to the end of file for each strain. 
Next, the strains can be concatenated into one pan-genome.
This can be done by executing the following commands:
```bash
tail -n +2 hs.grch37.fasta.non | sed "s/>.*//g" | cat - percent | tr -d '\n' > hs.grch37.txt
tail -n +2 hs.grch38.fasta.non | sed "s/>.*//g" | cat - dollar | tr -d '\n' > hs.grch38.txt
cat hs.grch37.txt hs.grch38.txt > hs.pangenome.txt
```
Where `hs.grch37.fasta.non` and `hs.grch38.fasta.non` are the result of the substitution of the N's for both strains and `hs.pangenome.txt` is the reference pan-genome used to build the index. 
`hs.grch37.txt` and `hs.grch38.txt` are intermediate files which can be removed once the pan-genome is obtained.
`dollar` and `percent` are text files containing a single character `$` and `%`, respectively. These files should be present in your directory before executing the above commands. 

For this reference pan-genome the implicit representation of the compressed de Bruijn graph needs to be built according to the build instructions described above.

### Reads
We sampled 100 000 reads from [an Illumina experiment dataset](ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz), which only contain `A`, `C`, `G` and `T` characters. One way to do this is to select the first 100 000 such reads. For example, using the following python3 script

```python
import sys

read_file = sys.argv[1]
write_file = sys.argv[2]
f = open(read_file, "r")
w = open(write_file, "w")

reads_written = 0

while(reads_written < 100000):
    fastq_record = [next(f) for x in range(4)]
    read = fastq_record[1]
    if(read.count('A') + read.count('C') + read.count('G') + read.count('T') != len(read.strip())):
        # contains a non-ACGT character
        continue
    for line in fastq_record:
        w.write(line)
    reads_written += 1

f.close()
w.close()


```

This python script takes an input name and requested output name. It samples 100 000 reads and writes them to a fastq file [output].fastq, it is assumed that the input file is a fastq file. Of course, sampling of the reads can be achieved in a myriad of ways.

The sampled dataset we used for the results in the paper is available [here](https://github.com/biointec/nexus/releases/download/v1.0.0/sampled_illumina_reads.fastq).

##  Strain-Fixed versus Strain-Free Pattern Matching

---
**NOTE**

These benchmarks are run for checkpoint sparseness factor 64. Make sure that during the build process of the pan-genome data structure, you include this sparseness factor as an option.

---
To reproduce the results of matching the reads to the linear reference text from Table 6, run Nexus's built-in version of Columba with maximal edit distance `[k]` (from the `build` folder) as follows:

Naive backtracking: 

```bash
./columba -e [k] -s 1 -p static -m editopt -ss naive basefile readsfile
```

Pigeonhole principle: 

```bash
./columba -e [k] -s 1 -p static -m editopt -ss custom ../search_schemes/pigeon/ basefile readsfile
```

Kucherov K + 1: 

```bash
./columba -e [k] -s 1 -p static -m editopt -ss custom ../search_schemes/kuch_k+1_adapted/ basefile readsfile
```

To reproduce the strain-fixed results from Table 6, run Nexus with maximal edit distance `[k]` (from the `build` folder) as follows:

Naive backtracking: 

```bash
./nexus -e [k] -s 1 -c 64 -p static -m editopt -ss naive basefile readsfile
```

Pigeonhole principle: 

```bash
./nexus -e [k] -s 1 -c 64 -p static -m editopt -ss custom ../search_schemes/pigeon/ basefile readsfile
```

Kucherov K + 1: 

```bash
./nexus -e [k] -s 1 -c 64 -p static -m editopt -ss custom ../search_schemes/kuch_k+1_adapted/ basefile readsfile
```
To reproduce the strain-free results from Table 6, run Nexus with maximal edit distance `[k]` (from the `build` folder) as follows:

Naive backtracking: 

```bash
./nexus -sfr -e [k] -s 1 -c 64 -f linear -p static -m editopt -ss naive basefile readsfile
```

Pigeonhole principle: 

```bash
./nexus -sfr -e [k] -s 1 -c 64 -f linear -p static -m editopt -ss custom ../search_schemes/pigeon/ basefile readsfile
```

Kucherov K + 1: 

```bash
./nexus -sfr -e [k] -s 1 -c 64 -f linear -p static -m editopt -ss custom ../search_schemes/kuch_k+1_adapted/ basefile readsfile
```

Runtime is reported as `Total duration`.
Search space size is reported as `Average no. index nodes`.
R/U can be calculated as (`Total no. reported node paths`)/(`Total no. unique node paths`).

##  Checkpoint Sparseness Space-Time Tradeoff

---
**NOTE**

These benchmarks are run for multiple different checkpoint sparseness factors. Make sure that during the build process of the pan-genome data structure, you include these sparseness factors as an option.

---

To reproduce the strain-free results from Figure 3, run Nexus with checkpoint sparseness factor `[cp]` (from the `build` folder) as follows:

```bash
./nexus -sfr -e 4 -s 1 -c [cp] -f linear -p static -m editopt -ss ../search_schemes/kuch_k+1_adapted/ basefile readsfile
```

Do this for `[cp] = {8,16,...,16384 and none}`.

Runtime is reported as `Total duration`.
The memory usage of node identifier mapping `IDmapping_right` can be found by inspecting the `basefile.right.map.[cp]` file size.