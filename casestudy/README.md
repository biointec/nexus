# Case Study

In this directory, we provide a case study that demonstrates how visualizing subgraphs can aid in extracting useful information from the pan-genome graph topology. Specifically, we study compensatory mutations corresponding to RRDR mutations causing antibiotic resistance. The pipeline implemented here is more of an ad hoc solution to a specific problem, rather than a general pipeline to be readily applied to other problems. Nevertheless, if there were any interest to extend these ideas into a more general pipeline, do not hesitate to contact the authors of this page: 
* Lore.Depuydt@Ugent.be 
* Jan.Fostier@Ugent.be

For the case study, a pan-genome of 341 M. tuberculosis strains must be built, which can be done as follows:

```bash
./nexusBuild -s 16 -c 128 -p MTuberculosisPanGenome 19
```
The preprocessed reference text can be downloaded [here](https://github.com/biointec/nexus/releases/download/v1.1.0/MTuberculosisPanGenome.txt). Annotation for the strains in the pan-genome can be found [here](https://github.com/biointec/nexus/releases/download/v1.1.0/MTuberculosisPanGenome.annotation.txt).

The case study can by run by executing the following command **inside the folder containing the index that was built earlier**:

```bash
./casestudy
```
As a result, a file `MTuberculosisPanGenome_Compensatory.tsv` is created that contains a set of candidate putative compensatory mutations corresponding to three RRDR mutations. 

