# Bacterial pan-genome data

## References

We obtained all complete genomes corresponding to eight species from the NCBI database:
* Escherichia coli: 3472 genomes
* Salmonella enterica: 1569 genomes
* Listeria monocytogenes: 335 genomes
* Pseudomonas aeruginosa: 863 genomes
* Bacillus subtilis: 333 genomes
* Limosilactobacillus fermentum: 62 genomes
* Enterococcus faecalis: 269 genomes
* Staphylococcus aureus: 1252 genomes

This was done using the following command:

```bash
ncbi-genome-download -s refseq -F fasta -l complete -g "Escherichia Coli,Salmonella enterica,Listeria monocytogenes,Pseudomonas aeruginosa,Bacillus subtilis,Limosilactobacillus fermentum,Enterococcus faecalis,Staphylococcus aureus" bacteria
```

A similar dataset was used to generate the [SPUMONI 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02958-1) results.

The resulting files were shuffled, generating the list available [here](metadata.tsv).

To create a pan-genome with X bacterial strains, we simply select the first X strains from this list.