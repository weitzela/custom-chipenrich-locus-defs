# custom-chipenrich-locus-defs  
This repository is designed to create custom locus definitions to use with ChipEnrich analysis, following the methods outlined in [[Welch, et al., 2014]](#1). It offers flexibiliity to tailor locus definitions across various genome releases and adjust the widths surrounding the TSS. **Disclaimer**: The rn7 genome was the primary focus during development and efforts were made to allow adaptability for other genomes. However, my limited experience in working with other genomes may present unexpected difficulties.  

**Required packages:**  
- tidyverse
- GenomicRanges
- [gUtils](https://github.com/mskilab-org/gUtils)
- plyranges
- biomaRt
- janitor
- chipenrich

**Example chipenrich command to use the custom locus definition:**  
```
chipenrich(
    peaks = peak_df, 
    locusdef = "rn7-locus_intergenic.txt", 
    genome = "rn6", 
    genesets = "GOBP"
)
```
*Note*: You must set the genome argument to one that exists within the chipenirch package itself. In this example, setting the argument to rn6 will return inaccurate TSS-related information in the resulting *peaks* dataframe because the `assign_peaks` function within the chipenrich package refers to built in rn6-TSS data. The pathway results do not face this problem.  

**Recreate Files:**  
```
Rscript R/create-locus-defs.R 
```

**Tips to Curate the Repository for Personal Use**
- Add a new `genome-info/[genome-version]_seq-info.RData` if you would like to define locus definitions for genes outside of chromosomes 1-20, or create locus definitions for a different genome version  
- Change the `genome_id` definition in the `R/create-locus-defs.R` file to match the one used in the seq-info filename  
- Update the list of possible entrez IDs, based on the pathway databases you plan to probe  

## References  
<a id="1">[1]</a> 
Welch RP, Lee C, Imbriano PM, Patil S, Weymouth TE, Smith RA, et al. ChIP-Enrich: gene set enrichment testing for ChIP-seq data. Nucleic Acids Res 2014;42:e105–e105. https://doi.org/10.1093/nar/gku463.
