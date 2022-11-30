# BioPlexAnalysis
Applications and downstream analysis of BioPlex PPI data

The [BioPlex project](https://bioplex.hms.harvard.edu/) uses affinity-purification mass spectrometry to profile protein-protein interactions (PPIs) in human cell lines.

To date, the BioPlex project has created two proteome-scale, cell-line-specific PPI networks. The first, BioPlex 3.0, results from affinity purification of 10,128 human proteins —- half the proteome —- in 293T cells and includes 118,162 interactions among 14,586 proteins. The second results from 5,522 immunoprecipitations in HCT116 cells and includes 70,966 interactions between 10,531 proteins.

For more information, please see:

* Huttlin et al. [The BioPlex network: a systematic exploration of the human interactome](https://doi.org/10.1016/j.cell.2015.06.043). *Cell*, 2015.
* Huttlin et al. [Architecture of the human interactome defines protein communities and disease networks](https://doi.org/10.1038/nature22366), *Nature*, 2017.
* Huttlin et al. [Dual proteome-scale networks reveal cell-specific remodeling of the human interactome](https://doi.org/10.1016/j.cell.2021.04.011), *Cell*, 2021.

The [BioPlex R package](https://github.com/ccb-hms/BioPlex)
implements access to the BioPlex protein-protein interaction networks and
related resources from within R. 
Besides protein-protein interaction networks for 293T and HCT116 cells,
this includes access to [CORUM](http://mips.helmholtz-muenchen.de/corum)
protein complex data, and transcriptome and proteome data for the two cell lines. 
             
Functionality focuses on importing these data resources and
storing them in dedicated Bioconductor data structures, as a foundation for
integrative downstream analysis of the data. This repository contains code for 
a set of downstream analyses and applications of BioPlexPPI data, bundled in the 
[BioPlexAnalysis R package](https://github.com/ccb-hms/BioPlexAnalysis)
and
[compiled analysis vignettes](https://ccb-hms.github.io/BioPlexAnalysis/).
