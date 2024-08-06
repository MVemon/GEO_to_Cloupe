# GEO_to_Cloupe

# Installation

Downloading directly through github

```
install.packages("devtools")
library(devtools)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")

install.packages("hdf5r")

remotes::install_github("10xGenomics/loupeR")

install_github("MVemon/GEO_to_Cloupe")
```



# Usage

`GEO_to_Cloupe()` function is used to directly download scRNA files from GEO databases and automatically converts them into cloupe files, which can be read through the 10X loupe browser. Can convert MTX and H5 files.

General pipeline involves converting the MTX and H5 scRNA files into Seurat object, followed by normalization and clustering, then lastly converted into a cloupe.

**Inputs/Parameters**

GEO_ID_List: (Required) Either a single GEO accession number (eg. "GSE123456") or a list of them (eg. c("GSE1111111", "GSE2222222")).

Downloaded: (default FALSE) Allows for processing scRNA files into cloupe if the files had been previously downloaded. Setting TRUE skips the downloading portion of the code.

Integrate: (default FALSE) Toggles whether you want to utilize Seurat's integrate function while merging multiple files together. This is a time consuming process that may or may not improve clustering.

Resolution: (default 0.5) Sets the sensitivity for clustering the scRNA UMAPs. Values between 0.0 to 1.0, with lower values indicating lower sensitivity (less number of clusters created, larger cluster sizes).

Merge: (default TRUE) Controls whether the different scRNA files are to be merged together and outputed as a single cloupe file, or create separate cloupe files for each dataset.
