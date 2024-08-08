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

# Package Version

Confirmed to work on R 4.3.2 and R 4.4.1 on both Windows and Mac.

GEOquery (>= 2.68.0)

stringr (>= 1.5.1)

Seurat (>= 5.0.0)

loupeR (>= 1.0.0)

hdf5r (>= 1.3.8)

devtools (>= 2.4.5)

BiocManager (>= 1.30.22)


# Usage

```
library(devtools)

library(GEOtoCloupe)

Seurat_lists <- GEO_to_Cloupe("GSE210187")
```

`GEO_to_Cloupe()` function is used to directly download scRNA files from GEO databases and automatically converts them into cloupe files, which can be read through the 10X loupe browser. Can convert MTX and H5 files.

General pipeline involves converting the MTX and H5 scRNA files into Seurat object, followed by normalization and clustering, then lastly converted into a cloupe.

**Inputs/Parameters**

GEO_ID_List: (Required) Either a single GEO accession number (eg. "GSE123456") or a list of them (eg. c("GSE1111111", "GSE2222222")).

```
Seurat_lists <- GEO_to_Cloupe(c("GSE110037", "GSE180416"))
```

Downloaded: (default FALSE) Allows for processing scRNA files into cloupe if the files had been previously downloaded. Setting TRUE skips the downloading portion of the code.

```
Seurat_lists <- GEO_to_Cloupe("GSE110037", Downloaded = TRUE)
```

Integrate: (default FALSE) Toggles whether you want to utilize Seurat's integrate function while merging multiple files together. This is a time consuming process that may or may not improve clustering.

```
Seurat_lists <- GEO_to_Cloupe("GSE110037", Integrate = TRUE)
```

Resolution: (default 0.1) Sets the sensitivity for clustering the scRNA UMAPs. Values between 0.0 to 1.0, with lower values indicating lower sensitivity (less number of clusters created, larger cluster sizes).

```
Seurat_lists <- GEO_to_Cloupe("GSE110037", Resolution = 0.5)
```

Merge: (default TRUE) Controls whether the different scRNA files are to be merged together and outputed as a single cloupe file, or create separate cloupe files for each dataset.

```
Seurat_lists <- GEO_to_Cloupe("GSE110037", Merge = FALSE)
```

Mitochondria: (default 30) Sets the threshold for mitochondrial QC metrics. From a scale of 0 to 100.

```
Seurat_lists <- GEO_to_Cloupe("GSE110037", Mitochondria = 15)
```

Combined Example:

```
Seurat_lists <- GEO_to_Cloupe("GSE110037", Merged = FALSE, Downloaded = TRUE, Resolution = 0.2, Mitochondria = 10)
```
