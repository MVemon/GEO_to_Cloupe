#' GEO scRNA Downloader and Cloupe Converter
#'
#' Automatically downloads scRNA files from GEO bases and converts them into cloupe file format
#' @param GEO_ID_List A list of strings that contain GEO accession codes
#' @param File_Format A string of the file format ending type
#' @return Returns a merged Seurat object of all the different scRNA datasets downloaded, while also exporting cloupe files into the directory "CloupeFiles"
#' @examples 
#' merged_Seurat <- GEO_to_Cloupe(c("GSE108677"), ".h5");
#' merged_Suerat <- (c("GSE108677", "GSE128033"), "mtx.gz" );
#' @export

GEO_to_Cloupe <- function(GEO_ID_List, File_Format){
  
  library(GEOquery)
  library(stringr)
  library(Seurat)
  library(loupeR)
  
  options(timeout = max(300, getOption("timeout")))

  Seurat_merged_list <- c()
  
  for (GEO_ID in GEO_ID_List){
    gse <- getGEO(GEO_ID, GSEMatrix = TRUE)
    
    GEO_accession_list <- c()
    
    for(i in 1: length(gse)){
      print(i)
      for (n in 1: length(gse[[i]]@phenoData@data[["geo_accession"]])){
        print(n)
        if (grepl(File_Format,getGEOSuppFiles(gse[[i]]@phenoData@data[["geo_accession"]][n], fetch_files = FALSE)[1])){
          GEO_accession_list <- append(GEO_accession_list, gse[[i]]@phenoData@data[["geo_accession"]][n])
        }
      }
    }
    
    Seurat_list <- c()
    
    for(GEO_accession in GEO_accession_list){
      print(GEO_accession)
      
      filePaths = getGEOSuppFiles(GEO_accession)
      
      fileNames = row.names(filePaths)
      
      #fileNames = getGEOSuppFiles(GEO_accession, fetch_files = FALSE)$fname #Use this for troubleshooting without having to download files again
      
      fixedNames <- str_split_i(fileNames, "_", -1)
      
      premetaNames <- str_split_i(fileNames[1], "/",-1)
      
      metaNames <- str_split(premetaNames, "_")
      
      #metaNames <- str_split(metaNames[1], end = -2)
      
      metaNames <- metaNames[[1]]
      
      metaNames_filtered <- metaNames[0:(length(metaNames)-1)]
      
      
      if (File_Format == "mtx.gz"){
        
        
        
        fixedNames <- sub(pattern = "genes.tsv.gz", replacement = "features.tsv.gz", x = fixedNames)
        
        file.rename(fileNames, paste0(GEO_accession, "/", fixedNames))
        
        Seurat_Object <- Read10X(
          GEO_accession,
          gene.column = 2,
          cell.column = 1,
          unique.features = TRUE,
          strip.suffix = FALSE
        )
        
      }
      
      else if ((File_Format == "h5")){
        
        Seurat_Object <- Read10X_h5(fileNames, use.names = TRUE)
        
      }
      
      
      Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = GEO_ID, min.features = 200)
      
      
      for(i in 1: length(metaNames_filtered)){
        
        colName = paste("variable",i)
        Seurat_Object <- AddMetaData(Seurat_Object, metaNames_filtered[i], col.name = colName)
        
      }
      
      Seurat_list <- append(Seurat_list, assign(GEO_accession, Seurat_Object))
      
    }
    
    Seurat_merged <- merge(Seurat_Object, y = Seurat_list[1:(length(Seurat_list)-1)], project = GEO_ID)
    
    Seurat_merged <- NormalizeData(Seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
    
    Seurat_merged <- FindVariableFeatures(Seurat_merged, selection.method = "vst", nfeatures = 2000)
    
    all.genes <- rownames(Seurat_merged)
    
    Seurat_merged <- ScaleData(Seurat_merged, features = all.genes)
    
    Seurat_merged <- RunPCA(Seurat_merged, features = VariableFeatures(object = Seurat_merged))
    
    Seurat_merged <- FindNeighbors(Seurat_merged, dims = 1:30)
    
    Seurat_merged <- FindClusters(Seurat_merged, resolution = 0.5)
    
    Seurat_merged <- RunUMAP(Seurat_merged, dims = 1:30)
    
    Seurat_merged_joined <- JoinLayers(Seurat_merged)

    dir.create("CloupeFiles")
    
    create_loupe_from_seurat(
      Seurat_merged_joined,
      output_dir = "CloupeFiles",
      output_name = paste0(GEO_ID,".",File_Format),
      dedup_clusters = FALSE,
      executable_path = NULL,
      force = TRUE)
    
    Seurat_merged_list <- append(Seurat_merged_list, assign(paste0(GEO_accession, "_merged"), Seurat_merged))

  }
  
  return(Seurat_merged_joined)
}

