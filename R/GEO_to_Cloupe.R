#' GEO scRNA Downloader and Cloupe Converter
#'
#' Automatically downloads scRNA files from GEO bases and converts them into cloupe file format
#' @param GEO_ID_List A list of strings that contain GEO accession codes
#' @param Downloaded A TRUE or FALSE that dictates whether the files have been downloaded prior. Useful for debugging and testing codes without having to re-download files.
#' @param Integrate A TRUE or FALSE that dictates whether to integrate the different scRNA datasets together
#' @param Resolution Set the value for UMAP clustering sensitivity from a value of 0.0 to 1.0 from less to more clusters
#' @param Merge A TRUE or FALSE that dictates whether you want to merge the multiple scRNA data 
#' @param Mitochondria Set threshold for mitochondria QC metrics. Values set between 0 to 100.
#' @return Returns a merged Seurat object of all the different scRNA datasets downloaded, while also exporting cloupe files into the directory "CloupeFiles"
#' @examples 
#' merged_Seurat <- GEO_to_Cloupe(c("GSE213338"));
#' merged_Suerat <- GEO_to_Cloupe(c("GSE213338", "GSE162733"), Downloaded = TRUE, Merge = FALSE);
#' @export

GEO_to_Cloupe <- function(GEO_ID_List, Downloaded = FALSE, Integrate = FALSE, Resolution = 0.5, Merge = TRUE, Mitochondria = 20){
  
  library(GEOquery)
  library(stringr)
  library(Seurat)
  library(loupeR)
  
  options(timeout = max(300, getOption("timeout")))

  Seurat_merged_list <- c()
  
  Seurat_separate_list <- c()
  
  
  for(GEO_ID in GEO_ID_List){
    
    Empty_samples <- TRUE
    
    gse <- getGEO(GEO_ID, GSEMatrix = TRUE)
    
    GEO_accession_list <- c()
    
    checkSamplesList <- getGEOSuppFiles(GEO_ID, fetch_files = FALSE)$fname
    
    for (suppFile in checkSamplesList){

      if(grepl("GSE", suppFile) && grepl(".tar", suppFile)){
        
        Empty_samples = FALSE
        
      }
      
    }
    
    if (!Empty_samples){
      
      file_type_check = getGEOSuppFiles(gse[[1]]@phenoData@data[["geo_accession"]][1], fetch_files = FALSE)$fname
      
      for(z in file_type_check){
        if(grepl("mtx.gz", z)){
          File_Format = "mtx.gz"
        } else if (grepl(".h5", z)){
          File_Format = "h5"
        }
      }
      
      
    } else if (Empty_samples){
      file_type_check = getGEOSuppFiles(GEO_ID, fetch_files = FALSE)$fname
      
      for(z in file_type_check){
        if(grepl("mtx.gz", z)){
          File_Format = "mtx.gz"
        } else if (grepl(".h5", z)){
          File_Format = "h5"
        }
      }
      
    }
    
    organism <- gse[[1]]@phenoData@data[["organism_ch1"]][1]
    
    if(organism == "Homo sapiens"){
      mitochondria <- "^MT-"
    } else if (organism == "Mus musculus"){
      mitochondria <- "^mt-"
    }

    if (!Empty_samples){
      
      for(i in 1: length(gse)){
        for (n in 1: length(gse[[i]]@phenoData@data[["geo_accession"]])){
          
          print(paste0("Checking Samples ",n, "/", length(gse[[i]]@phenoData@data[["geo_accession"]])))
          
          if (grepl(File_Format,getGEOSuppFiles(gse[[i]]@phenoData@data[["geo_accession"]][n], fetch_files = FALSE)[1])){
            GEO_accession_list <- append(GEO_accession_list, gse[[i]]@phenoData@data[["geo_accession"]][n])
          }
          
        }
      }
      
    } else if (Empty_samples){
      
      GEO_accession_list <- append(GEO_accession_list, GEO_ID)
      
    }
    
    Seurat_list <- c()
    
    for(GEO_accession in GEO_accession_list){
      
      print(paste0("Accessing files from ",GEO_accession))
      
      if (Downloaded) {
        
        fileNames = paste0(GEO_accession, "/",getGEOSuppFiles(GEO_accession, fetch_files = FALSE)$fname)
      
        } else {
        
        filePaths = getGEOSuppFiles(GEO_accession)
        
        fileNames = row.names(filePaths) 
      }
      
      fixedNames <- str_split_i(fileNames, "_", -1)
      
      premetaNames <- str_split_i(fileNames[1], "/",-1)
      
      metaNames <- str_split(premetaNames, "_")
      
      metaNames <- metaNames[[1]]
      
      metaNames_filtered <- metaNames[0:(length(metaNames)-1)]
      
      if (File_Format == "mtx.gz"){
        
        for(mtx_name in fileNames){
          
          if (grepl("mtx.gz", mtx_name)){
            
            mtx_file_name <-  mtx_name
            
            prefixNames <- str_split_i(mtx_file_name, "matrix.mtx.gz", 1)
            
            prefixNames <- str_split_i(prefixNames, "/", -1)
            
          }
          
        }
        
        fixedNames <- sub(pattern = "genes.tsv.gz", replacement = "features.tsv.gz", x = fixedNames)
        
        for (x in 1:length(fixedNames)){
          
          if (grepl(("matrix"), fixedNames[x])){
            fixedNames[x] <- "matrix.mtx.gz"
          }
          
          else if (grepl(("barcodes"), fixedNames[x])){
            fixedNames[x] <- "barcodes.tsv.gz"
            
          }
          
          else if (grepl(("features"), fixedNames[x])){
            fixedNames[x] <- "features.tsv.gz"
          }
          
        }
        
        
        file.rename(fileNames, paste0(GEO_accession, "/", fixedNames))
        
        print(paste0("Creating Seurat Object for ", GEO_accession))
        
        Seurat_Object <- Read10X(
          GEO_accession,
          gene.column = 2,
          cell.column = 1,
          unique.features = TRUE,
          strip.suffix = FALSE
        )
        
      } else if ((File_Format == "h5")){
        
        for(h5_name in fileNames){
          if (length(fileNames) == 1){
            
            h5_file_name <-  h5_name
            
            prefixNames <- str_split_i(h5_file_name, ".h5", 1)
            prefixNames <- str_split_i(prefixNames, "/", -1)
            
            
          } else if (grepl("filtered_feature_bc_matrix.h5", h5_name)){
            
            h5_file_name <-  h5_name
            
            prefixNames <- str_split_i(h5_file_name, "filtered_feature", 1)
            prefixNames <- str_split_i(prefixNames, "/", -1)
              
          } else if (grepl(".h5", h5_name)){
            
            h5_file_name <-  h5_name
            
            prefixNames <- str_split_i(h5_file_name, ".h5", -1)
            
            
          }
 
        }
        
        print(paste0("Creating Seurat Object for ", GEO_accession))
        
        Seurat_Object <- Read10X_h5(h5_file_name, use.names = TRUE)
        
      }
      
      # else if ((File_Format == "csv")){
      #   
      #   Seurat_Object <- read.csv(file = fileNames, header = TRUE, row.names = 1)
      #   Seurat_Object <- t(Seurat_Object)
      #   
      # }
      # 
      
      Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = GEO_ID, min.features = 50)
      
      
      for(i in 1: length(metaNames_filtered)){
        
        colName = paste("variable",i)
        Seurat_Object <- AddMetaData(Seurat_Object, metaNames_filtered[i], col.name = colName)
        
      }
      
      #Seurat_list <- append(Seurat_list, assign(GEO_accession, Seurat_Object))
      Seurat_list <- append(Seurat_list, Seurat_Object)
      
      
      # if(!Empty_samples){
      #   
      #   rm(list = ls(pattern = "^GSM"))
      #   
      # } else if(Empty_samples){
      #   
      #   rm(list = ls(pattern = "^GSE"))
      #   
      # }
      
      if (!Merge){
        
        print(paste0("Starting processing of ", GEO_accession, " Seurat object"))
        
        Seurat_Object[["percent.mt"]] <- PercentageFeatureSet(Seurat_Object, pattern = mitochondria)
        
        Seurat_Object <- subset(Seurat_Object, subset = percent.mt < Mitochondria)
        
        Seurat_Object <- NormalizeData(Seurat_Object, normalization.method = "LogNormalize", scale.factor = 10000)
        
        Seurat_Object <- FindVariableFeatures(Seurat_Object, selection.method = "vst", nfeatures = 2000)
        
        all.genes <- rownames(Seurat_Object)
        
        Seurat_Object <- ScaleData(Seurat_Object, features = all.genes)
        
        Seurat_Object <- RunPCA(Seurat_Object, features = VariableFeatures(object = Seurat_Object))
        
        Seurat_Object <- FindNeighbors(Seurat_Object, dims = 1:30)
        
        Seurat_Object <- FindClusters(Seurat_Object, resolution = Resolution)
        
        Seurat_Object <- RunUMAP(Seurat_Object, dims = 1:30)
        
        Seurat_Object <- JoinLayers(Seurat_Object)
        
        #Seurat_separate_list <- append(Seurat_separate_list, assign(paste0(GEO_accession), Seurat_Object))
        Seurat_separate_list <- append(Seurat_separate_list, Seurat_Object)
        
        
        dir.create("CloupeFiles", showWarnings = FALSE)
        
        print(paste0("Creating loupe file for ", GEO_accession))
        
        create_loupe_from_seurat(
          Seurat_Object,
          output_dir = "CloupeFiles",
          output_name = "TEST_MITO",
          dedup_clusters = FALSE,
          executable_path = NULL,
          force = TRUE)
        
      }
      
    }
    
    if (Merge){
      
      print(paste0("Starting merging and processing of ", GEO_ID))
      
      Seurat_merged <- merge(Seurat_Object, y = Seurat_list[1:(length(Seurat_list)-1)], project = GEO_ID)
      
      Seurat_merged[["percent.mt"]] <- PercentageFeatureSet(Seurat_merged, pattern = mitochondria)
      
      Seurat_merged <- subset(Seurat_merged, subset = percent.mt < Mitochondria)
      
      Seurat_merged <- NormalizeData(Seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
      
      Seurat_merged <- FindVariableFeatures(Seurat_merged, selection.method = "vst", nfeatures = 2000)
      
      all.genes <- rownames(Seurat_merged)
      
      Seurat_merged <- ScaleData(Seurat_merged, features = all.genes)
      
      Seurat_merged <- RunPCA(Seurat_merged, features = VariableFeatures(object = Seurat_merged))
      
      if (Integrate){
        
        print("Starting integration")
        
        Seurat_merged <- IntegrateLayers(object = Seurat_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
        
        Seurat_merged[["RNA"]] <- JoinLayers(Seurat_merged[["RNA"]])
        
        Seurat_merged <- FindNeighbors(Seurat_merged, reduction = "integrated.cca", dims = 1:30)
        
        Seurat_merged <- FindClusters(Seurat_merged, resolution = Resolution)
        
        Seurat_merged <- RunUMAP(Seurat_merged, dims = 1:30, reduction = "integrated.cca")
        
        dir.create("CloupeFiles", showWarnings = FALSE)
        
        print(paste0("Creating loupe file for ", GEO_ID))
        
        create_loupe_from_seurat(
          Seurat_merged,
          output_dir = "CloupeFiles",
          output_name = paste0(GEO_ID,".Integrated.",File_Format),
          dedup_clusters = FALSE,
          executable_path = NULL,
          force = TRUE)
        
      } else {
        
        Seurat_merged <- FindNeighbors(Seurat_merged, dims = 1:30)
        
        Seurat_merged <- FindClusters(Seurat_merged, resolution = Resolution)
        
        Seurat_merged <- RunUMAP(Seurat_merged, dims = 1:30)
        
        Seurat_merged <- JoinLayers(Seurat_merged)
        
        dir.create("CloupeFiles", showWarnings = FALSE)
        
        print(paste0("Creating loupe file for ", GEO_ID))
        
        create_loupe_from_seurat(
          Seurat_merged,
          output_dir = "CloupeFiles",
          output_name = paste0(GEO_ID,".",File_Format),
          dedup_clusters = FALSE,
          executable_path = NULL,
          force = TRUE)
        
      }
      
      # Seurat_merged_list <- append(Seurat_merged_list, assign(paste0(GEO_accession, "_merged"), Seurat_merged))
      Seurat_merged_list <- append(Seurat_merged_list, Seurat_merged)
      
    }
    
  }
  
  if(Merge){
    
    return(Seurat_merged_list)

  } else {
    
    return(Seurat_separate_list)
    
  }
  
}
