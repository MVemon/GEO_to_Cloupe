#' GEO scRNA Downloader and Cloupe Converter
#'
#' Automatically downloads scRNA files from GEO bases and converts them into cloupe file format
#' @param GEO_ID_List A list of strings that contain GEO accession codes
#' @param Downloaded A TRUE or FALSE that dictates whether the files have been downloaded prior. Useful for debugging and testing codes without having to re-download files.
#' @param Integrate A TRUE or FALSE that dictates whether to integrate the different scRNA datasets together
#' @param Resolution Set the value for UMAP clustering sensitivity from a value of 0.0 to 1.0 from less to more clusters
#' @param Merge A TRUE or FALSE that dictates whether you want to merge the multiple scRNA data 
#' @param Mitochondria Set threshold for mitochondria QC metrics. Values set between 0 to 100.
#' @param Dims Sets the dims values during clustering and Seurat processing
#' @param MetaData Set file name containing metadata information.
#' @return Returns a merged Seurat object of all the different scRNA datasets downloaded, while also exporting cloupe files into the directory "CloupeFiles"
#' @examples 
#' merged_Seurat <- GEO_to_Cloupe(c("GSE213338"));
#' merged_Suerat <- GEO_to_Cloupe(c("GSE213338", "GSE162733"), Downloaded = TRUE, Merge = FALSE);
#' @export

GEO_to_Cloupe <- function(GEO_ID_List, Downloaded = FALSE, Integrate = FALSE, Resolution = 0.10, Merge = TRUE, Mitochondria = 30, Dims = 1:30, MetaData = NULL){
  
  library(GEOquery)
  library(stringr)
  library(Seurat)
  library(loupeR)
  
  options(timeout = max(300, getOption("timeout")))
  
  #parallel = FALSE, parallel_workers = 4, maxSize = 850
  
  # if (parallel){
  #   
  #   library(future)
  #   
  #   plan("multicore", workers = parallel_workers)
  #   
  #   
  #   
  #   options(future.globals.maxSize= maxSize*1024^2)
  #   
  # }
  
  
  Seurat_merged_list <- c()
  
  Seurat_separate_list <- c()
  
  
  for(GEO_ID in GEO_ID_List){
    
    print(paste("Processing ",GEO_ID))
    
    Empty_samples <- TRUE
    
    gse <- getGEO(GEO_ID, GSEMatrix = TRUE)
    
    GEO_accession_list <- c()
    
    mouse_list <- c()
    
    human_list <- c()
    
    mixed_organism_dataset <- FALSE
    
    checkSamplesList <- getGEOSuppFiles(GEO_ID, fetch_files = FALSE)$fname
    
    for (suppFile in checkSamplesList){
      
      if(grepl("GSE", suppFile) && grepl(".tar", suppFile)){
        
        Empty_samples = FALSE
        
      }
      
    }
    
    # if (!Empty_samples){
    #   
    #   file_type_check = getGEOSuppFiles(gse[[1]]@phenoData@data[["geo_accession"]][1], fetch_files = FALSE)$fname
    #   
    #   for(z in file_type_check){
    #     if(grepl("mtx.gz", z)){
    #       File_Format = "mtx.gz"
    #     } else if (grepl(".h5", z)){
    #       File_Format = "h5"
    #     } else {
    #       File_Format = "ignore"
    #     }
    #   }
    #   
    #   
    # } else if (Empty_samples){
    #   file_type_check = getGEOSuppFiles(GEO_ID, fetch_files = FALSE)$fname
    #   
    #   for(z in file_type_check){
    #     if(grepl("mtx.gz", z)){
    #       File_Format = "mtx.gz"
    #     } else if (grepl(".h5", z)){
    #       File_Format = "h5"
    #     } else {
    #       File_Format = "ignore"
    #     }
    #   }
    #   
    # }
    
    organism <- gse[[1]]@phenoData@data[["organism_ch1"]][1]
    
    if(organism == "Homo sapiens"){
      mitochondria <- "^MT-"
    } else if (organism == "Mus musculus"){
      mitochondria <- "^mt-"
    }
    
    if (!Empty_samples){
      
      for(i in 1: length(gse)){
        
        organism_set <- gse[[i]]@phenoData@data[["organism_ch1"]][1]
        
        for (n in 1: length(gse[[i]]@phenoData@data[["geo_accession"]])){
          
          file_type_check = getGEOSuppFiles(gse[[i]]@phenoData@data[["geo_accession"]][n], fetch_files = FALSE)$fname
          
          for(z in file_type_check){
            if(grepl("mtx.gz", z)){
              File_Format = "mtx.gz"
            } else if (grepl(".h5", z)){
              File_Format = "h5"
            } else if (grepl("txt.gz", z)){
              File_Format = "txt.gz"
            } else {
              File_Format = "ignore"
            }
          }
          
          print(paste0("Checking Samples ",n, "/", length(gse[[i]]@phenoData@data[["geo_accession"]])))
          
          if (grepl(File_Format,getGEOSuppFiles(gse[[i]]@phenoData@data[["geo_accession"]][n], fetch_files = FALSE)[1])){
            
            if(organism_set == "Homo sapiens"){
              human_list <- append(human_list, gse[[i]]@phenoData@data[["geo_accession"]][n])
            } else if (organism_set == "Mus musculus"){
              mouse_list <- append(mouse_list, gse[[i]]@phenoData@data[["geo_accession"]][n])
              
            }
            
            GEO_accession_list <- append(GEO_accession_list, gse[[i]]@phenoData@data[["geo_accession"]][n])
            
            
          }
          
        }
      }
      
    } else if (Empty_samples){
      
      file_type_check = getGEOSuppFiles(GEO_ID, fetch_files = FALSE)$fname
      
      for(z in file_type_check){
        if(grepl("mtx.gz", z)){
          File_Format = "mtx.gz"
        } else if (grepl(".h5", z)){
          File_Format = "h5"
        } else if (grepl("txt.gz", z)){
          File_Format = "txt.gz"
        } else {
          File_Format = "ignore"
        }
      }
      
      GEO_accession_list <- append(GEO_accession_list, GEO_ID)
      
    }
    
    Seurat_list <- c()
    
    Mouse_Seurat_list <- c()
    
    Human_Seurat_list <- c()
    
    
    if(length(mouse_list) > 0 && length(human_list) > 0) {
      
      #genome_list <- list(mouse_list, human_list)
      
      mixed_organism_dataset = TRUE
      
    }
    
    for(GEO_accession in GEO_accession_list){
      
      print(paste0("Accessing files from ",GEO_accession))
      
      if (Downloaded) {
        
        fileNames = paste0(GEO_accession, "/",getGEOSuppFiles(GEO_accession, fetch_files = FALSE)$fname)
        
      } else {
        
        filePaths = getGEOSuppFiles(GEO_accession)
        
        fileNames = row.names(filePaths) 
      }
      
      #fixedNames <- str_split_i(fileNames, "_", -1)
      
      fixedNames <- str_split_i(fileNames, "/", -1)
      
      premetaNames <- str_split_i(fileNames[1], "/",-1)
      
      metaNames <- str_split(premetaNames, "_")
      
      metaNames <- metaNames[[1]]
      
      metaNames_filtered <- metaNames[0:(length(metaNames)-1)]
      
      file_type_check = getGEOSuppFiles(GEO_accession, fetch_files = FALSE)$fname

        for(z in file_type_check){
          if(grepl("mtx.gz", z)){
            File_Format = "mtx.gz"
          } else if (grepl(".h5", z)){
            File_Format = "h5"
          } else if (grepl("txt.gz", z)){
            File_Format = "txt.gz"
          } else {
            File_Format = "ignore"
          }
        }
      
      if (File_Format == "mtx.gz"){
        
        for(mtx_name in fixedNames){
          
          if (grepl("mtx.gz", mtx_name)){
            
            mtx_file_name <-  mtx_name
            
            prefixNames <- str_split_i(mtx_file_name, "matrix.mtx.gz", 1)
            
            prefixNames <- str_split_i(prefixNames, "/", -1)
            
          }
          
        }
        
        
        for (x in 1:length(fixedNames)){
          
          if (grepl(("matrix"), fileNames[x])){
            
            fixedNames[x] <- "matrix.mtx.gz"
            
          } else if (grepl(("barcodes"), fixedNames[x])){
            
            fixedNames[x] <- "barcodes.tsv.gz"
            
          } else if (grepl(("features"), fixedNames[x])){
            
            fixedNames[x] <- "features.tsv.gz"
            
          } else if (grepl(("genes"), fixedNames[x])){
            
            fixedNames[x] <- "genes.tsv.gz"
            
          }
          
        }
        
        fixedNames <- sub(pattern = "genes.tsv.gz", replacement = "features.tsv.gz", x = fixedNames)
        
        
        file.rename(fileNames, paste0(GEO_accession, "/", fixedNames))
        
        print(paste0("Creating Seurat Object for ", GEO_accession))
        
        df<-read.delim(paste0(GEO_accession, "/features.tsv.gz"),sep="\t")
        
        if (ncol(df) > 1){
          Seurat_Object <- Read10X(
            GEO_accession,
            gene.column = 2,
            cell.column = 1,
            unique.features = TRUE,
            strip.suffix = FALSE
          )
        } else if  (ncol(df ) == 1){
          Seurat_Object <- Read10X(
            GEO_accession,
            gene.column = 1,
            cell.column = 1,
            unique.features = TRUE,
            strip.suffix = FALSE
          )
        }
        
        
        
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
        
      } else if (File_Format == "txt.gz"){
        for(txt_name in fileNames){
          
          if (length(fileNames) == 1){
            
            txt_file_name <-  txt_name
            
            prefixNames <- str_split_i(txt_file_name, ".txt.gz", 1)
            prefixNames <- str_split_i(prefixNames, "/", -1)
            
            
          } else if (grepl("Counts", txt_name)){
            
            
            txt_file_name <-  txt_name
            
            prefixNames <- str_split_i(txt_file_name, "Counts", 1)
            prefixNames <- str_split_i(prefixNames, "/", -1)
            
            break
            
          } else if (grepl("counts", txt_name)){
            
            txt_file_name <-  txt_name
            
            prefixNames <- str_split_i(txt_file_name, "counts", 1)
            prefixNames <- str_split_i(prefixNames, "/", -1)
            
          } else if (grepl(".txt", txt_name)){
            

            txt_file_name <-  txt_name
            
            prefixNames <- str_split_i(txt_file_name, ".txt", -1)
            
            
          }
          
        }
        
        print(paste0("Creating Seurat Object for ", GEO_accession))
        
        
        Seurat_Object <- read.delim(txt_file_name, header = T, stringsAsFactors = F, row.names = 1)
        

      }
      
      # else if ((File_Format == "csv")){
      #   
      #   Seurat_Object <- read.csv(file = fileNames, header = TRUE, row.names = 1)
      #   Seurat_Object <- t(Seurat_Object)
      #   
      # }
      # 
      
      if (length(Seurat_Object) == 2){
        Seurat_Object <- Seurat_Object[["Gene Expression"]]
      }
      
      Seurat_Object <- CreateSeuratObject(counts = Seurat_Object, project = GEO_ID, min.features = 50)
      
      
      
      for(i in 1: length(metaNames_filtered)){
        
        colName = paste("variable",i)
        Seurat_Object <- AddMetaData(Seurat_Object, metaNames_filtered[i], col.name = colName)
        
        
      }
      
      if (!(is.null(MetaData))){
        
        print(paste0("Adding metadata to Seurat object from ", MetaData))
        
        metadata <- read.delim(paste0(GEO_accession,"/",MetaData), header = T, stringsAsFactors = F, row.names = 1)
        
        Seurat_Object <- AddMetaData(Seurat_Object, metadata)
        
      }
      
      #Seurat_list <- append(Seurat_list, assign(GEO_accession, Seurat_Object))
      
      if(mixed_organism_dataset){
        for(mouse_GEO in mouse_list){
          if(grepl(GEO_accession, mouse_GEO)){
            Mouse_Seurat_list <- append(Mouse_Seurat_list, Seurat_Object)
            
          }
        }
        
        for (human_GEO in human_list){
          if(grepl(GEO_accession, human_GEO)){
            Human_Seurat_list <- append(Human_Seurat_list, Seurat_Object)
            
          }
        }
        
        
      } else {
        
        Seurat_list <- append(Seurat_list, Seurat_Object)
        
      }
      
      
      
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
        
        Seurat_Object <- FindNeighbors(Seurat_Object, dims = Dims)
        
        Seurat_Object <- FindClusters(Seurat_Object, resolution = Resolution)
        
        Seurat_Object <- RunUMAP(Seurat_Object, dims = Dims)
        
        Seurat_Object <- RunTSNE(Seurat_Object, dims = Dims, check_duplicates = FALSE)
        
        
        Seurat_Object <- JoinLayers(Seurat_Object)
        
        #Seurat_separate_list <- append(Seurat_separate_list, assign(paste0(GEO_accession), Seurat_Object))
        Seurat_separate_list <- append(Seurat_separate_list, Seurat_Object)
        
        
        dir.create("CloupeFiles", showWarnings = FALSE)
        
        print(paste0("Creating loupe file for ", GEO_accession))
        
        create_loupe_from_seurat(
          Seurat_Object,
          output_dir = "CloupeFiles",
          output_name = paste0(prefixNames,".",File_Format),
          dedup_clusters = FALSE,
          executable_path = NULL,
          force = TRUE)
        
      }
      
    }
    
    if (Merge){
      
      if (!mixed_organism_dataset){
        
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
          
          Seurat_merged <- FindNeighbors(Seurat_merged, reduction = "integrated.cca", dims = Dims)
          
          Seurat_merged <- FindClusters(Seurat_merged, resolution = Resolution)
          
          Seurat_merged <- RunUMAP(Seurat_merged, dims = Dims, reduction = "integrated.cca")
          
          Seurat_merged <- RunTSNE(Seurat_merged, dims = Dims, reduction = "integrated.cca", check_duplicates = FALSE)
          
          dir.create("CloupeFiles", showWarnings = FALSE)
          
          print(paste0("Creating loupe file for ", GEO_ID))
          
          create_loupe_from_seurat(
            Seurat_merged,
            output_dir = "CloupeFiles",
            output_name = paste0(GEO_ID,".Integrated.",File_Format),
            dedup_clusters = FALSE,
            executable_path = NULL,
            force = TRUE)
          
        } else if (!Integrate){
          
          Seurat_merged <- FindNeighbors(Seurat_merged, dims = Dims)
          
          Seurat_merged <- FindClusters(Seurat_merged, resolution = Resolution)
          
          Seurat_merged <- RunUMAP(Seurat_merged, dims = Dims)
          
          Seurat_merged <- RunTSNE(Seurat_merged, dims = Dims, check_duplicates = FALSE)
          
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
        
      } else if (mixed_organism_dataset){
        
        print(paste0("Starting merging and processing of ", GEO_ID, " mouse data files"))
        
        last_Mouse_Seurat_object <- Mouse_Seurat_list[[length(Mouse_Seurat_list)]]
        
        Mouse_Seurat_merged <- merge(last_Mouse_Seurat_object, y = Mouse_Seurat_list[1:(length(Mouse_Seurat_list)-1)], project = GEO_ID)
        
        Mouse_Seurat_merged[["percent.mt"]] <- PercentageFeatureSet(Mouse_Seurat_merged, pattern = "^mt-")
        
        Mouse_Seurat_merged <- subset(Mouse_Seurat_merged, subset = percent.mt < Mitochondria)
        
        Mouse_Seurat_merged <- NormalizeData(Mouse_Seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
        
        Mouse_Seurat_merged <- FindVariableFeatures(Mouse_Seurat_merged, selection.method = "vst", nfeatures = 2000)
        
        all.genes <- rownames(Mouse_Seurat_merged)
        
        Mouse_Seurat_merged <- ScaleData(Mouse_Seurat_merged, features = all.genes)
        
        Mouse_Seurat_merged <- RunPCA(Mouse_Seurat_merged, features = VariableFeatures(object = Mouse_Seurat_merged))
        
        last_Human_Seurat_object <- Human_Seurat_list[[length(Human_Seurat_list)]]
        
        Human_Seurat_merged <- merge(last_Human_Seurat_object, y = Human_Seurat_list[1:(length(Human_Seurat_list)-1)], project = GEO_ID)
        
        Human_Seurat_merged[["percent.mt"]] <- PercentageFeatureSet(Human_Seurat_merged, pattern = "^MT-")
        
        Human_Seurat_merged <- subset(Human_Seurat_merged, subset = percent.mt < Mitochondria)
        
        Human_Seurat_merged <- NormalizeData(Human_Seurat_merged, normalization.method = "LogNormalize", scale.factor = 10000)
        
        Human_Seurat_merged <- FindVariableFeatures(Human_Seurat_merged, selection.method = "vst", nfeatures = 2000)
        
        all.genes <- rownames(Human_Seurat_merged)
        
        Human_Seurat_merged <- ScaleData(Human_Seurat_merged, features = all.genes)
        
        Human_Seurat_merged <- RunPCA(Human_Seurat_merged, features = VariableFeatures(object = Human_Seurat_merged))
        
        if (Integrate){
          
          print("Starting integration")
          
          Mouse_Seurat_merged <- IntegrateLayers(object = Mouse_Seurat_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
          
          Mouse_Seurat_merged[["RNA"]] <- JoinLayers(Mouse_Seurat_merged[["RNA"]])
          
          Mouse_Seurat_merged <- FindNeighbors(Mouse_Seurat_merged, reduction = "integrated.cca", dims = Dims)
          
          Mouse_Seurat_merged <- FindClusters(Mouse_Seurat_merged, resolution = Resolution)
          
          Mouse_Seurat_merged <- RunUMAP(Mouse_Seurat_merged, dims = Dims, reduction = "integrated.cca")
          
          Mouse_Seurat_merged <- RunTSNE(Mouse_Seurat_merged, dims = Dims, reduction = "integrated.cca", check_duplicates = FALSE)
          
          dir.create("CloupeFiles", showWarnings = FALSE)
          
          print(paste0("Creating loupe file for ", GEO_ID))
          
          create_loupe_from_seurat(
            Mouse_Seurat_merged,
            output_dir = "CloupeFiles",
            output_name = paste0(GEO_ID,".Mouse.Integrated.",File_Format),
            dedup_clusters = FALSE,
            executable_path = NULL,
            force = TRUE)
          
          Human_Seurat_merged <- IntegrateLayers(object = Human_Seurat_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
          
          Human_Seurat_merged[["RNA"]] <- JoinLayers(Human_Seurat_merged[["RNA"]])
          
          Human_Seurat_merged <- FindNeighbors(Human_Seurat_merged, reduction = "integrated.cca", dims = Dims)
          
          Human_Seurat_merged <- FindClusters(Human_Seurat_merged, resolution = Resolution)
          
          Human_Seurat_merged <- RunUMAP(Human_Seurat_merged, dims = Dims, reduction = "integrated.cca")
          
          Human_Seurat_merged <- RunTSNE(Human_Seurat_merged, dims = Dims, reduction = "integrated.cca", check_duplicates = FALSE)
          
          dir.create("CloupeFiles", showWarnings = FALSE)
          
          print(paste0("Creating loupe file for ", GEO_ID))
          
          create_loupe_from_seurat(
            Human_Seurat_merged,
            output_dir = "CloupeFiles",
            output_name = paste0(GEO_ID,".Human.Integrated.",File_Format),
            dedup_clusters = FALSE,
            executable_path = NULL,
            force = TRUE)
          
        } else if (!Integrate){
          
          Mouse_Seurat_merged <- FindNeighbors(Mouse_Seurat_merged, dims = Dims)
          
          Mouse_Seurat_merged <- FindClusters(Mouse_Seurat_merged, resolution = Resolution)
          
          Mouse_Seurat_merged <- RunUMAP(Mouse_Seurat_merged, dims = Dims)
          
          Mouse_Seurat_merged <- RunTSNE(Mouse_Seurat_merged, dims = Dims, check_duplicates = FALSE)
          
          Mouse_Seurat_merged <- JoinLayers(Mouse_Seurat_merged)
          
          dir.create("CloupeFiles", showWarnings = FALSE)
          
          print(paste0("Creating loupe file for ", GEO_ID))
          
          create_loupe_from_seurat(
            Mouse_Seurat_merged,
            output_dir = "CloupeFiles",
            output_name = paste0(GEO_ID,".mouse.",File_Format),
            dedup_clusters = FALSE,
            executable_path = NULL,
            force = TRUE)
          
          Human_Seurat_merged <- FindNeighbors(Human_Seurat_merged, dims = Dims)
          
          Human_Seurat_merged <- FindClusters(Human_Seurat_merged, resolution = Resolution)
          
          Human_Seurat_merged <- RunUMAP(Human_Seurat_merged, dims = Dims)
          
          Human_Seurat_merged <- RunTSNE(Human_Seurat_merged, dims = Dims, check_duplicates = FALSE)
          
          Human_Seurat_merged <- JoinLayers(Human_Seurat_merged)
          
          dir.create("CloupeFiles", showWarnings = FALSE)
          
          print(paste0("Creating loupe file for ", GEO_ID))
          
          create_loupe_from_seurat(
            Human_Seurat_merged,
            output_dir = "CloupeFiles",
            output_name = paste0(GEO_ID,".Human.",File_Format),
            dedup_clusters = FALSE,
            executable_path = NULL,
            force = TRUE)
          
        }
        
        Seurat_merged_list <- append(Seurat_merged_list, Mouse_Seurat_merged)
        
        Seurat_merged_list <- append(Seurat_merged_list, Human_Seurat_merged)
        
        
        
      }
      
    }
    
  }
  
  if(Merge){
    
    return(Seurat_merged_list)
    
  } else {
    
    return(Seurat_separate_list)
    
  }
  
}  
