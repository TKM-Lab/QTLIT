library(shiny)
qtl_result <- reactiveVal(NULL)
server <- function(input, output, session) {
  log <- reactiveVal("")
  fastqFiles <- reactiveVal(NULL)
  dir <- reactiveVal(NULL)
  vcf_file <- reactiveVal(NULL)
 dir_path <- "/srv/shiny-server/QTLIT"
#dir_path <- tempdir()  
#temp_dir <- file.path(dir_path,paste0("NIPB_", sample(100:999, 1))) 
  temp_dir <- file.path(dir_path, paste0("NIPB_", session$token))
  dir.create(temp_dir, recursive = TRUE)
  setwd(temp_dir)
  log_file <- file.path(getwd(), "progress.log")
  log_message <- function(msg) {
    write(msg, file = log_file, append = TRUE)
    log(paste0(log(), msg, "\n"))
  }
  output$log_output <- renderText({log() })
  
  # Quality Control Step
  observeEvent(input$run_qc, {
    log_message("Running QC...\n")
    
    if (!input$skip_qc && !is.null(input$fastqFiles)) {
      fastq_files <- input$fastqFiles$datapath
      qaSummary <- qa(fastq_files, type = "fastq", BPPARAM = SerialParam())
      
      QC_dir <- file.path(getwd(), "QC_results")
      if (!dir.exists(QC_dir)) dir.create(QC_dir, recursive = TRUE)
      
      report(qaSummary, type = "html", dest = QC_dir)
      log_message("QC report generated in QC_results folder.\n")
    } else if (input$skip_qc) {
      log_message("Skipping QC.\n")
    } else {
      log_message("Please upload FASTQ files for QC.\n")
    }
    log_message("QC Completed.\n")
  })
  #Alignment   
  observeEvent(input$generateBAM, {
    if (!is.null(input$fastqFiles) && !is.null(input$ref_genome_index)) {
      fastq_files <- input$fastqFiles$datapath
      original_names <- input$fastqFiles$name    
      print(original_names) 
      
      for (file in fastq_files) {
        if (!file.exists(file)) {
          log_message(paste("Error: File not found -", file, "\n"))
          return()
        }
      }
      
      # Rename the temporary files to original names
      renamed_files <- file.path(getwd(), original_names)
      file.copy(fastq_files, renamed_files, overwrite = TRUE)
      fastq_files <- renamed_files  
      fasta_path <- input$ref_genome_index$datapath
      genome_dir <- getwd()
      fastaFile <- FastaFile(fasta_path)
      log_message("Alignment started.\n")
      genome <- GmapGenome(fastaFile, genome_dir, create = TRUE)
      param <- GsnapParam(genome = genome, unique_only = TRUE, molecule = "DNA")
      session$userData$genome <- genome
      session$userData$param <- param
      
      log_message("Alignment Completed.\n")
      bam_dir <- file.path(getwd(), "BAM_files")
      if (!dir.exists(bam_dir)) dir.create(bam_dir, recursive = TRUE)
      
      sample_names <- gsub("_(R?[12]|1|2)\\.fastq(\\.gz)?$", "", original_names)
      fastq_df <- data.frame(sample = sample_names, file = fastq_files, stringsAsFactors = FALSE)
      paired_fastq_list <- split(fastq_files, sample_names)
      paired_fastq_list <- paired_fastq_list[sapply(paired_fastq_list, length) == 2]
      
      if (length(paired_fastq_list) == 0) {
        log_message("Error: No properly paired FASTQ files found.")
        return()
      }
      
      # Run GSNAP 
      for (sample in names(paired_fastq_list)) {
        fastq_pair <- paired_fastq_list[[sample]]
        bam_file <- file.path(bam_dir, paste0(sample))
        
        log_message(paste0("Running GSNAP for sample: ", sample, "\n"))
        print(paste("GSNAP input for", sample, ":", fastq_pair))
        
        alignments <- gsnap(fastq_pair[1], fastq_pair[2], param = session$userData$param, output = bam_file)
        
        for (file in fastq_pair) {
          if (file.exists(file)) {
            file.remove(file)
            print(paste("Deleted:", file))
          } else {
            print(paste("File not found:", file))
          }
        }
      }
      
      #BAM file renaming
      bam_files <- list.files(bam_dir, pattern = "\\.sam\\.bam$", full.names = TRUE)
      for (old_bam in bam_files) {
        new_bam <- gsub("\\.sam", "", old_bam)
        file.rename(old_bam, new_bam)
        
        old_bai <- paste0(old_bam, ".bai")
        new_bai <- paste0(new_bam, ".bai")
        
        if (!file.exists(old_bai)) {
          old_bai <- sub("\\.bam$", "", old_bam)
          new_bai <- sub("\\.bam$", "", new_bam)
        }
        
        if (file.exists(old_bai)) {
          file.rename(old_bai, new_bai)
        }
      }
      
      log_message("Completion of alignment. Download BAM files from BAM tab.\n")
    } else {
      log_message("Please upload FASTQ files and a reference genome.\n")
    }
  })
  
  # Variant Calling
  observeEvent(input$runVariantCalling, {
    log_message("Starting variant calling...\n")
    if (!is.null(input$bamFiles) && !is.null(input$ref_genome)) {
      bam_files <- input$bamFiles$datapath
      ref_fasta <- input$ref_genome$datapath
      
      # Rename temporary BAM files to original names
      original_bam_names <- input$bamFiles$name
      renamed_bam_files <- file.path(getwd(), original_bam_names)
      file.copy(bam_files, renamed_bam_files, overwrite = TRUE)
      bam_files <- renamed_bam_files  
      
      vcf_dir <- file.path(getwd(), "vcf_files")
      if (!dir.exists(vcf_dir)) dir.create(vcf_dir, recursive = TRUE)
      
      vcf_files <- c()
      sample_map_file <- file.path(vcf_dir, "sample_names.txt")
      write.table(data.frame(bam_files, original_bam_names), file = sample_map_file,
                  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
      
      # Run variant calling on each BAM file
      for (bam_file in bam_files) {
        sample_name <- gsub("\\.bam$", "", basename(bam_file))
        vcf_file <- file.path(vcf_dir, paste0(sample_name, ".vcf.gz"))
        
        command <- paste0("bcftools mpileup -Ou -f ", shQuote(ref_fasta), " ", shQuote(bam_file),
                          " | bcftools call -Ou -mv | bcftools filter -s LowQual -e 'QUAL<20 || DP>100' | bcftools view --types snps -m2 -M2 | bgzip > ", 
                          shQuote(vcf_file))
        system(command)
        system(paste("tabix -p vcf", shQuote(vcf_file)))
        
        reheader_vcf <- file.path(vcf_dir, paste0(sample_name, "_renamed.vcf.gz"))
        reheader_cmd <- paste("bcftools reheader -s", shQuote(sample_map_file), "-o", shQuote(reheader_vcf), shQuote(vcf_file))
        system(reheader_cmd)
        file.rename(reheader_vcf, vcf_file)
        
        vcf_files <- c(vcf_files, vcf_file)  # Store VCF file paths
      }
      
      # Index VCF files
      for (vcf_file in vcf_files) {
        system(paste("bcftools index", shQuote(vcf_file)))
      }
      
      # Merge VCF files if available
      if (length(vcf_files) > 0) {
        merged_vcf <- file.path(vcf_dir, "merged.vcf.gz")
        merge_command <- paste("bcftools merge -m id -O z -o", shQuote(merged_vcf), paste(vcf_files, collapse = " "))
        system(merge_command)
        system(paste("bcftools index", shQuote(merged_vcf)))
        
        log_message("Variant calling and merging completed.\n")
        
        vcf_path <- file.path(vcf_dir, "merged.vcf.gz")
        
        if (file.exists(vcf_path)) {
          log_message("Merged VCF file is ready for download and QTL analysis.")
          
          individual_vcfs <- list.files(vcf_dir, pattern = "\\.vcf(\\.gz)?$", full.names = TRUE)
          individual_vcfs <- individual_vcfs[individual_vcfs != vcf_path]  # Exclude merged file
          
          if (length(individual_vcfs) > 0) {
            file.remove(individual_vcfs)
          }
        }
        # Set merged VCF for further processing
        vcf_file(merged_vcf)
        
        # Proceed with downstream analysis
        vcf <- read.vcfR(merged_vcf)
        vcf_df <- as.data.frame(vcf@fix)
        
        output_dir <- "vcf_plots"
        if (!dir.exists(output_dir)) dir.create(output_dir)
        
        vcf_df$DP <- as.numeric(str_extract(vcf_df$INFO, "(?<=DP=)\\d+"))
        dp_plot <- ggplot(vcf_df, aes(x = DP)) + geom_histogram(binwidth = 5, fill = "blue", alpha = 0.6) +
          labs(title = "Depth of Coverage (DP) Distribution", x = "Read Depth", y = "Count") + theme_minimal()
        ggsave(file.path(output_dir, "DP_Distribution.png"), plot = dp_plot, width = 7, height = 5, dpi = 300)
        
        vcf_df$QUAL <- as.numeric(vcf_df$QUAL)
        qual_plot <- ggplot(vcf_df, aes(x = QUAL)) + geom_histogram(binwidth = 10, fill = "red", alpha = 0.6) +
          labs(title = "Variant Quality (QUAL) Distribution", x = "Quality Score", y = "Count") + theme_minimal()
        ggsave(file.path(output_dir, "Variant_Quality.png"), plot = qual_plot, width = 7, height = 5, dpi = 300)
        
        log_message("Variant calling and plotting completed.\n")
      } else {
        log_message("No VCF files generated. Cannot proceed with merging.\n")
      }
    } else {
      log_message("Please upload BAM and Reference Genome.\n")
    }
  })
  
  # Download Handlers
  output$download_qc <- downloadHandler(
    filename = "QC_results.zip",
    content = function(file) { zip::zipr(file, "QC_results") }
  )
  output$download_bam <- downloadHandler(
    filename = "BAM_files.zip",
    content = function(file) { zip::zipr(file, "BAM_files") }
  )
  output$download_vcf <- downloadHandler(
    filename = "VCF_plots.zip",
    content = function(file) { zip::zipr(file, "vcf_plots") }
  )
  output$downloadVCF <- downloadHandler(
    filename = function() {
      "merged.vcf.gz"
    },
    content = function(file) {
      vcf_path <- vcf_file()  
      
      if (!is.null(vcf_path) && file.exists(vcf_path)) {
        file.copy(vcf_path, file)
      } else {
        stop("Error: merged.vcf.gz file not found. Please run variant calling first.")
      }
    },
    contentType = "application/gzip"
  )
 
  # Linkage Mapping and QTL Analysis
  pheno_data <- reactive({
    req(input$pheno_file)
    read.csv(input$pheno_file$datapath, row.names = 1)
  })
  
  # Dynamic UI for phenotype column selection
  output$pheno_column_selector <- renderUI({
    req(pheno_data())
    selectInput("pheno_column", "Select Phenotype Column", choices = colnames(pheno_data()))
  })
  
  # Show phenotype data summary
  output$data_summary <- renderPrint({
    req(pheno_data())
    summary(pheno_data())
  })
  
  # QTL model radio button
  output$qtl_model_selector <- renderUI({
    radioButtons("qtl_model", "Select QTL Model", choices = c("scanone", "cim"))
  })
  
  # Observe the "Run Analysis" button
  observeEvent(input$run_analysis, {
    if (is.null(input$vcf_file)) {
      showNotification("Please upload a VCF file.", type = "error")
      return(NULL)
    }
    
    if (is.null(input$pheno_file)) {
      showNotification("Please upload a phenotype CSV file.", type = "error")
      return(NULL)
    }
    
    if (is.null(input$pheno_column)) {
      showNotification("Please select a phenotype column.", type = "error")
      return(NULL)
    }
    
    withProgress(message = 'Running QTL Analysis...', value = 0, {
      tryCatch({
        # Read and filter VCF
        GenoVCF <- readGenotypeTableFromPath(path = input$vcf_file$datapath)
        incProgress(0.1)
        
        filtered <- filterGenotypeTableSites(
          GenoVCF,
          siteMinCount = input$siteMinCount,
          siteMinAlleleFreq = input$siteMinAlleleFreq,
          siteMaxAlleleFreq = input$siteMaxAlleleFreq,
          minHeterozygous = input$minHeterozygous,
          maxHeterozygous = input$maxHeterozygous,
          removeMinorSNPStates = input$removeMinorSNPStates
        )
        
        exportGenotypeTable(filtered, file = "filtered_geno.vcf", format = "vcf")
        incProgress(0.2)
        
        tx <- readLines("filtered_geno.vcf")
        tx2 <- gsub("#CHROM", "CHROM", tx)
        writeLines(tx2, con = "cross.vcf")
        vcf <- read.table("cross.vcf", header = TRUE, stringsAsFactors = FALSE)
        vcf[, 10:ncol(vcf)] <- sapply(vcf[, 10:ncol(vcf)], function(x) substr(x, 1, 3))
        
        parent.A <- input$parent_A
        vcf[, 10:ncol(vcf)] <- t(apply(vcf[, 10:ncol(vcf)], 1, function(x) {
          ifelse(x == x[names(x) == parent.A], "A", 
                 ifelse(x == "0/1" | x == "1/0", "H", 
                        ifelse(x == "./.", "-", "B")))
        }))
        incProgress(0.3)
        
        marker_names <- vcf[, 3]
        marker_pos <- as.numeric(sub(".+_", "", marker_names)) / 1000000
        chromosomes <- sub("_.+", "", sub("S", "", marker_names))
        
        genotype_matrix <- vcf[, 10:ncol(vcf)]
        rownames(genotype_matrix) <- NULL
        genotype_matrix <- as.data.frame(t(genotype_matrix))
        
        genotype_matrix <- cbind(id = rownames(genotype_matrix), genotype_matrix)
        rownames(genotype_matrix) <- NULL
        
        marker_header <- c("id", as.character(marker_names))
        chr_header <- c("", chromosomes)
        pos_header <- c("", marker_pos)
        
        final_matrix <- rbind(marker_header, chr_header, pos_header, genotype_matrix)
        
        write.table(final_matrix, "cross.qtl.csv", quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE, na = "")
        incProgress(0.4)
        
        qtl.cross <- read.cross("csv", "", "cross.qtl.csv", crosstype = input$cross_type)
        
        pheno_df <- pheno_data()
        if (!("id" %in% colnames(pheno_df))) {
          pheno_df$id <- rownames(pheno_df)
        }
        pheno_df <- pheno_df[, c("id", setdiff(names(pheno_df), "id"))]
        numeric_phenos <- pheno_df[, sapply(pheno_df, is.numeric)]
        numeric_phenos <- numeric_phenos[, setdiff(names(numeric_phenos), "id")]
        qtl.cross <- add.phenos(qtl.cross, numeric_phenos)
        gt <- geno.table(qtl.cross)
        threshold <- 0.05
        distorted_markers <- rownames(gt[gt$P.value < threshold, ])
        length(distorted_markers)
        qtl.cross <- drop.markers(qtl.cross, distorted_markers)
        num_cols <- pheno_df[sapply(pheno_df, is.numeric)]
        
        cor_matrix <- cor(num_cols, use = "pairwise.complete.obs") #corr
        cor_data <- melt(cor_matrix)  # reshape to long format
        
        output$pheno_summary_plot <- renderPlot({
          pheno_col <- pheno_df[[input$pheno_column]]
          #p1 <- ggplot(pheno_df, aes(x = pheno_col)) + geom_histogram(binwidth = 1, fill = "skyblue") + labs(title = "Histogram")
          p1 <- ggplot(pheno_df, aes(y = pheno_col)) + geom_boxplot(fill = "orange") + labs(title = "Boxplot")
          missing_count <- sum(is.na(pheno_col))
          missing_percentage <- mean(is.na(pheno_col)) * 100
          p2 <- ggplot(pheno_df, aes(sample = pheno_col)) + stat_qq() + stat_qq_line(col = "red") +
            labs(title = paste("QQ Plot -", input$pheno_column),
                 subtitle = paste("Missing:", missing_count, "(", round(missing_percentage, 2), "%)")) +
            theme_minimal()
            
            # Violin plot
            p3 <- ggplot(pheno_df, aes(x = "", y = pheno_col)) +
              geom_violin(fill = "lightblue") +
              geom_boxplot(width = 0.1, fill = "orange", outlier.color = NA) +
              labs(title = paste("Violin Plot -", input$pheno_column), y = input$pheno_column) +
              theme_minimal()
            
            # Bar plot of mean ± SE
            mean_val <- mean(pheno_col, na.rm = TRUE)
            se_val <- sd(pheno_col, na.rm = TRUE) / sqrt(sum(!is.na(pheno_col)))
            
            summary_df <- data.frame(
              stat = "Mean",
              value = mean_val,
              se = se_val
            )
            
            p4 <- ggplot(summary_df, aes(x = stat, y = value)) +
              geom_bar(stat = "identity", fill = "steelblue") +
              geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.2) +
              labs(title = paste("Bar Plot of Mean ± SE -", input$pheno_column), y = input$pheno_column) +
              theme_minimal()
            
            # Correlation matrix heatmap
            p5 <- ggplot(cor_data, aes(Var1, Var2, fill = value)) +
              geom_tile(color = "white") +
              scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) +
              geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
              labs(title = "Correlation Heatmap", x = NULL, y = NULL) +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
          
            gridExtra::grid.arrange( p1, p2, p3, p4, nrow = 2, ncol = 3)
        })
        incProgress(0.5)
        
        suppressWarnings({
          qtl.cross <- mstmap.cross(qtl.cross, objective.fun = "ML", bychr = TRUE, id = "id",
                                    dist.fun = input$map_function, p.value = input$p_value, trace = FALSE,
                                    detectBadData = FALSE, return.imputed = FALSE)
        })
        
        # Linkage plot generation and save to PDF
        outfile <- input$outfile
        output$linkage_plot <- renderPlot({
          lmv.linkage.plot(qtl.cross, outfile)
        })
       
        output$download_linkage <- downloadHandler(
          filename = function() {
            input$outfile  
          },
          content = function(file) {
            pdf_path <- file.path(getwd(), input$outfile)
            if (file.exists(pdf_path)) {
              file.copy(pdf_path, file)
            } else {
              stop("The requested file does not exist.")
            }
          },
          contentType = "application/pdf"
        )
        
        selected_pheno_col <- which(colnames(pheno_df) == input$pheno_column)
        phenotype_type <- input$phenotype_type
        
        qtl.cross <- calc.genoprob(qtl.cross, map.function = input$map_function, step = 1)
        
        model_selected <- tolower(input$qtl_model)
        cat("QTL Model Selected:", model_selected, "\n")
        
        output_dir <- "qtl_analysis_output"
if (!dir.exists(output_dir)) dir.create(output_dir)

if (model_selected == "cim") {
  cat("Running CIM...\n")
  tryCatch({
    out <- cim(
      qtl.cross,
      pheno.col = selected_pheno_col,
      method = input$qtl_method,
      map.function = input$map_function,
      error.prob = input$error_prob,
      n.marcovar = 3,
      window = 10
    )
    qtl_result(out)

    # Save result as CSV
    write.csv(as.data.frame(out), file = file.path(output_dir, "qtl_result.csv"), row.names = FALSE)

  }, error = function(e) {
    showNotification("CIM failed.", type = "error")
    stop("CIM error: ", e$message)
  })
} else if (model_selected == "scanone") {
  cat("Running scanone...\n")
  tryCatch({
    out <- scanone(
      qtl.cross,
      pheno.col = selected_pheno_col,
      model = phenotype_type,
      method = input$qtl_method
    )
    qtl_result(out)

    # Save result as CSV
    write.csv(as.data.frame(out), file = file.path(output_dir, "qtl_result.csv"), row.names = FALSE)

  }, error = function(e) {
    showNotification("Scanone failed.", type = "error")
    stop("Scanone error: ", e$message)
  })

        } else {
          showNotification("Unknown QTL model selected. Please choose 'scanone' or 'cim'.", type = "error")
          stop("Invalid QTL model selected")
        }
        
        output$scanone_cim_plot <- renderPlot({
          p <- ntyped(qtl.cross, "mar")
          q <- ntyped(qtl.cross, what = c("ind", "mar"))
          r <- as.numeric(comparegeno(qtl.cross))
          
          plot1 <- ggplot(data.frame(x = 1:length(p), y = p),aes(x, y)) +
            geom_point(color = "darkgray") + labs(title = "Markers vs Individuals")
          
          plot2 <- ggplot(data.frame(x = 1:length(q), y = q),aes(x, y)) +
            geom_point(size= 1, color = "maroon") + labs(title = "Individuals vs Markers")
          
          plot3 <- ggplot(data.frame(x = 1:length(r), y = r),aes(x, y)) +
            geom_point(color = "blue") + labs(title = "Comparison of Genotypes")
          
         # plot4 <- ggplot(data.frame(x = 1:length(out$lod), y = out$lod),aes(x, y)) +
         #   geom_line(color = "black") + labs(title = "QTL Scanone CIM Plot")
          
          out$chr <- as.factor(out$chr) 
          out$cum_pos <- ave(out$pos, out$chr, FUN = function(x) seq_along(x)) 
          chr_lengths <- tapply(out$pos, out$chr, max)
          chr_offset <- c(0, cumsum(as.numeric(chr_lengths))[-length(chr_lengths)])
          names(chr_offset) <- levels(out$chr)
          out$cum_pos <- out$pos + chr_offset[out$chr]
          chr_ticks <- tapply(out$cum_pos, out$chr, mean)
          
          # Add labels for LOD > 3
          if (is.null(out$marker) || !("marker" %in% colnames(out))) {
            out$marker <- rownames(out)
          }
          out$marker_label <- ifelse(out$lod > 3, out$marker, NA)
          
          plot4<- ggplot(out, aes(x = cum_pos, y = out$lod, color = chr)) +
            geom_line(size = 1.5) +
            scale_x_continuous(breaks = chr_ticks, labels = names(chr_ticks)) +
            labs(x = "Chromosome", y = "LOD Score", title = "LOD Score Plot") +
            geom_hline(yintercept = 3, linetype = "dotted", color = "black", size = 1) +
            ggrepel::geom_text_repel(
              aes(label = marker_label),
              na.rm = TRUE,
              size = 3,
              max.overlaps = 20,
              show.legend = FALSE
            ) +
            theme_minimal() +
            theme(legend.position = "none")
          
          #plot5 <- ggplot(data.frame(
          #  x = 1:length(out$lod), 
          #  y = out$lod, 
          #  chr = factor(out$chr)
         # )) +
          #  geom_line(aes(x, y, color = chr), size = 1.2) +
          #  scale_color_viridis_d(option = "turbo")+
          #  geom_hline(yintercept = 3, linetype = "dotted", color = "black", size = 1) +
          #  labs(x = "Chromosome Position", y = "LOD Score", title = " LOD Score Plot") +
          #  theme_minimal() #+
            #guides(color = "none") 
          gridExtra::grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
          ggsave(file.path(output_dir, "plot1.pdf"), plot1)
          ggsave(file.path(output_dir, "plot2.pdf"), plot2)
          ggsave(file.path(output_dir, "plot3.pdf"), plot3)
          ggsave(file.path(output_dir, "plot4.pdf"), plot4)
          })
        
        output$download_qtl_zip <- downloadHandler(
          filename = function() {
            "qtl_analysis.zip"
          },
          content = function(file) {
            zip::zipr(
              zipfile = file,
              files = list.files("qtl_analysis_output", full.names = TRUE)
            )
          },
          contentType = "application/zip"
        )
        
        incProgress(0.8)
        
      }, error = function(e) {
        showNotification(paste("Error during analysis:", e$message), type = "error")
      })
    })
  })
}
