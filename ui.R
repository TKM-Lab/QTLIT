library(shiny)
options(mc.cores = detectCores() - 2)  # Use all but two core
options(shiny.maxRequestSize = Inf)

# UI 
custom_css <- "
  body {
    font-family: 'Tw Cen MT', sans-serif;
    background-color: #f4f4f4;
    color: #333;
  }
  .navbar {
    background-color:  #007958;
    height: 90px;
    border-bottom: 3px solid #FFBF00;
  }

  .navbar.qtl-active {
    background-color: #0047AB;  /* Change navbar color for QTL Analysis tab */
    border-bottom: 3px solid #6495ED;
  }

  .navbar-default .navbar-brand {
    color: #fff !important;
    font-weight: bold;
    font-size: 50px;
    text-align: center;
    position: absolute;
    left: 50%;
    transform:translateX(-50%);
  }
  
 .navbar .navbar-nav > li > a {
  color: #000 !important;
  font-weight: bold;
}
  .tab-title {
    font-size: 15px;
    font-weight: bold;
    color: #000;
  }

  .nav-tabs > li > a {
    color: #000;  /* Default: white for all tab titles */
    font-weight: bold;
  }

  .nav-tabs > li.active > a,
  .nav-tabs > li.active > a:focus,
  .nav-tabs > li.active > a:hover {
    color: #000;  /* Active tab: black */
    background-color: transparent;
    border: none;
  }

  .sidebar {
    background-color: #FFF5EE;
    border-right: 1px solid #ddd;
    padding: 15px;
  }

  .main-panel {
    background-color: #FFDEAD;
    padding: 20px;
    border-radius: 5px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  }

  .btn-primary {
    background-color: #FA5F55;
    border-color: #FFDEAD;
  }

  .btn-primary:hover {
    background-color: #0056b3;
    border-color: #004080;
  }

  .file-input {
    margin-bottom: 15px;
  }

  .plot-output {
    margin-top: 20px;
    border: 1px solid #ddd;
    border-radius: 5px;
    padding: 15px;
    background-color: #f9f9f9;
  }

  .log-output {
    background-color: #fff;
    border: 1px solid #ddd;
    border-radius: 5px;
    padding: 15px;
    margin-top: 20px;
    font-family: 'Courier New', monospace;
    font-size: 14px;
    color: #333;
  }

  /* Different background for second section */
  .second-section .main-panel {
    background-color: #badbe7;
  }
  
  .second-section .sidebar {
    background-color: #F0F8FF;
  }

  /* Different button color in second section */
  .second-section .btn-primary {
    background-color: #FFA07A;
    border-color: #FF7F50;
  }

  .second-section .btn-primary:hover {
    background-color: #0056b3;
    border-color: #004080;
  }

  /*head > link[rel='icon'] {  
    content: url('/home/tkm38/Downloads/passport size pic.jpg');
  }*/
"

ui <- navbarPage(
  title = "QTLIT", 
 position ="static-top", 
  header = tags$head(
    tags$style(HTML(custom_css)),
    # Optional external CSS
    # tags$link(rel = "stylesheet", type = "text/css", href = "file.css"),
    
    tags$script(HTML("
      function updateNavbarColor() {
        var activeTabText = $('.navbar-nav .active a').text().trim();
        if (activeTabText === 'QTL Analysis') {
          $('.navbar').addClass('qtl-active');
        } else {
          $('.navbar').removeClass('qtl-active');
        }
      }

      $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"]', function (e) {
        updateNavbarColor();
      });

      $(document).ready(function() {
        updateNavbarColor();
      });
    "))
  ),

  # Variant Calling Tab
  tabPanel(
    tags$div(class = "tab-title", "Variant Calling"),
    fluidPage(
      titlePanel("Variant Calling Analysis"),
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          tags$div(class = "file-input", 
                   fileInput("ref_genome_index", "Upload Reference Genome for alignment:", 
                             accept = c(".fasta", ".fa.gz"))),
          
          tags$div(class = "file-input", 
                   fileInput("fastqFiles", "Upload Paired-End FASTQ Files:", 
                             multiple = TRUE, accept = c(".fastq", ".fastq.gz"))),
          
          checkboxInput("skip_qc", "Skip Quality Control", FALSE),
          actionButton("run_qc", "Run Quality Control", class = "btn-primary"),
          actionButton("generateBAM", "Generate BAM Files", class = "btn-primary"),
          
          tags$div(class = "file-input", 
                   fileInput("bamFiles", "Upload BAM Files for Analysis:", 
                             multiple = TRUE, accept = c(".bam", ".bai"))),
          
          tags$div(class = "file-input", 
                   fileInput("ref_genome", "Upload Reference Genome for Variant Calling:", 
                             accept = ".fasta")),
          
          actionButton("runVariantCalling", "Run Variant Calling", class = "btn-primary")
        ),
        
        mainPanel(
          class = "main-panel",
          tabsetPanel(
            tabPanel("Logs", verbatimTextOutput("log_output")),
            tabPanel("QC Results", downloadButton("download_qc", "Download QC Results")),
            tabPanel("BAM Files", downloadButton("download_bam", "Download BAM Files")),
            tabPanel("Variant Calling", 
                     downloadButton("download_vcf", "Download VCF Results"),
                     downloadButton("downloadVCF", "Download Merged VCF"))
          )
        )
      )
    )
  ),

  # QTL Analysis Tab
  tabPanel(
    "QTL Analysis",
    tags$div(class = "second-section",
      fluidPage(
        titlePanel("Linkage Mapping and QTL Analysis"),
        sidebarLayout(
          sidebarPanel(
            class = "sidebar",
            tags$div(class = "file-input", 
                     fileInput("vcf_file", "Upload VCF File", accept = c(".vcf"))),
            
            tags$div(class = "file-input", 
                     fileInput("pheno_file", "Upload Pheno CSV", accept = c(".csv"))),
            
            textInput("parent_A", "Parent A for Alleles"),
            selectInput("phenotype_type", "Select Phenotype Type", 
                        choices = c("normal", "binary", "2part", "np")),
            
            uiOutput("pheno_column_selector"), 
            
            numericInput("siteMinCount", "Minimum Site Count", value = 0.0, min = 1, max = 10000),
            numericInput("siteMinAlleleFreq", "Minimum Allele Frequency", value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput("siteMaxAlleleFreq", "Maximum Allele Frequency", value = 1.0, min = 0, max = 1, step = 0.01),
            numericInput("minHeterozygous", "Minimum Heterozygous Proportion", value = 0.0, min = 0, max = 1, step = 0.01),
            numericInput("maxHeterozygous", "Maximum Heterozygous Proportion", value = 0.1, min = 0, max = 1, step = 0.01),
            
            checkboxInput("removeMinorSNPStates", "Remove Minor SNP States", value = FALSE),
            textInput("outfile", "Output File Name for Linkage Plot", "lmv_output.pdf"),
            
            radioButtons("map_function", "Map Function", choices = c("kosambi", "haldane")),
            selectInput("cross_type", "Cross Type", choices = c("riself", "dh", "bc", "bcsft")),
            numericInput("p_value", "Select P-value for Linkage Map", value = 1.0 , min= 0.000001, max = 1.0),
            
            radioButtons("qtl_model", "Select QTL Model", choices = c("scanone", "cim")),
            selectInput("qtl_method", "Select QTL Method", choices = c("em", "imp", "hk", "ehk")),
            numericInput("error_prob", "Error Probability for CIM", value = 0.001, min = 0.00001, max = 1),
            
            actionButton("run_analysis", "Run Analysis", class = "btn-primary")
          ),
          
          mainPanel(
            class = "main-panel",
            tabsetPanel(
              tabPanel("Data Summary", 
                       tags$div(class = "log-output", verbatimTextOutput("data_summary"))),
              
              tabPanel("Phenotype Summary", 
                       tags$div(class = "plot-output", plotOutput("pheno_summary_plot", height = "800px"))),
              
              tabPanel("Linkage Map", 
                       tags$div(class = "plot-output", plotOutput("linkage_plot")), 
                       downloadButton("download_linkage", "Download Linkage Map PDF")),
              
              tabPanel("Scanone/CIM Results", 
                       tags$div(class = "plot-output", plotOutput("scanone_cim_plot", height = "800px")),
                       downloadButton("download_qtl_zip", "Download All QTL Results"))
            )
          )
        )
      )
    )
  ),

  #Help and Download buttons
 tabPanel("Test Dataset",
    fluidPage(
      titlePanel("Test Dataset"),
      mainPanel(
        tags$h3("Sample Data"),
        tags$p("Click below to download the sample dataset for testing:"),
        tags$a(
          href = "sample_dataset.zip",
          "Download Sample Dataset",
          download = NA,
          target = "_blank"
        )
      )
    )
  ),
  
  tabPanel("Manual",
    fluidPage(
      titlePanel("User Manual"),
      mainPanel(
        tags$h3("User Guide (PDF)"),
        tags$p("Click below to download the user manual:"),
        tags$a(
          href = "user_manual.pdf",
          "Download User Manual (PDF)",
          download = NA,
          target = "_blank"
        )
      )
    )
  ),
  
  tabPanel("Contact",
    fluidPage(
      titlePanel("Contact Information"),
      mainPanel(
        tags$h3("Reach Out to Us"),
        tags$p("For support, questions, or feedback, contact:"),
        tags$ul(
          tags$li("Dr. Tapan Kumar Mondal"),
          tags$li("Email: ", tags$a(href="mailto:mondaltk@yahoo.com", "mondaltk@rediffmail.com")),
          tags$li("Institution: ICAR-NIPB, New Delhi"),
          tags$li("Phone: +91-99587-11064")
        ),
        tags$p("Website: ",
          tags$a(href = "https://www.nipb.res.in", target = "_blank", "https://www.nipb.res.in"))
      )
    )
  )
)
