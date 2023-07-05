#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("timecourse")


options(repos = BiocManager::repositories())

list.of.packages <- c("shiny", "shinythemes","shinyWidgets","reactable","tidyverse","plotly","cowplot","GGally","ggrepel","stringi","FactoMineR","factoextra","igraph","ggraph","tidygraph","pheatmap","readr", "ComplexUpset","eulerr","RColorBrewer","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#"BiocManager","timecourse",

library(shiny)
# library(reticulate)
library(shinythemes)
library(shinyWidgets)
library(reactable)
library(tidyverse)
library(plotly)
library(cowplot)
library(GGally)
library(ggrepel)
library(stringi)
library(FactoMineR)
library(factoextra)
library(igraph)
library(ggraph)
library(tidygraph)
library(pheatmap)
library(timecourse)
library(readr)
library(ComplexUpset)
library(eulerr)
library(RColorBrewer)


# use_condaenv("r-reticulate")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("timecourse")
#library(timecourse)


locale_input <- read.csv("org_shortcut_prot_ribo_withMouse2.csv")
# go.df<-read.csv("human_go_referencefile.csv")%>%
#   dplyr::mutate(GO=if_else(stringr::str_sub(GO,1,1)==" ",stringr::str_sub(GO,1,-1),GO))


goinput_start<-read.csv("human_mouse_go.csv")%>%
  dplyr::mutate(GO=if_else(stringr::str_sub(GO,1,1)==" ",stringr::str_sub(GO,2,-1),GO))%>%
  dplyr::filter(GO != "")

#go.df<-read.csv("C:/Harper/Side_analysis/go_parser/human_go_referencefile.csv")
int_df<-read.delim2("BioPlex_293T_Network_10K_Dec_2019.tsv")
bioplex_binary<-read.csv("bioplex_binary_interactions_final.csv")
#C:/Harper/Experiments/Exp0001_CellLines_Chaperone_Choice/BioPlex_293T_Network_10K_Dec_2019.tsv
corum_initial<-read.csv("corum_ref3.csv")
#C:/Harper/Experiments/Exp0003_MLN_Ctrl_Proteome/tmt_proteome_2way/
locale2<- read.csv("MH_organelle list_fromAlban_editedERproteins.csv")

# locale<-locale%>%
#   dplyr::select(Gene.Symbol,Compartment)%>%
#   dplyr::distinct()

locale2<-locale2%>%
  dplyr::select(Gene.Symbol,Compartment)%>%
  dplyr::distinct()


list()

# Define UI for application that draws a histogram
ui <- fluidPage(theme=shinytheme("united"),
  
  # Application title
  titlePanel("TMT Protein Abundance Analysis"),
  navbarPage(
    "ShinyApp Multiple Condition Comparisons",
    tabPanel("Data Import",
             sidebarPanel(
               tags$h3("Protein TMT Data Import:"),
               selectInput("txthumanmouse", "Is this human or mouse data?",
                           choices=list("human"="human","mouse"="mouse"),selected = "human"),
               fileInput('file1', 'Select the peptide-level TMT quant file',
                         accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
               fileInput('file2', 'Select the experimental design .csv file',
                         accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
               textInput("text1004","contaminants designation (case sensitive)","contaminant"),
               numericInput("txtNumPep","Number of peptides per protein required=",1,min=1),
               numericInput("num101","Sum signal:noise filter=",200,min=0),
               numericInput("num102","Isolation Specificity (0.5 is default)=",0.5,min=0),
               tags$h3("Normalize to (All default):"),
               textInput("txt402","Normalize to All (default) or One Protein", "All"),
               selectInput("txtMedianSumOptions", "Channel Sum or Median Normalization (Note: Protein Aggregation always summed):",
                           choices=list("Sum"="Sum","Median"="Median"),selected = "Sum"),
               plotOutput("plotCV"),
               tags$h3("Median CV Condition 1="),
               textOutput("medianC1"),
               tags$h3("Median CV Condition 2="),
               textOutput("medianC2"),
               tags$h3("Median CV Condition 3="),
               textOutput("medianC3"),
               tags$h3("Median CV Condition 4="),
               textOutput("medianC4"),
               downloadButton('labelCheckDownloadPlot', 'Download Label Check Plot')
             ),
             mainPanel(
               h1("TMT Data Overview"),
               
               h4("Pre-normalization"),
               plotOutput("plot1"),
               h4("Post-normalization"),
               plotOutput("plot2"),
               h4("Label Check"),
               plotOutput("plot3")               
             )
             
    ),
    tabPanel("Correlations",
             sidebarPanel(
               tags$h3("Number of proteins post-filtering:"),
               textOutput("totalproteins"),
               downloadButton('CorrDownloadPlot', 'Download Replicate Correlation Plot')
             ),
             mainPanel(
               h1("Replicate Correlation Plot"),
               h4("Plot"),
               plotOutput("plot4")
             )),
    tabPanel("PCA",
             sidebarPanel(
               downloadButton('PCADownloadPlot', 'Download PCA Plot'),
               selectInput("txtXpca", "PCA x-axis component:",
                           choices=list("PC1"="PC1","PC2"="PC2","PC3"="PC3","PC4"="PC4","PC5"="PC5",
                                        "PC6"="PC6","PC7"="PC7","PC8"="PC8","PC9"="PC9","PC10"="PC10"),selected = "PC1"
               ),
               selectInput("txtYpca", "PCA y-axis component:",
                           choices=list("PC1"="PC1","PC2"="PC2","PC3"="PC3","PC4"="PC4","PC5"="PC5",
                                        "PC6"="PC6","PC7"="PC7","PC8"="PC8","PC9"="PC9","PC10"="PC10"),selected = "PC2"
               ),
               selectInput("txtZpca", "Proteins that Drive this PC:",
                                  choices=list("PC1"="PC1","PC2"="PC2","PC3"="PC3","PC4"="PC4","PC5"="PC5",
                                               "PC6"="PC6","PC7"="PC7","PC8"="PC8","PC9"="PC9","PC10"="PC10"),selected = "PC1"
               ),
               plotOutput("plotPCA"),
               plotOutput("plotPCAadd")
             ),
             mainPanel(
               h1("PCA PC Drivers"),
               downloadButton('driversDownload', 'Download All PC Weights'),
               reactableOutput("PCAdrivers")
             )),
    tabPanel("Statistics",
             sidebarPanel(
               tags$h3("Individual Protein Abundance Plots"),
               selectInput("txtScale", "Row Scaling options",
                           choices=list("RowMean"="RowMean","Zscore"="Zscore"),selected = "Zscore"),
               selectInput("anovaorall", "ANOVA Significant or All Proteins for Heatmap",
                           choices=list("ANOVA"="ANOVA","All"="All"),selected = "ANOVA"),
               selectInput("txtcomparettest", "Comparison to One condition (Priority  #1) or All Possible",
                           choices=list("One"="One","All"="All"),selected = "All"),
               numericInput("num77","q-value cut-off",0.05,min=0,max=1),
               numericInput("num78","fold change cut-off (|abs|)",1.5,min=0),
               selectInput("txtclustercols", "Cluster columns?",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               numericInput("numclusters","Number of Clusters=",2,min=2),
               reactableOutput("numbercluster"),
               plotOutput("numclusterplot"),
               width=3
             ),
             mainPanel(
               tags$h3("If heatmap does not render, resize the window (will take a few minutes)"),
               downloadButton('plotheatmapprint', 'Download Heatmap'),
               downloadButton('printclustersplot', 'Download Clusters Plot'),
               plotOutput("heatmapplot"),
               plotOutput("plotclustertrace"),
               width=9
             )

    ),
    tabPanel("Results",
             mainPanel(
               h1("Volcano Plot - Student's T-test"),
               h3("Box selection of points will produce table of selected protein results"),
               downloadButton("alldownload", label="Download All Results"),
               reactableOutput("table16")
               #plotlyOutput("plot6"),
               #reactableOutput("brush")
               
             )
             
    ),
    tabPanel("Timecourse",
             sidebarPanel(
               h1("Timecourse analysis: Hotelling T2"),
               h3("Requires temporal data ordered by priority (earliest time point priority 1)"),
               selectInput("txttimecourse", "Is this timecourse data?",
                           choices=list("Yes"="Yes","No"="No"),selected = "Yes"),
               numericInput("numrepstimecourse","Number of replicates for each time point:",3,min=2,max=4),
               numericInput("numtimecourse1","Time point 0 for Hotelling:",1,min=1,max=8),
               numericInput("numtimecourse2","Comparison time point for Hotelling (based on priority):",2,min=1,max=8),
               textInput("text113","Gene Name (ALL CAPS)","GAPDH"),
               plotOutput("plottcprot")
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(
               downloadButton("printhotellingplot", label="Download Hotelling Plot"),
               downloadButton("plottcprotprint2", label="Download Protein Timecourse"),
               
               #numericInput("alphanum1","Alpha for volcano and correlation plot points",0.2,min=0, max=1),
               plotlyOutput("plothotelling"),
               reactableOutput("hotellingclick")
               
             )
             
    ),
    tabPanel("Timecourse Export",
             sidebarPanel(
               h1("Timecourse analysis: Hotelling T2"),
               h3("Requires temporal data ordered by priority (earliest time point priority 1)")
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(
               downloadButton("hotellingdownload", label="Download Hotelling Results"),
               #numericInput("alphanum1","Alpha for volcano and correlation plot points",0.2,min=0, max=1),
               reactableOutput("hotellingdf2")
               #plotlyOutput("plot6"),
               #reactableOutput("brush")
               
             )
             
    ),
    tabPanel("Volcano Plot",
             sidebarPanel(
               h1("Volcano Comparisons across different conditions: pairwise"),
               h3("Compares two conditions at a time"),
               numericInput("numtimecourse1vol","Baseline condition for Volcano:",1,min=1,max=8),
               numericInput("numtimecourse2vol","Comparison condtion (based on priority):",2,min=1,max=8),
               selectInput("txtinversion", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               numericInput("num55","Bottom and Top Ranked Proteins by Fold Change",5,min=1,max=40),
               plotOutput("rankfoldchange"),
               
               # numericInput("num3","q-value cut-off",0.05,min=0,max=1),
               # numericInput("num4","fold change cut-off (|abs|)",1.5,min=0),
               textInput("text11","Gene Name (ALL CAPS)","GAPDH"),
               plotOutput("plot5"),
               numericInput("numtimecourse3vol","x-axis baseline condition:",1,min=1,max=8),
               numericInput("numtimecourse4vol","x-axis experimental condition:",2,min=1,max=8),
               numericInput("numtimecourse5vol","y-axis baseline condition:",1,min=1,max=8),
               numericInput("numtimecourse6vol","y-axis experimental condition:",3,min=1,max=8),
               plotlyOutput("volcanoplotnorm"),
               width=5
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(
               downloadButton('printvolcano', 'Download Volcano Plot'),
               downloadButton('rankfoldchangeprint2', 'Download Top Up/Down FC Plot'),
               downloadButton('proteinBarplotDownloadPlot','Download Single Protein Barplot'),
               downloadButton('printcorrtwocond', 'Download Ratio Correlation Plot'),
               
               #numericInput("alphanum1","Alpha for volcano and correlation plot points",0.2,min=0, max=1),
               plotlyOutput("volcanoplot"),
               reactableOutput("volcanoclick"),
               
               width=7
               #plotlyOutput("plot6"),
               #reactableOutput("brush")
               
             )
             
    ),

    tabPanel("Volcano (labels)",
             sidebarPanel(
               h1("Labeled Volcano Plot"),
               numericInput("numtimecourse7vol","Time point 0 for Volcano:",1,min=1,max=8),
               numericInput("numtimecourse8vol","Comparison time point (based on priority):",2,min=1,max=8),
               numericInput("alphanum2","Alpha for volcano and correlation plot points",0.5,min=0, max=1)
               # ,
               # plotOutput("plotclustertrace")
             ),
             mainPanel(

               
               downloadButton('volcanoDownloadPlot', 'Download Volcano Plot (labels)'),
               plotOutput("plotVL")
             )

    ),
    
    tabPanel("Venn/Upset",
             sidebarPanel(
               h1("Venn Diagram"),
               selectInput("txtupdownoverlap", "Which regulated proteins for Venn:",
                           choices=list("up"="up","down"="down","both"="both"),selected = "both"),
               downloadButton('printploteuler', 'Download Venn Diagram'),
               plotOutput("ploteuler")
             ),
             mainPanel(
               h1("Upset Plot"),
               downloadButton('printplotupset', 'Download Upset Plot'),
               plotOutput("plotupset")
             )
             
    ),

    tabPanel("BioPlex Top",
             sidebarPanel(
               # h1("BioPlex Interactome"),
               # textInput("txtbioplex","Gene Name for Bioplex","GAPDH"),
               # plotOutput("bioplexExp"),
               # reactableOutput("tablebioplex"),
               h1("BioPlex Interactome Top Hits"),
               numericInput("numconbioplex3","Baseline condition:",1,min=1,max=8),
               numericInput("numconbioplex4","Comparison condtion (based on priority):",2,min=1,max=8),
               numericInput("num10001","minimium number of interactors",4,min=2),
               downloadButton("alldownloadbioplex", label="Download Bioplex Results"),
               reactableOutput("tableBinaryBioplex"),
               width = 5
             ),
             mainPanel(
               
               numericInput("numbioplex1","TopN most up regulated interactomes:",10,min=1),
               downloadButton("printbiolplextop", label="Download Top Up Bioplex Plot"),
               plotOutput("bioplextopN"),
               numericInput("numbioplex2","TopN most down regulated interactomes:",10,min=1),
               downloadButton("printbiolplexbottom", label="Download Top Down Bioplex Plot"),
               
               
               plotOutput("bioplexbottomN"),
               width = 7

             )

    ),
    tabPanel("BioPlex User",
             sidebarPanel(
               h1("BioPlex Interactome User Defined"),
               textInput("txtbioplex","Gene Name for Bioplex","GAPDH"),
               numericInput("numconbioplex1","Baseline condition:",1,min=1,max=8),
               numericInput("numconbioplex2","Comparison condtion (based on priority):",2,min=1,max=8),
               downloadButton("printbioplexExp", label="Download Single Condition Bioplex Plot"),
               plotOutput("bioplexExp"),
               downloadButton("printbioplexExpAll", label="Download Multiple Condition Bioplex Plot"),
               plotOutput("bioplexExpAll"),
               reactableOutput("tablebioplex"),
               h3("Protein Interactome Level: (baseline:condition above are different than last tab)"),
               numericInput("num10002","minimium number of interactors",4,min=2),
               reactableOutput("tableBinaryBioplex2"),
               width = 6
             ),
             mainPanel(
               # downloadButton("printbioplexNetwork", label="Download Network Plot"),
               plotOutput("bioplexNetwork"),
               width = 6

             )

    ),
    tabPanel("Corum",
             sidebarPanel(
               h1("Corum Complex Analysis"),
               downloadButton("alldownloadCorum",label = "Download All Corum Results"),
               h3("Complex Selection Parameter"),
               numericInput("numcon1","Baseline condition:",1,min=1,max=8),
               numericInput("numcon2","Comparison condtion (based on priority):",2,min=1,max=8),
               numericInput("numcomplexmin","Minimum number of identified complex members",4,min=2),
               numericInput("complexranknum","Complex Rank Number for Selection",1,min=1),
               h3("Complex-level Table"),
               reactableOutput("tablecorum"),
               h3("Top Regulated Complexes:"),
               downloadButton("printcorumExp", label="Download Top Up Corum Plot"),
               numericInput("numcorum","TopN most up regulated complexes:",10,min=1),
               plotOutput("corumExp"),
               downloadButton("printcorumExp2", label="Download Top Down Corum Plot"),
               numericInput("numcorum2","TopN most down regulated complexes:",10,min=1),
               plotOutput("corumExp2"),
               width = 6
             ),
             mainPanel(
               h3("Individual Complex Selection"),
               plotOutput("corumNetwork"),
               downloadButton("printcorumExp3", label="Download Protein (Single Comparison) Corum Plot"),
               plotOutput("corumExp3"),
               reactableOutput("tablecorumOneComplex"),
               downloadButton("printcorumExpallcond", label="Download Protein (Mulit-comparison) Corum Plot"),
               plotOutput("corumExpallcond"),
               width = 6

             )


    ),
    tabPanel("GO Clusters",
             sidebarPanel(
               
               
               numericInput("numcluster","cluster for GO enrichment",1,min=1),
               numericInput("numgo","q-value cut-off",0.05,min=0,max=1),
               numericInput("numtopn","topN number of GO terms (q-value ranked):",10,min=1,max=100),
               #downloadButton('plotGOall', 'Download GO scatter'),
               #downloadButton('TopGODownloadPlot', 'Download Top GO'),
               downloadButton("GOclusterdownload", label="Download Cluster GO Results"),
               downloadButton("printplotGO", label="Download All GO Plot (Cluster)"),
               plotlyOutput("plotGO")
             ),
             mainPanel(
               h3("Click on point to explore Proteins in GO term:"),
               downloadButton("printplotGOtopn2", label="Download TopN GO Plot (Cluster)"),
               plotlyOutput("plotGOtopn2"),
               reactableOutput("goclick")

             )

    ),
    tabPanel("GO Clusters 2",
             sidebarPanel(
               h3("Significance designations derived from user-defined cut-offs from the GO Clusters enrichment tab"),
               h3("GO terms with associated proteins:"),
               # downloadButton("alldownloadGO", label="Download All Results"),
               reactableOutput("tableGO"),
               width=5
             ),
             mainPanel(
               textInput("txtlistnum","List of GO ranks separated by commas ex: 1,14,24,16"),
               numericInput("numgo2","q-value cut-off",0.05,min=0,max=1),
               downloadButton("printplotGOuser", label="Download User-defined GO Plot (Cluster)"),
               plotOutput("plotGOuser"),
               h3("GO terms Ranked:"),
               reactableOutput("tableGO2"),
               width=7

             )


    ),
    tabPanel("GO Regulated",
             sidebarPanel(
               h3("Signficance designations from statistics section"),
               selectInput("txtupdown", "Which regulated proteins for GO:",
                           choices=list("up"="up","down"="down","both"="both"),selected = "up"),
               numericInput("numpriority1","Baseline condition for Volcano:",1,min=1,max=8),
               numericInput("numpriority2","Comparison condtion (based on priority):",2,min=1,max=8),
               numericInput("numgo1","q-value cut-off",0.05,min=0,max=1),
               numericInput("numtopn1","topN number of GO terms (q-value ranked):",10,min=1,max=100),
               downloadButton("GOregulateddownload", label="Download Regulated GO Results"),
               downloadButton("printplotGO2", label="Download All GO Plot (Regulated)"),
               plotlyOutput("plotGO2")
             ),
             mainPanel(
               h3("Click on point to explore Proteins in GO term:"),
               downloadButton("printplotGOtopn3", label="Download TopN GO Plot (Regulated)"),
               plotlyOutput("plotGOtopn3"),
               reactableOutput("goclick3")
               
             )
             
    ),
    tabPanel("GO Regulated 2",
             sidebarPanel(
               h3("Significance designations derived from user-defined cut-offs from the GO Regulated enrichment tab"),
               h3("GO terms with associated proteins:"),
               # downloadButton("alldownloadGO", label="Download All Results"),
               reactableOutput("tableGOreg"),
               width=5
             ),
             mainPanel(
               textInput("txtlistnum3","List of GO ranks separated by commas ex: 1,14,24,16"),
               numericInput("numgo3","q-value cut-off",0.05,min=0,max=1),
               downloadButton("printplotGOuserReg", label="Download User-defined GO Plot (Regulated)"),
               plotOutput("plotGOuserReg"),
               h3("GO terms Ranked:"),
               reactableOutput("tableGO2reg"),
               width=7
               
             )
             
             
    ),

    tabPanel("Localization v1",
             sidebarPanel(
               tags$h3("Localization Volcano"),
               selectInput("localization", "Localization for Volcano Plot:",
                           choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                        "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                        "Ribosome"="Ribosome",
                                        "Golgi.apparatus..160."="Golgi.apparatus..160." ,"Golgi.membrane..87."="Golgi.membrane..87.",
                                        "Golgi.membrane.associated..73."="Golgi.membrane.associated..73.",
                                        "Cis.Golgi..29."="Cis.Golgi..29.","Trans.Golgi..36."="Trans.Golgi..36."),selected = "Lysosome"),
               numericInput("numtimecourse1loc","Baseline condition (based on priority):",1,min=1,max=8),
               numericInput("numtimecourse2loc","Comparison condition (based on priority):",2,min=1,max=8),
               multiInput("locales", "Fold change across localizations (Select all that apply):",
                           choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                        "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                        "Ribosome"="Ribosome",
                                        "Golgi.apparatus..160."="Golgi.apparatus..160." ,"Golgi.membrane..87."="Golgi.membrane..87.",
                                        "Golgi.membrane.associated..73."="Golgi.membrane.associated..73.",
                                        "Cis.Golgi..29."="Cis.Golgi..29.","Trans.Golgi..36."="Trans.Golgi..36."),selected = "Lysosome"),
               selectInput("txtinversionloc", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               # ,
               downloadButton('downloadlocalizationdata', 'Download Localization Data'),
               downloadButton('printcorrcolorcompartment', 'Download Correlation Localization'),
               plotOutput("corrcolorcompartment"),
               numericInput("numtimecourse1loc2","x-axis baseline condition:",1,min=1,max=8),
               numericInput("numtimecourse2loc2","x-axis experimental condition:",2,min=1,max=8),
               numericInput("numtimecourse3loc2","y-axis baseline condition:",1,min=1,max=8),
               numericInput("numtimecourse4loc2","y-axis experimental condition:",3,min=1,max=8),
               downloadButton('printcorrcolorcompartmentmulti', 'Download Correlation Localization Ratio'),
               plotOutput("corrcolorcompartmentmulti")
             ),
             mainPanel(
               downloadButton('printplot7', 'Download Volcano Localization'),
               plotOutput("plot7"),
               numericInput("numlocale1","x-axis minimium",-2),
               numericInput("numlocale2","x-axis maximium",2),
               downloadButton('printplot8', 'Download Single Condition Foldchange for Localization'),
               plotOutput("plot8"),
               downloadButton('printplotmultilocale', 'Download Multi-condition Foldchange for Localization'),
               plotOutput("plotmultilocale")
             )

    ),
    tabPanel("Localization (labels) v1",
             sidebarPanel(
               tags$h3("Localization Volcano with Labels"),
               numericInput("numtimecourse1loc3","Baseline condition (based on priority):",1,min=1,max=8),
               numericInput("numtimecourse2loc3","Comparison condition (based on priority):",2,min=1,max=8),
               selectInput("txtinversionloc44", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               selectInput("localization2", "Localization for Volcano Plot:",
                           choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                        "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                        "Ribosome"="Ribosome",
                                        "Golgi.apparatus..160."="Golgi.apparatus..160." ,"Golgi.membrane..87."="Golgi.membrane..87.",
                                        "Golgi.membrane.associated..73."="Golgi.membrane.associated..73.",
                                        "Cis.Golgi..29."="Cis.Golgi..29.","Trans.Golgi..36."="Trans.Golgi..36."),selected = "Lysosome")
               # ,
               # downloadButton('volcanoLocaleDownloadPlot2', 'Download Volcano Locale Plot (Labels)')
             ),
             mainPanel(
               downloadButton('printplot44', 'Download Volcano Localization Labels'),
               plotOutput("plot44")
             )

    ),
    
    tabPanel("Localization v2",
             sidebarPanel(
               tags$h3("Localization Volcano"),
               selectInput("localization5", "Localization for Volcano Plot:",
                           choices=list("Actin binding proteins"="Actin binding proteins", 
                                        "Endosome"="Endosome",    "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "ER_sheet"="ER_Sheet","Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" , "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear (HA Frac)"="Nuclear (HA Frac)","Nuclear pore complex"="Nuclear pore complex", 
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome","Receptor"=="Receptor",
                                        "Ribosome"="Ribosome"),selected = "Lysosome"),
               numericInput("numtimecourse1loc5","Baseline condition (based on priority):",1,min=1,max=8),
               numericInput("numtimecourse2loc5","Comparison condition (based on priority):",2,min=1,max=8),
               multiInput("locales5", "Fold change across localizations (Select all that apply):",
                          choices=list("Actin binding proteins"="Actin binding proteins", 
                                       "Endosome"="Endosome",    "ER_high_curvature"= "ER_high_curvature",
                                       "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "ER_sheet"="ER_Sheet","Ergic/cisGolgi"="Ergic/cisGolgi",
                                       "Golgi"="Golgi" , "Lysosome"="Lysosome",
                                       "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                       "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                       "Nuclear (HA Frac)"="Nuclear (HA Frac)","Nuclear pore complex"="Nuclear pore complex", 
                                       "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome","Receptor"=="Receptor",
                                       "Ribosome"="Ribosome"),selected = "Lysosome"),
               selectInput("txtinversionloc5", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               downloadButton('printcorrcolorcompartmentloc5', 'Download Correlation Localization'),
               plotOutput("corrcolorcompartmentloc5"),
               numericInput("numtimecourse1loc6","x-axis baseline condition:",1,min=1,max=8),
               numericInput("numtimecourse2loc6","x-axis experimental condition:",2,min=1,max=8),
               numericInput("numtimecourse3loc6","y-axis baseline condition:",1,min=1,max=8),
               numericInput("numtimecourse4loc6","y-axis experimental condition:",3,min=1,max=8),
               downloadButton('printcorrcolorcompartmentmultiloc5', 'Download Correlation Localization Ratio'),
               plotOutput("corrcolorcompartmentmultiloc5")
             ),
             mainPanel(
               downloadButton('printplot7loc2', 'Download Volcano Localization'),
               plotOutput("plot7loc2"),
               numericInput("numlocale1loc2","x-axis minimium",-2),
               numericInput("numlocale2loc2","x-axis maximium",2),
               downloadButton('printplot8loc2', 'Download Single Condition Foldchange for Localization'),
               plotOutput("plot8loc2"),
               downloadButton('printplotmultilocaleloc2', 'Download Multi-condition Foldchange for Localization'),
               plotOutput("plotmultilocaleloc2")
             )
             
    ),
    tabPanel("Localization (labels) v2",
             sidebarPanel(
               tags$h3("Localization v2 Volcano with Labels"),
               numericInput("numtimecourse1loc8","Baseline condition (based on priority):",1,min=1,max=8),
               numericInput("numtimecourse2loc8","Comparison condition (based on priority):",2,min=1,max=8),
               selectInput("txtinversionloc88", "Reverese Condition Order (main volcano only):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               selectInput("localization8", "Localization for Volcano Plot:",
                           choices=list("Actin binding proteins"="Actin binding proteins", 
                                        "Endosome"="Endosome",    "ER_high_curvature"= "ER_high_curvature",
                                        "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "ER_sheet"="ER_Sheet","Ergic/cisGolgi"="Ergic/cisGolgi",
                                        "Golgi"="Golgi" , "Lysosome"="Lysosome",
                                        "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                        "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                        "Nuclear (HA Frac)"="Nuclear (HA Frac)","Nuclear pore complex"="Nuclear pore complex", 
                                        "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome","Receptor"=="Receptor",
                                        "Ribosome"="Ribosome"),selected = "Lysosome")
               # ,
               # downloadButton('volcanoLocaleDownloadPlot2', 'Download Volcano Locale Plot (Labels)')
             ),
             mainPanel(
               downloadButton('printplot88', 'Download Volcano Localization Labels'),
               plotOutput("plot88")
             )
             
    ),  
  
  tabPanel("User Select",
             sidebarPanel(
               tags$h3("Select proteins to plot"),
               fileInput('file4', 'Select Proteins of Interest Dataframe',
                         accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
               numericInput("numconduser1","Time point 0 for Volcano:",1,min=1,max=8),
               numericInput("numconduser2","Comparison time point (based on priority):",2,min=1,max=8),
               selectInput("txtinversionuser", "Reverese Condition Order (for volcano):",
                           choices=list("Yes"="Yes","No"="No"),selected = "No"),
               selectInput("txtlabeluser", "Label User-selected Proteins:",
                           choices=list("Yes"="Yes","No"="No"),selected = "Yes"),
               numericInput("alpha5","Alpha background points",0.2,min=0,max=1),
               numericInput("alpha6","Alpha user-defined points",0.9,min=0,max=1),
               numericInput("size1","Background point size",1,min=0,max=20),
               numericInput("size2","User-defined point size",3,min=0,max=20)
             ),
             mainPanel(
               downloadButton('printplotuser', 'Download User-defined Volcano'),
               plotOutput("plotuser"),
               downloadButton('printcategviolin', 'Download User-defined Violin Plot'),
               plotOutput("categviolin"),
               plotOutput("corruserplot")
             )
             
    ),  
  
  tabPanel("Volcano Maker",
           sidebarPanel(
             tags$h3("Select proteins to plot"),
             fileInput('file6', 'Select User defined input file',
                       accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
             # numericInput("numconduser1","Time point 0 for Volcano:",1,min=1,max=8),
             # numericInput("numconduser2","Comparison time point (based on priority):",2,min=1,max=8),
             # selectInput("txtinversionuser", "Reverese Condition Order (for volcano):",
             #             choices=list("Yes"="Yes","No"="No"),selected = "No"),
             selectInput("txtlabeluser2", "Label User-selected Proteins:",
                         choices=list("Yes"="Yes","No"="No"),selected = "Yes"),
             selectInput("txtcolorlabelonly", "Color only Label points(Yes) or by Regulation Designations (No):",
                         choices=list("Label"="Label","Significance"="Significance","Localization"="Localization"),selected = "Localization"),
             selectInput("localizationuser", "Localization for Volcano Plot:",
                         choices=list("Actin binding proteins"="Actin binding proteins", "Cytoplasm"="Cytoplasm",
                                      "Endosome"="Endosome",   "ER"= "ER", "ER_high_curvature"= "ER_high_curvature",
                                      "ER-Lumen"="ER-Lumen", "ER-Membrane"="ER-Membrane", "Ergic/cisGolgi"="Ergic/cisGolgi",
                                      "Golgi"="Golgi" ,"Large Protein Complex"="Large Protein Complex", "Lysosome"="Lysosome",
                                      "Mitochondria"="Mitochondria", "Mitochondrion-IM"="Mitochondrion-IM",
                                      "Mitochondrion-Mtx"="Mitochondrion-Mtx", "Mitochondrion-OM"="Mitochondrion-OM",
                                      "Nuclear pore complex"="Nuclear pore complex", "Nucleus"="Nucleus"  ,
                                      "Peroxisome"="Peroxisome", "Plasma membrane"="Plasma membrane","Proteasome"="Proteasome",
                                      "Ribosome"="Ribosome",
                                      "Golgi.apparatus..160."="Golgi.apparatus..160." ,"Golgi.membrane..87."="Golgi.membrane..87.",
                                      "Golgi.membrane.associated..73."="Golgi.membrane.associated..73.",
                                      "Cis.Golgi..29."="Cis.Golgi..29.","Trans.Golgi..36."="Trans.Golgi..36."),selected = "Lysosome"),
             selectInput("txtpvaluelog", "-Log10(p-value) or p-value as input:",
                         choices=list("NegLogarithm"="NegLogarithm","RawValue"="RawValue"),selected = "RawValue"),
             numericInput("alpha10","Alpha background points",0.5,min=0,max=1),
             numericInput("alpha11","Alpha user-defined points",0.9,min=0,max=1),
             numericInput("size10","Background point size",1,min=0,max=20),
             numericInput("size11","User-defined point size",1.5,min=0,max=20),
             numericInput("num10","q-value cut-off",0.05,min=0,max=1),
             numericInput("num11","fold change cut-off (|abs|)",1.5,min=0),
             textInput("textxaxis1","x-axis baseline name","CTL"),
             textInput("textxaxis2","x-axis experimental name","EXP"),
             textInput("textyaxis","p-value or q-value: place p or q","p")
             # selectInput("txtlabeluser2", "Label User-selected Proteins:",
             #             choices=list("Yes"="Yes","No"="No"),selected = "Yes"),
             
             
           ),
           mainPanel(
             downloadButton('printplotvolcanomaker', 'Download User-defined Volcano'),
             plotOutput("plotvolcanomaker")
             # downloadButton('printcategviolin', 'Download User-defined Violin Plot'),
             # plotOutput("categviolin"),
             # plotOutput("corruserplot")
           )
           
  )
  )

)

# Define server logic required to draw a histogram
options(shiny.maxRequestSize=1000*1024^2)
server <- function(input, output) {
  
  
  go.df <- reactive({
    go_input <- goinput_start%>%
      dplyr::filter(organism == input$txthumanmouse)%>%
      dplyr::select(-organism)
    return(go_input)
    # print(go_input[1:50,])
    # print(unique(go_input$organism))
    # print(unique(go_input$GO_cat))
      
  })
  
  corum <- reactive({
    print(head(corum_initial))
    corum_input <- corum_initial%>%
      dplyr::filter(Organism == input$txthumanmouse)%>%
      dplyr::select(-Organism)
    return(corum_input)
    
  })
  
  locale <- reactive({
    locale33<-locale_input %>%
      dplyr::select(Gene.Symbol,Compartment,Organism)%>%
      dplyr::filter(Organism == input$txthumanmouse)%>%
      dplyr::ungroup()%>%
      dplyr::select(-Organism)%>%
      dplyr::distinct()
    return(locale33)
  })
    
    
  
  df <- reactive({
    req(input$file1)
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    tbl <- read.csv(inFile$datapath)
    
    return(tbl)
  })

  exp_design <- reactive({
    req(input$file2)
    inFile1 <- input$file2
    if (is.null(inFile1))
      return(NULL)
    tbl2 <- read.csv(inFile1$datapath)
    return(tbl2)
  })

  

  
  df_slim<-reactive({
    df2<-df()%>%
      #Reference,Gene.Symbol,Annotation,
      dplyr::select(Peptide,Protein.ID,PA.Gene.Symbol,PA.Annotation,Isolation.Specificity,Sum.Sn,
                    X126.Raw.Intensity,X127n.Raw.Intensity,X127c.Raw.Intensity,X128n.Raw.Intensity,
                    X128c.Raw.Intensity,X129n.Raw.Intensity,X129c.Raw.Intensity,X130n.Raw.Intensity,
                    X130c.Raw.Intensity,X131n.Raw.Intensity,X131c.Raw.Intensity,X132n.Raw.Intensity,
                    X132c.Raw.Intensity,X133n.Raw.Intensity,X133c.Raw.Intensity,X134n.Raw.Intensity,
                    X134c.Raw.Intensity,X135n.Raw.Intensity)%>%
      dplyr::rename(Reference = Protein.ID,
                    Gene.Symbol = PA.Gene.Symbol,
                    Annotation = PA.Annotation)%>%
      dplyr::filter(Sum.Sn > input$num101 , Isolation.Specificity > input$num102)%>%
      dplyr::filter(str_detect(Reference,input$text1004)==F)%>%
      dplyr::distinct()%>%
      dplyr::group_by(Reference,Gene.Symbol,Annotation)%>%
      dplyr::mutate(numPep = n())%>%
      dplyr::filter(numPep > input$txtNumPep-1)%>%
      dplyr::ungroup()%>%
      dplyr::select(-Sum.Sn,-Isolation.Specificity,-Peptide,-numPep)%>%
      tidyr::unite("ProtID",Reference,Gene.Symbol,Annotation, sep="__X__")%>%
      tidyr::gather("channel","intensity",2:19)%>%
      dplyr::mutate(intensity = if_else(intensity == 0, 0.001,intensity))%>%
      dplyr::inner_join(.,exp_design(),by="channel")%>%
      dplyr::filter(condition != "empty")%>%
      dplyr::mutate(cond_rep = paste(condition, replicate,sep="_"))%>%
      dplyr::group_by(ProtID,condition,replicate,cond_rep,priority)%>%
      dplyr::summarise(intensity = sum(intensity))%>%
      dplyr::filter(str_detect(ProtID,"\\#\\#")==F)%>%
      dplyr::ungroup()
    
    df2$condition<-reorder(df2$condition,df2$priority,min)
    df2$cond_rep<-reorder(df2$cond_rep,df2$priority,min)

    return(df2)
  })
  
  norm_sum_final <- reactive({
    norm_sum <- df_slim()%>%
      dplyr::ungroup()%>%
      dplyr::group_by(condition,cond_rep,priority)%>%
      dplyr::summarise(sum_int = sum(intensity))
    
    med_med <- median(norm_sum$sum_int)
    max_val <- max(norm_sum$sum_int)
    
    norm_sum<-norm_sum%>%
      dplyr::mutate(correct_factor = med_med / sum_int,
             relative_intensity = sum_int / max_val)%>%
      dplyr::select(-sum_int)
    
    print(unique(norm_sum$condition))
    return(norm_sum)
  })
  
  
  inputlater1 <- reactive({
    
    if(input$txt402 != 'All'){
      if(input$txtMedianSumOptions == "Sum"){
      norm_sum2 <- df_slim()%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")%>%
        dplyr::filter(Gene.Symbol == input$txt402)%>%
        dplyr::ungroup()%>%
        dplyr::select(condition,cond_rep,priority,intensity)%>%
        dplyr::group_by(condition,cond_rep,priority)%>%
        dplyr::summarise(sum_int = sum(intensity))
      return(norm_sum2)
      }
      
      if(input$txtMedianSumOptions == "Median"){
        norm_sum2 <- df_slim()%>%
          tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")%>%
          dplyr::filter(Gene.Symbol == input$txt402)%>%
          dplyr::ungroup()%>%
          dplyr::select(condition,cond_rep,priority,intensity)%>%
          dplyr::group_by(condition,cond_rep,priority)%>%
          dplyr::summarise(sum_int = median(intensity))
        return(norm_sum2)
      }
    }
    else{
      if(input$txtMedianSumOptions == "Sum"){
      norm_sum2 <- df_slim()%>%
        dplyr::ungroup()%>%
        dplyr::group_by(condition,cond_rep,priority)%>%
        dplyr::summarise(sum_int = sum(intensity))
      return(norm_sum2)
      }
      if(input$txtMedianSumOptions == "Median"){
        norm_sum2 <- df_slim()%>%
          dplyr::ungroup()%>%
          dplyr::group_by(condition,cond_rep,priority)%>%
          dplyr::summarise(sum_int = median(intensity))
        return(norm_sum2)
      }
    }
    
  })
  
  norm_sum_final2<-reactive({
    
    med_med <- median(inputlater1()$sum_int)
    max_val <- max(inputlater1()$sum_int)
    
    norm_sum3<-inputlater1()%>%
      dplyr::mutate(correct_factor = med_med / sum_int,
                    relative_intensity = sum_int / max_val)%>%
      dplyr::select(-sum_int)
    
    return(norm_sum3)
  })
  
  df_final_res <- reactive({
    df_final<-full_join(df_slim(), norm_sum_final2(), by=c("condition","cond_rep","priority"))
    df_final<-df_final%>%
      dplyr::mutate(norm_int = intensity * correct_factor,
             log2_int = log2(norm_int))%>%
      dplyr::ungroup()%>%
      dplyr::mutate(log2_int=if_else(log2_int<=0,0,log2_int))
    return(df_final)
  })
  
  df_cv <- reactive({
    df_final_cv<-df_final_res()%>%
      dplyr::ungroup()%>%
      group_by(ProtID,condition,priority)%>%
      summarise(cv = sd(log2_int)/mean(log2_int))
    return(df_final_cv)
  })
  output$medianC1<-renderText({
    val<-unique(df_cv()$condition)[1]
    df_c1_cv<-df_cv()%>%
      dplyr::filter(condition == val,!is.na(cv))
    paste(unique(df_cv()$condition)[1],"=",round(median(df_c1_cv$cv,na.rm = T),4),sep = " ")
  })
  output$medianC2<-renderText({
    val2<-unique(df_cv()$condition)[2]
    df_c2_cv<-df_cv()%>%
      dplyr::filter(condition == val2,!is.na(cv))
    paste(unique(df_cv()$condition)[2],"=",round(median(df_c2_cv$cv,na.rm = T),4),sep = " ")
  })
  output$medianC3<-renderText({
    val3<-unique(df_cv()$condition)[3]
    df_c2_cv<-df_cv()%>%
      dplyr::filter(condition == val3,!is.na(cv))
    paste(unique(df_cv()$condition)[3],"=",round(median(df_c2_cv$cv,na.rm = T),4),sep = " ")
  })
  output$medianC4<-renderText({
    val4<-unique(df_cv()$condition)[4]
    df_c2_cv<-df_cv()%>%
      dplyr::filter(condition == val4,!is.na(cv))
    paste(unique(df_cv()$condition)[4],"=",round(median(df_c2_cv$cv,na.rm = T),4),sep = " ")
  })
  
  output$plotCV <- renderPlot({
    #reorder(condition,priority,min)
    ggplot(df_cv(),aes(condition,cv,color=condition))+geom_violin(draw_quantiles = c(0.25,0.5,0.75))+theme_classic()+
      labs(x="Condition",y="Coefficient of Variation",title="CV Plot")+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=10),plot.title = element_text(size=20, hjust = 0.5))+
      scale_colour_viridis_d(end=0.8)
  })
  
  output$plot1 <- renderPlot({
    ggplot(df_slim(),aes(log2(intensity),color=cond_rep))+geom_density()+theme_classic()+
      labs(x=bquote(""~Log[2]~"(Intensity)"),title="Pre-normalization")+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),plot.title = element_text(size=20, hjust = 0.5))
  })
  
  output$plot2 <- renderPlot({
    ggplot(df_final_res(),aes(log2(norm_int),color=cond_rep))+geom_density()+theme_classic()+
      labs(x=bquote(""~Log[2]~"(Intensity)"),title="Post-normalization")+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),plot.title = element_text(size=20, hjust = 0.5))
  })

  
  labelCheckPlotInput <- reactive({
    p1<-ggplot(norm_sum_final2(),aes(cond_rep,relative_intensity,fill=condition))+
      geom_col()+theme_classic()+scale_fill_viridis_d(end=0.8)+
      labs(x="",y="Relative Intensity",title="Ratio Check")+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),plot.title = element_text(size=20, hjust = 0.5),legend.position = "none")
  })
  
  output$plot3 <- renderPlot({
    #reorder(cond_rep,priority,min)
    ggplot(norm_sum_final2(),aes(cond_rep,relative_intensity,fill=condition))+
      geom_col()+theme_classic()+scale_fill_viridis_d(end=0.8)+
      labs(x="",y="Relative Intensity",title="Ratio Check")+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=10),plot.title = element_text(size=20, hjust = 0.5),legend.position = "none")
  })
  
  # output$labelCheckDownloadPlot <- downloadHandler(
  #   filename = function() { paste("labelCheckPlot", '.png', sep='') },
  #   content = function(file) {
  #     device <- function(..., width, height) grDevices::png(..., width = 14, height = 6.5, res = 300, units = "in")
  #     ggsave(file, plot = labelCheckPlotInput(), device = device)
  #   }
  # )
  
  output$labelCheckDownloadPlot <- downloadHandler(
    filename = function() { paste("labelCheckPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = labelCheckPlotInput(),height = 6.5,width=14)
    }
  )

  
  correlation_plot2<-reactive({
    correlation_plot <- df_final_res()%>%
      #dplyr::select(-intensity,-correct_factor,-relative_intensity,-norm_int)
      dplyr::select(ProtID , cond_rep , log2_int)%>%
      tidyr::spread(cond_rep,log2_int)
  })
  
  output$totalproteins <- renderText({
    paste(dim(correlation_plot2())[[1]]," total proteins analyzed")
  })
  
  output$plot4 <- renderPlot({
    ggpairs(correlation_plot2()[2:dim(correlation_plot2())[[2]]],lower = list(continuous = wrap("points",size=0.1)))+theme_classic()
  }, height = 1000, width = 1000 )
  
  
  
  
  
  
  
  output$plotPCA <- renderPlot({
    
    df_pca<-correlation_plot2()%>%
      #tibble::remove_rownames %>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(df_pca))
    
    df_pca_use<-data.frame(example_pca_test$x)%>%
      tibble::rownames_to_column("sample")

    PoV <- example_pca_test$sdev^2/sum(example_pca_test$sdev^2)
    pc_val1<-input$txtXpca
    pc_val2<-input$txtYpca

    pc_val1_num <- as.numeric(as.character(str_remove(pc_val1,"PC")))

    pc_val2_num <- as.numeric(as.character(str_remove(pc_val2,"PC")))
    
    pc_val1_num_plot<-pc_val1_num+1
    pc_val2_num_plot<-pc_val2_num+1

    var_pcVal1<-round(PoV[pc_val1_num]*100,4)
    var_pcVal2<-round(PoV[pc_val2_num]*100,4)
    
    df_pca_use<-df_pca_use%>%
               dplyr::mutate(loc_last = stringi::stri_locate_last_fixed(sample,"_")[2]-1,
                      condition = stringr::str_sub(sample, start=1,end=-3))

    
    #print(df_pca_use)

    ggplot(df_pca_use,aes(df_pca_use[,pc_val1_num_plot],df_pca_use[,pc_val2_num_plot],color=condition,label=sample))+
      geom_point(size=3,alpha=0.5)+
      geom_text_repel(max.overlaps = Inf,box.padding=1,show.legend=FALSE)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      labs(x=paste(input$txtXpca,": ", var_pcVal1,"%",sep=""),
           y=paste(input$txtYpca,": ", var_pcVal2,"%",sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.position = "none")
  })
  
  
  output$PCAdrivers <- renderReactable({
    correlation_plot<-correlation_plot2()%>%
      #tibble::remove_rownames %>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))

    var <- factoextra::get_pca_var(example_pca_test)
    contrib <- as.data.frame(var$contrib)

    contrib<-contrib%>%
      tibble::rownames_to_column("ProtID")%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    ##choose PC for drivers
    pc_val_drive<-input$txtZpca
    pc_val_drive1 <- as.numeric(as.character(str_remove(pc_val_drive,"PC")))
    pc_val_drive_sel<-all_of(pc_val_drive1) +3

    sel_column<-contrib%>%
      dplyr::select(1:3, all_of(pc_val_drive_sel))

    sel_column<-sel_column%>%
      plotly::arrange(desc(.[,4]))

    colnames(sel_column) <- sub("Dim\\.", "PC", colnames(sel_column))
    #sel_column<- sel_column[1:100,]
    reactable(sel_column, filterable = TRUE)
  })
  
  driverAll<- reactive({
    correlation_plot<-correlation_plot2()%>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))
    
    var <- factoextra::get_pca_var(example_pca_test)
    contrib <- as.data.frame(var$contrib)
    
    contrib<-contrib%>%
      tibble::rownames_to_column("ProtID")%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    ##choose PC for drivers
    # pc_val_drive<-input$txtZpca
    # pc_val_drive1 <- as.numeric(as.character(str_remove(pc_val_drive,"PC")))
    # pc_val_drive_sel<-pc_val_drive1 +3
    
    # sel_column<-contrib%>%
    #   dplyr::select(1:3, pc_val_drive_sel)
    
    sel_column<-contrib%>%
      plotly::arrange(desc(.[,4]))
    
    colnames(sel_column) <- sub("Dim\\.", "PC", colnames(sel_column))
    #sel_column<- sel_column[1:100,]
    return(sel_column)
  })
  
  
  
  output$driversDownload <- downloadHandler(
    filename = function() {
      paste("PCA_protein_weights_all", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(driverAll(), file, row.names = FALSE)
    }
  )
  
  
  output$plotPCAadd <-renderPlot({
    correlation_plot<-correlation_plot2()%>%
      #tibble::remove_rownames %>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))
    factoextra::fviz_eig(example_pca_test, addlabels = TRUE) + theme_classic()
  })
  
  
  pcaPlotInput <- reactive({
    correlation_plot<-correlation_plot2()%>%
      #tibble::remove_rownames %>%
      tibble::column_to_rownames(var="ProtID")
    example_pca_test<-prcomp(t(correlation_plot))

    df_pca_use<-data.frame(example_pca_test$x)%>%
      tibble::rownames_to_column("sample")

    PoV <- example_pca_test$sdev^2/sum(example_pca_test$sdev^2)
    pc_val1<-input$txtXpca
    pc_val2<-input$txtYpca

    pc_val1_num <- as.numeric(as.character(str_remove(pc_val1,"PC")))

    pc_val2_num <- as.numeric(as.character(str_remove(pc_val2,"PC")))

    pc_val1_num_plot<-pc_val1_num+1
    pc_val2_num_plot<-pc_val2_num+1

    var_pcVal1<-round(PoV[pc_val1_num]*100,4)
    var_pcVal2<-round(PoV[pc_val2_num]*100,4)

    df_pca_use<-df_pca_use%>%
      dplyr::mutate(loc_last = stringi::stri_locate_last_fixed(sample,"_")[2]-1,
                    condition = stringr::str_sub(sample, start=1,end=loc_last))

    p70<-ggplot(df_pca_use,aes(df_pca_use[,pc_val1_num_plot],df_pca_use[,pc_val2_num_plot],color=condition,label=sample))+
      geom_point(size=3,alpha=0.5)+
      geom_text_repel(max.overlaps = Inf,box.padding=1,show.legend=FALSE)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      labs(x=paste(input$txtXpca,": ", var_pcVal1,"%",sep=""),
           y=paste(input$txtYpca,": ", var_pcVal2,"%",sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size=14),legend.position = "none")

  })
  
  
  # output$PCADownloadPlot <- downloadHandler(
  #   filename = function() { paste("PCA_Plot", '.png', sep='') },
  #   content = function(file) {
  #     device <- function(..., width, height) grDevices::png(..., width = 10, height = 10, res = 300, units = "in")
  #     ggsave(file, plot = pcaPlotInput(), device = device)
  #   }
  # )
  
  output$PCADownloadPlot <- downloadHandler(
    filename = function() { paste("PCA_Plot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = pcaPlotInput(), width = 10, height = 10)
    }
  )
  
  corrPlotInput <- reactive({
    p72<-ggpairs(correlation_plot2()[2:dim(correlation_plot2())[[2]]],lower = list(continuous = wrap("points",size=0.1)))+theme_classic()
  })
  
  # output$CorrDownloadPlot <- downloadHandler(
  #   filename = function() { paste("Rep_Corr_Plot", '.png', sep='') },
  #   content = function(file) {
  #     device <- function(..., width, height) grDevices::png(..., width = 12, height = 12, res = 300, units = "in")
  #     ggsave(file, plot = corrPlotInput(), device = device)
  #   }
  # )
  
  output$CorrDownloadPlot <- downloadHandler(
    filename = function() { paste("Rep_Corr_Plot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrPlotInput(), width = 12, height = 12)
    }
  )
  
  

  ttest_df_initial1 <- reactive({

    ttest_df<-df_final_res()%>%
      dplyr::select(-intensity,-correct_factor,-relative_intensity,-norm_int)%>%
      dplyr::group_by(ProtID)%>%
      tidyr::nest()
    
    ttest_df<-ttest_df%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    if(input$txt402 != 'All'){
      ttest_df<-ttest_df%>%filter(Gene.Symbol != input$txt402)
    }
    
    return(ttest_df)
  })
  
  ttest_df_initial2 <- reactive({
    
    # t.test_fun<-function(data){
    #   df<-data.frame(data)
    #   c1<-df%>%dplyr::filter(condition == val1)
    #   c2<-df%>%dplyr::filter(condition == val2)
    #   test<-t.test(c1$log2_int,c2$log2_int, var.equal = T)
    #   return(test$p.value)
    # }
    
    t.test_fun<-function(data){
      df<-data.frame(data)
      # c1<-df%>%dplyr::filter(condition == val1)
      # c2<-df%>%dplyr::filter(condition == val2)
      # test<-t.test(c1$log2_int,c2$log2_int, var.equal = T)
      return(t.test(subset(df$log2_int,df$condition == val1),subset(df$log2_int,df$condition == val2), var.equal = T)$p.value)
    }
    
    
    # fc_fun<-function(data){
    #   df<-data.frame(data)
    #   c1<-df%>%dplyr::filter(condition == val1)
    #   c2<-df%>%dplyr::filter(condition == val2)
    #   test<-t.test(c1$log2_int,c2$log2_int, var.equal = T)
    #   return(test$estimate[[2]] - test$estimate[[1]]  )
    # }
    
    fc_fun<-function(data){
      df<-data.frame(data)
      test<-t.test(subset(df$log2_int,df$condition == val1),subset(df$log2_int,df$condition == val2), var.equal = T)
      return(test$estimate[[2]] - test$estimate[[1]]  )
    }
    
    
    
    
    anova_test<-function(data){
      df<-data.frame(data)
      test<-aov(log2_int ~ condition, data= df)
      return(summary(test)[[1]][["Pr(>F)"]][1])
    }
    
    # cond1_fun<-function(data){
    #   df<-data.frame(data)
    #   c1<-df%>%dplyr::filter(condition == value)
    #   return(mean(c1$log2_int)  )
    # }
    
    cond1_fun<-function(data){
      df<-data.frame(data)
      # c1<-df%>%dplyr::filter(condition == value)
      return(mean(subset(df$log2_int,df$condition == value))  )
    }
    
    if (length(unique(norm_sum_final()$condition)) >= 3) {
      print("three")
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      name1<-paste("p.val_",con2,"-",con1,sep="")
      qname1<-paste("q.val_",con2,"-",con1,sep="")
      fc1<-paste("log2FC_",con2,"-",con1,sep="")
      sig1<-paste("sig_",con2,"-",con1,sep="")
      value <- con1
      val1 <-  con1
      val2 <-  con2
      
      ttest_df<-ttest_df_initial1()%>%
        dplyr::mutate(p.val1 = data%>%purrr::map_dbl(t.test_fun),
                      !!fc1:= data%>%purrr::map_dbl(fc_fun),
                      !!con1:= data%>% purrr::map_dbl(cond1_fun))
      ttest_df$qval1<-p.adjust(ttest_df$p.val1,"fdr")
      
      ttest_df<-ttest_df%>%
        dplyr::mutate(!!sig1:= if_else(qval1 < input$num77,if_else(abs(.data[[fc1]])>log2(input$num78), if_else(.data[[fc1]] > 0, "up","down"),"n.s."),"n.s."))
      
      names(ttest_df)[names(ttest_df) == "p.val1"] <- name1
      names(ttest_df)[names(ttest_df) == "qval1"]  <- qname1
      
      name2<-paste("p.val_",con3,"-",con1,sep="")
      qname2<-paste("q.val_",con3,"-",con1,sep="")
      fc2<-paste("log2FC_",con3,"-",con1,sep="")
      sig2<-paste("sig_",con3,"-",con1,sep="")
      value <- con2
      val2 <-  con3
      
      ttest_df<-ttest_df%>%
        dplyr::mutate(p.val2 = data%>%purrr::map_dbl(t.test_fun),
                      !!fc2:= data%>%purrr::map_dbl(fc_fun),
                      !!con2:= data%>% purrr::map_dbl(cond1_fun))
      ttest_df$qval2<-p.adjust(ttest_df$p.val2,"fdr")
      
      ttest_df<-ttest_df%>%
        dplyr::mutate(!!sig2:= if_else(qval2 < input$num77,if_else(abs(.data[[fc2]])>log2(input$num78), if_else(.data[[fc2]] > 0, "up","down"),"n.s."),"n.s."))
      
      names(ttest_df)[names(ttest_df) == "p.val2"] <- name2
      names(ttest_df)[names(ttest_df) == "qval2"]  <- qname2
      
      name3<-paste("p.val_",con3,"-",con2,sep="")
      qname3<-paste("q.val_",con3,"-",con2,sep="")
      fc3<-paste("log2FC_",con3,"-",con2,sep="")
      sig3<-paste("sig_",con3,"-",con2,sep="")
      value <- con3
      val1 <-  con2

      
      ttest_df<-ttest_df%>%
        dplyr::mutate(!!con3:= data%>% purrr::map_dbl(cond1_fun))
      
      if (input$txtcomparettest == "All") {
        ttest_df<-ttest_df%>%
          dplyr::mutate(p.val3 = data%>%purrr::map_dbl(t.test_fun),
                        !!fc3:= data%>%purrr::map_dbl(fc_fun))
        ttest_df$qval3<-p.adjust(ttest_df$p.val3,"fdr")
        
        ttest_df<-ttest_df%>%
          dplyr::mutate(!!sig3:= if_else(qval3 < input$num77,if_else(abs(.data[[fc3]])>log2(input$num78), if_else(.data[[fc3]] > 0, "up","down"),"n.s."),"n.s."))
        
        names(ttest_df)[names(ttest_df) == "p.val3"] <- name3
        names(ttest_df)[names(ttest_df) == "qval3"]  <- qname3
      }
      
      
      ttest_df<-ttest_df%>%
        dplyr::mutate(p.val_anova= data %>% map_dbl(anova_test))
      
      ttest_df$q.val_anova <- p.adjust(ttest_df$p.val_anova,method="BH")
      
      
      ttest_df<-ttest_df%>%
        mutate(anova_sig = if_else(q.val_anova < input$num77, "sig","n.s."))
      
      if (length(unique(norm_sum_final()$condition)) >= 4) {
        print("four")
        con4<-as.character(unique(norm_sum_final()$condition)[4])
        
        name4<-paste("p.val_",con4,"-",con1,sep="")
        qname4<-paste("q.val_",con4,"-",con1,sep="")
        fc4<-paste("log2FC_",con4,"-",con1,sep="")
        sig4<-paste("sig_",con4,"-",con1,sep="")
        value <- con4
        val1 <-  con1
        val2 <-  con4
        ttest_df<-ttest_df%>%
          dplyr::mutate(p.val4 = data%>%purrr::map_dbl(t.test_fun),
                        !!fc4:= data%>%purrr::map_dbl(fc_fun),
                        !!con4:= data%>% purrr::map_dbl(cond1_fun))
        ttest_df$qval4<-p.adjust(ttest_df$p.val4,"fdr")
        
        ttest_df<-ttest_df%>%
          dplyr::mutate(!!sig4:= if_else(qval4 < input$num77,if_else(abs(.data[[fc4]])>log2(input$num78), if_else(.data[[fc4]] > 0, "up","down"),"n.s."),"n.s."))
        
        names(ttest_df)[names(ttest_df) == "p.val4"] <- name4
        names(ttest_df)[names(ttest_df) == "qval4"]  <- qname4
        
        

        if (input$txtcomparettest == "All") {
          name5<-paste("p.val_",con4,"-",con2,sep="")
          qname5<-paste("q.val_",con4,"-",con2,sep="")
          fc5<-paste("log2FC_",con4,"-",con2,sep="")
          sig5<-paste("sig_",con4,"-",con2,sep="")
          val1 <-  con2
          val2 <-  con4
          ttest_df<-ttest_df%>%
            dplyr::mutate(p.val5 = data%>%purrr::map_dbl(t.test_fun),
                          !!fc5:= data%>%purrr::map_dbl(fc_fun))
          ttest_df$qval5<-p.adjust(ttest_df$p.val5,"fdr")
          
          ttest_df<-ttest_df%>%
            dplyr::mutate(!!sig5:= if_else(qval5 < input$num77,if_else(abs(.data[[fc5]])>log2(input$num78), if_else(.data[[fc5]] > 0, "up","down"),"n.s."),"n.s."))
          
          names(ttest_df)[names(ttest_df) == "p.val5"] <- name5
          names(ttest_df)[names(ttest_df) == "qval5"]  <- qname5
          
          name6<-paste("p.val_",con4,"-",con3,sep="")
          qname6<-paste("q.val_",con4,"-",con3,sep="")
          fc6<-paste("log2FC_",con4,"-",con3,sep="")
          sig6<-paste("sig_",con4,"-",con3,sep="")
          val1 <-  con3
          val2 <-  con4
          ttest_df<-ttest_df%>%
            dplyr::mutate(p.val6 = data%>%purrr::map_dbl(t.test_fun),
                          !!fc6:= data%>%purrr::map_dbl(fc_fun))
          ttest_df$qval6<-p.adjust(ttest_df$p.val6,"fdr")
          
          ttest_df<-ttest_df%>%
            dplyr::mutate(!!sig6:= if_else(qval6 < input$num77,if_else(abs(.data[[fc6]])>log2(input$num78), if_else(.data[[fc6]] > 0, "up","down"),"n.s."),"n.s."))
          
          names(ttest_df)[names(ttest_df) == "p.val6"] <- name6
          names(ttest_df)[names(ttest_df) == "qval6"]  <- qname6
        }
        

        if (length(unique(norm_sum_final()$condition)) >= 5) {
          print("five")
          con5<-as.character(unique(norm_sum_final()$condition)[5])
          
          name7<-paste("p.val_",con5,"-",con1,sep="")
          qname7<-paste("q.val_",con5,"-",con1,sep="")
          fc7<-paste("log2FC_",con5,"-",con1,sep="")
          sig7<-paste("sig_",con5,"-",con1,sep="")
          value <- con5
          val1 <-  con1
          val2 <-  con5
          ttest_df<-ttest_df%>%
            dplyr::mutate(p.val7 = data%>%purrr::map_dbl(t.test_fun),
                          !!fc7:= data%>%purrr::map_dbl(fc_fun),
                          !!con5:= data%>% purrr::map_dbl(cond1_fun))
          
          ttest_df$qval7<-p.adjust(ttest_df$p.val7,"fdr")
          
          ttest_df<-ttest_df%>%
            dplyr::mutate(!!sig7:= if_else(qval7 < input$num77,if_else(abs(.data[[fc7]])>log2(input$num78), if_else(.data[[fc7]] > 0, "up","down"),"n.s."),"n.s."))
          
          names(ttest_df)[names(ttest_df) == "p.val7"] <- name7
          names(ttest_df)[names(ttest_df) == "qval7"]  <- qname7
          
          
          if (input$txtcomparettest == "All") {
            name8<-paste("p.val_",con5,"-",con2,sep="")
            qname8<-paste("q.val_",con5,"-",con2,sep="")
            fc8<-paste("log2FC_",con5,"-",con2,sep="")
            sig8<-paste("sig_",con5,"-",con2,sep="")
            val1 <-  con2
            val2 <-  con5
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val8 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc8:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval8<-p.adjust(ttest_df$p.val8,"fdr")
            
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig8:= if_else(qval8 < input$num77,if_else(abs(.data[[fc8]])>log2(input$num78), if_else(.data[[fc8]] > 0, "up","down"),"n.s."),"n.s."))
            
            names(ttest_df)[names(ttest_df) == "p.val8"] <- name8
            names(ttest_df)[names(ttest_df) == "qval8"]  <- qname8
            
            name9<-paste("p.val_",con5,"-",con3,sep="")
            qname9<-paste("q.val_",con5,"-",con3,sep="")
            fc9<-paste("log2FC_",con5,"-",con3,sep="")
            sig9<-paste("sig_",con5,"-",con3,sep="")
            val1 <-  con3
            val2 <-  con5
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val9 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc9:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval9<-p.adjust(ttest_df$p.val9,"fdr")
            
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig9:= if_else(qval9 < input$num77,if_else(abs(.data[[fc9]])>log2(input$num78), if_else(.data[[fc9]] > 0, "up","down"),"n.s."),"n.s."))
            
            names(ttest_df)[names(ttest_df) == "p.val9"] <- name9
            names(ttest_df)[names(ttest_df) == "qval9"]  <- qname9
            
            name10<-paste("p.val_",con5,"-",con4,sep="")
            qname10<-paste("q.val_",con5,"-",con4,sep="")
            fc10<-paste("log2FC_",con5,"-",con4,sep="")
            sig10<-paste("sig_",con5,"-",con4,sep="")
            val1 <-  con4
            val2 <-  con5
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val10 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc10:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval10<-p.adjust(ttest_df$p.val10,"fdr")
            
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig10:= if_else(qval10 < input$num77,if_else(abs(.data[[fc10]])>log2(input$num78), if_else(.data[[fc10]] > 0, "up","down"),"n.s."),"n.s."))
            
            names(ttest_df)[names(ttest_df) == "p.val10"] <- name10
            names(ttest_df)[names(ttest_df) == "qval10"]  <- qname10
          }
          

          
        
        
        if (length(unique(norm_sum_final()$condition)) >= 6) {
          print("six")
          con6<-as.character(unique(norm_sum_final()$condition)[6])
          
          name11<-paste("p.val_",con6,"-",con1,sep="")
          qname11<-paste("q.val_",con6,"-",con1,sep="")
          fc11<-paste("log2FC_",con6,"-",con1,sep="")
          sig11<-paste("sig_",con6,"-",con1,sep="")
          value <- con6
          val1 <-  con1
          val2 <-  con6
          ttest_df<-ttest_df%>%
            dplyr::mutate(p.val11 = data%>%purrr::map_dbl(t.test_fun),
                          !!fc11:= data%>%purrr::map_dbl(fc_fun),
                          !!con6:= data%>% purrr::map_dbl(cond1_fun))
          
          ttest_df$qval11<-p.adjust(ttest_df$p.val11,"fdr")
          
          ttest_df<-ttest_df%>%
            dplyr::mutate(!!sig11:= if_else(qval11 < input$num77,if_else(abs(.data[[fc11]])>log2(input$num78), if_else(.data[[fc11]] > 0, "up","down"),"n.s."),"n.s."))
          
          names(ttest_df)[names(ttest_df) == "p.val11"] <- name11
          names(ttest_df)[names(ttest_df) == "qval11"]  <- qname11
          
          
          if (input$txtcomparettest == "All") {
            name12<-paste("p.val_",con6,"-",con2,sep="")
            qname12<-paste("q.val_",con6,"-",con2,sep="")
            fc12<-paste("log2FC_",con6,"-",con2,sep="")
            sig12<-paste("sig_",con6,"-",con2,sep="")
            val1 <-  con2
            val2 <-  con6
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val12 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc12:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval12<-p.adjust(ttest_df$p.val12,"fdr")
            
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig12:= if_else(qval12 < input$num77,if_else(abs(.data[[fc12]])>log2(input$num78), if_else(.data[[fc12]] > 0, "up","down"),"n.s."),"n.s."))
            names(ttest_df)[names(ttest_df) == "p.val12"] <- name12
            names(ttest_df)[names(ttest_df) == "qval12"]  <- qname12
            
            name13<-paste("p.val_",con6,"-",con3,sep="")
            qname13<-paste("q.val_",con6,"-",con3,sep="")
            fc13<-paste("log2FC_",con6,"-",con3,sep="")
            sig13<-paste("sig_",con6,"-",con3,sep="")
            val1 <-  con3
            val2 <-  con6
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val13 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc13:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval13<-p.adjust(ttest_df$p.val13,"fdr")
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig13:= if_else(qval13 < input$num77,if_else(abs(.data[[fc13]])>log2(input$num78), if_else(.data[[fc13]] > 0, "up","down"),"n.s."),"n.s."))
            names(ttest_df)[names(ttest_df) == "p.val13"] <- name13
            names(ttest_df)[names(ttest_df) == "qval13"]  <- qname13
            
            name14<-paste("p.val_",con6,"-",con4,sep="")
            qname14<-paste("q.val_",con6,"-",con4,sep="")
            fc14<-paste("log2FC_",con6,"-",con4,sep="")
            sig14<-paste("sig_",con6,"-",con4,sep="")
            val1 <-  con4
            val2 <-  con6
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val14 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc14:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval14<-p.adjust(ttest_df$p.val14,"fdr")
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig14:= if_else(qval14 < input$num77,if_else(abs(.data[[fc14]])>log2(input$num78), if_else(.data[[fc14]] > 0, "up","down"),"n.s."),"n.s."))
            names(ttest_df)[names(ttest_df) == "p.val14"] <- name14
            names(ttest_df)[names(ttest_df) == "qval14"]  <- qname14
            
            name15<-paste("p.val_",con6,"-",con5,sep="")
            qname15<-paste("q.val_",con6,"-",con5,sep="")
            fc15<-paste("log2FC_",con6,"-",con5,sep="")
            sig15<-paste("sig_",con6,"-",con5,sep="")
            val1 <-  con5
            val2 <-  con6
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val15 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc15:= data%>%purrr::map_dbl(fc_fun))
            
            ttest_df$qval15<-p.adjust(ttest_df$p.val15,"fdr")
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig15:= if_else(qval15 < input$num77,if_else(abs(.data[[fc15]])>log2(input$num78), if_else(.data[[fc15]] > 0, "up","down"),"n.s."),"n.s."))
            names(ttest_df)[names(ttest_df) == "p.val15"] <- name15
            names(ttest_df)[names(ttest_df) == "qval15"]  <- qname15
          }
          
          if (length(unique(norm_sum_final()$condition)) >= 7) {
            print("seven")
            con7<-as.character(unique(norm_sum_final()$condition)[7])
            
            name16<-paste("p.val_",con7,"-",con1,sep="")
            qname16<-paste("q.val_",con7,"-",con1,sep="")
            fc16<-paste("log2FC_",con7,"-",con1,sep="")
            sig16<-paste("sig_",con7,"-",con1,sep="")
            value <- con7
            val1 <-  con1
            val2 <-  con7
            ttest_df<-ttest_df%>%
              dplyr::mutate(p.val16 = data%>%purrr::map_dbl(t.test_fun),
                            !!fc16:= data%>%purrr::map_dbl(fc_fun),
                            !!con7:= data%>% purrr::map_dbl(cond1_fun))
            
            ttest_df$qval16<-p.adjust(ttest_df$p.val16,"fdr")
            
            ttest_df<-ttest_df%>%
              dplyr::mutate(!!sig16:= if_else(qval16 < input$num77,if_else(abs(.data[[fc16]])>log2(input$num78), if_else(.data[[fc16]] > 0, "up","down"),"n.s."),"n.s."))
            
            names(ttest_df)[names(ttest_df) == "p.val16"] <- name16
            names(ttest_df)[names(ttest_df) == "qval16"]  <- qname16
            
            
            if (input$txtcomparettest == "All") {
              name17<-paste("p.val_",con7,"-",con2,sep="")
              qname17<-paste("q.val_",con7,"-",con2,sep="")
              fc17<-paste("log2FC_",con7,"-",con2,sep="")
              sig17<-paste("sig_",con7,"-",con2,sep="")
              val1 <-  con2
              val2 <-  con7
              ttest_df<-ttest_df%>%
                dplyr::mutate(p.val17 = data%>%purrr::map_dbl(t.test_fun),
                              !!fc17:= data%>%purrr::map_dbl(fc_fun))
              
              ttest_df$qval17<-p.adjust(ttest_df$p.val17,"fdr")
              
              ttest_df<-ttest_df%>%
                dplyr::mutate(!!sig17:= if_else(qval17 < input$num77,if_else(abs(.data[[fc17]])>log2(input$num78), if_else(.data[[fc17]] > 0, "up","down"),"n.s."),"n.s."))
              names(ttest_df)[names(ttest_df) == "p.val17"] <- name17
              names(ttest_df)[names(ttest_df) == "qval17"]  <- qname17
              
              name18<-paste("p.val_",con7,"-",con3,sep="")
              qname18<-paste("q.val_",con7,"-",con3,sep="")
              fc18<-paste("log2FC_",con7,"-",con3,sep="")
              sig18<-paste("sig_",con7,"-",con3,sep="")
              val1 <-  con3
              val2 <-  con7
              ttest_df<-ttest_df%>%
                dplyr::mutate(p.val18 = data%>%purrr::map_dbl(t.test_fun),
                              !!fc18:= data%>%purrr::map_dbl(fc_fun))
              
              ttest_df$qval18<-p.adjust(ttest_df$p.val18,"fdr")
              ttest_df<-ttest_df%>%
                dplyr::mutate(!!sig18:= if_else(qval18 < input$num77,if_else(abs(.data[[fc18]])>log2(input$num78), if_else(.data[[fc18]] > 0, "up","down"),"n.s."),"n.s."))
              names(ttest_df)[names(ttest_df) == "p.val18"] <- name18
              names(ttest_df)[names(ttest_df) == "qval18"]  <- qname18
              
              name19<-paste("p.val_",con7,"-",con4,sep="")
              qname19<-paste("q.val_",con7,"-",con4,sep="")
              fc19<-paste("log2FC_",con7,"-",con4,sep="")
              sig19<-paste("sig_",con7,"-",con4,sep="")
              val1 <-  con4
              val2 <-  con7
              ttest_df<-ttest_df%>%
                dplyr::mutate(p.val19 = data%>%purrr::map_dbl(t.test_fun),
                              !!fc19:= data%>%purrr::map_dbl(fc_fun))
              
              ttest_df$qval19<-p.adjust(ttest_df$p.val19,"fdr")
              ttest_df<-ttest_df%>%
                dplyr::mutate(!!sig19:= if_else(qval19 < input$num77,if_else(abs(.data[[fc19]])>log2(input$num78), if_else(.data[[fc19]] > 0, "up","down"),"n.s."),"n.s."))
              names(ttest_df)[names(ttest_df) == "p.val19"] <- name19
              names(ttest_df)[names(ttest_df) == "qval19"]  <- qname19
              
              name20<-paste("p.val_",con7,"-",con5,sep="")
              qname20<-paste("q.val_",con7,"-",con5,sep="")
              fc20<-paste("log2FC_",con7,"-",con5,sep="")
              sig20<-paste("sig_",con7,"-",con5,sep="")
              val1 <-  con5
              val2 <-  con7
              ttest_df<-ttest_df%>%
                dplyr::mutate(p.val20 = data%>%purrr::map_dbl(t.test_fun),
                              !!fc20:= data%>%purrr::map_dbl(fc_fun))
              
              ttest_df$qval20<-p.adjust(ttest_df$p.val20,"fdr")
              ttest_df<-ttest_df%>%
                dplyr::mutate(!!sig20:= if_else(qval20 < input$num77,if_else(abs(.data[[fc20]])>log2(input$num78), if_else(.data[[fc20]] > 0, "up","down"),"n.s."),"n.s."))
              names(ttest_df)[names(ttest_df) == "p.val20"] <- name20
              names(ttest_df)[names(ttest_df) == "qval20"]  <- qname20
              
              
              name21<-paste("p.val_",con7,"-",con6,sep="")
              qname21<-paste("q.val_",con7,"-",con6,sep="")
              fc21<-paste("log2FC_",con7,"-",con6,sep="")
              sig21<-paste("sig_",con7,"-",con6,sep="")
              val1 <-  con6
              val2 <-  con7
              ttest_df<-ttest_df%>%
                dplyr::mutate(p.val21 = data%>%purrr::map_dbl(t.test_fun),
                              !!fc21:= data%>%purrr::map_dbl(fc_fun))
              
              ttest_df$qval21<-p.adjust(ttest_df$p.val21,"fdr")
              ttest_df<-ttest_df%>%
                dplyr::mutate(!!sig21:= if_else(qval21 < input$num77,if_else(abs(.data[[fc21]])>log2(input$num78), if_else(.data[[fc21]] > 0, "up","down"),"n.s."),"n.s."))
              names(ttest_df)[names(ttest_df) == "p.val21"] <- name21
              names(ttest_df)[names(ttest_df) == "qval21"]  <- qname21
              
            }
            
            
            if (length(unique(norm_sum_final()$condition)) >= 8) {
              print("eight")
              con8<-as.character(unique(norm_sum_final()$condition)[8])
              
              name22<-paste("p.val_",con8,"-",con1,sep="")
              qname22<-paste("q.val_",con8,"-",con1,sep="")
              fc22<-paste("log2FC_",con8,"-",con1,sep="")
              sig22<-paste("sig_",con8,"-",con1,sep="")
              value <- con8
              val1 <-  con1
              val2 <-  con8
              ttest_df<-ttest_df%>%
                dplyr::mutate(p.val22 = data%>%purrr::map_dbl(t.test_fun),
                              !!fc22:= data%>%purrr::map_dbl(fc_fun),
                              !!con8:= data%>% purrr::map_dbl(cond1_fun))
              
              ttest_df$qval22<-p.adjust(ttest_df$p.val22,"fdr")
              
              ttest_df<-ttest_df%>%
                dplyr::mutate(!!sig22:= if_else(qval22 < input$num77,if_else(abs(.data[[fc22]])>log2(input$num78), if_else(.data[[fc22]] > 0, "up","down"),"n.s."),"n.s."))
              
              names(ttest_df)[names(ttest_df) == "p.val22"] <- name22
              names(ttest_df)[names(ttest_df) == "qval22"]  <- qname22
              
              
              if (input$txtcomparettest == "All") {
                name23<-paste("p.val_",con8,"-",con2,sep="")
                qname23<-paste("q.val_",con8,"-",con2,sep="")
                fc23<-paste("log2FC_",con8,"-",con2,sep="")
                sig23<-paste("sig_",con8,"-",con2,sep="")
                val1 <-  con2
                val2 <-  con8
                ttest_df<-ttest_df%>%
                  dplyr::mutate(p.val23 = data%>%purrr::map_dbl(t.test_fun),
                                !!fc23:= data%>%purrr::map_dbl(fc_fun))
                
                ttest_df$qval23<-p.adjust(ttest_df$p.val23,"fdr")
                
                ttest_df<-ttest_df%>%
                  dplyr::mutate(!!sig23:= if_else(qval23 < input$num77,if_else(abs(.data[[fc23]])>log2(input$num78), if_else(.data[[fc23]] > 0, "up","down"),"n.s."),"n.s."))
                names(ttest_df)[names(ttest_df) == "p.val23"] <- name23
                names(ttest_df)[names(ttest_df) == "qval23"]  <- qname23
                
                name24<-paste("p.val_",con8,"-",con3,sep="")
                qname24<-paste("q.val_",con8,"-",con3,sep="")
                fc24<-paste("log2FC_",con8,"-",con3,sep="")
                sig24<-paste("sig_",con8,"-",con3,sep="")
                val1 <-  con3
                val2 <-  con8
                ttest_df<-ttest_df%>%
                  dplyr::mutate(p.val24 = data%>%purrr::map_dbl(t.test_fun),
                                !!fc24:= data%>%purrr::map_dbl(fc_fun))
                
                ttest_df$qval24<-p.adjust(ttest_df$p.val24,"fdr")
                ttest_df<-ttest_df%>%
                  dplyr::mutate(!!sig24:= if_else(qval24 < input$num77,if_else(abs(.data[[fc24]])>log2(input$num78), if_else(.data[[fc24]] > 0, "up","down"),"n.s."),"n.s."))
                names(ttest_df)[names(ttest_df) == "p.val24"] <- name24
                names(ttest_df)[names(ttest_df) == "qval24"]  <- qname24
                
                name25<-paste("p.val_",con8,"-",con4,sep="")
                qname25<-paste("q.val_",con8,"-",con4,sep="")
                fc25<-paste("log2FC_",con8,"-",con4,sep="")
                sig25<-paste("sig_",con8,"-",con4,sep="")
                val1 <-  con4
                val2 <-  con8
                ttest_df<-ttest_df%>%
                  dplyr::mutate(p.val25 = data%>%purrr::map_dbl(t.test_fun),
                                !!fc25:= data%>%purrr::map_dbl(fc_fun))
                
                ttest_df$qval25<-p.adjust(ttest_df$p.val25,"fdr")
                ttest_df<-ttest_df%>%
                  dplyr::mutate(!!sig25:= if_else(qval25 < input$num77,if_else(abs(.data[[fc25]])>log2(input$num78), if_else(.data[[fc25]] > 0, "up","down"),"n.s."),"n.s."))
                names(ttest_df)[names(ttest_df) == "p.val25"] <- name25
                names(ttest_df)[names(ttest_df) == "qval25"]  <- qname25
                
                name26<-paste("p.val_",con8,"-",con5,sep="")
                qname26<-paste("q.val_",con8,"-",con5,sep="")
                fc26<-paste("log2FC_",con8,"-",con5,sep="")
                sig26<-paste("sig_",con8,"-",con5,sep="")
                val1 <-  con5
                val2 <-  con8
                ttest_df<-ttest_df%>%
                  dplyr::mutate(p.val26 = data%>%purrr::map_dbl(t.test_fun),
                                !!fc26:= data%>%purrr::map_dbl(fc_fun))
                
                ttest_df$qval26<-p.adjust(ttest_df$p.val26,"fdr")
                ttest_df<-ttest_df%>%
                  dplyr::mutate(!!sig26:= if_else(qval26 < input$num77,if_else(abs(.data[[fc26]])>log2(input$num78), if_else(.data[[fc26]] > 0, "up","down"),"n.s."),"n.s."))
                names(ttest_df)[names(ttest_df) == "p.val26"] <- name26
                names(ttest_df)[names(ttest_df) == "qval26"]  <- qname26
                
                
                name27<-paste("p.val_",con8,"-",con6,sep="")
                qname27<-paste("q.val_",con8,"-",con6,sep="")
                fc27<-paste("log2FC_",con8,"-",con6,sep="")
                sig27<-paste("sig_",con8,"-",con6,sep="")
                val1 <-  con6
                val2 <-  con8
                ttest_df<-ttest_df%>%
                  dplyr::mutate(p.val27 = data%>%purrr::map_dbl(t.test_fun),
                                !!fc27:= data%>%purrr::map_dbl(fc_fun))
                
                ttest_df$qval27<-p.adjust(ttest_df$p.val27,"fdr")
                ttest_df<-ttest_df%>%
                  dplyr::mutate(!!sig21:= if_else(qval27 < input$num77,if_else(abs(.data[[fc27]])>log2(input$num78), if_else(.data[[fc27]] > 0, "up","down"),"n.s."),"n.s."))
                names(ttest_df)[names(ttest_df) == "p.val27"] <- name27
                names(ttest_df)[names(ttest_df) == "qval27"]  <- qname27
                
                
                name28<-paste("p.val_",con8,"-",con7,sep="")
                qname28<-paste("q.val_",con8,"-",con7,sep="")
                fc28<-paste("log2FC_",con8,"-",con7,sep="")
                sig28<-paste("sig_",con8,"-",con7,sep="")
                val1 <-  con7
                val2 <-  con8
                ttest_df<-ttest_df%>%
                  dplyr::mutate(p.val28 = data%>%purrr::map_dbl(t.test_fun),
                                !!fc28:= data%>%purrr::map_dbl(fc_fun))
                
                ttest_df$qval28<-p.adjust(ttest_df$p.val28,"fdr")
                ttest_df<-ttest_df%>%
                  dplyr::mutate(!!sig27:= if_else(qval28 < input$num77,if_else(abs(.data[[fc28]])>log2(input$num78), if_else(.data[[fc28]] > 0, "up","down"),"n.s."),"n.s."))
                names(ttest_df)[names(ttest_df) == "p.val28"] <- name28
                names(ttest_df)[names(ttest_df) == "qval28"]  <- qname28
                
                
              }
              
              
              
            }
            
            
            
          }
          
          
          
        }
        }
        
      }
      
      
      
    }
    

    
    
    
    
    # ttest_df<-ttest_df_initial1%>%
    #   dplyr::mutate(p.val = data%>%purrr::map_dbl(t.test_fun),
    #          log2_foldchange = data%>%purrr::map_dbl(fc_fun),
    #          cond1 = data%>% purrr::map_dbl(cond1_fun))
    
    #ttest_df$q.val<- p.adjust(ttest_df$p.val, method="fdr")
    

    return(ttest_df)
  })
  
  
  # ttest_df_final<-reactive({
  #   ttest_df2<-ttest_df_initial()%>%
  #     dplyr::mutate(Significant = if_else(q.val < input$num3, if_else(abs(log2_foldchange)>log2(input$num4), if_else(log2_foldchange>0, "up","down"),"n.s."),"n.s."))
  #   return(ttest_df2)
  # })
  
  ttest_df_untidy_final<-reactive({
    ttest_df_untidy<- ttest_df_initial2()%>%
      tidyr::unnest(data)
    return(ttest_df_untidy)
  })
  
  ttest_slim<-reactive({
    # cond1_nam <- unique(norm_sum_final()$condition)[1]
    # cond2_nam <- unique(norm_sum_final()$condition)[2]
    
    ttest_df_slim<- ttest_df_initial2()%>%
      dplyr::select(-data)
    # %>%
    #   rename(!!cond1_nam:= cond1,
    #          !!cond2_nam:= cond2)
    return(ttest_df_slim)
  })
  
  heatmapinput <-  reactive({
    heatmap_input<-df_final_res()%>%
      dplyr::select(-intensity,-correct_factor,-relative_intensity,-norm_int)%>%
      dplyr::ungroup()%>%
      dplyr::group_by(ProtID)%>%
      dplyr::mutate(mean_row = mean(log2_int),
             sd_row = sd(log2_int),
             z_scale = (log2_int - mean_row)/sd_row,
             mean_scale = log2_int - mean_row)
    
    if (input$txtScale == "RowMean") {
      df_final2<-heatmap_input%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,cond_rep,mean_scale)%>%
        tidyr::spread(cond_rep,mean_scale)%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    } 
    if (input$txtScale == "Zscore") {
      df_final2<-heatmap_input%>%
        dplyr::ungroup()%>%
        dplyr::select(ProtID,cond_rep,z_scale)%>%
        tidyr::spread(cond_rep,z_scale)%>%
        tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    }
    
    
    if (input$anovaorall == "ANOVA") {
      sig<-ttest_slim()%>%dplyr::filter(anova_sig == "sig")
    }
    
    if (input$anovaorall == "All") {
      sig<-ttest_slim()
    }
    
    Refs<-sig$Reference
    
    #print(df_final2)
    
    df_final2<-df_final2%>%
      dplyr::filter(Reference %in% Refs)
    
    df_final2<-df_final2%>%
      tidyr::unite("ProtID",Reference,Gene.Symbol,Annotation, sep="__X__")
    
    # df3<-df_final2%>%
    #   dplyr::select(2:dim(df_final2)[[2]])
    df4<-df_final2%>%
      tibble::column_to_rownames("ProtID")
    
    # print(df4)
    
    return(df4)
    
    # map<-pheatmap(df4)
    # 
    # df4.clust <- cbind(df4, 
    #                    cluster = cutree(map$tree_row, 
    #                                     k = input$numclusters))
    # pheatmap(df4, cutree_rows = input$numclusters)

  })
  
  output$numclusterplot <- renderPlot({
    
    k.means.opt <- function(my.clusters) {
      k <- stats::kmeans(heatmapinput(), centers=my.clusters, iter.max = 500)
      if (is.finite(k$betweenss/k$totss)==T) {
        value<-k$betweenss/k$totss
      }
      return(value)
    }
    
    cluster.max <- 10
    ## make empty list to put the data
    k.means.opt.list <- vector("list", length=cluster.max) 
    clus_num<-c()
    max_var<-c()
    ## replicate function where, for each cluster from 1:cluster.max, replicates 1000 times and returns the (between.ss/tot.ss) ratio. Each of the [[elements]] of the list is 1000 ratios for each number of clusters. Also I want to time this step.
    for (i in 1:cluster.max) {
      k.means.opt.list[[i]] <- replicate(100, k.means.opt(i))
      clus_num <- c(clus_num,i)
      max_var<-c(max_var,max(k.means.opt.list[[i]]))
    }
    # print(clus_num)
    # print(max_var)
    ## now make a graph of all these replicates
    ## assign your x axis, number of clusters
    # clusters.x <- c(1:cluster.max)
    
    ## unlist your k.means.opt.list object to extract Ratio.BSS.TSS. Take the maximum Ratio.BSS.TSS of each 1000 kmeans replicates for each clutser, i (in this case, 1000 kmeans replicates for each cluster from 1 cluster (should be 0%) to 20 clusters (should be close to 100%)). If the number of clusters = number of data points, the "Ratio" will be undefined, which will show on the graph as 0% even though it's undefined.
    # k.means.data <- data.frame(Number.of.clusters= clusters.x, Ratio.BSS.TSS=unlist(lapply(k.means.opt.list, max)))
    
    k.means.data<-data.frame("cluster_num"=clus_num, "VarExp" = max_var)
    
    print(k.means.data)
    ## now plot these two things
    ggplot(k.means.data, aes(x = cluster_num, y = VarExp)) +
      geom_point(size=5,alpha=0.5) + geom_line()+
      theme_classic()+
      labs(y="BSS/TSS : Variance Explained B/w Groups", x="Cluster Number")+
      theme(axis.text = element_text(size = 18), axis.title = element_text(size=19))+ 
      scale_x_continuous(breaks=seq(0,10,2))
    
  
  })
  
  output$heatmapplot <- renderPlot({
    #cols <- colorRampPalette(brewer.pal(6,name="RdBu"))(12)
    # brks <- seq(-2,2,length.out=12)
    rg <- max(abs(heatmapinput()))
    
    #, cluster_cols = FALSE,color=cols,
    if(input$txtclustercols == "Yes"){
      p1<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercols == "No"){
      p1<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F, cluster_cols = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    p1
    
    
  })
  
  
  heatmapplotprint <- reactive({
    #cols <- colorRampPalette(brewer.pal(6,name="RdBu"))(12)
    # brks <- seq(-2,2,length.out=12)
    rg <- max(abs(heatmapinput()))
    
    #, cluster_cols = FALSE,color=cols,
    if(input$txtclustercols == "Yes"){
      p1<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F, breaks = seq(-rg, rg, length.out = 100))
    }
    
    if(input$txtclustercols == "No"){
      p1<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F, cluster_cols = FALSE, breaks = seq(-rg, rg, length.out = 100))
    }
    return(p1)
    
    
  })
  
  output$plotheatmapprint <- downloadHandler(
    filename = function() { paste("Heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = heatmapplotprint(), width = 10, height = 10)
    }
  )
  
  
  
  
  
  #,height=650

  clusterTraceDFinitial <- reactive({
    #, cluster_cols = FALSE
    
    if(input$txtclustercols == "Yes"){
      out<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F)
      heatmapinput2<-heatmapinput()%>%
        dplyr::select(colnames(heatmapinput()[,out$tree_col[["order"]]]))
      
      order_levels<-colnames(heatmapinput()[,out$tree_col[["order"]]])
    }
    
    if(input$txtclustercols == "No"){
      out<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "Heatmap",show_rownames=F, cluster_cols = FALSE)
      heatmapinput2<-heatmapinput()
    }
    #out<-pheatmap::pheatmap(heatmapinput(), cutree_rows = input$numclusters,main = "ANOVA Significant Proteins",show_rownames=F)
    
    
    
    # print(heatmapinput2)
    
    df4.clust <- cbind(heatmapinput2,
                       cluster = cutree(out$tree_row,
                                        k = input$numclusters))
    
    
    
    df4.clust2<-df4.clust%>%
      tibble::rownames_to_column("ProtID")%>%
      dplyr::group_by(ProtID,cluster)
    
    print(head(df4.clust2))
    
    return(df4.clust2)
    
  })
  
  hotellingdf<- reactive({
    if (input$txttimecourse == "Yes") {
      prot_values <- correlation_plot2()$ProtID
      
      df_timecourse <- correlation_plot2()%>%
        tibble::column_to_rownames(var = "ProtID")
      
      values <- c(names(df_timecourse))
      time.grp <- as.numeric(readr::parse_number(values))
      assay <- as.numeric(stringr::str_sub(values,-1,-1))
      if (input$numrepstimecourse == 3) {
        size <- rep(3, dim(df_timecourse)[[1]])
      }
      
      if (input$numrepstimecourse == 4) {
        size <- rep(4, dim(df_timecourse)[[1]])
      }
      
      if (input$numrepstimecourse == 2) {
        size <- rep(2, dim(df_timecourse)[[1]])
      }
      
      
      df_timecourse <- as.matrix(df_timecourse)
      
      hotelling <- timecourse::mb.long(df_timecourse, times=length(unique(time.grp)), reps=size, rep.grp=assay, time.grp=time.grp)
      hotelling_df<-data.frame("ProtID" = prot_values, "hotellingT2" = hotelling$HotellingT2)
    }
    
    if (input$txttimecourse == "No") {
      prot_values <- correlation_plot2()$ProtID
      hotelling_df<-data.frame("ProtID" = prot_values, "hotellingT2" = 0)
    }
    
    return(hotelling_df)
  })
  
  clustertimecourse<- reactive({
    
    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)
    
    HotelClusterDF <- dplyr::left_join(hotellingdf(),clusterdesignation,by="ProtID")%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
    
    HotelClusterDF <- HotelClusterDF%>%
      #dplyr::mutate(cluster = dplyr::if_else(is.na(cluster), "0", cluster))%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    all_table<-left_join(tableresults(),HotelClusterDF, by= c("Reference","Gene.Symbol","Annotation","cluster"))%>%
      dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    

    return(all_table%>%dplyr::arrange(desc(hotellingT2)))
  })
  
  
 
  
  
  output$hotellingdf2 <- renderReactable({
    #write.csv(clustertimecourse(),"C:/Harper/Side_analysis/Bobby_shiny/BH_timecourse_5cluster_HotellingT2.csv")
    reactable(clustertimecourse(),filterable = T)
  })
  
  
  output$hotellingdownload <- downloadHandler(
    filename = function() {
      paste("hotelling_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(clustertimecourse(), file, row.names = FALSE)
    }
  )

  #clustertimecourse

  output$plothotelling <- renderPlotly({

    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2])
    
    foldchangevalue2<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    photelling<-ggplot(clustertimecourse(), aes( .data[[foldchangevalue2]], log10(hotellingT2), color=factor(cluster),label=Gene.Symbol, key=ProtID))+
      geom_point(alpha=0.5,size=0.5)+theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x= paste(foldchangevalue2), y= "Log10(Hotelling T^2)")
    
    ggplotly(photelling)%>% layout(dragmode = "select")


  })
  
  
  plothotellingprint <- reactive({
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2])
    
    foldchangevalue2<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    photelling<-ggplot(clustertimecourse(), aes( .data[[foldchangevalue2]], log10(hotellingT2), color=factor(cluster)))+
      geom_point(alpha=0.5,size=0.5)+theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x= paste(foldchangevalue2), y= "Log10(Hotelling T^2)",color="Cluster")
    
    return(photelling)
    
    
  })
  
  
  output$printhotellingplot<- downloadHandler(
    filename = function() { paste("HotellingPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plothotellingprint(), width = 9, height = 5)
    }
  )
  
  
  
  
  output$hotellingclick <- renderReactable({
    d <- event_data("plotly_selected")
    req(d)
    df_sub3<-clustertimecourse()%>%
      # dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))%>%
      dplyr::filter(ProtID %in% d$key )%>%
      #dplyr::rename(!!unique(norm_sum_final()$condition)[1]:=cond1,!!unique(norm_sum_final()$condition)[2]:=cond2)%>%
      # dplyr::select(-data)%>%
      dplyr::ungroup()

    reactable(df_sub3, filterable = TRUE)
  })
  
  
  output$plottcprot <- renderPlot({
    req(input$text113)
    gene_name_input <- input$text113
    
    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)

    HotelClusterDF <- dplyr::left_join(hotellingdf(),clusterdesignation,by="ProtID")%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)

    HotelClusterDF <- HotelClusterDF%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    plot_protein<-left_join(ttest_df_untidy_final(),HotelClusterDF, by= c("Reference","Gene.Symbol","Annotation"))%>%
      dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
  print(plot_protein[1:3,])
  print(plot_protein$condition[1:8])
  
  # print(readr::parse_number(plot_protein$condition[1]))
    plot_protein<-plot_protein%>%
      dplyr::filter(Gene.Symbol==gene_name_input)%>%
      dplyr::mutate(condition = as.character(condition),
                    timepoint = readr::parse_number(condition))%>%
      dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))
    
    # plot_protein$timepoint <- gsub("[^0-9]","",plot_protein$condition)


    plot_protein_v2<-plot_protein%>%
      dplyr::ungroup()%>%
      dplyr::group_by(condition,timepoint,ProtID2,Reference)%>%
      dplyr::summarise(log2_int=median(log2_int))%>%
      dplyr::mutate(timepoint = as.numeric(as.character(timepoint)))
    
    plot_protein<-plot_protein%>%
      dplyr::mutate(timepoint=as.numeric(as.character(timepoint)))
    
    #,color=as.factor(timepoint)
    ggplot()+
      geom_line(data=plot_protein_v2,aes(x=timepoint,y=2^log2_int,color=Reference),alpha=0.5)+
      geom_point(data = plot_protein ,aes(x=timepoint,y=2^log2_int,color=Reference),size=3)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5),legend.position = "top")+
      labs(title = gene_name_input, x="Time",y="MS Intensity")
    
  })
  
  plottcprotprint <- reactive({
    req(input$text113)
    gene_name_input <- input$text113
    
    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)
    
    HotelClusterDF <- dplyr::left_join(hotellingdf(),clusterdesignation,by="ProtID")%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
    
    HotelClusterDF <- HotelClusterDF%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    plot_protein<-left_join(ttest_df_untidy_final(),HotelClusterDF, by= c("Reference","Gene.Symbol","Annotation"))%>%
      dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    
    # print(readr::parse_number(plot_protein$condition[1]))
    plot_protein<-plot_protein%>%
      dplyr::filter(Gene.Symbol==gene_name_input)%>%
      dplyr::mutate(condition = as.character(condition),
                    timepoint = readr::parse_number(condition))%>%
      dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))
    
    # plot_protein$timepoint <- gsub("[^0-9]","",plot_protein$condition)
    
    
    plot_protein_v2<-plot_protein%>%
      dplyr::ungroup()%>%
      dplyr::group_by(condition,timepoint,ProtID2,Reference)%>%
      dplyr::summarise(log2_int=median(log2_int))%>%
      dplyr::mutate(timepoint = as.numeric(as.character(timepoint)))
    
    plot_protein<-plot_protein%>%
      dplyr::mutate(timepoint=as.numeric(as.character(timepoint)))
    
    #,color=as.factor(timepoint)
    p1<-ggplot()+
      geom_line(data=plot_protein_v2,aes(x=timepoint,y=2^log2_int,color=Reference),alpha=0.5)+
      geom_point(data = plot_protein ,aes(x=timepoint,y=2^log2_int,color=Reference),size=3)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5),legend.position = "top")+
      labs(title = gene_name_input, x="Time",y="MS Intensity")
    return(p1)
    
  })
  

  
  output$plottcprotprint2 <- downloadHandler(
    filename = function() { paste("Timecourse_", input$text113 ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plottcprotprint(), width = 10, height = 10)
    }
  )
  
  
  
  output$volcanoplot <- renderPlotly({

    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    if (input$txtinversion == "No") {
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)")
    }
    
    if (input$txtinversion == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)")
    }
    
    ggplotly(volcano)%>% layout(dragmode = "select")
    
    
  })
  
  
  volcanoplotprint <- reactive({
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversion == "No") {
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)")
    }
    
    if (input$txtinversion == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      volcano<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)")
    }
    
    return(volcano)
    
    
  })
  
  
  output$printvolcano<- downloadHandler(
    filename = function() { paste("Volcano",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol]) ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = volcanoplotprint(), width = 9, height = 5)
    }
  )
  
  
  
  
  output$volcanoclick <- renderReactable({
    d <- event_data("plotly_selected")
    req(d)
    df_sub2<-tableresults()%>%
      dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))%>%
      dplyr::filter(ProtID %in% d$key )%>%
      #dplyr::rename(!!unique(norm_sum_final()$condition)[1]:=cond1,!!unique(norm_sum_final()$condition)[2]:=cond2)%>%
      # dplyr::select(-data)%>%
      dplyr::ungroup()
    
    reactable(df_sub2, filterable = TRUE)
  })

  output$volcanoplotnorm <- renderPlotly({
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse3vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse4vol])
    con3<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse5vol])
    con4<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse6vol])
    
    # foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    # qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    
    pair_plot<-tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    
    volcano2<-ggplot(pair_plot, aes( .data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]],label=Gene.Symbol, key=ProtID))+
      geom_point(alpha=0.3,size=0.5)+theme_classic()+
      # scale_color_viridis_d(end=0.8)+
      geom_abline(slope=1,linetype="dashed",size=0.5)+
      geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
      geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
      labs(x= paste(con2,"-",con1), y= paste(con4,"-",con3))
    
    ggplotly(volcano2)%>% layout(dragmode = "select")
    
  })
  
  volcanoplotnormprint <- reactive({
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse3vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse4vol])
    con3<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse5vol])
    con4<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse6vol])
    

    pair_plot<-tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_"))
    
    volcano2<-ggplot(pair_plot, aes( .data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]))+
      geom_point(alpha=0.3,size=0.5)+theme_classic()+
      geom_abline(slope=1,linetype="dashed",size=0.5)+
      geom_vline(xintercept = 0,linetype="dashed",size=0.5)+
      geom_hline(yintercept = 0,linetype="dashed",size=0.5)+
      labs(x= paste(con2,"-",con1), y= paste(con4,"-",con3))
    
    return(volcano2)
    
  })
  
  output$printcorrtwocond<- downloadHandler(
    filename = function() { paste("CorrPlot_x_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse4vol]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse3vol]),"_y_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse6vol]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse5vol]) ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = volcanoplotnormprint(), width = 9, height = 5)
    }
  )

  clusterTraceDF<-reactive({
    #, cluster_cols = FALSE
    
    if(input$txtclustercols == "Yes"){
      out<-pheatmap::pheatmap(heatmapinput(),   cutree_rows = input$numclusters,main = "ANOVA Significant Proteins",show_rownames=F)
      heatmapinput2<-heatmapinput()%>%
        dplyr::select(colnames(heatmapinput()[,out$tree_col[["order"]]]))
      order_levels<-colnames(heatmapinput()[,out$tree_col[["order"]]])
      
      length_dim<-dim(clusterTraceDFinitial())[[2]]-1
      
      df4.clust3<-clusterTraceDFinitial()%>%
        tidyr::gather("cond_rep","value",2:all_of(length_dim))
      
      df4.clust3$cond_rep <- factor(df4.clust3$cond_rep,levels = c(order_levels))
    }
    
    if(input$txtclustercols == "No"){
      out<-pheatmap::pheatmap(heatmapinput(),  cutree_rows = input$numclusters,main = "ANOVA Significant Proteins",show_rownames=F, cluster_cols = FALSE)
      heatmapinput2<-heatmapinput()
      
      length_dim<-dim(clusterTraceDFinitial())[[2]]-1
      
      df4.clust3<-clusterTraceDFinitial()%>%
        tidyr::gather("cond_rep","value",2:all_of(length_dim))
      
      df4.clust3$cond_rep <- factor(df4.clust3$cond_rep,levels = names(clusterTraceDFinitial())[-1])
  
    }
    
    
    return(df4.clust3)
    
  })

  output$plotclustertrace <- renderPlot({
    
    df_cluster<-clusterTraceDF()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")%>%
      dplyr::mutate(condition = stringr::str_sub(cond_rep,start=1,end=-3))
    
    df_cluster_med <- df_cluster%>%
      dplyr::group_by(cond_rep,cluster,condition)%>%
      dplyr::summarise(value = median(value))
      
    
    ggplot()+geom_violin(data=df_cluster,aes(cond_rep,value,color=condition),draw_quantiles = c(0.25,0.5,0.75))+geom_point(data=df_cluster_med,aes(cond_rep,value,color=condition))+
      theme_classic()+labs(x="",y="Scaled Intensity")+facet_wrap(vars(cluster),nrow=2)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept = 0,linetype="dashed")
  })
  
  
  plotclustertraceprint <- reactive({
    
    df_cluster<-clusterTraceDF()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")%>%
      dplyr::mutate(condition = stringr::str_sub(cond_rep,start=1,end=-3))
    
    df_cluster_med <- df_cluster%>%
      dplyr::group_by(cond_rep,cluster,condition)%>%
      dplyr::summarise(value = median(value))
    
    
    return(ggplot()+geom_violin(data=df_cluster,aes(cond_rep,value,color=condition),draw_quantiles = c(0.25,0.5,0.75))+geom_point(data=df_cluster_med,aes(cond_rep,value,color=condition))+
      theme_classic()+labs(x="",y="Scaled Intensity")+facet_wrap(vars(cluster),nrow=2)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept = 0,linetype="dashed"))
  })
  
  output$printclustersplot <- downloadHandler(
    filename = function() { paste("ClusterViolins", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotclustertraceprint(), width = 16, height = 10)
    }
  )
  
  

  output$numbercluster<-renderReactable({
    df_clustercount<-clusterTraceDF()%>%
      dplyr::select(ProtID,cluster)%>%
      dplyr::distinct()%>%
      dplyr::group_by(cluster)%>%
      dplyr::summarise(Num_Proteins_perCluster = n())
    reactable(df_clustercount)
  })

# output$plotvolcanoes<-renderPlot({
#   con1<-as.character(unique(norm_sum_final()$condition)[1])
#   con2<-as.character(unique(norm_sum_final()$condition)[2])
#   con3<-as.character(unique(norm_sum_final()$condition)[3])
#   
# })

  
  
  output$bioplexNetwork <-renderPlot({
    
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    binary_slim <- bioplex_binary%>%
      dplyr::select(Gene.Symbol_target,Gene.Symbol_interactor)%>%
      dplyr::rename(SymbolA = Gene.Symbol_target ,
                    SymbolB = Gene.Symbol_interactor)
    
    #style = list(fontFamily = 'Menlo',fontSize = '14px')
    
    # int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )%>%
    #   dplyr::bind_rows(.,to_bind) 
    # network<- igraph::graph_from_data_frame(int_df_ex,directed=F)  
    
    if (input$txthumanmouse == "mouse") {
      results<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::rename(name = Gene.Symbol)
    } 
    
    if (input$txthumanmouse == "human") {
      results<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::rename(name = Gene.Symbol)
    }
    
    
    
    # val1<-colnames(results)[8]
    # val2<-colnames(results)[9]
    results2<-results%>%
      # dplyr::rename(name = Gene.Symbol)%>%
      #,Log2_FC = log2_foldchange
      dplyr::mutate(med_two = (.data[[con1]] + .data[[con2]])/2)%>%
      dplyr::group_by(name)%>%
      dplyr::mutate(max_val = max(med_two))%>%
      dplyr::filter(med_two == max_val)%>%
      dplyr::ungroup()
    
    ref_name_df<-results2%>%
      dplyr::filter(name == input$txtbioplex)
    
    to_bind<-data.frame("SymbolA" = input$txtbioplex,"SymbolB" = input$txtbioplex)
    
    binary_select<-binary_slim%>%
      dplyr::filter(SymbolA == input$txtbioplex)%>%
      dplyr::bind_rows(.,to_bind)
    
    results3<-results2%>%
      filter(name %in% unique(binary_select$SymbolB))
    
    binary_select2<-binary_select%>%dplyr::filter(SymbolB %in% unique(results3$name))
    
    network<- igraph::graph_from_data_frame(binary_select2,directed=F)
    
    N<-tidygraph::as_tbl_graph(network)
    
    ggraph::set_graph_style()
    
    
    N%>%
      tidygraph::activate(nodes)%>%
      tidygraph::left_join(.,results3,by="name")%>%
      # tidygraph::activate(edges)%>%
      # tidygraph::inner_join(.,results2,by="name")%>%
      ggraph::ggraph(layout = "stress")+
      ggraph::geom_edge_fan(width=0.5,alpha=0.5)+
      ggraph::geom_node_point(aes(color=.data[[fc1]]),size=25,alpha=0.5)+
      ggraph::geom_node_text(aes(label=name),size=3)+ggraph::scale_color_viridis(end=0.8)
    #, repel = TRUE

    
  },width=800,height=800)
  
  


  
  
  
  
  output$bioplexExp<-renderPlot({
    req(input$txtbioplex)
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors")
    
    
    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }
    
    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }

    ggplot()+geom_boxplot(data=target_interactors,aes(plot_value,.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=target_interactors,aes(plot_value,.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Protein",y=paste(fc1))+
      theme(axis.text = element_text(size=18),axis.title = element_text(size=20),legend.title=element_text(size=17), 
            legend.text=element_text(size=15))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  bioplexExpprint <- reactive({
    req(input$txtbioplex)
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors")
    
    
    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }
    
    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }
    
    return(ggplot()+geom_boxplot(data=target_interactors,aes(plot_value,.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=target_interactors,aes(plot_value,.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Protein",y=paste(fc1))+
      theme(axis.text = element_text(size=18),axis.title = element_text(size=20),legend.title=element_text(size=17), 
            legend.text=element_text(size=15))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  
  output$printbioplexExp <- downloadHandler(
    filename = function() { paste("Bioplex_",input$txtbioplex,"_",as.character(unique(norm_sum_final()$condition)[input$numconbioplex2]),"-",as.character(unique(norm_sum_final()$condition)[input$numconbioplex1]) ,'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplexExpprint(), width = 12, height = 8)
    }
  )
  
  
  output$bioplexExpAll<-renderPlot({

    
    int_df_ex<-int_df%>%
      dplyr::select(SymbolA,SymbolB)%>%
      dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors",sep="")
    
    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }
    
    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }
    

    
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,target)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    print(unique(dftidy$condition))
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }
    
    

    #plotname=unique(dftidy$target)
    
    maxvalue<-abs(max(dftidy$log2_FC))+0.2
    
    ggplot(dftidy, aes(condition,log2_FC,color=condition,fill=condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Bioplex Foldchange over Conditions",title=paste(plot_value))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue))
    

  })
  
  
  bioplexExpAllprint <- reactive({
    
    
    int_df_ex<-int_df%>%
      dplyr::select(SymbolA,SymbolB)%>%
      dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors",sep="")
    
    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }
    
    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)%>%
        dplyr::mutate(target = paste(input$txtbioplex))
    }
    
    
    
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- target_interactors%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,target,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,target)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }
    
    
    
    #plotname=unique(dftidy$target)
    
    maxvalue<-abs(max(dftidy$log2_FC))+0.2
    
    return(ggplot(dftidy, aes(condition,log2_FC,color=condition,fill=condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Bioplex Foldchange over Conditions",title=paste(plot_value))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue)))
    
    
  })
  
  output$printbioplexExpAll <- downloadHandler(
    filename = function() { paste("Bioplex_",input$txtbioplex,"_multicomparison",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplexExpAllprint(), width = 16, height = 8)
    }
  )
  
  
  output$tablebioplex<-renderReactable({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    int_df_ex<-int_df%>%dplyr::select(SymbolA,SymbolB)%>%dplyr::filter(SymbolA == input$txtbioplex | SymbolB == input$txtbioplex )
    bioplex_interactors_target<- unique(append(int_df_ex$SymbolA,int_df_ex$SymbolB))
    plot_value <-paste(input$txtbioplex, " interactors")
    
    if (input$txthumanmouse == "human") {
      target_interactors<- tableresults()%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }
    if (input$txthumanmouse == "mouse") {
      target_interactors<- tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::filter(Gene.Symbol %in% bioplex_interactors_target)
    }
    
    

    reactable(target_interactors, filterable = TRUE)
  })
  
###based on gene.symbol now because of losing info with reference
  
  binary_bioplex_join<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    
    
    if (input$txthumanmouse == "human") {
      simplify_df<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }
    if (input$txthumanmouse == "mouse") {
      simplify_df<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }
    
    
    # %>%
    #   dplyr::rename(log2FC = log2_foldchange)
    simplify_df_target<- simplify_df%>%
      dplyr::rename(Gene.Symbol_target= Gene.Symbol)
    simplify_df_interactor<- simplify_df%>%
      dplyr::rename(Gene.Symbol_interactor= Gene.Symbol)
    bioplex_ref<-left_join(simplify_df_target,bioplex_binary,by="Gene.Symbol_target")
    bioplex_combine<- left_join(bioplex_ref,simplify_df_interactor,by="Gene.Symbol_interactor",suffix=c("_target","_interactor"))
    bioplex_combine2<-bioplex_combine%>%
      dplyr::filter(!is.na(.data[[qval1target]]),!is.na(.data[[qval1interactor]]))%>%
      dplyr::group_by(Gene.Symbol_target)%>%
      dplyr::mutate(interactomeFC=median(.data[[fc1interactor]]),
                    interactome_count=n())%>%
      dplyr::arrange(dplyr::desc(interactomeFC))
  })
  
  
  output$alldownloadbioplex <- downloadHandler(
    filename = function() {
      paste("all_biolpex_results_proteinLevel",as.character(unique(norm_sum_final()$condition)[input$numconbioplex2]),"-",as.character(unique(norm_sum_final()$condition)[input$numconbioplex1]), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(binary_bioplex_join()%>%dplyr::filter(interactome_count >= input$num10001), file, row.names = FALSE)
    }
  )
  
  #binary_bioplex_join
  
  output$tableBinaryBioplex<-renderReactable({
    req(input$num10001)
    bioplex_output<-binary_bioplex_join2()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10001)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)
    
    
    reactable(bioplex_output, filterable = TRUE)
  })
  
  output$tableBinaryBioplex2<-renderReactable({
    req(input$num10002)
    bioplex_output<-binary_bioplex_join()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10002)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)
    
    
    reactable(bioplex_output, filterable = TRUE)
  })
  
  bioplex_binary_final<-reactive({
    req(input$num10001)
    bioplex_output<-binary_bioplex_join()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10001)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)
    
    all_bioplex_info<-dplyr::inner_join(bioplex_output,binary_bioplex_join(),by=c("Gene.Symbol_target","interactomeFC","interactome_count"))
  })
  
  binary_bioplex_join2<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    
    
    if (input$txthumanmouse == "human") {
      simplify_df<-tableresults()%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }
    if (input$txthumanmouse == "mouse") {
      simplify_df<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
        dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    }
    
    
    # simplify_df<-tableresults()%>%
    #   tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
    #   dplyr::select(Gene.Symbol, .data[[qval1]],.data[[fc1]],.data[[sig1]])
    # %>%
    #   dplyr::rename(log2FC = log2_foldchange)
    simplify_df_target<- simplify_df%>%
      dplyr::rename(Gene.Symbol_target= Gene.Symbol)
    simplify_df_interactor<- simplify_df%>%
      dplyr::rename(Gene.Symbol_interactor= Gene.Symbol)
    bioplex_ref<-left_join(simplify_df_target,bioplex_binary,by="Gene.Symbol_target")
    bioplex_combine<- left_join(bioplex_ref,simplify_df_interactor,by="Gene.Symbol_interactor",suffix=c("_target","_interactor"))
    bioplex_combine2<-bioplex_combine%>%
      dplyr::filter(!is.na(.data[[qval1target]]),!is.na(.data[[qval1interactor]]))%>%
      dplyr::group_by(Gene.Symbol_target)%>%
      dplyr::mutate(interactomeFC=median(.data[[fc1interactor]]),
                    interactome_count=n())%>%
      dplyr::arrange(dplyr::desc(interactomeFC))
  })
  
  bioplex_binary_final2<-reactive({
    req(input$num10001)
    bioplex_output<-binary_bioplex_join2()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol_target,interactomeFC,interactome_count)%>%
      dplyr::distinct()%>%
      dplyr::filter(interactome_count >= input$num10001)%>%
      tibble::rownames_to_column("BioplexRank")%>%
      dplyr::mutate(BioplexRank=as.numeric(as.character(BioplexRank)))%>%
      dplyr::arrange(BioplexRank)
    
    all_bioplex_info<-dplyr::inner_join(bioplex_output,binary_bioplex_join2(),by=c("Gene.Symbol_target","interactomeFC","interactome_count"))
  })
  
    
  output$bioplextopN<-renderPlot({
    req(input$numbioplex1)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")
    
    select_bioplexdf<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank <= input$numbioplex1)
    
    print(colnames(select_bioplexdf))
    
    ggplot()+geom_boxplot(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
    })
  
  bioplextopNprint <- reactive({
    req(input$numbioplex1)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")
    
    select_bioplexdf<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank <= input$numbioplex1)
    

    
    return(ggplot()+geom_boxplot(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
    
  })
  
  

  output$bioplexbottomN<-renderPlot({
    req(input$numbioplex2)
    num_limit<-input$numbioplex2
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")
    
    
    shrink_group<-bioplex_binary_final2()%>%
      dplyr::arrange(desc(BioplexRank))%>%
      dplyr::select(BioplexRank)%>%
      dplyr::distinct()
    
    list_ranks<-c(shrink_group$BioplexRank[1:num_limit])
    
    select_bioplexdf2<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank %in% list_ranks)

    ggplot()+geom_boxplot(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
    
    
  })
  
  bioplexbottomNprint <- reactive({
    req(input$numbioplex2)
    num_limit<-input$numbioplex2
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconbioplex4])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qval1<-as.character(paste("q.val_",con2,"-",con1,sep=""))
    qval1target<-paste(qval1,"_target",sep="")
    qval1interactor <- paste(qval1,"_interactor",sep="")
    fc1interactor <- paste(fc1,"_interactor",sep="")
    sig1interactor <- paste(sig1,"_interactor",sep="")
    
    
    shrink_group<-bioplex_binary_final2()%>%
      dplyr::arrange(desc(BioplexRank))%>%
      dplyr::select(BioplexRank)%>%
      dplyr::distinct()
    
    list_ranks<-c(shrink_group$BioplexRank[1:num_limit])
    
    select_bioplexdf2<-bioplex_binary_final2()%>%
      dplyr::filter(BioplexRank %in% list_ranks)
    
    return(ggplot()+geom_boxplot(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]]),outlier.shape=NA)+
      geom_jitter(data=select_bioplexdf2,aes(reorder(Gene.Symbol_target, .data[[fc1interactor]], FUN=median),.data[[fc1interactor]],color=.data[[sig1interactor]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Protein Interactome",y=bquote(""~Log[2]~"(Fold Change Interactors)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
    
    
  })
  

  
  output$printbiolplextop <- downloadHandler(
    filename = function() { paste("TopN_UpBioplex",as.character(unique(norm_sum_final()$condition)[input$numconbioplex4]),"-",as.character(unique(norm_sum_final()$condition)[input$numconbioplex3]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplextopNprint(), width = 12, height = 8)
    }
  )
  
  output$printbiolplexbottom <- downloadHandler(
    filename = function() { paste("TopN_DownBioplex",as.character(unique(norm_sum_final()$condition)[input$numconbioplex4]),"-",as.character(unique(norm_sum_final()$condition)[input$numconbioplex3]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = bioplexbottomNprint(), width = 12, height = 8)
    }
  )
  
  corum_df<-reactive({
    
      con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
      con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
      fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
      medname <- paste("median_log2FC_",con2,"-",con1,sep="")

    corum_binder<-tableresults()%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
    corum_binder<-dplyr::left_join(corum(),corum_binder,by="Reference")
    corum_binder<-corum_binder%>%
      dplyr::group_by(ComplexName)%>%
      dplyr::mutate(total_complexMembers = n())
    corum_binder2<-corum_binder%>%
      dplyr::select(ComplexName,.data[[fc1]],total_complexMembers)%>%
      dplyr::filter(!is.na(.data[[fc1]]))%>%
      dplyr::group_by(ComplexName,total_complexMembers)%>%
      dplyr::summarise(complexMembersIDed=n(),
                !!medname:= median(.data[[fc1]]))%>%
      dplyr::filter(complexMembersIDed>input$numcomplexmin-1)%>%
      dplyr::arrange(desc(.data[[medname]]))%>%
      tibble::rownames_to_column("ComplexRank")%>%
      dplyr::mutate(ComplexRank=as.numeric(as.character(ComplexRank)))%>%
      dplyr::arrange(ComplexRank)
    # complete_corum_df<-dplyr::left_join(corum_binder,corum_binder2,by=c("ComplexName","total_complexMembers"))%>%
    #   dplyr::filter(complexMembersIDed>1)%>%
    #   dplyr::arrange(desc(median_log2FC))
    return(corum_binder2)
  })
  
  complete_corum_df2<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")
    
    corum_binder<-ttest_slim()%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
    corum_binder<-dplyr::left_join(corum(),corum_binder,by="Reference")
    corum_binder<-corum_binder%>%
      dplyr::group_by(ComplexName)%>%
      dplyr::mutate(total_complexMembers = n())
    corum_binder2<-corum_binder%>%
      dplyr::distinct()%>%
      dplyr::select(ComplexName,.data[[fc1]],total_complexMembers)%>%
      dplyr::filter(!is.na(.data[[fc1]]))%>%
      dplyr::group_by(ComplexName,total_complexMembers)%>%
      dplyr::summarise(complexMembersIDed=n(),
                       !!medname:= median(.data[[fc1]]))%>%
      dplyr::filter(complexMembersIDed>input$numcomplexmin-1)%>%
      dplyr::arrange(desc(.data[[medname]]))%>%
      tibble::rownames_to_column("ComplexRank")%>%
      dplyr::mutate(ComplexRank=as.numeric(as.character(ComplexRank)))%>%
      dplyr::arrange(ComplexRank)
    complete_corum_df<-dplyr::left_join(corum_binder,corum_binder2,by=c("ComplexName","total_complexMembers"))%>%
      dplyr::ungroup()%>%
      dplyr::filter(!is.na(Gene.Symbol),!is.na(complexMembersIDed))
    return(complete_corum_df)
  })
  
  
  output$corumNetwork <-renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")
    
    df_title<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank == input$complexranknum  )
    plot_title<-unique(df_title$ComplexName)
    df_for_network<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank == input$complexranknum  )%>%
      dplyr::mutate(complex="complex")%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol,complex)
    network<- igraph::graph_from_data_frame(df_for_network,directed=F)  
    
    N<-tidygraph::as_tbl_graph(network)
    
    ggraph::set_graph_style()
    
    df_for_network2<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank == input$complexranknum  )%>%
      dplyr::mutate(complex="complex")%>%
      dplyr::ungroup()%>%
      dplyr::rename(name = Gene.Symbol,
                    !!fc1:= .data[[fc1]])
    
    
    N%>%
      tidygraph::activate(nodes)%>%
      tidygraph::left_join(.,df_for_network2,by="name")%>%
      # tidygraph::activate(edges)%>%
      # tidygraph::inner_join(.,results2,by="name")%>%
      ggraph::ggraph(layout = "stress")+
      ggraph::geom_edge_fan(width=0.5,alpha=0.5)+
      ggraph::geom_node_point(aes(color=.data[[fc1]]),size=25,alpha=0.5)+
      ggraph::geom_node_text(aes(label=name),size=3)+ggraph::scale_color_viridis(end=0.8)+
      ggplot2::ggtitle(paste(plot_title))
    # plot(network,vertex.size=40,vertex.color="cadetblue1",edge.curved=0)
    # title(paste(plot_title),cex.main=1.2,col.main="black")
    # title()
  })
  
  output$corumExp<-renderPlot({
    req(input$numcorum)
    choosen_prot<-input$numcorum
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank<choosen_prot+1)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    # median_delta_value <- unique(select_complex$median_log2FC)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)

    ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  corumExpprint<-reactive({
    req(input$numcorum)
    choosen_prot<-input$numcorum
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank<choosen_prot+1)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    # median_delta_value <- unique(select_complex$median_log2FC)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)
    
    return(ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  
  
  output$corumExp2<-renderPlot({
    req(input$numcorum2)
    num_limit<-input$numcorum2
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    shrink_group<-complete_corum_df2()%>%
      dplyr::arrange(desc(ComplexRank))%>%
      dplyr::select(ComplexRank)%>%
      dplyr::distinct()
    
    list_ranks<-c(shrink_group$ComplexRank[1:num_limit])
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank %in% list_ranks)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    # median_delta_value <- unique(select_complex$median_log2FC)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)
    
    ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  corumExp2print<-reactive({
    req(input$numcorum2)
    num_limit<-input$numcorum2
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    shrink_group<-complete_corum_df2()%>%
      dplyr::arrange(desc(ComplexRank))%>%
      dplyr::select(ComplexRank)%>%
      dplyr::distinct()
    
    list_ranks<-c(shrink_group$ComplexRank[1:num_limit])
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank %in% list_ranks)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    # median_delta_value <- unique(select_complex$median_log2FC)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)
    
    return(ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),alpha=0.75,size=3)+
      theme_classic()+labs(x="Complex",y=bquote(""~Log[2]~"(Fold Change)"))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  
  output$printcorumExp <- downloadHandler(
    filename = function() { paste("TopN_UpCorum",as.character(unique(norm_sum_final()$condition)[input$numcon2]),"-",as.character(unique(norm_sum_final()$condition)[input$numcon1]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExpprint(), width = 12, height = 8)
    }
  )
  
  output$printcorumExp2 <- downloadHandler(
    filename = function() { paste("TopN_DownCorum",as.character(unique(norm_sum_final()$condition)[input$numcon2]),"-",as.character(unique(norm_sum_final()$condition)[input$numcon1]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExp2print(), width = 12, height = 8)
    }
  )
  
  
  
  output$corumExp3<-renderPlot({
    choosen_prot<-input$numcorum
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    # median_delta_value <- round(unique(select_complex$median_log2FC),3)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)
    
    ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Complex",y=paste(fc1))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip()
  })
  
  corumExp3print <- reactive({
    choosen_prot<-input$numcorum
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    # median_delta_value <- round(unique(select_complex$median_log2FC),3)
    total_comp<-unique(select_complex$total_complexMembers)
    id_comp <-unique(select_complex$complexMembersIDed)
    
    return(ggplot()+geom_boxplot(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]]),outlier.shape=NA)+
      geom_jitter(data=select_complex,aes(reorder(ComplexName_slim, .data[[fc1]], FUN=median),.data[[fc1]],color=.data[[sig1]]),size=5)+
      theme_classic()+labs(x="Complex",y=paste(fc1))+
      theme(axis.text = element_text(size=15),axis.title = element_text(size=18),legend.title=element_text(size=15),
            legend.text=element_text(size=13))+
      geom_hline(yintercept = 0, linetype="dashed")+
      scale_color_viridis_d(end=0.8)+coord_flip())
  })
  
  nameComplex <- reactive({
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    return(unique(select_complex$ComplexName_slim))
  })
  
  
  output$printcorumExp3 <- downloadHandler(
    filename = function() { paste("Corum_SingleCompare_",nameComplex(),"_",as.character(unique(norm_sum_final()$condition)[input$numcon2]),"-",as.character(unique(norm_sum_final()$condition)[input$numcon1]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExp3print(), width = 12, height = 8)
    }
  )
  

  
  output$corumExpallcond<-renderPlot({
    choosen_prot<-input$numcorum
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,ComplexName_slim)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    print(unique(dftidy$condition))
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }
    
    
    

    plotname=unique(dftidy$ComplexName_slim)
    
    maxvalue<-abs(max(dftidy$log2_FC))+0.2
    
    ggplot(dftidy, aes(condition,log2_FC,color=condition,fill=condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Complex Foldchange over Conditions",title=paste(plotname))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue))
    

  })
  
  output$alldownloadCorum <- downloadHandler(
    filename = function() {
      paste("all_corum_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(complete_corum_df2(), file, row.names = FALSE)
    }
  )
  
  
  corumExpallcondprint <- reactive({
    choosen_prot<-input$numcorum
    
    select_complex<-complete_corum_df2()%>%
      dplyr::filter(ComplexRank==input$complexranknum)%>%
      dplyr::mutate(ComplexName_slim = stringr::str_sub(ComplexName,1,50))
    
    
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- select_complex%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,ComplexName_slim,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,ComplexName_slim)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    print(unique(dftidy$condition))
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }


    plotname=unique(dftidy$ComplexName_slim)
    
    maxvalue<-abs(max(dftidy$log2_FC))+0.2
    
    return(ggplot(dftidy, aes(condition,log2_FC,color=condition,fill=condition))+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      geom_jitter(alpha=0.5,size=5)+theme_classic()+
      geom_hline(yintercept = 0,linetype="dashed")+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),subtitle="Complex Foldchange over Conditions",title=paste(plotname))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=20,hjust=0.5))+
      ylim(c(-maxvalue,maxvalue)))
    

  })
  
  output$printcorumExpallcond <- downloadHandler(
    filename = function() { paste("Corum_MultiComparison_",nameComplex(), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corumExpallcondprint(), width = 16, height = 8)
    }
  )
  
  output$tablecorum<-renderReactable({
    reactable(corum_df(), filterable = TRUE)
  })
  
  output$tablecorumOneComplex<-renderReactable({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
    fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
    qname1<-paste("q.val_",con2,"-",con1,sep="")
    medname <- paste("median_log2FC_",con2,"-",con1,sep="")
    
    df_out<-complete_corum_df2()%>%
      # dplyr::select(-type,-X,-Reference,-p.val,-Description)%>%
      dplyr::select(ComplexRank,ComplexName,.data[[medname]],total_complexMembers,complexMembersIDed,Gene.Symbol,.data[[sig1]],Annotation,.data[[fc1]],.data[[qname1]])%>%
      dplyr::filter(ComplexRank == input$complexranknum)
      
    
    reactable(df_out, filterable = TRUE)
  })
  
  # tablecorumOneComplexprint<- reactive({
  #   con1<-as.character(unique(norm_sum_final()$condition)[input$numcon1])
  #   con2<-as.character(unique(norm_sum_final()$condition)[input$numcon2])
  #   fc1<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
  #   sig1<-as.character(paste("sig_",con2,"-",con1,sep=""))
  #   qname1<-paste("q.val_",con2,"-",con1,sep="")
  #   medname <- paste("median_log2FC_",con2,"-",con1,sep="")
  #   
  #   df_out<-complete_corum_df2()%>%
  #     # dplyr::select(-type,-X,-Reference,-p.val,-Description)%>%
  #     dplyr::select(ComplexRank,ComplexName,.data[[medname]],total_complexMembers,complexMembersIDed,Gene.Symbol,.data[[sig1]],Annotation,.data[[fc1]],.data[[qname1]])
  #   # 
  #   # %>%
  #   #   dplyr::filter(ComplexRank == input$complexranknum)
  #   
  #   
  #   reactable(df_out, filterable = TRUE)
  # })
  
  go_df<-reactive({
    filt_ttest2<-tableresults()%>%
      dplyr::ungroup()%>%
      #tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
      dplyr::select(Reference,cluster)%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster)))%>%
      dplyr::mutate(set = if_else(cluster == input$numcluster, "sig","n.s."))
    
    
    
    df.bg <- filt_ttest2 %>% 
      plotly::distinct(Reference) %>% 
      dplyr::left_join(go.df()) %>% 
      dplyr::mutate(all = n_distinct(Reference)) %>% 
      dplyr::group_by(GO,GO_cat) %>% 
      dplyr::summarise("pos" = n_distinct(Reference),
                "neg" = all - pos) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(selection="background") %>% 
      plotly::ungroup()
    # print("first")
    # print(head(df.bg))
    go.df <- df.bg %>%
      dplyr::select(GO,GO_cat) %>%
      dplyr::left_join(go.df())
    # print("second")
    # print(head(go.df))
    df.fg <- filt_ttest2 %>%  
      dplyr::mutate(selection = "foreground") %>% 
      dplyr::distinct(Reference, selection, set) %>% 
      dplyr::right_join(go.df) %>% 
      dplyr::group_by(selection, set) %>% 
      dplyr::mutate(term_all = n_distinct(Reference)) %>% 
      dplyr::group_by(GO,GO_cat, selection, set) %>% 
      dplyr::summarise("pos" = n_distinct(Reference),
                "neg" = term_all - pos) %>% 
      plotly::ungroup() %>% 
      dplyr::filter(set == "sig")%>%
      dplyr::distinct()
    # print("third")
    # print(head(df.fg))
    df.nest <- df.bg %>% dplyr::left_join(df.fg %>% dplyr::distinct(GO,GO_cat, set)) %>% 
      dplyr::distinct() %>% 
      dplyr::full_join(df.fg) %>% 
      tidyr::drop_na(set) %>% 
      tidyr::drop_na(GO) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(GO,GO_cat, set) %>% 
      tidyr::nest() 
    # print("fourth")
    # print(head(df.nest))
    
    exact_pval_fun<-function(data){
      dataframe<-data.frame(data) %>% tibble::column_to_rownames(var = "selection")
      matrix<-as.matrix(dataframe) 
      exact_test_result<-fisher.test(matrix, alternative = "two.sided")
      return(exact_test_result$p.value)
    }
    
    ## function to perform fisher test and get enrichment
    enrich_fun<-function(data){
      dataframe<-data.frame(data) %>% 
        dplyr::mutate(frac = pos/neg) %>% 
        dplyr::select(-c(pos,neg)) %>% 
        tidyr::pivot_wider(names_from = selection, values_from = frac) %>% 
        dplyr::summarise(enrichment = foreground/background) 
      
      return(dataframe$enrichment[1])
    }
    
    
    fish.df <- df.nest %>% 
      dplyr::mutate(p.val = data %>% purrr::map_dbl(exact_pval_fun)) %>% #calc p val
      dplyr::mutate(enrich = data %>% purrr::map_dbl(enrich_fun)) %>%   #calc enrichment
      dplyr::select(-c(data)) %>%
      dplyr::filter(GO != "")%>%
      dplyr::arrange(p.val)
    # print("fifth")
    
    
    fish.df$p.val.adj = p.adjust(fish.df$p.val,method="fdr")
    
    fish.df<-fish.df%>%
      dplyr::arrange(p.val.adj)%>%    
      tibble::rownames_to_column("GORank")%>%
      dplyr::mutate(GORank=as.numeric(as.character(GORank)))
    # write.csv(fish.df,"C:/Harper/Side_analysis/go_output_test.csv")
    # print("sixth")
    # fish.df<-fish.df%>%
    #   dplyr::mutate(p.val.adj = p.adjust(p.val, method = "fdr"))
    
    fish.df.sig<-fish.df%>%
      dplyr::mutate(sig_GO = if_else(p.val.adj<input$numgo,"sig","n.s."))%>%
      dplyr::mutate(GO_id = str_sub(GO,-11,-2),
             GO_shrink1 = str_sub(GO,1,-13),
             GO_shrink = str_sub(GO_shrink1,1,70))
    
    fish.df.sig$GO_cat <- factor(fish.df.sig$GO_cat, levels = c("molecular_function","cellular_component","biological_process"))
    # print(head(fish.df.sig))
    # write.csv(fish.df,"C:/Harper/Side_analysis/go_output_test2.csv")
    return(fish.df.sig)
  })
  
  go_df2<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    filt_ttest<-tableresults()%>%
      dplyr::ungroup()%>%
      #tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")%>%
      dplyr::select(Reference,.data[[sig_name]])
    
    if (input$txtupdown == "up") {
      filt_ttest2 <- filt_ttest%>%
        dplyr::mutate(set = if_else(.data[[sig_name]] == "up", "sig","n.s."))
    }
    
    if (input$txtupdown == "down") {
      filt_ttest2 <- filt_ttest%>%
        dplyr::mutate(set = if_else(.data[[sig_name]] == "down", "sig","n.s."))
    }
    
    if (input$txtupdown == "both") {
      filt_ttest2 <- filt_ttest%>%
        dplyr::mutate(set = if_else(.data[[sig_name]] == "down", "sig",if_else(.data[[sig_name]] == "up","sig","n.s.")))
    }
    

    
    
    
    df.bg <- filt_ttest2 %>% 
      plotly::distinct(Reference) %>% 
      dplyr::left_join(go.df()) %>% 
      dplyr::mutate(all = n_distinct(Reference)) %>% 
      dplyr::group_by(GO,GO_cat) %>% 
      dplyr::summarise("pos" = n_distinct(Reference),
                       "neg" = all - pos) %>% 
      dplyr::distinct() %>% 
      dplyr::mutate(selection="background") %>% 
      plotly::ungroup()
    go.df <- df.bg %>%
      dplyr::select(GO,GO_cat) %>%
      dplyr::left_join(go.df())
    df.fg <- filt_ttest2 %>%  
      dplyr::mutate(selection = "foreground") %>% 
      dplyr::distinct(Reference, selection, set) %>% 
      dplyr::right_join(go.df) %>% 
      dplyr::group_by(selection, set) %>% 
      dplyr::mutate(term_all = n_distinct(Reference)) %>% 
      dplyr::group_by(GO,GO_cat, selection, set) %>% 
      dplyr::summarise("pos" = n_distinct(Reference),
                       "neg" = term_all - pos) %>% 
      plotly::ungroup() %>% 
      dplyr::filter(set == "sig")%>%
      dplyr::distinct()
    df.nest <- df.bg %>% dplyr::left_join(df.fg %>% dplyr::distinct(GO,GO_cat, set)) %>% 
      dplyr::distinct() %>% 
      dplyr::full_join(df.fg) %>% 
      tidyr::drop_na(set) %>% 
      tidyr::drop_na(GO) %>% 
      dplyr::distinct() %>% 
      dplyr::group_by(GO,GO_cat, set) %>% 
      tidyr::nest() 
    
    exact_pval_fun<-function(data){
      dataframe<-data.frame(data) %>% tibble::column_to_rownames(var = "selection")
      matrix<-as.matrix(dataframe) 
      exact_test_result<-fisher.test(matrix, alternative = "two.sided")
      return(exact_test_result$p.value)
    }
    
    ## function to perform fisher test and get enrichment
    enrich_fun<-function(data){
      dataframe<-data.frame(data) %>% 
        dplyr::mutate(frac = pos/neg) %>% 
        dplyr::select(-c(pos,neg)) %>% 
        tidyr::pivot_wider(names_from = selection, values_from = frac) %>% 
        dplyr::summarise(enrichment = foreground/background) 
      
      return(dataframe$enrichment[1])
    }
    
    
    fish.df <- df.nest %>% 
      dplyr::mutate(p.val = data %>% purrr::map_dbl(exact_pval_fun)) %>% #calc p val
      dplyr::mutate(enrich = data %>% purrr::map_dbl(enrich_fun)) %>%   #calc enrichment
      dplyr::select(-c(data)) %>%
      dplyr::filter(GO != "")%>%
      dplyr::arrange(p.val)
    
    
    fish.df$p.val.adj = p.adjust(fish.df$p.val,method="fdr")
    
    fish.df<-fish.df%>%
      dplyr::arrange(p.val.adj)%>%    
      tibble::rownames_to_column("GORank")%>%
      dplyr::mutate(GORank=as.numeric(as.character(GORank)))
    
    # fish.df<-fish.df%>%
    #   dplyr::mutate(p.val.adj = p.adjust(p.val, method = "fdr"))
    
    fish.df.sig<-fish.df%>%
      dplyr::mutate(sig_GO = if_else(p.val.adj<input$numgo,"sig","n.s."))%>%
      dplyr::mutate(GO_id = str_sub(GO,-11,-2),
                    GO_shrink1 = str_sub(GO,1,-13),
                    GO_shrink = str_sub(GO_shrink1,1,70))
    
    fish.df.sig$GO_cat <- factor(fish.df.sig$GO_cat, levels = c("molecular_function","cellular_component","biological_process"))
    return(fish.df.sig)
  })
  
  output$tableGO2<-renderReactable({
    reactable(goExport()%>%
                dplyr::ungroup()%>%
                dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj)%>%
                dplyr::distinct(), filterable = TRUE)
  })
  
  output$tableGO2reg<-renderReactable({
    reactable(goExport2()%>%
                dplyr::ungroup()%>%
                dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj)%>%
                dplyr::distinct(), filterable = TRUE)
  })
  
  # output$tableGO3<-renderReactable({
  #   reactable(goExport()%>%
  #               dplyr::ungroup()%>%
  #               dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj)%>%
  #               dplyr::distinct(), filterable = TRUE)
  # })
  
  
  goExport<-reactive({
    filt_ttest2<-tableresults()%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
    
    print(head(filt_ttest2))
    print(head(go.df()))
      #dplyr::select(Reference, Gene.Symbol,q.val,log2_foldchange,Annotation,Significant)
    annotate_go<-dplyr::inner_join(go.df()%>%dplyr::select(-Gene.names...primary..),filt_ttest2,by="Reference")
    go_output<-dplyr::inner_join(go_df(),annotate_go,by=c("GO","GO_cat"),suffix=c("_go","_protein"))
    return(go_output)
  })
  
  
  goExport2<-reactive({
    filt_ttest2<-tableresults()%>%
      tidyr::separate(Reference,into = c("type","Reference","Description"),sep="\\|")
    #dplyr::select(Reference, Gene.Symbol,q.val,log2_foldchange,Annotation,Significant)
  
    annotate_go<-dplyr::inner_join(go.df()%>%dplyr::select(-Gene.names...primary..),filt_ttest2,by="Reference")
    go_output<-dplyr::inner_join(go_df2(),annotate_go,by=c("GO","GO_cat"),suffix=c("_go","_protein"))
    return(go_output)
  })
  
  

  
  output$tableGO<-renderReactable({
    reactable(goExport()%>%dplyr::ungroup(),filterable = TRUE)
              # %>%
              #   dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj,Gene.Symbol,Reference,log2_foldchange,q.val,Significant,Annotation), filterable = TRUE)
  })
  
  output$GOclusterdownload <- downloadHandler(
    filename = function() {
      paste("All_GO_results_Cluster",input$numcluster, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(goExport(), file, row.names = FALSE)
    }
  )
  
  
  output$tableGOreg<-renderReactable({
    reactable(goExport2()%>%dplyr::ungroup(),filterable = TRUE)
    # %>%
    #   dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj,Gene.Symbol,Reference,log2_foldchange,q.val,Significant,Annotation), filterable = TRUE)
  })
  
  output$GOregulateddownload <- downloadHandler(
    filename = function() {
      paste("All_GO_results_Regulated_",as.character(unique(norm_sum_final()$condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final()$condition)[input$numpriority1]),"_",input$txtupdown, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(goExport2(), file, row.names = FALSE)
    }
  )
  
  output$plotGO <- renderPlotly({
    p10101<-ggplot(go_df(),  aes(enrich, -log10(p.val.adj), color=sig_GO,label=GO_id,label1=GO_shrink,key=GORank))+geom_point(alpha=0.75)+theme_classic()+
      geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+ 
      labs(x="Enrichment", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
    ggplotly(p10101)%>% layout(dragmode = "select")
    
  })
  
  plotGOprint <- reactive({
    #,label=GO_id,label1=GO_shrink,key=GORank
    return(ggplot(go_df(),  aes(enrich, -log10(p.val.adj), color=sig_GO))+geom_point(alpha=0.75)+theme_classic()+
      geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+ 
      labs(x="Enrichment", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18)))
  })
  
  output$printplotGO <- downloadHandler(
    filename = function() { paste("GO_plot_AllTerms_Cluster",input$numcluster, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotGOprint(), width =8, height = 8)
    }
  )
  
  output$plotGO2 <- renderPlotly({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    p10101<-ggplot(go_df2(),  aes(enrich, -log10(p.val.adj), color=sig_GO,label=GO_id,label1=GO_shrink,key=GORank))+geom_point(alpha=0.75)+theme_classic()+
      geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+ 
      labs(x="Enrichment", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
    ggplotly(p10101)%>% layout(dragmode = "select")
    
  })
  

  
  plotGO2print <- reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    #,label=GO_id,label1=GO_shrink,key=GORank
    return(ggplot(go_df2(),  aes(enrich, -log10(p.val.adj), color=sig_GO))+geom_point(alpha=0.75)+theme_classic()+
      geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+ 
      labs(x="Enrichment", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18)))
    # ggplotly(p10101)%>% layout(dragmode = "select")
    
  })
  
  output$printplotGO2 <- downloadHandler(
    filename = function() { paste("GO_plot_AllTerms_Regulated_",as.character(unique(norm_sum_final()$condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final()$condition)[input$numpriority1]),"_",input$txtupdown, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotGO2print(), width = 8, height = 8)
    }
  )
  
  
  output$plotGOuser <-renderPlot({
    req(input$txtlistnum)
    req(input$numgo2)

    list_index_ranks<-unlist(strsplit(input$txtlistnum, ","))

    ggplot(goExport()%>%filter(GORank %in% list_index_ranks),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+
      labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
      guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo2),linetype="dashed")


  })
  
  plotGOuserprint <- reactive({
    req(input$txtlistnum)
    req(input$numgo2)
    
    list_index_ranks<-unlist(strsplit(input$txtlistnum, ","))
    
    return(ggplot(goExport()%>%filter(GORank %in% list_index_ranks),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+
      labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
      guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo2),linetype="dashed"))
    
    
  })
  
  
  output$printplotGOuser <- downloadHandler(
    filename = function() { paste("GO_plot_UserTerms_Cluster",input$numcluster, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotGOuserprint(), width = 12, height = 8)
    }
  )
  
  output$plotGOuserReg <-renderPlot({
    req(input$txtlistnum3)
    req(input$numgo3)
    
    list_index_ranks2<-unlist(strsplit(input$txtlistnum3, ","))
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    ggplot(goExport2()%>%filter(GORank %in% list_index_ranks2),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+
      labs(x="GO terms", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
      guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo3),linetype="dashed")
    
    
  })
  
  plotGOuserRegprint <-reactive({
    req(input$txtlistnum3)
    req(input$numgo3)
    
    list_index_ranks2<-unlist(strsplit(input$txtlistnum3, ","))
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    return(ggplot(goExport2()%>%filter(GORank %in% list_index_ranks2),  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+
      labs(x="GO terms", y="-log10(q-value)",title = paste(sig_name,"_",input$txtupdown,sep = ""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))+
      guides(colour = guide_legend(override.aes = list(size=3)))+geom_hline(yintercept = -log10(input$numgo3),linetype="dashed"))
    
    
  })
  
  output$printplotGOuserReg <- downloadHandler(

    filename = function() { paste("GO_plot_UserTerms_Regulated_",as.character(unique(norm_sum_final()$condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final()$condition)[input$numpriority1]),"_",input$txtupdown, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotGOuserRegprint(), width = 12, height = 8)
    }
  )

  
  output$plotGOtopn2 <-renderPlotly({
    
    p333<-ggplot(go_df()[1:input$numtopn,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat,key=GORank))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+ 
      labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))
    # +
    #   guides(colour = guide_legend(override.aes = list(size=3)))
    ggplotly(p333)%>% layout(dragmode = "select")
    
  })
  
  plotGOtopn2print <-reactive({
    #,key=GORank
    return(ggplot(go_df()[1:input$numtopn,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+ 
      labs(x="GO terms", y="-log10(q-value)",title=paste("Cluster ",input$numcluster, sep=""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12)))
    # +
    #   guides(colour = guide_legend(override.aes = list(size=3)))
    # ggplotly(p333)%>% layout(dragmode = "select")
    
  })
  
  output$printplotGOtopn2 <- downloadHandler(
    filename = function() { paste("GO_plot_Top",input$numtopn,"_Cluster",input$numcluster,"_GOterms", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotGOtopn2print(), width = 12, height = 8)
    }
  )
  
  output$plotGOtopn3 <-renderPlotly({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    p333<-ggplot(go_df2()[1:input$numtopn1,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat,key=GORank))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+ 
      labs(x="GO terms", y="-log10(q-value)", title = paste(sig_name,"_",input$txtupdown,sep = ""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12))
    # +
    #   guides(colour = guide_legend(override.aes = list(size=3)))
    ggplotly(p333)%>% layout(dragmode = "select")
    
  })
  
  plotGOtopn3print <-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numpriority1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numpriority2])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    #,key=GORank
    return(ggplot(go_df2()[1:input$numtopn1,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
      geom_point()+theme_classic()+coord_flip()+
      scale_color_viridis_d(end=0.8)+ 
      labs(x="GO terms", y="-log10(q-value)", title = paste(sig_name,"_",input$txtupdown,sep = ""))+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 12)))
    # +
    #   guides(colour = guide_legend(override.aes = list(size=3)))
    # ggplotly(p333)%>% layout(dragmode = "select")
    
  })
  
  output$printplotGOtopn3 <- downloadHandler(
    
    filename = function() { paste("GO_plot_Top",input$numtopn1,"_Regulated_",as.character(unique(norm_sum_final()$condition)[input$numpriority2]),"-",as.character(unique(norm_sum_final()$condition)[input$numpriority1]),"_",input$txtupdown, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotGOtopn3print(), width = 12, height = 8)
    }
  )
  
  
  output$goclick <- renderReactable({
    d <- event_data("plotly_click")
    req(d)
    df_sub2<-goExport()%>%dplyr::filter(GORank %in% d$key )%>%
      #dplyr::rename(!!unique(norm_sum_final()$condition)[1]:=cond1,!!unique(norm_sum_final()$condition)[2]:=cond2)%>%
      # dplyr::select(-data)%>%
      dplyr::ungroup()
    # %>%
    #   dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj,Gene.Symbol,Reference,log2_foldchange,q.val,Significant)
    
    reactable(df_sub2, filterable = TRUE)
  })
  
  
  output$goclick3 <- renderReactable({
    d <- event_data("plotly_click")
    req(d)
    df_sub2<-goExport2()%>%dplyr::filter(GORank %in% d$key )%>%
      #dplyr::rename(!!unique(norm_sum_final()$condition)[1]:=cond1,!!unique(norm_sum_final()$condition)[2]:=cond2)%>%
      # dplyr::select(-data)%>%
      dplyr::ungroup()
    # %>%
    #   dplyr::select(GORank,GO,GO_cat,enrich,p.val.adj,Gene.Symbol,Reference,log2_foldchange,q.val,Significant)
    
    reactable(df_sub2, filterable = TRUE)
  })

  
  # topGOPlotInput <- reactive({
  #   p7233<-ggplot(go_df()[1:input$numtopn,],  aes(reorder(GO_shrink, p.val.adj, FUN=median),-log10(p.val.adj), size=enrich,color=GO_cat))+
  #     geom_point()+theme_classic()+coord_flip()+
  #     scale_color_viridis_d(end=0.8)+ 
  #     labs(x="GO terms", y="-log10(q-value)")+
  #     theme(axis.title = element_text(size=18),axis.text = element_text(size = 16))
  # })
  # 
  # 
  # output$TopGODownloadPlot <- downloadHandler(
  #   filename = function() { paste("TopGO_term_plot", '.pdf', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = topGOPlotInput(), width = 13, height = 10)
  #   }
  # )
  
  
  
  
  # output$alldownloadGO <- downloadHandler(
  #   filename = function() {
  #     paste("all_go_results", ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(go_df(), file, row.names = FALSE)
  #   }
  # )
  
  # GOPlotlyInput <- reactive({
  #   p774<-ggplot(go_df(),  aes(enrich, -log10(p.val.adj), color=sig_GO))+geom_point()+theme_classic()+
  #     geom_hline(yintercept = -log10(input$numgo))+scale_color_viridis_d(end=0.8)+ 
  #     labs(x="Enrichment", y="-log10(q-value)")+
  #     theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
  # })
  # 
  # output$plotGOall <- downloadHandler(
  #   filename = function() { paste("All_GO_Results_Plot", '.pdf', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = GOPlotlyInput(), width = 7, height = 7)
  #   }
  # )
  
  output$corrplotMeans<-renderPlotly({
    x_axis = unique(norm_sum_final()$condition)[1]
    y_axis = unique(norm_sum_final()$condition)[2]
    
    
    plot_df<-ttest_df_final()
    

    
    plot_df<-plot_df%>%
      dplyr::rename(!!unique(norm_sum_final()$condition)[1]:=cond1,!!unique(norm_sum_final()$condition)[2]:=cond2)%>%
      dplyr::select(-data)
    x_axis<-unique(norm_sum_final()$condition)[1]
    y_axis<-unique(norm_sum_final()$condition)[2]
    
    
    
    p4021<-ggplot(data=plot_df, aes(.data[[unique(norm_sum_final()$condition)[1]]],.data[[unique(norm_sum_final()$condition)[2]]],label=Gene.Symbol,color=Significant))+
      geom_point(alpha=input$alphanum1)+theme_classic()+
      geom_abline(slope = 1,color="turquoise")+
      labs(x=x_axis, y=y_axis)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))+scale_color_viridis_d(end=0.8)
    ggplotly(p4021)
  })
  
  
  
  output$rankfoldchange<-renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    sortedDF<-tableresults()%>%
      dplyr::select(Gene.Symbol,Reference,.data[[foldchangevalue]])%>%
      dplyr::distinct()%>%
      dplyr::arrange(desc(.data[[foldchangevalue]]))%>%
      tibble::rownames_to_column("colnum")%>%
      dplyr::mutate(colnum=as.numeric(as.character(colnum)))
    
    top_n<-sortedDF$Reference[1:input$num55]
    last_num<-max(sortedDF$colnum)
    last_num_min<-last_num-input$num55+1
    bottom_n<-sortedDF$Reference[last_num_min:last_num]
    
    
    ggplot(data=sortedDF, aes(colnum,.data[[foldchangevalue]]))+
      geom_point()+
      geom_text_repel(data=sortedDF%>%dplyr::mutate(label=dplyr::if_else(Reference %in% top_n,Gene.Symbol,dplyr::if_else(Reference %in% bottom_n,Gene.Symbol,""))),aes(label=label),max.overlaps = Inf,box.padding = 1)+
      theme_classic()+labs(x="Protein Rank",y=paste(foldchangevalue))
  })
  
  rankfoldchangeprint<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    sortedDF<-tableresults()%>%
      dplyr::select(Gene.Symbol,Reference,.data[[foldchangevalue]])%>%
      dplyr::distinct()%>%
      dplyr::arrange(desc(.data[[foldchangevalue]]))%>%
      tibble::rownames_to_column("colnum")%>%
      dplyr::mutate(colnum=as.numeric(as.character(colnum)))
    
    top_n<-sortedDF$Reference[1:input$num55]
    last_num<-max(sortedDF$colnum)
    last_num_min<-last_num-input$num55+1
    bottom_n<-sortedDF$Reference[last_num_min:last_num]
    
    
    p1<-ggplot(data=sortedDF, aes(colnum,.data[[foldchangevalue]]))+
      geom_point()+
      geom_text_repel(data=sortedDF%>%dplyr::mutate(label=dplyr::if_else(Reference %in% top_n,Gene.Symbol,dplyr::if_else(Reference %in% bottom_n,Gene.Symbol,""))),aes(label=label),max.overlaps = Inf,box.padding = 1)+
      theme_classic()+labs(x="Protein Rank",y=paste(foldchangevalue))
    return(p1)
  })
  
  output$rankfoldchangeprint2 <- downloadHandler(
    filename = function() { paste("Top_rankedUpDownFC_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = rankfoldchangeprint(), width = 13, height = 10)
    }
  )
  
  

  output$plot6 <- renderPlotly({
    x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    p1<-ggplot(data=ttest_df_final(), aes(log2_foldchange, -log2(q.val),color=Significant, label=Gene.Symbol,key=Reference))+
      geom_point(alpha=input$alphanum1)+theme_classic()+
      geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
      geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+ 
      labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
    
    
    ggplotly(p1) %>% layout(dragmode = "select")
    #,width = 1000,height=650
    
  })
  
  output$brush <- renderReactable({
    d <- event_data("plotly_selected")
    req(d)
    df_sub<-ttest_df_final()%>%dplyr::filter(Reference %in% d$key )%>%
      dplyr::rename(!!unique(norm_sum_final()$condition)[1]:=cond1,!!unique(norm_sum_final()$condition)[2]:=cond2)%>%
      dplyr::select(-data)
    
    reactable(df_sub, filterable = TRUE)
  })
  
  volcanoPlotInput <- reactive({
    x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    p771<-ggplot(data=ttest_df_final(), aes(log2_foldchange, -log2(q.val),color=Significant))+
      geom_point(alpha=input$alphanum1)+theme_classic()+
      geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
      geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+ 
      labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
  })
  
  # output$volcanoDownloadPlot <- downloadHandler(
  #   filename = function() { paste("Volcano_NoLabel_Plot", '.png', sep='') },
  #   content = function(file) {
  #     device <- function(..., width, height) grDevices::png(..., width = 12, height = 6.5, res = 300, units = "in")
  #     ggsave(file, plot = volcanoPlotInput(), device = device)
  #   }
  # )
  
  # output$volcanoDownloadPlot <- downloadHandler(
  #   filename = function() { paste("Volcano_NoLabel_Plot", '.pdf', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = volcanoPlotInput(), width = 12, height = 6.5)
  #   }
  # )
  
  
  output$ploteuler <- renderPlot({
    df<-tableresults()%>%
      dplyr::ungroup()%>%
      dplyr::select(contains("sig_"))
    
    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }

    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }
    
    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }
    
    a<-as.matrix(df2)
    
    a[is.na(a)] = F
    
    v = gsub("sig_","",colnames(a)) #Specific nomenclature
    fit = eulerr::euler(a)
    #brewer.pal(ncol(a),"Set3")
    plot(fit,
                   quantities = list(cex = 0.5),
                   fills = list(fill = brewer.pal(ncol(a),"Set3"), alpha = 0.8, cex=0.7),
                   # lty = 1:3,
                   adjust_labels=T,
                   labels = list(labels=v, font = 12, cex=1),
                   main = paste("Venn Diagram:",input$txtupdownoverlap))
    
    
  })
  
  ploteulerprint <- reactive({
    df<-tableresults()%>%
      dplyr::ungroup()%>%
      dplyr::select(contains("sig_"))
    
    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }
    
    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }
    
    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }
    
    a<-as.matrix(df2)
    a[is.na(a)] = F
    v = gsub("sig_","",colnames(a)) #Specific nomenclature
    fit = eulerr::euler(a)
    #brewer.pal(ncol(a),"Set3")
    return(plot(fit,
         quantities = list(cex = 0.5),
         fills = list(fill = brewer.pal(ncol(a),"Set3"), alpha = 0.8, cex=0.7),
         # lty = 1:3,
         adjust_labels=T,
         labels = list(labels=v, font = 12, cex=1),
         main = paste("Venn Diagram:",input$txtupdownoverlap)))
    
    
  })
  
  output$printploteuler <- downloadHandler(
    filename = function() { paste("VennDiagram_",input$txtupdownoverlap, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = ploteulerprint(), width = 8, height = 8)
    }
  )
  
  
  
  output$plotupset <-renderPlot({
    df<-tableresults()%>%
      dplyr::ungroup()%>%
      dplyr::select(contains("sig_"))
    
    colnames(df) <- gsub("sig_","",colnames(df)) 
    
    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }
    
    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }
    
    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }
    
    a<-as.matrix(df2)
    a = a[as.logical(rowSums(a)),] #Remove proteins without group

    plot3.3 = ComplexUpset::upset(
      as.data.frame(a),colnames(a), name = paste("Upset Plot:",input$txtupdownoverlap),
       min_size = 1
    )
    #group_by= "sets",
    
    plot3.3
  })
  
  plotupsetprint <- reactive({
    df<-tableresults()%>%
      dplyr::ungroup()%>%
      dplyr::select(contains("sig_"))
    
    colnames(df) <- gsub("sig_","",colnames(df)) 
    
    if (input$txtupdownoverlap == "both") {
      df2 <- (df == "up") | (df == "down")
    }
    
    if (input$txtupdownoverlap == "up") {
      df2 <- (df == "up")
    }
    
    if (input$txtupdownoverlap == "down") {
      df2 <- (df == "down")
    }
    
    a<-as.matrix(df2)
    a = a[as.logical(rowSums(a)),] #Remove proteins without group
    plot3.3 = upset(
      as.data.frame(a),colnames(a), name = paste("Upset Plot:",input$txtupdownoverlap),
      min_size = 1
    )
    
    return(plot3.3)
  })
  
  output$printplotupset <- downloadHandler(
    filename = function() { paste("UpsetPlot_",input$txtupdownoverlap, '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotupsetprint(), width = 12, height = 6.5)
    }
  )
  
  
  
  
  
  
  
  output$plotVL <- renderPlot({

    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse7vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse8vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    
    ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
      geom_point(alpha=input$alphanum2,size=0.5)+
      geom_text_repel(data=.%>%
                         mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(label=label),
                       max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
      geom_vline(xintercept=log2(input$num78),linetype="dashed")+
      geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
      labs(x= paste(foldchangevalue), y= "-Log10(q-value)")

  },width = 900,height=650)
  
  plotVLprint <- reactive({
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse7vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse8vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    
    p5<-ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[foldchangevalue]], -log10(.data[[qvaluename]]),color=.data[[sig_name]]))+
      geom_point(alpha=input$alphanum2,size=0.5)+
      geom_text_repel(data=.%>%
                        mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(label=label),
                      max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
      geom_vline(xintercept=log2(input$num78),linetype="dashed")+
      geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
      labs(x= paste(foldchangevalue), y= "-Log10(q-value)")
    
    return(p5)
    
  })
  
  output$volcanoDownloadPlot <- downloadHandler(
    filename = function() { paste("Volcano_Label_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse8vol]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse7vol]), '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotVLprint(), width = 12, height = 6.5)
    }
  )
  
  
  uniqueColProtein<-reactive({
    df_unique_val <-ttest_df_untidy_final()%>%
      dplyr::ungroup()%>%
      dplyr::select(Gene.Symbol)%>%
      plotly::distinct()
  })
  

  
  
  
  output$plot5 <- renderPlot({
    req(input$text11)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    
    gene_name_input <- input$text11
    plot_protein<-ttest_df_untidy_final()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)%>%
      dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))
    
    plot_protein_v2<-plot_protein%>%dplyr::group_by(condition,ProtID2)%>%dplyr::summarise(log2_int=median(log2_int))
    
    ggplot(plot_protein,aes(x=condition,2^log2_int,fill = condition,color=condition))+
      geom_col(data=plot_protein_v2,alpha=0.5)+
      geom_jitter(size=3)+
      theme_classic()+scale_fill_viridis_d(end=0.8)+scale_color_viridis_d(end=0.8)+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5))+
      labs(title = gene_name_input, x="",y="MS Intensity")+facet_wrap(~ProtID2,nrow=1,scales="free_y")

  })
  
  individProteinPlotInput <- reactive({
    req(input$text11)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1vol])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2vol])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    
    gene_name_input <- input$text11
    plot_protein<-ttest_df_untidy_final()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)%>%
      dplyr::mutate(ProtID2 = paste(Gene.Symbol,Reference,sep="_"))
    
    plot_protein_v2<-plot_protein%>%dplyr::group_by(condition,ProtID2)%>%dplyr::summarise(log2_int=median(log2_int))
    
    ggplot(plot_protein,aes(x=condition,2^log2_int,fill = condition,color=condition))+
      geom_col(data=plot_protein_v2,alpha=0.5)+
      geom_jitter(size=3)+
      theme_classic()+scale_fill_viridis_d(end=0.8)+scale_color_viridis_d(end=0.8)+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=24,hjust=0.5))+
      labs(title = gene_name_input, x="",y="MS Intensity")+facet_wrap(~ProtID2,nrow=1,scales="free_y")
  })
  
  
  output$proteinBarplotDownloadPlot <- downloadHandler(
    filename = function() { paste(input$text11,"_SingleProtein_Barplot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = individProteinPlotInput(), width = 10, height = 7)
    }
  )
  
  volcanoPlot2Input <- reactive({
    x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    p4000<-ggplot(data=ttest_df_final(), aes(log2_foldchange, -log2(q.val),color=Significant))+
      geom_point(alpha=input$alphanum1)+geom_text_repel(data=.%>%
                                     mutate(label=if_else(Significant != "n.s.",Gene.Symbol,"")),aes(label=label),
                                   max.overlaps = Inf,show.legend=FALSE,box.padding = 1)+theme_classic()+
      geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
      geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+ 
      labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18))
  })
  

  
  output$volcanoDownloadPlot2 <- downloadHandler(
    filename = function() { paste("Volcano_Label_Plot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = volcanoPlot2Input(), width = 12, height = 6.5)
    }
  )
  
  
  output$text12 <- renderText({
    gene_name_input <- input$text11
    plot_protein<-ttest_df_untidy_final()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)
    annotation_plot<-unique(plot_protein$Annotation)
    return(annotation_plot)
  })
  output$text13 <- renderText({
    gene_name_input <- input$text11
    plot_protein<-ttest_df_untidy_final()%>%
      dplyr::filter(Gene.Symbol==gene_name_input)
    ref_plot<-unique(plot_protein$Reference)
    return(ref_plot)
  })
  
  outputdata<-reactive({
    reps_df<-correlation_plot2()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    all_table<-left_join(ttest_slim(),reps_df, by= c("Reference","Gene.Symbol","Annotation"))
  })
    
    
  output$table1 <- renderReactable({
    cutoff<-input$num1
    fc_cutoff<-input$num2
    
    reps_df<-correlation_plot2()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    all_table<-left_join(ttest_slim(),reps_df, by= c("Reference","Gene.Symbol","Annotation"))
    
    slim_df_for_output<-all_table%>%
      dplyr::filter(q.val < cutoff, abs(log2_foldchange)>log2(fc_cutoff))
    reactable(slim_df_for_output, filterable = TRUE)
  })
  
  output$table16 <- renderReactable({
    # cutoff<-input$num1
    # fc_cutoff<-input$num2
    # 
    # results_df<-ttest_slim()%>%
    #   tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    reps_df<-correlation_plot2()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    HotelClusterDF <- dplyr::left_join(ttest_slim(),clusterdesignation,by=c("Reference","Gene.Symbol","Annotation"))%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
    
    all_table<-left_join(HotelClusterDF,reps_df, by= c("Reference","Gene.Symbol","Annotation"))
    # 
    # slim_df_for_output<-all_table%>%
    #   dplyr::filter(q.val < cutoff, abs(log2_foldchange)>log2(fc_cutoff))
    reactable(all_table, filterable = TRUE)
  })
  
  tableresults <- reactive({
    # cutoff<-input$num1
    # fc_cutoff<-input$num2
    # 
    reps_df<-correlation_plot2()%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")

    clusterdesignation<-clusterTraceDFinitial()%>%
      dplyr::ungroup()%>%
      dplyr::select(ProtID,cluster)%>%
      tidyr::separate(ProtID,c("Reference","Gene.Symbol","Annotation"), sep="__X__")
    
    HotelClusterDF <- dplyr::left_join(ttest_slim(),clusterdesignation,by=c("Reference","Gene.Symbol","Annotation"))%>%
      dplyr::mutate(cluster = as.numeric(as.character(cluster))) %>%
      # the good stuff here
      dplyr::mutate_if(is.numeric,dplyr::coalesce,0)
    
    all_table<-left_join(HotelClusterDF,reps_df, by= c("Reference","Gene.Symbol","Annotation"))
    # all_table<-left_join(ttest_slim(),reps_df, by= c("Reference","Gene.Symbol","Annotation"))
    # 
    # slim_df_for_output<-all_table%>%
    #   dplyr::filter(q.val < cutoff, abs(log2_foldchange)>log2(fc_cutoff))
    return(all_table)
  })
  
  
  
  output$alldownload <- downloadHandler(
    filename = function() {
      paste("all_ttest_results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(tableresults(), file, row.names = FALSE)
    }
  )
  ### localization start
  annotation_ttest<-reactive({
    
    if (input$txthumanmouse == "mouse") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol")
    }
    
    if (input$txthumanmouse == "human") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol")
    }
    # ann_combine<-tableresults()%>%
    #   dplyr::left_join(.,locale,by="Gene.Symbol")
    # 
    # if (input$txthumanmouse == "mouse") {
    #   results<-tableresults()%>%
    #     dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
    
    
    return(ann_combine)
  })
  
  annotation_ttest2<-reactive({
    if (input$txthumanmouse == "mouse") {
      ann_combine<-tableresults()%>%
        dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale2,by="Gene.Symbol")
    }
    
    if (input$txthumanmouse == "human") {
      ann_combine<-tableresults()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale2,by="Gene.Symbol")
    }
      
      
    # ann_combine<-tableresults()%>%
    #   dplyr::left_join(.,locale2,by="Gene.Symbol")

    
    
    return(ann_combine)
  })
  
  output$plot7 <- renderPlot({

    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversionloc == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    p1
  })
  
  plot7print <- reactive({
    
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversionloc == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    return(p1)
  })
  
  
  output$printplot7 <- downloadHandler(
    filename = function() { paste("Volcano_Localizationv1_",input$localization,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc]),"-", as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot7print(), width = 12, height = 6.5)
    }
  )
  
  
  
  output$plot7loc2 <- renderPlot({
    
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversionloc5 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc5 == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    p1
  })
  
  plot7loc2print <- reactive({
    
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversionloc5 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con2,"-",con1,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc5 == "Yes") {
      # foldchangevalue2<-as.character(paste("log2FC_",con1,"-",con2,sep=""))
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        #ggplot(tableresults()%>%dplyr::mutate(ProtID = paste(Reference,Gene.Symbol,sep="_")), aes( .data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=.data[[sig_name]],label=Gene.Symbol, key=ProtID))+
        geom_point(alpha=0.5,size=0.5)+theme_classic()+
        scale_color_viridis_d(end=0.8)+
        geom_hline(yintercept = -log10(input$num77),linetype="dashed")+
        geom_vline(xintercept=log2(input$num78),linetype="dashed")+
        geom_vline(xintercept=-log2(input$num78),linetype="dashed")+
        labs(x= paste("log2FC_",con1,"-",con2,sep=""), y= "-Log10(q-value)",title=input$localization5)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    return(p1)
  })
  
  output$printplot7loc2 <- downloadHandler(
    filename = function() { paste("Volcano_Localizationv2_",input$localization5,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5]),"-", as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot7loc2print(), width = 12, height = 6.5)
    }
  )
  # volcanoPlotLocaleInput <- reactive({
  #   x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
  #   
  #   p1400<-ggplot()+
  #     geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(log2_foldchange, -log2(q.val),color=Significant),color="black")+
  #     geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(log2_foldchange, -log2(q.val),color=Significant),color="red")+
  #     theme_classic()+
  #     geom_vline(xintercept = log2(input$num4))+geom_vline(xintercept = -log2(input$num4))+
  #     geom_hline(yintercept = -log2(input$num3))+scale_color_viridis_d(end=0.8)+ 
  #     labs(x=paste("log2(",x_axis,")"), y="-log2(q-value)",title=input$localization)+
  #     theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  # })
  
  output$corrcolorcompartment<-renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    # x_axis = unique(norm_sum_final()$condition)[1]
    # y_axis = unique(norm_sum_final()$condition)[2]

    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  
  # downloadlocalizationdata
  # 
  # input_for_localizationdataexport<- reactive({
  #   
  # })
  # annotation_ttest()
  
  
  output$downloadlocalizationdata <- downloadHandler(
    filename = function() {
      paste("All_results_stats_localizationAdded", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(annotation_ttest(), file, row.names = FALSE)
    }
  )
  
  corrcolorcompartmentprint<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    
    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  
  output$printcorrcolorcompartment <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv1_",input$localization,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc]),"vs", as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentprint(), width = 8, height = 8)
    }
  )
  
  
  output$corrcolorcompartmentloc5<-renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    # x_axis = unique(norm_sum_final()$condition)[1]
    # y_axis = unique(norm_sum_final()$condition)[2]
    
    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  corrcolorcompartmentloc5print<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5])
    
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    
    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con1]], .data[[con2]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con1]], .data[[con2]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=con1, y=con2,title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  output$printcorrcolorcompartmentloc5 <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv2_",input$localization5,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5]),"vs", as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentloc5print(), width = 8, height = 8)
    }
  )
  
  output$corrcolorcompartmentmulti<-renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc2])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc2])
    con3<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse3loc2])
    con4<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse4loc2])
    

    
    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  corrcolorcompartmentmultiprint<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc2])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc2])
    con3<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse3loc2])
    con4<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse4loc2])
    
    
    
    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  output$printcorrcolorcompartmentmulti <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv1_Ratio_",input$localization,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc2]),"_DIV_", as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc2]),"_vs_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse4loc2]),"_DIV_", as.character(unique(norm_sum_final()$condition)[input$numtimecourse3loc2]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentmultiprint(), width = 8, height = 8)
    }
  )
  
  output$corrcolorcompartmentmultiloc5<-renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc6])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc6])
    con3<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse3loc6])
    con4<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse4loc6])
    
    
    
    ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
  })
  
  corrcolorcompartmentmultiloc5print<-reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc6])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc6])
    con3<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse3loc6])
    con4<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse4loc6])
    
    
    
    return(ggplot()+
      geom_abline(slope=1,color="lightblue")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="black")+
      geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization5), aes(.data[[con2]] - .data[[con1]], .data[[con4]] - .data[[con3]]),color="turquoise")+
      theme_classic()+
      scale_color_viridis_d(end=0.8)+
      labs(x=paste(con2, "-",con1,sep=""), y= paste(con4,"-",con3,sep=""),title=input$localization5)+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
  })
  
  output$printcorrcolorcompartmentmultiloc5 <- downloadHandler(
    filename = function() { paste("Correlation_Localizationv2_Ratio_",input$localization5,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc6]),"_DIV_", as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc6]),"_vs_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse4loc6]),"_DIV_", as.character(unique(norm_sum_final()$condition)[input$numtimecourse3loc6]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = corrcolorcompartmentmultiloc5print(), width = 8, height = 8)
    }
  )
  
  
  # output$volcanoLocaleDownloadPlot <- downloadHandler(
  #   filename = function() { paste(input$localization,"_Volcano", '.pdf', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = volcanoPlotLocaleInput(), width = 12, height = 6.5)
  #   }
  # )
  
  
  
  output$plot8 <- renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    x_axis <- paste(con2, "-", con1,sep=" ")
    ggplot(data=annotation_ttest()%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1,input$numlocale2))
  })
  
  plot8print <- reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    x_axis <- paste(con2, "-", con1,sep=" ")
    return(ggplot(data=annotation_ttest()%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1,input$numlocale2)))
  })
  
  
  output$printplot8 <- downloadHandler(
    filename = function() { paste("Localev1_FC_MultipleLocale_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot8print(), width = 16, height = 12)
    }
  )
  
  output$plot8loc2 <- renderPlot({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    x_axis <- paste(con2, "-", con1,sep=" ")
    ggplot(data=annotation_ttest2()%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1loc2,input$numlocale2loc2))
  })
  
  plot8loc2print <- reactive({
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    
    x_axis <- paste(con2, "-", con1,sep=" ")
    return(ggplot(data=annotation_ttest2()%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(reorder(Compartment, .data[[foldchangevalue]], FUN = median),.data[[foldchangevalue]],color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
      # geom_jitter(alpha=0.1)+
      theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
      coord_flip()+
      ylim(c(input$numlocale1loc2,input$numlocale2loc2)))
  })
  
  output$printplot8loc2 <- downloadHandler(
    filename = function() { paste("Localev2_FC_MultipleLocale_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc5]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc5]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot8loc2print(), width = 16, height = 12)
    }
  )
  
  
  
  output$plotmultilocale <- renderPlot({
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    print(unique(dftidy$condition))
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }

    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }

    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }

    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }

    
    ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))

  })
  
  plotmultilocaleprint <- reactive({
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- annotation_ttest()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))

    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }
    
    
    return(ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
    
  })
  
  output$printplotmultilocale <- downloadHandler(
    filename = function() { paste("Localev1_FC_MultipleLocale_AlltoPriority1",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotmultilocaleprint(), width = 16, height = 8)
    }
  )
  
  output$plotmultilocaleloc2 <- renderPlot({
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    print(unique(dftidy$condition))
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }
    
    
    ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    
  })
  
  
  plotmultilocaleloc2print <- reactive({
    numcon<-length(unique(norm_sum_final()$condition))
    
    if (numcon == 3) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]])
    } 
    
    if (numcon  == 4) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]])
    } 
    
    if (numcon  == 5) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]])
    } 
    
    if (numcon  == 6) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]])
    } 
    
    if (numcon  == 7) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]])
    } 
    
    if (numcon  == 8) {
      con1<-as.character(unique(norm_sum_final()$condition)[1])
      con2<-as.character(unique(norm_sum_final()$condition)[2])
      con3<-as.character(unique(norm_sum_final()$condition)[3])
      con4<-as.character(unique(norm_sum_final()$condition)[4])
      con5<-as.character(unique(norm_sum_final()$condition)[5])
      con6<-as.character(unique(norm_sum_final()$condition)[6])
      con7<-as.character(unique(norm_sum_final()$condition)[7])
      con8<-as.character(unique(norm_sum_final()$condition)[8])
      fc1<-as.character(paste(con2,"-",con1,sep=""))
      fc2<-as.character(paste(con3,"-",con1,sep=""))
      fc3<-as.character(paste(con4,"-",con1,sep=""))
      fc4<-as.character(paste(con5,"-",con1,sep=""))
      fc5<-as.character(paste(con6,"-",con1,sep=""))
      fc6<-as.character(paste(con7,"-",con1,sep=""))
      fc7<-as.character(paste(con8,"-",con1,sep=""))
      
      df <- annotation_ttest2()%>%
        dplyr::mutate(!!fc1:=  .data[[con2]] - .data[[con1]],
                      !!fc2:=  .data[[con3]] - .data[[con1]],
                      !!fc3:= .data[[con4]] - .data[[con1]],
                      !!fc4:= .data[[con5]] - .data[[con1]],
                      !!fc5:= .data[[con6]] - .data[[con1]],
                      !!fc6:= .data[[con7]] - .data[[con1]],
                      !!fc7:= .data[[con8]] - .data[[con1]])%>%
        dplyr::ungroup()%>%
        dplyr::select(Reference, Gene.Symbol,Compartment,.data[[fc1]],.data[[fc2]],.data[[fc3]],.data[[fc4]],.data[[fc5]],.data[[fc6]],.data[[fc7]])
    } 
    
    
    
    size_df <- dim(df)[2]
    
    dftidy<-df%>%
      dplyr::group_by(Reference, Gene.Symbol,Compartment)%>%
      tidyr::gather("condition","log2_FC",4:size_df)%>%
      dplyr::ungroup()
    
    
    dftidy<-dftidy%>%
      dplyr::mutate(condition = str_remove(condition,pattern=paste("-",con1,sep="")))
    
    
    if (numcon  == 3) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3))
    }
    
    if (numcon  == 4) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4))
    }
    
    if (numcon  == 5) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5))
    }
    
    if (numcon  == 6) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6))
    }
    
    if (numcon  == 7) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7))
    }
    
    if (numcon  == 8) {
      dftidy$condition <- factor(dftidy$condition , levels=c(con2,con3,con4,con5,con6,con7,con8))
    }
    
    
    return(ggplot(dftidy%>%dplyr::filter(Compartment %in% input$locales5,!is.na(Compartment)), aes(condition,log2_FC,color=Compartment,fill=Compartment))+
      geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+
      geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+theme_classic()+
      geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
      labs(x="", y=paste("Log2 Foldchange to",con1),title="Foldchange by Localization")+
      theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5)))
    
  })
  
  
  output$printplotmultilocaleloc2 <- downloadHandler(
    filename = function() { paste("Localev2_FC_MultipleLocale_AlltoPriority1",'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotmultilocaleloc2print(), width = 16, height = 8)
    }
  )
  
  # multiLocaleInput <- reactive({
  #   x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
  #   p1700<-ggplot(data=annotation_ttest()%>%dplyr::filter(Compartment %in% input$locales,!is.na(Compartment)), aes(reorder(Compartment, log2_foldchange, FUN = median),log2_foldchange,color=Compartment,fill=Compartment))+
  #     geom_violin(size=0.35,draw_quantiles = c(0.25,  0.75),alpha=0)+geom_violin(size=1,draw_quantiles = c( 0.5),alpha=0.4)+
  #     #geom_jitter(alpha=0.25)+
  #     theme_classic()+
  #     geom_hline(yintercept = 0)+scale_color_viridis_d(end=0.8)+ scale_fill_viridis_d(end=0.8)+
  #     labs(x="Localization", y=paste("log2(",x_axis,")"),title="Fold change by localization")+
  #     theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+
  #     coord_flip()+
  #     ylim(c(input$numlocale1,input$numlocale2))
  # })
  
  # 
  # output$multiLocaleDownloadPlot <- downloadHandler(
  #   filename = function() { paste("MultiLocale_ViolinPlot", '.pdf', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = multiLocaleInput(), width = 16, height = 14)
  #   }
  # )
  # 
  
  output$plot44 <- renderPlot({
    # x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc3])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))

    
    if (input$txtinversionloc44 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization2), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=Significant),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=Significant),color="red")+
        geom_text_repel(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2)%>%
                          mutate(label=if_else(Compartment == input$localization2,Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc44 == "Yes") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization2), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        geom_text_repel(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2)%>%
                          mutate(label=if_else(Compartment == input$localization2,Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }

    p1

  },width = 1200,height=650)
  
  
  plot44print <- reactive({
    # x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc3])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc3])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    if (input$txtinversionloc44 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization2), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=Significant),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=Significant),color="red")+
        geom_text_repel(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2)%>%
                          mutate(label=if_else(Compartment == input$localization2,Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc44 == "Yes") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment != input$localization2), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]])),color="red")+
        geom_text_repel(data=annotation_ttest()%>%dplyr::filter(Compartment == input$localization2)%>%
                          mutate(label=if_else(Compartment == input$localization2,Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title=input$localization2)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    return(p1)
    
  })
  
  
  output$printplot44 <- downloadHandler(
    filename = function() { paste("Localizationv1_Volcano_Labels_",input$localization2,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc3]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc3]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot44print(), width = 12, height = 6.5)
    }
  )
  
  output$plot88 <- renderPlot({
    # x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc8])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc8])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    if (input$txtinversionloc88 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc88 == "Yes") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    p1
    
  },width = 1200,height=650)
  
  
  plot88print <- reactive({
    # x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc8])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc8])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    if (input$txtinversionloc88 == "No") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]])),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    if (input$txtinversionloc88 == "Yes") {
      p1<-ggplot()+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment != input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="black")+
        geom_point(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=Significant),color="red")+
        geom_text_repel(data=annotation_ttest2()%>%dplyr::filter(Compartment == input$localization8)%>%
                          mutate(label=if_else(Compartment == input$localization8,Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,color="red")+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title=input$localization8)+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))
    }
    
    return(p1)
    
  })
  
  output$printplot88 <- downloadHandler(
    filename = function() { paste("Localizationv2_Volcano_Labels_",input$localization8,"_",as.character(unique(norm_sum_final()$condition)[input$numtimecourse2loc8]),"-",as.character(unique(norm_sum_final()$condition)[input$numtimecourse1loc8]),'.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plot88print(), width = 12, height = 6.5)
    }
  )

  
  dfuser <- reactive({
    req(input$file4)
    inFile <- input$file4
    if (is.null(inFile))
      return(NULL)
    tbl <- read.csv(inFile$datapath)
    
    return(tbl)
  })
  
  userinput<-reactive({
    req(input$file4)
    ann_combine<-tableresults()%>%
      dplyr::left_join(.,dfuser(),by="Gene.Symbol")%>%
      dplyr::mutate(UserAnnotation = if_else(!is.na(UserAnnotation),UserAnnotation,"__"))
    return(ann_combine)
  })
  
  output$plotuser <- renderPlot({
    req(input$file4)
    req(input$size1)
    req(input$size2)
    req(input$alpha5)
    req(input$alpha6)
    # x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconduser2])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    p1
    
    
  })
  
  plotuserprint <- reactive({
    req(input$file4)
    req(input$size1)
    req(input$size2)
    req(input$alpha5)
    req(input$alpha6)
    # x_axis <- paste(unique(norm_sum_final()$condition)[2], "-", unique(norm_sum_final()$condition)[1],sep=" ")
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconduser2])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "none")%>%
                          mutate(label=if_else(UserAnnotation != "__",Gene.Symbol,"")),aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con2]] - .data[[con1]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con2,"-",con1,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot()+
        geom_point(data=userinput()%>%filter(UserAnnotation == "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha5,size=input$size1)+
        geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(.data[[con1]] - .data[[con2]], -log10(.data[[qvaluename]]),color=UserAnnotation),alpha=input$alpha6,size=input$size2)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num78))+geom_vline(xintercept = -log2(input$num78))+
        geom_hline(yintercept = -log10(input$num77))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC_",con1,"-",con2,sep=""), y="-log10(q-value)",title="User Input")+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),plot.title = element_text(size=24,hjust=0.5))+scale_color_viridis_d(end=0.8)
    }
    return(p1)
    
    
  })
  
  output$printplotuser <- downloadHandler(
    filename = function() { paste("UserDefined_VolcanoPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotuserprint(), width = 12, height = 8)
    }
  )
  
  
  
  
  output$categviolin <- renderPlot({
    
    req(input$file4)
    req(input$numconduser1)
    req(input$numconduser2)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconduser2])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }

    p1
    
    
  })
  
  categviolinprint <- reactive({
    
    req(input$file4)
    req(input$numconduser1)
    req(input$numconduser2)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconduser2])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "Yes") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
                          mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    if (input$txtinversionuser == "No" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con2]] - .data[[con1]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con2]] - .data[[con1]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    if (input$txtinversionuser == "Yes" &&  input$txtlabeluser == "No") {
      p1<-ggplot(userinput()%>%filter(UserAnnotation != "__"), aes(x=UserAnnotation,y=.data[[con1]] - .data[[con2]],color=UserAnnotation,fill=UserAnnotation))+
        geom_hline(yintercept = 0,linetype="dashed")+
        geom_jitter(size=0.75,alpha=0.5,position = position_jitter(seed = 1))+
        # geom_text_repel(data=userinput()%>%dplyr::filter(UserAnnotation != "__")%>%
        #                   mutate(label=if_else(.data[[sig_name]] != "n.s.",Gene.Symbol,"")),aes(UserAnnotation, .data[[con1]] - .data[[con2]],label=label,color=UserAnnotation),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5, position = position_jitter(seed = 1))+
        geom_violin(draw_quantiles = c(0.25,0.5, 0.75),alpha=0.3)+
        theme_classic()+scale_color_viridis_d(end=0.8)+scale_fill_viridis_d(end=0.8)+
        labs(x="", y=paste("Log2FC(",con1,"-",con2,")",sep=""),title="User Input Violin Plot")+
        theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
    }
    
    return(p1)
    
    
  })
  
  output$printcategviolin <- downloadHandler(
    filename = function() { paste("UserDefined_ViolinPlot", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = categviolinprint(), width = 16, height = 8)
    }
  )
  
  
  output$corruserplot <- renderPlot({
    
    req(input$file4)
    req(input$numconduser1)
    req(input$numconduser2)
    
    con1<-as.character(unique(norm_sum_final()$condition)[input$numconduser1])
    con2<-as.character(unique(norm_sum_final()$condition)[input$numconduser2])
    
    foldchangevalue<-as.character(paste("log2FC_",con2,"-",con1,sep=""))
    qvaluename<-paste("q.val_",con2,"-",con1,sep="")
    sig_name <- as.character(paste("sig_",con2,"-",con1,sep=""))
    

    ggplot()+
      # geom_hline(yintercept = 0,linetype="dashed")+
      #geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(x=.data[[con2]] - .data[[con1]],y=uservalue),size=0.25,alpha=0.2,color="grey")+
      geom_point(data=userinput()%>%filter(UserAnnotation != "__"), aes(x=.data[[con2]] - .data[[con1]],y=uservalue,color=UserAnnotation),size=1,alpha=0.5)+
      theme_classic()+scale_color_viridis_d(end=0.8)+
      labs(y="User-defined Values", x=paste("Log2FC(",con2,"-",con1,")",sep=""),title="User Correlation Plot")+
      theme(axis.title = element_text(size=16),axis.text = element_text(size = 14),plot.title = element_text(size=20,hjust=0.5), legend.title = element_blank())
      
    
  })
  
  dfuser2 <- reactive({
    req(input$file6)
    inFile <- input$file6
    if (is.null(inFile))
      return(NULL)
    tbl <- read.csv(inFile$datapath)
    
    return(tbl)
  })
  
  
  localizationvolcanomaker <-reactive({
    
    if (input$txthumanmouse == "mouse") {
      ann_combine<-dfuser2()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol")
    }
    
    if (input$txthumanmouse == "human") {
      ann_combine<-dfuser2()%>%
        # dplyr::mutate(Gene.Symbol = casefold(Gene.Symbol ,upper=TRUE))%>%
        dplyr::left_join(.,locale(),by="Gene.Symbol")%>%
        dplyr::mutate(Compartment = dplyr::if_else(is.na(Compartment), "__",Compartment))
    }
    return(ann_combine)
  })
  
  
  output$plotvolcanomaker<- renderPlot({
    req(input$file6)
    req(input$size10)
    req(input$size11)
    req(input$alpha10)
    req(input$alpha11)
    
    print(head(dfuser2()))
    
    userinput3<-dfuser2()%>%
      dplyr::mutate(sig = dplyr::if_else(abs(Log2FC)>log2(input$num11),dplyr::if_else(p.value<input$num10,dplyr::if_else(Log2FC>0,"up","down"),"n.s."),"n.s."))
    print(head(userinput3))
  ### need to finish the plotting aes
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    p1
  },width=1000,height=650)
  
  
  plotvolcanomakerprint<- reactive({
    req(input$file6)
    req(input$size10)
    req(input$size11)
    req(input$alpha10)
    req(input$alpha11)
    
    print(head(dfuser2()))
    
    userinput3<-dfuser2()%>%
      dplyr::mutate(sig = dplyr::if_else(abs(Log2FC)>log2(input$num11),dplyr::if_else(p.value<input$num10,dplyr::if_else(Log2FC>0,"up","down"),"n.s."),"n.s."))
    print(head(userinput3))
    ### need to finish the plotting aes
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, -log10(p.value),color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "RawValue") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, -log10(p.value)),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%dplyr::filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=sig),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        geom_text_repel(data=userinput3%>%filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color=labelplot),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Significance" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=sig),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=sig),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Label" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=userinput3%>%filter(labelplot == "no"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha10,size=input$size10)+
        geom_point(data=userinput3%>%filter(labelplot == "yes"), aes(Log2FC, p.value,color=labelplot),alpha=input$alpha11,size=input$size11)+
        # geom_text_repel(data=userinput3%>%filter(label == "yes")%>%
        #                   mutate(label=if_else(label == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color=sig),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+scale_color_viridis_d(end=0.8)+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())+scale_color_viridis_d(end=0.8)
    }
    
    
    if (  input$txtlabeluser2 == "Yes"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        geom_text_repel(data=localizationvolcanomaker()%>%dplyr::filter(labelplot == "yes")%>%
                          mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, p.value,label=label,color="red"),
                        max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    if (  input$txtlabeluser2 == "No"  & input$txtcolorlabelonly == "Localization" & input$txtpvaluelog == "NegLogarithm") {
      p1<-ggplot()+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment != input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha10,size=input$size10,color="black")+
        geom_point(data=localizationvolcanomaker()%>%dplyr::filter(Compartment == input$localizationuser), aes(Log2FC, p.value),alpha=input$alpha11,size=input$size11,color="red")+
        # geom_text_repel(data=userinput3%>%dplyr::filter(labelplot == "yes")%>%
        #                   mutate(label=dplyr::if_else(labelplot == "yes",Gene.Symbol,"")),aes(Log2FC, -log10(p.value),label=label,color="red"),
        #                 max.overlaps = Inf,show.legend=FALSE,box.padding = 1,size=5)+
        theme_classic()+
        geom_vline(xintercept = log2(input$num11))+geom_vline(xintercept = -log2(input$num11))+
        geom_hline(yintercept = -log10(input$num10))+ 
        labs(x=paste("log2FC(",input$textxaxis2,"-",input$textxaxis1,")",sep=""), y=paste("-log10(",input$textyaxis,"-value)",sep=""))+
        theme(axis.title = element_text(size=20),axis.text = element_text(size = 18),legend.title = element_blank())
    }
    
    return(p1)
  })
  
  output$printplotvolcanomaker <- downloadHandler(
    filename = function() { paste("UserDefined_Volcano", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, plot = plotvolcanomakerprint(), width = 12, height = 8)
    }
  )
  
  
}


# Run the application 
shinyApp(ui = ui, server = server)