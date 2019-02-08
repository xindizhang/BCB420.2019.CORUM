# `BCB420.2019.CORUM`

#### (CORUM data annotatation of human genes)

&nbsp;

###### [Xindi Zhang], Bioinformatics, University of Toronto, Canada. &lt;xindi.zhang@utoronto.ca&gt;

----

**If any of this information is ambiguous, inaccurate, outdated, or incomplete, please check the [most recent version](https://github.com/xindizhang/BCB420.2019.CORUM.git) of the package on GitHub and if the problem has not already been addressed, please [file an issue](https://github.com/xindizhang/BCB420.2019.CORUM.git/issues).**

----

## 1 About this package:
This pakage is designed using (hyginn/BCB420.2019.STRING) and (https://github.com/hyginn/rpt) as an template. The 2019.STRING package is written by Dr. Boris Steipe.

This package is designed to download network data from the [CORUM database](https://mips.helmholtz-muenchen.de/corum/), map annotated human protein complex data to [HGNC](https://www.genenames.org/), and provide examples for statistic computation of the databases.

The package serves dual duty, as an RStudio project, as well as an R package that can be installed. Package checks **pass without errors, warnings, or notes**.

&nbsp;

#### In this project ...

```text
 --BCB420.2019.CORUM/
   |__.gitignore
   |__.Rbuildignore
   |__BCB420.2019.CORUM.Rproj
   |__DESCRIPTION
   |__dev/
      |__toBrowser.R               # display .md files in your browser
   |__inst/
      |__extdata/
         |__symbolToCom.RData      # Tool for finding compleses encoded by a gene
   |__LICENSE
   |__NAMESPACE
   |__R/
      |__zzz.R                     # Welcoming file 
      |__script.R                  # A raw script for this work flow
   |__README.md                    # this file

```

&nbsp;

----

## 2 CORUM Data

CORUM is a database of annotated mammalian protein complexes. The main mammalian organisms are human(67%), mouse(15%), and rat(10%). There are data for 4274 protein complexes, which are encoded by 4473 genes (Giurgiu _et al._ 2019). All information is obtained from individual experiment. Each complex also has information crossed referenced from other datasets. All CORUM data is available under a CC-BY 4.0 license.




&nbsp;

#### 2.1 Data semantics

Available annotations and information are:

1. **Complex name and ID**
2. **UniPort annotation**: UniPort Id, Identification of subunits
3. **Gene Ontology (GO) information**: This is the functional annotation of the complex
4. **Protein complex purification method**: The experimental method used to purify the protein complexes
5. **Organism**
6. **PMID**
7. **Comment**: Information about disease and cellular function information

Annotation of splice variant complexes is included for functional properties since some complex isoform is tissue specific. This can be a great implication of the protein functions. 

In addition, a CORUM tool using IntAct database to predict protein- protein interactions within protein complexes and a visualization tool (Cytoscape written in javascript) of the interactions can be fouond at https://www.ebi.ac.uk/intact/. This tool is proved for reliable predictions(Giurgiu _et al._ 2019). 

CORUM integrates and cross-references data of mammalian protein complexes from experiment, Uniport, GO and other literature resouces. By demonstrating the importance of protein complex, it provide information for bioinformatics and biomedical research (especially cancer research). To build a larger datasets, input of other researchers are greatly welcoming. Please contact CORUM at andreas.ruepp@helmholtz-muenchen.de (Giurgiu _et al._ 2019).

&nbsp;

## 3 Data download and cleanup

To download the source data from CORUM:

1. Navigate to the [**CORUM** database](https://mips.helmholtz-muenchen.de/corum/#) and follow the link to the [download section](https://mips.helmholtz-muenchen.de/corum/#download).
2. Choose the type of complex data you wnat to download
3. This workflow uses complete complexes data. Download the following data file: 

* `allComplexes.txt.zip` (752kb)	complete complexes data;

4. Uncompress the file and place it into a sister directory of your working directory which is called `data`. (It should be reachable with `file.path("..", "data")`). 

&nbsp;

## 4 Mapping ENSEMBL IDs to HGNC symbols
CORUM database do have HGNC symbols for genes that encodes for each complex. However, these genes might be outdated, so the following workflow will update the outdated gene symbols in the datasets

&nbsp;

#### Preparations

To begin, we need to load the HGNC reference data:

```R
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame

```
&nbsp;

#### 4.1 Step one: Fiter for human data

&nbsp;

```R
  # Read data and check the information we have
  getwd()
  tmp <- readr::read_delim(file.path("../data", "allComplexes.txt"),
                           delim = "\t")
  nrow(tmp) # 4274
  # check datasets
  colnames(tmp)
  # [1] "ComplexID"                          
  # [2] "ComplexName"                        
  # [3] "Organism"                           
  # [4] "Synonyms"                           
  # [5] "Cell line"                          
  # [6] "subunits(UniProt IDs)"              
  # [7] "subunits(Entrez IDs)"               
  # [8] "Protein complex purification method"
  # [9] "GO ID"                              
  # [10] "GO description"                     
  # [11] "FunCat ID"                          
  # [12] "FunCat description"                 
  # [13] "subunits(Gene name)"                
  # [14] "Subunits comment"                   
  # [15] "PubMed ID"                          
  # [16] "Complex comment"                    
  # [17] "Disease comment"                    
  # [18] "SWISSPROT organism"                 
  # [19] "subunits(Gene name syn)"            
  # [20] "subunits(Protein name)"
  # filter human data
  human_data <- tmp[tmp$Organism == 'Human', ]
  human_data$Organism == 'Human'
  # Get basic information of the dataset 
  nrow(human_data)   # 2916
  # Check for missing variables
  sum(is.na(genes)) 
  sum(genes == "")  #0
  sum(genes == "N/A") #0

```

&nbsp;

#### 4.2  Step two: update the outdated symbols

The gene symbols encodes for one protein complex is quoted in one string and seperated by ; or ;;. To get a more accurate check of the missing data, I obtained a list of genes:

```R
  genes <- c()
  for (i in seq_len(nrow(human_data))){
    geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
    genes <- c(geneName, genes)
  }
  sel <- ( ! (genes %in% HGNC$sym)) 
  # we found that there are 81 unique genes that are outdated
  length(genes[ sel ] )  # 230
  length( unique(genes[ sel ])) # 81
  
  # Check which genes are outdated
  for (gene in unique(genes[ sel ]) ){
      iPrev <- grep(gene, HGNC$prev)[1]
      if (length(iPrev) == 1){
        print(HGNC$sym[iPrev])
      }
  }
  
  # Replace the oudated genes with the new version
  # These code uses Dr. Boris's BCB420.2019.STRING as a reference
  for (i in seq_len(nrow(human_data))){
    geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
    for (gene in geneName){
      iPrev <- grep(gene, HGNC$prev)[1]
      if ((gene != "") & (length(iPrev) == 1) & (! is.na(iPrev))) {
          newgene <- gsub(gene, HGNC$sym[iPrev], human_data$`subunits(Gene name)`[i])
          human_data$`subunits(Gene name)`[i] <- newgene
          count <- count + 1
      }
    }
  }
  ```
  
  &nbsp;
  
  
  #### 4.5 Final validation
  
  Validation: The validation shows no change which means that something is wrong. I need to check this and modify the code.
  
  ```R
  
    finalgenes <- c()
    for (i in seq_len(nrow(human_data))){
      geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
      finalgenes <- c(geneName, finalgenes)
    }
    loc <- which(finalgenes == "")
    finalgenes <- finalgenes[-loc]
    sum(finalgene == "N/A")
    sum(! is.na(finalgene)) * 100 / length(finalgene) 
    sel <- ( ! (finalgenes %in% HGNC$sym))
    # The performance is really bad, I will have to check it again
    length(finalgenes[ sel ] )  # 230
    length( unique(finalgenes[ sel ])) # 81


```

&nbsp;

# 5 Annotating gene sets with CORUM and inAct Data
Gereation of RData file which comtains a list with symbols and the protein complex the symbol encodes for. 
```R
elements <- 1
symbolList <- list()
for (i in seq_len(nrow(human_data))){
  geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
  for (gene in geneName){
    if (gene != ""){
      if (gene %in% names(symbolList)){
        sel <- which(names(symbolList) == gene)
        symbolList[[sel]] <- c(symbolList[[sel]], human_data$ComplexName[i])
      }else{
        symbolList[elements] <- c(human_data$ComplexName[i])
        names(symbolList)[elements] <- gene
        elements <- elements + 1
      }
    }
  }
}

# Save data
save(symbolList, file = file.path("inst", "extdata", "symbolToCom.RData"))

# Load data
load(file = file.path("inst", "extdata", "symbolToCom.RData"))

# To check what complex a gene encode for, use symbolList[gene symbol]
symbolList["CDKN1A"]
# [1] "CDKN1A"
# $`CDKN1A`
# [1] "Cell cycle kinase complex CDC2"
# [2] "Cell cycle kinase complex CDK2"
# [3] "Cell cycle kinase complex CDK4"

```

## 6 References

&nbsp;


* Giurgiu, M., Reinhard, J., Brauner, B., Dunger-Kaltenbach, I., Fobo, G., Frishman, G., Montrone, C., & Ruepp, A. (2019). CORUM: the comprehensive resource of mammalian protein complexes-2019. Nucleic acids research, D1, D559-D563.

&nbsp;

## 7 Acknowledgements

Thanks to Simon Kågedal's very useful [PubMed to APA reference tool](http://helgo.net/simon/pubmed/).

Thank you professor Boris Steipe for providing the R package templates, project templates and useful resources and tools for us. 


&nbsp;

&nbsp;

<!-- [END] -->
