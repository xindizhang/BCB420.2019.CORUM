# Read data
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
tmp[1,]
tmp$Synonyms[1:9]
tmp$`Cell line`[1:9]
tmp$`subunits(UniProt IDs)`[1:9]
tmp$`subunits(Entrez IDs)`[1:9]
tmp$`Protein complex purification method`[1:9]
tmp$`GO ID`[1:9]
tmp$`GO description`[1:5]
tmp$`FunCat ID`[1:5]
tmp$`subunits(Gene name)`[1:5]
tmp$`subunits(Gene name syn)`[1:5]
# filter human data
human_data <- tmp[tmp$Organism == 'Human', ]
human_data$Organism == 'Human'
nrow(human_data)   # 2916
sum(is.na(genes)) 
sum(genes == "") #0
sum(genes == "N/A") #0


# Check outdated data
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL)) 
head(HGNC)

genes <- c()
for (i in seq_len(nrow(human_data))){
  geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
  genes <- c(geneName, genes)
}
sel <- ( ! (genes %in% HGNC$sym)) 
length(genes[ sel ] )  # 230
length( unique(genes[ sel ])) # 81


# for (i in seq_len(nrow(human_data))){
#   print(i)
#   geneName <- strsplit(human_data$`subunits(Gene name)`[i], ";")
#   print(geneName)
#   for (j in length(geneName[[1]])){
#     if ((geneName[[1]][j] == "") & (geneName[[1]][j] == "N/A")){
#       geneName <- geneName[[1]][-j]
#     }
#   }
#   human_data$`subunits(Gene name)`[i] <- geneName
# }


# Check what genes are outdated
for (gene in unique(genes[ sel ]) ){
    iPrev <- grep(gene, HGNC$prev)[1]
    if (length(iPrev) == 1){
      print(HGNC$sym[iPrev])
    }
}

count <- 0
for (i in seq_len(nrow(human_data))){
  geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
  for (gene in geneName){
    # print(gene)
    iPrev <- grep(gene, HGNC$prev)[1]
    if ((gene != "") & (length(iPrev) == 1) & (! is.na(iPrev))) {
        # print(iPrev)
        # print(human_data$`subunits(Gene name)`[i])
        newgene <- gsub(gene, HGNC$sym[iPrev], human_data$`subunits(Gene name)`[i])
        human_data$`subunits(Gene name)`[i] <- newgene
        count <- count + 1
        # print(human_data$`subunits(Gene name)`[i])
      
    }
  }
}

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
length(finalgenes[ sel ] )  # 230
length( unique(finalgenes[ sel ])) # 81


# finalgene <- c()
# for (i in seq_len(nrow(human_data))){
#   for (j in length(human_data$`subunits(Gene name)`[i][[1]])){
#     finalgene <- c(human_data$`subunits(Gene name)`[i][[1]][j], finalgene)
#   }
# }
# length(finalgene)
# 
# for (i in seq_len(nrow(unkSym))) {
#   iPrev <- grep(unkSym$unk[i], HGNC$prev)[1] # take No. 1 if there are several
#   if (length(iPrev) == 1) {
#     unkSym$new[i] <- HGNC$sym[iPrev]
#   } else {
#     iSynonym <- which(grep(unkSym$unk[i], HGNC$synonym))[1]
#     if (length(iSynonym) == 1) {
#       unkSym$new[i] <- HGNC$sym[iSynonym]
#     }
#   }
# }
# 
# 
# 
# unlist(strsplit(geneNameEg, ";"))
# for (i in seq_len(nrow(human_data))){
#   geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
#   for (gene in geneName){
#     print(gene)
#     iPrev <- grep(gene, HGNC$prev)[1]
#     if (length(iPrev) == 1){
#       HGNC$sym[iPrev]
#     }
#   }
# }
#    

unique(human_data$`Protein complex purification method`)

df <- data.frame(symbols=character(), complexes=list())
de <- list(hello="hi", goodbye=list("bye", "fsd"))
df = rbind(df,de, stringsAsFactors=FALSE)

# Create a list containing a vector, a matrix and a list.
list_data <- list(c("Jan","Feb","Mar"), matrix(c(3,9,5,1,-2,8), nrow = 2),
                  list("green",12.3))

# Give names to the elements in the list.
names(list_data) <- c("1st Quarter", "A_Matrix", "A Inner list")

# Show the list.
print(list_data)
names(list_data)[1] <- "Hi"
 


elements <- 1
symbolList <- list()
for (i in seq_len(nrow(human_data))){
  geneName <- unlist(strsplit(human_data$`subunits(Gene name)`[i], ";"))
  for (gene in geneName){
    if (gene != ""){
      if (gene %in% names(symbolList)){
        sel <- which(names(symbolList) == gene)
        print(gene)
        print(symbolList[sel])
        symbolList[[sel]] <- c(symbolList[[sel]], human_data$ComplexName[i])
      }else{
        symbolList[elements] <- c(human_data$ComplexName[i])
        names(symbolList)[elements] <- gene
        elements <- elements + 1
      }
    }
  }
}
save(symbolList, file = file.path("inst", "extdata", "symbolToCom.RData"))
load(file = file.path("inst", "extdata", "symbolToCom.RData"))
symbolList["CDKN1A"]



