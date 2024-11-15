library(KEGGREST, quietly = TRUE)
library(tidyverse, quietly = TRUE)

# 返回信息很长，只取基因symbol.根据自己需要调整
symbolOnly <- function(x){
  items <- strsplit(x, ";", fixed = TRUE) %>% unlist()
  return(items[1])
}

# keggGet(x)[[1]]$GENE 数据基因名是个向量，其中奇数位置是 entrezgene_id 偶数位置是 symbol 
geneEntrez <- function(x){
  Sys.sleep(1)
  geneList <- keggGet(x)[[1]]$GENE
  if(!is.null(geneList)){
    listLength <- length(geneList)
    entrezList <- geneList[seq.int(from = 1, by = 2, length.out = listLength/2)]
    entrez <- stringr::str_c(entrezList, collapse = ",")
    return(entrez)
  }else{
    return(NA)
  }
  
}

# keggGet(x)[[1]]$GENE 数据基因名是个向量，其中奇数位置是 entrezgene_id 偶数位置是 symbol 
geneSymbol <- function(x){
  Sys.sleep(1)
  geneList <- keggGet(x)[[1]]$GENE
  if(!is.null(geneList)){
    listLength <- length(geneList)
    symbolList <- geneList[seq.int(from = 2, by = 2, length.out = listLength/2)] %>% map_chr(symbolOnly)
    symbol <- stringr::str_c(symbolList, collapse = ",")
    return(symbol)
  }else{
    return("")
  }
  
}

# 取得 hsaxxxxx 这种通路ID
pathwayID <- function(x){
  items <- strsplit(x, ":", fixed = TRUE) %>% unlist()
  return(items[2])
}

# 建议从这里开始读脚本。建议自己在交互模式下试一下用到的KEGGREST函数，看看返回数据的结构。
# 这是第一步，取得所有的KEGG通路列表
hsaList <- keggList("pathway", "hsa")
IDList <- names(hsaList)

# 将通路ID和通路名放在一个表格(tibble)里
hsaPathway <- tibble::tibble(pathway_id=IDList, pathway_name=hsaList)
head(hsaPathway, n=3) %>% print()

# 用前面定义函数，获得每个通路的基因，然后也放在表格里
pathwayFull <- hsaPathway %>% 
  dplyr::mutate(
    entrezgene_id=map_chr(pathway_id, geneEntrez), 
    hgnc_symbol=map_chr(pathway_id, geneSymbol)
  )

# 保存数据
write_tsv(pathwayFull, path="文章修订-20241114-01/data/KEGGREST.tsv")
dim(pathwayFull) %>% print()

# 会有通路没有基因，我的话只需要有基因的，所以把无基因的移除
pathwayWithGene <- dplyr::filter(pathwayFull, !is.na(entrezgene_id) & hgnc_symbol != "")
write_tsv(pathwayWithGene, path="文章修订-20241114-01/data/KEGGREST_WithGene.tsv")
dim(pathwayWithGene) %>% print()

# GMT
pathwayWithGene$pathway_name <- lapply(
  X = pathwayWithGene$pathway_name, function(x) gsub("\\s-\\sHomo\\ssapiens\\s\\(human\\)", "", x)
) %>% unlist()

keggmt <- "kegg_20241114.gmt"
apply(pathwayWithGene, 1, FUN = function(vec){
  sink(keggmt, append = TRUE)
  cat(paste(c(vec[2], vec[1], strsplit(vec[4], split = ",")[[1]]), collapse = "\t"))
  cat("\n")
  sink()
})


















