
# ---------------------------- maxquant ------------------------------
library(iq)
library(dplyr)
library(data.table)


#' @description
#' make protein list for maxLFQ, topN, meanInt and meadianPolish.
#' @param data.psm psm file in pd output
#' @param normalization bool. Whether do median normalization.
makeProteinListPD <- function(data.psm, normalization = T) {
  # >>> 合并pd中precursor abundance >>> 
  # 原则：去重总和
  psm <- data.psm[!data.psm$Master.Protein.Accessions == "", ]
  psm <- psm %>% 
    group_by(Annotated.Sequence, Charge, Spectrum.File) %>% 
    distinct(Precursor.Abundance, .keep_all = T) %>% 
    mutate(sum.precursor.abundance = sum(Precursor.Abundance, na.rm = T)) %>% 
    ungroup()
  
  psm <- psm %>% 
    select(Annotated.Sequence, Master.Protein.Accessions, Charge, Spectrum.File, sum.precursor.abundance) %>% 
    distinct(Annotated.Sequence, Master.Protein.Accessions, Charge, Spectrum.File, sum.precursor.abundance)
  psm$sum.precursor.abundance[psm$sum.precursor.abundance == 0] <- NA_real_
  
  ex <- unique(psm$Spectrum.File)
  if (normalization) {
    ex_median <- rep(NA, length(ex))
    names(ex_median) <- ex
    for (i in ex) {
      print(i)
      tmp <- subset(psm, Spectrum.File == i)
      ex_median[i] <- median(tmp$sum.precursor.abundance, na.rm = TRUE)
    }
    f <- mean(ex_median)/ex_median
    psm$sum.precursor.abundance.norm <- psm$sum.precursor.abundance * f[psm$Spectrum.File]
  }
  
  # create a protein list
  protein.list <- list()
  
  protein.group <- unique(psm$Master.Protein.Accessions)
  for (i in 1:length(protein.group)) {
    protein.group.select <- protein.group[i]
    
    a <- psm[psm$Master.Protein.Accessions == protein.group.select, ]
    b <- data.frame(cn = a$Spectrum.File, 
                    rn = paste(a$Annotated.Sequence, a$Charge, sep = "_"), 
                    quant = a$sum.precursor.abundance)
    b <- b[!is.na(b$quant),]
    
    m <- matrix(NA, nrow = length(unique(b$rn)), ncol = length(ex), 
                dimnames = list(unique(b$rn), ex))
    
    if (nrow(b) > 0) {
      for (j in 1:nrow(b)) {
        rn <- as.character(b$rn[j])
        cn <- as.character(b$cn[j])
        if (is.na(m[rn, cn])) {
          m[rn, cn] <- b$quant[j]
        } else {
          m[rn, cn] <- m[rn, cn] + b$quant[j]
        }
      }
      protein.list[[protein.group.select]] <- log2(m)
    } else {
      protein.list[[protein.group.select]] <- NA
    }
  }
  return(protein.list)
}

library(dplyr)
library(tidyr)
library(openxlsx2)
library(stringr)
library(shiny)
convertModifications <- function(mods) {
  if (mods == "") {
    return(mods)
  }
  # 提取所有修饰信息
  mod_list <- strsplit(mods, "; ")[[1]]
  
  # 使用正则表达式匹配位置和修饰类型
  mod_positions <- gsub("\\(([^)]+)\\)", "", mod_list)
  mod_types <- gsub(".*\\(([^)]+)\\).*", "\\1", mod_list)
  
  # 合并相同类型的修饰
  unique_mod_types <- unique(mod_types)
  result <- sapply(unique_mod_types, function(type) {
    positions <- mod_positions[mod_types == type]
    paste(length(positions), "x", type, " [", paste(positions, collapse = "; "), "]", sep = "")
  })
  
  return(paste(result, collapse = "; "))
}

processModification <- function(mod, constant) {
  # 使用正则表达式提取中括号中的内容
  bracket_content <- gsub(".*\\[(.*)\\].*", "\\1", mod)
  # 提取数字部分，并转换为整数
  numbers <- as.integer(unlist(regmatches(bracket_content, gregexpr("\\d+", bracket_content))))
  # 对每个数字加上常数
  if (length(numbers) == 0) {
    output_string <- mod
  } else {
    new_numbers <- numbers + constant
    modification <- unlist(regmatches(bracket_content, gregexpr("[A-Za-z]", bracket_content)))
    # 将新数字合并为字符串，并保留字母部分
    new_bracket_content <- paste(paste0(modification, new_numbers), collapse = ", ")
    # 替换原字符串中的中括号内容
    output_string <- gsub("\\[.*\\]", paste0("[", new_bracket_content, "]"), mod)
  }
  
  return(output_string)
}

addStartNumber <- function(input_string, constant) {
  modifications <- unlist(strsplit(input_string, "; "))
  # 对每个修饰项应用内部函数
  new_modifications <- sapply(modifications, processModification, constant = constant)
  # 将处理后的修饰项重新合并为一个字符串
  output_string <- paste(new_modifications, collapse = "; ")
  return(output_string)
}


makeResultListPD <- function(protein.table, data.protein.group, data.peptide.group, data.msms){
  # 定义表头
  protein.header1 <- c("First Protein ID", "IDs of protein", "Protein Descriptions", "Gene names", "Score", 
                       "Coverage (%)")
  protein.header2 <- c("peptides", "Unique peptides", "MS/MS count", 
                       "Mol.weight [kDa]")
  protein.table.header <- colnames(protein.table)[-1]
  peptide.header <- c("Peptide sequence", "Position", "Charges", "RT (min)", "Intensity", 
                      "Best MS/MS scan", "Missed cleavages", "theo.MH+", "PPM", "Modification")
  
  # >>> 初始化返回值 >>>
  result.table <- list(data.frame(t(c(protein.header1, protein.table.header, protein.header2))))
  rows.grouping <- list()
  
  data.msms <- data.msms[!data.msms$Master.Protein.Accessions == "", ]
  data.peptide.psm <- merge(data.peptide.group, 
                            data.msms %>% 
                              rowwise()%>% 
                              mutate(Annotated.Sequence = toupper(Annotated.Sequence), 
                                     Modifications = convertModifications(Modifications)), 
                            by = c("Master.Protein.Accessions", "Annotated.Sequence", "Modifications", "X..Missed.Cleavages", "Theo..MH...Da."))  # 可以删除Peptide group中一些错误识别的肽段
  
  
  if (length(protein.table.header) > 1) {
    pad.columns <- setNames(replicate(length(protein.table.header)-1, "", simplify = FALSE),
                            paste0("X", 12: (10 + length(protein.table.header))))
  } else {
    pad.columns <- NULL
  }
  
  row.n <- 2
  protein.group.names <- unique(data.msms$Master.Protein.Accessions)
  round.n <- length(protein.group.names)
  withProgress(message = "Extract info", {
    for (i in 1: round.n) {
      row.group.start <- row.n
      
      protein.group.name <- protein.group.names[i]
      first.protein.names <- strsplit(protein.group.name, split = "; ")[[1]]
      
      protein.info <- data.frame()
      j = 0
      while (nrow(protein.info) == 0) {
        j <- j+1
        if (j > length(first.protein.names)) {
          break
        }
        protein.info <- data.protein.group %>% 
          filter(Accession == first.protein.names[j])
      }
      intensity <- protein.table[protein.table$Accession == protein.group.name, ]
      
      if (nrow(protein.info) == 0 | nrow(intensity) == 0) {
        next
      }
      
      peptide.psm.tmp <- data.peptide.psm %>% 
        filter(Master.Protein.Accessions == protein.group.name)
      if (nrow(peptide.psm.tmp) == 0) {
        next
      }
      
      protein.group.tmp <- protein.info %>% 
        summarise(Accession = Accession, 
                  Accessions = paste(unique(unlist(strsplit(unique(peptide.psm.tmp$Protein.Accessions), "; "))), collapse = "; "),  #  合并所有protein accession
                  Description = Description, 
                  Genename = str_extract(Description, "(?<=GN=)\\w+"), 
                  Score = Score.Sequest.HT..Sequest.HT, 
                  Coverage = Coverage...., 
                  intensity[-1], 
                  Peptides = X..Peptides, 
                  Unique.Peptides = X..Unique.Peptides, 
                  MSMS.count = X..PSMs, 
                  Weight = MW..kDa.)
      
      peptide.tmp <- peptide.psm.tmp %>% 
        group_by(Annotated.Sequence, Modifications) %>% 
        reframe(Sequence = gsub("\\[.*?\\]|\\.", "", unique(Annotated.Sequence)), 
                Position = gsub(".*\\[(.*?)\\]", "\\1", strsplit(unique(Positions.in.Master.Proteins), "; ")[[1]][1]), 
                Charges = paste(unique(Charge), collapse = "; "), 
                RT = RT..min.[XCorr == max(XCorr)][1], 
                Intensity = Precursor.Abundance[XCorr == max(XCorr)], 
                Best.MSMS = paste(Spectrum.File[XCorr == max(XCorr)][1], First.Scan[XCorr == max(XCorr)][1], sep = "-"), 
                Missed.cleavage = unique(X..Missed.Cleavages), 
                Weight = unique(Theo..MH...Da.), 
                Ppm = DeltaM..ppm.[XCorr == max(XCorr)[1]], 
                Modification = addStartNumber(unique(Modifications), as.integer(strsplit(Position, "-")[[1]][1]) - 1)) %>% 
        select(-c(1, 2))
      
      protein.group.length <- length(c(protein.header1, protein.table.header, protein.header2))
      colnames(protein.group.tmp) <- paste0("X", 1: protein.group.length)
      row.n <- row.n + 1L
      rows.grouping[[i]] <- c(row.n, 0)
      
      # 添加肽段表头
      row.n <- row.n + 1
      
      # 添加肽段数据
      peptide.tmp <- peptide.tmp %>% 
        mutate(X1 = "", !!!pad.columns) %>% 
        select(X1, everything())
      colnames(peptide.tmp) <- paste0("X", 1: protein.group.length)
      
      
      row.n <- row.n + nrow(peptide.tmp)
      rows.grouping[[i]][2] <- row.n
      result.table <- c(result.table, list(protein.group.tmp, 
                                            data.table(t(c("", peptide.header, rep("", protein.group.length - 11)))), 
                                            peptide.tmp))
      incProgress(1 / round.n)
    }
    result.table <- rbindlist(result.table, use.names = F)
  })
  return(list(result.table = result.table, rows.grouping = rows.grouping))
}




# pd.psm <- readData("../LibrarySearchingCompare/data/new/dda/pd/MS20240511-HELA-DDA-160180200ng_PSMs.txt")
# pd.protein <- readData("../LibrarySearchingCompare/data/new/dda/pd/MS20240511-HELA-DDA-160180200ng_Proteins.txt")
# pd.peptide.group <- readData("../LibrarySearchingCompare/data/new/dda/pd/MS20240511-HELA-DDA-160180200ng_PeptideGroups.txt")
# protein.table <- makeProteinListPD(pd.psm)
# protein.table <- createProteinTable(protein.table, method = "MaxLFQ")
# result <- system.time({
#   res <- createResultTablePD(protein.table, pd.protein, pd.peptide.group, pd.psm)
#   groupRows(res$result.table, res$rows.grouping, "./pd_test.xlsx")
# })
# 
# result
