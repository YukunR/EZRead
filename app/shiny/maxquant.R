
# ---------------------------- maxquant ------------------------------
library(iq)
library(dplyr)
library(data.table)
library(tibble)

#' @description
#' extractLFQ is for extract lfq line in maxquant protein group file
#' @param data.protein protein group file in maxqaunt output
extractLFQ <- function(data.protein.group) {
  data.protein.group <- subset(data.protein.group, Reverse == "")    # remove reversed entries
  lfq <- grep("^LFQ", colnames(data.protein.group))
  lfq <- data.protein.group[, lfq] %>% 
    mutate(Accession = data.protein.group$Protein.IDs) %>% 
    select(c(Accession, everything()))
  return(lfq)
}


#' @description
#' make protein list for topN, meanInt and meadianPolish.
#' @param data.evidence evidence file
#' @param data.protein.group protein group file
#' @param normalization bool. Whether do median normalization.
makeProteinListMaxquant <- function(data.protein.group, data.evidence, normalization = T) {
  rownames(data.protein.group) <- data.protein.group$Protein.IDs
  data.protein.group <- subset(data.protein.group, Reverse == "")
  rownames(data.evidence) <- data.evidence$id
  data.evidence$Experiment <- as.character(data.evidence$Experiment)
  ex <- unique(data.evidence$Experiment)
  if (normalization) {
    ex_median <- rep(NA, length(ex))
    names(ex_median) <- ex
    for (i in ex) {
      tmp <- subset(data.evidence, Experiment == i)
      ex_median[i] <- median(tmp$Intensity, na.rm = TRUE)
    }
    f <- mean(ex_median) / ex_median
    data.evidence$Intensity <- data.evidence$Intensity * f[data.evidence$Experiment]
  }
  
  # create a protein list
  protein.list <- list()
  
  round.n <- nrow(data.protein.group)
  suppressWarnings({
    withProgress(message = "Make protein list", {
      for (i in 1: round.n) {
        tmp <- unlist(strsplit(as.character(data.protein.group$Evidence.IDs[i]), ";"))
        a <- data.evidence[tmp, ]
        b <- data.frame(cn = a$Experiment,
                        rn = paste(a$Modified.sequence, a$Charge, sep = "_"),
                        quant = a$Intensity) %>% 
          filter(!is.na(quant))
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
          protein.list[[rownames(data.protein.group)[i]]] <- log2(m)
        } else {
          protein.list[[rownames(data.protein.group)[i]]] <- NA
        }
        incProgress(1 / round.n)
      }

    })
  })
  return(protein.list)
}


# >>> format converting >>>
library(openxlsx2)
library(tidyr)

#' @description
#' create output table
#' @param protein.table intensity table, make by createProteinTable()
#' @param data.protein.group protein group file
#' @param data.peptide.group peptide group file
#' @param data.msms msms file
#' @param data.modification modification file
makeResultListMaxquant <- function(protein.table, data.protein.group, data.peptide.group, data.msms, data.modification) {
  # >>> 定义表头 >>>
  protein.header1 <- c("First Protein ID", "IDs of protein", "Potein Descriptions", "Gene names", "Score", 
                       "Coverage (%)")
  protein.header2 <- c("Peptides", "Unique peptides", "MS/MS count", 
                       "Mol.weight [kDa]", "Sequence length")
  protein.table.header <- colnames(protein.table)[-1]
  peptide.header <- c("Peptide sequence", "Position", "Charges", "RT (min)", "Intensity", 
                      "MS/MS count", "Best MS/MS scan", "Missed cleavages", "theo.MH+", "PPM", "Modification")
  # >>> 初始化返回值 >>>
  result.table <- list(data.frame(t(c(protein.header1, protein.table.header, protein.header2))))
  rows.grouping <- list()
  
  # >>> 数据预处理 >>>
  data.protein.group.select <- data.protein.group %>% 
    select(c(Protein.IDs, Protein.names, Gene.names, Score, 
             Sequence.coverage...., Peptides, Unique.peptides, 
             MS.MS.count, Mol..weight..kDa., Sequence.length, Peptide.IDs)) %>% 
    mutate(First.protein.ID = gsub(";.*$", "", Protein.IDs)) %>% 
    merge(protein.table, by.x = "Protein.IDs", by.y = "Accession") %>% 
    select(c(First.protein.ID, 1: 5, colnames(protein.table[-1]), everything()))
  
  row.names(data.peptide.group) <- data.peptide.group$id
  data.peptide.group.select <- data.peptide.group %>% 
    select(c(Sequence, Start.position, End.position, Charges, grep("Intensity\\.", colnames(data.peptide.group)), 
             MS.MS.Count, Best.MS.MS, Missed.cleavages, Mass,
             Oxidation..M..site.IDs)) %>% 
    mutate(Oxidation..M..site.IDs = ifelse(Oxidation..M..site.IDs == "", NA, Oxidation..M..site.IDs))
  
  row.names(data.modification) <- data.modification$id
  data.modification <- select(data.modification, c(Position.in.peptide, Amino.acid))
  
  row.names(data.msms) <- data.msms$id
  data.msms <- select(data.msms, c(Raw.file, Scan.number, Oxidation..M..site.IDs, Mass.error..ppm., Retention.time)) %>% 
    mutate(scan.number.combine = paste(as.character(Raw.file), as.character(Scan.number), sep = "-"))
  
  if (length(protein.table.header) > 1) {
    pad.columns <- setNames(replicate(length(protein.table.header)-1, "", simplify = FALSE),
                            paste0("X", 13: (11 + length(protein.table.header))))
  } else {
    pad.columns <- NULL
  }
  
  row.n <- 2
  
  # >>> 向表格中填入信息 >>>
  round.n <- nrow(data.protein.group.select)
  suppressWarnings({
    withProgress(message = "Extract info", {
      for (i in 1: round.n) {
        # >> 写Protein group信息 >>
        protein.group.tmp <- data.protein.group.select[i, ]
        
        # >> 提取peptide信息 >>
        peptides.id.tmp <- unlist(strsplit(as.character(protein.group.tmp$Peptide.IDs), ";"))
        peptides.tmp <- data.peptide.group.select[peptides.id.tmp, ] %>% 
          mutate(position = paste(as.character(Start.position), as.character(End.position), sep = "-"), 
                 modification = NA) %>% 
          select(Sequence, position, everything())
        data.msms.tmp <- data.msms[peptides.tmp$Best.MS.MS, ]
        peptides.tmp$Best.MS.MS <- data.msms.tmp$scan.number.combine
        peptides.tmp$ppm <- data.msms.tmp$Mass.error..ppm.
        peptides.tmp$rt <- data.msms.tmp$Retention.time
        
        # >> 提取ox modification >>
        if (sum(!is.na(peptides.tmp$Oxidation..M..site.IDs)) > 0) {
          for (j in 1: nrow(peptides.tmp)) {
            if (!is.na(peptides.tmp$Oxidation..M..site.IDs[j])) {
              modification.id.tmp <- unlist(strsplit(as.character(peptides.tmp$Oxidation..M..site.IDs[j]), ";"))
              modification.tmp <- data.modification[modification.id.tmp, ] %>% 
                mutate(modification = paste0(Amino.acid, peptides.tmp$Start.position[j] + Position.in.peptide - 1))
              modification.add <- paste0("Oxidation", "[", paste(modification.tmp$modification, collapse = "; "), "]")
            } else {
              modification.add <- NA
            }
            peptides.tmp$modification[j] <- modification.add
          }
        }
        
        # >>> write into result table >>>
        protein.group.length <- length(c(protein.header1, protein.table.header, protein.header2))
        protein.group.write <- data.protein.group.select[i, ][1: protein.group.length]
        colnames(protein.group.write) <- paste0("X", 1: protein.group.length)
        row.n <- row.n + 1L
        rows.grouping[[i]] <- c(row.n, 0)
        row.n <- row.n + 1L
        
        peptides.tmp <- select(peptides.tmp, -c(Start.position, End.position)) %>% 
          select(1: 3, rt, 4: 8, ppm, modification)
        peptides.tmp <- peptides.tmp %>% 
          mutate(X1 = "", !!!pad.columns) %>% 
          select(X1, everything())
        colnames(peptides.tmp) <- paste0("X", 1: protein.group.length)
        row.n <- row.n + nrow(peptides.tmp)
        rows.grouping[[i]][2] <- row.n
        result.table <- c(result.table, list(protein.group.write, 
                                             data.table(t(c("", peptide.header, rep("", protein.group.length - 12)))), 
                                             peptides.tmp))
        incProgress(1 / round.n)
      }
      result.table <- rbindlist(result.table, use.names = F)
    })
  })
  
  return(list(result.table = result.table, rows.grouping = rows.grouping))
}


# data.modification <- readData("../LibrarySearchingCompare/data/new/dda/maxquant/Oxidation (M)Sites.txt")
# data.msms <- readData("../LibrarySearchingCompare/data/new/dda/maxquant/msms.txt")
# data.protein.group <- readData("../LibrarySearchingCompare/data/new/dda/maxquant/proteinGroups.txt")
# data.peptide.group <- readData("../LibrarySearchingCompare/data/new/dda/maxquant/peptides.txt")
# protein.table <- extractLFQ(data.protein.group)
# result <- system.time({
#   res <- createResultTableMaxquant(protein.table, data.protein.group, data.peptide.group, data.msms, data.modification)
#   groupRows(res$result.table, res$rows.grouping, "./mq_test.xlsx")
# })
# 
# result
