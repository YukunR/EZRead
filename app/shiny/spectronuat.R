
# ---------------------------- spectronuat ------------------------------
library(iq)
library(dplyr)
library(data.table)
library(stringr)
library(openxlsx2)


#' @description
#' make protein list for maxLFQ, topN, meanInt and meadianPolish.
#' @param data.psm psm file in pd output
#' @param intensity.col character, one of PEP.MS1Quantity, PEP.MS2Quantity, ...
#' @param normalization bool. Whether do median normalization.
makeProteinListSP <- function(data.report, intensity.col, normalization = T) {
  # >>> quality control >>> 
  data.report <- data.report %>% 
    filter(!is.na(PG.Qvalue) & PG.Qvalue < 0.01 & 
             !is.na(EG.Qvalue) & EG.Qvalue < 0.01)
  norm.data <- preprocess(data.report, 
                          sample_id  = "R.FileName", 
                          secondary_id = c("EG.PrecursorId", "FG.Charge"), 
                          intensity_col = intensity.col, 
                          pdf_out = NULL, 
                          median_normalization = normalization)
  protein.list <- create_protein_list(norm.data)
  return(protein.list)
}

extract_modifications_and_positions <- function(start, sequence) {
  mods <- str_extract_all(sequence, "\\[.+?\\]")[[1]]  # 提取所有修饰
  mods_positions <- vector("list", length(mods))
  
  # 逐一处理每个修饰，删除其在序列中的第一次出现，以避免重复定位
  for (i in seq_along(mods)) {
    mod <- mods[i]
    pos <- str_locate(sequence, fixed(mod))[1]
    mods_positions[[i]] <- paste(mod, "at position", ifelse(pos == 1, start, start + pos - 2))
    sequence <- str_replace(sequence, fixed(mod), "")
  }
  # 转化为单一字符串，用分号分隔
  paste(unlist(mods_positions), collapse = "; ")
}


makeResultListSP <- function(protein.table, data.report) {
  protein.table <- as.data.table(protein.table)
  data.report <- as.data.table(data.report)
  protein.header1 <- c("First Protein ID", "IDs of protein", "Potein Descriptions", "Gene names", "Score", 
                       "Coverage (%)")
  protein.header2 <- c("Peptides", "Unique peptides", "MS/MS count", 
                       "Mol.weight [kDa]")
  protein.table.header <- colnames(protein.table)[-1]
  peptide.header <- c("Peptide sequence", "Position", "Charges", "RT (min)", "Intensity", 
                      "Best MS/MS scan", "Missed cleavages", "theo.MH+", "PPM", "Modification")
  
  # >>> 初始化返回值 >>>
  result.table <- list(data.frame(t(c(protein.header1, protein.table.header, protein.header2))))
  rows.grouping <- list()
  
  if (length(protein.table.header) > 1) {
    pad.columns <- setNames(replicate(length(protein.table.header)-1, "", simplify = FALSE),
                            paste0("X", 12: (10 + length(protein.table.header))))
  } else {
    pad.columns <- NULL
  }
  
  row.n <- 2
  round.n <- nrow(protein.table)
  suppressWarnings({
    withProgress(message = "Extract info", {
      for (i in 1: round.n) {
        # >> 计算pep信息
        protein.group.name <- protein.table[i, Accession]
        protein.group.tmp <- data.report[PG.ProteinGroups == protein.group.name, ]
        
        # 选择需要的列并去重
        peptide.tmp <- protein.group.tmp[, {
          first.pos <- as.numeric(strsplit(as.character(PEP.PeptidePosition), ";")[[1]][1])
          Positions <- paste(first.pos, first.pos + nchar(as.character(PEP.StrippedSequence)) - 1, sep = "-")
          
          # 直接在数据表中计算 Modifications
          Modifications <- extract_modifications_and_positions(first.pos, substr(EG.ModifiedSequence, 2, nchar(EG.ModifiedSequence)))
          
          # 计算其他列
          Intensity <- sum(FG.MS1Quantity)
          Charge <- paste(unique(FG.Charge), collapse = ";")
          bestIndex <- which.max(FG.MS1Quantity)
          
          Missedcleavages <- sum(ifelse(gregexpr("[KR](?!P|$)", PEP.StrippedSequence, perl = TRUE)[[1]] > 0,
                                        attr(gregexpr("[KR](?!P|$)", PEP.StrippedSequence, perl = TRUE)[[1]], "match.length"), 0))
          
          .(PEP.StrippedSequence = PEP.StrippedSequence[1],
            Positions = Positions,
            Charge = Charge,
            rt = EG.MeanApexRT[bestIndex],
            Intensity = Intensity,
            BestMSMSScan = paste(R.FileName[bestIndex], FG.MS2ApexScanIndex[bestIndex]),
            Missedcleavages = as.double(Missedcleavages),
            mass = max(as.numeric(FG.PrecMz) * as.numeric(FG.Charge)),
            ppm = "-",
            Modifications = paste(unique(Filter(function(x) nchar(str_trim(x)) > 0, unlist(Modifications))), collapse = "; ")
          )
        }, by = EG.ModifiedSequence] %>% unique()
        
        unique.peptide.tmp <- unique(protein.group.tmp[, .(FG.Id, PEP.IsProteinGroupSpecific)])
        # 处理 protein.group.tmp
        protein.group.write <- protein.group.tmp[1, ][
          , FirstProteinID := sapply(strsplit(as.character(PG.ProteinAccessions), ";"), `[`, 1)
        ][
          , Mol.weight := as.numeric(sapply(strsplit(as.character(PG.MolecularWeight), ";"), `[`, 1))
        ][
          , peptides := length(unique(unique.peptide.tmp$FG.Id))
        ][
          , unique.peptides := sum(unique.peptide.tmp$PEP.IsProteinGroupSpecific)
        ]
        
        # 选择和合并列
        protein.group.write <- cbind(protein.group.write, protein.table[i, ][, -1, with = FALSE])[
          , .(FirstProteinID, 
              PG.ProteinAccessions, 
              PG.ProteinDescriptions, 
              PG.Genes, 
              PG.Qvalue, 
              PG.Coverage, 
              .SD, 
              peptides, 
              unique.peptides, 
              R.Data.Points.per.Peak..MS2..EXT., 
              Mol.weight), 
          .SDcols = colnames(protein.table)[-1]
        ]
        
        peptide.tmp <- peptide.tmp[, !("EG.ModifiedSequence"), with = FALSE]
        
        protein.group.length <- length(c(protein.header1, protein.table.header, protein.header2))
        colnames(protein.group.write) <- paste0("X", 1: protein.group.length)
        row.n <- row.n + 1L
        rows.grouping[[i]] <- c(row.n, 0)
        row.n <- row.n + 1L
        peptide.tmp <- peptide.tmp %>% 
          mutate(X1 = "", !!!pad.columns) %>% 
          select(X1, everything()) 
        colnames(peptide.tmp) <- paste0("X", 1: protein.group.length)
        row.n <- row.n + nrow(peptide.tmp)
        rows.grouping[[i]][2] <- row.n
        result.table <- c(result.table, list(protein.group.write, 
                                             data.table(t(c("", peptide.header, rep("", protein.group.length - 11)))), 
                                             peptide.tmp))
        incProgress(1 / round.n)
      }
      result.table <- rbindlist(result.table, use.names = F)
    })
  })
  
  return(list(result.table = result.table, rows.grouping = rows.grouping))
}

# pep.long <- readData("../LibrarySearchingCompare/data/new/dia/20240621_122137_MS20240511-HELA-DIA-160ng_Report.tsv")
# protein.list.sp <- makeProteinListSP(pep.long, intensity.col = "FG.MS2Quantity")
# protein.table.sp <- createProteinTable(protein.list.sp, method = "MaxLFQ")
# 
# result1 <- system.time({
#   res <- createResultTableSP(protein.table.sp, pep.long)
# })
# 
# result1
# result2 <- system.time({
#   groupRows(res$result.table, res$rows.grouping, "./sp_test.xlsx")
# })
# result2
