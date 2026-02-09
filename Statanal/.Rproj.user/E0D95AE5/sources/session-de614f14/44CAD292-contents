# R functions


read_all_data <- function(dir) {
  files <- list.files(dir, pattern = "(plate)_[0-9]_summary.txt", full.names = TRUE)
  all_df <- lapply(files, read.table, header = TRUE, fill = F, sep = "\t",skip = 0)
  return(all_df)
}

get_param <- function(input, groupss, param, dix = T) {
  
  out <- data.frame()
  for (i in groupss) {
    tmp <- input %>% subset(Group == i ) %>% dplyr::select(Plate_ID, Isolate, param, Well)
    colnames(tmp) <- c(tmp[1,1], "Isolate", "gr")
    out <- rbind(out, tmp)
  }
  
  return(out)
}

permutation.test <- function(treatment, outcome, original, n){
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), median))
  }
  
  result=(1+sum(abs(distribution) >= abs(original)))/(n)
  return(list(result, distribution))
  
}

meandiff <- function(x,i) {
  E <- x[i,]
  return(mean(E$YGR067C_ko)-mean(E$BY4741))
}

nicefy_gene_results <- function(input, genelist) {
  temp_table <- input
  genes_in_input <- subset(genelist, (genelist$Gene_stable_ID %in% input$Gene_stable_ID))
  
  out_table <- merge(temp_table, genes_in_input, by.x = c("Gene_stable_ID"), by.y = c("Gene_stable_ID"), all = TRUE) %>% 
    relocate(Gene_name, .after = Gene_stable_ID)
  
  return(out_table)
}

megatest <- function(reslist, GS, experiments = 1, seed = 819230) {
  set.seed(seed)
  gup <- list()
  gdn <- list()
  for (i in seq(reslist)) {
    g <- FGSEA_wrap(reslist[[i]], GS, lims = c(5,200))
    gup[[i]] <- g$pAdjDistinctDirUp
    gdn[[i]] <- g$pAdjDistinctDirDn
  }
  
  n <- nrow(reslist[[1]])
  out <- list()
  for (i in seq(GS$gsc)) {
    out[[i]] <- data.frame(Method = c(rep("Fisher", experiments), 
                                      rep("Boschloo", experiments), 
                                      rep("FGSEA_up", experiments),
                                      rep("FGSEA_down", experiments)),
                           Experiment = rep(seq(experiments), 4),
                           log = rep(NA,4*experiments),
                           PDS = rep(NA,4*experiments))
    
    genelist <- GS$gsc[i]
    m <- length(genelist[[1]])
    for (j in seq(reslist)) {
      k <- reslist[[j]] %>% subset(padj < 0.05)
      q <- k %>% subset(Gene_stable_ID %in% genelist[[1]]) 
      
      k <- nrow(k)
      q <- nrow(q)
      
      print(paste0("Currently doing ", names(GS$gsc)[i], ". Dataset ", j, " out of ", length(reslist), "."))
      fish <- fisher.test(matrix(c(q, k-q, m-q, n-m-k+q), nr = 2), alternative = "g")$p.value
      bosch <- btest(k-q , n-m, q , m, tsmethod = "minlike")$p.value
      # fish <- runif(1)
      # bosch <- runif(1)
      
      r <- j  %% experiments + experiments * !(j %% experiments)
      c <-  floor((j-1)/experiments)+3
      
      out[[i]][r+experiments*0, c] <- fish
      out[[i]][r+experiments*1, c] <- bosch
      out[[i]][r+experiments*2, c] <- gup[[j]][i]
      out[[i]][r+experiments*3, c] <- gdn[[j]][i]
    }
  }
  names(out) <- names(GS$gsc)
  return(out)
}

FGSEA_wrap <- function(input, gscGO, lims = c(10, 200)) {
  tmp <- input %>% dplyr::filter(!is.na(padj) & !is.na(pvalue) & baseMean > 0)
  Pval <- tmp$pvalue
  Pval[Pval==0] <- min(Pval[Pval != 0])
  names(Pval) <- tmp$Gene_stable_ID
  FC <- tmp$log2FoldChange
  names(FC) <- tmp$Gene_stable_ID
  gsaRes <- runGSA(-log10(Pval)*sign(FC), geneSetStat="fgsea", gsc=gscGO, gsSizeLim=lims, nPerm=10000)
  return(gsaRes)
}

volcano_RNA <- function(input, plim = 0.05, flim = 0.58, qval = F) {
  
  # Let the first condition (often reference strain or log-phase) be A and second condition (often mutant strain or PDS-phase) be B
  # Given a linear model Y = K * x + m, where Y is the predicted peak area, x an indicator function [0 if A, 1 if B], m the intercept (i.e. predicted Y for A)
  # and K a vector of the slopes caused by B. Thus Y_A = K * 0 + m -> Y_A = m, and Y_B = K * 1 + m -> Y_B = K + m
  # Then the log2fold change (l2fc) is: log2(Y_B/Y_A) = log2((K+m)/m) = log2(K/m + 1) 
  tmp <- data.frame(feats = input[,1], pvals = input[,5+qval], l2fc = input[,3]) %>% mutate(diff = rep("NO"))
  colnames(tmp) <- c("feats", "pvals", "l2fc", "diff")
  
  # tmp <- tmp %>% mutate(diff = if_else((pvals <= plim) & (l2fc >= flim), "UP", "NO")) %>% mutate(diff = if_else((pvals <= plim) & (l2fc <= -flim), "DOWN", "NO"))
  tmp <- tmp %>% mutate(diff = case_when((pvals <= plim) & (l2fc >= flim) ~ "UP", (pvals <= plim) & (l2fc <= -flim) ~ "DOWN", .default = diff))
  
  groupcol <- c(UP = "#F8766D", NO = "grey70", DOWN = "#00B6EB")
  
  if (qval) {
    ytit <- "-log10(FDR)"
  } else {
    ytit <- "-log10(p-values)"
  }
  
  # pmax <- max(-log10(tmp$pvals))
  # fmax <- max(abs(tmp$l2fc))
  # 
  # print(pmax)
  # print(fmax)
  
  pmax <- 15
  fmax <- 7
  
  p <- ggplot(tmp, aes(x = l2fc, y = -log10(pvals), col = diff)) + geom_point() + geom_vline(xintercept = c(-flim, flim), linetype = "dashed") + geom_hline(yintercept = -log10(plim), linetype = "dashed") + xlim(c(-fmax, fmax)) + ylim(0, pmax) + scale_color_manual(values = groupcol) + theme_bw() + xlab("log2 Fold Change") + ylab(ytit) + theme(legend.position = "none")
  return(p)
}

myenrich <- function(input, p_thresh = 0.001) {
  tmp <- input %>% subset(lowest_p < p_thresh)
  
  tmp <- tmp %>% mutate(ups = str_count(Up_regulated, ",") + sign(str_count(Up_regulated, ",")), 
                        dns = str_count(Down_regulated, ",") + sign(str_count(Down_regulated, ","))) %>% mutate(tot = ups + dns) %>% mutate(ratio = (ups - dns)/tot)
  
  
  tmp <- tmp %>% mutate(class = if_else(ratio >= 0.6, "Upregulated", if_else(ratio <= -0.6, "Downregulated", "Mixed")))
  
  out <- ggplot(tmp, aes(x = Fold_Enrichment, y = reorder(Term_Description, -lowest_p), fill = -log10(lowest_p), size = tot)) + geom_point(pch = 21) +
    scale_fill_gradient2(
      mid = "white",
      high = "#A58AFF",
      midpoint = -log10(p_thresh),
      limits = c(-log10(p_thresh), 16)
    ) + xlim(0, max(tmp$Fold_Enrichment)+5) +
    labs(x = "Fold Enrichment", y = "", fill = expression(-log[10](p)), size = "# genes")
  
  out <- out + facet_grid(factor(class, levels=c("Upregulated", "Mixed", "Downregulated"))~. , scale = "free", space = "free") + theme_bw()
  return(out)
}

auto_scale <- function(expr) {
  m <- mean(expr)
  sd <- sd(expr)
  out <- (expr-m)/sd
  return(out)
}

volcano_MS <- function(input,plim = 0.05, flim = 0.58, qval = F) {
  
  # Let the first condition (often reference strain or log-phase) be A and second condition (often mutant strain or PDS-phase) be B
  # Given a linear model Y = K * x + m, where Y is the predicted peak area, x an indicator function [0 if A, 1 if B], m the intercept (i.e. predicted Y for A)
  # and K a vector of the slopes caused by B. Thus Y_A = K * 0 + m -> Y_A = m, and Y_B = K * 1 + m -> Y_B = K + m
  # Then the log2fold change (l2fc) is: log2(Y_B/Y_A) = log2((K+m)/m) = log2(K/m + 1) 
  tmp <- data.frame(feats = input[,1], pvals = input[,14+qval], l2fc = log2(input[,9]/input[,2]+1)) %>% mutate(diff = rep("Not significant"))
  
  # tmp <- tmp %>% mutate(diff = if_else((pvals <= plim) & (l2fc >= flim), "UP", "NO")) %>% mutate(diff = if_else((pvals <= plim) & (l2fc <= -flim), "DOWN", "NO"))
  tmp <- tmp %>% mutate(diff = case_when((pvals <= plim) & (l2fc >= flim) ~ "Upregulated", (pvals <= plim) & (l2fc <= -flim) ~ "Downregulated", .default = diff))
  
  colnames(tmp) <- c("feats", "pvals", "l2fc", "diff")
  groupcol <- c("Upregulated" = "#F8766D", "Not significant" = "grey70", "Downregulated" = "#00B6EB")
  
  if (qval) {
    ytit <- "-log10(FDR)"
  } else {
    ytit <- "-log10(p-values)"
  }
  
  # pmax <- max(-log10(tmp$pvals))
  # fmax <- max(abs(tmp$l2fc))
  # 
  # print(pmax)
  # print(fmax)
  
  pmax <- 4.5
  fmax <- 14
  
  p <- ggplot(tmp, aes(x = l2fc, y = -log10(pvals), col = diff)) + geom_point() + geom_vline(xintercept = c(-flim, flim), linetype = "dashed") + geom_hline(yintercept = -log10(plim), linetype = "dashed") + xlim(c(-fmax, fmax)) + ylim(0, pmax) + scale_color_manual(values = groupcol) + theme_bw() + xlab("log2 Fold Change") + ylab(ytit) + theme(legend.title = element_blank())
  return(p)
}

get_unique_metab_names <- function(input) {
  tmp <- input %>% subset(`Metabolite name`!="Unknown")%>% dplyr::select(`Metabolite name`, `m/z similarity`, l2fc, p)
  tmp <- tmp %>% mutate(`Metabolite name` = str_split_fixed(str_split_fixed(`Metabolite name`, ";",2)[,1], ": ",2)[,2])
  tmp <- tmp %>% arrange(`Metabolite name`)
  out <- data.frame(`Metabolite name` = character(), `m/z similarity` = character(), l2fc = double(), p = double())
  # print(typeof(tmp[2,2]))
  i <- 1
  while(i < nrow(tmp)-1) {
    if (tmp[i, 1] == tmp[i+1, 1]) {
      for (j in seq(i+1,(nrow(tmp)-1))) {
        if (tmp[i, 1] != tmp[j, 1]) {
          tmp_dup <- tmp[i:(j-1),] %>% arrange(desc(`m/z similarity`))
          if (tmp_dup[1,2] != tmp_dup[2,2]) {
            out <- rbind(out, tmp_dup[1,])
          } else {
            tmp_dup2 <- tmp_dup %>% subset(`m/z similarity` == max(`m/z similarity`))
            if (length(unique(sign(tmp_dup2$l2fc))) == 1) {
              out <- rbind(out, tmp_dup[1,])
            } else {
              print(paste0(tmp_dup[1,1], " was discarded due to uncertainty." ))
            }
          }
          i <- j
          break
          # print(tmp_dup)
          # print("")
        }
      }
    } else {
      out <- rbind(out, tmp[i,])
      i = i + 1
    }
  }
  return(out)
}

get_MAFS <- function(input, ref) {
  empteh <- rep("null", nrow(input))
  MAF <- data.frame(
    database_identifier = empteh,
    chemical_formula = empteh,
    smiles = empteh,
    inchi = empteh,
    metabolite_identification = empteh,
    mass_to_charge = empteh,
    fragmentation = empteh,
    modifications = empteh,
    charge = empteh,
    retention_time = empteh,
    taxid = empteh,
    species = empteh,
    database = empteh,
    database_version = empteh,
    reliability = empteh,
    uri = empteh,
    search_engine = empteh,
    search_engine_score = empteh,
    smallmolecule_abundance_sub = empteh,
    smallmolecule_abundance_stdev_sub = empteh,
    smallmolecule_abundance_std_error_sub = empteh,
    QC1	 = empteh,
    QC2 = empteh,
    R1_pre = empteh,
    R2_pre = empteh,
    R3_pre = empteh,
    Y1_pre = empteh,
    Y2_pre = empteh,
    Y3_pre = empteh,
    QC3 = empteh,
    R1_post = empteh,
    R2_post = empteh,
    R3_post = empteh,
    Y1_post = empteh,
    Y2_post = empteh,
    Y3_post = empteh,
    QC4 = empteh
  )
  
  MAF$chemical_formula <-  input$Formula
  MAF$smiles <- input$SMILES
  MAF$inchi <- input$INCHIKEY
  MAF$metabolite_identification <- input %>% dplyr::select(`Metabolite name`) %>% mutate(`Metabolite name` = if_else(`Metabolite name` == "Unknown", "null", str_split_fixed(str_split_fixed(`Metabolite name`, ";",2)[,1], ": ",2)[,2])) %>% pull()
  MAF$mass_to_charge <- input$`Average Mz`
  MAF$retention_time <-  input$`Average Rt(min)`
  MAF$taxid <- rep("NCBI:txid559292", nrow(MAF))
  MAF$species <- rep("Saccharomyces cerevisiae", nrow(MAF))
  
  MAF$QC1 <- input[,36]
  MAF$QC2 <- input[,37]
  MAF$R1_pre <- input[,38]
  MAF$R2_pre <- input[,39]
  MAF$R3_pre <- input[,40]
  MAF$Y1_pre <- input[,41]
  MAF$Y2_pre <- input[,42]
  MAF$Y3_pre <- input[,43]
  MAF$QC3 <- input[,44]
  MAF$R1_post <- input[,45]
  MAF$R2_post <- input[,46]
  MAF$R3_post <- input[,47]
  MAF$Y1_post <- input[,48]
  MAF$Y2_post <- input[,49]
  MAF$Y3_post <- input[,50]
  MAF$QC4 <- input[,51]
  
  MAF <- MAF %>% 
    left_join(ref %>% dplyr::select(SMILES, ChEBI), by = c("smiles" = "SMILES"), relationship = "many-to-many") %>% 
    mutate(database_identifier = paste0("CHEBI:",ChEBI)) %>% dplyr::select(-ChEBI) %>% 
    mutate(database_identifier = if_else(database_identifier == "CHEBI:NA", "null", database_identifier))
  
  MAF[MAF == "null"] <- ""
  MAF <- MAF %>% subset(metabolite_identification !="")
  return(MAF)
}

get_final_list <- function(input, refer) {
  tmp <- data.frame(Feature_ID = input[,1], l2fc = log2(input[,9]/input[,2]+1), p = input[,14], FDR = input[,15])
  tmp <- merge(refer$feature_data, tmp)
  posn <- nrow(refer$feature_data %>% subset(Split == "RP_POS"))
  tmp_pos <- tmp %>% subset(Split == "RP_POS")
  tmp_neg <- tmp %>% subset(Split == "RP_NEG") %>% mutate(Alignment_ID = Alignment_ID - posn)
  out <- list(tmp_pos, tmp_neg)
  names(out) <- c("POS", "NEG")
  return(out)
}