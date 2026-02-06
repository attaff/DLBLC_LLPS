library(limma)
library(dplyr)
library(ggplot2)
library(ggtext)
library(ggrepel)


########### 1. LLPS-related genelist construction ###########
GSE56315 <- read.delim("GSE56315_series_matrix.txt",skip = 90)
rownames(GSE56315) <- GSE56315[,1]
GSE56315 <- GSE56315[,-1]
GSE56315 <- na.omit(GSE56315)

group <- factor(c(rep("Tumor", 55), rep("Normal", 33)),
                levels = c("Normal", "Tumor")) 

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(Tumor - Normal, levels = design)
fit <- lmFit(GSE56315, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)

results <- topTable(fit2, 
                    number = Inf, 
                    adjust.method = "fdr",
                    coef = 1) %>%
  tibble::rownames_to_column("ProbeID")

sig_genes <- results %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC)))

probe_gpl570 <- read.delim("GPL570-55999.txt",skip = 16)

sig_genes$symbol <- probe_gpl570$Gene.Symbol[match(sig_genes$ProbeID,probe_gpl570$ID)]
sig_genes[sig_genes == "" | sig_genes == " "] <- NA
sig_genes <- na.omit(sig_genes)

results <- results %>%
  mutate(trend = case_when(
    logFC > 1 & adj.P.Val < 0.05 ~ "up",
    logFC < -1 & adj.P.Val < 0.05 ~ "down",
    TRUE ~ "non-sig"
  ))

library(ggrastr)
p <- ggplot(results,
            aes(x = logFC, y = -log10(adj.P.Val), colour = trend)) +
            ggrastr::geom_point_rast(alpha = 0.5, size = 1, raster.dpi = 600) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("up"     = "#C6295C",
                                 "down"   = "#2C6DB2",
                                 "non-sig"= "grey60")) +
  labs(x = "log<sub>2</sub>(Fold Change)",
       y = "-log<sub>10</sub>(adj.P.Val)",
       colour = "Trend") +
  theme_bw(14) + 
  theme(axis.text  = element_text(colour = "black"),
        axis.title = element_markdown(),
        panel.grid.minor = element_blank())
ggsave("Volcano_raster.pdf", p, width = 6, height = 6, device = cairo_pdf)


library(VennDiagram)
library(grid)
library(openxlsx)
deg_genes <- sig_genes$symbol
llpsmarker <- read.delim("LLPS_genelist_DrLLPS.txt")
llpsmarker <- llpsmarker[llpsmarker$Species == "Homo sapiens",]
llpsmarker[llpsmarker == ""] <- NA
llpsmarker <- na.omit(llpsmarker)#3611
gene_list <- list(
  LLPS_marker = llpsmarker$Gene.name,
  DEG = deg_genes
)
p <- venn.diagram(
  x = gene_list,
  scaled = FALSE,
  category.names = c("LLPS markers", "DEG"),
  fill = c("#6e48fb", "#ffa500"),
  cat.col = c("#6e48fb", "#ffa500"),
  alpha = 0.5,
  col = "black",
  lwd = 1.3,
  lty = 1,
  cex = 1.5,
  cat.cex = 1.3,
  cat.default.pos = "outer",
  cat.dist = 0.15,
  cat.pos = c(-90, 90),
  output = FALSE,
  filename = NULL
)
grid.newpage()
grid.draw(p)

llps_deg <- intersect(sig_genes$symbol,llpsmarker$Gene.name)  



########### 2. Data preparation for SUSUCC GSE11318 GSE10846 GSE53786 GSE56315  ###########
#————————————————————————————————————————————————SYSUCC cohort—————————————————————————————————————————————#
#library(openxlsx)
#clinical <- read.xlsx("baseline.xlsx")
#rownames(clinical) <- clinical[,1]
#expr <- read.xlsx("data_SYSUCC.xlsx",rowNames = T)
#expr <- expr[,colnames(expr) %in% clinical$ID]
#expr <- as.data.frame(t(expr))
#probe_gpl570 <- read.delim("GPL570-55999.txt",skip = 16)
#colnames(expr) <- probe_gpl570$Gene.Symbol[match(colnames(expr),probe_gpl570$ID)]
#expr <- expr[clinical$ID,]
#full_matrix <- cbind(clinical,expr)
#full_matrix <- full_matrix[,-1]
#saveRDS(full_matrix,"sysucc.RData")

#————————————————————————————————————————————————GSE10846 cohort—————————————————————————————————————————————#
# clinical imformation
file_path <- "GSE10846_series_matrix.txt"
all_lines <- readLines(file_path, n = 100)
selected_lines <- all_lines[30:83]
temp_file <- tempfile()
writeLines(selected_lines, temp_file)
data <- read.delim(temp_file, header = FALSE, sep = "\t")
unlink(temp_file)

clinical_data <- data[c(2,10,11,16,17,18,19,20,21,22),]
clinical_data <- clinical_data[,-1]
new_row_names <- c("sample_name", "sex","age","diagnosis","OS","OS.time","regimen","ECOG","stage","LDHratio")

extract_content <- function(x) {
  sub(".*:\\s*", "", x)
}
clinical_data[] <- lapply(clinical_data, function(column) {
  sapply(column, extract_content)
})
tdata <- as.data.frame(t(clinical_data))
tdata <- tdata %>%
  mutate(OS = ifelse(OS == "ALIVE", 0, 1))
# mRNA data preparation
gse10846 <- read.delim("GSE10846_series_matrix.txt",skip = 83)
probe_gpl570 <- read.delim("GPL570-55999.txt",skip = 16)
probe_gpl570 <- probe_gpl570[,c(1,11,12,15,16)]
expr10846 <- merge(probe_gpl570, gse10846, by.x = "ID", by.y = "ID_REF", all.x = TRUE)
expr10846 <- expr10846[,-c(1,4,5)]
# GSE10846 cohort, LLPS-related genelist#
llps_matrix <- expr10846[expr10846$Gene.Symbol %in% llps_variable,]
llps_matrix <- llps_matrix[!is.na(llps_matrix$Gene.Symbol) & llps_matrix$Gene.Symbol != "", ]
llps_matrix <- as.data.frame(llps_matrix)

df <- llps_matrix[,-1]
df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], as.numeric)
df_summary <- df %>%
  group_by(df$ENTREZ_GENE_ID) %>%
  summarise(across(.cols = 2:ncol(df), ~mean(.x, na.rm = TRUE)), .groups = 'drop')
df_summary <- df_summary %>%
  rename(ENTREZ_GENE_ID = `df$ENTREZ_GENE_ID` )

probe_gpl570 <- probe_gpl570%>%
  distinct(ENTREZ_GENE_ID, .keep_all = TRUE)

df_combined <- df_summary %>%
  left_join(probe_gpl570, by = "ENTREZ_GENE_ID") %>%
  select(Gene.Symbol, everything())
df_combined <- as.data.frame(df_combined)
rownames(df_combined) <- df_combined[,1]
df_combined <- df_combined[,-c(1,2,423:425)]

OS_surv <- tdata[,c(1,5,6)]
rownames(OS_surv) <- tdata[,1]
OS_surv <- OS_surv[,-1]
OS_surv <- as.data.frame(t(OS_surv)) 
OS_row <- OS_surv[1, ]
OS.time_row <- OS_surv[2, ]
merged_OS_row <- c(OS_row[names(df_combined)])
merged_OS.time_row <- c(OS.time_row[names(df_combined)])

llps_matrix <- rbind(merged_OS_row, merged_OS.time_row,df_combined)
rownames(llps_matrix) <- c("OS", "OS.time",rownames(df_combined) )

ML_matrix <- as.data.frame(t(llps_matrix))
ML_matrix[, 2:ncol(ML_matrix)] <- lapply(ML_matrix[, 2:ncol(ML_matrix)], as.numeric)
ML_matrix <- ML_matrix %>%
  filter(!is.na(OS.time)) %>%
  mutate_all(~ifelse(is.na(.), 0, .))
ML_matrix <- ML_matrix %>%
  filter(OS.time != 0)
ML_matrix$OS <- as.integer(ML_matrix$OS)
# saveRDS(ML_matrix,"F:\\LLPS_DLBCL\\gse10846_cohort.rdata")


########### 4. data preprocess ###########
gse_10846 <- readRDS("gse10846_cohort.rdata")
gse_10846 <- as.data.frame(gse_10846)
gse_10846$OS.time <- as.numeric(gse_10846$OS.time)
gse_10846$OS.time <- gse_10846$OS.time*365
gse_10846$OS <- as.integer(gse_10846$OS)
gse_10846[,c(4:ncol(gse_10846))] <- lapply(gse_10846[,c(4:ncol(gse_10846))],as.numeric)
gse_10846 <- na.omit(gse_10846)
gse_10846 <- gse_10846[gse_10846$OS.time> 30,]
##############################################################################
gse_11318 <- readRDS("gse11318_cohort.RData")
gse_11318 <- as.data.frame(gse_11318)
gse_11318$OS.time <- as.numeric(gse_11318$OS.time)
gse_11318$OS.time <- gse_11318$OS.time*365
gse_11318$OS <- as.integer(gse_11318$OS)
gse_11318[,c(4:ncol(gse_11318))] <- lapply(gse_11318[,c(4:ncol(gse_11318))],as.numeric)
gse_11318 <- na.omit(gse_11318)
gse_11318 <-gse_11318[gse_11318$OS.time> 30,]
##############################################################################
gse_53786 <- readRDS("gse53786_cohort.RData")
gse_53786 <- as.data.frame(gse_53786)
gse_53786$OS.time <- as.numeric(gse_53786$OS.time)
gse_53786$OS.time <- gse_53786$OS.time*365
gse_53786$OS <- as.integer(gse_53786$OS)
gse_53786[,c(4:ncol(gse_53786))] <- lapply(gse_53786[,c(4:ncol(gse_53786))],as.numeric)
gse_53786 <- na.omit(gse_53786)
gse_53786 <-gse_53786[gse_53786$OS.time> 30,]
##############################################################################
sysucc <- readRDS("sysucc.RData")
sysucc$ID <- rownames(sysucc)
sysucc <- sysucc[,c(54678,1:54677)]
sysucc <- sysucc[,c(1,3,2,4:ncol(sysucc))]
sysucc <-sysucc[sysucc$OS.time> 30,]
##############################################################################
trainlist=list(train= gse_10846, GSE11318= gse_11318,GSE53786= gse_53786,
               SYSUCC= sysucc)
## removeBatchEffect
library(limma)
do.call(rbind,lapply(trainlist, dim))
intersect <- Reduce(intersect, lapply(trainlist, colnames))

pheno_data <- purrr::map_dfr(
  names(trainlist), 
  function(dataset_name) {
    x <- trainlist[[dataset_name]] %>% as.data.frame()
    selected_cols <- intersect(colnames(x), intersect[1:3])
    x %>% 
      dplyr::select(all_of(selected_cols)) %>% 
      dplyr::mutate(dataset = dataset_name)
  }
)
rownames(pheno_data) <- paste(pheno_data$dataset, pheno_data$ID, sep = ".")

trans_data <- do.call(rbind,
                      lapply(trainlist, function(x) {
                        x <- x %>% as.data.frame() %>%
                          dplyr::select(.,all_of(c(intersect %>% .[4:length(.)])))
                        return(x)
                      })
) %>% as.data.frame()

exp = as.data.frame(trans_data)
ac <- data.frame(
  row.names = rownames(exp), 
  Group = pheno_data$dataset)
batch <- pheno_data$dataset
expr_limma <- removeBatchEffect(t(exp),batch=batch,batch2=NULL,
                                covariates=NULL)
expr_limma <- t(expr_limma) %>% as.data.frame()
expr_final <- pheno_data %>% rownames_to_column("tmp_ID") %>% select(., ID, OS.time, OS, tmp_ID, dataset) %>%
  inner_join(expr_limma %>% as.data.frame() %>%
               rownames_to_column("tmp_ID"), by = "tmp_ID") %>%
  select(., -tmp_ID)

trainlist <- split(expr_final, f = expr_final$dataset)



########### 5. LLPS-related model construction ###########
library(Mime1)

res <- Mime1::ML.Dev.Prog.Sig(
  train_data = trainlist$GSE11318,
  list_train_vali_Data = trainlist,
  unicox.filter.for.candi = T,
  unicox_p_cutoff = 0.001,
  candidate_genes = common_deg,
  mode = 'all',
  nodesize = 5,
  seed = 123
)
#--------------------------------------------------------------------------#
cindex_dis_all(res,validate_set = names(trainlist)[-1],order =names(trainlist),width = 0.35)
all.auc.1y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = trainlist$GSE11318,
                             inputmatrix.list = trainlist,mode = 'all',AUC_time = 1,
                             auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = trainlist$GSE11318,
                             inputmatrix.list = trainlist,mode = 'all',AUC_time = 3,
                             auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res(res.by.ML.Dev.Prog.Sig = res,train_data = trainlist$GSE11318,
                             inputmatrix.list = trainlist,mode = 'all',AUC_time = 5,
                             auc_cal_method="KM")
auc_dis_all(all.auc.5y,
            dataset = names(trainlist),
            validate_set= names(trainlist)[-1],
            order= names(trainlist),
            width = 0.3,
            year= 1)
library(survival)
library(survminer)
candidate_genes <- c("PELP1","DNM1L","TMOD2","C1orf198","PPIF","NR3C1")
exam_source <- trainlist$GSE11318
cox_llps <- exam_source[,colnames(exam_source) %in% c("OS","OS.time",candidate_genes)]
cox_llps <- na.omit(cox_llps)
llps_model <- coxph(Surv(OS.time,OS) ~ ., data = cox_llps)

# saveRDS(llps_model,"LLPS_DLBCL\\llps_model.RData")
coefficients <- as.data.frame(llps_model$coefficients)
fit <- summary(llps_model)
fit
cindex <- fit$concordance
cindex

predicted_scores <- predict(llps_model, type = "risk")
#coefficients <- coef(llps_model)
#coefficients <- as.data.frame(coefficients)

library(pROC)
roc_curve <- roc(cox_llps$OS, predicted_scores)
auc_value <- auc(roc_curve)
auc_value

coef_df <- data.frame(
  gene = rownames(fit$coefficients),
  coefficient = fit$coefficients[, "coef"],
  se = fit$coefficients[, "se(coef)"]
) %>%
  mutate(
    lower = coefficient - 1.96 * se,
    upper = coefficient + 1.96 * se,
    color_sign = coefficient > 0
  )

ggplot(coef_df, aes(x = reorder(gene, coefficient), y = coefficient)) +
  geom_point(aes(color = color_sign), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = color_sign), width = 0.2) +
  scale_color_manual(values = c("dodgerblue", "firebrick1"), guide = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Cox Forest Plot", x = "Genes", y = "Coefficient") +
  theme(plot.title = element_text(hjust = 0.5))

#----------------------------------------------------------------#
df_combined <- trainlist$GSE10846
predicted_survival <- predict(llps_model,newdata = df_combined,type = "lp")  # 使用线性预测值
roc_curves_exam <- roc(df_combined$OS,predicted_survival)
auc(roc_curves_exam)

library(timeROC)
str(df_combined)
df_combined$predicted_risk <- predict(llps_model, newdata = df_combined, type = "lp")

roc_1_year <- timeROC(T = df_combined$OS.time, delta = df_combined$OS, marker = df_combined$predicted_risk, 
                      cause = 1, times = 365, iid = TRUE)
roc_3_years <- timeROC(T = df_combined$OS.time, delta = df_combined$OS, marker = df_combined$predicted_risk, 
                       cause = 1, times = 1095, iid = TRUE)
roc_5_years <- timeROC(T = df_combined$OS.time, delta = df_combined$OS, marker = df_combined$predicted_risk, 
                       cause = 1, times = 1825, iid = TRUE)

par(bty = "l") 
plot(roc_1_year , time = 365, title = " ",col = "firebrick1", xlab = "1 - Specificity", ylab = "Sensitivity")
plot(roc_3_years, time = 1095, add = TRUE, col = "dodgerblue")
plot(roc_5_years, time = 1825, add = TRUE, col = "limegreen")
auc_1_year <- round(roc_1_year$AUC[2], 3)
ci_1_year <- confint(roc_1_year, level = 0.95,n.sim=2000)
ci_1_year_text <- paste(sprintf("%.3f", round(ci_1_year$CI_AUC[,"2.5%"]/100, 4)), "~", sprintf("%.3f", round(ci_1_year$CI_AUC[,"97.5%"]/100, 4)))

auc_3_years <- round(roc_3_years$AUC[2], 3)
ci_3_year <- confint(roc_3_years, parm=NULL, level = 0.95,n.sim=2000)
ci_3_year_text <- paste(sprintf("%.3f", round(ci_3_year$CI_AUC[,"2.5%"]/100, 4)), "~", sprintf("%.3f", round(ci_3_year$CI_AUC[,"97.5%"]/100, 4)))

auc_5_years <- round(roc_5_years$AUC[2], 3) 
ci_5_year <- confint(roc_5_years, parm=NULL, level = 0.95,n.sim=2000)
ci_5_year_text <- paste(sprintf("%.3f", round(ci_5_year$CI_AUC[,"2.5%"]/100, 4)), "~", sprintf("%.3f", round(ci_5_year$CI_AUC[,"97.5%"]/100, 4)))

legend("bottomright", 
       legend = c(paste("1 year\n(AUC =", auc_1_year, ", 95%CI:", ci_1_year_text,")"),
                  paste("3 years\n(AUC =", auc_3_years, ", 95%CI:", ci_3_year_text,")"),
                  paste("5 years\n(AUC =", auc_5_years, ", 95%CI:", ci_5_year_text,")")), 
       col = c("firebrick1", "dodgerblue", "limegreen"), 
       lwd = 2, cex = 0.4, text.font = 7)

my_list <- Map(function(mat, name) {
  first_col <- as.character(mat[, 1])
  new_rownames <- paste(name, first_col, sep = "_")
  rownames(mat) <- new_rownames
  return(mat)
}, trainlist, names(trainlist))

common_cols <- Reduce(intersect, lapply(my_list, colnames))
common_cols
result <- bind_rows(lapply(my_list, function(df) df[, common_cols]))
result$ID <- rownames(result)
saveRDS(result,"all_data.RDS")

#km
risk_matrix <- as.data.frame(trainlist$GSE11318) 
rownames(risk_matrix) <- risk_matrix[,1]
risk_matrix <- risk_matrix[,-1]
risk_matrix <- risk_matrix[,names(risk_matrix) %in% c("OS.time","OS",names(llps_model$coefficients))]

for (i in seq_along(names(llps_model$coefficients))) {
  risk_matrix[, i+2] <- risk_matrix[, i+2] * unname(llps_model$coefficients)[i]
}
risk_matrix <- risk_matrix %>% mutate(risk_score = rowSums(risk_matrix[,3:ncol(risk_matrix)]))

res.cut <- surv_cutpoint(risk_matrix, time = "OS.time", event = "OS",
                         variables = "risk_score")
summary(res.cut)
plot(res.cut, "risk_score", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)

res.cat$OS.time <- res.cat$OS.time / 365
fit <- survfit(Surv(OS.time, OS) ~ risk_score, data = res.cat)
ggsurvplot(fit, data = res.cat, pval = T)
bioCol=c("firebrick1","dodgerblue")
surPlot=ggsurvplot(fit, 
                   data= res.cat,
                   conf.int=F,
                   pval=T,
                   pval.size= 4,
                   
                   legend = "top",
                   font.legend= 10 ,
                   xlab="Time(years)",
                   xlim = c(0,7),
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.3,
                   risk.table.title = "event",
)
surPlot

risk_matrix_53786 <- as.data.frame(trainlist$GSE53786) 
rownames(risk_matrix_53786) <- risk_matrix_53786[,1]
risk_matrix_53786 <- risk_matrix_53786[,-1]
risk_matrix_53786$OS.time <- risk_matrix_53786$OS.time / 365
risk_matrix_53786 <- risk_matrix_53786[,names(risk_matrix_53786) %in% names(risk_matrix)]

for (i in seq_along(names(llps_model$coefficients))) {
  risk_matrix_53786[, i+2] <- risk_matrix_53786[, i+2] * unname(llps_model$coefficients)[i]
}
risk_matrix_53786 <- risk_matrix_53786 %>% mutate(risk_score = rowSums(risk_matrix_53786[,3:ncol(risk_matrix_53786)]))
risk_matrix_53786 <- mutate(risk_matrix_53786, group = ifelse(risk_score > -5.4970 ,"high_risk","low_risk"))

fit <- survfit(Surv(OS.time, OS) ~ risk_matrix_53786$group, data = risk_matrix_53786)
bioCol=c("firebrick1","dodgerblue")
surPlot=ggsurvplot(fit, 
                   data= risk_matrix_53786,
                   conf.int=F,
                   pval=T,
                   pval.size= 4,
                   legend = "top",
                   font.legend= 10 ,
                   xlab="Time(years)",
                   xlim = c(0,7),
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.3,
                   risk.table.title = "event",
)
surPlot

risk_matrix_10846 <- as.data.frame() 
rownames(risk_matrix_10846) <- risk_matrix_10846[,1]
risk_matrix_10846 <- risk_matrix_10846[,-1]
risk_matrix_10846$OS.time <- risk_matrix_10846$OS.time / 365
risk_matrix_10846 <- risk_matrix_10846[,names(risk_matrix_10846) %in% names(risk_matrix)]

for (i in seq_along(names(llps_model$coefficients))) {
  risk_matrix_10846[, i+2] <- risk_matrix_10846[, i+2] * unname(llps_model$coefficients)[i]
}
risk_matrix_10846 <- risk_matrix_10846 %>% mutate(risk_score = rowSums(risk_matrix_10846[,3:ncol(risk_matrix_10846)]))
risk_matrix_10846 <- mutate(risk_matrix_10846, group = ifelse(risk_score > -5.4970 ,"high_risk","low_risk"))

fit <- survfit(Surv(OS.time, OS) ~ risk_matrix_10846$group, data = risk_matrix_10846)
bioCol=c("firebrick1","dodgerblue")
surPlot=ggsurvplot(fit, 
                   data= risk_matrix_10846,
                   conf.int=F,
                   pval=T,
                   pval.size= 4,
                   legend = "top",
                   font.legend= 10 ,
                   xlab="Time(years)",
                   xlim = c(0,7),
                   break.time.by = 1,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.3,
                   risk.table.title = "event",
)
surPlot


matrix <- result
matrix_cluster <- matrix[,names(matrix) %in% candidate_genes]

for (i in seq_along(names(llps_model$coefficients))) {
  matrix_cluster[, i] <- matrix_cluster[, i] * unname(llps_model$coefficients)[i]
}
matrix_cluster <- matrix_cluster %>% mutate(risk_score = rowSums(matrix_cluster[,1:ncol(matrix_cluster)]))
matrix_cluster <- mutate(matrix_cluster, group = ifelse(risk_score > -5.4970 ,"high_risk","low_risk"))
################################################################################################
library(ggplot2)
matrix_cluster$OS <- matrix$OS[match(rownames(matrix_cluster),rownames(matrix))]
matrix_cluster <- matrix_cluster %>%
  mutate(risk_score_adjusted = risk_score + 5.497)

matrix_cluster<-matrix_cluster %>%
  arrange(risk_score_adjusted)


matrix_seperate <- matrix_cluster[grepl("^SYSUCC_", rownames(matrix_cluster)),]

ggplot(matrix_seperate, aes(x = seq_along(risk_score_adjusted), y = risk_score_adjusted, fill = factor(OS))) +
  geom_bar(stat = "identity", width = 0.8,position = "nudge") +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = "Risk score for every patient", y = "Risk Score", fill = "OS Status") +
  theme_minimal() +
  scale_fill_manual(values = c("0" = "#80b1d3", "1" = "#fb8072")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank())



########### 5. immune-filtration ###########
matrix_immune <- as.data.frame(t(matrix[,-c(2:4)]))  
names(matrix_immune) <- matrix_immune[1,]
matrix_immune <- matrix_immune[-1,]
matrix_immune <- matrix_immune %>% mutate_all(as.numeric)

matrix_immune_test <- 2^matrix_immune

library(IOBR)
cibersort <- deconvo_tme(eset = matrix_immune_test, method = "cibersort", arrays = T, perm = 500 )
epic <- deconvo_tme(eset = matrix_immune_test,method = "epic",arrays = T)
quantiseq <- deconvo_tme(eset = matrix_immune_test,method = "quantiseq",
                         arrays = T,tumor = T, scale_mrna = F)
MCPcounter <- deconvo_tme(eset = matrix_immune_test,method = "mcpcounter",tumor = T)
library(xCell)
xCELL <- xCellAnalysis(matrix_immune_test)
xCELL <- as.data.frame(t(xCELL))
names(xCELL) <- paste0(names(xCELL), "_xCELL")
estimate <- deconvo_tme(eset = matrix_immune_test,method = "estimate", tumor = T)

IPS <- deconvo_tme(eset = matrix_immune_test,method = "ips", tumor = T)

indications <- rep("dlbc", ncol(matrix_immune))
Timer <- deconvo_timer(eset = matrix_immune_test, project = NULL, indications = indications)

########################################
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

data_total <- cbind(cibersort[,-c(24:26)],epic[,-1],quantiseq[,-1],MCPcounter[,-1],xCELL, estimate[,-1],IPS[,-1],Timer[,-1])
data_total <- as.data.frame(t(data_total))
colnames(data_total) <- data_total[1,]
data_total <- dplyr::slice(data_total,-1)

data_total <- data_total %>% mutate_all(as.numeric)
data_total <- t(scale(t(data_total)))


highrisk_cols <- colnames(data_total)[matrix_cluster$group == "high_risk"]
lowrisk_cols  <- colnames(data_total)[matrix_cluster$group == "low_risk"]
data_total <- data_total[, c(highrisk_cols, lowrisk_cols)]

matrix_cluster <- matrix_cluster[colnames(data_total),]
grouplist <- matrix_cluster$group
grouplist


mat <- as.matrix(data_total)
pvals <- apply(mat, 1, function(x) {
  wilcox.test(x[group == "high_risk"], x[group == "low_risk"])$p.value
})
# logFC计算
logFC <- rowMeans(mat[, group == "high_risk", drop=FALSE]) -
  rowMeans(mat[, group == "low_risk", drop=FALSE])

get_star <- function(p) {
  if (p < 0.0001) return("****")
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("")
}
stars <- sapply(pvals, get_star)

sig_idx <- which(pvals < 0.05)
mat <- mat[sig_idx, ]
stars <- stars[sig_idx]
logFC <- logFC[sig_idx]

row_labels <- paste0(rownames(mat), " ", stars)
label_colors <- ifelse(logFC > 0, "red", "blue") 

methods <- sub(".*_", "", rownames(mat))
methods <- factor(methods, levels = unique(methods))

method_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(levels(methods))),
  levels(methods)
)

ann_row <- rowAnnotation(
  Method = methods,
  col = list(Method = method_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

risk_colors <- c(high_risk = "#D73027", low_risk = "#4575B4"
col_fun <- colorRamp2(c(-1, 0, 1), c("#336699", "white", "tomato"))
Heatmap(
  mat,
  name = "Infiltration",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = FALSE,
  row_names_side = "right", 
  row_names_rot = 0,
  row_names_gp = gpar(col = label_colors, fontsize = 6), 
  row_labels = row_labels,
  top_annotation = HeatmapAnnotation(Type = grouplist, col = list(Type = risk_colors)),
  column_split = grouplist,
  left_annotation = ann_row, 
  row_split = methods,  
  row_gap = unit(2, "mm"),  
  border = TRUE, 
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8),
  row_names_max_width = unit(10, "cm")
)

##################################################
### immune-related genes
genes <- c(
  "CD276","MS4A1", "CD19", "CD79B","CD160",
  "TNFRSF4","TNFRSF8","TNFRSF9","TNFSF9", "TNFSF4", 
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB", #HLA FAMILY
  "KIR2DL1","KIR2DL3","KIR2DL4","KIR3DL1","KIR3DL3",#KIR FAMILY
  "PDCD1", "CD274","TIGIT", "LAG3", "IDO1","CTLA4","ICOS","HAVCR2", # IMMUNE CHECK POINT FAMILY
  "SIGLEC7", "SIGLEC9", "SIGLEC10", # Siglec family
  "IL10", "IL10RA", "IL10RB","IL6", # Immunosuppressive cytokine and its receptors
  "VEGFA","VEGFB",
  "SIRPA","CD47",
  "CXCL8","CXCR2","CCL2","CCL5","CCR5"
)

setdiff(genes,names(matrix))

ICI_effect <- matrix[,names(matrix) %in% genes]
ICI_effect$Group <- matrix_cluster$group[match(rownames(ICI_effect), rownames(matrix_cluster))]

ICI_effect[, 1:length(genes)] <- lapply(ICI_effect[, 1:length(genes)], as.numeric)
ICI_effect$name <- rownames(ICI_effect)
library(reshape2)

ICI_long <- melt(ICI_effect, id.vars = c("name", "Group"), variable.name = "Gene", value.name = "Expression")
ICI_long$Gene <- factor(ICI_long$Gene, levels = genes) 
p <- ggplot(ICI_long, aes(x = Gene, y = Expression, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), 
               width = 0.6, alpha = 0.1, notch = T, size = 0.5) +
  labs(x = "Gene", y = "Expression", fill = "Risk Group", color = "Risk Group") +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    
    axis.title.x = element_text(size = 12, color = "black", face = "bold", vjust = 1.0),
    axis.title.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 1, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black"),
    plot.caption = element_text(hjust = 0.5, vjust = -1.0, color = "black", size = 9),
    
    legend.title = element_text(size = 11, face = "bold", color = "black"),
    legend.text = element_text(size = 10, color = "black"),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal"
  ) +
  scale_fill_manual(values = c("low_risk" = "#A3CDEA", "high_risk" = "#D9534F")) + 
  scale_color_manual(values = c("low_risk" = "#3B549D", "high_risk" = "#D3272B")) +
  geom_jitter(aes(fill = Group),
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), 
              alpha = 0.2,size = 0.25)
p <- p + stat_compare_means(aes(group = Group),
                            method = "wilcox.test", 
                            label = "p.signif",
                            hide.ns = F)+
  coord_flip()
print(p)

##############################################################################
##KEGG GO GSEA
set.seed(1234)

dat <- matrix[,c(4:ncol(matrix))]
dat <- dat[rownames(matrix_cluster),]
group_list <- as.factor(matrix_cluster$group)

design <- model.matrix(~ 0 + group_list)
design
colnames(design) <- levels(group_list)

dat <- t(dat)%>%as.data.frame()

contrast.matrix <- makeContrasts(high_vs_low = high_risk - low_risk, levels = design)
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, coef = "high_vs_low", number = Inf, adjust.method = "BH")
deg$change <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.585,
                     ifelse(deg$logFC > 0.585, "Up", "Down"),
                     "Stable")
diff_expr <- deg
dt <- deg
names(dt) <- c("log2FoldChange","AveExpr","t","pvalue","adjpVal","B","group")

p <- ggplot() +
  geom_point(data = dt,
             aes(x = log2FoldChange,
                 y = -log10(pvalue),
                 color = group),
             size = 1.6)
p
mycol <- c("#2DB2EB","grey90","#EB4232")
mytheme <- theme_bw() +
  theme(panel.grid=element_blank(),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12))
pvalue = 0.05
log2FC = 0.585

p1 <- p + 
  scale_colour_manual(name = '', values = alpha(mycol, 0.7)) +
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-2, 2, by = 0.1)) + 
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 40),
                     breaks = seq(0, 40, by = 10)) +
  geom_hline(yintercept = c(-log10(pvalue)), size = 0.6, color = "grey30", lty = 'dotdash') +
  geom_vline(xintercept = c(-log2FC, log2FC), size = 0.5, color = "grey30", lty = 'dashed') +
  mytheme
p1

p2 <- ggplot() +
  geom_point(data = dt,
             aes(x = log2FoldChange,
                 y = -log10(pvalue),
                 color = log2FoldChange, 
                 size = -log10(pvalue)),
             alpha = 0.7)
p2

c4a_gui() 
p3 <- p2 +
  scale_color_continuous_c4a_seq(palette = 'spectral',reverse = T) +
  scale_size_continuous(range = c(1,6)) +
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-2, 2, by = 0.5)) + 
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 40),
                     breaks = seq(0, 40, by = 10)) + 
  geom_hline(yintercept = c(-log10(pvalue)), size = 0.6, color = "grey30", lty = 'dotdash') +
  geom_vline(xintercept = c(-log2FC, log2FC), size = 0.5, color = "grey30", lty = 'dashed') +
  mytheme
p3


dt$gene <- rownames(dt)
top10_pos <- dt[order(dt$log2FoldChange, decreasing = TRUE), ][1:10, ]
top10_neg <- dt[order(dt$log2FoldChange, decreasing = FALSE), ][1:10, ]

top_genes <- rbind(top10_pos, top10_neg)

p3 + 
  geom_text_repel(
    data = top_genes,
    aes(x = log2FoldChange, y = -log10(pvalue), label = gene),
    size = 3, 
    box.padding = 0.5,
    max.overlaps = Inf,
    segment.color = "grey50"
  )

#######################################################################
### GO 
library(org.Hs.eg.db)
library(clusterProfiler) 
library(DOSE)
library(enrichplot) 

deg_up <- deg[deg$logFC > 0.585 & deg$adj.P.Val < 0.05, ]
deg_down <- deg[deg$logFC < -0.585 & deg$adj.P.Val < 0.05, ]

gene_up <- rownames(deg_up)
gene_down <- rownames(deg_down)

diff <- c(gene_up, gene_down)
diff_entrez <- bitr(diff,
                    fromType = 'SYMBOL',
                    toType = 'ENTREZID',
                    OrgDb = 'org.Hs.eg.db')

KEGG_diff <- enrichKEGG(gene = diff_entrez$ENTREZID,
                        organism = 'hsa',
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        minGSSize = 10,
                        maxGSSize = 500)

KEGG_diff <- setReadable(KEGG_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENTREZID')

dotplot(
  KEGG_diff, 
  x = 'GeneRatio', 
  color = 'p.adjust',
  showCategory = 20,
  font.size = 12,
  title = 'Top 20 of Pathway Enrichment',
  label_format = 50
)


GO_diff <- enrichGO(gene = diff_entrez$ENTREZID,
                        OrgDb = "org.Hs.eg.db",
                        ont = "ALL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        minGSSize = 10,
                        maxGSSize = 500)
GO_diff <- setReadable(GO_diff,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENTREZID')


library(ggalluvial)
library(ggstyle)
library(scales)
KEGG_filter <- KEGG_diff@result[KEGG_diff@result$p.adjust < 0.05,]

dflong = GO_diff %>% 
  separate_rows(geneID,sep = "/") %>% 
  mutate(Description = factor(Description)) %>% 
  mutate(geneID = factor(geneID))
dfsankey = to_lodes_form(dflong  %>%  select(c("geneID","Description")),
                         key = "x",
                         axes = c(1,2)) %>% 
  mutate(flowcolor = rep(dflong$Description,2))

n_genes <- length(unique(dfsankey$stratum))
stratum_colors <- rep(c("#3b4992", "#F36f43", "#80cba4", "#D62728", "#7d5cc6", 
                        "#fbda63", "#a30543", "#1f918d", "#e9f4a3", "#17BECF"), 
                      length.out = n_genes)

sankeyPlot <- ggplot(data = dfsankey,
       aes(x = x,
           stratum = stratum,
           alluvium = alluvium,
           y = 1)) +
  geom_flow(aes(fill = flowcolor), alpha = 0.3, width = 0, knot.pos = 0.1) +
  geom_stratum(aes(fill = stratum), width = 0.05, color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), 
            size = 3, hjust = 1, nudge_x = -0.03) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.2, 0, 0, 0)) +
  guides(fill = FALSE, color = FALSE) +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), units = "cm")
  ) +
  scale_fill_manual(values = hue_pal()(length(unique(dfsankey$flowcolor))), 
                    aesthetics = "fill", 
                    guide = FALSE) +
  scale_fill_manual(values = stratum_colors, 
                    aesthetics = "fill", 
                    guide = FALSE)

bubbleDf = GO_filter %>% 
  mutate(Term = factor(Description,levels =rev(KEGG_filter$Description))) %>%
  arrange(Term) %>%
  mutate(Term_num =cumsum(Count) - Count /2)

dot_plot <- ggplot(bubbleDf, aes(x = GeneRatio, y = Term_num,color= pvalue)) +
  geom_point(aes(size = Count)) +
  scale_y_continuous(expand =c(0,0),limits =c(0,sum(bubbleDf$Count,na.rm = T)))+
  scale_color_gradient(low ="red", high ="green") +
  scale_radius(range=c(2,5),name="Size") +
  guides(color =guide_colorbar(order =1), 
             size =guide_legend(order =2))+
  theme_bw() +
  labs(size ="Count", color ="PValue", y ="", x ="Ratio") +
  theme(axis.text.y =element_blank(),
        axis.ticks.y =element_blank(), 
        axis.title.y =element_blank(),
        plot.margin =unit(c(0,0,0,0),"inches"),
        panel.grid =element_blank())

sankeyPlot + dot_plot +
  plot_layout(widths =c(2,1))

##########################################################
### GSEA
library(GSEABase)
library(GSVA)
rt = diff_expr[order(diff_expr$logFC, decreasing=TRUE), ]
logFC = rt$logFC
names(logFC) = rownames(rt)

gmt = read.gmt("msigdb.v2023.2.Hs.symbols.gmt")

gmt = read.gmt("c2.all.v2025.1.Hs.symbols.gmt")


gsea = GSEA(logFC, TERM2GENE= gmt, pvalueCutoff= 1,
            minGSSize = 10, maxGSSize = 10000,
            pAdjustMethod = "BH")

gseaTab = as.data.frame(gsea)

gseaTab = gseaTab[gseaTab$p.adjust < 0.05, ]
gseaUp = gseaTab[gseaTab$NES > 0, ]
gseadown <- gseaTab[gseaTab$NES < 0,]

upTerm = c()
downTerm = c("FOROUTAN_INTEGRATED_TGFB_EMT_UP","GAVISH_3CA_METAPROGRAM_FIBROBLASTS_CAF_1",
             "ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
             "NABA_ECM_GLYCOPROTEINS")


gseaplot2(gsea, upTerm, title="Enriched in High", pvalue_table=TRUE,base_size = 8,ES_geom = "line")
gseaplot2(gsea, downTerm, title="Enriched in Low", pvalue_table=TRUE,base_size = 8)


library(GseaVis)
gseaRes <- gsea

term <- c("KEGG_DNA_REPLICATION","REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES_THROUGH_SYNTHESIS_DEPENDENT_STRAND_ANNEALING_SDSA",
          "CROONQUIST_IL6_DEPRIVATION_DN","REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX",
          "MANALO_HYPOXIA_DN")
term <- c("BOQUEST_STEM_CELL_UP","FOROUTAN_TGFB_EMT_UP","HAMAI_APOPTOSIS_VIA_TRAIL_UP",
          "KIM_GERMINAL_CENTER_T_HELPER_UP","HADDAD_T_LYMPHOCYTE_AND_NK_PROGENITOR_UP")

gseaNb(object = gseaRes,
       geneSetID = term,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.8,
       pvalY = 0.5,
       pFill = "transparent",
       addGene = T,
       markTopgene = T,
       geneCol = "black",
       )

##########################################################
### oncopredict
rt <- readRDS("./all_data.RDS") %>%
  .[,4:ncol(.)] %>%
  select(., -Group) %>%
  t()
max(rt)
min(rt)
expr_num <- log2(rt + 1)
GDSC2_Expr <- readRDS("DataFiles\\Training Data\\GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS("DataFiles\\Training Data\\GDSC2_Res.rds")
GDSC2_Res = exp(GDSC2_Res)

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = expr_num,
              batchCorrect = "eb",
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10,
              selection = 1,
              printOutput = TRUE,
              removeLowVaringGenesFrom = "homogenizeData")


res <- read.csv("DrugPredictions.csv",row.names = 1)
risk_matrix <- readRDS("F:\\LLPS_DLBCL\\all_data.RDS")
res$Group <- risk_matrix$Group[match(rownames(res),risk_matrix$ID)]
long_res_all <- pivot_longer(res, cols = -Group, names_to = "drugs", values_to = "IC50")

p_values <- sapply(unique(long_res_all$drugs), function(drug) {
  res <- long_res_all[long_res_all$drugs == drug, ]
  test_result <- wilcox.test(IC50 ~ Group, data = res)
  return(test_result$p.value)
}, USE.NAMES = TRUE)

p_values_df <- data.frame(drug = names(p_values), p_value = p_values)
p_values_sorted <- p_values_df[order(p_values_df$p_value), ]

significant_results <- p_values_sorted[p_values_sorted$p_value < 0.00000001, ] #-8e 
qualified_drugs <- significant_results$drug
long_qualified_clinicaluse <- long_res_all[long_res_all$drugs %in% qualified_drugs,]
long_qualified_clinicaluse$IC50 <- as.numeric(long_qualified_clinicaluse$IC50)

selected_drug <- "Sorafenib_1085"
long_test <- long_qualified_clinicaluse[long_qualified_clinicaluse$drugs == selected_drug,]

y_lower_limit <- 5
y_upper_limit <- 40

ggplot(long_test, aes(x = Group, y = IC50, fill = Group)) +
  geom_violin(trim = TRUE, color = "black",aes(fill = Group), alpha = 0.3, scale = "width",size = 1.5) +
  geom_boxplot(width = 0.6, aes(color = Group), color = "black",alpha = 0.75, 
               outlier.shape = NA, size = 1.8, notch=TRUE) +
  geom_jitter(shape = 21, size = 2, aes(fill = Group), color = "black", stroke = 0.8,
              position = position_jitter(0.15)) +
  stat_compare_means(aes(label = paste0("Wilcoxon, p = ", ..p.format..)),
                     method = "wilcox.test", label.x = 1.5, label.y = y_upper_limit) +
  scale_fill_manual(values = c("#9392BE", "#74B69F")) +
  scale_color_manual(values = c("#9392BE", "#74B69F")) +
  ylim(y_lower_limit, y_upper_limit) +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 1.0),  
        axis.line.y = element_line(color = "black", size = 1.0),  
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none",
        aspect.ratio = 1.5,
        panel.grid.major = element_line(color = "white"), 
        panel.grid.minor = element_blank(), )+  
  labs(x = "Expression of risk group", y = "Estimated IC50", 
       title = selected_drug)

#################################################################
llps_model <- readRDS("llps_model.RData")
coefficients <- llps_model$coefficients
coefficients

computing_score <- risk_matrix[,names(coefficients)]
for (i in seq_along(names(llps_model$coefficients))) {
  computing_score[, i] <- computing_score[, i] * unname(llps_model$coefficients)[i]
  }
computing_score <- computing_score %>% mutate(risk_score = rowSums(computing_score[,1:ncol(computing_score)]))
computing_score  <- mutate(computing_score , Group = ifelse(risk_score > -5.4970,"high_risk","low_risk"))


res$risk_score <- computing_score$risk_score[match(rownames(res),rownames(computing_score))]
drug_columns <- colnames(res)[1:198]
cor_results <- list()
# 计算相关性和 p 值
for (drug in drug_columns) {
  test_result <- cor.test(res[[drug]], res$risk_score)
  cor_results[[drug]] <- list(
    correlation = test_result$estimate,
    p.value = test_result$p.value
  )
}
cor_df <- do.call(rbind, lapply(cor_results, function(x) data.frame(correlation = x$correlation, p.value = x$p.value)))
cor_df$drug <- rownames(cor_df)
rownames(cor_df) <- NULL



library(ggplot2)
ggplot(cor_qualified, aes(x = reorder(drug, correlation), 
                          y = correlation, 
                          fill = correlation > 0)) +
  geom_col(aes(alpha = Log10P), width = 0.7) +   
  coord_flip() +
  geom_text(aes(label = paste0(drug, " (", round(correlation, 2), ")"),
                y = ifelse(correlation > 0, correlation - 0.05, correlation + 0.05),
                hjust = ifelse(correlation > 0, 1, 0)),
            color = "black", size = 4) +
  scale_fill_manual(values = c("FALSE" = "#e3716e", "TRUE" = "#2c7fb8"),
                    labels = c("Negative", "Positive"),
                    name = "Direction") +
  scale_alpha_continuous(range = c(0.3, 1), name = "-log10(p)") +
  geom_hline(yintercept = 0, color = 'black', size = 0.5, lty = 'dashed') +
  
  labs(x = NULL, y = "Correlation", 
       title = "Correlation of Drugs with Riskscore") +
  
  theme_bw() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 14),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
    panel.border = element_blank(),
    axis.line.x = element_line(size = 1, color = "black"),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

#####################################################################
expr_6genes <- risk_matrix[names(coefficients)]
expr_IC50 <- cbind(res[1:198],expr_6genes)
library(pheatmap)
drug_cols <- 1:198
gene_cols <- 199:ncol(expr_IC50)
results <- data.frame()
for (drug in colnames(expr_IC50)[drug_cols]) {
  for (gene in colnames(expr_IC50)[gene_cols]) {
    test <- cor.test(expr_IC50[[drug]], expr_IC50[[gene]], method = "spearman")
    results <- rbind(results, data.frame(
      Drug = drug,
      Gene = gene,
      Cor = test$estimate,
      Pval = test$p.value
    ))
  }
}
results_sig <- results %>% filter(Pval < 0.05 & Cor > 0.3)
ggplot(results_sig, aes(x = Gene, y = Drug, size = abs(Cor), color = -log10(Pval))) +
  geom_point(alpha = 0.9) +
  scale_color_viridis(option = "plasma", name = "-log10(P)") +
  scale_size(range = c(2, 8), name = "|Spearman R|") +
  labs(
    title = "Gene-Drug Correlation (p < 0.05)",
    x = "Gene",
    y = "Drug"
  ) +
  theme_minimal(base_size = 12) +
  coord_flip()+
  theme(
    panel.grid = element_blank(),             
    axis.line = element_line(color = "black"), 
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )
intersect(results_sig$Drug,clinical_use_drug)


load("TCGAdata_RNA\\test\\output_mRNA_lncRNA_expr\\TCGA-DLBC_mrna_expr_counts.rdata")
tcga_counts <- as.data.frame(t(mrna_expr_counts)) 
tcga_counts <- log2(tcga_counts)+1

llps_model <- readRDS("llps_model.RData")
coefficients <- llps_model$coefficients
coefficients

risk_matrix_tcga <- tcga_counts[,names(coefficients)]
for (i in seq_along(names(coefficients))) {
  risk_matrix_tcga[, i] <- risk_matrix_tcga[, i] * unname(coefficients)[i]
}
risk_matrix_tcga<- risk_matrix_tcga %>% mutate(risk_score = rowSums(risk_matrix_tcga))
rownames(risk_matrix_tcga) <- substr(rownames(risk_matrix_tcga), 1, 16)

tcga_surv <- read.delim("TCGA-DLBC.survival.tsv")
tcga_surv <- tcga_surv[,-3]
tcga_surv <- tcga_surv[tcga_surv$sample %in% rownames(risk_matrix_tcga),]
risk_matrix_tcga <- risk_matrix_tcga[rownames(risk_matrix_tcga) %in% tcga_surv$sample,]
risk_matrix_tcga <- risk_matrix_tcga[tcga_surv$sample,]
res_tcga <- cbind(tcga_surv,risk_matrix_tcga)
rownames(res_tcga) <- res_tcga[,1]
res_tcga <- res_tcga[,-1]
res_tcga$OS.time <- res_tcga$OS.time/365

risk_matrix_tcga$risk_group <- ifelse(risk_matrix_tcga$risk_score > -5.4970, "high", "low")

high_risk <- rownames(risk_matrix_tcga[risk_matrix_tcga$risk_group == "high", ])
low_risk <-  rownames(risk_matrix_tcga[risk_matrix_tcga$risk_group == "low", ])

library(maftools)
tcga_maf <- data.table::fread("dlbc_mutation.maf")
tcga_maf$sample <- substr(tcga_maf$Tumor_Sample_Barcode, 1, 16)
tcga_maf
kkkk <- intersect(tcga_maf$sample,low_risk)
llll <- intersect(tcga_maf$sample,high_risk)
##############################
primary.apl = system.file("extdata", "APL_primary.maf.gz", package = "maftools")
primary.apl = read.maf(maf = primary.apl)
tcga_maf_mafobj <- read.maf(maf = tcga_maf)
#########
##高低危分组
tcga_maf_low <- as.data.frame(tcga_maf[tcga_maf$sample %in% kkkk]) 
use_maf_low <- read.maf(maf = tcga_maf_low)
plotmafSummary(maf = use_maf_low, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)
oncoplot(maf = use_maf_low)
tmb_low <- tmb(maf = use_maf_low)

tcga_maf_high <- as.data.frame(tcga_maf[tcga_maf$sample %in% llll]) 
use_maf_high <- read.maf(maf = tcga_maf_high)
plotmafSummary(maf = use_maf_high, rmOutlier = T, addStat = 'median', dashboard = T, titvRaw = F)
oncoplot(maf = use_maf_high)
tmb_high <- tmb(maf = use_maf_high)


tmb_low <- as.data.frame(tmb_low)
tmb_high <- as.data.frame(tmb_high)
tmb_all <- rbind(tmb_low,tmb_high)
tmb_all$Tumor_Sample_Barcode <- substr(tmb_all$Tumor_Sample_Barcode, 1, 16)
rownames(tmb_all) <- tmb_all[,1]
tmb_all <- tmb_all[,-1]
rescat_45 <- risk_matrix_tcga[rownames(risk_matrix_tcga) %in% rownames(tmb_all),]
tmb_merge <- merge(tmb_all, rescat_45, by = "row.names", all = TRUE)

rownames(tmb_merge) <- tmb_merge[,1]
tmb_merge <- tmb_merge[,-1]


pt.vs.rt <- mafCompare(m1 = use_maf_high, m2 = use_maf_low, m1Name = 'High_risk', m2Name = 'Low_risk', minMut = 5) 
print(pt.vs.rt)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.5)
coBarplot(m1 = use_maf_high, m2 = use_maf_low, m1Name = 'High_risk', m2Name = 'Low_risk')

ggboxplot(tmb_merge, x = "risk_group", y = "total",
          fill = "risk_group", palette = c("#E69F00", "#56B4E9"),
          add = "jitter") +
  stat_compare_means(method = "wilcox.test",    # 非参数检验
                     label = "p.format",       # 只显示p值
                     label.y = max(tmb_merge$total) * 1.05) +
  labs(title = "TMB", x = NULL, y = "Total") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14))


res_tcga_45 <- risk_matrix_tcga[rownames(risk_matrix_tcga) %in% rownames(rescat_45),]
res_tcga_45 <- merge(res_tcga_45, tmb_all, by = "row.names", all = TRUE)
res_tcga_45$Group <- risk_matrix_tcga$risk_group[match(res_tcga_45$Row.names,rownames( risk_matrix_tcga))]

y <- res_tcga_45$total
x <- res_tcga_45$risk_score
corT=cor.test(x,y,method="spearman")
cor=corT$estimate
pValue=corT$p.value
cor_label <- paste0("R = ", round(cor, 3), ", p = ", format(pValue, digits = 3))
p1 <- ggplot(res_tcga_45, aes(x, y, color = Group)) + 
  xlab('risk score') +
  ylab('TMB') +
  geom_point() + 
  geom_smooth(method = "lm", formula = y ~ x, aes(group = 1), color = "midnightblue") + 
  theme_bw() +
  annotate("text", x = min(x), y = max(y), label = cor_label, size = 5, hjust = 0) + 
  scale_color_manual(values = c("firebrick1", "dodgerblue")) +
  guides(color = guide_legend(override.aes = list(linetype = "blank", shape = 16)))+
  theme(
    legend.position = c(0.08, 0.85),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.key.size = unit(1, 'lines'),
    legend.text = element_text(size = 10)
  )
library(ggExtra)
p2 = ggMarginal(p1, type = "density", xparams = list(fill = "firebrick1"), yparams = list(fill = "dodgerblue"))
print(p2)

###################################################
###nomogram
library(rms)
library(pROC)
library(timeROC)
library(ciTools)
#############################Training
#clinical imformation
file_path <- "GSE11318-GPL570_series_matrix.txt"
all_lines <- readLines(file_path, n = 100)
selected_lines <- all_lines[30:83]
temp_file <- tempfile()
writeLines(selected_lines, temp_file)
data <- read.delim(temp_file, header = FALSE, sep = "\t")
unlink(temp_file)

clinical_data <- data[c(2,10,11,16,17,18,19,20,21,22,23),]
clinical_data <- clinical_data[,-1]
new_row_names <- c("sample_name", "sex","age","diagnosis","OS","OS.time","regimen","ECOG","stage","LDHratio","extranodal_sites_number")

extract_content <- function(x) {
  sub(".*:\\s*", "", x)
}
clinical_data[] <- lapply(clinical_data, function(column) {
  sapply(column, extract_content)
})
tdata <- as.data.frame(t(clinical_data))
tdata <- tdata %>%
  mutate(OS = ifelse(OS == "ALIVE", 0, 1))
tdata <- tdata[,-c(4,7)]
tdata$sample_name <- paste0("GSE11318_", tdata$sample_name)
rownames(tdata) <- tdata[,1]
tdata <- tdata[,-1]
########################################################
llps_model <- readRDS("llps_model.RData")
coefficients <- llps_model$coefficients
coefficients

gse11318 <- grep("^GSE11318_", rownames(risk_matrix), value = FALSE)
gse11318 <- risk_matrix[gse11318, ]
gse11318 <- gse11318[,names(coefficients)]

for (i in seq_along(names(coefficients))) {
  gse11318[, i] <- gse11318[, i] * unname(coefficients)[i]
}
gse11318 <- gse11318 %>% mutate(risk_score = rowSums(gse11318))
tdata$riskscore <- gse11318$risk_score[match(rownames(tdata),rownames(gse11318))]
tdata <- tdata[!is.na(tdata$riskscore), ]
tdata$OS.time <- as.numeric(tdata$OS.time)
tdata$sex <- ifelse(is.na(tdata$sex), NA, ifelse(tdata$sex == "male", 0, 1))
tdata$sex <- as.integer(tdata$sex)
tdata$age <- as.numeric(tdata$age)
tdata$ECOG <- as.numeric(tdata$ECOG)
tdata$stage <- as.integer(tdata$stage)
tdata$extranodal_sites_number <- as.numeric(tdata$extranodal_sites_number)
tdata$LDHratio <- as.numeric(tdata$LDHratio)


univ_cox_results <- list()
vars <- c("age", "sex", "ECOG","stage","extranodal_sites_number","LDHratio","riskscore")
for (var in vars) {
  formula <- as.formula(paste("Surv(OS.time, OS) ~", var))
  cox_result <- coxph(formula, data = tdata)
  univ_cox_results[[var]] <- summary(cox_result)
}
significant_vars <- names(univ_cox_results)[sapply(univ_cox_results, function(x) x$coefficients[,"Pr(>|z|)"]) < 0.05]
significant_vars
univ_cox_summary <- do.call(rbind, lapply(univ_cox_results, function(x) {
  data.frame(
    coef  = x$coefficients[,"coef"],
    HR    = x$coefficients[,"exp(coef)"],
    se    = x$coefficients[,"se(coef)"],
    z     = x$coefficients[,"z"],
    p     = x$coefficients[,"Pr(>|z|)"],
    HR.95L = x$conf.int[,"lower .95"],
    HR.95U = x$conf.int[,"upper .95"]
  )
}))


dd <- datadist(tdata)
options(datadist = "dd")
mul_cox <- cph(Surv(OS.time, OS) ~ age + stage + riskscore, 
               data = tdata, x = TRUE, y = TRUE, surv = TRUE)

sur <- Survival(mul_cox)
sur1 <- function(x) sur(1, x)   # 1年生存
sur2 <- function(x) sur(2, x)   # 2年生存
sur3 <- function(x) sur(3, x)   # 3年生存
sur4 <- function(x) sur(4, x)   # 4年生存
sur5 <- function(x) sur(5, x)   # 5年生存
funlist <- list(sur1,sur3,sur5)
nom <- nomogram(mul_cox, 
                fun = funlist, 
                fun.at = c(0.1,0.3,0.5,0.7,0.9),
                funlabel = c("1 year survival",  
                             "3 years survival","5 years survival"),
                lp = FALSE)  

par(mar = c(5, 5, 4, 2))
plot(nom,
     xfrac = 0.3,
     lwd = 2,
     col.grid = "grey80",
     col.axis = "black",
     col.var = "black",
     cex.axis = 1.1,
     cex.var = 1.3,
     cex.fun = 1.2,
     lmgp = 0.3)

library(regplot)
mul_cox_2 <- coxph(Surv(OS.time,OS) ~ age + stage + LDHratio + riskscore, data = tdata) 
mul_cox_2
par(mar = c(5, 5, 4, 2)) 
regplot(mul_cox_2,
        observation = tdata[120, ],
        points = TRUE,
        plots = c("density", "no plot"),
        failtime = c(1, 2, 3, 4, 5),
        odds = FALSE,
        leftlabel = TRUE,
        prfail = TRUE,
        showP = TRUE,
        droplines = TRUE,
        colors = "steelblue",
        rank = "range",
        interval = "confidence",
        title = "Nomogram",
        cex.axis = 1.2,
        cex.var = 1.3,
        cex.points = 1.2,
        cex.title = 1.4,
        lwd = 2)

ggforest(mul_cox_2,
         data= tdata,
         main = "Hazard ratio", 
         cpositions = c(0.01, 0.1, 0.3),
         fontsize = 0.8,
         refLabel = "reference",
         noDigits = 3) 

par(bty = "l")
predicted_survival <- predict(mul_cox, type = "lp") 
roc_1_year <- timeROC(T = tdata$OS.time, delta = tdata$OS, marker = predicted_survival, 
                      cause = 1, times = 1, iid = TRUE)
roc_3_years <- timeROC(T = tdata$OS.time, delta = tdata$OS, marker = predicted_survival, 
                       cause = 1, times = 3, iid = TRUE)
roc_5_years <- timeROC(T = tdata$OS.time, delta = tdata$OS, marker = predicted_survival, 
                       cause = 1, times = 5, iid = TRUE)
par(bty = "l") 
plot(roc_1_year , time = 1, title = " ",col = "firebrick1", xlab = "1 - Specificity", ylab = "Sensitivity")
plot(roc_3_years, time = 3, add = TRUE, col = "dodgerblue")
plot(roc_5_years, time = 5, add = TRUE, col = "limegreen")
title(main = "Time dependent ROC curves of nomogram", line = 0.5)

auc_1_year <- round(roc_1_year$AUC[2], 3)
ci_1_year <- confint(roc_1_year, level = 0.95,n.sim=2000)
ci_1_year_text <- paste(sprintf("%.3f", round(ci_1_year$CI_AUC[,"2.5%"]/100, 4)), "~", sprintf("%.3f", round(ci_1_year$CI_AUC[,"97.5%"]/100, 4)))
auc_3_years <- round(roc_3_years$AUC[2], 3)
ci_3_year <- confint(roc_3_years, parm=NULL, level = 0.95,n.sim=2000)
ci_3_year_text <- paste(sprintf("%.3f", round(ci_3_year$CI_AUC[,"2.5%"]/100, 4)), "~", sprintf("%.3f", round(ci_3_year$CI_AUC[,"97.5%"]/100, 4)))
auc_5_years <- round(roc_5_years$AUC[2], 3)
ci_5_year <- confint(roc_5_years, parm=NULL, level = 0.95,n.sim=2000)
ci_5_year_text <- paste(sprintf("%.3f", round(ci_5_year$CI_AUC[,"2.5%"]/100, 4)), "~", sprintf("%.3f", round(ci_5_year$CI_AUC[,"97.5%"]/100, 4)))

legend("bottomright", 
       legend = c(paste("1 year\n(AUC =", auc_1_year, ", 95%CI:", ci_1_year_text,")"),
                  paste("3 years\n(AUC =", auc_3_years, ", 95%CI:", ci_3_year_text,")"),
                  paste("5 years\n(AUC =", auc_5_years, ", 95%CI:", ci_5_year_text,")")), 
       col = c("firebrick1", "dodgerblue", "limegreen"),
       lwd = 2, cex = 0.4, text.font = 7)
