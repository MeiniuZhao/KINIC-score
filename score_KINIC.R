#!/usr/bin/env Rscript
# KINIC-score v1.0.0：加载锁定RSF，按完整队列计算评分与固定截点分组。
# 输入：TSV，首列sample_id，其余列包含README列出的全部24个基因。
# 运行：Rscript score_KINIC.R input.tsv scores.tsv log2 [KINIC_RSF_model.rds]
# expression_scale：fpkm、log2或zscore；分别处理一个队列。

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Usage: Rscript score_KINIC.R input.tsv output.tsv {fpkm|log2|zscore} [model.rds]")
}
input_path <- args[1L]
output_path <- args[2L]
expression_scale <- args[3L]
if (!expression_scale %in% c("fpkm", "log2", "zscore")) {
  stop("expression_scale must be fpkm, log2 or zscore.")
}
this_script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
model_path <- if (length(args) == 4L) args[4L] else
  file.path(dirname(normalizePath(this_script)), "KINIC_RSF_model.rds")

if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
  stop("Install randomForestSRC 3.3.3 as described in README.md.")
}
if (as.character(utils::packageVersion("randomForestSRC")) != "3.3.3") {
  stop("This frozen model release uses randomForestSRC 3.3.3; see README.md.")
}
suppressPackageStartupMessages(library(randomForestSRC))
options(rf.cores = 1L)
set.seed(20261109L)
model <- readRDS(model_path)

# 精确匹配名称并按原模型顺序排列，缺失基因或非数值表达直接报错。
input <- utils::read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE,
                          colClasses = c(sample_id = "character"), na.strings = c("NA", ""))
if (names(input)[1L] != "sample_id" || anyDuplicated(names(input))) {
  stop("The first column must be sample_id, and column names must be unique.")
}
sample_ids <- as.character(input$sample_id)
if (!length(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) {
  stop("sample_id must contain unique, non-missing sample identifiers.")
}
missing_genes <- setdiff(model$gene_symbols, names(input))
if (length(missing_genes)) stop("Missing required genes: ", paste(missing_genes, collapse = ", "))
expression <- input[, model$gene_symbols, drop = FALSE]
if (!all(vapply(expression, is.numeric, logical(1)))) stop("Model gene columns must be numeric.")
expression <- as.matrix(expression)
if (!all(is.finite(expression))) stop("Model gene expression must be finite and complete.")

# 与原分析一致：FPKM先log2(x+1)，随后按目标队列各基因的均值/样本SD标准化。
# log2模式的输入应为已完成相应平台表达变换的矩阵；zscore模式直接使用输入。
if (expression_scale == "fpkm") {
  if (any(expression < 0)) stop("FPKM expression must be non-negative.")
  expression <- log2(expression + 1)
}
if (expression_scale != "zscore") {
  if (nrow(expression) < 2L) stop("Cohort standardization requires at least two samples.")
  cohort_mean <- colMeans(expression)
  cohort_sd <- apply(expression, 2L, stats::sd)
  if (any(!is.finite(cohort_sd) | cohort_sd <= 0)) {
    stop("Every model gene needs a finite, positive within-cohort sample SD.")
  }
  expression <- scale(expression, center = cohort_mean, scale = cohort_sd)
}
colnames(expression) <- model$model_features
rownames(expression) <- sample_ids

# 必须显式提供newdata：得到固定森林对这些表达的预测，而非训练集OOB值。
raw_score <- as.numeric(stats::predict(model$forest,
  newdata = as.data.frame(expression), na.action = "na.omit")$predicted)
score <- model$orientation_multiplier * raw_score
group <- ifelse(score > model$cutoff, "High", "Low")
output <- data.frame(sample_id = sample_ids, KINIC_score = sprintf("%.17g", score),
                     KINIC_group = group, check.names = FALSE)
utils::write.table(output, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Scored %d samples using %d genes; input_scale=%s; fixed_cutoff=%.17g\n",
            length(score), ncol(expression), expression_scale, model$cutoff))
