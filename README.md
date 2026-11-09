# KINIC-score

版本 / Release: **v1.0.0**

用于计算胶质瘤患者KINIC score的锁定24基因随机生存森林（RSF）。
This release provides the locked 24-gene random survival forest, scoring script and usage instructions for cohort-level KINIC scoring in diffuse glioma.

## 文件 / Files

- `KINIC_RSF_model.rds`：已拟合森林、输入基因顺序、评分方向和精确截点。
- `score_KINIC.R`：独立评分程序。
- `README.md`：输入、安装及运行说明。

## 软件 / Requirements

R >= 4.3.0；本次验证使用R 4.3.3及 **randomForestSRC 3.3.3**。无需安装Mime。
R >= 4.3.0 and randomForestSRC **3.3.3**. The frozen model was checked with R 4.3.3; Mime is not required for scoring.

```r
install.packages(c("data.tree", "DiagrammeR"), repos = "https://cloud.r-project.org")
install.packages("https://cran.r-project.org/src/contrib/Archive/randomForestSRC/randomForestSRC_3.3.3.tar.gz", repos = NULL, type = "source")
```

源码安装需要适合本机R的编译工具。The source installation requires the compiler toolchain for your R installation.

## 输入 / Input

输入为TSV：每行一个样本，首列名为`sample_id`，其余列包含下列24个基因的数值表达；多余列忽略，程序按模型顺序排列。样本ID与列名须唯一，模型基因不得缺失、重复或含非有限值。
Use a tab-separated table with samples in rows, `sample_id` as the first column, and numeric expression columns for all 24 genes below. Extra columns are ignored; the script reorders model genes. Keep unique sample IDs and exact gene symbols, with complete finite expression values.

```text
AMZ1, C2orf66, CKS2, CLVS1, EDN3, EN1, GRIA1, HIST1H2BJ, HOXA10, KLHDC8A, KSR2, LRRTM3, MN1, PRLHR, RCC1, SLC24A4, SLC30A7, SMC4, STK40, TRH, SYBU, TFPI, HOXD8, ISL2
```

HIST1H2BJ与C2orf66等名称保留训练时的标识，程序不自动替换基因别名。Symbols such as HIST1H2BJ and C2orf66 retain their training identifiers; aliases are not mapped automatically.

## 运行 / Run

将三个文件置于同一目录。每次只处理一个完整目标队列。Keep the three files together and process one target cohort per invocation.

```sh
# CGGA样式的原始RSEM FPKM：先log2(FPKM+1)，再按队列标准化
Rscript score_KINIC.R expression.tsv KINIC_scores.tsv fpkm

# 已完成平台表达变换的log2表达，例如原分析的Xena HiSeqV2矩阵
Rscript score_KINIC.R expression.tsv KINIC_scores.tsv log2

# 已按适当参考完成标准化的表达：直接输入模型，不再次标准化
Rscript score_KINIC.R standardized_expression.tsv KINIC_scores.tsv zscore
```

`fpkm`先计算log2(x+1)；`fpkm`与`log2`均使用当前队列逐基因均值和样本标准差（n−1）计算Z-score。`zscore`直接使用输入。
`fpkm` applies log2(x+1). Both `fpkm` and `log2` then center each gene by the current cohort mean and divide by its sample SD (n−1); `zscore` uses the supplied values directly.

原分析标准化分母对应最终纳入队列：CGGA325 313例、CGGA693 657例、TCGA 670例。复现原评分须保持相同样本集合与表达处理。不要把不同队列拼在一起标准化，或将已标准化的小样本子集再标准化。
The original normalization cohorts contained 313 CGGA325, 657 CGGA693 and 670 TCGA patients. Reproducing their scores requires the same included samples and expression processing. Normalize cohorts separately and do not re-standardize a subset of already standardized samples.

队列内标准化依赖该队列的表达分布。单个样本可使用已按指定参考标准化的`zscore`输入；单患者应用所需的固定表达参考和检测流程须另行确定。
Cohort normalization depends on the cohort expression distribution. A single sample can be supplied in `zscore` mode after normalization against an explicitly defined reference; single-patient application requires a fixed expression reference and assay workflow.

## 输出 / Output

输出为`sample_id`、`KINIC_score`、`KINIC_group`三列。Score为锁定森林的survival mortality预测，按训练阶段确定的方向输出；分数越高对应更高的模型预后风险。它不是生存概率。
The output columns are `sample_id`, `KINIC_score` and `KINIC_group`. The score is the forest survival-mortality prediction with the frozen training orientation; higher scores represent higher model prognostic risk. It is not a survival probability.

固定训练截点 / Fixed CGGA325 median cutoff: **101.69034285714299**。高于截点为High，其余为Low；不在新队列重选截点。
Values strictly above the cutoff are High; all other values are Low. The cutoff is not re-estimated in the target cohort.

## 核验 / Verification

已从项目冻结输入核对CGGA325、CGGA693和TCGA共1640例；导出森林与原模型评分一致，最大评分绝对差为1.05e-12，High/Low分组全部一致。原始表达预处理、独立CLI调用及重载模型均经核对，未重新训练模型。
The exported forest and standalone scoring program were checked against frozen scores for 1,640 patients across the discovery cohort and both external cohorts, including raw-expression preprocessing. All High/Low assignments matched; the model was not retrained.

公共数据来源 / Public data sources: [CGGA](https://www.cgga.org.cn/)；[UCSC Xena](https://xenabrowser.net/)。
固定版本 / Versioned release: https://github.com/MeiniuZhao/KINIC-score/releases/tag/v1.0.0

