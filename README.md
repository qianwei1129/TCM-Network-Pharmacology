# TCM-Network-Pharmacology

***中药网络药理学全过程***

## 网站和工具

### 成分/靶点获取网站

[TCMSP](https://old.tcmsp-e.com/tcmsp.php) [TCMBank](https://www.tcmbank.cn/) [HERB](http://herb.ac.cn/) [swisstarget](http://www.swisstargetprediction.ch/) [PubChem](https://pubchem.ncbi.nlm.nih.gov/)

### 疾病数据库

[DRUGBANK](https://www.drugbank.ca/) [Digsee](http://210.107.182.61/geneSearch/) [TTD](http://db.idrblab.net/ttd/) [Genecards](https://www.genecards.org/) [OMIM](https://omim.org/) [Disgenet](https://www.disgenet.org/) [Malacards](https://www.malacards.org/) [TCGA](https://portal.gdc.cancer.gov/)[GEO](https://www.ncbi.nlm.nih.gov/geo/) [HAGR](https://genomics.senescence.info/)[msigdb](https://www.gsea-msigdb.org/gsea/msigdb)

### 工具

[R](https://www.r-project.org/) [cytoscape](https://cytoscape.org/) [systemsDOCK](http://systemsdock.unit.oist.jp/iddp/home/index) [CPI-BIO](http://cpi.bio-x.cn/drar/) [PharmMapper](https://omictools.com/pharmmapper-tool)





## 1、黄芪六一汤有效成分的筛选

使用[TCMSP数据库](https://old.tcmsp-e.com/tcmsp.php)，选择检索种类为"Herb name"，检索黄芪。然后根据公认标准筛选口服生物利用度（OB≥30%）和药物相似性（DL≥0.18），共获得黄芪的20个作用化合物(\data\pre\chemicals\_hq.xlsx)和甘草的69个化合物。 然后根据smile值，使用[swisstarget数据库](http://www.swisstargetprediction.ch/)获取化合物作用靶点，其中部分化合物没有smile值则使用[TCMBank数据库](https://www.tcmbank.cn/)，或者[HERB](http://herb.ac.cn/)获取作用靶点





## 2、衰老相关基因

从[HAGR](https://genomics.senescence.info/)和[msigdb](https://www.gsea-msigdb.org/gsea/msigdb)数据获取细胞衰老相关基因，将两者取交集后构建基因蛋白互作网络


### HAGR数据库

该库本身提供了下载链接，我在下载后对其进行了清洗

### msigdb数据库
以"aging"作为关键词，Search Filters中collection设置为"all collections"，source species设置为"Homo sapiens",contributors设置为"all contributions"

在msigdb数据库中共得到56个衰老相关基因集





## 3、小细胞肝癌预后相关基因_TCGA数据库

```R
library("survival")
library("survminer")
```

使用[TCGA数据库](https://portal.gdc.cancer.gov/)通过TCGA-LIHC数据集(2023年1月,n=377)，
下载clinical（临床信息），点击Metadata（样本信息），cart（基因文件）[/data/TCGA数据库]获取转录组测序数据与临床数据  

接下来，基于拟合多因素Cox回归和lasso回归，构建预后效果评价分类模型  

以"stranded_first"作为特定基因的表达水平评价量，临床信息中选择vital_status：生存状态（例如，存活、死亡）和days_to_death：诊断后至死亡的天数，作为评价的临床指标（这边真的，这个数据集上传者职业素质真差，给出的临床数据好多错位的），
最后得到如下结果：[data/washed_data/TCGA/data_exp_clinical.xlsx]

最后整合得到如下表格：
| 病例ID | 基因名 | 生存时间 | 结局 |
| ------ | ------ | ------- | ---- |
| ID1    | GeneA  | 时间1   | 结果1 |
| ID2    | GeneB  | 时间2   | 结果2 |
| ID3    | GeneC  | 时间3   | 结果3 |

### 描述性分析
首先根据生存分析（或者lasso回归？）建立一个风险评分系统，对其进行分组

基于聚类分析，首先使用数据集群性评估，计算Hopkins统计量为0.7883237，表示数据高度倾向于聚类，是进行聚类分析的好候选  

然后选择PCA聚类分析，查找累计方差解释达到80%的位置，以选择主成分的个数，选择出3个最小主成分数



然后使用函数clusGap()来计算用于估计最优聚类数，所以选择
当聚类变量（k）设置为2时，组内相关性最高，组间相关性最低，
所以先用K均值聚类针对临床数据将数据分成两组，绘制Kaplan-Meier曲线









## 4、小细胞肝癌预后相关基因_GEO数据库

[GEO数据库](https://www.ncbi.nlm.nih.gov/geo/)筛选小细胞肝癌(HCC)预后相关基因





## 5、富集分析与分子验证

### 5.1 分子对接验证

### 5.2 分子动力学验证

### 5.3 关键靶点基因的生存分析
