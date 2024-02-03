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

使用[TCGA数据库](https://portal.gdc.cancer.gov/)通过TCGA-LIHC数据集(2023年1月,n=377)，下载clinical（临床信息），点击Metadata（样本信息），cart（基因文件）[/data/TCGA数据库]获取转录组测序数据与临床数据，并排除了临床数据不完整的患者队列

我们使用了Harmony软件包进行批次效应调整，

接下来，基于拟合多因素Cox回归和lasso回归，构建预后效果评价分类模型  

以"stranded_first"作为特定基因的表达水平评价量，临床信息中选择vital_status：生存状态（例如，存活、死亡）和days_to_death：诊断后至死亡的天数，作为评价的临床指标，
最后得到如下结果：[data/washed_data/TCGA/data_exp_clinical.xlsx]

最后整合得到如下表格：
| 病例ID | 基因名 | 生存时间 | 结局 |
| ------ | ------ | ------- | ---- |
| ID1    | GeneA  | 时间1   | 结果1 |
| ID2    | GeneB  | 时间2   | 结果2 |
| ID3    | GeneC  | 时间3   | 结果3 |

### 描述性分析
首先使用生存分析中的Cox比例风险模型建立一个风险评分系统，然后利用风险评分中位数对其进行分组，将其分为“高风险组”和”低风险组“，并绘制生存分析图（针对筛选出的基因）





## 4、KEGG富集分析
KEGG富集分析前40个结果涉及多个与癌症相关的信号通路，如“癌症通路”、“胰腺癌”、“前列腺癌”和“慢性髓性白血病”，这表明基因列表可能与癌症进程有关，或来源于癌症研究。

此外，还包含了多个与病毒感染相关的通路，如“卡波西肉瘤相关疱疹病毒感染”、“人类T细胞白血病病毒1型感染”和“乙型肝炎”，这些病毒感染可能与癌症相关，因为部分病毒被认为具有致癌性。

还有一些关键的信号通路，如“HIF-1信号通路”、“FoxO信号通路”和“PI3K-Akt信号通路”，这些通路对细胞生存、生长和代谢至关重要，通常在癌症中发生紊乱。

免疫应答和炎症相关的通路，如“TNF信号通路”和“C型凝集素受体信号通路”，这些通路在感染和癌症的背景下对免疫应答和炎症反应至关重要。

此外，“细胞周期”和“凋亡”等通路是细胞分裂和死亡的基础，它们的紊乱是癌症的标志之一。

“内分泌抵抗”和“EGFR酪氨酸激酶抑制剂抵抗”等通路表明研究关注了导致治疗抵抗的机制，这在癌症治疗中是一个重大挑战。

还有与代谢和其他疾病相关的通路，如“胰岛素抵抗”和“糖尿病并发症中的AGE-RAGE信号通路”，这表明代谢失调，也可能与癌症和其他疾病相关。


## 5、分子验证

### 5.1 分子对接验证

### 5.2 分子动力学验证

### 5.3 关键靶点基因的生存分析
