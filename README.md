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

## 2、肿瘤与衰老疾病共同作用靶点

使用[TCGA数据库](https://portal.gdc.cancer.gov/)和[GEO数据库](https://www.ncbi.nlm.nih.gov/geo/)筛选小细胞肝癌(HCC)预后相关基因，从[HAGR](https://genomics.senescence.info/)和[msigdb](https://www.gsea-msigdb.org/gsea/msigdb)数据获取衰老相关基因，将两者取交集后构建基因蛋白互作网络

## 3、富集分析与分子验证

### 3.1 分子对接验证

### 3.2 分子动力学验证

### 3.3 关键靶点基因的生存分析
