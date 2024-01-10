# TCM-Network-Pharmacology

***中药网络药理学全过程***

## 网站和工具

### 成分/靶点获取网站

[TCMSP](https://old.tcmsp-e.com/tcmsp.php) [TCMBank](https://www.tcmbank.cn/) [HERB](http://herb.ac.cn/) [swisstarget](http://www.swisstargetprediction.ch/)

### 工具

[R](https://www.r-project.org/) [cytoscape](https://cytoscape.org/)

## 黄芪六一汤靶点的获取

使用[TCMSP数据库](https://old.tcmsp-e.com/tcmsp.php)，选择检索种类为"Herb name"，检索黄芪。然后根据公认标准筛选口服生物利用度（OB≥30%）和药物相似性（DL≥0.18），共获得黄芪的20个作用化合物(\data\pre\chemicals\_hq.xlsx)和甘草的69个化合物。 然后根据smile值，使用[swisstarget数据库](http://www.swisstargetprediction.ch/)获取化合物作用靶点，其中部分化合物没有smile值则使用[TCMBank数据库](https://www.tcmbank.cn/)，h或者[HERB](http://herb.ac.cn/)获取作用靶点
