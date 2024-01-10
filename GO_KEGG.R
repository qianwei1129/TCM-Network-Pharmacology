#下载包#
install.packages("ggplot2")
install.packages("openxlsx")
#加载包#
library(ggplot2)
library(openxlsx)
#导入数据#
remove(list = ls()) #清除 Global Environment
getwd()  #查看当前工作路径
setwd("C:/Rdata/jc")  #设置需要的工作路径
list.files()  #查看当前工作目录下的文件
go_enrich = read.xlsx("enrich-gene.xlsx",sheet= "ONTOLOGY",sep=',') 
head(go_enrich)

#数据处理#
go_enrich$term <- paste(go_enrich$ID, go_enrich$Description, sep = ': ') #将ID与Description合并成新的一列
go_enrich$term <- factor(go_enrich$term, levels = go_enrich$term,ordered = T) #转成因子，防止重新排列

#纵向柱状图-根据ONTOLOGY类型绘制#
p1 <- ggplot(go_enrich,
             aes(x=term,y=Count, fill=ONTOLOGY)) + #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) +  #柱状图填充颜色
  coord_flip() +  #让柱状图变为纵向
  xlab("GO term") +  #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+  #设置标题
  theme_bw()
p1
#根据ONTOLOGY分类信息添加分组框#
p1+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

#横向柱状图-根据ONTOLOGY类型绘制#
p2 <- ggplot(go_enrich, 
             aes(x=term,y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values = c("#6666FF", "#33CC33", "#FF6666") ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
p2 
#根据ONTOLOGY分类信息添加分组框#
p2 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')

#纵向柱状图——-根据pvalue值绘制#
p3 <- ggplot(go_enrich,aes(y=term,x=Count,fill=pvalue))+  #x、y轴定义；根据pvalue填充颜色
  geom_bar(stat = "identity",width=0.8)+ #柱状图宽度设置
  scale_fill_gradient(low = "red",high ="blue" )+
  labs(title = "GO Terms Enrich",  #设置标题、x轴和Y轴名称
       x = "Gene number", 
       y = "GO Terms")+
  theme(axis.title.x = element_text(face = "bold",size = 16),
        axis.title.y = element_text(face = "bold",size = 16),
        legend.title = element_text(face = "bold",size = 16))+
  theme_bw()
p3
#根据ONTOLOGY分类信息添加分组框#
p3+facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')

#横向柱状图-根据pvalue值绘制#
p4 <- ggplot(go_enrich, 
             aes(x=term,y=Count, fill=pvalue)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_gradient(low = "red",high ="blue" ) + #柱状图填充颜色
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "GO Terms Enrich")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）
p4 
#根据GO富集分析分类信息添加分组框#
p4 + facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')

#气泡图#
P5 <- ggplot(go_enrich,
             aes(y=term,x=Count))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Ratio",y="GO term",title="GO Enrichment")+
  theme_bw()
P5
#根据ONTOLOGY分类信息添加分组框#
P5 + facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')