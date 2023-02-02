
##微生物网络鲁棒性评估
library(igraph)
library(ggplot2)
library(ggsci)

setwd("C:/BaiduNetdiskWorkspace/aquatic_microbial_ecology/network_imeta/robustness")


#读取网络邻接矩阵，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adj1 <- read.csv('n1_matrix0_1.csv', row.names = 1, check.names = FALSE)
adj2 <- read.csv('n2_matrix0_1.csv', row.names = 1, check.names = FALSE)
adj3 <- read.csv('n3_matrix0_1.csv', row.names = 1, check.names = FALSE)
adj4 <- read.csv('n4_matrix0_1.csv', row.names = 1, check.names = FALSE)
adj5 <- read.csv('n5_matrix0_1.csv', row.names = 1, check.names = FALSE)
adj6 <- read.csv('n6_matrix0_1.csv', row.names = 1, check.names = FALSE)






#计算自然连通度
#adj_matrix 是网络邻接矩阵
nc <- function(adj_matrix) {
  #获取 0-1 矩阵，1 表示节点间存在边，0 表示不存在边
  adj_matrix <- as.matrix(adj_matrix)
  adj_matrix[abs(adj_matrix) != 0] <- 1
  
  #矩阵的特征分解，获取特征值 λ
  lambda <- eigen(adj_matrix, only.values = TRUE)$values
  lambda <- sort(lambda, decreasing = TRUE)
  
  #计算“平均特征根”，获得自然连通度
  lambda_sum <- 0
  N = length(lambda)
  for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
  lambda_average <- log(lambda_sum/N, base = exp(1))
  lambda_average
}






#计算自然连通度
natural_connectivity1 <- nc(adj1)
natural_connectivity2 <- nc(adj2)
natural_connectivity3 <- nc(adj3)
natural_connectivity4 <- nc(adj4)
natural_connectivity5 <- nc(adj5)
natural_connectivity6 <- nc(adj6)




#转化为 igraph 邻接列表，计算节点平均度
g1 <- graph_from_adjacency_matrix(as.matrix(adj1), mode = 'undirected', diag = FALSE)
g2 <- graph_from_adjacency_matrix(as.matrix(adj2), mode = 'undirected', diag = FALSE)
g3 <- graph_from_adjacency_matrix(as.matrix(adj3), mode = 'undirected', diag = FALSE)
g4 <- graph_from_adjacency_matrix(as.matrix(adj4), mode = 'undirected', diag = FALSE)
g5 <- graph_from_adjacency_matrix(as.matrix(adj5), mode = 'undirected', diag = FALSE)
g6 <- graph_from_adjacency_matrix(as.matrix(adj6), mode = 'undirected', diag = FALSE)





average_degree1 <- mean(degree(g1))
average_degree2 <- mean(degree(g2))
average_degree3 <- mean(degree(g3))
average_degree4 <- mean(degree(g4))
average_degree5 <- mean(degree(g5))
average_degree6 <- mean(degree(g6))





for (i in 1:(nrow(adj1)-1)) {
  
  #在邻接矩阵中随机移除 i 个节点
  remove_node1 <- sample(1:nrow(adj1), i)
  adj1_remove <- adj1[-remove_node1,-remove_node1]
  
  #计算自然连通度
  natural_connectivity1 <- c(natural_connectivity1, nc(adj1_remove))
  
  #计算节点平均度
  g1 <- graph_from_adjacency_matrix(as.matrix(adj1_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree1 <- c(average_degree1, mean(degree(g1)))
  
}





for (i in 1:(nrow(adj2)-1)) {

  #在邻接矩阵中随机移除 i 个节点
  remove_node2 <- sample(1:nrow(adj2), i)
  adj2_remove <- adj2[-remove_node2,-remove_node2]

  #计算自然连通度
  natural_connectivity2 <- c(natural_connectivity2, nc(adj2_remove))

  #计算节点平均度
  g2 <- graph_from_adjacency_matrix(as.matrix(adj2_remove), mode = 'undirected', diag = FALSE)


  average_degree2 <- c(average_degree2, mean(degree(g2)))

}









for (i in 1:(nrow(adj3)-1)) {
  
  #在邻接矩阵中随机移除 i 个节点
  remove_node3 <- sample(1:nrow(adj3), i)
  adj3_remove <- adj3[-remove_node3,-remove_node3]
  
  #计算自然连通度
  natural_connectivity3 <- c(natural_connectivity3, nc(adj3_remove))
  
  #计算节点平均度
  g3 <- graph_from_adjacency_matrix(as.matrix(adj3_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree3 <- c(average_degree3, mean(degree(g3)))
  
}







for (i in 1:(nrow(adj4)-1)) {
  
  #在邻接矩阵中随机移除 i 个节点
  remove_node4 <- sample(1:nrow(adj4), i)
  adj4_remove <- adj4[-remove_node4,-remove_node4]
  
  #计算自然连通度
  natural_connectivity4 <- c(natural_connectivity4, nc(adj4_remove))
  
  #计算节点平均度
  g4 <- graph_from_adjacency_matrix(as.matrix(adj4_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree4 <- c(average_degree4, mean(degree(g4)))
  
}








for (i in 1:(nrow(adj5)-1)) {
  
  #在邻接矩阵中随机移除 i 个节点
  remove_node5 <- sample(1:nrow(adj5), i)
  adj5_remove <- adj5[-remove_node5,-remove_node5]
  
  #计算自然连通度
  natural_connectivity5 <- c(natural_connectivity5, nc(adj5_remove))
  
  #计算节点平均度
  g5 <- graph_from_adjacency_matrix(as.matrix(adj5_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree5 <- c(average_degree5, mean(degree(g5)))
  
}







for (i in 1:(nrow(adj6)-1)) {
  
  #在邻接矩阵中随机移除 i 个节点
  remove_node6 <- sample(1:nrow(adj6), i)
  adj6_remove <- adj6[-remove_node6,-remove_node6]
  
  #计算自然连通度
  natural_connectivity6 <- c(natural_connectivity6, nc(adj6_remove))
  
  #计算节点平均度
  g6 <- graph_from_adjacency_matrix(as.matrix(adj6_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree6 <- c(average_degree6, mean(degree(g6)))
  
}







####-------------随机去除一半节点






  #在邻接矩阵中随机移除 i 个节点
  remove_node1 <- sample(1:nrow(adj1), 0.1*nrow(adj1))
  adj1_remove <- adj1[-remove_node1,-remove_node1]
  
  #计算自然连通度
  natural_connectivity1 <- c(natural_connectivity1, nc(adj1_remove))
  
  #计算节点平均度
  g1 <- graph_from_adjacency_matrix(as.matrix(adj1_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree1 <- c(average_degree1, mean(degree(g1)))
  






  
  #在邻接矩阵中随机移除 i 个节点
  remove_node2 <- sample(1:nrow(adj2),  0.1*nrow(adj2))
  adj2_remove <- adj2[-remove_node2,-remove_node2]
  
  #计算自然连通度
  natural_connectivity2 <- c(natural_connectivity2, nc(adj2_remove))
  
  #计算节点平均度
  g2 <- graph_from_adjacency_matrix(as.matrix(adj2_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree2 <- c(average_degree2, mean(degree(g2)))
  










  #在邻接矩阵中随机移除 i 个节点
  remove_node3 <- sample(1:nrow(adj3),  0.1*nrow(adj3))
  adj3_remove <- adj3[-remove_node3,-remove_node3]
  
  #计算自然连通度
  natural_connectivity3 <- c(natural_connectivity3, nc(adj3_remove))
  
  #计算节点平均度
  g3 <- graph_from_adjacency_matrix(as.matrix(adj3_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree3 <- c(average_degree3, mean(degree(g3)))
  








  
  #在邻接矩阵中随机移除 i 个节点
  remove_node4 <- sample(1:nrow(adj4),  0.1*nrow(adj4))
  adj4_remove <- adj4[-remove_node4,-remove_node4]
  
  #计算自然连通度
  natural_connectivity4 <- c(natural_connectivity4, nc(adj4_remove))
  
  #计算节点平均度
  g4 <- graph_from_adjacency_matrix(as.matrix(adj4_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree4 <- c(average_degree4, mean(degree(g4)))
  









  
  #在邻接矩阵中随机移除 i 个节点
  remove_node5 <- sample(1:nrow(adj5), 0.1*nrow(adj5))
  adj5_remove <- adj5[-remove_node5,-remove_node5]
  
  #计算自然连通度
  natural_connectivity5 <- c(natural_connectivity5, nc(adj5_remove))
  
  #计算节点平均度
  g5 <- graph_from_adjacency_matrix(as.matrix(adj5_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree5 <- c(average_degree5, mean(degree(g5)))
  








  
  #在邻接矩阵中随机移除 i 个节点
  remove_node6 <- sample(1:nrow(adj6), 0.1*nrow(adj6))
  adj6_remove <- adj6[-remove_node6,-remove_node6]
  
  #计算自然连通度
  natural_connectivity6 <- c(natural_connectivity6, nc(adj6_remove))
  
  #计算节点平均度
  g6 <- graph_from_adjacency_matrix(as.matrix(adj6_remove), mode = 'undirected', diag = FALSE)
  
  
  average_degree6 <- c(average_degree6, mean(degree(g6)))
  



dat1 <- data.frame(remove_node = rep(c(0, 50, 40, 30, 20, 10), 2),
                     variable = c(rep('natural_connectivity', 6), rep('average_degree', 6)),
                     values = c(natural_connectivity1, average_degree1))


dat2 <- data.frame(remove_node = rep(c(0, 50, 40, 30, 20, 10), 2),
                   variable = c(rep('natural_connectivity', 6), rep('average_degree', 6)),
                   values = c(natural_connectivity2, average_degree2))



dat3 <- data.frame(remove_node = rep(c(0, 50, 40, 30, 20, 10), 2),
                   variable = c(rep('natural_connectivity', 6), rep('average_degree', 6)),
                   values = c(natural_connectivity3, average_degree3))




dat4 <- data.frame(remove_node = rep(c(0, 50, 40, 30, 20, 10), 2),
                   variable = c(rep('natural_connectivity', 6), rep('average_degree', 6)),
                   values = c(natural_connectivity4, average_degree4))




dat5 <- data.frame(remove_node = rep(c(0, 50, 40, 30, 20, 10), 2),
                   variable = c(rep('natural_connectivity', 6), rep('average_degree', 6)),
                   values = c(natural_connectivity5, average_degree5))



dat6 <- data.frame(remove_node = rep(c(0, 50, 40, 30, 20, 10), 2),
                   variable = c(rep('natural_connectivity', 6), rep('average_degree', 6)),
                   values = c(natural_connectivity6, average_degree6))


write.csv(dat1, 'dat11.csv', row.names = FALSE, quote = FALSE)
write.csv(dat2, 'dat21.csv', row.names = FALSE, quote = FALSE)
write.csv(dat3, 'dat31.csv', row.names = FALSE, quote = FALSE)
write.csv(dat4, 'dat41.csv', row.names = FALSE, quote = FALSE)
write.csv(dat5, 'dat51.csv', row.names = FALSE, quote = FALSE)
write.csv(dat6, 'dat61.csv', row.names = FALSE, quote = FALSE)


##########ggplot2 作图，拟合线--------------------------------------------------------------------------------



## 1. 抗毁性



dat1 <- data.frame(remove_node = rep(0:(nrow(adj1)-1), 2),
                  variable = c(rep('natural_connectivity', nrow(adj1)), rep('average_degree', nrow(adj1))),
                  values = c(natural_connectivity1, average_degree1))



dat2 <- data.frame(remove_node = rep(0:(nrow(adj2)-1), 2),
                   variable = c(rep('natural_connectivity', nrow(adj2)), rep('average_degree', nrow(adj2))),
                   values = c(natural_connectivity2, average_degree2))




dat3 <- data.frame(remove_node = rep(0:(nrow(adj3)-1), 2),
                   variable = c(rep('natural_connectivity', nrow(adj3)), rep('average_degree', nrow(adj3))),
                   values = c(natural_connectivity3, average_degree3))



dat4 <- data.frame(remove_node = rep(0:(nrow(adj4)-1), 2),
                   variable = c(rep('natural_connectivity', nrow(adj4)), rep('average_degree', nrow(adj4))),
                   values = c(natural_connectivity4, average_degree4))




dat5 <- data.frame(remove_node = rep(0:(nrow(adj5)-1), 2),
                   variable = c(rep('natural_connectivity', nrow(adj5)), rep('average_degree', nrow(adj5))),
                   values = c(natural_connectivity5, average_degree5))




dat6 <- data.frame(remove_node = rep(0:(nrow(adj6)-1), 2),
                   variable = c(rep('natural_connectivity', nrow(adj6)), rep('average_degree', nrow(adj6))),
                   values = c(natural_connectivity6, average_degree6))



# 
# write.csv(dat1, 'dat1.csv', row.names = FALSE, quote = FALSE)
# write.csv(dat2, 'dat2.csv', row.names = FALSE, quote = FALSE)
# write.csv(dat3, 'dat3.csv', row.names = FALSE, quote = FALSE)
# write.csv(dat4, 'dat4.csv', row.names = FALSE, quote = FALSE)
# write.csv(dat5, 'dat5.csv', row.names = FALSE, quote = FALSE)
# write.csv(dat6, 'dat6.csv', row.names = FALSE, quote = FALSE)





dat <- read.csv('dat.csv', header = T)

p1 = ggplot(dat, aes(x=nodes_removal_proportion, y=values, color = Salinity)) +
  geom_point(alpha = 0.6, size = 0.3) +
  scale_color_gradientn(colours = sal_col) +
  theme_light() + 
  facet_wrap(~variable, ncol = 1, scale = 'free') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.6,0.6))

p1





###  2.  linktype 百分比堆积柱状图

setwd("C:/Desktop/aquatic_microbiome/network")


data = read.csv("linktype.csv")

p2 = ggplot(data, aes(variable, value)) +
  
  geom_bar(aes(fill = linktype), stat = 'identity', position = 'fill', alpha = 1) +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = c("#B24745FF", "#0073C2FF")) +
  
  facet_wrap(~class, ncol = 1, scale = "free") + 
  
  theme_light() + guides(fill=F) +
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 

p2


###  3.  topological parameter

dat <- read.csv("indexes.csv", header = TRUE)


p3 <- ggplot(dat, aes(x= habitat, y= value))+ 
  
  geom_point(shape = 19, size = 3, alpha = 1, color = "#0073C2FF") + 
  
  scale_fill_npg(alpha = 1) + 
  
  theme_light() + theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(), 
                        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  
  facet_wrap(~index, scales = "free_y", ncol = 1)



p3








ggsave(plot = p1, "robuss.pdf", width = 60, height = 110, units = c("mm"), dpi = 300)

ggsave(plot = p2, "linktype.pdf", width = 60, height = 105, units = c("mm"), dpi = 300)

ggsave(plot = p3, "index.pdf", width = 60, height = 215, units = c("mm"), dpi = 300)


