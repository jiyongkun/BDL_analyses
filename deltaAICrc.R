#BDL_analyses.R
#Licence: GPL-2
#请将代码保存为扩展名为.R的文本文件
#并将其与输入树文件放在同一文件夹下
#运行该脚本程序前，必须先安装以下程序包：
	# Laser
	# ape
	# TreePar
	# TreeSim
	# subplex

#载入程序包
library("TreePar")	

#读取修剪和去除外群后的分歧时间树文件，nexus格式
tree <- read.nexus(file = "DivTime_NoOG.tre")

#估计输入树的实际物种形成和灭绝速率
#估算过程考虑了不完全类群采样的影响（0.76） 
tree.times <- branching.times(tree)
bd.s.o <- bd.shifts.optim(tree.times, sampling = c(0.76,1),1,0,0)
param <- bd.s.o[[2]][[1]]
k <- param[2] #extinction/speciation
d <- param[3] #diversification (speciation-extinction)
lambda <- d/(1-k) #获得物种形成速率
mu <- (d*k)/(1-k) #获得物种灭绝速率
cat("lambda = ", lambda, "\n")
cat("mu = ", mu, "\n")

#模拟加入缺少的物种
#模拟重复进行n次形成混合真实数据集		
n=1000
cat( "simulating missing species...\n")
for (i in 1:n){
    #模拟从该属的根节点时间开始
	xsim <- corsim(x=tree.times, lambda, mu, missing=20, told=6.5, tyoung=0)
    #对6个速率模型进行检验
    result <- fitdAICrc(x=xsim, modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate" ,"yule3rate"), ints=100)
    #由于fitdAICrc函数有BUG，yule3rate模型结果中只有st2，没有st1
    #因此单独对yule3rate模型进行检验
    result2 <- yule3rate(x=xsim, ints=100)
    #将结果写入逗号分隔的CSV文件
    #该文件可用EXCEL打开
    write.table(result, file="semi_empirical.csv", append=TRUE, sep=",", col.names=NA, qmethod="double")
    write.table(result2, file="semi_empirical_yule3rate.csv", append=TRUE, sep=",", col.names=NA, qmethod="double")
	cat(i, "\n")	
	}

#模拟n个速率恒定的时间树
#组成零分布数据集
n=1000
cat( "simulating trees...\n")
xsim <- sim.bd.taxa(n=85, numbsim = n, lambda, mu, frac = 1, complete = FALSE, stochsampling = FALSE)
cat( "fitdAICrc...\n")
#对6个速率模型进行检验
for (i in 1:n){
    result <- fitdAICrc(x=branching.times(xsim[[i]]), modelset=c("pureBirth", "bd", "DDX", "DDL", "yule2rate" ,"yule3rate"), ints=100)
    #结果写入CSV文件
    write.table(result, file = "NullDistribution.csv", append = TRUE, sep = ",", col.names = NA, qmethod = "double")
	cat(i, "\n")	
	}
cat("Finished!")
   


