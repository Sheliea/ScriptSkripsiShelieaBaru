data=read.delim("clipboard",header=T)
View(data)
str(data)

#install.packages("psych")
library("psych")
b<- corr.test(data[1:5])
b$r
b$t 
b$p 
b$se

sink("correlation.doc")
print(b)
sink()


pairs.panels(data[,-6],pch = 119, stars = T)
detach("package:psych",unload = TRUE)

#install.packages("PerformanceAnalytics")
#require("PerfomarmanceAnalytics")
#chart.Correlation(data[1:8], histrogam = TRUE, pch=100)
#detach("package:PerfomarmanceAnalytics", unload = TRUE)
