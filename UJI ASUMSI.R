dataa=read.delim("clipboard",header=TRUE)
View(dataa)
y=dataa$y
x1=dataa$x1
x2=dataa$x2
x3=dataa$x3
x4=dataa$x4
model=lm(y~x1+x2+x3+x4,data = dataa)
model
summary(model)

#Uji Heterokedastisitas
#Nonparametrik:Heterokedastisitas
#P-value > 0.05, maka gagal tolak H0 sehingga tidak ada hetero (parametrik)
#P-value < 0.05, maka tolak H0 sehingga ada hetero (non parametrik)
library(lmtest)
bptest(model)

#uji autokorelasi
#Nonparametrik:Autokorelasi
#Nilai D-W di bawah -2 artinya terdapat autokorelasi positif.
#Nilai D-W di antara -2 sampai +2 artinya tidak ada autokorelasi.
#Nilai D-W di atas +2 artinya terdapat autokorelasi negatif.
dwtest(model) 

#uji multikolnieritas
#Nonparametrik:Non multikolinieritas
#Nilai VIF < 10 maka tidak terjadi multikolinieritas
#Nilai VIF > 10 maka terjadi multikolinieritas
library(car)
vif(model) 

#uji normalitas
#Nonparametrik: Non Normalitas
#karena p-value > 0.05 maka data normal
#karena p-value < 0.05 maka data non normal
shapiro.test(residuals(model)) 
