#A = NNO
#B = N15NO
#C = 15NNO
#D = NN18O

#read in necessary packages for spline fitting 
library("splines")
library("MASS")
library("ConSpline")
library("cgam")

#first correcting the A data
#read in data
library("readxl")
caldata<- read_excel("Desktop/Dissertation Chapters/Calibration model/N2O analyzer calibration dataset_github.xlsx", sheet = "standards w indicators_D")
View(caldata)

#extracting observed N2O concentrations and true N2O concentrations
obs<-caldata$`[14N14N16O] = A` 
tru<-caldata$"TRUEA"
tru<-tru+0.0426 #the 0.0426 accounts for the trace N2O contamination in the zero-grade air

#plotting the true vs. observed [N2O]
plot(tru,obs)
abline(0,1)

#we want to do the polynomial fit on the log scale to better observe the differences in [N2O] and improve the fit 
logobs<-log(obs) 
logtru<-log(tru)
plot(logtru,obs)
plot(logtru,logobs)

#first we are going to correct all the A values over the entire range of standards 
m1<-lm(logobs~logtru+I(logtru^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
summary(m1) ##summary tells us that X and X^2 are both significantly different from zero 
lines(sort(logtru), m1$fitted.values[order(logtru)], lwd = 2, col = "red")
residm1<-m1$residuals #this plot of residuals tells us that model fit is less precise across entire A range
plot(logtru,residm1)

#solving quadratic equation to calculate A  
summary(m1) #quadratic equation formula: ax^2 + bx +c 
a1<-m1$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1$coefficients[1]-logobs

top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the estimated A values

#######################
#######################

#determined that correcting all A data across the entire calibration range is less precise, so we are
#now going to correct the A data in [NNO] concentration ranges, as indicated by plot of residuals from
#correcting all A data across the entire calibration range 

###########################
#####SMALL A CHUNK#########
#correcting A < or = to 1.49ppm N2O 
logobs.4<-logobs[which(logobs<=log(1.5))] #separating out the standards that are less than or equal to 1.49ppm 
logtru.4<-logtru[which(tru <=1.565)] #had to set it at 1.565 of non-logged true so the number of samples
#would match up for logobs vs. logtru (added in 0.0426 for contaminant) 
m1s<-lm(logobs.4~logtru.4+I(logtru.4^2)) #m1s = model 1 small NNO range
summary(m1s) ##summary tells us that X and X^2 are both significantly different from zero 
plot(logtru.4,logobs.4)
lines(sort(logtru.4), m1s$fitted.values[order(logtru.4)], lwd = 2, col = "red")
residm1s<-m1s$residuals #note that plots of residuals for each NNO range are now more homoscedastic 
plot(logtru.4,residm1s)

#solving quadratic equation to calculate A for 0-1.49ppm N2O   
summary(m1s) #quadratic equation formula: ax^2 + bx +c 
a1<-m1s$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1s$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1s$coefficients[1]-logobs.4

#solving the quadratic equation - when A < or = to 1.49ppm
top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the corrected small A values
Acorrsm<-exp(xhat1)
AcorrsmUSGS<-Acorrsm[c(17:20)] #separates USGS A vals
AcorrsmTANK<-Acorrsm[c(1:16)] #separates 500ppm tank A vals

############################
#####MEDIUM A CHUNK#########
#correcting A >1.49ppm and < or = to 33.12ppm N2O 
logobs33<-logobs[which(logobs>log(1.5)&logobs<=log(33.18))] #separating out the standards >1.49 and < or =33.12ppm N2O
logtru33<-log(tru[which(tru>1.565&tru<=33.185)]) #again, needed to tweak values a bit so numbers would match up
m1m<-lm(logobs33~logtru33+I(logtru33^2)) #m1m = model 1 medium NNO range
summary(m1m) ##summary tells us that X and X^2 are both significantly different from zero 
plot(logtru33,logobs33)
lines(sort(logtru33), m1m$fitted.values[order(logtru33)], lwd = 2, col = "red")
residm1m<-m1m$residuals
plot(logtru33,residm1m)

#solving quadratic equation to calculate A for >1.49ppm and < or = to 33.12ppm N2O  
summary(m1m) #quadratic equation formula: ax^2 + bx +c 
a1<-m1m$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1m$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1m$coefficients[1]-logobs33

#solving the quadratic equation when A >1.49ppm and < or = to 33.12ppm N2O 
top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the corrected medium A values
Acorrmed<-exp(xhat1)
AcorrmedUSGS<-Acorrmed[41:50] 
AcorrmedTANK<-Acorrmed[1:40]

##########################
#####LARGE A CHUNK#########
#correcting A >33.12 and < or = to 300ppm N2O 
logobs300<-logobs[which(logobs>log(33.18)&logobs<=max(logobs))] #separating out the standards >33.12ppm N2O
logtru300<-log(tru[which(tru>33.185&tru<=300.065)]) #again, kind of tweaked numbers to account for zero-air contaminant
m1l<-lm(logobs300~logtru300+I(logtru300^2)) #m1l = model 1 large NNO range
summary(m1l) ##summary tells us that X and X^2 are both significantly different from zero 
plot(logtru300,logobs300)
lines(sort(logtru300), m1l$fitted.values[order(logtru300)], lwd = 2, col = "red")
residm1l<-m1l$residuals
plot(logtru300,residm1l)

#solving quadratic equation to calculate A for >33.12 up to 300ppm N2O   
summary(m1l) #quadratic equation formula: ax^2 + bx +c 
a1<-m1l$coefficients[3] #a corresponds the x^2 term in quad. eqn., which is coefficient 4
b1<-m1l$coefficients[2] #b corresponds to bx term in quad. eqn., which is coefficient 3 
c1<-m1l$coefficients[1]-logobs300

#solving the quadratic equation when A >33.12 up to 300ppm N2O 
top <- -b1 + sqrt(b1^2-4*a1*c1)
bot <- 2*a1
xhat1<-top/bot
xhat1 #xhat1 is corrected log scale A values 
exp(xhat1) #exponentiate the logged A values to get the corrected large A values
Acorrlg<-exp(xhat1)
AcorrlgUSGS<-Acorrlg[c(17:21)] 
AcorrlgTANK<-Acorrlg[c(1:16)] 


######################################################################################################
######################################################################################################
######################################################################################################
#Now calibrating the D  values for both standards (USGS52 and 500ppm tank) in chunks corresponding to their low, medium, or high range.
#While these values typically align with their corresponding A range (e.g. smallO  tends to go with smallA),
#there are occasional discrepancies (e.g. a medO  value might correspond with a smallA value)

###########################
#####SMALL D CHUNK#########
Dobs<-caldata$`[14N14N18O] = D` #extracting D values from dataset 
AcorrUSGS<-c(AcorrsmUSGS,AcorrmedUSGS,AcorrlgUSGS) #corrected A values come from the polynomial fit for
#all the USGS-only data
AcorrTANK<-c(AcorrsmTANK,AcorrmedTANK,AcorrlgTANK) #corrected A values come from the polynomial fit for
#all the 500ppm tank-only data

#need to calculate the Dtru1 values for USGS and then the Dtru2 values for 500ppm tank and then stitch 
#those values back together when fitting m2 so it's a model for the small D chunk, BUT INCLUDING USGS
#AND the 500ppm tank data

#now separating out USGS only data b/c we'll use the known USGS alpha value as our correction 
usgs <- which(caldata$Source == "USGS52 stock")
Dobs1=Dobs[usgs] #D values USGS only 

#now calculating true D values from USGS data 
#we did this by rearranging the eqn: delta=((D/A-Rstd)/Rstd)*1000 b/c if we know the true delta val=40.64
#(as determined by Reston Stable Isotope lab), and we also calculated corrected A values, we can plug in 
#the known delta and the known A and thus solve for the concentration of D 
std<-0.0020052 #Rstandard for N2O 15N/14N ratio 
#a two-pool mixing model is used to account for the contaminant in zero-air in standards 
Dtru1<-(((40.64/1000)*std)+std)*(AcorrUSGS-0.0426)+(((0/1000)*std)+std)*(0.0426) #this is the rearranged eqn: D  = (((40.64/1000)*std)+std)*Acorr, but also
#accounting for zero-air contaminant concentration (0.0426) and assumed zero-air delta value (0)
Dtru1 #true D concentrations
plot(Dtru1,Dobs1)
#we are log transforming the D values for the same reasons we log transfored A - it improves 
#the polynomial fit when we use a log scale!
logDobs1<-log(Dobs1) 
logDtru1<-log(Dtru1)
plot(logDtru1,logDobs1)

#now separating out 500ppm data only b/c we'll need to use a different delta value, the estimated true 
#delta from the 500ppm tank (calculated elsewhere), as our known. 
isstock <- caldata$Source == "500ppm N2O stock"
stock <- which(caldata$Source == "500ppm N2O stock")
Dobs2=Dobs[stock] #D [] s 500ppm N2O only - these have not been corrected yet 
logDobs2<-log(Dobs2)   

#now calculating true D values from 500ppm tank data 
tankD<-41.56872 #this was determined previously 
std<-0.0020052
Dtru2<-(((tankD/1000)*std)+std)*(AcorrTANK-0.0426)+(((0/1000)*std)+std)*(0.0426)
Dtru2 #true D concentrations
plot(Dtru2,Dobs2)
#we are log transforming the D values for the same reasons we log transfored A - it improves 
#the polynomial fit when we use a log scale!
logDobs2<-log(Dobs2) 
logDtru2<-log(Dtru2)
plot(logDtru2,logDobs2)

#now that we have the true values for the USGS and 500ppm tanks, we will break the D [] s into their 
#corresponding low, medium, and high [] "chunks," and for each "chunk" we will fit a new version of m2
#(three chunks so three versions of m2 in total). After all of these models are fit, 
#we will stich the corrected values back together to see if the correction effectively eliminates
#bias in the data.

#fitting m2 for the small D chunk (D < or = 0.002987748) - this step uses the indicators in the 
#Excel file to extract just the small D data, including USGS and 500ppm tank values. 
DobsBoth=Dobs[which(caldata$IndicatorO2=="smallO")]
DtruBoth=c(Dtru2,Dtru1)
DtruBoth=DtruBoth[which(caldata$IndicatorO2=="smallO")]

logDobsBoth<-log(DobsBoth) 
logDtruBoth<-log(DtruBoth)
plot(logDtruBoth,logDobsBoth)

#m2s = model to correct small D concentrations 
m2s<-lm(logDobsBoth~ logDtruBoth + I(logDtruBoth^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logDtruBoth), m2s$fitted.values[order(logDtruBoth)], lwd = 2, col = "red")
resid<-m2s$residuals
sd(exp(m2s$fitted.values)/DtruBoth)
plot(logDtruBoth,resid)

#this is solving for coefficients for the small D chunk (500ppm tank and USGS data)
summary(m2s) #tells us if x and x^2 are significantly different from zero 
a2<-m2s$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2s$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2s$coefficients[1]-logDobsBoth # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale D concentrations

exp(xhat2) #these are the corrected non-log scale D concentrations 
Dcorrsm<-exp(xhat2)
DcorrsmTANK<-Dcorrsm[c(1:34)]
DcorrsmUSGS<-Dcorrsm[c(35:39)]

#Now need to calculate delta values for the 500ppm tank and USGS corrected data
#500ppm tank deltas
da2<-DcorrsmTANK/c(AcorrTANK[which(caldata$IndicatorO2=="smallO"&caldata$IndicatorO1=="smallAtank")],AcorrTANK[which(caldata$IndicatorO2=="smallO"&caldata$IndicatorO1=="medAtank")])
da2 
deltaD2<-(((da2-std))/std)*1000
deltaD2 #removed values 3 and 4 after calibration because residuals greater than 100, so assumed they're outliers 
#this (e.g. the mean, SD, plots of A vs. D ) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaD2)
sd(deltaD2)
plot(c(AcorrTANK[which(caldata$IndicatorO2=="smallO"&caldata$IndicatorO1=="smallAtank")],AcorrTANK[which(caldata$IndicatorO2=="smallO"&caldata$IndicatorO1=="medAtank")]),deltaD2)
plot(log(c(AcorrTANK[which(caldata$IndicatorO2=="smallO"&caldata$IndicatorO1=="smallAtank")],AcorrTANK[which(caldata$IndicatorO2=="smallO"&caldata$IndicatorO1=="medAtank")])),deltaD2)
deltaD2sm<-deltaD2

#USGS deltas
da1<-DcorrsmUSGS/c(AcorrsmUSGS,AcorrmedUSGS[c(1)])
da1 
deltaD1<-(((da1-std))/std)*1000
deltaD1
#this (e.g. the mean, SD, plots of A vs. D ) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaD1)
sd(deltaD1)
plot(c(AcorrsmUSGS,AcorrmedUSGS[c(1)]),deltaD1)
plot(log(c(AcorrsmUSGS,AcorrmedUSGS[c(1)])),deltaD1)
deltaD1sm<-deltaD1

############################
#####MEDIUM D CHUNK#########

#fitting m2 for the medium D chunk (D > 0.002987748 and < or = 0.066412224)
DobsBoth=Dobs[which(caldata$IndicatorO2=="medO")]
DtruBoth=c(Dtru2,Dtru1)
DtruBoth=DtruBoth[which(caldata$IndicatorO2=="medO")]

logDobsBoth<-log(DobsBoth) 
logDtruBoth<-log(DtruBoth)
plot(logDtruBoth,logDobsBoth)

#m2m = model to correct medium D concentrations 
m2m<-lm(logDobsBoth~ logDtruBoth + I(logDtruBoth^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logDtruBoth), m2m$fitted.values[order(logDtruBoth)], lwd = 2, col = "red")
resid<-m2m$residuals
sd(exp(m2m$fitted.values)/DtruBoth)
plot(logDtruBoth,resid)

#this is solving for coefficients for the med D chunk (500ppm tank AND USGS)
summary(m2m) #tells us if x and x^2 are significantly different from zero 
a2<-m2m$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2m$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2m$coefficients[1]-logDobsBoth # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale D concentrations

exp(xhat2) #these are the corrected non-log scale D concentrations 
Dcorrmed<-exp(xhat2)
DcorrmedTANK<-Dcorrmed[c(1:25)]
DcorrmedUSGS<-Dcorrmed[c(26:36)]


#Now need to calculate delta values for the 500ppm tank and USGS corrected data
#500ppm tank deltas
da2<-DcorrmedTANK/c(AcorrTANK[which(caldata$IndicatorO2=="medO"&caldata$IndicatorO1=="medAtank")],AcorrTANK[which(caldata$IndicatorO2=="medO"&caldata$IndicatorO1=="lgAtank")]) 
da2 
deltaD2<-(((da2-std))/std)*1000
deltaD2
#this (e.g. the mean, SD, plots of A vs. D ) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaD2)
sd(deltaD2)
plot(c(AcorrTANK[which(caldata$IndicatorO2=="medO"&caldata$IndicatorO1=="medAtank")],AcorrTANK[which(caldata$IndicatorO2=="medO"&caldata$IndicatorO1=="lgAtank")]),deltaD2)
plot(log(c(AcorrTANK[which(caldata$IndicatorO2=="medO"&caldata$IndicatorO1=="medAtank")],AcorrTANK[which(caldata$IndicatorO2=="medO"&caldata$IndicatorO1=="lgAtank")])),deltaD2)
deltaD2med<-deltaD2

#USGS deltas
da1<-DcorrmedUSGS/c(AcorrmedUSGS[c(2:10)],AcorrlgUSGS[c(1:2)])
da1 
deltaD1<-(((da1-std))/std)*1000
deltaD1
#this (e.g. the mean, SD, plots of A vs. D ) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaD1)
sd(deltaD1)
plot(c(AcorrmedUSGS[c(2:10)],AcorrlgUSGS[c(1:2)]),deltaD1)
plot(log(c(AcorrmedUSGS[c(2:10)],AcorrlgUSGS[c(1:2)])),deltaD1)
deltaD1med<-deltaD1

###########################
#####LARGE D CHUNK#########

#fitting m2 for the large D chunk (D  < or = 0.601560000, or really just D  > 0.066412224)
DobsBoth=Dobs[which(caldata$IndicatorO2=="lgO")]
DtruBoth=c(Dtru2,Dtru1)
DtruBoth=DtruBoth[which(caldata$IndicatorO2=="lgO")]

logDobsBoth<-log(DobsBoth) 
logDtruBoth<-log(DtruBoth)
plot(logDtruBoth,logDobsBoth)

#m2l = model to correct large D concentrations 
m2l<-lm(logDobsBoth~ logDtruBoth + I(logDtruBoth^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logDtruBoth), m2l$fitted.values[order(logDtruBoth)], lwd = 2, col = "red")
resid<-m2l$residuals
sd(exp(m2l$fitted.values)/DtruBoth)
plot(logDtruBoth,resid)

#this is solving for coefficients for the large D chunk (500ppm tank AND USGS)
summary(m2l) #tells us if x and x^2 are significantly different from zero 
a2<-m2l$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2l$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2l$coefficients[1]-logDobsBoth # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale D concentrations

exp(xhat2) #these are the corrected non-log scale D concentrations 
Dcorrlg<-exp(xhat2)
DcorrlgTANK<-Dcorrlg[c(1:13)]
DcorrlgUSGS<-Dcorrlg[c(14:16)] #to include the extra large value that's in the medium A bin 

#Now need to calculate delta values for the 500ppm tank and USGS corrected data
#500ppm tank deltas
da2<-DcorrlgTANK/AcorrTANK[which(caldata$IndicatorO2=="lgO"&caldata$IndicatorO1=="lgAtank")]
da2 
deltaD2<-(((da2-std))/std)*1000
deltaD2 #removed value 7 after calibration because residual greater than 100, so assumed it's an outlier
#this (e.g. the mean, SD, plots of A vs. D ) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaD2)
sd(deltaD2)
plot(AcorrlgTANK,deltaD2)
plot(log(AcorrlgTANK),deltaD2)
deltaD2lg<-deltaD2

#USGS deltas
da1<-DcorrlgUSGS/AcorrlgUSGS[c(3:5)]
da1 
deltaD1<-(((da1-std))/std)*1000
deltaD1
#this (e.g. the mean, SD, plots of A vs. D ) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaD1)
sd(deltaD1)
plot(c(AcorrlgUSGS[c(3:5)]),deltaD1)
plot(log(c(AcorrlgUSGS[c(3:5)])),deltaD1)
deltaD1lg<-deltaD1

######################################################################################################
######################################################################################################
######################################################################################################
#now setting up m3 with the three chunked versions of m2 stitched back together to see if any finer-scale
#degree of correction will be needed (this is like a "check")

#center USGS data first
cusgssm<-deltaD1sm-mean(deltaD1sm)
cusgsmed<-deltaD1med-mean(deltaD1med)
cusgslg<-deltaD1lg-mean(deltaD1lg)

cusgs<-c(cusgssm,cusgsmed,cusgslg)

#now centering the 500ppm N2O tank delta D values using the same method

ctanksm<-deltaD2sm-tankD
ctankmed<-deltaD2med-tankD
ctanklg<-deltaD2lg-tankD

ctank<-c(ctanksm,ctankmed,ctanklg)

#now need to make the centered delta D values correspond to their respective corrected log
#A values so I can refit the spline 
#CONCATONATING THE USGS AND TANK DATA  
x1<-c(AcorrTANK,AcorrUSGS) #these are all of the corrected A concentrations for USGS and the 500ppm tank
y<-c(ctank,cusgs) #these are all of the centered D delta values for USGS and the 500ppm tank
x2<-c(caldata$Indicator)  #this is an indicator term, 1=if tank, 0=if USGS

logx1<- log(x1) #need to log scale the "coarse" corrected A values 

myknots = quantile(logx1, probs = c(0.33, 0.66)) #per Zach's suggestion 
m3 <- lm( y ~ bs(logx1, knots = myknots))
summary(m3)
plot(log(AcorrTANK),ctank)
points(log(AcorrUSGS),cusgs, col="red")
lines(sort(logx1), m3$fitted.values[order(logx1)], lwd = 2, col = "red") #the flat line indicates that 
#there is no residual bias and additional correction is not necessary

#plotting residuals of centered data vs. size (aka corrected A) and vs. pressure 
residm5<-m3$residuals
plot(logx1,residm5) #this is residuals vs. the size effect (aka corrected A)
points(log(AcorrUSGS),cusgs, col="red")
lines(sort(logx1), m3$fitted.values[order(logx1)], lwd = 2, col = "red") #the flat line indicates that 
#there is no residual bias and additional correction is not necessary

#note that the tankD value calculated here differs slightly from the tankD value written in
#Stuchiner et al. (2020). This is because the final value used in our manuscript included ~20
#additional standard runs that we added into the dataset to compute the final, mean, d18O 
#value. The methodology for estimating the value still stands exactly the same. 
