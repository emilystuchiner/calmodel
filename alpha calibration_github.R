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
caldata<- read_excel("Desktop/Dissertation Chapters/Calibration model/N2O analyzer calibration dataset_github.xlsx", sheet = "standards w indicators_B&C")
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
#Now calibrating the B values for both standards (USGS52 and 500ppm tank) in chunks corresponding to their low, medium, or high range.
#While these values typically align with their corresponding A range (e.g. smallB tends to go with smallA),
#there are occasional discrepancies (e.g. a medB value might correspond with a smallA value)

###########################
#####SMALL B CHUNK#########
Bobs<-caldata$`[14N15N16O] = B` #extracting B values from dataset 
AcorrUSGS<-c(AcorrsmUSGS,AcorrmedUSGS,AcorrlgUSGS) #corrected A values come from the polynomial fit for
#all the USGS-only data
AcorrTANK<-c(AcorrsmTANK,AcorrmedTANK,AcorrlgTANK) #corrected A values come from the polynomial fit for
#all the 500ppm tank-only data

#need to calculate the Btru1 values for USGS and then the Btru2 values for 500ppm tank and then stitch 
#those values back together when fitting m2 so it's a model for the small B chunk, BUT INCLUDING USGS
#AND the 500ppm tank data

#now separating out USGS only data b/c we'll use the known USGS alpha value as our correction 
usgs <- which(caldata$Source == "USGS52 stock")
Bobs1=Bobs[usgs] #B values USGS only 

#now calculating true B values from USGS data 
#we did this by rearranging the eqn: delta=((B/A-Rstd)/Rstd)*1000 b/c if we know the true delta val=13.52
#(as determined by Reston Stable Isotope lab), and we also calculated corrected A values, we can plug in 
#the known delta and the known A and thus solve for the concentration of B
std<-0.0036765 #Rstandard for N2O 15N/14N ratio 
#a two-pool mixing model is used to account for the contaminant in zero-air in standards 
Btru1<-(((13.52/1000)*std)+std)*(AcorrUSGS-0.0426)+(((0/1000)*std)+std)*(0.0426) #this is the rearranged eqn: B = (((13.52/1000)*std)+std)*Acorr, but also
#accounting for zero-air contaminant concentration (0.0426) and assumed zero-air delta value (0)
Btru1 #true B concentrations
plot(Btru1,Bobs1)
#we are log transforming the B values for the same reasons we log transfored A - it improves 
#the polynomial fit when we use a log scale!
logBobs1<-log(Bobs1) 
logBtru1<-log(Btru1)
plot(logBtru1,logBobs1)

#now separating out 500ppm data only b/c we'll need to use a different delta value, the estimated true 
#delta from the 500ppm tank (calculated elsewhere), as our known. 
isstock <- caldata$Source == "500ppm N2O stock"
stock <- which(caldata$Source == "500ppm N2O stock")
Bobs2=Bobs[stock] #B []s 500ppm N2O only - these have not been corrected yet 
logBobs2<-log(Bobs2)   

#now calculating true B values from 500ppm tank data 
tankB<-20.20982 #this was determined previously 
std<-0.0036765
Btru2<-(((tankB/1000)*std)+std)*(AcorrTANK-0.0426)+(((0/1000)*std)+std)*(0.0426)
Btru2 #true B concentrations
plot(Btru2,Bobs2)
#we are log transforming the B values for the same reasons we log transfored A - it improves 
#the polynomial fit when we use a log scale!
logBobs2<-log(Bobs2) 
logBtru2<-log(Btru2)
plot(logBtru2,logBobs2)

#now that we have the true values for the USGS and 500ppm tanks, we will break the B []s into their 
#corresponding low, medium, and high [] "chunks," and for each "chunk" we will fit a new version of m2
#(three chunks so three versions of m2 in total). After all of these models are fit, 
#we will stich the corrected values back together to see if the correction effectively eliminates
#bias in the data.

#fitting m2 for the small B chunk (B < or = 0.005477985) - this step uses the indicators in the 
#Excel file to extract just the small B data, including USGS and 500ppm tank values. 
BobsBoth=Bobs[which(caldata$IndicatorB2=="smallB")]
BtruBoth=c(Btru2,Btru1)
BtruBoth=BtruBoth[which(caldata$IndicatorB2=="smallB")]

logBobsBoth<-log(BobsBoth) 
logBtruBoth<-log(BtruBoth)
plot(logBtruBoth,logBobsBoth)

#m2s = model to correct small B concentrations 
m2s<-lm(logBobsBoth~ logBtruBoth + I(logBtruBoth^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logBtruBoth), m2s$fitted.values[order(logBtruBoth)], lwd = 2, col = "red")
resid<-m2s$residuals
sd(exp(m2s$fitted.values)/BtruBoth)
plot(logBtruBoth,resid)

#this is solving for coefficients for the small B chunk (500ppm tank and USGS data)
summary(m2s) #tells us if x and x^2 are significantly different from zero 
a2<-m2s$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2s$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2s$coefficients[1]-logBobsBoth # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale B concentrations

exp(xhat2) #these are the corrected non-log scale B concentrations 
Bcorrsm<-exp(xhat2)
BcorrsmTANK<-Bcorrsm[c(1:11)]
BcorrsmUSGS<-Bcorrsm[c(12:14)]

#Now need to calculate delta values for the 500ppm tank and USGS corrected data
#500ppm tank deltas
ba2<-BcorrsmTANK/AcorrsmTANK[which(caldata$IndicatorB2=="smallB"&caldata$IndicatorB1=="smallAtank")] 
ba2 
deltaB2<-(((ba2-std))/std)*1000
deltaB2
#this (e.g. the mean, SD, plots of A vs. B) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaB2)
sd(deltaB2)
plot(AcorrsmTANK[which(caldata$IndicatorB2=="smallB"&caldata$IndicatorB1=="smallAtank")],deltaB2)
plot(log(AcorrsmTANK[which(caldata$IndicatorB2=="smallB"&caldata$IndicatorB1=="smallAtank")]),deltaB2)
deltaB2sm<-deltaB2

#USGS deltas
ba1<-BcorrsmUSGS/AcorrsmUSGS[c(1:3)]
ba1 
deltaB1<-(((ba1-std))/std)*1000
deltaB1
#this (e.g. the mean, SD, plots of A vs. B) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaB1)
sd(deltaB1)
plot(AcorrsmUSGS[c(1:3)],deltaB1)
plot(log(AcorrsmUSGS[c(1:3)]),deltaB1)
deltaB1sm<-deltaB1

############################
#####MEDIUM B CHUNK#########

#fitting m2 for the medium B chunk (B < or = 0.1217657)
BobsBoth=Bobs[which(caldata$IndicatorB2=="medB")]
BtruBoth=c(Btru2,Btru1)
BtruBoth=BtruBoth[which(caldata$IndicatorB2=="medB")]

logBobsBoth<-log(BobsBoth) 
logBtruBoth<-log(BtruBoth)
plot(logBtruBoth,logBobsBoth)

#m2m = model to correct medium B concentrations 
m2m<-lm(logBobsBoth~ logBtruBoth + I(logBtruBoth^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logBtruBoth), m2m$fitted.values[order(logBtruBoth)], lwd = 2, col = "red")
resid<-m2m$residuals
sd(exp(m2m$fitted.values)/BtruBoth)
plot(logBtruBoth,resid)

#this is solving for coefficients for the med B chunk (500ppm tank AND USGS)
summary(m2m) #tells us if x and x^2 are significantly different from zero 
a2<-m2m$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2m$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2m$coefficients[1]-logBobsBoth # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale B concentrations

exp(xhat2) #these are the corrected non-log scale B concentrations 
Bcorrmed<-exp(xhat2)
BcorrmedTANK<-Bcorrmed[c(1:45)]
BcorrmedUSGS<-Bcorrmed[c(46:55)] 

#Now need to calculate delta values for the 500ppm tank and USGS corrected data
#500ppm tank deltas
ba2<-BcorrmedTANK/c(AcorrsmTANK[which(caldata$IndicatorB2=="medB"&caldata$IndicatorB1=="smallAtank")],AcorrmedTANK) 
ba2 
deltaB2<-(((ba2-std))/std)*1000
deltaB2
#this (e.g. the mean, SD, plots of A vs. B) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaB2)
sd(deltaB2)
plot(c(AcorrsmTANK[which(caldata$IndicatorB2=="medB"&caldata$IndicatorB1=="smallAtank")],AcorrmedTANK),deltaB2)
plot(log(c(AcorrsmTANK[which(caldata$IndicatorB2=="medB"&caldata$IndicatorB1=="smallAtank")],AcorrmedTANK)),deltaB2)
deltaB2med<-deltaB2

#USGS deltas
ba1<-BcorrmedUSGS/c(AcorrsmUSGS[c(4)],AcorrmedUSGS[c(1:9)]) 
ba1 
deltaB1<-(((ba1-std))/std)*1000
deltaB1
#this (e.g. the mean, SD, plots of A vs. B) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaB1)
sd(deltaB1)
plot(c(AcorrsmUSGS[c(4)],AcorrmedUSGS[c(1:9)]),deltaB1)
plot(log(c(AcorrsmUSGS[c(4)],AcorrmedUSGS[c(1:9)])),deltaB1)
deltaB1med<-deltaB1

###########################
#####LARGE B CHUNK#########

#fitting m2 for the large B chunk (B < or = 1.10295, or really just B > 0.1217657)
BobsBoth=Bobs[which(caldata$IndicatorB2=="lgB")]
BtruBoth=c(Btru2,Btru1)
BtruBoth=BtruBoth[which(caldata$IndicatorB2=="lgB")]

logBobsBoth<-log(BobsBoth) 
logBtruBoth<-log(BtruBoth)
plot(logBtruBoth,logBobsBoth)

#m2l = model to correct large B concentrations 
m2l<-lm(logBobsBoth~ logBtruBoth + I(logBtruBoth^2)) #the squared term gives us the polynomial fit --> y=b0 +bX +bX^2
lines(sort(logBtruBoth), m2l$fitted.values[order(logBtruBoth)], lwd = 2, col = "red")
resid<-m2l$residuals
sd(exp(m2l$fitted.values)/BtruBoth)
plot(logBtruBoth,resid)

#this is solving for coefficients for the large B chunk (500ppm tank AND USGS)
summary(m2l) #tells us if x and x^2 are significantly different from zero 
a2<-m2l$coefficients[3] #a corresponds to bx^2 (because ax^2 in quadratic equation)
b2<-m2l$coefficients[2] #b corresponds to bx (because bx in quadratic equation)
c2<-m2l$coefficients[1]-logBobsBoth # c correspond to b0 (because c in quadratic equation)

#solving quadratic equation 
top <- -b2 + sqrt(b2^2-4*a2*c2)
bot <- 2*a2
xhat2<-top/bot
xhat2#xhat2 is corrected log scale B concentrations

exp(xhat2) #these are the corrected non-log scale B concentrations 
Bcorrlg<-exp(xhat2)
BcorrlgTANK<-Bcorrlg[c(1:16)] 
BcorrlgUSGS<-Bcorrlg[c(17:22)] #need to include the extra large value that's in the medium A bin 

#Now need to calculate delta values for the 500ppm tank and USGS corrected data
#500ppm tank deltas
ba2<-BcorrlgTANK/AcorrlgTANK
ba2 
deltaB2<-(((ba2-std))/std)*1000
deltaB2 #note: We removed values 1 and 10 after calibration b/c outliers
#this (e.g. the mean, SD, plots of A vs. B) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaB2)
sd(deltaB2)
plot(AcorrlgTANK,deltaB2)
plot(log(AcorrlgTANK),deltaB2)
deltaB2lg<-deltaB2

#USGS deltas
ba1<-BcorrlgUSGS/c(AcorrmedUSGS[c(10)],AcorrlgUSGS) 
ba1 
deltaB1<-(((ba1-std))/std)*1000
deltaB1 #note: We removed value 1 after calibration b/c outlier
#this (e.g. the mean, SD, plots of A vs. B) are just diagnostic type info - don't need to look at this, unless of interest. 
mean(deltaB1)
sd(deltaB1)
plot(c(AcorrmedUSGS[c(10)],AcorrlgUSGS),deltaB1)
plot(log(c(AcorrmedUSGS[c(10)],AcorrlgUSGS)),deltaB1)
deltaB1lg<-deltaB1

######################################################################################################
######################################################################################################
######################################################################################################
#now setting up m3 with the three chunked versions of m2 stitched back together to see if any finer-scale
#degree of correction will be needed (this is like a "check")

#center USGS data first
cusgssm<-deltaB1sm-mean(deltaB1sm)
cusgsmed<-deltaB1med-mean(deltaB1med)
cusgslg<-deltaB1lg-mean(deltaB1lg)

cusgs<-c(cusgssm,cusgsmed,cusgslg)

#now centering the 500ppm N2O tank delta B values using the same method

ctanksm<-deltaB2sm-tankB
ctankmed<-deltaB2med-tankB
ctanklg<-deltaB2lg-tankB

ctank<-c(ctanksm,ctankmed,ctanklg)

#now need to make the centered delta B values correspond to their respective corrected log
#A values so I can refit the spline 
#CONCATONATING THE USGS AND TANK DATA  
x1<-c(AcorrTANK,AcorrUSGS) #these are all of the corrected A concentrations for USGS and the 500ppm tank
y<-c(ctank,cusgs) #these are all of the centered B delta values for USGS and the 500ppm tank
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

