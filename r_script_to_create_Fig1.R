#R script used to create Figure 1 in ES Nano Praetorius et al. manuscript "Strategies for determining heteroaggregation attachment efficiencies of engineered nanoparticles in aquatic environments"

#this code is used to evaluate the relative importance of the 3 concurrent mechanisms (perikinetic, orthokinetic and differential sedimentation) 
#contributing to the overall collision rate between 2 aggregating particles of different sizes and densities. The size and density of particle i are fixed
#while size and density of particle j vary within a specified range. The fixed values and the ranges can be varied to explore different scenarios. 
#The values chosen here represent a metal oxide type engineered nanoparticle (i) heteroaggregating with a natural suspended matter (SPM) paticle (j)

#the collision rate constant for heteroaggregation is calculated according to Equation 2 in the paper following classical colloid theory
#This represents the most simple case of rectilinear collision mechanisms and Stokes Law is used to calculate settling
#In theory, this code can be extended to include curvilinear (or other) correction terms and/or to use different settling equations, 
#but these options are not implemented in the current code



### part 1: calculate collision rate constants


#packages required (install once before using code)
#install.packages("dplyr")    
#install.packages("plotly")

library(dplyr)    # alternatively, this also loads %>%
library(plotly)

#define constants
kb = 1.38064852e-23 #Boltzmann constant (in m2 kg s-2 K-1)
pi = 3.14159265
g = 9.80665  #gravitational acceleration (in m s-2)

#define variables (can be modified to simulate different conditions)
Tw = 293.0 #water temperature (in K)
eta = 1.002e-3 #viscosity (of water) at 293K (in kg m-1 s-1)
rhow = 1000.0 #density of water (in kg m-3) 
G = 10.0 #shear rate (in s-1)

#define particle size and density of fixed particle i (can be modified to explore different cases)
ri = 10e-9 #radius of particle i (in m)
rhoi = 4000.0 #density of particle i (in kg m-3) 


df1 = data.frame(rhoi = numeric(), rj = numeric(), kcollperi=numeric(), kcollshear = numeric(), kcollsed=numeric())
#dataframe to store inputs (densiy, size) and outputs (collision rate contants for perikinetic, orthokinetic and differential settling)


#calculate collision rate constants for each combination of size and density of particle j
#for (i in seq(1,4001, by = 10 )) sets loops through densities rhoj of 1000 to 5000 (in increments of 10).
#Note: smaller increments lead to better resolution but significantly longer computation times
for (i in seq(1,4001, by = 10 )) { 
  rhoj = 999.0 + i

  j=1.0
  #second loop (while loop) loops through sizes rj )from 1 nm to 10000 nm (10 micrometer)
  while(j<10001){
    rj=round((j)*1e-9,10) #variable radius of particle j in m (rounded to avoid numerical errors)
    
    #collision rate calculations according to Eq. 2 in paper
    kcollperi = 2*kb*Tw*(ri+rj)^2 / (3*eta*ri*rj) #perikinetic
    kcollshear = 4*G*(ri+rj)^3 / 3 #orthokinetic
    kcollsed = pi*(ri+rj)^2 * abs((2*g/(9*eta))*((rhoj-rhow)*(rj^2)-(rhoi-rhow)*(ri^2))) #differential sedimentation
    
    #store data from each run in dataframe
    z <- data.frame(rhoj, rj, kcollperi, kcollshear, kcollsed )
    df1 <- rbind(df1, z) 
    
    #note: here the size of the increments increases by a factor of 10, for each "order of magnitude bracket" 
    #this is done to achieve higher resolution at smaller sizes. Can be modified when looking at different size brackets
    if(j<10){
      j=j+0.1
    } else if (j<100) {
      j=j+1.0
    } else if (j<1000) {
      j=j+10.0
    }else { #more than 1000
      j=j+100.0
    }
  }
}

#write results into file: kcollspecific.dat
write.table(df1, "kcollspecific.dat", sep = " ", append = T, row.names=F, col.names =F )


### part 2: plot results in contour plots 

#read data from results file
kcoll_data <- read.table("kcollspecific.dat", header=FALSE, skip=1)

kcoll_data[,2] = kcoll_data[,2]*10^(6) #transform radius from m to micrometer
kcoll_data[,6] = kcoll_data[,3] + kcoll_data[,4] + kcoll_data[,5] 
#sum individual collision rate constants to obtain overall collision rate constants

#kcollsed is expected to contain 0 values (e.g. when particle i and j have same size and density)
#these will be plotted in a separate line. Since they cant be plotted in a log plot, 0 values are resaved as 10e-100
for (row in 1:nrow(kcoll_data)){
  if(kcoll_data[row,5]<10e-100){
    kcoll_data[row,5]=10e-100
}
}

#prepare graph properties
f1 <- list(
  size = 20
)

f2 <- list(
  size = 18
)

#x-axis
ax <- list(
  autotick = FALSE,
  ticks = "outside",
  ticklen = 5,
  tickwidth = 2,
  tickcolor = toRGB("black"),
  tickfont = f2,
  title = "radius particle <i>j</i> (&mu;m)",
  titlefont = f1,
  type = "log"
)

#y-axis
ay <- list(
  autotick = FALSE,
  ticks = "outside",
  tick0 = 1000,
  dtick = 1000,
  ticklen = 5,
  tickwidth = 2,
  tickcolor = toRGB("black"),
  tickfont = f2,
  title = "density particle <i>j</i> (kg/m<sup>3</sup>)",
  titlefont = f1
)


#plot perikinetic rate constants
p1 <- plot_ly(
  x = kcoll_data[,2], y = kcoll_data[,1], z = log10(kcoll_data[,3]),
  type = "contour", autocontour = F) %>% 
  colorbar(title = "log(rate \n constant \n in m<sup>3</sup>/s)", exponentformat = "E",   titlefont = f1, tickfont = f2
) %>%
  layout(
    title = "perikinetic aggregation",
    titlefont = f1,
    xaxis = ax,
    yaxis = ay 
  )
p1

#plot orthokinetic rate constants
p2 <- plot_ly(
  x = kcoll_data[,2], y = kcoll_data[,1], z = log10(kcoll_data[,4]),
  type = "contour")%>% 
  colorbar(title = "log(rate \n constant \n in m<sup>3</sup>/s)", exponentformat = "E",   titlefont = f1, tickfont = f2
  ) %>%
  layout(
    title = "orthokinetic aggregation",
    title = f1,
    xaxis = ax,
    yaxis = ay 
  )
p2

#plot differntial sedimentation rate constants without 0 values
p3 <- plot_ly(
  x = kcoll_data[,2], y = kcoll_data[,1], z = log10(kcoll_data[,5]),
  type = "contour", autocontour = F,   contours = list(
    start = -27,
    end = -14,
    size = 1
  )) %>% 
  colorbar(title = "log(rate \n constant \n in m<sup>3</sup>/s)", exponentformat = "E",   titlefont = f1, tickfont = f2
  ) %>% 
  layout(
    title = "differential sedimentation",
    title = f1,
    xaxis = ax,
    yaxis = ay 
  )
p3

#add trace of 0 values to differential sedimentation plot
df2 = data.frame(x_0 = numeric(), y_0 = numeric()) #dataframe to store 0 values
  
#caclulate the x and y value of the 'zero value curve'
  j=1
   while(j<10001){
    rj = (j)*1e-9  #step for nanoparticle size
    x_0=rj*10^(6)
    y_0=((rhoi-rhow)*(ri^2)/(rj^2))+rhow
   
    if(j<10){
      j=j+0.1
    } else if (j<100) {
      j=j+1
    } else if (j<1000) {
      j=j+10
    }else { #more than 1000
      j=j+100
    }
    k <- data.frame (x_0, y_0 )
    df2 <- rbind(df2, k) 
  }

x_curve = df2 %>% select(x_0)
y_curve = df2 %>% select(y_0)

#plot differntial sedimentation rate constants with curve of 0 values
p4 <- plot_ly(
  x = kcoll_data[,2], y = kcoll_data[,1], z = log10(kcoll_data[,5]),
  type = "contour", autocontour = F,   contours = list(
    start = -27,
    end = -14,
    size = 1
  )) %>% 
  colorbar(title = "log(rate \n constant \n in m<sup>3</sup>/s)", exponentformat = "E",   titlefont = f1, tickfont = f2
  ) %>% add_trace(x = x_curve[78:360,1], y= y_curve[78:360,1], type = "scatter", mode = "line", color = "orange") %>% 
  layout(
    title = "differential sedimentation",
    title = f1,
    xaxis = ax,
    yaxis = ay 
  )
p4

#plot overall rate constants
p5 <- plot_ly(
  x = kcoll_data[,2], y = kcoll_data[,1], z = log10(kcoll_data[,6]),
  type = "contour", autocontour = F   
    ) %>% 
  colorbar(title = "log(rate \n constant \n in m<sup>3</sup>/s)", exponentformat = "E",   titlefont = f1, tickfont = f2
  ) %>%
  layout(
    title = "overall collision rate",
    titlefont = f1,
    xaxis = ax,
    yaxis = ay 
  )
p5

#part 3: saving plots (optional)
#to save plots, install orca: https://github.com/plotly/orca
#and uncomment lines below

#save as png
#orca(p1, "peri.png")
#orca(p2, "ortho.png")
#orca(p3, "sedi_noZeros.png")
#orca(p4, "sedi_wZeros.png")
#orca(p5, "sum.png")


#save as svg
#orca(p1, "peri.svg")
#orca(p2, "ortho.svg")
#orca(p3, "sedi_noZeros.svg")
#orca(p4, "sedi_wZeros.svg")
#orca(p5, "sum.svg")
