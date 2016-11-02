library(ggplot2)
library(plyr)
X <- read.csv("/dos/Google_Drive/UNI_magistrale/2_anno/2_semestre_tesi/postprocessing/results/cdc1200/nacl/8800_160/1.0V/5_charge_results/tot_charge_run000_run056.dat",sep='\t',comment.char='#', header = FALSE)
Q_pos_3<-X[,3]
Q_pos_4<-X[,4]
num_data<-length(Q_pos_3)
time_data<- seq(0,num_data-1)*0.001
#plot(time_data,Q_pos_3,type='l')
cost1=min(Q_pos_3)
Q_1max<-max(Q_pos_3)
Q_2max<-max(Q_pos_4)
V_0<-1.0
C_1<-2.0*Q_1max/V_0
C_2<-2.0*Q_2max/V_0
#reduced units
t_0<-10^(-9)
R_0<-10^9
C_0<-1.60217662*10^(-19)
Q_0<-1.60217662*10^(-19)
a_red_conv<-(t_0/(R_0*C_0))
b_red_conv<-(t_0^2/(R_0^2*C_0^2))
c_red_conv<-(t_0^2*V_0/((R_0^2)*Q_0*C_0))
d_red_conv<-(t_0*V_0/(R_0*Q_0))

R_bulk<-4.060890

# Q1 <- function(t, R_l) {
#   a<-a_red_conv*((R_bulk+2*R_l)*C_1+(R_bulk+4*R_l)*C_2)/(R_l*(R_bulk+2*R_l)*C_1*C_2)*(10/1.602)
#   b<-b_red_conv*2.0/(R_l*(R_bulk+2*R_l)*C_1*C_2)
#   c_1<-c_red_conv*1.0/(R_l*(R_bulk+2*R_l)*C_2)
#   d<-d_red_conv*1.0/(R_bulk+2*R_l)
#   tau_1<-(2.0/(a+sqrt(a^2-4.0*b)))
#   tau_2<-(2.0/(a-sqrt(a^2-4.0*b)))
#   tildeA_1<-0.5*(1.0+((2.0*b*d-a*c_1)/(2*c_1*sqrt(a^2-4*b))))
#   tildeA_2<-0.5*(1.0-((2.0*b*d-a*c_1)/(2*c_1*sqrt(a^2-4*b))))
#   cost1 +Q_1max*(1-tildeA_1*exp(-t/tau_1) -tildeA_2*exp(-t/tau_2))
# }

Q1 <- function(t, tau_1,tau_2, tildeA) {
  cost1 +Q_1max*(1-tildeA*exp(-t/tau_1) -(1-tildeA)*exp(-t/tau_2))
}

df <- data.frame(time =time_data, charge = Q_pos_3, fit_charge=Q1(time_data,2,2,0.5))
ggplot(df, aes(x=time))+
  geom_line(aes(y=charge),color='red')+
  geom_line(aes(y=fit_charge),color='blue')+
  guides(color=FALSE)



nls(time_data ~ Q1(time_data))


