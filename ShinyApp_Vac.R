## R shiny app to simulate the impact of vaccination:
## how population susceptibility & disease prevelance change with vaccination rate
## Model structure assumes a 2-age-group (children & adults) construct to test vaccination in children
## Disease dynamics follow an SIR form, i.e. assuming life-long immunity after infection/vaccination
## 2/24/18, by Wan Yang

library(shiny)
library(deSolve)
SIR2ageGrsVac <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
  
    dSC = nu*(1-p) - SC * (betaCC * IC + betaCA * IA) - muC * SC - lC * SC
    dIC = SC * (betaCC * IC + betaCA * IA) - gamma * IC - muC * IC - lC * IC
    dSA = lC * SC - SA * (betaAC * IC + betaAA * IA) - muA * SA
    dIA = lC * IC + SA * (betaAC * IC + betaAA * IA) - gamma * IA - muA * IA
    
    list(c(dSC, dIC, dSA, dIA))
  })
}

# parameters taken from Keeling & Rohani 2008
BETA=matrix(c(1100,110,110,220),2,2)
lC=0.066667; # 1/15; 
muC=0; muA=0.0166667; #1/(75-15); # mean life expectancy=75yr; subtract the 15 yrs in youth with 0 death rate
nC=muA/(lC+muA); nA=1-nC;
nu=(lC+muA)*nC;

IC0=.00005; SC0=0.05;
IA0=.00001; SA0=0.03;
state=c(SC=SC0,IC=IC0,SA=SA0,IA=IA0)
## run the model for 100 years to test the model
times=seq(1,200,by=7/365)


ui <- fluidPage(h2('Impact of Vaccination'),
                
                # outputs
                sidebarLayout(
                  sidebarPanel(width=4,
                               sliderInput(inputId = 'vacc',h6('Vaccination rate (%)'),value=50,min=0,max=100,step=1),
                               h6('Assume 100% Vaccine effectiveness'),
                               h6('Simulation based on a 2 age-group model with the following parameters:'),
                               h6('R0 = 15'),
                               h6('Infectious period: 14 days'),
                               h6('Children: 20% of population'),
                               h6('Adults: 80% of population')),
                              
                              mainPanel(
                                plotOutput(outputId = 'plots',width = '100%', height = "550px")
                                )
                              )
                
                )

server <- function(input, output){
  
  output$plots=renderPlot({
    R0=15; D=14/365;
    p=input$vacc/100
    # scale BETA to get the right R0
    r0=eigen(diag(c(.2,.8),2,2)%*%BETA)$values[1]*D  #  eigen(diag(c(.2,.8),2,2)%*%BETA0)$values[1]*D 
    betaCC.rt=BETA[1,1]*R0/r0;
    betaCA.rt=BETA[1,2]*R0/r0;
    betaAC.rt=BETA[2,1]*R0/r0;
    betaAA.rt=BETA[2,2]*R0/r0;
    
    parametersVac=c(betaCC=betaCC.rt,betaCA=betaCA.rt, betaAC=betaAC.rt, betaAA=betaAA.rt,
                    gamma=1/D,lC=lC,muC=muC,muA=muA,p=p)
    ## run the model for 100 years to test the model
    simVac=ode(times = times, func = SIR2ageGrsVac,y = state, parms = parametersVac)
    fSVac=simVac[,c('SC','SA')]/matrix(c(nC,nA),nrow(simVac),2,byrow=T)
    fIVac=simVac[,c('IC','IA')]/matrix(c(nC,nA),nrow(simVac),2,byrow=T)
    Stot=rowSums(simVac[,c('SC','SA')]);
    Itot=rowSums(simVac[,c('IC','IA')])
    
    # plot results
    par(mfrow=c(2,1),mar=c(3,3,1,1),cex=1.2,mgp=c(1.6,.5,0),cex.axis=.9)
    matplot(fSVac,xaxt='n',ylab='Proportion of class susceptible, Si/ni',xlab='Time (years)',
            main=paste0('Vaccination rate = ',p*100,'%'),type='l',col=c('blue','red'),lty=1,lwd=1.5)
    axis(1,at=seq(1,365/7*200,by=365/7),lab=1:200)
    legend('right',c('Children','Adults'),col=c('blue','red'),lty=1,lwd=1.5,cex=.9, bty='n')
    mtext(paste0('Population immunity at end of simulation: ',100-round(tail(Stot,1)*100,0),'%'),side=3,line=-1.1,outer=F)
    
    matplot(fIVac,xaxt='n',ylab='Proportion of class prevalence, Ii/ni',xlab='Time (years)',
            main=paste0('Vaccination rate = ',p*100,'%'),type='l',col=c('blue','red'),lty=1,lwd=1.5)
    axis(1,at=seq(1,365/7*200,by=365/7),lab=1:200)
    legend('right',c('Children','Adults'),col=c('blue','red'),lty=1,lwd=1.5,cex=1, bty='n')
    mtext(paste0('Prevalence at end of simulation: ',round(tail(Itot,1)*100,4),'%'),side=3,line=-1.1,outer=F)
    
  })
  
}

shinyApp(ui=ui, server = server)