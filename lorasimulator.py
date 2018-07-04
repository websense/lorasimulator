##   Lora Physical Layer Simulator Python code
##   Copyright (C) 2017  Benjamin Dix-Matthews

##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.


## To run: load this file into an interactive Python session using
## exec(open("./LoraSimulator.py").read())
## then call testParameters or parameterOptimisation 

import math

import matplotlib.pyplot as plt
import matplotlib  

# Plane Earth channel model
# NPE constant for terrain: NPE=4 standard, >4 for trees or buildings
# for best results fit a site-specific value to experimental data
# maxPathLoss transmit_dBm minus receiver sensitivity floor
# hTransmit, hReceive height in meters of transmit and receive antennas (respectively)
# return maxDist in kilometers
def planeEarth(maxPathLoss,hTransmit,hReceive,NPE):
    maxDist=pow(10,((maxPathLoss+20*math.log10(hTransmit*hReceive))/(10*NPE)))/1000
    return maxDist

# Spreading Factor (SF) can be 6-12 for 
# Coding Rate (CR) can be 1-4 
# Power can be -2 to 20 dBm 
# Bandwidth (BW) can be 125, 250 or 500 kHz (represented as 1, 2 or 3)
# parameters for Plane Earch channel model:
# NPE typically 4, higher for buildings or trees
# hTransmit, hReceive antenna heights in meters typically =1 meter
def testParameters(spreadingFactor,codingRate,programmedPower,BWidth,NPE,hTransmit,hReceive):
    PS = 32                    #Payload size in Bytes
    SF = spreadingFactor         #Spread factor
    H = 1                        #Explicit header (H=0 when the Explicit header is enabled)
    DE = 0                       #Low DR optimize (DE=0 when Low data rate optimisation is disabled)
    CR = codingRate              #Coding rate (from 1 to 4)
    PreSymb = 6                  #Preamble symbols
    BW = {1:125,2:250,3:500}[BWidth]        #Bandwidth in kHz
    transmit_dBm = programmedPower  #from -2 to 20 dB
    ACK = 3          #number of acknowledgment bytes
    #Determine Timing
    Tsym = float(pow(2, SF)) / (BW * 1000.0)*1000
    Tpreamble =(PreSymb+4.25)*Tsym 
    payloadSymbNb = 8 + max(math.ceil(float(8*PS-4*SF+28+16-20*H)/(4*(SF-2*DE)))*(CR+4),0)
    Tpayload = payloadSymbNb * Tsym
    Tpacket = Tpreamble + Tpayload
    #DReq = 1000*SF/Tsym*(4/(CR+4)) #not used
    #Determine Power Usage
    transmitPower={-2:22,-1:22,0:22,1:23,2:24,3:24,4:24,5:25,6:25,7:25,8:25,9:26,10:31,11:32,12:34,13:35,14:44,15:82,16:85,17:90,18:105,19:115,20:125}
    CAD_Rx={125:10.8,250:11.6,500:13}
    CAD_proc={125:5.6,250:6.5,500:8}
    sleepPower=0.1  
    #the current considering only the time on air of the original packet
    #disregards Ack and CADs
    energy_uC = float(Tpacket)*transmitPower[transmit_dBm]/(PS*1000.0) 
    #the current in micro joules considering only the time on air of the original packet
    #disregards Ack and CAD
    energy_uJ=Tpacket*transmitPower[transmit_dBm]*3.3/(1000)
    IDDR_L = (pow(2,SF)+32)/BW
    IDDC_L = pow(2,SF)*SF/1750
    totalCAD = IDDR_L+IDDC_L
    powerCAD = CAD_Rx[BW]*IDDR_L + CAD_proc[BW]*IDDC_L
    ackSymbNb =(PreSymb+4.25)*2+8+max(math.ceil((8*ACK-4*SF+28+16-20*H)/(4*(SF-2*DE)))*(CR+4),0)
    ackTime=ackSymbNb*Tsym
    # Receiver sensitivity floor
    if(BWidth==1):
        RFS={6:-121,7:-124,8:-127,9:-130,10:-133,11:-135,12:-137}
    elif(BWidth==2):
        RFS={6:-118,7:-122,8:-125,9:-128,10:-130,11:-132,12:-135}
    else:
        RFS={6:-111,7:-116,8:-119,9:-122,10:-125,11:-128,12:-129}
    RS=RFS[SF]   
    #To change the transmission channel alter the code in the following section
    #It is currently using the Plane Earth model
    maxPathLoss=transmit_dBm-RS
    MaxDist=planeEarth(maxPathLoss,hTransmit,hReceive,NPE)
    return [MaxDist,energy_uC,spreadingFactor,codingRate,BW,programmedPower]

## examples    
## print (testParameters(12,2,14,1,4,1,1) )   
## testParameters(6,2,14,1,4,1,1)    

## energy_uCMax, energy_uCMin, energy budget range as transmit charge per byte
## distMax,distMin, required transmission range in kilometers
## distWeighting,powerWeighting, weights for the utility function, default 1,1
## plotOn boolean whether to show graph of search space results
## returns currentBest vector: 
##         MaxDist,energy_uC,spreadingFactor,codingRate,BW,programmedPower
def parameterOptimisation(energy_uCMax,energy_uCMin,distMax,distMin,NPE,hTransmit,hReceive,distWeighting,powerWeighting,plotOn):
    transPower=0
    SF=6
    CR=1 #was 2
    BW=1
    testResults=[]
    power=[]
    dist=[]
    while transPower<21:
        while SF<13:
            #in this model coding rate only increases the length of packet
            #so the least FEC will always be chosen
            #
            #while CR<5:
            while BW<4:
                testResults.append(testParameters(SF,CR,transPower,BW,NPE,hTransmit,hReceive))
                BW+=1
            #CR+=1
            SF+=1
            BW=1
            #CR=1
        transPower+=1
        SF=6
    if(plotOn):
        plt.scatter([row[1] for row in testResults],[row[0] for row in testResults],alpha=0.4)
        plt.title('Plane Earth Parameter Optimisation')
        plt.grid(True)
        plt.xlabel('Transmit Charge per Byte (uC)')
        plt.ylabel('Transmit Distance (km)')
        plt.gca().set_xscale('log')
        plt.gca().set_autoscale_on(True)
        #plt.axis([0.005,10, 0,6]) #deprecated
    
    currentScore=None
    currentBest=[]
    for row in testResults:
        if((row[1]<energy_uCMax)&(row[1]>energy_uCMin)&(row[0]>distMin)&(row[0]<distMax)):
            score=row[0]*distWeighting-row[1]*powerWeighting
            if currentScore is None:
                currentScore=score
            if(score>currentScore):            
                currentBest=row
                currentScore=score
    #if error, just return empty best param list
    #if(len(currentBest)==0):
        #raise ValueError("There are no possibilities within this region")
    if(plotOn):
        #plt.hold(True) #deprecated
        #Determine the Optimum parameter in RED
        plt.scatter(currentBest[1],currentBest[0],c='red',marker='o')   
        #create a box for the feasible region
        #Rectangle((xl,yl),width=w,height=h)
        rect = matplotlib.patches.Rectangle((energy_uCMin,distMin),(energy_uCMax-energy_uCMin),(distMax-distMin),linewidth=1,edgecolor='r',facecolor='none')
        plt.gca().add_patch(rect)
        plt.show()
    return currentBest
    #return testResults #alternative - return the full results array for further processing
  
## examples
## parameterOptimisation(2,0.1,3,0,4,1,1,1,1,True)  
    
