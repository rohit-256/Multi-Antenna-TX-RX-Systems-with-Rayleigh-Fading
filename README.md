main1.m can simulate a 1 transmit L receive system using Maximal Ratio Combining (MRC) or a 2 transmit L receive system using Alamouti coding. 
A narrowband multipath fading channel is considered and the receiver is assumed to have perfect channel state information (CSI).

main1(scheme,L,modln,M)
scheme - 'MRC' or 'ALM'
L - number of receiver antennas
modln - 'MPSK' or 'MQAM'
M - number of points in the constellation

#####################################################3

main2.m considers the case where the receiver is designed based on the channel gain obtained at time t,
but the transmission happens at time t+kTs, where Ts is the symbol duration and k is some integer.
The receiver moves with a speed of v m/sec and the carrier frequency is fc Hz. The channels vary
with time as per Jakes’s correlation model

main2(scheme,L,modln,M,v,fc)
scheme - 'MRC' or 'ALM'
L - number of receiver antennas
modln - 'MPSK' or 'MQAM'
M - number of points in the constellation
v - velocity in m/s
fc - carrier frequency in Hz
