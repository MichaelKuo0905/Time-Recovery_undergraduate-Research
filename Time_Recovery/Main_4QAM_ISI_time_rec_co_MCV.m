%***********************************************************************
%             Simulate PAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
clear all; close all;

%---------------- preparation part -------------------------
N=1000;               % Message length
trainLen=800;              % Training Length
OverSamp=2;                % Number of samples between each symbol
ml=2;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
l=50;                            % 1/2 length of pulse shape (in symbols)


%---------------- Channel and Clock Offset Effect ---------------
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset 
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%----------------- Start Calculation -----------------------
nloop=1;                 % Number of simulation loops
noe=0;                     % Number of error data
nod=0;                     % Number of TX data
%-------------       Data Generation -------------------------
current_data=randn(1,N*ml)>0;
data_seq=Mapper(current_data,ml);
m=data_seq;
    

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*ml);                % upsample the data by placing...
sup(1:OverSamp:N*2)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,OverSamp,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter
%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.002;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=x_deltap-x_deltam; % numerical derivative
  
  qx=quantalph(xs(i));  % quantize to alphabet
  tau=tau+mu*dx*(abs(qx)-abs(xs(i)));              % alg update (energy)
  tnow=tnow+OverSamp; tausave(i)=tau;      % save for plotting
        

end

figure(1) 
plot(x,'o')        
ylabel('offset estimates'), xlabel('iterations')
figure(2)
plot(xs,'o')        
ylabel('offset estimates'), xlabel('iterations')
figure(3)
subplot(2,1,1), plot(real(xs),'b.')        % plot constellation diagram
title('constellation diagram');
subplot(2,1,2), plot(imag(xs),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
figure(4)
 plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')