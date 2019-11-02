%***********************************************************************
%             Simulate PAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
clear all; close all;

%---------------- preparation part -------------------------
N=1000;               % Message length
trainLen=800;              % Training Length
OverSamp=4;                % Number of samples between each symbol
% ApproxPeriods=6;           % SRRC filter one-sided span in no. of symbols
% Tb=5*10^(-3);              % Baud interval in seconds
% Ts=Tb/OverSamp;            % Sampling period
% fc=200;                    % Carrier frequency (Hz)
ml=2;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
% N=10000;                         % number of data points
% M=2;                             % oversampling factor
% beta=0.3;                        % rolloff parameter for srrc
 l=50;                            % 1/2 length of pulse shape (in symbols)
% ApproxPeriods=6;                 % SRRC filter one-sided span in no. of symbols
% 
% d=randn(1,N*2)>0;                % random data sequence


%---------------- Channel and Clock Offset Effect ---------------
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset 
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%----------------- Start Calculation -----------------------
nloop=1;                 % Number of simulation loops
noe=0;                     % Number of error data
nod=0;                     % Number of TX data
    
for iii=1:nloop
    %------------- Data Generation -------------------------
    current_data=randn(1,N*ml)>0;
    data_seq=Mapper(current_data,ml);  
     m=data_seq;
%     train=randn(1,trainLen*ml)>0;
%     pre_data=randn(1,100*ml)>0;
%     next_data=randn(1,100*ml)>0;
%     
%     %------------- QPSK Modulation -------------------------
%     data_seq=1/sqrt(2)*((2*current_data(1:2:end)-1)+j*(2*current_data(2:2:end)-1));
%     train_seq=1/sqrt(2)*((2*train(1:2:end)-1)+j*(2*train(2:2:end)-1));
%     pre_data_seq=1/sqrt(2)*((2*pre_data(1:2:end)-1)+j*(2*pre_data(2:2:end)-1));
%     next_data_seq=1/sqrt(2)*((2*next_data(1:2:end)-1)+j*(2*next_data(2:2:end)-1));
%     m=[pre_data_seq train_seq data_seq next_data_seq]; 
% %     m=[ data_seq ];
 end

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*ml);                % upsample the data by placing...
sup(1:OverSamp:N*4)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,OverSamp,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1;  xs=zeros(1,N); xs2=zeros(1,N);% initialize variables
tau(1)=0;
tausave=zeros(1,N); tausave(1)=tau; i=0;
tausave2=zeros(1,N);
tau2=0;
mu1=0.01;                           % algorithm stepsize
mu2=0.002;
delta=0.1;                          % time for derivative

while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau+tau2,l);% interp at tnow+tau+tau2
  xs2(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  
  x_deltap2=interpsinc(x,tnow+tau+tau2+delta,l);  % value to right
  x_deltam2=interpsinc(x,tnow+tau+tau2-delta,l);  % value to left
  
  dx2=(abs(x_deltap2)-abs(x_deltam2))/delta;             % numerical derivative
  tau2=tau2-mu2*dx2*abs(xs(i))^5;            
  
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau-mu1*dx*abs(xs2(i))^5;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
  tausave2(i)=tau2;
 
end
figure(1)
subplot(2,1,1), plot(real(xs),'b.')        % plot constellation diagram
title('constellation diagram');
subplot(2,1,2), plot(imag(xs),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
figure(2)
 plot(tausave2(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

figure(3) 
plot(x,'o')        
ylabel('offset estimates'), xlabel('iterations')
figure(4)
plot(xs,'o')        
ylabel('offset estimates'), xlabel('iterations')
% % % %     %**************************** QPSK Demodulation *****************************
%     demodata=zeros(1,N*ml);
%     demodata(1:2:end)=real(xs)>0;
%     demodata(2:2:end)=imag(xs)>0;
%     
%     %************************** Bit Error Rate (BER) ****************************
%     
%     noe2=sum(abs(current_data-demodata));  % sum: built in function
%     nod2=length(N);  % length: built in function
%     noe=noe+noe2;
%     nod=nod+nod2;

 





% %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
% tnow=l*M+1; tau=0; xs=zeros(1,N);   % initialize variables
% tausave=zeros(1,N); tausave(1)=tau; i=0;
% mu=0.05;                            % algorithm stepsize
% delta=0.1;                          % time for derivative
% while tnow<length(x)-l*M            % run iteration
%   i=i+1;
%   xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
%   dx=x_deltap-x_deltam;             % numerical derivative
%   tau=tau+mu*dx*xs(i);              % alg update (energy)
%   tnow=tnow+M; tausave(i)=tau;      % save for plotting
% end
% 
% %---------------- Plot Results -------------------------
% subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
% title('constellation diagram');
% ylabel('estimated symbol values')
% subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
% ylabel('offset estimates'), xlabel('iterations')

