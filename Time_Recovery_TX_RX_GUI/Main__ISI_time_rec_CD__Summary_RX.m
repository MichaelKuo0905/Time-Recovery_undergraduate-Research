%***********************************************************************
%             Simulate all kinds of  Signal with Clock Drift
%             different Methods Used for Clock Recovery 
%***********************************************************************
% clear all; close all;

%---------------- preparation part -------------------------
N=10000;                   % Message length
% trainLen=800;            % Training Length
OverSamp=4;                % Number of samples between each symbol
ml=4;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;                 % SRRC filter rolloff factor
l=50;                      % 1/2 length of pulse shape (in symbols)
Delay=101;
Q_Delay_Tb=1;              % OQPSK
load train;
load  packet_length;
adaptive_method=3;        % 1:Minimizing cluster variance; 2:Maximizing Output Power; 3: Absolute  Sampler(one degree Output Power)    
% %*********************E4406A Instrument Setting****************
%    [sig_bb,sig_IQ,fs]=E4406A_control;
%  save sig_bb;            
  load sig_bb;                  % 8PSK ISI
%  load sig_bb2;                               
%  load sig_bb1;                  
%  load sig_bb3                 
%     load sig_bb4                  
%     load sig_bb4_1                 
RX_Input=resample(sig_bb,300e4,15e6);
%---------------- find  packet position----------------------
a=conv(RX_Input,train(100*OverSamp:-1:1));
figure(1)
plot(abs(a));
[c,b]=sort(abs(a),2,'descend');
%-------  in case of find untire packet ----------------------------
for i=1
    if  b(1,i) > length(RX_Input)-N*OverSamp
        i=i+1;
    end
end
peak_position=b(1,i);
cor=train.*conj(RX_Input(peak_position-length(train)+1:peak_position));
ph_fix=angle(sum(cor)/length(train));                         %find phase               
RX_Input=RX_Input(b(1,i)+1:b(1,i)+packet_length);             % change to one packet
% RX_Input=RX_Input(21845:61840);
RX_Output=RX_Input/max(abs(RX_Input));
figure(2)
plot(RX_Output,'g');
RX_Output=RX_Output*exp(j*ph_fix);
figure(3),plot(RX_Output,'g')
x=RX_Output;
%---------------- Resample to change clock period for clock linear drift -------------
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
if  Q_Delay_Tb=1;              %OQPSK Signal
  x1=[zeros(1,OverSamp/2) real(z(1:end-OverSamp/2))]; % 實部虛部對齊
  x2=imag(z);
  x=x1+j*x2;                                        
else
   x=z;                          % relabel signal
end
%------------------------adaptive method-----------------------------------
switch(adaptive_method)   
 case 1
%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
%-------------------------Real PART--------------------------------------
tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);   % initialize variables
tausave1=zeros(1,N); tausave1(1)=tau1; i=0;
mu=0.3;                            % algorithm stepsize
delta=0.3;                          % time for derivative
while tnow1<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs1(i)=interpsinc(x1,tnow1+tau1,l,beta);   % interp at tnow+tau
  x_deltap1=interpsinc(x1,tnow1+tau1+delta,l,beta);  % value to right
  x_deltam1=interpsinc(x1,tnow1+tau1-delta,l,beta);  % value to left
  dx1=x_deltap1-x_deltam1; % numerical derivative
  
  qx1=quantalph1(xs1(i)*sqrt(10),ml);           % quantize to alphabet
  tau1=tau1+mu*dx1*(qx1-xs1(i));              % alg update (energy)
  tnow1=tnow1+OverSamp; tausave1(i)=tau1;      % save for plotting
end
%-----------------------------ImagPART--------------------------------------
tnow2=l*OverSamp+1; tau2=0; xs2=zeros(1,N);   % initialize variables
tausave2=zeros(1,N); tausave2(1)=tau2; i=0;
mu=0.3;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow2<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs2(i)=interpsinc(x2,tnow2+tau2,l,beta);   % interp at tnow+tau
  x_deltap2=interpsinc(x2,tnow2+tau2+delta,l,beta);  % value to right
  x_deltam2=interpsinc(x2,tnow2+tau2-delta,l,beta);  % value to left
  dx2=x_deltap2-x_deltam2; % numerical derivative
  
  qx2=quantalph1(xs2(i)*sqrt(10),ml);  % quantize to alphabet
  tau2=tau2+mu*dx2*(qx2-xs2(i));              % alg update (energy)
  tnow2=tnow2+OverSamp; tausave2(i)=tau2;      % save for plotting
end
xs=xs1+j*xs2;
 case 2
 %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.2;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
 

end
 %---------------- Clock Recovery Using the Absolute Sampler Method ---------------   
 case 3
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.2;                            % algorithm stepsize  4QAM(ISI):0.02     
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau+mu*dx;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
end
end 
xs=xs(2000:8000);
figure(4)
subplot(2,1,1), plot(real(xs),'b.')     % plot constellation diagram
title('constellation diagram');
subplot(2,1,2), plot(imag(xs),'b.')      % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
figure(5)
 plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

figure(6) 
plot(x,'o')        
ylabel('offset estimates'), xlabel('iterations')
figure(7)
plot(xs,'o')        
ylabel('offset estimates'), xlabel('iterations')
