 N=str2double(get(handles.NumberSymbol,'string'));
 ml=(get(handles.BitperSymbol,'value'));
 beta=str2double(get(handles.SRRCRolloffFactor,'string'));
 OverSamp=(get(handles.OverSamp,'value'));
 Q_Delay_Tb=(get(handles.Q_Delay_Tb,'value'));
 clock_drift=(get(handles.ClockDrift,'value'));
 l=str2double(get(handles.DelayPulseShape,'string'));
 adaptive_method=(get(handles.AdaptiveMethod,'value'));
 mu=str2double(get(handles.StepSize,'string'));
 average_power=str2double(get(handles.AveragePower,'string'));
%  fc=str2double(get(handles.Fc,'string'));
%  fs=str2double(get(handles.Fs,'string'));
%% ***************************************************************************************************
%%                                                           1. Overall Simulation Setup
Q_Delay_Tb=Q_Delay_Tb-1;
clock_drift=clock_drift-1;
load train
load packet_length
%---------------- Channel and Clock Offset Effect ---------------
% %*********************E4406A Instrument Setting****************
%    [sig_bb4_1,sig_IQ,fs]=E4406A_control;
%       save sig_bb4_1;                       
%   load sig_bb;                  %QPSK     
   load sig_bb_1;                %OQPSK ISI                                  
%   load sig_bb2;                 %8PSK  ISI         
%  load sig_bb3                   %16QAM ISI
%     load sig_bb4                %64QAM Ideal                
%      load sig_bb4_1              %64QAM  ISI                
RX_Input=resample(sig_bb_1,300e4,15e6);
%---------------- find  packet position----------------------
a=conv(RX_Input,train(100*OverSamp:-1:1));
axes(handles.axes1);
plot(abs(a));xlabel('Time'), ylabel('Amplitude')
[c,b]=sort(abs(a),2,'descend');
%-------  in case of find untire packet ----------------------------
for i=1
    if  b(1,i) > length(RX_Input)-N*OverSamp
        i=i+1;
    end
end
peak_position=b(1,i);
set(handles.PeakPosition,'String',num2str(peak_position))
cor=train.*conj(RX_Input(peak_position-length(train)+1:peak_position));
ph_fix=angle(sum(cor)/length(train));                         %find phase  
set(handles.PhaseOffset,'String',num2str(ph_fix))
RX_Input=RX_Input(b(1,i)+1:b(1,i)+packet_length);             % change to one packet
% RX_Input=RX_Input(21845:61840);
RX_Output=RX_Input/max(abs(RX_Input));
% figure(2)
% plot(RX_Output,'g');
RX_Output=RX_Output*exp(1i*ph_fix);
% figure(3),plot(RX_Output,'g')
x=RX_Output;
%---------------- Resample to change clock period for clock linear drift -------------
if clock_drift==1
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset

if  Q_Delay_Tb==1;              %OQPSK Signal
  x1=[zeros(1,OverSamp/2) real(z(1:end-OverSamp/2))]; % 實部虛部對齊
  x2=imag(z);
  x=x1+j*x2;                                        
else
   x=z;                          % relabel signal
end
else
    if  Q_Delay_Tb==1;              %OQPSK Signal
  x1=[zeros(1,OverSamp/2) real(x(1:end-OverSamp/2))]; % 實部虛部對齊
  x2=imag(x);
  x=x1+j*x2;                                        
else
   x=x;                          % relabel signal
    end
end
%---------------eye diagram----------------------------------------------
n=real(x);
neye=5;
c=floor(length(n)/(neye*OverSamp));
xp=n(4:end);  %6:OQPSK  4:otherwise            % dont plot transients at start
q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
axes(handles.axes2);
plot(q)
axis([1,5,-0.5,0.5])
ylabel('Amplitude'), xlabel('symbols')
title('Eye diagram for sinc pulse shape')

%------------------------adaptive method-----------------------------------
switch(adaptive_method)   
 case 1
 x1=real(x);
 x2=imag(x);
    
%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
%-------------------------Real PART--------------------------------------
tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);   % initialize variables
tausave1=zeros(1,N); tausave1(1)=tau1; i=0;
% mu=0.3;                            % algorithm stepsize
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
% mu=0.3;                            % algorithm stepsize
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
% mu=0.2;                            % algorithm stepsize
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
% mu=0.2;                            % algorithm stepsize  4QAM(ISI):0.02     
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
%---------------eye diagram----------------------------------------------
b=real(xs(3000:8000));
neye=5;
c=floor(length(b)/neye);
xp=b(2:end);  % dont plot transients at start
q=reshape(xp,neye,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
axes(handles.axes3)
plot(q);ylabel('Amplitude'),xlabel('symbols')



xs=xs(3000:8000);
axes(handles.axes4)
plot(real(xs),'b.')     % plot constellation diagram
 title('constellation diagram');
ylabel('Amplitude')
axes(handles.axes5)
 plot(imag(xs),'b.')      % plot constellation diagram
ylabel('Amplitude')
 axes(handles.axes6)
 plot(tausave(1:i-2))        % plot trajectory of tau
 ylabel('offset estimates'), xlabel('iterations')

axes(handles.axes7)
plot(x,'o')        
ylabel('Amplitude'), xlabel('iterations')
axes(handles.axes8)
plot(xs,'o')        
ylabel('Amplitude'), xlabel('iterations')
