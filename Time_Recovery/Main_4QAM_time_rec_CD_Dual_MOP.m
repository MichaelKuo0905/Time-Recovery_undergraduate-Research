%***********************************************************************
%             Simulate QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
clear all; close all;

%---------------- preparation part -------------------------
N=1000;               % Message length
trainLen=800;              % Training Length
OverSamp=4;                % Number of samples between each symbol
ml=2;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
l=60;                            % 1/2 length of pulse shape (in symbols)


%---------------- Channel and Clock Offset Effect ---------------
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%----------------- Start Calculation -----------------------
% nloop=1;                 % Number of simulation loops
% noe=0;                     % Number of error data
% nod=0;                     % Number of TX data
current_data=randn(1,N*ml)>0;
data_seq=Mapper(current_data,ml);  
m=data_seq;
% for iii=1:nloop
%     %------------- Data Generation -------------------------
%     current_data=randn(1,N*ml)>0;
% %     train=randn(1,trainLen*ml)>0;
% %     pre_data=randn(1,100*ml)>0;
% %     next_data=randn(1,100*ml)>0;
%     
%     %------------- QPSK Modulation -------------------------
%     data_seq=1/sqrt(2)*((2*current_data(1:2:end)-1)+j*(2*current_data(2:2:end)-1));
% %     train_seq=1/sqrt(2)*((2*train(1:2:end)-1)+j*(2*train(2:2:end)-1));
% %     pre_data_seq=1/sqrt(2)*((2*pre_data(1:2:end)-1)+j*(2*pre_data(2:2:end)-1));
% %     next_data_seq=1/sqrt(2)*((2*next_data(1:2:end)-1)+j*(2*next_data(2:2:end)-1));
% %     m=[pre_data_seq train_seq data_seq next_data_seq]; 
%     m=[ data_seq ];
% end

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*ml);                % upsample the data by placing...
sup(1:OverSamp:N*4)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,OverSamp,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

%---------------- Resample to change clock period for clock linear drift -------------
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal
%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1;  xs=zeros(1,N); xs2=zeros(1,N);% initialize variables
tau(1)=0;
tausave1=zeros(1,N); tausave1(1)=tau; i=0;
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
  tausave1(i)=tau;      % save for plotting
  tausave2(i)=tau2;
 
end

% %---------------eye diagram----------------------------------------------
% n=real(x);
% neye=5;
% c=floor(length(n)/(neye*OverSamp));
% xp=n(1:end);  % dont plot transients at start
% q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% figure(5)
%  plot(q)
% axis([0,10,-3,3])
% title('Eye diagram for sinc pulse shape')
% %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
% tnow=l*OverSamp+1; 
% tau1=0; xs1=zeros(1,N);   % initialize variables
% tau1save=zeros(1,N); tau1save(1)=tau1; i=0;
% tau2=0; xs2=zeros(1,N);   % initialize variables
% tau2save=zeros(1,N); tau2save(1)=tau2;
% 
% mu1=0.01; mu2=0.001;                        % algorithm stepsizes
% delta1=0.1;                          % time for derivative
% delta2=0.1;                          % time for derivative
% while tnow<length(x)-l*OverSamp            % run iteration
%   i=i+1;
%   xs1(i)=interpsinc(x,tnow+tau1,l);   % interp at tnow+tau
%   x1_deltap=interpsinc(x,tnow+tau1+delta1,l);  % value to right
%   x1_deltam=interpsinc(x,tnow+tau1-delta1,l);  % value to left
%   dx1=(abs(x1_deltap)-abs(x1_deltam))/delta1;             % numerical derivative
%   tau1=tau1+mu1*dx1*abs(xs1(i))^5;              % alg update (energy)
%   tau1save(i)=tau1;      % save for plotting
%   
%   xs2(i)=interpsinc(x,tnow+tau1+tau2,l);   % interp at tnow+tau
%   x2_deltap=interpsinc(x,tnow+tau1+tau2+delta2,l);  % value to right
%   x2_deltam=interpsinc(x,tnow+tau1+tau2-delta2,l);  % value to left
%   dx2=(abs(x2_deltap)-abs(x2_deltam))/delta2;             % numerical derivative
%   tau2=tau2+mu2*dx2*abs(xs2(i))^5;              % alg update (energy)
%   tnow=tnow+OverSamp; tau2save(i)=tau2;      % save for plotting
% end

% %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
% tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
% tausave=zeros(1,N); tausave(1)=tau; i=0;
% mu=0.002;                            % algorithm stepsize
% delta=0.1;                          % time for derivative
% while tnow<length(x)-l*OverSamp            % run iteration
%   i=i+1;
%   xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
%   dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
%   tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
%   tnow=tnow+OverSamp; 
%   tausave(i)=tau;      % save for plotting
%  
% 
% end
figure(1)
subplot(2,1,1), plot(real(xs2),'b.')        % plot constellation diagram
title('constellation diagram');
subplot(2,1,2), plot(imag(xs2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
figure(2)
subplot(2,1,1), plot(tausave1(1:i-2))
subplot(2,1,2), plot(tausave2(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

% figure(2)
%  plot(tausave(1:i-2))        % plot trajectory of tau
% ylabel('offset estimates'), xlabel('iterations')

figure(3) 
plot(x,'o')        
ylabel('offset estimates'), xlabel('iterations')
figure(4)
plot(xs2,'o')        
ylabel('offset estimates'), xlabel('iterations')
% %---------------eye diagram----------------------------------------------
% b=real(xs2);
% neye=5;
% c=floor(length(b)/(neye*OverSamp));
% xp=b(1:end);  % dont plot transients at start
% q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% figure(6)
%  plot(q)
% axis([0,10,-2,2])
% title('Eye diagram for sinc pulse shape')
% 

% % % %     %**************************** QPSK Demodulation *****************************
%     demodata=length(xs);
%     demodata(1:2:end)=real(xs)>0;
%     demodata(2:2:end)=imag(xs)>0;
%     
%     %************************** Bit Error Rate (BER) ****************************
%     
%     noe2=sum(abs(current_data-demodata));  % sum: built in function
%     nod2=length(N);  % length: built in function
%     noe=noe+noe2;
%     nod=nod+nod2;


