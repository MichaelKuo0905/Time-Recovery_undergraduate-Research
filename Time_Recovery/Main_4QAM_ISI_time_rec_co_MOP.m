%***********************************************************************
%             Simulate QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
% clear all; close all;

%---------------- preparation part -------------------------
N=10000;               % Message length
% trainLen=800;              % Training Length
OverSamp=2;               % Number of samples between each symbol
ml=2;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
l=50;                            % 1/2 length of pulse shape (in symbols)

%---------------- Channel and Clock Offset Effect ---------------
chan=[1 0.7 0 0];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%-----------------    data generation -----------------------
% nloop=1;                 % Number of simulation loops
% noe=0;                     % Number of error data
% nod=0;                     % Number of TX data
current_data=randn(1,N*ml)>0;
data_seq=Mapper(current_data,ml);  
m=data_seq;

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*ml);                % upsample the data by placing...
sup(1:OverSamp:N*OverSamp)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,OverSamp,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter
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


%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.0002;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
 

end
figure(1)
subplot(2,1,1), plot(real(xs),'b.')        % plot constellation diagram
title('constellation diagram');
subplot(2,1,2), plot(imag(xs),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
figure(2)
 plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

figure(3) 
plot(x,'o')        
ylabel('offset estimates'), xlabel('iterations')
figure(4)
plot(xs,'o')        
ylabel('offset estimates'), xlabel('iterations')
% %---------------eye diagram----------------------------------------------
% b=real(xs);
% neye=5;
% c=floor(length(b)/(neye*OverSamp));
% xp=b(1:end);  % dont plot transients at start
% q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% figure(6)
%  plot(q)
% axis([0,10,-2,2])
% title('Eye diagram for sinc pulse shape')

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

 





