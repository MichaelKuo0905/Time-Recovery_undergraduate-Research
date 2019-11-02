%***********************************************************************
%             Simulate QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
clear all; close all;

%---------------- preparation part -------------------------
N=10000;               % Symbol length
trainLen=800;              % Training Length
OverSamp=4;                % Number of samples between each symbol
ml=4;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
l=50;                            % 1/2 length of pulse shape (in symbols)


%---------------- Channel and Clock Offset Effect ---------------
chan=[1 0.3 0.3 0];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%----------------- Start Calculation -----------------------
current_data=randn(1,N*ml)>0;
% save current_data
data_seq=Mapper(current_data,ml);  
m=data_seq;

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*OverSamp);                % upsample the data by placing...
sup(1:OverSamp:N*OverSamp)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,OverSamp,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter
train=sign(randn(1,1000));
x=[train x]*exp(j*pi/6);
a=conv(x,train(1000:-1:1));
[c,b]=sort(abs(a),2,'descend');
peak_position=b(1,1);
theata=angle(train(1,1000).*conj(x(1,peak_position)));
x=x*exp(j*theata);
x=x(1000+1:end);

%---------------- Resample to change clock period for clock linear drift -------------
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal
%---------------eye diagram----------------------------------------------
n=real(x);
neye=5;
c=floor(length(n)/(neye*OverSamp));
xp=n(4:end);  % dont plot transients at start
q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
figure(5)
 plot(q)
axis([0,10,-3,3])
ylabel('Amplitude'), xlabel('symbols')
title('Eye diagram for sinc pulse shape')

%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N*OverSamp);   % initialize variables
tausave=zeros(1,N*OverSamp); tausave(1)=tau; i=0;
mu=0.02;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp     % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau+mu*dx;              % alg update (energy)
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
plot(xs(2000:6000),'o')
axis([-2,2,-2,2])
ylabel('offset estimates'), xlabel('iterations')
%---------------eye diagram----------------------------------------------
b=real(xs(3000:8000));
neye=5;
c=floor(length(b)/neye);
xp=b(2:end);  % dont plot transients at start
q=reshape(xp,neye,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
figure(6)
 plot(q)
 ylabel('Amplitude'), xlabel('symbols')
title('Eye diagram for sinc pulse shape')
figure(7)
plotspec(xs,10000);


 





