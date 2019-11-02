%***********************************************************************
%             Simulate QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
% clear all; close all;

%---------------- preparation part -------------------------
N=10000;               % Message length
% trainLen=800;              % Training Length
OverSamp=4;                % Number of samples between each symbol
ml=4;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
l=50;                            % 1/2 length of pulse shape (in symbols)
average_power=sqrt(10);

%---------------- Channel and Clock Offset Effect ---------------
chan=[1 0.7 0 0];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset 
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%-------------       Data Generation -------------------------
current_data=randn(1,N*ml)>0;
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
% x=x/max(abs(x));
% x=x*exp(j*pi/6);
 x=x/max(abs(x));
% save x x;
%---------------- Resample to change clock period for clock linear drift -------------
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal

x1=real(x);
x2=imag(x);
% %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
%-------------------------Real PART--------------------------------------
tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);    % initialize variables
tausave=zeros(1,N); tausave(1)=tau1; i=0;
 mu=0.9;                            % algorithm stepsize
delta=0.3;                          % time for derivative
while tnow1<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs1(i)=interpsinc(x1,tnow1+tau1,l,beta);   % interp at tnow+tau
  x_deltap1=interpsinc(x1,tnow1+tau1+delta,l,beta);  % value to right
  x_deltam1=interpsinc(x1,tnow1+tau1-delta,l,beta);  % value to left
  dx1=x_deltap1-x_deltam1; % numerical derivative
  
  qx1=quantalph1(xs1(i)*average_power,ml);           % quantize to alphabet
  tau1=tau1+mu*dx1*(qx1-xs1(i));              % alg update (energy)
  tnow1=tnow1+OverSamp; tausave(i)=tau1;      % save for plotting
end
%-----------------------------ImagPART--------------------------------------
tnow2=l*OverSamp+1; tau2=0; xs2=zeros(1,N);   % initialize variables
tausave2=zeros(1,N); tausave2(1)=tau2; i=0;
 mu=0.9;                            % algorithm stepsize
delta=0.3;                          % time for derivative
while tnow2<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs2(i)=interpsinc(x2,tnow2+tau2,l,beta);   % interp at tnow+tau
  x_deltap2=interpsinc(x2,tnow2+tau2+delta,l,beta);  % value to right
  x_deltam2=interpsinc(x2,tnow2+tau2-delta,l,beta);  % value to left
  dx2=x_deltap2-x_deltam2; % numerical derivative
  
  qx2=quantalph1(xs2(i)*average_power,ml);  % quantize to alphabet
  tau2=tau2+mu*dx2*(qx2-xs2(i));              % alg update (energy)
  tnow2=tnow2+OverSamp; tausave2(i)=tau2;      % save for plotting
end
 xs=xs1+j*xs2;

% %-------------------------Real PART--------------------------------------
% tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);   % initialize variables
% tausave1=zeros(1,N); tausave1(1)=tau1; i=0;
% mu=0.5;                            % algorithm stepsize
% delta=0.3;                          % time for derivative
% while tnow1<length(x)-2*l*OverSamp % run iteration
%     
%   i=i+1;
%   xs1(i)=interpsinc(x1,tnow1+tau1,l,beta);   % interp at tnow+tau
%   x_deltap1=interpsinc(x1,tnow1+tau1+delta,l,beta);  % value to right
%   x_deltam1=interpsinc(x1,tnow1+tau1-delta,l,beta);  % value to left
%   dx1=x_deltap1-x_deltam1; % numerical derivative
%   
%   qx1=quantalph1(xs1(i)*average_power,ml);           % quantize to alphabet
%   tau1=tau1+mu*dx1*(qx1-xs1(i));              % alg update (energy)
%   tnow1=tnow1+OverSamp; tausave1(i)=tau1;      % save for plotting
% end
% %-----------------------------ImagPART--------------------------------------
% tnow2=l*OverSamp+1; tau2=0; xs2=zeros(1,N);   % initialize variables
% tausave2=zeros(1,N); tausave2(1)=tau2; i=0;
% mu=0.5;                            % algorithm stepsize
% delta=0.3;                          % time for derivative
% while tnow2<length(x)-2*l*OverSamp % run iteration
%     
%   i=i+1;
%   xs2(i)=interpsinc(x2,tnow2+tau2,l,beta);   % interp at tnow+tau
%   x_deltap2=interpsinc(x2,tnow2+tau2+delta,l,beta);  % value to right
%   x_deltam2=interpsinc(x2,tnow2+tau2-delta,l,beta);  % value to left
%   dx2=x_deltap2-x_deltam2; % numerical derivative
%   
%   qx2=quantalph1(xs2(i)*average_power,ml);        % quantize to alphabet
%   tau2=tau2+mu*dx2*(qx2-xs2(i));               % alg update (energy)
%   tnow2=tnow2+OverSamp; tausave2(i)=tau2;       % save for plotting
% end
% xs=xs1*max(abs(xs1))+j*xs2*max(abs(xs2));
%-----------------plot----------------------------------------------------
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
plot(xs(5000:8000),'o')
% axis([-0.6,0.6,-0.6,0.6])
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
% figure(7)
% plotspec(xs,10000);

% figure(1) 
% plot(x,'o')        
% ylabel('offset estimates'), xlabel('iterations')
% figure(2)
% plot(xs(6000:10000),'o')        
% ylabel('offset estimates'), xlabel('iterations')
% figure(3)
% subplot(2,1,1), plot(xs1,'b.')        % plot constellation diagram
% title('constellation diagram');
% subplot(2,1,2), plot(xs2,'b.');       % plot constellation diagram
% title('constellation diagram');
% ylabel('estimated symbol values')
% figure(4)
%  plot(tausave1(1:i-2))        % plot trajectory of tau
% ylabel('offset estimates'), xlabel('iterations')