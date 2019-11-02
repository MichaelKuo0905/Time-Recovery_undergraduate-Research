%****************************************************************************************
%             Simulate PAM Signal with Clock Linear Drift
%             Dual Maximizing Output Power Method Used for Clock Recovery 
%****************************************************************************************
clear all; close all;

%---------------- PAM Signal Generator -------------------------
N=10000;                         % number of data points
M=2;                             % oversampling factor
beta=0.5;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
d=randn(1,N*2)>0;                % random data sequence
d_p=reshape(d,N,2);              % 2 bits for one symbol
table=[-3 3 -1 1];               % PAM candidate symbols (gray coding)
idx=d_p*[2;1]+1;                 % symbols index
m=table(idx);                    % 4-level PAM signal of length N

%---------------- Channel and Clock Offset Effect ---------------
chan=[1];                        % T/m "channel"
toffset=-0.2;                      % initial timing offset
pulshap=srrc(l,beta,M,toffset);  % srrc pulse shape with timing offset
hh=conv(pulshap,chan);           % ... and pulse shape

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*M);                % upsample the data by placing...
sup(1:M:N*M)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,M,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

%---------------- Resample to change clock period for clock linear drift -------------
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal

%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*M+1; 
tau1=0; xs1=zeros(1,N);   % initialize variables
tau1save=zeros(1,N); tau1save(1)=tau1; i=0;
tau2=0; xs2=zeros(1,N);   % initialize variables
tau2save=zeros(1,N); tau2save(1)=tau2;

mu1=0.01; mu2=0.001;                        % algorithm stepsizes
delta1=0.1;                          % time for derivative
delta2=0.1;                          % time for derivative
while tnow<length(x)-l*M            % run iteration
  i=i+1;
  xs1(i)=interpsinc(x,tnow+tau1,l);   % interp at tnow+tau
  x1_deltap=interpsinc(x,tnow+tau1+delta1,l);  % value to right
  x1_deltam=interpsinc(x,tnow+tau1-delta1,l);  % value to left
  dx1=x1_deltap-x1_deltam;             % numerical derivative
  tau1=tau1+mu1*dx1*xs1(i);              % alg update (energy)
  tau1save(i)=tau1;      % save for plotting
  
  xs2(i)=interpsinc(x,tnow+tau1+tau2,l);   % interp at tnow+tau
  x2_deltap=interpsinc(x,tnow+tau1+tau2+delta2,l);  % value to right
  x2_deltam=interpsinc(x,tnow+tau1+tau2-delta2,l);  % value to left
  dx2=x2_deltap-x2_deltam;             % numerical derivative
  tau2=tau2+mu2*dx2*xs2(i);              % alg update (energy)
  tnow=tnow+M; tau2save(i)=tau2;      % save for plotting
end

%---------------- Plot Results -------------------------
figure(1)
subplot(2,1,1), plot(xs1(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tau1save(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

figure(2)
subplot(2,1,1), plot(xs2(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tau2save(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

