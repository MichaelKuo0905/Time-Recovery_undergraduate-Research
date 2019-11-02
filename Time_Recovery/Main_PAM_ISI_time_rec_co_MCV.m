%************************************************************************************
%             Simulate PAM Signal with ISI and Clock Offset
%             Minimizing Cluster Variance Method Used for Clock Recovery 
%************************************************************************************
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
chan=[1 0.3];            % T/m "channel"  ISI Channel
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,M,toffset);  % srrc pulse shape with timing offset
hh=conv(pulshap,chan);           % ... and pulse shape

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*M);                % upsample the data by placing...
sup(1:M:N*M)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,M,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
tnow=l*M+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.002;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*M          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,3,-1,1]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+M; tausave(i)=tau;      % save for plotting
end

%---------------- Plot Results -------------------------
figure(1)
plot(x(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
figure(2)
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')


