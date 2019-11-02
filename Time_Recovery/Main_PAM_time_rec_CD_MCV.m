%***********************************************************************************
%             Simulate PAM Signal with Clock Linear Drift
%             Minimizing Cluster Variance Method Used for Clock Recovery 
%***********************************************************************************
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
chan=[1 0.7];                        % T/m "channel"
toffset=-1;                    % initial timing offset
pulshap=srrc(l,beta,M,toffset);  % srrc pulse shape with timing offset
hh=conv(pulshap,chan);           % ... and pulse shape

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*M);                % upsample the data by placing...
sup(1:M:N*M)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,M,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter
% %---------------eye diagram----------------------------------------------
% n=real(x);
% neye=5;
% c=floor(length(n)/(neye*M));
% xp=n(1:end);  % dont plot transients at start
% q=reshape(xp,neye*M,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% figure(1)
%  plot(q)
% axis([0,10,-4,4])
% title('Eye diagram for sinc pulse shape')

%---------------- Resample to change clock period for clock linear drift -------------
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset
x=z;                          % relabel signal
% %---------------eye diagram----------------------------------------------
% n=real(x);
% neye=5;
% c=floor(length(n)/(neye*M));
% xp=n(1:end);  % dont plot transients at start
% q=reshape(xp,neye*M,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% figure(2)
%  plot(q)
% axis([0,10,-4,4])
% title('Eye diagram for sinc pulse shape')

%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
tnow=l*M+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.02;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-2*l*M          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l);   % interp value at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l); % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l); % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  qx=quantalph(xs(i),[-3,-1,1,3]);  % quantize to alphabet
  tau=tau+mu*dx*(qx-xs(i));         % alg update: DD
  tnow=tnow+M; tausave(i)=tau;      % save for plotting
end

%---------------- Plot Results -------------------------
figure(3)
subplot(2,1,1), plot(xs(1:i-2),'b.')        % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')
% %---------------eye diagram----------------------------------------------
% n=real(xs);
% neye=15;
% c=floor(length(n)/(neye*M));
% xp=n(1:end);  % dont plot transients at start
% q=reshape(xp,neye*M,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% figure(4)
%  plot(q)
% axis([0,10,-4,4])
% title('Eye diagram for sinc pulse shape')


