%***********************************************************************
%             Simulate PAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
clear all; close all;

%---------------- PAM Signal Generator -------------------------
N=10000;                         % number of data points
M=2;                             % oversampling factor
beta=0.3;                        % rolloff parameter for srrc
l=50;                            % 1/2 length of pulse shape (in symbols)
d=randn(1,N*2)>0;                % random data sequence
d_p=reshape(d,N,2);              % 2 bits for one symbol
table=[-3 3 -1 1];               % PAM candidate symbols (gray coding)
idx=d_p*[2;1]+1;                 % symbols index
m=table(idx);                    % 4-level PAM signal of length N
Delay=101;
%---------------- Channel and Clock Offset Effect ---------------
chan=[1];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap=srrc(l,beta,M,toffset);  % srrc pulse shape with timing offset
hh=conv(pulshap,chan);           % ... and pulse shape

%---------------- RX Signal -------------------------------------
sup=zeros(1,N*M);                % upsample the data by placing...
sup(1:M:N*M)=m;                  % ... M-1 zeros between each data point
r=conv(hh,sup);                  % ... to get received signal
matchfilt=srrc(l,beta,M,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter
TX_Output=[zeros(1,Delay) x zeros(1,Delay)];

%*************** TX RF IQI Effect ***************************
%------- I/Q Imbalance Effect ---------------
I=real(TX_Output);
Q=imag(TX_Output);

% %------- Phase Imbalance Model --------------
% I_temp=I+Q*sin(PhaseImbalance);
% Q_temp=Q*cos(PhaseImbalance);
% 
% %------- Amplitude Imbalance Model ----------
% I_out=I_Gain*I_temp;
% Q_out=Q_Gain*Q_temp;
% RF_Output=I_out+j*Q_out;
RF_Output=I+j*Q;
RF_Output=RF_Output/max(abs(RF_Output));

% %***************Channel Effect and AWGN Effect***************************
% %--------------------Attennation Calculation ----------------------------
% spow=real(Sample_Output*Sample_Output')/N;
% attn=sqrt(0.5*spow/ml*10^(-ebn0/10));
% 
% %-------------------Rayleigh fading channel---------------------------
% fading=1;
% RX_fad=fading(1)*RF_Output;
% 
% %-------------------- AWGN -------------------
% [ich5,qch5]=comb(real(RX_fad),imag(RX_fad),attn);
% RX_Input=ich5+j*qch5;

% ---------------------- For E4438C RF TX Setting ----------------------
fc = 0.7e9;                  %   Unit : MHz
fs=15e4;                       
[status] = E4438C_control(RF_Output, fs, fc,'192.168.0.8')
