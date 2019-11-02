%***********************************************************************
%             Simulate QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
% clear all; close all;

%---------------- preparation part -------------------------
N=75000;               % Message length
% trainLen=800;              % Training Length
OverSamp=4;                % Number of samples between each symbol
ml=2;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;               % SRRC filter rolloff factor
L=50;                            % 1/2 length of pulse shape (in symbols)
Delay=101;

%---------------- Channel and Clock Offset Effect ---------------
chan=[1 0 0 0];                        % T/m "channel"
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(L,beta,OverSamp,toffset);% srrc pulse shape with timing offset
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
matchfilt=srrc(L,beta,OverSamp,0);      % matched filter = srrc pulse shape
x=conv(r,matchfilt);             % convolve signal with matched filter

TX_Output=[zeros(1,Delay) x zeros(1,Delay)];


%*************** TX RF IQI Effect ***************************
%------- I/Q Imbalance Effect ---------------
% I=real(TX_Output);
% Q=imag(TX_Output);

% %------- Phase Imbalance Model --------------
% I_temp=I+Q*sin(PhaseImbalance);
% Q_temp=Q*cos(PhaseImbalance);
% 
% %------- Amplitude Imbalance Model ----------
% I_out=I_Gain*I_temp;
% Q_out=Q_Gain*Q_temp;
% RF_Output=I_out+j*Q_out;
RF_Output=x;
RF_Output=RF_Output/max(abs(RF_Output));
% save RF_Output;

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