%***********************************************************************
%             Simulate QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
% clear all; close all;

%---------------- preparation part -------------------------
N=2200;                   % Message length
% trainLen=800;            % Training Length
OverSamp=4;                % Number of samples between each symbol
ml=2;                      % Number of modulation levels (BPSK: ml=1, QPSK: ml=2, 16QAM: ml=4)
beta=0.25;                 % SRRC filter rolloff factor
L=50;                      % 1/2 length of pulse shape (in symbols)
Delay=101;

% load RF_Output;
% RF_Input=RF_Output;
% x=RF_Input;
% %*********************E4406A Instrument Setting****************
[sig_bb,sig_IQ,fs]=E4406A_control;
% save sig_bb
% load sig_bb
RX_Input=resample(sig_bb,15e4,15e6);
RX_Output=RX_Input/max(abs(RX_Input));
RX_Output=RX_Output/max(abs(RX_Output));
figure(6),plot(RX_Output,'g')
% x=RX_Output;
%*************Receiver Processing***************
I_out1=real(RX_Output(1:2:end)); 
Q_out1=imag(RX_Output(1:2:end)); 
I_out2=real(RX_Output(2:2:end)); 
Q_out2=imag(RX_Output(2:2:end)); 

val1=sum(I_out1(1:500).^2+Q_out1(1:500).^2);
val2=sum(I_out2(1:500).^2+Q_out2(1:500).^2);
if (val1>val2)
    I_out=I_out1;
    Q_out=Q_out1;
else
    I_out=I_out2;
    Q_out=Q_out2;
end

sample=fix(Delay/2)+2*L;

% %------- Least Square Method ----------------
% G=[I_Gain I_Gain*sin(PhaseImbalance);0 Q_Gain*cos(PhaseImbalance)];
% tmp=inv(G)*[I_out;Q_out];
% I_LS=tmp(1,:); Q_LS=tmp(2,:);
% 
% %------ EVM versus LS calculation ---------
% Ideal_Symbol=Input;
% LS_Symbol=I_LS+j*Q_LS;
% 
% ErrorVector=Ideal_Symbol-LS_Symbol(sample+1:sample+N);
% Error=mean(ErrorVector.*conj(ErrorVector));
% 
% EVM_LS=10*log10(Error);
% 
% %------ EVM versus IQI calculation ---------
% Ideal_Symbol=Input;
% IQI_Symbol=I_out+j*Q_out;
% ErrorVector=Ideal_Symbol-IQI_Symbol(sample+1:sample+N);
% Error=mean(ErrorVector.*conj(ErrorVector));
% 
% EVM_IQI=10*log10(Error);

%*******************************************
%------- Adaptive DM Method ----------------
f1=1; f2=0; f1_buf=[]; f2_buf=[];
f3=0; f4=1; f3_buf=[]; f4_buf=[];
muI=0.05; muQ=0.05;
% for i=1:562
for i=1:2250
    rI=I_out(i); 
    rQ=Q_out(i); 
    yI=f1*rI+f2*rQ;
    yQ=f3*rI+f4*rQ;
    eI=(1-yI^2)*yI;
    eQ=(1-yQ^2)*yQ;    
    f1=f1+muI*eI*rI;
    f2=f2+muI*eI*rQ;
    f3=f3+muQ*eQ*rI;
    f4=f4+muQ*eQ*rQ;
    f1_buf=[f1_buf f1];
    f2_buf=[f2_buf f2];
    f3_buf=[f3_buf f3];
    f4_buf=[f4_buf f4];
end
F_DM=[f1 f2;f3 f4];
tmp=F_DM*[I_out;Q_out];
I_DM=tmp(1,:); Q_DM=tmp(2,:);
x=I_DM+j*Q_DM;


%---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=L*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
mu=0.02;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-L*OverSamp          % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,L);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,L);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,L);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
 

end
%  for i = 1:1500
% %   i=i+1;
%   xs(i)=interpsinc(x,tnow+tau,L);   % interp at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,L);  % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,L);  % value to left
%   dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
%   tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
%   tnow=tnow+OverSamp; 
%   tausave(i)=tau;      % save for plotting
%  end
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
