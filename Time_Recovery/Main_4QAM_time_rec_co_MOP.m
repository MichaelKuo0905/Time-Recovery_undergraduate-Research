%***********************************************************************
%             Simulate 16QAM Signal with Clock Offset
%             Maximizing Output Power Method Used for Clock Recovery 
%***********************************************************************
clear all; close all;
   % Program 3-21
% qam16
%
% Simulation program to realize 16QAM transmission system
%
% Programmed by J. H. Deng
%

%******************** preparation part *************************************

sr=256000.0; % Symbol rate
ml=4;        % ml:Number of modulation levels (BPSK:ml=1, QPSK:ml=2, 16QAM:ml=4)
br=sr .* ml; % Bit rate
nd = 500;    % Number of symbols that simulates in each loop
ebn0=100;      % Eb/N0
IPOINT=8;    % Number of oversamples
l=50;        %1/2 length of pulse shape
%********************** Filter initialization   **************************

OVR=8;			
L=3;          % SRRC filter one-sided span in no. of symbols
rolloff=0.5;  % SRRC filter rolloff factor 
toffset=0.3;
rcpul=rcosine(l,OVR,'fir',rolloff,toffset);  % SRRC pulse shaping

%************************** START CALCULATION *******************************

nloop=1000;  % Number of simulation loops

noe = 0;    % Number of error data
nod = 0;    % Number of transmitted data

for iii=1:nloop
    
%*************************** Data generation ********************************

    
	data=rand(1,nd*ml)>0.5;

%*************************** 16QAM Modulation ********************************

	[ich,qch]=qammod(data,1,nd,ml);
     ich=[zeros(1,L) ich];    % Phesudo the initial values being zeros
     qch=[zeros(1,L) qch];    % Phesudo the initial values being zeros

	[ich1,qch1]= compoversamp(ich,qch,length(ich),IPOINT); 
	[ich2,qch2]= compconv(ich1,qch1,rcpul); 

%**************************** Attenuation Calculation ***********************
	
%     spow=sum(ich2.*ich2+qch2.*qch2)/nd;
% 	attn=0.5*spow*sr/br*10.^(-ebn0/10);
% 	attn=sqrt(attn);

%********************** Fading channel **********************

  % Generated data are fed into a fading simulator
  %  [ifade,qfade]=sefade(ich2,qch2,itau,dlvl,th1,n0,itnd1,now1,length(ich2),tstp,fd,flat);
  
  % Updata fading counter
  %  itnd1 = itnd1+ itnd0;

%********************* Add White Gaussian Noise (AWGN) **********************
	matchfilt=rcosine(l,OVR,'fir',rolloff,0);
%     [ich3,qch3]= comb(ich2,qch2,attn);% add white gaussian noise
	[ich4,qch4]= compconv(ich2,qch2,matchfilt);

	sampl=2*L*IPOINT+1; % note
	ich5 = ich4(sampl:IPOINT:length(ich4));
	qch5 = qch4(sampl:IPOINT:length(ich4));

    ich5 = ich5(L+1:nd+L);
    qch5 = qch5(L+1:nd+L);
    
%**************************** 16QAM Demodulation *****************************
	
    [demodata]=qamdemod(ich5,qch5,1,nd,ml);

%******************** Bit Error Rate (BER) ****************************
	
    noe2=sum(abs(data-demodata));
	nod2=length(data);
	noe=noe+noe2;
	nod=nod+nod2;

	fprintf('%d\t%e\n',iii,noe2/nod2);
end % for iii=1:nloop    

%********************** Output result ***************************

ber = noe/nod;
fprintf('%d\t%d\t%d\t%e\n',ebn0,noe,nod,noe/nod);
fid = fopen('BERqam.dat','a');
fprintf(fid,'%d\t%e\t%f\t%f\t\n',ebn0,noe/nod,noe,nod);
fclose(fid);

%******************** end of file ***************************

% % Hardware 
% data_rom=[data zeros(1,2048-length(data))];
% aa=reshape(data,2,length(data)/2);
% bb=[2 1]*aa+1;
% I_ch_idx=bb(1:2:end);
% Q_ch_idx=bb(2:2:end);
% 
% figure(1)
% subplot(211)
% plot(ich(1:100))
% title('ich output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(ich_HW(1:100))
% title('ich output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(2)
% subplot(211)
% plot(ich1(1:400))
% title('ich1 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(ich1_HW(1:400))
% title('ich1 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(3)
% subplot(211)
% plot(ich2(1:500))
% title('ich2 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(ich2_HW(1:500))
% title('ich2 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(4)
% subplot(211)
% plot(qch(1:100))
% title('qch output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(qch_HW(1:100))
% title('qch output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(5)
% subplot(211)
% plot(qch1(1:400))
% title('qch1 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(qch1_HW(1:400))
% title('qch1 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(6)
% subplot(211)
% plot(qch2(1:500))
% title('qch2 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(qch2_HW(1:500))
% title('qch2 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(7)
% subplot(211)
% plot(ich4(1:300))
% title('ich4 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(ich4_HW(1:300))
% title('ich4 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(8)
% subplot(211)
% plot(ich5(1:100))
% title('ich5 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(ich5_HW(1:100))
% title('ich5 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(9)
% subplot(211)
% plot(qch4(1:300))
% title('qch4 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(qch4_HW(1:300))
% title('qch4 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(10)
% subplot(211)
% plot(qch5(1:100))
% title('qch5 output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(qch5_HW(1:100))
% title('qch5 output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(11)
% subplot(211)
% plot(demodata(1:200))
% title('demodata output of 16QAM modulation after MATLAB software simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(demodata_HW(1:200))
% title('demodata output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% 
% figure(12)
% subplot(211)
% plot(noe_HW)
% title('Noe output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
% subplot(212)
% plot(nod_HW)
% title('Nod output of 16QAM modulation after System Gen hardware simulation','fontsize',18)
% xlabel('Number of Samples','fontsize',18);
% ylabel('Amplitude','fontsize',18);
 
    
    
    
    
    
    
    