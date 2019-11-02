function [sig_bb,sig_IQ,fs]=E4406A_control

% ===========================================
%   Agilent E4406A VSA Control and Signal Acqusition
%   Programmer: Jeng-Kuang Hwang
%   Last update: Nov. 19, 2013
%   The IP right of this program is of YZU CSP LAB. Do not distribute!!
%===========================================

%% ========= Connect the E4406A to MATLAB==============
t = tcpip('192.168.0.6', 5025);
set(t,'InputBufferSize',4000000)      % Set instrument buffer size 3000
set(t,'ByteOrder','littleEndian')   % Set the instrument's byte order is little-endian
fopen(t)                            % connect to the instrument IP
fprintf(t,'*IDN?');
idn = fscanf(t)

%======== Setting the acqusition parameters============
fprintf(t,'%s',[':INST:SEL BASIC' 13 10]);                % Set the analyzer in Basic mode
fprintf(t,'%s',[':FREQ:CENTER 700000000Hz' 13 10]);             % Set center frequency
fprintf(t,'%s',[':WAV:BWID:RES 10MHz' 13 10]);            % Set RBW
fprintf(t,'%s',[':WAV:SWE:TIME 30ms' 13 10]);             % Set sweep time
fprintf(t,'%s',[':FORM REAL,32' 13 10]);              % Set the returned data type to REAL, 32
fprintf(t,'%s',[':FORM:BORD SWAP' 13 10]);             % Set the Byte order to Swap
fprintf(t,'%s',[':WAV:BAND:TYPE FLAT' 13 10]);             % Set filter to Flattop
% 
% fprintf(t,'%s',[':CAL:TCOR ON' 13 10]);             % Turn time corrections on
% fprintf(t,'%s',['POW:RF:ATT 0' 13 10]);             % Set attenuation
% fprintf(t,'%s',[':WAV:ACQ:PACK auto' 13 10]);             % Set data packing if necessary (Only needed for RBWs between 1.2 and 7.5 MHz where you need
%      

%======= Read the binary data from the E4406A===========
fprintf(t,'%s',[':INIT:CONT 1' 13 10]);             % Set analyzer to single sweep
fprintf(t,'%s',[':READ:WAV0?' 13 10]);             % Returns the I/Q data

% more dynamic range, although this may not be recommended)
head=char(fread(t,1));
if head~='#'
    disp('Read data error: beginning char is not a #  ');
    return;
end
no_char=str2num(char(fread(t,1)));
bytes_no=str2num(char(fread(t,no_char).'));
sig_length=bytes_no/(2*4);
sig_IQ=fread(t,[2 sig_length],'float32');
sig_bb=sig_IQ(1,:)+j*sig_IQ(2,:);      % complex envelope

%======== read the sampling time==========
fprintf(t,'%s',[':READ:WAV1?' 13 10]);
no_char=get(t,'BytesAvailable');
head2=char(fread(t,5)).';
Ts=fread(t,1,'float32');
fs=1/Ts                % sampling rate

% % % %====== Plotting the result=========
% % % figure(1);
% % % plot(sig_IQ(1,:),sig_IQ(2,:));
% % % title('The IQ plot of the captured E4406A raw data');
% % % figure(2);
% % % subplot(211);
% % % tx=(0:sig_length-1)*Ts;
% % % plot(tx,sig_IQ(1,:),tx,sig_IQ(2,:));
% % % title('The BB waveforms of the captured IQ data');
% % % xlabel('Time in sec');
% % % subplot(212);
% % % fx=((0:sig_length-1)/sig_length-0.5)*fs;
% % % semilogy(fx,fftshift(abs(fft(sig_bb))));
% % % title('The BB spectrum of the captured IQ data');
% % % xlabel('Freq in Hz');

fclose(t);
return;