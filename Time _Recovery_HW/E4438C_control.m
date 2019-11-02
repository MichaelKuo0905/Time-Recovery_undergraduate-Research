function [status]=E4438C_control(IQData, Fs, fc,IP)

% % % define a marker matrix and activate a marker to indicate the beginning of the waveform
Markers = zeros(2,length(IQData)); % fill marker array with zero, i.e no markers set
Markers(1,1) = 1; % set marker to first point of playback

%% Open a connection session with the signal generator. 
% make a new connection to theAgilent MXG/PSG over the GPIB interface
% io = agt_newconnection('gpib',0,19);
% IP = '192.168.0.8';
io = agt_newconnection('tcpip',IP);

% verify that communication with the Agilent MXG/PSG has been established
[status, status_description, query_result] = agt_query(io,'*idn?');
if (status < 0) ,return; end

% set the carrier frequency and power level on the Agilent MXG/PSG using the Agilent
% Waveform Download Assistant
fc = num2str(fc);                               % Central frquency(MHz)
[status, status_description] = agt_sendcommand(io,['SOURce:FREQuency' ' ' fc]);
[status, status_description] = agt_sendcommand(io, 'POWer -30');
[status, status_description] = agt_sendcommand(io, 'Mod ON');
% % % define the ARB sample clock for playback
sampclk = Fs;

%% Download the I/Q data
% download the iq waveform to the PSG baseband generator for playback
[status, status_description] = agt_waveformload(io, IQData, 'LTE_UL', sampclk, 'play', 'no_normscale',Markers);

% turn on Modulation output 
[status, status_description ] = agt_sendcommand( io, ':OUTPut:MODulation:STATe ON' )

% turn on RF output power
[status, status_description ] = agt_sendcommand( io, 'OUTPut:STATe ON' )

% % plot(real(IQData))

%% Release memory
agt_closeAllSessions;