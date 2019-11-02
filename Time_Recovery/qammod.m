% Program 3-23
% qammod.m
%
% This function is used for Gray coding of 16QAM modulation  
%
% programmed by R.Funada and H.Harada
%

function [iout,qout]=qammod(paradata,para,nd,ml)

%****************** variables *************************
% paradata : input data (para-by-nd matrix)
% iout :output Ich data
% qout :output Qch data
% para   : Number of paralell channels
% nd : Number of data
% ml : Number of modulation levels
% (QPSK ->2  16QAM -> 4)
% *****************************************************

% The constellation power

k=sqrt(10);
%iv=[3 1 -3 -1];
iv=[-3 -1 3 1];  % note
               
aa=reshape(paradata,2,length(paradata)/2);
bb=[2 1]*aa+1;
I_ch_idx=bb(1:2:end);
Q_ch_idx=bb(2:2:end);

iout=iv(I_ch_idx)/k;
qout=iv(Q_ch_idx)/k;


%******************** end of file ***************************

      
 
