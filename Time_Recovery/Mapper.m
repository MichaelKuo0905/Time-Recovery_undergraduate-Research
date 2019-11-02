function [ OutputSymbols ] = Mapper( InputBits,BitsPerSymbol )
% %BitsPerSymbol:1,2,4,6,3-->BPSK,QPSK,16QAM,64QAM,MPSK
% N=10000;
% InputBits=randn(1,N)>0.5;
% BitsPerSymbol=3;
persistent BPSK_LUT QPSK_LUT QAM16_LUT QAM64_LUT MPSK_LUT 
if(isempty(BPSK_LUT))
    BPSK_LUT=[-1;1];
    QPSK_LUT=[-1;1]/sqrt(2);
    QAM16_LUT=[-3;-1;3;1]/sqrt(10);
    QAM64_LUT=[-7;-5;-1;-3;7;5;1;3]/sqrt(42);
    MPSK_LUT=[-1/sqrt(2)-j/sqrt(2);-1;+j;-1/sqrt(2)+j/sqrt(2);-j;1/sqrt(2)-j/sqrt(2);1/sqrt(2)+j/sqrt(2);1];
%     PSK16_LUT=[]
end
NumberOfSymbols=floor(length(InputBits)/BitsPerSymbol);
OutputSymbols=zeros(1,NumberOfSymbols);
for i=1:NumberOfSymbols
    Start=1+(i-1)*BitsPerSymbol;
    Stop=Start+BitsPerSymbol-1;
    BitGroup=InputBits(1,Start:Stop);
    switch(BitsPerSymbol)
        case 1
            Symbol=BPSK_LUT(BitGroup(1,1)+1,1);
        case 2
            Symbol=QPSK_LUT(BitGroup(1,1)+1,1)+...
                j*QPSK_LUT(BitGroup(1,2)+1,1);
        case 4
            Symbol=QAM16_LUT(BitGroup(1,1)*2+BitGroup(1,2)+1,1)+...
                j*QAM16_LUT(BitGroup(1,3)*2+BitGroup(1,4)+1,1);
        case 6
            Symbol=QAM64_LUT(BitGroup(1,1)*4+BitGroup(1,2)*2+BitGroup(1,3)+1,1)+...
                j*QAM64_LUT(BitGroup(1,4)*4+BitGroup(1,5)*2+BitGroup(1,6)+1,1);
        case 3 
             Symbol=MPSK_LUT(BitGroup(1,1)*4+BitGroup(1,2)*2+BitGroup(1,3)+1,1);
    end
    OutputSymbols(1,i)=Symbol;
end
