function [ OutputSymbols ] = Mapper( InputBits,BitsPerSymbol )
%BitsPerSymbol:1,2,4,6-->BPSK,QPSK,16QAM,64QAM
persistent BPSK_LUT QPSK_LUT QAM16_LUT QAM64_LUT
if(isempty(BPSK_LUT))
    BPSK_LUT=[-1;1];
    QPSK_LUT=[-1;1]/sqrt(2);
    QAM16_LUT=[-3;-1;1;3]/sqrt(10);
    QAM64_LUT=[-7;-5;-3;-1;1;3;5;7]/sqrt(42);
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
    end
    OutputSymbols(1,i)=Symbol;
end