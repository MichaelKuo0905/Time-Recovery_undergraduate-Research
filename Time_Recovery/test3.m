% y=quantalph(x,alphabet)
%
% quantize the input signal x to the alphabet
% using nearest neighbor method
% input x - vector to be quantized
%       alphabet - vector of discrete values that y can take on
%                  sorted in ascending order
% output y - quantized vector

x=[0.5+1i, 0.8-0.5i];
y1=zeros([1,length(x)]);
y2=zeros([1,length(x)]);
for i=1:length(x)    
if  real(x(i))>0
    y1(i)=1;
else 
 y1(i)=-1;
end
end
for  i=1:length(x)    
if  imag(x(i))>0
    y2(i)=1;
else y2(i)=-1;

end
end
x=y2+j*y2;
% demodata((1:para),(1:ml:ml*nd-1))=idata((1:para),(1:nd))>=0; % note
% demodata((1:para),(2:ml:ml*nd))=qdata((1:para),(1:nd))>=0; % note
