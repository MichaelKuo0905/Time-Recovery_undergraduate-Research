% y=quantalph(x,alphabet)
%
% quantize the input signal x to the alphabet
% using nearest neighbor method
% input x - vector to be quantized
%       alphabet - vector of discrete values that y can take on
%                  sorted in ascending order
% output y - quantized vector

% function y=quantalph(x,alphabet)
%  y1=zeros([1,length(x)]);
% y2=zeros([1,length(x)]);
% y=zeros([1,length(x)]);
% for i=1:length(x)    
% if  real(x(i))>0.5
%     y1(i)=1;
% else 
%  y1(i)=-1;
% end
% end
% for  i=1:length(x)    
% if  imag(x(i))>0.5
%     y2(i)=1;
% else y2(i)=-1;
% 
% end
% end
% y=y1+j*y2;


function y=quantalph(x,alphabet)
alphabet=alphabet(:);
x=x(:);
alpha=alphabet(:,ones(size(x)))';
dist=(x(:,ones(size(alphabet)))-alpha).^2;
[v,i]=min(dist,[],2);
y=alphabet(i);
% 

