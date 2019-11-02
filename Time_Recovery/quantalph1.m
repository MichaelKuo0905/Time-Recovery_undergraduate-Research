% y=quantalph(x,alphabet) 
% quantize the input signal x to the alphabet
% using nearest neighbor method
% input x - vector to be quantized
% alphabet - vector of discrete values that y can take on
%  sorted in ascending order
% output y - quantized vector

function y=quantalph1(x,ml)
 y1=zeros([1,length(x)]);
% ml=4;
% x=zeros(1,5);
% x=[-2.5,1.5,-0.1,0.1,1]
switch(ml) 
    case 2
        for i=1:length(x)    
        if  real(x(i))>0
         y1(i)=1/sqrt(2);
        else 
        y1(i)=-1/sqrt(2);
        end
        end
    case 4
           for i=1:length(x)
            if  real(x(i))>0
                if   real(x(i))>2
                    y1(i)=3/sqrt(10);
                else  y1(i)=1/sqrt(10);
                end
            else   
                 if real(x(i))<-2
                     y1(i)=-3/sqrt(10);
                 else
                     y1(i)=-1/sqrt(10);
                    
                 end
            end
           end
    case 3 
           for i=1:length(x)
              if real(x(i))>0
                  if real(x(i))>0.85
                      y1(i)=1;
                  else
                   if real(x(i))>0.35
                       y1(i)=1/sqrt(2);
                   else
                       y1(i)=0;
                   end
                   
                  end
              else
                  if real(x(i))<-0.85
                       y1(i)=-1;
                  else
                      if  real(x(i))<-0.35
                          y1(i)=-1/sqrt(2);
                      else
                          y1(i)=0;
                      end
                  end
              end
           end
                   
               
    case 6
            for i=1:length(x)
              if real(x(i))>0
                 if real(x(i))>6
                     y1(i)=7/sqrt(42);
                 else 
                     if real(x(i))>4
                         y1(i)=5/sqrt(42);
                     else 
                         if real(x(i))>2
                             y1(i)=3/sqrt(42);
                         else
                             y1(i)=1/sqrt(42);
                         end
                     end
                 end
              else  
                    if real(x(i))<-6      
                     y1(i)=-7/sqrt(42);
                 else 
                     if real(x(i))<-4
                         y1(i)=-5/sqrt(42);
                     else 
                         if real(x(i))<-2
                             y1(i)=-3/sqrt(42);
                         else
                             y1(i)=-1/sqrt(42);
                         end   
                         
                     end
                    end
                        
              end     
            end
end
        y=y1;
% for  i=1:length(x)    
% if  imag(x(i))>0.5
%     y2(i)=1/sqrt(2);
% else y2(i)=-1/sqrt(2);
% 
% end
% end


% 
% function y=quantalph(x,alphabet)
% alphabet=alphabet(:);
% x=x(:);
% alpha=alphabet(:,ones(size(x)))';
% dist=(x(:,ones(size(alphabet)))-alpha).^2;
% [v,i]=min(dist,[],2);
% y=alphabet(i);
% 

