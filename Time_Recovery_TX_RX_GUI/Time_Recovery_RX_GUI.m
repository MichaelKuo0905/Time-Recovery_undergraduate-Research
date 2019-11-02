function varargout = Time_Recovery_RX_GUI(varargin)
% TIME_RECOVERY_RX_GUI MATLAB code for Time_Recovery_RX_GUI.fig
%      TIME_RECOVERY_RX_GUI, by itself, creates a new TIME_RECOVERY_RX_GUI or raises the existing
%      singleton*.
%
%      H = TIME_RECOVERY_RX_GUI returns the handle to a new TIME_RECOVERY_RX_GUI or the handle to
%      the existing singleton*.
%
%      TIME_RECOVERY_RX_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIME_RECOVERY_RX_GUI.M with the given input arguments.
%
%      TIME_RECOVERY_RX_GUI('Property','Value',...) creates a new TIME_RECOVERY_RX_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Time_Recovery_RX_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Time_Recovery_RX_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Time_Recovery_RX_GUI

% Last Modified by GUIDE v2.5 09-Aug-2016 13:38:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Time_Recovery_RX_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Time_Recovery_RX_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before Time_Recovery_RX_GUI is made visible.
function Time_Recovery_RX_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Time_Recovery_RX_GUI (see VARARGIN)

% Choose default command line output for Time_Recovery_RX_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Time_Recovery_RX_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Time_Recovery_RX_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
N=str2double(get(handles.NumberSymbol,'string'));
ml=(get(handles.BitperSymbol,'value'));
beta=str2double(get(handles.SRRCRolloffFactor,'string'));
 OverSamp=str2double(get(handles.OverSamp,'string'));
 Q_Delay_Tb=(get(handles.Q_Delay_Tb,'value'));
 clock_drift=(get(handles.ClockDrift,'value'));
 l=str2double(get(handles.DelayPulseShape,'string'));
 adaptive_method=(get(handles.AdaptiveMethod,'value'));
 mu=str2double(get(handles.StepSize,'string'));
 average_power=(get(handles.AveragePower,'value'));
%% ***************************************************************************************************
%%                                                           1. Overall Simulation Setup
Q_Delay_Tb=Q_Delay_Tb-1;
clock_drift=clock_drift-1;
% save Q_Delay_Tb Q_Delay_Tb
% save clock_drift clock_drift
% save adaptive_method adaptive_method
% save mu mu
% save average_power average_power
% save N N

load train
load packet_length
switch average_power
    case 1
        average_power=sqrt(2);
    case  2 
        average_power=sqrt(10);
    case 3 
        average_power=sqrt(42);
end

switch ml
    case 1;
    ml=ml+1;         
    case 2;
    ml=ml+1;
   case 3;
    ml=ml+1;
   case 4;
    ml=ml+2;
end
  %*********************E4406A Instrument Setting****************
             [sig_bb,sig_IQ,fs]=E4406A_control;
 %                  save sig_bb3_1 sig_bb3_1;                       
%         load sig_bb;                  %QPSK ISI   
%         load sig_bb1;                  %QPSK ideal   
 %       load sig_bb_1;                %OQPSK ISI                                  
%       load sig_bb2;                  % 8PSK ISI                        
%            load sig_bb3                      %16QAM Ideal               
%              load sig_bb4                    %64QAM   mistake                    
%           load sig_bb51                          
 RX_Input=resample(sig_bb,300e4,15e6);
%---------------- find  packet position----------------------
a=conv(RX_Input,train(100*OverSamp:-1:1));
% figure(1)
 axes(handles.axes11);
plot(abs(a));
[c,b]=sort(abs(a),2,'descend');
%-------  in case of find untire packet ----------------------------
 i=1;
while b(1,i) > length(RX_Input)-N*OverSamp
     
        i=i+1;
    end

peak_position=b(1,i);
set(handles.PeakPosition,'String',num2str(peak_position));
cor=train.*conj(RX_Input(peak_position-length(train)+1:peak_position));
ph_fix=angle(sum(cor)/length(train));                  %find phase   
set(handles.PhaseOffset,'String',num2str(ph_fix))
RX_Input=RX_Input(b(1,i)+1:b(1,i)+packet_length);             % change to one packet
% RX_Input=RX_Input(21845:61840);
RX_Output=RX_Input/max(abs(RX_Input));
% figure(2)
% plot(RX_Output,'g');
 RX_Output=RX_Output*exp(j*ph_fix);
% figure(3),plot(RX_Output,'g')
 x=RX_Output;
%---------------- Resample to change clock period for clock linear drift -------------
if clock_drift==1
fac=1.0001; z=zeros(size(x)); % percent change in period
t=l+1:fac:length(x)-2*l;      % vector of new times
for i=l+1:length(t)           % resample x at new rate
  z(i)=interpsinc(x,t(i),l);  % to create received signal
end                           % with period offset

if  Q_Delay_Tb==1;              %OQPSK Signal
  x1=[zeros(1,OverSamp/2) real(z(1:end-OverSamp/2))]; % 實部虛部對齊
  x2=imag(z);
  x=x1+j*x2;                                        
else
   x=z;                          % relabel signal
end
else
    if  Q_Delay_Tb==1;              %OQPSK Signal
  x1=[zeros(1,OverSamp/2) real(x(1:end-OverSamp/2))]; % 實部虛部對齊
  x2=imag(x);
  x=x1+j*x2;                                        
else
   x=x;                          % relabel signal
    end
end
% figure(2)
% plot(x,'o')
% %---------------eye diagram----------------------------------------------
n=real(x(30000:35000));
neye=5;
c=floor(length(n)/(neye*OverSamp));
xp=n(2:end);  %6:OQPSK  4:otherwise            % dont plot transients at start
q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
 axes(handles.axes2);
plot(q)
axis([1,20,-1,1])
ylabel('Amplitude'), xlabel('symbols')
title('Eye diagram for sinc pulse shape')

%------------------------adaptive method-----------------------------------
switch(adaptive_method)   
 case 1
x1=real(x);
x2=imag(x);
    
%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
%-------------------------Real PART--------------------------------------
tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau1; i=0;
% mu=0.9;                            % algorithm stepsize
delta=0.3;                          % time for derivative
while tnow1<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs1(i)=interpsinc(x1,tnow1+tau1,l,beta);   % interp at tnow+tau
  x_deltap1=interpsinc(x1,tnow1+tau1+delta,l,beta);  % value to right
  x_deltam1=interpsinc(x1,tnow1+tau1-delta,l,beta);  % value to left
  dx1=x_deltap1-x_deltam1; % numerical derivative
  
  qx1=quantalph1(xs1(i)*sqrt(42),ml);           % quantize to alphabet
  tau1=tau1+mu*dx1*(qx1-xs1(i));              % alg update (energy)
  tnow1=tnow1+OverSamp; tausave(i)=tau1;      % save for plotting
end
%-----------------------------ImagPART--------------------------------------
tnow2=l*OverSamp+1; tau2=0; xs2=zeros(1,N);   % initialize variables
tausave2=zeros(1,N); tausave2(1)=tau2; i=0;
% mu=0.9;                            % algorithm stepsize
delta=0.3;                          % time for derivative
while tnow2<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs2(i)=interpsinc(x2,tnow2+tau2,l,beta);   % interp at tnow+tau
  x_deltap2=interpsinc(x2,tnow2+tau2+delta,l,beta);  % value to right
  x_deltam2=interpsinc(x2,tnow2+tau2-delta,l,beta);  % value to left
  dx2=x_deltap2-x_deltam2; % numerical derivative
  
  qx2=quantalph1(xs2(i)*sqrt(42),ml);  % quantize to alphabet
  tau2=tau2+mu*dx2*(qx2-xs2(i));              % alg update (energy)
  tnow2=tnow2+OverSamp; tausave2(i)=tau2;      % save for plotting
end
xs=xs1+j*xs2;
 case 2
 %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
% mu=0.2;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
 

end
 %---------------- Clock Recovery Using the Absolute Sampler Method ---------------   
 case 3
tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
tausave=zeros(1,N); tausave(1)=tau; i=0;
% mu=0.2;                            % algorithm stepsize  4QAM(ISI):0.02     
delta=0.1;                          % time for derivative
while tnow<length(x)-l*OverSamp            % run iteration
  i=i+1;
  xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
  x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
  x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
  dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
  tau=tau+mu*dx;              % alg update (energy)
  tnow=tnow+OverSamp; 
  tausave(i)=tau;      % save for plotting
end
end 
%---------------eye diagram----------------------------------------------
b=real(xs(3000:8000));
neye=5;
c=floor(length(b)/neye);
xp=b(2:end);  % dont plot transients at start
q=reshape(xp,neye,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
axes(handles.axes3)
plot(q);ylabel('Amplitude'),xlabel('symbols')



 xs=xs(3000:8000);
axes(handles.axes4)
plot(real(xs),'b.')     % plot constellation diagram
 title('constellation diagram');
ylabel('Amplitude')
axes(handles.axes5)
 plot(imag(xs),'b.')      % plot constellation diagram
ylabel('Amplitude')
 axes(handles.axes6)
 plot(tausave(1:i-2))        % plot trajectory of tau
 ylabel('offset estimates'), xlabel('iterations')

axes(handles.axes7)
plot(x,'o')        
ylabel('Amplitude'), xlabel('iterations')
axes(handles.axes8)
plot(xs,'o')        
ylabel('Amplitude'), xlabel('iterations')

%  N=str2double(get(handles.NumberSymbol,'string'));
%  ml=(get(handles.BitperSymbol,'value'));
%  beta=str2double(get(handles.SRRCRolloffFactor,'string'));
%  OverSamp=str2double(get(handles.OverSamp,'string'));
%  Q_Delay_Tb=(get(handles.Q_Delay_Tb,'value'));
%  clock_drift=(get(handles.ClockDrift,'value'));
%  l=str2double(get(handles.DelayPulseShape,'string'));
%  adaptive_method=(get(handles.AdaptiveMethod,'value'));
%  mu=str2double(get(handles.StepSize,'string'));
%  average_power=(get(handles.AveragePower,'value'));
% %  fc=str2double(get(handles.Fc,'string'));
% %  fs=str2double(get(handles.Fs,'string'));
% %% ***************************************************************************************************
% %%                                                           1. Overall Simulation Setup
% Q_Delay_Tb=Q_Delay_Tb-1;
% clock_drift=clock_drift-1;
% % save Q_Delay_Tb Q_Delay_Tb
% % save clock_drift clock_drift
% % save adaptive_method adaptive_method
% % save mu mu
% % save average_power average_power
% % save N N
% 
% load train
% load packet_length
% switch average_power
%     case 1
%         average_power=sqrt(2);
%     case  2 
%         average_power=sqrt(10);
%     case 3 
%         average_power=sqrt(42);
% end
% 
% switch ml
%     case 1;
%     ml=ml+1;         
%     case 2;
%     ml=ml+1;
%    case 3;
%     ml=ml+1;
%    case 4;
%     ml=ml+2;
%     end
% %    save average_power average_power     
% %---------------- Channel and Clock Offset Effect ---------------
% % %*********************E4406A Instrument Setting****************
% %             [sig_bb5,sig_IQ,fs]=E4406A_control;
% %                  save sig_bb3_1 sig_bb3_1;                       
% %       load sig_bb;                  %QPSK ISI   
% %       load sig_bb1;                  %QPSK ideal   
% %       load sig_bb_1;                %OQPSK ISI                                  
% %     load sig_bb2;                  % 8PSK ISI                        
%       load sig_bb3                      %16QAM Ideal               
% %         load sig_bb5                    %64QAM                      
% %       load sig_bb4                           
% RX_Input=resample(sig_bb3,300e4,15e6);
% %---------------- find  packet position----------------------
% a=conv(RX_Input,train(100*OverSamp:-1:1));
%  axes(handles.axes11);
%  plot(abs(a));xlabel('Time'), ylabel('Amplitude')
%  [c,b]=sort(abs(a),2,'descend');
% %-------  in case of find untire packet ----------------------------
% for i=1
%     if  b(1,i) > length(RX_Input)-N*OverSamp
%         i=i+1;
%     end
% end
% peak_position=b(1,i);
% set(handles.PeakPosition,'String',num2str(peak_position));
% cor=train.*conj(RX_Input(peak_position-length(train)+1:peak_position));
% ph_fix=angle(sum(cor)/length(train));                         %find phase  
% set(handles.PhaseOffset,'String',num2str(ph_fix))
% RX_Input=RX_Input(b(1,i)+1:b(1,i)+packet_length);             % change to one packet
% % RX_Input=RX_Input(21845:61840);
% RX_Output=RX_Input/max(abs(RX_Input));
% % figure(2)
% % plot(RX_Output,'g');
% RX_Output=RX_Output*exp(1i*ph_fix);
% % figure(3),plot(RX_Output,'g')
% x=RX_Output;
% %---------------- Resample to change clock period for clock linear drift -------------
% if clock_drift==1
% fac=1.0001; z=zeros(size(x)); % percent change in period
% t=l+1:fac:length(x)-2*l;      % vector of new times
% for i=l+1:length(t)           % resample x at new rate
%   z(i)=interpsinc(x,t(i),l);  % to create received signal
% end                           % with period offset
% 
% if  Q_Delay_Tb==1;              %OQPSK Signal
%   x1=[zeros(1,OverSamp/2) real(z(1:end-OverSamp/2))]; % 實部虛部對齊
%   x2=imag(z);
%   x=x1+j*x2;                                        
% else
%    x=z;                          % relabel signal
% end
% else
%     if  Q_Delay_Tb==1;              %OQPSK Signal
%   x1=[zeros(1,OverSamp/2) real(x(1:end-OverSamp/2))]; % 實部虛部對齊
%   x2=imag(x);
%   x=x1+j*x2;                                        
% else
%    x=x;                          % relabel signal
%     end
% end
% save x x;
% % %---------------eye diagram----------------------------------------------
% n=real(x(30000:35000));
% neye=5;
% c=floor(length(n)/(neye*OverSamp));
% xp=n(2:end);  %6:OQPSK  4:otherwise            % dont plot transients at start
% q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% axes(handles.axes2);
% plot(q)
% axis([1,20,-1,1])
% ylabel('Amplitude'), xlabel('symbols')
% title('Eye diagram for sinc pulse shape')
% 
% %------------------------adaptive method-----------------------------------
% switch(adaptive_method)   
%  case 1
%  x1=real(x);
%  x2=imag(x);
%     
% %---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
% %-------------------------Real PART--------------------------------------
% tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);   % initialize variables
% tausave=zeros(1,N); tausave(1)=tau1; i=0;
% % mu=0.3;                            % algorithm stepsize
% delta=0.3;                          % time for derivative
% while tnow1<length(x)-2*l*OverSamp % run iteration
%     
%   i=i+1;
%   xs1(i)=interpsinc(x1,tnow1+tau1,l,beta);   % interp at tnow+tau
%   x_deltap1=interpsinc(x1,tnow1+tau1+delta,l,beta);  % value to right
%   x_deltam1=interpsinc(x1,tnow1+tau1-delta,l,beta);  % value to left
%   dx1=x_deltap1-x_deltam1; % numerical derivative
%   
%   qx1=quantalph1(xs1(i)*average_power,ml);           % quantize to alphabet
%   tau1=tau1+mu*dx1*(qx1-xs1(i));              % alg update (energy)
%   tnow1=tnow1+OverSamp; tausave(i)=tau1;      % save for plotting
% end
% %-----------------------------ImagPART--------------------------------------
% tnow2=l*OverSamp+1; tau2=0; xs2=zeros(1,N);   % initialize variables
% tausave2=zeros(1,N); tausave2(1)=tau2; i=0;
% % mu=0.3;                            % algorithm stepsize
% delta=0.3;                          % time for derivative
% while tnow2<length(x)-2*l*OverSamp % run iteration
%     
%   i=i+1;
%   xs2(i)=interpsinc(x2,tnow2+tau2,l,beta);   % interp at tnow+tau
%   x_deltap2=interpsinc(x2,tnow2+tau2+delta,l,beta);  % value to right
%   x_deltam2=interpsinc(x2,tnow2+tau2-delta,l,beta);  % value to left
%   dx2=x_deltap2-x_deltam2; % numerical derivative
%   
%   qx2=quantalph1(xs2(i)*average_power,ml);  % quantize to alphabet
%   tau2=tau2+mu*dx2*(qx2-xs2(i));              % alg update (energy)
%   tnow2=tnow2+OverSamp; tausave2(i)=tau2;      % save for plotting
% end
% xs=xs1+j*xs2;
%  case 2
%  %---------------- Clock Recovery Using the Maximizing Output Power Method ---------------
% tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
% tausave=zeros(1,N); tausave(1)=tau; i=0; 
% % mu=0.2;                            % algorithm stepsize
% delta=0.1;                          % time for derivative
% while tnow<length(x)-l*OverSamp            % run iteration
%   i=i+1;
%   xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
%   dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
%   tau=tau-mu*dx*abs(xs(i))^5;              % alg update (energy)
%   tnow=tnow+OverSamp; 
%   tausave(i)=tau;      % save for plotting
%  
% 
% end
%  %---------------- Clock Recovery Using the Absolute Sampler Method ---------------   
%  case 3
% tnow=l*OverSamp+1; tau=0; xs=zeros(1,N);   % initialize variables
% tausave=zeros(1,N); tausave(1)=tau; i=0;
% % mu=0.2;                            % algorithm stepsize  4QAM(ISI):0.02  64QAM=0.15   
% delta=0.1;                          % time for derivative
% while tnow<length(x)-l*OverSamp            % run iteration
%   i=i+1;
%   xs(i)=interpsinc(x,tnow+tau,l,beta);   % interp at tnow+tau
%   x_deltap=interpsinc(x,tnow+tau+delta,l,beta);  % value to right
%   x_deltam=interpsinc(x,tnow+tau-delta,l,beta);  % value to left
%   dx=(abs(x_deltap)-abs(x_deltam))/delta;             % numerical derivative
%   tau=tau+mu*dx;              % alg update (energy)
%   tnow=tnow+OverSamp; 
%   tausave(i)=tau;      % save for plotting
% end
% end 
% %---------------eye diagram----------------------------------------------
% b=real(xs(3000:8000));
% neye=5;
% c=floor(length(b)/neye);
% xp=b(2:end);  % dont plot transients at start
% q=reshape(xp,neye,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
% axes(handles.axes3)
% plot(q);ylabel('Amplitude'),xlabel('symbols')
% 
% 
% 
% xs=xs(3000:8000);
% axes(handles.axes4)
% plot(real(xs),'b.')     % plot constellation diagram
%  title('constellation diagram');
% ylabel('Amplitude')
% axes(handles.axes5)
%  plot(imag(xs),'b.')      % plot constellation diagram
% ylabel('Amplitude')
%  axes(handles.axes6)
%  plot(tausave(1:i-2))        % plot trajectory of tau
%  ylabel('offset estimates'), xlabel('iterations')
% 
% axes(handles.axes7)
% plot(x,'o')        
% ylabel('Amplitude'), xlabel('iterations')
% axes(handles.axes8)
% plot(xs,'o')        
% ylabel('Amplitude'), xlabel('iterations')



function StepSize_Callback(hObject, eventdata, handles)
% hObject    handle to StepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StepSize as text
%        str2double(get(hObject,'String')) returns contents of StepSize as a double


% --- Executes during object creation, after setting all properties.
function StepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StepSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in AdaptiveMethod.
function AdaptiveMethod_Callback(hObject, eventdata, handles)
% hObject    handle to AdaptiveMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns AdaptiveMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AdaptiveMethod


% --- Executes during object creation, after setting all properties.
function AdaptiveMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AdaptiveMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ClockDrift.
function ClockDrift_Callback(hObject, eventdata, handles)
% hObject    handle to ClockDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ClockDrift contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClockDrift


% --- Executes during object creation, after setting all properties.
function ClockDrift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClockDrift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NumberSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to NumberSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumberSymbol as text
%        str2double(get(hObject,'String')) returns contents of NumberSymbol as a double


% --- Executes during object creation, after setting all properties.
function NumberSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumberSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BitperSymbol.
function BitperSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to BitperSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BitperSymbol contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BitperSymbol


% --- Executes during object creation, after setting all properties.
function BitperSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BitperSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function OverSamp_Callback(hObject, eventdata, handles)
% hObject    handle to OverSamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OverSamp as text
%        str2double(get(hObject,'String')) returns contents of OverSamp as a double


% --- Executes during object creation, after setting all properties.
function OverSamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OverSamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SRRCRolloffFactor_Callback(hObject, eventdata, handles)
% hObject    handle to SRRCRolloffFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SRRCRolloffFactor as text
%        str2double(get(hObject,'String')) returns contents of SRRCRolloffFactor as a double


% --- Executes during object creation, after setting all properties.
function SRRCRolloffFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SRRCRolloffFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Q_Delay_Tb.
function Q_Delay_Tb_Callback(hObject, eventdata, handles)
% hObject    handle to Q_Delay_Tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Q_Delay_Tb contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Q_Delay_Tb


% --- Executes during object creation, after setting all properties.
function Q_Delay_Tb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q_Delay_Tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DelayPulseShape_Callback(hObject, eventdata, handles)
% hObject    handle to DelayPulseShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DelayPulseShape as text
%        str2double(get(hObject,'String')) returns contents of DelayPulseShape as a double


% --- Executes during object creation, after setting all properties.
function DelayPulseShape_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DelayPulseShape (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AveragePower_Callback(hObject, eventdata, handles)
% hObject    handle to AveragePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AveragePower as text
%        str2double(get(hObject,'String')) returns contents of AveragePower as a double


% --- Executes during object creation, after setting all properties.
function AveragePower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AveragePower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function PhaseOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PhaseOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text151_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text151 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function PeakPosition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeakPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function PeakPosition_Callback(hObject, eventdata, handles)
% hObject    handle to PeakPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PeakPosition as text
%        str2double(get(hObject,'String')) returns contents of PeakPosition as a double
