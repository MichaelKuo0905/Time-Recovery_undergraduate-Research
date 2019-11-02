
function varargout = untitled1(varargin)
% UNTITLED1 MATLAB code for untitled1.fig
%      UNTITLED1, by itself, creates a new UNTITLED1 or raises the existing
%      singleton*.
%
%      H = UNTITLED1 returns the handle to a new UNTITLED1 or the handle to
%      the existing singleton*.
%
%      UNTITLED1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED1.M with the given input arguments.
%
%      UNTITLED1('Property','Value',...) creates a new UNTITLED1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled1

% Last Modified by GUIDE v2.5 29-Dec-2015 19:29:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled1_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled1_OutputFcn, ...
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


% --- Executes just before untitled1 is made visible.
function untitled1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled1 (see VARARGIN)

% Choose default command line output for untitled1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = untitled1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
backgroundImage = importdata('1.bmp');
 axes(handles.axes8);
image(backgroundImage);
axis off



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function OFDMSymbol_Callback(hObject, eventdata, handles)
% hObject    handle to OFDMSymbol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OFDMSymbol as text
%        str2double(get(hObject,'String')) returns contents of OFDMSymbol as a double


% --- Executes during object creation, after setting all properties.
function OFDMSymbol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OFDMSymbol (see GCBO)
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


% --- Executes on selection change in IPOINT.
function IPOINT_Callback(hObject, eventdata, handles)
% hObject    handle to IPOINT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IPOINT contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IPOINT


% --- Executes during object creation, after setting all properties.
function IPOINT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IPOINT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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



function OFDMBW_Callback(hObject, eventdata, handles)
% hObject    handle to OFDMBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of OFDMBW as text
%        str2double(get(hObject,'String')) returns contents of OFDMBW as a double


% --- Executes during object creation, after setting all properties.
function OFDMBW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OFDMBW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function FallingEdgePosition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FallingEdgePosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function CoarseCFOOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoarseCFOOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function FineCFOOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FineCFOOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function LTPeakPosition_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LTPeakPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function EVM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function Fc_Callback(hObject, eventdata, handles)
% hObject    handle to Fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fc as text
%        str2double(get(hObject,'String')) returns contents of Fc as a double


% --- Executes during object creation, after setting all properties.
function Fc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bw_Callback(hObject, eventdata, handles)
% hObject    handle to Bw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bw as text
%        str2double(get(hObject,'String')) returns contents of Bw as a double


% --- Executes during object creation, after setting all properties.
function Bw_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fs_Callback(hObject, eventdata, handles)
% hObject    handle to Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fs as text
%        str2double(get(hObject,'String')) returns contents of Fs as a double


% --- Executes during object creation, after setting all properties.
function Fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Gain_Callback(hObject, eventdata, handles)
% hObject    handle to Gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Gain as text
%        str2double(get(hObject,'String')) returns contents of Gain as a double


% --- Executes during object creation, after setting all properties.
function Gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 N=str2double(get(handles.NumberSymbol,'string'));
 ml=(get(handles.BitperSymbol,'value'));
 beta=str2double(get(handles.SRRCRolloffFactor,'string'));
 OverSamp=(get(handles.OverSamp,'value'));
 Q_Delay_Tb=(get(handles.Q_Delay_Tb,'value'));
 clock_drift=(get(handles.ClockDrift,'value'));
 l=str2double(get(handles.DelayPulseShape,'string'));
 adaptive_method=(get(handles.AdaptiveMethod,'value'));
 mu=str2double(get(handles.StepSize,'string'));
 fc=str2double(get(handles.Fc,'string'));
 fs=str2double(get(handles.Fs,'string'));
%% ***************************************************************************************************
%%                                                           1. Overall Simulation Setup
%---------------- Channel and Clock Offset Effect ---------------
% %*********************E4406A Instrument Setting****************
%     [sig_bb_1,sig_IQ,fs]=E4406A_control;
%   save sig_bb_1;            
%   load sig_bb;                  % 8PSK ISI
%  load sig_bb2;                               
%   load sig_bb1;                   % OQPSK ISI               
%  load sig_bb3                 
%     load sig_bb4                  
%     load sig_bb4_1                 
RX_Input=resample(sig_bb_1,fs,15e6);
%---------------- find  packet position----------------------
a=conv(RX_Input,train(100*OverSamp:-1:1));
axes(handles.axes1);
plot(abs(a));xlabel('Time'); ylabel('Amplitude')
[c,b]=sort(abs(a),2,'descend');
%-------  in case of find untire packet ----------------------------
for i=1
    if  b(1,i) > length(RX_Input)-N*OverSamp
        i=i+1;
    end
end
peak_position=b(1,i);
set(handles.peak_position,'String',num2str(PeakPosition))
cor=train.*conj(RX_Input(peak_position-length(train)+1:peak_position));
ph_fix=angle(sum(cor)/length(train));                         %find phase  
set(handles.ph_fis,'String',num2str(PhaseOffset))
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
%---------------eye diagram----------------------------------------------
n=real(x);
neye=5;
c=floor(length(n)/(neye*OverSamp));
xp=n(4:end);  % dont plot transients at start
q=reshape(xp,neye*OverSamp,c);      % plot in clusters of size 5*Mt=(1:198)/50+1;
axes(handles.axes2)
plot(q)
axis([0,10,-3,3])
ylabel('Amplitude'), xlabel('symbols')
title('Eye diagram for sinc pulse shape')

%------------------------adaptive method-----------------------------------
switch(adaptive_method)   
 case 1
%---------------- Clock Recovery Using the Minimizing Cluster Variance Method ---------------
%-------------------------Real PART--------------------------------------
tnow1=l*OverSamp+1; tau1=0; xs1=zeros(1,N);   % initialize variables
tausave1=zeros(1,N); tausave1(1)=tau1; i=0;
% mu=0.3;                            % algorithm stepsize
delta=0.3;                          % time for derivative
while tnow1<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs1(i)=interpsinc(x1,tnow1+tau1,l,beta);   % interp at tnow+tau
  x_deltap1=interpsinc(x1,tnow1+tau1+delta,l,beta);  % value to right
  x_deltam1=interpsinc(x1,tnow1+tau1-delta,l,beta);  % value to left
  dx1=x_deltap1-x_deltam1; % numerical derivative
  
  qx1=quantalph1(xs1(i)*sqrt(10),ml);           % quantize to alphabet
  tau1=tau1+mu*dx1*(qx1-xs1(i));              % alg update (energy)
  tnow1=tnow1+OverSamp; tausave1(i)=tau1;      % save for plotting
end
%-----------------------------ImagPART--------------------------------------
tnow2=l*OverSamp+1; tau2=0; xs2=zeros(1,N);   % initialize variables
tausave2=zeros(1,N); tausave2(1)=tau2; i=0;
% mu=0.3;                            % algorithm stepsize
delta=0.1;                          % time for derivative
while tnow2<length(x)-2*l*OverSamp % run iteration
    
  i=i+1;
  xs2(i)=interpsinc(x2,tnow2+tau2,l,beta);   % interp at tnow+tau
  x_deltap2=interpsinc(x2,tnow2+tau2+delta,l,beta);  % value to right
  x_deltam2=interpsinc(x2,tnow2+tau2-delta,l,beta);  % value to left
  dx2=x_deltap2-x_deltam2; % numerical derivative
  
  qx2=quantalph1(xs2(i)*sqrt(10),ml);  % quantize to alphabet
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
plot(q)
 ylabel('Amplitude'), xlabel('symbols')
title('Eye diagram for sinc pulse shape')


xs=xs(3000:8000);
axes(handles.axes4)
subplot(2,1,1), plot(real(xs),'b.')     % plot constellation diagram
title('constellation diagram');
subplot(2,1,2), plot(imag(xs),'b.')      % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
axes(handles.axes5)
 plot(tausave(1:i-2))        % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

axes(handles.axes6)
plot(x,'o')        
ylabel('Amplitude'), xlabel('iterations')
axes(handles.axes7)
plot(xs,'o')        
ylabel('Amplitude'), xlabel('iterations')
% --- Executes on selection change in CorrectCFO.
function CorrectCFO_Callback(hObject, eventdata, handles)
% hObject    handle to CorrectCFO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns CorrectCFO contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CorrectCFO


% --- Executes during object creation, after setting all properties.
function CorrectCFO_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CorrectCFO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in UseMaxCombing.
function UseMaxCombing_Callback(hObject, eventdata, handles)
% hObject    handle to UseMaxCombing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns UseMaxCombing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from UseMaxCombing


% --- Executes during object creation, after setting all properties.
function UseMaxCombing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UseMaxCombing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SampleAdvance_Callback(hObject, eventdata, handles)
% hObject    handle to SampleAdvance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SampleAdvance as text
%        str2double(get(hObject,'String')) returns contents of SampleAdvance as a double


% --- Executes during object creation, after setting all properties.
function SampleAdvance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SampleAdvance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
