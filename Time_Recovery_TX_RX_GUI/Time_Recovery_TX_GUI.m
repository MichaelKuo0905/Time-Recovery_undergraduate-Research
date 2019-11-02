function varargout = Time_Recovery_TX_GUI(varargin)
% TIME_RECOVERY_TX_GUI MATLAB code for Time_Recovery_TX_GUI.fig
%      TIME_RECOVERY_TX_GUI, by itself, creates a new TIME_RECOVERY_TX_GUI or raises the existing
%      singleton*.
%
%      H = TIME_RECOVERY_TX_GUI returns the handle to a new TIME_RECOVERY_TX_GUI or the handle to
%      the existing singleton*.
%
%      TIME_RECOVERY_TX_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIME_RECOVERY_TX_GUI.M with the given input arguments.
%
%      TIME_RECOVERY_TX_GUI('Property','Value',...) creates a new TIME_RECOVERY_TX_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Time_Recovery_TX_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Time_Recovery_TX_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Time_Recovery_TX_GUI

% Last Modified by GUIDE v2.5 11-Aug-2016 23:12:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Time_Recovery_TX_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Time_Recovery_TX_GUI_OutputFcn, ...
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


% --- Executes just before Time_Recovery_TX_GUI is made visible.
function Time_Recovery_TX_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Time_Recovery_TX_GUI (see VARARGIN)

% Choose default command line output for Time_Recovery_TX_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%icon
% UIWAIT makes Time_Recovery_TX_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Time_Recovery_TX_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
 backgroundImage = importdata('1.bmp');
 axes(handles.axes4);
 image(backgroundImage);
 axis off


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


% --- Executes on selection change in SRRCRolloffFactor.
function SRRCRolloffFactor_Callback(hObject, eventdata, handles)
% hObject    handle to SRRCRolloffFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SRRCRolloffFactor contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SRRCRolloffFactor


% --- Executes during object creation, after setting all properties.
function SRRCRolloffFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SRRCRolloffFactor (see GCBO)
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



function PropagationDelay_Callback(hObject, eventdata, handles)
% hObject    handle to PropagationDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PropagationDelay as text
%        str2double(get(hObject,'String')) returns contents of PropagationDelay as a double


% --- Executes during object creation, after setting all properties.
function PropagationDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PropagationDelay (see GCBO)
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


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function TimeOffset_Callback(hObject, eventdata, handles)
% hObject    handle to TimeOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TimeOffset as text
%        str2double(get(hObject,'String')) returns contents of TimeOffset as a double


% --- Executes during object creation, after setting all properties.
function TimeOffset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Multipath.
function Multipath_Callback(hObject, eventdata, handles)
% hObject    handle to Multipath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Multipath contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Multipath


% --- Executes during object creation, after setting all properties.
function Multipath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Multipath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 N=str2double(get(handles.NumberSymbol,'string'));
 ml=(get(handles.BitperSymbol,'value'));
 beta=str2double(get(handles.SRRCRolloffFactor,'string'));
 OverSamp=str2double(get(handles.OverSamp,'string'));
 Q_Delay_Tb=(get(handles.Q_Delay_Tb,'value'));
 Multipath=(get(handles.Multipath,'value'));
 Delay=str2double(get(handles.PropagationDelay,'string'));
 l=str2double(get(handles.DelayPulseShape,'string'));
 toffset=str2double(get(handles.TimeOffset,'string'));   
fc=str2double(get(handles.Fc,'string'));
fs=str2double(get(handles.Fs,'string'));
%% ***************************************************************************************************
%%                                                           1. Overall Simulation Setup
%---------------- Channel and Clock Offset Effect ---------------
Multipath=Multipath-1;
Q_Delay_Tb=Q_Delay_Tb-1;
switch ml
    case 1;
    ml=ml+1;         
    case 2;
    ml=ml+1;
   case 3;
    ml=ml+1;
    otherwise
    ml=ml+2;
    end

if Multipath==1
chan=[1 0.3 0 0];  % T/m "channel"
else chan=[1 0 0 0];
end
 save chan chan
toffset=-0.3;                    % initial timing offset
pulshap_tx=srrc(l,beta,OverSamp,toffset);% srrc pulse shape with timing offset
hh=conv(pulshap_tx,chan);           % ... and pulse shape
%-----------------    data generation -----------------------
% nloop=1;                 % Number of simulation loops
% noe=0;                     % Number of error data
% nod=0;                     % Number of TX data
current_data=randn(1,N*ml)>0;
% save current_data(ml=2);
data_seq=Mapper(current_data,ml);  
m=data_seq;
% save m m;
trainmessage=sign(randn(1,100));


%---------------- TX Signal -------------------------------------
% trainlen=100;
% train=zeros(1,trainlen*OverSamp);
% train(1:OverSamp:trainlen*OverSamp)=trainmessage;      % train sequence oversample
%  save train train;
load train;
sup=zeros(1,N*OverSamp);                   % upsample the data by placing...
sup(1:OverSamp:N*OverSamp)=m;                  % ... M-1 zeros between each data point
if Q_Delay_Tb==1
 ich1=[real(sup) zeros(1,OverSamp/2)];    % note Delay 1Tb
 qch1=[zeros(1,OverSamp/2) imag(sup)];    % note Delay 1Tb
 sup=ich1+j*qch1;
end
sup=[train sup];
r=conv(hh,sup);                     % ... to get received signal
matchfilt=srrc(l,beta,OverSamp,0);  % matched filter = srrc pulse shape
x=conv(r,matchfilt);                % convolve signal with matched filter
packet_length=length(x)-length(train(100*OverSamp:-1:1));
save packet_length packet_length;
%    x=[train x];                        % add train sequence

TX_Output=[zeros(1,Delay) x zeros(1,Delay)];
RF_Output=TX_Output/max(abs(TX_Output));
%  save RF_Output;

% ---------------------- For E4438C RF TX Setting ----------------------
% fc = 0.7e9;                  %   Unit : MHz
% fs=300e4;                       
[status] = E4438C_control(RF_Output, fs, fc,'192.168.0.8')         
ich_TX=real(RF_Output);
qch_TX=imag(RF_Output);

axes(handles.axes2);
plot(ich_TX);xlabel('Time'); ylabel('Amplitude') 
axes(handles.axes3);
plot(qch_TX);xlabel('Time'); ylabel('Amplitude') 
Ts=1/fs; 
% Ts=1/fs2;
N=length(x);                               % length of the signal x
% t=Ts*(1:N);                                % define a time vector
ssf=(ceil(-N/2):ceil(N/2)-1)/(Ts*N);       % frequency vector
fx=fft(x(1:N));                            % do DFT/FFT
fxs=fftshift(fx);                          % shift it for plotting
axes(handles.axes1);
plot(ssf,20*log10(abs(fxs)));  % plot magnitude spectrum (dB)
% plotspec(real(RF_Output),Ts) % plot magnitude spectrum (dB)
% plot(ssf,abs(fxs))         % plot magnitude spectrum
xlabel('Frequency'); ylabel('Magnitude (dB)')  % label the axes



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


% --- Executes during object creation, after setting all properties.
