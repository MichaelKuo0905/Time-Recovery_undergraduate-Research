function varargout = Time_Recovery_RX_GUI(varargin)
%TIME_RECOVERY_RX_GUI M-file for Time_Recovery_RX_GUI.fig
%      TIME_RECOVERY_RX_GUI, by itself, creates a new TIME_RECOVERY_RX_GUI or raises the existing
%      singleton*.
%
%      H = TIME_RECOVERY_RX_GUI returns the handle to a new TIME_RECOVERY_RX_GUI or the handle to
%      the existing singleton*.
%
%      TIME_RECOVERY_RX_GUI('Property','Value',...) creates a new TIME_RECOVERY_RX_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Time_Recovery_RX_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TIME_RECOVERY_RX_GUI('CALLBACK') and TIME_RECOVERY_RX_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TIME_RECOVERY_RX_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Time_Recovery_RX_GUI

% Last Modified by GUIDE v2.5 09-Aug-2016 09:49:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Time_Recovery_RX_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Time_Recovery_RX_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
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
