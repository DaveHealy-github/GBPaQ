function varargout = GBPaQ(varargin)
% GBPAQ MATLAB code for GBPaQ.fig
%      GBPAQ, by itself, creates a new GBPAQ or raises the existing
%      singleton*.
%
%      H = GBPAQ returns the handle to a new GBPAQ or the handle to
%      the existing singleton*.
%
%      GBPAQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GBPAQ.M with the given input arguments.
%
%      GBPAQ('Property','Value',...) creates a new GBPAQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GBPaQ_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GBPaQ_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%% Copyright
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.


% Edit the above text to modify the response to help GBPaQ

% Last Modified by GUIDE v2.5 30-Sep-2021 13:11:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;

gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GBPaQ_OpeningFcn, ...
                   'gui_OutputFcn',  @GBPaQ_OutputFcn, ...
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


%   my initialisation can go here... 
% 
% --- Executes just before GBPaQ is made visible.
function GBPaQ_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GBPaQ (see VARARGIN)

% Choose default command line output for GBPaQ
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GBPaQ wait for user response (see UIRESUME)
% uiwait(handles.figure_main);
axes(handles.axes_logo) ;
imLogo = imread('GBPaQlogo.jpeg') ; 
imshow(imLogo) ; 

axes(handles.axes_icon) ;
imIcon = imread('GBPaQicon.jpeg') ; 
imshow(imIcon) ; 

axes(handles.axes_tracemap) ; 

set(handles.popupmenu_histolengthbins,'Value',2) ;
set(handles.popupmenu_histoanglebins,'Value',2) ;
set(handles.popupmenu_rosebins,'Value',2) ;

handles.GBPaQversion = '0.9' ; 
disp(['GBPaQ version ', handles.GBPaQversion]) ; 
% Update handles structure
guidata(hObject, handles);

v = ver('MATLAB') ; 
iRelease = regexp(v.Release, 'R.....') ; 
nReleaseYear = str2num(v.Release(iRelease+1:iRelease+4)) ; 
if nReleaseYear < 2015 
    h = msgbox('Warning: GBPaQ needs MATLAB release of R2015a, or later.', 'Important Warning!') ;  
    uiwait(h) ; 
end ; 

%Reposition each panel to same location as panel 1
set(handles.uibgTabLengths,'position',get(handles.uibgTabMaps,'position'));
set(handles.uibgTabAngles,'position',get(handles.uibgTabMaps,'position'));

set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 

% --- Outputs from this function are returned to the command line.
function varargout = GBPaQ_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_intensitymap.
function checkbox_intensitymap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_intensitymap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_intensitymap
if get(hObject, 'Value') || get(handles.checkbox_densitymap, 'Value')
    set(handles.checkbox_showcircles, 'Enable', 'on') ; 
    set(handles.edit_numscancircles, 'Enable', 'on') ; 
    set(handles.label_numscancircles, 'Enable', 'on') ; 
else 
    set(handles.checkbox_showcircles, 'Enable', 'off') ; 
    set(handles.checkbox_showcircles, 'Value', 0) ; 
    set(handles.edit_numscancircles, 'Enable', 'off') ; 
    set(handles.edit_numscancircles, 'Value', 0) ; 
    set(handles.label_numscancircles, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_densitymap.
function checkbox_densitymap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_densitymap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_densitymap
if get(hObject, 'Value') || get(handles.checkbox_intensitymap, 'Value')
    set(handles.checkbox_showcircles, 'Enable', 'on') ; 
    set(handles.edit_numscancircles, 'Enable', 'on') ; 
    set(handles.label_numscancircles, 'Enable', 'on') ; 
else 
    set(handles.checkbox_showcircles, 'Enable', 'off') ; 
    set(handles.checkbox_showcircles, 'Value', 0) ; 
    set(handles.edit_numscancircles, 'Enable', 'off') ; 
    set(handles.edit_numscancircles, 'Value', 0) ; 
    set(handles.label_numscancircles, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_tracemap.
function checkbox_tracemap_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tracemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tracemap
if get(hObject, 'Value') 
    set(handles.checkbox_shownodes, 'Enable', 'on') ; 
else 
    set(handles.checkbox_shownodes, 'Enable', 'off') ; 
    set(handles.checkbox_shownodes, 'Value', 0) ; 
end ; 
    
% --- Executes on button press in checkbox_shownodes.
function checkbox_shownodes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_shownodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_shownodes


% --- Executes on button press in checkbox_histolength.
function checkbox_histolength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histolength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histolength
if get(hObject, 'Value')
    set(handles.text_histolengthbins, 'Enable', 'on') ; 
    set(handles.popupmenu_histolengthbins, 'Enable', 'on') ; 
else 
    set(handles.text_histolengthbins, 'Enable', 'off') ; 
    set(handles.popupmenu_histolengthbins, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_logloglength.
function checkbox_logloglength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_logloglength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_logloglength
if get(hObject, 'Value')
    set(handles.checkbox_censor, 'Enable', 'on') ; 
    set(handles.checkbox_mle, 'Enable', 'on') ; 
else    
    set(handles.checkbox_censor, 'Enable', 'off') ; 
    set(handles.checkbox_mle, 'Enable', 'off') ; 
    set(handles.checkbox_censor, 'Value', 0) ; 
    set(handles.checkbox_mle, 'Value', 0) ; 
    set(handles.edit_lc, 'Enable', 'off') ; 
    set(handles.edit_uc, 'Enable', 'off') ; 
    set(handles.label_lc, 'Enable', 'off') ; 
    set(handles.label_uc, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_mle.
function checkbox_mle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mle
if get(hObject, 'Value')
    set(handles.edit_lc, 'Enable', 'on') ; 
    set(handles.edit_uc, 'Enable', 'on') ; 
    set(handles.label_lc, 'Enable', 'on') ; 
    set(handles.label_uc, 'Enable', 'on') ; 
else
    set(handles.edit_lc, 'Enable', 'off') ; 
    set(handles.edit_uc, 'Enable', 'off') ; 
    set(handles.label_lc, 'Enable', 'off') ; 
    set(handles.label_uc, 'Enable', 'off') ; 
end ; 

% --- Executes on button press in checkbox_histoangle.
function checkbox_histoangle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histoangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histoangle
if get(hObject, 'Value')
    set(handles.text_histoanglebins, 'Enable', 'on') ; 
    set(handles.popupmenu_histoanglebins, 'Enable', 'on') ; 
else 
    set(handles.text_histoanglebins, 'Enable', 'off') ; 
    set(handles.popupmenu_histoanglebins, 'Enable', 'off') ; 
end ; 


% --- Executes on button press in checkbox_rose.
function checkbox_rose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rose
if get(hObject, 'Value')
    set(handles.text_rosebins, 'Enable', 'on') ; 
    set(handles.popupmenu_rosebins, 'Enable', 'on') ; 
    set(handles.edit_degfromnorth, 'Enable', 'on') ; 
    set(handles.label_degfromnorth, 'Enable', 'on') ; 
    set(handles.checkbox_roselengthweighted, 'Enable', 'on') ; 
    set(handles.checkbox_showrosemean, 'Enable', 'on') ; 
    set(handles.checkbox_rosecolour, 'Enable', 'on') ; 
else 
    set(handles.text_rosebins, 'Enable', 'off') ; 
    set(handles.popupmenu_rosebins, 'Enable', 'off') ; 
    set(handles.edit_degfromnorth, 'Enable', 'off') ; 
    set(handles.label_degfromnorth, 'Enable', 'off') ; 
    set(handles.checkbox_roselengthweighted, 'Enable', 'off') ; 
    set(handles.checkbox_showrosemean, 'Enable', 'off') ; 
    set(handles.checkbox_rosecolour, 'Enable', 'off') ; 
    set(handles.checkbox_roselengthweighted, 'Value', false) ; 
    set(handles.checkbox_showrosemean, 'Value', false) ; 
    set(handles.checkbox_rosecolour, 'Value', false) ; 
end ; 


% --- Executes on button press in checkbox_crossplot.
function checkbox_crossplot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_crossplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_crossplot


% --- Executes on button press in checkbox_showcircles.
function checkbox_showcircles_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showcircles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showcircles


% --- Executes on button press in checkbox_censor.
function checkbox_censor_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_censor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_censor



function edit_houghpeaks_Callback(hObject, eventdata, handles)
% hObject    handle to edit_houghpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_houghpeaks as text
%        str2double(get(hObject,'String')) returns contents of edit_houghpeaks as a double


% --- Executes during object creation, after setting all properties.
function edit_houghpeaks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_houghpeaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_houghthreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edit_houghthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_houghthreshold as text
%        str2double(get(hObject,'String')) returns contents of edit_houghthreshold as a double


% --- Executes during object creation, after setting all properties.
function edit_houghthreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_houghthreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fillgap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fillgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fillgap as text
%        str2double(get(hObject,'String')) returns contents of edit_fillgap as a double


% --- Executes during object creation, after setting all properties.
function edit_fillgap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fillgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minlength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minlength as text
%        str2double(get(hObject,'String')) returns contents of edit_minlength as a double


% --- Executes during object creation, after setting all properties.
function edit_minlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_degfromnorth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_degfromnorth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_degfromnorth as text
%        str2double(get(hObject,'String')) returns contents of edit_degfromnorth as a double


% --- Executes during object creation, after setting all properties.
function edit_degfromnorth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_degfromnorth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_scaling_Callback(hObject, eventdata, handles)
% hObject    handle to edit_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_scaling as text
%        str2double(get(hObject,'String')) returns contents of edit_scaling as a double


% --- Executes during object creation, after setting all properties.
function edit_scaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numscancircles_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numscancircles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numscancircles as text
%        str2double(get(hObject,'String')) returns contents of edit_numscancircles as a double


% --- Executes during object creation, after setting all properties.
function edit_numscancircles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numscancircles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_exit.
function pushbutton_exit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all ; 
close(gcf) ; 


function edit_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filename as text
%        str2double(get(hObject,'String')) returns contents of edit_filename as a double
% handles.filename = uigetfile({'*.txt'; '*.jpeg'; '*.tiff'}, 'Select an input file for GBPaQ2D') ; 
% if length(handles.filename) > 0 
%     set(hObject, 'String', handles.filename) ; 
% else 
%     set(hObject, 'String', '(no file selected)') ; 
% end ; 

% --- Executes during object creation, after setting all properties.
function edit_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white') ;
end


% --- Executes on button press in pushbutton_browse.
function pushbutton_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes_tracemap) ; 
handles.axes_tracemap.YDir = 'normal' ; 

%   from v1.8 - all platforms get access to .svg input files 
sFileTypes = {'*.txt'; '*.svg'; '*.jpeg'; '*.jpg'; '*.tiff'; '*.tif'} ; 

[ handles.selfile, handles.selpath ] = uigetfile(sFileTypes, 'Select input file(s) for GBPaQ2D') ; 
assignin('base','user_file_name',handles.selfile)             %save that variable to the workspace to use it later in the file guiPrint.m

if ~isempty(handles.selfile)
    
    sFileName = char(handles.selfile) ; 

    if ~isempty(regexp(sFileName, 'colour', 'ONCE')) 
        iHex = regexp(sFileName, '[A-F,0-9]{6}', 'ONCE') ; 
        %   check each file has 'HHHHHH' in the name i.e. a valid hex string for colour 
        if isempty(iHex) 
            hError = errordlg('For a colour file, all filenames must contain a 6 character hexadecimal colour code', 'Input error', 'modal') ; 
            uicontrol(handles.edit_filename) ; 
            flagError = true ; 
            return ; 
        else 
            iRGB = strfind(sFileName, 'colour') ; 
            sRGB = sFileName(iRGB+6:iRGB+11) ; 
            sColour = hexRGB(sRGB) ; 
            handles.cstrColours = cellstr(sColour) ; 
        end ; 

    else
        %   MATLAB default blue 
        sColour = '[0, 0, 1]' ; 
        handles.cstrColours = cellstr(sColour) ; 
    end ; 

    handles.iCurrentColour = 1 ; 

    set(handles.edit_filename, 'String', sFileName) ; 
    set(handles.pushbutton_preview, 'Enable', 'on') ; 

    if ~isempty(strfind(upper(sFileName), '.TIF')) || ...
       ~isempty(strfind(upper(sFileName), '.TIFF')) || ...
       ~isempty(strfind(upper(sFileName), '.JPEG')) || ...
       ~isempty(strfind(upper(sFileName), '.JPG'))  
        set(handles.radiobutton_image, 'Value', 1) ;
        set(handles.radiobutton_node, 'Value', 0) ;
        set(handles.edit_houghpeaks, 'Enable', 'on') ; 
        set(handles.edit_houghthreshold, 'Enable', 'on') ; 
        set(handles.edit_fillgap, 'Enable', 'on') ; 
        set(handles.edit_minlength, 'Enable', 'on') ; 
    else 
        set(handles.radiobutton_node, 'Value', 1) ;
        set(handles.radiobutton_image, 'Value', 0) ;
        set(handles.edit_houghpeaks, 'Enable', 'off') ; 
        set(handles.edit_houghthreshold, 'Enable', 'off') ; 
        set(handles.edit_fillgap, 'Enable', 'off') ; 
        set(handles.edit_minlength, 'Enable', 'off') ; 
    end ; 
        
else
    
    set(handles.edit_filename, 'String', '(no file selected)') ; 
    set(handles.pushbutton_preview, 'Enable', 'off') ; 
    
end ; 
    
set(handles.pushbutton_flipx, 'Enable', 'off') ; 
set(handles.pushbutton_flipy, 'Enable', 'off') ; 
set(handles.pushbutton_run, 'Enable', 'off') ; 
set(handles.text_message, 'String', 'Click Preview to view the file contents.') ;                                 

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_preview.
function pushbutton_preview_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cla(handles.axes_tracemap) ; 

flagError = false ; 

%   image file validation 
if get(handles.radiobutton_image, 'Value')
    
    sValue = get(handles.edit_houghpeaks, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) < 100 || str2double(sValue) > 5000
        hError = errordlg('Hough peaks must be a whole number, 100-5000', 'Input error', 'modal') ; 
        uicontrol(handles.edit_houghpeaks) ; 
        flagError = true ; 
        return ; 
    end ; 

    sValue = get(handles.edit_houghthreshold, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) <= 0.0 || str2double(sValue) > 1.0  
        hError = errordlg('Hough threshold must be a number > 0.0, and <= 1.0', 'Input error', 'modal') ; 
        uicontrol(handles.edit_houghthreshold) ; 
        flagError = true ; 
        return ; 
    end ; 

    sValue = get(handles.edit_fillgap, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) <= 1.0    
        hError = errordlg('Merge gap must be a number of pixels, > 1 and < max. image size', 'Input error', 'modal') ; 
        uicontrol(handles.edit_fillgap) ; 
        flagError = true ; 
        return ; 
    end ; 

    sValue = get(handles.edit_minlength, 'String') ; 
    if isnan(str2double(sValue)) || str2double(sValue) <= 1   
        hError = errordlg('Discard length must be a number of pixels, > 1 and < max. image size', 'Input error', 'modal') ; 
        uicontrol(handles.edit_minlength) ; 
        flagError = true ; 
        return ; 
    end ; 
    
    Iinfo = imfinfo(get(handles.edit_filename, 'String')) ; 
    if Iinfo.BitDepth ~= 8 || ~strcmp('grayscale', Iinfo.ColorType)
        disp('Image file format error - check BitDepth (8) and ColorType(''grayscale'')') ; 
        disp(Iinfo) ; 
        hError = errordlg({'Selected image file is not 8-bit & grayscale';'Check Command Window for file info, BitDepth and ColorType'}, ... 
                                'Input error', 'modal') ; 
        uicontrol(handles.edit_filename) ; 
        flagError = true ; 
        return ; 
    end ; 

    handles.selmultifile = { } ; 

%   node file validation 
else 
    
    sFilename = [ char(handles.selpath), char(handles.selfile) ] ; 
    
    %   convert .svg file to .txt file 
    if ~isempty(strfind(sFilename, '.svg'))  
        
        nConvert = convertSVG2txt_colour(sFilename) ; 
%         disp(nConvert) ; 
        
        if nConvert > 0 
            
            %   now handle possible multiple files 
            if nConvert > 1 
                
                sFilenameHead = sFilename(1:strfind(sFilename, '.svg')-1) ; 
                infoConvertedFiles = dir([sFilenameHead, '_colour*_converted.txt']) ;
                handles.selmultifile = { infoConvertedFiles.name } ;
                sColour = '' ;
                for i = 1:nConvert 
                    
                    sThisFileName = char(handles.selmultifile(i)) ; 
                    iHex = regexp(sThisFileName, '[A-F,0-9]{6}', 'ONCE') ; 
                    
%                     disp(sThisFileName) ; 
%                     disp(iHex) ; 
                     
                    %   check each file has 'HHHHHH' in the name i.e. a valid hex string for colour 
                    if isempty(iHex) 
                        hError = errordlg('For a colour file, all filenames must contain a 6 character hexadecimal colour code', 'Input error', 'modal') ; 
                        uicontrol(handles.edit_filename) ; 
                        flagError = true ; 
                        return ; 
                    else 
                        iRGB = strfind(sThisFileName, 'colour') ; 
                        sRGB = sThisFileName(iRGB+6:iRGB+11) ; 
                        sColour = [ sColour ; hexRGB(sRGB) ]  ; 
                    end ; 

                end ; 
                
                handles.cstrColours = cellstr(sColour) ; 
                    
            else
                
                handles.selmultifile = { } ; 
                handles.selfile = strrep(handles.selfile, '.svg', '_converted.txt') ; 
                
            end ; 
            
        else
            
            hError = errordlg('Error converting .svg file to .txt; check .svg file is version 1.1', 'Input error', 'modal') ; 
            uicontrol(handles.edit_filename) ; 
            flagError = true ; 
            
        end ; 
        
    else
        
        handles.selmultifile = { } ; 

    end ; 
    
end ; 

sValue = get(handles.edit_scaling, 'String') ; 
if ~isempty(sValue)
    if isnan(str2double(sValue)) || str2double(sValue) <= 0.0 
        hError = errordlg('Scaling can be blank (i.e. scale as is) or a positive number', 'Input error', 'modal') ; 
        uicontrol(handles.edit_scaling) ; 
        flagError = true ; 
        return ; 
    end ; 
end ; 

if ~flagError 

    sFilename = [ char(handles.selpath), char(handles.selfile) ] ; 
    %   get the current RGB colour string 
    sColour = char(handles.cstrColours(handles.iCurrentColour)) ; 

    if isempty(get(handles.edit_scaling, 'String'))
        nPixels = 0 ; 
    else
        nPixels = str2double(get(handles.edit_scaling, 'String')) ; 
    end ; 
    
    if get(handles.radiobutton_image, 'Value') 
        
        nHoughPeaks = str2double(get(handles.edit_houghpeaks, 'String')) ;
        dHoughThreshold = str2double(get(handles.edit_houghthreshold, 'String')) ;
        dfillgap = str2double(get(handles.edit_fillgap, 'String')) ;
        dminlength = str2double(get(handles.edit_minlength, 'String')) ;

        set(handles.text_message, 'String', 'Running Hough transform on image...') ;  
        [ handles.traces, handles.xmin, handles.ymin, handles.xmax, handles.ymax ] = GBPaQimage(sFilename, nPixels, ...
                                                            nHoughPeaks, dHoughThreshold, dfillgap, dminlength, handles.axes_tracemap) ; 
    else
        
        set(handles.text_message, 'String', 'Reading node file...') ;                                 

        %   allow for multiple files, one for each colour 
        if ~isempty(handles.selmultifile)

            handles.traces = struct([]) ; 
            handles.nTraces = 0 ; 
            handles.nSegments = 0 ; 
            handles.nNodes = 0 ;
            handles.xmin = 3e38 ; 
            handles.ymin = 3e38 ; 
            handles.xmax = -3e38 ;
            handles.ymax = -3e38 ; 
            hold on ; 
            for i = 1:length(handles.selmultifile) 
                
                sFilename = [ char(handles.selpath), char(handles.selmultifile(i)) ] ;
                sColour = char(handles.cstrColours(i)) ; 
                
                [ traces, xmin, ymin, xmax, ymax ] = GBPaQlist(sFilename, nPixels, handles.axes_tracemap, sColour) ; 
                
                handles.traces = [ handles.traces, traces ] ; 

                if xmin < handles.xmin
                    handles.xmin = xmin ; 
                end ; 
                if ymin < handles.ymin 
                    handles.ymin = ymin ; 
                end ; 
                if xmax > handles.xmax
                    handles.xmax = xmax ; 
                end ; 
                if ymax > handles.ymax 
                    handles.ymax = ymax ; 
                end ; 
                
            end ; 
            hold off ; 
           
        else
            
            hold on ; 
            [ handles.traces, handles.xmin, handles.ymin, handles.xmax, handles.ymax ] = GBPaQlist(sFilename, nPixels, handles.axes_tracemap, sColour) ; 
            hold off ; 
            
        end ; 
        
    end ; 
    
    set(handles.text_message, 'String', 'Collecting fracture trace data...') ;     

    setappdata(handles.pushbutton_run, 'Traces', handles.traces) ;                                                 
    setappdata(handles.pushbutton_run, 'xMax', handles.xmax) ;                                                 
    setappdata(handles.pushbutton_run, 'yMax', handles.ymax) ;                                                 
    setappdata(handles.pushbutton_run, 'xMin', handles.xmin) ;                                                 
    setappdata(handles.pushbutton_run, 'yMin', handles.ymin) ;                                                 

    handles.nTraces = length(handles.traces) ; 
    handles.nSegments = sum([handles.traces(:).nSegments]) ; 
    handles.nNodes = sum([handles.traces(:).nNodes]) ; 

    set(handles.label_stats, 'String', {['Min. X coordinate: ', num2str(handles.xmin)], ...
                                        ['Min. Y coordinate: ', num2str(handles.ymin)], ...
                                        ['Max. X coordinate: ', num2str(handles.xmax)], ...
                                        ['Max. Y coordinate: ', num2str(handles.ymax)], ...
                                        ['Number of traces: ', num2str(handles.nTraces)], ...
                                        ['Number of segments: ', num2str(handles.nSegments)], ...
                                        ['Number of nodes: ', num2str(handles.nNodes)]}) ;

    set(handles.text_message, 'String', 'Ready. Click Run to generate maps and graphs.') ;     
  
    set(handles.pushbutton_run, 'Enable', 'on') ; 
    set(handles.pbMapsTab, 'Enable', 'on') ; 
    set(handles.pbLengthsTab, 'Enable', 'on') ; 
    set(handles.pbAnglesTab, 'Enable', 'on') ; 
    set(handles.pushbutton_flipx, 'Enable', 'on') ; 
    set(handles.pushbutton_flipy, 'Enable', 'on') ; 
%    uicontrol(handles.pushbutton_run) ; 
    uicontrol(handles.pbMapsTab) ; 
    set(handles.uibgTabMaps, 'Visible', 'on') ; 
    
    % Update handles structure
    guidata(hObject, handles);
    
end ; 

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

flagError = false ; 

sValue = get(handles.edit_numscancircles, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) < 2
    hError = errordlg('Number of scan circles must be a positive whole number (e.g. 10-20)', 'Input error', 'modal') ; 
    uicontrol(handles.edit_numscancircles) ; 
    flagError = true ; 
    return ; 
end ; 

sValue = get(handles.edit_degfromnorth, 'String') ; 
if isnan(str2double(sValue)) || str2double(sValue) > 360.0 || str2double(sValue) < -360.0
    hError = errordlg('Degrees from North must be a number <= 360 (use negative for CCW)', 'Input error', 'modal') ; 
    uicontrol(handles.edit_degfromnorth) ; 
    flagError = true ; 
    return ; 
end ; 

if ~flagError 
    
    traces = getappdata(hObject, 'Traces') ; 
    xmin = getappdata(hObject, 'xMin') ; 
    ymin = getappdata(hObject, 'yMin') ; 
    xmax = getappdata(hObject, 'xMax') ; 
    ymax = getappdata(hObject, 'yMax') ; 

    %   get the current RGB colour string 
    if ~isempty(handles.selmultifile)
        sColour = '[ 0 0 1 ]' ; 
    else
        sColour = char(handles.cstrColours(handles.iCurrentColour)) ; 
    end ; 
    
    if length(traces) > 0 

        %   get values from text boxes 
        if isempty(get(handles.edit_scaling, 'String'))
            nPixels = 0 ; 
        else
            nPixels = str2double(get(handles.edit_scaling, 'String')) ; 
        end ; 
        
        if strcmp(handles.axes_tracemap.YDir, 'reverse')
            flag_reverseY = true ; 
        else
            flag_reverseY = false ; 
        end ;
        
        if strcmp(handles.axes_tracemap.XDir, 'reverse')
            flag_reverseX = true ; 
        else
            flag_reverseX = false ; 
        end ;
        
        nNorth = str2double(get(handles.edit_degfromnorth, 'String')) ;
        nMLE_lc = str2double(get(handles.edit_lc, 'String')) ;  
        nMLE_uc = str2double(get(handles.edit_uc, 'String')) ;  
        
        %   get values from drop down lists 
        list = get(handles.popupmenu_histolengthbins, 'String') ; 
        nHistoLengthBins = str2double(list{get(handles.popupmenu_histolengthbins, 'Value')}) ;
        
        list = get(handles.popupmenu_histoanglebins, 'String') ; 
        nHistoAngleBins = str2double(list{get(handles.popupmenu_histoanglebins, 'Value')}) ;
        
        list = get(handles.popupmenu_rosebins, 'String') ; 
        nRoseBins = str2double(list{get(handles.popupmenu_rosebins, 'Value')}) ; 
    
        %   call the various maps & graphs 
        flag_tracemap = get(handles.checkbox_tracemap, 'Value') ; 
        flag_shownodes = get(handles.checkbox_shownodes, 'Value') ; 
        flag_tracesbylength = get(handles.checkbox_tracemaplength, 'Value') ; 
        flag_segmentsbylength = get(handles.checkbox_segmentmaplength, 'Value') ; 
        flag_segmentsbystrike = get(handles.checkbox_segmentmapstrike, 'Value') ; 
        
        flag_tracemap_any = false ; 
        if sum([flag_tracemap, flag_tracesbylength, ...
                flag_segmentsbylength, flag_segmentsbystrike]) > 0 
            flag_tracemap_any = true ; 
        end ; 
        
        flag_intensitymap = get(handles.checkbox_intensitymap, 'Value') ; 
        flag_densitymap = get(handles.checkbox_densitymap, 'Value') ; 
        flag_showcircles = get(handles.checkbox_showcircles, 'Value') ; 
        nCircles = round(str2double(get(handles.edit_numscancircles, 'String'))) ; 

        flag_histolength = get(handles.checkbox_histolength, 'Value') ; 
        flag_logloglength = get(handles.checkbox_logloglength, 'Value') ;
        flag_crossplot = get(handles.checkbox_crossplot, 'Value') ;
        flag_censor = get(handles.checkbox_censor, 'Value') ; 
        flag_mle = get(handles.checkbox_mle, 'Value') ;
        flag_seglenvario = get(handles.checkbox_seglenvariogram, 'Value') ; 

        flag_histoangle = get(handles.checkbox_histoangle, 'Value') ;
        flag_roseangle = get(handles.checkbox_rose, 'Value') ;
        flag_rosemean = get(handles.checkbox_showrosemean, 'Value') ; 
        flag_roselengthweighted = get(handles.checkbox_roselengthweighted, 'Value') ; 
        flag_rosecolour = get(handles.checkbox_rosecolour, 'Value') ; 
        
        flag_cracktensor = get(handles.checkbox_cracktensor, 'Value') ;

        if length(handles.cstrColours) > 1 
            flag_multicolour = 1 ; 
        else 
            flag_multicolour = 0 ; 
        end ;

        if flag_tracemap_any 
            GBPaQtracemap(traces, nPixels, nNorth, ...
                    xmin, ymin, xmax, ymax, ...
                    flag_shownodes, flag_reverseY, flag_reverseX, sColour, ...
                    flag_multicolour, flag_tracemap, ...
                    flag_tracesbylength, flag_segmentsbylength, flag_segmentsbystrike) ; 
        end ; 

        flag_length = sum([flag_histolength, flag_logloglength, flag_blocksize, flag_seglenvario]) ; 
        if flag_length
            GBPaQlength(traces, nPixels, nNorth, xmin, ymin, xmax, ymax, ...
                               nHistoLengthBins, flag_histolength, flag_logloglength, ...
                               flag_crossplot, flag_mle, flag_censor, ...
                               flag_blocksize, flag_reverseY, flag_reverseX, ...
                               sColour, nMLE_lc, nMLE_uc, flag_seglenvario) ; 
        end ; 

        flag_angle = sum([flag_histoangle, flag_roseangle, flag_cracktensor]) ; 
        if flag_angle
            GBPaQangle(traces, nNorth, xmin, ymin, xmax, ymax, nHistoAngleBins, nRoseBins, ...
                                    flag_histoangle, flag_roseangle, ...
                                    flag_reverseY, flag_reverseX, ...
                                    flag_cracktensor, ...
                                    flag_roselengthweighted, flag_rosemean, ...
                                    flag_rosecolour, sColour) ; 
        end ; 

        flag_pattern = sum([flag_intensitymap, flag_densitymap]) ; 
        if flag_pattern
            GBPaQpattern(traces, nPixels, xmin, ymin, xmax, ymax, ...
                                    flag_intensitymap, flag_densitymap, flag_showcircles, nCircles, flag_reverseY, flag_reverseX, sColour) ; 
        end ; 

        set(handles.text_message, 'String', 'All done. Check current folder for plot figure .tif files.') ;   
        uicontrol(handles.pushbutton_run) ; 

    else 
        
        errordlg('No traces found; check the input file is OK.') ; 
        uicontrol(handles.pushbutton_browse) ; 

    end ; 
   
end ; 

% --- Executes on button press in radiobutton_image.
function radiobutton_image_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_image
set(handles.edit_houghpeaks, 'Enable', 'on') ; 
set(handles.edit_houghthreshold, 'Enable', 'on') ; 
set(handles.edit_fillgap, 'Enable', 'on') ; 
set(handles.edit_minlength, 'Enable', 'on') ; 


% --- Executes on button press in radiobutton_node.
function radiobutton_node_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_node (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_node
set(handles.edit_houghpeaks, 'Enable', 'off') ; 
set(handles.edit_houghthreshold, 'Enable', 'off') ; 
set(handles.edit_fillgap, 'Enable', 'off') ; 
set(handles.edit_minlength, 'Enable', 'off') ; 


% --- Executes during object creation, after setting all properties.
function figure_main_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in popupmenu_histolengthbins.
function popupmenu_histolengthbins_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_histolengthbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_histolengthbins contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_histolengthbins


% --- Executes during object creation, after setting all properties.
function popupmenu_histolengthbins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_histolengthbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_histoanglebins.
function popupmenu_histoanglebins_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_histoanglebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_histoanglebins contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_histoanglebins


% --- Executes during object creation, after setting all properties.
function popupmenu_histoanglebins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_histoanglebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_rosebins.
function popupmenu_rosebins_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_rosebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_rosebins contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_rosebins


% --- Executes during object creation, after setting all properties.
function popupmenu_rosebins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_rosebins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_flipy.
function pushbutton_flipy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_flipy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.axes_tracemap.YDir, 'reverse')
    handles.axes_tracemap.YDir = 'normal' ;
else
    handles.axes_tracemap.YDir = 'reverse' ;
end ;     
set(handles.text_message, 'String', 'Y-axis flipped.') ;   

% --- Executes on button press in pushbutton_flipx.
function pushbutton_flipx_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_flipx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.axes_tracemap.XDir, 'reverse')
    handles.axes_tracemap.XDir = 'normal' ;
else
    handles.axes_tracemap.XDir = 'reverse' ;
end ;     
set(handles.text_message, 'String', 'X-axis flipped.') ;   


% --- Executes on button press in checkbox_cracktensor.
function checkbox_cracktensor_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_cracktensor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_cracktensor


% --- Executes on button press in checkbox_roselengthweighted.
function checkbox_roselengthweighted_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_roselengthweighted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_roselengthweighted


% --- Executes on button press in checkbox_showrosemean.
function checkbox_showrosemean_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showrosemean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showrosemean


% --- Executes on button press in checkbox_blocksize.
function checkbox_blocksize_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_blocksize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_blocksize



function edit_lc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lc as text
%        str2double(get(hObject,'String')) returns contents of edit_lc as a double


% --- Executes during object creation, after setting all properties.
function edit_lc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_uc_Callback(hObject, eventdata, handles)
% hObject    handle to edit_uc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_uc as text
%        str2double(get(hObject,'String')) returns contents of edit_uc as a double


% --- Executes during object creation, after setting all properties.
function edit_uc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_uc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_rosecolour.
function checkbox_rosecolour_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rosecolour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rosecolour


% --- Executes on button press in checkbox_seglenvariogram.
function checkbox_seglenvariogram_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_seglenvariogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_seglenvariogram


% --- Executes on button press in pbMapsTab.
function pbMapsTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbMapsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'on') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in pbLengthsTab.
function pbLengthsTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbLengthsTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'on') ; 
set(handles.uibgTabAngles, 'Visible', 'off') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in pbAnglesTab.
function pbAnglesTab_Callback(hObject, eventdata, handles)
% hObject    handle to pbAnglesTab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uibgTabMaps, 'Visible', 'off') ; 
set(handles.uibgTabLengths, 'Visible', 'off') ; 
set(handles.uibgTabAngles, 'Visible', 'on') ; 
set(handles.uibgTabFluid, 'Visible', 'off') ; 
set(handles.uibgTabWavelets, 'Visible', 'off') ; 


% --- Executes on button press in checkbox_tracemaplength.
function checkbox_tracemaplength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_tracemaplength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_tracemaplength


% --- Executes on button press in checkbox_segmentmaplength.
function checkbox_segmentmaplength_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentmaplength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentmaplength


% --- Executes on button press in checkbox_segmentmapstrike.
function checkbox_segmentmapstrike_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_segmentmapstrike (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_segmentmapstrike


function edit_DegString_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DegString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DegString as text
%        str2double(get(hObject,'String')) returns contents of edit_DegString as a double


% --- Executes during object creation, after setting all properties.
function edit_DegString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DegString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
