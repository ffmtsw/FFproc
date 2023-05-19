%% ===================================================================== %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%     MATLAB-SCRIPT FOR EXTRACTING INFORMATION FROM *.TBL(txt)-files      %
%                     (for MTU-5A and V8 systems)                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

% Selecting TBL-files (You can choose multiple files at the same time)
[file,path] = uigetfile({'*.txt','Text files (*.txt)';...
                         '*.*','All files (*,*)'},...
                         'Select input file',...
                         'multiselect','on');
                     
% Creating path...   
file = cellstr(file)'; 
FILE = file;
file_path = strcat(path,file);
FILE = strrep(FILE,'.txt','');

N_file = length(file);                          % Number of files
Properties_table = cell(N_file,40);             % Size of properties table
infmt = 'yyyy/MM/dd HH:mm:ss';                  % Format for date and time
infmt2 = 'yyyy/MM/dd HH:mm:ss'; 

for i = 1:N_file
    
    clear instring
    fileID = fopen(file_path{i});
    split_path = regexp(file_path{1},'\','split');
    instring = textscan(fileID,'%s');               % Reading each element of the file
    
    instring = strrep(instring{1,1},',','');        % Deleting ','
    
    %% FINDING PARAMETERS
    
    FOLDER = split_path{1,end-1};                                   % Folder containind TBL File
    HW = instring{find(strcmp(instring,'HW'))+7,1};                 % Box model
    SNUM = instring{find(strcmp(instring,'SNUM'))+7,1};             % Box Serial Number
    SITE = instring{find(strcmp(instring,'SITE'))+7,1};             % Site
    STIM_date = instring{find(strcmp(instring,'STIM'))+7,1};        % Programed Start-Up date
    STIM_hour = instring{find(strcmp(instring,'STIM'))+8,1};        % Programed Start-Up time
    ETIM_date = instring{find(strcmp(instring,'ETIM'))+7,1};        % Programed End-Time date
    ETIM_hour = instring{find(strcmp(instring,'ETIM'))+8,1};        % Programed End-Time time
    HTIM_date = instring{find(strcmp(instring,'HTIM'))+7,1};        % Programmed HF Start-Up Logging date
    HTIM_hour = instring{find(strcmp(instring,'HTIM'))+8,1};        % Programmed HF Start-Up Logging time
    ETMH_date = instring{find(strcmp(instring,'ETMH'))+7,1};        % Programmed HF End-Time Logging date
    ETMH_hour = instring{find(strcmp(instring,'ETMH'))+8,1};        % Programmed HF End-Time Logging time
    FTIM_date = instring{find(strcmp(instring,'FTIM'))+7,1};        % Start-Up date
    FTIM_hour = instring{find(strcmp(instring,'FTIM'))+8,1};        % Start-Up time
    LTIM_date = instring{find(strcmp(instring,'LTIM'))+7,1};        % End-Time date
    LTIM_hour = instring{find(strcmp(instring,'LTIM'))+8,1};        % End-Time time
    EXLN = instring{find(strcmp(instring,'EXLN'))+7,1};             % X-Dipole longitude (m)
    EYLN = instring{find(strcmp(instring,'EYLN'))+7,1};             % Y-Dipole longitude (m)
    HXSN = instring{find(strcmp(instring,'HXSN'))+7,1};             % Hx Sensor
    HYSN = instring{find(strcmp(instring,'HYSN'))+7,1};             % Hy Sensor
    HZSN = instring{find(strcmp(instring,'HZSN'))+7,1};             % Hz Sensor
    EGN = instring{find(strcmp(instring,'EGN'))+7,1};               % E-Chan Gain
    HGN = instring{find(strcmp(instring,'HGN'))+7,1};               % H-Chan Gain
    LFRQ = instring{find(strcmp(instring,'LFRQ'))+7,1};             % Power-Line Filter
    TOTL = instring{find(strcmp(instring,'TOTL'))+7,1};             % Total Records
    LATG = instring{find(strcmp(instring,'LATG'))+7,1};             % Latitude
    LNGG = instring{find(strcmp(instring,'LNGG'))+7,1};             % Longitude
    
    %% CALCULATING MEASURED TIMES
    STIM = datetime([STIM_date,' ',STIM_hour],'inputformat',infmt);
    ETIM = datetime([ETIM_date,' ',ETIM_hour],'inputformat',infmt);
    % Elapsed time for Programmed measured time
    PT_dur = ETIM - STIM;
        
    HTIM = datetime([HTIM_date,' ',HTIM_hour],'inputformat',infmt);
    ETMH = datetime([ETMH_date,' ',ETMH_hour],'inputformat',infmt);
    % Elapsed time for High Frequencies recording
    HF_dur = ETMH - HTIM;
    
    FTIM = datetime([FTIM_date,' ',FTIM_hour],'inputformat',infmt);
    LTIM = datetime([LTIM_date,' ',LTIM_hour],'inputformat',infmt);
    % Elapsed time for Measured time
    MEAS_dur = LTIM - FTIM;
    
    %% TIME IN MEXICO
    STIM_MX = STIM - hours(6);
    ETIM_MX = ETIM - hours(6);
    HTIM_MX = HTIM - hours(6);
    ETMH_MX = ETMH - hours(6);
    FTIM_MX = FTIM - hours(6);
    LTIM_MX = LTIM - hours(6);
    
    %% TRANSFORMING GEOGRAPHIC COORDINATES TO UTM COORDINATES
    
    LATG = insertAfter(LATG,2,'°');
    LATG = insertAfter(LATG,10,'''');
    LNGG = insertAfter(LNGG,3,'°');
    LNGG = insertAfter(LNGG,11,'''');
    ELEV = instring{find(strcmp(instring,'ELEV'))+7,1};
    lon = strrep(LNGG,'W','');
    lon = strrep(lon,'''','');
    lon_g = extractBefore(lon,'°');
    lon_m = extractAfter(lon,'°');
    lon = -1*(str2num(lon_g)+ str2num(lon_m)/60);
    lat = strrep(LATG,'N','');
    lat = strrep(lat,'''','');
    lat_g = extractBefore(lat,'°');
    lat_m = extractAfter(lat,'°');
    lat = str2num(lat_g)+ str2num(lat_m)/60;
    [X,Y] = deg2utm(lat,lon);                                       % UTM Coordinates
    
    % SAVING INFORMATION
    
    Properties_table(i,:) = strrep({FILE{i,1},FOLDER,...
                               HW,SNUM,SITE,...
                               STIM_date,STIM_hour,ETIM_date,ETIM_hour,datestr(STIM_MX),datestr(ETIM_MX),datestr(PT_dur,'HH:MM:SS'),...
                               HTIM_date,HTIM_hour,ETMH_date,ETMH_hour,datestr(HTIM_MX),datestr(ETMH_MX),datestr(HF_dur,'HH:MM:SS'),...
                               FTIM_date,FTIM_hour,LTIM_date,LTIM_hour,datestr(FTIM_MX),datestr(LTIM_MX),datestr(MEAS_dur,'HH:MM:SS'),...
                               EXLN,EYLN,HXSN,HYSN,HZSN,...
                               EGN,HGN,LFRQ,TOTL,...
                               LATG,LNGG,num2str(X),num2str(Y),ELEV},',','');
                           
     fclose(fileID);

end

% Table Header
Properties_header = {'FILE','FOLDER',...
                     'HW','SNUM','SITE',...
                     'STIM_date','STIM_hour','ETIM_date','ETIM_hour','STIM_MX','ETIM_MX','PT_dur',...
                     'HTIM_date','HTIM_hour','ETMH_date','ETMH_hour','HTIM_MX','ETMH_MX','HF_dur'...
                     'FTIM_date','FTIM_hour','LTIM_date','LTIM_hour','FTIM_MX','LTIM_MX','MEAS_dur'...
                     'EXLN','EYLN','HXSN','HYSN','HZSN',...
                     'EGN','HGN','LFRQ','TOTL',...
                     'LATG','LNGG','X','Y','ELEV'};

%% TABLE (Results)
Properties = [Properties_header;Properties_table];

%% ---------------------------- END PROGRAM ---------------------------- %%