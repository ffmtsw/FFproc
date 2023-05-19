% Creates a *.PFC file for use by Phoenix software syscal_v6.exe
% Frequency spacing should be fine enough to cover notches
% Order of frequency is arbitrary but must be monotonous.

% version 1.1 / 20feb2019 / aj
% version 1.2 / 12022020 / cc      Improved to create PFC directory, to be
%                                  used in CreateCTS_.

function create_pfc

    % Create PFC directory
    [user,~] = getuser;
    cd(fullfile(user,'Documents','FFMT','FFproc','Calibration Files','Phoenix Geophysics'))
    path = pwd;
    
    path = uigetdir(path,'Directory in which be created the new PFC folder');
    cd(path)
    mkdir PFC
    cd PFC

    % Output filename
    fileout = {'PFC_2.pfc';'PFC_3.pfc';'PFC_4.pfc';'PFC_5.pfc'};
    % nb of frequencies between 1e1 Hz and 1e4 Hz (maximum frequency, may be changed)
    Nf1 = 1001;
    % nb of frequencies between 1e-4 Hz and 1e1 Hz
    Nf2 = 101;   
    % nb equals band nb of *.TSn file, will be added to fileout
    nb = 2:5;                     

    freq1 = logspace(4,1,Nf1);
    freq2 = logspace(1,-4,Nf2);
    freq = [freq1 freq2];
    freq = unique(freq);
    Nf = length(freq);

    for i = 1:length(nb)
        fileid = fopen(fileout{i},'w');  % output file name
        text = 'FldType, 1, \n';
        fprintf(fileid,text);
        for j = 1:Nf
            text = ['Frequency, ' num2str(nb(i)) ', ' num2str(freq(j)) ' \n'];
            fprintf(fileid,text);
        end

        disp(['** ',num2str(Nf),' frequencies saved in ',fileout{i},' **'])
        fclose(fileid);
    end

end