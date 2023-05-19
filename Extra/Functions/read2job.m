% read2job reads all the sites/stations recorded duuring a survey/campaign
% and automatically creates all the jobs to run single-site processing
% using FFMT - FFproc. Results are saved in the AUTOJOBS folder.

% version 1.0 / 20jan2023 / cc script to create jobs for:
%                              - Metronix ADU07-e
%                              - Phoenix MTU5-A
%                              - Phoenix MTU5-C/Phoenix MTU8-A

function read2job(varargin)

    % Input device selection
    if isempty(varargin)
        answer = menu('Select input device','Metronix ADU-07e', ...
                                            'Phoenix MTU5-A', ...
                                            'Phoenix MTU5-C');
        devopts = {'Metronix ADU-07e','Phoenix MTU5-A','Phoenix MTU5-C'};
        device = devopts{answer};
    else
        device = varargin{1};
    end

    % Select root directory containing all the measurements
    root = uigetdir(pwd,['Select folder containing ',device,' measurements']);
    cd(root);
    info = dir(root);
    % Remove folders with /. and /..
    info(strcmp(extractfield(info,'name'),'.') | strcmp(extractfield(info,'name'),'..')) = [];

    % Frequencies to be excluded. Should be entered in the same line with
    % space separated
    prompt = {'Frequencies to be excluded (Hz):'};
    dlgtitle = ['Load ',device];
    dims = [1 100];
    definput = {'60'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    if iscell(answer)
        singf = str2num(answer{1});
    else
        singf = str2num(answer);
    end
    phsingf = 3*ones(1,numel(singf));
    mult = singf;
    phmult = phsingf;
    
    if isempty(singf) | singf == 0
        exclude = false;
    else
        exclude = true;
    end   

    
    % Waitbar 
    wb = waitbar(0,['Create job for site: ',info(1).name],'Name','FFproc - Creating jobs');

    % Failed sites - Errors occurred when creating job
    fail = false(numel(info),1);
    for i = 1:numel(info)
        try
            % Update waitbar
            waitbar(i/numel(info),wb,['Create job for site: ',info(i).name],'Interpreter','None');

            % Device selection
            switch device
                case 'Metronix ADU07-e'
                    job = read2job_ATS(fullfile(root,info(i).name));
                case 'Phoenix MTU5-A'
                    job = read2job_MTU5A(fullfile(root,info(i).name));
                case 'Phoenix MTU5-C'
                    job = read2job_MTU5C(fullfile(root,info(i).name));
            end

            % Adding excluded frequencies
            for j = 1:numel(job)
                job(j).exclude = exclude;
                job(j).singf = singf;
                job(j).phsing = phsingf;
                job(j).mult = mult;
                job(j).phmult = phmult;
            end

            % Saving file
            jobdir = fullfile(root,'AUTOJOBS');
            if ~isfolder(jobdir)
                mkdir(jobdir)
            end
            save(fullfile(jobdir,[info(i).name,'.mat']),'job','-mat')

        catch            
            fail(i) = true;
        end

    end

    % Close and delete waitbar
    close(wb)
    delete(wb)
    
    % Create job list with all the jobs (single site processing)
    listfiles = extractfield(info,'name');
    listfiles(fail) = [];
    catjobs(jobdir,listfiles)

end