function catjobs(varargin)

    if isempty(varargin)
        % Selecting job files to be concatenated
        [file,path] = uigetfile('*.mat','Select jobs to be concatenated:','MultiSelect','on');
    else
        path = varargin{1};
        file = varargin{2};
    end

    % Loading and concatenating each file
    jobcat = [];
    for j = 1:numel(file)
        load(fullfile(path,file{j}))
        jobcat = [jobcat,job];
    end

    job = jobcat;

    % Saving file
    if isempty(varargin)
        [file,path] = uiputfile('_JOB.mat','Save Job file:');
        cd(path)
        save(fullfile(path,file),'job','-mat')
    else
        file = 'AUTOJOB.mat';
        save(fullfile(path,file),'job','-mat')
    end

end