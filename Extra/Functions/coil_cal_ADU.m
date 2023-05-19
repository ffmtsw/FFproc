% function coil_cal_ADU reads serial number/name and sampling rate from 
% egts time series, calculation calibration data for each sensor. The
% output is a cell array with 3 cells containing information from each
% sensor in a matrix (usually 3 columns).
% version 1.0 / ????????? / AG
% version 2.0 / 09apr2020 / cc   output was squezzed into a cell array
%                                instead of having 3 different variables as
%                                an input and output
% version 2.1 / 27jan2021 / cc   user ID is added to retrieve calibration
%                                files from Documents Directory

function coil_cal = coil_cal_ADU(coil,sr)

    % Identifying Operative System
    [user,~] = getuser;
    
    % Converting string typed in calibration section to numbers
    coil = str2double(cellstr(coil));
    
    % Compare coil SN
    for i = 1:numel(coil)
        try
            if coil(1)<102 && coil(1)>92 || ...
               coil(2)<102 && coil(2)>92 || ...
               coil(3)<102 && coil(3)>92
                if coil(i) < 100
                    name = ['HC.','0',num2str(coil(i))];
                else
                    name = ['HC.',num2str(coil(i))];
                end
                fid = fopen(fullfile(user,'Documents','FFMT','FFproc','Calibration Files','Metronix',name));
                fseek(fid,432,'bof');
                cal = fscanf(fid,'%f %f %f',[3 inf]);
                cal = cal';
                [~,fi] = sort(cal(:,1));
                coil_cal{i} = cal(fi,:);
                coil_cal{i}(:,2) = coil_cal{i}(:,2)*800;
                fclose(fid);
                clear fid cal name
            else
                if sr < 513
                    chopper = 'on';
                else
                    chopper = 'off';
                end
                name = [chopper,'.',num2str(coil(i))];
                fid = fopen(fullfile(user,'Documents','FFMT','FFproc','Calibration Files','Metronix',name));
                cal = fscanf(fid,'%f %f %f',[3 inf]);
                cal = cal';
                [~,fi] = sort(cal(:,1));
                coil_cal{i} = cal(fi,:);
                fclose(fid);
                clear fid name

                % Scale factor for amplitudes
                if coil(1) == 130 || coil(1) == 115 ||...
                   coil(2) == 130 || coil(2) == 115 ||...
                   coil(3) == 130 || coil(3) == 115
                    coil_cal{i}(:,2) = coil_cal{i}(:,2).*coil_cal{i}(:,1)*640;
                else
                    coil_cal{i}(:,2) = coil_cal{i}(:,2).*coil_cal{i}(:,1)*1000;
                end
            end
        catch
            coil_cal{i}(:,1) = coil_cal{i-1}(:,1);
            coil_cal{i}(:,2) = ones(size(coil_cal{i-1}(:,1),1),1);
            coil_cal{i}(:,3) = ones(size(coil_cal{i-1}(:,1),1),1);
        end
    end
end