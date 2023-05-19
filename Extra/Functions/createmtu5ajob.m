clear, close, clc

path = uigetdir(pwd,'Select input folder:');
cd(path)
info = dir(fullfile(path,'*.TS*'));

for i = 1:numel(info)

    jobnew = job_create;
    ts = ts_input(path,info(i).name,false,false);
    a = find(strcmp(ts.device,cat(1,param.device)));

    jobnew.name = {ts.name};
    jobnew.ffts = ts;
    jobnew.sr = ts.sr;   
    jobnew.maxf = param(a).(strrep(lower(ts.ext),'.','')).maxf;
    jobnew.minf = param(a).(strrep(lower(ts.ext),'.','')).minf;
    jobnew.ftarg = 6;
    jobnew.osc = param(a).(strrep(lower(ts.ext),'.','')).osc;
    jobnew.exclude = true;
    jobnew.singf = 50;
    jobnew.phsing = 4;
    jobnew.mult = 50;
    jobnew.phmult = 4;
    jobnew.bivar = false;
    jobnew.ev = true;
    jobnew.model = 'Advanced';
    jobnew.reg = 'LSQ';
    jobnew.outlier = true;
    jobnew.alpha = 0.99;
    jobnew.zrank = 'EVPCC';
    jobnew.trank = 'EVPCC';
    jobnew.ranking = 25;
    jobnew.eind = 'EV4';
    jobnew.evpow = 1;
    jobnew.rpow = 1;
    jobnew.polyfit = true;
    jobnew.opt = true;
    jobnew.raw = false;
    jobnew.fc = false;
    jobnew.rnk = true;
    jobnew.ncores = feature('numcores') - 1;

    job(i) = jobnew;

end