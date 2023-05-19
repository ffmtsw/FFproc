function job = jobconfig(job,ffts,device,type)

    switch device
        case 'Metronix ADU-07e'
            switch ffts.sr
                case 40960
                    osc = 2000;
                    maxf = 15000;
                    minf = 6000;
                case 32768
                    osc = 2000;
                    maxf = 12000;
                    minf = 400;
                case 8192
                    osc = 2000;
                    maxf = 3500;
                    minf = 600;
                case 4096
                    osc = 2000;
                    maxf = 1750;
                    minf = 200;
                case 2048
                    osc = 2000;
                    maxf = 800;
                    minf = 60;
                case 1024
                    osc = 2000;
                    maxf = 400;
                    minf = 6;
                case 512
                    osc = 500;
                    maxf = 200;
                    minf = 0.1;
                case 256
                    osc = 500;
                    maxf = 100;
                    minf = 0.025;
                case 128
                    osc = 500;
                    maxf = 60;
                    minf = 0.001;
                case 64
                    osc = 200;
                    maxf = 20;
                    minf = 0.001;
                case 2
                    osc = 50;
                    maxf = 0.1;
                    minf = 0.001;
            end
        case 'Phoenix MTU5-A'
            switch type
                case 'MT'
                    switch ffts.sr
                        case 2400
                            osc = 1000;
                            maxf = 600;
                            minf = 40;
                        case 150
                            osc = 500;
                            maxf = 60;
                            minf = 4;
                        case 15
                            osc = 100;
                            maxf = 6;
                            minf = 0.0005;
                    end
                case 'AMT'
                    switch ffts.sr
                        case 24000
                            osc = 2000;
                            maxf = 10000;
                            minf = 400;        
                        case 2400
                            osc = 1000;
                            maxf = 600;
                            minf = 40;
                        case 150
                            osc = 500;
                            maxf = 60;
                            minf = 0.4;
                    end
                case 'BAMT'
                    switch ffts.sr
                        case 24000
                            osc = 2000;
                            maxf = 10000;
                            minf = 400;        
                        case 2400
                            osc = 1000;
                            maxf = 600;
                            minf = 40;
                        case 150
                            osc = 500;
                            maxf = 60;
                            minf = 0.0005;        
                    end
                case 'BBMT'
                    switch ffts.sr
                        case 2400
                            osc = 1000;
                            maxf = 600;
                            minf = 40;
                        case 150
                            osc = 500;
                            maxf = 60;
                            minf = 4;
                        case 15
                            osc = 100;
                            maxf = 6;
                            minf = 0.0005;     
                    end
            end    
        case 'Phoenix MTU5-C'
            switch ffts.sr
                case 24000
                    osc = 2000;
                    maxf = 10000;
                    minf = 40;   
                case 150
                    osc = 250;
                    maxf = 60;
                    minf = 0.0005;   
            end
    end

    job.name = {ffts.name};
    job.ffts = ffts;
    job.sr = ffts.sr;
    job.maxf = maxf;
    job.minf = minf;
    job.ftarg = 6;
    job.osc = osc;

    job.bivar = false;
    job.ev = true;
    job.model = 'Advanced';
    job.reg = 'Robust';
    job.outlier = false;
    job.remmethod = 'N/A';
    job.zrank = 'EV';
    job.trank = 'EV';
    job.ranking = 25;
    job.eind = 'EI3';
    job.evpow = 1;
    job.rpow = 0;
    job.polyfit = true;
    job.opt = true;
    job.raw = false;
    job.fc = false;
    job.rnk = true;
    job.ncores = feature('numcores') - 1;
    
end