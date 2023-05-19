function cal_inst = calibration_Artificial(job,tseg)
          
    sr = job.sr;            
    tint = tseg;
    nxt = round(tint*sr);
    for i = 1:5
        cal{i} = 1;
    end
    cal_inst = cal;

end