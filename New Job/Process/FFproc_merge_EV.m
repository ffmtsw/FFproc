function EV = FFproc_merge_EV(TF,EigVal)

    EV = struct;
    for i = 1:numel(TF)
        for j = 1:numel(TF{1,i})
            EV(i).EVal(j,:) = TF{i}(j).EVal;
        end
        EV(i).EI0 = EigVal{i}.EI0;
        EV(i).EI1 = EigVal{i}.EI1;
        EV(i).EI2 = EigVal{i}.EI2;
        EV(i).EI3 = EigVal{i}.EI3;
        EV(i).EI4 = EigVal{i}.EI4;
        EV(i).EI = EigVal{i}.EI;
        EV(i).EI_sort = EigVal{i}.EI_sort;
        EV(i).EI_ind = EigVal{i}.EI_ind;
        EV(i).evpow = EigVal{i}.evpow;
        EV(i).rpow = EigVal{i}.rpow;
        EV(i).Zxx_pcorr = EigVal{i}.Zxx_pcorr;
        EV(i).Zxy_pcorr = EigVal{i}.Zxy_pcorr;
        EV(i).Zyx_pcorr = EigVal{i}.Zyx_pcorr;
        EV(i).Zyy_pcorr = EigVal{i}.Zyy_pcorr;
        EV(i).txz_pcorr = EigVal{i}.txz_pcorr;
        EV(i).tyz_pcorr = EigVal{i}.tyz_pcorr;
        EV(i).EIPC_Zxx =  EigVal{i}.EIPC_Zxx;
        EV(i).EIPC_Zxx_ind = EigVal{i}.EIPC_Zxx_ind;
        EV(i).EIPC_Zxy =  EigVal{i}.EIPC_Zxy;
        EV(i).EIPC_Zxy_ind =  EigVal{i}.EIPC_Zxy_ind;
        EV(i).EIPC_Zyx = EigVal{i}.EIPC_Zyx;
        EV(i).EIPC_Zyx_ind =  EigVal{i}.EIPC_Zyx_ind;
        EV(i).EIPC_Zyy =  EigVal{i}.EIPC_Zyy;
        EV(i).EIPC_Zyy_ind = EigVal{i}.EIPC_Zyy_ind;
        EV(i).EIPC_txz =  EigVal{i}.EIPC_txz;
        EV(i).EIPC_txz_ind = EigVal{i}.EIPC_txz_ind;
        EV(i).EIPC_tyz = EigVal{i}.EIPC_tyz;
        EV(i).EIPC_tyz_ind = EigVal{i}.EIPC_tyz_ind;
        EV(i).segments = EigVal{i}.segments;
    end    
end