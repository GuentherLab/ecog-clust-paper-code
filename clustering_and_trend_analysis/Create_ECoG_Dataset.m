function ECoGALL = Create_ECoG_Dataset(inpFiles, SETTING, ECoG_Fs, Audio_Fs)

    ECoGALL = cell(1,length(inpFiles));
    for ds = 1:length(inpFiles)
        
        disp (inpFiles{ds})
        load(inpFiles{ds});
        ECoG = [];
        for ch = 1:256
            if SETTING ==1
                tempVar = eval(sprintf('LFPx_RZ2_chn%03d',ch));
                a = double(tempVar.dat');
                ECoG.Fs(ch) = mean(tempVar.fs);

            elseif SETTING ==2
                tempVar = eval(sprintf('LFPx.LFPx%d',ch));
                a = double(tempVar');
                ECoG.Fs(ch) = ECoG_Fs; 

            elseif SETTING ==3
                tempVar = eval(sprintf('Data.LFPx%d.data',ch));
                a = double(tempVar');
                ECoG.Fs(ch) = ECoG_Fs;
                
            elseif SETTING == 4
                tempVar = data(ch+3).dat;
                a = double(tempVar');
                ECoG.Fs(ch) = ECoG_Fs;
                
            end


            ECoG.data(ch,:)=resample(a, 1000, ECoG.Fs(ch))';
        end
        
        % for ch = 1:32
        %     tempVar = eval(sprintf('PDes_RA16_chn%03d',ch));
        %     a = double(tempVar.dat');
        %     ECoG.Fs(ch) = mean(tempVar.fs);
        %     ECoG.data(ch,:)=resample(a,1000,ECoG.Fs(ch))';
        % end
        
        if SETTING ==1
            tempVar =  eval('Inpt_RZ2_chn001');
            ECoG.Speaker = tempVar.dat';
            ECoG.Fs(size(ECoG.data,1)+1) = mean(tempVar.fs);
            tempVar =  eval('Inpt_RZ2_chn002');
            ECoG.Trigger = tempVar.dat';
            ECoG.Fs(size(ECoG.data,1)+2) = mean(tempVar.fs);
            tempVar =  eval('Inpt_RZ2_chn003');
            ECoG.Mic = tempVar.dat';
            ECoG.Fs(size(ECoG.data,1)+3) = mean(tempVar.fs);
        elseif SETTING ==2
            tempVar =  eval('Inpt.Inpt2');
            ECoG.Speaker = tempVar';
            ECoG.Fs(size(ECoG.data,1)+1) = Audio_Fs;
            tempVar =  eval('Inpt.Inpt1');
            ECoG.Trigger = tempVar';
            ECoG.Fs(size(ECoG.data,1)+2) = Audio_Fs;
            tempVar =  eval('Inpt.Inpt3');
            ECoG.Mic = tempVar';
            ECoG.Fs(size(ECoG.data,1)+3) = Audio_Fs;

        elseif SETTING ==3
            tempVar =  eval('Data.Inpt2.data');
            ECoG.Speaker = tempVar';
            ECoG.Fs(size(ECoG.data,1)+1) = Audio_Fs;
            tempVar =  eval('Data.Inpt1.data');
            ECoG.Trigger = tempVar';
            ECoG.Fs(size(ECoG.data,1)+2) = Audio_Fs;
            tempVar =  eval('Data.Inpt3.data');
            ECoG.Mic = tempVar';
            ECoG.Fs(size(ECoG.data,1)+3) = Audio_Fs;

        elseif SETTING == 4
            ECoG.Speaker = data(1).dat;
            ECoG.Fs(end+1) = Audio_Fs; % data(1).fs;
            ECoG.Trigger = data(2).dat;
            ECoG.Fs(end+1) = Audio_Fs; % data(2).fs;
            ECoG.Mic = data(3).dat;
            ECoG.Fs(end+1) = Audio_Fs; % data(3).fs;
            
        end
        
        ECoGALL{1,ds} = ECoG;
    end
    
end




