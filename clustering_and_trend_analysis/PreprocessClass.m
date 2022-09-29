classdef PreprocessClass
    
    properties
        data
        num_chans
        num_samples
        Fs  
        pre_samples
        post_samples
        chan_ids
    end
    
    methods
        
        function obj = PreprocessClass(data, chan_ids, Fs)
            obj.data = data;
            obj.num_chans = size(data, 1);
            obj.num_samples = size(data, 2);
            obj.Fs = Fs;
            obj.pre_samples = Fs/2;
            obj.post_samples = Fs/2; % Includes onset sample
            obj.chan_ids = chan_ids;
            
            if size(data,1) ~= length(chan_ids)
                error('Not correctly labelled');
            end
        end
        
        function data_hilbert = HilbertTransform(obj, freqs, num_freq_bins)
           
            error('This isn"t the correct spacing');
            if nargin < 2
                fprintf('Defaulting to High Gamma Analysis\n');
                num_freq_bins = 8;
                freqs = [70, 150];
            end
            min_freq = freqs(1);
            max_freq = freqs(2);

            CentralFreq = logspace(log10(min_freq), log10(max_freq), num_freq_bins);
            bw = diff(CentralFreq)/2;
            lf = [bw(1) bw];
            hf = [bw bw(end)];

            Hd = cell(1, num_freq_bins);
            for i = 1:num_freq_bins
                Fst1 = CentralFreq(i) - lf(i) - 1;
                Fp1 = CentralFreq(i) - lf(i);
                Fp2 = CentralFreq(i) + hf(i);
                Fst2 = CentralFreq(i) + hf(i) + 1;
                Ast1 = 60;
                Ap = 1;
                Ast2 = 60;
                H = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', Fst1, Fp1, Fp2, Fst2, Ast1, Ap, Ast2, obj.Fs);
                Hd{i} = design(H, 'butter');
            end

            data_hilbert = nan(size(obj.data,1), obj.num_samples);
            for ch = 1:size(obj.data, 1)

                tempSig = [zeros(1, 10000) obj.data(ch,:) zeros(1, 10000)];

                hilbert_tempVar = nan(num_freq_bins, obj.num_samples);
                for i = 1:num_freq_bins,

                    tempVar= filtfilt(Hd{i}.sosMatrix, Hd{i}.ScaleValues, tempSig);
                    tempVar1 = tempVar(10001:end-10000);
                    hilbert_tempVar(i,:) = abs(hilbert(tempVar1));

                end
                data_hilbert(ch,:) = mean(hilbert_tempVar, 1);
            end
            
        end % HilbertTransform

        function x = z_score(obj, x, baseline_inds)

            % Z score
            x = (x - repmat(nanmean(x(:, baseline_inds, :), 2), 1, size(x,2))) ...
                ./ (repmat(nanstd(x(:, baseline_inds, :), [], 2), 1, size(x,2)));

        end
        
        function x = z_score_with_bl_data(obj, x, baseline_data)
            
            % Z score
            x = (x - repmat(nanmean(baseline_data, 2), 1, size(x,2))) ...
                ./ (repmat(nanstd(baseline_data, [], 2), 1, size(x,2)));
            
        end
            

    end % methods
end

