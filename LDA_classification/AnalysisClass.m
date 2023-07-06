classdef AnalysisClass < handle
    
    properties
        data
        num_chans
        num_samples
        Fs
        fig        
        event_sample
        erp = [];
        erp_sd = [];
        chan_ids
    end
    
    methods
        function obj = AnalysisClass(data, chan_ids, Fs)
            obj.data = data;
            obj.fig = struct([]);
            obj.num_chans = size(data, 1);
            obj.num_samples = size(data, 2);
            obj.Fs = Fs;
            obj.event_sample = Fs/2;
            obj.chan_ids = chan_ids;
            
            if size(data, 1) ~= length(chan_ids)
                error('These should be the same');
            end
        end
        
        function fig_ = SetupGridPlot(obj, chan_inds, n_rows, n_cols, order, labels)
            
            if size(chan_inds, 2) > size(chan_inds, 1),
                chan_inds = chan_inds';
            end
            
            if nargin <= 3,
                % TODO: Make these a function of chan_inds
                n_rows = ceil(sqrt(length(chan_inds)));
                n_cols = n_rows;
                order = 'across';
                nLabels = n_cols*n_rows;
                labels = [chan_inds; nan(nLabels-length(chan_inds),1)];
            elseif nargin <= 5,
                nLabels = n_cols*n_rows;
                labels = [chan_inds; nan(nLabels-length(chan_inds),1)];
            end
            
            addpath(genpath('bounLine'));
            
            fig_handle = figure('Position', [1,26,1680,919], 'color', [1 1 1]);
            
            % Starts top left and goes down left side
            down_1st = reshape(labels, n_rows, n_cols);
            % Starts top left and goes across top 
            across_1st = reshape(labels, n_cols, n_rows)';
            switch order,
                case 'across',
                    fig_index = across_1st';
                case 'down',
                    fig_index = down_1st';
            end
            
            % Setup subplots
            edge_space = 0.07;
            all_plt_space_width = 1-edge_space*2;
            plt_width = all_plt_space_width/n_cols*.8;
            plt_height = all_plt_space_width/n_rows*.8;
            
            col_spacing = (all_plt_space_width-plt_width)/(n_cols-1);
            col_offset = edge_space;
            row_spacing = (all_plt_space_width-plt_height)/(n_rows-1);
            row_offset = (1-edge_space-plt_height);
            
            aH_map = across_1st';
            aH = zeros(size(aH_map));
            for plt=1:numel(aH),
                
                plt_x_start = col_offset + rem(plt-1, n_cols) * col_spacing;
                plt_y_end = row_offset - floor((plt-1) / n_cols) * row_spacing;

                aH(plt) = axes('Units', 'normalized', 'Position', ...
                             [plt_x_start, plt_y_end, plt_width plt_height], 'Visible', 'off');
            end
            
            fig_ = struct('handle', fig_handle, 'plot_ref', fig_index, 'plot_axes_handles', aH, 'axes_ref', aH_map, 'n_rows', n_rows, 'n_cols', n_cols);

        end
        
        function [erp, erp_sd, obj] = GetERP(obj, percent)
                
            if ~isempty(obj.erp),
                erp = obj.erp;
            else
                
                % Percentage bounds to take ERP over
                if nargin < 2,
                    percent = 10;
                end

                % Calculate ERP
                erp =  trimmean(obj.data, percent, 3);
                erp_sd = nanstd(obj.data, [], 3)/ sqrt(size(obj.data, 3));

                % Smooth the trajectories
                for ch = 1:obj.num_chans,
                    erp(ch,:) = smooth(erp(ch,:), .1, 'lowess');
                end
                
                obj.erp = erp;
                obj.erp_sd = erp_sd;
                
            end
            
        end
            
        function [fig] = PlotERP(obj, varargin)
            
            % Inputs
            ip = inputParser;
            
            defaultLabels = 1:obj.num_chans;
            addOptional(ip, 'chan_labels', defaultLabels);
            
            parse(ip, varargin{:});
            
            chan_labels = ip.Results.chan_labels;
            
            [erp, erp_sd] = obj.GetERP();
            
            fig_ = obj.SetupGridPlot((1:obj.num_chans)');

            for  chan = fig_.plot_ref(:)',

                if isnan(chan),
                    continue;
                end
                
                aH = fig_.plot_axes_handles(chan);
                axes(aH);
                
                box(aH, 'off');
                hold(aH, 'all');
                line([obj.event_sample obj.event_sample], [-10 10], 'LineWidth', 3, 'Color', [1 1 1]/3);

                if ~isempty(erp(chan,:)),

                    [l1, p1] = boundedline(1:obj.num_samples, erp(chan,:), ...
                        erp_sd(chan,:), 'r'); % FIXME: Make 2nd ERP standard error
                    set(l1, 'LineWidth', 2);
                    outlinebounds(l1, p1);
                    hold on;

                end

                set(gca, 'YTick', [0 2 4 6], 'YTickLabel', {'','','',''}, ....
                    'XTick', [0 obj.event_sample obj.num_samples], ...
                    'XTickLabel', {'','','','',''}, 'TickLength', [0.025 0.025]);
                
                minLim = -.5 + nanmin(erp(:));
                maxLim =  .5 + nanmax(erp(:));
                axis([0 obj.num_samples minLim maxLim]) ;
                title(aH, num2str(chan_labels(chan)), 'FontSize', 14, 'FontWeight', 'bold', 'Visible', 'on');
            end

            fig = fig_;

        end
        
        function PlotTrialReponseImage(obj, varargin)
            
            % Inputs
            ip = inputParser;
            
            defaultLabels = 1:obj.num_chans;
            addOptional(ip, 'chan_labels', defaultLabels);
            
            parse(ip, varargin{:});
            
            chan_labels = ip.Results.chan_labels;
            
            fig_ = obj.SetupGridPlot((1:obj.num_chans)');
            
%             erp = obj.GetERP();
%             minLim = -.5 + min(erp(:));
%             maxLim =  .5 + max(erp(:));
            fprintf('Hack: Fixing the bounds...\n');
            minLim = -5;
            maxLim =  5;
            
            for  plt = 1:obj.num_chans,
                
                aH = fig_.plot_axes_handles(plt);
                axes(aH);

                tempData = squeeze(obj.data(fig_.chan(plt),:,:))';
                imagesc(tempData,[minLim maxLim]);
                
                set(gca, 'YTick', [], 'XTick', [0 obj.event_sample obj.num_samples], ...
                    'XTickLabel', {'','',''}, 'TickLength', [0.025 0.025]);
                
                title(aH, num2str(chan_labels(plt)),'FontSize',14, 'FontWeight','bold', 'Visible', 'on');
            end

        end
        
    end
    
end

