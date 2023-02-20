function [labels] = LabelSignal(sig, fs)
    
    sections = [...
        Pattern.BRADYPNEA, ...
        Pattern.EUPNEA, ...
        Pattern.HYPOPNEA, ...
        Pattern.HYPERPNEA, ...
        Pattern.TACHYPNEA, ...
        Pattern.EUPNEA, ...
        Pattern.KUSSMAUL, ...
        Pattern.CHEYNE_STOKES, ...
        Pattern.TACHYPNEA, ...
        Pattern.BIOT, ...
        Pattern.BIOT, ...
        Pattern.APNEA, ...
    ];

    colors = {
        [0 0.4470 0.7410], ...
        [0.8500 0.3250 0.0980], ...
        [0.9290 0.6940 0.1250], ...
        [0.4940 0.1840 0.5560], ...
        [0.4660 0.6740 0.1880], ...
        [0.3010 0.7450 0.9330], ...
        [0.6350 0.0780 0.1840], ...
        'red', ...
        'green', ...
        'blue', ...
        'cyan', ...
        'magenta', ...
    };

    labels = uint8(ones(length(sig),1))*Pattern.UNDEFINED;
    
    %% label data
    ShowSignal(sig, fs, 'FullScreen', true, 'PlotParams', 'k');
    t = 0:1/fs:(length(sig)-1)*fs;
    limits = ylim;
    hold on;
    for i=1:length(sections)
        label = sections(i);
        title(label.string);

        success = false;
        while ~success
            [t_lim, ~] = ginput(2);
            ids = find(t>=t_lim(1) & t<=t_lim(2));

            h_fig = figure();
            plot(t(ids), sig(ids));
            answ = questdlg('Segmentation okay?', 'Question', 'No', 'Yes', 'Yes');
            if strcmpi(answ, 'Yes')
                success = true;
            end

            close(h_fig);
        end  
        
        labels(ids) = label;

        plot(t_lim(1)*ones(2,1),limits, 'Color', colors{i});
        plot(t_lim(2)*ones(2,1),limits, 'Color', colors{i});
        plot(t(ids), sig(ids), 'Color', colors{i});
        text(mean(t_lim), limits(2), label.string, ...
            'VerticalAlignment', 'top',...
            'HorizontalAlignment', 'center',...
            'Color', colors{i});
    end
    close(gcf);
end