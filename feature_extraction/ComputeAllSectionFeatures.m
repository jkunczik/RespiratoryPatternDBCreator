function ComputeAllSectionFeatures(dataset_dir, config)
%COMPUTEFEATURES Computes the features for all section of the given data set

    load(fullfile(dataset_dir, 'sections.mat'));

    sig_names = fieldnames(sections);
    sig_names(contains(sig_names, 'REF')) = [];
    for i=1:length(sig_names)
        sec_names = fieldnames(sections.(sig_names{i}));

        h_fig = figure;
        title([sig_names{i} ' dataset']);
        xlabel('Respiration Rate [BPM]');
        ylabel('Respiration Rate Variance [BPM]');
        zlabel('Amplitude [n.u]')
        view(0,0);
        legend();
        
        for j=1:length(sec_names)
            sections = ComputeSectionFeatures(sections, sig_names{i}, sec_names{j}, config.fs);
            plot_ds_section(sections, sig_names{i}, sec_names{j}, h_fig);
        end 
    end

    sections = add_ref_features(sections);
    
    mkdir(config.dirs.aug_ds_dir)
    save(fullfile(dataset_dir, 'sections.mat'), 'sections');
end



function [fig] = plot_ds_section(ds_section, modality, label, fig)
     if nargin >= 4
        figure(fig);
     else
        fig = figure();
     end

     label_str = label;
     replace(label_str, '_', ' ');

     hold on;
     scatter3( ...
         [ds_section.(modality).(label)(:).rr_med], ...
         [ds_section.(modality).(label)(:).rr_var], ...
         [ds_section.(modality).(label)(:).a_med], ...
         'DisplayName', label ...
     );

     if nargin < 4
         title('DataSet')
         xlabel('Respiration Rate [BPM]');
         ylabel('Respiration Rate Variance [BPM]');
         zlabel('Amplitude [n.u]')
         view(0,0);
         legend();
     end
end

function [sections] = add_ref_features(sections)
    sec_names = fieldnames(sections.REF);
    for k=1:length(sec_names)
        for i=1:length(sections.REF.(sec_names{k}))
            rr = sections.REF.(sec_names{k})(i).rr;
            a_rr = sections.REF.(sec_names{k})(i).a_rr;
            sections.REF.(sec_names{k})(i).rr_med = nanmedian(rr(rr ~= 0));
            sections.REF.(sec_names{k})(i).rr_var = nanvar(rr);
            sections.REF.(sec_names{k})(i).a_med = nanmedian(a_rr(rr ~= 0));
            sections.REF.(sec_names{k})(i).a_var = nanvar(a_rr);
             if isnan(sections.REF.(sec_names{k})(i).rr_med)
                sections.REF.(sec_names{k})(i).rr_med = 0;
             end
             if isnan(sections.REF.(sec_names{k})(i).rr_var)
                sections.REF.(sec_names{k})(i).rr_var = 0;
             end
             if isnan(sections.REF.(sec_names{k})(i).a_med)
                sections.REF.(sec_names{k})(i).a_med = 0;
             end
              if isnan(sections.REF.(sec_names{k})(i).a_var)
                sections.REF.(sec_names{k})(i).a_var = 0;
             end
        end
    end
end

