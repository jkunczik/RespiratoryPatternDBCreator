function [ds_section] = ComputeSectionFeatures(ds_section, modality, label, fs)
     for i=1:length(ds_section.(modality).(label))
         if ( ~isfield(ds_section.(modality).(label)(i), 'rr') || ...
             isempty(ds_section.(modality).(label)(i).rr))
                    [rr, a_rr, var_rr, var_a_rr, w_rr] =...
                  BreathFeatures(ds_section.(modality).(label)(i).sig, fs);

             %time parameters
             ds_section.(modality).(label)(i).w_rr = int16(w_rr*1e3);
             ds_section.(modality).(label)(i).rr = int16(rr*1e3);
             ds_section.(modality).(label)(i).a_rr = int16(a_rr*1e3);
             ds_section.(modality).(label)(i).var_rr = int16(var_rr*1e2);
             ds_section.(modality).(label)(i).var_a_rr = int16(var_a_rr*1e5);
         end

         if ~isfield(ds_section.(modality).(label)(i), 'rr_med') || ...
             isempty(ds_section.(modality).(label)(i).rr_med)

             %statistical parameters
             rr = double(ds_section.(modality).(label)(i).rr)/1e3;
             a_rr = double(ds_section.(modality).(label)(i).a_rr)/1e3;
    
             ds_section.(modality).(label)(i).rr_med = nanmedian(rr(rr ~= 0));
             ds_section.(modality).(label)(i).rr_var = nanvar(rr);
             ds_section.(modality).(label)(i).a_med  = nanmedian(a_rr(rr ~= 0));
             ds_section.(modality).(label)(i).a_var  = nanvar(a_rr);
    
             if isnan(ds_section.(modality).(label)(i).rr_med)
                ds_section.(modality).(label)(i).rr_med = 0;
             end
             if isnan(ds_section.(modality).(label)(i).rr_var)
                ds_section.(modality).(label)(i).rr_var = 0;
             end
             if isnan(ds_section.(modality).(label)(i).a_med)
                ds_section.(modality).(label)(i).a_med = 0;
             end
              if isnan(ds_section.(modality).(label)(i).a_var)
                ds_section.(modality).(label)(i).a_var = 0;
              end
         end
     end
end