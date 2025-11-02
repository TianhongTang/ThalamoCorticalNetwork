idx_file = ['../GLM_data/Sl_Thalamus_Muscimol.mat'];
load(idx_file, "idxUncI_posi", "idxUncN_posi", "idxUncI_nega", "idxUncN_nega", ...
    "unitID_all");
idx_all = [idxUncI_posi, idxUncN_posi, idxUncI_nega, idxUncN_nega];
for session = 1:10
    session_suffix = '-001_';

    raster_file = ['../GLM_data/MuscimolPre_full/raster_MuscimolPre_full_', int2str(session), ...
        '_0.mat'];
    load(raster_file, "cell_area", "cell_id", "session_name_full");
    N = length(cell_area);
    cell_type = zeros(1, N);
    for i=1:N
        if strcmp(cell_area{i},"Thalamus")
            cell_id_full = [session_name_full(1:8), session_suffix, cell_id{i}];
            idx_in_table = find(strcmp(unitID_all, cell_id_full));
            if ~isempty(idx_in_table)
                % calc value
                % posi_I = 101, posi_N = 102, posi_IN=103
                % nega_I = 110, nega_N = 120, nega_IN=130
                cell_type(i) = idxUncI_posi(idx_in_table)+...
                    idxUncN_posi(idx_in_table)*2+...
                    idxUncI_nega(idx_in_table)*10+...
                    idxUncN_nega(idx_in_table)*20+100;
            end
        end
    end
    type_file = ['../GLM_data/MuscimolPre_full/celltype_MuscimolPre_full_', int2str(session), ...
        '.mat'];
    save(type_file, "cell_type");
end