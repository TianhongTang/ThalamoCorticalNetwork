% Generate simulated data for the GLM model
% Use inferred post data
% Validate possible model structures
% UseGPU = gpuDeviceCount > 0;

clear;
%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main

session_types = {'Muscimol'};
session_nums = {10, 5}; % number of sessions for each type
states = {'RestOpen', 'RestClose', 'Task'};
for session_type_idx = 1:length(session_types)
    session_type = session_types{session_type_idx};
    session_num = session_nums{session_type_idx};
    for session_idx = 1:session_num
        session_name_full = sprintf('%s_Session%d', session_type, session_idx);
        for state_idx = 1:length(states)
            state = states{state_idx};
            simulation_from_post(state, session_idx);
        end
    end
end