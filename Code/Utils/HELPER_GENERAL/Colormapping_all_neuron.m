function Colormapping_all_neuron()
clear all
close all

folder_name = 'C:\Documents and Settings\tachibanay\desktop\Huey analysis';
data_dir = [folder_name, '\Rex_analyze\Rex_PSTH_mat\'];
result_dir = [folder_name, '\Rex_analyze\Rex_PSTH_mat\'];
figure_dir = [folder_name, '\Rex_analyze\Rex_PSTH_Fig'];

%%%%sort setting%%%%%
sort_timing_mode = 2;%0=original;1*=fix;2*=cue;3*=sac;4*=rew;
sort_condition_mode = 2;%0=original;*1=LS;*2=LB;*3=RB;*4=RS;

sort_FIX_starttime = 100;
sort_FIX_endtime = 400;
sort_FIX_duration = sort_FIX_endtime-sort_FIX_starttime;
sort_CUE_starttime = 100;
sort_CUE_endtime = 400;
sort_CUE_duration = sort_CUE_endtime-sort_CUE_starttime;
sort_SAC_starttime = 0;
sort_SAC_endtime = 300;
sort_SAC_duration = sort_SAC_endtime-sort_SAC_starttime;
sort_REW_starttime = 0;
sort_REW_endtime = 500;
sort_REW_duration = sort_REW_endtime-sort_REW_starttime;

%%%%initial setting%%%%%
window_duration = 100;
window_step = 20;
max_cutoff = 2;

analyze_base_starttime = -1300;%unchanged
analyze_base_endtime = -300;%unchanged
analyze_base_duration = analyze_base_endtime - analyze_base_starttime;

analyze_FIX_starttime = -200;
analyze_FIX_endtime = 400;
analyze_FIX_duration = analyze_FIX_endtime - analyze_FIX_starttime;
analyze_FIX_Timecourse = [analyze_FIX_starttime:window_step:analyze_FIX_endtime];

analyze_CUE_starttime = -300;
analyze_CUE_endtime = 900;
analyze_CUE_duration = analyze_CUE_endtime - analyze_CUE_starttime;
analyze_CUE_Timecourse = [analyze_CUE_starttime:window_step:analyze_CUE_endtime];

analyze_SAC_starttime = -300;
analyze_SAC_endtime = 300;
analyze_SAC_duration = analyze_SAC_endtime - analyze_SAC_starttime;
analyze_SAC_Timecourse = [analyze_SAC_starttime:window_step:analyze_SAC_endtime];

analyze_REW_starttime = -200;
analyze_REW_endtime = 1000;
analyze_REW_duration = analyze_REW_endtime - analyze_REW_starttime;
analyze_REW_Timecourse = [analyze_REW_starttime:window_step:analyze_REW_endtime];

subplot_x = 0.0002;
subplot_y_height = 0.15;
subplot_spacex = 0.015;

%%%%data pickup%%%%%
cd(data_dir)

[srtlist] = fuf (data_dir,'detail');
srtlist=srtlist';
X=1;
while (X<=length(srtlist))
    XX=strfind(srtlist{1,X},'m');
    if (~isempty(XX))
        Empty(X)=1;
    else
        Empty(X)=0;
    end
    X=X+1;
end
Y=find(Empty==0);
srtlist(Y)=[];
srtlist=srtlist';

ZZ=1;
YY=strfind(srtlist{ZZ,1},'\');
YY_max=max(YY);
while ZZ<=length(srtlist)
    srtlist{ZZ,1}(1:YY_max)=[];
    ZZ =ZZ+1;
end

ZZ=1;
while ZZ<=length(srtlist)
    YY=strfind(srtlist{ZZ,1},'.mat');
    YY_max=max(YY);
    if~isempty(YY_max)
        srtlist{ZZ,1}(YY_max:YY_max+3)=[];
    end
    ZZ =ZZ+1;
end

keyboard%1000 is arbitrary value

List=1;
List_analyze_LS_FIX = inf*ones(length(srtlist),1000);
List_analyze_LS_CUE = inf*ones(length(srtlist),1000);
List_analyze_LS_SAC = inf*ones(length(srtlist),1000);
List_analyze_LS_REW = inf*ones(length(srtlist),1000);
List_analyze_RB_FIX = inf*ones(length(srtlist),1000);
List_analyze_RB_CUE = inf*ones(length(srtlist),1000);
List_analyze_RB_SAC = inf*ones(length(srtlist),1000);
List_analyze_RB_REW = inf*ones(length(srtlist),1000);
List_analyze_LB_FIX = inf*ones(length(srtlist),1000);
List_analyze_LB_CUE = inf*ones(length(srtlist),1000);
List_analyze_LB_SAC = inf*ones(length(srtlist),1000);
List_analyze_LB_REW = inf*ones(length(srtlist),1000);
List_analyze_RS_FIX = inf*ones(length(srtlist),1000);
List_analyze_RS_CUE = inf*ones(length(srtlist),1000);
List_analyze_RS_SAC = inf*ones(length(srtlist),1000);
List_analyze_RS_REW = inf*ones(length(srtlist),1000);
neuron_num = length(srtlist);

while (List<=neuron_num)
    indicator(List,neuron_num);
    load_data_name = srtlist{List,1};
    eval(['load ' load_data_name]);

    %%%%analyze%%%%%
    %%FIX%%
    Sum_cor_LS_FIX = cor_LS_FIX_Spike_list;
    Sum_cor_LS_FIX =Sum_cor_LS_FIX(:);
    i = 0;
    while i<= ((analyze_FIX_endtime-analyze_FIX_starttime)/window_step)
        analyze_LS_FIX(i+1) = length(find((analyze_FIX_starttime+(i)*window_step) < Sum_cor_LS_FIX &...
            Sum_cor_LS_FIX <= (analyze_FIX_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LS_FIX_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LS_FIX = analyze_LS_FIX/mean(cor_LS_BASE_Spike);
    analyze_LS_FIX_max_finder = find(analyze_LS_FIX>max_cutoff);
    analyze_LS_FIX(analyze_LS_FIX_max_finder) = max_cutoff;

    Sum_cor_RB_FIX = cor_RB_FIX_Spike_list;
    Sum_cor_RB_FIX = Sum_cor_RB_FIX(:);
    i = 0;
    while i<= ((analyze_FIX_endtime-analyze_FIX_starttime)/window_step)
        analyze_RB_FIX(i+1) = length(find((analyze_FIX_starttime+(i)*window_step) < Sum_cor_RB_FIX &...
            Sum_cor_RB_FIX <= (analyze_FIX_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RB_FIX_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RB_FIX = analyze_RB_FIX/mean(cor_RB_BASE_Spike);
    analyze_RB_FIX_max_finder = find(analyze_RB_FIX>max_cutoff);
    analyze_RB_FIX(analyze_RB_FIX_max_finder) = max_cutoff;

    Sum_cor_LB_FIX = cor_LB_FIX_Spike_list;
    Sum_cor_LB_FIX = Sum_cor_LB_FIX(:);
    i = 0;
    while i<= ((analyze_FIX_endtime-analyze_FIX_starttime)/window_step)
        analyze_LB_FIX(i+1) = length(find((analyze_FIX_starttime+(i)*window_step) < Sum_cor_LB_FIX &...
            Sum_cor_LB_FIX <= (analyze_FIX_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LB_FIX_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LB_FIX = analyze_LB_FIX/mean(cor_LB_BASE_Spike);
    analyze_LB_FIX_max_finder = find(analyze_LB_FIX>max_cutoff);
    analyze_LB_FIX(analyze_LB_FIX_max_finder) = max_cutoff;

    Sum_cor_RS_FIX = cor_RS_FIX_Spike_list;
    Sum_cor_RS_FIX = Sum_cor_RS_FIX(:);
    i = 0;
    while i<= ((analyze_FIX_endtime-analyze_FIX_starttime)/window_step)
        analyze_RS_FIX(i+1) = length(find((analyze_FIX_starttime+(i)*window_step) < Sum_cor_RS_FIX &...
            Sum_cor_RS_FIX <= (analyze_FIX_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RS_FIX_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RS_FIX = analyze_RS_FIX/mean(cor_RS_BASE_Spike);
    analyze_RS_FIX_max_finder = find(analyze_RS_FIX>max_cutoff);
    analyze_RS_FIX(analyze_RS_FIX_max_finder) = max_cutoff;

    %%CUE%%
    Sum_cor_LS_CUE = cor_LS_CUE_Spike_list;
    Sum_cor_LS_CUE = Sum_cor_LS_CUE(:);
    i = 0;
    while i<= ((analyze_CUE_endtime-analyze_CUE_starttime)/window_step)
        analyze_LS_CUE(i+1) = length(find((analyze_CUE_starttime+(i)*window_step) < Sum_cor_LS_CUE &...
            Sum_cor_LS_CUE <= (analyze_CUE_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LS_CUE_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LS_CUE = analyze_LS_CUE/mean(cor_LS_BASE_Spike);
    analyze_LS_CUE_max_finder = find(analyze_LS_CUE>max_cutoff);
    analyze_LS_CUE(analyze_LS_CUE_max_finder) = max_cutoff;

    Sum_cor_RB_CUE = cor_RB_CUE_Spike_list;
    Sum_cor_RB_CUE = Sum_cor_RB_CUE(:);
    i = 0;
    while i<= ((analyze_CUE_endtime-analyze_CUE_starttime)/window_step)
        analyze_RB_CUE(i+1) = length(find((analyze_CUE_starttime+(i)*window_step) < Sum_cor_RB_CUE &...
            Sum_cor_RB_CUE <= (analyze_CUE_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RB_CUE_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RB_CUE = analyze_RB_CUE/mean(cor_RB_BASE_Spike);
    analyze_RB_CUE_max_finder = find(analyze_RB_CUE>max_cutoff);
    analyze_RB_CUE(analyze_RB_CUE_max_finder) = max_cutoff;

    Sum_cor_LB_CUE = cor_LB_CUE_Spike_list;
    Sum_cor_LB_CUE = Sum_cor_LB_CUE(:);
    i = 0;
    while i<= ((analyze_CUE_endtime-analyze_CUE_starttime)/window_step)
        analyze_LB_CUE(i+1) = length(find((analyze_CUE_starttime+(i)*window_step) < Sum_cor_LB_CUE &...
            Sum_cor_LB_CUE <= (analyze_CUE_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LB_CUE_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LB_CUE = analyze_LB_CUE/mean(cor_LB_BASE_Spike);
    analyze_LB_CUE_max_finder = find(analyze_LB_CUE>max_cutoff);
    analyze_LB_CUE(analyze_LB_CUE_max_finder) = max_cutoff;

    Sum_cor_RS_CUE = cor_RS_CUE_Spike_list;
    Sum_cor_RS_CUE = Sum_cor_RS_CUE(:);
    i = 0;
    while i<= ((analyze_CUE_endtime-analyze_CUE_starttime)/window_step)
        analyze_RS_CUE(i+1) = length(find((analyze_CUE_starttime+(i)*window_step) < Sum_cor_RS_CUE &...
            Sum_cor_RS_CUE <= (analyze_CUE_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RS_CUE_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RS_CUE = analyze_RS_CUE/mean(cor_RS_BASE_Spike);
    analyze_RS_CUE_max_finder = find(analyze_RS_CUE>max_cutoff);
    analyze_RS_CUE(analyze_RS_CUE_max_finder) = max_cutoff;

    %%SAC%%
    Sum_cor_LS_SAC = cor_LS_SAC_Spike_list;
    Sum_cor_LS_SAC = Sum_cor_LS_SAC(:);
    i = 0;
    while i<= ((analyze_SAC_endtime-analyze_SAC_starttime)/window_step)
        analyze_LS_SAC(i+1) = length(find((analyze_SAC_starttime+(i)*window_step) < Sum_cor_LS_SAC &...
            Sum_cor_LS_SAC <= (analyze_SAC_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LS_SAC_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LS_SAC = analyze_LS_SAC/mean(cor_LS_BASE_Spike);
    analyze_LS_SAC_max_finder = find(analyze_LS_SAC>max_cutoff);
    analyze_LS_SAC(analyze_LS_SAC_max_finder) = max_cutoff;

    Sum_cor_RB_SAC = cor_RB_SAC_Spike_list;
    Sum_cor_RB_SAC = Sum_cor_RB_SAC(:);
    i = 0;
    while i<= ((analyze_SAC_endtime-analyze_SAC_starttime)/window_step)
        analyze_RB_SAC(i+1) = length(find((analyze_SAC_starttime+(i)*window_step) < Sum_cor_RB_SAC &...
            Sum_cor_RB_SAC <= (analyze_SAC_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RB_SAC_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RB_SAC = analyze_RB_SAC/mean(cor_RB_BASE_Spike);
    analyze_RB_SAC_max_finder = find(analyze_RB_SAC>max_cutoff);
    analyze_RB_SAC(analyze_RB_SAC_max_finder) = max_cutoff;

    Sum_cor_LB_SAC = cor_LB_SAC_Spike_list;
    Sum_cor_LB_SAC = Sum_cor_LB_SAC(:);
    i = 0;
    while i<= ((analyze_SAC_endtime-analyze_SAC_starttime)/window_step)
        analyze_LB_SAC(i+1) = length(find((analyze_SAC_starttime+(i)*window_step) < Sum_cor_LB_SAC &...
            Sum_cor_LB_SAC <= (analyze_SAC_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LB_SAC_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LB_SAC = analyze_LB_SAC/mean(cor_LB_BASE_Spike);
    analyze_LB_SAC_max_finder = find(analyze_LB_SAC>max_cutoff);
    analyze_LB_SAC(analyze_LB_SAC_max_finder) = max_cutoff;

    Sum_cor_RS_SAC = cor_RS_SAC_Spike_list;
    Sum_cor_RS_SAC = Sum_cor_RS_SAC(:);
    i = 0;
    while i<= ((analyze_SAC_endtime-analyze_SAC_starttime)/window_step)
        analyze_RS_SAC(i+1) = length(find((analyze_SAC_starttime+(i)*window_step) < Sum_cor_RS_SAC &...
            Sum_cor_RS_SAC <= (analyze_SAC_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RS_SAC_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RS_SAC = analyze_RS_SAC/mean(cor_RS_BASE_Spike);
    analyze_RS_SAC_max_finder = find(analyze_RS_SAC>max_cutoff);
    analyze_RS_SAC(analyze_RS_SAC_max_finder) = max_cutoff;

    %%REW%%
    Sum_cor_LS_REW = cor_LS_REW_Spike_list;
    Sum_cor_LS_REW = Sum_cor_LS_REW(:);
    i = 0;
    while i<= ((analyze_REW_endtime-analyze_REW_starttime)/window_step)
        analyze_LS_REW(i+1) = length(find((analyze_REW_starttime+(i)*window_step) < Sum_cor_LS_REW &...
            Sum_cor_LS_REW <= (analyze_REW_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LS_REW_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LS_REW = analyze_LS_REW/mean(cor_LS_BASE_Spike);
    analyze_LS_REW_max_finder = find(analyze_LS_REW>max_cutoff);
    analyze_LS_REW(analyze_LS_REW_max_finder) = max_cutoff;

    Sum_cor_RB_REW = cor_RB_REW_Spike_list;
    Sum_cor_RB_REW = Sum_cor_RB_REW(:);
    i = 0;
    while i<= ((analyze_REW_endtime-analyze_REW_starttime)/window_step)
        analyze_RB_REW(i+1) = length(find((analyze_REW_starttime+(i)*window_step) < Sum_cor_RB_REW &...
            Sum_cor_RB_REW <= (analyze_REW_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RB_REW_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RB_REW = analyze_RB_REW/mean(cor_RB_BASE_Spike);
    analyze_RB_REW_max_finder = find(analyze_RB_REW>max_cutoff);
    analyze_RB_REW(analyze_RB_REW_max_finder) = max_cutoff;

    Sum_cor_LB_REW = cor_LB_REW_Spike_list;
    Sum_cor_LB_REW = Sum_cor_LB_REW(:);
    i = 0;
    while i<= ((analyze_REW_endtime-analyze_REW_starttime)/window_step)
        analyze_LB_REW(i+1) = length(find((analyze_REW_starttime+(i)*window_step) < Sum_cor_LB_REW &...
            Sum_cor_LB_REW <= (analyze_REW_starttime+(i)*window_step+window_duration)))*1000/(size(cor_LB_REW_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_LB_REW = analyze_LB_REW/mean(cor_LB_BASE_Spike);
    analyze_LB_REW_max_finder = find(analyze_LB_REW>max_cutoff);
    analyze_LB_REW(analyze_LB_REW_max_finder) = max_cutoff;

    Sum_cor_RS_REW = cor_RS_REW_Spike_list;
    Sum_cor_RS_REW = Sum_cor_RS_REW(:);
    i = 0;
    while i<= ((analyze_REW_endtime-analyze_REW_starttime)/window_step)
        analyze_RS_REW(i+1) = length(find((analyze_REW_starttime+(i)*window_step) < Sum_cor_RS_REW &...
            Sum_cor_RS_REW <= (analyze_REW_starttime+(i)*window_step+window_duration)))*1000/(size(cor_RS_REW_Spike_list,1)*window_duration);
        i = i + 1;
    end
    analyze_RS_REW = analyze_RS_REW/mean(cor_RS_BASE_Spike);
    analyze_RS_REW_max_finder = find(analyze_RS_REW>max_cutoff);
    analyze_RS_REW(analyze_RS_REW_max_finder) = max_cutoff;

    List_analyze_LS_FIX(List,1:size(analyze_LS_FIX,2)) = analyze_LS_FIX;
    List_analyze_LS_CUE(List,1:size(analyze_LS_CUE,2)) = analyze_LS_CUE;
    List_analyze_LS_SAC(List,1:size(analyze_LS_SAC,2)) = analyze_LS_SAC;
    List_analyze_LS_REW(List,1:size(analyze_LS_REW,2)) = analyze_LS_REW;
    List_analyze_RB_FIX(List,1:size(analyze_RB_FIX,2)) = analyze_RB_FIX;
    List_analyze_RB_CUE(List,1:size(analyze_RB_CUE,2)) = analyze_RB_CUE;
    List_analyze_RB_SAC(List,1:size(analyze_RB_SAC,2)) = analyze_RB_SAC;
    List_analyze_RB_REW(List,1:size(analyze_RB_REW,2)) = analyze_RB_REW;
    List_analyze_LB_FIX(List,1:size(analyze_LB_FIX,2)) = analyze_LB_FIX;
    List_analyze_LB_CUE(List,1:size(analyze_LB_CUE,2)) = analyze_LB_CUE;
    List_analyze_LB_SAC(List,1:size(analyze_LB_SAC,2)) = analyze_LB_SAC;
    List_analyze_LB_REW(List,1:size(analyze_LB_REW,2)) = analyze_LB_REW;
    List_analyze_RS_FIX(List,1:size(analyze_RS_FIX,2)) = analyze_RS_FIX;
    List_analyze_RS_CUE(List,1:size(analyze_RS_CUE,2)) = analyze_RS_CUE;
    List_analyze_RS_SAC(List,1:size(analyze_RS_SAC,2)) = analyze_RS_SAC;
    List_analyze_RS_REW(List,1:size(analyze_RS_REW,2)) = analyze_RS_REW;

    List = List+1;
end

maxX_analyze_LS_FIX = size(analyze_LS_FIX,2);
zero_analyze_LS_FIX = abs(analyze_FIX_starttime)/analyze_FIX_duration;
zeroX_analyze_LS_FIX = (maxX_analyze_LS_FIX-1)*zero_analyze_LS_FIX+0.5;
maxX_analyze_LS_CUE = size(analyze_LS_CUE,2);
zero_analyze_LS_CUE = abs(analyze_CUE_starttime)/analyze_CUE_duration;
zeroX_analyze_LS_CUE = (maxX_analyze_LS_CUE-1)*zero_analyze_LS_CUE+0.5;
maxX_analyze_LS_SAC = size(analyze_LS_SAC,2);
zero_analyze_LS_SAC = abs(analyze_SAC_starttime)/analyze_SAC_duration;
zeroX_analyze_LS_SAC = (maxX_analyze_LS_SAC-1)*zero_analyze_LS_SAC+0.5;
maxX_analyze_LS_REW = size(analyze_LS_REW,2);
zero_analyze_LS_REW = abs(analyze_REW_starttime)/analyze_REW_duration;
zeroX_analyze_LS_REW = (maxX_analyze_LS_REW-1)*zero_analyze_LS_REW+0.5;

maxX_analyze_LB_FIX = size(analyze_LB_FIX,2);
zero_analyze_LB_FIX = abs(analyze_FIX_starttime)/analyze_FIX_duration;
zeroX_analyze_LB_FIX = (maxX_analyze_LB_FIX-1)*zero_analyze_LB_FIX+0.5;
maxX_analyze_LB_CUE = size(analyze_LB_CUE,2);
zero_analyze_LB_CUE = abs(analyze_CUE_starttime)/analyze_CUE_duration;
zeroX_analyze_LB_CUE = (maxX_analyze_LB_CUE-1)*zero_analyze_LB_CUE+0.5;
maxX_analyze_LB_SAC = size(analyze_LB_SAC,2);
zero_analyze_LB_SAC = abs(analyze_SAC_starttime)/analyze_SAC_duration;
zeroX_analyze_LB_SAC = (maxX_analyze_LB_SAC-1)*zero_analyze_LB_SAC+0.5;
maxX_analyze_LB_REW = size(analyze_LB_REW,2);
zero_analyze_LB_REW = abs(analyze_REW_starttime)/analyze_REW_duration;
zeroX_analyze_LB_REW = (maxX_analyze_LB_REW-1)*zero_analyze_LB_REW+0.5;

maxX_analyze_RB_FIX = size(analyze_RB_FIX,2);
zero_analyze_RB_FIX = abs(analyze_FIX_starttime)/analyze_FIX_duration;
zeroX_analyze_RB_FIX = (maxX_analyze_RB_FIX-1)*zero_analyze_RB_FIX+0.5;
maxX_analyze_RB_CUE = size(analyze_RB_CUE,2);
zero_analyze_RB_CUE = abs(analyze_CUE_starttime)/analyze_CUE_duration;
zeroX_analyze_RB_CUE = (maxX_analyze_RB_CUE-1)*zero_analyze_RB_CUE+0.5;
maxX_analyze_RB_SAC = size(analyze_RB_SAC,2);
zero_analyze_RB_SAC = abs(analyze_SAC_starttime)/analyze_SAC_duration;
zeroX_analyze_RB_SAC = (maxX_analyze_RB_SAC-1)*zero_analyze_RB_SAC+0.5;
maxX_analyze_RB_REW = size(analyze_RB_REW,2);
zero_analyze_RB_REW = abs(analyze_REW_starttime)/analyze_REW_duration;
zeroX_analyze_RB_REW = (maxX_analyze_RB_REW-1)*zero_analyze_RB_REW+0.5;

maxX_analyze_RS_FIX = size(analyze_RS_FIX,2);
zero_analyze_RS_FIX = abs(analyze_FIX_starttime)/analyze_FIX_duration;
zeroX_analyze_RS_FIX = (maxX_analyze_RS_FIX-1)*zero_analyze_RS_FIX+0.5;
maxX_analyze_RS_CUE = size(analyze_RS_CUE,2);
zero_analyze_RS_CUE = abs(analyze_CUE_starttime)/analyze_CUE_duration;
zeroX_analyze_RS_CUE = (maxX_analyze_RS_CUE-1)*zero_analyze_RS_CUE+0.5;
maxX_analyze_RS_SAC = size(analyze_RS_SAC,2);
zero_analyze_RS_SAC = abs(analyze_SAC_starttime)/analyze_SAC_duration;
zeroX_analyze_RS_SAC = (maxX_analyze_RS_SAC-1)*zero_analyze_RS_SAC+0.5;
maxX_analyze_RS_REW = size(analyze_RS_REW,2);
zero_analyze_RS_REW = abs(analyze_REW_starttime)/analyze_REW_duration;
zeroX_analyze_RS_REW = (maxX_analyze_RS_REW-1)*zero_analyze_RS_REW+0.5;

realign_FIX = zeroX_analyze_LS_FIX + (sort_FIX_starttime/window_step)+0.5;
realign_CUE = zeroX_analyze_LS_CUE + (sort_CUE_starttime/window_step)+0.5;
realign_SAC = zeroX_analyze_LS_SAC + (sort_SAC_starttime/window_step)+0.5;
realign_REW = zeroX_analyze_LS_REW + (sort_REW_starttime/window_step)+0.5;

%%%%%%%%%%%data sorting%%%%%%%%%%%%
realign_sum = (1:1:neuron_num)';

if (sort_timing_mode == 1)
    if (sort_condition_mode == 1)
        realign_sum = List_analyze_LS_FIX(:,realign_FIX:realign_FIX+(sort_FIX_duration/window_step)-1);
    elseif (sort_condition_mode == 2)
        realign_sum = List_analyze_LB_FIX(:,realign_FIX:realign_FIX+(sort_FIX_duration/window_step)-1);
    elseif (sort_condition_mode == 3)
        realign_sum = List_analyze_RB_FIX(:,realign_FIX:realign_FIX+(sort_FIX_duration/window_step)-1);
    elseif (sort_condition_mode == 4)
        realign_sum = List_analyze_RS_FIX(:,realign_FIX:realign_FIX+(sort_FIX_duration/window_step)-1);
    end
end

if (sort_timing_mode == 2)
    if (sort_condition_mode == 1)
        realign_sum = List_analyze_LS_CUE(:,realign_CUE:realign_CUE+(sort_CUE_duration/window_step)-1);
    elseif (sort_condition_mode == 2)
        realign_sum = List_analyze_LB_CUE(:,realign_CUE:realign_CUE+(sort_CUE_duration/window_step)-1);
    elseif (sort_condition_mode == 3)
        realign_sum = List_analyze_RB_CUE(:,realign_CUE:realign_CUE+(sort_CUE_duration/window_step)-1);
    elseif (sort_condition_mode == 4)
        realign_sum = List_analyze_RS_CUE(:,realign_CUE:realign_CUE+(sort_CUE_duration/window_step)-1);
    end
end

if (sort_timing_mode == 3)
    if (sort_condition_mode == 1)
        realign_sum = List_analyze_LS_SAC(:,realign_SAC:realign_SAC+(sort_SAC_duration/window_step)-1);
    elseif (sort_condition_mode == 2)
        realign_sum = List_analyze_LB_SAC(:,realign_SAC:realign_SAC+(sort_SAC_duration/window_step)-1);
    elseif (sort_condition_mode == 3)
        realign_sum = List_analyze_RB_SAC(:,realign_SAC:realign_SAC+(sort_SAC_duration/window_step)-1);
    elseif (sort_condition_mode == 4)
        realign_sum = List_analyze_RS_SAC(:,realign_SAC:realign_SAC+(sort_SAC_duration/window_step)-1);
    end
end

if (sort_timing_mode == 4)
    if (sort_condition_mode == 1)
        realign_sum = List_analyze_LS_REW(:,realign_REW:realign_REW+(sort_REW_duration/window_step)-1);
    elseif (sort_condition_mode == 2)
        realign_sum = List_analyze_LB_REW(:,realign_REW:realign_REW+(sort_REW_duration/window_step)-1);
    elseif (sort_condition_mode == 3)
        realign_sum = List_analyze_RB_REW(:,realign_REW:realign_REW+(sort_REW_duration/window_step)-1);
    elseif (sort_condition_mode == 4)
        realign_sum = List_analyze_RS_REW(:,realign_REW:realign_REW+(sort_REW_duration/window_step)-1);
    end
end

realign_sum = realign_sum';
realign_sum = (sum(realign_sum)');
[B, sort_order]=sort(realign_sum,'descend');

List_analyze_LS_FIX = List_analyze_LS_FIX(sort_order,:);
List_analyze_LS_CUE = List_analyze_LS_CUE(sort_order,:);
List_analyze_LS_SAC = List_analyze_LS_SAC(sort_order,:);
List_analyze_LS_REW = List_analyze_LS_REW(sort_order,:);
List_analyze_LB_FIX = List_analyze_LB_FIX(sort_order,:);
List_analyze_LB_CUE = List_analyze_LB_CUE(sort_order,:);
List_analyze_LB_SAC = List_analyze_LB_SAC(sort_order,:);
List_analyze_LB_REW = List_analyze_LB_REW(sort_order,:);
List_analyze_RB_FIX = List_analyze_RB_FIX(sort_order,:);
List_analyze_RB_CUE = List_analyze_RB_CUE(sort_order,:);
List_analyze_RB_SAC = List_analyze_RB_SAC(sort_order,:);
List_analyze_RB_REW = List_analyze_RB_REW(sort_order,:);
List_analyze_RS_FIX = List_analyze_RS_FIX(sort_order,:);
List_analyze_RS_CUE = List_analyze_RS_CUE(sort_order,:);
List_analyze_RS_SAC = List_analyze_RS_SAC(sort_order,:);
List_analyze_RS_REW = List_analyze_RS_REW(sort_order,:);

%%%%%%%%%%%figure%%%%%%%%%%%%
dify=0.23;
difx=0.44;
subplot_1x = 0.1;         subplot_1y = 0.775;
subplot_2x = subplot_1x + subplot_x*analyze_FIX_duration + subplot_spacex;  subplot_2y = subplot_1y;
subplot_3x = subplot_2x + subplot_x*analyze_CUE_duration + subplot_spacex;  subplot_3y = subplot_1y;
subplot_4x = subplot_3x + subplot_x*analyze_SAC_duration + subplot_spacex;  subplot_4y = subplot_1y;
%subplot_5x = subplot_4x + subplot_x*analyze_REW_duration + subplot_spacex;  subplot_5y = subplot_1y;

subplot_5x = subplot_1x;  subplot_5y = subplot_1y-dify;
subplot_6x = subplot_2x;  subplot_6y = subplot_5y;
subplot_7x = subplot_3x;  subplot_7y = subplot_5y;
subplot_8x = subplot_4x;  subplot_8y = subplot_5y;

subplot_9x = subplot_1x;  subplot_9y = subplot_1y-2*dify;
subplot_10x = subplot_2x;  subplot_10y = subplot_9y;
subplot_11x = subplot_3x;  subplot_11y = subplot_9y;
subplot_12x = subplot_4x;  subplot_12y = subplot_9y;

subplot_13x = subplot_1x;  subplot_13y = subplot_1y-3*dify;
subplot_14x = subplot_2x;  subplot_14y = subplot_13y;
subplot_15x = subplot_3x;  subplot_15y = subplot_13y;
subplot_16x = subplot_4x;  subplot_16y = subplot_13y;

figure('position',[50,200,504,712]);
set(gcf,'PaperPosition',[0.25 0.25 8.00 10.50]);

%%%%%%%%%%%Left Small%%%%%%%%%%%%
subplot('Position',[subplot_1x subplot_1y subplot_x*analyze_FIX_duration subplot_y_height]);
List_analyze_LS_FIX(:,maxX_analyze_LS_FIX+1) = 0;
List_analyze_LS_FIX(:,maxX_analyze_LS_FIX+2) = 2;
List_analyze_LS_FIX(:,maxX_analyze_LS_FIX+3:1000) = [];
imagesc(List_analyze_LS_FIX)
hold on
plot([zeroX_analyze_LS_FIX zeroX_analyze_LS_FIX], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LS_FIX-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LS_FIX maxX_analyze_LS_FIX-0.5],'XTickLabel',{num2str(analyze_FIX_starttime),'FIX',num2str(analyze_FIX_endtime)});
set(gca,'YTick',[0.5 neuron_num+0.5],'YTickLabel',[0 neuron_num]);

subplot('Position',[subplot_2x subplot_2y subplot_x*analyze_CUE_duration subplot_y_height]);
List_analyze_LS_CUE(:,maxX_analyze_LS_CUE+1) = 0;
List_analyze_LS_CUE(:,maxX_analyze_LS_CUE+2) = 2;
List_analyze_LS_CUE(:,maxX_analyze_LS_CUE+3:1000) = [];
imagesc(List_analyze_LS_CUE)
hold on
plot([zeroX_analyze_LS_CUE zeroX_analyze_LS_CUE], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LS_CUE-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LS_CUE maxX_analyze_LS_CUE-0.5],'XTickLabel',{'','CUE',num2str(analyze_CUE_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_3x subplot_3y subplot_x*analyze_SAC_duration subplot_y_height]);
List_analyze_LS_SAC(:,maxX_analyze_LS_SAC+1) = 0;
List_analyze_LS_SAC(:,maxX_analyze_LS_SAC+2) = 2;
List_analyze_LS_SAC(:,maxX_analyze_LS_SAC+3:1000) = [];
imagesc(List_analyze_LS_SAC)
hold on
plot([zeroX_analyze_LS_SAC zeroX_analyze_LS_SAC], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LS_SAC-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LS_SAC maxX_analyze_LS_SAC-0.5],'XTickLabel',{'','SAC',num2str(analyze_SAC_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_4x subplot_4y subplot_x*analyze_REW_duration subplot_y_height]);
List_analyze_LS_REW(:,maxX_analyze_LS_REW+1) = 0;
List_analyze_LS_REW(:,maxX_analyze_LS_REW+2) = 2;
List_analyze_LS_REW(:,maxX_analyze_LS_REW+3:1000) = [];
imagesc(List_analyze_LS_REW)
hold on
plot([zeroX_analyze_LS_REW zeroX_analyze_LS_REW], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LS_REW-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LS_REW maxX_analyze_LS_REW-0.5],'XTickLabel',{'','REW',num2str(analyze_REW_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

%%%%%%%%%%%Left Big%%%%%%%%%%%%
subplot('Position',[subplot_5x subplot_5y subplot_x*analyze_FIX_duration subplot_y_height]);
List_analyze_LB_FIX(:,maxX_analyze_LB_FIX+1) = 0;
List_analyze_LB_FIX(:,maxX_analyze_LB_FIX+2) = 2;
List_analyze_LB_FIX(:,maxX_analyze_LB_FIX+3:1000) = [];
imagesc(List_analyze_LB_FIX)
hold on
plot([zeroX_analyze_LB_FIX zeroX_analyze_LB_FIX], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LB_FIX-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LB_FIX maxX_analyze_LB_FIX-0.5],'XTickLabel',{num2str(analyze_FIX_starttime),'FIX',num2str(analyze_FIX_endtime)});
set(gca,'YTick',[0.5 neuron_num+0.5],'YTickLabel',[0 neuron_num]);

subplot('Position',[subplot_6x subplot_6y subplot_x*analyze_CUE_duration subplot_y_height]);
List_analyze_LB_CUE(:,maxX_analyze_LB_CUE+1) = 0;
List_analyze_LB_CUE(:,maxX_analyze_LB_CUE+2) = 2;
List_analyze_LB_CUE(:,maxX_analyze_LB_CUE+3:1000) = [];
imagesc(List_analyze_LB_CUE)
hold on
plot([zeroX_analyze_LB_CUE zeroX_analyze_LB_CUE], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LB_CUE-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LB_CUE maxX_analyze_LB_CUE-0.5],'XTickLabel',{'','CUE',num2str(analyze_CUE_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_7x subplot_7y subplot_x*analyze_SAC_duration subplot_y_height]);
List_analyze_LB_SAC(:,maxX_analyze_LB_SAC+1) = 0;
List_analyze_LB_SAC(:,maxX_analyze_LB_SAC+2) = 2;
List_analyze_LB_SAC(:,maxX_analyze_LB_SAC+3:1000) = [];
imagesc(List_analyze_LB_SAC)
hold on
plot([zeroX_analyze_LB_SAC zeroX_analyze_LB_SAC], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LB_SAC-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LB_SAC maxX_analyze_LB_SAC-0.5],'XTickLabel',{'','SAC',num2str(analyze_SAC_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_8x subplot_8y subplot_x*analyze_REW_duration subplot_y_height]);
List_analyze_LB_REW(:,maxX_analyze_LB_REW+1) = 0;
List_analyze_LB_REW(:,maxX_analyze_LB_REW+2) = 2;
List_analyze_LB_REW(:,maxX_analyze_LB_REW+3:1000) = [];
imagesc(List_analyze_LB_REW)
hold on
plot([zeroX_analyze_LB_REW zeroX_analyze_LB_REW], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_LB_REW-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_LB_REW maxX_analyze_LB_REW-0.5],'XTickLabel',{'','REW',num2str(analyze_REW_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

%%%%%%%%%%%Right Big%%%%%%%%%%%%
subplot('Position',[subplot_9x subplot_9y subplot_x*analyze_FIX_duration subplot_y_height]);
List_analyze_RB_FIX(:,maxX_analyze_RB_FIX+1) = 0;
List_analyze_RB_FIX(:,maxX_analyze_RB_FIX+2) = 2;
List_analyze_RB_FIX(:,maxX_analyze_RB_FIX+3:1000) = [];
imagesc(List_analyze_RB_FIX)
hold on
plot([zeroX_analyze_RB_FIX zeroX_analyze_RB_FIX], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RB_FIX-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RB_FIX maxX_analyze_RB_FIX-0.5],'XTickLabel',{num2str(analyze_FIX_starttime),'FIX',num2str(analyze_FIX_endtime)});
set(gca,'YTick',[0.5 neuron_num+0.5],'YTickLabel',[0 neuron_num]);

subplot('Position',[subplot_10x subplot_10y subplot_x*analyze_CUE_duration subplot_y_height]);
List_analyze_RB_CUE(:,maxX_analyze_RB_CUE+1) = 0;
List_analyze_RB_CUE(:,maxX_analyze_RB_CUE+2) = 2;
List_analyze_RB_CUE(:,maxX_analyze_RB_CUE+3:1000) = [];
imagesc(List_analyze_RB_CUE)
hold on
plot([zeroX_analyze_RB_CUE zeroX_analyze_RB_CUE], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RB_CUE-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RB_CUE maxX_analyze_RB_CUE-0.5],'XTickLabel',{'','CUE',num2str(analyze_CUE_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_11x subplot_11y subplot_x*analyze_SAC_duration subplot_y_height]);
List_analyze_RB_SAC(:,maxX_analyze_RB_SAC+1) = 0;
List_analyze_RB_SAC(:,maxX_analyze_RB_SAC+2) = 2;
List_analyze_RB_SAC(:,maxX_analyze_RB_SAC+3:1000) = [];
imagesc(List_analyze_RB_SAC)
hold on
plot([zeroX_analyze_RB_SAC zeroX_analyze_RB_SAC], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RB_SAC-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RB_SAC maxX_analyze_RB_SAC-0.5],'XTickLabel',{'','SAC',num2str(analyze_SAC_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_12x subplot_12y subplot_x*analyze_REW_duration subplot_y_height]);
List_analyze_RB_REW(:,maxX_analyze_RB_REW+1) = 0;
List_analyze_RB_REW(:,maxX_analyze_RB_REW+2) = 2;
List_analyze_RB_REW(:,maxX_analyze_RB_REW+3:1000) = [];
imagesc(List_analyze_RB_REW)
hold on
plot([zeroX_analyze_RB_REW zeroX_analyze_RB_REW], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RB_REW-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RB_REW maxX_analyze_RB_REW-0.5],'XTickLabel',{'','REW',num2str(analyze_REW_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

%%%%%%%%%%%Right Small%%%%%%%%%%%%
subplot('Position',[subplot_13x subplot_13y subplot_x*analyze_FIX_duration subplot_y_height]);
List_analyze_RS_FIX(:,maxX_analyze_RS_FIX+1) = 0;
List_analyze_RS_FIX(:,maxX_analyze_RS_FIX+2) = 2;
List_analyze_RS_FIX(:,maxX_analyze_RS_FIX+3:1000) = [];
imagesc(List_analyze_RS_FIX)
hold on
plot([zeroX_analyze_RS_FIX zeroX_analyze_RS_FIX], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RS_FIX-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RS_FIX maxX_analyze_RS_FIX-0.5],'XTickLabel',{num2str(analyze_FIX_starttime),'FIX',num2str(analyze_FIX_endtime)});
set(gca,'YTick',[0.5 neuron_num+0.5],'YTickLabel',[0 neuron_num]);

subplot('Position',[subplot_14x subplot_14y subplot_x*analyze_CUE_duration subplot_y_height]);
List_analyze_RS_CUE(:,maxX_analyze_RS_CUE+1) = 0;
List_analyze_RS_CUE(:,maxX_analyze_RS_CUE+2) = 2;
List_analyze_RS_CUE(:,maxX_analyze_RS_CUE+3:1000) = [];
imagesc(List_analyze_RS_CUE)
hold on
plot([zeroX_analyze_RS_CUE zeroX_analyze_RS_CUE], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RS_CUE-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RS_CUE maxX_analyze_RS_CUE-0.5],'XTickLabel',{'','CUE',num2str(analyze_CUE_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_15x subplot_15y subplot_x*analyze_SAC_duration subplot_y_height]);
List_analyze_RS_SAC(:,maxX_analyze_RS_SAC+1) = 0;
List_analyze_RS_SAC(:,maxX_analyze_RS_SAC+2) = 2;
List_analyze_RS_SAC(:,maxX_analyze_RS_SAC+3:1000) = [];
imagesc(List_analyze_RS_SAC)
hold on
plot([zeroX_analyze_RS_SAC zeroX_analyze_RS_SAC], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RS_SAC-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RS_SAC maxX_analyze_RS_SAC-0.5],'XTickLabel',{'','SAC',num2str(analyze_SAC_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);

subplot('Position',[subplot_16x subplot_16y subplot_x*analyze_REW_duration subplot_y_height]);
List_analyze_RS_REW(:,maxX_analyze_RS_REW+1) = 0;
List_analyze_RS_REW(:,maxX_analyze_RS_REW+2) = 2;
List_analyze_RS_REW(:,maxX_analyze_RS_REW+3:1000) = [];
imagesc(List_analyze_RS_REW)
hold on
plot([zeroX_analyze_RS_REW zeroX_analyze_RS_REW], [0 10000],'-','color',[0 0 0],'LineWidth',1)
axis([0.5 maxX_analyze_RS_REW-0.5 0.5 neuron_num+0.5])
set(gca,'XTick',[0.5 zeroX_analyze_RS_REW maxX_analyze_RS_REW-0.5],'XTickLabel',{'','REW',num2str(analyze_REW_endtime)});
set(gca,'YTick',[],'YTickLabel',[]);












