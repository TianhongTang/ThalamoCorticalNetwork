% Attempt to test significance of number of neurons in each group of info
% coding, in animals pooled and in each animal
% 2020-11-20 coded by ESBM


% reward only, aversive only, salience, value, neither rew/ave significant
type_names = {'rewo','aveo','sali','valu','none'};

monk_names = {'bat','sab','all'};

monk_ntotal = [303; 333];
type_n = [26 17 11 0 ; 97 12 4 3];
type_n = [type_n (monk_ntotal - sum(type_n,2))]; % add 'nonsig' cells

type_p_under_null_of_no_coding = [0.05*(1-0.05) 0.05*(1-0.05) (0.05*0.05*0.5) (0.05*0.05*0.5)];

% add 'all' monk
type_n(end+1,:) = sum(type_n,1);
monk_ntotal(end+1,1) = sum(monk_ntotal);
type_p_under_null_of_no_coding(end+1) = (1-0.05)*(1-0.05);


% add more categories
type_names{end+1} = 'rew ';
type_n(:,end+1) = sum(type_n(:,[1 3 4]),2);
type_p_under_null_of_no_coding(end+1) = sum(type_p_under_null_of_no_coding(:,[1 3 4]));

% add more categories
type_names{end+1} = 'ave ';
type_n(:,end+1) = sum(type_n(:,[2 3 4]),2);
type_p_under_null_of_no_coding(end+1) = sum(type_p_under_null_of_no_coding(:,[2 3 4]));

% do right-tailed test b/c there cannot 'truly' be fewer significant neurons than chance
psig = binopvalue(type_n,repmat(monk_ntotal,1,size(type_n,2)),repmat(type_p_under_null_of_no_coding,numel(monk_names),1),'right');


% put results in tables for easy viewing in matlab
table_args = {'VariableNames',type_names,'RowNames',monk_names};

tmp = mat2cell(type_n,3,ones(size(type_n,2),1)); 
number_of_units = table(tmp{:},table_args{:})

tmp = mat2cell(type_n ./ monk_ntotal,3,ones(size(type_n,2),1)); 
prob_units = table(tmp{:},table_args{:})

tmp = mat2cell(type_p_under_null_of_no_coding,1,ones(size(type_n,2),1)); 
prob_units_under_null = table(tmp{:},'VariableNames',type_names)

tmp = mat2cell(psig,3,ones(size(psig,2),1)); 
binomial_test_pvalues = table(tmp{:},table_args{:})

