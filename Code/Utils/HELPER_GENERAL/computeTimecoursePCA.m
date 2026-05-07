function [pcadata, e_vector, e_value] = computeTimecoursePCA(data,newdims)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Principal Component Analysis (PCA)
%
% ver 1.00 16-12-1999  Yusuke MURAYAMA
% modified 17-02-2000  David  LEOPOLD
% modified 22-06-2000  for fMRI data DAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


w_ts = length(data(1,:));
mean_vec = zeros(1, w_ts);
cvar_mat = zeros( w_ts,  w_ts);
tmp_timecourse = zeros(1,  w_ts);
tmp_mat = zeros( w_ts,  w_ts);
[n_data, n] = size(data);

wb = waitbar(0, 'Computing Principal Components');
for i=1:n_data
  tmp_timecourse = data(i,:);
  mean_vec = mean_vec + tmp_timecourse;
  tmp_mat = (tmp_timecourse') * tmp_timecourse;
  cvar_mat = cvar_mat + tmp_mat;
  waitbar(i/n_data,wb);
end
close(wb);

mean_vec = mean_vec / n_data;
cvar_mat = cvar_mat / n_data - (mean_vec') * mean_vec;

% Now, get eigenvectors and eigenvalues
[e_vector, e_value] = eig(cvar_mat);

for i=1:w_ts
   e_value(1,i) = e_value(i,i);
end
e_value = e_value(1,:);
[e_value, e_value_idx] = sort(e_value);
e_value = fliplr(e_value);
e_value_idx = fliplr(e_value_idx);
tmp_mat = e_vector(:, e_value_idx);
e_vector = tmp_mat;

% plots variance and cumulative variance.
norm_var = e_value / sum(e_value);
cumu_var = cumsum(norm_var);

% compute values on new axis
if nargin < 2
  newdims = size(data,2);
end
new_axs = e_vector(:, 1:newdims);
pcadata = data * new_axs;




