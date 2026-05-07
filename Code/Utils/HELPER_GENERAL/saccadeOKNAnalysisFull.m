function sacc2D = saccadeOKNAnalysisFull(em_h,em_v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two one-dimensional saccade analyses, followed by a 2-D stage
% where the final info is computed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsacc  = saccadeOKNAnalysis1D(em_h);
vsacc  = saccadeOKNAnalysis1D(em_v);
if ~isempty(hsacc) & ~isempty(vsacc)
   sacc2D = saccadeOKNAnalysis2D(hsacc, vsacc);
else
   sacc2D = [];
end


