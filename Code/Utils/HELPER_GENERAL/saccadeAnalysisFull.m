function sacc2D = saccadeAnalysisFull(em_h,em_v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Two one-dimensional saccade analyses, followed by a 2-D stage
% where the final info is computed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hsacc  = saccadeAnalysis1D(em_h);
vsacc  = saccadeAnalysis1D(em_v);
if ~isempty(hsacc) & ~isempty(vsacc)
   sacc2D = saccadeAnalysis2D(hsacc, vsacc);
else
   sacc2D = [];
end


