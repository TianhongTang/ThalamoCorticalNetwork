function [tf, nenames] = structcomp(s1,s2)
% STRUCTCOMP   Compare structures and/or handle graphics objects.
%     [tf, nenames] = structcomp(s1,s2)
%     STRUCTCOMP compares the two structures and/or graphics objects s1 and
%     s2, field by field (or property by property). There are three cases:
%     1. The structures have non-identical fields (order does not matter).
%        In this case, tf = false and nenames = {}.
%     2. The structures have identical fields with identical values. In
%        this case, tf = true and nenames = {}.
%     3. The structures have identical fields, but some values are
%        different. In this case, tf = false, and nenames is a cell array
%        containing the names of the non-equal fields.

if ~isstruct(s1) && ishandle(s1)
    s1 = get(s1);
end
if ~isstruct(s2) && ishandle(s2)
    s2 = get(s2);
end

fnames = fieldnames(s1);
fnames2 = fieldnames(s2);

nenames = {};

if length(fnames)~=length(fnames2) || ~all(strcmp(fnames,fnames2))
    tf = false;
    return
end

for k=1:length(fnames)
    v1 = s1.(fnames{k});
    v2 = s2.(fnames{k});
    if ~strcmp(class(v1),class(v2)) ...
            || length(v1)~=length(v2) ...
            || ~all(v1(:)==v2(:))
        nenames = [nenames fnames{k}]; %#ok<AGROW>
    end
end

tf = isempty(nenames);