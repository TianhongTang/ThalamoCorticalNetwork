function ob = set(ob,data)
% ARCHDAT/SET  Set the value of an archdat reference
% Usage:
%   ad = set(ad, data)
% The SET method is provided to modify the archdat data, but it should be
% used with caution. In general, archdat objects are best used to point to
% static data that is no longer being changed (e.g., the end result of a
% major computation). If the SET method must be used, there are two
% "gotchas" to be aware of. First, the archdat object automatically writes
% the new data to disk when SET is called; multiple calls to SET will
% therefore be likely to be very slow. Second, if multiple archdat objects
% point to the same archived data, changing the value of the data with SET
% should change the value reported by GET on all objects. However, if any
% of the objects are currently ready, they will continue to report the old
% value.

classname = classConst('classname');

if ~isscalar(ob)
    error([classname ':nonScalarSet'], ...
        'Can only set data on scalar archobj''s.');
end

if isnull(ob)
    error([classname ':setNull'], 'Cannot set data on a null archobj.');
end

checkArchive(ob.file,ob.var);

if ob.inmem
    ob.data = data;
end

s.(ob.var) = data;
save(ob.file, '-struct', 's', '-append');