function val = classConst(name)

switch name
    case 'classname'
        val = 'archdat';
    otherwise
        error([classConst('classname') ':classConst'], ...
            'Undefined class constant "%s".', name);
end