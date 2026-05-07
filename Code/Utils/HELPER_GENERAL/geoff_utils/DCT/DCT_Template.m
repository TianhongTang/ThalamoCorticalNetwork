function job = DCT_Template(arg1, arg2, etc, sched)
% This is a template for writing functions that take advantage of the
% Distributed Computing Toolbox. This function does not itself do anything,
% and should not be called. Instead, copy this file and modify it to your
% needs.
% First, change the input arguments above, as desired. I suggest keeping
% 'sched' as an input argument.
% Next, go through the code below. I have tried to point out places where
% you should modify code as follows:
%-MODIFY------------------------------------------------------------------%
%(some code)
%-------------------------------------------------------------------------%
% When you are ready to call your code, execute the following at the MATLAB
% prompt:
% >> sched = findResource('scheduler','type','local');
% >> job = DCT_Template(arg1, arg2, etc, sched);
% (changing the name of the DCT_Template function and the argument list as
% necessary, of course)
% To check the status of your job, use the command:
% >> findTask(job)
% Finally, to retrieve your output arguments, use the following:
% >> outputs = getAllOutputArguments(job);
% >> out1 = outputs(:,1);
% >> out2 = outputs(:,2);
% etc. In this case, out1 and out2 will be cell arrays representing the
% first and second output arguments from the job (respectively), and the
% k'th element of each is the result of the k'th iteration of the loop.
% 
% This template makes use of a few helper functions intended to make it
% easier to pass shared input arguments to all the workers.
% addLargeSharedArgument and addGlobal modify the job's jobData property.
% Therefore you should either use these functions instead of modifying the
% jobData property directly, or be cautious in setting the jobData property
% yourself. In particular, they require that the jobData property is a
% structure with the fields "bigArgs" and "globals". You can add additional
% fields if you wish, but removing these fields will generate an error in
% the workers.

%-MODIFY------------------------------------------------------------------%
% Insert your preprocessing code here.
%-------------------------------------------------------------------------%
job = createJob(sched);
%-MODIFY------------------------------------------------------------------%
% Change the following to the number of output arguments from each loop
% iteration:
numLoopOutputs = 1;
% This will be the horizontal size of the outputs from
% getAllOutputArguments when your job has finished running.
%-------------------------------------------------------------------------%

%-MODIFY------------------------------------------------------------------%
% Replace "bigArg1 bigArg2 etc" with any large input arguments that all
% tasks will share:
addLargeSharedArgument bigArg1 bigArg2 etc;
%-------------------------------------------------------------------------%

%-MODIFY------------------------------------------------------------------%
% If you would like to make global variables available to the workers,
% uncomment the following line. If you would like to make only specific
% globals available, replace "-all" with a list of the globals you would
% like to use:
%addGlobal -all

% Note that this will make the global available in the workers' global
% workspaces. However, these workspaces are NOT shared. This means that if
% your loop iteration code modifies the value of a global variable, other
% workers will NOT see that change, nor will that change be made available
% in your MATLAB environment. The above is useful only for making globals
% available for reading access.
%-------------------------------------------------------------------------%

%-MODIFY------------------------------------------------------------------%
% Replace numLoops with the expression for the number of iterations:
for k=1:numLoops
%-------------------------------------------------------------------------%
%-MODIFY------------------------------------------------------------------%
% Replace loopArgs(k) with a cell array of the iteration-specific input
% arguments (e.g., "{x(k), y, z(:,k)}")
    createTask(job, @loopIteration, numLoopOutputs, loopArgs(k));
%-------------------------------------------------------------------------%
end
submit(job);
end % DCT_Template

%-MODIFY------------------------------------------------------------------%
% Set your input and output arguments as desired:
function [out1, out2, etc] = loopIteration(in1, in2, etc)
%-------------------------------------------------------------------------%
getSharedJobData;
%-MODIFY------------------------------------------------------------------%
% Your loop code goes here. Any large arguments specified in the
% addLargeSharedArgument command (above) are assigned by getSharedJobData,
% as are any globals specifed in the addGlobal command. Therefore you can
% use them as you would if this code really were in a loop, rather than in
% a separate function.
%-------------------------------------------------------------------------%
end % loopIteration


%% Helper functions (No need to modify):
function addLargeSharedArgument(varargin)
if ~all(cellfun(@ischar,varargin))
    error('DCT_Template:addLargeSharedArgument:badArgument', ...
        'Variable names must be strings.');
end
job = evalin('caller','job');
jobData = job.jobData;
if isempty(jobData)
    jobData = struct('bigArgs',struct,'globals',struct);
end
for k=1:numel(varargin)
    jobData.bigArgs.(varargin{k}) = evalin('caller',varargin{k});
end
job.jobData = jobData;
end % addLargeSharedArgument

function addGlobal(varargin)
if ~all(cellfun(@ischar,varargin))
    error('DCT_Template:addGlobal:badArgument', ...
        'Variable names must be strings.');
end
if nargin==1 && strcmp(varargin{1},'-all')
    globals = evalin('caller','who(''global'')');
else
    globals = varargin;
end
job = evalin('caller','job');
jobData = job.jobData;
if isempty(jobData)
    jobData = struct('bigArgs',struct,'globals',struct);
end
for k=1:length(numel(globals))
    jobData.globals.(globals{k}) = eval(globals{k});
end
job.jobData = jobData;
end % addGlobal

function getSharedJobData
job = getCurrentJob;
jobData = job.jobData;
if isempty(jobData)
    return
end
bigArgs = jobData.bigArgs;
bigArgNames = fieldnames(bigArgs);
globals = jobData.globals;
globalNames = fieldnames(globals);
for k=1:numel(bigArgNames)
    assignin('caller',bigArgNames{k},bigArgs.(bigArgNames{k}));
end
for k=1:numel(globalNames)
    evalstr = sprintf('%s = globals.(globalNames{k}); global %s;', ...
        globalNames{k}, globalNames{k});
    eval(evalstr);
    evalin('caller',['global ' globalNames{k}]);
end
end % getSharedJobData