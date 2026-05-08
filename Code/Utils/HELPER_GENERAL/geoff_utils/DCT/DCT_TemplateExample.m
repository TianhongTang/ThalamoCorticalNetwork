function job = DCT_TemplateExample(sig, sigLevels, sched)
% This is an example of the use of the DCT_Template m-file.
% This example takes an input signal "sig", which can be a vector or
% matrix, and a vector of signal-to-noise ratios "sigLevels". For each
% given signal-to-noise ratio, it computes a noisy signal at the given
% ratio. It does this in a distributed fashion, by treating the problem of
% computing the various noisy signals as a job, and each noisy signal as a
% task.
% To see this example in action, do the following on a machine with the DCT
% installed:
% >> sched = findResource('scheduler','Type','local');
% >> sig = peaks(100);
% >> snrs = [1 2 5 10 20 50];
% >> job = DCT_TemplateExample(sig,snrs,sched);
% Then wait until the job is finished. You can make MATLAB do this for you
% with the command:
% >> waitForState(job)
% You will regain control of the prompt after a few seconds, when the job
% is finished. Then, capture your output and display it:
% >> output = getAllOutputArguments(job);
% >> destroy(job);
% >> figure;
% >> for k=1:6, subplot(2,3,k); imagesc(output{k}); end

% Preprocessing:
sigVariance = var(sig(:));   % Signal variance
numLoops = numel(sigLevels); % Number of iterations needed

job = createJob(sched);

numLoopOutputs = 1; % This example only uses one output argument

% Here, sig is our "large" input argument that every iteration of the loop
% sees, so it is better to share it than pass it individually. sigVariance
% is small (a double scalar), but it is still shared, so we can pass it in
% a shared manner if we wish.
addLargeSharedArgument sig sigVariance;

% This example uses no global variables. If it did, I would want to declare
% them global somewhere above this, and then uncomment the following line.
%addGlobal -all

for k=1:numLoops
    % Note that the single input argument here, sigLevels(k), MUST be
    % passed in a cell array.
    createTask(job, @loopIteration, numLoopOutputs, {sigLevels(k)});
end
submit(job);
end % DCT_Template

function out = loopIteration(sigLevel)
% getSharedJobData assigns the variables we passed in
% addLargeSharedArgument, so after the following line, sig and sigVariance
% are assigned:
getSharedJobData;

% Compute noise:
noise = randn(size(sig));
noiseVariance = var(noise(:));

% Add noise at appropriate S/N.
out = sig + noise*sigVariance/(noiseVariance*sigLevel);
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