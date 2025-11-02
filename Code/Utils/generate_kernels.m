%% 

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
% gaussian function
gaus = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

%% kernel generate 
% Each kernel set can contain multiple 'conn_kernel's (connection) 
% and 'PS_kernel's (post-spike). All kernels must have same length,
% add zeros to align all kernels.

% kernel file: "conn_kernels", "PS_kernels",
% "n_conn_kernel", "n_PS_kernel", "kernel_len"

% exponential decay kernel
tau1=10; % synaptic integration time constant (in ms)
T1=3*tau1; % cutoff on the sum for the neuron total input
tau_all1=0:T1;
tt_start=1+T1;

kernel = exp(-tau_all1./tau1);
kernel(1) = 0; % remove simultanuous corr
kernel = kernel/sum(kernel); % Normalize

conn_kernels = {kernel};
PS_kernels = {};
kernel_len=T1+1;

n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save_folder = fullfile(root, 'Data', 'Working', 'kernel');
check_path(save_folder);
save(fullfile(save_folder, 'kernel_expDecay10.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

% zero-delay kernel
kernel = ones(1);
conn_kernels = {kernel};
PS_kernels = {};
kernel_len = 1;

n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_zeroDelay.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

% gaussian kernel

% multi-kernel group
% expMulti200
T = 200;
t=0:T;
kernel_len = T+1;
tau1=10;
tau2=50;
n_conn_kernel = 2;
n_PS_kernel = 2;

k1 = exp(-t/tau1);
k2 = t.*exp(-t/tau2);
k1 = k1/sum(k1);
k2 = k2/sum(k2);
conn_kernels = {k1, k2};
k1(1)=0;
k1 = k1/sum(k1);
PS_kernels = {k1, k2};

save(fullfile(save_folder, 'kernel_expMulti200.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");



% steps 50
conn_kernels = cell(1, 10);
PS_kernels = cell(1, 10);
for i=1:10
    kernel = zeros(50);
    kernel((i*5-4):(i*5))=1;
    kernel = kernel/sum(kernel); % Normalize

    conn_kernels{i} = kernel;
    PS_kernels{i} = kernel;
end
kernel_len = 50;

n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_steps50.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

% steps 25
conn_kernels = cell(1, 5);
PS_kernels = cell(1, 5);
for i=1:5
    kernel = zeros(25);
    kernel((i*5-4):(i*5))=1;
    kernel = kernel/sum(kernel); % Normalize

    conn_kernels{i} = kernel;
    PS_kernels{i} = kernel;
end
kernel_len = 25;

n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_steps25.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

%% ----------------expGauss60
T = 60;
t = 0:T;
kernel_len = T+1;
tau = 5;
center = 30;
sigma = 5;
k1 = exp(-t/tau);
k2 = gaus(t, center, sigma, 1, 0);
k1 = k1/sum(k1);
k2 = k2/sum(k2);
conn_kernels = {k1, k2};
k1(1)=0;
k1 = k1/sum(k1);
PS_kernels = {k1, k2};

% plot kernels
figure("Visible", "off");
hold on;
plot(t, k1, 'r');
plot(t, k2, 'b');
hold off;
xlabel('Time (ms)');
ylabel('Amplitude');
title('expGauss60');
legend('Exponential', 'Gaussian');
saveas(gcf, fullfile(save_folder, 'kernel_expGauss60.png'));

% save kernel
n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_expGauss60.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

%% ----------------exp5Gauss5C## groups
T = 60;
t = 0:T;
tau = 5;
sigma = 5;
kernel_len = T+1;
for center = [20, 30, 40]
    k1 = exp(-t/tau);
    k2 = gaus(t, center, sigma, 1, 0);
    k1 = k1/sum(k1); % how to normalize?
    k2 = k2/sum(k2);
    conn_kernels = {k1, k2};
    k1(1)=0;
    k1 = k1/sum(k1);
    PS_kernels = {k1, k2};

    % plot kernels
    figure("Visible", "off");
    hold on;
    plot(t, k1, 'r');
    plot(t, k2, 'b');
    hold off;
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title(['exp5Gauss5C', num2str(center)]);
    legend('Exponential', 'Gaussian');
    saveas(gcf, fullfile(save_folder, ['kernel_exp5Gauss5C', num2str(center), '.png']));

    % save kernel
    n_conn_kernel=length(conn_kernels);
    n_PS_kernel=length(PS_kernels);
    save(fullfile(save_folder, ['kernel_exp5Gauss5C', num2str(center), '.mat']), "conn_kernels", "PS_kernels", ...
        "n_conn_kernel", "n_PS_kernel", "kernel_len");

%% ----------------Delta: 5ms exp, 10ms gauss centered at 40, 40ms gauss centered at 120
T = 200;
t = 0:T;
kernel_len = T+1;
tau = 5;
sigma1 = 5;
sigma2 = 20;
center1 = 40;
center2 = 120;
k1 = exp(-t/tau);
k2 = gaus(t, center1, sigma1, 1, 0);
k3 = gaus(t, center2, sigma2, 1, 0);
k1 = k1/sum(k1);
k2 = k2/sum(k2);
k3 = k3/sum(k3);
conn_kernels = {k1, k2, k3};
% plot kernels
figure("Visible", "off");
hold on;
plot(t, k1, 'r');
plot(t, k2, 'b');
plot(t, k3, 'g');
hold off;
xlabel('Time (ms)');
ylabel('Amplitude');
title('Delta');
legend('Exponential 5ms', 'Gaussian 40±5ms', 'Gaussian 120±20ms');
saveas(gcf, fullfile(save_folder, 'kernel_Delta.png'));

k1(1)=0;
k1 = k1/sum(k1);
PS_kernels = {k1, k2, k3};

% save kernel
n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_Delta.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

%% ----------------DeltaPure: 5ms exp, 10ms gauss centered at 40, 40ms gauss centered at 120, no self-connection
T = 200;
t = 0:T;
kernel_len = T+1;
tau = 5;
sigma1 = 5;
sigma2 = 20;
center1 = 40;
center2 = 120;
k1 = exp(-t/tau);
k2 = gaus(t, center1, sigma1, 1, 0);
k3 = gaus(t, center2, sigma2, 1, 0);
k1 = k1/sum(k1);
k2 = k2/sum(k2);
k3 = k3/sum(k3);
conn_kernels = {k1, k2, k3};
% plot kernels
figure("Visible", "off");
hold on;
plot(t, k1, 'r');
plot(t, k2, 'b');
plot(t, k3, 'g');
hold off;
xlabel('Time (ms)');
ylabel('Amplitude');
title('Kernel: Delta');
legend('Exponential 5ms', 'Gaussian 40±5ms', 'Gaussian 120±20ms');
saveas(gcf, fullfile(save_folder, 'kernel_DeltaPure.png'));

PS_kernels = {};

% save kernel
n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_DeltaPure.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");
end

%% ----------------DoublePure: 5ms exp, 10ms gauss centered at 40, no self-connection
T = 200;
t = 0:T;
kernel_len = T+1;
tau = 5;
sigma1 = 5;
center1 = 40;
k1 = exp(-t/tau);
k2 = gaus(t, center1, sigma1, 1, 0);
k1 = k1/sum(k1);
k2 = k2/sum(k2);
conn_kernels = {k1, k2};

% plot kernels
figure("Visible", "off");
hold on;
plot(t, k1, 'r');
plot(t, k2, 'b');
hold off;
xlabel('Time (ms)');
ylabel('Amplitude');
title('Kernel: Double');
legend('Exponential 5ms', 'Gaussian 40±5ms');
saveas(gcf, fullfile(save_folder, 'kernel_DoublePure.png'));

PS_kernels = {};

% save kernel
n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_DoublePure.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

%% ----------------SinglePure: 5ms exp only, no self-connection
T = 200;
t = 0:T;
kernel_len = T+1;
tau = 5;
k1 = exp(-t/tau);
k1 = k1/sum(k1);
conn_kernels = {k1};

% plot kernels
figure("Visible", "off");
hold on;
plot(t, k1, 'r');
hold off;
xlabel('Time (ms)');
ylabel('Amplitude');
title('Kernel: Single');
legend('Exponential 5ms');
saveas(gcf, fullfile(save_folder, 'kernel_SinglePure.png'));

PS_kernels = {};

% save kernel
n_conn_kernel=length(conn_kernels);
n_PS_kernel=length(PS_kernels);
save(fullfile(save_folder, 'kernel_SinglePure.mat'), "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");