
[pvalcorr1,corr1] = permutation_pair_test_fast(a , b ,20000,'rankcorr')
[pvalcorr2,corr2] = permutation_pair_test_fast(c , d ,20000,'rankcorr')


Np = length(a);
Nb = 20000;
bootcorr1 = zeros(Nb,1);
bootcorr2 = zeros(Nb,1);
for B = 1:Nb
    bootsample = randi(Np,1,Np);
    [g, bootcorr1(B)] = permutation_pair_test_fast(a(bootsample),b(bootsample) ,1 ,'rankcorr');
    [g,bootcorr2(B)] = permutation_pair_test_fast(c(bootsample),d(bootsample) ,1 ,'rankcorr');
    clear g
end

alpha = 0.05; % probability coverage - 0.05 for 95% CI
hi = floor((1-alpha/2)*Nb+.5);
lo = floor((alpha/2)*Nb+.5);

% for each correlation
boot1sort = sort(bootcorr1);
boot2sort = sort(bootcorr2);
boot1ci = [boot1sort(lo) boot1sort(hi)];
boot2ci = [boot2sort(lo) boot2sort(hi)];

% for the difference between correlations
bootdiff = bootcorr1 - bootcorr2;
bootdiffsort = sort(bootdiff);
diffci = [bootdiffsort(lo) bootdiffsort(hi)];


clc;

bootval=bootdiff;
nullval=0;

p=min(1,2*min(mean(bootval<=nullval),mean(bootval>=nullval)))


pvalue = mean(bootdiffsort<0);
pvalue = 2*min(pvalue,1-pvalue)



% Bieniek, M.M., Frei, L.S. & Rousselet, G.A. (2013) Early ERPs to faces: aging, luminance, and individual differences. Frontiers in psychology, 4, 268.
% 
% Nieuwenhuis, S., Forstmann, B.U. & Wagenmakers, E.J. (2011) Erroneous analyses of interactions in neuroscience: a problem of significance. Nat Neurosci, 14, 1105-1107.
% 
% Pernet, C.R., Wilcox, R. & Rousselet, G.A. (2012) Robust correlation analyses: false positive and power validation using a new open source matlab toolbox. Front Psychol, 3, 606.
% 
% Wilcox, R.R. (2009) Comparing Pearson Correlations: Dealing with Heteroscedasticity and Nonnormality. Communications in Statistics-Simulation and Computation, 38, 2220-2234.
% 
% Wilcox, R.R. (2012) Introduction to robust estimation and hypothesis testing. Academic Press, San Diego, CA.
% 
% Wilcox, R.R. (2016) Comparing dependent robust correlations. Brit J Math Stat Psy, 69, 215-224.
% 
