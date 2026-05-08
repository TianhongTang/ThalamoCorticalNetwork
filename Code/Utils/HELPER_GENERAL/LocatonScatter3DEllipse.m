
%%%%% 1 - we need to fit a real sphere to confidence intervals i plot and find which neurons are IN or OUT 
%%%%% 2 - we need to make the confidence intervals real 95% / right now they are scaled see comment below
%%%%% 3 - we need to use this and use Slicer-Matlab interface that I am
%%%%% working on to find regions within CL

%%%%% note to self: we need to include Dan Moran in authorship or
%%%%% acknowledgements later

clear; close all;
load('/Users/Ilya/Desktop/Plot_3D_all_session_Lemmy_neuron.mat')
THR=0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figuren;

NoveltyNEURONSnotsig=Neuron3D(find(Neuron3D(:,4)>=THR),1:3); %%%%%%%% later we can make these for excited only or something. we will see how the data go.
xall=NoveltyNEURONSnotsig(:,1);
yall=NoveltyNEURONSnotsig(:,2);
zall=NoveltyNEURONSnotsig(:,3);
scatter3(xall,yall,zall,20,'k','filled'); hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoveltyNEURONS=Neuron3D(find(Neuron3D(:,4)<THR),1:3); 

A(:,1)=NoveltyNEURONS(:,1);
A(:,2)=NoveltyNEURONS(:,2);
A(:,3)=NoveltyNEURONS(:,3);

scatter3(A(:,1),A(:,2),A(:,3),100,'r','filled');

[C, SCORE, Variance, TSQUARED, EXPLAINED, MU] = pca(A,'Centered','off');

%% Find the range of each column of SCORE - column 1 is
%% PC1, column 2 is PC2 etc.

R = range(SCORE)/2

%%%%%%%%%%%%%%%%%%%%%% Scale each principal component vector by its range (replace with 95% CI)

X = [-C(1,1)*R(1) -C(1,2)*R(2) -C(1,3)*R(3) ; C(1,1)*R(1) C(1,2)*R(2) C(1,3)*R(3)]
    
Y = [-C(2,1)*R(1) -C(2,2)*R(2) -C(2,3)*R(3) ; C(2,1)*R(1) C(2,2)*R(2) C(2,3)*R(3)]

Z = [-C(3,1)*R(1) -C(3,2)*R(2) -C(3,3)*R(3) ; C(3,1)*R(1) C(3,2)*R(2) C(3,3)*R(3)]

MU=mean(A); %%% shift by mean
X2 = X+MU(1);
Y2 = Y+MU(2);
Z2 = Z+MU(3);
 
line(X2,Y2,Z2,'LineWidth',5); 

xlabel('ML')
ylabel('AP')
zlabel('Depth')


