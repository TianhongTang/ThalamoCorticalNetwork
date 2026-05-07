Col20 = 20*ones(size(clusterAnalysis.dim1, 1),1);

input1 = [Col20, repmat(clusterAnalysis.dim1, 1, 100), -Col20] ;
input2 = [Col20, repmat(clusterAnalysis.dim2, 1, 100), -Col20] ;
input3 = [Col20, repmat(clusterAnalysis.dim3, 1, 100), -Col20] ;
input4 = [Col20, repmat(clusterAnalysis.dim4, 1, 100), -Col20] ;
input5 = [Col20, repmat(clusterAnalysis.dim6, 1, 100), -Col20] ;
input6 = [Col20, repmat(clusterAnalysis.dim7, 1, 100), -Col20] ;
input7 = [Col20, repmat(clusterAnalysis.dim8, 1, 100), -Col20] ;

inSDFs = { input1, input2,input3,input4,input5,input6,input7};
inTimes = { [0:101], [0:101], [0:101], [0:101], [0:101], [0:101], [0:101]};       
myEpocs ={ [1:50], [1:50], [1:50], [1:50], [1:50], [1:50], [1:50]};
myEpocInds =[ 1, 2, 3, 4, 5, 6, 7];
minN = 10;

[sortIDs, idxDist, raw, respSumStruct, respLink] = consensusCluster(inSDFs, inTimes,'-e',myEpocs,'-ei',myEpocInds,'-r',0,'-mn',minN);