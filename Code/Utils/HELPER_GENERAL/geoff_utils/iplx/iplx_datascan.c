
#include "mex.h"

#ifndef NULL
#define NULL 0
#endif

#define WS_IN    prhs[0]

#define BL_OUT   plhs[0]
#define DATA_OUT plhs[1]

#define SHORTWORDHDRSIZE 8
#define NUMWFOFF         6
#define NUMWWFOFF        7


static void nextBlock(short **streamPtrPtr)
{ /* Advance the stream pointer to the next block in the data stream: */
    /* The next block is located SHORTWORDHDRSIZE + (the number of words
     * in the block) words after the current block. */
    short numWF,  /* number of waveforms in this block */
          numWWF; /* number of words per waveform */
    numWF  = *((*streamPtrPtr)+NUMWFOFF);
    numWWF = *((*streamPtrPtr)+NUMWWFOFF);
    /* The number of words in the block is numWF*numWWF */
    if (numWF < 0 || numWWF < 0)
        mexErrMsgTxt("Data integrity failure.");
    else
        *streamPtrPtr = (*streamPtrPtr)+SHORTWORDHDRSIZE+numWF*numWWF;
}

static unsigned int countBlocks(short wordstream[], unsigned int streamLen)
{ /* Count the number of data blocks in the stream */
    short *streamPtr = wordstream;
    int numBlocks = 0;
    /* Keep going until the end of the stream is reached: */
    while (streamPtr < (wordstream+streamLen))
    {
        numBlocks++; /* Count the block */
        nextBlock(&streamPtr); /* Go to the next block */
    }
    return numBlocks;
}

static void blockLocs(short wordstream[], unsigned int numBlocks, 
                      unsigned int indices[])
{ /* Return the indices, in short words, of each block: */
    short *streamPtr = wordstream;
    unsigned int k;
    for (k = 0; k < numBlocks; ++k)
    {
        indices[k] = streamPtr-wordstream + 1; /* Get the current index */
        nextBlock(&streamPtr); /* Go to the next block */
    }
}

static void blockData(short wordstream[], unsigned int numBlocks,
                      unsigned int indices[], mxArray *dataCell)
{ /* Put the data from each block into a cell array: */
    short *streamPtr = wordstream;
    mxArray *dataVec; /* Cell array element pointer */
    short *dataPtr; /* Pointer to the data in the cell array element */
    short numWF,  /* number of waveforms in a block */
          numWWF, /* number of words per waveform */
          numWB;  /* number of words in a block */
    unsigned int k;
    for (k = 0; k < numBlocks; ++k)
    {
        indices[k] = streamPtr-wordstream + 1; /* Get the current index */
        numWF  = *(streamPtr+NUMWFOFF);
        numWWF = *(streamPtr+NUMWWFOFF);
        numWB  = numWF*numWWF;
        /* Create a new MATLAB vector */
        dataVec = mxCreateNumericMatrix(numWB, 1, mxINT16_CLASS, mxREAL);
        dataPtr = mxGetData(dataVec);
        if (numWB) /* If there is data in this block, */
        {
            /* copy the data from the stream to the MATLAB vector. */
            memcpy(dataPtr,streamPtr+SHORTWORDHDRSIZE,sizeof(short)*numWB);
        }
        /* Put the new vector into the cell array. */
        mxSetCell(dataCell,k,dataVec); /* This uses 0-based indexing. Grrr. */
        /* Go to the next block. */
        streamPtr = streamPtr+SHORTWORDHDRSIZE+numWB;
    }
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
    short *wordstream; /* Pointer to the data stream */
    unsigned int *indices; /* Locations of data blocks, in short word indices */
    unsigned int streamLen; /* Length of the data stream */
    unsigned int numBlocks; /* Number of data blocks in the stream */
    
    if (nrhs != 1) { 
	mexErrMsgTxt("One input argument required."); 
    } else if (nlhs > 2) {
	mexErrMsgTxt("Too many output arguments."); 
    }
    if (!mxIsInt16(WS_IN))
    {
        mexErrMsgTxt("Wordstream must be int16 type.");
    }
    
    wordstream = mxGetData(WS_IN);
    streamLen = mxGetNumberOfElements(WS_IN);
    
    /* Get the total number of blocks in the data stream so we can
     * initialize our variables: */
    numBlocks = countBlocks(wordstream,streamLen);
    
    /* Initialize block location output argument */
    BL_OUT = mxCreateNumericMatrix(numBlocks,1, mxUINT32_CLASS, mxREAL);
    indices = mxGetData(BL_OUT);
    if (nlhs < 2)
    { /* If there is only one output argument, use the simpler blockLocs */
        blockLocs(wordstream,numBlocks,indices);
    }
    else
    { /* Otherwise, use blockData to get the data in the blocks as well. */
        /* Initialize the data output cell array: */
        DATA_OUT = mxCreateCellMatrix(numBlocks,1);
        blockData(wordstream,numBlocks,indices,DATA_OUT);
    }
    return;
}
