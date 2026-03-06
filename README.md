# ThalamoCorticalNetwork


### Working data structures

#### raster: Rasterized spike data.
Parameters: 
    'session' [struct]
        Specifies cell group. Cells in the same session should be the same.

        'animal_name'
        'session_name'
        'session_idx'
        'date': Experiment date. Format: 'MMDDYYYY'. 
        'injection'


    'state' [struct]
        Specifies state in a session (pre/post injection, resting/task, alignment...). 

        'prepost': 

Data fields:
    'meta': Meta data.
        'N': Neuron number.
        'trial_num': Trial number.

    'rasters' [cell] (1, trial_num) 
        Cell array of rasters. Each cell is a (N, trial_len) binary matrix. The matrix is rasterized spikes in one trial. Time step = 1ms.
    
    'firing_rates': Firing rate of each trial. Unit: Hz.
        size: (N, trial_num)
    
