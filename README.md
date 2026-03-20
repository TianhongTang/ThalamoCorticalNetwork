[TOC] 

# ThalamoCorticalNetwork


## Abstract


## Code Structure


## Data Structure


### General File Format

Each data file, except `metadata.mat`, is stored in a `[datatype]_##.mat` file containing two structs named `meta` and `data`. The file name consists of the type of the data followed by the identification index constructed from its meta data.

| Field | Description |
|------|------|
| meta | Identification and description data. Also storaged in `metadata.mat` |
| data | Large data matrices and variables. |

**Data shape notation**: All 1-dimentional data with shape `(data_len)` is stored as row vectors with shape `(1, data_len)`.


### `metadata.mat` - Meta data file
Each field is a struct array for a data type. The field name is the type of the data. The struct array contains all meta data in the `.metadata` field and the full file path of the data file.

| Field | Description |
|------|------|
| raster | |
|  |  |

<!--------------------------------->
### `raster_##.mat` - Rasterized spikes

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `Pre` or Post-injection `Post`. |
| state | string | State identifier. |
| area | string | Dataset brain area identifier. Single area `ACC/VLPFC/Thalamus` or merged area type `Full/Cortex` |
| align | string | Alignment type. For unaligned data, this field is `'None'`. |
| session_idx | int | Session index. |
| date | string | Recording date. Format: `MMDDYYYY`. |
| N | int | Neuron number. |
| dt | double | Time bin size in seconds. |
| trial_num | int | Total trial number. |
| trial_len | int | Sorted `data.trial_len`. |
| file_name | string | File name. |
| max_len | int | Maximum trial length. |
| min_len | int | Minimum trial length. |
| total_len | int | Total trial length. |
| resting_dur_threshold (optional) | double | Minimum duration threshold for resting state trials. |
| align_kernel (optional) | string | Kernel used to align data. |
| align_kernel_len (optional) | string | Kernel length. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| rasters | cell(int) | `{trial_num}(N, trial_len)` | Raster matrices. Each cell contains a binary matrix indicating whether the neuron fired a spike within a time bin. |
| spikes | cell(int) | `{trial_num}(spike_num)` | Raw spike timings relative to trial start in ms. |
| trial_len | int | `(trial_num)` | Time bin number for each trial. |
| cell_id | string | `(N)` | Neuron identifier. |
| cell_area | string | `(N)` | Brain area for each neuron. |
| channel | double | `(N)` | Neuron location relative to the probe. 64 channels per probe. |
| cuetype | cell(double) | `{trial_num}(cue_num)` | Cue type identifiers for each trial. Each cell contains all cue identifiers in the trial. For non-task sessions, cell contains an empty array. |
| firing_rates | cell(double) | `{trial_num}(N)` | Mean firing rate for each trial. Each cell contains firing rates of neurons in Hz.|


<!--------------------------------->
### `border_##.mat` - Brain area border

#### meta
| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| area | string | Dataset brain area identifier. Merged area type `Full/Cortex` |
| align | string | Alignment type. For aligned data only. |
| session_idx | int | Session index. |
| N | int | Number of neurons. |
| area_num | int | Number of brain areas. |
| file_name | string | File name. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| borders | int | `(area_num+1)` | First `area_num` elements: Neuron indices for each area's first neuron. Last element: `N+1`. |


<!--------------------------------->
### `spike_##.mat` - Array dataset spike data
Extracted spike data from the array dataset.

#### meta
| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal identifier. |
| date | string | Recording date. Format: `MMDDYYYY`. |
| session_name_raw | string | Raw data file identifier. |
|  |  |  |
|  |  |  |
|  |  |  |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
|  |  |  |  |
|  |  |  |  |
|  |  |  |  |
|  |  |  |  |


<!--------------------------------->
### `sortidx_##.mat` - Neuron sorting indices

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| state | string | State identifier. |
| area | string | Dataset brain area identifier. Merged area type `Full/Cortex` |
| align | string | Alignment type. For aligned data only. |
| session_idx | int | Session index. |
| criterion | string | Sorting criterion. |
| file_name | string | File name. |
| N | int | Number of neurons. |
| kernel (optional) | string | GLM model kernel name. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| sort_idx | int | (N) | Neuron indice in the original data. |


<!--------------------------------->
### `shuffled_##.mat` - Shuffled raster file

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| state | string | State identifier. |
| area | string | Dataset brain area identifier. Merged area type `Full/Cortex` |
| align | string | Alignment type. For aligned data only. |
| session_idx | int | Session index. |
| shuffle_mode | string | Shuffle mode. `None`: No shuffling. `Within trial`: Shuffle time bins within each trial. `Across trial`: Shuffle time bins in all trials. |
| shuffle_idx | int | Shuffle index. 0 is default for no shuffle. |
| shuffle_seed | int | Random seed used in shuffling. |
| file_name | string | File name. |
| N | int | Number of neurons. |
| dt | double | Time bin size in seconds. |
| trial_num | int | Number of trials. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| rasters | cell(int) | `{trial_num}(N, trial_len)` | Shuffled raster matrices. Each cell contains a binary matrix indicating whether the neuron fired a spike within a time bin. |
| trial_len | int | `(trial_num)` | Time bin number for each trial. |
| firing_rates | cell(double) | `{trial_num}(N)` | Mean firing rate for each trial. Each cell contains firing rates of neurons in Hz.|

<!--------------------------------->
### `crossval_##.mat` - Splitted raster groups for cross validation
Rasters are splitted and assigned to folds for cross validation. For task sessions, each trial is assigned to a fold. For resting sessions, rasters are splitted into equal-length segments. Each segment is assigned to a splitted fold. Assignment details are stored in the field `data.assignments`.

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| state | string | State identifier. |
| area | string | Dataset brain area identifier. Merged area type `Full` or `Cortex` |
| align | string | Alignment type. For aligned data only. |
| session_idx | int | Session index. |
| shuffle_idx | int | Shuffle index. 0 is default for no shuffle. |
| file_name | string | File name. |
| N | int | Number of neurons. |
| fold_num | int | Number of splitted folds. |
| assignment_num | int | Number of assignments. |
| fold_total_len | int | Total length of each fold. |
| fold_trial_lens | cell(int) | Same as `data.fold_trial_lens`. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| fold_rasters | cell(cell)(int) | {fold_num}{segment_num}(N, fold_trial_len) | Raster segments. |
| fold_trial_lens | cell(int) | {fold_num}(segment_num) | Duration of raster segments. |
| assignments | cell(struct) | {assignment_num} | Details of fold assignment. Fields: <br>`trial_index`: Index for the source trial. <br>`segment_index`: Index of the segment in the source trial. <br>`fold`: Fold index assigned to. <br> `length`: Segment length. |
|  |  |  |  |


<!--------------------------------->
### `GLMdata_##.mat` - GLM training data
Convolved and concatenated rasters.

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| state | string | State identifier. |
| area | string | Dataset brain area identifier. Merged area type `Full` or `Cortex` |
| align | string | Alignment type. For aligned data only. |
| session_idx | int | Session index. |
| shuffle_idx | int | Shuffle index. 0 is default for no shuffle. |
| kernel_name | string | Kernels used in GLM. |
| file_name | string | File name. |
| N | int | Number of neurons. |
| fold_num | int | Number of cross validation folds. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| folds | cell(struct) | {fold_num} | Concatenated and convolved rasters of each cross validation fold. Fields: <br> `B`: `int`, Concatenated raster length. <br> `raster`: `int(N, B)`, Raster of the fold. <br> `predjs_conn`: `double(N, B, n_conn_kernel)`, raster convolved by connection kernels. <br> `predjs_PS`: `double(N, B, n_PS_kernel)`, raster convolved by post-spike kernels. |
| kernel | struct | - | Kernel used in GLM. See `Kernel_##.mat`. |

### `GLM_##.mat` - Trained GLM parameters

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal name. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| state | string | State identifier. |
| area | string | Dataset brain area identifier. Merged area type `Full` or `Cortex` |
| align | string | Alignment type. For aligned data only. |
| session_idx | int | Session index. |
| shuffle_idx | int | Shuffle index. 0 is default for no shuffle. |
| kernel_name | string | Kernels used in GLM. |
| reg_name | struct | Regularization used in training. |
| epoch | int | Epochs trained. |
| fold_num | int | Number of cross validation folds. |
| fold_idx | int | Left out fold index. 0 is default for no cross validation models. |
| file_name | string | File name. |
| N | int | Number of neurons. |
| N_filtered | int | Number of neurons after filtering out no-spike and ####(other criterions?)#### neurons. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| model_par | double | `(N, 1 + n_PS_kernel + N*n_conn_kernel)` | All parameters of the GLM. Consists of three parts: <br> First column: `h_i`, baseline activity of neuron i. <br> Next `n_PS_kernel` columns: `P_ik`, k-th post-spike kernel weight of neuron i. <br> Next `N*n_conn_kernel` columns: `J_ijk`, k-th connection kernel weight from neuron j to neuron i.|
| model_err | struct | - | Standard error of `model_par`. Computed by the square root of inverse hessian. Fields: <br> `minuslogL`: Ignore regularization. <br> `total`: Include regularization. |
| train_loss | struct | - | Training loss of the last epoch. Fields: `minuslogL`, `rag`, `total`. |
| test_loss | struct | - | Testing loss of the last epoch. Fields: `minuslogL`, `rag`, `total`. |
| filter | int | (N) | Binary filter for valid neurons. |
| reg | struct | - | Regularization. Fields: <br> `name`: Regularization name. <br> `l1`: L1 regularization. <br> `l2`: L2 regularization. |
| kernel | struct | - | Kernels used in GLM. See `kernel_##.mat`. |


### `kernel_##.mat` - GLM kernel

#### meta

| Field | Type | Description |
|------|------|-------------|
| kernel_name | string | Kernel name. |
| n_conn_kernel | int | Number of connection kernels. |
| n_PS_kernel | int | Number of post-spike kernels. |
| kernel_len | int | Kernel length. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| conn_kernels | cell(double) | {n_conn_kernel}(kernel_len) | Weight of connection kernels. |
| PS_kernels | cell(double) | {n_conn_kernel}(kernel_len) | Weight of post-spike kernels. |


### `_##.mat` - file

#### meta

| Field | Type | Description |
|------|------|-------------|
|  |  |  |
|  |  |  |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
|  |  |  |  |
|  |  |  |  |
|  |  |  |  |
|  |  |  |  |


### `_##.mat` - file

#### meta

| Field | Type | Description |
|------|------|-------------|
|  |  |  |
|  |  |  |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
|  |  |  |  |
|  |  |  |  |
|  |  |  |  |
|  |  |  |  |