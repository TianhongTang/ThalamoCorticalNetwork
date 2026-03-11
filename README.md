# ThalamoCorticalNetwork

## Abstract


## Data File Structure


### General File Format

Each data file, except `metadata.mat`, is stored as a `[datatype]_##.mat` file containing two structs named `meta` and `data`. The file name consists of the type of the data followed by the identification index constructed from its meta data.

| Field | Description |
|------|------|
| meta | Meta data of the data. Also storaged in `metadata.mat`|
| data | Large variables containing the data. |

**Data shape notation**: All 1-dimentional data with shape `(data_len)` is stored as row vectors with shape `(1, data_len)`.


### `metadata.mat` - Meta data file: 
Each field is a struct array for a data type. The field name is the type of the data. The struct array contains all meta data in the `.meta` field and the full file path of the data file.

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
| session_name | string | Session full name as identifier. |
| session_idx | int | Session index. |
| date | string | Recording date. Format: `MMDDYYYY`. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection `pre` or Post-injection `post`. |
| N | int | Neuron number. |
| dt | double | Time bin size in seconds. |
| area | string | Dataset brain area identifier. Single area `ACC/VLPFC/Thalamus` or merged area type `Full/Cortex` |
| align | string | Alignment type. For aligned data only. |
| trial_num | int | Total trial number. |
| state | string | State identifier. |
| file_name | string | File name. |
|  |  |  |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| rasters | cell(int) | `{trial_num}(N, trial_len)` | Raster matrices. Each cell contains a binary matrix indicating whether the neuron fired a spike within a time bin. |
| trial_len | int | `(trial_num)` | Time bin number for each trial. |
| cell_id | string | `(N)` | Neuron identifier. |
| cell_area | string | `(N)` | Brain area for each neuron. |
| channel | double | `(N)` | Neuron location relative to the probe. 64 channels per probe. |
| cuetype | cell(double) | `{trial_num}(cue_num)` | Cue type identifiers for each trial. Each cell contains all cue identifiers in the trial. For non-task sessions, cell contains an empty array. |
| firing_rates | cell(double) | `{trial_num}(N)` | Mean firing rate for each trial. Each cell contains firing rates of neurons in Hz.|
|  |  |  |  |
|  |  |  |  |


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
| session_name | string | Session full name as identifier. |
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
| kernel | string | GLM model kernel name. |
| criterion | string | Sorting criterion. |
| file_name | string | File name. |
| N | int | Number of neurons. |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| sort_idx | int | (N) | Neuron indice in the original data. |


<!--------------------------------->
### file: `_##.mat`

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


### file: `_##.mat`

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