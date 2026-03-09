# ThalamoCorticalNetwork

## Abstract


## Data File Structure


### General File Format

Each data file, except `metadata.mat`, is stored as a `[datatype]_##.mat` file containing two structs named `meta` and `data`. The file name consists of the type of the data followed by identification index constructed from its meta data.

| Field | Description |
|------|------|
| meta | Meta data of the data. Also storaged in `metadata.mat`|
| data | Large variables containing the data. |


### Meta data file: `metadata.mat`
Each field is a struct array for a data type. The field name is the type of the data. The struct array contains all meta data in the `.meta` field and the full file path of the data file.

| Field | Description |
|------|------|
| raster | |
|  |  |

### Raster file: `raster_##.mat`

#### meta

| Field | Type | Description |
|------|------|-------------|
| animal_name | string | Animal identifier. |
| session_name | string | Session name. |
| session_idx | int | Session index. |
| date | string | Recording date. Format: `MMDDYYYY`. |
| injection | string | Injection type. `Saline`, `Muscimol` or `No injection`. |
| prepost | string | Pre-injection(`pre`) or Post-injection(`post`). |
| N | int | Number of neurons |
| timestep | int | Time bin size |
| bin_size | float | Time bin size in ms |
|  |  |  |

#### data

| Field | Type | Shape | Description |
|------|------|------|-------------|
| rasters | cell | `{1, trial_num}` | Raster matrices. Each cell contains a (N, trial_len) binary matrix. |
| trial_len | int | `(1, trial_num)` | Time bin number of each trial. |
|  |  |  |
|  |  |  |