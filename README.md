[![DOI](https://zenodo.org/badge/280235977.svg)](https://zenodo.org/badge/latestdoi/280235977)

# Snakemake workflow: clinical DICOMs to BIDS

## Description

Snakemake workflow to convert a clinical dicom directory into BIDS structure.

## Requirements

* dcm2niix (v1.0.20200427)
* python requirements (defined in `workflow/envs/mapping.yaml`):
    * dcmstack>=0.7.0
    * dicognito>=0.11.0
    * heudiconv>=0.8.0
    * pandas>=0.24.2
    * pydicom>=1.0.2
    * setuptools>=39.2.0
    * snakemake>=5.23.0

## Input directory structure

The input directory with dicoms should be setup as follows:

```sh
data/
├── dicoms/
│   └── <subject>/
│        ├── <sequence>/<dicom_files.dcm>
│        ├── <sequence>/<dicom_files.dcm>
│        └── <sequence>/<dicom_files.dcm>
└── output/
```

* `data` the main directory, which stores input subject directories
* `dicoms` directory that stores the source DICOM files
   * `<subject>` is the identifier for the subject in the form `sub-001`, `sub-002` etc.
   * `<sequence>` is the directory for a specific imaging sequence and can be given any name
* `output` directory will store all outputs from the pipeline

## Clinical event dates

One main feature of this clinical pipeline is that the final output can be stored around a clinical event. For instance if the clinical event was an operation, which included preop, periop, and postop imaging, then the output directory would contain three session folders:
* **ses-pre:** for imaging data acquired prior to the clinical event
* **ses-peri:** for imaging data acquired on the same day as the clinical event (i.e. intraoperative imaging)
* **ses-post:** for imaging data acquired after the clinical event

To divide the imaging data based on a clinical event, the event date is required for each subject. The date should be defined in a tab seperated text file named `clinical_events.tsv`, which has the clinical event date defined for all subjects. 

Here is an example of what this file should look like:

| subject  | event_date   |
|:---------|:-----------|
| sub-001<img width="100"/>  | 2014_09_28<img width="100"/> |
| sub-002  | 2018_03_26 |
| sub-003  | n/a |

If the event date for a subject is unknown or the subject, place `n/a` under **event_date** for that subject. In this case, the output BIDS directory will contain only the `pre` session for that subject (with all the imaging data stored there):

```
output/bids/sub-P001/
  └── ses-pre/anat/...
```

## Configuration

### Modify the configuration file

#### File paths

If you are running this pipeline locally, edit the `config/config.yaml` file to include the paths to the following:

<center>

|Variable   |Description        |
|:----------|:------------------|
| `dicom_dir`<img width="200"/> | full path to where the input dicom directory is stored   |
| `out_dir`   | full path to where the pipeline should output the data   |
| `clinical_event_file` **[optional]** | path to the `clinical_events.tsv` file with subject clinical event dates |
| `heuristic`  | heudiconv template file to sort/name dicoms according to BIDS standard |
| `dcm_config`  | used to modify the input parameters for dcm2niix |
| `anonymize` | whether to anonymize the dicom files prior to storing in tar archives (default = True)|

</center>

#### Clinical event session determination

How one medical centre acquires imaging for/around a clinical event will differ from another centre. Additionally, patients may experience multiple clinical events, which increases the complexity of ensuring the imaging studies are stored in the appropriate session directory for each clinical event. Within the dicom2bids pipeline, the user can modify how `pre`, `peri` and `post` sessions are defined. Within the `config/config.yaml`, the following settings should be used to customize:

<center>

|Variable   |Description        |
|:----------|:------------------|
| `peri` [default: 0]<img width="280"/> | number of days ± around clinical event that will deemed `peri` |
| `dur_multi_event` [default: -30]| in the event of multiple clinical events (i.e. repeat surgery), the max num days, before the subsequent clinical event, that imaging studies will be deemed `pre`<sup>1</sup> |
| `override_peri` [default: True]  | in the event an imaging study is deemed `peri`, but the heuristic file detects a post clinical event flag, it will be moved to `post` <sup>2</sup> |

</center>

* <sup>1</sup> the default value is quantified as 30 days before the subsequent clinical event, any imaging acquisitions during that period will be stored in the `pre` session
* <sup>2</sup> this may occur if the clinical event and post event imaging occur on the same day.

### DICOM sort rules

DICOM files will first be sorted based on modality and acquisition date. The default DICOM sorting heuristic should work on any dataset but may need some adjustment. The tar archives will then be converted to nifit format and sorted based on image sequence name. This heurisitc will most likely need to be modified based.

#### dicom2tar sort rule

The default sort rule can be found in [workflow/scripts/dicom2tar/sort_rules.py](workflow/scripts/dicom2tar/sort_rules.py) and is the function **sort_rule_clinical**. This sort rule separates the DICOM files based on image type (MRI/CT/fluoro) as well as acquisition date. Initially, image acquisitions occurring on different days are stored in separate Tar archives. The tar archives are then given sequentially numbering based on acquisition date.

#### tar2bids sort rule

The sort rule for **HeuDiConv** is called a **heuristic**. Refer to this detailed description of the [heuristic and how to write your own](https://heudiconv.readthedocs.io/en/latest/heuristics.html). The default heuristic file within the dicom2bids workflow can be found in [workflow/scripts/heudiconv/clinical_imaging.py](workflow/scripts/heudiconv/clinical_imaging.py). 

## Running locally

1. Install Snakemake using [Python](https://www.python.org/downloads/):

    ```sh
    python -m pip install snakemake
    ```

    For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

2. Clone a copy of this repository to your system:

    ```sh
    git clone https://github.com/greydongilmore/clinical_dicom2bids_smk.git
    ```

3. Install the Python dependencies by opening a terminal, changing to the project directory and running:

    ```sh
    python -m pip install -r requirements.txt
    ```

### Local run

All the following commands should be run in the root of the project directory.

1. Prior to running you can do a dry run:

   ```sh
   snakemake -n
   ```

2. To locally run the pipeline, run the following:

   ```sh
   snakemake -j $N
   ```

   where `$N` specifies the number of cores to use.


## Description of the pipeline

### Repository structure

The repository has the following scheme:
```
├── README.md
├── workflow
│   ├── rules
│   │   └── dicom2bids.smk
│   ├── envs
│   │   └── mapping.yaml
│   ├── scripts
│   |   ├── dicom2tar               # sorts and stores the dicoms into Tarballs
|   │   |   ├── clinical_helpers.py
|   |   |   ├── DicomSorter.py
|   |   |   ├── main.py
|   │   |   └── sort_rules.py
|   |   ├── heudiconv               # heuristic file for clinical imaging
|   │   |   └── clinical_imaging.py 
|   |   └── post_tar2bids           # refactors heudiconv output into final BIDS structure
|   │       └── clean_sessions.py
|   └── Snakefile
├── config
│   └── config.yaml
└── data                            # contains test input dicoms
```

### Rule 01: dicom2tar

<center>

| Variable  | Description                              |
|:----------|:-----------------------------------------|
|Overview   | Sorts and stores the dicom files into Tar archives |
|Input      | MRI/CT dicoms |
|output     | MRI/CT Tar archives|

</center>

The dicom2tar pipeline is modified from the [dicom2tar](https://github.com/khanlab/dicom2tar) master branch (version date:16/07/2020). The code has been modified to fit the dicom2bids workflow.  

First, the dicom2tar pipeline will create tar archives with the acquisition date embedded in the archive name:

```
output/
  └── tars/
        ├── P185_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

The final dicom2tar output will sort the tar archives based on date and assign sequential sessions numbers, which are embedded in the archive name:

```
output/
  └── tars/
        ├── P185_001_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_002_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_003_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_004_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_005_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

### Rule 02: tar2bids

<center>

| Variable  | Description                              |
|:----------|:-----------------------------------------|
|Overview   | Converts the Tar archives into BIDS compliant format |
|Input      | MRI/CT dicom Tar archives |
|output     | MRI/CT nifti files stored in BIDS format|

</center>

The output from the rule tar2bids will be:

```
output/
  └── bids_tmp/
        └── sub-P185/
              ├── ses-001/anat/...
              ├── ses-002/anat/...
              ├── ses-003/anat/...
              ├── ses-004/anat/...
              └── ses-005/anat/...
```

### Rule 03: cleanSessions

<center>

| Variable  | Description                              |
|:----------|:-----------------------------------------|
|Overview   | Reorganizes the HeuDiConv BIDS output into the final session(s) |
|Input      | BIDS directory with default session numbering |
|output     | BIDS directory with clinically relevant session naming |

</center>

The output from the rule cleanSessions will be:

```
output/
  └── bids/
        └── sub-P185/
              ├── ses-peri/anat/...
              ├── ses-post/anat/...
              └── ses-pre/anat/...
```
