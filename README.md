# Snakemake workflow: clinical DICOMs to BIDS

## Description

Snakemake workflow to convert a clinical dicom directory into BIDS structure.

## Requirements

* dcm2niix (v1.0.20200427)
* python requirements (defined in `workflow/envs/mapping.yaml`):
    * pydicom>=1.0.2
    * setuptools>=39.2.0
    * extractCMRRPhysio>=0.1.1
    * dcmstack>=0.7.0
    * pandas>=0.24.2
    * heudiconv>=0.8.0

## Input directory structure

The input directory with dicoms should be setup as follows:
```sh

sourcedata/
└── <subject>/
    ├── <sequence>/<dicom_files.dcm>
    ├── <sequence>/<dicom_files.dcm>
    └── <sequence>/<dicom_files.dcm>
```
* `<subject>` is the identifier for the subject in the form `sub-001`, `sub-002` etc.
* `<sequence>` is the directory for a specific imaging sequence and can be given any name

## Operation dates

One main feature of this clinical pipeline is that the final output is stored in three session directories for each subject:
* **ses-presurg:** for imaging data aquired prior to the subjects surgery date
* **ses-perisurg:** for imaging data aquired on the same day of the subjects surgery (generally MRI/CT with sterotactic frame)
* **ses-postsurg:** for imaging data aquired the immediate day after the subjects surgery onwards

To divide the imaging data into the three sessions, the operation date is required foreach subject. This date is attempted to be gleaned automatically assuming some type of intraoperative imaging is performed and the **SeriesDescription** and/or **StudyDescription** DICOM header tag includes an identifer that it was aquired intra-operatively. 

If the operation date cannot be determined automatically you may provide a path to an `or_dates.tsv` file that has the operation dates for all subjects. 

Here is an example of what this file should look like:

| subject  | or_date   |
|:---------|:-----------|
| sub-001<img width="100"/>  | 2014_09_28<img width="100"/> |
| sub-002  | 2018_03_26 |
| sub-003  | n/a |

If the operation date for a subject is unknown or the subject did not undergo surgery, you can put `n/a` under **or_date** for that subject. In this case, the output BIDS directory will contain only the `presurg` session for that subject (with all the imaging data stored there):

```
output/bids/sub-P001/
  └── ses-presurg/anat/...
```

## Setting up

### Step 1: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

```sh
conda create -c bioconda -c conda-forge -n snakemake snakemake
```

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 2: Clone a copy of this repository

To clone a local copy of this repository to your system:

```sh
git clone https://github.com/greydongilmore/sampleClinicalWorkflow.git
```

### Step 3: Modify the configuration file

#### file paths

Edit the `config/config.yaml` file to include the paths to the following:

<center>

|Variable   |Description        |
|:----------|:------------------|
| `dicom_dir`<img width="200"/> | full path to where the input dicom directory is stored   |
| `out_dir`   | full path to where the pipeline should output the data   |
| `heuristic`  | heudiconv template file to sort/name dicoms according to BIDS standard |
| `or_dates_file` **[optional]** | path to the `or_dates.tsv` file with subject surgery dates |

</center>

#### session determination

How one medical center performs aquisition of imaging data around the surgery date will most likely differ from another medical center. Further, in many clinical cases a patient will undergo surgery more than once. Within this pipeline you are able to modify how `presurgery`, `perisurgery` and `postsurgery` are defined. Within the `config/config.yaml` you will notice the following settings under **session_calc**: 
<center>

|Variable   |Description        |
|:----------|:------------------|
| `periop` [default: 0]<img width="280"/> | number of days ± around surgery day that will deemed `perisurg` |
| `dur_multi_surg` [default: -30]| in the event of multi patient surgery, the max num days before the follow-up surgery imaging data will be deemed `presurg`<sup>1</sup> |
| `override_periop` [default: True]  | in the event an imaging study is deemed `perisurg` but it contains an **electrode** flag it will be moved to `postsurg` <sup>2</sup> |

</center>

* <sup>1</sup> the default value is quantified as 30 days before the subsequent surgery any imaging aquisitions will be sorted into the `presurg` session
* <sup>2</sup> this may occur if the surgery occurs in the morning and a post-op imaging study is performed to localize implanted electrodes later the same day


### Step 4: DICOM sort rules

Depending on your dicom dataset you may need to create a sorting rule for your data. There are two points in the pipeline where the imaging headers are parsed:
* within **dicom2tar** to create the Tarball archives
* within **tar2bids** to create the BIDS filename conventions

#### dicom2tar sort rule

The default sort rule can be found in [workflow/scripts/dicom2tar/sort_rules.py](workflow/scripts/dicom2tar/sort_rules.py) and is the function **sort_rule_clinical**. This sort rule seperates the DICOM files based on image type (MRI/CT/fluoro) as well as scan date. Scan sessions occuring on different days are stored in different Tarball archives, even if they are the same image type. Depending on the clinical scanner used at your center and the information stored within the DICOM header tags you may need to add/substract from this sort rule.

#### tar2bids sort rule

The sort rule for **heudiconv** is called a **heuristic**. They provide a detailed description of the [heuristic and how to write your own](https://heudiconv.readthedocs.io/en/latest/heuristics.html). The default heuristic file within this pipeline can be found in [workflow/scripts/heudiconv/clinical_imaging.py](workflow/scripts/heudiconv/clinical_imaging.py). 

### Step 5: Run the pipeline

All the following commands should be run in the root of the directory.

1. Active the conda environment by running:

```sh
conda activate snakemake
```

2. Prior to running you can do a dry run:

```sh
snakemake --use-conda -n
```

3. To locally run the pipeline, run the following:

```sh
snakemake --use-conda -j $N
```

where `$N` specifices the number of cores to use.


## Description of the pipeline

### Repository structure

The repository has the following scheme:
```
├── README.md
├── workflow
│   ├── rules
|   │   └── dicom2bids.smk
│   ├── envs
|   │   └── mapping.yaml
│   ├── scripts
|   |   ├── dicom2tar               # sorts and stores the dicoms into Tarballs
|   │   |   ├── clinical_helpers.py
|   |   |   ├── DicomSorter.py
|   |   |   ├── main.py
|   │   |   └── sort_rules.py
|   |   ├── heudiconv               # heuristic file for clinical imaging
|   │   |   ├── clinical_imaging.py 
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
|Overview   | Sorts and stores the dicom files into Tarball files |
|Input      | MRI/CT dicoms |
|output     | MRI/CT Tarballs|

</center>

This workflow is modified from the [dicom2tar](https://github.com/khanlab/dicom2tar) master branch (version date:16/07/2020). The code has been modfied to fit the current pipeline.  

The first pass will provide an intermediate output, Tarball files with the associated scan date in the name:
```
output/
  └── tars/
        ├── P185_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

The final output will sort the Tarball files by date and assign sequential numbers to the files:
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
|Overview   | Converts the Tarball archives into BIDS compliant format |
|Input      | MRI/CT dicom Tarball archives |
|output     | MRI/CT nifti files stored in BIDS layout|

</center>

output will be:
```
output/bids/sub-P185/
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
|Overview   | Restructures the BIDS session naming to more meaningful session types |
|Input      | BIDS directory with default session numbering |
|output     | BIDS directory with session naming according to surgery date|

</center>

output will be:
```
output/bids_final/sub-P185/
  ├── ses-perisurg/anat/...
  ├── ses-postsurg/anat/...
  └── ses-presurg/anat/...
```
