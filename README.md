# Snakemake workflow: clinical DICOMs to BIDS

## Description

Snakemake workflow to convert a clinical dicom directory into BIDS structure.

## Requirements

* dcm2niix (v1.0.20200427)
* python requirements in `requirements.txt`:
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

## Setting up

### Step 1: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 2: Clone a copy of this repository

To clone a local copy of this repository to your system:

```sh
git clone https://github.com/greydongilmore/sampleClinicalWorkflow.git
```

#### Modify the configuration file

Edit the `config/config.yaml` file to include the paths to the following:

<center>

|Variable   |Description        |
|:----------|:------------------|
| **dicom_dir** | full path to where the input dicom directory is stored   |
| **out_dir**   | full path to where the pipeline should output the data   |
| **heuristic**  | heudiconv template file to sort/name dicoms according to BIDS standard |

</center>

### Step 3: Run the pipeline

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
|   │   ├── dicom2tar.smk
|   |   ├── tar2bids.smk
|   │   └── cleanSessions.smk
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
|   |   ├── post_tar2bids           # refactors heudiconv output into final BIDS structure
|   │   |   ├── clean_sessions.py
|   │   └── process_complete.py     # writes a log file when process completes
|   └── Snakefile
├── config
│   └── config.yaml
└── data                            # contains test input dicoms
```

### dicom2tar

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

## tar2bids

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

## post tar2bids cleanup

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


