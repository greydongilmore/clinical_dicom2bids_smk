
# Sample Workflow

structure of repository:
```
data/           # contains input dicoms and working dirs for tars and bids
dicom2tar/      # parts of the dicom2tar master branch (with modified clinical sort rule) 
heudiconv/      # heuristic file for clinical imaging
post_tar2bids/  # python script to refactor heudiconv output into final structure
```

input dicoms have been anonymized, example dicom directory is:
```
data/dicoms
  └── sub-180/
        └── <dicom_dirs>
```

## requirements

* dcm2niix (v1.0.20200427)
* python requirements in `requirements.txt`: heudiconv (v0.8.0)
    * pydicom>=1.0.2
    * setuptools>=39.2.0
    * extractCMRRPhysio>=0.1.1
    * dcmstack>=0.7.0
    * pandas>=0.24.2
    * heudiconv>=0.8.0

## dicom2tar

forked dicom2tar master branch version (date:16/07/2020), packed here with the most recent `sort_rule.py`. There are changes to clinical rule not currently on the dicom2tar master branch, once finalized I will open a PR. 

```
cd */sampleClinicalWorkflow

python /dicom2tar/main.py data/dicoms data/tars --clinical_scans
```

the intermediate output, with no session identifiers, will be:
```
data/
  └── tars/
        ├── P185_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

within the dicom2tar `main.py`, the function `tarSessions` is imported from `clinical_helpers.py`. The Tarballs are sorted by date and sequential session numbering is added to the Tarball scheme. The output will be:
```
data/
  └── tars/
        ├── P185_001_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_002_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_003_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_004_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_005_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

## tar2bids

the session identifiers are parsed using the following `infotoids` (defined within the `clinical_imaging.py` heuristic):
```python
def infotoids(seqinfos, outdir):
	subject = get_unique(seqinfos, 'example_dcm_file').split('_')[0]
	session = get_unique(seqinfos, 'example_dcm_file').split('_')[1]
	
	ids = {
    'locator': '',
    'session': session,
    'subject': subject,
	}
				
	return ids
``` 
running with heudiconv (v 0.8.0) directly for now:

```
heudiconv --files data/tars -o data/bids -f heudiconv/clinical_imaging.py -c dcm2niix -b
```

output will be:
```
data/bids/sub-P185/
  ├── ses-001/anat/...
  ├── ses-002/anat/...
  ├── ses-003/anat/...
  ├── ses-004/anat/...
  └── ses-005/anat/...
```

## post tar2bids cleanup

[WIP] the code is not well documented.

requires two inputs:
* **dicom_dir:** need to extract OR date from source dicoms (uses the intraoperative fluoro to determine OR date)
* **bids_dir:** will copy original bids directory to new directory **bids_final**

```
python post_tar2bids/clean_sessions.py --dicom_dir data/dicoms --bids_dir data/bids 
```

output will be:
```
data/bids_final/sub-P185/
  ├── ses-perisurg/anat/...
  ├── ses-postsurg/anat/...
  └── ses-presurg/anat/...
```

