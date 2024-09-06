#!/bin/bash
/home/greydon/.local/bin/heudiconv --files /home/greydon/Documents/data/DBS/sourcedata/tars/P266/P266_001_2022_11_29_20221129_061626_MN220021224.36909769_MR.tar -o /home/greydon/Documents/data/DBS/bids_tmp -f resources/heuristics/clinical_imaging_old.py -c dcm2niix --dcmconfig resources/dcm_config.json -g all
