
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/greydon/.local/lib/python3.8/site-packages', '/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/workflow/scripts/dicom2tar']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95\xa4\x06\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c1/media/veracrypt7/dicom_dbs_backup/dicoms/sub-121\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x05dicom\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cJ/media/veracrypt6/projects/stealthMRI/working_dir/out/sourcedata/tars/P121\x94a}\x94(h\x0c}\x94\x8c\x03tar\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh)h&ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8cE/media/veracrypt6/projects/stealthMRI/working_dir/clinical_events.tsv\x94\x8c:/media/veracrypt6/projects/stealthMRI/working_dir/out/logs\x94e}\x94(h\x0c}\x94(\x8c\x0fclinical_events\x94K\x00N\x86\x94\x8c\x07log_dir\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh<h8h>h9ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94\x8c\x03121\x94a}\x94(h\x0c}\x94\x8c\x07subject\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94b\x8c\x07subject\x94hMub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhcK\x01heK\x01ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\tdicom_dir\x94\x8c)/media/veracrypt7/dicom_dbs_backup/dicoms\x94\x8c\x07out_dir\x94\x8c5/media/veracrypt6/projects/stealthMRI/working_dir/out\x94\x8c\x13clinical_event_file\x94h8\x8c\x10participants_tsv\x94N\x8c\theuristic\x94\x8c(resources/heuristics/clinical_imaging.py\x94\x8c\ndcm_config\x94\x8c\x19resources/dcm_config.json\x94\x8c\tanonymize\x94\x89\x8c\x0csession_calc\x94}\x94(\x8c\x04peri\x94K\x00\x8c\x0fdur_multi_event\x94J\xe2\xff\xff\xff\x8c\roverride_peri\x94\x88u\x8c\tsub_group\x94\x8c\x07patient\x94\x8c\x0esubject_prefix\x94\x8c\x01P\x94u\x8c\x04rule\x94\x8c\tdicom2tar\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cQ/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/workflow/scripts/dicom2tar\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/workflow/scripts/dicom2tar/main.py';
######## snakemake preamble end #########
#!/usr/bin/env python
'''
sort or tar CFMM' data with DicomSorter


Update: 2019-11-21

Add DICOM tag option arguments to handel DICOM files missing(or empty): 
StudyDesciption, StudyDate, PatientName, StudyID, StudyInstanceUID,
SeriesNumber,InstanceNumber, and SOPInstanceUID. These tags are 
necessary to CFMM sort rule


Author: YingLi Lu
Email:  yinglilu@gmail.com
Date:   2018-05-22

Note:
    Tested on windows 10/ubuntu 16.04, python 2.7.14

'''
import os
import logging
from dicognito.anonymizer import Anonymizer
import sort_rules
import DicomSorter
import pydicom
import clinical_helpers as ch

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s -%(message)s')

class Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def main():
    '''
    use DicomSorter sort or tar CFMM's dicom data

    input:
        dicom_dir: folder contains dicom files(and/or compressed files:.zip/.tgz/.tar.gz/.tar.bz2)
        output_dir: output sorted or tar files to this folder
    '''

    args = Namespace(dicom_dir=snakemake.input.dicom, output_dir=snakemake.output.tar, clinical_scans=True,
        clinical_events = snakemake.params.clinical_events, log_dir=snakemake.params.log_dir, config=snakemake.config)

    logger = logging.getLogger(__name__)

    if not os.path.exists(args.dicom_dir):
        logger.error("{} not exist!".format(args.dicom_dir))
        return False
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if not os.path.exists(args.log_dir):
        os.makedirs(args.log_dir)

    if snakemake.config['anonymize']:
        print('De-identifying imaging data for {}\n'.format(os.path.split(args.dicom_dir)[-1]))
        anonymizer = Anonymizer()
        for root, folders, files in os.walk(os.path.join(args.dicom_dir)):
            for file in files:
                fullpath = os.path.abspath(os.path.join(root,file))
                old_time={}
                with pydicom.dcmread(fullpath, force=True) as dataset:
                    old_time={'StudyDate': dataset.StudyDate if 'StudyDate' in dataset else None,
                              'SeriesDate': dataset.SeriesDate if 'SeriesDate' in dataset else None,
                              'StudyTime': dataset.StudyTime if 'StudyTime' in dataset else None,
                              'SeriesTime': dataset.SeriesTime if 'SeriesTime' in dataset else None,
                              'ContentDate': dataset.ContentDate if 'ContentDate' in dataset else None,
                              'ContentTime': dataset.ContentTime if 'ContentTime' in dataset else None,
                              'AcquisitionDate': dataset.AcquisitionDate if 'AcquisitionDate' in dataset else None
                              }
                    anonymizer.anonymize(dataset)
                    for key,val in old_time.items():
                        if val is not None:
                            dataset[key].value=val
                    dataset.save_as(fullpath)
                    
            
    ######
    # CFMM sort rule
    ######
    try:
        if not args.clinical_scans:

            with DicomSorter.DicomSorter(sort_rules.sort_rule_CFMM, args) as d:
                # #######
                # # sort
                # #######
                # sorted_dirs = d.sort()
                # # logging
                # for item in sorted_dirs:
                #     logger.info("sorted directory created: {}".format(item))

                #######
                # tar
                #######
                # pi/project/study_date/patient/studyID_and_hash_studyInstanceUID
                tar_full_filenames = d.tar(5)
                # logging
                for item in tar_full_filenames:
                    logger.info("tar file created: {}".format(item))

            # ######
            # # demo sort rule
            # ######
            # with DicomSorter.DicomSorter(args.dicom_dir, sort_rules.sort_rule_demo, output_dir) as d:
            #     # sort
            #     sorted_dirs = d.sort()
            #     #logging
            #     for item in sorted_dirs:
            #         logger.info("sorted directory created: {}".format(item))

            #     # tar
            #     # patient_name/study_date/series_number/new_filename.dcm
            #     tar_full_filenames = d.tar(2)
            #     # logging
            #     for item in tar_full_filenames:
            #         logger.info("tar file created: {}".format(item))

        else:
            ######
            # Clinical sort rule
            ######
            logger.info("These are clinical scans.")
            
            with DicomSorter.DicomSorter(sort_rules.sort_rule_clinical, args) as d:
                # if os.path.exists(os.path.join(args.output_dir, 'errorInfo.tsv')):
                #     os.remove(os.path.join(args.output_dir, 'errorInfo.tsv'))
                # if os.path.exists(os.path.join(args.output_dir, 'or_dates.tsv')):
                #     os.remove(os.path.join(args.output_dir, 'or_dates.tsv'))
                # tar
                # study_date/patient/modality/series_number/new_filename.dcm
                tar_full_filenames = d.tar(4)

                # logging
                for item in tar_full_filenames:
                    logger.info("tar file created: {}".format(item))

            ch.tarSession(args)

    except Exception as e:
        logger.exception(e)

if __name__ == "__main__":

    main()
