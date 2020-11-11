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

    args = Namespace(dicom_dir=snakemake.input.dicom, output_dir=snakemake.output.tar, clinical_scans=True, clinical_events = snakemake.params.clinical_events, log_dir=snakemake.params.log_dir)

    logger = logging.getLogger(__name__)

    if not os.path.exists(args.dicom_dir):
        logger.error("{} not exist!".format(args.dicom_dir))
        return False
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
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
