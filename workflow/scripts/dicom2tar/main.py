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
import sys
import os
import logging
import argparse

import sort_rules
import DicomSorter

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

    args = Namespace(dicom_dir=snakemake.input[0], output_dir=snakemake.output[0], clinical_scans=True, get_or_dates = snakemake.params[0])

    logger = logging.getLogger(__name__)

    if not os.path.exists(args.dicom_dir):
        logger.error("{} not exist!".format(args.dicom_dir))
        return False
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
                
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
