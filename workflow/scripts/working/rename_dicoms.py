#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import shutil
import string
import argparse
import pydicom as dicom
import traceback
import sys
import re
import datetime


class DicomAggregator( object ):
	"""
	This class is from the brainvisa software
	see: https://github.com/brainvisa/aims-free/blob/79e7a21772bb8fd3a3c6ca6892d09c98049e50f4/pyaims/python/soma/dicom/aggregator.py
	"""
	def __init__( self ):
		self._dicom_sources = []
		self.modalities = []
		self.frameOfReferenceUID = {}

	def add_dicom_sources( self, source ):
		if os.path.isdir( source ):
			for root, dirs, files in os.walk( source ):
				for name in files:
					self._dicom_sources.append( os.path.join( root, name ) )
		
		self._dicom_sources = list( set( self._dicom_sources ) )
	
	def aggregate( self ):
		self._aggregate_sources = {}
		serie_sop_uids = {}
		positions_and_orientations = []
		for f in self._dicom_sources:
			try:
				ds = dicom.read_file( f , stop_before_pixels=True)
			except Exception as e:
				print(e)
				continue
			
			if len(self.modalities) != 0 and \
			   ds.Modality not in self.modalities:
				continue
			
			try:
				serie_uid = ds.SeriesInstanceUID
			except Exception as e:
				print(e)
				continue
			
			if serie_uid not in self._aggregate_sources:
				self._aggregate_sources[ serie_uid ] = []
				serie_sop_uids[ serie_uid ] = []
			
			
			position_and_orientation=[]
			for tag in ("ImagePositionPatient","ImageOrientationPatient","FrameReferenceTime"):
				position_and_orientation.append(ds.get(tag) if tag in ds else None)
			
			sop_uid = ds.SOPInstanceUID
			if not sop_uid in serie_sop_uids[ serie_uid ] and \
			   not tuple(position_and_orientation) in positions_and_orientations:
				serie_sop_uids[ serie_uid ].append( sop_uid )
				self._aggregate_sources[ serie_uid ].append( f )
				
				if all(x is not None for x in position_and_orientation):
					positions_and_orientations.append( position_and_orientation )


def createTree(serie, patientDirName, sessionPatientDirNames):
	"""
	This function creates directories to tidy all the DICOM elements in a serie/acquisition.
	It creates a folder for the patient, with the name of the patient, and a 
	subfolder for the serie/acquisition, named with the modality and the description.

	:param serie: list of all the path DICOM of a serie
	:return: 
	"""

	def clean_path(path):
		return re.sub(r'[^a-zA-Z0-9.-]', '_', '{0}'.format(path))

	if len(serie) == 0:
		return

	first = serie[0]
	
	# Get DICOM information to build file tree
	ds = dicom.read_file(first, stop_before_pixels=True)

	dsAttrs={}
	for tag in ("PatientName","Modality","SeriesDescription","SeriesInstanceUID","StudyDate","StudyTime"):
		dsAttrs[tag]=clean_path(ds.get(tag)) if tag in ds else ''
	
	try:
		
		# Inspect
		inc = 1
		if os.path.exists(patientDirName) and \
		   patientDirName not in sessionPatientDirNames:
			inc = 2
			while os.path.exists(patientDirName + "_" + str(inc)) and \
					patientDirName + "_" + str(inc) not in sessionPatientDirNames:
				inc += 1

		if inc > 1:
			patientDirName += "_" + str(inc)

		if patientDirName not in sessionPatientDirNames:
			os.makedirs(patientDirName)
		sessionPatientDirNames.append(patientDirName)

	except:
		print("Exception in user code:")
		print('-'*60)
		traceback.print_exc(file=sys.stdout)
		print('-'*60)

	try:
		if dsAttrs["StudyDate"] != '':
			dsAttrs["StudyDate"]=datetime.datetime.strptime(dsAttrs["StudyDate"], '%Y%m%d').strftime('%Y-%m-%d')
		if dsAttrs["StudyTime"] != '':
			dsAttrs["StudyDate"]=dsAttrs["StudyDate"]+datetime.datetime.strptime(dsAttrs["StudyTime"].split('.')[0], '%H%M%S').strftime('_%H-%M-%S')
		
		# Create serie/acquisition directory
		serieDir = os.path.join(patientDirName, "_".join([dsAttrs["Modality"], dsAttrs["SeriesDescription"], dsAttrs["StudyDate"]]))
		
		# Inspect if the serie directory exist and change its name in case
		inc2 = 1
		if os.path.exists(serieDir):
			inc2 = 2
			while os.path.exists(serieDir + "_" + str(inc2)):
				inc2 += 1

		if inc2 > 1:
			serieDir += "_" + str(inc2)
		
		os.makedirs(serieDir)

		# Copy DICOM files into serie directory
		for f in serie:
			baseDestFile = os.path.join(serieDir, os.path.basename(f))
			destFile = baseDestFile
			num = 2
			while os.path.exists(destFile):
				destFile = baseDestFile + "_" + str(num)
				num += 1
			
			if not destFile.endswith('.dcm'):
				destFile=os.path.splitext(destFile)[0]+'.dcm'
			
			shutil.copyfile(f, destFile)
		#os.chmod(serieDir, 0o444)
	except:
		print("Exception in user code:")
		print('-'*60)
		traceback.print_exc(file=sys.stdout)
		print('-'*60)

def sorted_nicely(lst):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	sorted_lst = sorted(lst, key = alphanum_key)
	
	return sorted_lst


#%%

debug = False

if debug:
	class dotdict(dict):
		"""dot.notation access to dictionary attributes"""
		__getattr__ = dict.get
		__setattr__ = dict.__setitem__
		__delattr__ = dict.__delitem__
	
	class Namespace:
		def __init__(self, **kwargs):
			self.__dict__.update(kwargs)
	
	input=dotdict({
				'dicom_dir': '/home/greydon/Documents/data/SNSX_7T/sourcedata/Khan_NeuroAnalytics_20230127_2023_01_27_SNSX3_P110_1.A07D0585',
				})
	
	output=dotdict({
				'out_dir':'/home/greydon/Documents/data/SNSX_7T/sourcedata/out',
				})
	
	snakemake = Namespace(output=output, input=input)
	

sourcedata=r'/home/greydon/Documents/data/single/sourcedata/dicoms'
sub_num='sub-014'


for sub_num in sorted_nicely([x for x in os.listdir(sourcedata) if x.startswith('sub')]):
	
	sub_dir = os.path.join(sourcedata,sub_num)
	old_folders = [os.path.join(sub_dir,x) for x in os.listdir(sub_dir) if os.path.isdir(os.path.join(sub_dir,x)) and 'sorted_dicoms' not in x]
	sub_dir_out = os.path.join(sourcedata,sub_num,'sorted_dicoms')
	
	if not os.path.exists(sub_dir_out):
		# Aggregate input files
		aggregator = DicomAggregator()
		aggregator.add_dicom_sources(sub_dir)
		aggregator.aggregate()
		
		sessionPatientDirNames = []
		# Browse aggregated DICOM series
		for serie in aggregator._aggregate_sources.values():
			createTree(serie, sub_dir_out, sessionPatientDirNames)
		
		print(f"Finished sorting dicoms for {sub_num}")
		
		for ifold in old_folders:
			shutil.rmtree(ifold)



