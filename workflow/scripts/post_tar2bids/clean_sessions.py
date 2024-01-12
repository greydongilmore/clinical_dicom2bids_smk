#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 18:03:02 2020

@author: greydon
"""
import os
import datetime
import pandas as pd
import numpy as np
import re
from collections import OrderedDict
import shutil
import glob

class Namespace:
	def __init__(self, **kwargs):
		self.__dict__.update(kwargs)

def sorted_nicely( l ):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key = alphanum_key)

def copytree(src, dst, symlinks=False, ignore=None):
	for item in os.listdir(src):
		s = os.path.join(src, item)
		d = os.path.join(dst, item)
		if os.path.isdir(s):
			shutil.copytree(s, d, symlinks, ignore)
		else:
			shutil.copy2(s, d)
			
def make_bids_filename(subject_id, session_id, suffix, prefix, task=None, acq=None, ce=None, rec=None, 
	dir=None, mod=None, echo=None, hemi=None, space=None,res=None,den=None,label=None,part=None,desc=None,run=None):
	if isinstance(session_id, str):
		if 'ses' in session_id:
			session_id = session_id.split('-')[1]

	order = OrderedDict([('ses', session_id),
						('task', task),
						('acq', acq),
						('ce', ce),
						('rec', rec),
						('dir', dir),
						('mod', mod),
						('echo', echo),
						('hemi', hemi),
						('space', space),
						('res', res),
						('den', den),
						('label', label),
						('part', part),
						('desc', desc),
						('run', run)])
	filename = []
	if subject_id is not None:
		filename.append(subject_id)
	for key, val in order.items():
		if val is not None:
			filename.append('%s-%s' % (key, val))

	if isinstance(suffix, str):
		filename.append(suffix)

	filename = '_'.join(filename)
	if isinstance(prefix, str):
		filename = os.path.join(prefix, filename)
		
	return filename

def make_bids_folders(subject_id, session_id, kind, output_path, make_dir, overwrite):
	path = []
	path.append(subject_id)
		
	if isinstance(session_id, str):
		if 'ses' not in session_id:
			path.append('ses-%s' % session_id)
		else:
			path.append(session_id)
			
	if isinstance(kind, str):
		path.append(kind)
	
	path = os.path.join(*path)  
	path = os.path.join(output_path, path)
		
	if make_dir == True:
		if not os.path.exists(path):
			os.makedirs(path)
		elif overwrite:
			shutil.rmtree(path)
			os.makedirs(path)
			
	return path

import time
def date_parser(string_list):
	return [time.ctime(float(x)) for x in string_list]

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
	
	isub='sub-P109'
	data_dir=r'/home/greydon/Documents/data/SEEG/derivatives/atlasreg'
	
	input=dotdict({
				'tmp_dir': f'/home/greydon/Documents/datasets/DBS/bids_tmp/sub-P309',
				})
	
	config={
		'bids_dir':'/home/greydon/Documents/datasets/DBS',
		'remove_temp':True,
		'sub_group':'patient',
		'ses_calc':{
			'peri':0,
			'dur_multi_event':-30,
			'override_peri':True
			}
		}
	
	snakemake = Namespace(input=input,config=config)

clinical_event_file=os.path.join(snakemake.config['bids_dir'], 'events.tsv')

final_dir = os.path.join(snakemake.config['bids_dir'], 'bids')
if not os.path.exists(final_dir):
	os.mkdir(final_dir)

#os.remove(snakemake.input.touch_tar2bids)
sub_num=''.join([x for x in os.path.basename(snakemake.input.tmp_dir) if x.isnumeric()])
isub = os.path.basename(snakemake.input.tmp_dir)

print('Converting subject {} ...'.format(os.path.basename(snakemake.input.tmp_dir)))

subject_event = []
if os.path.exists(clinical_event_file):
	event_dates = pd.read_csv(clinical_event_file, sep='\t',dtype = str)
	subject_event = [datetime.datetime.strptime(x, '%Y_%m_%d') for x in [y for y in event_dates[event_dates['subject']==sub_num]['event_date'].values] if x is not np.nan]

orig_sessions = sorted_nicely([x for x in os.listdir(snakemake.input.tmp_dir) if os.path.isdir(os.path.join(snakemake.input.tmp_dir, x)) and 'ses' in x])

#this loop just assigns each session to 'pre','peri','post'
#for now, only need to consider the day ignore time
sessionDates = {'ses_num':[],'session':[]}
for ises in orig_sessions:
	if not subject_event:
		try:
			sessionDates['ses_num']._append(ises)
			sessionDates['session']._append('pre')
		except:
			sessionDates['ses_num'].append(ises)
			sessionDates['session'].append('pre')
	else:
		scans_tsv = [x for x in os.listdir(os.path.join(snakemake.input.tmp_dir, ises)) if x.endswith('scans.tsv')]
		scans_data = pd.read_table(os.path.join(snakemake.input.tmp_dir, ises, scans_tsv[0]))
		if isinstance(scans_data['acq_time'].values[0],np.float64):
			sessionDates['ses_num'].append(ises)
			sessionDates['session'].append('pre')
		else:
			#idate = datetime.datetime.strptime(scans_data['acq_time'].values[0], '%Y-%m-%dT%H:%M:%S')
			idate = datetime.datetime.strptime(scans_data['acq_time'].values[0].split('T')[0], '%Y-%m-%d')
			dateAdded = False
			for ievent in subject_event:
				if dateAdded:
					if 'post' in sessionDates['session'][-1]:
						if abs((idate-ievent).days) <= snakemake.config['ses_calc']['peri']:
							post_scan=False
							if snakemake.config['ses_calc']['override_peri']:
								for root, folders, files in os.walk(os.path.join(snakemake.input.tmp_dir, ises)):
									for file in files:
										if 'electrode' in file.lower():
											post_scan = True

							if post_scan:
								sessionDates['session'][-1]='post'
							else:
								sessionDates['session'][-1]='peri'

						elif snakemake.config['ses_calc']['dur_multi_event'] < (idate-ievent).days < 0:
							sessionDates['session'][-1]='pre'
				else:
					sessionDates['ses_num'].append(ises)
					if (idate-ievent).days > 0:
						sessionDates['session'].append('post')
					elif abs((idate-ievent).days) <= snakemake.config['ses_calc']['peri']:
						post_scan=False
						if snakemake.config['ses_calc']['override_peri']:
							for root, folders, files in os.walk(os.path.join(snakemake.input.tmp_dir, ises)):
								for file in files:
									if 'electrode' in file.lower():
										post_scan = True
						
						if post_scan:
							sessionDates['session'].append('post')
						else:
							sessionDates['session'].append('peri')
							
					elif (idate-ievent).days < 0:
						sessionDates['session'].append('pre')
						
					dateAdded=True
	
sessionDates = pd.DataFrame.from_dict(sessionDates)

for ilabel in sessionDates.session.unique():
	sessions = sessionDates[sessionDates['session']==ilabel]['ses_num'].values
	scans_tsv_new = []
	for ises in sessions:
		sub_code_path = make_bids_folders(isub.split('-')[1], ilabel, 'info', os.path.join(final_dir,'.heudiconv'), True, False)
		scans_tsv = [x for x in os.listdir(os.path.join(snakemake.input.tmp_dir, ises)) if x.endswith('scans.tsv')]
		scans_data = pd.read_csv(os.path.join(snakemake.input.tmp_dir, ises, scans_tsv[0]),sep='\t',parse_dates=True,date_format='%Y-%m-%dT%H:%M:%S')
		scans_data['acq_time']=[pd.Timestamp(d).round(freq='S') for d in scans_data['acq_time']]
	
		for ifile in scans_data['filename']:
			
			sub_path = make_bids_folders(isub, ilabel, ifile.split('/')[0], final_dir, True, False)
			ises_path = make_bids_folders(isub, ises, ifile.split('/')[0], os.path.dirname(snakemake.input.tmp_dir), False, False)
			
			key_dict={
				'task': [],
				'acq': [],
				'ce':[],
				'rec':[],
				'dir':[],
				'mod':[],
				'echo':[],
				'hemi':[],
				'space':[],
				'res':[],
				'den':[],
				'desc': [],
				'label':[],
				'part':[],
				'run':[],
			}
			
			
			
			for ikey in key_dict.keys():
				key_dict[ikey]=ifile.split(f'{ikey}-')[1].split('_')[0] if f'{ikey}-' in ifile else None
			
			key_dict['suffix']=ifile.split('_')[-1].split('.nii')[0]
			key_dict['prefix']=sub_path
			
			new_file_novel = make_bids_filename(isub, 'ses-'+ilabel, **key_dict)
			
			while os.path.exists(new_file_novel+'.nii.gz'):
				run_num=int(new_file_novel.split('run-')[1].split('_')[0])
				key_dict['run']=str(run_num+1).zfill(2)
				new_file_novel = make_bids_filename(isub, 'ses-'+ilabel, **key_dict)
			
			
			for iallfiles in glob.glob(os.path.join(ises_path, ifile.split('/')[-1].split('.nii')[0])+ '*'):
				old_ending='_'+new_file_novel.split('_')[-1]
				new_ending='_'+os.path.basename(iallfiles).split('_')[-1]
				shutil.copy(iallfiles, new_file_novel.replace(old_ending,new_ending))
				os.chmod(new_file_novel.replace(old_ending,new_ending), 0o777)
			
			data_tmp=scans_data[scans_data['filename']==ifile].reset_index(drop=True)
			data_tmp.loc[0,'filename']='/'.join([ifile.split('/')[0], os.path.basename(new_file_novel+'.nii.gz')])
			scans_tsv_new.append(data_tmp)
		
		copytree(os.path.join(os.path.dirname(snakemake.input.tmp_dir), '.heudiconv', isub.split('-')[1], ises,'info'), sub_code_path)
	
	if glob.glob(os.path.join(os.path.dirname(snakemake.input.tmp_dir),'*scans.json')):
		scans_json = make_bids_filename(isub, 'ses-' + ilabel, 'scans.json', os.path.dirname(sub_path))
		shutil.copy(glob.glob(os.path.join(os.path.dirname(snakemake.input.tmp_dir),'*scans.json'))[0], scans_json)
	
	scans_file = make_bids_filename(isub, 'ses-'+ilabel, 'scans.tsv', os.path.dirname(sub_path))
	scans_tsv_new = pd.concat(scans_tsv_new)
	scans_tsv_new.sort_values(["acq_time"]).to_csv(scans_file, sep='\t', index=False, na_rep='n/a',date_format='%Y-%m-%dT%H:%M:%S')

# Check to see if this is the last subject complete, copy main BIDS files if so
check_status = [x for x in os.listdir(os.path.join(snakemake.config['bids_dir'],'bids_tmp')) if os.path.isdir(os.path.join(snakemake.config['bids_dir'], 'bids_tmp', x)) and not x.startswith('.')]
if len(check_status)==1:
	bids_files = [x for x in os.listdir(os.path.join(snakemake.config['bids_dir'], 'bids_tmp')) if os.path.isfile(os.path.join(snakemake.config['bids_dir'], 'bids_tmp', x))]
	for ifile in bids_files:
		if ifile == 'participants.tsv':
			if os.path.exists(os.path.join(final_dir, ifile)):
				patient_tsv_old = pd.read_csv(os.path.join(snakemake.config['bids_dir'], 'bids_tmp', 'participants.tsv'), sep='\t')
				patient_tsv = pd.read_csv(os.path.join(final_dir, 'participants.tsv'), sep='\t')
				try:
					patient_tsv = patient_tsv._append(patient_tsv_old).reset_index(drop=True)
				except:
					patient_tsv = patient_tsv.append(patient_tsv_old).reset_index(drop=True)
			else:
				patient_tsv = pd.read_csv(os.path.join(snakemake.config['bids_dir'], 'bids_tmp', 'participants.tsv'), sep='\t')

			patient_tsv=patient_tsv.drop_duplicates(subset='participant_id', keep="last")
			patient_tsv = patient_tsv.sort_values(by=['participant_id']).reset_index(drop=True)
			patient_tsv['group'] = patient_tsv['group'].replace('control',snakemake.config['sub_group'])
			patient_tsv.to_csv(os.path.join(final_dir, ifile), sep='\t', index=False, na_rep='n/a', float_format='%.0f')
		else:
			if not os.path.exists(os.path.join(final_dir, ifile)):
				shutil.copy(os.path.join(snakemake.config['bids_dir'], 'bids_tmp', ifile), os.path.join(final_dir, ifile))
	if snakemake.config['remove_temp']:
		shutil.rmtree(os.path.join(snakemake.config['bids_dir'], 'bids_tmp'))
else:
	if snakemake.config['remove_temp']:
		shutil.rmtree(os.path.join(snakemake.config['bids_dir'], 'bids_tmp',isub))
		
