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
			
def make_bids_filename(subject_id, session_id, task, acq, desc, run, suffix, prefix):
	if isinstance(session_id, str):
		if 'ses' in session_id:
			session_id = session_id.split('-')[1]
			
	order = OrderedDict([('ses', session_id if session_id is not None else None),
						 ('task', task if task is not None else None),
						 ('acq', acq if acq is not None else None),
						 ('desc', desc if desc is not None else None),
						 ('run', run if run is not None else None)])

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

def main():
	output_dir=os.path.dirname(os.path.dirname(snakemake.input.touch_tar2bids))

	final_dir = os.path.join(output_dir, 'bids')
	if not os.path.exists(final_dir):
		os.mkdir(final_dir)
	
	#os.remove(snakemake.input.touch_tar2bids)

	print('Converting subject {} ...'.format(os.path.basename(snakemake.params.bids_fold)))
	subject_event = []
	if snakemake.params.clinical_events:
		event_dates = pd.read_csv(snakemake.params.clinical_events, sep='\t')
		subject_event = [datetime.datetime.strptime(x, '%Y_%m_%d') for x in [y for y in event_dates[event_dates['subject']==os.path.basename(snakemake.params.bids_fold).split('-')[1]]['event_date'].values] if x is not np.nan]
	
	orig_sessions = sorted_nicely([x for x in os.listdir(snakemake.params.bids_fold) if os.path.isdir(os.path.join(snakemake.params.bids_fold, x)) and 'ses' in x])
	
	sessionDates = {'ses_num':[],'session':[]}
	for ises in orig_sessions:
		if not subject_event:
			sessionDates['ses_num'].append(ises)
			sessionDates['session'].append('pre')
		else:
			scans_tsv = [x for x in os.listdir(os.path.join(snakemake.params.bids_fold, ises)) if x.endswith('scans.tsv')]
			scans_data = pd.read_table(os.path.join(snakemake.params.bids_fold, ises, scans_tsv[0]))
			idate = datetime.datetime.strptime(scans_data['acq_time'].values[0].split('T')[0], '%Y-%m-%d')
			dateAdded = False
			for ievent in subject_event:
				if dateAdded:
					if 'post' in sessionDates['session'][-1]:
						if abs((idate-ievent).days) <= snakemake.params.ses_calc['peri']:
							post_scan=False
							if snakemake.params.ses_calc['override_peri']:
								for root, folders, files in os.walk(os.path.join(snakemake.params.bids_fold, ises)):
									for file in files:
										if 'electrode' in file.lower():
											post_scan = True

							if post_scan:
								sessionDates['session'][-1]='post'
							else:
								sessionDates['session'][-1]='peri'

						elif snakemake.params.ses_calc['dur_multi_event'] < (idate-ievent).days < 0:
							sessionDates['session'][-1]='pre'
				else:
					sessionDates['ses_num'].append(ises)
					if (idate-ievent).days > 0:
						sessionDates['session'].append('post')
					elif abs((idate-ievent).days) <= snakemake.params.ses_calc['peri']:
						post_scan=False
						if snakemake.params.ses_calc['override_peri']:
							for root, folders, files in os.walk(os.path.join(snakemake.params.bids_fold, ises)):
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
	isub = os.path.basename(snakemake.params.bids_fold)
	for ilabel in sessionDates.session.unique():
		sessions = sessionDates[sessionDates['session']==ilabel]['ses_num'].values
		scans_tsv_new = []
		for ises in sessions:
			scans_tsv = [x for x in os.listdir(os.path.join(snakemake.params.bids_fold, ises)) if x.endswith('scans.tsv')]
			scans_data = pd.read_table(os.path.join(snakemake.params.bids_fold, ises, scans_tsv[0]))
			scan_type = [x for x in os.listdir(os.path.join(snakemake.params.bids_fold, ises)) if os.path.isdir(os.path.join(snakemake.params.bids_fold, ises, x))]
			for iscan in scan_type:
				sub_path = make_bids_folders(isub, ilabel, iscan, final_dir, True, False)
				files = [x for x in os.listdir(os.path.join(snakemake.params.bids_fold, ises, iscan)) if os.path.isfile(os.path.join(snakemake.params.bids_fold, ises, iscan, x))]
				
				for ifile in files:

					key_dict={
						'task': [],
						'acq': [],
						'desc': [],
						'run':[]
					}

					key_dict['suffix']=ifile.split('_')[-1]

					for ikey in key_dict.keys():
						key_dict[ikey]=ifile.split(f'{ikey}-')[1].split('_')[0] if f'{ikey}-' in ifile else None

					key_dict['suffix']=ifile.split('_')[-1]
					key_dict['prefix']=sub_path

					new_file = make_bids_filename(isub, 'ses-'+ilabel, **key_dict)
					
					shutil.copyfile(os.path.join(snakemake.params.bids_fold, ises, iscan, ifile), new_file)
										
					if iscan+'/'+ifile in scans_data['filename'].values:
						name_idx = [i for i,x in enumerate(scans_data['filename'].values) if x == iscan+'/'+ifile][0]
						data_temp = scans_data.iloc[name_idx,:].to_dict()
						data_temp['filename']=iscan+'/'+os.path.basename(new_file)
						
						scans_tsv_new.append(data_temp)
			
			sub_code_path = make_bids_folders(isub.split('-')[1], ilabel, 'info', os.path.join(final_dir,'.heudiconv'), True, False)
			copytree(os.path.join(os.path.dirname(snakemake.params.bids_fold), '.heudiconv', isub.split('-')[1], ises,'info'), sub_code_path)
			
		scans_file = make_bids_filename(isub, 'ses-' + ilabel, None, None, None, None, 'scans.json', os.path.dirname(sub_path))
		scans_json = [x for x in os.listdir(os.path.join(snakemake.params.bids_fold, ises)) if x.endswith('scans.json')]
		if scans_json:
			shutil.copyfile2(os.path.join(snakemake.params.bids_fold, ises, scans_json[0]), scans_file)
		
		scans_file = make_bids_filename(isub, 'ses-'+ilabel, None, None, None, None, 'scans.tsv', os.path.dirname(sub_path))
		scans_tsv_new = pd.DataFrame(scans_tsv_new)
		scans_tsv_new.to_csv(scans_file, sep='\t', index=False, na_rep='n/a', line_terminator="")
	
	# Check to see if this is the last subject complete, copy main BIDS files if so
	check_status = [x for x in os.listdir(os.path.join(output_dir,'bids_tmp')) if os.path.isdir(os.path.join(output_dir, 'bids_tmp', x)) and not x.startswith('.')]
	if len(check_status)==(snakemake.params.num_subs):
		bids_files = [x for x in os.listdir(os.path.join(output_dir, 'bids_tmp')) if os.path.isfile(os.path.join(output_dir, 'bids_tmp', x))]
		for ifile in bids_files:
			if ifile == 'participants.tsv':
				if os.path.exists(os.path.join(final_dir, ifile)):
					patient_tsv_old = pd.read_csv(os.path.join(output_dir, 'bids_tmp', 'participants.tsv'), sep='\t')
					patient_tsv = pd.read_csv(os.path.join(final_dir, 'participants.tsv'), sep='\t')
					patient_tsv = patient_tsv.append(patient_tsv_old).reset_index(drop=True)
				else:
					patient_tsv = pd.read_csv(os.path.join(output_dir, 'bids_tmp', 'participants.tsv'), sep='\t')

				patient_tsv = patient_tsv.sort_values(by=['participant_id']).reset_index(drop=True)
				patient_tsv['group'] = patient_tsv['group'].replace('control',snakemake.params.sub_group)
				patient_tsv.to_csv(os.path.join(final_dir, ifile), sep='\t', index=False, na_rep='n/a', line_terminator="")
			else:
				if not os.path.exists(os.path.join(final_dir, ifile)):
					shutil.copyfile(os.path.join(output_dir, 'bids_tmp', ifile), os.path.join(final_dir, ifile))

		shutil.rmtree(os.path.join(output_dir, 'bids_tmp'))
	else:
		shutil.rmtree(os.path.join(output_dir, 'bids_tmp',isub))
		
if __name__ == "__main__":

	main()
			
