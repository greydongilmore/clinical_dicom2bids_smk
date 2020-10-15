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
			
def make_bids_filename(subject_id, session_id, task_id, acq_id, run_num, suffix, prefix):
	if isinstance(session_id, str):
		if 'ses' in session_id:
			session_id = session_id.split('-')[1]
			
	order = OrderedDict([('ses', session_id if session_id is not None else None),
						 ('task', task_id if task_id is not None else None),
						 ('acq', acq_id if acq_id is not None else None),
						 ('run', run_num if run_num is not None else None)])

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
# config = {
# 	'periop': 0,
# 	'dur_multi_surg': -30,
# 	'override_periop': True
# 	}
# args = Namespace(tar2bids_done = '/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/data/output/sub-P197_tar2bids.done', 
# 				 output_dir=os.path.dirname('/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/data/output/sub-P197_tar2bids.done'), 
# 				 or_dates_file='',
# 				 bids_fold='/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/data/output/bids/sub-P197',
# 				 num_subs=3, 
# 				 ses_calc=config)

def main():

	args = Namespace(tar2bids_done = snakemake.input[0], output_dir=os.path.dirname(snakemake.input[0]), or_dates_file=snakemake.params[0], 
				  bids_fold=snakemake.params[1],num_subs=snakemake.params[2], ses_calc=snakemake.params[3])

	if not args.or_dates_file:
		or_dates = pd.read_csv(os.path.join(args.output_dir, 'or_dates.tsv'), sep='\t')
	else:
		or_dates = pd.read_csv(args.or_dates_file, sep='\t')
		
	final_dir = os.path.join(args.output_dir, 'bids')
	if not os.path.exists(final_dir):
		os.mkdir(final_dir)

	# Check to see if this is the last subject complete, copy main BIDS files if so
	check_status = [x for x in os.listdir(args.output_dir) if x.endswith('_dicom2bids.done')]
	if len(check_status)==(args.num_subs)-1:
		bids_files = [x for x in os.listdir(os.path.join(args.output_dir, 'bids_tmp')) if os.path.isfile(os.path.join(args.output_dir, 'bids_tmp', x))]
		for ifile in bids_files:
			if ifile == 'participants.tsv':
				patient_tsv = pd.read_csv(os.path.join(args.output_dir, 'bids_tmp', 'participants.tsv'), sep='\t')
				patient_tsv = patient_tsv.sort_values(by=['participant_id']).reset_index(drop=True)
				patient_tsv.to_csv(os.path.join(final_dir, ifile), sep='\t', index=False, na_rep='n/a', line_terminator="")
			else:
				shutil.copyfile(os.path.join(args.output_dir, 'bids_tmp', ifile), os.path.join(final_dir, ifile))
	
	os.remove(args.tar2bids_done)

	print('Converting subject {} ...'.format(os.path.basename(args.bids_fold)))
	subject_or = [datetime.datetime.strptime(x, '%Y_%m_%d') for x in [y for y in or_dates[or_dates['subject']==os.path.basename(args.bids_fold).split('-')[1]]['or_date'].values] if x is not np.nan]
	orig_sessions = sorted_nicely([x for x in os.listdir(args.bids_fold) if os.path.isdir(os.path.join(args.bids_fold, x)) and 'ses' in x])
	
	sessionDates = {'ses_num':[],'session':[]}
	for ises in orig_sessions:
		if not subject_or:
			sessionDates['ses_num'].append(ises)
			sessionDates['session'].append('presurg')
		else:
			scans_tsv = [x for x in os.listdir(os.path.join(args.bids_fold, ises)) if x.endswith('scans.tsv')]
			scans_data = pd.read_table(os.path.join(args.bids_fold, ises, scans_tsv[0]))
			idate = datetime.datetime.strptime(scans_data['acq_time'].values[0].split('T')[0], '%Y-%m-%d')
			dateAdded = False
			for ior in subject_or:
				if dateAdded:
					if 'postsurg' in sessionDates['session'][-1]:
						if abs((idate-ior).days) <= args.ses_calc['periop']:
							postop_scan=False
							if args.ses_calc['override_periop']:
								for root, folders, files in os.walk(os.path.join(args.bids_fold, ises)):
									for file in files:
										if 'electrode' in file.lower():
											postop_scan = True

							if postop_scan:
								sessionDates['session'][-1]='postsurg'
							else:
								sessionDates['session'][-1]='perisurg'

						elif args.ses_calc['dur_multi_surg'] < (idate-ior).days < 0:
							sessionDates['session'][-1]='presurg'
				else:
					sessionDates['ses_num'].append(ises)
					if (idate-ior).days > 0:
						sessionDates['session'].append('postsurg')
					elif abs((idate-ior).days) <= args.ses_calc['periop']:
						postop_scan=False
						if args.ses_calc['override_periop']:
							for root, folders, files in os.walk(os.path.join(args.bids_fold, ises)):
								for file in files:
									if 'electrode' in file.lower():
										postop_scan = True

						if postop_scan:
							sessionDates['session'].append('postsurg')
						else:
							sessionDates['session'].append('perisurg')
							
					elif (idate-ior).days < 0:
						sessionDates['session'].append('presurg')
						
					dateAdded=True
		
	sessionDates = pd.DataFrame.from_dict(sessionDates)
	isub = os.path.basename(args.bids_fold)
	for ilabel in sessionDates.session.unique():
		sessions = sessionDates[sessionDates['session']==ilabel]['ses_num'].values
		scans_tsv_new = []
		for ises in sessions:
			scans_tsv = [x for x in os.listdir(os.path.join(args.bids_fold, ises)) if x.endswith('scans.tsv')]
			scans_data = pd.read_table(os.path.join(args.bids_fold, ises, scans_tsv[0]))
			scan_type = [x for x in os.listdir(os.path.join(args.bids_fold, ises)) if os.path.isdir(os.path.join(args.bids_fold, ises, x))]
			for iscan in scan_type:
				sub_path = make_bids_folders(isub, ilabel, iscan, final_dir, True, False)
				files = [x for x in os.listdir(os.path.join(args.bids_fold, ises, iscan)) if os.path.isfile(os.path.join(args.bids_fold, ises, iscan, x))]
				
				for ifile in files:
					acq_id = ifile.split('acq-')[1].split('_')[0] if 'acq-' in ifile else None
					task_id = ifile.split('task-')[1].split('_')[0] if 'task-' in ifile else None
					if acq_id is not None:
						if task_id is not None:
							number_files = [x for x in os.listdir(sub_path) if x.endswith(ifile.split('_')[-1]) and 'acq-'+acq_id in x and 'task-'+task_id in x]
						else:
							number_files = [x for x in os.listdir(sub_path) if x.endswith(ifile.split('_')[-1]) and 'acq-'+acq_id in x and 'task-' not in x]
					else:
						if task_id is not None:
							number_files = [x for x in os.listdir(sub_path) if x.endswith(ifile.split('_')[-1]) and 'task-'+task_id in x and 'acq-' not in x]
						else:
							number_files = [x for x in os.listdir(sub_path) if x.endswith(ifile.split('_')[-1]) and 'acq-' not in x and 'task-' not in x]
							
					new_file = make_bids_filename(isub, 'ses-'+ilabel, task_id, acq_id, str(len(number_files)+1).zfill(2), ifile.split('_')[-1], sub_path)
						
					shutil.copyfile(os.path.join(args.bids_fold, ises, iscan, ifile), new_file)
					
					if iscan+'/'+ifile in scans_data['filename'].values:
						name_idx = [i for i,x in enumerate(scans_data['filename'].values) if x == iscan+'/'+ifile][0]
						data_temp = scans_data.iloc[name_idx,:].to_dict()
						data_temp['filename']=iscan+'/'+os.path.basename(new_file)
						
						scans_tsv_new.append(data_temp)
			
			sub_code_path = make_bids_folders(isub.split('-')[1], ilabel, 'info', os.path.join(final_dir,'.heudiconv'), True, False)
			copytree(os.path.join(os.path.dirname(args.bids_fold), '.heudiconv', isub.split('-')[1], ises,'info'), sub_code_path)
			
		scans_file = make_bids_filename(isub, 'ses-'+ilabel, None, None, None, 'scans.json', os.path.dirname(sub_path))
		scans_json = [x for x in os.listdir(os.path.join(args.bids_fold, ises)) if x.endswith('scans.json')]
		shutil.copyfile(os.path.join(args.bids_fold, ises, scans_json[0]), scans_file)
		
		scans_file = make_bids_filename(isub, 'ses-'+ilabel, None, None, None, 'scans.tsv', os.path.dirname(sub_path))
		scans_tsv_new = pd.DataFrame(scans_tsv_new)
		scans_tsv_new.to_csv(scans_file, sep='\t', index=False, na_rep='n/a', line_terminator="")
				
if __name__ == "__main__":

	main()
			
