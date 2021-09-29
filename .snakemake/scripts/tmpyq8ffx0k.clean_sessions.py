
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/home/greydon/.local/lib/python3.8/site-packages', '/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/workflow/scripts/post_tar2bids']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95=\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8cQ/media/veracrypt6/projects/stealthMRI/working_dir/out/logs/sub-P197_tar2bids.done\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x0etouch_tar2bids\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x18)}\x94\x8c\x05_name\x94h\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8cS/media/veracrypt6/projects/stealthMRI/working_dir/out/logs/sub-P197_dicom2bids.done\x94a}\x94(h\x0c}\x94\x8c\x10touch_dicom2bids\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bh)h&ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8cE/media/veracrypt6/projects/stealthMRI/working_dir/clinical_events.tsv\x94\x8cG/media/veracrypt6/projects/stealthMRI/working_dir/out/bids_tmp/sub-P197\x94K\xd9}\x94(\x8c\x04peri\x94K\x00\x8c\x0fdur_multi_event\x94J\xe2\xff\xff\xff\x8c\roverride_peri\x94\x88u\x8c\x07patient\x94e}\x94(h\x0c}\x94(\x8c\x0fclinical_events\x94K\x00N\x86\x94\x8c\tbids_fold\x94K\x01N\x86\x94\x8c\x08num_subs\x94K\x02N\x86\x94\x8c\x08ses_calc\x94K\x03N\x86\x94\x8c\tsub_group\x94K\x04N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhAh8hCh9hEK\xd9hGh:hIh>ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94\x8c\x03197\x94a}\x94(h\x0c}\x94\x8c\x07subject\x94K\x00N\x86\x94sh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94b\x8c\x07subject\x94hXub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bhnK\x01hpK\x01ub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x12sNt\x94bh\x13h\x16h\x18\x85\x94R\x94(h\x18)}\x94h\x1ch\x13sNt\x94bub\x8c\x06config\x94}\x94(\x8c\tdicom_dir\x94\x8c)/media/veracrypt7/dicom_dbs_backup/dicoms\x94\x8c\x07out_dir\x94\x8c5/media/veracrypt6/projects/stealthMRI/working_dir/out\x94\x8c\x13clinical_event_file\x94h8\x8c\x10participants_tsv\x94N\x8c\theuristic\x94\x8c(resources/heuristics/clinical_imaging.py\x94\x8c\ndcm_config\x94\x8c\x19resources/dcm_config.json\x94\x8c\tanonymize\x94\x89\x8c\x0csession_calc\x94h:\x8c\tsub_group\x94h>\x8c\x0esubject_prefix\x94\x8c\x01P\x94u\x8c\x04rule\x94\x8c\rcleanSessions\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cU/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/workflow/scripts/post_tar2bids\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/home/greydon/Documents/GitHub/clinical_dicom2bids_smk/workflow/scripts/post_tar2bids/clean_sessions.py';
######## snakemake preamble end #########
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

def main():
	output_dir=os.path.dirname(os.path.dirname(snakemake.input.touch_tar2bids))

	final_dir = os.path.join(output_dir, 'bids')
	if not os.path.exists(final_dir):
		os.mkdir(final_dir)
	
	sub_num=''.join(filter(str.isdigit, os.path.basename(snakemake.params.bids_fold)))
	
	print('Converting subject {} ...'.format(os.path.basename(snakemake.params.bids_fold)))
	subject_event = []
	if snakemake.params.clinical_events:
		event_dates = pd.read_csv(snakemake.params.clinical_events, sep='\t')
		subject_event = [datetime.datetime.strptime(x, '%Y_%m_%d') for x in [y for y in event_dates[event_dates['subject']==int(sub_num)]['event_date'].values] if x is not np.nan]
	
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
						
					shutil.copyfile(os.path.join(snakemake.params.bids_fold, ises, iscan, ifile), new_file)
					
					if iscan+'/'+ifile in scans_data['filename'].values:
						name_idx = [i for i,x in enumerate(scans_data['filename'].values) if x == iscan+'/'+ifile][0]
						data_temp = scans_data.iloc[name_idx,:].to_dict()
						data_temp['filename']=iscan+'/'+os.path.basename(new_file)
						
						scans_tsv_new.append(data_temp)
			
			sub_code_path = make_bids_folders(isub.split('-')[1], ilabel, 'info', os.path.join(final_dir,'.heudiconv'), True, False)
			copytree(os.path.join(os.path.dirname(snakemake.params.bids_fold), '.heudiconv', isub.split('-')[1], ises,'info'), sub_code_path)
			
		scans_file = make_bids_filename(isub, 'ses-'+ilabel, None, None, None, 'scans.json', os.path.dirname(sub_path))
		scans_json = [x for x in os.listdir(os.path.join(snakemake.params.bids_fold, ises)) if x.endswith('scans.json')]
		shutil.copyfile(os.path.join(snakemake.params.bids_fold, ises, scans_json[0]), scans_file)
		
		scans_file = make_bids_filename(isub, 'ses-'+ilabel, None, None, None, 'scans.tsv', os.path.dirname(sub_path))
		scans_tsv_new = pd.DataFrame(scans_tsv_new)
		scans_tsv_new.to_csv(scans_file, sep='\t', index=False, na_rep='n/a', line_terminator="")
	
	# Check to see if this is the last subject complete, copy main BIDS files if so
	check_status = [x for x in os.listdir(os.path.join(output_dir,'logs')) if x.endswith('_dicom2bids.done')]
	if len(check_status)==(snakemake.params.num_subs)-1:
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
		
if __name__ == "__main__":

	main()
			
