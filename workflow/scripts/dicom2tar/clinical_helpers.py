# -*- coding: utf-8 -*-
"""
Created on Mon May 13 14:03:24 2019

@author: User
"""
import os
import numpy as np
import pandas as pd
from tempfile import mkdtemp
import  tarfile
import shutil
import csv
import datetime

def tarSession(args):
	output_dir = os.path.dirname(os.path.dirname(args.output_dir))
	tar_files = [f for f in os.listdir(args.output_dir) if f.endswith('.tar')]
	subjects = np.unique([x.split('_')[0] for x in tar_files])


	if os.path.exists(os.path.join(output_dir, 'clinical_events.tsv')) and not args.clinical_events:
		event_dates = pd.read_csv(os.path.join(output_dir, 'clinical_events.tsv'), sep='\t')
		event_dates = or_dates.sort_values(by=['subject']).reset_index(drop=True)
		event_dates.to_csv(os.path.join(output_dir, 'clinical_events.tsv'), sep='\t', index=False, na_rep='n/a', line_terminator="")
		
	for isub in subjects:
		tar_files_sub = [f for f in os.listdir(args.output_dir) if f.startswith(isub)]
		tar_files_sub_check = tar_files_sub[0].split('_')
		if len(tar_files_sub_check) > 7:
			print('Incorrect .tar filname for subject ' + isub)
			continue
		else:
			#TODO: fix for multiple OR dates
			# subject_or = [datetime.datetime.strptime(x, '%Y_%m_%d') for x in [y for y in or_dates[or_dates['subject']==isub]['or_date'].values] if x is not np.nan]
			# if len(subject_or)>1:
			#     subject_or = [subject_or[0]] # only date first OR date for now
			tar_files_dates = [datetime.datetime.strptime(x.split('_')[4], '%Y%m%d') for x in tar_files_sub]
			date_sort_idx = np.argsort(tar_files_dates)
			
			tar_files_sub = [tar_files_sub[i] for i in date_sort_idx]
			tar_files_dates = [tar_files_dates[i] for i in date_sort_idx]
			
			sessionDates = {'date':[],'session':[]}
			session = 1
			for idate in tar_files_dates:
				sessionDates['date'].append(idate.strftime('%Y%m%d'))
				sessionDates['session'].append(str(session).zfill(3))
				session += 1
					
			sessionDates = pd.DataFrame.from_dict(sessionDates)
			
			for i in range(len(tar_files_sub)):
				
				session = sessionDates.loc[i,'session']
				tmpdir = mkdtemp(prefix='DCM')
				
				tf = tarfile.open(os.path.join(args.output_dir, tar_files_sub[i]))
				tmembers = tf.getmembers()
				for tm in tmembers:
					tm.mode = 0o700
				tf_content = [m.name for m in tmembers if m.isfile()]
				fullTarPath = [os.path.join(tmpdir, f) for f in tf_content]
				tf.extractall(path=tmpdir, members=tmembers)
				tf.close()
				
				newName = os.path.join(args.output_dir, '_'.join([tar_files_sub[i].split('_')[0], session]+ tar_files_sub[i].split('_')[1:]))
				
				with tarfile.open(newName, "a") as tar:
					for filename in fullTarPath:
						arcNameOld = os.path.normpath(filename).split(os.sep)[-6:]
						firstLevel = '_'.join([arcNameOld[0].split('_')[0], session]+ arcNameOld[0].split('_')[1:])
						lastLevel = '_'.join([arcNameOld[-1].split('_')[0], session]+ arcNameOld[-1].split('_')[1:])
						arcname = '/'.join([firstLevel]+arcNameOld[1:-1] + [lastLevel])
						tar.add(filename, arcname=arcname)
				
				shutil.rmtree(tmpdir)
				os.remove(os.path.join(args.output_dir, tar_files_sub[i]))
				
				print('Finished subject ' + isub + ' session ' + session)