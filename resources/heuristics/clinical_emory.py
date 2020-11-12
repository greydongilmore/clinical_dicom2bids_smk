# -*- coding: utf-8 -*-

def create_key(template, outtype=('nii'), annotation_classes=None):
	if template is None or not template:
		raise ValueError('Template must be a valid format string')
	return (template, outtype, annotation_classes)

def get_unique(seqinfos, attr):
	"""Given a list of seqinfos, which must have come from a single study
	get specific attr, which must be unique across all of the entries
	If not -- fail!
	"""
	values = set(getattr(si, attr) for si in seqinfos)
#    assert (len(values) == 1)
	return values.pop()

def infotoids(seqinfos, outdir):
	
	subject = get_unique(seqinfos, 'example_dcm_file').split('_')[0]
	session = get_unique(seqinfos, 'example_dcm_file').split('_')[1]
	
	ids = {
	'locator': '',
	'session': session,
	'subject': subject,
	}
				
	return ids

def infotodict(seqinfo):
	"""Heuristic evaluator for determining which runs belong where

	allowed template fields - follow python string module:

	item: index within category
	subject: participant id
	seqitem: run number during scanning
	subindex: sub index within group
	"""
	
	#Anat
	t1w = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_run-{item:02d}_T1w')
	t1w_acq = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T1w')
	t2w = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T2w')
	
	#CT
	ct = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_run-{item:02d}_ct')
	
	#Diffusion
	dwi = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_run-{item:02d}_dwi')

	info = {t1w:[],
			t1w_acq:[],
			t2w:[],
			ct:[],
			dwi:[]
			}
	
	for idx,s in enumerate(seqinfo):
		if any(substring in s.series_description.upper() for substring in {'STEALTH','3D','STEREO'}) and not any(substring in s.series_description.upper() for substring in {'IR_FSPGR', 'FSPGR'}):
			if any(substring in s.series_description.upper() for substring in {'IR_FSPGR', 'FSPGR', 'IR-FSPGR'}):
				info[t1w_acq].append({'item': s.series_id, 'acq': 'FSPGR'})
			elif 'MPRAGE' in s.series_description.upper():
				info[t1w_acq].append({'item': s.series_id, 'acq': 'mprage'})
			else:
				info[t1w].append({'item': s.series_id})
		elif any(substring in s.series_description.upper() for substring in {'MPRAGE','MPR'}):
			if 'AX' in s.series_description.upper():
				info[t1w_acq].append({'item': s.series_id, 'acq': 'mprage2D'})
			else:
				info[t1w_acq].append({'item': s.series_id, 'acq': 'mprage3D'})
				
		elif any(substring in s.series_description.upper() for substring in {'IR_FSPGR', 'FSPGR','IR-FSPGR'}):
			info[t1w_acq].append({'item': s.series_id, 'acq': 'fspgr'})
		
		elif any(substring in s.series_description.upper() for substring in {'DWI', 'DTI', 'DIFFUSION'}):
			info[dwi].append({'item': s.series_id})

		elif any(substring in s.series_description.upper() for substring in {'AX', 'COR','SAG'}) and s.TR != -1 and s.TE != -1:
			if ('AX' in s.series_description.upper()):
				orientation = 'tra'
			elif ('COR' in s.series_description.upper()):
				orientation = 'cor'
			elif ('SAG' in s.series_description.upper()):
				orientation = 'sag'
			
		elif s.TR == -1 and s.TE == -1 or any(substring in s.series_description.upper() for substring in {'AX BRAIN THIN','AX BONE THIN'}):
			info[ct].append({'item': s.series_id})
	
	return info
