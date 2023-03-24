# -*- coding: utf-8 -*-
import regex as re

diff_dic={
	'ADC':'ADC',
	'TRACEW':'TRACEW',
	'WEIGHTED TRACE':'TRACEW',
	'COLFA':'COLFA',
	'FA':'FA',
	'TENSOR':'TENSOR',
	'TENSOR_B0':'TENSOR',
}

def create_key(template, outtype=('nii.gz'), annotation_classes=None):
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

# Utility:  Check a list of regexes for truthyness
def regex_search_label(regexes, label):
	if any(regex.search(label) for regex in regexes):
		out_list=[regex.search(label)[0] for regex in regexes if regex.search(label) is not None]
		return [x for x in out_list if x !='']
	else:
		return []

def is_diffusion_derived(label):
	regexes = [
		re.compile('ADC$', re.IGNORECASE),
		re.compile('TRACEW$', re.IGNORECASE),
		re.compile('WEIGHTED TRACE$', re.IGNORECASE),
		re.compile('COLFA$', re.IGNORECASE),
		re.compile('FA$', re.IGNORECASE),
		re.compile('TENSOR[_B0]*', re.IGNORECASE)
		]
	return regex_search_label(regexes, label)

def is_swi_derived(label):
	regexes = [
		re.compile('Mag$', re.IGNORECASE),
		re.compile('Pha$', re.IGNORECASE),
		re.compile('mIP$', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

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
	t1w_pd = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_PD')
	t1w_flair = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_FLAIR')
	fa = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_run-{item:02d}_angio')
	fa_acq = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_angio')
	
	#fmap
	fmap = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_run-{item:02d}_epi')
	fmap_acq = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_epi')

	#Diffusion
	dwi = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_run-{item:02d}_dwi')
	dwi_acq = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_dwi')
	
	#CT
	ct = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_run-{item:02d}_ct')
	ct_acq = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_ct')
	ct_acq_desc = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_acq-{acq}_desc-{desc}_run-{item:02d}_ct')
	
	#pet
	pet = create_key('{bids_subject_session_dir}/pet/{bids_subject_session_prefix}_task-rest_run-{item:02d}_pet')
	pet_acq = create_key('{bids_subject_session_dir}/pet/{bids_subject_session_prefix}_task-rest_acq-{acq}_run-{item:02d}_pet')
	pet_task = create_key('{bids_subject_session_dir}/pet/{bids_subject_session_prefix}_task-{task}_acq-{acq}_run-{item:02d}_pet')

	#MIPS
	mips_sag = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-SWI_rec-DIS2D_run-{item:02d}_sagMIP')
	mips_cor = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-SWI_rec-DIS2D_run-{item:02d}_corMIP')
	mips_tra = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-SWI_rec-DIS2D_run-{item:02d}_traMIP')

	#GRE (Susc3D)
	swi_gre = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-SWI_run-{item:02d}_GRE')
	swi_gre_elec = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-ElectrodeSWI_run-{item:02d}_GRE')
	mag_gre = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_part-mag_run-{item:02d}_GRE')
	phase_gre = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_part-phase_run-{item:02d}_GRE')

	info = {t1w:[],
			t1w_acq:[],
			t2w:[],
			t1w_pd:[],
			t1w_flair:[],
			fmap:[],
			fmap_acq:[],
			dwi:[],
			dwi_acq:[],
			fa:[],
			fa_acq:[],
			ct:[],
			ct_acq:[],
			ct_acq_desc:[],
			pet:[],
			pet_acq:[],
			pet_task:[],
			mips_sag:[],
			mips_cor:[],
			mips_tra:[],
			swi_gre:[],
			swi_gre_elec:[],
			mag_gre:[],
			phase_gre:[],
	}
	
	for idx, s in enumerate(seqinfo):
		if any(substring in s.study_description.upper() for substring in {'MR'}):
			postop = False
			if 'SAR' in s.series_description.upper() or any(x in s.protocol_name.upper() for x in {'SAFE', 'STIMULATOR', 'STIM SAFE', 'POST', 'POST OP','POST-OP'}):
				if not any(x.upper() in s.protocol_name.upper() for x in {'POST STROKE','GAD','+C','STEALTH POST','MPRAGE POST'}):
					postop = True
				
			if any(substring in s.series_description.upper() for substring in {'STEALTH','3D','STEREO','STEREO-INCLUDE','1.5 MM ANATOMY','MPRAGE'}) and not any(substring in s.series_description.upper() for substring in {'IR_FSPGR', 'FSPGR'}):
				if any(substring in s.series_description.upper() for substring in {'IR_FSPGR', 'FSPGR', 'IR-FSPGR','3D T1 BRAVO'}):
					if postop:
						info[t1w_acq].append({'item': s.series_id, 'acq': 'ElectrodeFSPGR'})
					else:
						info[t1w_acq].append({'item': s.series_id, 'acq': 'FSPGR'})
				elif 'MPGR' in s.series_description.upper():
					if postop:
						info[t1w_acq].append({'item': s.series_id, 'acq': 'ElectrodeMPGR'})
					elif 'AX' in s.series_description.upper():
						info[t1w_acq].append({'item': s.series_id, 'acq': 'MPGR2D'})
					else:
						info[t1w_acq].append({'item': s.series_id, 'acq': 'MPGR3D'})
				else:
					if postop:
						if 'T2' in s.series_description.upper():
							info[t2w].append({'item': s.series_id, 'acq': 'Electrode3D'})
						else:
							info[t1w_acq].append({'item': s.series_id, 'acq': 'Electrode3D'})
					else:
						info[t1w].append({'item': s.series_id})
					
			elif 'EPI' in s.series_description.upper():
				if postop:
					info[fmap_acq].append({'item': s.series_id, 'acq': 'Electrode'})
				elif 'T2' in s.series_description.upper():
					info[fmap].append({'item': s.series_id})
					
			elif 'MPGR' in s.series_description.upper():
				if postop:
					info[t1w_acq].append({'item': s.series_id, 'acq': 'ElectrodeMPGR'})
				elif 'AX' in s.series_description.upper():
					info[t1w_acq].append({'item': s.series_id, 'acq': 'MPGR2D'})
				else:
					info[t1w_acq].append({'item': s.series_id, 'acq': 'MPGR3D'})
			
			elif 'RAGE' in s.series_description.upper():
				if postop:
					info[t1w_acq].append({'item': s.series_id, 'acq': 'ElectrodeMPRAGE'})
				elif 'AX' in s.series_description.upper():
					info[t1w_acq].append({'item': s.series_id, 'acq': 'MPRAGE2D'})
				else:
					info[t1w_acq].append({'item': s.series_id, 'acq': 'MPRAGE'})

			elif any(substring in s.series_description.upper() for substring in {'IR_FSPGR', 'FSPGR','IR-FSPGR'}):
				if postop:
					info[t1w_acq].append({'item': s.series_id, 'acq': 'ElectrodeFSPGR'})
				else:
					info[t1w_acq].append({'item': s.series_id, 'acq': 'FSPGR'})
			
			elif any(substring in s.series_description.upper() for substring in {'SWI'}):
				if len(is_swi_derived(s.series_description))>=1:
					str_derv=is_swi_derived(s.series_description)[0].upper()
					if str_derv == 'MAG':
						info[mag_gre].append({'item': s.series_id})
					elif str_derv == 'PHA':
						info[phase_gre].append({'item': s.series_id})
					elif str_derv == 'MIP':
						if 'AX' in s.series_description.upper():
							info[mips_tra].append({'item': s.series_id})
						elif 'COR' in s.series_description.upper():
							info[mips_cor].append({'item': s.series_id})
						elif 'SAG' in s.series_description.upper():
							info[mips_sag].append({'item': s.series_id})
				else:
					if postop:
						info[swi_gre_elec].append({'item': s.series_id})
					else:
						info[swi_gre].append({'item': s.series_id})

			elif any(substring.upper() in s.series_description.upper() for substring in {'DWI', 'DTI', 'DIFFUSION','DWI-DTI'}):
				if len(is_diffusion_derived(s.series_description))>=1:
					str_derv=diff_dic[is_diffusion_derived(s.series_description)[0].upper()]
					info[dwi_acq].append({'item': s.series_id,'acq': str_derv if not postop else 'Electrode'+str_derv})
				else:
					if postop:
						info[dwi_acq].append({'item': s.series_id, 'acq': 'Electrode'})
					else:
						info[dwi].append({'item': s.series_id})

			elif any(substring in s.series_description.upper() for substring in {'AX', 'COR','SAG'}) and ('3D' not in s.series_description.upper()):
				if ('AX' in s.series_description.upper()):
					orientation = 'Tra'
				elif ('COR' in s.series_description.upper()):
					orientation = 'Cor'
				elif ('SAG' in s.series_description.upper()):
					orientation = 'Sag'
				
				if postop:
					if any(substring in s.series_description.upper() for substring in {'T2', '2D'}):
						info[t2w].append({'item': s.series_id, 'acq': 'Electrode' + orientation})
					elif any(substring in s.series_description.upper() for substring in {'PD'}):
						info[t1w_pd].append({'item': s.series_id, 'acq': 'Electrode' + orientation})
					elif any(substring in s.series_description.upper() for substring in {'FLAIR'}):
						info[t1w_flair].append({'item': s.series_id, 'acq': 'Electrode' + orientation})
				else:
					if any(substring in s.series_description.upper() for substring in {'T2','2D'}):
						info[t2w].append({'item': s.series_id, 'acq': orientation})
					elif any(substring in s.series_description.upper() for substring in {'PD'}):
						info[t1w_pd].append({'item': s.series_id, 'acq': orientation})
					elif any(substring in s.series_description.upper() for substring in {'FLAIR'}):
						info[t1w_flair].append({'item': s.series_id, 'acq': orientation})
					elif any(substring in s.series_description.upper() for substring in {'SSFSE'}):
						info[t1w_acq].append({'item': s.series_id, 'acq': 'SSFSE' + orientation})

			elif any(substring in s.series_description.upper() for substring in {'FRACTIONAL', 'ANSIO', 'ANSIO.'}):
				if postop:
					info[fa_acq].append({'item': s.series_id, 'acq': 'Electrode'})
				else:
					info[fa].append({'item': s.series_id})
		
		
		elif any(substring in s.study_description.upper() for substring in {'CT'}) and all(x not in s.series_description.upper() for x in ('SCOUT','SUMMARY')):
			electrode_list = {'OVER', 'UNDER', 'ELECTRODE', 'ROUTINE', 'F_U_HEAD', 'F/U_HEAD', 'ER_HEAD', 'POST OP','POSTOP'}
			frame_list = {'STEROTACTIC', 'STEREOTACTIC','STEREOTACTIC FRAME', 'STEALTH', 'CTA_COW','Axial 1.200 CE'}
			
			if any(substring in s.protocol_name.upper() for substring in electrode_list):
				if any(x.upper() in s.series_description.upper() for x in ('BONE','SEMAR')):
					info[ct_acq_desc].append({'item': s.series_id, 'acq': 'Electrode', 'desc':'BONE'})
				else:
					info[ct_acq].append({'item': s.series_id, 'acq': 'Electrode'})
			elif any(x.upper() in s.protocol_name.upper() for x in frame_list):
				if any(x.upper() in s.series_description.upper() for x in ('BONE','SEMAR')):
					info[ct_acq_desc].append({'item': s.series_id, 'acq': 'Frame', 'desc':'BONE'})
				else:
					info[ct_acq].append({'item': s.series_id, 'acq': 'Frame'})
				
		elif any(substring in s.study_description.upper() for substring in {'PET'}) or any(substring in s.series_description.upper() for substring in {'PET CORR'}):
			if any(substring in s.series_description.upper() for substring in {'ITERATIVE','RECON','FBP','MIP','PET CORR'}):
				if '3D_FBP' not in s.series_description.upper():
					info[pet].append({'item': s.series_id})
				
	return info
