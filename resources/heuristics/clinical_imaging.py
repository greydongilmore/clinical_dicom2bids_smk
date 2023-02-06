# -*- coding: utf-8 -*-
import regex as re

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

orien_dic={
	'ax':'Tra',
	'cor':'Cor',
	'sag':'Sag',
	'2d':'Tra',
	'3d':'3D',
	'spc':'SPACE',
	'space':'SPACE',
}

diff_dic={
	'ADC':'ADC',
	'TRACEW':'TRACEW',
	'WEIGHTED TRACE':'TRACEW',
	'COLFA':'COLFA',
	'FA':'FA',
	'TENSOR':'TENSOR',
	'TENSOR_B0':'TENSOR',
}

# Utility:  Check a list of regexes for truthyness
def regex_search_label(regexes, label):
	if any(regex.search(label) for regex in regexes):
		out_list=[regex.search(label)[0] for regex in regexes if regex.search(label) is not None]
		return [x for x in out_list if x !='']
	else:
		return []

def is_2d(label):
	regexes = [
		re.compile('(?:\s+|^)ax[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)cor[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)sag[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)2d[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_3d(label):
	regexes = [
		re.compile('(?:\s+|^)3d[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)spc[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)space[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return [x.strip() for x in regex_search_label(regexes, label)]

def is_postop(label):
	regexes = [
		re.compile('stimulator', re.IGNORECASE),
		re.compile("[stim\s]*safe", re.IGNORECASE),
		re.compile('post([ ]{0,1}[\s-]?)op', re.IGNORECASE),
		re.compile('sar', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_post_ignore(label):
	regexes = [
		re.compile('post stroke', re.IGNORECASE),
		re.compile('post grappa', re.IGNORECASE),
		re.compile('post gad', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_stealth(label):
	regexes = [
		re.compile('stealth', re.IGNORECASE),
		re.compile('[stereo\s-include]*', re.IGNORECASE),
		re.compile('1.5 mm anatomy', re.IGNORECASE),
		re.compile("ax[\s]3d[\s]mprage", re.IGNORECASE),
		re.compile("axial[\s]mprage", re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_post_contrast(label):
	regexes = [
		re.compile('gad', re.IGNORECASE),
		re.compile('post[\s]gad', re.IGNORECASE),
		re.compile('contrast', re.IGNORECASE),
		re.compile('(?:\s+|^)+c[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_fspgr(label):
	regexes = [
		re.compile('[IR\s-_]*fspgr', re.IGNORECASE),
		re.compile('axial ir', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_gradient_echo(label):
	regexes = [
		re.compile('(?:\s+|^)mpgr[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)mprage[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

# Proton Density
def is_proton_density(label):
	regexes = [
		re.compile('(?:\s+|^)pd[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('pd[\s]tse', re.IGNORECASE),
		re.compile('(?=.*proton)(?=.*density)', re.IGNORECASE)
		]
	return regex_search_label(regexes, label)

def is_flair(label):
	regexes = [
		re.compile('(?:\s+|^)flair[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_t2(label):
	regexes = [
		re.compile('(?:\s+|^)t2[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('t2[\s]tse', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_epi(label):
	regexes = [
		re.compile('(?:\s+|^)epi[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)

def is_fa(label):
	regexes = [
		re.compile('(?:\s+|^)FRACTIONAL[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		re.compile('(?:\s+|^)ANSIO[!,?."\"]*(?=\s+|$)', re.IGNORECASE),
		]
	return regex_search_label(regexes, label)



# Diffusion
def is_diffusion(label):
    regexes = [
        re.compile('dti', re.IGNORECASE),
        re.compile('dwi', re.IGNORECASE),
        re.compile('diff_', re.IGNORECASE),
        re.compile('diffusion', re.IGNORECASE),
        re.compile('(?=.*diff)(?=.*dir)', re.IGNORECASE),
        re.compile('dwi[\s-\s]dti', re.IGNORECASE),
        re.compile('t2-weighted trace', re.IGNORECASE)
        ]
    return regex_search_label(regexes, label)

# Diffusion - Derived
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
	t2w = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_run-{item:02d}_T2w')
	t2w_acq = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T2w')
	t1w_pd = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_PD')
	flair = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_run-{item:02d}_FLAIR')
	flair_acq = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_FLAIR')
	fa = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_run-{item:02d}_angio')
	fa_acq = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_angio')
	
	#fmap
	fmap = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_run-{item:02d}_epi')
	fmap_acq = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_epi')

	#Diffusion
	dwi = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_run-{item:02d}_dwi')
	dwi_acq = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_dwi')# fractional anisotropy,trace-weighted image,apparent diffusion coefficient,color fractional anisotropy

	#CT
	ct = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_run-{item:02d}_ct')
	ct_acq = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_ct')
	ct_acq_desc = create_key('{bids_subject_session_dir}/ct/{bids_subject_session_prefix}_acq-{acq}_desc-{desc}_run-{item:02d}_ct')
	
	#pet
	pet = create_key('{bids_subject_session_dir}/pet/{bids_subject_session_prefix}_task-rest_run-{item:02d}_pet')
	pet_acq = create_key('{bids_subject_session_dir}/pet/{bids_subject_session_prefix}_task-rest_acq-{acq}_run-{item:02d}_pet')
	pet_task = create_key('{bids_subject_session_dir}/pet/{bids_subject_session_prefix}_task-{task}_acq-{acq}_run-{item:02d}_pet')

	info = {t1w:[],
			t1w_acq:[],
			t2w:[],
			t2w_acq:[],
			t1w_pd:[],
			flair:[],
			flair_acq:[],
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
			pet_task:[]}
	
	for idx, s in enumerate(seqinfo):
		if any(substring in s.study_description.upper() for substring in {'MR'}):
			
			tmp={'item': s.series_id}

			postop=[]
			if len(is_post_ignore(s.series_description))==0:
				postop=is_postop(s.series_description)
			
			orientation=is_2d(s.series_description)
			dimensions=is_3d(s.series_description)

			acq_label=[]
			if len(orientation) > 0:
				acq_label=orien_dic[orientation[0].lower()]

			if len(dimensions)>0:
				acq_label = (len(dimensions)==1 and dimensions[0]) or (any('space' in x.lower() for x in dimensions) and 'space'.upper()) or dimensions[0]

			if len(is_stealth(s.series_description))>=1 and 'fspgr' not in s.series_description.lower() or len(is_post_contrast(s.series_description))>=1  and len(orientation) > 0:
				
				if 'mpgr' in s.series_description.lower():
					if len(postop) >=1:
						info[t1w_acq].append({**tmp, 'acq': 'ElectrodeMPGR3D'})
					else:
						info[t1w].append(tmp)

				elif len(is_flair(s.series_description))>=1:
					if len(postop) >=1:
						info[flair_acq].append({**tmp, 'acq': 'ElectrodeFLAIR'})
					else:
						info[flair].append(tmp)

				elif 't2' in s.series_description.lower():
					if len(postop) >=1:
						info[t2w_acq].append({**tmp,'acq': 'Electrode'})
					else:
						info[t2w].append(tmp)
				else:
					if len(postop) >=1:
						info[t1w_acq].append({**tmp, 'acq': 'Electrode3D'})
					else:
						info[t1w].append(tmp)

			elif len(is_fspgr(s.series_description))>=1  and len(orientation) >0:
				if len(postop) >=1:
					info[t1w_acq].append({**tmp,'acq': 'Electrode'})
				else:
					info[t1w].append(tmp)

			elif len(is_epi(s.series_description))>=1 and len(orientation) > 0:
				if len(postop) >=1:
					info[fmap_acq].append({**tmp, 'acq': 'Electrode'})
				else:
					info[fmap].append(tmp)
					
			elif len(is_gradient_echo(s.series_description))>=1  and len(orientation) > 0:
				if len(postop) >=1:
					info[t1w_acq].append({**tmp, 'acq': 'Electrode'})
				else:
					info[t1w].append(tmp)

			if len(is_fa(s.series_description))>=1  and len(orientation) > 0:
				if len(postop) >=1:
					info[fa_acq].append({**tmp, 'acq': 'Electrode'})
				else:
					info[fa].append(tmp)

			elif len(is_diffusion(s.series_description))>=1:
				if len(is_diffusion_derived(s.series_description))>=1:
					str_derv=diff_dic[is_diffusion_derived(s.series_description)[0].upper()]
					info[dwi_acq].append({**tmp,'acq': str_derv if len(postop) >=1 else 'Electrode'+str_derv})
				else:
					info[dwi].append(tmp)
			
			elif len(is_t2(s.series_description))>=1 and len(is_proton_density(s.series_description))==0 and len(is_flair(s.series_description))==0:
				if len(postop) >=1:
					info[t2w_acq].append({**tmp, 'acq': 'Electrode'})
				elif len(acq_label) >=1:
					info[t2w_acq].append({**tmp, 'acq': acq_label})
				else:
					info[t2w].append(tmp)
			
			elif len(is_proton_density(s.series_description))>=1:
				if len(postop) >=1:
					info[t1w_pd].append({**tmp, 'acq': 'Electrode'})
				elif len(acq_label) >=1:
					info[t1w_pd].append({**tmp, 'acq': acq_label})
				else:
					info[t1w_pd].append({**tmp, 'acq': 'PD'})

			elif len(is_flair(s.series_description))>=1:
				if len(postop) >=1:
					info[flair_acq].append({**tmp, 'acq': 'Electrode'})
				elif len(acq_label) >=1:
					info[flair_acq].append({**tmp, 'acq': acq_label})
				else:
					info[flair].append(tmp)

			elif 'ssfse' in  s.series_description.lower()[0] !='':
				if len(acq_label) >=1:
					info[t1w_acq].append({**tmp, 'acq': 'SSFSE'+acq_label})
				else:
					info[t1w_acq].append({**tmp, 'acq': 'SSFSE'})

		elif any(substring in s.study_description.upper() for substring in {'PET'}) or any(substring in s.series_description.upper() for substring in {'PET CORR'}):
			if any(substring in s.series_description.upper() for substring in {'ITERATIVE','RECON','FBP','MIP','PET CORR'}):
				if '3D_FBP' not in s.series_description.upper():
					info[pet].append({'item': s.series_id})
		else:
			#---CT SCANS
			ct_scan = False
			if s.study_description =='':
				if any(x.upper() in s.series_description.upper() for x in {'NO ANGLE'}):
					ct_scan = True
			elif any(x.upper() in s.study_description.upper() for x in {'CT'}):
				ct_scan = True
			elif any(x.upper() in s.series_description.upper() for x in {'H31S','PF HEAD','H37S','HEAD','AXIAL'}):
				ct_scan = True
			
			if ct_scan:
				electrode_list = {'OVER', 'UNDER', 'ELECTRODE', 'F_U_HEAD', 'F/U_HEAD', 'ER_HEAD', 'POST', 'POST OP'}
				frame_list = {'STEROTACTIC', 'STEREOTACTIC','FRAME', 'FRAME', 'STEALTH', 'CTA_COW','Axial 1.200 CE'}
				
				if all(x not in s.series_description.upper() for x in ('SCOUT','SUMMARY')):
					if any(x.upper() in s.protocol_name.upper() for x in electrode_list):
						if any(x.upper() in s.series_description.upper() for x in ('BONE','SEMAR')):
							info[ct_acq_desc].append({'item': s.series_id, 'acq': 'Electrode', 'desc':'BONE'})
						else:
							info[ct_acq].append({'item': s.series_id, 'acq': 'Electrode'})
					elif any(x.upper() in s.protocol_name.upper() for x in frame_list):
						if any(x.upper() in s.series_description.upper() for x in ('BONE','SEMAR')):
							info[ct_acq_desc].append({'item': s.series_id, 'acq': 'Frame', 'desc':'BONE'})
						else:
							info[ct_acq].append({'item': s.series_id, 'acq': 'Frame'})
					else:
						info[ct].append({'item': s.series_id})
				
	return info
