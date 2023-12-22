def get_dicom_dir(wildcards):
	subji = wildcards.subject.strip(config['subject_prefix'])
	if config['anonymize']:
		for root, folders, files in walk(join(config['dicom_dir'],'sub-' + subji)):
			for file in files:
				fileN = '_'.join([basename(root),file+'.dcm']) if not file.endswith('.dcm') else '_'.join([basename(root),file])
				
				if not exists(join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + subji)):
					makedirs(join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + subji))

				copyfile(abspath(join(root,file)), join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + subji, fileN))

		return join(config['out_dir'], 'sourcedata', 'dicoms','sub-' + subji)
	else:
		return join(config['dicom_dir'], 'sub-' + subji)

rule dicom2tar:
	input:
		dicom = get_dicom_dir
	output:
		tar = directory(join(config['out_dir'], 'sourcedata', 'tars', subject_id))
	params:
		clinical_events=config['clinical_event_file'],
		log_dir=join(config['out_dir'],'logs'),
		prefix=config['subject_prefix'],
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script:
		"../scripts/dicom2tar/main.py"

rule tar2bids:
	input:
		tar = join(config['out_dir'], 'sourcedata', 'tars', subject_id),
	params:
		heuristic_file = config['heuristic'],
		bids = directory(join(config['out_dir'], 'bids_tmp')),
		dcm_config=config['dcm_config']
	output:
		tmp_dir = directory(join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id))
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	shell:
		'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix --dcmconfig {params.dcm_config} -b'

if config['fastsurfer']['run'] or config['fmriprep']['run']:
	rule cleanSessions:
		input:
			#touch_tar2bids=join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done"),
			tmp_dir=rules.tar2bids.output.tmp_dir,
		output:
			#touch_dicom2bids=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_cleanSessions.done")),
			tmp_dir = directory(join(config['out_dir'], 'bids', 'sub-' + subject_id)),
			t1w_file= bids(root=join(config['out_dir'], 'bids'), subject=subject_id, datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),
		params:
			clinical_events=config['clinical_event_file'],
			bids_fold = join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id),
			num_subs = len(subjects),
			ses_calc = config['session_calc'],
			sub_group = config['sub_group']
		#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
		script:
			"../scripts/post_tar2bids/clean_sessions.py"
else:
	rule cleanSessions:
		input:
			#touch_tar2bids=join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done"),
			tmp_dir=rules.tar2bids.output.tmp_dir,
		output:
			tmp_dir = directory(join(config['out_dir'], 'bids', 'sub-' + subject_id)),
		params:
			clinical_events=config['clinical_event_file'],
			bids_fold = join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id),
			num_subs = len(subjects),
			ses_calc = config['session_calc'],
			sub_group = config['sub_group']
		#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
		script:
			"../scripts/post_tar2bids/clean_sessions.py"

final_outputs.extend(expand(join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id), subject=subjects))
final_outputs.extend(expand(join(config['out_dir'], 'bids', 'sub-' + subject_id), subject=subjects))

if config['fastsurfer']['run'] or config['fmriprep']['run']:
	final_outputs.extend(expand(bids(root=join(config['out_dir'], 'bids'), subject=subject_id, datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'), subject=subjects))
