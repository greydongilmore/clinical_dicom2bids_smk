def get_dicom_dir(wildcards):
	if config['anonymize']:
		for root, folders, files in walk(join(config['dicom_dir'],'sub-' + wildcards.subject)):
			for file in files:
				fileN = '_'.join([basename(root),file+'.dcm']) if not file.endswith('.dcm') else '_'.join([basename(root),file])
				
				if not exists(join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + wildcards.subject)):
					makedirs(join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + wildcards.subject))

				copyfile(abspath(join(root,file)), join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + wildcards.subject, fileN))

		return join(config['out_dir'], 'sourcedata', 'dicoms','sub-' + wildcards.subject)
	else:
		return join(config['dicom_dir'],'sub-' + wildcards.subject)

rule dicom2tar:
	input:
		dicom = get_dicom_dir
	output:
		tar = directory(join(config['out_dir'], 'sourcedata', 'tars', subject_id))
	params:
		clinical_events=config['clinical_event_file'],
		log_dir=join(config['out_dir'],'logs'),
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
		bids_fold = directory(join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id)),
		touch_tar2bids=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done")),
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	shell:
		'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix --dcmconfig {params.dcm_config} -b'

if config['sort_sessions']:
	rule cleanSessions:
		input:
			touch_tar2bids=join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done"),
			bids_fold = join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id),
		output:
			touch_cleanSessions=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_cleanSessions.done")),
			t1w_vol=bids(root=join(config['out_dir'], 'bids'), subject=subject_id, datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),
		params:
			clinical_events=config['clinical_event_file'],
			num_subs = len(subjects),
			ses_calc = config['session_calc'],
			sub_group = config['sub_group']
		#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
		script:
			"../scripts/post_tar2bids/clean_sessions.py"
else:
	rule cleanSessions:
		input:
			touch_tar2bids=join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done")
		output:
			touch_cleanSessions=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_cleanSessions.done")),
			t1w_vol=bids(root=join(config['out_dir'], 'bids'), subject=subject_id, datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),

final_outputs.extend(expand(rules.tar2bids.output.touch_tar2bids, subject=subjects))
final_outputs.extend(expand(rules.cleanSessions.output.touch_cleanSessions, subject=subjects))
