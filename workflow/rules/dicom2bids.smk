def get_dicom_dir(wildcards):
	if config['anonymize']:
		for root, folders, files in walk(path.join(config['dicom_dir'],'sub-' + wildcards.subject)):
			for file in files:
				fileN = '_'.join([path.basename(root),file+'.dcm']) if not file.endswith('.dcm') else '_'.join([path.basename(root),file])
				
				if not path.exists(path.join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + wildcards.subject)):
					makedirs(path.join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + wildcards.subject))

				copyfile(path.abspath(path.join(root,file)), path.join(config['out_dir'], 'sourcedata', 'dicoms', 'sub-' + wildcards.subject, fileN))

		return path.join(config['out_dir'], 'sourcedata', 'dicoms','sub-' + wildcards.subject)
	else:
		return path.join(config['dicom_dir'],'sub-' + wildcards.subject)

rule dicom2tar:
	input:
		dicom = get_dicom_dir
	output:
		tar = directory(path.join(config['out_dir'], 'sourcedata', 'tars', subject_id))
	params:
		clinical_events=config['clinical_event_file'],
		log_dir=path.join(config['out_dir'],'logs'),
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script:
		"../scripts/dicom2tar/main.py"

rule tar2bids:
	input:
		tar = path.join(config['out_dir'], 'sourcedata', 'tars', subject_id),
	params:
		heuristic_file = config['heuristic'],
		bids = directory(path.join(config['out_dir'], 'bids_tmp')),
		dcm_config=config['dcm_config']
	output:
		touch_tar2bids=touch(path.join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done")),
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	shell:
		'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix --dcmconfig {params.dcm_config} -b'

rule cleanSessions:
	input:
		touch_tar2bids=path.join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done"),
	output:
		touch_dicom2bids=touch(path.join(config['out_dir'], 'logs', 'sub-' + subject_id + "_dicom2bids.done")),
	params:
		clinical_events=config['clinical_event_file'],
		bids_fold = path.join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id),
		num_subs = len(subjects),
		ses_calc = config['session_calc'],
		sub_group = config['sub_group']
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script:
		"../scripts/post_tar2bids/clean_sessions.py"