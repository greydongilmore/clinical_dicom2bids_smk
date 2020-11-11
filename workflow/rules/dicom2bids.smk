def get_dicom_dir(wildcards):
	if config['anonymize']:
		for root, folders, files in walk(join(config['dicom_dir'],'sub-' + wildcards.subject)):
			for file in files:
				filePath = abspath(join(root,file)).split('sub-' + wildcards.subject+os.sep)[-1]
				if not isdir(dirname(join(config['out_dir'], 'sourcedata', 'dicoms','sub-' + wildcards.subject,filePath))):
					makedirs(dirname(join(config['out_dir'], 'sourcedata', 'dicoms','sub-' + wildcards.subject,filePath)))
				copyfile(abspath(join(root,file)), join(config['out_dir'], 'sourcedata', 'dicoms','sub-' + wildcards.subject,filePath))

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
		touch_tar2bids=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done")),
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	shell:
		'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix --dcmconfig {params.dcm_config} -b'

rule cleanSessions:
	input:
		touch_tar2bids=join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done"),
	output:
		touch_dicom2bids=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_dicom2bids.done")),
	params:
		clinical_events=config['clinical_event_file'],
		bids_fold = join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id),
		num_subs = len(subjects),
		ses_calc = config['session_calc'],
		sub_group = config['sub_group']
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script:
		"../scripts/post_tar2bids/clean_sessions.py"