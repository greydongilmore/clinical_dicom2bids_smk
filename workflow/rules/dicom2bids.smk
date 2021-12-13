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
		touch_dicom2bids=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_cleanSessions.done")),
	params:
		clinical_events=config['clinical_event_file'],
		bids_fold = join(config['out_dir'], 'bids_tmp', 'sub-' + subject_id),
		num_subs = len(subjects),
		ses_calc = config['session_calc'],
		sub_group = config['sub_group']
	#container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script:
		"../scripts/post_tar2bids/clean_sessions.py"

final_outputs.extend(expand(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_tar2bids.done"), subject=subjects))
final_outputs.extend(expand(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_cleanSessions.done"), subject=subjects))
