rule dicom2tar:
	input:
		dicom = join(dicom_dir,'sub-' + '{subject}')
	output:
		tar = directory(join(out_dir, 'tars', 'P' +'{subject}'))
	params:
		config['or_dates_file']
	container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script: 
		"../scripts/dicom2tar/main.py"

rule tar2bids:
	input:
		tar = join(out_dir, 'tars', 'P' + '{subject}'),
	params:
		heuristic_file = heuristic_file,
		bids = directory(join(out_dir, 'bids_tmp')),
		dcm_config=config['dcm_config']
	output:
		touch(join(out_dir, 'sub-P' + '{subject}' + "_tar2bids.done")),
	container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	shell:
		'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix --dcmconfig {params.dcm_config} -b'

rule cleanSessions:
	input:
		join(out_dir, 'sub-P' + '{subject}' + "_tar2bids.done"),
	output:
		touch(join(out_dir, 'sub-P' + '{subject}' + "_dicom2bids.done")),
	params:
		or_dates_file = config['or_dates_file'],
		bids_fold = join(out_dir, 'bids_tmp', 'sub-P' + '{subject}'),
		num_subjects = len(subjects),
		ses_calcs = config['session_calc'],
	container: 'docker://greydongilmore/dicom2bids-clinical:latest'
	script:
		"../scripts/post_tar2bids/clean_sessions.py"