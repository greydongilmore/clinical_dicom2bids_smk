rule dicom2tar:
	input:
		dicom = join(dicom_dir,'sub-' + '{subject}')
	output:
		tar = directory(join(out_dir, 'tars', 'P' +'{subject}'))
	params:
		script = "workflow/scripts/dicom2tar/main.py"
	shell:
		'python3 {params.script} {input.dicom} {output.tar} --clinical_scans'