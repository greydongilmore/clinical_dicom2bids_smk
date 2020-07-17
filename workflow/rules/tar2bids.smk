rule tar2bids:
	input:
		tar = join(out_dir, 'tars', 'P' + '{subject}'),
	params:
		heuristic_file = heuristic_file,
		bids = directory(join(out_dir, 'bids')),
		script = "workflow/scripts/process_complete.py"
	output:
		#bids = directory(join(out_dir, 'bids')),
		out_file = join(out_dir, 'logs', 'sub-P'+'{subject}'+ "_bids_complete.txt")
	shell:
	#	'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix -b;'
		'heudiconv --files {input.tar} -o {params.bids} -f {params.heuristic_file} -c dcm2niix -b;python3 {params.script} --out_file {output.out_file}'
