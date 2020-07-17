rule cleanSessions:
	input:
		bids_file = join(out_dir, 'logs', 'sub-P'+'{subject}'+ "_bids_complete.txt"),
	output:
		out_file = join(out_dir, 'logs', 'sub-P'+'{subject}'+ "_complete.txt")
	params:
		script = "workflow/scripts/post_tar2bids/clean_sessions.py",
		final_script = "workflow/scripts/process_complete.py",
		bids_fold = join(out_dir, 'bids', 'sub-P' + '{subject}'),
	shell:
		'python3 {params.script} --output_dir {out_dir} --bids_fold {params.bids_fold};python3 {params.final_script} --out_file {output.out_file}'