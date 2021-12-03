def get_pre_t1_filename(wildcards):
    files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='anat', session='pre', run='*', suffix='T1w.nii.gz'))
    if len(files) <=1:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
    else:
        file=files[-1]
    return file

if config['fastsurfer']['seg_only']:
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = join(config['fastsurfer']['home'], 'run_fastsurfer.sh'),
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
        output:
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id)),
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        shell:
            '{params.fastsurfer_run} --t1 {input.t1} --sd {output.fastsurfer_out} --sid {params.sid} --seg_with_cc_only --order {params.order} \
            --py {params.py} --parallel'
else:
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = join(config['fastsurfer']['home'], 'run_fastsurfer.sh'),
            sid = config['fastsurfer']['sid'],
            threads = config['fastsurfer']['threads'],
            batch = config['fastsurfer']['batch'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
        output:
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id)),
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        shell:
            '{params.fastsurfer_run} --t1 {input.t1} --sd {output.fastsurfer_out} --sid {params.sid} --vol_segstats --order {params.order} \
            --batch {params.batch} --py {params.py} --parallel'

#rule fmriprep_seg:
#    input: 
#    	sub = get_pre_t1_filename,
#        t1 = bids(root=join(config['out_dir'], 'derivatives', 'atlasreg'), subject=subject_id, desc='n4', suffix='T1w.nii.gz')
#    params:
#        bids_dir = join(config['out_dir'], 'bids'),
#        out_dir = join(config['out_dir'], 'derivatives'),
#        license = config['fmriprep_license'],
#        sub = subject_id,
#        fmriprep=join(dirname(workflow.basedir),config['singularity']['fmriprep']),
#        bids_filter=config['fmriprep_bids_filter'],
#    output:
#        touch_fmriprep=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fmriprep.done"))
#    #container: config['singularity']['fmriprep']
#    shell:
#        'singularity exec --bind {params.bids_dir}:{params.bids_dir},{params.out_dir}:{params.out_dir} --writable-tmpfs {params.fmriprep} fmriprep --skip_bids_validation --bids-filter-file {params.bids_filter} {params.bids_dir} {params.out_dir} participant --participant-label {params.sub} --anat-only --fs-license-file {params.license}'
