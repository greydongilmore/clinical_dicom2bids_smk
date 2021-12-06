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
            fastsurfer_run = config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
        output:
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id)),
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        threads:6
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&{params.fastsurfer_run}/run_fastsurfer.sh --t1 {input.t1} --sd {output.fastsurfer_out} --sid {params.sid} --order {params.order} \
            --py {params.py}"
else:
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run =config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            threads = config['fastsurfer']['threads'],
            batch = config['fastsurfer']['batch'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
        output:
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id)),
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        threads:6
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&{params.fastsurfer_run}/run_fastsurfer.sh --t1 {input.t1} --sd {output.fastsurfer_out} --sid {params.sid} --order {params.order} \
            --py {params.py} --parallel"

final_outputs.extend(expand(rules.fastsurfer_seg.output.touch_fastsurfer, subject=subjects))
