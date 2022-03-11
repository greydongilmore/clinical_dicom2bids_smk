def get_pre_t1_filename(wildcards):
    if config['noncontrast_t1']['acq']:
        files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='*', suffix='T1w.nii.gz'))

        if len(files) <=1:
            file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
            if len(file)==0:
                file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
        else:
            files.sort(key=lambda f: int(re.sub('\D', '', f)))
            file=files[-1]
    else:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
    
    print(f'Pre T1w non-contrast file: {basename(file[0])}')
    return file

if config['fastsurfer']['seg_only']:
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            threads = config['fastsurfer']['threads'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
        output:
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id)),
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        #threads:config['fastsurfer']['threads']
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 {params.fastsurfer_run}/run_fastsurfer.sh \
--t1 {input.t1} --sd {output.fastsurfer_out} --sid {params.sid} --order {params.order} --py {params.py} --run_viewagg_on cpu"
else:
    rule fastsurfer_seg:
        input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run =config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            threads = config['fastsurfer']['threads'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
        output:
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', config['fastsurfer']['sid'], 'sub-' + subject_id)),
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        threads:config['fastsurfer']['threads']
        shell:
            "export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 &&export FASTSURFER_HOME={params.fastsurfer_run} &&{params.fastsurfer_run}/run_fastsurfer.sh \
--t1 {input.t1} --sd {output.fastsurfer_out} --sid {params.sid} --order {params.order} --py {params.py} --threads {params.threads} --batch {params.batch} --run_viewagg_on cpu"

final_outputs.extend(expand(rules.fastsurfer_seg.output.touch_fastsurfer, subject=subjects))
