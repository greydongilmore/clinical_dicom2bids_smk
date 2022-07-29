
def get_pre_t1_filename(wildcards):
    if config['noncontrast_t1']['acq']:
        files=glob(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+f'{wildcards.subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='*', suffix='T1w.nii.gz'))
        if len(files) <=1:
            file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
            if not exists(file[0]):
                file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', acq=config['noncontrast_t1']['acq'], run='02', suffix='T1w.nii.gz'),subject=wildcards.subject)
                if not exists(file[0]):
                    file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
        else:
            files.sort(key=lambda f: int(re.sub('\D', '', f)))
            file=files[-1]
    else:
        file=expand(bids(root=join(config['out_dir'], 'bids'), subject=config['subject_prefix']+'{subject}', datatype='anat', session='pre', run='01', suffix='T1w.nii.gz'),subject=wildcards.subject)
    
    print(f'Pre T1w non-contrast file: {basename(file[0])}')
    return file

if config['fastsurfer']['seg_only']:
    rule fastsurfer_seg_only:
         input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            threads = config['fastsurfer']['threads'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', 'fastsurfer')),
            subjid=expand('sub-' + subject_id,subject=subjects),
        output:
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        #threads:config['fastsurfer']['threads']
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 {params.fastsurfer_run}/run_fastsurfer.sh \
--t1 {input.t1} --sd {params.fastsurfer_out} --sid {params.subjid} --order {params.order} --py {params.py} --run_viewagg_on cpu --fsaparc --parallel --surfreg"

    final_outputs.extend(expand(rules.fastsurfer_seg_only.output.touch_fastsurfer, subject=subjects))
else:
    rule fastsurfer_all:
         input: 
            t1 = get_pre_t1_filename,
        params:
            fastsurfer_run = config['fastsurfer']['home'],
            sid = config['fastsurfer']['sid'],
            batch = config['fastsurfer']['batch'],
            threads = config['fastsurfer']['threads'],
            order = config['fastsurfer']['order'],
            py = config['fastsurfer']['py'],
            fastsurfer_out = directory(join(config['out_dir'], 'derivatives', 'fastsurfer')),
            subjid=expand('sub-' + subject_id,subject=subjects),
        output:
            touch_fastsurfer=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fastsurfer.done")),
        #threads:config['fastsurfer']['threads']
        shell:
            "export FASTSURFER_HOME={params.fastsurfer_run} &&PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:4096 {params.fastsurfer_run}/run_fastsurfer.sh \
--t1 {input.t1} --sd {params.fastsurfer_out} --sid {params.subjid} --order {params.order} --py {params.py} --run_viewagg_on cpu --fsaparc --parallel --surfreg"

    final_outputs.extend(expand(rules.fastsurfer_all.output.touch_fastsurfer, subject=subjects))


