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

rule fmriprep_seg:
    input: 
        t1 = get_pre_t1_filename,
    params:
        bids_dir = join(config['out_dir'], 'bids'),
        out_dir = join(config['out_dir'], 'derivatives', 'fmriprep'),
        license = config['fmriprep']['fmriprep_license'],
        license_name = basename(config['fmriprep']['fmriprep_license']),
        sub = subject_id,
        fmriprep_img=join(dirname(workflow.basedir), config['singularity']['fmriprep']),
        bids_filter=config['fmriprep']['bids_filter'],
    output:
        touch_fmriprep=touch(join(config['out_dir'], 'logs', 'sub-' + subject_id + "_fmriprep.done"))
    container: config["singularity"]["fmriprep"]
    shell:
        'export SINGULARITYENV_FS_LICENSE=$HOME/.freesurfer.txt&&\
singularity run --cleanenv \
--bind {params.bids_dir}:/tmp/input \
--bind {params.out_dir}:/tmp/output \
--bind {params.license}:/tmp/{params.license_name} \
{params.fmriprep_img} /tmp/input  /tmp/output participant --skip_bids_validation \
--participant_label {params.sub} --anat-only \
--fs-license-file /tmp/{params.license_name} \
--bids-filter-file {params.bids_filter}'

final_outputs.extend(expand(rules.fmriprep_seg.output, subject=subjects))
