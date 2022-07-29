# Variables
hemi = ["lh", "rh"]
surf_suffix = ["pial", "white", "inflated", "thickness"]

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

rule fs_surf_to_gii:
    """ Convert surfaces (pial, white, etc.) to gifti (surf.gii) """
    input: join(config['out_dir'], 'derivatives', config['fastsurfer']['sid']) + ''.join(['/sub-',subject_id])+ "/surf/{hemi}.{surf_suffix}"
    params:
        struct = lambda wildcards: "CORTEX_LEFT" if wildcards.hemi == "lh" else "CORTEX_RIGHT"
    output: join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + "/surf/{hemi,(lh|rh)}.{surf_suffix,(pial|white|inflated|sphere.reg)}.surf.gii"
    shell:
        "mris_convert {input} {output} && "
        "wb_command -set-structure {output} {params.struct}"

rule compute_cortical_thickness:
    """ Compute thickness from resampled surface """
    input: 
        pial = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/surf/{{hemi}}.pial.surf.gii",
        white = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/surf/{{hemi}}.white.surf.gii",
    output:
        gii = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/metric/{{hemi}}.thickness.shape.gii",
    shell:
        "wb_command -surface-to-surface-3d-distance {input.pial} {input.white} {output}"    

rule gen_depth_surfaces:
    """ Generate surfaces for different depths """
    input:         
        pial = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) +f"/surf/{{hemi}}.pial.surf.gii",
        white = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/surf/{{hemi}}.white.surf.gii",
    output:
        depth = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/surf/{{hemi,(lh|rh)}}.surf.gii"
    shell:
        "wb_command -surface-cortex-layer {input.white} {input.pial} {output.depth}"

rule sample_depth_surfaces:
    """ Sample values at different depths """
    input:
        t1 = get_pre_t1_filename,
        depth = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/surf/{{hemi}}.surf.gii",
    params:
        sample_method = "trilinear"
    output:
        depth = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id]) + f"/metric/{{hemi}}.T1w.shape.gii",
    shell:
        "wb_command -volume-to-surface-mapping {input.t1} {input.depth} {output.depth} -{params.sample_method}"

rule qc_surf:
    """
    Create visualization to QC generated white and pial surfaces
    """
    input:
        scene_template = "resources/surf_qc_template.scene",
        lh_pial = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id])+ "/surf/lh.pial.surf.gii",
        lh_white = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id])+ "/surf/lh.white.surf.gii",
        rh_pial = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id])+ "/surf/rh.pial.surf.gii",
        rh_white = join(config['out_dir'], 'derivatives', "gifti") + ''.join(['/sub-',subject_id])+ "/surf/rh.white.surf.gii",
        t1 = get_pre_t1_filename
    output: 
        report = join(config['out_dir'], 'derivatives', 'gifti') + ''.join(['/sub-',subject_id]) + '/qc/sub-' + subject_id +'_surfqc.svg',
    script: "../scripts/viz_surfqc.py"

final_outputs.extend(expand(join(config['out_dir'], 'derivatives', 'gifti') + ''.join(['/sub-',subject_id]) + '/qc/sub-' + subject_id +'_surfqc.svg',subject=subjects))