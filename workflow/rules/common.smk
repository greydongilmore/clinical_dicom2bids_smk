def bids(root=None, datatype=None, prefix=None, suffix=None, subject=None, session=None,include_subject_dir=True,include_session_dir=True,**entities):
    """Helper function for generating bids paths for snakemake workflows

    File path is of the form:

    [root]/[sub-{subject}]/[ses-{session]/[prefix]_[sub-{subject}]_[ses-{session}]_[{key}-{val}_ ... ]_[suffix]

    root -- root folder to include in the path (e.g. 'results'))
    datatype -- folder to include after sub-/ses- (e.g. anat, dwi )
    prefix -- string to prepend to the file name (typically not defined, unless you want tpl-{tpl}, or a datatype)
    suffix -- bids suffix including extension (e.g. 'T1w.nii.gz')
    subject -- subject to use, for folder and filename
    session -- session to use, for folder and filename
    include_subject_dir -- whether to include the sub-{subject} folder if subject defined (default: True)
    include_session_dir -- whether to include the ses-{session} folder if session defined (default: True)
    **entities -- dictionary of bids entities (e.g. space=T1w for space-T1w)

    Returns: bids-like file path

    Example:

        Below is a rule using bids naming for input and output:

        rule proc_img:
            input: 'sub-{subject}_T1w.nii.gz'
            output: 'sub-{subject}_space-snsx32_desc-preproc_T1w.nii.gz'

        With bids() you can instead use:

         rule proc_img:
            input: bids(subject='{subject}',suffix='T1w.nii.gz')
            output: bids(subject='{subject}',space='snsx32',desc='preproc',suffix='T1w.nii.gz')

        Note that here we are not actually using "functions as inputs" in snakemake, which would require
        a function definition with wildcards as the argument, and restrict to input/params, but bids()
        is being used simply to return a string.

        Also note that space, desc and suffix are NOT wildcards here, only {subject} is.
        This makes it easy to combine wildcards and non-wildcards with bids-like naming.

        However, you can still use bids() in a lambda function, this is especially useful if your wildcards
        are named the same as bids entities (e.g. {subject}, {session}, {task} etc..):
 
        rule proc_img:
            input: lambda wildcards: bids(**wildcards,suffix='T1w.nii.gz')
            output: bids(subject='{subject}',space='snsx32',desc='preproc',suffix='T1w.nii.gz')

        Or another example where you may have many bids-like wildcards used in your workflow:
        
        rule denoise_func:
            input: lambda wildcards: bids(**wildcards, suffix='bold.nii.gz')
            output: bids(subject='{subject}',session='{session}',task='{task}',acq='{acq}',desc='denoise',suffix='bold.nii.gz')

        In this example, all the wildcards will be determined from the output and passed on to bids() for inputs.
        The output filename will have a 'desc-denoise' flag added to it.


        Also note that even if you supply entities in a different order, the
        entities will be ordered based on the OrderedDict defined here.
        If you entities not known are provided, they will be just be placed
        at the end (before the suffix), the the order provided.

        Note: For maximum flexibility all arguments are optional (if none are specified, will return empty string)

    -- some code adapted from mne-bids
      https://mne.tools/mne-bids/stable/_modules/mne_bids/utils.html


    """

    from collections import OrderedDict
    from os.path import join
    

    #replace underscores in keys (needed to that users can use reserved keywords by appending a _)
    entities = { k.replace('_', ''): v for k, v in entities.items() }
       
 
   
    #strict ordering of bids entities is specified here:
    order = OrderedDict([('task', None),
                         ('acq', None),
                         ('ce', None),
                         ('rec', None),
                         ('dir', None),
                         ('run', None),
                         ('mod', None),
                         ('echo', None),
                         ('hemi', None),
                         ('space', None),
                         ('res', None),
                         ('den', None),
                         ('label', None),
                         ('desc', None)])

    #now add in entities (this preserves ordering above)
    for key, val in entities.items():
        order[key] = val

    #initialize lists for filename and folder
    # will append to these, then '_'.join() os.path.join() respectively
    filename = []
    folder = []

    #root directory
    if isinstance(root,str):
        folder.append(root)

    #if prefix is defined, put it before other anything else
    if isinstance(prefix, str):
        filename.append(prefix)

    #if subject defined then append to file and folder
    if isinstance(subject,str):
        if include_subject_dir is True:
            folder.append(f'sub-{subject}')
        filename.append(f'sub-{subject}')

    #if session defined then append to file and folder
    if isinstance(session,str):
        if include_session_dir is True:
            folder.append(f'ses-{session}')
        filename.append(f'ses-{session}')
    
    if isinstance(datatype,str):
        folder.append(datatype)
    
    #add the entities
    for key, val in order.items():
        if val is not None:
            filename.append(f'{key}-{val}')

    #if suffix is defined, append it
    if isinstance(suffix, str):
        filename.append(suffix)


    if len(filename) == 0:    
        return ''

    #now, join up the lists:
    filename = '_'.join(filename)

    if len(folder)>0:
        filename = join(*folder,filename)
    
    return filename


