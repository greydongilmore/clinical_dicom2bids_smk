# Snakemake workflow: clinical DICOMs to BIDS

## Description

Snakemake workflow to convert a clinical dicom directory into BIDS structure.

## Requirements

* dcm2niix (v1.0.20200427)
* python requirements (defined in `workflow/envs/mapping.yaml`):
    * snakemake>=5.23.0
    * pydicom>=1.0.2
    * setuptools>=39.2.0
    * extractCMRRPhysio>=0.1.1
    * dcmstack>=0.7.0
    * pandas>=0.24.2
    * heudiconv>=0.8.0

## Input directory structure

The input directory with dicoms should be setup as follows:

```sh

data/
├── dicoms/
│   └── <subject>/
│        ├── <sequence>/<dicom_files.dcm>
│        ├── <sequence>/<dicom_files.dcm>
│        └── <sequence>/<dicom_files.dcm>
└── output/

```

* `data` the main working directory, the path to this directory is passed to the `docker run` command
* `dicoms` directory that stores the source DICOM files
   * `<subject>` is the identifier for the subject in the form `sub-001`, `sub-002` etc.
   * `<sequence>` is the directory for a specific imaging sequence and can be given any name
* `output` directory will store all outputs from the pipeline

## Operation dates

One main feature of this clinical pipeline is that the final output is stored in three session directories for each subject:
* **ses-presurg:** for imaging data acquired prior to the subjects surgery date
* **ses-perisurg:** for imaging data acquired on the same day of the subjects surgery (generally MRI/CT with sterotactic frame)
* **ses-postsurg:** for imaging data acquired the immediate day after the subjects surgery onwards

To divide the imaging data into the three sessions, the operation date is required for each subject. This date is attempted to be gleaned automatically assuming some type of intra-operative imaging is performed and the **SeriesDescription** and/or **StudyDescription** DICOM header tag includes an identifier that it was acquired intra-operatively. 

If the operation date cannot be determined automatically you may provide a path to an `or_dates.tsv` file that has the operation dates for all subjects. 

Here is an example of what this file should look like:

| subject  | or_date   |
|:---------|:-----------|
| sub-001<img width="100"/>  | 2014_09_28<img width="100"/> |
| sub-002  | 2018_03_26 |
| sub-003  | n/a |

If the operation date for a subject is unknown or the subject did not undergo surgery, you can put `n/a` under **or_date** for that subject. In this case, the output BIDS directory will contain only the `presurg` session for that subject (with all the imaging data stored there):

```
output/bids/sub-P001/
  └── ses-presurg/anat/...
```

## Configuration

### Modify the configuration file

#### file paths

If you are running this pipeline locally, edit the `config/config.yaml` file to include the paths to the following:

<center>

|Variable   |Description        |
|:----------|:------------------|
| `dicom_dir`<img width="200"/> | full path to where the input dicom directory is stored   |
| `out_dir`   | full path to where the pipeline should output the data   |
| `heuristic`  | heudiconv template file to sort/name dicoms according to BIDS standard |
| `or_dates_file` **[optional]** | path to the `or_dates.tsv` file with subject surgery dates |

</center>

If you modify the `config_docker.yml` file, the  docker image will need to be rebuilt. Preferably leave this configuration as set and pass input arguments to the `docker run` command to modify data paths (see Docker section below).

#### session determination

How one medical center performs aquisition of imaging data around the surgery date will most likely differ from another medical center. Further, in many clinical cases a patient will undergo surgery more than once. Within this pipeline you are able to modify how `presurgery`, `perisurgery` and `postsurgery` are defined. Within the `config/config.yaml` you will notice the following settings under **session_calc**: 
<center>

|Variable   |Description        |
|:----------|:------------------|
| `periop` [default: 0]<img width="280"/> | number of days ± around surgery day that will deemed `perisurg` |
| `dur_multi_surg` [default: -30]| in the event of multi patient surgery, the max num days before the follow-up surgery imaging data will be deemed `presurg`<sup>1</sup> |
| `override_periop` [default: True]  | in the event an imaging study is deemed `perisurg` but it contains an **electrode** flag it will be moved to `postsurg` <sup>2</sup> |

</center>

* <sup>1</sup> the default value is quantified as 30 days before the subsequent surgery any imaging acquisitions will be sorted into the `presurg` session
* <sup>2</sup> this may occur if the surgery occurs in the morning and a post-op imaging study is performed to localize implanted electrodes later the same day

### DICOM sort rules

Depending on your dicom dataset you may need to create a sorting rule for your data. There are two points in the pipeline where the imaging headers are parsed:
* within **dicom2tar** to create the Tarball archives
* within **tar2bids** to create the BIDS filename conventions

#### dicom2tar sort rule

The default sort rule can be found in [workflow/scripts/dicom2tar/sort_rules.py](workflow/scripts/dicom2tar/sort_rules.py) and is the function **sort_rule_clinical**. This sort rule seperates the DICOM files based on image type (MRI/CT/fluoro) as well as scan date. Scan sessions occuring on different days are stored in different Tarball archives, even if they are the same image type. Depending on the clinical scanner used at your center and the information stored within the DICOM header tags you may need to add/substract from this sort rule.

#### tar2bids sort rule

The sort rule for **heudiconv** is called a **heuristic**. They provide a detailed description of the [heuristic and how to write your own](https://heudiconv.readthedocs.io/en/latest/heuristics.html). The default heuristic file within this pipeline can be found in [workflow/scripts/heudiconv/clinical_imaging.py](workflow/scripts/heudiconv/clinical_imaging.py). 

## Run with Docker

To run the the pipeline in Docker, you will first need to install Docker.

### Docker Installation

#### Linux

1. Update the `apt` package index and install packages to allow `apt` to use a repository over HTTPS:

    ```sh
    sudo apt-get install -y apt-transport-https ca-certificates curl gnupg-agent software-properties-common
    ```

2. Add Docker’s official GPG key:

    ```sh
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    ```

3. Use the following command to set up the stable repository:

    ```sh
    sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
       $(lsb_release -cs) stable"
    ```

4. Install the latest version of Docker Engine and containerd:

    ```sh
    sudo apt-get install docker-ce docker-ce-cli containerd.io
    ```

5. Verify that Docker Engine is installed correctly by running the hello-world image:

    ```sh
    sudo docker run hello-world
    ```

6. **Docker-compose** is a tool for defining and running multi-container Docker applications. With Compose, you use a YAML file to configure your application’s services. Then, with a single command, you create and start all the services from your configuration. Run this command to download the current stable release of Docker Compose:

    ```sh
    sudo curl -L "https://github.com/docker/compose/releases/download/1.26.2/docker-compose-$(uname -s)-$(uname -m)" \
        -o /usr/local/bin/docker-compose
    ```

7. Apply executable permissions to the binary:

    ```sh
    sudo chmod +x /usr/local/bin/docker-compose
    ```

8. **[optional]** To run Docker without `sudo`:

    ```sh
    # Create the Docker group
    sudo groupadd docker

    #Add your user to the docker group
    sudo usermod -aG docker $USER
    ```

9. You will need to log out and log back in so that your group membership is re-evaluated.

#### MacOS

1. Download the <a href="https://download.docker.com/mac/stable/Docker.dmg">`Docker.dmg`</a> file.

2. Double-click `Docker.dmg` to open the installer, then drag the Docker icon to the Applications folder.

3. Double-click `Docker.app` in the Applications folder to start Docker.


#### Windows

1. Download the <a href="https://download.docker.com/win/stable/Docker%20Desktop%20Installer.exe">`Docker.dmg`</a> file.

2. Double-click `Docker Desktop Installer.exe` to run the installer.

3. When prompted, ensure the **Enable Hyper-V Windows Features** option is selected on the Configuration page.

4. Follow the instructions on the installation wizard to authorize the installer and proceed with the install.

5. When the installation is successful, click **Close** to complete the installation process.

6. If your admin account is different to your user account, you must add the user to the **docker-users** group. Run **Computer Management** as an administrator and navigate to **Local Users and Groups > Groups > docker-users**. Right-click to add the user to the group. Log out and log back in for the changes to take effect.

### Docker Setup

1. Clone this repository to your system:

    ```sh
    git clone https://github.com/greydongilmore/clinical_dicom2bids_smk.git
    ```

2. To setup the container in your system, run the following command in the root directory:

   ```sh
   # build the composition in `docker-compose.yml`
   docker-compose build
   ```

   ```sh
   # run the container
   docker-compose up -d
   ```
  * `build` builds the image
  * `up` builds the image if the image do not exist and starts the container
    * `--build` if you add this input argument then images will be built even when not needed (i.e. previous built container exists)
  * `d` starts the container in detached mode so it will run in the background

3. If you run `docker images ls` in the terminal you will also notice the built image `d2b-clinical_image`. Now if you type `docker ps -a`, you should see the corresponding container named `d2b-clinical`.

### Docker run

1. Now you can run the built image using the following default command:

    ```sh
    # Start the instance
    docker run -it --rm d2b-clinical_image bash
    ```

2. Inside the container the main directory will be located at `/d2b-clinical`. By default the data directory is the one found in this repository (it is copied to the Docker image when building). Inside the container, the test DICOMs are stored in `/data/dicoms` and the default output will be `/data/output`. You can test the Docker instance by running the following (in the root of `/d2b-clinical`):

    ```sh
    #Run the Snakemake pipline
    snakemake -j4
    ```

3. To link your dataset to the container you need to map the data path within the `docker run` command. Ensure your input path has the structure outline above (in the example the path supplied here would be `*/data):

    ```sh
    docker run -v <path-to-data-dir>:/data -it --rm d2b-clinical_image bash
    ```

4. Now when you run Snakemake the data ouput will appear at your specified path `output` directory.

## Run Locally

1. Install Snakemake using [Python](https://www.python.org/downloads/):

    ```sh
    python -m pip install snakemake
    ```

    For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

2. Clone a copy of this repository to your system:

    ```sh
    git clone https://github.com/greydongilmore/clinical_dicom2bids_smk.git
    ```

3. Install the Python dependencies by opening a terminal, changing to the project directory and running:

    ```sh
    python -m pip install -r requirements.txt
    ```

### Local run

All the following commands should be run in the root of the project directory.

1. Prior to running you can do a dry run:

   ```sh
   snakemake -n
   ```

2. To locally run the pipeline, run the following:

   ```sh
   snakemake -j $N
   ```

   where `$N` specifies the number of cores to use.


## Description of the pipeline

### Repository structure

The repository has the following scheme:
```
├── README.md
├── workflow
│   ├── rules
|   │   └── dicom2bids.smk
│   ├── envs
|   │   └── mapping.yaml
│   ├── scripts
|   |   ├── dicom2tar               # sorts and stores the dicoms into Tarballs
|   │   |   ├── clinical_helpers.py
|   |   |   ├── DicomSorter.py
|   |   |   ├── main.py
|   │   |   └── sort_rules.py
|   |   ├── heudiconv               # heuristic file for clinical imaging
|   │   |   ├── clinical_imaging.py 
|   |   └── post_tar2bids           # refactors heudiconv output into final BIDS structure
|   │       └── clean_sessions.py
|   └── Snakefile
├── config
│   └── config.yaml
└── data                            # contains test input dicoms
```

### Rule 01: dicom2tar

<center>

| Variable  | Description                              |
|:----------|:-----------------------------------------|
|Overview   | Sorts and stores the dicom files into Tarball files |
|Input      | MRI/CT dicoms |
|output     | MRI/CT Tarballs|

</center>

This workflow is modified from the [dicom2tar](https://github.com/khanlab/dicom2tar) master branch (version date:16/07/2020). The code has been modified to fit the current pipeline.  

The first pass will provide an intermediate output, Tarball files with the associated scan date in the name:
```
output/
  └── tars/
        ├── P185_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

The final output will sort the Tarball files by date and assign sequential numbers to the files:
```
output/
  └── tars/
        ├── P185_001_2017_09_15_20170915_I5U57IAQF6Q0.46F51E3C_MR.tar
        ├── P185_002_2017_11_10_20171110_NW1NTJVCBOJH.F06BAD6C_MR.tar
        ├── P185_003_2018_03_26_20180326_76EGBNLGE0PA.06FE3A9D_CT.tar
        ├── P185_004_2018_03_26_20180326_K6S2I5FI1MB9.92C5EAF3_CT.tar
        └── P185_005_2018_03_27_20180327_9WK33LJUKNJP.F61DC193_CT.tar
```

### Rule 02: tar2bids

<center>

| Variable  | Description                              |
|:----------|:-----------------------------------------|
|Overview   | Converts the Tarball archives into BIDS compliant format |
|Input      | MRI/CT dicom Tarball archives |
|output     | MRI/CT nifti files stored in BIDS layout|

</center>

output will be:
```
output/bids/sub-P185/
  ├── ses-001/anat/...
  ├── ses-002/anat/...
  ├── ses-003/anat/...
  ├── ses-004/anat/...
  └── ses-005/anat/...
```

### Rule 03: post tar2bids cleanup

<center>

| Variable  | Description                              |
|:----------|:-----------------------------------------|
|Overview   | Restructures the BIDS session naming to more meaningful session types |
|Input      | BIDS directory with default session numbering |
|output     | BIDS directory with session naming according to surgery date|

</center>

output will be:
```
output/bids_final/sub-P185/
  ├── ses-perisurg/anat/...
  ├── ses-postsurg/anat/...
  └── ses-presurg/anat/...
```
