# Introduction to Using Singularity & Dockerhub
In scientific settings, Singularity is often preferred due to its user-friendly security and integration with HPC cluster environments. Docker, while versatile, often requires higher administrative control, making it less ideal for multi-user scientific systems. Singularity also allows direct use of Docker images from dockerhub, providing greater flexibility for researchers.

Here, we focus on how to use Singularity with two bioinformatics tools: FastQC and MultiQC, commonly used for quality control in sequencing data analysis.

**Contents**
- [Getting started](#getting-started)
- [Pulling FastQC and MultiQC Images](#pulling-fastqc-and-multiqc-images)
- [Executing Commands Using Singularity](#executing-commands-using-singularity)
- [Pull and execute together](#pull-and-execute-together)
- [Building a sif from a definition file](#building-a-sif-from-a-definition-file)
  - [Basic multiQC Singularity definition file](#basic-multiqc-singularity-definition-file)

## Getting started
Firstly, [install singularity](https://docs.sylabs.io/guides/4.0/user-guide/quick_start.html). If on a HPC cluster your systems administrator will need to do this. A good encouragement for them is saying that after this, you'll never ask them to install something for you again!

## Pulling FastQC and MultiQC Images
To pull FastQC and MultiQC Docker images and convert them into Singularity images, run the following commands:

```
singularity pull docker://biocontainers/fastqc:v0.11.9_cv8
singularity pull docker://ewels/multiqc:latest
```

These commands will save the images as .sif files in your current directory, such as fastqc_v0.11.9_cv8.sif and multiqc_latest.sif. You can then refer to them in your pipeline

## Executing Commands Using Singularity
After pulling the images, you can run FastQC and MultiQC locally with Singularity using the exec command. For instance, to run FastQC:

```
singularity exec fastqc_v0.11.9_cv8.sif fastqc sample.fastq
```

We can also use the bind function to specify folders on your system to be visable inside the container. Here we bind the current working directory to the same location in the Singularity container and execute MultiQC in that directory MultiQC.
```
singularity exec --bind `pwd`:`pwd` --pwd `pwd` multiqc_latest.sif multiqc .
```

These exec commands let you execute the respective program as if they were installed locally, using the containers you pulled from DockerHub.

Note that Singularity doesn't store downloaded and available containers like docker does (viewed with ```docker images```), and you should direct to the downloaded ```.sif``` files.

## Pull and execute together
You can directly execute commands from a Docker image without explicitly pulling it first. Singularity will automatically pull the image and run the command in a single step. 

To run FastQC directly from its Docker image without pulling it first:
```
singularity exec docker://biocontainers/fastqc:v0.11.9_cv8 fastqc sample.fastq
```

Similarly, for MultiQC:
```
singularity exec --bind `pwd`:`pwd` --pwd `pwd` docker://ewels/multiqc:latest multiqc .
```
In both cases, Singularity pulls the Docker image, converts it to a Singularity image, and executes the specified command, all in one go. This can be particularly useful for one-off tasks or for testing different versions of a tool without keeping the image locally.

## Building a sif from a definition file
You may also have a text definition file that specifies how to build the container. This is a limited file that uses ubuntu as the base, and installs the two programs.

### Basic multiQC Singularity definition file
```
Bootstrap: docker
From: ubuntu:latest

%post
    # Update and install Python
    apt-get update && apt-get install -y python3 python3-pip

    # Install MultiQC
    pip3 install multiqc

%runscript
    echo "Run MultiQC with 'multiqc'"
```

Save this content into a file named, for example, MultiQC_Singularity.def. Then, you can build the Singularity image with:
```
sudo singularity build MultiQC.sif MultiQC_Singularity.def
```

Note the requirement of sudo for building the image, so you may need to build the image locally rather than in your HPC environment.

[def]: #introduction-to-using-singularity--dockerhub
