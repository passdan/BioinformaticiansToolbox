# Queuing with SLURM
So far everything we have done has been directly on the commandline. You type your code or program, press enter and wait for it to complete. However this has some issues, particularly for long processes and unstable internet connections (wi-fi goes down, your job stops!), and when many people are trying to use a server at the same time.

Most server clusters have a queuing system and we'll look at (one of?) the most popular one: **SLURM**. This allows you to put your process into a file, submit it and then it will join a queue and start when there is processing space available without you needing to be online.

## Checking the slurm queue status
Two important commands will let you see the status of the queue and the options:
```
$ squeue

```

```
$ sinfo -O "%C %P"
```

## Making your script into a slurm script
### Loop script
Lets take a simple QC script:
```
#!/bin/bash
mkdir -P trim_fastq

for file in *R1.fastq.gz
do
    R1=$(basename $file | cut -f1 -d.)
    base=$(echo $R1 | sed 's/_R1$//')

    fastp   -i fastq/${base}_R1.fastq.gz \
            -I fastq/${base}_R2.fastq.gz \
            -o trim_fastq/${base}_trim_R1.fastq.gz\
            -O trim_fastq/${base}_trim_R2.fastq.gz \
            -w 2
done
```

### Loop script in a slurm
Now lets make that script into a slurm script which we can submit. For the most simple version we just need to create a header for the loop and then can paste the same script. 

Note the ```${SLURM_CPUS_PER_TASK}``` parameter, which takes it's value from the header

```
#!/bin/bash
#SBATCH --job-name=fastqc-job   # Job identifier
#SBATCH --partition=normal      # Name of queue
#SBATCH --output=%j.out         # Output file
#SBATCH --error=%j.err          # Error file
#SBATCH --cpus-per-task=4       # Request cores per node
#SBATCH --mem-per-cpu=2400M     # Request RAM per core

mkdir -P trim_fastq

for file in *R1.fastq.gz
do
    R1=$(basename $file | cut -f1 -d.)
    base=$(echo $R1 | sed 's/_R1$//')

    fastp   -i fastq/${base}_R1.fastq.gz \
            -I fastq/${base}_R2.fastq.gz \
            -o trim_fastq/${base}_trim_R1.fastq.gz\
            -O trim_fastq/${base}_trim_R2.fastq.gz \
            -w ${SLURM_CPUS_PER_TASK}
done
```
This is a minimum example, there are many more optional parameters such as "email me when job completes".

Also you may need to load modules or docker/singularity containers, or give full paths to your files in all cases.

### Array slurm script
What we have above works fine to make a process run in a queue system, however we can improve. That will do all steps in a series, one after the other like a normal loop.

However we can instead create an array of jobs that start independently and use the processing power more efficiently:
```
#!/bin/bash
#SBATCH --job-name=fastqc-job   # Job identifier
#SBATCH --partition=normal      # Name of queue
#SBATCH --output=%j.out         # Output file
#SBATCH --error=%j.err          # Error file
#SBATCH --cpus-per-task=4       # Request cores per node
#SBATCH --mem-per-cpu=2400M     # Request RAM per core
#SBATCH --array=1-30%4          # Do 30 samples, 4 at a time

mkdir -P trim_fastq

# Note the ${SLURM_ARRAY_TASK_ID} which containes a range of files
file=$(ls fastq/*_R1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

R1=$(basename $file | cut -f1 -d.)
base=$(echo $R1 | sed 's/_R1$//')

fastp   -i fastq/${base}_R1.fastq.gz \
        -I fastq/${base}_R2.fastq.gz \
        -o trim_fastq/${base}_trim_R1.fastq.gz\
        -O trim_fastq/${base}_trim_R2.fastq.gz \
        -w ${SLURM_CPUS_PER_TASK}
```
### Running your slurm script
We have our script, now we need to submit it to the queue. Super simple, and it will return the jobID:
```
$ sbatch my_fastqc_slurm.sh
512738
```
Made a mistake?? Cancel it with scancel & the jobID:
```
$ scancel 512738
```

That should cover the main aspects of using a queuing system, SLURM or another! There are lots more parameters and I recommend reading the full SLURM documentation. Your institute will likely also have some of their own rules for you to follow such as maximum time limits!

Good luck!