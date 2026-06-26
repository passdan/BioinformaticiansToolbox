## Refresher on Commandline Bioinformatics

N.b. This is **not** a from-scratch introduction! It assumes you've used the terminal before but might be rusty. If you've never used Linux at all, work through the fuller `Introduction_to_Linux.md` first.

By the end you'll have moved around the filesystem, pulled code and data off the internet, run a couple of real bioinformatics tools, and wrapped one of them in a loop.

## Before you start

A few reminders that prevent 90% of early errors:

- The leading `$` in examples just marks a command line — you don't type it.
- **Tab-complete everything.** Start a file or folder name, press `Tab`, and let the shell finish it. If it won't complete, it doesn't exist (usually a typo or you're in the wrong folder).
- `Ctrl + C` cancels a stuck command; `q` quits a manual (`man`) or `less` view.
- Never put spaces in file or folder names: use `_` or `-`.
- Lost? `pwd` tells you where you are; `ls` shows what's around you.

Throughout, we run tools (`fastqc`, `fastp`) directly by name. In the main sessions we'll run them inside conda, or singularity containers (don't worry about that yet!).

## 1. Moving around the filesystem

The everyday commands, in one place:

| Command | What it does |
| -- | -- |
| `pwd` | **p**rint **w**orking **d**irectory — where am I? |
| `ls` | list files here (`ls -l` for detail, `ls -lh` for human-readable sizes) |
| `cd` | change directory (`cd ..` up one, `cd ~` home, `cd -` back to previous) |
| `mkdir` | make a new directory (`mkdir -p a/b/c` makes the whole path) |
| `cp` | copy a file (`cp -r` for a whole folder) |

```
$ pwd                      # where am I?
$ ls -lh                   # what's here, with sizes
$ mkdir -p refresher       # make a working folder
$ cd refresher             # move into it
$ cp ~/Share/somefile.txt .   # copy a file into "here" (the . means current folder)
```

### EXERCISES

1. Use `pwd` to confirm where you are, then `ls -lh` to see what's in your home directory.
2. Make a fresh working directory called `refresher` and `cd` into it. Confirm with `pwd` that you're inside it.
3. From inside `refresher`, make a subfolder called `raw_data` (you'll use it shortly).

## 2. Getting code and data: git and wget

Two ways to bring things onto the server from the internet.

**`git clone`** copies an entire GitHub repository (code, scripts, small data) to your machine:

```
$ git clone https://github.com/<user>/<repo>.git
```

This creates a folder named after the repo. It's how you'll grab published pipelines and scripts — including the Nextflow ones later.

**`wget`** downloads a single file from a URL — handy for data that's too big to keep in a git repo:

```
$ wget https://example.com/path/to/file.fastq.gz
```

For this session we'll download a small set of practice reads as a single compressed archive, then unpack it with `tar`:

```
$ cd ~/refresher
$ wget http://<location>/fastqs.tar   # download the archive
$ tar -xf refresher_fastqs.tar -C raw_data     # extract into raw_data/
$ ls raw_data
sample1.fastq.gz  sample2.fastq.gz  sample3.fastq.gz  sample4.fastq.gz  sample5.fastq.gz
```

`tar -xzf` means e**x**tract, un**z**ip (gzip), from this **f**ile; `-C raw_data` puts the results in that folder.

### EXERCISES

There is a tar zipped file containing 5 raw Chip-seq fastq files at the http://github.com/passdan/BioinformaticiansToolbox/fastqs. Copy the URL for the file in that folder, then:

1. Download the fastq archive from with `wget` and extract it into your `raw_data` folder.
2. Use `ls -lh raw_data` to confirm you have **five** `.fastq.gz` files. How big is each one?

## 3. Running a command-line program

A bioinformatics tool is just another command: the program name, some flags (aka parameters), and the file(s) to work on.

**FastQC** produces a quality report for a sequencing file:

```
$ fastqc raw_data/sample1.fastq.gz
```

This writes `sample1_fastqc.html` (open it in a browser) and a `.zip` of the details, right next to the input.

**fastp** trims low-quality bases and adapters, and also writes its own QC report:

```
$ fastp \
    -i raw_data/sample1.fastq.gz \
    -o trimmed/sample1.trim.fastq.gz
```

Here `-i` is the **i**nput and `-o` the **o**utput. fastp won't create the output folder for you, so make it first with `mkdir -p trimmed`.

> Tip: almost every tool has built-in help i.e. `fastqc -h` or `fastp -h`,  which lists every flag. When a command fails, the first suspects are a typo, a missing output folder, or pointing at a file that isn't there (check with `ls`).

### EXERCISES
1. Ensure you have fastqc and fastp available. Depending on your system this could be via modules, or alternatively download them directly using wget from [**here** (fastqc)](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip) and [**here** (fastp)](http://opengene.org/fastp/fastp). Be sure to unzip the folders after they download
2. Run `fastqc` on one of the samples in the `raw_data` folder and confirm the `.html` report appears.
3. Make a `trimmed` folder, then run `fastp` on the sample as above. Read the summary fastp prints to the screen. How many reads passed the filter?
4. Run `fastqc` again, this time on your **trimmed** file, and compare it to the raw one.

## 4. Writing a simple loop

Running `fastp` by hand on one sample is fine. Doing it five times (or five hundred) by hand is not. A `for` loop runs the same command over every file that matches a pattern.

Create a script with `nano` (e.g. `nano run_fastp.sh`) and write:

```
#!/bin/bash

mkdir -p trimmed

for sample in raw_data/*.fastq.gz
do
    name=$(basename $sample .fastq.gz)  
    fastp -i $sample -o trimmed/${name}.trim.fastq.gz
done
```

Breaking it down:

- `for sample in raw_data/*.fastq.gz` — the `*` glob matches all five files; the loop runs once per file, with `$sample` holding the current one.
- `basename $sample .fastq.gz` strips the folder and the extension so we can build a tidy output name.
- everything between `do` and `done` runs each time.

Make it executable and run it:

```
$ chmod +x run_fastp.sh
$ ./run_fastp.sh
```

### EXERCISES

1. Write and run the loop above. Confirm with `ls trimmed` that you have five trimmed files.
2. Add a `fastqc` line inside the loop so each trimmed file is also QC'd in the same pass.
3. The loop printed a line per sample — change the `echo` message to also report which file it's writing to.

<details>
  <summary>Example loop with both steps (don't look until you've tried!)</summary>

```
#!/bin/bash

mkdir -p trimmed fastqc

for sample in raw_data/*.fastq.gz
do
    name=$(basename $sample .fastq.gz)
    echo "Processing $name -> trimmed/${name}.trim.fastq.gz"
    fastp -i $sample -o trimmed/${name}.trim.fastq.gz
    fastqc trimmed/${name}.trim.fastq.gz -o fastqc
done
```
</details>

<details>
  <summary>[Extension] Looping over an explicit list</summary>

Instead of a glob you can loop over a named list — useful when your samples don't share a tidy pattern. You'll see this style in the course's processing scripts:

```
#!/bin/bash
list=("sample1" "sample2" "sample3" "sample4" "sample5")

for name in ${list[@]}
do
    fastp -i raw_data/${name}.fastq.gz -o trimmed/${name}.trim.fastq.gz
done
```
</details>

## Ready for Nextflow

Look again at what that loop does: **run one tool over many files, naming the outputs from the inputs.** That is the single most common pattern in bioinformatics 

In the Nextflow session you'll see the same idea expressed differently:

- the **glob** (`raw_data/*.fastq.gz`) becomes a **channel** of files,
- the **command** in your loop body becomes a **process**,
- and instead of running one-after-another, Nextflow runs them **in parallel** and tracks every output so that you can add samples, resume after a failure, and scale to a cluster without rewriting anything.
