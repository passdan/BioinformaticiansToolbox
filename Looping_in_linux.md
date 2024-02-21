# Loops – When processing lots of data

A ‘for loop’ in linux allows you to run the same command repeatedly, with slight change i.e. for lots of different files, or changing a parameter for each loop. 

The ```$``` character says that what follows is a variable, and the ```{i}``` is the name of the variable. All the commands in the loop need to be surrounded by ‘do’ and ‘done’

### LOOPS USING NUMERICAL VARIABLES
To run a command repeatedly but changing the parameter k each time
```
for i in {1..25}
do
   myCommand.py -k ${i} myFile.fasta
done
```

### LOOPS USING STRINGS (LISTS) AS VARIABLES
1. To run a command repeatedly but changing the input file each time:
```
for i in *fastq
do
   myCommand.py -k 7 ${i}
done
```

2. Alternative method creating a list of IDs:
```
# List of sample 'root' names
list=("sample1" "sample2" "sample3")

for i in ${list[@]}
do
  myCommand.py ${i}_1.fastq ${i}_2.fastq ....
done
```

3. Alternative reading a list of files from a text file:
```
## Samples read in from a text file
for i in $(cat listOfSamples.txt)
do
  myCommand.py ${i}_1.fastq ${i}_2.fastq ....
done
```

Similarly, without needing to create a list file we can list all of the R1 files. Here is something more complex which takes all of the forward files in a pair of fastqs, and extracts just the name:
```
## Samples read from the terminal
for file in $(ls *R1.fastq)
do
  R1=$(basename $file | cut -f1 -d.)
  base=$(echo $R1 | sed 's/_R1$//')
  myCommand.py ${base}_1.fastq ${base}_2.fastq ....
done
```

### EXERCISES
#### Exercise 1
1. Make a new script using nano or vi and run the following script:
```
#!/bin/bash

for i in {10..1}
do
   echo "T-minus: $i"
done
echo 'Blastoff!'
```

2. Remember to make the script executable once made! 
3. Think about: How is this script different to the countdown one we made earlier?


And just in case you forget for vi (this is to exit without saving!):
![exit vi](images/exit-vi.png)

#### Exercise 2
1. Write a script that will run fastqc on each of the four fastq files in the folder Day1/looping but only writing the command once.

### [Extension]
2. Once the fastqc loop in step one is working, add fastp trimming to your script!
3. Once both steps are working, repeat fastqc on the trimmed data, and run multiqc on all the outputs.

It might seem simple now but will be really useful when working with lots of files! (Example 3 (list of files) is the version that will scale up best to lots of files).