
##What is NINJA-OPS?

---


NINJA is Not Just Another - OTU Picking Solution. Making use of BWT-enabled aligners, NINJA-OPS quickly generates an OTU map
and table from a given fasta file of sequence reads. NINJA-OPS also allows for convenient quality control on your data, such
as fast reverse complementing, base pair trimming, and a specialized denoising transformation.

NINJA-OPS can transform an entire MiSeq run into a QIIME-formatted biom table in under 10 minutes on your laptop, achieving higher accuracy and more exact matches than USEARCH.

Moreover, NINJA is entirely free and open source. 

The latest NINJA version can always be found on Github (https://github.com/GabeAl/NINJA-OPS) or at http://ninja-ops.ninja.


##Version
---

1.2
Removed extraneous code, added some interface features, simplified denoising.  

1.1 (November 23rd, 2015)  
Improved Stability  

1.0 (July 17th, 2015)  
Initial submitted version  

##Dependencies
---


Python 2.7 for wrapper (no additional libraries, such as numpy, are necessary)  
Bowtie2



##Installation
---



0. To use the wrapper, make sure you have the a version of Python 2.7 installed on your computer. Python can be found at https://www.python.org/downloads/
1. Download and unzip the ninja package.
2. Download and unzip bowtie2 ( http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ ). If you wish to use a different sequence aligner, you cannot run the wrapper ninja.py. Instructions below.
3. Copy the executable "bowtie-2-align-s" and paste it to the main folder of the ninja package. Alternatively, you can add the bowtie-2-align-s executable (or all Bowtie2 executables) to your system PATH variable.

IMPORTANT: Do not copy "bowtie-2-align-s" to any other folder (for instance, bt2db). NINJA-OPS looks for bowtie2 only in the same folder ninja.py is located in, or on the system executable path.

##Instructions
---

###Normal usage
Simply run "python ninja.py" in your command line to put your data through the entire NINJA pipeline. 

Sample commands are as follows:

```
# Takes a fasta file as input (-i) and outputs its OTU map, OTU table and a list of sequences culled by NINJA to the default output folder ninja_output in your current working directory
python ninja.py -i seqs.fna 

# Takes a fasta file as input (-i), reverse complements all reads within it (-r) and outputs its OTU map, OTU table and a list of sequences culled by NINJA to the folder seqs_output (-o)
python ninja.py -i seqs.fna -o seqs_output -r

# Takes a fasta file as input (-i), trims all sequences down to 200 base pairs (-t 200), denoises data by discarding all reads that appear less than 3 times and all kmers that appear less than 5 times (-d 3.0001) and outputs to folder seqs_output (-o)
python ninja.py -i seqs.fna -t 200 -d 3.005 -o seqs_output

# Takes a fasta file as input (-i) and runs bowtie2 with maximum sensitivity (-m max)
python ninja.py -i seqs.fna -m max

# Takes a fasta file as input (-i) and runs bowtie2 with maximum speed (-m fast)
python ninja.py -i seqs.fna -m fast

# Runs NINJA-OPS on paired-end data,
# Trimming from the end of the forward read to 150 base pairs
# and the beginning of the reverse read to 125 base pairs
# Allowing up to 600 bases length in the alignment from the start of the forward
# to the end of the reverse read
python ninja.py -i forward.fna,reverse.fna -o ninja -I 600 -t 150 -T 125
```


###Custom usage
You can also run each step of the ninja pipeline individually. The steps are detailed below.

NOTE: When using Bowtie2, you must point Bowtie2 to the database files by prefix (such as "Ninja97" if your files are Ninja97blabla.bla) with absolute full path. Otherwise, bowtie will fail with internal error #2. Also, any typo in the commandlines given to it will fail with error #2.


0. (Optional if bowtie2 is included in this package) Ensure bowtie2 is installed and working with the included Ninja97 DB.
1. (Optional) Compile the c files with gcc. E.g. (you may have to replace * with each real ninja c filename):  
	```gcc -m64 -O3 -flto ninja_*.c -o ninja_*```
1b. OR use the included binaries
2. Run ninja_filter on your raw fasta reads. E.g. for an input "seqs.fna":  
	```./ninja_filter seqs.fna my_prefix```
You can include e.g, "170 RC" after that to trim at 170 bases and perform reverse complementing (doing it here is much faster than elsewhere).  
You can append D X.Y to the NINJA_filter commandline, where X is whole read (duplicate) filtering level and Y is the rarity threshold for k-mer filtering.
3. Run your filtered sequences through bowtie2 or other aligner with headerless SAM output. E.g. for 4 threads (p option) and 97% id (one plus the third parameter in --score-min) on seqs_filtered.fna with the db in directory /bt2db:  
	```bowtie2-align-s --no-head -x /bt2db/Ninja97 -S alignments.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.03" -k 1 --norc -p 4 -f seqs_filtered.fna```
	
For presets instead, specify e.g. ```--preset fast```

or ninjaMAX instead, use this combo: ```--no-head -x Ninja97 -S outTrial.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.03" --norc -p 4 -f filtered.fa --very-sensitive -D 40 -R 4 -N 0 -L 10 -i "S,1,0.50"```
	
4. Run ninja_parse on the alignments, e.g. with the included 97 taxonomy mapping files masterDB97.map and taxmap97:  
	```./ninja_parse my_prefix alignments.txt masterDB97.map taxmap97.txt```
	
For legacy table output, add ``--legacy`` to the end of the command above
(Also, to log which sequences failed to align, add LOG at the very end of both ninja_filter and ninja_parse_filtered)


Contact
---

Created by Gabe Al Ghalith, Henry Ward, Emmanuel Montassier, and Dan Knights in the Knights Lab at the University of Minnesota (algh0022@umn.edu). Questions, comments and concerns can also be directed to Dan Knights (dknights@umn.edu).


Licensing
---
NINJA-OPS is licensed under the ISC license (included in package). Bowtie2 is distributed under the GPLv3 license. 
