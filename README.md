
## What is NINJA-OPS?

---


NINJA is Not Just Another - OTU Picking Solution. Making use of BWT-enabled aligners, NINJA-OPS quickly generates an OTU map
and table from a given fasta file of sequence reads. NINJA-OPS also allows for convenient quality control on your data, such
as fast reverse complementing, base pair trimming, and a specialized denoising transformation.

NINJA-OPS can transform an entire MiSeq run into a QIIME-formatted biom table in under 10 minutes on your laptop, achieving higher accuracy and more exact matches than USEARCH.

Moreover, NINJA is entirely free and open source. 

The latest NINJA version can always be found on Github (https://github.com/GabeAl/NINJA-OPS) or at http://ninja-ops.ninja.


## Citing NINJA-OPS
Al-Ghalith GA, Montassier E, Ward HN, Knights D. NINJA-OPS: Fast Accurate Marker Gene Alignment Using Concatenated Ribosomes. PLoS Computational Biology. 2016 Jan;12(1).


## Version
---
1.5 (March 20, 2016)
New OTU compaction step using compressed concatesome (compressasome), enabled by default in wrapper. This pools tied OTU assignments into a single arbitrary OTU ID instead of maintaining a distribution across equally-likely OTU candidates. Additionally, NINJA-OPS -m max (NINJA-MAX) settings are now the default speed settings. 

1.4 (March 13, 2016)
Updated readme, restructured file directory, included output OTU map in full output mode, rewrote ninja_prep for better database generation, chimera detection bugfixes

1.3 (January 29, 2016)  
Added explicit paired-end support, filter reasons to log file, streamlined BIOM format, improved Mac stability, alternative ambiguity handling schemes, updated commandline interface, and more.

1.2 (December 7, 2015)  
Removed extraneous code, added some interface features, simplified denoising.  

1.1 (November 23rd, 2015)  
Improved Stability  

1.0 (July 17th, 2015)  
Initial submitted version  

## Dependencies
---


Bowtie2

Python 2.7 for wrapper



## Installation
---



0. To use the wrapper, you need to have a version of Python 2.7 installed on your computer. If you have a recent version of the Mac OS or Linux/Unix OS then you will already have Python 2.7. You can test this by opening the terminal and typing "python". It will print the version. Then type "exit()" to quit. If you have Windows you may need to install Python. Python can be found at https://www.python.org/downloads/
1. Download and unzip the ninja package.
2. Download and unzip bowtie2 ( http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ ). If you wish to use a different sequence aligner, you cannot run the wrapper ninja.py. Instructions below.
3. Copy the executable "bowtie2-align-s" and paste it to the main folder of the ninja package. Alternatively, you can add the bowtie2-align-s executable (or all Bowtie2 executables) to your system PATH variable.

IMPORTANT: Do not copy "bowtie2-align-s" to any other folder (for instance, bt2db). NINJA-OPS looks for bowtie2 only in the same folder ninja.py is located in, or on the system executable path.

### Instructions
---

## Normal usage
Simply run `python /path/to/your/ninja/directory/bin/ninja.py` in your command line to put your data through the entire NINJA pipeline. You will have to replace `/path/to/your/ninja/directory/bin/` with the proper filepath to your `ninja.py` file, which is in the `bin` directory in the NINJA-OPS folder.

To download and use more databases, see the readme inside the 'databases' subdirectory. 

Sample commands are as follows:

```
# Takes a fasta file as input (-i) and outputs its OTU map, OTU table and a
# list of sequences culled by NINJA to the default output folder
# ninja_output in your current working directory
python /path/to/ninja.py -i seqs.fna 

# Takes a fasta file as input (-i), reverse complements all reads within it
# (-r) and outputs its OTU map, OTU table and a list of sequences culled by
# NINJA to the folder ninja (-o)
python /path/to/ninja.py -i seqs.fna -o ninja -r

# Takes a fasta file as input (-i), trims all sequences down to 200 base pairs
# (-t 200), denoises data by discarding all reads that appear fewer than 3
# times and outputs to folder ninja (-o)
python /path/to/ninja.py -i seqs.fna -t 200 -d 3 -o ninja

# Takes a fasta file as input (-i) and runs bowtie2 with maximum sensitivity (-m max)
python /path/to/ninja.py -i seqs.fna -m max

# Takes a fasta file as input (-i) and runs bowtie2 with maximum speed (-m fast)
python /path/to/ninja.py -i seqs.fna -m fast

# Runs NINJA-OPS on paired-end data,
# Trimming from the end of the forward read to 150 base pairs
# and the beginning of the reverse read to 125 base pairs
# Allowing up to 600 bases length in the alignment from the start of the forward
# to the end of the reverse read
python /path/to/ninja.py -i forward.fna,reverse.fna -o ninja -I 600 -t 150 -T 125
```


## Custom usage
You can also run each step of the ninja pipeline individually. The steps are detailed below.

NOTE: When using Bowtie2, you must point Bowtie2 to the database files by prefix (such as "Ninja97" if your files are Ninja97blabla.bla) with absolute full path. Otherwise, bowtie will fail with internal error #2. Also, any typo in the commandlines given to it will fail with error #2.


0. (Optional if bowtie2 is included in this package) Ensure bowtie2 is installed and working with the included Ninja97 DB.
1. (Optional) Compile the c files with gcc. E.g. (you may have to replace * with each real ninja c filename):  
	```gcc -m64 -std=gnu99 -O3 ninja_*.c -o ninja_*```
(OR use the included binaries)
2. Run ninja_filter on your raw fasta reads. E.g. for an input "seqs.fna":  
	```./ninja_filter seqs.fna my_prefix```
You can include e.g, "170 RC" after that to trim at 170 bases and perform reverse complementing (doing it here is much faster than elsewhere).  
You can append D X.Y to the NINJA_filter commandline, where X is whole read (duplicate) filtering level and Y is the rarity threshold for k-mer filtering.
3. Run your filtered sequences through bowtie2 or other aligner with headerless SAM output. E.g. for 4 threads (p option) and 97% id (one plus the third parameter in --score-min) on seqs_filtered.fna with the db in directory /bt2db:  
	```bowtie2-align-s --no-head -x /bt2db/Ninja97 -S alignments.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.03" -k 1 --norc -p 4 -f seqs_filtered.fna```
For presets instead, specify e.g. ```--preset fast```
or NINJA-MAX instead, use this combo: ```--no-head -x Ninja97 -S alignments.sam --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.03" --norc -p 4 -f filtered.fa --very-sensitive -D 40 -R 4 -N 0 -L 10 -i "S,1,0.50"```
4. (Optional) Run ninja_compact on the alignments, e.g. with the included compressasome:  
	```./ninja_compact alignments.sam greengenes97.tcf new_alignments.sam```
	
5. Run ninja_parse on the alignments, e.g. with the included 97 taxonomy mapping files masterDB97.map and taxmap97:  
	```./ninja_parse my_prefix new_alignments.sam masterDB97.map taxmap97.txt```
	
For legacy table output, add ``--legacy`` to the end of the command above
(Also, to log which sequences failed to align, add LOG at the very end of both ninja_filter and ninja_parse_filtered)


## Building a database
You can build your own marker gene database for use with NINJA-OPS. To do so, use ninja_prep in the /bin directory. 
The format of the input FASTA file must respect the formatting considerations outlined in formats_readme.txt in the root directory of NINJA-OPS. Other formats may work, but use at your own risk.

Here is an example workflow for creating a new database using the 99% greengenes 13_8 OTUs, where the input reference FASTA file is called ```99_otus.fasta```:

```ninja_prep 99_otus.fasta greengenes99```

This produces a concatesome FASTA called ```greengenes99.fa``` and a NINJA-OPS db called ```greengenes99.db```. Then build a bowtie2 database for the concatesome using ```bowtie2-build-s``` like this:

```bowtie2-build-s greengenes99.fa greengenes99```

(On Linux or Mac, if bowtie is not in your PATH, you may have to copy the bowtie2-build-s binary to the folder with your sequences and run it as ./bowtie2-build-s)

To include taxonomic annotation in your OTU tables, create a tab-delimited text file with two columns, the first containing the OTU ids, and the second containing the corresponding taxonomy in string form. The taxonomy should be semicolon-delimited if compartments are desired in BIOM format. No sorting is required as of v1.4. The taxonomy file should named ```greengenes99.taxonomy```. The taxonomy file is optional.

You may now place all files starting with the name ```greengenes99``` (or whatever your database name is) except for the concatesome file into a folder named ```greengenes99``` within the "databases" directory of NINJA-OPS. You may delete the concatesome file. The folder should contain (with the ```.taxonomy``` file optional):

```
greengenes99.1.bt2
greengenes99.2.bt2
greengenes99.3.bt2
greengenes99.4.bt2
greengenes99.db
greengenes99.rev.1.bt2
greengenes99.rev.2.bt2
greengenes99.taxonomy
greengenes99.tcf
```

Now NINJA-OPS may be invoked using the custom database like this:

```ninja.py -i input.fna -b greengenes99```


### Contact
---

Created by Gabe Al-Ghalith (algh0022@umn.edu), Emmanuel Montassier, Henry Ward, and Dan Knights* (dknights@umn.edu) in the Knights Lab at the University of Minnesota.


## Licensing
---
NINJA-OPS is licensed under the ISC license (included in package). Bowtie2 is distributed under the GPLv3 license. 

