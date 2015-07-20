
What is NINJA-OPS?

- - - - - - - - - 


NINJA is Not Just Another - OTU Picking Solution. Making use of BWT-enabled aligners, NINJA quickly generates an OTU map
and table from a given fasta file of sequence reads. NINJA also allows for convenient quality control on your data, such
as fast reverse complementing, base pair trimming, and a specialized denoising transformation. 

Moreover, NINJA is entirely free and open source. 

The latest NINJA version can always be found on Github (https://github.umn.edu/algh0022/NINJA-OPS).



Version

- - - - - - - - - 


1.0 (July 17th, 2015)



Dependencies
- - - - - - - - - 


Python 2.7 for wrapper (no additional libraries, such as numpy, are necessary)
Bowtie2



Installation

- - - - - - - - - 



0) To use the wrapper, make sure you have the a version of Python 2.7 installed on your computer. Python can be found at https://www.python.org/downloads/
1) Download and unzip the ninja package.
2) Download and unzip bowtie2 ( http://sourceforge.net/projects/bowtie-bio/files/bowtie2/ ). If you wish to use a different sequence aligner, you cannot run the wrapper ninja.py. Instructions below.
3) Take the executable "bowtie-2-align-s" and copy it to the main folder ninja_package.

IMPORTANT: Do NOT copy "bowtie-2-align-s" to any other folder (for instance, bt2db). Ninja.py looks for bowtie only in the same folder ninja.py is located in.



Instructions (wrapper)

- - - - - - - - - 


With the NINJA wrapper, simply run "python ninja.py" in your command line to put your data through the entire NINJA pipeline. 

Sample commands are as follows:



# Takes a fasta file as input (-i) and outputs its OTU map, OTU table and a list of sequences culled by NINJA to the default output folder ninja_output in your current working directory

> python ninja.py -i seqs.fna 



# Takes a fasta file as input (-i), reverse complements it (-r) and outputs its OTU map, OTU table and a list of sequences culled by NINJA to the folder seqs_output (-o)

> python ninja.py -i seqs.fna -o seqs_output -r



# Takes a fasta file as input (-i), trims all sequences down to 200 base pairs (-t 200), denoises data by discarding all reads that appear less than 3 times and all kmers that appear less than 1 percent of the time (-d 3.0001) and outputs to folder seqs_output (-o)

> python ninja.py -i seqs.fna -t 200 -d 3.0001



# Takes a fasta file as input (-i) and runs bowtie2 with maximum sensitivity (-m max)

> ninja.pyc -i seqs.fna -m max



# Takes a fasta file as input (-i) and runs bowtie2 with maximum speed (-m fast)

> ninja.pyc -i seqs.fna -m fast



Without the ninja wrapper, you have to run each step of the ninja pipeline individually. The steps are detailed below.



Instructions (without wrapper)

- - - - - - - - - 


NOTE: When using Bowtie2, you must point Bowtie2 to the database files by prefix (such as "Ninja97" if your files are Ninja97blabla.bla) with absolute full path. Otherwise, bowtie will fail with internal error #2. Also, any typo in the commandlines given to it will fail with error #2.



0. (Optional if bowtie2 is included in this package) Ensure bowtie2 is installed and working with the included Ninja97 DB.



1. (Optional) Compile the c files with gcc. E.g. (you may have to replace * with each real ninja c filename): 
	gcc -m64 -Ofast -flto ninja_*.c -o ninja_*



2. Run ninja_filter on your raw fasta reads. E.g. for an input "seqs.fna":
	./ninja_filter seqs.fna seqs_filtered.fna seqs.db
	
You can include e.g, "170 RC" after that to trim at 170 bases and perform reverse complementing (doing it here is much faster than elsewhere)
You can append D X.Y to the NINJA_filter commandline, where X is whole read (duplicate) filtering level and Y is the rarity threshold for k-mer filtering


3. Run your filtered sequences through bowtie2 or other aligner. E.g. for 4 threads (p option) and 97% id (one plus the third parameter in --score-min) on seqs_filtered.fna with the db in directory /bt2db:
	bowtie2-align-s --no-head --no-unal -x /bt2db/Ninja97 -S alignments.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.03" -k 1 --norc -p 4 -f seqs_filtered.fna
	
For presets instead, specify e.g. --preset fast

For ninjaMAX instead, use this combo: --no-head --no-unal -x Ninja97 -S outTrial.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.03" --norc -p 4 -f filtered.fa --very-sensitive -D 40 -R 4 -N 0 -L 10 -i "S,1,0.50"
	


4. Run ninja_parse on the alignments, e.g. with the included 97 taxonomy mapping files masterDB97.map and taxmap97:
	./ninja_parse seqs.db alignments.txt masterDB97.map taxmap97.txt otutable.biom
	
For legacy table, add --legacy and keep the file as .txt



Contact
- - - - - - - - - 


Created by Gabe Al Ghalith in the Knights Lab at the University of Minnesota (algh0022@umn.edu). Questions, comments and concerns can also be directed to Dan Knights (dknights@umn.edu). NINJA wrapper created by Henry Ward (henry.n.ward@lawrence.edu).



Licensing
- - - - - - - - - 


NINJA-OPS is licensed under the ISC license (included in package). Bowtie2 is distributed under the GPLv3 license. 