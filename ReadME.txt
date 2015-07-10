Dependency: Bowtie2 (in vanilla configuration; any BWT aligner that outputs per-read valid alignments to headerless SAM will work)

You must also point Bowtie2, if you use it, to the database files, by prefix (such as "Ninja97" if your files are Ninja97blabla.bla), with absolute full path. Otherwise it will fail with internal error #2. Also, any typo in the commandlines given to it will fail with Error #2.



How to use NINJA (without wrapper):
0. Ensure bowtie2 is installed and working with the included Ninja97 DB.
1. Compile the c files with gcc. E.g. (you may have to replace * with each real ninja c filename): 
	gcc -m64 -Ofast -flto ninja_*.c -o ninja_*
2. Run ninja_filter on your raw fasta reads. E.g. for an input "seqs.fna":
	./ninja_filter seqs.fna seqs_filtered.fna seqs.db
	You can include e.g, "170 RC" after that to trim at 170 bases and perform reverse complementing (doing it here is much faster than elsewhere)
3. Run your filtered sequences through bowtie2 or other aligner. E.g. for 4 threads (p option) and 97% id (one plus the third parameter in --score-min) on seqs_filtered.fna with the db in directory /bt2db:
	bowtie2-align-s --no-head --no-unal -x /bt2db/Ninja97 -S alignments.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.03" -k 1 --norc -p 4 -f seqs_filtered.fna
	For presets instead, specify e.g. --preset fast
	For ninjaMAX instead, use this combo: --no-head --no-unal -x Ninja97 -S outTrial.txt --np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-0.03" --norc -p 4 -f filtered.fa --very-sensitive -D 40 -R 4 -N 0 -L 10 -i "S,1,0.50"
	For filtering, append D X.Y to the NINJA_filter commandline, where X is whole read (duplicate) filtering level and Y is the rarity threshold for k-mer filtering
4. Run ninja_parse on the alignments, e.g. with the included 97 taxonomy mapping files masterDB97.map and taxmap97:
	./ninja_parse seqs.db alignments.txt masterDB97.map taxmap97.txt otutable.biom
	For legacy table, add --legacy and keep the file as .txt
