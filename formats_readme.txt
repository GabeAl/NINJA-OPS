FORMATS README
--------------
There are some rules regarding formatting in NINJA. All FASTA files must be formatted with a header line, followed by
the newline character (\n), followed by the sequence in uppercase or lowercase without any additional newlines or line
breaks, followed by the next header starting with ">", and so on until the end of the file, which terminates in a newline
(\n). 
For example:
>Header
SEQUENCE WITHOUT LINE BREAKS
>Header2
SEQUENCE2 WITHOUT LINE BREAKS
...
[blank end line]

Header formatting instructions -- DATABASE:
Database fasta file restrictions are less stringent. Try your file and report if it does not work.

Header formatting instructions -- QUERIES/READS:
Headers are formatted according to QIIME protocol. That is, minimally, in the format SAMPLENAME_SEQNUMBER [extra info].
For example, the first sequence header may look like the following:
>100217.1246514_0 A222Y121112:1:1101:10003:15638/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0
where 100217.1246514 is a sample, and 0 is the index corresponding to the first read from the sequencing run.

Other non-wrapped fasta files will work; however, the nubmer of samples may be cast to 1 if this header type is not used.