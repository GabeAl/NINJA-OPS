Silva 119 database
==================

This index was building using SILVA version 119 and bowtie2-build version 2.2.4
The following steps were taken:

Retrieve SILVA (QIIME compatible) database:

$ wget http://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_119_provisional_release.zip
$ unzip Silva_119_provisional_release.zip

Index database:

$ ninja_prep_linux Silva119_release/rep_set/97/Silva_119_rep_set97.fna silva119_prep.fa silva119.db --no-spacers 
$ bowtie2-build silva119_prep.fa silva119

Move the taxonomy to working folder:

$ mv Silva119_release/taxonomy/97/taxonomy_97_7_levels.txt silva119.taxonomy
