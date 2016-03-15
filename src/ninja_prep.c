/* NINJA-OPS: NINJA Is Not Just Another - OTU Picking Solution
   Database preparation program. 
   http://ninja-ops.ninja
   This program generates the databases for use with a BWT-enabled aligner.
   
   Compilation information (GCC):
   Ascribes to std=gnu89 multi-platform
   Flags: -m64 -O3 -fwhole-program ninja_prep.c
*/
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#define LINELEN 1000000
#define SEP "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
#define NINJA_VER "1.4.1"
#define PRINT_USAGE() \
{\
	printf( "\nNINJA Is Not Just Another - OTU Picking Solution v" NINJA_VER "\n");\
	printf( "Database preparation program. Usage:\n"); \
	printf( "ninja_prep in_refs.fna out_concatesome.fa out_concatesome.db [options]\n"); \
	printf("\nINPUT PARAMETERS:\n");\
	printf( "in_refs.fa: the references with which a concatesome is to be generated.\n"); \
	printf( "    This feeds into your aligner's database builder: eg bowtie2-build-s.\n"); \
	printf( "\n" "OUTPUT PARAMETERS:\n");\
	printf( "out_concatesome.fa: Concatenated refs for making a BW database.\n");\
	printf( "out_concatesome.db: NINJA-OPS' index for the resulting database.\n");\
	printf( "\nOPTIONS:\n");\
	printf( "--no-spacer: Do not insert spacer between references\n");\
	exit(1);\
}
int main( int argc, char *argv[] )
{
	FILE *in_seqs, *out_DB, *out_map;
	char *inmode = "rb", *outmode = "wb"; 
	int doSpacer = 1, argS = 1; 
	if ( argc < 4 || argc > 5 ) PRINT_USAGE()
	in_seqs = fopen(argv[argS++], inmode);
	out_DB = fopen(argv[argS++], outmode);
	out_map = fopen(argv[argS++], outmode);
	if (argc == 5 && !strcmp(argv[--argc],"--no-spacer")) 
		doSpacer = 0, puts("Not using internal spacers.");
	if (!in_seqs || !out_DB || !out_map) {
		fputs("Can't open input/output file(s)! Check parameter order.\n",stderr);
		exit(1);
	}
	char *spacer = doSpacer ? SEP : "\0";
	char *lineO = malloc(LINELEN+1), *line=lineO; 
	if (!(line = fgets(line,LINELEN,in_seqs))) 
		{ fputs("Bad input.\n",stderr); exit(2); }
	uint64_t len = 0, stLen, spLen = strlen(spacer), oldLen = 0; 
	stLen = strlen(line);
	if (line[stLen-1] == '\n') --stLen;
	if (line[stLen-1] == '\r') --stLen;
	line[stLen] = 0;
	fprintf(out_map,"%s\t",line+1);
	fprintf(out_DB,">NINJA-OPS_140\n%s",spacer);
	while (line = fgets(line,LINELEN,in_seqs)) {
		stLen = strlen(line);
		if (line[stLen-1] == '\n') --stLen;
		if (line[stLen-1] == '\r') --stLen;
		line[stLen] = 0;
		if (*line == '>') { // new header
			if (len == oldLen) continue;
			fprintf(out_DB,"%s%s",spacer,spacer);
			fprintf(out_map,"%llu\n%s\t",oldLen,line+1);
			oldLen = (len += spLen + spLen);
			continue;
		}
		len += fprintf(out_DB,"%s",line); // within sequence
	}
	fprintf(out_DB,"%s\n",spacer); // endcap spacer
	fprintf(out_map,"%llu\n",oldLen); // finalize map
	free(lineO);
	puts("Finished.");
	return 0;
}