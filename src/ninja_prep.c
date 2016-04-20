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
#define BUFFER_LEN 100000000
#define SEP "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
#define NINJA_VER "1.5.1"
#define PRINT_USAGE() \
{\
	printf( "\nNINJA Is Not Just Another - OTU Picking Solution v" NINJA_VER "\n"); \
	printf( "Database preparation program. Usage:\n"); \
	printf( "ninja_prep in_refs.fna out_PREFIX [options]\n"); \
	printf("\nINPUT PARAMETERS:\n"); \
	printf( "in_refs.fna: the references with which a concatesome is to be generated.\n"); \
	printf( "    This feeds into your aligner's database builder (bowtie2-build-s).\n"); \
	printf( "\nOUTPUT PARAMETERS:\n"); \
	printf( "out_PREFIX: the name used to make the following output files: \n"); \
	printf( "    PREFIX.fa  (Concatasome for making a BW database)\n"); \
	printf( "    PREFIX.db  (NINJA-OPS' index for the resulting database)\n"); \
	printf( "    PREFIX.tcf (Compressasome for OTU compaction program)\n"); \
	printf( "\nOPTIONS (must be used in order):\n"); \
	printf( "--spacers: Insert ambiguous spacers between references.\n"); \
	printf( "--no-comp: Repress compressasome (PREFIX.tcf) creation.\n"); \
	exit(1); \
}
typedef uint32_t PAYLOAD;
static const PAYLOAD BAD_PAYLOAD = (PAYLOAD)-1;
typedef struct TrieNode TrieNode; 
struct TrieNode {
	uint8_t letter;
	uint32_t depth;
	PAYLOAD content;
	TrieNode *parent;
	TrieNode *Child[256];
};
TrieNode * TN_new() { return calloc(1,sizeof(TrieNode)); }
static inline void TN_add(TrieNode *tn, uint8_t *string, uint32_t len, PAYLOAD x) {
	uint32_t i = 0; for (; i < len; ++i) {
		if (tn->Child[string[i]]) tn = tn->Child[string[i]];
		else break; 
	}
	for (; i < len; ++i) {
		TrieNode *n = calloc(1,sizeof(*n));
		n->letter = string[i];
		n->depth = i+1;
		n->content = BAD_PAYLOAD;
		n->parent = tn;
		tn->Child[string[i]] = n;
		tn = tn->Child[string[i]];
	}
	tn->content = x;
}
static inline PAYLOAD TN_query(TrieNode *tn, uint8_t *string, uint32_t len) {
	for (uint32_t i = 0; i < len; ++i) {
		if (tn->Child[string[i]]) 
			tn = tn->Child[string[i]];
		else return BAD_PAYLOAD; 
	}
	return tn->content;
}
size_t compress(char *in, size_t len, uint16_t **outP) { // make TCF
	TrieNode *tn = TN_new();
	char *b = malloc(1);
	unsigned nxt = 0; for ( ; nxt < 256; ++nxt) 
		*b = nxt, TN_add(tn,b,1,nxt);
	size_t outsz = 1000;
	uint16_t *out = malloc(outsz*sizeof(*out));
	uint8_t *buffer = malloc(BUFFER_LEN+1);
	size_t loc = 0, outix = 0;
	for (size_t i = 0; i < len; ++i) {
		buffer[loc++] = in[i];
#ifdef DEBUG
		if (loc >= BUFFER_LEN) {puts("Error. Exceeded <buffer>."); exit(3);}
#endif
		if ( TN_query(tn,buffer,loc) == BAD_PAYLOAD ) {
			if (nxt < UINT16_MAX) TN_add(tn,buffer,loc, nxt++); 
			if (outix == outsz) out = realloc(out,(outsz *= 2)*sizeof(*out));
			out[outix++] = TN_query(tn,buffer,loc-1);
			*buffer = in[i], loc = 1;
		}
	}
	out[outix++] = TN_query(tn,buffer,loc);
	out = realloc(out,outix*sizeof(*out));
	*outP = out;
	free(b); free(buffer);
	return outix;
}
int main( int argc, char *argv[] )
{
	FILE *in_seqs, *out_DB, *out_map;
	char *inmode = "rb", *outmode = "wb"; 
	int doSpacer = 0, doComp = 1; 
	if ( argc < 3 || argc > 5 ) PRINT_USAGE()
	char *name = malloc(2048), *prefix = argv[2];
	in_seqs = fopen(argv[1], inmode);
	sprintf(name,"%s.fa",prefix), out_DB = fopen(name, outmode);
	sprintf(name,"%s.db",prefix), out_map = fopen(name, outmode);
	if (argc > 3 && !strcmp(argv[argc-1],"--no-comp")) 
		--argc, doComp = 0, puts("Repressing compressasome.");
	if (argc > 3 && !strcmp(argv[argc-1],"--spacers")) 
		--argc, doSpacer = 1, puts("Using internal concatesome spacers.");
	if (!in_seqs || !out_DB || !out_map) {
		fputs("Can't open input/output file(s)! Check parameter order.\n",stderr);
		exit(1);
	}
	char *spacer = doSpacer ? SEP : "\0"; // empty spacer
	char *lineO = malloc(LINELEN+1), *line=lineO; 
	if (!(line = fgets(line,LINELEN,in_seqs))) 
		{ fputs("Bad input.\n",stderr); exit(2); }
	uint64_t len = 0, stLen, spLen = strlen(spacer), oldLen = 0; 
	stLen = strlen(line);
	if (line[stLen-1] == '\n') --stLen;
	if (line[stLen-1] == '\r') --stLen;
	line[stLen] = 0;
	fprintf(out_map,"%s\t",line+1);
	fprintf(out_DB,">NINJA-OPS_%s\n%s",NINJA_VER,spacer);
	while (line = fgets(line,LINELEN,in_seqs)) { // convert fasta
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
	free(lineO); fclose(out_DB); fclose(out_map);
	if (doComp) { // create compressasome
		sprintf(name,"%s.tcf",prefix); FILE *comp_db = fopen(name,outmode);
		sprintf(name,"%s.fa",prefix), out_DB = fopen(name,inmode);
		if (!comp_db || !out_DB) {
			fputs("Failed to create compressasome.\n",stderr); exit(1); }
		fseek(out_DB,0,SEEK_END); size_t sz = ftell(out_DB); rewind(out_DB);
		char *db_txt = malloc(sz+1);
		size_t db_len = fread(db_txt, sz, 1, out_DB);
		db_txt[sz] = 0;
		fclose(out_DB);
		uint16_t *Compressed;
		size_t compLen = compress(db_txt,sz,&Compressed);
		fwrite(Compressed,sizeof(*Compressed),compLen,comp_db);
	}
	puts("Finished.");
	return 0;
}