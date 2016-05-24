/* NINJA-OPS: NINJA Is Not Just Another - OTU Picking Solution
   Alignment parser and OTU table construction program
   http://ninja-ops.ninja
   This program parses the alignments and generates reports and an OTU table.
   
   Compilation information (GCC):
   Ascribes to std=gnu99 multi-platform
   Flags: -m64 -O3 -ffast-math -std=gnu99 -fwhole-program ninja_prep.c
   
   This program assumes the user or wrapper script has run an aligner (bowtie2):
   bowtie2-align-s --no-head -p4 -x greengenes97 -f reads.fna -S align.sam
   [--np 0 --mp "1,1" --rdg "0,1" --rfg "0,1" --score-min "L,0,-.03" --norc -k 1]
   
   Compilation information/directives (GCC):
   Use:     -m64 -Ofast -fwhole-program parse.c 
   Debug:   -m64 -D DEBUG -D PROFILE -ggdb
   More -D: LOGMATCHES (produces parseLog.txt and map_seqid_reps.txt)
            CHTOL=X (chimera tolerance. Allows reads across refs, X bases)
*/
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define NINJA_VER "1.5.1"
#ifndef CHTOL
	#define CHTOL 0
#endif
#define SHOW_USAGE() {\
	printf( "\nNINJA Is Not Just Another - OTU Picking Solution v" NINJA_VER "\n");\
	printf( "Alignment parser and OTU table construction program. Usage:\n");\
	printf( "ninja_parse in_PREFIX in_aligns.sam in_NINJA.db [in_taxa.txt] [--legacy] [LOG]\n" );\
	printf("\nINPUT PARAMETERS:\n");\
	printf( "in_PREFIX: use the same prefix as with ninja_filter.\n"); \
	printf( "in_aligns: the bowtie2 (headerless, match-only) short read alignment\n");\
	printf( "in_NINJA.db: the (included) sequence index -> OTU database\n");\
	printf( "in_taxa.txt (optional): the (included) sorted OTU -> taxonomy table\n");\
	printf( "\n" "OUTPUT PARAMETERS:\n"); \
	printf( "--legacy: optional, create legacy tab-delimited OTU table (txt).\n");\
	printf( "LOG: Output failures, OTU mappings. ninja_filter with LOG is required.\n"); \
	printf( "Note: OTU table will be output as PREFIX.biom\n"); \
	exit(1); \
}

/** 
 *  Utility functions
 */
 
 // String comparison
inline int xcmp(str1, str2) register const char *str1, *str2; {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}
inline size_t crBST(char *key, size_t sz, char **String) {
	char **p = String;
	while (sz) {
		size_t w = sz >> 1; 
		char *ref_s = *(p+w+1), *key_s = key;
		while (*ref_s == *key_s++) if (!*ref_s++) return p+w+1-String; 
		if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	char *ref_s = *p, *key_s = key;
	while (*ref_s == *key_s++) if (!*ref_s++) return p - String;
	return -1;
	//return p - String
}

// Wide binary search: uncovers range of result rather than single result
unsigned long int uWBS(unsigned long *ixList, unsigned long key, unsigned long range) {
	// wide binary search index list for correct OTU 
	unsigned long middle, low = 0, high = range;
	while (low <= high) {
		middle = low + ((high - low) >> 1);
		if (key > ixList[middle]) low = middle + 1;
		else if (key < ixList[middle]) high = middle - 1;
		else break; 
	}
	if (ixList[middle] > key) --middle;
	return middle;
}

// Opens tax map from file handle with the given delimiter and outputs numbers to 
// Col2 array and strings to Col1. File is closed after running.
unsigned long parse_str_lu_map(FILE * fp, char delim, char *** Col1, unsigned long ** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned long ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); // *ixStrp = ixStr; 
	if (!ixStr) { puts("PSM out of memory.\n"); return 0; }
	fread(ixStr, 1, ixSize, fp); fclose(fp);
	ixStr[ixSize] = 0;
	unsigned long ilines = 0; char *iix = ixStr;
	while (*iix) ilines+=*iix++=='\n';
	
	// Grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines * sizeof(char *));
	*Col2 = malloc(ilines * sizeof(unsigned long)); 
	if (!*Col1 || !*Col2) { puts("PSL out of memory.\n"); return 0; }
	char **Col1p = *Col1; unsigned long *Col2p = *Col2;
	
	int which = 1; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(10000), *bufp = buffer;
	iix = ixStr;
	char c; while ((c=*iix++)) {
		if (c=='\n' || c==delim) { // commit, switch buffer pointer
			memset(bufp,'\0',1);
			if (!which) *Col2p++ = (unsigned long)atol(buffer);
			else {
				if (*(bufp-1)=='\r') memset(bufp-1,'\0',1);
				int amt; *Col1p = malloc(amt = strlen(buffer)+1);
				if (!*Col1p) return 0; 
				strncpy(*Col1p++, buffer, amt);
			}
			which ^= 1; bufp = buffer; 
		}
		else *bufp++ = c;
	}
	free(ixStr);
	return ilines;
}

// Similar to the above, but both columns are strings
unsigned long parse_str_str_map(FILE * fp, char delim, char *** Col1, char *** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned long ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); // *ixStrp = ixStr; 
	if (!ixStr) { puts("PSM out of memory.\n"); return 0; }
	fread(ixStr, 1, ixSize, fp); fclose(fp);
	ixStr[ixSize] = 0;
	unsigned long ilines = 0; char *iix = ixStr;
	while (*iix) ilines+=*iix++=='\n';
	
	// Grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines*sizeof(char*));
	*Col2 = malloc(ilines*sizeof(char*)); 
	if (!*Col1 || !*Col2) { puts("PSS out of memory.\n"); return 0; }
	char **Col1p = *Col1, **Col2p = *Col2;
	
	int which = 0; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(10000), *bufp = buffer;
	iix = ixStr;
	char c; while ((c=*iix++)) {
		if (c=='\n' || c==delim) { // commit, switch buffer pointer
			memset(bufp,'\0',1);
			if (!which) {
				if (*(bufp-1)=='\r') memset(bufp-1,'\0',1);
				int amt; *Col1p = malloc(amt = strlen(buffer)+1);
				if (!*Col1p) return 0; 
				strncpy(*Col1p++, buffer, amt);
			}
			else {
				if (*(bufp-1)=='\r') memset(bufp-1,'\0',1);
				int amt; *Col2p = malloc(amt = strlen(buffer)+1);
				if (!*Col2p) return 0; 
				strncpy(*Col2p++, buffer, amt);
			}
			
			which ^= 1; bufp = buffer; 
		}
		else *bufp++ = c;
	}
	free(ixStr);
	return ilines;
}

// Opens a file handle (from fopen) and outputs the string to the given address. 
// Returns number of lines. 
unsigned long parse_strings(FILE *fp, char *** Strings) {
	fseek(fp, 0, SEEK_END); unsigned long ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *Dump = malloc(ixSize + 1); 
	if (!Dump) { puts("PS out of memory.\n"); return 0; }
	fread(Dump, 1, ixSize, fp); fclose(fp);
	Dump[ixSize] = 0;
	unsigned long ilines = 0, curLen = 1000; 
	*Strings = malloc(curLen * sizeof(char **));
	char *iix = Dump - 1, *bufP, *Buffer, // = malloc(1000000), *bufP = Buffer,
		**StrP = *Strings;
	while (*++iix) {
		bufP = iix;
		while (*++iix != '\n'); //skips empty lines
		*StrP = malloc(iix - bufP + 1);
		Buffer = *StrP++; 
		do *Buffer++ = *bufP; while (++bufP < iix);
		memset(Buffer,'\0',1);
		if (++ilines == curLen) {
			unsigned long offset = StrP - *Strings;
			*Strings = realloc(*Strings, (curLen *= 2) * sizeof(char **));
			if (!*Strings) { printf("out of D-memory.\n"); return 1; }
			StrP = *Strings + offset;
		}
	}
	*Strings = realloc(*Strings, ilines * sizeof( char **));
	free(Dump);
	return ilines;
}

#ifndef min 
	#define min(a, b) ((a)<=(b) ? (a) : (b)) 
#endif
#define ch(i) *(**(a+i) + depth) 
#define med3(ia, ib, ic) med3func(a, ia, ib, ic, depth)
#define CUTOFF 10
#define MEDCUT 50
// Swaps two characters in a vector
inline void swap(char ***a, int i, int j) 
	{ char **t = *(a+i); *(a+i) = *(a+j); *(a+j) = t; }
inline void vecswap(char ***a, int i, int j, int n) 
	{ while (n-->0) swap(a, i++, j++); }
// Returns median of ints, used in twrqs
inline int med3func(char ***a, int ia, int ib, int ic, int depth) {
	int va, vb, vc;
	if ((va=ch(ia)) == (vb=ch(ib))) return ia;
	if ((vc=ch(ic)) == va || vc == vb) return ic;
	return va < vb ?
		(vb < vc ? ib : (va < vc ? ic : ia ) ) : 
		(vb > vc ? ib : (va < vc ? ia : ic ) ); 
} 
// Insertion sort delegated to by twrqs
inline void inssort(char ***a, int n, int depth) {
	int i, j;
	for (i = 1; i < n; i++) for (j = i; j > 0; j--) {
		if (strcmp(**(a+j-1)+depth, **(a+j)+depth) <= 0) break;
		swap(a, j, j-1);
	} 
}  
// 3-way Radix Quicksort 
void twrqs(char ***a, unsigned n, int depth) {
	if (n < CUTOFF) { inssort(a, n, depth); return; }
	unsigned pl = 0, pm = n >> 1, d;
	int le, lt, gt, ge, r, v, pn = n-1;
	// if large enough, get median of median
	if (n > MEDCUT) {
		d = n >> 3;
		pl = med3(pl, pl+d, pl + (d << 1));
		pm = med3(pm-d, pm, pm+d);
		pn = med3(pn - (d << 1), pn-d, pn);
	}
	pm = med3(pl, pm, pn);
	swap(a, 0, pm);
	v = ch(0); // grab first letter
	for (le = 1; le < n && ch(le) == v; le++);  
	if (le == n) {
		if (v != 0) twrqs(a, n, depth+1);
		return;
	}
	lt = le; gt = ge = n-1;
	// core QS module; partition the data recursively
	for (;;) {
		for ( ; lt <= gt && ch(lt) <= v; lt++)
			if (ch(lt) == v) swap(a, le++, lt);
		for ( ; lt <= gt && ch(gt) >= v; gt--)
			if (ch(gt) == v) swap(a, gt, ge--);
		if (lt > gt) break;
		swap(a, lt++, gt--);
	}
		r = min(le, lt-le); 
		vecswap(a, 0, lt-r, r); 
		r = min(ge-gt, n-ge-1);
		vecswap(a, lt, n-r, r);
		twrqs(a, lt-le, depth);
		if (v != 0) twrqs(a + lt-le, le + n-ge-1, depth+1);
		twrqs(a + n-(ge-gt), ge-gt, depth); 
}

/** 
 *  Main parsing algorithm.
 *  INPUT     prefix:           db file output from ninja_filter
 *            alignmentsFile:   alignment file output from bowtie2
 *            masterDBFile:     master db file packaged with ninja
 *            taxMapFile:       reference taxonomy map packaged with NINJA-OPS
 **/
int main ( int argc, char *argv[] )
{
	// Starts timer for profiling
	clock_t start;
	start = clock();
	// Checks supported number of args specified
	if ( argc < 4 || argc > 7 ) SHOW_USAGE();
    int legacy = 0, doLog = 0;
	if (!strcmp(argv[argc-1],"LOG")) {
		printf("Logging enabled.\n");
		doLog = 1;
		--argc;
	}
	if (!strcmp(argv[argc-1],"--legacy")) {
		printf("Legacy output mode toggled.\n");
		legacy = 1;
		--argc;
	}
	// Controls whether taxonomy will be considered
    int doTaxmap = (argc == 5) ?: 0;
	// Assumes argv[n] are filenames to open, in the order specified.
	char *prefixStr = argv[1];
	char *db_sx = ".db", *dp_sx = "_dupes.txt", 
		 *tab_sx = legacy? "_otutable.txt" : "_otutable.biom",
		 *lf_sx = "_fail.txt", *lp_sx = "_pass.log", *om_sx = "_otumap.txt";
	char *inputDB = calloc(1,1+strlen(prefixStr)+strlen(db_sx)),
		 *inputDP = calloc(1,1+strlen(prefixStr)+strlen(dp_sx)),
		 *outTable= calloc(1,1+strlen(prefixStr)+strlen(tab_sx)),
		 *outLogF = calloc(1,1+strlen(prefixStr)+strlen(lf_sx)),
		 *outLogP = calloc(1,1+strlen(prefixStr)+strlen(lp_sx)),
		 *outLogM = calloc(1,1+strlen(prefixStr)+strlen(om_sx));
	strcpy(inputDB,prefixStr); strcpy(inputDB+strlen(prefixStr),db_sx);
	strcpy(inputDP,prefixStr); strcpy(inputDP+strlen(prefixStr),dp_sx);
	strcpy(outTable,prefixStr); strcpy(outTable+strlen(prefixStr),tab_sx);
	strcpy(outLogF,prefixStr); strcpy(outLogF+strlen(prefixStr),lf_sx);
	strcpy(outLogP,prefixStr); strcpy(outLogP+strlen(prefixStr),lp_sx);
	strcpy(outLogM,prefixStr); strcpy(outLogM+strlen(prefixStr),om_sx);
	
	FILE *mfp = fopen( inputDB, "rb" );
	FILE *ifp = fopen( argv[2], "rb" ), *ifi = fopen( argv[3], "rb"), 
		 *ofp = fopen( outTable, "wb"), 
		 *itx = 0, *logFail = 0, *logPass = 0, *logMap = 0, *inDupes = 0;
	if (doTaxmap) itx = fopen( argv[4], "rb"); // otu-tax map provided
	printf("Opened %s for OTU Table writing\n", outTable);
	if (doLog) {
		puts("Logging enabled.");
		logFail = fopen(outLogF, "wb");
		inDupes = fopen(inputDP,"rb");
		logPass = fopen(outLogP, "wb");
		logMap = fopen(outLogM, "wb");
		if (!logFail || !inDupes || !logPass || !logMap) {
			puts("Couldn't initialize logging."); 
			puts("Note: LOG is also required during filter with same prefix.");
			exit(2); 
		}
	}
	if ( !mfp || !ifp || !ifi || !ofp || (doTaxmap && itx == 0)) {
		fprintf(stderr, "Could not open one or more files.\n");
		SHOW_USAGE();
	}
	// Parse the loaded files: sample file is "pre-parsed" in raw string form
	char **OtuList, **OtuMap_otusPre, **OtuMap_taxaPre, **SampDBdump, **AllSamps;
	unsigned long *ixList, 
		ilines = parse_str_lu_map(ifi, '\t', &OtuList, &ixList),
		blines = doTaxmap ? parse_str_str_map(itx, '\t', &OtuMap_otusPre, &OtuMap_taxaPre) : 0;
	char **OtuMap_otus = malloc(blines*sizeof(*OtuMap_otus)), 
		 **OtuMap_taxa = malloc(blines*sizeof(*OtuMap_taxa));
	printf("Total OTUs available to pick from: %lu\n", ilines);
	
	if (doTaxmap) {
		if (!OtuMap_otus || !OtuMap_taxa) {fputs("Error: memory.\n",stderr); exit(3);}
		char ***OtuMap_otusP = malloc(blines*sizeof(*OtuMap_otusP)); 
		if (!OtuMap_otusP) {fputs("Error: memory.\n",stderr); exit(3);}
		for (size_t i = 0; i < blines; ++i) OtuMap_otusP[i] = OtuMap_otusPre+i; 
		//inline int strPcmp(const void *a, const void *b) { return strcmp(**(char ***)a,**(char ***)b); }
		//qsort(OtuMap_otusP, blines, sizeof(*OtuMap_otusP), strPcmp);
		twrqs(OtuMap_otusP,blines,0);
		for (size_t i = 0; i < blines; ++i) OtuMap_otus[i] = *OtuMap_otusP[i], 
			OtuMap_taxa[i] = OtuMap_taxaPre[OtuMap_otusP[i]-OtuMap_otusPre];
		free(OtuMap_otusPre); free(OtuMap_taxaPre); free(OtuMap_otusP); 
		//for (size_t i=0; i<blines; ++i) printf("%s\t%s\n",OtuMap_otus[i],OtuMap_taxa[i]);
	}
	
	// Parse samples from pre-parsed sample file
	unsigned long slines = parse_strings(mfp, &SampDBdump);
	unsigned long flines = 0;
	if (doLog) {
		flines = parse_strings(inDupes, &AllSamps);
		unsigned long z = 0; for (; z < flines; ++z) {
			char *entry = AllSamps[z] - 1;
			while (*++entry) if (*entry == '\t') *entry = '\n';
		}
	}
	if (!ilines || (doTaxmap && !blines) || !slines) 
		{ printf("Unparsable: ilines %lu, blines %lu, slines %lu.\n",ilines, blines, slines); return 1; }
#ifdef PROFILE
	printf("->Time for list parse: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	unsigned long numSamps = atol(*SampDBdump++), numReads = slines - numSamps - 1; // max reads possible
	char **Seq2samp = SampDBdump + numSamps;
	printf("Number of unique samples: %lu, max reads: %lu\n", numSamps, numReads);
    // Processes SAM file generated from bowtie
	fseeko(ifp, 0, SEEK_END); size_t fsize = ftello(ifp); fseeko(ifp, 0, SEEK_SET); 

	// Create OTU table counts matrix
	unsigned long *OtuTable = calloc(numSamps * ilines, sizeof (unsigned long)), otuIX;
	if (!OtuTable) {printf("Couldn't allocate OTU table.\n"); return 1;}
	unsigned long tabs, rix, six, cnt, alignPos, totalReadsCnt = 0, chimeras = 0;
	//char *cix = string - 1, *startS, 
	char *curSamp; 
// Lots each match and duplication count to a separate, fixed filename in the directory
#ifdef LOGMATCHES
	FILE *log = fopen("parseLog.txt", "wb");
	FILE *log2 = fopen("map_seqid_reps.txt", "wb");
#endif
	// New parser loop
	char *line = malloc(UINT16_MAX+1), *lineO = line;
	while (line = fgets(line,UINT16_MAX+1,ifp)) {
		char *cix = line, *startS = line;
		tabs = 0; do if (*++cix == '\t') ++tabs; while (tabs < 2);
		if (*++cix == '*') {
			if (doLog) fprintf(logFail,"%s",AllSamps[atol(startS)]);
			while (*++cix != '\n')
			; continue; // skip to next
		}
		do if (*++cix == '\t') ++tabs; while (tabs < 3); // one more tab
		alignPos = atol(cix);
		rix = atol(startS);
		curSamp = *(Seq2samp + rix); // look up rix
		
		otuIX = uWBS(ixList, alignPos, ilines); 
		// slays chimeras
		do if (*++cix == '\t') ++tabs; while (tabs < 9); // proper SAM format required!
		char *beginSeq = cix + 1;
		while (*++cix != '\t');
		if (otuIX < ilines && ((cix-beginSeq) > *(ixList+otuIX+1) - alignPos + CHTOL)) { 
			if (doLog) fprintf(logFail,"%s",AllSamps[rix]);
			++chimeras;
			while (*++cix != '\n')
			; continue; 
		}
		if (doLog) {
			char *dupes = AllSamps[rix];
			while (*dupes) {
				char *begin = dupes;
				while (*dupes != '\n') ++dupes;
				*dupes = '\0';
				fprintf(logPass,"%s\t%s\n", begin, OtuList[otuIX]);
				*dupes++ = '\t';
			}
		}
#ifdef LOGMATCHES
		fprintf(log,"%lu\t%s\t%s\n", rix, OtuList[otuIX],!doTaxmap ? "" : *(OtuMap_taxa + crBST(*(OtuList + otuIX), blines-1, OtuMap_otus))); //deleteme
#endif
		unsigned long amt = 0;
		do { //startS new scope
			startS = curSamp;
			while (*++curSamp != ':');
			six = atol(startS);
			startS = ++curSamp; // + 1;
			while (*++curSamp != ':');
			cnt = atol(startS);
			amt += cnt;
			*(OtuTable + otuIX * numSamps + six) += cnt;
			totalReadsCnt += cnt;
		} while (*++curSamp);
#ifdef LOGMATCHES
		fprintf(log2, "%lu\t%lu\n", rix, amt); 
#endif
	}
	free(lineO);
#ifdef PROFILE
	printf("->Time for matrix generation: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	printf("Total chimeric alignments: %lu\n",chimeras);
	printf("Total reads expanded: %lu\n", totalReadsCnt);
	
	// Legacy table format output handler
	if (legacy) {
		// Prints sample and taxonomy header
		fprintf(ofp, "#OTU ID");
		char **SampP = SampDBdump - 1;
		char **OtuP = OtuList - 1;
		cnt = numSamps; do fprintf(ofp,"\t%s",*++SampP); while (--cnt);
		if (doTaxmap) fprintf(ofp,"\ttaxonomy");
		unsigned long i, *row, *rowP;
		// Print taxonomy values
		for (i = 0; i < ilines; i++) {
			// Screens if any numbers in this row
			row = OtuTable + i*numSamps; rowP = row;
			cnt = numSamps; do if (*rowP++) break; while (--cnt);
			++OtuP;
			if (cnt) {
				fprintf(ofp, "\n%s", *OtuP);
				rowP = row; cnt = numSamps; 
				do fprintf(ofp, "\t%lu", *rowP++); while (--cnt);
				if (doTaxmap) {
					size_t mp = crBST(OtuList[i],blines-1,OtuMap_otus);
					fprintf(ofp,"\t%s", mp != (size_t)-1 ? 
						*(OtuMap_taxa + mp) : "UNSUPPORTED_MAPPING"); 
				}
			}
		}
		fputc('\n',ofp);
	}
	// BIOM 1.0 (sparse) format output handler
	else { 
		unsigned long *imap = malloc(ilines*sizeof(*imap));
		unsigned long i, j, curIx = 0; // *row, *rowP;
		for (i = 0; i < ilines; i++) { 
			imap[i] = -1;
			for (j = 0; j < numSamps; ++j) 
				if (OtuTable[i*numSamps + j]) {imap[i] = curIx++; break;}
		}
		time_t t = time(NULL);
		struct tm tm = *localtime(&t);
		fprintf(ofp, "{\n\"id\":null,\n\"format\": \"Biological Observation Matrix 1.0.0\",\n"
		"\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",\n"
		"\"type\": \"OTU table\",\n\"generated_by\": \"NINJA v" NINJA_VER "\",\n"
		"\"date\": \"%d-%s%d-%s%dT%s%d:%s%d:%s%d\",\n", tm.tm_year + 1900, tm.tm_mon + 1 < 10 ? "0" : "", 
		tm.tm_mon + 1, tm.tm_mday < 10 ? "0" : "", tm.tm_mday, 
		tm.tm_hour < 10 ? "0" : "", tm.tm_hour, tm.tm_min < 10 ? "0" : "", tm.tm_min,
		tm.tm_sec < 10 ? "0" : "", tm.tm_sec);

		// Write the rows in the biom format
		fprintf(ofp, "\"rows\":[");
		char **OtuP = OtuList - 1; size_t pcnt = 0;
		int ix = ilines; do { // display the "row" lines
			if (imap[ilines - ix] == -1) {++OtuP; continue;}
			fprintf(ofp,"\n\t{\"id\":\"%s\", \"metadata\":",*++OtuP);
			++pcnt;
			if (doTaxmap) {
				fprintf(ofp,"{\"taxonomy\":[");
				char *taxon, *tP;
				size_t mp = crBST(OtuList[ilines-ix],blines-1,OtuMap_otus);
				if (mp != (size_t)-1) {
					taxon = *(OtuMap_taxa + mp), tP = taxon;
					for (;;) {
						while (*taxon && *taxon != ';') ++taxon;
						if (!*taxon) break;
						*taxon = 0;
						fprintf(ofp, "\"%s\", ", tP);
						*taxon++ = ';';       // delimiter exclusion
						tP = taxon;
						if (*tP == ' ') ++tP; // space exclusion
					}
				} else tP = "UNSUPPORTED_MAPPING";
				fprintf(ofp,"\"%s\"]}}", tP);
			}
			else fprintf(ofp,"null}");
			fprintf(ofp,",");
		} while (--ix);
		fseek(ofp, -1, SEEK_CUR);
		fprintf(ofp, "\n],\n");
		
		// Write sparse columns
		char **SampP = SampDBdump - 1;
		fprintf(ofp, "\"columns\": [");
		ix = numSamps; do 
			fprintf(ofp, "\n\t{\"id\":\"%s\", \"metadata\":null},", *++SampP);
		while (--ix);
		fseek(ofp, -1, SEEK_CUR);
		fprintf(ofp, "\n],\n");
		
		// Write structure, data
		fprintf(ofp,"\"matrix_type\": \"sparse\",\n\"matrix_element_type\": \"int\",\n"
			"\"shape\": [%lu, %lu],\n", pcnt, numSamps);
		fprintf(ofp,"\"data\":[");
		
		// loop thru all points
		for (i = 0; i < ilines; i++) for (j = 0; j < numSamps; j++) {
			if (OtuTable[i * numSamps + j]) 
				fprintf(ofp, "[%lu,%lu,%lu],\n\t", imap[i], j, OtuTable[i * numSamps + j]);
		}
		fseek(ofp, -3, SEEK_CUR);
		fprintf(ofp,"]\n}");
		free(imap);
	}
	if (doLog) {
		fclose(logPass); fclose(logFail);
		if (doTaxmap) free(OtuMap_otus), free(OtuMap_taxa);
		free(OtuTable);
		logPass = fopen(outLogP, "rb");
		if (!logPass) {fprintf(stderr,"Couldn't open logPass\n"); exit(1);}
		char **queriesPre, **refsPre;
		unsigned long xlines = 
			parse_str_str_map(logPass, '\t', &queriesPre, &refsPre);
		char ***refsP = malloc(xlines*sizeof(*refsP)); 
		for (size_t i = 0; i < xlines; ++i) refsP[i] = refsPre+i; 
		twrqs(refsP,xlines,0);
		char **queries = malloc(xlines*sizeof(*queries)), **refs = malloc(xlines*sizeof(*refs));
		for (size_t i = 0; i < xlines; ++i) refs[i] = *refsP[i], 
			queries[i] = queriesPre[refsP[i]-refsPre];
		free(refsPre); free(queriesPre); free(refsP);
		// For each same ref, spit out all queries with tab after it
		char *lastRef = *refs;
		fprintf(logMap,"%s\t%s\t",lastRef,*queries);
		for (size_t i = 1; i < xlines; ++i) {
			if (!strcmp(lastRef,refs[i]))
				fprintf(logMap,"%s\t",queries[i]);
			else fprintf(logMap,"\n%s\t%s\t",(lastRef = refs[i]),queries[i]);
		}
		fputc('\n',logMap);
	}
#ifdef PROFILE
	printf("->Time for filtering and table write: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	printf("Run complete.\n"); 
	return 0;
}
