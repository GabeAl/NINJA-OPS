/* NINJA OTU Picker: NINJA Is Not Just Another OTU Picker: parser program
   OTU table parser, mapper, generator. Uses B-W LF mapping search to
   quickly align each read to the appended DB of all OTU sequences.
   This program assumes the user or companion script has run bowtie2:
   [bowtie2-build-s -o3 ninjaDB.fasta Ninja97] or use included DB
   bowtie2-align-s --no-head --no-unal -o3 -p4 -f reads.fna -x Ninja97 -S align.txt
   [--mp "1,1" --rdg "1,1" --rfg "1,1" --score-min "L,0,-.03" -k 1 --norc --fast]
   
   Compilation information (GCC):
   Ascribes to std=gnu and doesn't require C99 support or UNIX intrinsics.
   Use:     -m64 -Ofast -fwhole-program parse.c -fprofile-generate [-fprofile-use]
   Profile: -m64 -Ofast -D PROFILE -flto parse.c
   Debug:   -m64 -D DEBUG -D PROFILE -ggdb
   More -D: LOGMATCHES (produces parseLog.txt and map_seqid_reps.txt)
*/
   
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define ARRSZ 16
#define NINJA_VER "1.3"
#define SHOW_USAGE() {\
	printf( "\nNINJA Is Not Just Another - OTU Picking Solution v" NINJA_VER ": Parser. Usage:\n");\
	printf( "ninja_parse in_PREFIX in_aligns.sam in_NINJA.map [in_taxa.txt] [--legacy] [LOG]\n" );\
	printf("\nINPUT PARAMETERS:\n");\
	printf( "in_PREFIX: use the same prefix as with ninja_filter.\n"); \
	printf( "in_aligns: the bowtie2 (headerless, match-only) short read alignment\n");\
	printf( "in_NINJA.map: the (included) sequence index -> OTU database\n");\
	printf( "in_taxa.txt (optional): the (included) sorted OTU -> taxonomy table\n");\
	printf( "\n" "OUTPUT PARAMETERS:\n"); \
	printf( "--legacy: optional, create legacy tab-delimited OTU table (txt).\n");\
	printf( "LOG: optional, outputs failures. ninja_filter with LOG is required.\n"); \
	printf( "Note: OTU table will be output as PREFIX.biom\n"); \
	exit(1); \
}

/** 
 *  Utility functions
 */
 
 // String comparison functions. The first returns -1, 0, and 1 only as <, ==, >
 // Second reports magnitude of difference from point of divergence
inline int ycmp(st1p, st2p) register const char *st1p, *st2p; { 
	while (*st1p && *st2p && *st1p++ == *st2p++); 
	return *--st1p<*--st2p?-1:(*st1p>*st2p); 
}
inline int xcmp(str1, str2) register const char *str1, *str2; {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}

// Resizes array without realloc: unsigned (top) and string (bottom)
inline int res_unsigned_arr(unsigned **arr, unsigned oldS, unsigned newS) {
	unsigned *newArray = malloc(newS * sizeof(unsigned int)), 
		*nap = newArray, *oap = *arr;
	if (!newArray) { printf("Out of memory.\n"); return 0; }
	int j=-1; while (++j<=oldS) *nap++ = *oap++;
	unsigned * addr = *arr; free(*arr); 
	return (*arr = newArray) - addr;
}
inline int res_str_arr(char ***arr, unsigned oldS, unsigned newS) {
	char **newArray = malloc(newS * sizeof(char *)),
		**nap = newArray, **oap = *arr;
	if (!newArray) { printf("Out of memory.\n"); return 0; }
	int j=-1; while (++j<=oldS) *nap++ = *oap++;
	char ** addr = *arr; free(*arr); 
	return (*arr = newArray) - addr;
}

// Bubbling comparison of characters
int cmpS(const void *v1, const void *v2) {
	const char i1 = **(const char **)v1;
	const char i2 = **(const char **)v2;
	return i1<i2?-1:(i1>i2);
}
// Bubbling comparison of integers
int cpx(const void *v1, const void *v2) {
	const int i1 = *(const int *)v1;
	const int i2 = *(const int *)v2;
	return i1<i2?-1:(i1 > i2);
}
// Wrapper for xcmp
int cmp(const void *v1, const void *v2) {
	const char *i1 = **(const char ***)v1;
	const char *i2 = **(const char ***)v2;
	return xcmp(i1,i2);
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
// Col1 array and strings to Col2. File is closed after running.
unsigned long parse_string_map(FILE * fp, char delim, unsigned long ** Col1, char *** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned long ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); // *ixStrp = ixStr; 
	if (!ixStr) { puts("PSM out of memory.\n"); return 0; }
	fread(ixStr, ixSize, 1, fp); fclose(fp);
	unsigned long ilines = 0; char *iix = ixStr - 1, iptr;
	while ((iptr = *++iix)) iptr != '\n' ?: ++ilines, (void)0;
	
	// Grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines * sizeof(unsigned long)); 
	*Col2 = malloc(ilines * sizeof(char *));
	if (!*Col1 || !*Col2) { puts("PSM out of memory.\n"); return 0; }
	unsigned long *Col1p = *Col1; char **Col2p = *Col2;
		
	int which = 0; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(1024), *bufp = buffer;
	iix = ixStr;
	char c; while ((c=*iix++)) {
		if (c=='\n' || c==delim) { // commit, switch buffer pointer
			memset(bufp,'\0',1);
			if (!which) *Col1p++ = atol(buffer);
			else {
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

// Opens int map from file handle with the given delimiter and first col numbers to 
// Col1 array and second column numbers to Col2. File is closed after running.
unsigned long parse_unsigned_map(FILE * fp, char delim, unsigned long ** Col1, unsigned long ** Col2) {
	// treat the file to a full read-in and caching
	fseek(fp, 0, SEEK_END); unsigned long ixSize = ftell(fp); fseek(fp, 0, SEEK_SET); 
	char *ixStr = malloc(ixSize + 1); 
	if (!ixStr) { puts("PUM out of memory.\n"); return 0; }
	fread(ixStr, ixSize, 1, fp); fclose(fp);
	unsigned long ilines = 0; char *iix = ixStr - sizeof(char), iptr;
	while ((iptr = *++iix)) iptr != '\n' ?: ++ilines;
	
	// grab values in 'Col1[delim]Col2[\n]' format
	*Col1 = malloc(ilines * sizeof(unsigned long)); 
	*Col2 = malloc(ilines * sizeof(unsigned long));
	if (!*Col1 || !*Col2) { puts("PUM out of memory.\n"); return 0; }
	unsigned long *Col1p = *Col1, *Col2p = *Col2;
		 
	int which = 0; // here 0 is the Col1 field and !0 is the Col2 field
	char *buffer = malloc(32), *bufp = buffer;
	iix = ixStr;
	char c; while ((c=*iix++)) {
		if (c=='\n' || c==delim) { // commit, switch buffer pointer
			memset(bufp,'\0',1);
			which ? (*Col2p++ = atol(buffer)) : (*Col1p++ = atol(buffer));
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
	fread(Dump, ixSize, 1, fp); fclose(fp);
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
		 *lf_sx = "_fail.txt", *lp_sx = "_pass.log";
	char *inputDB = calloc(1,1+strlen(prefixStr)+strlen(db_sx)),
		 *inputDP = calloc(1,1+strlen(prefixStr)+strlen(dp_sx)),
		 *outTable= calloc(1,1+strlen(prefixStr)+strlen(tab_sx)),
		 *outLogF = calloc(1,1+strlen(prefixStr)+strlen(lf_sx)),
		 *outLogP = calloc(1,1+strlen(prefixStr)+strlen(lp_sx));
	strcpy(inputDB,prefixStr); strcpy(inputDB+strlen(prefixStr),db_sx);
	strcpy(inputDP,prefixStr); strcpy(inputDP+strlen(prefixStr),dp_sx);
	strcpy(outTable,prefixStr); strcpy(outTable+strlen(prefixStr),tab_sx);
	strcpy(outLogF,prefixStr); strcpy(outLogF+strlen(prefixStr),lf_sx);
	strcpy(outLogP,prefixStr); strcpy(outLogP+strlen(prefixStr),lp_sx);
	
	FILE *mfp = fopen( inputDB, "rb" );
	FILE *ifp = fopen( argv[2], "rb" ), *ifi = fopen( argv[3], "rb"), 
		 *ofp = fopen( outTable, "wb"), 
		 *itx = 0, *logFail = 0, *logPass = 0, *inDupes = 0;
	if (doTaxmap) itx = fopen( argv[4], "rb"); // otu-tax map provided
	printf("Opened %s for OTU Table writing\n", outTable);
	if (doLog) {
		logFail = fopen(outLogF, "wb");
		inDupes = fopen(inputDP,"rb");
		//logPass = fopen(outLogP, "wb");
		if (!logFail || !inDupes) { // || !logPass) 
			puts("Couldn't create output logs."); 
			puts("LOG is also required during filter with same prefix.");
			exit(2); 
		}
	}
	if ( !mfp || !ifp || !ifi || !ofp || (doTaxmap && itx == 0)) {
		fprintf(stderr, "Could not open one or more files.\n");
		SHOW_USAGE();
	}
	// Parse the loaded files: sample file is "pre-parsed" in raw string form
	unsigned long *OtuList, *ixList,
		ilines = parse_unsigned_map(ifi, ',', &ixList, &OtuList);
	char **OtuMap_taxa, **SampDBdump, **AllSamps; 
	unsigned long *OtuMap_otus, blines = doTaxmap ? 
		parse_string_map(itx, '\t', &OtuMap_otus, &OtuMap_taxa) : 0;
	printf("Total OTUs available: %lu\n", ilines);
	
	// Parse samples from pre-parsed sample file
	unsigned long slines = parse_strings(mfp, &SampDBdump);
	unsigned long flines = 0;
	if (doLog) {
		//printf("Trying to read inDupes: %s\n",inputDP);
		flines = parse_strings(inDupes, &AllSamps);
		unsigned long z = 0; for (; z < flines; ++z) {
			char *entry = AllSamps[z] - 1;
			while (*++entry) if (*entry == '\t') *entry = '\n';
		}
	}
	//for (int i = 0; i < flines; ++i) fprintf(logFail,"%s",AllSamps[i]);
	if (!ilines || (doTaxmap && !blines) || !slines) 
		{ printf("Unparsable: ilines %lu, blines %lu, slines %lu.\n",ilines, blines, slines); return 1; }
#ifdef PROFILE
	printf("->Time for list parse: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	unsigned long numSamps = atol(*SampDBdump++), numReads = slines - numSamps - 1; // max reads possible
	char **Seq2samp = SampDBdump + numSamps;
	printf("Number of unique samples: %lu, max reads: %lu\n", numSamps, numReads);
    // Processes SAM file generated from bowtie
	fseek(ifp, 0, SEEK_END); size_t fsize = ftell(ifp); fseek(ifp, 0, SEEK_SET); 
	//size_t sz = lseek(fileno(fp), 0, SEEK_END); rewind(fp);
	char *string = malloc(fsize + 1); // Allocates a block of memory (direct bytes)
	if (string == NULL) {
		fprintf(stderr, "Insufficient memory for caching input file.\n");
		return 1;
	}
	fread(string, fsize, 1, ifp); //read into string: elems of fsize bytes, 1 elem, using ifp pointer
	fclose(ifp); // close the file 
	
#ifdef PROFILE
	printf("->Time for read-in: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// Create OTU table counts matrix
	unsigned long *OtuTable = calloc(numSamps * ilines, sizeof (unsigned long)), otuIX;
	if (!OtuTable) {printf("couldn't allocate OTU table.\n"); return 1;}
	unsigned long tabs, rix, six, cnt, alignPos, totalReadsCnt = 0;
	char *cix = string - 1, *startS, *curSamp; 
// Lots each match and duplication count to a separate, fixed filename in the directory
#ifdef LOGMATCHES
	FILE *log = fopen("parseLog.txt", "wb");
	FILE *log2 = fopen("map_seqid_reps.txt", "wb");
#endif

	unsigned long lcounter = 0;
    //puts("I'm not a seg fault!");
	while (*++cix) { // != '\t' && *cix != '\n') { // work through the sample string
        if (++lcounter >= numReads) break; 
		startS = cix;
		// lookahead to the read map region
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
		
		otuIX = uWBS(ixList, alignPos, ilines); //*(OtuList + otuIX) is actual otu
		
		// slays chimeras
		do if (*++cix == '\t') ++tabs; while (tabs < 9); // proper SAM format required!
		char *beginSeq = cix + 1;
		while (*++cix != '\t');
		if (otuIX < ilines && ((cix-beginSeq) > (*(ixList+otuIX+1) - alignPos))) { 
			if (doLog) fprintf(logFail,"%s",AllSamps[rix]);
			while (*++cix != '\n')
			; continue; 
		}

#ifdef LOGMATCHES
		fprintf(log,"%lu\t%lu\t%s\n", rix, OtuList[otuIX],!doTaxmap ? "" : *(OtuMap_taxa + uWBS(OtuMap_otus, *(OtuList + otuIX), blines))); //deleteme
#endif
		unsigned long amt = 0;
		//fprintf(ofp,"-- This read maps to otu %u, (%s)\n", OtuList[otuIX], OtuMap_taxa[uWBS(OtuMap_otus,OtuList[otuIX],blines)]); //*(OtuList + otu));
		//printf("rix=%u, curSamp=%u, alignPos=%u,cix=%u\n",rix,curSamp,alignPos,cix);
		do { //startS new scope
			//++amt;
			startS = curSamp;
			while (*++curSamp != ':');
			six = atol(startS);
			startS = ++curSamp; // + 1;
			while (*++curSamp != ':');
			cnt = atol(startS);
			amt += cnt;
			//fprintf(ofp,"-- Found %u copies in sample id %u (%s).\n", cnt, six, SampDBdump[six] );
			*(OtuTable + otuIX * numSamps + six) += cnt;
			totalReadsCnt += cnt;
		} while (*++curSamp);
#ifdef LOGMATCHES
		fprintf(log2, "%lu\t%lu\n", rix, amt); 
#endif
		while (*++cix != '\n');
	}
endGame:
#ifdef PROFILE
	printf("->Time for matrix generation: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	printf("Total reads expanded: %lu\n", totalReadsCnt);
	// Legacy table format output handler
	if (legacy) {
		// Prints sample and taxonomy header
		fprintf(ofp, "#OTU ID");
		char **SampP = SampDBdump - 1;
		unsigned long *OtuP = OtuList - 1;
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
				fprintf(ofp, "\n%lu", *OtuP);
				rowP = row; cnt = numSamps; 
				do fprintf(ofp, "\t%lu", *rowP++); while (--cnt);
				if (doTaxmap) 
					fprintf(ofp,"\t%s", 
						*(OtuMap_taxa + uWBS(OtuMap_otus, *(OtuList + i), blines)));
				}
			}
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
			fprintf(ofp, "{\n\"id\":null,\n\"format\": \"NINJA-BIOM " NINJA_VER " (BIOM 1.0)\",\n"
			"\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\",\n"
			"\"type\": \"OTU table\",\n\"generated_by\": \"NINJA " NINJA_VER "\",\n"
			"\"date\": \"%d-%d-%dT%s%d:%s%d:%s%d\",\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, 
			tm.tm_hour < 10 ? "0" : "", tm.tm_hour, tm.tm_min < 10 ? "0" : "", tm.tm_min,
			tm.tm_sec < 10 ? "0" : "", tm.tm_sec);

			// Write the rows in the biom format
			fprintf(ofp, "\"rows\":[");
			unsigned long *OtuP = OtuList - 1;
			int ix = ilines; do { // display the "row" lines
				if (imap[ilines - ix] == -1) {++OtuP; continue;}
				fprintf(ofp,"\n\t{\"id\":\"%lu\", \"metadata\":",*++OtuP);
				if (doTaxmap) {
					fprintf(ofp,"{\"taxonomy\":[");
					char *taxon = *(OtuMap_taxa + uWBS(OtuMap_otus, *OtuP, blines)) - 2, *tP;
					int i; for (i = 0; i < 6; i++) { // there are 6 ;'s to consider for 7 taxa
						tP = taxon + 2;
						while (*++taxon != ';')
							; *taxon = '\0';
						fprintf(ofp, "\"%s\", ", tP);
						*taxon = ';';
					}
					fprintf(ofp,"\"%s\"]}}", taxon + 2);
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
				"\"shape\": [%lu, %lu],\n", ilines, numSamps);
			fprintf(ofp,"\"data\":[");
			
			// loop thru all points
			for (i = 0; i < ilines; i++) for (j = 0; j < numSamps; j++) {
				if (OtuTable[i * numSamps + j]) 
					fprintf(ofp, "[%lu,%lu,%lu],\n\t", imap[i], j, OtuTable[i * numSamps + j]);
			}
			fseek(ofp, -3, SEEK_CUR);
			fprintf(ofp,"]\n}");
		}
		
#ifdef PROFILE
	printf("->Time for filtering and table write: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	printf("Run complete.\n"); 
	return 0;
}
