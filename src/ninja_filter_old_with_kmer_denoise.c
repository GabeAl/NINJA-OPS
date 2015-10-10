/* NINJA OTU Picker: NINJA Is Not Just Another OTU Picker -- filter program
   Knights Lab (www.knightslab.org/ninja)
   This program generates the databases necessary for OTU mapping with bowtie2.
   Commented and revised by Henry Ward (henry.n.ward@lawrence.edu)
   
   Compilation information (GCC):
   Ascribes to std=gnu and doesn't require C99 support or UNIX intrinsics.
   Flags: -m64 -O3 -flto ninja_filter.c -o ninja_filter
   Optional Definitions (use -D DEF[=SETTING] in gcc at compile time to use): 
    WORDTYPE=(C intrinsic type), KMER=(length of kmer), DEBUG (print debug info), 
    PROFILE (print speed info), LOGK (logs kmers)
*/
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#define PRINT_USAGE() \
{\
	printf( "\nNINJA Is Not Just Another OTU Picker: filter program. Usage:\n");\
	printf( "ninja_filter in_reads.fa out_filtered.fa out_sampDB.db [<trim>] [RC] [D [x]]\n" );\
	printf("\nINPUT PARAMETERS:\n");\
	printf( "in_reads.fa: the reads you wish to process\n");\
	printf( "\n" "OUTPUT PARAMETERS:\n");\
	printf( "out_filtered.fa: the filtered fasta to feed to the aligner (BT2)\n");\
	printf( "out_sampDB: bookkeeping DB required by ninja_parse\n");\
	printf( "<trim> (optional, numeric): specify the number of bases to keep\n");\
	printf( "[RC] (optional): if \"RC\" is specified, reverse-complement seqs\n");\
	printf( "[D] <x.y> (optional): Denoise [duplicates x, kmer rarity y]\n");\
	return 1;\
}

/**
 *  Utility macros
 */
#ifndef min 
#define min(a, b) ((a)<=(b) ? (a) : (b)) 
#endif
#define ch(i) *(**(a+i) + depth) 
#define med3(ia, ib, ic) med3func(a, ia, ib, ic, depth)
#define CUTOFF 10
#define MEDCUT 50
// Denoising complement utility macros
#define C2Xb(x) x=='A' ? A_ID : x=='C' ? C_ID : x=='G' ? G_ID : x=='T' ? T_ID : 0
#define X2Cb(x)  x==A_ID ? 'A' : x==C_ID ? 'C' : x==G_ID ? 'G' : x==T_ID ? 'T' : 'N'
#define C2X(x) x=='A' ? 0 : x=='C' ? 1 : x=='G' ? 2 : x=='T' ? 3 : 0
#define X2C(x) x==0 ? 'A' : x==1 ? 'C' : x==2 ? 'G' : x==3 ? 'T' : 'N' 

// Returns complement of a base pair
#define CMPT(x) \
	x=='A' ? 'T' : x=='C' ? 'G' : x=='G' ? 'C' : x=='T' ? 'A' :  \
	x=='a' ? 'T' : x=='c' ? 'G' : x=='g' ? 'C' : x=='t' ? 'A' : x 
// Returns transcribed complement of a base pair
#define CMPTR(x) \
	x=='A' ? 'U' : x=='C' ? 'G' : x=='G' ? 'C' : x=='U' ? 'A' :  \
	x=='a' ? 'U' : x=='c' ? 'G' : x=='g' ? 'C' : x=='u' ? 'A' : x 

// Set up the word size shifting for denoising
#ifndef WORDTYPE
#define WORDTYPE uint32_t
#endif
#ifndef KMER
#define KMER 8
#endif
unsigned WORDSIZE=KMER;
unsigned WSHFT = (((sizeof(WORDTYPE)*4) << 1) - 2), 
	RSHFT = sizeof(WORDTYPE) * 8 - KMER*2; // Later auto-set as well to ((WORDSIZE << 1) - 2)
size_t MAXLEN; 
WORDTYPE A_ID, C_ID, G_ID, T_ID; // Will be defined later in PREP_WORDS
WORDTYPE *C2XbL;

// Prepares word lookup matrix, kmer bounds, bitshift profiles
#define PREP_WORDS() \
{ \
	if (KMER > sizeof(WORDTYPE) * 4) {printf("invalid K for wordsize\n"); return 1; } \
	MAXLEN = (size_t)1 << (WORDSIZE * 2); \
	WSHFT = (((sizeof(WORDTYPE)*4) << 1) - 2); \
	RSHFT = sizeof(WORDTYPE) * 8 - KMER*2; \
	A_ID = (WORDTYPE)0; C_ID = (WORDTYPE)1 << WSHFT; \
	G_ID = (WORDTYPE)2 << WSHFT; T_ID = (WORDTYPE)3 << WSHFT; \
	C2XbL = calloc(128,sizeof(WORDTYPE));\
	C2XbL['a'] = A_ID; C2XbL['A'] = A_ID; \
	C2XbL['c'] = C_ID; C2XbL['C'] = C_ID; \
	C2XbL['g'] = G_ID; C2XbL['G'] = G_ID; \
	C2XbL['t'] = T_ID; C2XbL['T'] = T_ID; \
}

/** 
 *  Utility functions
 */
// Swaps two characters in a vector
inline void swap(char ***a, int i, int j) 
	{ char **t = *(a+i); *(a+i) = *(a+j); *(a+j) = t; }
inline void vecswap(char ***a, int i, int j, int n) 
	{ while (n-->0) swap(a, i++, j++); }

// Specialized qcomp comparator functions for pointer indexing
int xcmp(register const char *str1, register const char *str2) {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return (*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}
int ycmp(register const char *str1, register const char *str2) {
	while (*str1 == *str2++) if (!*str1++) return 0; 
	return *str1 && *(str2 - 1); //(*(const unsigned char *)str1 - *(const unsigned char *)(str2 - 1));
}
inline int ucmp(int1, int2) register const unsigned int *int1, *int2; {
	return **(unsigned **)int1 - **(unsigned **)int2; 
}
inline int comp(int1, int2) register const void *int1, *int2; {
	return *(unsigned *)int1 - *(unsigned *)int2; 
}
int cmp64u(int1, int2) register const void *int1, *int2; {
	return *(uint64_t *)int1 < *(uint64_t *)int2 ? -1 : 
		*(uint64_t *)int1 > *(uint64_t *)int2; 
}
// Specialized inline character-based binary string search
inline size_t crBST(char *key, size_t sz, char **String) {
	char **p = String; //, *ref_s, *key_s; 
	while (sz) {
		size_t w = sz >> 1; //ref_p = p + w + 1;
		char *ref_s = *(p+w+1), *key_s = key;
		
		while (*ref_s == *key_s++) if (!*ref_s++) return p+w+1-String; 
		if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	return p - String;
}

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
		if (xcmp(**(a+j-1)+depth, **(a+j)+depth) <= 0) break;
		swap(a, j, j-1);
	} 
}  

// Counting sort, for bins of k-mers in denoising
typedef uint64_t type;   
type * countSrt(type *array, size_t length) {
	type *ptr = array, max = *array;
	size_t ctr = length;
	while (--ctr) if (*++ptr > max) max = *ptr;
	
	type *Countarr = calloc(max,sizeof(type)), *caP = Countarr - 1;
	ptr = array - 1; ctr = length;
	do ++Countarr[*++ptr]; while (--ctr);
	type *new = malloc(length * sizeof(type));
	ptr = new - 1; ctr = max + 1;
	//if (!ctr) return 0;
	do {
		type reps = *++caP;
		while (reps--) *++ptr = caP - Countarr;
	} while (--ctr);
	free(Countarr);
	return new;
}

// 3-way Radix Quicksort (Dobbs-optimized)
void twrqs(char ***a, unsigned n, int depth) {
	if (n < CUTOFF) { inssort(a, n, depth); return; }
	unsigned pl = 0, pm = n >> 1, d;
	int le, lt, gt, ge, r, v, pn = n-1;
	// If large enough, get median of median
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
	// Core QS module. Partitions the data recursively
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

// Print large numbers of arbitrary type (WORDTYPE)
int print_bignum(WORDTYPE n) {
  if (n == 0)  return printf("0\n");

  char str[40] = {0}; // log10(1 << 128) + '\0'
  char *s = str + sizeof(str) - 1; // Start at the end
  while (n != 0) {
    if (s == str) return -1; // Never happens

    *--s = "0123456789"[n % 10]; // Save last digit
    n /= 10;                     // Drop last digit
  }
  return printf("%s", s);
}

// Generates, permutes, and displays a test word from type WORDTYPE
void testWord() {
	char *testSeq = "ACGTTGACAACCCT\0", *testSeqP = testSeq;
 	WORDTYPE hoho= 0; 
	do {
		printf("Loop (%d): %u. ", (unsigned)(testSeqP - testSeq), hoho);
		hoho += C2Xb(*testSeqP);
		printf("After '%c' (%u): %u.\n", *testSeqP, C2X(*testSeqP), hoho);
		// Replay what's inside hoho, letter by letter.
		size_t x=WORDSIZE; WORDTYPE hohoTemp = hoho; char* hi = calloc(WORDSIZE+1,1); //"XXXXXXXX\0"; 
		memset(hi,'X',WORDSIZE);
		while (x--) {
			WORDTYPE temp = hohoTemp >> WSHFT;
			printf("Temp: %u ('%c')\n", temp, X2C(temp));
			hi[x] = X2C(temp);
			hohoTemp <<= 2;
		}  
		printf("hi: %s\n", hi);
		hoho >>= 2;
		printf("After right shift: ");
		print_bignum(hoho); 
		printf("\n");
		
	} while (*++testSeqP);
}

/** 
 * Main filter algorithm. Parameters are as follows 
 * INPUT	in_reads.fa:			sequence reads to be processed
 * OUTPUT	out_filtered.fa:		filtered fasta to feed to the aligner (Bowtie2)
 *			out_sampDB:				bookkeeping DB required by ninja_parse as an input parameter
 *			<trim> (optional):		number of base pairs to keep in each sequence, trims ends of sequences
 *			[RC] (optional):		reverse-complement sequences
 *			[D <x.y>] (optional):	denoise at x duplicates and kmer rarity y
 * 
 */
int main( int argc, char *argv[] )
{
	// Requires 3 parameters for execution. Else, displays below info
	if ( argc < 4 || argc > 8 ) PRINT_USAGE();
	PREP_WORDS(); // Prepares for kmer denoising by creating WORDTYPE objects
	// Flags for optional arguments
	int doRC = 0, trim = 0; 
	double filt_i = 0.f;
	unsigned copyNumThres = 0; 
	// Denoises at default intensity
	if (argc > 4 && !strcmp(argv[argc-1],"D")) { 
		filt_i = .0001f;
		printf("Performing NINJA statistical denoising of DEFAULT intensity %f (%.3f%%)\n", filt_i,filt_i*100.f);
		--argc;
	}
	// Denoises at specified intensity in the form x.y
	else if (argc > 4 && !strcmp(argv[argc-2],"D")) {
		filt_i = atof(argv[argc-1]);
		if (filt_i < 0 || filt_i > 100) printf("Invalid denoising intensity (expect #REPS.%%%%).\n"); 
		else {
			if (filt_i >= 1.f) {
				copyNumThres = filt_i;
				filt_i -= copyNumThres;
				printf("Performing NINJA replicon-denoising of constructed reps: %u\n", copyNumThres);
			}
			if (filt_i) { // Use the decimal remainder as kmer denoising
				printf("Performing NINJA statistical denoising of intensity %f (%.3f%%)\n", filt_i,filt_i*100.f);
				if (copyNumThres) ++copyNumThres;
			}
		}
		argc -= 2;
	}
	// Flags for reverse complementing
	if (argc > 4 && !strcmp(argv[argc-1],"RC")) {
		printf("Reverse complementing the sequences.\n");
		doRC = 1; --argc;
	}
	// Flags for truncation after specified base
	if (argc == 5) { 
		trim = atoi(argv[argc-1]);
		printf("Trimming input sequences to %d bases.\n", trim);
	}
	char *inmode = "rb", *outmode = "wb"; 
	clock_t start; double cpu_time_used; start = clock();
	
	// Takes in a file with all the user's short reads for later sorting/filtering
	// Opens all three required parameters as files. Displays specific error upon failure.
	FILE *seqs = fopen(argv[1], inmode);
	if (!seqs) { printf("Error: cannot open input reads fasta.\n"); return 1; }
	FILE *filtered = fopen(argv[2], outmode);
	if (!filtered) {printf("could not open filtered reads output file.\n"); return 1; }
	FILE *sampDBfile = fopen(argv[3], outmode);
	if (!sampDBfile) {printf("could not write to sample DB.\n"); return 1; }
	
#ifdef DEBUG
	if (filt_i) 
		printf("WORDSIZE: %I64u, KMER: %u, MAXLEN: %I64u, WSHFT: %u, RSHFT: %u\n",WORDSIZE,KMER,MAXLEN,WSHFT,RSHFT); 
#endif
	
	// Finds length of input sequence FASTA file, in bytes
	fseek(seqs, 0, SEEK_END); size_t fsize = ftell(seqs); fseek(seqs, 0, SEEK_SET); 

	// Reads seq file as a string, stores in "string." Displays error if file too large
	char *string = malloc(fsize + 1), *stringP = string - 1; 
	if (!string) { printf("Insufficient memory for short read read-in.\n"); return 1; }
	fread(string, fsize, 1, seqs); // Reads into string: elems of fsize bytes, 1 elem, using ifp pointer
	fclose(seqs); // Closes the file
	
	// Max read length: 100KB. Max sample length: 1KB. 
	// Outputs specific error if limits exceeded. 
	size_t numEntries = 0, curLen = 1000;  // Floating tally
	int RLEN = 100000, SLEN = 1024; // Common sense max bounds on read and sample lengths
	char *seqBuf = malloc(RLEN), *seqBufP = seqBuf; // Assumption: no read longer than 100KB
	char **seqArr = malloc(curLen * sizeof(char *)), **seqArrP = seqArr,
	*smpBuf = malloc(SLEN), *smpBufP = smpBuf, // Assumption: no sample name > 1KB
	**smpArr = malloc(curLen * sizeof(char *)), **smpArrP = smpArr,
	*thisSeq, *thisSmp;
	if (!smpBuf || !smpArr || !seqArr) 
		{ printf("Cannot allocate post-run memory.\n"); return 1; }
	
	// Main parsing loop. Also performs kmer generation and read trimming
	size_t smpCharNum = 0, bufCharNum = 0, inSeq = 0, seqChars = 0, k = 0;
	// Sets up a kmer dictionary for all possible kmers
	uint64_t *words = calloc(filt_i ? MAXLEN : 1, sizeof(uint64_t)), *wordPtr = words;
	size_t sum = 0;
	if (!words) {printf("out of word allocation memory.\n"); return 1;}
	WORDTYPE word = 0;  // Uses kmer as index into "words" dictionary
	
	stringP = string - 1; 
	while (*++stringP) { 
		// Skips '>'
		if (*stringP == '>' || *stringP == '-') { continue; } 
		else if (*stringP == '\n') { 
			// Prepares to parse sequence
			if (!inSeq) { 
				if (++numEntries == curLen) { // Resizes all floating arrays
					size_t offset = seqArrP - seqArr; 
					seqArr = realloc(seqArr, (curLen *= 2) *sizeof(char *));
					smpArr = realloc(smpArr, curLen *sizeof(char *));
					if (!seqArr || !smpArr) 
						{ printf("out of resizing memory.\n"); return 1; }
					seqArrP = seqArr + offset;
					smpArrP = smpArr + offset;
				}
				*smpArrP = malloc(smpCharNum + 1);
#ifdef DEBUG
				if (!*smpArrP) { printf("sample memory depleted!\n"); return 1; }
#endif
				thisSmp = *smpArrP++;
				smpBufP = smpBuf;
				do *thisSmp++ = *smpBufP++; while (--smpCharNum);
				memset(thisSmp, '\0', 1);
				seqBufP = seqBuf;
				bufCharNum = 0;
			}
			// Prepares to parse header
			else { 
				*seqArrP = malloc(bufCharNum + 1);
#ifdef DEBUG
				if (!*seqArrP) { printf("sequence memory depleted!\n"); return 1; }
#endif
				thisSeq = *seqArrP++;
				seqBufP = seqBuf;
				do *thisSeq++ = *seqBufP++; while (--bufCharNum);
				memset(thisSeq, '\0', 1);
				smpBufP = smpBuf;
				smpCharNum = 0;
				if (filt_i) {
					word = 0; k = 0;
				}
			}
			inSeq ^= 1;  // Toggles between header and sequence for storing
		}
		// Parses the sequence
		else if (inSeq) { 
			*seqBufP++ = *stringP;
			++bufCharNum;
			// Generates kmer at current letter
			if (filt_i) {
				word += C2XbL[*stringP];
				if (++k >= WORDSIZE) { 
					++*(words + (word >> RSHFT)); // Increments word count at ix=word
					++sum; // Increments total sum of words considered
				}
				word >>= 2u; // Left-shift bits in current word for kmer
			}
#ifdef DEBUG
			if (bufCharNum > RLEN) {printf("sequence buffer overflow\n"); return 1; }
#endif
			// Trims sequences if user specified
			if (trim && bufCharNum >= trim) { while (*++stringP != '\n'); --stringP; }
		}
		else { // Finds sample tag in current header
			if (*stringP == '_') { // Fast-forwards to end of line after '_'
				while (*++stringP != '\n'); 
				--stringP; 
				continue;
			} 
			*smpBufP++ = *stringP;
			++smpCharNum;
#ifdef DEBUG
			if (smpCharNum > SLEN) {printf("sample name buffer overflow\n"); return 1; }
#endif
		}
	}
	// Reallocate the active vs paged memory on sequence and sample arrays
	seqArr = realloc(seqArr, numEntries * sizeof(char *));
	smpArr = realloc(smpArr, numEntries * sizeof(char *));
	
#ifdef PROFILE
	printf("->Short read parse: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	printf("Total short reads: %lu\n", numEntries);
	// Prepares for denoising
	uint64_t *sortWrd;
	size_t len = MAXLEN, rareCutoff, rareIX;
	double rareThres = filt_i ? filt_i : -1.f;
	unsigned nThres;
	
	// Considers kmer component of denoising
	if (filt_i) {
		wordPtr = words;
		sortWrd = malloc(MAXLEN *sizeof(uint64_t));
		uint64_t *sortWrdPtr = sortWrd;
		size_t lenC = MAXLEN; 
		do *sortWrdPtr++ = *wordPtr++; while (--lenC);
		qsort(sortWrd, MAXLEN, sizeof(uint64_t), cmp64u);
		//sortWrd = countSrt(words,MAXLEN);
		
#ifdef PROFILE
		printf("->SortIX: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
		wordPtr = sortWrd-1; // Resets pointer head to sorted array
		rareCutoff = rareThres * sum;
		rareIX = 0; sum = 0;
		// Traverses distribution until reaches rarity threshold
		while ((sum += *++wordPtr) < rareCutoff) ++rareIX; 
		if (!rareIX) printf("Warning: k=%u may result in spurious denoising.\n",KMER);
		nThres = sortWrd[rareIX];
#ifdef DEBUG
		printf("rareCutoff=%u, rareThres=%f, rareDupes=%u, %%Coverage=%f\n",rareCutoff, rareThres, nThres,(double)(rareIX+1)/MAXLEN);
#endif
		if (nThres == 1) {
			printf("Chosen K-mer denoising level has no effect.\n");
			filt_i = 0;
			--copyNumThres;
		}
		
	}
	// Sorts the sequences and runs deduplification
	// 1) Allocates memory and checks if memory available
	char ***parray = malloc(numEntries * sizeof(char **)), ***pp = parray,
	***smparray = malloc(numEntries * sizeof(char **)), ***smpp = smparray,
	**SmpDD = malloc(numEntries * sizeof(char *)), **ddp = SmpDD;
	if (!parray || !smparray || !SmpDD) 
		{ printf("Out of post-memory: parray.\n"); return 1; }
		
	// 2) Stores addresses of strings in auxiliary array for pointer sorting
	seqArrP = seqArr; smpArrP = smpArr;
	size_t nE = numEntries; do { *pp++ = seqArrP++; *smpp++ = smpArrP++; } while (--nE); 
	
	// 3) Sorts arrays on pointers generated above
	twrqs(parray, numEntries, 0);
	twrqs(smparray, numEntries, 0);
	
	// 4) De-dupes the sorted sample list and makes a parallel counts array
	unsigned copies = 1, dupes = 0; 
	nE = numEntries; smpp = smparray;
	while (--nE) {
		if (!xcmp(**smpp, **(smpp + 1))) ++dupes;
		else *ddp++ = **smpp;
		++smpp;
	}
	unsigned numUniq = numEntries - dupes, nU = numUniq;
	if (ddp - SmpDD <= numUniq) *ddp = **smpp; // endcap
	printf("Unique samples: %d\n", numUniq);
	// Writes into sample database file in format: sample number \n all samples
	fprintf(sampDBfile, "%u\n", numUniq);
	SmpDD = realloc(SmpDD, numUniq * sizeof(char *)); ddp = SmpDD; 
	nU = numUniq; do fprintf(sampDBfile, "%s\n", *ddp++); while (--nU);
	
#ifdef PROFILE
	printf("->Short read sample prep: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	// Creates counts array of integers parallel to the unique samples array
	unsigned *Counts = calloc(numUniq, sizeof (unsigned int)), *cpp = Counts;
	if (!Counts) {printf("unable to allocate counts\n"); return 1;}
	copies = 1; dupes = 0; 
	nE = numEntries; 
	pp = parray; 
	smpBuf = malloc(SLEN*100); smpBufP = smpBuf; // Assumption: no sample string > SLEN * 100
	if (!smpBuf) {printf("Out of memory: concatSample\n"); return 1; }
	unsigned rix = 0;
	int six;
	unsigned totCopies = 0;
	seqBufP = seqBuf + RLEN - 2; // Re-uses seqBuf and seqBufP
	memset(seqBufP,'\0',1); // termi-null

#ifdef LOGK
	FILE *log = fopen("bad_k_log.txt", outmode);
	FILE *log_dn = fopen("log_dn.txt", outmode);
#endif

	// Stores occurances of words meeting rarity threshold
	unsigned wBad = 0; 
	int (*cmpF)(register const char *, register const char *) = copyNumThres ? &ycmp : &xcmp;
	//int *cmpF = copyNumThres ? ycmp(const char *, const char *) : xcmp(const char *, const char *);
	if (!copyNumThres) copyNumThres = filt_i ? -1 : 1;
	size_t numUniq_mx = numUniq - 1;
	
	// Performs deduplication and denoising using sorted data
	//double propBad = 0;
	nE = numEntries; while (nE--) {
		++*(Counts + (six=crBST(*(smpArr + (*pp - seqArr)), numUniq_mx, SmpDD))); // specific count
		if (nE && !cmpF(**pp, **(pp + 1))) { // Dupe! 
			++copies; ++dupes; 
		} 
		else { // Writes the last sequence and its number of copies
			smpBufP = smpBuf; 
			six = 0;
			char *read = **pp;
			unsigned rLen = 0;
			double propBad = 1.f;
			if (filt_i) { // && copies <= copyNumThres) {
				
				k = 0; word = 0; wBad = 0;
#ifdef LOGK
				fprintf(log,"\n");
#endif
				unsigned lastBad = 0;
				do {
					++rLen;
					word += C2XbL[*read];
					if (++k >= WORDSIZE) {
						if (*(words + (word >> RSHFT)) < nThres) {
							if (rLen > lastBad) {
								++wBad;
								lastBad = rLen + WORDSIZE - 1;
#ifdef LOGK
								fprintf(log,"0 ");
							} else fprintf(log,"* ");
						} else fprintf(log,"1 ");
#else
							}
						}
#endif
					}
					word >>= 2u; 
				} while (*++read);
				propBad = wBad; // (double)wBad/rLen;
#ifdef LOGK
				fprintf(log,"\nCopies: %u, Bad=%u, prop=%f", copies,wBad,propBad);
#endif
			} 
			if (copies >= copyNumThres || (!propBad && 
				(copyNumThres==-1 || (copies > (int)copyNumThres - 2)))) { 
#ifdef LOGK
				fprintf(log_dn,"%u: %u\n", rix, copies);
#endif
				nU = numUniq; cpp = Counts; do {
					if (*cpp) {
						smpBufP += sprintf(smpBufP, "%u:%u:", six, *cpp);
						*cpp = 0;
					}
					++six; ++cpp;
				} while (--nU);
				fprintf(sampDBfile, "%s\n",smpBuf);
				read = **pp;
				// Takes the reverse complement of the read, if user specified
				if (doRC) { 
					do *--seqBufP = CMPT(*read); while (*++read);
					read = seqBufP;
					seqBufP = seqBuf + RLEN - 2;
				}
				fprintf(filtered, ">%u\n%s\n", rix++, read);
				totCopies += copies;
			}
			else {
				memset(Counts,0,numUniq*sizeof (unsigned int));
#ifdef LOGK
				fprintf(log_dn,"%BAD u: %u\n",rix,copies);
#endif
			}
			copies = 1;
		}
		++pp;
	}
	// Program ran successfully. Outputs summaries to user
	printf("Dupes: %d, total copies: %d (%f%%)\n", dupes, totCopies, (double)totCopies/numEntries*100);
	printf("Optimized fasta written.\n");

#ifdef PROFILE
	printf("->Read prep and Fasta write: %f\n", ((double) (clock() - start)) / CLOCKS_PER_SEC); 
	start = clock();
#endif
	return 0;
}