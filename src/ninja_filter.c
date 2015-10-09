#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define PRINT_USAGE() \
{\
	printf( "\nNINJA Reparse: ninja_minifilter. Usage:\n");\
	printf( "ninja_filter in_reads.fa out_filtered.fa out_sampDB.db [<trim>] [RC] [D [x]]\n" );\
	printf("\nINPUT PARAMETERS:\n");\
	printf( "in_reads.fa: the reads you wish to process\n");\
	printf( "\n" "OUTPUT PARAMETERS:\n");\
	printf( "out_filtered.fa: the filtered fasta to feed to the aligner (BT2)\n");\
	printf( "out_sampDB: bookkeeping DB required by ninja_parse\n");\
	printf( "<trim> (optional, numeric): specify the number of bases to keep\n");\
	printf( "[RC] (optional): if \"RC\" is specified, reverse-complement seqs\n");\
	printf( "[D] <x.y> (optional): Denoise [duplicates x, kmer rarity y]\n");\
	printf( "(Note: .y k-mer filtering is NOT YET IMPLEMENTED in minifilter.\n");\
	return 2;\
}

//#include "fgets2.c"
#define LINELEN UINT16_MAX
#ifndef PACKSIZE
	#define PACKSIZE 32
#endif
#if PACKSIZE==64
	#define WTYPE __uint128_t
	#define SEQPACKS 1024
	#define RSHFT 126
#elif PACKSIZE==32
	#define WTYPE uint64_t
	#define SEQPACKS 2048
	#define RSHFT 62
#elif PACKSIZE==16
	#define WTYPE uint32_t
	#define SEQPACKS 4096
	#define RSHFT 30
#elif PACKSIZE==8
	#define WTYPE uint16_t
	#define SEQPACKS 8192
	#define RSHFT 14
#elif PACKSIZE==4
	#define WTYPE uint8_t
	#define SEQPACKS 16384
	#define RSHFT 6
#endif
//#define SEQPACKS LINELEN/PACKSIZE
//#define RSHFT (PACKSIZE*2)-2

typedef struct __attribute__ ((__packed__)) {
	WTYPE word;
	uint32_t ix;
	uint16_t length;
} SortBlock2;

#ifdef USE_QSORT
#include "qsort.h"
void SB2_qsort(SortBlock2 *arr, unsigned n) {
	#define SB2_LT(a,b) ((a->word < b->word) || \
		(a->word == b->word && a->length < b->length))
	QSORT(SortBlock2, arr, n, SB2_LT);
}
#endif

WTYPE *C2Xb;
char *X2C = "ACGTNNNNNNNNNNNNNNNN";
char *X2C_RC = "TGCANNNNNNNNNNNNNNNN";

inline void num2word(WTYPE num, char * word) {
	int go = 0; for (; go < PACKSIZE; go++) {
		WTYPE temp = (WTYPE)num >> RSHFT;
		word[go] = X2C[temp];
		num <<= 2;
	}
}

inline void num2wordRC(WTYPE num, char * word) {
	int go = PACKSIZE-1; for (; go > -1; go--) {
		WTYPE temp = (WTYPE)num >> RSHFT;
		word[go] = X2C_RC[temp];
		num <<= 2;
	}
}

inline char * decodeStringX(WTYPE * Seq, uint16_t length, char *word, char *newString) {
	unsigned clumps = length/PACKSIZE;
	if (PACKSIZE*clumps < length) ++clumps;
	int z = 0; for (; z < clumps-1; z++) 
		num2word(Seq[z],newString + z*PACKSIZE);
	num2word(Seq[clumps-1],newString+z*PACKSIZE);
	newString[length] = 0;
	return newString;
}

inline char * decodeStringXRC(WTYPE * Seq, uint16_t length, char *word, char *newString) {
	newString[length] = 0;
	unsigned clumps = length/PACKSIZE;
	if (PACKSIZE*clumps < length) ++clumps;
	int z = clumps-2; for (; z > -1; z--) 
		num2wordRC(Seq[z],newString + length - (z+1) *PACKSIZE);
	num2wordRC(Seq[clumps-1],word);
	register int fold = length % PACKSIZE; if (!fold) fold = PACKSIZE;
	memcpy(newString,word+PACKSIZE-fold, fold); 
	return newString;
}

// strict comparator
int xcmp(WTYPE *Seq1, WTYPE *Seq2, uint16_t len1, uint16_t len2) {
	unsigned length = len1 < len2 ? len1 : len2; //len1 is min
	register unsigned clumps = (unsigned)length/PACKSIZE;
	if (PACKSIZE*clumps < length) ++clumps;
	int z = 0; for (; z < clumps; ++z) if (Seq1[z]!=Seq2[z])
		return Seq1[z] < Seq2[z] ? -1 : 1;
	return len1 < len2 ? -1 : len1 > len2;
}

// pre-sorted compactor
int ycmp(WTYPE *Seq1, WTYPE *Seq2, uint16_t len1, uint16_t len2) {
	if (len1 > len2) return 1; // lexicographic guarantee
	int clumps = (unsigned)len1/PACKSIZE;
	if (PACKSIZE*clumps < len1) ++clumps;
	int z = 0; for (; z < clumps-1; ++z) if (Seq1[z]!=Seq2[z])
		return 1; 
	// Can differ by length in last clump
	if (Seq1[z] == Seq2[z]) return 0;
	if (Seq1[z] > Seq2[z]) return 1; // seq2 must be superset
	unsigned shift = len1 % PACKSIZE;
	if (shift) shift = (PACKSIZE - shift) * 2;
	return (Seq1[z] >> shift) != (Seq2[z] >> shift);
}

// pre-sorted filter
int zcmp(WTYPE *Seq1, WTYPE *Seq2, uint16_t len1, uint16_t len2) {
	if (len1 != len2) return 1;
	register unsigned clumps = (unsigned)len1/PACKSIZE;
	if (PACKSIZE*clumps < len1) ++clumps;
	int z = 0; for (; z < clumps; ++z) if (Seq1[z]!=Seq2[z])
		return 1;
	return 0;
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

inline size_t crBST(char *key, size_t sz, char **String) {
	char **p = String;
	while (sz) {
		size_t w = sz >> 1; 
		char *ref_s = *(p+w+1), *key_s = key;
		
		while (*ref_s == *key_s++) if (!*ref_s++) return p+w+1-String; 
		if (*ref_s < *(key_s-1)) { p+=w+1; sz-=w+1; }
		else sz = w;
	}
	return p - String;
}

int SB2Cmp(blk1, blk2) register const void *blk1, *blk2; {
	if (((SortBlock2 *)blk1)->word < ((SortBlock2 *)blk2)->word) return -1;
	if (((SortBlock2 *)blk1)->word > ((SortBlock2 *)blk2)->word) return 1;
	if (((SortBlock2 *)blk1)->length == ((SortBlock2 *)blk2)->length) return 0;
	if (((SortBlock2 *)blk1)->length < ((SortBlock2 *)blk2)->length) return -1;
	return 1;
}

void superSort2(uint32_t *SeqIX, WTYPE **base, uint16_t *Lengths, 
int depth, size_t beginRange, size_t endRange) {
	size_t n = endRange - beginRange; // endRange is one after last index
	SortBlock2 *BinPtrs = malloc(n * sizeof(SortBlock2)); 
	if (!BinPtrs) {puts("Error-MemoryBinPtrs"); return;}
	size_t depthSize = (depth+1) * PACKSIZE;
	size_t i = beginRange; for (; i < endRange; ++i) 
		BinPtrs[i-beginRange] = (SortBlock2){base[SeqIX[i]][depth],SeqIX[i],
			Lengths[SeqIX[i]] <= depthSize ? Lengths[SeqIX[i]] : 0};
	#ifdef USE_QSORT
		SB2_qsort(BinPtrs,n);
	#else
		qsort(BinPtrs, n, sizeof(*BinPtrs), SB2Cmp);
	#endif
	for (i=beginRange; i < endRange; ++i) 
		SeqIX[i] = BinPtrs[i-beginRange].ix; 
	free(BinPtrs);
	
	#define CASCADE_MERGE() \
	if (i != lastUniq + 1) { \
		/* Merge swapping indices for truncated pairs */ \
		size_t z = lastUniq; for (; z < i; ++z) { \
			if (Lengths[SeqIX[z]] <= depthSize) { \
				 if (z > lastUniq) { \
					/* swap this ix with the ix at lastUniq++ */ \
					uint32_t temp = SeqIX[z]; \
					SeqIX[z] = SeqIX[lastUniq]; \
					SeqIX[lastUniq] = temp; \
				}  \
				++lastUniq; \
			} \
		} \
		/* Spawn a new sort on the remainder */ \
		if (lastUniq < i-1) \
			superSort2(SeqIX, base, Lengths, depth+1, lastUniq, i); \
	}
	
	// Check for duplicates; for each set, move truncations to top
	WTYPE curElem = base[SeqIX[beginRange]][depth]; 
	size_t lastUniq = beginRange;
	for (i=beginRange + 1; i < endRange; ++i) {
		if (base[SeqIX[i]][depth] != curElem) {
			CASCADE_MERGE();
			curElem = base[SeqIX[i]][depth];
			lastUniq = i;
		}
	}
	CASCADE_MERGE(); // end cap
}
	
int main( int argc, char *argv[] ) {
	clock_t start; double cpu_time_used; start = clock(); // profiler
	// Debugging statements
	printf("type size=%u, shift=%u, pack=%u\n", sizeof(WTYPE), RSHFT, PACKSIZE);
	printf("max int size=%u/%u\n",sizeof(unsigned),sizeof(uint64_t));
	printf("Size of SortBlock2=%u\n",sizeof(SortBlock2));
	//return 0;
	if ( argc < 4 || argc > 8 ) PRINT_USAGE();
	char *inputFilename = argv[1], *outputFasta = argv[2], *outputDB = argv[3];
	FILE *fp = fopen(inputFilename, "rb");
	if (fp == NULL) { puts("Invalid input"); return 2; }
	FILE *off = fopen(outputFasta, "wb"), *ofd = fopen(outputDB,"wb");
	if (!off || !ofd) { puts("Invalid output file(s)"); return 2; }
	size_t trim = UINT16_MAX;
	int doRC = 0; 
	double filt_i = 0.f; unsigned copyNumThres = 0; // denoisers (k-filter not yet implemented)
	// Denoises at default intensity
	if (argc > 4 && !strcmp(argv[argc-1],"D")) { 
		//filt_i = .0001f; // k-filtering not yet implemented!
		//printf("Performing NINJA statistical denoising of DEFAULT intensity %f (%.3f%%)\n", filt_i,filt_i*100.f);
		puts("k-denoise not implemented yet.");
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
	//(WTYPE *Seq1, WTYPE *Seq2, uint16_t len1, uint16_t len2)
	int (*cmpF)(WTYPE *, WTYPE *, uint16_t, uint16_t) = 
		copyNumThres ? &ycmp : &zcmp;
	if (!copyNumThres) copyNumThres = 1 ; // filt_i ? -1 : 1; // uncomment when .y implemented
	if (argc > 4 && !strcmp(argv[argc-1],"RC")) {
		printf("Reverse complementing the sequences.\n");
		doRC = 1; --argc;
	}
	// Flags for truncation after specified base
	if (argc == 5) { 
		trim = atoi(argv[argc-1]);
		printf("Trimming input sequences to %d bases.\n", trim);
	}
	
	C2Xb = calloc(128,sizeof(WTYPE));
	C2Xb['a'] = 0; C2Xb['A'] = 0; 
	C2Xb['c'] = 1; C2Xb['C'] = 1; 
	C2Xb['g'] = 2; C2Xb['G'] = 2;
	C2Xb['t'] = 3; C2Xb['T'] = 3;
	//ctx_t *ctx = init_fgets_sse2 (LINELEN*32); // delete
	//next_t *ne; // delete
	
	size_t numElem = 1000, ns=0;
	char **Samples = malloc(numElem*sizeof(char *));
	WTYPE **ReadsX = malloc(numElem*sizeof(WTYPE *));
	uint16_t *Sizes = calloc(numElem,sizeof(uint16_t));
	char *line = malloc(LINELEN + 1); // read up to  65k
	
	// new (testing)
	/* char *buffer = malloc(100*LINELEN+1); // read up to 6.5mb
	buffer[100*LINELEN] = 0;
	size_t numRead = 0; */
	while (line = fgets(line,LINELEN,fp)) { // old
	//while (numRead = fread(buffer,1,100*LINELEN,fp)) { // new (testing)
		if (ns == numElem) {
			numElem *= 2;
			Samples = realloc(Samples,numElem * sizeof(char *));
			ReadsX = realloc(ReadsX, numElem * sizeof(WTYPE *));
			Sizes = realloc(Sizes, numElem*sizeof(uint16_t));
			if (!Samples || !ReadsX || !Sizes) {puts("Error in resize"); return 0;}
			memset(Sizes+numElem/2 + 1,0,(numElem/2-1)*sizeof(uint16_t));
		}
		// copy in the sample name up to _ or null minus 1
		char *src = line + 1;
		
		while (*src != '_' && *src != ' ' && *src != '\n') ++src; 
		Samples[ns] = malloc(src - line);
		if (!Samples[ns]) {puts("Not enough Samples[ns] mem"); return 0;}
		
		char *dest = Samples[ns]; //char *src = line + 1;
		char *beginSample = line + 1; while (beginSample < src) 
			*dest++ = *beginSample++;
		*dest = 0;

		// copy in the encoded sequence
		if (!(line = fgets(line,LINELEN,fp))) 
			{ puts("Error reading file."); return 0; }
		src = line;
		
		register size_t length = strlen(src);
		if (src[length-1] == '\n') --length; // lop off newline(s)
		if (src[length-1] == '\r') --length; 
		if (trim < length) length = trim;
		
		size_t numPacks = length/PACKSIZE;
		if (numPacks * PACKSIZE < length) ++numPacks;
		
		Sizes[ns] = length; 
		ReadsX[ns] = malloc(numPacks*sizeof(WTYPE));
		if (!ReadsX[ns]) {puts("Bad ReadsX[ns] mem"); return 1; }
		
		WTYPE *thisPack = ReadsX[ns];
		WTYPE clump; int k;
		while (length--) {
			k = 1; clump = C2Xb[*src++];
			while (k < PACKSIZE && length) {
				clump <<= 2u;
				clump += C2Xb[*src++];
				++k; --length;
			}
			if (k != PACKSIZE) clump <<= (2* (PACKSIZE-k));
			*thisPack++ = clump;
		}
		++ns;
	}
	fclose(fp);
	free(line);
	
	// Shrink data structures for more memory
	Samples = realloc(Samples,ns * sizeof(char *));
	ReadsX = realloc(ReadsX, ns * sizeof(WTYPE *));
	Sizes = realloc(Sizes, ns * sizeof(uint16_t));

	//printf("Num of sequences: %u\n",ns);
	if (ns > UINT32_MAX) {puts("Too many sequences (>4 bil)."); return 4;}
	printf("Total reads considered: %u\n",ns);
#ifdef PROFILE
	printf("->Short read parse: %f\n", 
		((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	
	// Create index structure for sequences read (in 32-bit)
	uint32_t *SeqIX = malloc(sizeof(uint32_t) * ns);
	size_t k = 0; 
	for (; k < ns; ++k) SeqIX[k] = k;
	superSort2(SeqIX, ReadsX, Sizes, 0,0,ns);
	printf("\nDONE SORTING. \n");
	
	char ***smpSrt = malloc(ns * sizeof(char **)),
		**SmpDD = malloc(ns * sizeof(char *));
	if (!smpSrt || !SmpDD) 
		{ printf("Out of post-memory: parray.\n"); return 3; }
	for (k=0; k < ns; ++k) smpSrt[k] = &Samples[k];
	twrqs(smpSrt, ns, 0);
	*SmpDD = **smpSrt; // store first sample
	unsigned x = 1; for (k=1; k < ns; ++k) 
		if (strcmp(*smpSrt[k-1],*smpSrt[k])) SmpDD[x++] = *smpSrt[k];
	free(smpSrt);
	SmpDD = realloc(SmpDD,sizeof(char*)*x);
	printf("%d Samples found.\n",x);
	fprintf(ofd, "%u\n", x);
	for (k=0; k < x; ++k) fprintf(ofd,"%s\n",SmpDD[k]);
	
#ifdef PROFILE
	printf("->Short read sample prep: %f\n", 
		((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// Create counts array of integers parallel to the unique samples array
	unsigned *Counts = calloc(x, sizeof(unsigned));
	if (!Counts) {puts("unable to allocate counts"); return 3;}
	#define WRITE_SUPPORTED_DUPE() {\
		six = 0; \
		if (copies >= copyNumThres) { \
			int y = 0; for (; y < x; ++y) \
				if (Counts[y]) fprintf(ofd,"%u:%u:",y,Counts[y]), Counts[y] = 0; \
			fprintf(ofd,"\n"); \
			if (doRC) fprintf(off,">%u\n%s\n",rix++, decodeStringXRC(ReadsX[prevIX], \
				Sizes[prevIX],word,string)); \
			else fprintf(off,">%u\n%s\n", rix++, decodeStringX(ReadsX[prevIX], \
				Sizes[prevIX],word,string)); \
		} \
		else memset(Counts,0,x*sizeof(unsigned)); \
		copies = 1; \
	}
	
	unsigned copies = 1, dupes = 0, six, rix=0;
	++Counts[crBST(Samples[*SeqIX],x-1,SmpDD)]; // add first count (counts as a copy)
	char *string = malloc(UINT16_MAX), *word = calloc(PACKSIZE+1,1);
	unsigned prevIX, thisIX;
	for (k=1; k < ns; ++k) {
		prevIX = SeqIX[k-1]; thisIX = SeqIX[k];
		++Counts[crBST(Samples[thisIX],x-1,SmpDD)];
		//++*(Counts + crBST(*(Samples+*(SeqIX+k)),x-1,SmpDD));
		if (cmpF(ReadsX[prevIX],ReadsX[thisIX],Sizes[prevIX], Sizes[thisIX])) 
			WRITE_SUPPORTED_DUPE()
		else ++copies, ++dupes;
	}
	prevIX = thisIX;
	WRITE_SUPPORTED_DUPE();
	
#ifdef PROFILE
	printf("->Mapping and file writing: %f\n", 
		((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// todo: free more, process more
	free (SeqIX);
	free (string);
}