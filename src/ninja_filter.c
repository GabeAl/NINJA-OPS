#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define NINJA_VER "1.2"
#define PRINT_USAGE() \
{\
	printf( "\nNINJA Is Not Just Another - OTU Picking Solution v" NINJA_VER ": Filter. Usage:\n");\
	printf( "ninja_filter in_reads.fna out_PREFIX [<trim>] [RC] [D [x[.y]]] [LOG]\n" );\
	printf("\nINPUT PARAMETERS:\n");\
	printf( "in_reads.fa: the reads you wish to process\n");\
	printf( "\n" "OUTPUT PARAMETERS:\n");\
	printf( "out_PREFIX: prefix for all output files produced\n");\
	printf( "<trim> (optional, numeric): specify the number of bases to keep\n");\
	printf( "[RC] (optional): if \"RC\" is specified, reverse-complement seqs\n");\
	printf( "[D] <x.y> (optional): Denoise [duplicates x, kmer duplicates/100 y]\n");\
	printf( "Note: using .y discards reads with k-mers < y*1000 duplicates.\n");\
	printf( "LOG: optional, outputs failures.\n"); \
	return 2;\
}

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
char WORDTEMP[PACKSIZE+1] = {0};


typedef struct 
#if PACKSIZE<64
__attribute__ ((__packed__)) 
#endif
{
	WTYPE word;
	uint32_t ix;
	uint16_t length;
} SortBlock2;

typedef struct 
//#if PACKSIZE<64
__attribute__ ((__packed__)) 
//#endif
{
	WTYPE word;
	uint64_t count;
} KMer;

typedef struct KMerX KMerX;
struct 
//#if PACKSIZE<64
__attribute__ ((__packed__)) 
//#endif
KMerX {
	WTYPE word;
	uint64_t count;
	KMerX *left, *right;
	//size_t left, right;
};

// Explicit thread memory management
KMerX ***KBANK = 0; 
size_t KBANK_MAXK = 10000, KBANK_INITBINS = 100;
size_t *KBANK_BIN =0, *KBANK_BINCNT = 0, *KBANK_IX = 0;

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

inline char * num2word(WTYPE num, char * word) {
	int go = 0; for (; go < PACKSIZE; go++) {
		WTYPE temp = (WTYPE)num >> RSHFT;
		word[go] = X2C[temp];
		num <<= 2;
	}
	return word;
}

inline char * num2wordRC(WTYPE num, char * word) {
	int go = PACKSIZE-1; for (; go > -1; go--) {
		WTYPE temp = (WTYPE)num >> RSHFT;
		word[go] = X2C_RC[temp];
		num <<= 2;
	}
	return word;
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
	char *ref_s = *p, *key_s = key;
	while (*ref_s == *key_s++) if (!*ref_s++) return p - String;
	return -1;
	//return p - String; // replace last 3 lines for unsave ver
}
#define arr_t KMer
inline int uBST(arr_t *array, size_t sz, WTYPE key, size_t *ix){
	arr_t *p = array;
	while (sz) {
		size_t w = sz >> 1;
		if((p + w + 1)->word <= key) { p += w + 1; sz -= w + 1; } 
		else sz = w;
	}
	//return *p==key ? p - array : -1;
	*ix = p - array; return p->word == key;
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

void isv(KMer *array, size_t sz) {
    for (size_t x=1; x<sz; ++x) {
        for (size_t y=x; y && array[y-1].word > array[y].word; y--) {
            KMer temp = array[y];
            array[y] = array[y-1];
            array[y-1] = temp;
        }
    }
} 
void radix_sortC(KMer *array, size_t sz, int depth) {
    int x,y;
	KMer temp, value;
    uint64_t last[256] = { 0 }, pointer[256];

    for (x=0; x<sz; ++x) {
        ++last[(array[x].word >> depth) & 0xFF];
    }

    pointer[0] = 0;
    for (x=1; x<256; ++x) {
        pointer[x] = last[x-1];
        last[x] += last[x-1];
    }
    for (x=0; x<256; ++x) {
        while (pointer[x] != last[x]) { // last == start in other
            value = array[pointer[x]];
            y = (value.word >> depth) & 0xFF;
            while (x != y) {
                temp = array[pointer[y]];
                array[pointer[y]++] = value;
                value = temp;
                y = (value.word >> depth) & 0xFF;
            }
            array[pointer[x]++] = value;
        }
    }

    if (depth > 0) {
        depth -= 8;
        for (x=0; x<256; ++x) {
            y = x > 0 ? pointer[x] - pointer[x-1] : pointer[0];
            if (y > 256) 
                radix_sortC(array + pointer[x] - y, y, depth);
        }
    }
}

int wCmp(w1, w2) register const void *w1, *w2; {
	if (*(WTYPE *)w1 < *(WTYPE *)w2) return -1;
	return *(WTYPE *)w1 > *(WTYPE *)w2;
}

int kCmp(blk1, blk2) register const void *blk1, *blk2; {
	if (((KMer *)blk1)->word < ((KMer *)blk2)->word) return -1;
	return ((KMer *)blk1)->word > ((KMer *)blk2)->word;
}

inline size_t kmerDedupe(KMer ** Kmers, size_t num, size_t oldSz) {
	// New paradigm. Sort new chunk only, merge-sort new and old,
	// then resize the array with new bounds
	KMer * newBatch = *Kmers + oldSz, *oldBatch = *Kmers;
	size_t newSz = num - oldSz;
	#ifdef USE_QSORT
		#define KMER_LT(a,b) (a->word < b->word)
		QSORT(KMer, newBatch, newSz, KMER_LT);
	#else
		radix_sortC(newBatch, newSz, sizeof(WTYPE)*8 - 8);
		isv(newBatch, newSz);
	#endif
	
	// now the new batch is sorted. Start the merge
	// loop thru new batch, stick into old
	KMer *temp = malloc(num*sizeof(KMer));
	if (!temp) { puts("memory error in kmerDedupe"); return 0; }
	
	size_t new_ix=0, old_ix=0, k=0; // sum=0;
	for (; new_ix < newSz && old_ix < oldSz; k++) {
		if (newBatch[new_ix].word < oldBatch[old_ix].word)
			temp[k] = newBatch[new_ix++];
		else temp[k] = oldBatch[old_ix++];
	}
	// most likely one of the two will need appending on the end
	for (; new_ix < newSz; k++) temp[k] = newBatch[new_ix++];
	for (; old_ix < oldSz; k++) temp[k] = oldBatch[old_ix++];
	
	k = 0; 
	for (size_t i=1; i < num; ++i) {
		if (temp[i-1].word != temp[i].word) {
			temp[++k] = temp[i];
		}
		else temp[k].count += temp[i].count;
	}
	free(*Kmers); *Kmers=temp;
	return k + 1;
}
int resizeKmers(KMer ** AllKmers, size_t sz, size_t old, size_t *numK) { 
	// "old" contains previously sorted size. Can be 0
	*numK = kmerDedupe(AllKmers,*numK, old); 
	*AllKmers = realloc(*AllKmers,sz * sizeof(KMer)); //KMer
	if (!*AllKmers) { puts("Error in k-resize"); return 0; } 
	return 1;
} 

inline KMerX * xalloc(int thread, WTYPE word) { // KBANK,KBANK_INITBINS,KBANK_MAXK,KBANK_BIN,KBANK_IX
	KMerX *Kptr = KBANK[thread][KBANK_BIN[thread]] + KBANK_IX[thread];
	*Kptr = (KMerX){word,1,0,0};
	if (++KBANK_IX[thread] == KBANK_MAXK) { // reset the ix, increment bin
		KBANK_IX[thread] = 0;
		if (++KBANK_BIN[thread] == KBANK_BINCNT[thread]) { // resize bin array
			KBANK[thread] = realloc(KBANK[thread],
				sizeof(*KBANK[thread])*(KBANK_BINCNT[thread]*=2));
			if (!KBANK[thread]) { puts("ERROR: xalloc 1"); exit(3); }
			for (size_t x=KBANK_BINCNT[thread]/2; x<KBANK_BINCNT[thread]; ++x) {
				KBANK[thread][x] = malloc(KBANK_MAXK*sizeof(*KBANK[thread][x]));
				if (!KBANK[thread][x]) { puts("ERROR: xalloc 2"); exit(3); }
			}
		}
	}
	return Kptr;
}

void rexalloc(int threads) {
	for (int i = 0; i < threads; ++i) {
		KBANK[i] = realloc(KBANK[i],sizeof(*KBANK[i])*KBANK_BIN[i]);
		KBANK[i][KBANK_BIN[i]] = realloc(KBANK[i][KBANK_BIN[i]],
			sizeof(*KBANK[i][KBANK_BIN[i]]) * KBANK_IX[i]);
	}
}

// returns whether new node was created; a counter
int deTree(KMerX *tree, WTYPE word) { 
	do {
		if (word > tree->word) { // go right
			if (!tree->right) {
				KMerX *new = malloc(sizeof *new);
				*new = (KMerX){word,1,0,0};
				tree->right = new;
				return 1;
			}
			tree = tree->right; 
		}
		else if (word < tree->word) { // go left
			if (!tree->left) {
				KMerX *new = malloc(sizeof *new);
				*new = (KMerX){word,1,0,0};
				tree->left = new;
				return 1;
			}
			tree = tree->left;
		}
	} while (word != tree->word);
	++tree->count;
	return 0;
}

// returns whether new node was created; a counter
int xeTree(KMerX *tree, WTYPE word, int T) { 
	do {
		if (word > tree->word) { // go right
			if (!tree->right) {
				tree->right = xalloc(T,word);
				return 1;
			}
			tree = tree->right; 
		}
		else if (word < tree->word) { // go left
			if (!tree->left) {
				tree->left = xalloc(T,word);
				return 1;
			}
			tree = tree->left;
		}
	} while (word != tree->word);
	++tree->count;
	return 0;
}

// for repopulating an existing tree
void reTree(KMerX *tree, KMerX *node) { 
	for (;;) {
		if (node->word > tree->word) { // go right
			if (!tree->right) {
				node->left = 0; node->right = 0;
				tree->right = node;
				return;
			}
			tree = tree->right; 
		}
		else { // go left
			if (!tree->left) {
				node->left = 0; node->right = 0;
				tree->left = node;
				return;
			}
			tree = tree->left;
		}
	}
}

// for merging existing trees (returns if new node added)
int meNode(KMerX *tree, KMerX *node) { 
	do {
		if (node->word > tree->word) { // go right
			if (!tree->right) {
				node->left = 0; node->right = 0;
				tree->right = node;
				return 1; // node->count;
			}
			tree = tree->right; 
		}
		else if (node->word < tree->word) { // go left
			if (!tree->left) {
				node->left = 0; node->right = 0;
				tree->left = node;
				return 1; // node->count;
			}
			tree = tree->left;
		}
	} while (node->word != tree->word);
	tree->count += node->count;
	return 0;
}
// find in tree
size_t fiTree(KMerX *tree, WTYPE word) {
	do {
		if (word > tree->word) { // go right
			if (!tree->right) return 0;
			tree = tree->right; 
		}
		else if (word < tree->word) { // go left
			if (!tree->left) return 0;
			tree = tree->left;
		}
	} while (word != tree->word);
	return tree->count;
}
// get in tree (known existence)
size_t giTree(KMerX *tree, WTYPE word) {
	do {
		if (word > tree->word) tree = tree->right; 
		else if (word < tree->word) tree = tree->left;
	} while (word != tree->word);
	return tree->count;
}
// merge trees
void meTree(KMerX *tree, KMerX *tree2, size_t *totals) {
	KMerX *left = tree2->left, *right = tree2->right;
	*totals += meNode(tree,tree2);
	if (left) meTree(tree, left, totals);
	if (right) meTree(tree, right, totals);
}

void traceBalance(KMerX *tree, KMerX **array, size_t *ix) {
	
	if (tree->left) traceBalance(tree->left, array, ix);
	array[(*ix)++] = tree; // if on top, DFS. If mid, IOS, if bot: LFS
	if (tree->right) traceBalance(tree->right, array, ix);
	
}
void traceCnt(KMerX *tree, size_t *ix) {
	
	if (tree->left) traceCnt(tree->left, ix);
	++*ix;
	if (tree->right) traceCnt(tree->right, ix);
	
}
void buildBalanceL(KMerX *tree, KMerX **array, size_t sz);
void buildBalanceR(KMerX *tree, KMerX **array, size_t sz);
#define BUILDBALANCE() \
	if (!sz) { \
		CHILD = *array; \
		CHILD->left = 0; \
		CHILD->right = 0; \
		return; \
	} \
	size_t ix = sz >> 1; \
	CHILD = array[ix]; \
	if (ix) buildBalanceL(CHILD,array,ix-1); \
	else CHILD->left = 0; \
	buildBalanceR(CHILD,array+(ix+1), sz-(ix+1));

// set a branch of the given tree, and recurse with that branch as root
void buildBalanceL(KMerX *tree, KMerX **array, size_t sz) {
	#define CHILD tree->left
		BUILDBALANCE()
	#undef CHILD
}
void buildBalanceR(KMerX *tree, KMerX **array, size_t sz) {
	#define CHILD tree->right
		BUILDBALANCE()
	#undef CHILD
}

void aufbau(KMerX *tree, KMerX **array, size_t sz) {
	if (!sz) { reTree(tree,*array); return; }
	size_t ix = sz >> 1;
	reTree(tree,array[ix]);
	if (ix) aufbau(tree,array,ix-1);
	aufbau(tree,array+(ix+1),sz-(ix+1));
}
void freeTree(KMerX *tree) {
	if (tree->left) freeTree(tree->left);
	if (tree->right) freeTree(tree->right);
	free(tree);
}
void traceTree(KMerX *tree) {
	if (tree->left) traceTree(tree->left);
	printf("%s\t%I64u\n",num2word(tree->word,WORDTEMP),tree->count);
	if (tree->right) traceTree(tree->right);
}
void traceTreeDetail(KMerX *tree, int depth) {
	printf("%d\t%s\t%I64u\n",depth, num2word(tree->word,WORDTEMP),tree->count);
	if (tree->left) traceTreeDetail(tree->left, depth+1);
	if (tree->right) traceTreeDetail(tree->right, depth+1);
}
size_t buildDepth(KMerX *node, int depth, int *depthMax, size_t *depthTot,
size_t *count) {
	if (node->left) buildDepth(node->left,depth+1,depthMax,depthTot,count);
	if (depth > *depthMax) *depthMax = depth;
	++(*count);
	(*depthTot) += depth;
	//printf("%s\t%I64u\n",num2word(tree->word,WORDTEMP),tree->count);
	if (node->right) buildDepth(node->right,depth+1,depthMax,depthTot,count);
}
void reportAvMaxDepth(KMerX *tree) {
	int depthMax = 0;
	size_t count = 0, depthTot=0;
	buildDepth(tree,1,&depthMax, &depthTot, &count);
	double depthAv = (double)(depthTot)/count;
	printf("Total nodes = %lu. Max depth=%d, Avg=%f\n",count,depthMax,depthAv);
}

int tfs_cmp(const void *a, const void *b) {
	KMerX *b1 = *(KMerX **)a, *b2 = *(KMerX **)b;
	return (b1->count > b2->count) ? -1 : (b1->count < b2->count);
}
void treeFreqSort(KMerX **arr, size_t n) {
#ifdef USE_QSORT
	#define NODEFREQGT(a,b) ((*a)->count > (*b)->count)
	QSORT(KMerX*, arr, n, NODEFREQGT);
#else
	qsort(arr, n, sizeof(*arr), tfs_cmp);
#endif
}
int tns_cmp(const void *a, const void *b) {
	KMerX *b1 = *(KMerX **)a, *b2 = *(KMerX **)b;
	return (b1->word < b2->word) ? -1 : (b1->word > b2->word);
}
void treeNameSort(KMerX **arr, size_t n) {
#ifdef USE_QSORT
	#define NODENAMELT(a,b) ((*a)->word < (*b)->word)
	QSORT(KMerX*, arr, n, NODENAMELT);
#else
	qsort(arr, n, sizeof(*arr), tns_cmp);
#endif
}

KMerX * balanceTree(KMerX *tree, size_t sz, size_t totalCount) {
	// set limits
	#define MAX_NODES 1000000
	#define TOP_SHIFT 5
	if (sz > MAX_NODES) return tree;
	size_t ix = 0;
	KMerX **array = malloc(sizeof(*array) * (sz+1));
	traceBalance(tree, array, &ix);
	
	// experimental intervention: frequency-first tree construction
	treeFreqSort(array, sz+1);
	
	// Adaptive threshold determination
	size_t limit = 0;
	size_t thres = 0, cutoff = totalCount >> TOP_SHIFT;
	while ((thres += array[limit++]->count) < cutoff);
	//while (++limit <= sz && array[0]->count/array[limit]->count < 4); 
	//printf("limit = %lld\n",limit);
	tree = array[0];
	tree->left = 0; tree->right=0;
	
	if (limit > 2) {
		for (size_t i = 1; i < limit; ++i) reTree(tree,array[i]);
		
		// balance the top
		KMerX **top = malloc(sizeof(*top)* limit);
		for (size_t i = 0; i < limit; ++i) top[i] = array[i];
		treeNameSort(top, limit);
		ix = (limit-1)/2;
		tree = top[ix];
		buildBalanceL(tree, top, ix-1);
		buildBalanceR(tree, top + (ix+1), limit - 1 - (ix+1));
		free (top);
	} else limit = 1;
	size_t limit2 = (sz+1);
	
	// Add in the rest
	//for (size_t i = limit; i < limit2; ++i) reTree(tree,array[i]);
	
	// Add in the rest v2 (stacade)
	//size_t limit3 = (sz+1); limit2 = (sz+1)/2; // comment out to enable bottomBalance
	if (limit2 > limit) {
		int L = 1;
		for (long long i = limit; i < limit2-1; ++i) {
			reTree(tree,array[i+L]);
			L = -L;
		}
		if (L==-1) reTree(tree,array[limit2-2]);
		else reTree(tree,array[limit2-1]);
	}
	
	free(array);
	return tree;
}

inline void clumpParachute(KMerX **Roots, WTYPE *Clumps, size_t *NumsInserted,
size_t *TotalCounts, size_t *BalanceThreshes, size_t length) {
	//printf("here we go...\n");
	#pragma omp parallel for schedule(dynamic,1000) 
	for (int i = 0; i < length; ++i) {
		int tid = 0;
		#ifdef _OPENMP
		tid = omp_get_thread_num();
		#endif
		++TotalCounts[tid];
		
		if (!Roots[tid]->count) { 
			*Roots[tid] = (KMerX){Clumps[i],1,0,0}, NumsInserted[tid]=1; 
			continue;
		}
		
		NumsInserted[tid] += xeTree(Roots[tid],Clumps[i],tid); 
		if (NumsInserted[tid] >= BalanceThreshes[tid]) { 
			//printf("Balancing tree (tid %d) at %lu\n",tid,NumsInserted[tid]); 
			Roots[tid] = balanceTree(Roots[tid],NumsInserted[tid]-1, TotalCounts[tid]); 
			/* reportAvMaxDepth(root); */
			/* if (balanceThres==65535) {traceTreeDetail(root,0); 
			return 0;  } */ 
			BalanceThreshes[tid]=(BalanceThreshes[tid]+1)*2-1;  
		
		}
	}
}

KMerX* mergeParachutes(KMerX **Roots, int T, size_t *NumsInserted, 
size_t *TotalCounts, size_t *numInserted, size_t *totalCount) {
	*totalCount = *TotalCounts; // *TotalCounts;
	*numInserted = *NumsInserted;
	for (size_t i = 1; i < T; ++i) {
		*totalCount += TotalCounts[i];
		meTree(*Roots,Roots[i],numInserted);
	}
	//traceCnt(*Roots,numInserted);
	return *Roots;
}

size_t findRarestK(KMerX *tree, WTYPE *seq, uint16_t length) {
	size_t numPacks = length/PACKSIZE;
	//if (numPacks * PACKSIZE < length) ++numPacks;
	size_t min = (size_t)-1, cur;
	for (int i = 0; i < numPacks; ++i) {
		cur = giTree(tree, seq[i]);
		//printf("%llu, ", cur);
		if (cur < min) min = cur;
	}
	
	/* if (numPacks * PACKSIZE < length) { // handle endpiece
		unsigned mod = length % PACKSIZE; // guarantee: never 0
		// rightshift by 2xmodulo
		WTYPE prev = (seq[numPacks-1] << (2*mod)) + (seq[numPacks] >> (2*(PACKSIZE-mod)));
		cur = fiTree(tree, prev);
		if (cur < min) min = cur;
	} */
	//printf("\n");
	return min;
}

inline size_t findRarestK2(KMerX *tree, WTYPE *seq, uint16_t length) {
	size_t min = giTree(tree,*seq), cur;
	unsigned offset = 0, basePack = 1;
	for (int i = 1, b=length-PACKSIZE+1; i < b; ++i) {
		if (++offset == PACKSIZE) 
			cur = giTree(tree, seq[basePack]), 
			++basePack, offset = 0;
		else cur = giTree(tree, (*(seq+basePack-1) << (offset << 1)) + 
				(*(seq+basePack) >> ((PACKSIZE-offset) << 1)));
		//printf("%llu [%s=%u:%u], ", cur,num2word(this,WORDTEMP),offset,basePack);
		//printf("[%llu] ",cur);
		if (cur < min) min = cur;
	}
	//printf("MIN=%llu\n",min);
	return min;
}

/* inline size_t findRarestK3(KMerX **Trees, WTYPE *seq, uint16_t length, int T) {
	int limit = length - PACKSIZE + 1;
	size_t *Curs = malloc(sizeof(*Curs) * (T*limit));
	#pragma omp parallel for schedule(static,1) 
	for (int z=0; z<T; ++z) {
		int tid = 0;
		#ifdef _OPENMP
		tid = omp_get_thread_num();
		#endif
		KMerX *tree = Trees[z];
		size_t *Cur = Curs + z*limit;
		*Cur = fiTree(tree,*seq);
		
		unsigned offset = 0, basePack = 1;
		for (int i = 1, b=length-PACKSIZE+1; i < b; ++i) {
			if (++offset == PACKSIZE) 
				Cur[i] = fiTree(tree, seq[basePack]), 
				++basePack, offset = 0;
			else Cur[i] = fiTree(tree, (*(seq+basePack-1) << (offset << 1)) + 
					(*(seq+basePack) >> ((PACKSIZE-offset) << 1)));
			//printf("%llu [%s=%u:%u], ", cur,num2word(this,WORDTEMP),offset,basePack);
			//printf("[%llu] ",cur);
			//if (cur < Mins[z]) Mins[z] = cur;
		}
	}
	//printf("MIN=%llu\n",min);
	size_t min = (size_t)-1, sum;
	for (int j = 0; j < limit; ++j) {
		
		//size_t *Cur = Curs + j*limit;
		sum = Curs[j]; for (int k = 1; k < T; ++k) sum += Curs[k*limit + j];
		if (sum < min) min = sum;
		
	}
	free(Curs);
	return min;
} */

int main( int argc, char *argv[] ) {
	clock_t start; double cpu_time_used; start = clock(); // profiler
	// Debugging statements
	#ifdef DEBUG
		printf("type size=%u, shift=%u, pack=%u\n", sizeof(WTYPE), RSHFT, PACKSIZE);
		printf("max int size=%u/%u\n",sizeof(unsigned),sizeof(uint64_t));
		printf("Size of SortBlock2=%u\n",sizeof(SortBlock2));
	#endif
	if ( argc < 3 || argc > 9 ) PRINT_USAGE();
	char *inputFilename = argv[1], *prefixStr = argv[2];
	char *fasta_sx = "_filt.fa", *db_sx = ".db", *dp_sx = "_dupes.txt",
		 *filt_sx = "_filtered.txt";
	char *outputFasta = calloc(1,1+strlen(prefixStr)+strlen(fasta_sx)), 
		 *outputDB = calloc(1,1+strlen(prefixStr)+strlen(db_sx)),
		 *outputDP = calloc(1,1+strlen(prefixStr)+strlen(dp_sx)),
		 *outputFL = calloc(1,1+strlen(prefixStr)+strlen(filt_sx));
	strcpy(outputFasta,prefixStr); strcpy(outputFasta+strlen(prefixStr),fasta_sx);
	strcpy(outputDB,prefixStr); strcpy(outputDB+strlen(prefixStr),db_sx);
	strcpy(outputDP,prefixStr); strcpy(outputDP+strlen(prefixStr),dp_sx);
	strcpy(outputFL,prefixStr); strcpy(outputFL+strlen(prefixStr),filt_sx);
	/* printf("strings (prefix %s):\n%s\n%s\n%s\n",prefixStr, outputFasta,outputDB,outputDP);
	fopen(outputFasta,"wb");
	fopen(outputDB,"wb");
	fopen(outputDP,"wb");
	return 0; */
	FILE *fp = fopen(inputFilename, "rb");
	if (fp == NULL) { puts("Invalid input FASTA"); return 2; }
	FILE *off = fopen(outputFasta, "wb"), *ofd = fopen(outputDB,"wb");
	if (!off || !ofd) { puts("Invalid output prefix"); return 2; }
	FILE *ofdp=0, *ofdpF = 0;
	size_t trim = UINT16_MAX;
	int doRC = 0, doLog = 0; 
	double filt_i = 0.f; int copyNumThres = 0; // denoisers (k-filter not yet implemented)
	int numThreads = 1;
	if (strcmp(argv[argc-1],"ST")) { // enable MT if not "ST"
		#ifdef _OPENMP
			numThreads = omp_get_max_threads();
		#endif
	} else --argc;
	if (argc > 3 && !strcmp(argv[argc-1],"LOG")) { // enable MT if not "ST"
		doLog = 1;
		ofdp=fopen(outputDP,"wb");
		ofdpF=fopen(outputFL,"wb");
		if (!ofdp || !ofdpF) { puts("Invalid output prefix"); exit(2); }
		puts("Log writing enabled.");
		--argc;
	}
	#ifdef _OPENMP
		omp_set_num_threads(numThreads);
	#else
		numThreads = 1;
	#endif
	// Denoises at default intensity
	if (argc > 3 && !strcmp(argv[argc-1],"D")) { 
		filt_i = 2.f; 
		printf("Performing NINJA k-mer denoising at DEFAULT intensity: %.0f k-mers\n",filt_i);
		--argc;
	}
	// Denoises at specified intensity in the form x.y
	else if (argc > 3 && !strcmp(argv[argc-2],"D")) {
		filt_i = atof(argv[argc-1]);
		if (filt_i < 0) printf("Invalid denoising intensity (expect #REPS[.###Kmers]).\n"); 
		else {
			if (filt_i >= 1.f) {
				copyNumThres = filt_i;
				filt_i -= copyNumThres;
				printf("Performing NINJA replicon-denoising at %u compacted reads.\n", copyNumThres);
			}
			if (filt_i) { // Use the decimal remainder as kmer denoising
				printf("Performing NINJA k-mer denoising at %.0f k-mers\n", filt_i*1000.f);
				filt_i *= 1000;
				if (copyNumThres) ++copyNumThres;
			}
		}
		argc -= 2;
	}
	int (*cmpF)(WTYPE *, WTYPE *, uint16_t, uint16_t) = 
		copyNumThres ? &ycmp : &zcmp; //zcmp can replace xcmp
	if (!copyNumThres) copyNumThres = filt_i ? -1 : 1; 
	if (argc > 3 && !strcmp(argv[argc-1],"RC")) {
		printf("Reverse complementing the sequences.\n");
		doRC = 1; --argc;
	}
	// Flags for truncation after specified base
	if (argc == 4) { 
		trim = atoi(argv[argc-1]);
		printf("Trimming input sequences to %d bases.\n", trim);
	}
	
	C2Xb = calloc(128,sizeof(WTYPE));
	C2Xb['a'] = 0; C2Xb['A'] = 0; 
	C2Xb['c'] = 1; C2Xb['C'] = 1; 
	C2Xb['g'] = 2; C2Xb['G'] = 2;
	C2Xb['t'] = 3; C2Xb['T'] = 3;
	
	size_t numElem = 1000, ns=0;
	char **Samples = malloc(numElem*sizeof(char *));
	char **SeqIDs = doLog ? malloc(numElem*sizeof(*SeqIDs)) : 0;
	WTYPE **ReadsX = malloc(numElem*sizeof(WTYPE *));
	uint16_t *Sizes = calloc(numElem,sizeof(uint16_t));
	char *line = malloc(LINELEN + 1); // read up to  65k
	
	// prepare array and initial array size for kmer filtering
	/* KMer *AllKmers = 0; size_t kArraySize = 100000, 
		kQuantum= kArraySize, numK = 0, oldNumK = 0, kIX; */

	// MT versions of k-denoisers
	size_t queuedClumps = 0, fireThres = 1000000;
	#define BAL_THRES 255
	KMerX **Roots = 0;
	WTYPE *Clumps = 0; 
	size_t *NumsInserted=0, *TotalCounts=0, *BalanceThreshes=0;
	//KMerX **KBANK = 0; // GLOBAL VARIABLES
	// size_t KBANK_MAXK = 1000, *KBANK_BIN =0, *KBANK_IX = 0;
	
	if (filt_i) {
		printf("Number of threads for k-mer denoise: %d\n",numThreads);
		//AllKmers = malloc(kArraySize*sizeof(KMer));
		
		Roots = malloc(numThreads*sizeof(*Roots));
		Clumps = malloc(fireThres*sizeof(*Clumps));
		NumsInserted = calloc(numThreads,sizeof(*NumsInserted));
		TotalCounts = calloc(numThreads,sizeof(*TotalCounts));
		BalanceThreshes = malloc(numThreads*sizeof(*BalanceThreshes));
		KBANK = malloc(numThreads*sizeof(*KBANK));
		KBANK_BIN = calloc(numThreads,sizeof(*KBANK_BIN));
		KBANK_BINCNT = malloc(numThreads*sizeof(*KBANK_BINCNT));
		KBANK_IX = calloc(numThreads,sizeof(*KBANK_IX));
		
		for (int i = 0; i < numThreads; ++i) {
			Roots[i] = malloc(sizeof(*Roots[i]));
			*Roots[i] = (KMerX){0,0,0,0};
			BalanceThreshes[i] = BAL_THRES;
			
			KBANK[i] = malloc(KBANK_INITBINS*sizeof(*KBANK[i]));
			KBANK_BINCNT[i] = KBANK_INITBINS;
			for (int j=0; j < KBANK_INITBINS; ++j) {
				KBANK[i][j] = malloc(KBANK_MAXK*sizeof(*KBANK[i][j])); // init this bin's kmers
				if (!KBANK[i][j]) {puts("error: xalloc 0"); exit(3); }
			}
		}
	}
	
	while (line = fgets(line,LINELEN,fp)) { 
		if (ns == numElem) {
			numElem *= 2;
			Samples = realloc(Samples,numElem * sizeof(char *));
			ReadsX = realloc(ReadsX, numElem * sizeof(WTYPE *));
			Sizes = realloc(Sizes, numElem*sizeof(uint16_t));
			if (!Samples || !ReadsX || !Sizes) {puts("Error in resize"); return 3;}
			memset(Sizes+numElem/2 + 1,0,(numElem/2-1)*sizeof(uint16_t));
			if (doLog) {
				SeqIDs = realloc(SeqIDs, numElem * sizeof(*SeqIDs));
				if (!SeqIDs) {puts("Error in resize"); return 3;}
			}
		}
		// copy in the sample name up to _ or null minus 1
		char *src = line + 1;
		
		while (*src != '_' && *src != ' ' && *src != '\n') ++src; 
		if (doLog) { // also trace until whitespace for sample id
			char *seqID = src;
			while (*seqID != ' ' && *seqID != '\n') ++seqID;
			SeqIDs[ns] = malloc(seqID - src + 1);
			if (!SeqIDs[ns]) {puts("Out of memory for SeqIDs"); return 3;}
			char *d = SeqIDs[ns];
			char *b = src; while (b < seqID) *d++ = *b++;
			*d = 0;
		}
		Samples[ns] = malloc(src - line);
		if (!Samples[ns]) {puts("Not enough Samples[ns] mem"); return 3;}
		
		char *dest = Samples[ns]; 
		char *beginSample = line + 1; while (beginSample < src) 
			*dest++ = *beginSample++;
		*dest = 0;

		// copy in the encoded sequence
		if (!(line = fgets(line,LINELEN,fp))) 
			{ puts("Error reading file."); return 2; }
		src = line;
		
		register size_t length = strlen(src);
		if (src[length-1] == '\n') --length; // lop off newline(s)
		if (src[length-1] == '\r') --length; // supports every platform!
		if (trim < length) length = trim;
		if (length >= UINT16_MAX) { 
			printf("Warning: truncating read %llu.\n",ns);
			length = UINT16_MAX - 1;
		}
		size_t numPacks = length/PACKSIZE;
		if (numPacks * PACKSIZE < length) ++numPacks;
		
		Sizes[ns] = length; 
		ReadsX[ns] = malloc(numPacks*sizeof(WTYPE));
		if (!ReadsX[ns]) {puts("Bad ReadsX[ns] mem"); return 3; }
		
		WTYPE *thisPack = ReadsX[ns];
		WTYPE clump = C2Xb[*src++]; 
		
		#define GENERATE_WORD_PRE() \
		for (int k = 1, z = 2; k < length; ++k, ++z) { \
			clump <<= 2u; \
			clump += C2Xb[*src++]; \
			if (z == PACKSIZE) *thisPack++ = clump, z = 0; 
		#define GENERATE_KMER() \
		if (k + 2 > PACKSIZE) { \
			Clumps[queuedClumps++] = clump; \
			if (queuedClumps == fireThres) { \
				clumpParachute(Roots,Clumps,NumsInserted, \
					TotalCounts,BalanceThreshes,fireThres); \
					queuedClumps = 0; \
			} \
		} 
		#define GENERATE_WORD_POST() } 
		
		if (filt_i) 
			GENERATE_WORD_PRE()
			GENERATE_KMER()
			GENERATE_WORD_POST()
		else
			GENERATE_WORD_PRE()
			GENERATE_WORD_POST()
		numPacks *= PACKSIZE;
		if (numPacks > length) *thisPack++ = clump << ((numPacks - length) << 1);
		++ns;
	}
	KMerX *master = 0;
	if (filt_i) { 
		if (queuedClumps) clumpParachute(Roots,Clumps,NumsInserted, 
			TotalCounts,BalanceThreshes,queuedClumps); 
		size_t numInserted = 0, totalCount = 0;
		
		if (numThreads > 1) {
			master = mergeParachutes(Roots, numThreads, 
				NumsInserted, TotalCounts,&numInserted, &totalCount);
			//master = balanceTree(master,numInserted-1, totalCount);
			//master = quickBalance(master,numInserted-1);
		}
		else {
			master = *Roots; 
			numInserted = *NumsInserted; totalCount = *TotalCounts;
		}
		//traceTreeDetail(*Roots,0);
		printf("Distinct K-mers found: %lu, Total k-mers: %llu\n",numInserted,totalCount);
		#ifdef DEBUG
			reportAvMaxDepth(master);
		#endif
		rexalloc(numThreads);
	}
	
	fclose(fp);
	free(line);
	// Shrink data structures for more memory
	Samples = realloc(Samples,ns * sizeof(char *));
	ReadsX = realloc(ReadsX, ns * sizeof(WTYPE *));
	Sizes = realloc(Sizes, ns * sizeof(uint16_t));
	/* if (doLog) SeqIDs = realloc(SeqIDs, ns *sizeof(*SeqIDs));
		for (unsigned i = 0; i < ns; ++i) 
		fprintf(ofdp,"%s\n",SeqIDs[i]);
	exit(3); */
	printf("Num of sequences: %u\n",ns);
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
	int64_t i_copyThres = filt_i ? copyNumThres-1 : INT64_MAX;
	size_t filt_n = filt_i;
#ifdef DEBUG
	printf("copyNumThres=%d, copyThres=%llu, filt_i=%f [%u]\n", copyNumThres,i_copyThres,filt_i,filt_n);
#endif
	#define WRITE_SUPPORTED_DUPE() {\
		if (copies >= copyNumThres || (copies >= i_copyThres && \
			findRarestK(master, ReadsX[prevIX], Sizes[prevIX]) >= filt_n)) { \
			/* printf("\nfound rarest K=%llu\n",findRarestK2(master, ReadsX[prevIX], Sizes[prevIX])); */ \
			if (doLog) { \
				/* while (++lastLogged <= k) */ \
				for (unsigned w = lastLogged; w < k; ++w) \
					++committed, fprintf(ofdp,"%s%s\t",Samples[SeqIX[w]],SeqIDs[SeqIX[w]]); \
				fprintf(ofdp,"\n"); \
				lastLogged = k; \
			} \
			for (int y = 0; y < x; ++y) \
				if (Counts[y]) fprintf(ofd,"%u:%u:",y,Counts[y]), Counts[y] = 0; \
			fprintf(ofd,"\n"); \
			if (doRC) fprintf(off,">%u\n%s\n",rix++, decodeStringXRC(ReadsX[prevIX], \
				Sizes[prevIX],word,string)); \
			else fprintf(off,">%u\n%s\n", rix++, decodeStringX(ReadsX[prevIX], \
				Sizes[prevIX],word,string)); \
		} \
		else { \
			if (doLog) { \
				for (unsigned w = lastLogged; w < k; ++w) \
					++rejected, fprintf(ofdpF,"%s%s\n",Samples[SeqIX[w]],SeqIDs[SeqIX[w]]); \
				lastLogged = k; \
			} \
			memset(Counts,0,x*sizeof(unsigned)); \
		} \
		copies = 1; \
	}
	size_t committed = 0, rejected = 0;
	unsigned copies = 1, dupes = 0, rix=0;
	char *string = malloc(UINT16_MAX), *word = calloc(PACKSIZE+1,1);
	unsigned prevIX, thisIX, lastLogged = 0;
	for (k=1; k < ns; ++k) {
		prevIX = SeqIX[k-1]; thisIX = SeqIX[k];
		++Counts[crBST(Samples[prevIX],x-1,SmpDD)];
		if (cmpF(ReadsX[prevIX],ReadsX[thisIX],Sizes[prevIX], Sizes[thisIX])) 
			WRITE_SUPPORTED_DUPE()
		else { ++copies; ++dupes; }
	}
	prevIX = thisIX;
	++Counts[crBST(Samples[prevIX],x-1,SmpDD)]; // add last count
	WRITE_SUPPORTED_DUPE()
	if (doLog) printf("Number rejected = %llu, committed = %llu\n",
		rejected,committed);
#ifdef PROFILE
	printf("->Mapping and file writing: %f\n", 
		((double) (clock() - start)) / CLOCKS_PER_SEC); start = clock();
#endif
	// todo: free more, add k-mer filtering
	free (SeqIX);
	free (string);
	return 0;
}