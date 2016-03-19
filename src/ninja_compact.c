#define _FILE_OFFSET_BITS 64
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#define LINELEN 1000000
#define PRINT_USAGE() { \
	printf( "\nNINJA-OPS Coelescescence. v1.5.0a\n"); \
	printf( "Usage: ninja_coelescescence in_aligns.sam in_concatesome.fa out_newAligns.sam\n");\
	exit(1);\
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

int ycmp(register const char *str1, register const char *str2) { 
	while (*str1 == *str2++) if (!*str1++) return 0;  
	return *str1 && *(str2 - 1); 
} 

int main( int argc, char *argv[] )
{
	FILE *in, *db, *out;
	if ( argc < 4 || argc > 4 ) PRINT_USAGE()
	in = fopen(argv[1], "rb");
	db = fopen(argv[2], "rb");
	out = fopen(argv[3], "wb");
	if (!in || !out || !db) {
		fputs("Can't open input/output file(s)!.\n",stderr);
		exit(1);
	}
	// parse sam
	size_t arr_sz = 1000, ns = 0;
	unsigned long long *Sam_ptrs = malloc(arr_sz*sizeof(*Sam_ptrs));
	uint16_t *Sam_lens = malloc(arr_sz*sizeof(*Sam_lens));
	char *lineO = malloc(LINELEN+1), *line; 
	while (line = fgets(lineO,LINELEN,in)) {
		if (ns == arr_sz) {
			arr_sz *= 2;
			Sam_ptrs = realloc(Sam_ptrs, arr_sz*sizeof(*Sam_ptrs));
			Sam_lens = realloc(Sam_lens, arr_sz*sizeof(*Sam_lens));
			if (!Sam_ptrs || !Sam_lens) {
				fputs("Error: Out of memory.\n",stderr); exit(1); }
		}
		int tabs = 0; while (tabs < 3) if (*line++=='\t') ++tabs; // at pos
		if (*line == '0') continue; // sequence didn't match
		Sam_ptrs[ns] = (size_t)atoll(line); 
		while (tabs < 9) if (*line++=='\t') ++tabs; // at seq
		char *startSeq = line;
		while (tabs < 10) if (*line++=='\t') ++tabs;
		Sam_lens[ns] = line - startSeq - 1;
		++ns;
	}
	
	// parse DB
	fseek(db,0,SEEK_END); size_t sz = ftell(db); rewind(db);
	char *db_txt = malloc(sz+1);
	size_t db_len = fread(db_txt, sz, 1, db);
	db_txt[sz] = 0;
	char *db_seq = strchr(db_txt,'\n'); // has ">ACGCGA..."
	fclose(db);
	
	//construct refs
	char **Refs = malloc(ns*sizeof(*Refs)), 
		***RefsP = malloc(ns*sizeof(*RefsP));
	for (size_t i = 0; i < ns; ++i) {
		Refs[i] = malloc(Sam_lens[i]+1);
		Refs[i] = strncpy(Refs[i],db_seq+Sam_ptrs[i],Sam_lens[i]);
		Refs[i][Sam_lens[i]] = 0;
		RefsP[i] = Refs + i;
		//fprintf(out,"%llu: %s\n",Sam_ptrs[i],Refs[i]);
	}
	free(db_txt);
	twrqs(RefsP,ns,0);
	//for (size_t i = 0; i < ns; ++i) fprintf(out,"%llu: %s\n",Sam_ptrs[RefsP[i]-Refs],*RefsP[i]);
	for (size_t i = 1; i < ns; ++i) {
		if (!ycmp(*RefsP[i-1],*RefsP[i])) Sam_ptrs[RefsP[i]-Refs] = Sam_ptrs[RefsP[i-1]-Refs];
	}
	//for (size_t i = 0; i < ns; ++i) fprintf(out,"%llu: %s\n",Sam_ptrs[RefsP[i]-Refs],*RefsP[i]);
	
	// Reparse and replace
	rewind(in); ns = 0;
	while (line = fgets(lineO,LINELEN,in)) {
		char *begin = line;
		int tabs = 0; while (tabs < 3) if (*line++=='\t') ++tabs; // at pos
		if (*line == '0') {
			fputs(begin,out);
			continue; 
		}
		*line = 0;
		//fputs(begin,out); //line till here
		fprintf(out,"%s%llu",begin,Sam_ptrs[ns]);
		while (tabs < 4) if (*++line=='\t') ++tabs; // at seq
		fputs(line,out);
		++ns;
	}
	free(lineO);
	return 0;
}