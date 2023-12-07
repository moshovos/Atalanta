// A. Moshovos 2021
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PROBS 16
#define PROB_BITS 10
#define DEPTH_MAX 2

typedef struct pte_t {
	int vmin;
	int off;
	int abits; // symbol encoding bits
	int obits; // offset encoding bits
	int vcnt; // how many values
} pte_t;

static const int pcnt = PROBS;
int vmax = 0;
pte_t pnew[PROBS+1];
pte_t pbest[PROBS+1];

static int verbose = 1;

static int
lg(int i) {
  if (i <= 1) return 0;
  return 31 - __builtin_clz(i-1)+1;
// i |= (i >> 1); i |= (i >> 2); i |= (i >> 4); i |= (i >> 8); i |= (i >> 16); return i - (i >> 1);
}

static void
pt_off_set(pte_t *ptp) {
  int i;

  for (i = 0; i < PROBS; i++) {
    int dist = ptp[i+1].vmin - ptp[i].vmin;
    ptp[i].off = lg(dist);
  }
}

static void
pt_init(pte_t *ptp, int vmax) {
  int i;

  int vstep = vmax / PROBS;
  for (i = 0; i <= PROBS; i++)
    ptp[i].vmin = i * vstep;

  pt_off_set(ptp);
}

static void
pt_print_final(pte_t *ptp) {
  int i;
  int tbits;
  int vcnt = 0;
  int abits = 0;
  int obits = 0;

  tbits = 0;
  fprintf(stderr, "PT_FINAL: vmin off abits obits ab100 ob100 tb100 vcnt frac\n");
  for (i = 0; i <=PROBS; i++) {
     tbits += ptp[i].abits + ptp[i].obits;
     abits += ptp[i].abits;
     obits += ptp[i].obits;
     vcnt = ptp[i].vcnt;
  }
  
  for (i = 0; i <=PROBS; i++) {
    
    fprintf(stderr, "[%3d, %2d] : %10d %10d %.3f %.3f %.3f %10d %.3f\n", ptp[i].vmin, ptp[i].off, ptp[i].abits, ptp[i].obits, ptp[i].abits/(double) tbits, ptp[i].obits/(double) tbits, (ptp[i].abits+ptp[i].obits)/(double)tbits, ptp[i].vcnt, ptp[i].vcnt/(double) vcnt);
  }
  fprintf(stderr, "\n");
}

static void
pt_print(pte_t *ptp) {
  int i;

  fprintf(stderr, "PT_INIT: ");
  for (i = 0; i <=PROBS; i++)
    fprintf(stderr, "[%d, %d (%d)] ", ptp[i].vmin, ptp[i].off, ptp[i].vcnt);
  fprintf(stderr, "\n");
}

static float 
entropy_precision (float f) {
  //float t = f;
  f *= (1 << PROB_BITS);
  f = round(f);
  f /= (1 << PROB_BITS);
  //fprintf (stderr, "ep %f --> %f\n", t, f);
  if (f == 0) return 0;
  return log2(f);
}

static float
pt_encoded_size(int *hist, pte_t *ptp) {
   int i, v, ptotal, ototal;
   float btotal;
   int rcnt[PROBS+1];

   ptotal = 0;
   ototal = 0;
   pt_off_set(ptp);
   //printf ("PTENCODEDSIZE: ");
   //pt_print(ptp);

   for (i = 0; i < PROBS; i++) {
      rcnt[i] = 0;
      ptp[i].obits = 0;
      for (v = ptp[i].vmin; v < ptp[i+1].vmin; v++) {
         rcnt[i] += hist[v];        
         ototal += hist[v] * ptp[i].off;
         ptotal += hist[v];
	 ptp[i].obits += hist[v] * ptp[i].off;
         
	 //fprintf(stderr, "V %d r %d o %d p %d\n", i, rcnt[i], ototal, ptotal);
      }
      ptp[i].vcnt = rcnt[i];
   }

   btotal = 0;
   for (i = 0; i < PROBS; i++) {
	   //fprintf(stderr, "EP <%d> %d:\n", i, rcnt[i]);
      btotal += rcnt[i] * entropy_precision(((float) rcnt[i])/ (float) ptotal);
      ptp[i].abits = -(rcnt[i] * entropy_precision((rcnt[i])/ (float) ptotal) ); 
   }
   //fprintf(stderr, "o %d p %d b %f\n", ototal, ptotal, btotal);
   return ototal - btotal;
}

static void inline pt_copy(pte_t *d, pte_t *s) {
  int i;
  for (i = 0; i <= PROBS; i++) d[i] = s[i];
}

float
search_try(int *hist, pte_t *trial_in, float *score_best, pte_t *ptbest, int depth, int around) {
   pte_t trial[PROBS+1];
   int i;
   int c;

   pt_copy (trial, trial_in);

   for (c = 1; c < PROBS; c++) {
      if ((around >= 0) && (abs(c - around) != 1)) continue;
        while (trial[c].vmin > trial[c-1].vmin) {
          trial[c].vmin -= 1;
          if (depth < DEPTH_MAX) search_try(hist, trial, score_best, ptbest, depth+1, c);
          else {
            float score_new = pt_encoded_size(hist, trial);
	    //pt_print(trial);
            if (score_new < *score_best) {
                pt_copy(ptbest, trial);
                *score_best = score_new;
		if (verbose == 1) { fprintf(stderr, "PTBEST: "); pt_print(ptbest); }
            }
          }
        }
        while (trial[c].vmin < trial[c+1].vmin) {
          trial[c].vmin += 1;
          if (depth < DEPTH_MAX) search_try(hist, trial, score_best, ptbest, depth+1, c);
          else {
            float score_new = pt_encoded_size(hist, trial);
	    //pt_print(trial);
            if (score_new < *score_best) {
                pt_copy(ptbest, trial);
                *score_best = score_new;
		if (verbose == 1) { fprintf(stderr, "PTBEST: "); pt_print(ptbest); }
            }
          }
        }
   }
}

#include <time.h>

clock_t start, end;
double cpu_time_used;

void
search(int bits, int *hist, pte_t *ptab, int in_verbose) { 
    int i;
    float score_best, prev_best, score_raw;

    //start = clock();
    vmax = 1<<bits;

    // for (i=0; i<256; i++) fprintf(stderr, "%d->%d,",i,lg(i));
    //verbose = 0;
    pt_init(pbest, vmax);
    score_raw = score_best = pt_encoded_size(hist, pbest);
    if (verbose == 1) pt_print(pbest);
    
    while (1) {
       score_best = pt_encoded_size(hist, pbest);
       for (i = 0; i <= PROBS; i++) pnew[i] = pbest[i];
       prev_best = score_best;
       if (verbose == 1) printf ("ENCODED: %f\n", score_best);
       search_try(hist, pnew, &score_best, pbest, 2, -2);
       if (verbose == 1) printf ("ENCODED: %f\n", score_best);
       if (score_best/prev_best > 0.99) break;
    }
    //end = clock();
    //cpu_time_used = ((double) (end - start));
    //fprintf(stderr, "CSEARCH: %f clocks %f clocks/sec\n", cpu_time_used, CLOCKS_PER_SEC);


   if (verbose > 1) pt_print_final(pbest);
   if (verbose > 1) printf ("ENCODED: %f\n", score_best);
   for (i = 0; i< PROBS; i++) ptab[i] = pbest[i];
   if (verbose) fflush(stderr);

    return;

    for (i = 0; i < (1<<bits); i++)
    {
       printf ("%d %d\n", i, hist[i]);
    }
}
