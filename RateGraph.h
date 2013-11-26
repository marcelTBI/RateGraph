#ifndef __BHGBUILD_H
#define __BHGBUILD_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#include <vector>
#include <queue>
#include <string>
#include <map>
#include <set>
#include <algorithm>

extern "C" {
  #include "findpath.h"
}

#include "RNAutils.h"
#include "RNAstruc.h"

using namespace std;

enum type1 {INTER_CLUSTER, REPRESENT, CRIT_EDGE, NEW_FOUND};
const char type1_str[][14] = {"INTER_CLUSTER", "REPRESENT", "CRIT_EDGE", "NEW_FOUND"};
const int type1_len = 4;

class DSU {

public:
  // sequence
  char *seq;
  short *s0;
  short *s1;

   // temperature
  double _kT;

  // structures
  vector<RNAlocmin> LM;  // contains memory (after split - should be in both programs)
  int number_lm;  // number of lm in the beginning
  int gl_maxen;

  // -- linkCP
    // vertex sets
    map<RNAstruc, int> vertex_l;  // points to number in LM
    //map<RNAstruc, int> vertex_s;  // points to number in saddles

    vector<RNAsaddle> saddles; // contains memory

    // edge sets
    set<edgeLL> edges_l; // LM to LM
    set<edgeSS> edges_s; // saddle to saddle
    set<edgeLS> edges_ls; // i is LM, j is saddle

    // edges for graph search
    vector< set<edgeLL> > edgesV_l;

  // -- Height-first Search
  // mapping for not full matrices
  vector<int> mapping;
  vector<int> mapping_rev;

private:
  DSU() {};

public:
  DSU(FILE *input); // read seq + structs from input (barriers/RNAlocmin output)
  ~DSU();

public:
  // big ones

  // visualisation
  void VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug);
  vector<SimplePath> ConstructAllPaths(int source, int dest, int max_length, int threshold);
  void ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold);
  void PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual, bool print_energies); // print dot file to filename, dot_prog - use dot or neato?; print - print dot output to file_print, visual - use tree for visualisation

  // print text
  void PrintLinkCP(FILE *output = stdout, bool fix = true);
  void PrintLM(FILE *output = stdout, bool fix = true);
  void PrintSaddles(FILE *output = stdout, bool fix = true);
  void PrintBarr(FILE *output = stdout);

  // print files:
  void PrintMatrix(char *filename, bool full, char type); // print matrices (E - energy, D - distance, G - graph distance)
  void PrintRates(char *filename, bool full, double temp, char mode);

  // sort results and fix connections, output TRUE if recomputed connections
  bool SortFix();

  // small
  int Size() {return LM.size();}

  // graph techniques
  // -- Height-first Search  -- returns energy barriers to get i-th minima from start minima + distances in graph
  vector<std::pair<int, int> > HeightSearch(int start, vector< set<edgeLL> > &edgesV_l);

  // prints optimal path between start and stop to filename.
  void GetPath(int start, int stop,  vector< set<edgeLL> > &edgesV_l, char *filename, int maxkeep);
  void GetPath(int start, int stop, int maxkeep = 0);

  // evaluation
  void EHeights(FILE *heights, bool full);
  void ERank(FILE *rank, bool barrier, bool out_conns = false);

  void Histo(FILE *histo);
};

struct Rate {
  double rate;
  int saddle_num;

  Rate operator+(const Rate second) {
    Rate res(rate + second.rate, saddle_num);
    if (rate < second.rate) res.saddle_num = second.saddle_num;
    return res;
  }

  Rate(int r, int s) {
    rate = r;
    saddle_num = s;
  }
  Rate() {
    rate = 0;
    saddle_num = -1;
  }
};
struct pq_rate {
  int conns;
  int lm;
  int energy;

  pq_rate(int c, int l, int e) {
    conns = c;
    lm = l;
    energy = e;
  }

  bool operator<(const pq_rate &scnd) const {
    if (conns == scnd.conns) {
      return energy < scnd.energy;
    } else return conns > scnd.conns;
  }
};

void SwapMinima(int dim_rates, double *rate_matrix, int a, int b);
void MxRShorten(double **shorten, int fulldim, int gdim);

class RateGraph
{
  // rates
  vector < map<int, double> > rates;

  // removed:
  set <int> removed;
  map <int, int> lm_to_pos;
  vector<int> pos_to_lm;

  // int
  int edge_count;

  // priority queue for removal
  priority_queue<pq_rate> to_rem;

  // rate_matrix for shur removal
  double *rate_matrix;
  unsigned int dim_rates;  // dimension of the matrix rate_matrix

  // local minima
  vector<RNAlocmin> lms;
  set<RNAstruc> filter;

  // redirection of local minima

public:
  RateGraph(DSU &dsu, double temp);
  ~RateGraph();

  int ReadFilter(char *filename);
  int ConstructQueue(char order, int number_remove);

  void PrintDot(FILE *dot);
  void PrintDot(char *filename, bool to_eps);

  int RemoveOne(int rem_lm);
  int RemoveX(int x, int stop_fraction = 50, bool reeval = false);
  int RemoveShur(int x, int step);

  void PrintRates(char *filename);
  void PrintRates(FILE *filname);

  void PrintOutput(FILE *output);
};
#endif
