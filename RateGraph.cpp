#include <stdio.h>
#include <string.h>
#include <stdlib.h>

extern "C" {
  #include "utils.h"
  #include "fold_vars.h"
  #include "pair_mat.h"

  #include "fold.h"
  #include "findpath.h"
  #include "move_set.h"
  #include "mxccm.h"
}

#include "RateGraph.h"
#include "hash_util.h"

#include <algorithm>

using namespace std;

inline double rate(int en_from, int en_to, double _kT) {
  return 1.0*exp(-(en_to-en_from)/100.0/_kT);
}

RateGraph::RateGraph(DSU &dsu, double temp)
{
  double _kT = 0.00198717*(273.15 + temp);

  rates.resize(dsu.Size());

  for (unsigned int i=0; i<dsu.edgesV_l.size(); i++) {
    for (auto j=dsu.edgesV_l[i].begin(); j!=dsu.edgesV_l[i].end(); j++) {
      // extract info
      int en = j->en;
      int lm1 = j->i;
      int lm2 = j->j;
      double rate1 = rate(dsu.LM[lm1].energy, en, _kT);
      double rate2 = rate(dsu.LM[lm2].energy, en, _kT);

      // construct Rate-s
      //Rate rt1 (rate1, j->saddle);
      //Rate rt2 (rate2, j->saddle);

      // add 'em
      rates[lm1][lm2] = rate1;
      rates[lm2][lm1] = rate2;
      edge_count++;
    }
  }


  // copy local minima:
  lms.resize(dsu.LM.size());
  for (unsigned int i=0; i<dsu.LM.size(); i++) {
    lms[i] = dsu.LM[i];
    lms[i].structure = allocopy(dsu.LM[i].structure);
    lms[i].str_ch = NULL;
    lms[i].recompute_str();
  }

  // maps:
  pos_to_lm.resize(lms.size());
  for (unsigned int i=0; i<lms.size(); i++) {
    pos_to_lm[i] = i;
    lm_to_pos[i] = i;
  }

  // null it.
  rate_matrix = NULL;
}

RateGraph::~RateGraph()
{
  if (rate_matrix) free(rate_matrix);
  for (auto it=filter.begin(); it!=filter.end(); it++) {
    free(it->structure);
    if (it->str_ch) free(it->str_ch);
  }
  for (unsigned int i=0; i<lms.size(); i++) {
    free(lms[i].structure);
    if (lms[i].str_ch) free(lms[i].str_ch);
  }
}

int RateGraph::ConstructQueue(char order, int number_remove)
{
  switch (order) {
  case 'C':
    // make a priority queue with function as a function of comparison
    for (unsigned int i=0; i<rates.size(); i++) {
      to_rem.push(pq_rate(rates[i].size(), i, lms[i].energy));
    }
    break;
  case 'E':
    // the same, but insert there only high energy ones. assume sorted order.
    for (unsigned int i=rates.size()-1; i>=rates.size()-number_remove; i--) {
      to_rem.push(pq_rate(rates[i].size(), i, lms[i].energy));
    }
    break;
  case 'F':
    // make a priority queue with function as a function of comparison
    for (unsigned int i=0; i<rates.size(); i++) {
      if (filter.count(lms[i])==0) {
        to_rem.push(pq_rate(rates[i].size(), i, lms[i].energy));
      }
    }
    break;
  default: return 1;
  }
  return 0;
}

int RateGraph::RemoveOne(int rem_lm)
{
  int count = 0;

  // calculate sum
  double sum = 0.0;
  for (auto it=rates[rem_lm].begin(); it!=rates[rem_lm].end(); it++) {
    sum += it->second;
  }

  // go through doubles
  for (auto it=rates[rem_lm].begin(); it!=rates[rem_lm].end(); it++) {

    double rR_1 = it->second;
    int lm1 = it->first;
    auto rat_it = rates[lm1].find(rem_lm);
    double r1_R = rat_it->second;
    rates[lm1].erase(rat_it);
    count--;

    auto it2 = it; it2++;
    for (; it2!=rates[rem_lm].end(); it2++) {
      // we have rates rem_lm -> lm1; rem_lm> lm2
      // get others
      double rR_2 = it2->second;
      int lm2 = it2->first;
      auto rat_it = rates[lm2].find(rem_lm);
      double r2_R = rat_it->second;
      //rates[lm2].erase(rat_it);

      //join 'em
      double r1_2 = r1_R*rR_2/sum;
      double r2_1 = r2_R*rR_1/sum;

      // add them to rates:
      auto f12 = rates[lm1].find(lm2);
      auto f21 = rates[lm2].find(lm1);
      if (f12!=rates[lm1].end()) {
        f12->second += r1_2;
      } else {
        count++;
        rates[lm1][lm2] = r1_2;
      }
      if (f21!=rates[lm2].end()) {
        f21->second += r2_1;
      } else {
        count++;
        rates[lm2][lm1] = r2_1;
      }
    }
  }

  // erase edges:
  count -= (int)rates[rem_lm].size();
  rates[rem_lm].clear();

  // return how many edges inserted
  return count;
}

int RateGraph::RemoveX(int x, int stop_fraction, bool reeval)
{
  int count = 0;
  while (count!=x) {

    // get lowest
    pq_rate pqr = to_rem.top(); to_rem.pop();

    // stopping criterion:
    if (pqr.conns * stop_fraction > (int)rates.size()-count) {
      to_rem.push(pqr);
      break;
    }

    // check it:
    if (reeval && pqr.conns != (int)rates[pqr.lm].size()) {
      pqr.conns = (int)rates[pqr.lm].size();
      to_rem.push(pqr);
      continue;
    }

    // debug:
    if (count%10000 == 0 || (pqr.conns>100 && count%1000==0) || (pqr.conns>150 && count%100==0)) fprintf(stderr, "removing node %8d (%8d remaining) (conns=%4d, energy=%6.2f)\n", count, (int)rates.size()-count, pqr.conns, pqr.energy/100.0);

    /*// debug
    char name[100];
    sprintf(name, "test%d.dot", (int)rates.size()-count);
    PrintDot(name, true);
*/
    // remove it:
    int new_edges = RemoveOne(pqr.lm);
    removed.insert(pqr.lm);
    edge_count += new_edges;
    count++;
  }

  /*// debug
  char name[100];
  sprintf(name, "test%d.dot", (int)rates.size()-count);
  PrintDot(name, true);
*/
  // reduce size of rates vector:
  // update maps:
  auto it_rem = removed.begin();
  for(unsigned int i=0, j=0; i<rates.size(); i++) {
    if (it_rem==removed.end() || *it_rem != (int)i) {
      lm_to_pos[i] = j;
      j++;
    } else {
      if (it_rem!=removed.end()) it_rem++;
    }
  }
  pos_to_lm.resize(lm_to_pos.size());
  for (auto it=lm_to_pos.begin(); it!=lm_to_pos.end(); it++) {
    pos_to_lm[it->second] = it->first;
  }

  // rewrite rates vector:
  for (unsigned int i=0; i<rates.size(); i++) {
    map<int, double> tmp;
    for (auto it=rates[i].begin(); it!=rates[i].end(); it++) {
      tmp.insert(make_pair(lm_to_pos[it->first], it->second));
    }
    rates[i] = tmp;
  }
  // and remove empty lines
  int new_i = 0;
  for(unsigned int i=0; i<rates.size()-removed.size(); i++, new_i++) {
    while(removed.count(new_i)) new_i++;
    rates[i] = rates[new_i];
  }
  rates.resize(rates.size()-removed.size());

  /*// print the conversion
  for (unsigned int i=0; i<pos_to_lm.size(); i++) {
    fprintf(stderr, "%8d -> %8d\n", pos_to_lm[i]+1, i+1);
  }*/

  // return how many removed
  return count;
}

void SwapMinima(int dim_rates, double *rate_matrix, int a, int b)
{
  for (int i=0; i<dim_rates; i++) {
    //fprintf(stderr, "swapping %5d %5d %d (%d)\n", a, b, dim_rates, i );
    // change lines:
    swap(rate_matrix[a*dim_rates+i], rate_matrix[b*dim_rates+i]);
  }
  for (int i=0; i<dim_rates; i++) {
    // change columns
    swap(rate_matrix[i*dim_rates+a], rate_matrix[i*dim_rates+b]);
  }
}

void MxRShorten(double **shorten, int fulldim, int gdim)
{
  //does: shortened = GG - GB*BB^(-1)*BG, where matrix tmp_rates is split as:
  //tmp_rates = (GG | GB)
  //            (BG | BB)
  // GG has dimension gdim*gdim; tmp_rates fulldim*fulldim
  // create matrices:

  int bdim = fulldim - gdim;
  int i,j;

  //double *gg = (double *)calloc(gdim*gdim,sizeof(double));
  double *bg = (double *)calloc(bdim*gdim,sizeof(double));
  double *bb = (double *)calloc(bdim*bdim,sizeof(double));
  double *gb = (double *)calloc(gdim*bdim,sizeof(double));

  // first we need to fix the diagonal entries tmp_rates[i][i] = sum_j tmp_rates[i][j]
  double *tmp_rates = *shorten;
  for (i = 0; i < fulldim; i++) tmp_rates[fulldim*i+i] = 0.0;
  for (i = 0; i < fulldim; i++) {
    double tmp = 0.00;
    // calculate row sum
    for(j = 0; j < fulldim; j++)  tmp += tmp_rates[fulldim*i+j];
    tmp_rates[fulldim*i+i] = -tmp;
  }

  // fill the matrices: (row = i; column = j)
  for (i=0; i<bdim; i++) {
    for (j=0; j<gdim; j++) {
      bg[gdim*i+j] = tmp_rates[fulldim*(i+gdim)+j];
    }
  }

  for (i=0; i<gdim; i++) {
    for (j=0; j<bdim; j++) {
      gb[bdim*i+j] = tmp_rates[fulldim*i+j+gdim];
    }
  }

  for (i=0; i<bdim; i++) {
    for (j=0; j<bdim; j++) {
      bb[bdim*i+j] = tmp_rates[fulldim*(i+gdim)+j+gdim];
    }
  }

  /*MxFPrintD(tmp_rates, "Q", my_dim, my_dim, stderr);
  MxFPrintD(gg, "GG", dim, dim, stderr);
  MxFPrintD(bg, "BG", bdim, dim, stderr);
  MxFPrintD(gb, "GB", dim, bdim, stderr);
  MxFPrintD(bb, "BB", bdim, bdim, stderr);
*/
  // result2 = gb*bb^(-1)*bg
  minv(bb, bdim);
  //MxFPrintD(bb, "BBinv", bdim, bdim, stderr);
  double *result = (double *)calloc(gdim*bdim,sizeof(double));
  mmul_singular(result, gb, bb, gdim, bdim, bdim, 0);
  //MxFPrintD(result, "gb*bb-1", dim, bdim, stderr);
  double *result2 = (double *)calloc(gdim*gdim,sizeof(double));
  mmul_singular(result2, result, bg, gdim, bdim, gdim, 1);

  //if (opt.want_verbose) MxFPrintD(result2, "gb*bb-1*bg", gdim, gdim, stderr);

  // result2 = gg - result2
  for (i=0; i<gdim; i++) {
    for (j=0; j<gdim; j++) {
      result2[gdim*i+j] = tmp_rates[fulldim*i+j] - result2[gdim*i+j];
    }
  }

  //MxFPrintD(result2, "matrix after shortening", dim ,dim, stderr);
  free(result);
  free(*shorten);
  free(gb);
  free(bg);
  free(bb);
  *shorten = result2;
}

struct sort_struc {
  int number;
  bool to_remove;

  sort_struc(int number, bool to_remove) {
    this->number = number;
    this->to_remove = to_remove;
  }

  const bool operator<(const sort_struc &second) const {
    if (to_remove == second.to_remove) {
      return number < second.number;
    } else {
      return to_remove < second.to_remove;
    }
  }
};

int RateGraph::RemoveShur(int x, int step)
{
  //fprintf(stderr, "calling Shur %5d %5d\n", x, step );
  // get set to remove:
  /*set<int> to_remove;
  while ((int)to_remove.size() != x) {
    pq_rate tmp = to_rem.top(); to_rem.pop();
    to_remove.insert(lm_to_pos[tmp.lm]);
  }*/

  //build rate matrix:
  dim_rates = rates.size();
  if (rate_matrix) free(rate_matrix);
  rate_matrix = (double *) malloc(sizeof(double)*dim_rates*dim_rates);
  for (unsigned int i=0; i<dim_rates*dim_rates; i++) rate_matrix[i] = 0.0;

  for (unsigned int i=0; i<dim_rates; i++) {
    double sum = 0.0;
    for (auto a=rates[i].begin(); a!=rates[i].end(); a++) {
      rate_matrix[i*dim_rates + a->first] = a->second;
      sum += a->second;
    }
    rate_matrix[i*dim_rates + i] = -sum;
  }

  //PrintRates("before_reorder.rat");

  // now resort matrix to have the minima to remove on bottom.
  vector<sort_struc> tmp_sort;
  for (unsigned int i=0; i<dim_rates; i++) {
    tmp_sort.push_back(sort_struc(i, false));
  }
  int to_remove = 0;
  while (to_remove != x) {
    pq_rate tmp = to_rem.top(); to_rem.pop();
    tmp_sort[lm_to_pos[tmp.lm]].to_remove = true;
    to_remove++;
  }
  sort(tmp_sort.begin(), tmp_sort.end());
  // make inverse:
  vector<int> tmp_inv(tmp_sort.size());
  for (unsigned int i=0; i<tmp_sort.size(); i++) {
    tmp_inv[tmp_sort[i].number] = i;
  }

  // sort it by swapping (in place ordering):
  for (unsigned int i=0; i<dim_rates; i++) {
    while ((int)i!=tmp_inv[i]) {
      SwapMinima(dim_rates, rate_matrix, i, tmp_inv[i]);

      swap(pos_to_lm[i], pos_to_lm[tmp_inv[i]]);
      swap(lm_to_pos[pos_to_lm[tmp_inv[i]]], lm_to_pos[pos_to_lm[i]]);

      swap(tmp_inv[i], tmp_inv[tmp_inv[i]]);
    }
  }

/*
  // sort it temporary:
  {
    double *tmp_rates = (double*) malloc(dim_rates*dim_rates*sizeof(double));

    for (int i=0; i<dim_rates; i++) {
      for (int j=0; j<dim_rates; j++) {
        tmp_rates[i*dim_rates+j] = rate_matrix[tmp_sort[i].number*dim_rates+tmp_sort[j].number];
      }
    }
    swap(rate matrix, tmp_rates);
    free(tmp_rates);
  }

  // adjust mapping:
  for (int i=0; i<dim_rates; i++) {
    swap(pos_to_lm[i], pos_to_lm[tmp_sort[i].number]);
    swap(lm_to_pos[pos_to_lm[*it]], lm_to_pos[pos_to_lm[i]]);
  }*/

  //PrintRates("after_reorder.rat");

  //fprintf(stderr, "calling Shur %5d %5d\n", x, step );

  // finally run removal on that matrix
  int res_dim = dim_rates-x;
  while ((int)dim_rates > res_dim) {
    // reduce of:
    int reduce = min(step, (int)dim_rates-res_dim);
    // now Shur it!
    MxRShorten(&rate_matrix, dim_rates, dim_rates-reduce);
    dim_rates -= reduce;
    // report:
    fprintf(stderr, "reduced dimension from %8d to %8d\n", dim_rates+reduce, dim_rates);
  }

  // add those, that we have removed
  for (unsigned int i=res_dim; i<dim_rates+x; i++) {
    removed.insert(pos_to_lm[i]);
    lm_to_pos.erase(pos_to_lm[i]);
  }
  pos_to_lm.resize(res_dim);

  // lastly, rewrite the "rates" data structure:
  rates.clear();
  rates.resize(res_dim);
  for (int i=0; i<res_dim; i++) {
    for (int j=0; j<res_dim; j++) {
      if (rate_matrix[i*res_dim+j]>0 && i!=j) {
        rates[i][j] = rate_matrix[i*res_dim+j];
      }
    }
  }

  return x;
}

void RateGraph::PrintRates(char *filename)
{
  FILE *rate_file = fopen(filename, "w");
  if (rate_file) {
    PrintRates(rate_file);
    fclose(rate_file);
  }
}

void RateGraph::PrintRates(FILE *filname)
{
  if (!rate_matrix) {
      //build rate matrix:
    dim_rates = rates.size();
    rate_matrix = (double *) malloc(sizeof(double)*dim_rates*dim_rates);
    for (unsigned int i=0; i<dim_rates*dim_rates; i++) rate_matrix[i] = 0.0;

    for (unsigned int i=0; i<dim_rates; i++) {
      double sum = 0.0;
      for (auto a=rates[i].begin(); a!=rates[i].end(); a++) {
        rate_matrix[i*dim_rates + a->first] = a->second;
        sum += a->second;
      }
      rate_matrix[i*dim_rates + i] = -sum;
    }
    dim_rates = dim_rates;
  }

  // just print the rate matrix:
  for (unsigned int i=0; i<dim_rates; i++) {
    for (unsigned int j=0; j<dim_rates; j++) {
        fprintf(filname, "%11.5g ", rate_matrix[dim_rates*i+j]);
    }
    fprintf(filname, "\n");
  }
}

void RateGraph::PrintDot(char *filename, bool to_eps)
{
  FILE *dot = fopen(filename, "w");
  if (dot) {
    PrintDot(dot);
    fclose(dot);
  }

  if (to_eps) {
    // start neato/dot:
    char syst[200];
    char filename2[200];
    strcpy(filename2, filename);
    filename2[strlen(filename2)-4]='\0';
    sprintf(syst, "%s -Tps < %s > %s.eps", "dot", filename, filename2);
    system(syst);
  }
}

void RateGraph::PrintDot(FILE *dot)
{
  fprintf(dot, "Digraph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
  //nodes LM:
  for (unsigned int i=0; i<rates.size(); i++) {
    fprintf(dot, "\"%d\" [label=\"%d\"]\n", i+1, pos_to_lm[i]+1/*, lms[pos_to_lm[i]].energy*/);
  }
  fprintf(dot, "\n");

  for (unsigned int i=0; i<rates.size(); i++) {
    for (auto it=rates[i].begin(); it!=rates[i].end(); it++) {
      // edges l-l
      fprintf(dot, "\"%d\" -> \"%d\" [label=\"%6.4g\"]\n", i+1, (it->first)+1, it->second);
    }
  }
  fprintf(dot, "\n}\n");
}

int RateGraph::ReadFilter(char *filename)
{
  FILE *file_filt;
  file_filt = fopen(filename, "r");
  if (file_filt) {
    char *line;
    line = my_getline(file_filt);
    while (line) {
      bool empty_line = false;
      short *tmp = NULL;
      float energy_tmp;


      char *p;
      for (int i=0; i<4; i++) {
        p = strtok(i==0?line:NULL, " ");
        switch (i) {
        case 0:
        case 1:
          if (!p) {empty_line = true; break;}
          if (isStruct(p)) {
            if (tmp) free(tmp); // only one struct per line!
            tmp = make_pair_table(p);
            break;
          }
        case 2:
          if (tmp) sscanf(p, "%f", &energy_tmp);
          break;
        }
        if (empty_line) break;
      }
      if (tmp) {
        RNAstruc filtered;
        filtered.energy = energy_tmp;
        filtered.structure = tmp;
        filtered.recompute_str();
        filter.insert(filtered);
      }
      line = my_getline(file_filt);
    }

    fclose(file_filt);
  }
  return filter.size();
}

void RateGraph::PrintOutput(FILE *output)
{
  for (unsigned int i=0; i<pos_to_lm.size(); i++) {
    fprintf(output, "%5d\n", pos_to_lm[i]+1);
  }
}

DSU::DSU(FILE *input) {

  // NULL::
  seq = NULL;
  s0 = NULL;
  s1 = NULL;
  gl_maxen = INT_MIN;

  if (!input) return ;

  // read seq
  char *line;
  int num = 0;

  line = my_getline(input);
  char *seq2 = strtok(line, " ");
  if (!isSeq(seq2)) {
    free(line);
    line = my_getline(input);
    seq2 = strtok(line, " ");
    if (!isSeq(seq2)) {
      free(line);
      return ;
    }
  }
  seq = (char*) malloc((strlen(seq2)+1)*sizeof(char));
  strcpy(seq, seq2);
  free(line);

  //init
  make_pair_matrix();
  s0 = encode_sequence(seq, 0);
  s1 = encode_sequence(seq, 1);

  // read .dsu file

  //read structs
  line = my_getline(input);
  bool saddle_reading = false;
  while (line) {
    bool empty_line = false;
    short *tmp = NULL;
    float energy_tmp;
    int type;

    char *p;
    for (int i=0; i<4; i++) {
      p = strtok(i==0?line:NULL, " ");
      switch (i) {
      case 1:
        if (!p) {empty_line = true; break;}
        if (isStruct(p)) {
          if (tmp) free(tmp); // only one struct per line!
          tmp = make_pair_table(p);
          break;
        }
      case 2:
        sscanf(p, "%f", &energy_tmp);
        break;
      case 3:
        for (int i=0; i<4; i++) {
          if (!saddle_reading) {
            if (strcmp(p, LMtype_string[i])==0) {
              type = i;
              break;
            }
          } else {
            if (strcmp(p, SDtype_string[i])==0) {
              type = i;
              break;
            }
          }
        }
        break;
      }
      if (empty_line) break;
    }

    if (empty_line) {
      saddle_reading = true;
      free(line);
      line = my_getline(input);

      edgesV_l.resize(LM.size());
      continue;
    }

    vector<int> LM_c;
    vector<int> saddle_c;
    int number;
    if (saddle_reading) {
      while ((p=strtok(NULL, " "))) {
        if (p[strlen(p)-1]=='S') {
          sscanf(p, "%dS", &number);
          saddle_c.push_back(number-1);
        } else {
          sscanf(p, "%d", &number);
          LM_c.push_back(number-1);
        }
      }
    }

    if (saddle_reading && LM_c.size()<2) {
      fprintf(stderr, "File reading error -- too few LM connected by saddle %d\n", (int)saddles.size()+1);
      exit(EXIT_FAILURE);
    }

    // add info:
    if (tmp) {

      if (!saddle_reading) {
        RNAlocmin struc;
        struc.structure = tmp;
        struc.str_ch = pt_to_char(struc.structure);
        struc.energy = en_fltoi(energy_tmp);
        struc.type = (LMtype)type;

        // insert new LM
        vertex_l[struc] = LM.size();
        LM.push_back(struc);

        gl_maxen = max(gl_maxen, struc.energy);
      } else {
        RNAsaddle struc(LM_c[0], LM_c[1], (SDtype)type);
        struc.structure = tmp;
        struc.str_ch = pt_to_char(struc.structure);
        struc.energy = en_fltoi(energy_tmp);
        // insert new sadle
        saddles.push_back(struc);

        // edge l-l
        for (int i=0; i<(int)LM_c.size(); i++) {
          for (int j=i+1; j<(int)LM_c.size(); j++) {
            edgeLL e(LM_c[i], LM_c[j], struc.energy, saddles.size()-1);
            edges_l.insert(e);
            edgesV_l[LM_c[i]].insert(e);
            edgesV_l[LM_c[j]].insert(e);
          }
        }
        // edge l-s
        for (int i=0; i<(int)LM_c.size(); i++) {
          edgeLS e2(LM_c[i], saddles.size()-1);
          edges_ls.insert(e2);
        }

        // edges s-s
        for (int i=0; i<(int)saddle_c.size(); i++) {
          edgeSS e(saddle_c[i], saddles.size()-1);
          edges_s.insert(e);
        }

      }
    } else {
      free(tmp);
    }

    free(line);
    line = my_getline(input);
    num++;
  }

  // sort them!
  SortFix();
	number_lm = (int)LM.size();
	//printf("------------------------------------------------------------%d\n", number_lm);
}

DSU::~DSU() {
  for (unsigned int i=0; i<LM.size(); i++) {
    if (LM[i].structure) free(LM[i].structure);
    if (LM[i].str_ch) free(LM[i].str_ch);
  }
  if (seq) free(seq);
  if (s0) free(s0);
  if (s1) free(s1);
  LM.clear();

  for (unsigned int i=0; i<saddles.size(); i++) {
    if (saddles[i].structure) free(saddles[i].structure);
    if (saddles[i].str_ch) free(saddles[i].str_ch);
  }
}

void DSU::PrintDot(char *filename, bool dot_prog, bool print, char *file_print, bool visual, bool print_energies)
{
  // landmap not supported yet
  // start find_union stuff
  UF_set connected;

  int color = 180;
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
    //nodes LM:
    for (unsigned int i=0; i<LM.size(); i++) {
      char energy[20] = "";
      if (print_energies) sprintf(energy, "\\n%.2f", LM[i].energy/100.0);
      switch (LM[i].type) {
        case NORMAL:
        case NORM_CF: fprintf(dot, "\"%d\" [label=\"%d%s\"]\n", i+1, i+1, energy); break;
        case EE_DSU: fprintf(dot, "\"%d\" [label=\"%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(0, 0, 255), rgb(0, 0, 255)); break;
        case EE_COMP: fprintf(dot, "\"%d\" [label=\"%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(255, 0, 0), rgb(255, 0, 0)); break;
      }
    }
    fprintf(dot, "\n");

    // visualisation option (not finished -- currently it prints out only without saddles)
    if (visual) {
      connected.enlarge_parent(LM.size());
      set<edgeLL, edgeLL_compen> tmp;
      tmp.insert(edges_l.begin(), edges_l.end());
      for (set<edgeLL>::iterator it=tmp.begin(); it!=tmp.end(); it++) {
        if (!connected.joint(it->i, it->j)) {
          fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0);
          connected.union_set(it->i, it->j);
        }
      }
      fprintf(dot, "\n");

    } else {

      //nodes saddle:
      for (unsigned int i=0; i<saddles.size(); i++) {
        char energy[20] = "";
        if (print_energies) sprintf(energy, "\\n%.2f", saddles[i].energy/100.0);
        switch (saddles[i].type) {
          /*case DIRECT: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(color, color, color), rgb(color, color, color)); break;
          case LDIRECT: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(color+30, color+30, color+30), rgb(color+30, color+30, color+30)); break;
          case NOT_SURE: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(color-30, color-30, color-30), rgb(color-30, color-30, color-30)); break;
          case COMP: fprintf(dot, "\"S%d\" [label=\"S%d (%d %d)\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, saddles[i].lm1+1, saddles[i].lm2+1, rgb(255, color, color), rgb(255, color, color)); break;*/

          case DIRECT: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(color, color, color), rgb(color, color, color)); break;
          case LDIRECT: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(color+30, color+30, color+30), rgb(color+30, color+30, color+30)); break;
          case NOT_SURE: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(color-30, color-30, color-30), rgb(color-30, color-30, color-30)); break;
          case COMP: fprintf(dot, "\"S%d\" [label=\"S%d%s\", color=\"%s\", fontcolor=\"%s\"]\n", i+1, i+1, energy, rgb(255, color, color), rgb(255, color, color)); break;
        }
      }
      fprintf(dot, "\n");
      // edges l-l
      for (set<edgeLL>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
        bool component = (saddles[it->saddle].type == COMP);
        fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\", color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, it->en/100.0, (component?rgb(255, 0, 0):rgb(0, 0, 0)), (component?rgb(255, 0, 0):rgb(0, 0, 0)));
      }
      fprintf(dot, "\n");
      // edges l-s
      for (set<edgeLS>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
        bool component = (saddles[it->j].type == COMP);
        fprintf(dot, "\"%d\" -- \"S%d\" [color=\"%s\", fontcolor=\"%s\"]\n", (it->i)+1, (it->j)+1, (component?rgb(255, color, color):rgb(color, color, color)), (component?rgb(255, color, color):rgb(color, color, color)));
      }

      fprintf(dot, "\n");
      // edges s-s
      for (set<edgeSS>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
        fprintf(dot, "\"S%d\" -- \"S%d\" [color=\"%s\", fontcolor=\"%s\"]\n",(it->i)+1, (it->j)+1, rgb(color, color, color),rgb(color, color, color));
      }
    }
    fprintf(dot, "\n}\n");
  }

  fclose(dot);

  // start neato/dot:
  if (dot && print && file_print) {
    char syst[200];
    sprintf(syst, "%s -Tps < %s > %s", (dot_prog ? "dot" : "neato"), filename, file_print);
    system(syst);
    //printf("%s returned %d\n", syst, res);
  }
}

static int EN_BARRIERS = true;

struct pq_path {
  int lm;
  int dist;
  int en_barr;

  int value() {
    return (EN_BARRIERS ? en_barr : dist);
  }

  bool operator<(const pq_path &second) const {
    if (EN_BARRIERS) {
      if (en_barr != second.en_barr) {
        return en_barr<second.en_barr;
      }
    }
    if (dist == second.dist) {
      return lm<second.lm;
    }
    return dist<second.dist;
  }

  pq_path(int lm, int dist, int en_barr) {
    this->lm = lm;
    this->dist = dist;
    this->en_barr = en_barr;
  }
};

void DSU::VisPath(int src, int dest, bool en_barriers, int max_length, bool dot_prog, bool debug)
{
  EN_BARRIERS = en_barriers;

  char filename[50];
  char file_print[50];
  sprintf(filename, "path%d_%d.dot", src+1, dest+1);
  sprintf(file_print, "path%d_%d.eps", src+1, dest+1);

  float color = 0.5;
  //open file
  FILE *dot;
  dot = fopen(filename, "w");
  if (dot) {
    fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");

    // actual nodes:
      // forward pass:
    vector<int> LM_tmp(LM.size(), INT_MAX);
    priority_queue <pq_path> pq_tmp;
    pq_path to_insert(src, 0, LM[src].energy);
    LM_tmp[src] = to_insert.value();
    pq_tmp.push(to_insert);

    bool found_dst = false;
    int maximum = INT_MAX;
    while (!pq_tmp.empty()) {
      pq_path point = pq_tmp.top();
      pq_tmp.pop();

      if (point.lm == dest) {
        maximum = min(point.value(), maximum);
        found_dst = true;
        if (!EN_BARRIERS) {
          break;
        }
      }

      if (LM_tmp[point.lm] < point.value()) continue; // we have found better and this is obsolete

      for (set<edgeLL>::iterator it=edgesV_l[point.lm].begin(); it!=edgesV_l[point.lm].end(); it++) {
        int en_barr = it->en;
        int goesTo = it->goesTo(point.lm);
        if (found_dst && (maximum < en_barr)) continue; // we dont want higher energy (really dirty programming :/)

        // insert new elements into pq (if they are better)
        pq_path to_insert(goesTo, point.dist+1, max(point.en_barr, en_barr));
        if (LM_tmp[goesTo] <= to_insert.value()) continue;
        pq_tmp.push(to_insert);
        LM_tmp[goesTo] = to_insert.value();
      }
    }

    // output
    set<int> LM_out;
    LM_out.insert(dest);
    LM_out.insert(src);
    set<edgeLL> edge_out;

      // backward pass -- find all paths with same dist/same en_barrier
    if (EN_BARRIERS) {
      // collect all paths (NP-complete)
      vector<SimplePath> paths = ConstructAllPaths(src, dest, max_length, LM_tmp[dest]);

      int en_barr;
      if (paths.size()>0) {
        en_barr = paths[0].max_energy;

        for (unsigned int i=0; i<paths.size(); i++) {
          // stopping condition:
          if (paths[i].max_energy != en_barr) break;

          // print them?
          if (debug) paths[i].Print(true);

          // output them
          for (unsigned int j=0; j<paths[i].points.size(); j++) {
            LM_out.insert(paths[i].points[j]);
            if (j>0) {
              edgeLL e(paths[i].points[j-1], paths[i].points[j], paths[i].energies[j-1], -1);
              edge_out.insert(e);
            }
          }
        }
      }

    } else {
      // find the shortest path
      int max = LM_tmp[dest];
      if (max == INT_MAX) {
        fprintf(stderr, "WARNING: Cannot reach %d from %d!\n", src+1, dest+1);
      } else {
        queue<pq_path> que;
        que.push(pq_path(dest, max, INT_MAX));

        while (!que.empty()) {
          pq_path point = que.front();
          que.pop();

          for (set<edgeLL>::iterator it=edgesV_l[point.lm].begin(); it!=edgesV_l[point.lm].end(); it++) {
            int goesTo = it->goesTo(point.lm);

            // if this point is on shortest path:
            if (LM_tmp[goesTo] == point.dist-1) {
              LM_out.insert(goesTo);
              edge_out.insert(*it);

              if (point.dist>1) que.push(pq_path(goesTo, point.dist-1, INT_MAX));
            } // else nothing
          }
        }
      }
    }

    // nodes
    for (set<int>::iterator it=LM_out.begin(); it!=LM_out.end(); it++) {
      if ((*it) == dest || (*it) == src) {
        fprintf(dot, "\"%d\" [label=\"%d\"]\n", (*it)+1, (*it)+1);
      } else fprintf(dot, "\"%d\" [label=\"%d\", color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", (*it)+1, (*it)+1, color, color);
    }
    fprintf(dot, "\n");

    // edges:
    for (set<edgeLL>::iterator it=edge_out.begin(); it!=edge_out.end(); it++) {
      fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\", color=\"0.0 0.0 %.1f\", fontcolor=\"0.0 0.0 %.1f\"]\n", (it->i)+1, (it->j)+1, it->en/100.0, color, color);
    }

    fprintf(dot, "}\n");
  }

  fclose(dot);

  // start neato/dot:
  char syst[200];
  sprintf(syst, "%s -Tps < %s > %s", (dot_prog ? "dot" : "neato"), filename, file_print);
  system(syst);
  //printf("%s returned %d", syst, res);
}

vector<vector<std::pair<int, int> > > matrix;
bool generated = false;

void DSU::PrintMatrix(char *filename, bool full, char type)
{
  if (!full && mapping.size() == 0) {
    mapping_rev.resize(LM.size(), 0);
    fprintf(stderr, "Mapping: \n");
    for (int i=0; i<(int)LM.size(); i++) {
      if (LM[i].type==NORMAL) {
        mapping.push_back(i);
        mapping_rev[i] = mapping.size()-1;
        fprintf(stderr, "%4d[%4d] ", i+1, (int)mapping.size());
      }
    }
    fprintf(stderr, "(%d mapped to %d) \n", (int)LM.size(), (int)mapping.size());
  }

  int size = (full?LM.size():mapping.size());
  FILE *energies;
  energies = fopen(filename, "w");
  if (energies) {
    // create matrix
    if (type != 'D' && !generated) {
      matrix.clear();
      for (int i=0; i<size; i++) {
        int to_search = full?i:mapping[i];
        matrix.push_back(HeightSearch(to_search, edgesV_l));
        // discard those we dont want
        if (!full) {
          for (int j=0; j<(int)mapping.size(); j++) {
            matrix[i][j] = matrix[i][mapping[j]];
          }
        }
        if ((int)matrix[i].size()!=size) {
          matrix[i].resize(size);
        }
      }
      // resolve i==i
      for (unsigned int i=0; i<matrix.size(); i++) {
        matrix[i][i] = make_pair(LM[i].energy, 0);
      }
      generated = true;
    }

    // print
    for (unsigned int i=0; i<matrix.size(); i++) {
      fprintf(energies, "%6d %s ", mapping[i]+1, LM[mapping[i]].str_ch);
      for (unsigned int j=0; j<matrix[i].size(); j++) {
        switch (type) {
        case 'E': fprintf(energies, "%6.2f ", matrix[i][j].first/100.0); break;
        case 'D': fprintf(energies, "%5d ", HammingDist(LM[i].structure, LM[j].structure)); break;
        case 'G': fprintf(energies, "%5d ", matrix[i][j].second); break;
        }
      }
      fprintf(energies, "\n");
    }
  }
  fclose(energies);
}

void DSU::PrintRates(char *filename, bool full, double temp, char mode)
{
  int size = (full?LM.size():number_lm);
  FILE *rates;
  rates = fopen(filename, "w");
  if (rates) {
    mode_rates mod;
    switch (mode) {
      case 'V': mod = VERTEX_CONTR; break;
      case 'E': mod = EDGE_CONTR_MIN; break;
      case 'M': mod = EDGE_CONTR_MAX; break;
      case 'F': mod = NO_CONTR; break;
      case 'S': mod = VERTEX_CONTR_SUM; break;
      default: mod = VERTEX_CONTR;
    }
    Graph graph(edges_l, LM, mod);
    fprintf(stderr, "graph created...%d LM\n", (int)LM.size());
    if (mod!=NO_CONTR) {
      for(int i=LM.size()-1; i>=size; i--) {
        /*char filename[20];
        char filename_eps[20];
        sprintf(filename, "smth%d.dot", i);
        sprintf(filename_eps, "smth%d.eps", i);
        graph.PrintDot(filename, true, true, filename_eps);*/

        // call appropriate contraction
        if (mod == VERTEX_CONTR || mod == VERTEX_CONTR_SUM) {
          graph.RemoveLastPoint();
        } else if (mod == EDGE_CONTR_MAX || mod == EDGE_CONTR_MIN) {
          while (!graph.RemoveLowestEdge()) ;
        }
        if (i%100 ==0) fprintf(stderr, "removed point %d\n", i);
      }
      char filename[20];
		  char filename_eps[20];
		  sprintf(filename, "reduced%c.dot", mode);
		  sprintf(filename_eps, "reduced%c.eps", mode);
		  graph.PrintDot(filename, true, true, filename_eps);
    }
    graph.PrintRates(rates, temp);
  }
  fclose(rates);
}

vector<SimplePath> DSU::ConstructAllPaths(int source, int dest, int max_length, int threshold)
{
  vector<SimplePath> paths;

  SimplePath path;
  path.AddLast(source, INT_MIN);

  // construct all path recursively
  ConstructPath(paths, path, dest, max_length-1, threshold);

  // score them
  for (unsigned int i=0; i<paths.size(); i++) {
    paths[i].Score();
  }

  // sort
  sort(paths.begin(), paths.end());

  return paths;
}

void DSU::ConstructPath(vector<SimplePath> &paths, SimplePath &path, int dest, int max_length, int threshold)
{
  int num = path.points[path.points.size()-1];

  if (num == dest) {
    paths.push_back(path);
    //if (paths.size()%100==1) printf("found paths: %d\n", (int)paths.size());
    return ;
  }

  if (max_length == 0) return;

  // all edges
  for (set<edgeLL>::iterator it=edgesV_l[num].begin(); it!=edgesV_l[num].end(); it++) {
    int goesTo = it->goesTo(num);
    if (!path.ContainsNode(goesTo) && it->en <= threshold) {
      path.AddLast(goesTo, it->en);
      ConstructPath(paths, path, dest, max_length-1, threshold);
      path.RemoveLast();
    }
  }
}

void DSU::PrintLinkCP(FILE *output, bool fix)
{
  PrintLM(output, true);
  PrintSaddles(output, false);
}

void DSU::PrintLM(FILE *output, bool fix)
{
  if (fix) SortFix();
  fprintf(output, "      %s\n", seq);

  //printf("Local minima (%4d):\n", (int)LM.size());
  for (unsigned int i=0; i<LM.size(); i++) {
    fprintf(output, "%4d  %s %7.2f %8s\n", i+1, LM[i].str_ch, LM[i].energy/100.0, LMtype_string[LM[i].type]);
  }
}

void DSU::PrintBarr(FILE *output)
{
  fprintf(output, "     %s\n", seq);

  //get heights
  vector<vector<std::pair<int, int> > > res(number_lm);
  for (int i=0; i<number_lm; i++) {
    res[i] = HeightSearch(i, edgesV_l);
  }

  // get minimal height:
  vector<int> saddle_en(number_lm);
  vector<int> father(number_lm);
  for (int i=number_lm-1; i>=1; i--) {
    int minim = INT_MAX;
    int index;
    for (int j=i-1; j>=0; j--) {
      if (res[i][j].first <= minim) {
        minim = res[i][j].first;
        index = j;
      }
    }
    saddle_en[i] = minim;
    father[i] = index+1;
  }

  saddle_en[0] = LM[0].energy+0.1;
  father[0] = 0;


  for (int i=0; i<number_lm; i++) {
    fprintf(output, "%4d %s %6.2f %4d %6.2f\n", i+1, LM[i].str_ch, LM[i].energy/100.0, father[i], (saddle_en[i]-LM[i].energy)/100.0);
  }
}


void DSU::PrintSaddles(FILE *output, bool fix)
{
  if (fix) SortFix();
  //printf("Saddles (%4d):\n", (int)saddles.size());
  fprintf(output, "\n");
  // collect saddle info
  vector<set<int> > saddle_connLM (saddles.size());
  vector<set<int> > saddle_connSadd (saddles.size());
  for (set<edgeLS>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
    saddle_connLM[it->j].insert(it->i);
  }
  for (set<edgeSS>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
    saddle_connSadd[it->j].insert(it->i);
    saddle_connSadd[it->i].insert(it->j);
  }
  for (unsigned int i=0; i<saddles.size(); i++) {
    fprintf(output, "%4dS %s %7.2f %8s", i+1, saddles[i].str_ch, saddles[i].energy/100.0, SDtype_string[saddles[i].type]);
    for (set<int>::iterator it=saddle_connLM[i].begin(); it!=saddle_connLM[i].end(); it++) fprintf(output, " %4d", (*it)+1);
    for (set<int>::iterator it=saddle_connSadd[i].begin(); it!=saddle_connSadd[i].end(); it++) fprintf(output, " %3dS", (*it)+1);
    fprintf(output, "\n");
  }
}

bool DSU::SortFix()
{
  //return ;
  vector<int> mappingLM(LM.size());
  vector<int> mappingSD(saddles.size());
  {
    // create mapping (and sort LMs)
    vector<std::pair<RNAlocmin, int> > temp(LM.size());
    for (unsigned int i=0; i<LM.size(); i++) {
      temp[i] = make_pair(LM[i], i);
    }
    sort(temp.begin(), temp.end());
    for (unsigned int i=0; i<temp.size(); i++) {
      LM[i]=temp[i].first;
      mappingLM[temp[i].second]=i;
      vertex_l[temp[i].first]=i;
    }
  }
  {
    // create mapping (and sort saddles)
    vector<std::pair<RNAsaddle, int> > temp;
    for (unsigned int i=0; i<saddles.size(); i++) {
      temp.push_back(make_pair(saddles[i], i));
    }
    sort(temp.begin(), temp.end());
    for (unsigned int i=0; i<temp.size(); i++) {
      saddles[i]=temp[i].first;
      saddles[i].lm1 = mappingLM[saddles[i].lm1];
      saddles[i].lm2 = mappingLM[saddles[i].lm2];
      if (saddles[i].lm1 > saddles[i].lm2) swap(saddles[i].lm1, saddles[i].lm2);
      mappingSD[temp[i].second]=i;
    }
  }

  // now process all edges:
  set<edgeLL> edges_l_t;
  for (set<edgeLL>::iterator it=edges_l.begin(); it!=edges_l.end(); it++) {
    edges_l_t.insert(edgeLL(mappingLM[it->i], mappingLM[it->j], saddles[mappingSD[it->saddle]].energy, mappingSD[it->saddle]));
  }
  edges_l = edges_l_t;

  set<edgeLS> edges_ls_t;
  for (set<edgeLS>::iterator it=edges_ls.begin(); it!=edges_ls.end(); it++) {
    edges_ls_t.insert(edgeLS(mappingLM[it->i], mappingSD[it->j]));
  }
  edges_ls = edges_ls_t;

  set<edgeSS> edges_s_t;
  for (set<edgeSS>::iterator it=edges_s.begin(); it!=edges_s.end(); it++) {
    edges_s_t.insert(edgeSS(mappingSD[it->i], mappingSD[it->j]));
  }

  // process edgesV_l:
  vector<set<edgeLL> > edges_temp(edgesV_l.size());
  for (unsigned int i=0; i<edgesV_l.size(); i++) {
    for (set<edgeLL>::iterator it=edgesV_l[i].begin(); it!=edgesV_l[i].end(); it++) {
      edges_temp[i].insert(edgeLL(mappingLM[it->i], mappingLM[it->j], saddles[mappingSD[it->saddle]].energy, mappingSD[it->saddle]));
    }
  }
  for (unsigned int i=0; i<edgesV_l.size(); i++) {
    edgesV_l[mappingLM[i]]=edges_temp[i];
  }

  return false;
}
