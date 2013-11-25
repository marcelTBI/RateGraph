#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "RateGraph.h"


extern "C" {
  #include "RateGraph_cmdline.h"
}

using namespace std;

int main(int argc, char **argv)
{
  clock_t time = clock();
  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "Argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  // print-all
  if (args_info.print_all_flag) {
    args_info.barr_file_given = 1; free(args_info.barr_file_arg); args_info.barr_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.barr_file_arg, "barrier.lm");
    args_info.energy_file_given = 1; free(args_info.energy_file_arg); args_info.energy_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.energy_file_arg, "enbarr.eb");
    args_info.dist_file_given = 1; free(args_info.dist_file_arg); args_info.dist_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.dist_file_arg, "bpdistances.dist");
    args_info.gdist_file_given = 1; free(args_info.gdist_file_arg); args_info.gdist_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.gdist_file_arg, "graphdistances.gdist");
    args_info.rates_file_given = 1; free(args_info.rates_file_arg); args_info.rates_file_arg = (char*) malloc(50*sizeof(char));  strcpy(args_info.rates_file_arg, "rates.rat");
  }

  // code
    // BHGbuilder
  DSU dsu(stdin);
  fprintf(stderr, "reading input took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

  dsu.PrintDot(args_info.dot_file_arg, args_info.dot_flag, args_info.graph_file_given, args_info.graph_file_arg, args_info.tree_visualise_flag, args_info.dot_energies_flag);
  fprintf(stderr, "printing dot file took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

  // Shur removal
  if (args_info.max_given || args_info.filter_given) {
    RateGraph rg(dsu, args_info.rates_temp_arg);
    if (args_info.filter_given) {
      rg.ReadFilter(args_info.filter_arg);
      args_info.ordering_arg[0] = 'F';
    }
    rg.ConstructQueue(args_info.ordering_arg[0], args_info.max_arg);
    fprintf(stderr, "creating graph took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    int x = rg.RemoveX(dsu.Size()-args_info.max_arg, args_info.fraction_arg);
    fprintf(stderr, "removal of %d lm took %.2f secs.\n", x, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
    x = rg.RemoveShur(dsu.Size()-args_info.max_arg-x, args_info.Shur_step_arg);
    fprintf(stderr, "removal of %d lm took %.2f secs.\n", x, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

    if (args_info.rates_file_given) {
      FILE *file = fopen(args_info.rates_file_arg, "w");
      if (file){
        rg.PrintRates(file);
        fprintf(stderr, "printing rates into %s took %.2f secs.\n", args_info.rates_file_arg, (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
        fclose(file);
      }
    }
  }

  // barriers-like output
  if (args_info.barr_file_given) {
    FILE *file_h;
    file_h = fopen(args_info.barr_file_arg, "w");
    dsu.PrintBarr(file_h);
    fclose(file_h);
    fprintf(stderr, "printing barrier output took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
  }


  //print optimal path
  for (int i=0; i<(int)args_info.get_path_given; i++) {
    int a, b;
    if (sscanf(args_info.get_path_arg[i], "%d=%d", &a, &b)!=2) {
      fprintf(stderr, "WARNING: wrong use of --get-path option (%s)\n", args_info.get_path_arg[i]);
    } else {
      if (a<=0 || b<=0) {
        fprintf(stderr, "WARNING: non-positive number in --get-path (%s)\n", args_info.get_path_arg[i]);
      } else {
        a--;
        b--;
        if (a>=dsu.Size() || b>=dsu.Size()) {
          fprintf(stderr, "WARNING: visualisation number(s) exceeds number of minima (%d) (%s)\n", dsu.Size(), args_info.get_path_arg[i]);
        } else dsu.GetPath(a, b, 1000);
      }
    }
  }

  // print energy matrix
  if (args_info.energy_file_given) {
    dsu.PrintMatrix(args_info.energy_file_arg, args_info.print_full_flag, 'E');
  }

  // print dist matrix
  if (args_info.dist_file_given) {
    dsu.PrintMatrix(args_info.dist_file_arg, args_info.print_full_flag, 'D');
  }

  // print graph distance matrix
  if (args_info.gdist_file_given) {
    dsu.PrintMatrix(args_info.gdist_file_arg, args_info.print_full_flag, 'G');
  }

  fprintf(stderr, "printing matrices took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();

  // print rates matrix
  /*if (args_info.rates_file_given) {
    for (int i=0; i<max(1, (int)args_info.rates_mode_given); i++) {
      char filename[strlen(args_info.rates_file_arg)+2];
      strcpy(filename, args_info.rates_file_arg);
      filename[strlen(args_info.rates_file_arg)]=args_info.rates_mode_arg[i][0];
      filename[strlen(args_info.rates_file_arg)+1]='\0';
      //fprintf(stderr, filename);
      dsu.PrintRates(filename, args_info.print_full_flag, args_info.rates_temp_arg, args_info.rates_mode_arg[i][0]);
    }
    fprintf(stderr, "printing rates(%d) took %.2f secs.\n", max(1, (int)args_info.rates_mode_given), (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
  }*/

  // visualisation
  for (unsigned int i=0; i<args_info.visualise_given; i++) {
    int a, b;
    if (sscanf(args_info.visualise_arg[i], "%d=%d", &a, &b)!=2) {
      fprintf(stderr, "WARNING: wrong use of visualisation option (%s)\n", args_info.visualise_arg[i]);
    } else {
      if (a<=0 || b<=0) {
        fprintf(stderr, "WARNING: non-positive number in visualisation (%s)\n", args_info.visualise_arg[i]);
      } else {
        a--;
        b--;
        if (a>=dsu.Size() || b>=dsu.Size()) {
          fprintf(stderr, "WARNING: visualisation number(s) exceeds number of minima (%d) (%s)\n", dsu.Size(), args_info.visualise_arg[i]);
        } else dsu.VisPath(a, b, !args_info.vis_dist_flag, args_info.vis_length_arg, args_info.dot_flag, args_info.debug_flag);
      }
    }
  }

  // evaluation
  if (args_info.energy_heights_given) {
    FILE *file_h;
    file_h = fopen(args_info.energy_heights_arg, "w");
    bool full = true;
    dsu.EHeights(file_h, full);
    fclose(file_h);
  }

  // evaluation
  if (args_info.energy_rank_given) {
    FILE *file_h;
    file_h = fopen(args_info.energy_rank_arg, "w");
    dsu.ERank(file_h, false, true);
    fclose(file_h);
  }

  // evaluation
  if (args_info.energy_barrier_given) {
    FILE *file_h;
    file_h = fopen(args_info.energy_barrier_arg, "w");
    dsu.ERank(file_h, true);
    fclose(file_h);
  }

  // gaph analyze:
  if (args_info.analyze_graph_flag) {
    FILE *histo;
    histo = fopen("histogram.txt", "w");
    dsu.Histo(histo);
    fclose(histo);
  }

  fprintf(stderr, "RateGraph exitting succesfully!\n");

  cmdline_parser_free(&args_info);
  fprintf(stderr, "rest took %.2f secs.\n", (clock()-time)/(double)CLOCKS_PER_SEC); time = clock();
}
