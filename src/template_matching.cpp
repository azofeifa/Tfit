
#include "load.h"
#include "model.h"
#include <limits>
#include <iostream>
#include <algorithm>
#include "template_matching.h"
#include <fstream>
#include <random>
#include "omp.h"
#include <cmath>
#include "BIC.h"
#include "FDR.h"
using namespace std;

double nINF = -exp(1000);
double INF  = exp(1000);

vector<vector<double>> merge(vector<vector<double>> X, double window) {
   int N = X.size(), t = 0;
   vector<vector<double> > nX;
   while (t < N) {
      double start = X[t][0] - window, stop = X[t][1], n = 0;
      vector<double> S = {0.0, 0.0, 0.0};
      while (t < N and stop + window > X[t][0] - window ) {
         stop = X[t][1];
         S[0] += X[t][2], S[1] += X[t][3], S[2] += X[t][4];
         t += 1, n += 1;
      }
      vector<double> row = {start , stop + window, S[0] / n, S[1] / n, S[2] / n };
      nX.push_back(row);
   }
   return nX;
}


vector<vector<double>> bubble_sort3(vector<vector<double>> X) { //sort vector of vectors by second
   bool changed = true;
   if (X.size() < 2) {
      return X;
   }
   while (changed) {
      changed = false;
      for (int i = 0; i < X.size() - 1; i++  )  {
         if (X[i][2] < X[i + 1][2]) {
            vector<double> copy   = X[i];
            X[i]          = X[i + 1];
            X[i + 1]          = copy;
            changed = true;
         }
      }
   }
   return X;
}

int sample_centers(vector<double> centers, double p) {
   random_device rd;
   mt19937 mt(rd());
   default_random_engine generator;
   uniform_int_distribution<int> distribution(0, centers.size() - 1);
   int i   = distribution(mt);
   return i;//required in model.o (ugh...)
}

void BIC_template(segment * data,  double * BIC_values, double * densities, double * densities_r, double window,
                  double sigma, double lambda, double foot_print, double pi, double w) {
   double vl;
   int NN        = int(data->XN);
   int threads   = omp_get_max_threads();
   int counts    = NN / threads;

   #pragma omp parallel num_threads(threads)
   {
      int tid   = omp_get_thread_num();
      int start   = tid * counts;
      int stop  = (tid + 1) * (counts);
      if (tid + 1 == threads) {
         stop  = NN;
      }
      int j = start, k = start;
      double N_pos = 0, N_neg = 0;
      double total_density;
      for (int i = start; i < stop; i++) {
         while ((j < data->XN) and ((data->X[0][j] - data->X[0][i]) < -window)) {
            N_pos -= data->X[1][j];
            N_neg -= data->X[2][j];
            j++;
         }
         while ((k < data->XN) and ((data->X[0][k] - data->X[0][i]) < window)) {
            N_pos += data->X[1][k];
            N_neg += data->X[2][k];
            k++;
         }

         if (k < data->XN  and j < data->XN and k != j ) {
            total_density   = (N_pos / (data->X[0][k] - data->X[0][j])) + (N_neg / (data->X[0][k] - data->X[0][j]));
            densities[i]    = N_pos ;
            densities_r[i]  = N_neg ;

            double mu   = (data->X[0][k] + data->X[0][j]) / 2.;

            BIC_values[i]   = BIC3(data->X,  j,  k,  i, N_pos,  N_neg,
                                   sigma, lambda, foot_print, pi, w);
         } else {
            BIC_values[i]   = 0;
            densities[i]  = 0;
            densities_r[i]  = 0;
         }
      }
   }
}

bool check_hit(double a, double b, double c, double x, double y, double z) {
   if (a > x and b > y and c > z) {
      return true;
   }
   return false;
}

double run_global_template_matching(vector<segment*> segments,
                                    string out_dir,  params * P, slice_ratio SC) {

   double CTT                    = 5; //filters for low coverage regions

   double ns                     = stod(P->p["-ns"]);
   double window                 = stod(P->p["-pad"]) / ns;
   double sigma, lambda, foot_print, pi, w;
   double ct                     = stod(P->p["-bct"]);
   sigma   = stod(P->p["-sigma"]) / ns , lambda = ns / stod(P->p["-lambda"]);
   foot_print = stod(P->p["-foot_print"]) / ns , pi = stod(P->p["-pi"]), w = stod(P->p["-w"]);

   bool SCORES     = not P->p["-scores"].empty();

   ofstream FHW_scores;

   if (SCORES) {
      FHW_scores.open(P->p["-scores"]);
   }
   string annotation;
   int prev, prev_start, stop;
   int N;
   int start, center;

   for (int i = 0; i < segments.size(); i++) {
      double * BIC_values   = new double[int(segments[i]->XN)];
      double * densities    = new double[int(segments[i]->XN)];
      double * densities_r  = new double[int(segments[i]->XN)];

      double l    =  segments[i]->maxX - segments[i]->minX;
      double ef     = segments[i]->fN * ( 2 * (window * ns) * 0.05  / (l * ns ));
      double er     = segments[i]->rN * ( 2 * (window * ns) * 0.05 / (l * ns ));
      double stdf   = sqrt(ef * (1 - (  2 * (window * ns) * 0.05 / (l * ns )  ) )  );
      double stdr   = sqrt(er * (1 - (  2 * (window * ns) * 0.05 / (l * ns ) ) )  );
      BIC_template(segments[i],  BIC_values, densities, densities_r, window, sigma, lambda, foot_print, pi, w);
      double start = -1, rN = 0.0 , rF = 0.0, rR = 0.0, rB = 0.0;
      vector<vector<double>> HITS;
      for (int j = 1; j < segments[i]->XN - 1; j++) {
         bool HIT = check_hit(BIC_values[j], densities[j], 
            densities_r[j], SC.threshold, ef + CTT * stdf, er + CTT * stdr  );
         if (SCORES) {
            double vl   = BIC_values[j];
            if (std::isnan(double(vl)) or std::isinf(double(vl))) {
               vl    = 0;
            }
            int DENS    = densities[j] + densities_r[j] ;
            FHW_scores << segments[i]->chrom << "\t" << to_string(int(segments[i]->X[0][j - 1]*ns + segments[i]->start)) << "\t";
            FHW_scores << to_string(int(segments[i]->X[0][j]*ns + segments[i]->start )) << "\t" << to_string(vl) + "\t" + to_string(densities[j]) + "\t" + to_string(densities_r[j]) + "\t" +  to_string(int(HIT)) << endl;
         }
         if ( HIT ) {
            if (start < 0) {
               start = segments[i]->X[0][j - 1] * ns + segments[i]->start;
            }
            start += 1, rN += 1 , rF += densities[j], rR += densities_r[j], rB += log10( SC.pvalue(BIC_values[j]) + pow(10, -20)) ;
         }
         if (not HIT and start > 0 ) {
            vector<double> row = {start , segments[i]->X[0][j - 1]*ns + segments[i]->start, rB / rN , rF / rN, rR / rN  };
            HITS.push_back(row);
            start = -1, rN = 0.0 , rF = 0.0, rR = 0.0, rB = 0.0;
         }
      }
      for (int j = 0; j < HITS.size(); j++) {
         segments[i]->bidirectional_bounds.push_back(HITS[j]);
      }
      segments[i]->bidirectional_bounds   = merge(segments[i]->bidirectional_bounds, window * 0.5);
   }
   return 1.0;
}















