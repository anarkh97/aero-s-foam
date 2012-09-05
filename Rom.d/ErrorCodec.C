#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>

using namespace std;

void print_help();
void print_syntax();

int main (int argc, char *argv[]) {

  if (argc == 1) {
    print_help();
    return EXIT_SUCCESS;
  }

  if (argc == 2) {
    print_syntax();
    return EXIT_FAILURE;
  }

  //open input file streams
  ifstream truth_file (argv[1]);
  ifstream comp_file (argv[2]);

  string header_buffer;
  int num_nodes, length2, num_time_steps;
  double time1, time2;
  double a1, b1, c1, a2, b2, c2;
  double dummy1, dummy2, dummy3, sumx, sumy, sumz, sumx2, sumy2, sumz2;
  double cum_normx, cum_normy, cum_normz, normalize_factorx, normalize_factory, normalize_factorz;
  double relative_errorx, relative_errory, relative_errorz;

  // check to see of both files were succefully opened
  if(truth_file.is_open() && comp_file.is_open()) {

    //get header line and length of displacement vector 
    getline(truth_file,header_buffer);
    truth_file>>num_nodes;
    getline(comp_file,header_buffer);
    comp_file>>length2;

    // check to see if displacement vectors are same length
    if(num_nodes != length2){
      cout << "incompatible files" << endl;
      return EXIT_FAILURE;}

    //first: loop over all timesteps
    cum_normx = 0; cum_normy = 0; cum_normz = 0;
    num_time_steps = 0;
    while(truth_file>>time1) {
      num_time_steps += 1;
      comp_file>>time2;

//      cout << "time stamp 1 = " << time1 << endl;
//      cout << "time stamp 2 = " << time2 << endl;

      if(time1 != time2) {
        cout << "time stamps do not match" << endl;
        return EXIT_FAILURE;}

      //begin computation for L2-norm for timestep, 'num_time_step'
      sumx = 0; sumy = 0; sumz = 0; sumx2 = 0; sumy2 = 0; sumz2 = 0;
      // second: loop of nodes
      for(int counter=0;counter<num_nodes;counter++) {
        // third: read in all dofs
	  truth_file>>a1; truth_file>>b1; truth_file>>c1;
 	  comp_file>>a2; comp_file>>b2; comp_file>>c2;

	  sumx += pow((a1-a2),2);      
          sumy += pow((b1-b2),2);
          sumz += pow((c1-c2),2);

          sumx2 += pow(a1,2);
          sumy2 += pow(b1,2);
          sumz2 += pow(c1,2);
     }

     //sum 2-norm for timestep,"num_time_step" (sum(n=1->n_t)[|v^n - v^n_I|]
     cum_normx += pow(sumx,0.5);
     cum_normy += pow(sumy,0.5);
     cum_normz += pow(sumz,0.5);

     normalize_factorx += pow(sumx2,0.5);
     normalize_factory += pow(sumy2,0.5);
     normalize_factorz += pow(sumz2,0.5);


     if(!truth_file)
        break;
    }

    //divide cummulative 2-norm by number of time steps
     relative_errorx = cum_normx/(normalize_factorx);
     relative_errory = cum_normy/(normalize_factory);
     relative_errorz = cum_normz/(normalize_factorz);


//     cout << "cummulative 2-norm = " << cum_norm << endl;
//     cout << "number of time steps = " << num_time_steps << endl;
//     cout << "normalize factor = " << normalize_factor << endl;
     cout << "*** relative error: x ***  = " << relative_errorx*100 << "%"<< endl;
     cout << "*** relative error: y ***  = " << relative_errory*100 << "%"<< endl;
     cout << "*** relative error: z ***  = " << relative_errorz*100 << "%"<< endl;

  }
}

void print_syntax() {
  std::printf("Syntax: relerror truthfile comparisonfile\n");
}

void print_help() {
  std::printf("Relative error computation executable\n");
  print_syntax();;
}
