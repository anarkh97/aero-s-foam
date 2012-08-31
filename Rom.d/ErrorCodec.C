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
  double dummy1, dummy2, dummy3, sum, sum2, cum_norm, normalize_factor;
  double relative_error;

  // check to see of both files were succefully opened
  if(truth_file.is_open() && comp_file.is_open()) {

    //get header line and length of displacement vector 
    getline(truth_file,header_buffer);
    truth_file>>num_nodes;
    getline(comp_file,header_buffer);
    comp_file>>length2;

    //TODO: figure out how to allocate for number of time steps
    double matrix1[num_nodes][3][500];

    //outpout displacement vector lengths to screen
    //cout << num_nodes << endl;
    //cout << length2 << endl;

    // check to see if displacement vectors are same length
    if(num_nodes != length2){
      cout << "incompatible files" << endl;
      return EXIT_FAILURE;}

    //first: loop over all timesteps
    cum_norm = 0; num_time_steps = 0;
    while(truth_file>>time1) {
      num_time_steps += 1;
      comp_file>>time2;

//      cout << "time stamp 1 = " << time1 << endl;
//      cout << "time stamp 2 = " << time2 << endl;

      if(time1 != time2) {
        cout << "time stamps do not match" << endl;
        return EXIT_FAILURE;}

      //begin computation for L2-norm for timestep, 'num_time_step'
      sum = 0; sum2 = 0;
      // second: loop of nodes
      for(int counter=0;counter<num_nodes;counter++) {
        // third: loop over dofs
	for(int i=0;i<3;i++) {
	  truth_file>>a1;comp_file>>a2;
          matrix1[counter][i][num_time_steps-1] = a1;
	  sum += pow((a1-a2),2);      
          sum2 += pow(a1,2);
	}
     }

     //sum 2-norm for timestep,"num_time_step" (sum(n=1->n_t)[|v^n - v^n_I|]
     cum_norm += pow(sum,0.5);
     normalize_factor += pow(sum2,0.5);

     if(!truth_file)
        break;
    }

    /*/compute max(for i,j elements of 1,..,nt) ||v[i]_I - v[j]_I||_2
    normalize_factor = 0;
    for(int i = 0;i<num_time_steps;i++) {
      for(int j = i+1;j<num_time_steps;j++) {
        dummy3 = 0;
        for(int counter=0;counter<num_nodes;counter++){
	  for(int dim=0;dim<3;dim++){
           dummy3 += pow((matrix1[counter][dim][i] - matrix1[counter][dim][j]),2);
       	  }
        }
        dummy3 = pow(dummy3,0.5);
        if(dummy3 > normalize_factor){
        normalize_factor = dummy3;
        }
      }
    }*/

    //divide cummulative 2-norm by number of time steps
    //relative_error = cum_norm/(num_time_steps*normalize_factor);
     relative_error = cum_norm/(normalize_factor);

     cout << "cummulative 2-norm = " << cum_norm << endl;
     cout << "number of time steps = " << num_time_steps << endl;
     cout << "normalize factor = " << normalize_factor << endl;
     cout << "*** relative error ***  = " << relative_error*100 << "%"<< endl;

  }

}

void print_syntax() {
  std::printf("Syntax: relerror truthfile comparisonfile\n");
}

void print_help() {
  std::printf("Relative error computation executable\n");
  print_syntax();;
}
