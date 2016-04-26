#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>

using namespace std;

void print_help();
void print_syntax();

int main (int argc, char *argv[]) {

  if (argc == 1) {
    print_help();
    return EXIT_SUCCESS;
  }

  if (!(argc == 3 || 4)) {
    print_syntax();
    return EXIT_FAILURE;
  }

  // open input file streams
  ifstream truth_file (argv[1]);
  ifstream comp_file (argv[2]);

  string header_buffer;
  int num_nodes, length2;
  double time1, time2, tFinal;
  double a1, a2;
  double sumx, sumx2;
  double cum_normx, normalize_factorx;
  double relative_errorx;
  double maxabs = 0, maxrel = 0;
  bool getTime1 = true, getTime2 = true;
  // check to see of both files were successfully opened
  if(truth_file.is_open() && comp_file.is_open()) {

    if( argc == 3 ) {
      std::cout << "calculate error up to time: ";
      std::cin >> tFinal;
      std::cout << std::endl;
    } else if (argc == 4) {
      std::string tFinalString(argv[3]);
 
      tFinal = atof(tFinalString.c_str()); 
    }

    // get header line and length of result vector 
    getline(truth_file, header_buffer);
    truth_file >> num_nodes;
    getline(comp_file, header_buffer);
    comp_file >> length2;

    // check to see if result vectors are same length
    if(num_nodes != length2) {
      cout << "incompatible files" << endl;
      return EXIT_FAILURE;
    }

    // initialize variables
    sumx = 0; sumx2 = 0;
    
    // begin Froebenius norm computation
    // first: loop over all timesteps
    int tcounter = 1;
    while(true) {

      if(getTime1) {
       truth_file >> time1;
       if(truth_file.eof() || time1 > tFinal) {
         break;
       }
     }

      if(getTime2) {
        comp_file >> time2;
        if(comp_file.eof() || time2 > tFinal) {
          break;
        }
      }

      printf("\r time stamp %d = %f",tcounter,time1);
      tcounter++;
      // second: loop over nodes
      for(int counter = 0; counter < num_nodes; counter++) {

        // if the timestamps are the same then read data from both files
        if(time1 == time2) {

          // third: read in all dofs
          truth_file >> a1;
          comp_file >> a2;
          getTime1 = true;
          getTime2 = true;
	  
          sumx += pow((a1-a2),2);      
          sumx2 += pow(a1,2);

          maxabs = max(maxabs,fabs(a1-a2));
          if(a1 != 0) maxrel = max(maxrel,fabs(a1-a2)/fabs(a1));
        }
        else {

          if(time1 < time2) {
            truth_file >> a1;
            if(counter == 0) {
              std::cout << "\nskipping time step " << time1 << " in truthfile" << std::endl;
              getTime1 = true;
              getTime2 = false;
            }
          }

          if(time1 > time2) {
            comp_file >> a2;
            if(counter == 0) {
              std::cout << "\nskipping time step " << time2 << " in comparisonfile" << std::endl;
              getTime1 = false;
              getTime2 = true;
            }
          }
        }
      }

      if(!truth_file)
        break;
    }

    // square root of differences
    cum_normx = pow(sumx,0.5);

    // square root of absolute
    normalize_factorx = pow(sumx2,0.5);

    relative_errorx = cum_normx/(normalize_factorx);

    cout << "*** absolute L2 error = " << cum_normx << endl;
    cout << "*** relative L2 error = " << relative_errorx*100 << "%" << endl;
    cout << "*** absolute L1 error = " << maxabs << endl;
    cout << "*** relative L1 error = " << maxrel*100 << "%" << endl;

  }
}

void print_syntax() {
  std::printf("Syntax: relerr truthfile comparisonfile\n");
}

void print_help() {
  std::printf("Relative error computation executable\n");
  print_syntax();;
}
