/* 
 * -----------------------------------------------------------------------------
 * -----                               main.cpp                            -----
 * -----                          REED-SOLOMON CODES                       -----
 * -----------------------------------------------------------------------------
 *
 * File Description:
 *   This is the simulation file file for `reedSolomon` encoder/decoder
 *
 * Assumptions:
 *   None
 *
 * References:
 *   - http://downloads.bbc.co.uk/rd/pubs/whp/whp-pdf-files/WHP031.pdf
 *
 * Revision History
 *   Jun 02, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

#include "reedSolomon.h"

const int num_data = 1000;          // number of different channel samples
const int num_trials_per_pt = 1000; // number of trials at each data point

////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////
/*
 * plotResults
 * Description:
 *   Plot a set of (x,y) value pairs on a graph
 *
 * Arguments: 
 *   - xData = array of x-axis values
 *   - yData = array of y-axis values
 *   - dataSize = size of both data arrays.
 *   - xLabel = x-axis label
 *   - yLabel = y-axis label
 *   - title = plot title
 *
 * Return: 
 *   None
 *
 * Assumptions:
 *   - both data arrays are he same size, dataSize
 *
 * Operation:
 *   - Ensure gnuplot exists on the machine. If it doesn't exit with a
 *     notice to user.
 *   - Call the `gnuplot` program with x,y pairs coming from corresponding
 *     elements in xData and yData
 *  
 * Revision History
 *   Jun 02, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
void plotResults(double* xData, double* yData, int dataSize, 
                 const char *xLabel, const char * yLabel, const char *title) 
{  
  FILE *gnuplotPipe,*tempDataFile;  
  const char *tempDataFileName;  
  double x,y;  
  int i;  
  tempDataFileName = "tempData";  
  gnuplotPipe = popen("gnuplot","w");  
  if (gnuplotPipe) {  
    // turn off legend
    fprintf(gnuplotPipe, "unset key;");
    // set title
    fprintf(gnuplotPipe, "set title \"%s\";", title);
    // set x-axis label
    fprintf(gnuplotPipe, "set xlabel \"%s\";", xLabel);
    // set y-axis label
    fprintf(gnuplotPipe, "set ylabel \"%s\";", yLabel);
    
    fprintf(gnuplotPipe,"plot \"%s\" with points\n",tempDataFileName);  
    fflush(gnuplotPipe);  
    tempDataFile = fopen(tempDataFileName,"w");  
    for (i=0; i <= dataSize; i++) {  
      x = xData[i];  
      y = yData[i];              
      fprintf(tempDataFile,"%lf %lf\n",x,y);          
    }          
    fclose(tempDataFile);          
    printf("press enter to continue...");          
    getchar();          
    remove(tempDataFileName);          
    fprintf(gnuplotPipe,"exit \n");      
  } else {          
    printf("gnuplot not found...");      
  }  
}


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
/*
 * main
 * Description:
 *   main loop
 *
 * Arguments: 
 *   none
 *
 * Return: 
 *   Program exit status
 *
 * Operation:
 *   - initialize rand()'s seed. I do this because I like consistency in results
 *
 *   - pick default m,t values or prompt user to provide these values.
 *  
 *   - For the purpose of this simulation, I determined the performance of the
 *     performance of the Reed-Solomon codes for two different choices of params
 *     + m=7, t=60
 *     + m-7, t=30
 *
 *   - error check whatever values of m,t are to be used
 *     + ensure m is less than number of bits in a system int. We want our
 *       Galois Field elements to always fit in an int.
 *     + ensure k (== n-2t == 2^m -1 -2t) > 0
 *
 *   - The simulation procedure is as follows:
 *     + I fix the number of data points to be used `num_data_pts` and then the
 *       Probability of channel error, `Ps`, values will be incremented by a
 *       delta of 1/num_data. So each data point will have a unique `Ps`
 *       with `Ps` values range from [0.0, 1/(`num_data`-1))
 *      
 *     + I fix the number of trials at each data point, `num_trials_per_pt`.
 *       Remember that a data point is associated with a specific `Ps` value. So
 *       each batch of trials at a data point will work with the same `Ps` value
 *      
 *     + Each trial has a random mesage generated, encoded and passed through a
 *       channel with error probability `Ps`. The resulting received message is
 *       decoded and checked for correctness (by comparing to the original
 *       encoded codeword). The number of failed decodings is counted and taken
 *       as a percentage of `num_trials_per_pt`. This gives the error rate for
 *       that Ps.
 *
 *       ++ The more trials we have at a data point, the more confident we can
 *          be about the error rate associated with that data point (or that 
 *          `Ps`, depending on how you read it).
 *
 *     + This data is plotted appropriately on the screen.
 *
 *     + From the two plots it is evident that when t is smaller, the error
 *       rate gets to 1 from a lower Ps than expected if t were larger. 
 *       ++ This makes sense as we would expect that when the error correcting 
 *          ability of the Reed-Solom code is decreased we would only be 
 *          increasing its Error rate at any given Ps.
 *  
 * Revision History
 *   Jun 02, 2011    Nnoduka Eruchalu    Initial Revision
 *   Mar 16, 2014    Nnoduka Eruchalu    Cleaned up comments
 */
int main(void)
{
  srand(time(0));            // initialize random seed 
  
  unsigned int m, t;
  int k;
  
  // pick default m,t values
  m = 7;                         // probably want values < 16 for top speed
  t= 60; // 30                   // remember n = 2^m-1, so pick t accordingly
  
  // prompt user to provide m and t values 
  /*  cout << "enter m: ";
      cin >> m;
      cout << "enter t: ";
      cin >> t;*/
  
  // crude error checking
  // check that m <= number of bits in an int
  if(m > sizeof(int)*8)
    {
      cout << "m (== " << m << ") has to be <= int bit count of "
           << sizeof(int)*8 << endl;
      exit(0);
    }
  
  // check that k > 0
  k = (pow(2,m) - 1 - 2*t);
  if(k <= 0)
    {
      cout << "k (== n-2t == 2^m -1 -2t) = " << k << "is negative!!" << endl;
      exit(0);
    }
  
  // delta between consecutive channel probability of error (Ps) values
  double delta = 1.0/num_data;
  // array holding Ps value at each data point
  double *Ps = new double[num_data];
  // array holding error rate at each data point
  double *Error_Rate = new double[num_data]; 
  
  // initial data point will have Probability of error == 0
  double Pss = 0.0;
  
  // loop through data points, being sure to increment Ps by delta on each run
  for(int i = 0; i < num_data; i++, Pss+=delta)
    {
      // keep users informed on progress status
      cout << i << endl;
      
      // perform a number of trials of generating 
      int num_errors = 0;
      for(int j=0; j< num_trials_per_pt; j++)
        {
          reedSolomon rs(m, t); // create the reed solomon object with m and t
          rs.gen_rand_msg();    // generate a random  message
          rs.encode();          // encode the given message
          rs.sim_channel(Pss);     // pass encoded message through channel
          rs.decode();          // now decoded received message
          //rs.print_params();    // print some good stuff
          bool correctly_decoded = rs.compare();  // check if decoder worked
            
          if(!correctly_decoded) num_errors++;
        } 
      
      // log Ps value and error rate for this data point
      Ps[i] = Pss;
      Error_Rate[i] = (double) num_errors/num_trials_per_pt;    
    }
  
  // remind user what parameters were used in simulation
  cout << "m: " << m << "\tt: " << t << endl;
  
  string title = "Error Rate vs Ps with m=";
  title += to_string(m);
  title += ", t=";
  title += to_string(t);
    
  // plot results
  plotResults(Ps, Error_Rate, num_data, "Ps", "Error Rate", title.c_str());
  
  delete [] Ps;                 // cleanup
  delete [] Error_Rate;
  return 1;
}
