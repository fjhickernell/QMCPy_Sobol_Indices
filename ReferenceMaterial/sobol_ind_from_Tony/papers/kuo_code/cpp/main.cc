#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "BooleanMatrix.hh"
#include "tools.hh"
#include <ctime>


struct generator{
	unsigned d; // Dimension generated
	unsigned s; // Primitive polynomial degree
	unsigned a; // Primitive polynomial coefficients
	std::vector <int> m; // Initial directional numbers
	std::vector< std::vector<bool> > C; // Generating matrix in Z2
};

int main(){
  std::clock_t timer; // initializing a clock type
  timer = std::clock(); // starting time of clock

  unsigned D = 19U; // Number of dimensions to work with
  unsigned digits = 24U; // Number of digits and 2^digits points
  std::vector<generator> geners(D);

  // Setting dimension 1 which is defined independently
  std::vector<bool> vector_zeros(digits,false);
  for (unsigned i=0U; i<digits; i++){
	vector_zeros[i] = true;
	for (unsigned k=0U; k<D; k++) geners[k].C.push_back(vector_zeros); // Setting the identity matrix C to all dimensions by default
	vector_zeros[i] = false;
  }
  geners[0].d = 1U; geners[0].s = 0U;
  
  // Reading Kuo and Joe data file joe-kuo-old.1111 new-joe-kuo-6.21201
  std::ifstream infile("../new-joe-kuo-6.21201",std::ios::in);
  if (!infile) {
    std::cout << "Input file containing direction numbers cannot be found!\n";
    //exit(1);
  }
  char buffer[1000];
  infile.getline(buffer,1000,'\n');
  std::vector<int> m_bin;
  int next_m;
  int s_m;
  for (unsigned j = 1U; j<D; j++) {
    // Read in parameters from file 
    infile >> geners[j].d >> geners[j].s >> geners[j].a;
    for (unsigned i=0U; i<geners[j].s; i++) geners[j].m.push_back(0U); // Preallocate m
    for (unsigned i=0U; i<geners[j].s; i++) infile >> geners[j].m[i]; // Get m from file
    // Build matrix C for each dimension j
    for (unsigned i=1U; i<digits; i++){
	if (i<geners[j].s){
	        m_bin.clear();
        	m_bin = bin_index(geners[j].m[i]);
		for (unsigned k=1U; k<m_bin.size(); k++) geners[j].C[i-m_bin[k]][i] = true;
	}
	else{
		m_bin.clear();
		m_bin = bin_index(geners[j].a);
		s_m = geners[j].m.size();
		next_m = (unsigned) (pow(2,geners[j].s)*geners[j].m[s_m-geners[j].s]) ^ geners[j].m[s_m-geners[j].s];
		for (unsigned k=0U; k<m_bin.size(); k++) next_m ^= (unsigned) pow(2,geners[j].s-m_bin[k]-1)*geners[j].m[s_m-geners[j].s+m_bin[k]+1];
		geners[j].m.push_back(next_m);
		m_bin.clear();
                m_bin = bin_index(geners[j].m[i]);
                for (unsigned k=1U; k<m_bin.size(); k++) geners[j].C[i-m_bin[k]][i] = true;
	}
    }
  }
   infile.close();

   // Computing t_values table
	for (unsigned dim = 1U; dim<D; dim++){
		for (unsigned j = 0U; j<dim; j++) std::cout << t_value(geners[j].C,geners[dim].C,digits) << " ";
	std::cout << std::endl;
	}

   // Final total elapsed time
  std::cout << (std::clock() - timer)/((double)CLOCKS_PER_SEC) << " seconds." << std::endl;
}
