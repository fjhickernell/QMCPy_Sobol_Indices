#ifndef TOOLS_HH
#define TOOLS_HH

#include <vector>
#include "BooleanMatrix.hh"

std::vector<int> bin_index(int m){
	// This method computes the binary decomposition indices that are 1
	std::vector<int> m_bin;
	int index = 0;
	if (m != 0){
		while (m != 0){
			if (2*(m >> 1) != m) m_bin.push_back(index);
			m >>= 1;
			index++;
		}
	}
	return m_bin;
};

int t_value(std::vector< std::vector<bool> > Cj, std::vector< std::vector<bool> > Cd,int m){
// Gives the t-value between C_j and C_d when m is fixed
std::vector< std::vector<bool> > aux;
BooleanMatrix M;
int rd;
    for (int t = 0U; t<m; t++){
        if (t == m - 1U) return t;
        for (int rj = 1U; rj <= m - t; ++rj){
            rd = m - t - rj;
	    aux.clear();
	    for (unsigned k = 0U; k<rj+rd; k++){
		if (k<rj) aux.push_back(Cj[k]);
		else aux.push_back(Cd[k-rj]);
	    }
//		for (unsigned pp = 0U; pp<rj+rd; pp++){
//			for (unsigned pp2 = 0U; pp2<m; pp2++){
//				std::cout << aux[pp][pp2] << " ";
//			}
//			std::cout << std::endl;
//		}
//	    std::cout << std::endl;
	    M.set(aux,rj+rd,m);
//	    M.print();
//	    std::cout << M.getRank() << std::endl;
//	    std::cout << "-----------------------" << std::endl;
	    aux.clear();
            if (M.getRank() < rj+rd) break; // M is not full rank, we increase t
            else if (rj == m - t) return t;
        }
    }
};

unsigned power2(unsigned b){
	unsigned solution = 1U;
	for (unsigned k=0U; k<b; k++) solution *= 2U;
	return solution;
};


#endif // TOOLS_HH
