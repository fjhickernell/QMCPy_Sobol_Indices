#ifndef BOOLEANMATRIX_HH
#define BOOLEANMATRIX_HH

#include <vector>

class BooleanMatrix{
    /*This class contains a method gaussElimination()
     *that diagonalizes mat and returns the rank in Z2
     */

    std::vector< std::vector<bool> > mat; //boolean matrix
    int n, m;           //size of matrix nxm
    int rank;           //rank of the matrix

    public:

    /*Constructor
     * Required Parameters:
     * M ==> boolean matrix
     * n ==> number of rows
     * m ==> number of columns
     */

    BooleanMatrix(){
    this -> n = 0;
    this -> m = 0;
    }

    BooleanMatrix(std::vector< std::vector<bool> > M, int n, int m){
        this -> n = n;
        this -> m = m;
	this -> mat = M;
        gaussElimination();
    }
/*    template <size_t size_m>
    BooleanMatrix(bool M[][size_m], int n, int m){
        this -> n = n;
        this -> m = m;
        for (int i = 0; i < n; i++){
            std::vector<bool> row(m);
            for (int j = 0; j < m; j++) row[j] = M[i][j];
            mat.push_back(row);         
        }
        gaussElimination();
    }
*/

   /* Resets any BooleanMatrix to new values */	   
    void set(std::vector< std::vector<bool> > M, int n, int m){
        this -> n = n;
        this -> m = m;
        this -> mat = M;
        gaussElimination();
    }

    /* Does Gauss Elimination with partial pivoting on the matrix */
    /* It diagonalizes matrix mat and computes the rank */
    /* We work on n>=m. Otherwise, we work with the transpose */
     bool gaussElimination(){
	if (n >= m){
		rank = m;
	        for (int i = 0; i < n; i++){
        	    if (!mat[i][i]){
                	int j;
                	for (j = i+1; j < n && !mat[j][i]; j++);
                	if (j < n){
                    		for (int k = i; k < m; k++){
                       			bool t = mat[i][k];
                        		mat[i][k] = mat[j][k];
                        		mat[j][k] = t;
                   	 	}
            		}
		    }
            	    for (int j = i+1; j < n; j++){
                	if (mat[j][i]){
                    	for (int k = i; k < m; k++)
                        	mat[j][k] = mat[j][k] - mat[i][k];
                	}
            	    }
        	}
		for (int j = 0; j < m; j++){if (!mat[j][j]) rank--;};
    	}
	else{
                rank = n;
                for (int i = 0; i < m; i++){
                    if (!mat[i][i]){
                        int j;
                        for (j = i+1; j < m && !mat[i][j]; j++);
                        if (j < m){
                                for (int k = i; k < n; k++){
                                        bool t = mat[k][i];
                                        mat[k][i] = mat[k][j];
                                        mat[k][j] = t;
                                }
                        }
                    }
                    for (int j = i+1; j < m; j++){
                        if (mat[i][j]){
                        for (int k = i; k < n; k++)
                                mat[k][j] = mat[k][j] - mat[k][i];
                        }
                    }
                }
                for (int j = 0; j < n; j++){if (!mat[j][j]) rank--;};
	}
   }


    void print(){
            for (unsigned pp = 0U; pp < n; pp++){
                for (unsigned pp2 = 0U; pp2 < m; pp2++){
                        std::cout << (unsigned) mat[pp][pp2] << " ";
                }
	    std::cout << std::endl;
            }
    }

    /* Get the row rank of the boolean matrix
     * If you require the rank of the matrix, make sure that n >= m.
     * i.e. if n < m, call the constructor over the transpose.
     */
    int getRank(){
        return rank;
    }
};

#endif // BOLEANMATRIX_HH
