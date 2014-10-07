#ifndef GULO_DISCRETERATEMATRIX_HPP
#define GULO_DISCRETERATEMATRIX_HPP
 
#include <iostream>
#include <cassert>
#include <Eigen/Dense>
 #include <unsupported/Eigen/MatrixFunctions>

 
namespace Gulo
{
    class DiscreteRateMatrix
    {
    public:
        
        friend std::ostream& operator<< (std::ostream& , const DiscreteRateMatrix&);
    
        DiscreteRateMatrix(int numStates, double rate = 1) :
            n(numStates),
            rates(n,n),
            eval(n,n),
            evec(n,n),
            dirty(true)
        {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    evec(i,j) = 0;
                    eval(i,j) = 0;
                    if (i == j) {
                        rates(i,j) = -(rate*(n-1));
                    } else {
                        rates(i,j) = rate;
                    }
                }
            }
        }
    
        ~DiscreteRateMatrix() 
        {
        }
    
        int numStates() const 
        {
            return n;
        }
    
        void setRate(int i, int j, double r) 
        {
            assert(i != j); 
            if (rates(i,j) != r) {
                rates(i,j) = r;
                double diag = 0;
                for (int k = 0; k < n; ++k) {
                    if (i != k) {
                        diag -= rates(i,k);
                    }
                }
                rates(i,i) = diag;
                dirty = true;
            }
        }
        
        double rate(int i, int j) const
        {
            return rates(i,j);
        }
        
        Eigen::MatrixXd probabilityMatrix(double time) const
        {
          return (rates*time).exp();
        }
    
    private:
    
        int n;
        Eigen::MatrixXd rates;
        mutable Eigen::MatrixXd eval;
        mutable Eigen::MatrixXd evec;
        mutable bool dirty;
    };
    
    std::ostream& operator<< (std::ostream& stream, const DiscreteRateMatrix& matrix)
    {
        stream << matrix.rates;
        return stream;
    }
}

#endif // GULO_DISCRETERATESMATRIX_HPP