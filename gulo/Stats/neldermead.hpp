#ifndef GULO_NELDERMEAD_HPP
#define GULO_NELDERMEAD_HPP

#include <cmath>
#include <algorithm>
#include <functional>

// uses adaptive values for alpha, beta, gamma and delta from:
// Gao and Han (2010) Implementing the Nelder-Mead simplex algorithm with adaptive centreOfGravity()()

namespace Gulo
{
    class NelderMead
    {
        struct ScoredVector
        {
            ScoredVector(const std::vector<double> & vector, double funcval) 
            : values(vector),
            score(funcval)
            {
            }
            
            std::vector<double> values;
            double score;
        };
        
    public:
        
        NelderMead(const std::function<double(const std::vector<double> &)> & f, const std::vector<double> & guess, double prec = 0.0001)
        : dim(guess.size()),
        func(f),
        alpha(1),
        beta(1.0 + 2.0/dim),
        gamma(0.75 - 1.0/(2.0*dim)),
        delta(1.0 - 1.0/dim),
        precision(prec)
        {       
            insert(guess);
            for (size_t i = 0; i < dim; ++i) {
                std::vector<double> x = guess;
                for (size_t j = 0; j < dim; ++j) {
                    double hj = (guess[j] == 0) ? 0.00025 : 0.05;
                    double ej = (i == j) ? 1 : 0;
                    x[j] += hj * ej;
                }
                insert(x);
            }
        }
        
        NelderMead(const std::function<double(const std::vector<double> &)> & f, const std::vector<double> & guess, const std::vector<double> scale, double prec = 0.0001)
        : dim(guess.size()),
        func(f),
        alpha(1),
        beta(1.0 + 2.0/dim),
        gamma(0.75 - 1.0/(2.0*dim)),
        delta(1.0 - 1.0/dim),
        precision(prec)
        {
            insert(guess);
            for (size_t i = 0; i < dim; ++i) {
                std::vector<double> x = guess;
                for (size_t j = 0; j < dim; ++j) {
                    double hj = (guess[j] == 0) ? 0.00025 : 0.05;
                    double ej = (i == j) ? scale[i] : 0;
                    x[j] += hj * ej;
                }
                insert(x);
            }
        }
        
        double optimize()
        {
            size_t count = 0;
            while(!finished()) {
                ++count;

                std::sort(vectors.begin(), vectors.end(), [](const ScoredVector & first, const ScoredVector & second)->bool{return first.score < second.score;});
                auto x0 = centreOfGravity();
                 std::cout << func(x0) << " " << func(vectors[0].values) << " " << maxDistance() << std::endl;
//                 for (int i = 0; i < vectors[0].values.size(); ++i) {
//                     std::cout <<exp( x0[i]) << " ";
//                 }
//                 std::cout << std::endl;
                
                /////////////
                // reflect //
                /////////////
                std::vector<double> xr(dim, 0);
                for (size_t i = 0; i < dim; ++i) {
                     xr[i] = x0[i] + alpha * (x0[i] - vectors[dim].values[i]);
                }
                double fxr = func(xr);
                if (vectors[0].score <= fxr && fxr < vectors[dim-1].score) {
                    vectors[dim] = ScoredVector(xr, fxr);   
                    continue;
                }
                
                ////////////
                // expand //
                ////////////
                if (fxr < vectors[0].score) {
                    std::vector<double> xe(dim, 0);
                    for (size_t i = 0; i < dim; ++i) {
                        xe[i] = x0[i] + beta * (x0[i] - vectors[dim].values[i]);
                    }
                    double fxe = func(xe);
                    if (fxe < fxr) {
                        vectors[dim] = ScoredVector(xe, fxe);  
                        continue;
                    } else {
                        vectors[dim] = ScoredVector(xr, fxr);
                        continue;
                    }
                } 
                
                //////////////////////
                // outside contract //
                //////////////////////
                std::vector<double> xc(dim, 0);
                for (size_t i = 0; i < dim; ++i) {
                    xc[i] = x0[i] + gamma * (vectors[dim].values[i] - x0[i]);
                }
                double fxc = func(xc);
                if (fxc < vectors[dim].score) {
                    vectors[dim] = ScoredVector(xc, fxc);
                    continue;
                }
                
                /////////////////////
                // inside contract //
                /////////////////////
                for (size_t i = 0; i < dim; ++i) {
                    xc[i] = x0[i] - gamma * (vectors[dim].values[i] - x0[i]);
                }
                fxc = func(xc);
                if (fxc < vectors[dim].score) {
                    vectors[dim] = ScoredVector(xc, fxc);
                    continue;
                }
                
                ////////////
                // reduce //
                ////////////
                for (size_t i = 1; i <= dim; ++i) {
                    for (size_t j = 0; j < dim; ++j) {
                        vectors[i].values[j] = vectors[0].values[j] + delta * (vectors[i].values[j] - vectors[0].values[j]);
                    }
                }
                
                
            }
            
            std::sort(vectors.begin(), vectors.end(), [](const ScoredVector & first, const ScoredVector & second)->bool{return first.score < second.score;});
            return func(vectors[0].values);
        }
        
        std::vector<double> centreOfGravity() const
        {
            std::vector<double> cog(dim, 0);
            for (size_t i = 0; i < dim; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                    cog[j] += vectors[i].values[j]/dim;
                }
            }
            
            return cog;
        }
        
        std::vector<double> parameters() const
        {
            return vectors[0].values;
        }
        
    private:
        
        void insert(const std::vector<double> & values)
        {
            vectors.push_back(ScoredVector(values, func(values)));
        }
        
        double maxDistance()
        {
            double result = 0;
            for (size_t i=0; i<dim+1; ++i) {
                for (size_t j=i+1; j<dim+1; ++j) {
                    double distance = 0;
                    for (size_t a = 0; a < dim; ++a) {
                        distance += pow(vectors[i].values[a] - vectors[j].values[a], 2);
                    }
                    if (sqrt(distance) > result) {
                        result = sqrt(distance);
                    }
                }
            }
            return result;
        }

        bool finished() const
        {
            for (size_t i=0; i<dim+1; ++i) {
                for (size_t j=i+1; j<dim+1; ++j) {
                    double distance = 0;
                    for (size_t a = 0; a < dim; ++a) {
                        distance += pow(vectors[i].values[a] - vectors[j].values[a], 2);
                    }
                    if (sqrt(distance) > precision) {
                        return false;
                    }
                }
            }
            return true;
        }
        
        size_t dim;
        const std::function<double(const std::vector<double> &)> & func;
        double precision;
        std::vector<ScoredVector> vectors;
        double alpha;
        double beta;
        double gamma;
        double delta;
    };
}

#endif // GULO_NELDERMEAD_HPP