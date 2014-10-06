#ifndef GULO_NORMALDISTRIBUTION_HPP
#define GULO_NORMALDISTRIBUTION_HPP

#include <cassert>
#include <cmath>
#include <Eigen/Dense>

namespace Gulo
{
    template <int N>
    class NormalDistribution
    {
        static_assert(N>0, "Invalid dimensionality");
        
    public:
        
        NormalDistribution(double mean = 0, double variance = 1)
        {
            U.fill(mean);
            E.fill(0);
            for (size_t j = 0; j < N; ++j) {
                E(j,j) = variance;
            }
        }
        
        NormalDistribution(const Eigen::Matrix<double, N, 1> & u, const Eigen::Matrix<double, N, N> & e)
        : U(u),
          E(e)
        {
        }
        
        ~NormalDistribution()
        {
        }
        
        bool operator ==(const NormalDistribution<N> & other) const
        {
            return U == other.U && E == other.E;
        }
        
        bool operator != (const NormalDistribution<N> & other) const
        {
            return U != other.U || E != other.E;
        }
        
        size_t dimensions() const
        {
            return N;
        }
        
        double mean(size_t i) const
        {
            assert (i < N);
            return U(i);
        }
        
        double variance(size_t i) const
        {
            assert (i < N);
            return E(i,i);
        }
        
        double rho(size_t i, size_t j) const 
        {
            assert (i < N);
            assert (j < N);
            return (E(i,i) == 0 || E(j,j) == 0) ? 0 : E(i,j)/(sqrt(E(i,i))*sqrt(E(j,j)));
        }
        
        const Eigen::Matrix<double, N, 1> & meanVector() const
        {
            return U;
        }
        
        const Eigen::Matrix<double, N, N> & covarianceMatrix() const
        {
            return E;
        }
            
        double logPdf(const Eigen::Matrix<double, N, 1> & x) const
        {
            return -log(sqrt((M_PI * 2 * E).determinant())) - (0.5 * (x - U).transpose() * E.inverse() * (x-U))(0,0);
        }
        
        void setMean(size_t i, double u) 
        {
            assert (i < N);
            U(i,0) = u;
        }
        
        void setVariance(size_t i, double v)
        {
            assert (i < N);
            assert(v >= 0);
            double oldsd = sqrt(E(i,i));
            double newsd = sqrt(v);
            for (size_t j = 0; j < N; ++j) {
                if (j != i) {
                    double sdj = sqrt(E(j,j));
                    E(i,j) = E(j,i) = newsd * sdj * ((oldsd == 0 || sdj == 0) ? 0 : E(i,j)/(oldsd*sdj));
                } else {
                    E(i,j) = v;
                }
            }
        }
        
        void setRho(size_t i, size_t j, double r) 
        {
            assert (i < N);
            assert (j < N);
            assert(i != j);
            assert(r >= -1 && r <= 1);
            E(i,j) = E(j,i) = sqrt(E(i,i))*sqrt(E(j,j))*r;
        }
        
        NormalDistribution<N> operator +(const NormalDistribution<N> & other) const
        {
            NormalDistribution<N> d = *this;
            d += other;
            return d;
        }
        
        NormalDistribution<N> & operator +=(const NormalDistribution<N> & other) 
        {
            U += other.U;
            E += other.E;
            return *this;
        }
        
        NormalDistribution<N> operator *(const NormalDistribution<N> & other) const
        {
            auto ec = (E.inverse() + other.E.inverse()).inverse();
            auto uc = ec*(other.E.inverse()*other.U + E.inverse()*U);
            return NormalDistribution<N>(uc, ec);
        }
        
        NormalDistribution<N> & operator *=(const NormalDistribution<N> & other)
        {
            auto ec = (E.inverse() + other.E.inverse()).inverse();
            U = ec*(other.E.inverse()*other.U + E.inverse()*U);
            E = ec;
            return *this;
        }
        
        NormalDistribution<N> operator *(double d) const
        {
            return NormalDistribution<N>(U * d, E * d);
        }
        
        NormalDistribution<N> & operator *=(double d)
        {
            U *= d;
            E *= d;
            return *this;
        }
        
    private:
        
        Eigen::Matrix<double, N, 1> U;
        Eigen::Matrix<double, N, N> E;
    };
    
    template <>
    class NormalDistribution<-1>
    {
    public:
        
        NormalDistribution(size_t i = 1, double mean=0, double variance=1)
        : U(i),
          E(i,i)
        {
            assert(i >= 1);
            U.fill(mean);
            E.fill(0);
            for (size_t j = 0; j < i; ++j) {
                E(j,j) = variance;
            }
        }
        
        NormalDistribution(const Eigen::VectorXd & u, const Eigen::MatrixXd & e)
        : U(u),
          E(e)
        {
            assert(U.size() >= 1);
            assert(U.size() == E.rows());
            assert(U.size() == E.cols());
        }
        
        ~NormalDistribution()
        {
        }
        
        bool operator ==(const NormalDistribution<-1> & other) const
        {
            return U == other.U && E == other.E;
        }
        
        bool operator != (const NormalDistribution<-1> & other) const
        {
            return U != other.U || E != other.E;
        }
          
        size_t dimensions() const
        {
            return U.size();
        }
        
        double mean(size_t i) const
        {
            assert (i < dimensions());
            return U(i);
        }
        
        double variance(size_t i) const
        {
            assert (i < dimensions());
            return E(i,i);
        }
        
        double rho(size_t i, size_t j) const 
        {
            assert (i < dimensions());
            assert (j < dimensions());
            return (E(i,i) == 0 || E(j,j) == 0) ? 0 : E(i,j)/(sqrt(E(i,i))*sqrt(E(j,j)));
        }
        
        const Eigen::VectorXd & meanVector() const
        {
            return U;
        }
        
        const Eigen::MatrixXd & covarianceMatrix() const
        {
            return E;
        }
        
        double logPdf(const Eigen::VectorXd & x) const
        {
            assert(x.size() == U.size());
            return -log(sqrt((M_PI * 2 * E).determinant())) - (0.5 * (x - U).transpose() * E.inverse() * (x-U))(0,0);
        }
        
        void setMean(size_t i, double u) 
        {
            assert (i < dimensions());
            U(i) = u;
        }
        
        void setVariance(size_t i, double v)
        {
            assert (i < dimensions());
            assert(v >= 0);
            double oldsd = sqrt(E(i,i));
            double newsd = sqrt(v);
            for (size_t j = 0; j < dimensions(); ++j) {
                if (j != i) {
                    double sdj = sqrt(E(j,j));
                    E(i,j) = E(j,i) = newsd * sdj * ((oldsd == 0 || sdj == 0) ? 0 : E(i,j)/(oldsd*sdj));
                } else {
                    E(i,j) = v;
                }
            }
        }
        
        void setRho(size_t i, size_t j, double r) 
        {
            assert (i < dimensions());
            assert (j < dimensions());
            assert(i != j);
            assert(r >= -1 && r <= 1);
            E(i,j) = E(j,i) = sqrt(E(i,i))*sqrt(E(j,j))*r;
        }
        
        NormalDistribution<-1> operator +(const NormalDistribution<-1> & other) const
        {
            assert(dimensions() == other.dimensions());
            NormalDistribution<-1> d = *this;
            d += other;
            return d;
        }
        
        NormalDistribution<-1> & operator +=(const NormalDistribution<-1> & other) 
        {
            assert(dimensions() == other.dimensions());
            U += other.U;
            E += other.E;
            return *this;
        }
        
        NormalDistribution<-1> operator *(const NormalDistribution<-1> & other) const
        {
            assert(dimensions() == other.dimensions());
            auto ec = (E.inverse() + other.E.inverse()).inverse();
            auto uc = ec*(other.E.inverse()*other.U + E.inverse()*U);
            return NormalDistribution<-1>(uc, ec);
        }
        
        NormalDistribution<-1> & operator *=(const NormalDistribution<-1> & other)
        {
            assert(dimensions() == other.dimensions());
            auto ec = (E.inverse() + other.E.inverse()).inverse();
            U = ec*(other.E.inverse()*other.U + E.inverse()*U);
            E = ec;
            return *this;
        }
        
        NormalDistribution<-1> operator *(double d) const
        {
            return NormalDistribution<-1>(U * d, E * d);
        }
        
        NormalDistribution<-1> & operator *=(double d)
        {
            U *= d;
            E *= d;
            return *this;
        }
        
        NormalDistribution<-1> conditionalDistribution(const std::vector<size_t> & indices, const std::vector<double> & values, bool useOrigSize = false) const
        {
            assert(indices.size() > 0 && indices.size() <= dimensions());
            assert(values.size() == indices.size());
            auto dim = dimensions();
            
            if (indices.size() == dimensions()) {
                NormalDistribution<-1> result(dimensions());
                for (size_t i = 0; i < dimensions(); ++i) {
                    result.setMean(i, values[i]);
                    result.setVariance(i, 0);
                }
                return result;
            }
            
            // partition
            std::vector<size_t> reverseOrder(dim);
            std::vector<size_t> order;
            order.reserve(dim);
            for (size_t i = 0; i < dim; ++i) {
                if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
                    order.push_back(i);
                    reverseOrder[i] = order.size()-1;
                }
            }
            for (size_t i = 0; i < indices.size(); ++i) {
                assert(indices[i] < dim);
                order.push_back(indices[i]);
                reverseOrder[indices[i]] = order.size()-1;
            }
            Eigen::VectorXd UU(dim);
            Eigen::MatrixXd EE(dim, dim);
            for (size_t i = 0; i < order.size(); ++i) {
                UU(i) = U(order[i],0);
                for (size_t j = 0; j < order.size(); ++j) {
                    EE(i,j) = E(order[i],order[j]);
                }
            }   
            
            // conditional distribution of unfixed variables
            size_t q = dim - indices.size();
            size_t nmq = indices.size();
            Eigen::MatrixXd a(nmq, 1);
            for (size_t i = 0; i < values.size(); ++i) {
                a(i,0) = values[i];
            }
            auto E12E22inv = EE.block(0,q,q,nmq)*EE.block(q,q,nmq,nmq).inverse();
            auto nU = UU.block(0,0,q,1) + E12E22inv * (a - UU.block(q,0,nmq,1));    // mean vector of the q-dimensional multinormal distribution of unfixed values
            auto nE = EE.block(0,0,q,q) - E12E22inv * EE.block(q,0,nmq,q);          // covariance matrix of the q-dimensional multinormal distribution of unfixed values
            
            if (!useOrigSize) {
                return NormalDistribution<-1>(nU, nE);
            }
            
            // regenerate distribution of original size        
            NormalDistribution<-1> d(dim);
            size_t pos = 0;
            size_t pos2 = 0;
            for (size_t i = 0; i < dim; ++i) { // this is going to be really slow...
                if (reverseOrder[i] >= q) {
                    d.setMean(indices[pos], values[pos]);
                    d.setVariance(indices[pos], 0);
                    ++pos;
                } else {
                    d.setMean(i, nU(pos2, 0));
                    d.setVariance(i, nE(pos2, pos2));
                    ++pos2;
                }
            }
            for (size_t i = 0; i < q; ++i) {
                for (size_t j = i+1; j < q; ++j) {
                    d.setRho(order[i], order[j], (nE(i,i) == 0 || nE(j,j) == 0) ? 0 : nE(i,j)/(sqrt(nE(i,i))*sqrt(nE(j,j))));
                }
            }
            return d;
        }
    
    private:
        
        Eigen::VectorXd U;
        Eigen::MatrixXd E;
    };
    
    template <>
    class NormalDistribution<1>
    {
    public:
        
        NormalDistribution(double mean = 0, double variance = 1)
        {
            assert(variance >= 0);
            U(0,0) = mean;
            E(0,0) = variance;
        }
        
        NormalDistribution(const Eigen::Matrix<double,1,1> & u, const Eigen::Matrix<double,1,1> & e)
        : U(u),
          E(e)
        {
        }
        
        ~NormalDistribution()
        {
        }
        
        bool operator ==(const NormalDistribution<1> & other) const
        {
            return U == other.U && E == other.E;
        }
        
        bool operator != (const NormalDistribution<1> & other) const
        {
            return U != other.U || E != other.E;
        }
        
        size_t dimensions() const
        {
            return 1;
        }
        
        double mean(size_t i = 0) const
        {
            assert (i < 1);
            return U(0,0);
        }
        
        double variance(size_t i = 0) const
        {
            assert (i < 1);
            return E(0,0);
        }
        
        double rho(size_t i = 0, size_t j = 0) const 
        {
            assert (i < 1);
            assert (j < 1);
            return 0;
        }
        
        const Eigen::Matrix<double, 1, 1> & meanVector() const
        {
            return U;
        }
        
        const Eigen::Matrix<double, 1, 1> & covarianceMatrix() const
        {
            return E;
        }
        
        double logPdf(double x) const
        {
            return -pow(x-U(0,0),2)/(2*E(0,0)) - 0.91893853320467274178 - log(sqrt(E(0,0)));
        }
        
        double logPdf(const Eigen::Matrix<double, 1, 1> & x) const
        {
            return logPdf(x(0,0));
        }
        
        void setMean(double u)
        {
            U(0,0) = u;
        }
        
        void setMean(size_t i, double u) 
        {
            assert (i == 0);
            U(0,0) = u;
        }
        
        void setVariance(size_t i, double v)
        {
            assert(i == 0);
            assert(v >= 0);
            E(0,0) = v;
        }
        
        void setVariance(double v)
        {
            assert( v >= 0);
            E(0,0) = v;
        }
        
        void setRho(size_t i, size_t j, double r) 
        {
        }
        
        NormalDistribution<1> operator +(const NormalDistribution<1> & other) const
        {
            return NormalDistribution<1>(U + other.U, E + other.E);
        }
        
        NormalDistribution<1> & operator +=(const NormalDistribution<1> & other) 
        {
            U += other.U;
            E += other.E;
            return *this;
        }
        
        NormalDistribution<1> operator *(const NormalDistribution<1> & other) const
        {
             return NormalDistribution<1>((other.mean() * variance() + mean() * other.variance())/(variance() + other.variance())
                                        , (variance() * other.variance())/(variance() + other.variance()));
        }
        
        NormalDistribution<1> & operator *=(const NormalDistribution<1> & other)
        {
            U(0,0) = (other.mean() * variance() + mean() * other.variance())/(variance() + other.variance());
            E(0,0) = (variance() * other.variance())/(variance() + other.variance());
            return *this;
        }
        
        NormalDistribution<1> operator *(double d) const
        {
            return NormalDistribution<1>(U * d, E * d);
        }
        
        NormalDistribution<1> & operator *=(double d)
        {
            U *= d;
            E *= d;
            return *this;
        }
        
    private:
        
        Eigen::Matrix<double, 1, 1> U;
        Eigen::Matrix<double, 1, 1> E;
    };
}

#endif // GULO_NORMALDISTRIBUTION_HPP
