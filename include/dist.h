#ifndef _INCLUDE_UTILS_DIST_H
#define _INCLUDE_UTILS_DIST_H
#include <cmath>
#include <random>
#include <cassert>
namespace Distribution
{
    // modify from https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.html
    double normal_01_cdf_inv( double p );
    double r8poly_value_horner( int m, double const c[], double x );
    class TruncatedNormal
    {
        // ctor
        public:
            TruncatedNormal(double mu, double sigma, double lb, double ub);
        // public func
        public:
            double operator()();

        // internal func
        protected:
            inline double normal_01_cdf(double x)
            {
                return std::erfc( -x / std::sqrt(2) ) / 2;
            }

        // internal data member
        protected:
            double mu_;
            double sigma_;
            double lb_;
            double ub_;
            std::mt19937 gen_{std::random_device{}()};
            std::uniform_real_distribution<> uniform_dis_{0, 1};

            double alpha_cdf_;
            double beta_cdf_;
    };
}
#endif
