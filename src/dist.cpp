#include <fmt/core.h>
#include "dist.h"
namespace Distribution
{
    TruncatedNormal::TruncatedNormal(double mu, double sigma, double lb, double ub):
        mu_(mu), sigma_(sigma), lb_(lb), ub_(ub)
    {
        auto alpha = (lb_ - mu_) / sigma_;
        auto beta = (ub_ - mu_) / sigma_;

        alpha_cdf_ = normal_01_cdf(alpha);
        beta_cdf_ = normal_01_cdf(beta);

        assert(lb_ < ub_ &&
                fmt::format("lower bound: {} must < upper bound: {}",
                    lb_, ub_).c_str());
        assert(lb_ < mu_ && mu_ < ub_ &&
                fmt::format("lower bound: {} must < mu: {}, and "
                    "mu: {} must > upper bound: {}",
                    lb_, mu_, mu_, ub_).c_str());
    }
    double TruncatedNormal::operator()()
    {
        auto u = uniform_dis_(gen_);
        auto xi_cdf = alpha_cdf_ + u * (beta_cdf_ - alpha_cdf_);
        auto xi = normal_01_cdf_inv(xi_cdf);
        
        auto x = mu_ + sigma_ * xi;
        return x;
    }
    double normal_01_cdf_inv( double p )
    {
        double const a[8] = {
            3.3871328727963666080,     1.3314166789178437745E+2,
            1.9715909503065514427E+3,  1.3731693765509461125E+4,
            4.5921953931549871457E+4,  6.7265770927008700853E+4,
            3.3430575583588128105E+4,  2.5090809287301226727E+3 };
        double const b[8] = {
            1.0,                       4.2313330701600911252E+1,
            6.8718700749205790830E+2,  5.3941960214247511077E+3,
            2.1213794301586595867E+4,  3.9307895800092710610E+4,
            2.8729085735721942674E+4,  5.2264952788528545610E+3 };
        double const c[8] = {
            1.42343711074968357734,     4.63033784615654529590,
            5.76949722146069140550,     3.64784832476320460504,
            1.27045825245236838258,     2.41780725177450611770E-1,
            2.27238449892691845833E-2,  7.74545014278341407640E-4 };
        double const const1 = 0.180625;
        double const const2 = 1.6;
        double const d[8] = {
            1.0,                        2.05319162663775882187,
            1.67638483018380384940,     6.89767334985100004550E-1,
            1.48103976427480074590E-1,  1.51986665636164571966E-2,
            5.47593808499534494600E-4,  1.05075007164441684324E-9 };
        double const e[8] = {
            6.65790464350110377720,     5.46378491116411436990,
            1.78482653991729133580,     2.96560571828504891230E-1,
            2.65321895265761230930E-2,  1.24266094738807843860E-3,
            2.71155556874348757815E-5,  2.01033439929228813265E-7 };
        double const f[8] = {
            1.0,                        5.99832206555887937690E-1,
            1.36929880922735805310E-1,  1.48753612908506148525E-2,
            7.86869131145613259100E-4,  1.84631831751005468180E-5,
            1.42151175831644588870E-7,  2.04426310338993978564E-15 };
        double q;
        double r;
        double split1 = 0.425;
        double split2 = 5.0;
        double value;

        if ( p <= 0.0 )
        {
            value = - HUGE_VAL;
            return value;
        }

        if ( 1.0 <= p )
        {
            value = HUGE_VAL;
            return value;
        }

        q = p - 0.5;

        if ( fabs ( q ) <= split1 )
        {
            r = const1 - q * q;
            value = q * r8poly_value_horner ( 7, a, r ) 
                      / r8poly_value_horner ( 7, b, r );
        }
        else
        {
            if ( q < 0.0 )
            {
                r = p;
            }
            else
            {
                r = 1.0 - p;
            }

            if ( r <= 0.0 )
            {
                value = HUGE_VAL;
            }
            else
            {
                r = std::sqrt ( - std::log ( r ) );

                if ( r <= split2 )
                {
                    r = r - const2;
                    value = r8poly_value_horner ( 7, c, r ) 
                      / r8poly_value_horner ( 7, d, r );
                }
                else
                {
                    r = r - split2;
                    value = r8poly_value_horner ( 7, e, r ) 
                       / r8poly_value_horner ( 7, f, r );
                }
            }

            if ( q < 0.0 )
            {
                value = - value;
            }

        }

        return value;
    }
    double r8poly_value_horner( int m, double const c[], double x )
    {
        int i;
        double value;

        value = c[m];

        for ( i = m - 1; 0 <= i; i-- )
        {
            value = value * x + c[i];
        }

        return value;
    }
}
