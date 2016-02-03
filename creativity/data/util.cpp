#include "creativity/data/util.hpp"

namespace creativity { namespace data {

std::string double_str(double d, unsigned precision) {
    double round_trip;
    for (unsigned prec : {precision-2, precision-1}) {
        std::stringstream ss;
        ss.precision(prec);
        ss << d;
        ss >> round_trip;
        if (round_trip == d) { ss << d; return ss.str(); }
    }
    std::stringstream ss;
    ss.precision(precision);
    ss << d;
    return ss.str();
}

double quantile(const std::vector<double> &vals, double prob) {
    return quantile(vals.begin(), vals.end(), prob);
}

double quantile(const Eigen::Ref<const Eigen::VectorXd> &vals, double prob) {
    if (prob < 0 or prob > 1) throw std::logic_error("Requested quantile probability is invalid");
    if (vals.size() == 0) return std::numeric_limits<double>::quiet_NaN();
    double index = prob * (vals.size()-1);
    unsigned below = std::floor(index), above = std::ceil(index);
    if (below == above) return vals[above];
    return (above - index) * vals[below] + (index - below) * vals[above];
}

double variance(const std::vector<double> &vals, double mean) {
    if (vals.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(mean)) {
        mean = 0;
        for (const auto &v : vals) mean += v;
        mean /= vals.size();
    }
    // Keep a sum of the demeaned values: mathematically they should equal 0, but numerical error
    // can creep in, and so this is essentially collecting the numerical error so that it can be
    // corrected for.
    double demeaned_squares = 0, demeaned_sum = 0;
    for (const auto &v : vals) {
        double d = v - mean;
        demeaned_squares += d*d;
        demeaned_sum += d;
    }
    return (demeaned_squares - demeaned_sum*demeaned_sum / vals.size()) / (vals.size()-1);
}

}}
