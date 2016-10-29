#include <creativity/data/TreatmentFilter.hpp>
#include <memory>

namespace creativity { namespace data {

using namespace Eigen;

TreatmentFilter::TreatmentFilter(const Treatment &source,
        std::function<bool(const Properties&)> filter,
        std::function<bool(bool pre, bool piracy, bool policy, bool short_run)> stage_filter) {
    std::vector<unsigned> stages;
    stages.reserve(source.rowsPerSimulation());
    unsigned stage_i = 0;
    // Data is always added in this order: pre, piracy.SR, piracy.LR, policy.SR, policy.LR
    has_pre_ = false; // Everything else defaults to false already
    if (source.hasPre()) {
        if (stage_filter(true, false, false, false)) {
            has_pre_ = true;
            stages.push_back(stage_i);
        }
        stage_i++;
    }
    if (source.hasPiracySR()) {
        if (stage_filter(false, true, false, true)) {
            has_piracy_sr_ = true;
            stages.push_back(stage_i);
        }
        stage_i++;
    }
    if (source.hasPiracy()) {
        if (stage_filter(false, true, false, false)) {
            has_piracy_ = true;
            stages.push_back(stage_i);
        }
        stage_i++;
    }
    if (source.hasPolicySR()) {
        if (stage_filter(false, false, true, true)) {
            has_policy_sr_ = true;
            stages.push_back(stage_i);
        }
        stage_i++;
    }
    if (source.hasPolicy()) {
        if (stage_filter(false, false, true, false)) {
            has_policy_ = true;
            stages.push_back(stage_i);
        }
        stage_i++;
    }
    rows_per_sim_ = stages.size();
    if (rows_per_sim_ == 0) throw std::logic_error("TreatmentFilter: stage_filter admitted no stages!");

    const unsigned rowincr = std::min<unsigned>(1000, source.data().rows()/source.rowsPerSimulation()) * rows_per_sim_;
    data_ = MatrixXd(rowincr, source.data().cols());
    unsigned datapos = 0;

    for (unsigned i = 0; i < source.data().rows(); i += source.rowsPerSimulation()) {
        Properties props;
        props.source = source.sourceFiles()[i / source.rowsPerSimulation()];
        int j = 0;
        if (source.hasPre())      props.pre.reset(new StageProperties(source, i + j++));
        if (source.hasPiracySR()) props.piracy_SR.reset(new StageProperties(source, i + j++));
        if (source.hasPiracy())   props.piracy.reset(new StageProperties(source, i + j++));
        if (source.hasPolicySR()) props.policy_SR.reset(new StageProperties(source, i + j++));
        if (source.hasPolicy())   props.policy.reset(new StageProperties(source, i + j++));

        if (filter(props)) {
            if (datapos >= data_.rows()) data_.conservativeResize(data_.rows() + rowincr, NoChange);
            for (unsigned k : stages) {
                data_.row(datapos++) = source.data().row(i + k);
            }
            source_.push_back(props.source);
        }
    }

    data_.conservativeResize(datapos, NoChange);

    data_column_ = source.columns();
}

double TreatmentFilter::StageProperties::value(const std::string &field) const {
    return source_.data()(row_, source_.columns().at(field));
}

}}
