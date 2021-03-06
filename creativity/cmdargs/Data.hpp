#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <eris/types.hpp>
#include <limits>
#include <string>
#include <vector>

namespace boost { namespace program_options { class variables_map; } }

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-data arguments (for generating data from crstate files). */
class Data : public CmdArgs {
    public:
        /// Constructor for data arguments; takes no arguments
        Data();

        /// Struct for storing --verify-* options
        struct {
            /// If non-zero, this verifies that given crstate data files have the same piracy_begins parameter
            eris::eris_time_t piracy_begins = 0;
            /// If non-zero, verify the total number of simulation periods
            eris::eris_time_t periods = 0;
            /// If non-zero, verify the policy_begins of given simulations
            eris::eris_time_t policy_begins = 0;
        } verify;

        /// Struct for storing --skip-* options
        struct {
            /// Skip piracy data if true
            bool piracy = false;
            /// Skip policy phase data if true
            bool policy = false;
            /// Skip short-run effects
            bool short_run = false;
        } skip;

        /// If true, produce output in human-readable format; otherwise CSV output.
        bool human_readable = false;

        /// If true, suppress the CSV header (doesn't apply when `human_readable` is true).
        bool no_csv_header = false;

        /// If true, *only* output the CSV header (no files may be specified).
        bool only_csv_header = false;

        /// The number of periods to use to calculate average values
        unsigned int data_periods = 25;

        /// The precision of double values.  The default is full double precision.
        unsigned int double_precision = std::numeric_limits<double>::max_digits10;

        /// Number of threads to use
        unsigned int threads = 0;

        /// Number of files to preload from disk into memory
        unsigned int preload = 0;

        /// The input files to load data from
        std::vector<std::string> input;

        /// If true, decompress xz files into memory (note that decompression into memory happens
        /// with `preload > 0` regardless of this option).
        bool memory_xz = false;

        /** If memory_xz is false, decompress files into temporary files in this directory instead
         * of into the same directory as the input file.
         */
        std::string tmpdir;

        /// Overridden to add " CRSTATE [CRSTATE ...]"
        virtual std::string usage() const override;

        /// Overridden to add info about input files
        virtual std::string help() const override;

        /// Overridden to add " -- simulation data collector"
        virtual std::string versionSuffix() const override;

    protected:
        /// Adds data collector command-line options into the option descriptions
        virtual void addOptions() override;

        /** Overridden to change default precision to 8 if `--human-readable` is specified (and
         * `--precision` isn't).
         */
        virtual void postParse(boost::program_options::variables_map &vars) override;

};

}}
