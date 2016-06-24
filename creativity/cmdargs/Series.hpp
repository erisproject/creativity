#pragma once
#include "creativity/cmdargs/CmdArgs.hpp"
#include <eris/types.hpp>
#include <limits>
#include <string>
#include <vector>

namespace creativity { namespace cmdargs {

/** CmdArgs subclass for creativity-series arguments (for generating series values from crstate files). */
class Series : public CmdArgs {
    public:
        /// Constructor for series arguments; takes no arguments
        Series();

        /** comma-separated list of series to display
         */
        std::string series = "net_u,reader_spending,books_bought,books_pirated,books_written,book_quality,book_profit,book_author_effort,book_author_scale,book_p0";

        /// A plea for help about possible series values
        bool help_series = false;

        /** The number of periods.  If default (0), read from the first file.  Files must have at
         * least this number of periods.
         */
        eris::eris_time_t periods = 0;

        /// Whether files are allowed to contain more than `periods` periods.
        bool allow_unused_periods = false;

        /** The piracy begins period.  All simulations must have the same value (unless it is larger
         * than `periods`).
         *
         * The default (0) reads the value from the first file loaded.
         */
        eris::eris_time_t piracy_begins = 0;

        /** The public sharing begins period.  All simulations must have the same value (unless it
         * is larger than `periods`).
         *
         * The default (0) reads the value from the first file loaded.
         */
        eris::eris_time_t public_sharing_begins = 0;

        /// If true, just warn instead of aborting for files that can't be read or contain invalid period values.
        bool ignore_errors = false;

        /// The precision of double values.  The default is full double precision.
        unsigned int double_precision = std::numeric_limits<double>::max_digits10;

        /// Number of threads to use
        unsigned int threads = 0;

        /// The input files to load data from
        std::vector<std::string> input;

        /// If true, decompress xz files into memory
        bool memory_xz = false;

        /** If memory_xz is false, decompress files into temporary files in this directory instead
         * of into the same directory as the input file.
         */
        std::string tmpdir;

        /** The output directory; each variable VAR will be written to a "series-VAR.csv" file in
         * the directory.  The directory will be created if it doesn't yet exist.
         */
        std::string output_dir;

        /// Overridden to add " FILE [FILE ...]"
        virtual std::string usage() const override;

        /// Overridden to add info about input files
        virtual std::string help() const override;

        /// Overridden to add " -- simulation data series generator"
        virtual std::string versionSuffix() const override;

    protected:
        /// Adds series command-line options into the option descriptions
        virtual void addOptions() override;

        /** Overridden to make sure an output directory is given.
         */
        virtual void postParse(boost::program_options::variables_map &vars) override;
};

}}
