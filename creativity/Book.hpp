#include <eris/WrappedPositional.hpp>
#include <eris/Good.hpp>

namespace creativity {

class BookBase : public eris::Good::Discrete {
    public:
        BookBase(double quality) : quality_(quality) {}

        double quality() const { return quality_; }
    private:
        double quality_;
};

typedef eris::WrappedPositional<BookBase> Book;

}

