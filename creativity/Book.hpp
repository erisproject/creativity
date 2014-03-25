#pragma once
#include <eris/WrappedPositional.hpp>
#include <eris/Good.hpp>

using namespace eris;

namespace creativity {

class Reader; // forward-declaration

class Book : public Positional<Good::Discrete> {
    public:
        Book(Position p, SharedMember<Reader> r)
            : Positional<Good::Discrete>(p), author_(r) {}

        /// Returns the age of the book, in simulation periods.
        unsigned long age() {
            return simulation()->t() - created_;
        }

        virtual void added() override {
            created_ = simulation()->t();
        }

        SharedMember<Reader> author() {
            return simAgent<Reader>(author_);
        }

    private:
        unsigned long created_;
        eris_id_t author_;
};

}

