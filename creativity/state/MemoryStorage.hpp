#pragma once
#include "creativity/state/Storage.hpp"
#include <vector>

namespace creativity { namespace state {

/** Class for in-memory storage (using an underlying std::vector).  All States are stored in memory
 * (via std::shared_ptr) for immediate retrieval.
 */
class MemoryStorage : public Storage {
    public:
        MemoryStorage() = delete;

        /// Creates an empty MemoryStorage object, using the given creativity settings reference.
        MemoryStorage(CreativitySettings &set);

        /** Creates a MemoryStorage object by coping the States and settings of the given Storage
         * object into new in-memory states.  Note that this copies std::shared_ptr<State> from the
         * given object; no deep copying is performed.
         */
        MemoryStorage(const Storage &copy);

        /** Returns a shared pointer to the stored State data.  Note that since the State object is
         * stored, this method returns a shared_pointer with `.use_count() > 1`.
         */
        virtual std::shared_ptr<const State> operator[](size_t i) const override;

        /// Returns the number of states currently stored.
        virtual size_t size() const override;

        /// Reserves the requested capacity in the underlying std::vector
        virtual void reserve(size_t capacity) override;

        /// Does nothing; settings are not stored outside the CreativitySettings reference.
        virtual void updateSettings() override {}

    protected:
        /// Adds a new State pointer to this storage container.
        virtual void push_back_(std::shared_ptr<const State> &&state) override;

    private:
        std::vector<std::shared_ptr<const State>> states_;
};

}}
