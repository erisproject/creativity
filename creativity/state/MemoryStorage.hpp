#pragma once
#include "creativity/state/StorageBackend.hpp"
#include <eris/types.hpp>
#include <cstddef>
#include <memory>
#include <vector>

namespace creativity { struct CreativitySettings; }

namespace creativity { namespace state {

class Storage;

/** Class for in-memory storage (using an underlying std::vector).  All States are stored in memory
 * (via std::shared_ptr) for immediate retrieval.
 */
class MemoryStorage final : public StorageBackend {
    public:
        /** Creates a new, blank MemoryStorage object. */
        MemoryStorage() = default;

        /** Creates a MemoryStorage object by coping the States and settings of the given Storage
         * object into new in-memory states.  Note that this copies std::shared_ptr<State> from the
         * given object; no deep copying is performed.
         */
        MemoryStorage(const Storage &copy);

        /// Returns the number of states currently stored.
        virtual size_t size() const override;

        /// Reserves the requested capacity in the underlying std::vector
        virtual void reserve(size_t capacity) override;

        /// Does nothing; settings are not stored by this class.
        virtual void writeSettings(const CreativitySettings&) override;

        /// Does nothing; settings are not stored by this class.
        virtual void readSettings(CreativitySettings&) const override;

        /// Returns the given index.
        virtual std::shared_ptr<const State> load(eris::eris_time_t t) const override;

        /// Adds the given state to the stored states_ directly.
        virtual void enqueue(std::shared_ptr<const State> &&s) override;

    private:
        std::vector<std::shared_ptr<const State>> states_;
};

}}
