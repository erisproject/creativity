#pragma once
#include <cstdint>
#include <string>
#include <eris/serialize/serializer.hpp>

namespace creativity {

/** Class for interpreting and/or constructing encoded policy values.  The encoded value is stored
 * in the Creativity `settings.policy` parameter.  That value should be generated and interpreted
 * only through the abstractions provided by this class.
 */
class Policy final {
public:
    /// Default construction results in an all-disabled policy.
    constexpr Policy() : Policy(0) {};

    /** Constructs a Policy wrapper from an encoded policy value. This constructor intentionally
     * allows implicit conversion from a uint32_t previously obtained from a Policy object.
     */
    constexpr Policy(uint32_t policy_code) : code_{policy_code} {}

    /** Constructs a Policy wrapper from booleans for each policy type.
     *
     * \param public_sharing if true, constructs a policy with public sharing enabled
     * \param public_voting if true, constructs a policy with public sharing with voting
     * enabled
     * \param catch_pirates if true, constructs a policy with catching pirates enabled
     */
    constexpr Policy(
            bool public_sharing,
            bool public_voting,
            bool catch_pirates) : code_{
        public_sharing * POLICY_PUBLIC_SHARING +
        public_voting * POLICY_PUBLIC_VOTING +
        catch_pirates * POLICY_CATCH_PIRATES}
    {}

    /** Constructs a Policy wrapper from a comma-separated list of policies.  The following values
     * are permitted:
     *
     * - "public-sharing" (or "public")
     * - "public-voting" (or "voting" or "vote")
     * - "catch-pirates" (or "catch")
     * - "none" (or "") -- this value is silently ignored
     *
     * Any other value will result in a runtime_error exception being thrown.
     */
    explicit Policy(const std::string &p);

    /// Returns a constexpr Policy value for public sharing
    constexpr static Policy PublicSharing() { return POLICY_PUBLIC_SHARING; }
    /// Returns a constexpr Policy value for public voting
    constexpr static Policy PublicVoting()  { return POLICY_PUBLIC_VOTING; }
    /// Returns a constexpr Policy value for catch pirates
    constexpr static Policy CatchPirates()  { return POLICY_CATCH_PIRATES; }
    /// Returns a constexpr Policy value with all (known) policies enabled
    constexpr static Policy All() { return Policy(true, true, true); }

    /// Returns true if this policy has basic public sharing enabled
    constexpr bool publicSharing() const { return code_ & POLICY_PUBLIC_SHARING; }
    /// Returns true if this policy has public sharing with voting enabled
    constexpr bool publicVoting() const { return code_ & POLICY_PUBLIC_VOTING; }
    /// Returns true if this policy has catching pirates enabled
    constexpr bool catchPirates() const { return code_ & POLICY_CATCH_PIRATES; }

    /** Returns true if the policy has some unknown policy enabled (that is, some policy code that
     * isn't constructible from the known policies).
     */
    constexpr bool unknown() const { return code_ & ~(POLICY_PUBLIC_SHARING | POLICY_PUBLIC_VOTING | POLICY_CATCH_PIRATES); }

    /** Returns true if all (known) policies are enabled. */
    constexpr bool all() const { return publicSharing() && publicVoting() && catchPirates(); }

    /** Implicit conversion to bool: returns true if any known policy is enabled *and* there are no
     * unknown policy bits set in the policy code.
     */
    constexpr operator bool() const { return (publicSharing() || publicVoting() || catchPirates()) && !unknown(); }

    /** Implicit conversion to uint32_t policy code. */
    constexpr operator uint32_t() const { return code_; }

    /** Returns true if the two policy objects represent the same policy choice. */
    constexpr bool operator==(const Policy &other) const { return code_ == other.code_; }

    /** Negation of ==. */
    constexpr bool operator!=(const Policy &other) const { return code_ != other.code_; }

    /** Modifies this Policy value by enabling any policies enabled in the given Policy.  Currently
     * enabled policies remain enabled, even if not enabled in the target.
     *
     * Unknown policies in the target are ignored.
     */
    Policy& operator+=(const Policy &add);

    /** Modifies this Policy value by disabling any policies that are enabled in the given Policy.
     * Unknown policies are ignored.
     */
    Policy& operator-=(const Policy &remove);

    /** Returns a new Policy object with all policies enabled that are enabled in either argument.
     */
    constexpr Policy operator+(const Policy &plus) const {
        return Policy(
                publicSharing() || plus.publicSharing(),
                publicVoting()  || plus.publicVoting(),
                catchPirates()  || plus.catchPirates()
                );
    }

private:
    uint32_t code_;

    friend class eris::serialize::serializer<Policy>;
    friend class eris::serialize::serializer<const Policy>;

    // The bits for the known policies:
    constexpr static uint32_t
        POLICY_PUBLIC_SHARING = 1 << 0,
        POLICY_CATCH_PIRATES  = 1 << 1,
        POLICY_PUBLIC_VOTING  = 1 << 2;
};

}

namespace eris { namespace serialize {

/// Specialization of eris::serializer for a Policy reference
template <> class serializer<creativity::Policy> : public serializer<uint32_t> {
public:
    /// Constructs a serializer that serializes the policy's internal uint32_t policy code
    explicit serializer(creativity::Policy &p) : serializer<uint32_t>(p.code_) {}
};
/// Specialization of eris::serializer for a const Policy reference
template <> class serializer<const creativity::Policy> : public serializer<const uint32_t> {
public:
    /// Constructs a serializer that serializes the policy's internal uint32_t policy code
    explicit serializer(const creativity::Policy &p) : serializer<const uint32_t>(p.code_) {}
};

}}
