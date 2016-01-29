#pragma once
#include <memory>
#include <ostream>
#include <cairomm/context.h>
namespace creativity { namespace data { namespace graph {

/// Simple class holding red, green, blue, and alpha values.  Each should be a value in [0,1].
class RGBA final {
    public:
        RGBA() = delete;

        /** Creates a RGBA colour object with a red, green, blue, and (optional) alpha value.
         *
         * \param r the red channel value, from 0 to 1.  0 is off, 1 is full intensity.
         * \param g the red channel value, from 0 to 1.  0 is off, 1 is full intensity.
         * \param b the red channel value, from 0 to 1.  0 is off, 1 is full intensity.
         * \param a the alpha channel value, from 0 to 1.  0 is completely transparent, 1 is
         * completely opaque.  If omitted, defaults to 1 (fully opaque).
         * \throws std::out_of_range if any of the given values are outside the [0,1] range.
         */
        constexpr RGBA(const double &r, const double &g, const double &b, const double &a = 1.0) :
            red  {r < 0 or r > 1 ? throw std::out_of_range("") : r},
            green{g < 0 or g > 1 ? throw std::out_of_range("") : g},
            blue {b < 0 or b > 1 ? throw std::out_of_range("") : b},
            alpha{a < 0 or a > 1 ? throw std::out_of_range("") : a} {}

        /** Creates a RGBA object with a gray value and optional alpha value.  The gray intensity is
         * used for the three red, green, and blue values.
         *
         * \throws std::out_of_range if any of the given values are outside the [0,1] range.
         */
        explicit constexpr RGBA(const double &i, const double &a = 1.0) : RGBA(i,i,i,a) {}

        /// Copy assignment operator
        RGBA& operator=(const RGBA &copy) { r_(copy.red); g_(copy.green); b_(copy.blue); a_(copy.alpha); return *this; }

        /// Compares two RGBA objects: they are considered equal if all components are equal.
        constexpr bool operator==(const RGBA& other) const { return equalsRGB(other) and alpha == other.alpha; }

        /// Negation of ==
        constexpr bool operator!=(const RGBA& other) const { return not (*this == other); }

        /** Returns true if the red, green, and blue channels are equal.  The alpha values do not have
         * to match.
         */
        constexpr bool equalsRGB(const RGBA& other) const { return red == other.red and blue == other.blue and green == other.green; }

        /** Returns true if the object is opaque, that is, has alpha = 1. */
        constexpr bool opaque() const { return alpha == 1; }

        /** Returns true if the object is completely transparent, that is, has alpha = 0. */
        constexpr bool transparent() const { return alpha == 0; }

        /** Using an RGBA object as a bool returns true unless the object is fully transparent.  In
         * other words, `if (color)` is equivalent to `if (not color.transparent())`.
         */
        constexpr operator bool() const { return alpha > 0; }

        /** "Multiplies" two RGBA objects together via alpha compositing "over" filter.  The result
         * is the colour that would result if the left hand side argument was painted, then the
         * right-hand-side argument was painted over top.  If the RHS argument is opaque or LHS is
         * transparent, the result is simply a copy of the RHS object; if RHS transparent, the
         * result is a copy of the LHS object; otherwise the result is a new, combined color.  The
         * final opacity of the object will be between the larger opacity value of the two objects
         * and 1 (i.e. overlaying results in increasing opacity).
         */
        RGBA operator*(const RGBA& other) const { return RGBA(*this) *= other; }

        /** Multiplies this colour by another colour.  See RGBA::operator* for details. */
        RGBA& operator*=(const RGBA& other) {
            if (other.opaque() or transparent()) *this = other;
            else if (not other.transparent()) {
                double myweight = alpha - alpha*other.alpha;
                double newalpha = other.alpha + myweight;
                r_((myweight*red   + other.alpha*other.red  ) / newalpha);
                g_((myweight*green + other.alpha*other.green) / newalpha);
                b_((myweight*blue  + other.alpha*other.blue ) / newalpha);
                a_(newalpha);
            }
            return *this;
        }

        /// Applies this colour to the given Cairo::Context as the source colour
        void applyTo(Cairo::RefPtr<Cairo::Context> ctx) const { ctx->set_source_rgba(red, green, blue, alpha); }

        /** Creates the struct from an object with get_red(), get_green(), get_blue(), and get_alpha()
         * methods.  In particular, this is designed to be usable with Gdk::RGBA (without actually having
         * to add a dependency on gdkmm).
         */
        template <class T, typename = typename std::enable_if<
            std::is_same<double, decltype(T::get_red())>::value and
            std::is_same<double, decltype(T::get_green())>::value and
            std::is_same<double, decltype(T::get_blue())>::value and
            std::is_same<double, decltype(T::get_alpha())>::value>::type>
        RGBA(const T& rgba) : RGBA(rgba.get_red(), rgba.get_blue(), rgba.get_green(), rgba.get_alpha()) {}

        const double red, ///< The red channel value
                     green, ///< The green channel value
                     blue, ///< The blue channel value
                     alpha; ///< The alpha channel value (1 == opaque, 0 == transparent)

        /** Sents this RGBA to an output stream in the form "RGBA(r,g,b,a)" where r, g, b, and a are
         * replaced with their [0,1] values.
         */
        friend std::ostream& operator<<(std::ostream &os, const RGBA &colour) {
            os << "RGBA(" << colour.red << "," << colour.green << "," << colour.blue << "," << colour.alpha << ")";
            return os;
        }

    private:
        // Private access to set new values.  Values outside [0,1] are truncated.
        void r_(double r) { const_cast<double&>(red)   = r < 0 ? 0 : r > 1 ? 1 : r; }
        void g_(double g) { const_cast<double&>(green) = g < 0 ? 0 : g > 1 ? 1 : g; }
        void b_(double b) { const_cast<double&>(blue)  = b < 0 ? 0 : b > 1 ? 1 : b; }
        void a_(double a) { const_cast<double&>(alpha) = a < 0 ? 0 : a > 1 ? 1 : a; }
};

constexpr RGBA
    Transparent{0, 0}, ///< Pre-declared Transparent colour; the colour is actually black, but with alpha set to 0
    Black{0}, ///< Pre-declared black colour, with full opacity
    White{1}; ///< Pre-declared white colour, with full opacity

/** Style class for a drawn line. */
struct LineStyle {
    /** Constructs a line style directly from an Gdk::RGBA object. */
    LineStyle(RGBA rgba, double thickness = 1.0, double length = std::numeric_limits<double>::quiet_NaN())
        : colour{std::move(rgba)}, thickness{std::move(thickness)}, length{std::move(length)}
    {}

    /** Applies the line style to the given Cairo::Context */
    void applyTo(Cairo::RefPtr<Cairo::Context> ctx) const { colour.applyTo(ctx); ctx->set_line_width(thickness); }
    /// The colour of the line
    RGBA colour;
    /// The thickness of the line
    double thickness;
    /** Length of the line, if applicable.  This is used for things like tick marks where
     * `thickness` is the line width and this is the line length.  In other cases, this value is
     * unused.
     */
    double length;
};

/** Style class for a filled region. */
class FillStyle {
    public:
        /** Constructs a FillStyle using just a fill colour: the borders will be omitted (i.e.
         * transparent and 0 width). */
        FillStyle(RGBA fill) : FillStyle(std::move(fill), LineStyle(Transparent, 0)) {};
        /** Constructs a FillStyle using a fill colour and line style for the borders. */
        FillStyle(RGBA fill, LineStyle border) : fill_colour{std::move(fill)}, border{std::move(border)} {}
        /// Constructs a FillStyle using just a LineStyle; the fill will be omitted (i.e. transparent).
        FillStyle(LineStyle border) : FillStyle(Transparent, border) {}
        /// Move constructor
        FillStyle(FillStyle&&) = default;
        /// Copy constructor
        FillStyle(const FillStyle&) = default;
        /// Default move assignment
        FillStyle& operator=(FillStyle&&) = default;
        /// Default copy assignment
        FillStyle& operator=(const FillStyle&) = default;

        /// The fill colour for the region
        RGBA fill_colour;

        /// The border colour for the region
        LineStyle border;
};

/** Style class for a filled, rectangular region.  In addition to the fill colour and border of a
 * FillStyle, this also allows any of the 4 borders to be disabled.
 */
class RectangleStyle : public FillStyle {
    public:
        /// Inherit constructors from FillStyle
        using FillStyle::FillStyle;

        /// Default move constructor
        RectangleStyle(RectangleStyle&&) = default;
        /// Default copy constructor
        RectangleStyle(const RectangleStyle&) = default;
        /// Default move assignment
        RectangleStyle& operator=(RectangleStyle&&) = default;
        /// Default copy assignment
        RectangleStyle& operator=(const RectangleStyle&) = default;

        /// Initialize using a FillStyle
        RectangleStyle(FillStyle f) : FillStyle(std::move(f)) {}

        ///@{
        /// True if the named border should be drawn, false if not. Only applies when drawing rectangles.
        bool border_top = true, border_right = true, border_bottom = true, border_left = true;
        ///@}

        /// Returns the border width if the top border is to be drawn, 0 otherwise.
        double borderT() const { return border_top ? border.thickness : 0; }
        /// Returns the border width if the right border is to be drawn, 0 otherwise.
        double borderR() const { return border_right ? border.thickness : 0; }
        /// Returns the border width if the bottom border is to be drawn, 0 otherwise.
        double borderB() const { return border_bottom ? border.thickness : 0; }
        /// Returns the border width if the left border is to be drawn, 0 otherwise.
        double borderL() const { return border_left ? border.thickness : 0; }
        /// Returns the total drawn horizontal border width, i.e. borderR() + borderL()
        double borderLR() const { return borderL() + borderR(); }
        /// Returns the total drawn vertical border width, i.e. borderT() + borderB()
        double borderTB() const { return borderT() + borderB(); }
};

}}}
