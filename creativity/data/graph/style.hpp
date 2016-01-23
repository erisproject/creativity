#pragma once
#include <memory>
#include <cairomm/context.h>
namespace creativity { namespace data { namespace graph {

/// Simple struct holding red, green, blue, and alpha values.  Each should be a value in [0,1].
struct RGBA {
    RGBA() = delete;
    /// Creates the struct with a red, green, blue, and (optional) alpha value.
    RGBA(double r, double g, double b, double a = 1.0) : red{r}, green{g}, blue{b}, alpha{a} {
        if (r<0 or r>1 or g<0 or g>1 or b<0 or b>1 or a<0 or a>1)
            throw std::invalid_argument("Invalid RGBA value: all values must be >= 0 and <= 1");
    }
    /// Compares two RGBA objects: they are considered equal if all components are equal.
    bool operator==(const RGBA& other) {
        return equalsRGB(other) and alpha == other.alpha;
    }
    /// Negation of ==
    bool operator!=(const RGBA& other) { return not (*this == other); }
    /** Returns true if the red, green, and blue channels are equal.  The alpha values do not have
     * to match.
     */
    bool equalsRGB(const RGBA& other) {
        return red == other.red and blue == other.blue and green == other.green;
    }
    /// Applies the colour to the given Cairo::Context
    void applyTo(Cairo::RefPtr<Cairo::Context> ctx) const { ctx->set_source_rgba(red, green, blue, alpha); }
    /** Creates the struct from an object with get_red(), get_green(), get_blue(), and get_alpha()
     * methods.  In particular, this is designed to be used with Gdk::RGBA (but we don't want to
     * have to depend on gdkmm, so this is left as a templated constructor).
     */
    template <class T, typename = typename std::enable_if<
        std::is_same<double, decltype(T::get_red())>::value and
        std::is_same<double, decltype(T::get_green())>::value and
        std::is_same<double, decltype(T::get_blue())>::value and
        std::is_same<double, decltype(T::get_alpha())>::value>::type>
    RGBA(const T& rgba) : RGBA(rgba.get_red(), rgba.get_blue(), rgba.get_green(), rgba.get_alpha()) {}

    double red, ///< The red channel value
           green, ///< The green channel value
           blue, ///< The blue channel value
           alpha; ///< The alpha channel value (1 == opaque, 0 == transparent)
};

extern const RGBA
    Transparent, ///< Pre-declared Transparent colour; the colour is actually black, but with alpha set to 0
    Black, ///< Pre-declared black colour, with full opacity
    White; ///< Pre-declared white colour, with full opacity

/** Style class for a drawn line. */
struct LineStyle {
    /** Constructs a line style directly from an Gdk::RGBA object. */
    LineStyle(RGBA rgba, double thickness = 1.0)
        : colour{std::move(rgba)}, thickness{std::move(thickness)}
    {}

    /** Applies the line style to the given Cairo::Context */
    void applyTo(Cairo::RefPtr<Cairo::Context> ctx) const { colour.applyTo(ctx); ctx->set_line_width(thickness); }
    /// The colour of the line
    RGBA colour;
    /// The thickness of the line
    double thickness;
};

/** Style class for a filled region. */
struct FillStyle {
    /** Constructs a FillStyle using just a fill colour: the borders will be omitted (i.e.
     * transparent and 0 width). */
    FillStyle(RGBA fill) : FillStyle(std::move(fill), LineStyle(Transparent, 0)) {};
    /** Constructs a FillStyle using a fill colour and line style for the borders. */
    FillStyle(RGBA fill, LineStyle border) : fill_colour{std::move(fill)}, border{std::move(border)} {}
    /// Constructs a FillStyle using just a LineStyle; the fill will be omitted (i.e. transparent).
    FillStyle(LineStyle border) : FillStyle(Transparent, border) {}
    /// The fill colour for the region
    RGBA fill_colour;

    /// The border colour for the region
    LineStyle border;
};

}}}
