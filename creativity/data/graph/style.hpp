#include <gdkmm/rgba.h>
#include <memory>
namespace creativity { namespace data { namespace graph {

/// Pre-declared transparent colour; the colour is actually black, but has opacity set to 0.
extern const Gdk::RGBA transparent;

/** Style class for a drawn line. */
struct LineStyle {
    /** Constructs a line style directly from an Gdk::RGBA object. */
    LineStyle(Gdk::RGBA rgba, double thickness = 1.0)
        : colour{std::move(rgba)}, thickness{std::move(thickness)}
    {}

    /** Constructs a line style from a colour string; this can be a name such as "red"; a hex
     * sequence of 3, 6, 9, or 12 hex characters beginning with '#' (for example, `#f00`, `#ff0000`,
     * `#ffff00000000` all represent red); a string in the form `rgb(r,g,b)` where `r`, `g`, `b` can
     * be either 0--255 integer values or numeric percentages between 0% and 100%; or a string in
     * the form `rgba(r,g,b,a)` where the additional `a` is a value in [0, 1] indicating the opacity
     * level.
     *
     * \param colour the colour name (passed to Gdk::RGBA).
     */
    LineStyle(std::string colour, double thickness = 1.0) : LineStyle(Gdk::RGBA(colour), thickness) {}

    /// The colour of the line
    Gdk::RGBA colour;
    /// The thickness of the line
    double thickness;
};

/** Style class for a filled region. */
struct FillStyle {
    /** Constructs a FillStyle using just a fill colour: the borders will be omitted (i.e.
     * transparent). */
    FillStyle(Gdk::RGBA fill) : FillStyle(std::move(fill), transparent) {};
    /** Same as above, but takes a string to pass to the Gdk::RGBA constructor. */
    FillStyle(std::string fill) : FillStyle(Gdk::RGBA(fill)) {}
    /** Constructs a FillStyle using a fill colour and line style for the borders. */
    FillStyle(Gdk::RGBA fill, LineStyle border) : fill{std::move(fill)}, border{std::move(border)} {}
    /** Same as above, but takes a string to pass to the Gdk::RGBA constructor. */
    FillStyle(std::string fill, LineStyle border) : FillStyle(Gdk::RGBA(fill), std::move(border)) {}
    /// Constructs a FillStyle using just a LineStyle; the fill will be omitted (i.e. transparent).
    FillStyle(LineStyle border) : FillStyle(transparent, border) {}
    /// The fill colour for the region
    Gdk::RGBA fill;
    /// The border colour for the region
    LineStyle border;
};

}}}
