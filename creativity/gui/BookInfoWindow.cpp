#include "creativity/gui/BookInfoWindow.hpp"
#include "creativity/gui/GUI.hpp"

namespace creativity { namespace gui {

using namespace eris;
using namespace creativity::state;

BookInfoWindow::BookInfoWindow(
        std::shared_ptr<const state::State> state,
        std::shared_ptr<Gtk::Window> main_window,
        eris_id_t book
        )
    : InfoWindow(std::move(main_window), book, "Book details (" + std::to_string(book) + ")")
{
    auto &grid_status = new_grid();

    data_append(grid_status, "id", "ID");
    data_append(grid_status, "position", "Position");
    data_append(grid_status, "market", "Market");
    data_append(grid_status, "price", "Price");
    data_append(grid_status, "quality", "Quality");
    data_append(grid_status, "created", "Period written");
    data_append(grid_status, "age", "Age");
    data_append(grid_status, "revenue", "Revenue (lifetime)");
    data_append(grid_status, "revenueLast", "Revenue (current)");
    data_append(grid_status, "salesPrivate", "Copies sold (lifetime)");
    data_append(grid_status, "salesPublic", "Public market copies sold (lifetime)");
    data_append(grid_status, "salesLast", "Copies sold (current)");
    data_append(grid_status, "pirated", "Copies pirated (lifetime)");
    data_append(grid_status, "piratedLast", "Copies pirated (current)");
    data_append(grid_status, "copies", "Copies in world");
    data_append(grid_status, "author", "Author ID");
    override_background_color(Gdk::RGBA{"white"});

    add(grid_status);

    refresh(state);

    show_all();
}

void BookInfoWindow::refresh(std::shared_ptr<const State> state) {
    if (state->t == t_) return;

    if (state->books.count(id)) {
        set_title("Book details (" + std::to_string(id) + ")");
        auto &b = state->books.at(id);
        updateValue("id", b.id);
        updateValue("author", b.author);
        updateValue("position", GUI::pos_to_string(b.position));
        updateValue("market", b.market_private ? "yes (private)" : b.market_public() ? "yes (public)" : "no");
        updateValue("price", b.price);
        updateValue("quality", b.quality);
        updateValue("created", b.created);
        updateValue("age", state->t - b.created);
        updateValue("revenue", b.revenue_lifetime);
        updateValue("revenueLast", b.revenue);
        updateValue("salesPrivate", b.sales_lifetime_private);
        updateValue("salesPublic", b.sales_lifetime_public);
        updateValue("salesLast", b.sales);
        updateValue("pirated", b.pirated_lifetime);
        updateValue("piratedLast", b.pirated);
        updateValue("copies", b.copies_lifetime());
    }
    else {
        // If the user navigates back in time to a period where the book doesn't exist, set
        // everything to N/A and change the title.
        set_title("Book details (" + std::to_string(id) + "): not yet created!");
        for (auto &key : {"id", "author", "position", "market", "price", "quality", "age", "created",
                "revenue", "revenueLast", "salesPrivate", "salesPublic", "salesLast", "pirated", "piratedLast", "copies"}) {
            updateValue(key, "N/A");
        }
    }

    if (initial_refresh_) initial_refresh_ = false;
}

}}
