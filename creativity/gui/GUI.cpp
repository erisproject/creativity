#include "creativity/gui/GUI.hpp"
#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/BookStore.hpp"
#include "creativity/gui/ReaderInfoWindow.hpp"
#include "creativity/gui/BookInfoWindow.hpp"
#include "creativity/gui/GraphArea.hpp"
#include "creativity/gui/InfoWindow.hpp"
#include "creativity/state/State.hpp"
#include "creativity/state/Storage.hpp"
#include "creativity/Creativity.hpp"
#include "creativity/config.hpp"
#include <eris/Random.hpp>
#include <boost/filesystem/path.hpp>
#include <cairomm/matrix.h>
#include <cairomm/pattern.h>
#include <cairomm/refptr.h>
#include <gdkmm/cursor.h>
#include <gdkmm/device.h>
#include <gdkmm/rgba.h>
#include <gdkmm/window.h>
#include <giomm/application.h>
#include <glibmm/dispatcher.h>
#include <glibmm/fileutils.h>
#include <glibmm/main.h>
#include <glibmm/refptr.h>
#include <glibmm/signalproxy.h>
#include <glibmm/ustring.h>
#include <gtkmm/application.h>
#include <gtkmm/box.h>
#include <gtkmm/builder.h>
#include <gtkmm/button.h>
#include <gtkmm/checkbutton.h>
#include <gtkmm/colorbutton.h>
#include <gtkmm/combobox.h>
#include <gtkmm/comboboxtext.h>
#include <gtkmm/dialog.h>
#include <gtkmm/entry.h>
#include <gtkmm/enums.h>
#include <gtkmm/expander.h>
#include <gtkmm/filechooser.h>
#include <gtkmm/filechooserdialog.h>
#include <gtkmm/filefilter.h>
#include <gtkmm/frame.h>
#include <gtkmm/grid.h>
#include <gtkmm/label.h>
#include <gtkmm/main.h>
#include <gtkmm/messagedialog.h>
#include <gtkmm/notebook.h>
#include <gtkmm/radiobutton.h>
#include <gtkmm/scale.h>
#include <gtkmm/scrolledwindow.h>
#include <gtkmm/spinbutton.h>
#include <gtkmm/stock.h>
#include <gtkmm/treeiter.h>
#include <gtkmm/treemodel.h>
#include <gtkmm/treemodelcolumn.h>
#include <gtkmm/treepath.h>
#include <gtkmm/treeselection.h>
#include <gtkmm/treeview.h>
#include <gtkmm/treeviewcolumn.h> // IWYU pragma: keep
#include <gtkmm/viewport.h>
#include <gtkmm/widget.h>
#include <gtkmm/window.h>
#include <sigc++/connection.h>
#include <cstddef>
#include <exception>
#include <iterator>
#include <set>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>

using namespace eris;
using namespace creativity::state;
using namespace std::placeholders;

namespace creativity { namespace gui {

GUI::GUI(std::shared_ptr<Creativity> creativity,
        std::function<void(Parameter)> configure,
        std::function<void()> initialize,
        std::function<void(unsigned int count)> run,
        std::function<void()> stop,
        std::function<void()> resume,
        std::function<void()> step,
        std::function<void()> quit)
    :
        creativity_{std::move(creativity)},
        on_configure_{std::move(configure)},
        on_initialize_{std::move(initialize)},
        on_run_{std::move(run)},
        on_stop_{std::move(stop)},
        on_resume_{std::move(resume)},
        on_step_{std::move(step)},
        on_quit_{std::move(quit)}
{}

double GUI::sb(const std::string &widget_name) {
    return widget<Gtk::SpinButton>(widget_name)->get_value();
}
int GUI::sb_int(const std::string &widget_name) {
    return widget<Gtk::SpinButton>(widget_name)->get_value_as_int();
}

void GUI::start(const cmdargs::GUI &args) {
    if (gui_thread_.joinable())
        throw std::runtime_error("GUI thread can only be started once!");

    app_ = Gtk::Application::create("ca.imaginary.eris.creativity", Gio::APPLICATION_NON_UNIQUE);

    builder_ = Gtk::Builder::create();

    std::list<std::string> datadirs;
    for (const auto &dir : DATADIR) {
        datadirs.push_back(dir);
    }
    datadirs.push_back(".");
    char *envdatadir = getenv("CREATIVITY_DATADIR");
    if (envdatadir) {
        std::string stddatadir(envdatadir);
        if (stddatadir != "") datadirs.push_front(stddatadir);
    }

    bool gui_glade_loaded = false;
    for (auto &dir : datadirs) {
        if (Glib::file_test(dir + "/gui.glade", Glib::FILE_TEST_IS_REGULAR)) {
            builder_->add_from_file(dir + "/gui.glade");
            gui_glade_loaded = true;
        }
    }
    if (not gui_glade_loaded) {
        throw Glib::FileError(Glib::FileError::Code::NO_SUCH_ENTITY, "Unable to find gui.glade; try setting CREATIVITY_DATADIR to the directory containing gui.glade");
    }

    // Have to set the initial seed value *before* starting the GUI thread, because we want the
    // seed() call to affect the main thread, not the GUI thread.
    widget<Gtk::Entry>("set_seed")->set_text(std::to_string(eris::Random::seed()));

    // Update the number of periods now, so that it has the right value when the GUI comes up
    widget<Gtk::SpinButton>("set_periods")->set_value(args.periods);

    gui_thread_ = std::thread(&GUI::thr_run, this, args);

    // Wait for the thread to finish its startup before returning
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this]() -> bool { return thread_running_; });
    lock.unlock();
}

GUI::~GUI() {
    mutex_.lock();
    if (dispatcher_) {
        signal_queue_.emplace_back(Signal::Type::quit);
        dispatcher_->emit();
    }
    mutex_.unlock();
    if (gui_thread_.joinable())
        gui_thread_.join();
}

void GUI::thr_run(const cmdargs::GUI &args) {
    main_window_ = std::shared_ptr<Gtk::Window>(widget<Gtk::Window>("window1"));

    main_window_title_ = main_window_->get_title() + " v" + std::to_string(VERSION[0]) + "." + std::to_string(VERSION[1]) + "." + std::to_string(VERSION[2]);
    main_window_->set_title(main_window_title_);

    auto disable_on_load = {"fr_agents", "fr_run", "fr_sim"};
    // When the "load" radio button is activated, disable all the new simulation parameter fields
    widget<Gtk::RadioButton>("radio_load")->signal_clicked().connect([this,disable_on_load] {
        for (const auto &fr : disable_on_load) widget<Gtk::Frame>(fr)->set_sensitive(false);
    });
    // Undo the above when "New simulation" is activated
    widget<Gtk::RadioButton>("radio_new")->signal_clicked().connect([this,disable_on_load] {
            for (const auto &fr : disable_on_load) widget<Gtk::Frame>(fr)->set_sensitive(true);
    });
    // Clicking the load file chooser button directly first fires the load radio, triggering the above
    widget<Gtk::Button>("btn_load")->signal_clicked().connect([this] {
        widget<Gtk::RadioButton>("radio_load")->set_active(true);

        Gtk::FileChooserDialog fdlg(*main_window_, "Select simulation data output file", Gtk::FILE_CHOOSER_ACTION_OPEN);
        fdlg.set_modal(true);
        fdlg.set_urgency_hint(true);
        fdlg.set_skip_taskbar_hint(true);
        fdlg.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
        fdlg.add_button(Gtk::Stock::OPEN, Gtk::RESPONSE_ACCEPT);
        fdlg.set_filter(fileFilter());
        if (load_ != "")
            fdlg.set_filename(load_);
        else
            fdlg.set_current_folder(".");

        auto result = fdlg.run();
        if (result == Gtk::RESPONSE_ACCEPT) {
            load_ = fdlg.get_filename();
            loadSim();
        }
        else {
            load_ = "";
        }
    });

    // Clicking the save file chooser likewise fires the save radio, then also starts up the file
    // chooser dialog.
    widget<Gtk::Button>("btn_choose_save")->signal_clicked().connect([this] {
        widget<Gtk::RadioButton>("radio_save")->set_active(true);
        Gtk::FileChooserDialog fdlg(*main_window_, "Select simulation data output file", Gtk::FILE_CHOOSER_ACTION_SAVE);
        fdlg.set_do_overwrite_confirmation(true);
        fdlg.set_modal(true);
        fdlg.set_urgency_hint(true);
        fdlg.set_skip_taskbar_hint(true);
        fdlg.add_button(Gtk::Stock::CANCEL,Gtk::RESPONSE_CANCEL);
        fdlg.add_button(Gtk::Stock::SAVE,Gtk::RESPONSE_ACCEPT);
        fdlg.set_filter(fileFilter());
        if (save_ != "")
            fdlg.set_filename(save_);
        else {
            fdlg.set_current_folder(".");
            fdlg.set_current_name("creativity-" + widget<Gtk::Entry>("set_seed")->get_text() + ".crstate");
        }

        auto result = fdlg.run();
        if (result == Gtk::RESPONSE_ACCEPT) {
            save_ = fdlg.get_filename();
            saveSim();
        }
        else {
            save_ = "";
        }
    });

    widget<Gtk::Button>("btn_start")->signal_clicked().connect([this] {
        initializeSim();
        runSim();
    });

    widget<Gtk::Button>("btn_init")->signal_clicked().connect([this] {
        initializeSim();
    });

    widget<Gtk::Button>("btn_run")->signal_clicked().connect([this] {
        runSim();
    });

    widget<Gtk::Button>("btn_pause")->signal_clicked().connect([this] {
        queueEvent(Event::Type::stop);
    });

    widget<Gtk::Button>("btn_resume")->signal_clicked().connect([this] {
        queueEvent(Event::Type::resume);
    });

    widget<Gtk::Button>("btn_step")->signal_clicked().connect([this] {
        queueEvent(Event::Type::step);
    });

    auto thrbox = widget<Gtk::ComboBoxText>("combo_threads");
    // The gui setup has two entries: no threads, and 1 thread.  The latter is pointless (except for
    // eris debugging), so remove it, but first check to see if it was selected; if it was, we'll
    // update the default selection to the maximum number of threads; otherwise we'll leave it on no
    // threads.
    thrbox->remove_text(1);
    unsigned int last = 1;
    bool added_requested_threads = args.threads == 0;
    // List common number of threads by going up by increasing amounts so that the increment is
    // always greater than 25% and at most 50%, so we get:
    // 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, ...
    // But stop that, of course, at the maximum number of hardware threads supported.
    for (unsigned int i = 2, incr = 1; i <= std::thread::hardware_concurrency(); i += incr) {
        // If the user requested some oddball number that we skip with the next one, we need to add it
        if (not added_requested_threads and i >= args.threads) {
            if (i > args.threads) {
                // We're about to skip over it, so add it
                thrbox->append(std::to_string(args.threads), std::to_string(args.threads) +
                        (args.threads > 1 ? " threads" : " thread"));
            }
            // Otherwise it's equal, which means we're about to add it
            added_requested_threads = true;
        }

        thrbox->append(std::to_string(i), std::to_string(i) + " threads");
        if (i >= 4*incr) incr *= 2;
        last = i;
    }
    // We might also need to add the hardware value and the requested value; the former if it is
    // some oddball amount not caught above; the latter similarly, or if it exceeds the number of
    // hardware threads.
    unsigned int more[] = {std::thread::hardware_concurrency(), args.threads};
    if (more[1] < more[0]) std::swap(more[0], more[1]);

    for (auto &t : more) {
        if (t > last) {
            std::string num = std::to_string(t);
            thrbox->append(num, num + " threads");
            last = t;
        }
    }

    // Select the requested (or default) number of threads
    thrbox->set_active_id(std::to_string(args.threads));

    thrbox->signal_changed().connect([this] {
        Parameter th;
        th.param = ParamType::threads; th.ul = std::stoul(widget<Gtk::ComboBoxText>("combo_threads")->get_active_id());
        queueEvent(th);
    });

    widget<Gtk::Entry>("set_seed")->signal_icon_press().connect([this](Gtk::EntryIconPosition icon_position, const GdkEventButton*) -> void {
        if (icon_position == Gtk::EntryIconPosition::ENTRY_ICON_SECONDARY)
            widget<Gtk::Entry>("set_seed")->set_text(std::to_string(std::random_device{}()));
    });

    widget<Gtk::Scale>("scale_state")->signal_value_changed().connect([this] {
        unsigned long t = std::lround(widget<Gtk::Scale>("scale_state")->get_value());
        thr_set_state(t);
    });

    Gtk::Viewport *vis;
    builder_->get_widget("view_vis", vis);

    // Lock access to everything until we're set up.  We'll unlock and notify on cv_ once the
    // mainloop is running.
    std::unique_lock<std::mutex> lock{mutex_};

    graph_ = std::unique_ptr<GraphArea>(new GraphArea(*this));

    hand_ = Gdk::Cursor::create(main_window_->get_display(), Gdk::HAND1);

    motion_handler_ = [this](GdkEventMotion *event) -> bool {
        double x = event->x, y = event->y;
        auto c2g = graph_->canvas_to_graph();
        c2g.transform_point(x, y);
        rt_point clicked(x, y);

        double thresh_x = 5.0, thresh_y = -5.0;
        c2g.transform_distance(thresh_x, thresh_y);

        // When the mouse moves over the image, change the cursor to a hand any time the cursor
        // is within 5 pixels (on both dimensions) of a graph item (book, reader, etc.).
        auto nearest = thr_nearest(clicked);

        if (not nearest.empty() and std::fabs(nearest[0].first.get<0>() - x) <= thresh_x
                and std::fabs(nearest[0].first.get<1>() - y) <= thresh_y)
            main_window_->get_window()->set_cursor(hand_);
        else
            main_window_->get_window()->set_cursor();

        return false;
    };
    graph_->add_events(Gdk::POINTER_MOTION_MASK);
    motion_handler_conn_ = graph_->signal_motion_notify_event().connect(motion_handler_);

    graph_->add_events(Gdk::EventMask::BUTTON_PRESS_MASK);
    graph_->signal_button_press_event().connect([this](GdkEventButton *event) -> bool {
        double x = event->x, y = event->y;
        auto c2g = graph_->canvas_to_graph();
        c2g.transform_point(x, y);
        rt_point clicked(x, y);

        // Figure out the closest Reader and Book on the canvas relative to where the click happened
        auto nearest = thr_nearest(clicked);

        // If the nearest point isn't within 5 screen pixels (in each dimension) of the click,
        // ignore the click.
        double thresh_x = 5.0, thresh_y = -5.0;
        c2g.transform_distance(thresh_x, thresh_y);

        if (nearest.empty() or std::fabs(nearest[0].first.get<0>() - x) > thresh_x
                or std::fabs(nearest[0].first.get<1>() - y) > thresh_y)
            return false;

        thr_info_dialog(nearest[0].second);
        return true;
    });

    graph_->add_events(Gdk::EventMask::SCROLL_MASK);
    // Gtk 3.16.1 reversed the scroll direction of the mousewheel on horizontal scales; having the
    // scale above the graph scroll in the opposite direction as scrolling on the graph itself is
    // really weird, so go with the flow for scrolling on the graph itself
    const GdkScrollDirection scroll_back =
        (::gtk_get_major_version() > 3 or (::gtk_get_major_version() == 3 and
             (::gtk_get_minor_version() > 16 or (::gtk_get_minor_version() == 16 and ::gtk_get_micro_version() >= 1))))
        ? /* >= 3.16.1 */ GdkScrollDirection::GDK_SCROLL_DOWN
        : /* < 3.16.1 */ GdkScrollDirection::GDK_SCROLL_UP;
    const GdkScrollDirection scroll_fwd = scroll_back == GdkScrollDirection::GDK_SCROLL_DOWN
        ? GdkScrollDirection::GDK_SCROLL_UP
        : GdkScrollDirection::GDK_SCROLL_DOWN;

    graph_->signal_scroll_event().connect([this,scroll_back,scroll_fwd](GdkEventScroll *event) -> bool {
        if (event->direction == scroll_back) {
            if (state_curr_ > 0) thr_set_state(state_curr_ - 1);
            return true;
        }
        else if (event->direction == scroll_fwd) {
            if (state_curr_ + 1 < state_num_) thr_set_state(state_curr_ + 1);
            return true;
        }
        return false;
    });


    // Set up handlers and defaults for the visualization settings
#define GUI_SETUP_VIS_COLOUR(FIELD) {\
        auto colour_button = widget<Gtk::ColorButton>("colour_" #FIELD); \
        double r, g, b, a; Gdk::RGBA rgba; \
        graph_->design.colour.FIELD->get_rgba(r, g, b, a); \
        rgba.set_red(r); rgba.set_green(g); rgba.set_blue(b); rgba.set_alpha(a); \
        colour_button->set_rgba(rgba); \
        colour_button->signal_color_set().connect([this,colour_button] { \
            auto colour = colour_button->get_rgba(); \
            graph_->design.colour.FIELD = Cairo::SolidPattern::create_rgba(colour.get_red(), colour.get_green(), colour.get_blue(), colour.get_alpha()); \
            graph_->resetCache(); \
            graph_->queue_draw(); \
        }); \
    }
#define GUI_SETUP_VIS_SETTING_(FIELD, AFFECTS_RTREE) {\
        auto enable_button = widget<Gtk::CheckButton>("enable_" #FIELD); \
        enable_button->set_active(graph_->design.enabled.FIELD); \
        enable_button->signal_toggled().connect([this,enable_button] { \
            bool was = graph_->design.enabled.FIELD; \
            bool now = enable_button->get_active(); \
            if (was != now) { \
                graph_->design.enabled.FIELD = now; \
                graph_->resetCache(); \
                graph_->queue_draw(); \
                if (AFFECTS_RTREE) thr_reset_rtrees(); \
            } \
        }); \
    } \
    GUI_SETUP_VIS_COLOUR(FIELD)
#define GUI_SETUP_VIS_SETTING(FIELD) GUI_SETUP_VIS_SETTING_(FIELD, false)
#define GUI_SETUP_VIS_SETTING_AFFECTS_RTREE(FIELD) GUI_SETUP_VIS_SETTING_(FIELD, true)

    GUI_SETUP_VIS_SETTING_AFFECTS_RTREE(reader)
    GUI_SETUP_VIS_SETTING_AFFECTS_RTREE(book_live)
    GUI_SETUP_VIS_SETTING_AFFECTS_RTREE(book_dead)
    GUI_SETUP_VIS_SETTING_AFFECTS_RTREE(book_public)
    GUI_SETUP_VIS_SETTING(friendship)
    GUI_SETUP_VIS_SETTING(movement)
    GUI_SETUP_VIS_SETTING(author_live)
    GUI_SETUP_VIS_SETTING(author_dead)
    GUI_SETUP_VIS_SETTING(author_public)
    GUI_SETUP_VIS_SETTING(reading)
    GUI_SETUP_VIS_SETTING(utility_gain)
    GUI_SETUP_VIS_SETTING(utility_loss)
    GUI_SETUP_VIS_SETTING(axes)
    GUI_SETUP_VIS_COLOUR(background)
#undef GUI_SETUP_VIS_SETTING
#undef GUI_SETUP_VIS_COLOUR

    // Copy parameters (which will be the defaults) into the Agent Attributes settings
    thr_update_parameters();

    // Update the attributes tooltips: only the SpinButtons in the glade file have tooltips, so copy
    // those tooltips into the associated Label as well (which is just left of the spinbutton in the
    // grid)
    auto attribs_grid = widget<Gtk::Grid>("grid_sim_params");
    bool done_copying_tooltips = false;
    for (int c = 1; not done_copying_tooltips; c += 2) {
        for (int r = 1; ; r++) {
            auto *widget = attribs_grid->get_child_at(c, r), *left = attribs_grid->get_child_at(c-1, r);
            if (not widget or not left) {
                if (r == 1) done_copying_tooltips = 1;
                break;
            }
            // Don't copy the tooltip left if there is no tooltip, or the left element already has a
            // tooltip
            auto tooltip = widget->get_tooltip_markup();
            if (not tooltip.empty() and left->get_tooltip_markup().empty())
                left->set_tooltip_markup(tooltip);
        }
    }


    rdr_win_ = widget<Gtk::ScrolledWindow>("win_rdr");
    bk_win_ = widget<Gtk::ScrolledWindow>("win_bk");

    dispatcher_ = std::unique_ptr<Glib::Dispatcher>(new Glib::Dispatcher);
    dispatcher_->connect([this] { thr_signal(); });

    vis->add(*graph_);
    graph_->show();

    Glib::signal_idle().connect_once([this,&lock,&args] {
            if (not args.input.empty()) {
                load_ = args.input;
            }
            if (not args.output.empty()) {
                save_ = args.output;
                widget<Gtk::Label>("lbl_save")->set_text(save_);
                main_window_->set_title(main_window_title_ + " [" + boost::filesystem::path(save_).filename().string() + "]");
            }
            bool init = false, run = false;
            if (args.start) {
                init = run = true;
            }
            else if (args.initialize) {
                init = true;
            }

            // Have to do this *after* doing everything we want with args: GUI::start() waits on
            // this lock, and our args reference isn't guaranteed to survive once start() returns.
            // The actual calls to initializeSim() etc. can't be above, however, because they send a
            // signal which would deadlock before releasing the thread startup lock:
            thread_running_ = true;
            lock.unlock();
            cv_.notify_all();

            if (not load_.empty())
                loadSim();
            if (init) {
                initializeSim();
                if (run)
                    runSim();
            }

    });

    app_->run(*main_window_);

    queueEvent(Event::Type::quit);
    lock.lock();
    dispatcher_.reset();
}

void GUI::thr_set_state(unsigned long t) {
    if (t == state_curr_ or state_num_ == 0) return;

    std::shared_ptr<const State> state;
    {
        size_t storage_size;
        {
            auto st = creativity_->storage();
            state = (*st.first)[t]; // Will throw if `t` is invalid
            storage_size = st.first->size();
        }

        // Enlarge if necessary
        if (rdr_models_.size() < storage_size) rdr_models_.resize(storage_size);
        if (rdr_trees_.size() < storage_size) rdr_trees_.resize(storage_size);
        if (bk_models_.size() < storage_size) bk_models_.resize(storage_size);
        if (bk_trees_.size() < storage_size) bk_trees_.resize(storage_size);
        if (rtrees_.size() < storage_size) rtrees_.resize(storage_size);
    }

    // Values may be null pointers--if so, create new model stores
    if (not rdr_models_[t]) rdr_models_[t] = ReaderStore::create(state);
    if (not bk_models_[t]) bk_models_[t] = BookStore::create(state);

    Gtk::TreeModel::Path rdr_select, bk_select;

    if (state_curr_ == (unsigned long) -1) {
        // This is the first state, replacing the window initialization fake state.
        //
        // Apply default sort orders (readers go in ID ascending, books go by age ascending).  Do
        // this here rather that {Reader,Book}Store::create because the sort order will be copied to
        // new state transitions; if the user changes the sort order, it's that new order rather
        // than this default order that we want to apply.
        rdr_models_[t]->set_sort_column(rdr_models_[t]->columns.id, Gtk::SortType::SORT_ASCENDING);
        bk_models_[t]->set_sort_column(bk_models_[t]->columns->age, Gtk::SortType::SORT_ASCENDING);
    }
    else {
        // Transitioning from one (actual) state to another: preserve sort order and reader/book
        // selection from the old treeview&model to the newly selected treeview&model.
        int old_sort_col, new_sort_col;
        Gtk::SortType old_sort_order, new_sort_order;

        // Readers:
        rdr_models_[state_curr_]->get_sort_column_id(old_sort_col, old_sort_order);
        rdr_models_[t]->get_sort_column_id(new_sort_col, new_sort_order);
        if (old_sort_col != new_sort_col or old_sort_order != new_sort_order)
            rdr_models_[t]->set_sort_column(old_sort_col, old_sort_order);

        // Books:
        bk_models_[state_curr_]->get_sort_column_id(old_sort_col, old_sort_order);
        bk_models_[t]->get_sort_column_id(new_sort_col, new_sort_order);
        if (old_sort_col != new_sort_col or old_sort_order != new_sort_order)
            bk_models_[t]->set_sort_column(old_sort_col, old_sort_order);

        // Figure out if a reader/book is selected so that we can also select it in the new model
        auto rdr_sel_iter = rdr_trees_[state_curr_]->get_selection()->get_selected();
        if (rdr_sel_iter) {
            auto selected_id = rdr_models_[state_curr_]->member(rdr_sel_iter).id;
            rdr_select = rdr_models_[t]->find(selected_id, rdr_sel_iter);
        }

        auto bk_sel_iter = bk_trees_[state_curr_]->get_selection()->get_selected();
        if (bk_sel_iter) {
            auto selected_id = bk_models_[state_curr_]->member(bk_sel_iter).id;
            bk_select = bk_models_[t]->find(selected_id, bk_sel_iter);
        }

    }

    // Make sure we have properly set up Gtk::TreeView references and not null pointers
    if (not rdr_trees_[t]) {
        auto *tv = new Gtk::TreeView;
        rdr_models_[t]->appendColumnsTo(*tv);
        tv->set_fixed_height_mode(true);
        tv->signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
            thr_info_dialog(rdr_models_[state_curr_]->member(path).id);
        });
        tv->set_model(rdr_models_[t]);
        tv->show();
        rdr_trees_[t].reset(tv);
    }
    if (not bk_trees_[t]) {
        auto *tv = new Gtk::TreeView;
        bk_models_[t]->appendColumnsTo(*tv);
        tv->set_fixed_height_mode(true);
        tv->signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
            thr_info_dialog(bk_models_[state_curr_]->member(path).id);
        });
        tv->set_model(bk_models_[t]);
        tv->show();
        bk_trees_[t].reset(tv);
    }


    // Apply the previous selection to the new treeviews
    if (rdr_select.empty())
        rdr_trees_[t]->get_selection()->unselect_all();
    else {
        rdr_trees_[t]->get_selection()->select(rdr_select);
        rdr_trees_[t]->scroll_to_row(rdr_select);
    }

    if (bk_select.empty())
        bk_trees_[t]->get_selection()->unselect_all();
    else {
        bk_trees_[t]->get_selection()->select(bk_select);
        bk_trees_[t]->scroll_to_row(bk_select);
    }

    // Replace the current treeview in the window with the new one
    rdr_win_->remove();
    rdr_win_->add(*rdr_trees_[t]);

    bk_win_->remove();
    bk_win_->add(*bk_trees_[t]);

    if (not rtrees_[t]) {
        // Stick all points into an rtree so that we can quickly find the nearest one (needed, in
        // particular, for fast mouseovers).
        //
        // Note that we store graph positions, not canvas positions: canvas positions change when
        // the graph size changes (i.e. the user resizes the window).  Rather than invalidating and
        // needing to rebuild the rtrees when that happens, we store graph positions and translate
        // to graph positions as needed in the mouseover/click handlers.
        //
        // Also note that only *visible* points get added here (according to
        // graph_->design.enabled.*)
        rtrees_[t] = std::unique_ptr<RTree>(new RTree);
        thr_init_rtree(*rtrees_[t], state);
    }

    // Now go through any open reader/book info dialog windows: delete any that have been closed,
    // and refresh the information on any still open.
    std::vector<eris_id_t> del;
    // NB: can't do this in a single pass without the extra erase() lookups: the iteration order of
    // unordered_map after an .erase() isn't guaranteed to be the same until C++14 (C++11 only
    // guarantees that iterators remain valid, but not necessarily in the same order).
    for (auto &w : info_windows_) {
        if (w.second->get_visible())
            w.second->refresh(state);
        else
            del.push_back(w.first);
    }
    for (auto &d : del) {
        info_windows_.erase(d);
    }

    state_curr_ = t;
    // This must be *after* state_curr_ gets updated, otherwise it'll trigger a recursive call
    widget<Gtk::Scale>("scale_state")->set_value(t);

    widget<Gtk::Label>("lbl_tab_agents")->set_text("Agents (" + std::to_string(state->readers.size()) + ")");
    widget<Gtk::Label>("lbl_tab_books")->set_text("Books (" + std::to_string(state->books.size()) + ")");

    // Update the summary tab fields
#define GUI_SUMMARY_CALC(NAME, VALUE) do { std::ostringstream os; os << VALUE; widget<Gtk::Label>("summary_"#NAME)->set_text(os.str()); } while (0)
#define GUI_SUMMARY(NAME) GUI_SUMMARY_CALC(NAME, NAME)
    GUI_SUMMARY(t);
    GUI_SUMMARY_CALC(agents, state->readers.size());
    double n = (double) state->readers.size(), u = 0, ulife = 0, book_price = 0, book_revenue = 0, book_quality = 0;
    int authors = 0, last_wrote = 0, books_new = 0, books_market = 0, books_w_sales = 0, books_pirated = 0, book_sales = 0, book_pirated = 0,
        readers_bought = 0, readers_pirated = 0, readers_either = 0, library_bought = 0, library_pirated = 0;
    for (auto &rp : state->readers) {
        auto &r = rp.second;
        u += r.u;
        ulife += r.u_lifetime;
        if (r.wrote.size() > 0) { authors++; last_wrote += t - state->books.at(*r.wrote.crbegin()).created; }
        if (r.library_purchased_new > 0) { readers_bought++; }
        if (r.library_pirated_new > 0) { readers_pirated++; }
        if (r.library_purchased_new > 0 or r.library_pirated_new > 0) { readers_either++; }
        library_bought += r.library_purchased;
        library_pirated += r.library_pirated;
    }
    for (auto &bp : state->books) {
        auto &b = bp.second;
        if (b.created == t) { books_new++; }
        // FIXME: what about distinguishing public/private sales?
        if (b.market_private) { books_market++; book_sales += b.sales; book_price += b.price; book_revenue += b.revenue; book_quality += b.quality; }
        if (b.sales > 0) { books_w_sales++; }
        if (b.pirated > 0) { books_pirated++; book_pirated += b.pirated; }
    }
    GUI_SUMMARY_CALC(u, u / n);
    GUI_SUMMARY_CALC(ulife, ulife / n);
    GUI_SUMMARY(authors);
    GUI_SUMMARY_CALC(last_wrote, last_wrote / (double) authors);
    GUI_SUMMARY(books_new);
    GUI_SUMMARY(books_market);
    GUI_SUMMARY(books_w_sales);
    GUI_SUMMARY(books_pirated);
    GUI_SUMMARY_CALC(book_sales, book_sales / (double) books_market);
    GUI_SUMMARY_CALC(book_price, book_price / books_market);
    GUI_SUMMARY_CALC(book_revenue, book_revenue / books_market);
    GUI_SUMMARY_CALC(book_quality, book_quality / books_market);
    GUI_SUMMARY_CALC(book_pirated, book_pirated / (double) books_pirated);
    GUI_SUMMARY(readers_bought);
    GUI_SUMMARY(readers_pirated);
    GUI_SUMMARY(readers_either);
    GUI_SUMMARY_CALC(library_bought, library_bought / n);
    GUI_SUMMARY_CALC(library_pirated, library_pirated / n);
#undef GUI_SUMMARY
#undef GUI_SUMMARY_CALC

    if (graph_->get_is_drawable()) graph_->queue_draw();
}

void GUI::thr_reset_rtrees() {
    rtrees_.clear();
    if (state_num_ == 0) return; // No states to worry about!
    auto st = creativity_->storage();
    rtrees_.resize(st.first->size());
    auto state = (*st.first)[state_curr_];

    rtrees_[state_curr_] = std::unique_ptr<RTree>(new RTree);
    thr_init_rtree(*rtrees_[state_curr_], state);
}

void GUI::thr_init_rtree(RTree &rt, const std::shared_ptr<const State> &state) const {
    if (graph_->design.enabled.reader) {
        for (auto &r : state->readers) {
            auto gp = graph_->graph_position(r.second.position);
            rt.insert(std::make_pair(rt_point{gp[0], gp[1]}, r.second.id));
        }
    }
    if (graph_->design.enabled.book_live or graph_->design.enabled.book_dead) {
        for (auto &b : state->books) {
            if (b.second.market_private ? graph_->design.enabled.book_live :
                    b.second.market_public() ? graph_->design.enabled.book_public : graph_->design.enabled.book_dead) {
                auto gp = graph_->graph_position(b.second.position);
                rt.insert(std::make_pair(rt_point{gp[0], gp[1]}, b.second.id));
            }
        }
    }
}

std::vector<GUI::rt_val> GUI::thr_nearest(const rt_point &point, int n) {
    std::vector<rt_val> nearest;
    if (state_num_ > 0)
        rtrees_[state_curr_]->query(boost::geometry::index::nearest(point, n), std::back_inserter(nearest));
    return nearest;
}

decltype(Gtk::FileFilter::create()) GUI::fileFilter() const {
    if (not ff_) {
        ff_ = Gtk::FileFilter::create();
        ff_->set_name("Creativity simulation state files");
        ff_->add_pattern("*.crstate");
    }
    return ff_;
}

void GUI::thr_info_dialog(eris_id_t member_id) {
    auto already_open = info_windows_.find(member_id);
    if (already_open != info_windows_.end()) {
        // If the window is already active, just present it again
        already_open->second->present();
    }
    else {
        std::shared_ptr<const State> state((*creativity_->storage().first)[state_curr_]);

        if (state->readers.count(member_id))
            info_windows_.emplace(
                member_id, std::unique_ptr<InfoWindow>(
                    new ReaderInfoWindow(state, main_window_, member_id, std::bind(&GUI::thr_info_dialog, this, _1))
                )
            );
        else if (state->books.count(member_id))
            info_windows_.emplace(
                member_id, std::unique_ptr<InfoWindow>(
                    new BookInfoWindow(state, main_window_, member_id)
                )
            );
        else
            throw std::out_of_range("thr_info_dialog: requested member id does not exist");

    }
}

void GUI::thr_update_parameters() {
#define SET_SB(PARAMETER) widget<Gtk::SpinButton>("set_" #PARAMETER)->set_value(creativity_->parameters.PARAMETER)
    SET_SB(dimensions);
    SET_SB(readers);
    widget<Gtk::SpinButton>("set_density")->set_value(creativity_->densityFromBoundary());
    SET_SB(reader_step_mean);
    SET_SB(book_distance_mean);
    SET_SB(book_quality_sd);
    SET_SB(reader_creation_shape);
    SET_SB(reader_creation_scale_min);
    SET_SB(reader_creation_scale_max);
    SET_SB(creation_time);
    SET_SB(creation_fixed);

#define SET_INIT_SB(PARAMETER) widget<Gtk::SpinButton>("set_init_" #PARAMETER)->set_value(creativity_->parameters.initial.PARAMETER)
    SET_INIT_SB(prob_write);
    SET_INIT_SB(q_min);
    SET_INIT_SB(q_max);
    SET_INIT_SB(p_min);
    SET_INIT_SB(p_max);
    SET_INIT_SB(prob_keep);
    SET_INIT_SB(keep_price);
#undef SET_INIT_SB

    SET_SB(income);
    SET_SB(cost_market);
    SET_SB(cost_unit);
    SET_SB(cost_piracy);

    SET_SB(prior_scale);
    SET_SB(prior_scale_burnin);
    SET_SB(burnin_periods);
    SET_SB(prediction_draws);

    SET_SB(piracy_begins);
    widget<Gtk::SpinButton>("set_piracy_link_proportion")->set_value(creativity_->parameters.piracy_link_proportion * 100.0);
    SET_SB(prior_scale_piracy);

    SET_SB(public_sharing_begins);
    SET_SB(public_sharing_tax);
    SET_SB(prior_scale_public_sharing);
#undef SET_SB
}

void GUI::thr_signal() {
    std::unique_lock<std::mutex> lock(mutex_);
    // Keep track of the *last* state change we see, and whether or not we see a new_states signal
    // (we want to skip everything except the last).
    std::stringstream errors;
    bool init = false;
    Signal last_state, last_progress, last_new_states;
    for (Signal &s : signal_queue_) {
        switch (s.type) {
            case Signal::Type::none:
                throw std::logic_error("Thread received invalid signal with type=none");
                break;
            case Signal::Type::quit:
                app_->quit();
                return;
                break;
            case Signal::Type::new_states:
                last_new_states = s;
                break;
            case Signal::Type::initialized:
                init = true;
            case Signal::Type::running:
            case Signal::Type::stopped:
                last_state = s;
                break;
            case Signal::Type::progress:
                last_progress = s;
                break;
            case Signal::Type::error:
                errors << s.message << "\n";
                break;
        }
    }
    signal_queue_.clear();
    lock.unlock();

    if (init) {
        widget<Gtk::Notebook>("nb_tabs")->set_current_page(1);

        // Disable spin buttons:
        for (auto &field : widget<Gtk::Grid>("grid_sim_params")->get_children()) {
            if (auto *sb = dynamic_cast<Gtk::SpinButton*>(field)) {
                sb->set_sensitive(false);
            }
        }

        // Update the agent attributes settings with the settings from the creativity object.
        // These may be different from the defaults, particularly when opening a file.
        thr_update_parameters();

        widget<Gtk::Button>("btn_init")->set_sensitive(false);
        widget<Gtk::Button>("btn_start")->set_sensitive(false);

        auto l = widget<Gtk::Label>("label_summary_t_desc");
        std::string d = l->get_text();
        auto pos = d.find("###");
        if (pos != d.npos) l->set_markup(d.replace(pos, 3, std::to_string(creativity_->parameters.piracy_begins)));

        std::ostringstream incstr; incstr << creativity_->parameters.income;
        for (auto &w : {"label_summary_u_desc", "label_summary_ulife_desc"}) {
            l = widget<Gtk::Label>(w);
            d = l->get_text();
            pos = d.find("###");
            if (pos != d.npos) l->set_markup(d.replace(pos, 3, incstr.str()));
        }
    }

    if (last_state.type == Signal::Type::initialized) {
        widget<Gtk::Button>("btn_run")->set_visible(true);
        widget<Gtk::Button>("btn_run")->set_sensitive(true);
        widget<Gtk::Button>("btn_resume")->set_visible(false);
        widget<Gtk::Button>("btn_pause")->set_visible(false);
        widget<Gtk::Button>("btn_step")->set_sensitive(true);
    }
    else if (last_state.type == Signal::Type::running) {
        widget<Gtk::Button>("btn_run")->set_visible(false);
        widget<Gtk::Button>("btn_run")->set_sensitive(true);
        widget<Gtk::Button>("btn_resume")->set_visible(false);
        widget<Gtk::Button>("btn_pause")->set_visible(true);
        widget<Gtk::Button>("btn_step")->set_sensitive(false);

        widget<Gtk::ComboBox>("combo_threads")->set_sensitive(false);
    }
    else if (last_state.type == Signal::Type::stopped) {
        // Turn off pause and turn on either play or resume buttons
        widget<Gtk::Button>("btn_run")->set_visible(!last_state.boolean);
        widget<Gtk::Button>("btn_resume")->set_visible(last_state.boolean);
        widget<Gtk::Button>("btn_pause")->set_visible(false);
        widget<Gtk::Button>("btn_step")->set_sensitive(true);

        widget<Gtk::ComboBox>("combo_threads")->set_sensitive(true);
    }
    if (last_progress.type == Signal::Type::progress) {
        // Progress bar update
        // .uls contains: t, end
        // .doubles has: speed
        auto scale = widget<Gtk::Scale>("scale_state");
        scale->set_range(0, last_progress.uls[1]);
        if (not piracy_tick_added_ and last_progress.uls[1] >= creativity_->parameters.piracy_begins) {
            widget<Gtk::Scale>("scale_state")->add_mark(creativity_->parameters.piracy_begins, Gtk::POS_BOTTOM, "");
            piracy_tick_added_ = true;
        }
        if (not public_tick_added_ and last_progress.uls[1] >= creativity_->parameters.public_sharing_begins) {
            widget<Gtk::Scale>("scale_state")->add_mark(creativity_->parameters.public_sharing_begins, Gtk::POS_BOTTOM, "");
            public_tick_added_ = true;
        }

        widget<Gtk::Label>("lbl_total")->set_text(std::to_string(last_progress.uls[1]));
        scale->set_fill_level(last_progress.uls[0]);
        scale->set_value(state_curr_);
        std::ostringstream speed;
        speed << std::setw(7) << std::showpoint << last_progress.doubles[0];

        widget<Gtk::Label>("status_speed")->set_text(speed.str());
    }

    std::string errstr = errors.str();
    if (not errstr.empty()) {
        Gtk::MessageDialog dlg_error(*main_window_, "An error has occured!", false, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_CLOSE, true);
        dlg_error.set_secondary_text(errstr);
        dlg_error.run();
        dlg_error.hide();
    }

    if (last_new_states.type == Signal::Type::new_states) {
        auto old_state_num = state_num_;
        state_num_ = creativity_->storage().first->size();

        if (state_num_ > old_state_num) {
            // - 1 here because we don't really count 0 as a stage
            //widget<Gtk::Label>("lbl_total")->set_text(std::to_string(state_num_-1));

            widget<Gtk::Scale>("scale_state")->set_fill_level(state_num_-1);

            // If the new_state signal doesn't carry a state instruction, and the user was already
            // on the last state, switch to the new last state:
            if (last_new_states.uls.empty() and (old_state_num == 0 or state_curr_ == old_state_num - 1))
                thr_set_state(state_num_-1);
        }

        if (not last_new_states.uls.empty()) {
            thr_set_state(std::min(state_num_-1, last_new_states.uls[0]));
        }
    }
}

void GUI::queueSignal(Signal &&s) {
    std::unique_lock<std::mutex> lock(mutex_);
    if (not dispatcher_) /// The thread has already gone away, ignore the signal
        return;
    bool empty = signal_queue_.empty();
    signal_queue_.push_back(std::move(s));
    if (empty) {
        // We only need to both signalling if the queue was empty (otherwise whatever added to the
        // queue last time already signalled, and the thread hasn't yet responded to that signal).
        dispatcher_->emit();
    }
}

// Load a simulation from a saved file:
void GUI::loadSim() {
    std::vector<Parameter> params;

    if (load_ == "") throw std::runtime_error("loadSim() called without load_ set to file to load");

    Parameter load;
    load.param = ParamType::load;
    load.ptr = &load_;
    queueEvent(load);

    widget<Gtk::Label>("lbl_load")->set_text(load_);

    main_window_->set_title(main_window_title_ + " [" + boost::filesystem::path(load_).filename().string() + "]");

    // Disable and hide the simulation controls (at the top of the window), since loading is read-only
    auto controls = widget<Gtk::Box>("box_controls");
    controls->set_sensitive(false);
    controls->set_visible(false);

    // Hide the speed indicator
    for (const auto &l : {"lbl_speed", "status_speed", "lbl_speed_units"})
        widget<Gtk::Label>(l)->set_visible(false);

    // Disable all the settings items that aren't modifiable anymore (basically everything except
    // the visualization settings and save button):
    for (const auto &disable : {"fr_new", "fr_agents", "fr_sim", "fr_run"})
        widget<Gtk::Frame>(disable)->set_sensitive(false);

    // Completely disable and hide the expanders for settings that are totally irrelevant
    for (const auto &kill : {"ex_sim", "ex_run"}) {
        auto ex = widget<Gtk::Expander>(kill);
        ex->set_expanded(false);
        ex->set_sensitive(false);
        ex->set_visible(false);
    }
    // Open the visualization settings (it's the only thing that's actually useful on the page)
    widget<Gtk::Expander>("ex_vis")->set_expanded(true);
}

// Save a simulation to a file.  This copies the current simulation states (if any) and keeps the
// file open for any future states.
void GUI::saveSim() {
    if (save_ == "") throw std::runtime_error("saveSim() called without save_ set to file to write");

    Parameter save;
    save.param = ParamType::save_as;
    save.ptr = &save_;
    queueEvent(save);

    widget<Gtk::Label>("lbl_save")->set_text(save_);

    main_window_->set_title(main_window_title_ + " [" + boost::filesystem::path(save_).filename().string() + "]");

    widget<Gtk::Button>("btn_load")->set_sensitive(false);
    widget<Gtk::RadioButton>("radio_load")->set_sensitive(false);
    widget<Gtk::RadioButton>("radio_memory")->set_sensitive(false);
}

// Initialize a new simulation:
void GUI::initializeSim() {
    // Set the parameters directly
    auto &set = creativity_->set();
#define COPY_SB_D(PARAMETER) set.PARAMETER = sb("set_"#PARAMETER)
#define COPY_SB_I(PARAMETER) set.PARAMETER = lround(sb("set_"#PARAMETER))
    COPY_SB_I(dimensions);
    COPY_SB_I(readers);
    set.boundary = Creativity::boundaryFromDensity(set.readers, set.dimensions, sb("set_density"));
    COPY_SB_I(piracy_begins);
    COPY_SB_D(piracy_link_proportion) * 0.01; // From percentage
    COPY_SB_D(book_distance_mean);
    COPY_SB_D(book_quality_sd);
    COPY_SB_D(reader_step_mean);
    COPY_SB_D(reader_creation_shape);
    COPY_SB_D(reader_creation_scale_min);
    COPY_SB_D(reader_creation_scale_max);
    COPY_SB_D(income);
    COPY_SB_D(cost_market);
    COPY_SB_D(cost_unit);
    COPY_SB_D(cost_piracy);
    COPY_SB_D(prior_scale);
    COPY_SB_D(prior_scale_piracy);
    COPY_SB_D(prior_scale_burnin);
    COPY_SB_I(burnin_periods);
    COPY_SB_I(prediction_draws);
    COPY_SB_I(creation_time);
    COPY_SB_I(creation_fixed);
#define COPY_SB_INIT_D(PARAMETER) set.initial.PARAMETER = sb("set_init_"#PARAMETER)
    COPY_SB_INIT_D(prob_write);
    COPY_SB_INIT_D(q_min);
    COPY_SB_INIT_D(q_max);
    COPY_SB_INIT_D(p_min);
    COPY_SB_INIT_D(p_max);
    COPY_SB_INIT_D(prob_keep);
    COPY_SB_INIT_D(keep_price);
#undef COPY_SB_D
#undef SET_INIT_SB

    // Other non-simulation options need to be queued via a configure event:
    Parameter p;
    p.param = ParamType::seed; p.ul = std::stoul(widget<Gtk::Entry>("set_seed")->get_text());
    queueEvent(p);
    p.param = ParamType::threads; p.ul = std::stoul(widget<Gtk::ComboBoxText>("combo_threads")->get_active_id());
    queueEvent(p);

    queueEvent(Event::Type::initialize);
}

void GUI::runSim() {
    queueEvent(Event::Type::run, sb_int("set_periods"));
}

void GUI::newStates(unsigned long switch_to) {
    queueSignal(
            (switch_to == (unsigned long)-1)
            ? Signal::Type::new_states
            : Signal(Signal::Type::new_states, switch_to)
            );
}

void GUI::initialized() {
    queueSignal(Signal::Type::new_states);
    queueSignal(Signal::Type::initialized);
}

void GUI::running() { queueSignal(Signal::Type::running); }

void GUI::progress(unsigned long t, unsigned long end, double speed) {
    Signal s{Signal::Type::progress};
    s.uls = {t, end};
    s.doubles = {speed};
    queueSignal(std::move(s));
}

void GUI::stopped(bool manual) { queueSignal({Signal::Type::stopped, manual}); }

void GUI::error(std::string message) { queueSignal({ Signal::Type::error, message }); }

GUI::Event::Event(Event::Type t) : type{t} {}
GUI::Event::Event(const Parameter &p) : type{Event::Type::configure}, parameter(p) {}
GUI::Event::Event(Event::Type t, unsigned long ulval): type{t}, ul{ulval} {}
GUI::Event::operator bool() { return type != Type::none; }

void GUI::checkEvents() {
    std::unique_lock<std::mutex> lock(mutex_);
    processEvents_(lock);
}

void GUI::waitEvents() {
    // Just like the above except that we wait until there is at least one item in the queue:
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this]() -> bool { return not event_queue_.empty(); });
    processEvents_(lock);
}

void GUI::waitForEvent(Event::Type t) {
    bool done = false;
    while (not done) {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this]() -> bool { return not event_queue_.empty(); });
        for (auto &e : event_queue_) {
            if (e.type == t) { done = true; break; }
        }
        processEvents_(lock);
    }
}

// Moves out the current queue elements, unlocks the lock, then processes the queue elements.
void GUI::processEvents_(std::unique_lock<decltype(mutex_)> &lock) {
    decltype(event_queue_) dequeued;
    dequeued.swap(event_queue_);
    lock.unlock();

    for (auto &e : dequeued) handleEvent(e);
}


void GUI::handleEvent(const Event &event) {
    std::vector<std::string> setup_errors;
    try {
        switch (event.type) {
            case Event::Type::configure:
                if (on_configure_) {
                    try {
                        on_configure_(event.parameter);
                    }
                    catch (const std::exception &e) {
                        setup_errors.push_back(e.what());
                    }
                }
                break;
            case Event::Type::initialize:
                if (on_initialize_) on_initialize_();
                break;
            case Event::Type::run:
                if (on_run_) on_run_(event.ul);
                break;
            case Event::Type::stop:
                if (on_stop_) on_stop_();
                break;
            case Event::Type::resume:
                if (on_resume_) on_resume_();
                break;
            case Event::Type::step:
                if (on_step_) on_step_();
                break;
            case Event::Type::quit:
                if (on_quit_) on_quit_();
                break;
            case Event::Type::none:
                break; // no op
        }
    } catch (std::exception &e) {
        queueSignal({ Signal::Type::error, std::string{"An exception has occured: "} + e.what() });
    }

    for (auto &e : setup_errors) {
        queueSignal({ Signal::Type::error, e });
    }
}

std::string GUI::pos_to_string(const eris::Position &pos) {
    std::ostringstream opos;
    opos << "(" << std::setprecision(3) << std::showpoint;
    for (size_t i = 0; i < pos.dimensions; i++) {
        if (i > 0) opos << ",";
        opos << pos[i];
    }
    opos << ")";
    return opos.str();
}


} }
