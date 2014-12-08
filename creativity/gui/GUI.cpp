#include "creativity/gui/GUI.hpp"
#include "creativity/gui/ReaderStore.hpp"
#include "creativity/gui/BookStore.hpp"
#include "creativity/Creativity.hpp"
#include "creativity/config.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <eris/Random.hpp>
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace eris;
using namespace creativity::state;
using namespace std::placeholders;

namespace creativity { namespace gui {

GUI::GUI(std::shared_ptr<Creativity> creativity,
        std::function<void(Parameter)> setup,
        std::function<void(unsigned int count)> run,
        std::function<void()> stop,
        std::function<void()> resume,
        std::function<void()> step,
        std::function<void()> quit)
    :
        creativity_{std::move(creativity)},
        on_setup_{std::move(setup)},
        on_run_{std::move(run)},
        on_stop_{std::move(stop)},
        on_resume_{std::move(resume)},
        on_step_{std::move(step)},
        on_quit_{std::move(quit)}
{}

GUI::Exception::Exception(const std::string &what) : std::runtime_error(what) {}

double GUI::sb(const std::string &widget_name) {
    return widget<Gtk::SpinButton>(widget_name)->get_value();
}
int GUI::sb_int(const std::string &widget_name) {
    return widget<Gtk::SpinButton>(widget_name)->get_value_as_int();
}

void GUI::start(int argc, char *argv[]) {
    if (gui_thread_.joinable())
        throw std::runtime_error("GUI thread can only be started once!");

    app_ = Gtk::Application::create(argc, argv, "ca.imaginary.eris.creativity",
            Gio::APPLICATION_NON_UNIQUE | Gio::APPLICATION_HANDLES_OPEN);
    builder_ = Gtk::Builder::create();

    std::list<std::string> datadirs;
    datadirs.push_front(DATADIR);
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

    gui_thread_ = std::thread(&GUI::thr_run, this);

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

void GUI::thr_run() {
    main_window_ = std::shared_ptr<Gtk::Window>(widget<Gtk::Window>("window1"));

    auto disable_on_load = {"fr_agents", "fr_save", "fr_run", "fr_sim"};
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

    // Set up a signal handler that is invoked if the application is started up with a filename
    app_->signal_open().connect_notify([this](const std::vector<Glib::RefPtr<Gio::File>> &files, const Glib::ustring&) -> void {
        if (files.size() != 1) throw std::runtime_error("Unable to open multiple files!");
        load_ = files[0]->get_path();
        app_->add_window(*main_window_);
        main_window_->show_all();
    }, true);

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
        }
        else {
            save_ = "";
        }
        widget<Gtk::Label>("lbl_save")->set_text(save_ == "" ? "(no file selected)" : save_);
    });

    widget<Gtk::Button>("btn_start")->signal_clicked().connect([this] {
        setupSim();
        runSim();
    });

    widget<Gtk::Button>("btn_init")->signal_clicked().connect([this] {
        setupSim();
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
    auto active = thrbox->get_active_row_number();
    // The gui setup has two entries: no threads, and 1 thread.  The latter is pointless (except for
    // eris debugging), so remove it, but first check to see if it was selected; if it was, we'll
    // update the default selection to the maximum number of threads; otherwise we'll leave it on no
    // threads.
    thrbox->remove_text(1);
    unsigned int last = 1;
    for (unsigned int i = 2, incr = 1; i <= std::thread::hardware_concurrency(); i += incr) {
        thrbox->append(std::to_string(i), std::to_string(i) + " threads");
        if (i >= 4*incr) incr *= 2;
        last = i;
    }
    if (last < std::thread::hardware_concurrency()) {
        // In case the processor has some weird number of max concurrency (such as 5,
        // 7, 9-11, 13-15, 17-23, etc.) that doesn't get added above, add it.
        std::string num = std::to_string(std::thread::hardware_concurrency());
        thrbox->append(num, num + " threads");
    }
    // If the GUI's default was the 1 thread option, update it to max threads
    if (active > 0 and std::thread::hardware_concurrency() > 1)
        thrbox->set_active_id(std::to_string(std::thread::hardware_concurrency()));

    thrbox->signal_changed().connect([this] {
        int threads = widget<Gtk::ComboBoxText>("combo_threads")->get_active_row_number();

        std::vector<Parameter> params;
        Parameter p;
        p.param = ParamType::threads;
        p.ul = threads;
        queueEvent(Event::Type::setup, std::move(params));
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
    graph_->signal_scroll_event().connect([this](GdkEventScroll *event) -> bool {
        if (event->direction == GdkScrollDirection::GDK_SCROLL_UP) {
            if (state_curr_ > 0) thr_set_state(state_curr_ - 1);
            return true;
        }
        if (event->direction == GdkScrollDirection::GDK_SCROLL_DOWN) {
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
    GUI_SETUP_VIS_SETTING(friendship)
    GUI_SETUP_VIS_SETTING(movement)
    GUI_SETUP_VIS_SETTING(author_live)
    GUI_SETUP_VIS_SETTING(author_dead)
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
    auto attribs_grid = widget<Gtk::Grid>("grid_sim_attribs");
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

    Glib::signal_idle().connect_once([this,&lock] {
            thread_running_ = true;
            lock.unlock();
            cv_.notify_all();
            if (not load_.empty()) loadSim();
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
        bk_models_[t]->set_sort_column(bk_models_[t]->columns.age, Gtk::SortType::SORT_ASCENDING);
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
        if (w.second.get_visible())
            w.second.refresh(state);
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
        for (auto &r : state->readers)
            rt.insert(std::make_pair(rt_point{r.second.position[0], r.second.position[1]}, r.second.id));
    }
    if (graph_->design.enabled.book_live or graph_->design.enabled.book_dead) {
        for (auto &b : state->books) {
            if (b.second.market ? graph_->design.enabled.book_live : graph_->design.enabled.book_dead)
                rt.insert(std::make_pair(rt_point{b.second.position[0], b.second.position[1]}, b.second.id));
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
        already_open->second.present();
    }
    else {
        std::shared_ptr<const State> state;
        {
            state = (*creativity_->storage().first)[state_curr_];
        }

        if (state->readers.count(member_id)) {
            info_windows_.emplace(std::piecewise_construct,
                    std::tuple<eris_id_t>{member_id},
                    std::tuple<decltype(state), decltype(main_window_), eris_id_t, std::function<void(eris_id_t)>>{
                        state, main_window_, member_id, std::bind(&GUI::thr_info_dialog, this, _1)});
        }
        else if (state->books.count(member_id)) {
            info_windows_.emplace(std::piecewise_construct,
                    std::tuple<eris_id_t>{member_id},
                    std::tuple<decltype(state), decltype(main_window_), eris_id_t>{state, main_window_, member_id});
        }
        else {
            throw std::out_of_range("thr_info_dialog: requested member id does not exist");
        }

    }
}

void GUI::thr_update_parameters() {
#define SET_SB(PARAMETER) widget<Gtk::SpinButton>("set_" #PARAMETER)->set_value(creativity_->parameters.PARAMETER)
    SET_SB(readers);
    SET_SB(density);
    SET_SB(book_distance_sd);
    SET_SB(book_quality_sd);
    SET_SB(reader_step_sd);
    SET_SB(reader_creation_shape);
    SET_SB(reader_creation_scale_min);
    SET_SB(reader_creation_scale_max);
    SET_SB(piracy_begins);
    SET_SB(income);
    SET_SB(cost_fixed);
    SET_SB(cost_unit);
    SET_SB(cost_piracy);
    SET_SB(prior_weight);
    SET_SB(prior_weight_piracy);
    widget<Gtk::SpinButton>("set_piracy_link_proportion")->set_value(creativity_->parameters.piracy_link_proportion * 100.0);
#define SET_INIT_SB(PARAMETER) widget<Gtk::SpinButton>("set_init_" #PARAMETER)->set_value(creativity_->parameters.initial.PARAMETER)
    SET_INIT_SB(prob_write);
    SET_INIT_SB(q_min);
    SET_INIT_SB(q_max);
    SET_INIT_SB(p_min);
    SET_INIT_SB(p_max);
    SET_INIT_SB(prob_keep);
    SET_INIT_SB(keep_price);
#undef SET_SB
#undef SET_INIT_SB
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

        // FIXME: This really doesn't look good, and doesn't work properly until the scale has a
        // range of at least piracy_begins (until then it appears at the beginning).  For now,
        // leave it off, but it might be worthwhile adding in later.
        //widget<Gtk::Scale>("scale_state")->add_mark(creativity_->parameters.piracy_begins, Gtk::POS_BOTTOM, "Piracy");

        // Disable spin buttons:
        for (auto &widg : {"set_readers", "set_density", "set_piracy_link_proportion",
                "set_book_distance_sd", "set_book_quality_sd", "set_piracy_begins", "set_income",
                "set_cost_fixed", "set_cost_unit", "set_cost_piracy", "set_init_prob_write",
                "set_init_q_min", "set_init_q_max", "set_init_p_min", "set_init_p_max",
                "set_init_prob_keep", "set_init_keep_price"})
            widget<Gtk::SpinButton>(widg)->set_sensitive(false);

        // Update the agent attributes settings with the settings from the creativity object.
        // These may be different from the defaults, particularly when opening a file.
        thr_update_parameters();

        widget<Gtk::Button>("btn_init")->set_sensitive(false);
        widget<Gtk::Button>("btn_start")->set_sensitive(false);
        widget<Gtk::Button>("btn_run")->set_visible(true);
        widget<Gtk::Button>("btn_run")->set_sensitive(true);
        widget<Gtk::Button>("btn_resume")->set_visible(false);
        widget<Gtk::Button>("btn_pause")->set_visible(false);
        widget<Gtk::Button>("btn_step")->set_sensitive(true);
    }
    else if (last_state.type == Signal::Type::running) {
        widget<Gtk::Notebook>("nb_tabs")->set_current_page(1);

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
        widget<Gtk::Label>("lbl_total")->set_text(std::to_string(last_progress.uls[1]));
        scale->set_fill_level(last_progress.uls[0]);
        scale->set_value(state_curr_);
/*        auto progress = widget<Gtk::ProgressBar>("progress_stage");
        progress->set_text(std::to_string(last_progress.uls[0]) + " / " + std::to_string(last_progress.uls[1]));
        progress->set_fraction((double) last_progress.uls[0] / (double) last_progress.uls[1]);*/
        std::ostringstream speed;
        speed << std::setw(7) << std::showpoint << last_progress.doubles[0];

        widget<Gtk::Label>("status_speed")->set_text(speed.str());
    }

    std::string errstr = errors.str();
    if (not errstr.empty()) {
        Gtk::MessageDialog dlg_error(*main_window_, "An error has occured!", false, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_CLOSE, true);
        dlg_error.set_secondary_text(errors.str());
        dlg_error.run();
        dlg_error.hide();
    }

    if (last_new_states.type == Signal::Type::new_states) {
        auto old_state_num = state_num_;
        state_num_ = creativity_->storage().first->size();

        if (state_num_ > old_state_num) {
            // - 1 here because we don't really count 0 as a stage
            widget<Gtk::Label>("lbl_total")->set_text(std::to_string(state_num_-1));

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
    params.push_back(load);

    queueEvent(Event::Type::setup, std::move(params));

    widget<Gtk::Label>("lbl_load")->set_text(load_);

    // Disable and hide the simulation controls (at the top of the window), since loading is read-only
    auto controls = widget<Gtk::Box>("box_controls");
    controls->set_sensitive(false);
    controls->set_visible(false);

    // Hide the speed indicator
    for (const auto &l : {"lbl_speed", "status_speed", "lbl_speed_units"})
        widget<Gtk::Label>(l)->set_visible(false);

    // Disable all the settings items that aren't modifiable anymore (basically everything except
    // the visualization settings):
    for (const auto &disable : {"fr_new", "fr_agents", "fr_sim", "fr_save", "fr_run"})
        widget<Gtk::Frame>(disable)->set_sensitive(false);

    // Completely disable and hide the expanders for settings that are totally irrelevant
    for (const auto &kill : {"ex_sim", "ex_save", "ex_run"}) {
        auto ex = widget<Gtk::Expander>(kill);
        ex->set_expanded(false);
        ex->set_sensitive(false);
        ex->set_visible(false);
    }
    // Open the visualization settings (it's the only thing that's actually useful on the page)
    widget<Gtk::Expander>("ex_vis")->set_expanded(true);
}

// Start a new simulation:
void GUI::setupSim() {
    // Set the parameters directly
    auto &set = creativity_->set();
    set.readers = lround(sb("set_readers"));
    set.use_density = true;
    set.piracy_begins = lround(sb("set_piracy_begins"));
#define COPY_SB_D(PARAMETER) set.PARAMETER = sb("set_"#PARAMETER)
    COPY_SB_D(density);
    COPY_SB_D(piracy_link_proportion) * 0.01; // From percentage
    COPY_SB_D(book_distance_sd);
    COPY_SB_D(book_quality_sd);
    COPY_SB_D(reader_step_sd);
    COPY_SB_D(reader_creation_shape);
    COPY_SB_D(reader_creation_scale_min);
    COPY_SB_D(reader_creation_scale_max);
    COPY_SB_D(income);
    COPY_SB_D(cost_fixed);
    COPY_SB_D(cost_unit);
    COPY_SB_D(cost_piracy);
    COPY_SB_D(prior_weight);
    COPY_SB_D(prior_weight_piracy);
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
    std::vector<Parameter> params;
    Parameter p;
    p.param = ParamType::seed; p.ul = std::stoul(widget<Gtk::Entry>("set_seed")->get_text()); params.push_back(p);
    int threads = widget<Gtk::ComboBoxText>("combo_threads")->get_active_row_number();
    if (threads < 0) threads = 0; // -1 means no item selected (shouldn't be possible, but just in case)
    p.param = ParamType::threads; p.ul = (unsigned long)threads; params.push_back(p);

    if (save_ != "") {
        p.param = ParamType::save_as; p.ptr = &save_; params.push_back(p);
    }

    queueEvent(Event::Type::setup, std::move(params));
}

void GUI::runSim() {
    queueEvent(Event::Type::run, sb_int("set_periods"));
}

void GUI::newStates(unsigned long switch_to) {
    queueSignal((switch_to == (unsigned long)-1)
            ? Signal::Type::new_states
            : Signal(Signal::Type::new_states, switch_to));
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
GUI::Event::Event(Event::Type t, std::vector<Parameter> &&p)
    : type{t}, parameters{std::move(p)} {}
GUI::Event::Event(Event::Type t, unsigned long ulval)
    : type{t}, ul{ulval} {}
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
            case Event::Type::setup:
                if (on_setup_ and not event.parameters.empty()) {
                    on_setup_({ParamType::begin, {0}});
                    for (auto &p : event.parameters) {
                        try {
                            on_setup_(p);
                        }
                        catch (std::exception &e) {
                            setup_errors.push_back(e.what());
                        }
                    }
                    on_setup_({ setup_errors.empty() ? ParamType::finished : ParamType::erred, {0} });
                }
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
            default: // Ignore anything else
                break;
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
