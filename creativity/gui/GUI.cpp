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

    // When the "load" radio button is activated, disable all the new simulation parameter fields
    widget<Gtk::RadioButton>("radio_load")->signal_clicked().connect([this] {
        widget<Gtk::Frame>("fr_agents")->set_sensitive(false);
        widget<Gtk::Frame>("fr_data")->set_sensitive(false);
        widget<Gtk::Frame>("fr_run")->set_sensitive(false);
    });
    // Undo the above when "New simulation" is activated
    widget<Gtk::RadioButton>("radio_new")->signal_clicked().connect([this] {
        widget<Gtk::Frame>("fr_agents")->set_sensitive(true);
        widget<Gtk::Frame>("fr_data")->set_sensitive(true);
        widget<Gtk::Frame>("fr_run")->set_sensitive(true);
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
    int default_threads = 0;
    for (unsigned int i = 2; i <= std::thread::hardware_concurrency(); i++) {
        thrbox->append(std::to_string(i) + " threads");
        default_threads = i;
    }
    thrbox->set_active(default_threads);
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

    widget<Gtk::ComboBox>("combo_state")->signal_changed().connect([this] {
        thr_set_state(widget<Gtk::ComboBox>("combo_state")->get_active_row_number());
    });

    widget<Gtk::Button>("btn_prev")->signal_clicked().connect([this] {
        if (state_curr_ > 0) thr_set_state(state_curr_ - 1);
    });
    widget<Gtk::Button>("btn_next")->signal_clicked().connect([this] {
        if (state_curr_ + 1 < state_num_) thr_set_state(state_curr_ + 1);
    });
    widget<Gtk::Button>("btn_last")->signal_clicked().connect([this] {
        thr_set_state(state_num_ - 1);
    });

    Gtk::Viewport *vis;
    builder_->get_widget("view_vis", vis);

    // Lock access to everything until we're set up.  We'll unlock and notify on cv_ once the
    // mainloop is running.
    std::unique_lock<std::mutex> lock{mutex_};

    graph_ = std::unique_ptr<GraphArea>(new GraphArea(*this));
    graph_->add_events(Gdk::EventMask::BUTTON_PRESS_MASK);

    hand_ = Gdk::Cursor::create(main_window_->get_display(), Gdk::HAND1);

    motion_handler_ = [this](GdkEventMotion *event) -> bool {
        rt_point clicked(event->x, event->y);

        // When the mouse moves over the image, change the cursor to a hand any time the cursor
        // is within 5 pixels of a graph item (book, reader, etc.)
        auto nearest = thr_nearest(clicked);

        if (not nearest.empty() and boost::geometry::distance(nearest[0].first, clicked) <= 5)
            main_window_->get_window()->set_cursor(hand_);
        else
            main_window_->get_window()->set_cursor();

        return false;
    };
    graph_->add_events(Gdk::POINTER_MOTION_MASK);
    motion_handler_conn_ = graph_->signal_motion_notify_event().connect(motion_handler_);

    graph_->signal_button_press_event().connect([this](GdkEventButton *event) -> bool {
        rt_point clicked(event->x, event->y);

        // Figure out the closest Reader and Book on the canvas relative to where the click happened
        auto nearest = thr_nearest(clicked);

        if (nearest.empty() or boost::geometry::distance(nearest[0].first, clicked) > 5)
            return false;

        thr_info_dialog(nearest[0].second);
        return true;
    });

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

    if (rtrees_[t].empty()) {
        // Stick all points into an rtree so that we can quickly find the nearest one (needed, in
        // particular, for fast mouseovers).
        auto &rt = rtrees_[t];
        auto g2c = graph_->graph_to_canvas();
        for (auto &r : state->readers) {
            double x = r.second.position[0], y = r.second.position[1];
            g2c.transform_point(x, y);
            rt.insert(std::make_pair(rt_point{x, y}, r.second.id));
        }
        for (auto &b : state->books) {
            double x = b.second.position[0], y = b.second.position[1];
            g2c.transform_point(x, y);
            rt.insert(std::make_pair(rt_point{x, y}, b.second.id));
        }
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
    widget<Gtk::ComboBox>("combo_state")->set_active(t);
    widget<Gtk::Label>("lbl_tab_agents")->set_text("Agents (" + std::to_string(state->readers.size()) + ")");
    widget<Gtk::Label>("lbl_tab_books")->set_text("Books (" + std::to_string(state->books.size()) + ")");

    widget<Gtk::Button>("btn_prev")->set_sensitive(state_curr_ > 0);
    widget<Gtk::Button>("btn_next")->set_sensitive(state_curr_ + 1 < state_num_);
    widget<Gtk::Button>("btn_last")->set_sensitive(state_curr_ + 1 < state_num_);

    if (graph_->get_is_drawable()) graph_->queue_draw();
}

std::vector<GUI::rt_val> GUI::thr_nearest(const rt_point &point, int n) {
    std::vector<rt_val> nearest;
    if (state_num_ > 0)
        rtrees_[state_curr_].query(boost::geometry::index::nearest(point, n), std::back_inserter(nearest));
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
        for (auto &widg : {"set_dimensions", "set_readers", "set_density", "set_book_sd", "set_quality_draw_sd", "set_cost_fixed", "set_cost_unit"})
            widget<Gtk::SpinButton>(widg)->set_sensitive(false);

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
        auto progress = widget<Gtk::ProgressBar>("progress_stage");
        progress->set_text(std::to_string(last_progress.uls[0]) + " / " + std::to_string(last_progress.uls[1]));
        progress->set_fraction((double) last_progress.uls[0] / (double) last_progress.uls[1]);
        std::ostringstream speed;
        speed << std::setw(7) << std::showpoint << last_progress.doubles[0];

        widget<Gtk::Label>("status_speed")->set_text(speed.str());
    }

    std::string errstr = errors.str();
    if (not errstr.empty()) {
        auto dlg_error = widget<Gtk::MessageDialog>("dlg_error");
        dlg_error->set_secondary_text(errors.str());
        dlg_error->run();
        dlg_error->hide();
    }

    if (last_new_states.type == Signal::Type::new_states) {
        auto old_state_num = state_num_;
        state_num_ = creativity_->storage().first->size();

        if (state_num_ > old_state_num) {
            // - 1 here because we don't really count 0 as a stage
            widget<Gtk::Label>("lbl_total")->set_text(std::to_string(state_num_-1));

            // In case the next/last buttons are disable, reactivate them
            widget<Gtk::Button>("btn_next")->set_sensitive(true);
            widget<Gtk::Button>("btn_last")->set_sensitive(true);

            auto periods = widget<Gtk::ComboBoxText>("combo_state");
            for (unsigned long t = old_state_num; t < state_num_; t++) {
                periods->append(std::to_string(t));
            }

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

void GUI::loadSim() {
    std::vector<Parameter> params;
    Parameter p;
    int threads = widget<Gtk::ComboBoxText>("combo_threads")->get_active_row_number();
    if (threads < 0) threads = 0; // -1 means no item selected (shouldn't be possible, but just in case)
    p.param = ParamType::threads; p.ul = (unsigned long)threads; params.push_back(p);

    if (load_ == "") throw std::runtime_error("loadSim() called without load_ set to file to load");

    p.param = ParamType::load; p.ptr = &load_; params.push_back(p);

    queueEvent(Event::Type::setup, std::move(params));

    // Disable all the simulation controls, since loading is read-only
    widget<Gtk::Box>("box_settings")->set_sensitive(false);
    widget<Gtk::Box>("box_controls")->set_sensitive(false);
}

void GUI::setupSim() {
    std::vector<Parameter> params;
    Parameter p;
    p.param = ParamType::readers; p.ul = sb_int("set_readers"); params.push_back(p);
    p.param = ParamType::dimensions; p.ul = sb_int("set_dimensions"); params.push_back(p);
    p.param = ParamType::density; p.dbl = sb("set_density"); params.push_back(p);
    p.param = ParamType::book_sd; p.dbl = sb("set_book_sd"); params.push_back(p);
    p.param = ParamType::quality_draw_sd; p.dbl = sb("set_quality_draw_sd"); params.push_back(p);
    p.param = ParamType::cost_fixed; p.dbl = sb("set_cost_fixed"); params.push_back(p);
    p.param = ParamType::cost_unit; p.dbl = sb("set_cost_unit"); params.push_back(p);
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
        queueSignal({ Signal::Type::error, e.what() });
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
