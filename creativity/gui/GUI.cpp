#include "creativity/gui/GUI.hpp"
#include "creativity/config.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace eris;

namespace creativity { namespace gui {

GUI::GUI(std::shared_ptr<Simulation> sim,
        std::function<void(Parameter)> setup,
        std::function<void(unsigned int count)> run,
        std::function<void()> stop,
        std::function<void()> resume,
        std::function<void()> step,
        std::function<void()> quit)
    :
        sim_{sim},
        on_setup_{std::move(setup)},
        on_run_{std::move(run)},
        on_stop_{std::move(stop)},
        on_resume_{std::move(resume)},
        on_step_{std::move(step)},
        on_quit_{std::move(quit)}
{}

GUI::Exception::Exception(const std::string &what) : std::runtime_error(what) {}

template<class T>
T* GUI::widget(const std::string &widget_name) {
    T* widget;
    builder_->get_widget(widget_name, widget);
    return widget;
}

double GUI::sb(const std::string &widget_name) {
    return widget<Gtk::SpinButton>(widget_name)->get_value();
}
int GUI::sb_int(const std::string &widget_name) {
    return widget<Gtk::SpinButton>(widget_name)->get_value_as_int();
}

void GUI::start(int argc, char *argv[]) {
    if (gui_thread_.joinable())
        throw std::runtime_error("GUI thread can only be started once!");

    app_ = Gtk::Application::create(argc, argv, "ca.imaginary.test.cairo-drawing",
            Gio::APPLICATION_NON_UNIQUE);
    builder_ = Gtk::Builder::create();

    std::string datadir(DATADIR);
    char *envdatadir = getenv("CREATIVITY_DATADIR");
    if (envdatadir) {
        std::string stddatadir(envdatadir);
        if (stddatadir != "") datadir = stddatadir;
    }

    builder_->add_from_file(datadir + "/gui.glade"); // May throw exception

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

    auto thrbox = widget<Gtk::ComboBox>("combo_threads");
    auto thrlist = Glib::RefPtr<Gtk::ListStore>::cast_static(thrbox->get_model());
    for (unsigned int i = 2; i <= std::thread::hardware_concurrency(); i++) {
        auto iter = thrlist->append();
        iter->set_value(0, i);
        iter->set_value(1, std::to_string(i) + " threads");
    }
    //thrbox->set_active(*(thrlist->children().rbegin()));
    thrbox->set_active(*(thrlist->children().begin()));

    Gtk::Viewport *vis;
    builder_->get_widget("view_vis", vis);

    // Lock access to everything until we're set up.  We'll unlock and notify on cv_ once the
    // mainloop is running.
    std::unique_lock<std::mutex> lock{mutex_};

    graph_ = std::unique_ptr<GraphArea>(new GraphArea{10., 10., -10., -10., sim_, *this});
    graph_->add_events(Gdk::EventMask::BUTTON_PRESS_MASK);

    hand_ = Gdk::Cursor::create(main_window_->get_display(), Gdk::HAND1);

    motion_handler_ = [this](GdkEventMotion *event) -> bool {
        // When the mouse moves over the image, change the cursor to a hand any time the cursor
        // is within 5 pixels of a graph item (book, reader, etc.)
        rt_point clicked(event->x, event->y);

        std::vector<rt_val> nearest;
        rtree_.query(boost::geometry::index::nearest(clicked, 1), std::back_inserter(nearest));
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
        std::vector<rt_val> nearest;
        rtree_.query(boost::geometry::index::nearest(clicked, 1), std::back_inserter(nearest));
        if (nearest.empty() or boost::geometry::distance(nearest[0].first, clicked) > 5)
            return false;

        thr_info_dialog(nearest[0].second);
        return true;
    });

    rdr_model_ = ReaderStore::create(sim_);
    rdr_tree_ = std::unique_ptr<Gtk::TreeView>(new Gtk::TreeView);
    rdr_tree_->set_model(rdr_model_);
    rdr_model_->appendColumnsTo(*rdr_tree_);
    rdr_tree_->set_fixed_height_mode(true);
    rdr_model_->set_sort_column(rdr_model_->columns.id, Gtk::SortType::SORT_ASCENDING);
    rdr_tree_->signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
            thr_info_dialog(rdr_model_->reader(path));
    });

    ERIS_DBG("");
    widget<Gtk::ScrolledWindow>("win_rdr")->add(*rdr_tree_);
    ERIS_DBG("");
    rdr_tree_->show();
    ERIS_DBG("");

    bk_model_ = BookStore::create(sim_);
    bk_tree_ = std::unique_ptr<Gtk::TreeView>(new Gtk::TreeView);
    bk_tree_->set_model(bk_model_);
    bk_model_->appendColumnsTo(*bk_tree_);
    bk_tree_->set_fixed_height_mode(true);
    bk_model_->set_sort_column(bk_model_->columns.id, Gtk::SortType::SORT_DESCENDING);
    bk_tree_->signal_row_activated().connect([this] (const Gtk::TreeModel::Path &path, Gtk::TreeViewColumn*) -> void {
            thr_info_dialog(bk_model_->book(path));
    });

    ERIS_DBG("");
    widget<Gtk::ScrolledWindow>("win_bk")->add(*bk_tree_);
    ERIS_DBG("");
    bk_tree_->show();
    ERIS_DBG("");


    dispatcher_ = std::unique_ptr<Glib::Dispatcher>(new Glib::Dispatcher);
    dispatcher_->connect([this] { thr_signal(); });

    vis->add(*graph_);
    graph_->show();

    Glib::signal_idle().connect_once([this,&lock] {
            thread_running_ = true;
            lock.unlock();
            cv_.notify_all();
            });

    app_->run(*main_window_);

    queueEvent(Event::Type::quit);
    lock.lock();
    dispatcher_.reset();
}

void GUI::thr_update_readers() {
    rdr_model_->resync();
    /*
    ERIS_DBG("starting");

    auto new_readers = sim_->agents<Reader>([this] (const Reader &r) -> bool { return r.id() > rdr_biggest_id_; });

    // Save the current sort column and order
    int sort_col;
    Gtk::SortType sort_type;
    rdr_model_->get_sort_column_id(sort_col, sort_type);
    
    // Freeze tree to prevent each row modification triggering notifications
    rdr_tree_->freeze_child_notify();

    // Detach model from the tree
    rdr_tree_->unset_model();

    // Turn off sorting
    rdr_model_->set_sort_column(Gtk::TreeSortable::DEFAULT_UNSORTED_COLUMN_ID, sort_type);

    ERIS_DBG("Sorting off; updating");
    // Updating existing rows, removing if any reader isn't in the simulation anymore
    for (auto row : rdr_model_->children()) {
        eris_id_t id = row[rdr_cols_->id];
        if (not sim_->hasAgent(id)) {
            ERIS_DBG("FIXME: remove agent!");
            throw std::runtime_error("FIXME: remove agent!");
        }
        else {
            rdr_cols_->updateRow(row, sim_->agent<Reader>(id));
        }
    }


    ERIS_DBG("adding new ones");

    // Re-add all readers (whose info has almost always changed), potentially including new ones
    for (auto &r : sim_->agents<Reader>()) {
        rdr_cols_->appendRow(rdr_model_, r);
    }
    ERIS_DBG("resorting");

    // Turn sorting back on, which will resort everything again
    //rdr_model_->set_sort_column(sort_col, sort_type);

    ERIS_DBG("set_model");
    // Reassociate the model with the tree, and thaw the tree
    rdr_tree_->set_model(rdr_model_);
    ERIS_DBG("thawing");
    rdr_tree_->thaw_child_notify();

    ERIS_DBG("finished");
    */
}

void GUI::thr_update_books() {
    bk_model_->resync();
}

void GUI::thr_info_dialog(SharedMember<Member> member) {
    auto already_open = info_windows_.find(member->id());
    if (already_open != info_windows_.end()) {
        // If the window is already active, just present it again
        already_open->second.present();
    }
    else {
        // Otherwise we need to create a new window
        SharedMember<Reader> reader;
        SharedMember<Book> book;
        try {
            reader = member;
        }
        catch (std::bad_cast&) {
            // Not a reader: must be a book (if this throws again, don't catch it)
            book = member;
        }

        if (reader) {
            info_windows_.emplace(std::piecewise_construct,
                    std::tuple<eris_id_t>{reader},
                    std::tuple<decltype(reader), decltype(main_window_)>{reader, main_window_});
        }
        else {
            info_windows_.emplace(std::piecewise_construct,
                    std::tuple<eris_id_t>{book},
                    std::tuple<decltype(book), decltype(main_window_)>{book, main_window_});
        }
    }
}

void GUI::thr_signal() {
    std::unique_lock<std::mutex> lock(mutex_);
    // Keep track of the *last* state change we see, and whether or not we see a redraw (we want to
    // skip everything except the last).
    std::stringstream errors;
    Signal last_state{Signal::Type::quit}, last_progress{Signal::Type::quit}, last_redraw{Signal::Type::quit};
    for (Signal &s : signal_queue_) {
        switch (s.type) {
            case Signal::Type::quit:
                app_->quit();
                return;
                break;
            case Signal::Type::redraw:
                last_redraw = s;
                break;
            case Signal::Type::initialized:
            case Signal::Type::running:
            case Signal::Type::stopped:
                last_state = s;
                break;
            case Signal::Type::progress:
                last_progress = s;
                break;
            case Signal::Type::error:
                ERIS_DBG("Thr received error signal");
                errors << s.message << "\n";
                break;
        }
    }
    signal_queue_.clear();
    lock.unlock();

    if (last_state.type == Signal::Type::initialized) {
        widget<Gtk::Notebook>("nb_tabs")->set_current_page(1);

        // Disable spin buttons:
        for (auto &widg : {"set_dimensions", "set_readers", "set_book_sd", "set_quality_draw_sd", "set_cost_fixed", "set_cost_unit"})
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

        // Disconnect motion handler while running
        if (motion_handler_conn_) motion_handler_conn_.disconnect();
    }
    else if (last_state.type == Signal::Type::stopped) {
        // Turn off pause and turn on either play or resume buttons
        widget<Gtk::Button>("btn_run")->set_visible(!last_state.boolean);
        widget<Gtk::Button>("btn_resume")->set_visible(last_state.boolean);
        widget<Gtk::Button>("btn_pause")->set_visible(false);
        widget<Gtk::Button>("btn_step")->set_sensitive(true);

        // Stick all points into an rtree so that we can quickly find the nearest one (needed, in
        // particular, for fast mouseovers).
        rtree_.clear();
        auto g2c = graph_->graph_to_canvas();
        for (auto &r : sim_->agents<Reader>()) {
            double x = r->position()[0], y = r->position()[1];
            g2c.transform_point(x, y);
            rtree_.insert(std::make_pair(rt_point{x, y}, SharedMember<Member>{r}));
        }
        for (auto &b : sim_->goods<Book>()) {
            double x = b->position()[0], y = b->position()[1];
            g2c.transform_point(x, y);
            rtree_.insert(std::make_pair(rt_point{x, y}, SharedMember<Member>{b}));
        }

        // Reenable mouseover handler
        if (not motion_handler_conn_)
            motion_handler_conn_ = graph_->signal_motion_notify_event().connect(motion_handler_);
    }
    if (last_progress.type == Signal::Type::progress) {
        // .uls contains: t, end, readers, books
        // .doubles has: speed
        auto progress = widget<Gtk::ProgressBar>("progress_stage");
        progress->set_text(std::to_string(last_progress.uls[0]) + " / " + std::to_string(last_progress.uls[1]));
        progress->set_fraction((double) last_progress.uls[0] / (double) last_progress.uls[1]);
        widget<Gtk::Label>("status_readers")->set_text(std::to_string(last_progress.uls[2]));
        widget<Gtk::Label>("status_books")->set_text(std::to_string(last_progress.uls[3]));
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

    if (last_redraw.type == Signal::Type::redraw) {
        thr_update_readers();
        thr_update_books();
        if (graph_->get_is_drawable()) {
            graph_->queue_draw();
        }
        else {
            ERIS_DBG("Fake redraw (not currently drawable) completed.");
            // Not currently drawable (perhaps not on visualization tab), so send back a fake redraw
            // event in case the caller is going to wait for a redraw to complete.
            queueEvent(Event::Type::redraw_complete);
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

void GUI::setupSim() {
    std::vector<Parameter> params;
    Parameter p;
    p.param = ParamType::readers; p.ul = sb_int("set_readers"); params.push_back(p);
    p.param = ParamType::dimensions; p.ul = sb_int("set_dimensions"); params.push_back(p);
    p.param = ParamType::book_sd; p.dbl = sb("set_book_sd"); params.push_back(p);
    p.param = ParamType::quality_draw_sd; p.dbl = sb("set_quality_draw_sd"); params.push_back(p);
    p.param = ParamType::cost_fixed; p.dbl = sb("set_cost_fixed"); params.push_back(p);
    p.param = ParamType::cost_unit; p.dbl = sb("set_cost_unit"); params.push_back(p);
    p.param = ParamType::speed_limit; p.dbl = sb("set_speed"); params.push_back(p);
    p.param = ParamType::redraw; p.dur_ms = std::chrono::milliseconds{sb_int("set_redraw")}; params.push_back(p);
    unsigned int threads;
    widget<Gtk::ComboBox>("combo_threads")->get_active()->get_value(0, threads);
    p.param = ParamType::threads; p.ul = (unsigned long)threads; params.push_back(p);

    queueEvent(Event::Type::setup, std::move(params));
}

void GUI::runSim() {
    queueEvent(Event::Type::run, sb_int("set_periods"));
}

void GUI::redraw(bool sync) {
    queueSignal(Signal::Type::redraw);
    if (sync) waitForEvent(Event::Type::redraw_complete);
}

void GUI::initialized() { queueSignal(Signal::Type::initialized); }

void GUI::running() { queueSignal(Signal::Type::running); }

void GUI::progress(const unsigned long &end, const double &speed) {
    Signal s{Signal::Type::progress};
    s.uls = {sim_->t(), end, sim_->countAgents<Reader>(), sim_->countGoods<Book>()};
    s.doubles = {speed};
    queueSignal(std::move(s));
}

void GUI::stopped(bool manual) { queueSignal({Signal::Type::stopped, manual}); }

void GUI::error(std::string message) { queueSignal({ Signal::Type::error, message }); }

GUI::Event::Event() : none{true} {}
GUI::Event::Event(Event::Type t) : type{t} {}
GUI::Event::Event(Event::Type t, std::vector<Parameter> &&p)
    : type{t}, parameters{std::move(p)} {}
GUI::Event::Event(Event::Type t, const unsigned long &ulval)
    : type{t}, ul{ulval} {}
GUI::Event::operator bool() { return not none; }

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
