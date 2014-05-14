#include "creativity/GUI.hpp"
#include "creativity/config.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include "creativity/BookMarket.hpp"
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace eris;

namespace creativity {

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
    auto *main = widget<Gtk::Window>("window1");

    auto *run = widget<Gtk::Button>("btn_start");

    run->signal_clicked().connect([this] {
        ERIS_DBG("run clicked");
        std::vector<Parameter> params;
        Parameter p;
        p.param = ParamType::readers; p.ul = sb_int("set_readers"); params.push_back(p);
        p.param = ParamType::dimensions; p.ul = sb_int("set_dimensions"); params.push_back(p);
        p.param = ParamType::book_sd; p.dbl = sb("set_book_sd"); params.push_back(p);
        p.param = ParamType::quality_draw_sd; p.dbl = sb("set_quality_draw_sd"); params.push_back(p);
        p.param = ParamType::cost_fixed; p.dbl = sb("set_cost_fixed"); params.push_back(p);
        p.param = ParamType::cost_unit; p.dbl = sb("set_cost_unit"); params.push_back(p);
        p.param = ParamType::speed_limit; p.dur_ms = std::chrono::milliseconds{sb_int("set_speed")}; params.push_back(p);
        p.param = ParamType::redraw; p.dur_ms = std::chrono::milliseconds{sb_int("set_redraw")}; params.push_back(p);
        unsigned int threads;
        widget<Gtk::ComboBox>("combo_threads")->get_active()->get_value(0, threads);
        p.param = ParamType::threads; p.ul = (unsigned long)threads; params.push_back(p);

        queueEvent(Event::Type::setup, std::move(params));
        queueEvent(Event::Type::run, sb_int("set_periods"));
        ERIS_DBG("run events sent");
    });

    widget<Gtk::Button>("btn_run")->signal_clicked().connect([this] {
        queueEvent(Event::Type::run, sb_int("set_periods"));
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

    graph_ = std::unique_ptr<GUIGraphArea>(new GUIGraphArea{10., 10., -10., -10., sim_, *this});
    main->add_events(Gdk::EventMask::BUTTON_PRESS_MASK);

    // FIXME: change mouse to clickable when over readers/books?

    graph_->signal_button_press_event().connect([this](GdkEventButton *event) -> bool {
        // Calculate the *graph* position where the user clicked
        auto c2g = graph_->graph_to_canvas();
        c2g.invert();

        double x = event->x, y = event->y;
        c2g.translate(x, y);

        Position pos{{x, y}};
        // Figure out the closest thing (Reader or Book) to where the click happened
        SharedMember<Reader> reader;
        SharedMember<Book> book;

        for (auto &r : sim_->agents<Reader>())
            if (not reader or r->distance(pos) < reader->distance(pos))
                reader = r;
        for (auto &r : sim_->goods<Book>())
            if (not book or r->distance(pos) < book->distance(pos))
                book = r;

        double dist = std::min(reader->distance(pos), book->distance(pos));

        // If not within 5 pixels of the closest thing, ignore.
        if (dist > 5) return false;

        if (reader->distance(pos) <= book->distance(pos)) {
            // If a reader, open up and fill out the reader info dialog
            auto dlg_reader = widget<Gtk::Dialog>("dlg_readerinfo");
            widget<Gtk::Label>("vrd_id")->set_text(std::to_string(reader->id()));
            std::ostringstream pos;
            pos << "(" << std::setw(7) << std::showpoint << reader->position()[0] << "," << reader->position()[1] << ")";
            widget<Gtk::Label>("vrd_pos")->set_text(pos.str());
            widget<Gtk::Label>("vrd_u")->set_text(std::to_string(reader->u()));
            widget<Gtk::Label>("vrd_ulife")->set_text(std::to_string(reader->uLifetime()));
            widget<Gtk::Label>("vrd_libsize")->set_text(std::to_string(reader->library().size()));
            widget<Gtk::Label>("vrd_libnew")->set_text(std::to_string(reader->newBooks().size()));
            widget<Gtk::Label>("vrd_books")->set_text(std::to_string(reader->wrote().size()));
            auto last_book = sim_->good<Book>(reader->wrote().back());
            widget<Gtk::Label>("vrd_lastbookage")->set_text(std::to_string(last_book->age()));
            dlg_reader->run();
            dlg_reader->hide();
        }
        else {
            // If a book, open up and fill out the book info dialog
            auto dlg_book = widget<Gtk::Dialog>("dlg_bookinfo");
            widget<Gtk::Label>("vbk_id")->set_text(std::to_string(book->id()));
            std::ostringstream pos;
            pos << "(" << std::setw(7) << std::showpoint << book->position()[0] << "," << book->position()[1] << ")";
            widget<Gtk::Label>("vbk_pos")->set_text(pos.str());
            widget<Gtk::Label>("vbk_mkt")->set_text(book->hasMarket()
                    ? std::to_string(book->market()->id())
                    : "(not on market)"
            );
            widget<Gtk::Label>("vbk_p")->set_text(book->hasMarket()
                    ? std::to_string(book->market()->price())
                    : "(not on market)"
            );
            widget<Gtk::Label>("vbk_age")->set_text(std::to_string(book->age()));
            widget<Gtk::Label>("vbk_sales")->set_text(std::to_string(book->lifeSales()));
            unsigned long copies = 0;
            for (auto &r : sim_->agents<Reader>()) {
                if (r->library().count(book))
                    copies++;
            }
            widget<Gtk::Label>("vbk_copies")->set_text(std::to_string(copies));
            widget<Gtk::Label>("vbk_author")->set_text(std::to_string(book->author()));
            dlg_book->run();
            dlg_book->hide();
        }
        return true;
    });

    dispatcher_ = std::unique_ptr<Glib::Dispatcher>(new Glib::Dispatcher);
    dispatcher_->connect([this] { thr_signal(); });

    vis->add(*graph_);
    graph_->show();

    Glib::signal_idle().connect_once([this,&lock] {
            thread_running_ = true;
            lock.unlock();
            cv_.notify_all();
            });

    app_->run(*main);

    queueEvent(Event::Type::quit);
    lock.lock();
    dispatcher_.reset();
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
            case Signal::Type::running:
            case Signal::Type::stopped:
                last_state = s;
                break;
            case Signal::Type::progress:
                last_progress = s;
                break;
            case Signal::Type::error:
                std::cerr << "Thr received error signal\n";
                errors << s.message << "\n";
                break;
        }
    }
    signal_queue_.clear();
    lock.unlock();

    if (last_state.type == Signal::Type::running) {
        std::cerr << "running signal received\n";
        widget<Gtk::Notebook>("nb_tabs")->set_current_page(1);

        // Disable spin buttons:
        for (auto &widg : {"set_dimensions", "set_readers", "set_book_sd", "set_quality_draw_sd", "set_cost_fixed", "set_cost_unit", "set_speed"})
            widget<Gtk::SpinButton>(widg)->set_sensitive(false);

        widget<Gtk::Button>("btn_start")->set_sensitive(false);
        widget<Gtk::Button>("btn_run")->set_visible(false);
        widget<Gtk::Button>("btn_run")->set_sensitive(true);
        widget<Gtk::Button>("btn_resume")->set_visible(false);
        widget<Gtk::Button>("btn_pause")->set_visible(true);
        widget<Gtk::Button>("btn_step")->set_sensitive(false);
    }
    else if (last_state.type == Signal::Type::stopped) {
        // Turn off pause and turn on either play or resume buttons
        widget<Gtk::Button>("btn_run")->set_visible(last_state.boolean);
        widget<Gtk::Button>("btn_resume")->set_visible(!last_state.boolean);
        widget<Gtk::Button>("btn_pause")->set_visible(false);
        widget<Gtk::Button>("btn_step")->set_sensitive(true);
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
        std::cerr << "queueing a redraw\n";
        if (graph_->get_is_drawable()) {
            graph_->queue_draw();
        }
        else {
            std::cerr << "...but not visible, so skipping\n";
            // Not currently drawable (perhaps not on visualization tab), so send back a fake redraw
            // event in case the caller is going to wait for a redraw to complete.
            queueEvent(Event::Type::redraw);
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

void GUI::redraw(bool sync) {
    queueSignal(Signal::Type::redraw);
    if (sync) waitForEvent(Event::Type::redraw);
}

void GUI::running() { queueSignal(Signal::Type::running); }

void GUI::progress(const unsigned long &end, const double &speed) {
    Signal s{Signal::Type::progress};
    s.uls = {sim_->t(), end, sim_->countAgents<Reader>(), sim_->countGoods<Book>()};
    s.doubles = {speed};
    queueSignal(std::move(s));
}

void GUI::stopped(bool done) { queueSignal({Signal::Type::stopped, done}); }

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
            std::cerr << "saw event\n";
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

}
