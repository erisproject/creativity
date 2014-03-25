#include "creativity/GUI.hpp"
#include "creativity/config.hpp"
#include "creativity/Reader.hpp"
#include "creativity/Book.hpp"
#include <iostream>
#include <cstdlib>
#include <iomanip>

namespace creativity {

GUI::GUI(std::shared_ptr<eris::Simulation> sim,
        std::function<void(Parameters)> setup,
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
        Parameters p;
        p.readers = sb_int("set_readers");
        p.dimensions = sb_int("set_dimensions");
        p.prob_writer = sb("set_create_prob");
        p.book_sd = sb("set_sd");
        p.speed_limit = std::chrono::milliseconds{sb_int("set_speed")};
        p.redraw = std::chrono::milliseconds{sb_int("set_redraw")};

        for (auto &widg : {"set_readers", "set_dimensions", "set_create_prob", "set_sd", "set_speed"})
            widget<Gtk::SpinButton>(widg)->set_sensitive(false);

        queueEvent(Event::Type::setup, std::move(p));
        queueEvent(Event::Type::run, sb_int("set_periods"));
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

    Gtk::Viewport *vis;
    builder_->get_widget("view_vis", vis);

    // Lock access to everything until we're set up.  We'll unlock and notify on cv_ once the
    // mainloop is running.
    std::unique_lock<std::mutex> lock{mutex_};

    graph_ = std::unique_ptr<GUIGraphArea>(new GUIGraphArea{11., 11., -11., -11., *this, sim_});

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
                // FIXME:
                // - Show error popup
                break;
        }
    }
    signal_queue_.clear();
    lock.unlock();

    if (last_state.type == Signal::Type::running) {
        widget<Gtk::Notebook>("nb_tabs")->set_current_page(1);

        // - Replace Play button with Pause button
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
        auto progress = widget<Gtk::ProgressBar>("progressbar1");
        progress->set_text(std::to_string(last_progress.uls[0]) + " / " + std::to_string(last_progress.uls[1]));
        progress->set_fraction((double) last_progress.uls[0] / (double) last_progress.uls[1]);
        widget<Gtk::Label>("status_readers")->set_text(std::to_string(last_progress.uls[2]));
        widget<Gtk::Label>("status_books")->set_text(std::to_string(last_progress.uls[3]));
        std::ostringstream speed;
        speed << std::setw(7) << std::showpoint << last_progress.doubles[0];

        widget<Gtk::Label>("status_speed")->set_text(speed.str());
    }
    if (last_redraw.type == Signal::Type::redraw) {
        graph_->queue_draw();
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

void GUI::redraw() { queueSignal(Signal::Type::redraw); }

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
GUI::Event::Event(Event::Type t, Parameters &&p)
    : type{t}, parameters{std::make_shared<Parameters>(std::move(p))} {}
GUI::Event::Event(Event::Type t, const unsigned long &ulval)
    : type{t}, ul{ulval} {}
GUI::Event::operator bool() { return not none; }

void GUI::checkEvents() {
    std::unique_lock<std::mutex> lock(mutex_);
    for (auto &e : event_queue_) handleEvent(e);
    event_queue_.clear();
}

void GUI::waitEvents() {
    // Just like the above except that we wait until there is at least one item in the queue:
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this]() -> bool { return not event_queue_.empty(); });

    for (auto &e : event_queue_) handleEvent(e);
    event_queue_.clear();
}

void GUI::handleEvent(const Event &event) {
    try {
        switch (event.type) {
            case Event::Type::setup:
                if (on_setup_) on_setup_(*event.parameters);
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
        }
    } catch (std::exception e) {
        // FIXME
        std::cerr << "Caught exception (" << e.what() << ") from handler.  FIXME!\n";
    }
}

}
