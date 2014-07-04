#pragma once
#include <memory>
#include <gtkmm/treemodel.h>
#include <gtkmm/treeview.h>
#include <eris/Simulation.hpp>
#include "creativity/Reader.hpp"

namespace creativity { namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Reader information in the list of readers in
 * the main GUI window.
 */
class ReaderStore : public Gtk::TreeModel, Glib::Object {
    public:
        ReaderStore() = delete;

        /** Interface class between a simulation's Readers and a Gtk::TreeView.  This internally
         * stores a vector of Readers which can be updated (if needed) by calling the update method.
         *
         * This exposes 9 columns:
         * - ID
         * - x position
         * - y position
         * - current utility
         * - lifetime utility
         * - books owned
         * - new books
         * - books written
         * - age of most recently written book (simulation age if no books written)
         */
        static Glib::RefPtr<ReaderStore> create(std::shared_ptr<eris::Simulation> sim);

        /** Sychronizes the list of readers with the stored Simulation. */
        void resync();

        class Columns : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<eris::eris_id_t> id;
                Gtk::TreeModelColumn<double> posX, posY, u, uLifetime;
                Gtk::TreeModelColumn<std::string> posstr;
                Gtk::TreeModelColumn<size_t> booksOwned, booksNew, booksWritten, bookLastAge;

                Columns() {
                    add(id); add(posX); add(posY); add(posstr); add(u); add(uLifetime);
                    add(booksOwned); add(booksNew); add(booksWritten); add(bookLastAge);
                }
        };

        Columns columns;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        void appendColumnsTo(Gtk::TreeView &v) const;

    protected:
        /// Protected constructor; this object should be constructed using create().
        ReaderStore(std::shared_ptr<eris::Simulation> &&sim);

        /** Returns Gtk::TreeModel flags (specifically, the LIST_ONLY flag). */
        virtual Gtk::TreeModelFlags get_flags_vfunc() const override;
        /** Returns 9, the number of reader data elements mapped into model columns. */
        virtual int get_n_columns_vfunc() const override;
        /** Returns the column type of the given position.  See the list of virtual columns in the
         * class description.
         *
         * \sa ReaderStore
         */
        virtual GType get_column_type_vfunc(int index) const override;

        virtual bool get_iter_vfunc(const Path &path, iterator& iter) const override;
        virtual bool iter_next_vfunc(const iterator &iter, iterator &iter_next) const override;
        virtual bool iter_children_vfunc(const iterator &parent, iterator &iter) const override;
        virtual bool iter_parent_vfunc(const iterator &child, iterator &iter) const override;
        virtual bool iter_nth_child_vfunc(const iterator &parent, int n, iterator &iter) const override;
        virtual bool iter_nth_root_child_vfunc(int n, iterator &iter) const override;
        virtual bool iter_has_child_vfunc(const iterator &iter) const override;
        virtual int iter_n_children_vfunc(const iterator &iter) const override;
        virtual int iter_n_root_children_vfunc() const override;
//        virtual void ref_node_vfunc(const iterator &iter) const override;
//        virtual void unref_node_vfunc(const iterator &iter) const override;
        virtual Path get_path_vfunc(const iterator &iter) const override;
        virtual void get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const override;

    private:
        // The maximum reader eris_id_t currently in readers_, used to identify new readers during resync():
        eris::eris_id_t max_id_;
        std::shared_ptr<eris::Simulation> sim_;
        std::vector<eris::SharedMember<Reader>> readers_;
        // Tracks model changes by being incremented whenever such a change occurs
        int stamp_ = 1;

        /*
        // sort_by_id isn't needed: the default sort (std::less) works fine
        static bool sort_by_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->position()[0] < b->position()[0]; }
        static bool sort_by_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->position()[1] < b->position()[1]; }
        // First x, then y for ties
        static bool sort_by_pos(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
            auto ax = a->position()[0], bx = b->position()[0];
            return ax == bx ? a->position()[1] < b->position()[1] : ax < bx;
        }
        static bool sort_by_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->u() < b->u(); }
        static bool sort_by_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->uLifetime() < b->uLifetime(); }
        static bool sort_by_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->library().size() < b->library().size(); }
        static bool sort_by_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->newBooks().size() < b->newBooks().size(); }
        static bool sort_by_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) { return a->wrote().size() < b->wrote().size(); }
        static bool sort_by_bookLatestAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b) {
            return (a->wrote().empty() ? std::numeric_limits<unsigned long>::max() : a->wrote().back()->age())
                 < (b->wrote().empty() ? std::numeric_limits<unsigned long>::max() : b->wrote().back()->age());
        }
        */

        template <typename T, typename = typename std::enable_if<std::is_base_of<Gtk::TreeModelColumnBase, T>::value>>
        void appendCol(Gtk::TreeView &v, const std::string &label, T &col, int width, bool sortable = true) const {
            v.append_column(label, col);
            auto *c = v.get_column(v.get_n_columns()-1);
            c->set_sizing(Gtk::TREE_VIEW_COLUMN_FIXED);
            c->set_fixed_width(width);
            if (sortable)
                c->set_sort_column(col);
        }

};

}}
