#pragma once
#include "creativity/gui/MemberStore.hpp"
#include "creativity/Reader.hpp"

namespace creativity { namespace gui {

/** Gtk::TreeModel::ColumnRecord subclass for handling Reader information in the list of readers in
 * the main GUI window.
 */
class ReaderStore : public MemberStore<Reader>, Glib::Object {
    public:
        /** Interface class between a simulation's Readers and a Gtk::TreeView.  This internally
         * stores a vector of Readers which can be updated (if needed) by calling the update method.
         *
         * This exposes 9 columns:
         * - ID
         * - x position
         * - y position
         * - position string (made from x position and y position)
         * - current utility
         * - lifetime utility
         * - books owned
         * - new books
         * - books written
         * - age of most recently written book (simulation age if no books written)
         */
        static Glib::RefPtr<ReaderStore> create(std::shared_ptr<eris::Simulation> sim);

        /** Sychronizes the list of readers with the stored Simulation. Typically called after a
         * simulation period runs. */
        virtual std::vector<eris::SharedMember<Reader>> resync_add() override;

        /** ColumnRecord object for a ReaderStore.  This object contains the columns for this Book
         * model.  This should not be used directly, but rather accessed via the public `columns`
         * member.
         */
        class ColRec : public Gtk::TreeModel::ColumnRecord {
            public:
                Gtk::TreeModelColumn<eris::eris_id_t> id;
                Gtk::TreeModelColumn<double> posX, posY, u, uLifetime;
                Gtk::TreeModelColumn<std::string> posstr;
                Gtk::TreeModelColumn<size_t> booksOwned, booksNew, booksWritten, lastBookAge;

            private:
                ColRec() {
                    add(id); add(posX); add(posY); add(posstr); add(u); add(uLifetime);
                    add(booksOwned); add(booksNew); add(booksWritten); add(lastBookAge);
                }
                friend class ReaderStore;
        };

        /** The columns of this ReaderStore.  For example, to access the utility column, use
         * `readerstore.columns.u`.
         */
        ColRec columns;

        /** Takes a Gtk::TreeView and adds this object's columns to it. */
        void appendColumnsTo(Gtk::TreeView &v) const override;

    protected:
        /// Protected constructor; this object should be constructed using create().
        ReaderStore(std::shared_ptr<eris::Simulation> &&sim);

        /** Returns the column type of the given position.  This is typically invoked via
         * get_column_type, itself given a column member of the `.columns` ColRec object.
         *
         * \sa ReaderStore::ColRec
         */
        virtual GType get_column_type_vfunc(int index) const override;

        /** Returns `obj.columns.size()`, the number of model columns. */
        virtual int get_n_columns_vfunc() const override;

        /** Accesses a column value.
         *
         * \param iter a valid iterator referencing the row to access
         * \param column the index of the column to access
         * \param value a Glib::Value<TYPE> object (where TYPE is the appropriate type for the
         * requested `column`) in which to store the value.
         */
        virtual void get_value_vfunc(const iterator &iter, int column, Glib::ValueBase &value) const override;

        // TreeSortable overrides:

        /** Sets the model sort column and sort order.  If the sort_column and order differ from the
         * current values, the model data is resorted.  If resorting occurs, the sort_column_changed
         * and rows_reordered signals will fire.
         *
         * This uses a stable sort: elements that are equal under the new sort column will preserve
         * their current ordering.  Note that this means that sorting by the current column but in
         * the opposite order will *not* reverse the ordering of equal-value elements.
         *
         * \param sort_column_id the index of the new sort column
         * \param order the new sort order (Gtk::SORT_ASCENDING or Gtk::SORT_DESCENDING).
         */
        virtual void set_sort_column_id_vfunc(int sort_column_id, Gtk::SortType order) override;


    private:
        // The various comparison functions; one of these gets passed to std::stable_sort.
        static bool less_id(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_id(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_posX(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_posY(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_posstr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_posstr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_uCurr(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_uLife(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_booksOwned(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_booksNew(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_booksWritten(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool less_lastBookAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);
        static bool greater_lastBookAge(const eris::SharedMember<Reader> &a, const eris::SharedMember<Reader> &b);

};

}}
