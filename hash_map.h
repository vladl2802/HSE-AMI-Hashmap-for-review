//
// Created by Vlad on 22.01.2023.
//

#ifndef HASHTABLE_HASH_MAP_H
#define HASHTABLE_HASH_MAP_H

#include <iostream>
#include <list>
#include <random>
#include <unordered_map>
#include <algorithm>
#include <cassert>

const int SIZES[] = {3, 7, 23, 47, 107, 227, 467, 983, 2027, 4079, 8423, 17027, 34319, 68639, 137279, 274679, 549503,
                     1099103, 2198447, 4396907, 8793839, 17587763, 35175527, 70351223};
const unsigned int SIZES_LEN = 23;
const unsigned int BUCKETS_COUNT = 8;

template<class KeyType, class ValueType, class Hash = std::hash<KeyType>, class KeyEqual = std::equal_to<KeyType>>
class HashMap {
    //using KeyType = int;
    //using ValueType = int;
    using SizeType = unsigned int;
    //using Hash = std::hash<KeyType>;
    //using KeyEqual = std::equal_to<KeyType>;
    using Element = std::pair<const KeyType, ValueType>;
    const bool CHECK_VALID_AFTER = false;

public:
    class iterator {
        friend class HashMap;

        SizeType cur_;
        std::list<SizeType>::iterator pos_;
        const HashMap* parent_;

    public:
        explicit iterator(const HashMap* parent = nullptr);
        iterator(SizeType cur, std::list<SizeType>::iterator pos, const HashMap* parent);

        std::pair<const KeyType, ValueType>& operator*() const;
        std::pair<const KeyType, ValueType>* operator->();

        iterator& operator++();
        iterator operator++(int);

        bool operator==(const iterator& other) const;
        bool operator!=(const iterator& other) const;
    };

    class const_iterator {
        friend class HashMap;

        SizeType cur_;
        std::list<SizeType>::iterator pos_;
        const HashMap* parent_;

    public:
        explicit const_iterator(const HashMap* parent = nullptr);
        const_iterator(SizeType cur, std::list<SizeType>::iterator pos, const HashMap* parent);
        explicit const_iterator(iterator it);

        const std::pair<const KeyType, ValueType>& operator*() const;
        const std::pair<const KeyType, ValueType>* operator->();

        const_iterator& operator++();
        const_iterator operator++(int);

        bool operator==(const const_iterator& other) const;
        bool operator!=(const const_iterator& other) const;
    };
private:
    
    struct TableBucket {
        bool free;
        std::list<SizeType>::iterator pos; // position in jumps list
        SizeType next; // null if equal to _buckets
        SizeType prev; // null if equal to _buckets
        // SizeType end;  // null if equal to _buckets
        Element* elem;

        explicit TableBucket(bool free = true,
                              std::list<SizeType>::iterator pos = static_cast<std::list<SizeType>::iterator>(nullptr),
                              SizeType next = 0, SizeType prev = 0, Element* elem = nullptr);
        ~TableBucket();
    };

    friend class iterator;
    friend class const_iterator;

    SizeType buckets_;
    Hash hash_;
    KeyEqual key_equal_;
    TableBucket* table_;
    std::list<SizeType> jumps_;
    SizeType size_;
    SizeType last_free_;
    iterator begin_;
    iterator end_;

    SizeType real_hash(const KeyType& key) const;
    std::pair<iterator, bool> unchecked_insert(std::pair<KeyType, ValueType>&& elem);
    bool check_overload();
    bool check_valid();
public:
    explicit HashMap(const Hash &hash = Hash(), const KeyEqual &key_equal = KeyEqual());

    template<class order_iterator>
    HashMap(order_iterator begin, order_iterator end, const Hash &hash = Hash(),
            const KeyEqual &key_equal = KeyEqual());
    HashMap(std::initializer_list<std::pair<KeyType, ValueType>> il, const Hash &hash = Hash(),
            const KeyEqual &key_equal = KeyEqual());
    HashMap(const HashMap& other);
    HashMap& operator=(const HashMap& other);
    //HashMap& operator=(HashMap&& other);
    //HashMap(HashMap&& other) noexcept = default;
    ~HashMap();

    [[nodiscard]] bool empty() const;
    [[nodiscard]] SizeType size() const;

    [[nodiscard]] Hash hash_function() const;
    [[nodiscard]] KeyEqual key_equal_function() const;

    std::pair<iterator, bool> insert(std::pair<KeyType, ValueType>&& elem);
    bool erase(const KeyType& key);

    iterator begin();
    [[nodiscard]] const_iterator begin() const;
    iterator end();
    [[nodiscard]] const_iterator end() const;

    [[nodiscard]] const_iterator find(const KeyType& key) const;
    iterator find(const KeyType& key);

    ValueType& operator[](const KeyType& key);
    ValueType& operator[](KeyType&& key);
    [[nodiscard]] const ValueType& at(const KeyType& key) const;

    void reserve(SizeType new_size);
    void clear();
};

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::TableBucket::TableBucket(bool free, std::list<SizeType>::iterator pos,
                                                                        SizeType next, SizeType prev, Element* elem)
        : free(free), pos(pos), next(next), prev(prev), elem(elem) {
    //std::cerr << "TableBucket constructor\n";
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::TableBucket::~TableBucket() {
    delete elem;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::const_iterator(const HashMap *parent) : cur_(0),
                                                                                                     pos_(nullptr),
                                                                                                     parent_(parent) {}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::const_iterator(HashMap::SizeType cur,
                                                                            std::list<HashMap::SizeType>::iterator pos,
                                                                            const HashMap *parent)  : cur_(cur),
                                                                                                      pos_(pos),
                                                                                                      parent_(parent) {}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::const_iterator(HashMap::iterator it) : cur_(it.cur_),
                                                                                                    pos_(it.pos_),
                                                                                                    parent_(it.parent_) {}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
const std::pair<const KeyType, ValueType> &
HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::operator*() const {
    return *parent_->table_[cur_].elem;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
const std::pair<const KeyType, ValueType> *HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::operator->() {
    return parent_->table_[cur_].elem;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator &
HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::operator++() {
    if (parent_->table_[cur_].next == parent_->buckets_) {
        ++pos_;
        cur_ = *pos_;
    } else {
        cur_ = parent_->table_[cur_].next;
    }
    return *this;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::operator++(int) {
    auto temp = *this;
    ++*this;
    return temp;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::operator==(
        const HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator &other) const {
    return cur_ == other.cur_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator::operator!=(
        const HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator &other) const {
    return cur_ != other.cur_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::iterator(const HashMap *parent) : cur_(0), pos_(nullptr),
                                                                                         parent_(parent) {}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::iterator(SizeType cur, std::list<SizeType>::iterator pos,
                                                                const HashMap *parent) : cur_(cur), pos_(pos),
                                                                                         parent_(parent) {}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
std::pair<const KeyType, ValueType>& HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::operator*() const {
    return *parent_->table_[cur_].elem;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
std::pair<const KeyType, ValueType>* HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::operator->() {
    return parent_->table_[cur_].elem;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator &
HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::operator++() {
    if (parent_->table_[cur_].next == parent_->buckets_) {
        ++pos_;
        cur_ = *pos_;
    } else {
        cur_ = parent_->table_[cur_].next;
    }
    return *this;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::operator++(int) {
    auto temp = *this;
    ++*this;
    return temp;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::operator==(
        const HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator &other) const {
    return cur_ == other.cur_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator::operator!=(
        const HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator &other) const {
    return cur_ != other.cur_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::HashMap(const Hash &hash, const KeyEqual &key_equal) :
        buckets_(SIZES[0]), hash_(hash), key_equal_(key_equal), table_(new TableBucket[buckets_ + 1]), size_(0),
        last_free_(buckets_ - 1) {
    auto it = jumps_.insert(jumps_.end(), buckets_); // push the ned
    begin_ = {buckets_, it, this};
    end_ = {buckets_, it, this};
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
template<class order_iterator>
HashMap<KeyType, ValueType, Hash, KeyEqual>::HashMap(order_iterator begin, order_iterator end, const Hash &hash,
                                                     const KeyEqual &key_equal) :
        buckets_(SIZES[0]), hash_(hash), key_equal_(key_equal), table_(new TableBucket[buckets_ + 1]), size_(0),
        last_free_(buckets_ - 1) {
    auto it = jumps_.insert(jumps_.end(), buckets_); // push the ned
    begin_ = {buckets_, it, this};
    end_ = {buckets_, it, this};

    for (auto cur = begin; cur != end; ++cur) {
        insert(*cur);
    }
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::HashMap(std::initializer_list<std::pair<KeyType, ValueType>> il,
                                                     const Hash &hash, const KeyEqual &key_equal) :
        buckets_(SIZES[0]), hash_(hash), key_equal_(key_equal), table_(new TableBucket[buckets_ + 1]), size_(0),
        last_free_(buckets_ - 1) {
    auto it = jumps_.insert(jumps_.end(), buckets_); // push the end
    begin_ = {buckets_, it, this};
    end_ = {buckets_, it, this};

    for (auto tmp : il) {
        insert(std::move(tmp));
    }
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::HashMap(const HashMap<KeyType, ValueType, Hash, KeyEqual> &other) :
        buckets_(other.buckets_), hash_(other.hash_), key_equal_(other.key_equal_), jumps_(other.jumps_),
        size_(other.size_), last_free_(buckets_ - 1) {
    begin_ = {*jumps_.begin(), jumps_.begin(), this};
    end_ = {buckets_, --jumps_.end(), this};

    table_ = new TableBucket[buckets_ + 1];
    for (SizeType i = 0; i < buckets_ + 1; ++i) {
        table_[i] = other.table_[i];
        if (other.table_[i].elem == nullptr) table_[i].elem = nullptr;
        else table_[i].elem = new Element(*other.table_[i].elem);
    }
}

/*template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::HashMap(HashMap<KeyType, ValueType, Hash, KeyEqual> &&other) noexcept {
    buckets_ = std::move(other.buckets_);
    hash_ = std::move(other.hash);
    key_equal_ = std::move(other.key_equal_);
    table_ = std::move(other.table_);
    jumps_ = std::move(other.jumps_);
    size_ = std::move(other.size_);
    last_free_ = std::move(other.last_free_);
    begin_ = std::move(other.begin_);
    end_ = std::move(other.end_);
}*/

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual>::~HashMap() {
    delete[] table_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
HashMap<KeyType, ValueType, Hash, KeyEqual> &HashMap<KeyType, ValueType, Hash, KeyEqual>::operator=(const HashMap &other) {
    if (this == &other) return *this;
    delete[] table_;

    buckets_ = other.buckets_;
    hash_ = other.hash_;
    key_equal_ = other.key_equal_;
    jumps_ = other.jumps_;
    size_ = other.size_;
    last_free_ = buckets_ - 1;

    begin_ = {*jumps_.begin(), jumps_.begin(), this};
    end_ = {buckets_, jumps_.end(), this};

    table_ = new TableBucket[buckets_ + 1];
    for (SizeType i = 0; i < buckets_ + 1; ++i) {
        table_[i] = other.table_[i];
        if (other.table_[i].elem == nullptr) table_[i].elem = nullptr;
        else table_[i].elem = new Element(*other.table_[i].elem);
    }
    return *this;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::SizeType
HashMap<KeyType, ValueType, Hash, KeyEqual>::real_hash(const KeyType &key) const {
    SizeType tmp = hash_(key) % buckets_;
    //if (tmp < 0) tmp += buckets_;
    return tmp;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::check_overload() {
    return size_ * 10 >= buckets_ * 9;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
std::pair<typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator, bool>
HashMap<KeyType, ValueType, Hash, KeyEqual>::unchecked_insert(std::pair<KeyType, ValueType>&& elem) {
    auto &key = elem.first;
    auto hash = real_hash(key);
    SizeType pos, prev;
    std::list<SizeType>::iterator jumps_it;
    if (table_[hash].free) {
        // lucky day, just start chain with hash from here
        jumps_it = jumps_.insert(--jumps_.end(), hash);
        prev = buckets_;
        pos = hash;
    } else {
        // find any empty slot in table to fit in element
        while (!table_[last_free_].free) {
            if (last_free_ == 0) last_free_ = buckets_;
            last_free_--;
        }
        // check chaining in hash
        if (table_[hash].prev != buckets_) {
            // chain go through hash, that's not nice, need to move this element away
            SizeType t_prev = table_[hash].prev, t_next = table_[hash].next;
            auto& tmp = table_[last_free_];
            tmp.free = false;
            tmp.pos = table_[hash].pos;
            tmp.next = t_next;
            tmp.prev = t_prev;
            std::swap(tmp.elem, table_[hash].elem);
            table_[t_prev].next = last_free_;
            if (t_next != buckets_) table_[t_next].prev = last_free_;
            // and start chain with hash from here
            jumps_it = jumps_.insert(--jumps_.end(), hash);
            prev = buckets_;
            pos = hash;
        } else {
            // chaining is ok, try to find key in chain
            SizeType cur = hash;
            while (true) {
                if (table_[cur].elem->first == key) {
                    return {{0, static_cast<std::list<SizeType>::iterator>(nullptr), this}, false};
                }
                if (table_[cur].next == buckets_) break;
                cur = table_[cur].next;
            }
            // if nothing found then find the use end of chain and connect element
            table_[cur].next = last_free_;
            jumps_it = static_cast<std::list<SizeType>::iterator>(nullptr);
            pos = last_free_;
            prev = cur;
        }
    }
    auto& tmp = table_[pos];
    tmp.free = false;
    tmp.pos = jumps_it;
    tmp.next = buckets_;
    tmp.prev = prev;
    delete tmp.elem;
    tmp.elem = new Element(std::forward<std::pair<KeyType, ValueType>>(elem));
    auto result = iterator(pos, jumps_it, this);
    if (size_ == 0) begin_ = result; // first inserted element become begin
    size_++;
    return {result, true};
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::SizeType
HashMap<KeyType, ValueType, Hash, KeyEqual>::size() const {
    return size_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::empty() const {
    return size_ == 0;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
Hash HashMap<KeyType, ValueType, Hash, KeyEqual>::hash_function() const {
    return hash_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
KeyEqual HashMap<KeyType, ValueType, Hash, KeyEqual>::key_equal_function() const {
    return key_equal_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
std::pair<typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator, bool>
HashMap<KeyType, ValueType, Hash, KeyEqual>::insert(std::pair<KeyType, ValueType>&& elem) {
    if (CHECK_VALID_AFTER) {
        assert(check_valid());
    }
    if (check_overload()) reserve(buckets_ * 2);
    return unchecked_insert(std::forward<std::pair<KeyType, ValueType>>(elem));
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::erase(const KeyType &key) {
    if (CHECK_VALID_AFTER) {
        assert(check_valid());
    }

    auto hash = real_hash(key);
    // first mean that no chain go in hash. second mean that if chain go through hash, then no chain start in hash
    if (table_[hash].free || table_[hash].prev != buckets_) return false;
    // if something start in hash, just go through chain and check it all
    SizeType cur = hash;
    while (cur != buckets_) {
        if (table_[cur].elem->first == key) {
            // found elem to delete
            SizeType pos = buckets_;
            if (table_[cur].prev == buckets_) {
                // that's begin of chain, it can't be moved, thus if next element exists move it to begin of chain else delete chain
                if (table_[cur].next == buckets_) {
                    // chain should be deleted
                    auto it = jumps_.erase(table_[cur].pos);
                    pos = cur;

                    // if this chain was begin, then we need to change begin
                    if (begin_.cur_ == cur) {
                        begin_ = {*it, it, this};
                    }
                } else {
                    // next element should be moved to cur
                    auto next = table_[cur].next;
                    pos = next;
                    std::swap(table_[cur].elem, table_[next].elem);
                    table_[cur].next = table_[next].next;
                    if (next != buckets_) {
                        auto next_next = table_[next].next;
                        table_[next_next].prev = cur;
                    }
                }
            } else {
                // that's middle of chain, just connect left and right parts together
                auto prev = table_[cur].prev, next = table_[cur].next;
                table_[prev].next = next;
                if (next != buckets_) table_[next].prev = prev;
                pos = cur;
            }
            // element was found
            table_[pos].free = true;
            delete table_[pos].elem;
            table_[pos].elem = nullptr;
            size_--;
            // and deleted
            return true;
        }
        cur = table_[cur].next;
    }
    // element wasn't found
    return false;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::begin() {
    return begin_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::begin() const {
    return const_iterator(begin_);
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::end() {
    return end_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::end() const {
    return const_iterator(end_);
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::const_iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::find(const KeyType &key) const {
    auto hash = real_hash(key);
    // first mean that no chain go in hash. second mean that if chain go through hash, then no chain start in hash
    if (table_[hash].free || table_[hash].prev != buckets_) return const_iterator(end_);
    // if something start in hash, just go through chain and check it all
    auto cur = hash;
    while (cur != buckets_) {
        if (table_[cur].elem->first == key) {
            return {cur, table_[hash].pos, this};
        }
        cur = table_[cur].next;
    }
    return const_iterator(end_);
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
typename HashMap<KeyType, ValueType, Hash, KeyEqual>::iterator
HashMap<KeyType, ValueType, Hash, KeyEqual>::find(const KeyType &key) {
    auto hash = real_hash(key);
    // first mean that no chain go in hash. second mean that if chain go through hash, then no chain start in hash
    if (table_[hash].free || table_[hash].prev != buckets_) return end_;
    // if something start in hash, just go through chain and check it all
    auto cur = hash;
    while (cur != buckets_) {
        if (table_[cur].elem->first == key) {
            return {cur, table_[hash].pos, this};
        }
        cur = table_[cur].next;
    }
    return end_;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
ValueType &HashMap<KeyType, ValueType, Hash, KeyEqual>::operator[](const KeyType& key) {
    auto it = find(key);
    if (it == end_) {
        auto res = insert({key, ValueType()});
        it = res.first;
    }
    return it->second;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
ValueType &HashMap<KeyType, ValueType, Hash, KeyEqual>::operator[](KeyType &&key) {
    auto it = find(key);
    if (it == end_) {
        auto res = insert({std::forward<KeyType>(key), ValueType()});
        it = res.first;
    }
    return it->second;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
const ValueType &HashMap<KeyType, ValueType, Hash, KeyEqual>::at(const KeyType &key) const {
    auto it = find(key);
    if (it == end()) throw std::out_of_range("Stupid kiddo");
    return it->second;
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
void HashMap<KeyType, ValueType, Hash, KeyEqual>::reserve(SizeType new_size) {
    if (new_size <= buckets_) return;

    auto best = std::lower_bound(SIZES, SIZES + SIZES_LEN, new_size) - SIZES;
    new_size = SIZES[best];

    auto temp_table = new TableBucket[new_size + 1];
    std::swap(temp_table, table_);
    std::swap(new_size, buckets_);
    jumps_.clear();
    auto it = jumps_.insert(jumps_.end(), buckets_); // return the end
    begin_ = {buckets_, it, this};
    end_ = {buckets_, it, this};
    last_free_ = buckets_ - 1;
    size_ = 0;
    for (SizeType i = 0; i < new_size; ++i) {
        if (temp_table[i].free || temp_table[i].prev != new_size) continue;
        auto cur = i;
        while (cur != new_size) {
            unchecked_insert(std::move(*temp_table[cur].elem));
            cur = temp_table[cur].next;
        }
    }
    delete[] temp_table;

    if (CHECK_VALID_AFTER) {
        assert(check_valid());
    }
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
void HashMap<KeyType, ValueType, Hash, KeyEqual>::clear() {
    // no shrinking of buckets
    for (const auto& it : jumps_) {
        auto cur = it;
        while (cur != buckets_) {
            table_[cur].free = true;
            delete table_[cur].elem;
            table_[cur].elem = nullptr;
            cur = table_[cur].next;
        }
    }
    jumps_.clear();
    auto it = jumps_.insert(jumps_.end(), buckets_); // return the end
    begin_ = {buckets_, it, this};
    end_ = {buckets_, it, this};
    last_free_ = buckets_ - 1;
    size_ = 0;

    if (CHECK_VALID_AFTER) {
        assert(check_valid());
    }
}

template<class KeyType, class ValueType, class Hash, class KeyEqual>
bool HashMap<KeyType, ValueType, Hash, KeyEqual>::check_valid() {
    SizeType cnt = 0, tmp = 0;
    for (SizeType i = 0; i < buckets_; ++i) {
        if (table_[i].free) {
            // if element marked free then other param can be invalid except elem, this should be nullptr (and deleted before)
            if (table_[i].elem != nullptr) {
                std::cerr << "elem is free but not nullptr\n";
                return false;
            }
        } else {
            tmp++;
            // else param must be valid
            // firstly check, that elem not nullptr
            if (table_[i].elem == nullptr) {
                std::cerr << "elem isn't free but nullptr\n";
                return false;
            }
            // elem is begin of chain if prev = buckets_ therefore pos should be not nullptr
            if (table_[i].prev == buckets_) {
                auto key = table_[i].elem->first;
                if (real_hash(key) != i) {
                    std::cerr << "elem is begin of chain, but hash is wrong\n";
                    return false;
                }
                cnt++;
            }
            if (table_[i].prev == buckets_ && table_[i].pos == static_cast<std::list<SizeType>::iterator>(nullptr)) {
                std::cerr << "it's begin of chain but no iterator to list\n";
                return false;
            }
            // check valid of next and prev indexes
            if (table_[i].next != buckets_) {
                if (table_[table_[i].next].prev != i) {
                    std::cerr << "next or prev index is invalid\n";
                    return false;
                }
            }
            if (table_[i].prev != buckets_) {
                if (table_[table_[i].prev].next != i) {
                    std::cerr << "next or prev index is invalid\n";
                    return false;
                }
            }
        }
    }
    // cnt - is count of chains, check is jumps_ contains all chains
    if (cnt + 1 != jumps_.size()) {
        std::cerr << "jumps size not match chains count\n";
        return false;
    }
    // tmp - is count of elements, check is size valid
    if (tmp != size_) {
        std::cerr << "size is invalid\n";
        return false;
    }
    return true;
}

#endif //HASHTABLE_HASH_MAP_H