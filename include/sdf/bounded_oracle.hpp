#pragma once
#include "sdf/generalization.hpp"
#include <cstring>
#include <list>
#include <unordered_map>
namespace sdf {
class BoundedOracle final:public Oracle {
    struct Key{std::array<std::uint64_t,3> bits;bool operator==(const Key& k)const{return bits==k.bits;}};
    struct Hash{std::size_t operator()(const Key& k)const{std::size_t h=0;for(auto v:k.bits)h^=std::hash<std::uint64_t>{}(v)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);return h;}};
    GeneralizationOracle base_;std::size_t capacity_;mutable std::list<std::pair<Key,double>> lru_;mutable std::unordered_map<Key,decltype(lru_)::iterator,Hash> index_;
public:
    mutable std::size_t logical{},physical{},hits{},evictions{};
    explicit BoundedOracle(GeometrySpec spec,std::size_t capacity=262144):base_(std::move(spec)),capacity_(capacity){if(!capacity_)throw std::invalid_argument("Zero LRU capacity");index_.reserve(capacity_);}
    double distance(Vec3 p)const override{++logical;Key k;double v[3]={p.x,p.y,p.z};std::memcpy(k.bits.data(),v,24);auto it=index_.find(k);if(it!=index_.end()){++hits;lru_.splice(lru_.begin(),lru_,it->second);return it->second->second;}++physical;auto d=base_.distance(p);if(index_.size()==capacity_){index_.erase(lru_.back().first);lru_.pop_back();++evictions;}lru_.push_front({k,d});index_.emplace(k,lru_.begin());return d;}
    Vec3 surface(std::mt19937_64& rng)const override{return base_.surface(rng);}
    std::string name()const override{return base_.name();}
    void clear(){index_.clear();lru_.clear();logical=physical=hits=evictions=0;}
    std::size_t entries()const{return index_.size();}
};
}
