#include "sdf/compact_surface.hpp"
namespace sdf {
CompactExtraction extract_compact_surface(const CompactMeshField &field,unsigned cells) {
  if(cells<2||cells>512)throw std::invalid_argument("Invalid compact extraction resolution");
  const unsigned n=cells+1;const auto count=std::size_t(n)*n*n;
  CompactExtraction result;result.lattice.resize(count);std::vector<unsigned char> active(count);
  auto point=[&](std::size_t i){return Box{}.point({double(i%n)/cells,double((i/n)%n)/cells,double(i/(std::size_t(n)*n))/cells});};
  for(std::size_t i=0;i<count;++i){auto side=field.classify(point(i));++result.classification_sites;
    if(!side.resolved||!std::isfinite(side.distance)){++result.unresolved_sites;continue;}
    result.lattice[i]=side.boundary?0.:side.distance;
  }
  if(result.unresolved_sites)return result;
  for(unsigned z=0;z<cells;++z)for(unsigned y=0;y<cells;++y)for(unsigned x=0;x<cells;++x) {
    std::size_t ids[8];bool inside=false,outside=false;
    for(unsigned c=0;c<8;++c){auto id=(std::size_t(z+bool(c&4))*n+y+bool(c&2))*n+x+bool(c&1);ids[c]=id;inside|=result.lattice[id]<0;outside|=result.lattice[id]>=0;}
    if(inside&&outside)for(auto id:ids)active[id]=1;
  }
  for(std::size_t i=0;i<count;++i)if(active[i]) {
    auto answer=field.evaluate(point(i));++result.active_distance_sites;
    if(!answer.resolved||!std::isfinite(answer.distance)){++result.unresolved_sites;continue;}
    if((answer.distance<0)!=(result.lattice[i]<0))throw std::runtime_error("Classification/full-distance sign disagreement");
    result.lattice[i]=answer.distance;
  }
  if(result.unresolved_sites)return result;
  // extract_surface calls its field exactly once per lattice site in index
  // order, then keeps the standard zero perturbation and mesh traversal.
  std::size_t index=0;result.surface=extract_surface([&](Vec3){return result.lattice.at(index++);},cells);
  if(index!=count)throw std::runtime_error("Extraction lattice consumption changed");
  return result;
}
}
