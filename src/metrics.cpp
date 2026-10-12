#include "sdf/metrics.hpp"
#include <chrono>
#include <numeric>
#include <limits>

namespace sdf {
Metrics evaluate(const Distance& f,const Samples& samples,double epsilon) {
    if(samples.empty() || epsilon<0 || !std::isfinite(epsilon)) throw std::invalid_argument("Invalid evaluation samples/epsilon");
    Metrics m; m.count=samples.size(); std::vector<double> errors; errors.reserve(samples.size());
    std::size_t wrong=0,inside_correct=0,outside_correct=0;
    for(const auto& s:samples) {
        const double pred=f(s.point),e=std::abs(pred-s.distance);
        if(!std::isfinite(e)) throw std::runtime_error("Non-finite evaluation distance");
        errors.push_back(e); m.rmse+=e*e; m.mae+=e;
        if(std::abs(s.distance)>epsilon) {
            ++m.sign_count;
            if(s.distance<0) { ++m.inside_count; if(pred<0) ++inside_correct; else ++wrong; }
            else { ++m.outside_count; if(pred>0) ++outside_correct; else ++wrong; }
        }
    }
    std::sort(errors.begin(),errors.end());
    const auto tail=std::max<std::size_t>(1,static_cast<std::size_t>(std::ceil(0.01*errors.size())));
    m.rmse=std::sqrt(m.rmse/m.count); m.mae/=m.count;
    m.p99=errors[static_cast<std::size_t>(std::ceil(0.99*errors.size()))-1];
    m.cvar99=std::accumulate(errors.end()-tail,errors.end(),0.0)/tail;
    m.max_error=errors.back();
    m.sign_error=m.sign_count?double(wrong)/m.sign_count:std::numeric_limits<double>::quiet_NaN();
    m.inside_recall=m.inside_count?double(inside_correct)/m.inside_count:std::numeric_limits<double>::quiet_NaN();
    m.outside_recall=m.outside_count?double(outside_correct)/m.outside_count:std::numeric_limits<double>::quiet_NaN();
    return m;
}
Timing time_queries(const Distance& f,const Samples& samples,unsigned repeats) {
    if(samples.empty() || repeats<3) throw std::invalid_argument("Timing needs samples and at least three repeats");
    volatile double sink=0;
    double checksum=0; for(const auto& s:samples) checksum+=f(s.point); sink=checksum;
    std::vector<double> times;
    for(unsigned r=0;r<repeats;++r) {
        const auto start=std::chrono::steady_clock::now(); double sum=0;
        for(const auto& s:samples) sum+=f(s.point);
        const auto end=std::chrono::steady_clock::now(); sink=sum;
        times.push_back(std::chrono::duration<double,std::nano>(end-start).count()/samples.size()); checksum+=sum;
    }
    (void)sink; std::sort(times.begin(),times.end());
    return {times[times.size()/2],times[static_cast<std::size_t>(std::ceil(0.95*times.size()))-1],checksum};
}
} // namespace sdf
