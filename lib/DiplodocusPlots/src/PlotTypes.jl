abstract type PlotType end

struct Static <: PlotType end
struct Animated <: PlotType end
struct Interactive <: PlotType end