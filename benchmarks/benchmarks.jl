# # Benchmarks for GeometryOps

#=
=#

# These are the geospatial packages we'll use.  LibGEOS is the gold standard as far as 
# geographic packages go, and GeometryBasics is useful to construct various geometries.
# GeoJSON is the format of choice for single-file geometry, so it's useful to load in various geometries
# as well.
import GeometryOps as GO, 
    GeoInterface as GI, 
    GeometryBasics as GB, 
    LibGEOS as LG
import GeoJSON
# In order to benchmark, we'll actually use the new [Chairmarks.jl](https://github.com/lilithhafner/Chairmarks.jl), 
# since it's significantly faster than BenchmarkTools.  To keep benchmarks organized, though, we'll still use BenchmarkTools' 
# `BenchmarkGroup` structure.
using Chairmarks
import BenchmarkTools: BenchmarkGroup
# We use CairoMakie to visualize our results!
using CairoMakie, MakieThemes, GeoInterfaceMakie
# Finally, we import some general utility packages:
using Statistics

# We also set up some utility functions for later on.
"""
Returns LibGEOS and GeometryOps' native geometries, 
so as not to disadvantage either package.
"""
lg_and_go(geometry) = (GI.convert(LG, geometry), GO.tuples(geometry))

"This is the main benchmark suite, from which all other suites flow."
SUITE = BenchmarkGroup()

#=
## Centroid and area

Centroids and areas have very similar calculation algorithms, since they all iterate over points.

Thus, geometries which challenge one algorithm will challenge the other.

We'll start by defining some geometry:
=#

circle_area_suite = SUITE["area"]["circle"] = BenchmarkGroup(["title = Area of a circle", "subtitle = Regular circle"])

n_points_values = round.(Int, exp10.(LinRange(log10(10), log10(100_000), 5)))

@time for n_points in n_points_values
    circle = GI.Wrappers.Polygon([[reverse(sincos(θ)) for θ in LinRange(0, 2π, n_points)]])
    closed_circle = GO.ClosedRing()(circle)
    lg_circle, go_circle = lg_and_go(closed_circle)
    circle_area_suite["GeometryOps"][n_points] = @be GO.area($go_circle) seconds=1
    circle_area_suite["LibGEOS"][n_points] = @be LG.area($lg_circle) seconds=1
end

plot_trials(circle_area_suite)


# ## Vancouver watershed benchmarks
#=

Vancouver Island has ~1,300 watersheds.  LibGEOS uses this exact data
in their tests, so we'll use it in ours as well!

We'll start by loading the data, and then we'll use it to benchmark various operations.

=#

# The CRS for this file is EPSG:3005, or as a PROJ string,
# `"+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"`
# TODO: this doesn't work with WellKnownGeometry.  Why?
wkt_geoms = LG.readgeom.(eachline("/Users/anshul/Downloads/watersheds.wkt"), (LG.WKTReader(LG.get_global_context()),))
vancouver_polygons = GI.getgeom.(wkt_geoms, 1); #.|> GO.tuples;

import SortTileRecursiveTree as STR
tree = STR.STRtree(vancouver_polygons)
query_result = STR.query(tree, GI.extent(vancouver_polygons[1]))

GO.intersects.((vancouver_polygons[1],), vancouver_polygons[query_result])

all_intersected = falses(length(vancouver_polygons))
accumulator = deepcopy(vancouver_polygons[1])
all_intersected[1] = true
i = 1
# query_result = STR.query(tree, GI.extent(accumulator))
# for idx in query_result
#     println("Assessing $idx")
#     if !all_intersected[idx] && LG.intersects(vancouver_polygons[idx], accumulator)
#         println("Assimilating $idx")
#         result = LG.union(vancouver_polygons[idx], accumulator#=; target = GI.PolygonTrait()=#)
#         # @show length(result)
#         accumulator = result#[1]
#         all_intersected[idx] = true
#     end
# end
display(poly(vancouver_polygons[all_intersected]; color = rand(RGBf, sum(all_intersected))))
display(poly(accumulator))
@time while !(all(all_intersected)) && i < length(vancouver_polygons)
    println("Began iteration $i")
    query_result = STR.query(tree, GI.extent(accumulator))
    for idx in query_result
        if !(all_intersected[idx] || !(LG.intersects(vancouver_polygons[idx], accumulator)))
            println("Assimilating $idx")
            result = LG.union(vancouver_polygons[idx], accumulator#=; target = GI.PolygonTrait()=#)
            # @show length(result)
            accumulator = result#[1]
            all_intersected[idx] = true
        end
    end
    display(poly(vancouver_polygons[all_intersected]; color = rand(RGBf, sum(all_intersected))))
    println("Finished iteration $i")
    # wait_for_key("Press any key to continue to the next iteration.")
    i += 1
end 

fig = Figure()
ax = Axis(fig[1, 1]; title = "STRTree query for polygons", aspect = DataAspect())
for (i, result_index) in enumerate(query_result)
    poly!(ax, vancouver_polygons[result_index]; color = Makie.wong_colors()[i], label = "$result_index")
end
Legend(fig[1, 2], ax)
fig



# TODO: 
# - Vancouver watersheds:
#    - Intersection on pre-buffered geometry
#    - Polygon union by reduction (perhaps pre-sort by border order, so we don't end up with useless polygons)
#    - Queries using STRTree.jl
#    - Potentially using a prepared geometry based approach to multithreaded reductive joining
#    - Implement multipolygon joining.  How?  Query intersection or touching for each individual geometry,
#      and implement a 


## Segmentization

