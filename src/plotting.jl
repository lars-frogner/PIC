meshgrid(x, y) = (repeat(x, outer = length(y)), repeat(y, inner = length(x)))

@userplot ScalarHeatMap
@recipe function f(p::ScalarHeatMap)
    field, = p.args
    @series begin
        seriestype := :heatmap
        aspect_ratio := :equal
        transpose(field)
    end
end

@userplot VectorFieldQuiver
@recipe function f(p::VectorFieldQuiver)
    x, y, vectorfield = p.args
    n_arrows = 20
    skip_x = ceil(Int, length(x) / n_arrows)
    skip_y = ceil(Int, length(y) / n_arrows)
    vx_mesh = vectorfield[1][1:skip_x:end, 1:skip_y:end][:]
    vy_mesh = vectorfield[2][1:skip_x:end, 1:skip_y:end][:]
    x_mesh, y_mesh = meshgrid(x[1:skip_x:end], y[1:skip_y:end])
    Lx = x[end] - x[1]
    max_v = maximum(sqrt.(vx_mesh .^ 2 .+ vy_mesh .^ 2))
    scale = 0.2 * Lx / max_v
    headlength = 0.3
    headwidth = 0.2
    @series begin
        seriestype := :quiver
        aspect_ratio := :equal
        xlims := (x[1], x[end])
        ylims := (y[1], y[end])
        arrow := arrow(:closed, :head, headlength, headwidth)
        quiver := (vx_mesh * scale, vy_mesh * scale)
        x_mesh, y_mesh
    end
end
