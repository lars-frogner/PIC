struct Interpolator
    interpolation::Interpolations.BSplineInterpolation

    function Interpolator(field::Array{<:Any,N}) where {N,D}
        new(Interpolations.interpolate(field, BSpline(Linear())))
    end
end

function interpolate(interpolator::Interpolator, continuous_indices)
    interpolator.interpolation(continuous_indices...)
end
