const Float = Float64

SVectorF{N} = SVector{N,Float}
SVectorU{N} = SVector{N,UInt}

cross(a, b) = [
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1],
]
