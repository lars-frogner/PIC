b₀(ξ) = abs(ξ) < 0.5 ? 1 : 0
b₁(ξ) = (ξ >= -1 && ξ <= 1) ? (1 - abs(ξ)) : 0
b₂(ξ) = (ξ >= -1.5 && ξ <= 1.5) ?
    ((ξ >= -0.5 && ξ <= 0.5) ? (0.75 - ξ^2) : (0.5 * (1.5 - abs(ξ))^2)) : 0
