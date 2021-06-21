mutable struct ComplexNums
    r::Vector{Float64}
    x::Vector{Float64}
    z::Vector{Complex}
    g::Vector{Float64}
    b::Vector{Float64}
    y::Vector{Complex}

    function ComplexNums(r::Vector{Float64},
                         x::Vector{Float64},
                         z::Vector{Complex{T}},
                         g::Vector{Float64},
                         b::Vector{Float64},
                         y::Vector{Complex{T}}) where T
        new(r,x,z,g,b,y)
    end
end

function ComplexNums(i::Vector{Float64},
                     j::Vector{Float64};
                     type = "impedance")
    if type == "impedance"
        return ImpedanceComplexNums(i, j)
    end
    
    return AdmitanceComplexNums(i, j)
end

function ImpedanceComplexNums(r::Vector{Float64},
                              x::Vector{Float64})
    z = Complex.(r, x)
    y = z.^-1
    g = real.(z)
    b = imag.(z)
    
    return ComplexNums(r,x,z,g,b,y)
end

function AdmitanceComplexNums(g::Vector{Float64},
                              b::Vector{Float64})
    y = Complex.(g, b)
    z = y.^-1
    r = real.(y)
    x = imag.(y)
    return ComplexNums(r,x,z,g,b,y)
end


mutable struct Buses
    B::Vector{Int64}
    buses::Vector{Int64}
    capacitor_bank::Vector

    function Buses(buses::Vector{Int64}, 
                   capacitor_bank::Vector) where T
        B = collect(1:length(buses))
        new(B, buses, capacitor_bank)
    end
end

mutable struct Lines
    L::Vector{Int64}
    lines::Vector{Int64}
    from::Vector{Int64}
    to::Vector{Int64}
    data::ComplexNums

    function Lines(lines::Vector{Int64},
                   from::Vector{Int64},
                   to::Vector{Int64},
                   i::Vector{T},
                   j::Vector{T};
                   type = "impedance") where T
        L = collect(1:length(lines))
        data = ComplexNums(i, j, type = type)
        new(
            L,
            lines,
            from,
            to,
            data      
        )
    end
end

mutable struct Input
    buses::Buses
    lines::Lines

    function Input(buses::Vector{Int64},
                   lines::Vector{Int64},
                   from::Vector{Int64},
                   to::Vector{Int64},
                   i_line::Vector{T},
                   j_line::Vector{T};
                   capacitor_bank::Vector = Any[],
                   type_line = "impedance") where T
        
        new(
            Buses(buses, capacitor_bank),
            Lines(lines, from, to, i_line,j_line, type = type_line)
        )
    end
end

mutable struct AlgorithmResults
    V::Matrix
    θ::Matrix
    Δ::Matrix
    gap::Vector
end

mutable struct Results
    P::Matrix
    Q::Matrix
    S::Matrix
    V::Vector
    θ::Vector
    algorithm_results::AlgorithmResults
end