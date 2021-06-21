function capacitor_bank_in_ybus(input)
    @unpack buses = input
    @unpack B = buses
    n_bus = length(B)
    modify_ybus = Array{Any, 2}(undef, n_bus, n_bus)
    modify_ybus[:,:] .= 0.0
    if !isempty(buses.capacitor_bank)
        for b in 1:n_bus
            modify_ybus[b,b] = buses.capacitor_bank[b]
        end
    end
    return modify_ybus
end

function ybus(input::Input)
    @unpack lines, buses = input
    @unpack B = buses
    @unpack L, from, to, data = lines
    @unpack y = data

    n_bus = length(B)
    n_lines = length(L)
    YBUS = zeros(Number, n_bus, n_bus)
    for l in L
        YBUS[from[l], to[l]] += - y[l]
        YBUS[to[l], from[l]] += - y[l]
        YBUS[from[l], from[l]] += y[l]
        YBUS[to[l], to[l]] += y[l]
    end
    if !isempty(buses.capacitor_bank)
        YBUS = YBUS .+ capacitor_bank_in_ybus(input)
    end
    return YBUS
end