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
    
    return YBUS
end