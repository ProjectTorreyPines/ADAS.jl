
function make_struct(structname::Symbol, fieldnames::Vector{Symbol}, ; k=2)
    types_list = get_types_list(k)
    @assert length(fieldnames) < length(values(types_list))
    for f in fieldnames
        @assert f ∉ values(types_list)
    end
    fields = [(f, types_list[i]) for (i, f) in enumerate(fieldnames)]
    header = expr_curly(structname, [f[2] for f in fields])
    blk = Expr(:block)
    for f in fields
        push!(blk.args, Expr(:(::), f[1], f[2]))
    end
    return Expr(:struct, false, header, blk)
end
params = Dict(:a => 1.0, :b => [1.0, 2.0], :c => "abc")
eval(make_struct(:MyStruct, collect(keys(params))))
mystruct = MyStruct(collect(values(params))...)
show(mystruct)
