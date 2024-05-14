#=
Author: Jerome Guterl (guterlj@fusion.gat.com)
Company: General Atomics
ADAS.jl (c) 2024
=#

function get_block_attr(block_header::String, attr::String; outputformat=Int64)
    rg = Regex("$attr\\s*=\\s*([\\d]+)")
    m = match(rg, block_header)
    if m === nothing
        return nothing
    else
        return parse.(outputformat, m.captures[1])
    end
end

function split_blocks(lines::Vector{String})
    blocks = Vector{ADASBlock}()
    iprevious = 1
    for i in 1:length(lines)
        if startswith(strip(lines[i]), "--") || startswith(strip(lines[i]), "C-")
            push!(blocks, ADASBlock(lines[max(1, iprevious - 1)], lines[iprevious:i-1]))
            iprevious = i + 1
        end

    end
    @debug "# of blocks = $(length(blocks))"
    return blocks
end

function get_comments(lines::Vector{String})
    comments = Vector{String}()
    for i in 1:length(lines)
        if startswith(lines[i], "C")
            push!(comments, lines[i])
        end
    end
    return comments
end

function get_units(comments::Vector{String})
    for c in comments
        m = match(Regex("IN UNITS OF"), c) # todo: capture directly units through regex
        if m !== nothing
            return string(split(comments[3], "IN UNITS OF")[end])
        end
    end
    return string("")
end