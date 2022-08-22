const ASCII_X_ORDER = ForwardOrdered()
const ASCII_Y_ORDER = ReverseOrdered()

has_layers(::Type{ASCIIFile}) = false

function _open(f, ::Type{ASCIIFile}, filename::AbstractString; write = false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    _open(f, ASCIIFile, ASCIIparams(filename; write))
end
_open(f, ::Type{ASCIIFile}, params::ASCIIparams; kw...) = f(params)

Array(params::ASCIIparams) = _read_ascii(filename(params); lazy = false)[1]

struct ASCIIparams{T,F}
    filename::F
    params::Dict{Symbol, Real}
    write::Bool
end

filename(p::ASCIIparams) = p.filename
params(p::ASCIIparams) = p.params

function ASCIIparams(filename::AbstractString; write = false)
    meta = _read_ascii(filename; lazy = true)
    ASCIIparams{Float64, typeof(filename)}(filename, meta, write)
end

# Utils
###################################################

"""
    _read_ascii

Reads an ASCII file. Parameters are parsed according to the [AAIGrid](https://gdal.org/drivers/raster/aaigrid.html) format.
"""
function _read_ascii(filename::AbstractString; lazy = false)

    output = open(filename, "r") do file
        nc = parse(Int, match(r"ncols (.+)", readline(file)).captures[1])
        nr = parse(Int, match(r"nrows (.+)", readline(file)).captures[1])
        xll = parse(Float64, match(r"xllcorner (.+)", readline(file)).captures[1])
        yll = parse(Float64, match(r"yllcorner (.+)", readline(file)).captures[1])
        dx = parse(Float64, match(r"dx (.+)", readline(file)).captures[1])
        dy = parse(Float64, match(r"dy (.+)", readline(file)).captures[1])
        NA = parse(Float64, match(r"NODATA_value (.+)", readline(file)).captures[1])

        params = Dict(:xll => xll, :yll => yll, :dx => dx, :dy => dy, :NA => NA)

        if !lazy
            out = Array{Float64}(undef, nr, nc)

            for row in nr:-1:1
                out[row, :] = parse.(Float64, split(readline(file), " ")[2:end]) # data lines start with a space
            end
            output = (out, params)
        else
            output = params
        end
    end

    return output
end