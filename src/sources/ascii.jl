const ASCII_X_ORDER = ForwardOrdered()
const ASCII_Y_ORDER = ReverseOrdered()

has_layers(::Type{ASCIIFile}) = false

defaultcrs(::Type{ASCIIFile}) = EPSG(4326)
defaultmappedcrs(::Type{ASCIIFile}) = EPSG(4326)
struct ASCIIparams{T,F}
    filename::F
    params::Dict{Symbol, Real}
    write::Bool
end

function _open(f, ::Type{ASCIIFile}, filename::AbstractString; write = false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    _open(f, ASCIIFile, ASCIIparams(filename; write); kw...)
end
_open(f, ::Type{ASCIIFile}, params::ASCIIparams; kw...) = f(params)

Array(params::ASCIIparams) = _read_ascii(filename(params); lazy = false)[1]

filename(p::ASCIIparams) = p.filename
params(p::ASCIIparams) = p.params

Base.size(p::ASCIIparams) = (params(p)[:nc], params(p)[:nr])
Base.eltype(p::ASCIIparams{T}) where T = T

function ASCIIparams(filename::AbstractString; write = false)
    meta = _read_ascii(filename; lazy = true)
    ASCIIparams{Float64, typeof(filename)}(filename, meta, write)
end

# DimensionalData methods
######################################################
function DD.dims(ascp::ASCIIparams, crs=nothing, mappedcrs=nothing)
    crs = crs isa Nothing ? ASCII_DEFAULT_CRS : crs

    nc, nr = size(ascp)

    pars = params(ascp)

    xbounds = (pars[:xll], pars[:xll] + pars[:dx] * nc)
    ybounds = (pars[:yll], pars[:yll] + pars[:dy] * nr)

    # Always intervals
    xspan = (xbounds[2] - xbounds[1]) / nc
    yspan = (ybounds[2] - ybounds[1]) / nr

    # Not fully implemented yet
    xy_metadata = Metadata{ASCIIFile}(Dict())

    xindex = LinRange(xbounds[1], xbounds[2] - xspan, nc)
    yindex = LinRange(ybounds[2] - yspan, ybounds[1], nr)

    xlookup = Projected(xindex;
        order=ASCII_X_ORDER,
        span=Regular(xspan),
        sampling=Intervals(Start()),
        metadata=xy_metadata,
        crs=crs,
        mappedcrs=mappedcrs,
        dim=X()
    )
    ylookup = Projected(yindex;
        order= ReverseOrdered(),
        span=Regular(yspan),
        sampling=Intervals(Start()),
        metadata=xy_metadata,
        crs=crs,
        mappedcrs=mappedcrs,
        dim=Y()
    )

    x = X(xlookup)
    y = Y(ylookup)

    return x,y

end

missingval(ascp::ASCIIparams) = params(ascp)[:NA]

DD.metadata(ascp::ASCIIparams) = Metadata{ASCIIFile}()

# Array
######################################################
function FileArray(ascp::ASCIIparams, filename = filename(ascp); kw...)
    size_ = size(ascp)
    eachchunk = DiskArrays.GridChunks(size_, size_)
    haschunks = DiskArrays.Unchunked()
    T = eltype(ascp)
    N = length(size_)
    FileArray{ASCIIFile, T, N}(filename, size_; eachchunk, haschunks, kw...)
end

# Base methods
######################################################
# maybe not needed for the moment idk

# AbstrackRasterStack methods
function Base.open(f::Function, A::FileArray{ASCIIFile}, key...; write = A.write)
    _open(ASCIIFile, filename(A); key=key(A), write, kw...) do dat
        f(RasterDiskArray{ASCIIFile}(dat, DA.eachchunk(A), DA.haschunks(A)))
    end
end

# Utils
######################################################

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

        params = Dict(:nr => nr, :nc => nc, :xll => xll, :yll => yll, :dx => dx, :dy => dy, :NA => NA)

        if !lazy
            # for use as Raster.data, the coordinates must be swapped
            # compared to a "normal" matrix
            out = Array{Float64}(undef, nc, nr)

            for row in 1:nr
                out[:, row] = parse.(Float64, split(readline(file), " ")[2:end]) # data lines start with a space
            end
            output = (out, params)
        else
            output = params
        end
    end

    return output
end