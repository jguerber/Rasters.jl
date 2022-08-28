const ASCII_X_ORDER = ForwardOrdered()
const ASCII_Y_ORDER = ReverseOrdered()

has_layers(::Type{ASCIIfile}) = false

defaultcrs(::Type{ASCIIfile}) = EPSG(4326)
defaultmappedcrs(::Type{ASCIIfile}) = EPSG(4326)
struct ASCIIparams{T,F}
    filename::F
    params::Dict{Symbol, Real}
    write::Bool
end

function _open(f, ::Type{ASCIIfile}, filename::AbstractString; write = false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    _open(f, ASCIIfile, ASCIIparams(filename; write); kw...)
end
_open(f, ::Type{ASCIIfile}, params::ASCIIparams; kw...) = f(params)

Array(ascp::ASCIIparams) = _io_ascii(Array, ascp; write = ascp.write)

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
    xy_metadata = Metadata{ASCIIfile}(Dict())

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

DD.metadata(ascp::ASCIIparams) = Metadata{ASCIIfile}()
DD.ndims(ascp::ASCIIparams) = 2

# Array
######################################################
function FileArray(ascp::ASCIIparams, filename = filename(ascp); kw...)
    size_ = size(ascp)
    eachchunk = DiskArrays.GridChunks(size_, size_)
    haschunks = DiskArrays.Unchunked()
    T = eltype(ascp)
    N = length(size_)
    FileArray{ASCIIfile, T, N}(filename, size_; eachchunk, haschunks, kw...)
end

# Base methods
######################################################
function Base.write(filename::String, ::Type{ASCIIfile}, A::AbstractRaster{T,2}) where T
    _write_ascii(filename, A)
end

function Base.write(filename::String, ::Type{ASCIIfile}, A::AbstractRaster{T,3}) where T
    # TODO: this is really clumsy, improve behavior
    if hasdim(A, Band)
        _write_ascii(filename, A[Band(1)])
    else
        throw("Unsupported thrid dimension when attempting to write a .asc array. Please drop third dimension.")
    end
end

function Base.write(ascp::ASCIIparams, dat::AbstractArray{T,2}) where T
    pars = params(acsp)
    _write_ascii(filename(ascp), pars, dat)
end

# AbstrackRasterStack methods
function Base.open(f::Function, A::FileArray{ASCIIfile}, key...; write = A.write)
    _open(ASCIIfile, filename(A); write) do dat
        mappedA = f(Array(dat))
        RasterDiskArray{ASCIIfile}(mappedA, DA.eachchunk(A), DA.haschunks(A))
    end
end

# Utils
######################################################

function _io_ascii(f, ascp::ASCIIparams; write = false)
    A = f(_read_ascii(filename(ascp); lazy = false)[1])
    if write
        println(A[1,1])
        w = _write_ascii(filename(ascp), params(ascp), A)
    end
    return A
end

"""
    _read_ascii

Reads an ASCII file. Parameters are parsed according to the [AAIGrid](https://gdal.org/drivers/raster/aaigrid.html) format.
"""
function _read_ascii(filename::AbstractString; lazy = false)
    isfile(filename) || throw(ArgumentError("File $filename does not exist"))
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

            for row in 1:nr # columns of a (nc, nr) array
                out[:, row] = parse.(Float64, split(readline(file), " ")[2:end]) # data lines start with a space
            end
            output = (out, params)
        else
            output = params
        end
    end

    return output
end

"""
    _write_ascii
"""
function _write_ascii(filename::String, A::AbstractRaster)
    # Collect Parameters
    x, y = DD.dims(A)

    ncols, nrows = length(x), length(y)
    xll = min(bounds(x)...)
    yll = min(bounds(y)...)

    dx = abs((bounds(x)[2] - bounds(x)[1])/ncols)
    dy = abs((bounds(y)[2] - bounds(y)[1])/nrows)
    nodatavalue = missingval(A)
    if nodatavalue isa Nothing
        nodatavalue = -9999
    end
    # Write
    _write_ascii(filename, ncols, nrows, xll, yll, dx, dy, nodatavalue, Array(A))
end

function _write_ascii(filename, pars::Dict{Symbol, Real}, dat::AbstractArray)
    _write_ascii(filename, pars[:nc], pars[:nr], pars[:xll], pars[:yll], pars[:dx], pars[:dy], pars[:NA], dat)
end

function _write_ascii(filename, ncols, nrows, xll, yll, dx, dy, nodatavalue, dat)
    size(dat) == (nrows, ncols) || throw(ArgumentError("$nrows rows and $ncols cols incompatible with array of size $(size(dat))"))
    # Write
    open(filename, "w") do f
        write(f,
            """
            ncols        $(string(ncols))
            nrows        $(string(nrows))
            xllcorner    $(string(xll))
            yllcorner    $(string(yll))
            dx           $(string(dx))
            dy           $(string(dy))
            NODATA_value  $(string(nodatavalue))
            """
        )
        for col in 1:nrows # ascii format is column by column
            write(f, " " * join(dat[:, col], " ") * "\n")
        end
    end
    return filename
end