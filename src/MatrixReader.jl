import Base: eof, close
import Statistics: mean
export MatrixWriter, close, writeMatrix, MatrixReader, value, eof, advance!, mean, cor, cov

mutable struct MatrixWriter
    fileStream
    zerocount::Int
    maxnum::Int ##Somewhat of unused number
    expsize::Int
    datatype::DataType
end

function MatrixWriter(fileName::String, expsize::Int64, datatype::DataType)
    f = open(fileName, "w")
    assert(datatype in [UInt8, UInt16])
    if datatype == UInt8 #UInt8 mode stores up to 100
        mw = MatrixWriter(f, 1, 100, expsize, datatype)
    elseif datatype == UInt16 #UInt16 mode stores up to 60000
        mw = MatrixWriter(f, 1, 60000, expsize, datatype)
    end
    writeHeader(mw)
    mw
end

function writeHeader(mw::MatrixWriter)
    #Write other information.
    write(mw.fileStream, UInt16(mw.expsize)) 
    write(mw.fileStream, UInt16(mw.maxnum))
    
    #Indicate datatype.
    if mw.datatype == UInt8
        write(mw.fileStream, UInt16(0)) 
    elseif mw.datatype == UInt16
        write(mw.fileStream, UInt16(1))
    end
end

close(mw::MatrixWriter) = close(mw.fileStream)

function writeMatrix(mw::MatrixWriter, matrix::Array{Int64,2})
    
    skipnum = typemax(mw.datatype)-mw.maxnum
    for j in 1:size(matrix)[1]
        for i in 1:size(matrix)[2]
            if matrix[j, i] == 0
                mw.zerocount += 1
            elseif matrix[j, i] != 0
                
                # Figure out how much entry you have skipped and write it out.
                if mw.zerocount != 1
                    temp = mw.zerocount-1
                    _multi = floor(Int, temp/skipnum)
                    _const = temp%skipnum
                    
                    for _ in 1:_multi
                        write(mw.fileStream, mw.datatype(typemax(mw.datatype)))
                    end
                    
                    if _const != 0
                        write(mw.fileStream, mw.datatype(mw.maxnum+_const))
                    end
                    
                end
                
                # Write non-zero number.
                if matrix[j, i] < mw.maxnum+1 
                    # Write actual number
                    write(mw.fileStream, mw.datatype(matrix[j, i]))
                else
                    println("$(matrix[j, i]) is too big. Writing $(mw.maxnum) instead")
                    write(mw.fileStream, mw.datatype(mw.maxnum))
                end
                
                # Start new zero count
                mw.zerocount = 1
            end
        end
    end
end

mutable struct MatrixReader
    # For basic matrix reading
    fileStream
    expsize::Int
    maxnum::Int
    binsize::Int
    datatype::DataType
    blocksize::Int
    data::Array{Int64, 2}
    offset::Int
    
    # For buffering bites
    bufferpointer::Int
    buffersize::Int
    buffer
end

function MatrixReader(fileName::String, blocksize; buffsize=10000000)
    f = open(fileName)
    
    # These are read as header of file.
    _expsize = Int(read!(f, Array{UInt16}(undef, 1))[1])
    _maxnum = Int(read!(f, Array{UInt16}(undef, 1))[1])
    _datatype = Int(read!(f, Array{UInt16}(undef, 1))[1])
    
    if _datatype == 0
        dt = UInt8
    elseif _datatype == 1
        dt = UInt16
    end
    
    _buffsize = buffsize # This is 10Mb of buffer size for default
    mr = MatrixReader(f, _expsize, _maxnum, _maxnum, dt, blocksize, zeros(Int64, (blocksize, _expsize)), 0, _buffsize+1, _buffsize, zeros(dt, _buffsize))
end

close(mr::MatrixReader) = close(mr.fileStream)
value(mr::MatrixReader) = mr.data
eof(mr::MatrixReader) = eof(mr.fileStream) && mr.bufferpointer > length(mr.buffer)

function advance!(mr::MatrixReader)
    temp = zeros(Int64, (mr.blocksize, mr.expsize))
    pointer = mr.offset
    while !eof(mr) && pointer < mr.blocksize*mr.expsize
        #Reload buffer
        if mr.bufferpointer > length(mr.buffer)
            if mr.datatype == UInt8
                mr.buffer = read(mr.fileStream, mr.buffersize)
            elseif mr.datatype == UInt16
                tempbuffer = read(mr.fileStream, mr.buffersize*2)
                assert(length(tempbuffer)%2==0)
                mr.buffer = reinterpret(UInt16, tempbuffer)
            end
            mr.bufferpointer = 1
        end
        
        # Read in from buffer
        value = mr.buffer[mr.bufferpointer]
        mr.bufferpointer += 1
        
        # Fill matrix
        if value < mr.maxnum+1
            temp[ceil(Int, (pointer+1)/mr.expsize), pointer%mr.expsize+1] = value
            pointer += 1
        else
            pointer += value-mr.maxnum
        end
    end
    mr.offset = pointer - mr.blocksize*mr.expsize
    mr.data = temp
end

function mean(mr::MatrixReader; verbose=0)
    total = zeros(1, mr.expsize)
    advance!(mr)
    
    count = 0
    while !eof(mr)
        total += sum(mr.data, 1)
        count += mr.blocksize
        advance!(mr)
        if verbose>0 && count % (verbose*mr.blocksize) == 0
            println(count)
        end
    end
    (total/count)[1, 1:end]
end

function cor(mr::MatrixReader, means; verbose=0)
    # loop through all the chunks
    XtX = zeros(Float64, mr.expsize, mr.expsize)
    advance!(mr)
    
    count = 0
    while !eof(mr)
        temp = mr.data
        _x = temp'.-means
        BLAS.syrk!('U', 'N', 1.0, _x, 1.0, XtX)
        advance!(mr)
        if verbose>0 && count % (verbose) == 0
            println(count, ":", size(XtX))
        end
        count += 1
    end
    
    # Converting top half covariance to correlation.
    uppercov = XtX
    uppercor = cov2cor!(uppercov)
    cor = uppercor+uppercor'-eye(size(uppercor)[1])
    cor
end

function cov(mr::MatrixReader, means; verbose=0)
    # loop through all the chunks
    XtX = zeros(Float64, mr.expsize, mr.expsize)
    advance!(mr)
    
    count = 0
    while !eof(mr)
        temp = mr.data
        _x = temp'.-means
        BLAS.syrk!('U', 'N', 1.0, _x, 1.0, XtX)
        advance!(mr)
        if verbose>0 && count % (verbose) == 0
            println(count, ":", size(XtX))
        end
        count += 1
    end
    
    cov = XtX+XtX'
    for i in 1:size(XtX)[1]
       cov[i, i] ./= 2 
    end
    cov
end
