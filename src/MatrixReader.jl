import Base: eof, close, mean
export MatrixWriter, close, writeMatrix, MatrixReader, value, eof, advance!, mean, cor, cov

type MatrixWriter
    fileStream
    zerocount::Int
    binsize::Int
    expsize::Int
end

function MatrixWriter(fileName::String, expsize)
    f = open(fileName, "w")
    mw = MatrixWriter(f, 1, 100, expsize)
    writeHeader(mw)
    mw
end

function writeHeader(mw::MatrixWriter)
   write(mw.fileStream, UInt16(mw.expsize)) 
   write(mw.fileStream, UInt16(mw.binsize)) 
end

close(mw::MatrixWriter) = close(mw.fileStream)

function writeMatrix(mw::MatrixWriter, matrix::Array{Int64,2})
    for j in 1:size(matrix)[1]
        for i in 1:size(matrix)[2]
            #println(matrix[j, i],":",zerocount)
            if matrix[j, i] == 0
                mw.zerocount += 1
            elseif matrix[j, i] != 0
                if mw.zerocount != 1
                    temp = mw.zerocount-1
                    _multi = floor(Int, temp/155)
                    _const = temp%155
                    
                    for _ in 1:_multi
                        write(mw.fileStream, UInt8(255))
                    end
                    
                    if _const != 0
                        write(mw.fileStream, UInt8(100+_const))
                    end
                    
                end 
                assert(matrix[j, i] < 101 )
                write(mw.fileStream, UInt8(matrix[j, i]))
                mw.zerocount = 1
            end
        end
    end
end

type MatrixReader
    # For basic matrix reading
    fileStream
    expsize::Int
    binsize::Int
    blocksize::Int
    data::Array{Int64, 2}
    offset::Int
    
    # For buffering bites
    bufferpointer::Int
    buffersize::Int
    buffer::Array{UInt8,1}
end

function MatrixReader(fileName::String, blocksize; buffsize=5000000)
    f = open(fileName)
    
    # These are read as header of file.
    _expsize = read(f, UInt16, 1)[1]
    _binsize = read(f, UInt16, 1)[1] 
    
    _buffsize = buffsize # This is 10Mb of buffer size for default
    mr = MatrixReader(f, _expsize, _binsize, blocksize, zeros(Int64, (blocksize, _expsize)), 0, _buffsize+1, _buffsize, zeros(UInt8, _buffsize))
    #advance!(mr)
    #mr
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
            mr.buffer = read(mr.fileStream, mr.buffersize)
            mr.bufferpointer = 1
        end
        
        # Read in from buffer
        value = mr.buffer[mr.bufferpointer]
        mr.bufferpointer += 1
        
        # Fill matrix
        if value < 101
            temp[ceil(Int, (pointer+1)/mr.expsize), pointer%mr.expsize+1] = value
            pointer += 1
        else
            pointer += value-100
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