include("../src/MatrixReader.jl")

function fakedata(fakecount, blocksize, expsize; maxval=100)
    testdata = zeros(Int, (blocksize, expsize))
    for _ in 1:fakecount
        i = rand(1:blocksize)
        j = rand(1:expsize)
        v = rand(1:maxval)
        testdata[i,j] = v
    end
    testdata
end

# Testing with 10000 random data with random sparcity, dimensions.

println("For UInt8 version")

for itr in 1:10000
    expsize = rand(20:100)
    
    ######################
    # Write in fake data #
    ######################
    mw = MatrixWriter("test.data", expsize, UInt8)

    # data1
    t1 = fakedata(rand(1:200), 100, expsize, maxval=100)
    writeMatrix(mw, t1)
    # data2
    t2 = fakedata(rand(1:200), 100, expsize, maxval=100)
    writeMatrix(mw, t2)
    close(mw)
    
    #####################
    # Read in fake data #
    #####################
    mr = MatrixReader("test.data", 25, buffsize=rand(10:50))
    new = advance!(mr)
    m1 = mr.data
    new = advance!(mr)
    m2 = mr.data
    new = advance!(mr)
    m3 = mr.data
    new = advance!(mr)
    m4 = mr.data
    new = advance!(mr)
    m5 = mr.data
    new = advance!(mr)
    m6 = mr.data
    new = advance!(mr)
    m7 = mr.data
    new = advance!(mr)
    m8 = mr.data
    close(mr)

    #################
    # Test equality #
    #################
    x1 = vcat(t1, t2)
    new = vcat(m1, m2, m3, m4, m5, m6, m7, m8)
    @assert new == x1
    if itr%1000 == 0
        println(itr," test passed. ", size(x1), ":", sum(x1))
    end
end

println("For UInt16 version")

for itr in 1:10000
    expsize = rand(20:100)
    
    ######################
    # Write in fake data #
    ######################
    mw = MatrixWriter("test.data", expsize, UInt16)

    # data1
    t1 = fakedata(rand(1:200), 100, expsize, maxval=60000)
    writeMatrix(mw, t1)
    # data2
    t2 = fakedata(rand(1:200), 100, expsize, maxval=60000)
    writeMatrix(mw, t2)
    close(mw)
    
    #####################
    # Read in fake data #
    #####################
    mr = MatrixReader("test.data", 25, buffsize=rand(10:50))
    new = advance!(mr)
    m1 = mr.data
    new = advance!(mr)
    m2 = mr.data
    new = advance!(mr)
    m3 = mr.data
    new = advance!(mr)
    m4 = mr.data
    new = advance!(mr)
    m5 = mr.data
    new = advance!(mr)
    m6 = mr.data
    new = advance!(mr)
    m7 = mr.data
    new = advance!(mr)
    m8 = mr.data
    close(mr)

    #################
    # Test equality #
    #################
    x1 = vcat(t1, t2)
    new = vcat(m1, m2, m3, m4, m5, m6, m7, m8)
    @assert new == x1
    if itr%1000 == 0
        println(itr," test passed. ", size(x1), ":", sum(x1))
    end
end
