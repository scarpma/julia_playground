using Base.Threads
@threads for i in 0:100
    println("thread id ", threadid()," gives ", i)
end
println("finished")

