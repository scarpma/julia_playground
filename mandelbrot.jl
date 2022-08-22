function mandelbrot(a, n)
    z = 0
    for i=1:n
        z = z^2 + a
    end
    return z
end

spacing = 0.005
y=range(-1.,1.,step=spacing)
x=range(-2.,0.5,step=spacing)
M = Array{Int}(undef,(length(x),length(y)))

for i=1:length(y)
    for j=1:length(x)
        M[j,i] = abs(mandelbrot(complex(x[j], y[i]), 50)) < 2 ? 1 : 0
    end
end

using Plots
heatmap(M', color = :greys)
savefig("my_plot.pdf")
