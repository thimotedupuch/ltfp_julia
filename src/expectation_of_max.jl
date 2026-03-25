using CairoMakie
using Statistics

# ---- First panel ----
n = 10_000
nrep = 1000

# definition of cumulative max function
cummax(x) = accumulate(max, x)
cummax(A; dims=1) = mapslices(x -> accumulate(max, x), A; dims=dims)



X = randn(n, nrep)

fig = Figure()

ax1 = Axis(fig[1, 1])
lines!(ax1, X[:, 1], color = :blue)              # sample
lines!(ax1, cummax(X[:, 1]), color = :red)       # cumulative max
ax1.xlabel = "n"
ax1.ylabel = "single draw"

# ---- Second panel ----
n = 10_000
nrep = 10

X = randn(n, nrep)

ax2 = Axis(fig[1, 2])
xvals = sqrt.(log.(1:n))

for j in 1:nrep
    lines!(ax2, xvals, cummax(X[:, j]), color = :blue)
end

ax2.xlabel = "(log n)^{1/2}"
ax2.ylabel = "Samples of max"

# ---- Third panel ----
ax3 = Axis(fig[1, 3])

cm = cummax(X, dims = 1)             # cumulative max along rows
cm = Array(cm)                       # ensure standard array

m = mean(cm, dims = 2)
s = std(cm, dims = 2)

lines!(ax3, xvals, vec(m), color = :blue)
lines!(ax3, xvals, vec(m .+ s), linestyle = :dot, color = :blue)
lines!(ax3, xvals, vec(m .- s), linestyle = :dot, color = :blue)

ax3.xlabel = "(log n)^{1/2}"
ax3.ylabel = "Expectation of max"

save("expectation_of_max.svg", fig)
